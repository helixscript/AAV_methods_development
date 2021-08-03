tmpFile <- function(){ paste0(paste0(stringi::stri_rand_strings(30, 1, '[A-Za-z0-9]'), collapse = ''), '.tmp') }

ppNum = function(n) format(n, big.mark = ",", scientific = FALSE, trim = TRUE)


blastReads <- function(reads, blastn, db, minEval){
  f <- tmpFile()
  writeFasta(reads, file = f)
  o <- tmpFile()
  system(paste0(blastn, ' -word_size 7 -evalue ', minEval, ' -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand" -query ', f, ' -db ', db, ' -out ', o))
  b <- data.frame()
  if(file.info(o)$size > 0){
    b <- read.table(o, sep = '\t', header = FALSE)
    names(b) <- c('qname', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'sstrand')
    b$alignmentLength <- b$qend - b$qstart + 1
  }     
  invisible(file.remove(f))
  invisible(file.remove(o))
  return(b)
}

blast2rearangements <- function(x){
  if(nrow(x) == 0) return(data.frame())
  x <- subset(x, alignmentLength >= 15 & evalue <= 1e-5 & gapopen == 0 & pident >= 98)
  if(nrow(x) == 0) return(data.frame())
  
  library(parallel)
  n <- 25
  cluster <- makeCluster(n)
  
  # Here we create a spliting variable across the blast data frame
  # being careful not to split read ids into different chunks.
  a <- floor(n_distinct(x$qname) / n)
  b <- 1
  
  message('step 1')
  library(data.table)
  z <- dplyr::select(x,  qname, qstart, qend, sstart, send, evalue) %>% group_split(qname)
  z <- as.data.table(bind_rows(mapply(function(x, n){ x$n <- n; x }, z, dplyr::ntile(1:length(z), n), SIMPLIFY = FALSE)))
  
  message('step 2')
  r <- bind_rows(parLapply(cluster, split(z, z$n), function(b){
  #r <- rbindlist(lapply(split(z, z$n), function(b){
    library(dplyr)
    library(IRanges)
    library(data.table)
    
    rbindlist(lapply(split(b, b$qname), function(b2){
   
      #### if(b2$qname %in% readIDs.list) browser()
         
      # Bin start and end positions to help mitigate sequencing and aligner error.
      # Binned positions are shifted to the lower end of 5 NT intervals.
      breaks <- seq(1, max(b2$qstart), by = 5)
      b2$qstart_binned <- as.integer(as.character(cut(b2$qstart, breaks = c(breaks, Inf), include.lowest = TRUE, labels = breaks)))
     
      breaks <- seq(1, max(b2$qend), by = 5)
      b2$qend_binned <- as.integer(as.character(cut(b2$qend, breaks = c(breaks, Inf), include.lowest = TRUE, labels = breaks)))
      
      # Sort BLAST results by query start position and evalue (low to high).
      b2 <- arrange(b2, qstart_binned, evalue)
      b2$strand <- ifelse(b2$send < b2$sstart, '-', '+')
      
      # Alignment to the negative starnd will result in the subject end to come before the start.
      # Switch it back so that they are sequential.
      b2$sstart2 <- ifelse(b2$sstart > b2$send, b2$send, b2$sstart)
      b2$send2   <- ifelse(b2$sstart > b2$send, b2$sstart, b2$send)
      
      # Shrink the ranges to prevent closely spaced ranges from being assembeled.
      # These positions have already been binned a bit to handle alignment chatter.
      b2$qstart_binned <- b2$qstart_binned + 3
      b2$qend_binned   <- b2$qend_binned - 3
      
      # Create IRanges
      ir <- IRanges(start = b2$qstart_binned, end = b2$qend_binned)
      if(length(ir) == 0) return(data.frame())
      
      # Name the ranges with the binned query positions followed by the actual subject positions.
      names(ir) <- paste0(b2$qstart, '..', b2$qend, '[', b2$sstart2, b2$strand, b2$send2, ']')
      
      o <- ir[1]
      invisible(lapply(split(ir, 1:length(ir)), function(a){
        if(all(! countOverlaps(o, a) > 0)){
          o <<- c(o, a)
        }
      }))
      
      if(length(o) == 0) return(data.frame())
      
      # Undo the range shrinkage.
      #n1 <- as.integer(stringr::str_extract(unlist(lapply(strsplit(names(o), '\\['), '[', 1)), '^\\d+')) - 3
      #n2 <- as.integer(stringr::str_extract(unlist(lapply(strsplit(names(o), '\\['), '[', 1)), '\\d+$')) + 3
      #names(o) <- paste0(n1, '..', n2, stringr::str_extract(names(o), '\\[.+\\]'))
      
      data.frame(readID = b2$qname[1], rearrangement = paste0(names(o), collapse = ';'))
    }))
  }))
  
  stopCluster(cluster)
  r
}

