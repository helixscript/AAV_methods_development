library(yaml)
library(ShortRead)
library(tidyverse)
library(parallel)
library(lubridate)
options(stringsAsFactors = FALSE)
blastn <- '/home/everett/ext/blast+/bin/blastn'
blastDB <- '/home/everett/projects/AAV_methods_development/dbs/pCCVC_TTRm8003'

# To run AAVengeR in RStudio, replace the following three lines with
# this line, edit to point to the intended configuration file.

config <- read_yaml('/home/everett/projects/AAV_methods_development/210715_M04734_0298_000000000-DCNLN/config.yml')

source(file.path(config$softwareDir, 'AAVengeR.lib.R'))

config$startTime <- ymd_hms(format(Sys.time(), "%y-%m-%d %H:%M:%S"))

# configFile <- commandArgs(trailingOnly = TRUE)
# if(! file.exists(configFile)) stop('Error -- configuration file not found.')
# config  <- read_yaml(configFile)


# Prepare run session.
if(dir.exists(config$outputDir)) stop('Error -- output directory already exists.')
dir.create(config$outputDir)

invisible(sapply(c('tmp', 'readIDs', 'readsRemoved', 'seqChunks', 'sampleReads', 'logs', 'fragReads'), 
                 function(x) dir.create(file.path(config$outputDir, x))))
write(capture.output(sessionInfo()), file = file.path(config$outputDir, 'sessionInfo.txt'))
config$logFile <- file.path(config$outputDir, 'logs', 'log')
write(date(), file = config$logFile)


# Read in the sample data and add default column values if not present.
samples <- read_delim(config$sampleConfigFile, delim  = '\t', col_names = TRUE, col_types = cols())


# Reverse compliment index1 sequences if requested.
if(config$indexReads.rc) samples$index1.seq <- as.character(reverseComplement(DNAStringSet(samples$index1.seq)))



if(! 'subject' %in% names(samples))   samples$subject <- 'subject'
if(! 'replicate' %in% names(samples)) samples$replicate <- 1
if(any(grepl('~|\\|', paste(samples$subject, samples$sample, samples$replicate)))) stop('Error -- tildas (~) are reserved characters and can not be used in the subject, sample, or replicate sample configuration columns.')


# Create unique sample identifiers -- subject~sample~replicate
samples$uniqueSample <- paste0(samples$subject, '~', samples$sample, '~', samples$replicate)
if(any(duplicated(samples$uniqueSample))) stop('Error -- all subject, sample, replicate id combinations are not unique.')


cluster <- makeCluster(config$demultiplexing.CPUs)
clusterExport(cluster, c('config', 'samples'))

# Quality trim virus reads and break reads.
invisible(parLapply(cluster, 
                    list(c(config$adriftReadsFile,  config$sequence.chunk.size, 'adriftReads',  file.path(config$outputDir, 'seqChunks')),
                         c(config$anchorReadsFile,  config$sequence.chunk.size, 'anchorReads',  file.path(config$outputDir, 'seqChunks'))), 
                    function(x){
                      library(ShortRead)
                      source(file.path(config$softwareDir, 'AAVengeR.lib.R'))
                      qualTrimReads(x[[1]], x[[2]], x[[3]], x[[4]])
                    }))


# Collate trimmed reads.
system(paste('cat', file.path(config$outputDir, 'seqChunks', 'anchor*'), ' > ', file.path(config$outputDir, 'trimmedAnchorReads.fastq')))
system(paste('cat', file.path(config$outputDir, 'seqChunks', 'adrift*'), ' > ', file.path(config$outputDir, 'trimmedAdriftReads.fastq')))
system(paste('rm',  file.path(config$outputDir, 'seqChunks', '*')))


# Convert reads to DNAStringSets and sync reads.
logMsg(config, 'Converting index reads to DNA strings.', config$logFile)
index1Reads <- readFastq(config$index1ReadsFile)
index1Reads <- Reduce('append', parLapply(cluster, split(index1Reads, ntile(1:length(index1Reads), config$demultiplexing.CPUs)), 
                                          function(x){source(file.path(config$softwareDir, 'AAVengeR.lib.R')); shortRead2DNAstringSet(x)}))


# Correct Golay encoded barcodes if requested.
if(config$correctGolayIndexReads){
  logMsg(config, 'Correcting golay bar code reads.', config$logFile)
  index1Reads <- Reduce('append', parLapply(cluster, split(index1Reads, ntile(1:length(index1Reads), config$demultiplexing.CPUs)), golayCorrection))
}


# Convert quality trimmed reads to xstring DNA objects.
logMsg(config, 'Converting anchor reads to DNA strings.', config$logFile)
anchorReads <- readFastq(file.path(config$outputDir, 'trimmedAnchorReads.fastq'))
anchorReads <- Reduce('append', parLapply(cluster, split(anchorReads, ntile(1:length(anchorReads), config$demultiplexing.CPUs)), 
                                          function(x){source(file.path(config$softwareDir, 'AAVengeR.lib.R')); shortRead2DNAstringSet(x)}))

logMsg(config, 'Converting break reads to DNA strings.', config$logFile)
adriftReads <- readFastq(file.path(config$outputDir, 'trimmedAdriftReads.fastq'))
adriftReads <- Reduce('append', parLapply(cluster, split(adriftReads, ntile(1:length(adriftReads), config$demultiplexing.CPUs)), 
                                          function(x){source(file.path(config$softwareDir, 'AAVengeR.lib.R')); shortRead2DNAstringSet(x)}))

invisible(file.remove(file.path(config$outputDir, 'trimmedAnchorReads.fastq')))
invisible(file.remove(file.path(config$outputDir, 'trimmedAdriftReads.fastq')))


# Synchronize reads.
logMsg(config, 'Synchronizing trimmed reads.', config$logFile)
reads <- syncReads(index1Reads, anchorReads, adriftReads)
index1Reads <- reads[[1]];  anchorReads  <- reads[[2]];  adriftReads  <- reads[[3]]
rm(reads)
gc()


# Split the trimmed reads into chunks for parallel processing.
logMsg(config, 'Distributing read data into chunks ...', config$logFile)
chunkNum <- 1
d <- tibble(i = ntile(1:length(index1Reads), config$demultiplexing.CPUs), n = 1:length(index1Reads))
invisible(lapply(split(d, d$i), function(x){
  index1Reads <- index1Reads[min(x$n):max(x$n)]
  anchorReads  <- anchorReads[min(x$n):max(x$n)]
  adriftReads  <- adriftReads[min(x$n):max(x$n)]
  save(index1Reads, anchorReads, adriftReads, file = file.path(config$outputDir, 'seqChunks', chunkNum))
  chunkNum <<- chunkNum + 1
}))
logMsg(config, paste0(chunkNum-1, ' data chunks created.'), config$logFile)

# Clean up and free up memory. 
rm(d, chunkNum, index1Reads, anchorReads, adriftReads)
gc()


save(list = ls(all=TRUE), file = file.path(config$outputDir, 'savePoint1.RData'))
logMsg(config, 'Starting sample chunk threads, resetting timer.', config$logFile)
config$startTime <- ymd_hms(format(Sys.time(), "%y-%m-%d %H:%M:%S"))

if(! dir.exists(file.path(config$outputDir, 'logs', 'cutadapt'))) dir.create(file.path(config$outputDir, 'logs', 'cutadapt'))
if(! dir.exists(file.path(config$outputDir, 'tmp', 'cutadapt'))) dir.create(file.path(config$outputDir, 'tmp', 'cutadapt'))


invisible(parLapply(cluster, list.files(file.path(config$outputDir, 'seqChunks'), full.names = TRUE), function(f){
#invisible(lapply(list.files(file.path(config$outputDir, 'seqChunks'), full.names = TRUE), function(f){
  library(ShortRead)
  library(tidyverse)
  source(file.path(config$softwareDir, 'AAVengeR.lib.R'))
  
  load(f)
  
  # Capture the chunk identifier.
  chunk.n <- unlist(str_match_all(f, '(\\d+)$'))[2]
  logFile <- file.path(config$outputDir, 'logs', paste0('seqChunk_', chunk.n, '.log'))
  
  
  # Loop through samples in sample data file to demultiplex and apply read specific filters.
  invisible(lapply(1:nrow(samples), function(r){
    r <- samples[r,]
    
    # Create barcode demultiplexing vectors.
    v1 <- vcountPattern(r$index1.seq, index1Reads, max.mismatch = config$index1Reads.maxMismatch) > 0
    logMsg(config, paste0('Chunk ', chunk.n, ': (', r$uniqueSample, ') ', sum(v1), ' reads pass barcode filter.'), logFile)
    
    log.report <- tibble(sample = r$uniqueSample, demultiplexedIndex1Reads = sum(v1))
    
    # Create break read linker barcode demultiplexing vector.
    v2 <- rep(TRUE, length(adriftReads))
    if('adriftReads.linkerBarcode.maxMismatch' %in% names(config)){
      ### browser()
      testSeq <- substr(r$adriftRead.linker.seq, r$adriftRead.linkerBarcode.start, r$adriftRead.linkerBarcode.end)
      v2 <- vcountPattern(testSeq, subseq(adriftReads, r$adriftRead.linkerBarcode.start, r$adriftRead.linkerBarcode.end), max.mismatch = config$adriftReads.linkerBarcode.maxMismatch) > 0
      logMsg(config, paste0('Chunk ', chunk.n, ': (', r$uniqueSample, ') ', sum(v2), ' reads pass linker code filter.'), logFile)
      log.report$demultiplexedLinkerReads <- sum(v2)
    } else {
      log.report$demultiplexedLinkerReads <- NA
    }
    
    
    # Test to see if any reads demultiplex to this row of the sample table and then subset reads to this sample.
    i <- base::intersect(which(v1), which(v2))
    if(length(i) == 0){
      logMsg(config, paste0('Chunk ', chunk.n, ': (', r$uniqueSample, ') No reads demultiplexed.'), logFile)
      log.report$demultiplexedReads <- 0
      write.table(log.report, sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE, file = file.path(config$outputDir, 'tmp', paste0(r$uniqueSample, '.', chunk.n, '.logReport')))
      return()
    } else {
      reads <- syncReads(index1Reads[i], anchorReads[i], adriftReads[i])
      index1Reads <- reads[[1]];  anchorReads  <- reads[[2]]; adriftReads  <- reads[[3]]
      log.report$demultiplexedReads <- length(index1Reads)
    }
    
    
    # Test the start of anchor reads (static or blast search options only)
    # Not compatible with anchorReads.captureLTRseq.method = 'lentiViralHMM'
    v1 <- rep(TRUE, length(anchorReads))
    if('anchorReads.startTest.maxMismatch' %in% names(config) & config$anchorReads.captureLTRseq.method != 'lentiViralHMM'){
      v1 <- Reduce('|', 
                   lapply(unlist(strsplit(r$anchorRead.identification, ',')), function(x){ 
                     testSeq <- substr(unlist(strsplit(x, ':'))[2], 1,  config$anchorReads.startTest.length)
                     vcountPattern(testSeq, 
                                   subseq(anchorReads, 1, config$anchorReads.startTest.length), 
                                   max.mismatch = config$anchorReads.startTest.maxMismatch) > 0
                   }))
      
      
      logMsg(config, paste0('Chunk ', chunk.n, ': (', r$uniqueSample, ') ', sum(!v1), ' reads removed by virus read start test filter (', config$anchorReads.startTest.length, ' NTs).'), logFile)
      log.report$readsPassingAnchorStartTest <- sum(v1)
    } else {
      log.report$readsPassingAnchorStartTest <- NA
    }
    
    if(! any(v1)){
      logMsg(config, paste0('Chunk ', chunk.n, ': (', r$uniqueSample, ') ', sum(!v1), ' All reads removed by virus read start test filter.'), logFile)
      write.table(log.report, sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE, file = file.path(config$outputDir, 'tmp', paste0(r$uniqueSample, '.', chunk.n, '.logReport')))
      return()
    } else {
      anchorReads <- anchorReads[v1]
    }
    
    
    # Create the anchor read over-read sequences from the last 10 NT of linker sequences.
    anchorReadOverReadSeq <- substr(r$adriftRead.linker.seq, nchar(r$adriftRead.linker.seq) - 10, nchar(r$adriftRead.linker.seq))
    anchorReadOverReadSeq <- Biostrings::DNAString(anchorReadOverReadSeq)
    anchorReadOverReadSeq <- as.character(Biostrings::reverseComplement(anchorReadOverReadSeq))
    
    logMsg(config, paste0('Chunk ', chunk.n, ': (', r$uniqueSample, ') ', ' anchor read over read sequence determined to be ', 
                          anchorReadOverReadSeq, ' from parsing provided linker sequence.'), logFile)
    
    # Add sample name and set scale...
    # createReadLengthPlots(anchorReads, adriftReads, paste0('ReadLengths_chunk_preAdapterTrim.', r$uniqueSample, '.', chunk.n, '.pdf'))
    
    
    # Trim virus read over-read sequences by looking for the reverse complement of end of the common linker sequence.
    anchorReads <- trimOverReadSeq(anchorReads, anchorReadOverReadSeq, logFile = paste0('seqChunk_', chunk.n, '_', r$uniqueSample, '.anchorReads.trimadapt.log'))
    
    
    # Determine the length of the study sequence and then select reads >= the study length then truncate 
    # all reads to the study length.
    
    n <- max(nchar(unlist(lapply(strsplit(unlist(strsplit(r$anchorRead.identification, ',')), ':'), '[[', 2))))
    anchorReads <- anchorReads[width(anchorReads) >= n]
    
    if(length(anchorReads) == 0){
      logMsg(config, paste0('Chunk ', chunk.n, ': (', r$uniqueSample, ') No reads remaining after over-read trimming and length check.'), logFile)
      write.table(log.report, sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE, file = file.path(config$outputDir, 'tmp', paste0(r$uniqueSample, '.', chunk.n, '.logReport')))
      return()
    }
    
    anchorReads <- subseq(anchorReads, 1, n)
    
    # Write out final reads. Add sample names to read IDs.
    names(adriftReads) <- paste0(names(adriftReads), '|', r$uniqueSample)
    names(anchorReads) <- paste0(names(anchorReads), '|', r$uniqueSample)
    
    writeFasta(adriftReads,  file = file.path(config$outputDir, 'tmp', paste0(r$uniqueSample, '.adriftReads.', chunk.n, '.fasta')))
    writeFasta(anchorReads,  file = file.path(config$outputDir, 'tmp', paste0(r$uniqueSample, '.anchorReads.', chunk.n, '.fasta')))
    
    write.table(log.report, sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE, file = file.path(config$outputDir, 'tmp', paste0(r$uniqueSample, '.', chunk.n, '.logReport')))
    
    logMsg(config, paste0('Chunk ', chunk.n, ': (', r$uniqueSample, ') completed with ', length(anchorReads), ' reads.'), logFile)
  }))
  
  logMsg(config, paste0('Read data chunk ', chunk.n, ' completed.'), file.path(config$outputDir, 'logs', 'log'))
}))

stopCluster(cluster)


# Collect all the logs from the different computational nodes and create a single report.
logReport <- bind_rows(lapply(list.files(file.path(config$outputDir, 'tmp'), pattern = '*.logReport$', full.names = TRUE), function(f){
  read.table(f, header = TRUE, sep = '\t')
}))


logReport <- bind_rows(lapply(split(logReport, logReport$sample), function(x){
  o <- data.frame(lapply(2:length(x), function(y){
    if(all(is.na(x[,y]))){
      return(NA)
    } else {
      return(sum(x[,y], na.rm = TRUE))
    }
  }))
  
  names(o) <- names(x)[2:length(x)]
  bind_cols(data.frame(sample = x[1,1]), o)
}))



# Collate demultiplexed reads using file name snibets from tmp directory to sampleReads directory.
collateSampleReads('anchorReads')
collateSampleReads('adriftReads')




#--------------------



blastReads <- function(f, db, analysis.width = 200){
  o <- tmpFile()
  system(paste0(blastn, ' -num_threads 30 -word_size 7 -evalue 50 -outfmt 6 -query ', f, ' -db ', db, ' -out ', o))
  b <- data.frame()
  if(file.info(o)$size > 0){
    b <- read.table(o, sep = '\t', header = FALSE)
    names(b) <- c('qname', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore')
    b$alignmentLength <- b$qend - b$qstart + 1
    b$pcoverage <- (b$alignmentLength / analysis.width) * 100
  }     
  invisible(file.remove(o))
  return(b)
}

blast2rearangements <- function(x){
  if(nrow(x) == 0) return(data.frame())
  x <- subset(x, alignmentLength >= 15 & evalue <= 1e-10)
  if(nrow(x) == 0) return(data.frame())
  
  library(parallel)
  n <- 40
  cluster <- makeCluster(n)
  
  # Here we create a spliting variable across the blast data frame
  # being careful not to split read ids into different chunks.
  a <- floor(n_distinct(x$qname) / n)
  b <- 1
  
  o <- bind_rows(lapply(split(x, x$qname), function(x2){
    x2$n <- b
    b <<- b + 1
    if(b > a) b <<- 1
    x2
  }))
  
  r <- bind_rows(parLapply(cluster, split(o, o$n), function(b){
    library(dplyr)
    library(IRanges)
    
    bind_rows(lapply(split(b, b$qname), function(b2){
      b2 <- arrange(b2, qstart)
      b2$strand <- ifelse(b2$send < b2$sstart, '-', '+')
      
      b2$sstart2 <- ifelse(b2$sstart > b2$send, b2$send, b2$sstart)
      b2$send2   <- ifelse(b2$sstart > b2$send, b2$sstart, b2$send)
      
      # Shrink the ranges to prevent closely spaced ranges from being assembeled.
      b2$qstart <- b2$qstart + 3
      b2$qend   <- b2$qend - 3
      
      ir <- IRanges(start = b2$qstart, end = b2$qend)
      if(length(ir) == 0) return(data.frame())
      
      names(ir) <- paste0(b2$qstart, '..', b2$qend, '[', b2$sstart2, b2$strand, b2$send2, ']')
      
      o <- ir[1]
      invisible(lapply(split(ir, 1:length(ir)), function(a){
        if(all(! countOverlaps(o, a) > 0)){
          o <<- c(o, a)
        }
      }))
      
      if(length(o) == 0) return(data.frame())
      
      # Undo the range shrinkage.
      n1 <- as.integer(stringr::str_extract(unlist(lapply(strsplit(names(o), '\\['), '[', 1)), '^\\d+')) - 3
      n2 <- as.integer(stringr::str_extract(unlist(lapply(strsplit(names(o), '\\['), '[', 1)), '\\d+$')) + 3
      names(o) <- paste0(n1, '..', n2, stringr::str_extract(names(o), '\\[.+\\]'))
      
      data.frame(readID = b2$qname[1], rearrangement = paste0(names(o), collapse = ';'))
    }))
  }))
  
  stopCluster(cluster)
  r
}


o <- bind_rows(lapply(list.files(file.path(config$outputDir, 'sampleReads'), pattern='anchor', full.names = TRUE), function(f){
  message(f, '\b')
  blast2rearangements(blastReads(f, blastDB))
}))


o$a <- unlist(lapply(o$rearrangement, function(x){
         n <- as.integer(unlist(str_extract_all(x, '\\d+')))
         ifelse(n[2] - n[1] >= 185, 0, 1)
       }))

# Count the number of vector recombinations
o$b1 <- unlist(lapply(o$rearrangement, function(x){
         length(unlist(strsplit(x, ';'))) - 1
       }))

# 0 if no vector recombinations, 1 if 1 or more vector recombinations
o$b2 <- unlist(lapply(o$rearrangement, function(x){
         ifelse(length(unlist(strsplit(x, ';'))) == 1, 0, 1) 
         }))

o$c <- unlist(lapply(o$readID, function(x){
        sub('~\\d+$', '', unlist(strsplit(x, '\\|'))[2])
}))

group_by(o, c) %>%
  summarise(reads = n(), percentRecomb = sprintf('%.2f%%', (sum(b2) / n())*100), percentOdd = sprintf('%.2f%%', (sum(a) / n())*100)) %>%
  ungroup()


# Reads which have a truncated piece of vector and then something else (non-vector)
subset(o, a == 1 & b == 1)

R2 <- readFastq(config$anchorReadsFile)
ids <- sub('\\s+.+$', '', as.character(R2@id))
R2 <- R2@sread
names(R2) <- ids

# subset(o, a == 1 & b1 == 0) %>% select(c, readID, rearrangement)
# R2[names(R2) == 'M04734:298:000000000-DCNLN:1:1102:25779:7162']




