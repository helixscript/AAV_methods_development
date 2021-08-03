library(ShortRead)
library(dplyr)
library(parallel)
library(stringr)
options(stringsAsFactors = FALSE)
nCPUs <- 25
requiredCharacterizedReadLengthPercentage <- 90
maxReadAlignmentStartPos <- 5
minFragWidth <- 100
minBlastEval <- 1e-5
plasmidLength <- 11976
blastn <- '/home/everett/ext/blast+/bin/blastn'
blastDB <- '/home/everett/projects/AAV_methods_development/dbs/pCCVC_TTRm8003'
source('lib.R')


# Read in PEAR merged paired-end reads which were filtered for fragments with widths >= minFragWidth
water <- readFastq('pearSeqs/water_control_S4.assembled.cutadapt.fastq')
plasmid <- readFastq('pearSeqs/8003_plasmid_S1.assembled.cutadapt.fastq')
vector <- readFastq('pearSeqs/8003_vector_noklenow_S3.assembled.cutadapt.fastq')
vectorKlenow <- readFastq('pearSeqs/8003_vector_plusklenow_S2.assembled.cutadapt.fastq')


# Remove fragments below the specified width, fragments may fall below the width defined 
# in the PEAR mergers after trimming with cutadapt.
water <- water[width(water) >= minFragWidth]
plasmid <- plasmid[width(plasmid) >= minFragWidth]
vector <- vector[width(vector) >= minFragWidth] 
vectorKlenow <- vectorKlenow[width(vectorKlenow) >= minFragWidth]


# Convert ShortRead objects to named DNAString objects.
water.ids <- sub('\\s+.+$', '', as.character(water@id))
water <- water@sread
names(water) <- paste0(water.ids, '|water')

plasmid.ids <- sub('\\s+.+$', '', as.character(plasmid@id))
plasmid <- plasmid@sread
names(plasmid) <- paste0(plasmid.ids, '|plasmid')

vector.ids <- sub('\\s+.+$', '', as.character(vector@id))
vector <- vector@sread
names(vector) <- paste0(vector.ids, '|vectorNoKlenow')

vectorKlenow.ids <- sub('\\s+.+$', '', as.character(vectorKlenow@id))
vectorKlenow <- vectorKlenow@sread
names(vectorKlenow) <- paste0(vectorKlenow.ids, '|vectorPlusKlenow')


# Remove duplicated reads.
water <- water[! duplicated(water)]
plasmid <- plasmid[! duplicated(plasmid)]
vector <- vector[! duplicated(vector)]
vectorKlenow <- vectorKlenow[! duplicated(vectorKlenow)]


# Combine reads now that their read ids can be used to separate them downstream.
reads <- Reduce('append', list(water, plasmid, vector, vectorKlenow))
readWidths <- tibble(readID = names(reads), width = width(reads))


# Blast reads against the vector sequence database.
cluster <- makeCluster(nCPUs)
clusterExport(cluster, c('blastReads', 'tmpFile', 'blastn', 'blastDB', 'minBlastEval'))

b <- bind_rows(parLapply(cluster, split(reads, ntile(1:length(reads), nCPUs)), function(f){
       library(ShortRead)
       blastReads(f, blastn, blastDB, minBlastEval)
     }))

saveRDS(b, file = 'b.rds')


r <- blast2rearangements(b)
saveRDS(r, file = 'r.rds')


# Determine the first and last alignment positions of the read models.
r$firstAlnPos <- unlist(lapply(r$rearrangement, function(x){
  as.integer(str_extract(x, '^\\d+'))
}))

r$lastAlnPos <- unlist(lapply(r$rearrangement, function(x){
                  as.integer(rev(stringr::str_match_all(x, '\\.\\.(\\d+)')[[1]][,2])[1])
                }))

message('Removing ', sprintf("%.2f%%", (sum(r$firstAlnPos > maxReadAlignmentStartPos) / nrow(r))*100), 
        ' of reads because their alignments start beyond position ', maxReadAlignmentStartPos)

r <- r[r$firstAlnPos <= maxReadAlignmentStartPos,]



# Add read widths to read data.frame.
r <- left_join(r, readWidths, by = 'readID')
r$coverage <- (((r$lastAlnPos - r$firstAlnPos)+1) / r$width)*100


# Remove reafs with incomplete characterization of their ends.
message(round((nrow(subset(r, coverage < requiredCharacterizedReadLengthPercentage)) / nrow(r))*100, digits = 2), 
        '% reads removed because of incomplete characterization')
r <- subset(r, coverage >= requiredCharacterizedReadLengthPercentage)


# Add additional metadata.
r$rearranged <- grepl(';', r$rearrangement)
r$sample <- unlist(lapply(strsplit(r$readID, '\\|'), '[', 2))


# Create a summary table of % recobmined reads.
readSummaryTable <- group_by(r, sample) %>%
  summarise(uniqueReads = ppNum(n()), percentRearranged = sprintf("%.2f%%", (sum(rearranged)/n())*100)) %>%
  ungroup() %>%
  data.frame()



r$rearrangementPos <- unlist(lapply(r$rearrangement, function(x){
  if(! grepl(';', x)) return(NA)
  #browser()
  o <- unlist(stringr::str_extract_all(x, '\\[[\\d+\\-\\+]+\\]'))
  a <- unlist(str_extract_all(o[1], '\\d+'))[ifelse(grepl('\\+', o[1]), 2, 1)]
  b <- NA
  c <- unlist(str_extract_all(o[length(o)], '\\d+'))[ifelse(grepl('\\+', o[length(o)]), 1, 2)]

  if(length(o) > 2){
    b <- unlist(lapply(o[2:(length(o)-1)], function(e){
      unlist(str_extract_all(e, '\\d+'))
    }))
    
  } else {
    b <- NA
  }
  
  d <- c(a, b, c)  
  d <- as.integer(unique(d[! is.na(d)]))
  d <- sapply(d, function(x) ifelse(x > plasmidLength, x - plasmidLength, x))
  paste0(d, collapse = ',')
}))


r2 <- subset(r, rearranged == TRUE)
p <- bind_rows(lapply(split(r2, r2$sample), function(x){
       data.frame(sample = x$sample[1], 
                  pos = unlist(lapply(x$rearrangementPos, function(x2){
                          as.integer(unlist(strsplit(x2, ',')))
                         })))
      }))
  

# bins = 240 selected for 50 NT bins
plots <- lapply(split(p, p$sample), function(x){
           d <- data.frame(x1=c(1, 146, 376, 498, 4880, 4926), 
                           x2=c(141, 370, 483, 4871, 4925, 5066), 
                           y1=rep(-10, 6), 
                           y2=rep(-100,6), 
                           f = c("5' ITR", "Promoter", "SynIntron", "Transgene", "Poly-A", "3' ITR"))
           d$f <- factor(d$f, levels = unique(d$f))

           ggplot(x, aes(pos)) + 
           theme_bw()+ 
           ggtitle(x$sample[1]) +
           geom_histogram(bins = 240) +
           scale_fill_manual(name = 'Vector component', values = c('red', 'green3', 'gray50', 'gold', 'dodgerblue1', 'red4')) +
           geom_rect(data=d, inherit.aes = FALSE, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2, fill=f), color="black") +
           labs(x = 'Plasmid position', y = 'Read junctures') + 
           theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 panel.background = element_blank(), axis.line = element_line(colour = "black"))
          })






p <- subset(r, sample == 'plasmid' & rearranged == TRUE)
z <- subset(r, sample == 'vectorPlusKlenow' & rearranged == TRUE)

openxlsx::write.xlsx(p, file = 'plasmid_reads.xlsx')
