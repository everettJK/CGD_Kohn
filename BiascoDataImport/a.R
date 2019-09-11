library(ShortRead)
library(stringr)
library(dplyr)
options(stringsAsFactors = FALSE)

PR0119           <- read.table('PR0119.tsv', sep = '\t')
PR0119$cellType  <- 'Nutrophils'
PR0119$timePoint <- unname(unlist(lapply(sapply(PR0119$V1, strsplit, ' '), '[[', 5)))
PR0119$patient   <- unname(unlist(lapply(sapply(PR0119$V1, strsplit, ' '), '[[', 1)))
PR0119$enz       <- unname(unlist(lapply(sapply(PR0119$V1, strsplit, ' '), '[[', 4)))
PR0119$runID     <- 'PR0119-83522440'
PR0119           <- select(PR0119, patient, cellType, enz, timePoint, runID, V2, V3)


PR0227           <- read.table('PR0227.tsv', sep = '\t')
PR0227$cellType  <- unname(unlist(lapply(sapply(PR0227$V1, strsplit, ' '), '[[', 2)))
PR0227$enz       <- unname(unlist(lapply(sapply(PR0227$V1, strsplit, ' '), '[[', 3)))
PR0227$timePoint <- unname(unlist(lapply(sapply(PR0227$V1, strsplit, ' '), '[[', 4)))
PR0227$patient   <- unname(unlist(lapply(sapply(PR0227$V1, strsplit, ' '), '[[', 1)))
PR0227$runID     <- 'PR0227-133032946'
PR0227           <- select(PR0227, patient, cellType, enz, timePoint, runID, V2, V3)
PR0227$cellType  <- sub('T', 'T-Cells', PR0227$cellType)
PR0227$cellType  <- sub('B', 'B-Cells', PR0227$cellType)
PR0227$cellType  <- sub('NK', 'NK Cells', PR0227$cellType)
PR0227$cellType  <- sub('Mono', 'Monocytes', PR0227$cellType)

LTR.barcodes <- read.table('LTR.barcodes', sep = '\t')
names(LTR.barcodes) <- c('id', 'seq')

breakPoint.barcodes <- read.table('breakPoint.barcodes', sep = '\t')
names(breakPoint.barcodes) <- c('id', 'seq')

breakPoint.barcodes$id <- gsub('\\-LC', '.LC', breakPoint.barcodes$id)

d <- bind_rows(PR0119, PR0227)
d$bcSeq <- paste0(LTR.barcodes[match(d$V2, LTR.barcodes$id),]$seq, breakPoint.barcodes[match(d$V3, breakPoint.barcodes$id),]$seq)
d$timePoint <- paste0('m', stringr::str_extract(d$timePoint, '\\d+'))
d$daysPostTrtmnt <- as.integer(stringr::str_extract(d$timePoint, '\\d+'))*30

samplesToReg <- d[! duplicated(paste(d$patient, d$timePoint, d$cellType)),]
samplesToReg$GTSP <- paste0('GTSP', 3101:(3101+nrow(samplesToReg)-1))

samplesToReg <- rowwise(samplesToReg) %>%
  mutate(sql = sprintf("insert into gtsp (SpecimenAccNum, SamplePatientCode, Patient, CellType, SampleCellType, Trial, DaysPostTrtmnt, Timepoint, PI, `PI Email`, Gender, refGenome) values ('%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s');",
                       GTSP, patient, patient, cellType, cellType, 'Biasco early CGD', daysPostTrtmnt, timePoint, 'Luca Biasco', 'l.biasco@ucl.ac.uk', 'm', 'hg38')) %>%
  ungroup()

write(samplesToReg$sql, file = 'sampleCreation.sql')

d$GTSP <- samplesToReg[match(paste(d$patient, d$cellType, d$timePoint), paste(samplesToReg$patient, samplesToReg$cellType, samplesToReg$timePoint)),]$GTSP
d$alias <- paste0(d$GTSP, '-', d$enz)

d$linkerSequence <- 'GTAAGGNNNNNNAGATCTGGAATGAACTGGCC'
d$gender         <- 'm'
d$primer         <- 'GAAAATC'
d$ltrBit         <- 'TCTAGCA'
d$largeLTRFrag   <- 'TGCTAGAGATTTTCCACACTGACTAAAAGGGT'
d$vectorSeq      <- 'vector_Biasco.fa'
d$refGenome      <- 'hg38'


# Manually check that each sample has three digests all with the sa
group_by(d, runID, patient, cellType, timePoint) %>%
summarise(cuts = paste0(enz, collapse = ', '),
          barcodes = paste0(unique(bcSeq), collapse = ', ')) %>%
ungroup()


filter(d, runID == 'PR0119-83522440') %>% 
group_by(patient, cellType, timePoint) %>%
  slice(1) %>%
  mutate(alias = paste0(GTSP, '-1')) %>%
  ungroup() %>%
select(alias, linkerSequence, bcSeq, gender, primer, ltrBit, largeLTRFrag, vectorSeq, refGenome) %>%
write.csv(row.names = FALSE, quote = FALSE, file = 'data-PR0119-83522440-sampleInfo.csv')

filter(d, runID == 'PR0227-133032946') %>% 
  group_by(patient, cellType, timePoint) %>%
  slice(1) %>%
  mutate(alias = paste0(GTSP, '-1')) %>%
  ungroup() %>%
  select(alias, linkerSequence, bcSeq, gender, primer, ltrBit, largeLTRFrag, vectorSeq, refGenome) %>%
  write.csv(row.names = FALSE, quote = FALSE, file = 'data-PR0227-133032946-sampleInfo.csv')

system('cp *.csv ~/task.sequencingReg.sampleData/sampleInfoFiles/')


#--------------------------------------------------------------------------------------------------

R1 <- readFastq('PR0119-83522440/Undetermined_S0_L001_R1_001.fastq.org.gz')
R2 <- readFastq('PR0119-83522440/Undetermined_S0_L001_R2_001.fastq.org.gz')

R1.bcs <- substr(sread(R1), 13, 20)
R2.bcs <- substr(sread(R2), 13, 20)
bcs <- paste0(R1.bcs, R2.bcs)
readIDs <- as.character(ShortRead::id(R1))
sprintf("%.2f%%", (sum(R1.bcs %in% LTR.barcodes$seq) / length(R1.bcs))*100)
sprintf("%.2f%%", (sum(R2.bcs %in% breakPoint.barcodes$seq) / length(R2.bcs))*100)

bcOutput <- paste0('@', readIDs, '\n', bcs, '\n+\n', 'IIIIIIIIIIIIIIII')
if(file.exists('PR0119-83522440/Undetermined_S0_L001_I1_001.fastq.gz')) file.remove('PR0119-83522440/Undetermined_S0_L001_I1_001.fastq.gz')
write(bcOutput, file = 'PR0119-83522440/Undetermined_S0_L001_I1_001.fastq')
system('gzip PR0119-83522440/Undetermined_S0_L001_I1_001.fastq')


# Narrow the reads to remove diversity and barcode NTs and write out the sequences.
# Swap R1 and R2 labels to match expected INSPIIRED input.
if(file.exists('PR0119-83522440/Undetermined_S0_L001_R2_001.fastq.gz')) file.remove('PR0119-83522440/Undetermined_S0_L001_R2_001.fastq.gz')
writeFastq(narrow(R1, 39, width(R1)), compress = TRUE, 
           file = 'PR0119-83522440/Undetermined_S0_L001_R2_001.fastq.gz')

if(file.exists('PR0119-83522440/Undetermined_S0_L001_R1_001.fastq.gz')) file.remove('PR0119-83522440/Undetermined_S0_L001_R1_001.fastq.gz')
writeFastq(narrow(R2, 56, width(R2)), compress = TRUE, 
           file = 'PR0119-83522440/Undetermined_S0_L001_R1_001.fastq.gz')
                        
system('scp PR0119-83522440/Undetermined_S0_L001_I1_001.fastq.gz PR0119-83522440/Undetermined_S0_L001_R1_001.fastq.gz PR0119-83522440/Undetermined_S0_L001_R2_001.fastq.gz microb120:/media/sequencing/Illumina/PR0119-83522440/Data/Intensities/BaseCalls/')


#                                               GAAAATC
#                                                      TCTAGCA
# R1: LTR read {12 div}{8 bc}{ACCCTTTTAGTCAGTGTGGAAAATCTCTAGCA}
# R1: select NTs 39 (12 + 8 + 18) + 1 ->
#                                                                        GTAAGGNNNNNNAGATCTGGAATGAACTGGCC
# R2: break point read {12 div}{8 bc}{GATCTGAATTCAGTGGCACAGCAGTTAGGNNNNNNGTAAGGNNNNNNAGATCTGGAATGAACTGGCC}
# R2 select NTs 56 (12 + 8 + 35) + 1 ->



# R1.trimmed <- trimTailw(R1, 2, '?', 5)
# R2.trimmed <- trimTailw(R2, 2, '?', 5)
# R1.trimmed <- R1.trimmed[width(R1.trimmed) >= 100]
# R2.trimmed <- R2.trimmed[width(R2.trimmed) >= 100]
# R1.trimmed <- subseq(sread(R1.trimmed), 1, 100)
# R2.trimmed <- subseq(sread(R2.trimmed), 1, 100)
# 
# m1 <- vmatchPattern('GATCTGAATTCAGTGGCACAGCAGTTAGGNNNNNNGTAAGGNNNNNNAGATCTGGAATGAACTGGCC',
#                     R2.trimmed, fixed = FALSE, max.mismatch = 3)
# 
# m2 <- vmatchPattern('GATCTGAATTCAGTGGCACAGCAGTTAGG',
#                     R2.trimmed, fixed = FALSE, max.mismatch = 2)
# 
# m1.s <- unlist(lapply(startIndex(m1), length))
# m2.s <- unlist(lapply(startIndex(m2), length))
# 
# m1.s <- ifelse(m1.s > 0, 1, 0)
# m2.s <- ifelse(m2.s > 0, 1, 0)
# 
# sum(m1.s)/length(m1.s)
# sum(m2.s)/length(m2.s)


#--------------------------------------------------------------------------------------------------


R1 <- readFastq('PR0227-133032946/Undetermined_S0_L001_R1_001.fastq.org.gz')
R2 <- readFastq('PR0227-133032946/Undetermined_S0_L001_R2_001.fastq.org.gz')

R1.bcs <- substr(sread(R1), 13, 20)
R2.bcs <- substr(sread(R2), 13, 20)
bcs <- paste0(R1.bcs, R2.bcs)
readIDs <- as.character(ShortRead::id(R1))
sprintf("%.2f%%", (sum(R1.bcs %in% LTR.barcodes$seq) / length(R1.bcs))*100)
sprintf("%.2f%%", (sum(R2.bcs %in% breakPoint.barcodes$seq) / length(R2.bcs))*100)

bcOutput <- paste0('@', readIDs, '\n', bcs, '\n+\n', 'IIIIIIIIIIIIIIII')

if(file.exists('PR0227-133032946/Undetermined_S0_L001_I1_001.fastq.gz')) file.remove('PR0227-133032946/Undetermined_S0_L001_I1_001.fastq.gz')
write(bcOutput, file = 'PR0227-133032946/Undetermined_S0_L001_I1_001.fastq')
system('gzip PR0227-133032946/Undetermined_S0_L001_I1_001.fastq')

if(file.exists('PR0227-133032946/Undetermined_S0_L001_R2_001.fastq.gz')) file.remove('PR0227-133032946/Undetermined_S0_L001_R2_001.fastq.gz')
writeFastq(narrow(R1, 39, width(R1)), compress = TRUE, 
           file = 'PR0227-133032946/Undetermined_S0_L001_R2_001.fastq.gz')

if(file.exists('PR0227-133032946/Undetermined_S0_L001_R1_001.fastq.gz')) file.remove('PR0227-133032946/Undetermined_S0_L001_R1_001.fastq.gz')
writeFastq(narrow(R2, 56, width(R2)), compress = TRUE, 
           file = 'PR0227-133032946/Undetermined_S0_L001_R1_001.fastq.gz')


system('scp PR0227-133032946/Undetermined_S0_L001_I1_001.fastq.gz PR0227-133032946/Undetermined_S0_L001_R1_001.fastq.gz PR0227-133032946/Undetermined_S0_L001_R2_001.fastq.gz microb120:/media/sequencing/Illumina/PR0227-133032946/Data/Intensities/BaseCalls/')


# R1.trimmed <- trimTailw(R1, 2, '?', 5)
# R2.trimmed <- trimTailw(R2, 2, '?', 5)
# R1.trimmed <- R1.trimmed[width(R1.trimmed) >= 100]
# R2.trimmed <- R2.trimmed[width(R2.trimmed) >= 100]
# R1.trimmed <- subseq(sread(R1.trimmed), 1, 100)
# R2.trimmed <- subseq(sread(R2.trimmed), 1, 100)
# 
# m1 <- vmatchPattern('GATCTGAATTCAGTGGCACAGCAGTTAGGNNNNNNGTAAGGNNNNNNAGATCTGGAATGAACTGGCC',
#                     R2.trimmed, fixed = FALSE, max.mismatch = 3)
# 
# m2 <- vmatchPattern('GATCTGAATTCAGTGGCACAGCAGTTAGG',
#                     R2.trimmed, fixed = FALSE, max.mismatch = 2)
# 
# m1.s <- unlist(lapply(startIndex(m1), length))
# m2.s <- unlist(lapply(startIndex(m2), length))
# 
# sum(m1.s)/length(m1.s)
# sum(m2.s)/length(m2.s)

# R1: LTR read {12 div}{8 bc}{ACCCTTTTAGTCAGTGTGGAAAATCTCTAGCA}
# R2: break point read {12 div}{8 bc}{GATCTGAATTCAGTGGCACAGCAGTTAGGNNNNNNGTAAGGNNNNNNAGATCTGGAATGAACTGGCC}
