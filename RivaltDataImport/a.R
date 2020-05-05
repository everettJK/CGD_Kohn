library(ShortRead)
library(dplyr)
options(stringsAsFactors = FALSE)

# LTR read (Rivat R1)
# Diversity    Barcode  LTR                              Genomic ...  
# AGATTGCACCTT AAACATCG ACCCTTTTAGTCAGTGTGGAAAATCTCTAGCA GGGTTCCCTAGTTAGCCAGAGAGCTGGCCAGTTCATTCCAGA...TTCAGATCTCCGTCTTTCTCATAGAATAAGATCGGAAGAGCACACGTCTGGACTCCAGTCACATCTCGTATGCCGTCTTCTGCTTGAAAAAAA
# CATTATTCTTGT AACAACCA ACCCTTTTAGTCAGTGTGGAAAATCTCTAGCA ATCTCGAACTCCTGGCCTCAGGTGATCCACCCGGGCCAGTTC...GCCTTACAATTCCCCTAACTGCTGTGCCACTGAATTCAGATCTCTCGGTTCAAATTCAACTAAGATCGGAAGAGCACACGTCTGAACTCCAGT
# TTGCGTTTTTTC AACAACCA ACCCTTTTAGTCAGTGTGGAAAATCTCTAGTA GGTTCCCTAGTTAGCCAGAGAGCTGGCCAGTTCATTCCAGAT...ATTCAGATCTCTCGGTTAAAATAAAAAGCAGATCGGAAGAGCACACGTCTGAACTCCAGTCACATCTCGTATGCCGTCTTCTGCTTGAAAAAA
# TTTGGGTGTGTA AAACATCG ACCCTTTTAGTCAGTGTGGAAAATCTCTAGCA GTGGGTTCCCTAGTTAGCCAGAGAGCTCGGCCAGTTCATTCC...TCAGATCTGGTTGTTACAATAACTTTTAGATCGGAAGAGCACACGTCTGAACTCCAGTCACATCTCGTATGCCGTCTTCTGCTTGAAAAAAAA
# 
# 
# Break point read  (Rivat R2)
# Diversity    Barcode  Linker1                       NNNNNN GTAAGG NNNNNN Linker2              Genomic...
# TGTTTTGATGAA AACAACCA GATCTGAATTCAGTGGCACAGCAGTTAGG TGCGTA GTAAGG ATATCA AGATCTGGAATGAACTGGCC CGCCCAA...AAGCTGCTAGAGATTTTCCACACTGACTAAAAGGGTCGATGTTTGAAAACAGAATAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGAT
# TAACTTTTTGTT AACAACCA GATCTGAATTCAGTGGCACAGCAGTTAGG TTCCCT GTAAGG CGCATC AGATCTGGAATGAACTGGCC CGTGGAC...TGCTTGCTAGAGATTTTCCACACTGACTAAAAGGGTTGGTTGTTGAGAATAAAAAGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGAT
# TGTGATTTGTTG AAACATCG GATCTGAATTCAGTGGCACAGCAGTTAGG TAGTCG GTAAGG GGAGGC AGATCTGGAATGAACTGGCC CGCCACT...AAAAGGGTTGGTTGTTAACCAACATAAAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTAAAAAAA
# TTATTCTCGACT AACTCACC GATCTGAATTCAGTGGCACAGCAGTTAGG TACGTG GTAAGG TGAAGG AGATCTGGGATGAACTGGCC GCTCTCT...CCCACTGCTTCTAGAGATTTTCCACACTGACTAAAAGGGTCGATGTTTTACCTTATTATAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

PR0245 <- read.table('PR0245.tmp.tsv', sep = '\t', header = FALSE)
PR0245$cellType  <- 'Neutrophils'
PR0245$timePoint <- PR0245$V4
PR0245$patient   <- PR0245$V2
PR0245$enz       <- PR0245$V5
PR0245$runID     <- 'PR0245-137703573'
PR0245$LTRcode   <- PR0245$V6
PR0245$Breakcode <- PR0245$V7
PR0245           <- select(PR0245, patient, cellType, enz, timePoint, runID, LTRcode, Breakcode)

# Copied to local dir for export.
LTR.barcodes <- read.table('../BiascoDataImport/LTR.barcodes', sep = '\t')
names(LTR.barcodes) <- c('id', 'seq')

breakPoint.barcodes <- read.table('../BiascoDataImport/breakPoint.barcodes', sep = '\t')
names(breakPoint.barcodes) <- c('id', 'seq')

breakPoint.barcodes$id <- gsub('\\-LC', '.LC', breakPoint.barcodes$id)


# Creation of data base records for imported samples.
#--------------------------------------------------------------------------------------------------
d <- PR0245
d$patient <- paste0('p', d$patient)
d$bcSeq <- paste0(LTR.barcodes[match(d$LTRcode, LTR.barcodes$id),]$seq, breakPoint.barcodes[match(d$Breakcode, breakPoint.barcodes$id),]$seq)
d$timePoint <- paste0('m', stringr::str_extract(d$timePoint, '\\d+'))
d$daysPostTrtmnt <- as.integer(stringr::str_extract(d$timePoint, '\\d+'))*30

samplesToReg <- d[! duplicated(paste(d$patient, d$timePoint, d$cellType)),]
samplesToReg$GTSP <- paste0('GTSP', 3228:(3228+nrow(samplesToReg)-1))

samplesToReg <- rowwise(samplesToReg) %>%
  mutate(sql = sprintf("insert into gtsp (SpecimenAccNum, SamplePatientCode, Patient, CellType, SampleCellType, Trial, DaysPostTrtmnt, Timepoint, PI, `PI Email`, Gender, refGenome) values ('%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s');",
                       GTSP, patient, patient, cellType, cellType, 'Rivalt CGD', daysPostTrtmnt, timePoint, 'Luca Biasco', 'l.biasco@ucl.ac.uk', 'm', 'hg38')) %>%
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


# Manually check that digests
group_by(d, runID, patient, cellType, timePoint, bcSeq) %>%
  summarise(cuts = paste0(enz, collapse = ', '),
            barcodes = paste0(unique(bcSeq), collapse = ', ')) %>%
  ungroup()


group_by(d, patient, cellType, timePoint, bcSeq) %>%
slice(1) %>%
mutate(alias = paste0(GTSP, '-1')) %>%
ungroup() %>%
select(alias, linkerSequence, bcSeq, gender, primer, ltrBit, largeLTRFrag, vectorSeq, refGenome) %>%
write.csv(row.names = FALSE, quote = FALSE, file = 'data-PR0245-137703573-sampleInfo.csv')

### system('cp *.csv /home/everett/data/BushmanGeneTherapy/sampleInfoFiles/')


#------------

R1 <- readFastq('PR0245-137703573/Undetermined_S0_L001_R1_001.fastq.org.gz')
R2 <- readFastq('PR0245-137703573/Undetermined_S0_L001_R2_001.fastq.org.gz')

R1.bcs <- substr(sread(R1), 13, 20)
R2.bcs <- substr(sread(R2), 13, 20)
bcs <- paste0(R1.bcs, R2.bcs)
readIDs <- as.character(ShortRead::id(R1))
sprintf("%.2f%%", (sum(R1.bcs %in% LTR.barcodes$seq) / length(R1.bcs))*100)
sprintf("%.2f%%", (sum(R2.bcs %in% breakPoint.barcodes$seq) / length(R2.bcs))*100)


bcOutput <- paste0('@', readIDs, '\n', bcs, '\n+\n', 'IIIIIIIIIIIIIIII')
if(file.exists('PR0245-137703573/Undetermined_S0_L001_I1_001.fastq.gz')) file.remove('PR0245-137703573/Undetermined_S0_L001_I1_001.fastq.gz')
fileConn <-file('PR0245-137703573/Undetermined_S0_L001_I1_001.fastq')
writeLines(bcOutput, fileConn)
close(fileConn)
system('gzip PR0245-137703573/Undetermined_S0_L001_I1_001.fastq')


# fasta2fastq <- function(x, code = '?'){
#   r <- readFastq(x)
#   r@quality <- FastqQuality(sapply(width(r), function(x) paste0(rep(code, x), collapse = '')))
#   file.remove(x)
#   writeFastq(r, file = x, compress = FALSE)
# }


# Narrow the reads to remove diversity and barcode NTs and write out the sequences.
# Swap R1 and R2 labels to match expected INSPIIRED input.
if(file.exists('PR0245-137703573/Undetermined_S0_L001_R2_001.fastq.gz')) file.remove('PR0245-137703573/Undetermined_S0_L001_R2_001.fastq.gz')
writeFastq(narrow(R1, 39, width(R1)), compress = TRUE, 
           file = 'PR0245-137703573/Undetermined_S0_L001_R2_001.fastq.gz')

if(file.exists('PR0245-137703573/Undetermined_S0_L001_R1_001.fastq.gz')) file.remove('PR0245-137703573/Undetermined_S0_L001_R1_001.fastq.gz')
writeFastq(narrow(R2, 56, width(R2)), compress = TRUE, 
           file = 'PR0245-137703573/Undetermined_S0_L001_R1_001.fastq.gz')


# Test
i <- readFastq('PR0245-137703573/Undetermined_S0_L001_I1_001.fastq.gz')
table(as.character(i@id) == as.character(R1@id))

r1 <- readFastq('PR0245-137703573/Undetermined_S0_L001_R2_001.fastq.gz')
table(as.character(i@id) == as.character(r1@id))


# Export reads to different server for intSiteCaller processing.
system('scp PR0245-137703573/Undetermined_S0_L001_I1_001.fastq.gz PR0245-137703573/Undetermined_S0_L001_R1_001.fastq.gz PR0245-137703573/Undetermined_S0_L001_R2_001.fastq.gz microb120:/media/sequencing/Illumina/PR0245-137703573/Data/Intensities/BaseCalls/')









