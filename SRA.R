library(tidyverse)
library(RMySQL)
load(file = 'data/data.RData')
processedRunPath <- '/data/internal/geneTherapy/processedRuns'

dbConn1 <- DBI::dbConnect(MySQL(), group = 'intsites_miseq')
dbConn2 <- DBI::dbConnect(MySQL(), group = 'specimen_management')

# Create a local copy of the sample db table.
sampleData <- DBI::dbGetQuery(dbConn2, "select * from gtsp")

# Retrieve GTSP ids
query <- paste('select distinct miseqid from samples where ', paste0(' sampleName like "%', unique(intSites$GTSP), '%"', collapse = ' OR '))
runIDs <- unique(unname(unlist(DBI::dbGetQuery(dbConn1, query))))

f <- tibble(path = system(paste('find ', paste0(file.path(processedRunPath, runIDs, 'Data', 'demultiplexedReps'), collapse = ' '), ' -name "GTSP*"'), intern = TRUE),
            sample = unlist(str_extract_all(path, 'GTSP\\d+\\-\\d+'))) %>%
     group_by(sample) %>%
     mutate(nReads = n()) %>%
     ungroup() %>%
     mutate(GTSP = sub('\\-\\d+$', '', sample))

if(any(f$nReads != 2)) stop('Error -- one or more samples did not have 2 reads files.')
if(! all(intSites$GTSP %in% f$GTSP)) stop('Error -- all of the analysis samples are not in the retrieved sample table.')

# Remove samples which did not result in intSite calls.
f <- subset(f, GTSP %in% intSites$GTSP)

f$patient   <- toupper(sampleData[match(f$GTSP, sampleData$SpecimenAccNum),]$Patient)
f$cellType  <- toupper(sampleData[match(f$GTSP, sampleData$SpecimenAccNum),]$CellType)
f$timePoint <- toupper(sub('_', '.', sampleData[match(f$GTSP, sampleData$SpecimenAccNum),]$Timepoint))

# Metadata updates to match analysis.
f$patient   <- sub('^P', '', f$patient)
f$cellType  <- gsub('^CD34$', 'BM CD34+', f$cellType, ignore.case = TRUE)
f$cellType  <- gsub('CD3\\+\\scells', 'T cells', f$cellType, ignore.case = TRUE)
f$timePoint <- sub('M12', 'Y1',   f$timePoint)
f$timePoint <- sub('M18', 'Y1.5', f$timePoint)
f$timePoint <- sub('M24', 'Y2',   f$timePoint)
f$timePoint <- sub('M32', 'Y2.5', f$timePoint)

if(any(unname(unlist(is.na(f))))) stop('Error -- one or more data table entries are missing.')

# Crate a local copy of the reads for upload.
if(! dir.exists('readArchive')) dir.create('readArchive')
system(paste('cp ', paste0(f$path, collapse = ' '), 'readArchive/'))


# Create a local read archive for upload.
invisible(sapply(f$path, function(x) system(paste('cp ', x, 'readArchive/'))))

