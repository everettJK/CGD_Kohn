library(gt23)
library(RMySQL)
library(tidyverse)
library(GenomicRanges)
library(RColorBrewer)
library(grDevices)
library(wordcloud)
library(scales)
library(reldist)
library(vegan)
library(png)
library(gridExtra)
library(ggpubr)
library(parallel)
library(tidyr)
library(xlsx)
options(useFancyQuotes = FALSE, stringsAsFactors = FALSE)
#/home/opt/R-3.4.0/bin/Rscript


# Retrieve sample data and intSites.
patients <- c("pB101", "pN102", "pU103", "pN104", "pN105", "pLW", "pAK", "pPM")
dbConn   <- DBI::dbConnect(RMySQL::MySQL(), group = 'specimen_management')
samples  <- DBI::dbGetQuery(dbConn, paste('select * from gtsp where Trial = "CGD" and  Patient in (', paste0(sQuote(patients), collapse = ',')  ,')'))
DBI::dbDisconnect(dbConn)

intSites <- getDBgenomicFragments(samples$SpecimenAccNum, 'specimen_management', 'intsites_miseq') %>%
            stdIntSiteFragments() %>%
            collapseReplicatesCalcAbunds() %>%
            annotateIntSites() %>%
            data.frame()

intSites.reps <- getDBgenomicFragments(samples$SpecimenAccNum, 'specimen_management', 'intsites_miseq') %>%
                 stdIntSiteFragments()

# Metadata updates.
intSites$cellType  <- gsub('^CD34$', 'BM CD34+', intSites$cellType, ignore.case = TRUE)
intSites$cellType  <- gsub('CD3\\+\\scells', 'T cells', intSites$cellType, ignore.case = TRUE)
intSites$patient   <- sub('^p', '', intSites$patient)


convertName <- function(x){
  nameConversion <- data.frame(internal=c('LW',        'AK',         'PM',        'B101',       'N102',      'U103',       'N104',      'N105'),
                               external=c('Patient 3', 'Patient 5',  'Patient 9', 'Patient 2',  'Patient 4', 'Patient 6',  'Patient 7', 'Patient 8'))
  nameConversion[match(tolower(x), tolower(nameConversion$internal)),]$external
}

intSites$patient <- as.character(convertName(intSites$patient))


# Create cellsPerSample column in order to select samples with enough inferred cells to justify using relative abundances.
intSites <- dplyr::group_by(intSites, GTSP) %>%
            dplyr::mutate(cellsPerSample = sum(estAbund)) %>%
            dplyr::ungroup()


# Create labeled nearest features.
intSites$labeledNearestFeature <- 
  intSites %>%
  mutate(labeledNearestFeature = paste0(nearestFeature, ' ')) %>% 
  mutate(labeledNearestFeature = ifelse(inFeature, paste0(labeledNearestFeature, '*'), labeledNearestFeature)) %>% 
  mutate(labeledNearestFeature = ifelse(abs(nearestOncoFeatureDist) <= 50000, paste0(labeledNearestFeature, '~'), labeledNearestFeature)) %>%
  mutate(labeledNearestFeature = ifelse(abs(nearestlymphomaFeatureDist) <= 50000, paste0(labeledNearestFeature, '!'), labeledNearestFeature)) %>%
  select(labeledNearestFeature) %>% unlist() %>% unname() 


# Select last nearest gene in comma delimit list if delimited.
# intSites$labeledNearestFeature <- unlist(lapply(strsplit(intSites$labeledNearestFeature, ','), function(x){ x[[length(x)]] }))


intSites$timePoint <- sub('M12', 'Y1',   intSites$timePoint)
intSites$timePoint <- sub('M18', 'Y1.5', intSites$timePoint)
intSites$timePoint <- sub('M24', 'Y2',   intSites$timePoint)
intSites$timePoint <- sub('M32', 'Y2.5', intSites$timePoint)

save(list = ls(all.names = TRUE), file = 'data/data.RData', envir = .GlobalEnv, compress = TRUE, compression_level = 9)



# Create summary table.
#--------------------------------------------------------------------------------------------------

calculateUC50 <- function(abund){
  stopifnot(is.vector(abund) & is.numeric(abund))
  abund <- abund[order(abund)]
  accum <- sapply(1:length(abund), function(i){sum(abund[1:i])})
  length(accum[accum >= sum(abund)/2])
}

ppNum <- function(n){ if(is.na(n)) return(''); format(n,big.mark=",", scientific=FALSE, trim=TRUE) }


summaryTable <- 
  group_by(intSites, GTSP) %>%
  summarise(timePointDays = timePointDays[1],
            Patient       = patient[1],
            Timepoint     = timePoint[1],
            CellType      = cellType[1],
            TotalReads    = ppNum(sum(reads)),
            InferredCells = ppNum(sum(estAbund)),
            UniqueSites   = ppNum(n_distinct(posid)),
            Gini          = sprintf("%.3f", gini(estAbund)),
            Chao1         = ppNum(round(estimateR(estAbund, index='chao')[2], 0)),
            Shannon       = sprintf("%.2f", diversity(estAbund)),
            Pielou        = sprintf("%.3f", diversity(estAbund)/log(specnumber(estAbund))),
            UC50          = ppNum(calculateUC50(estAbund))) %>%
  dplyr::ungroup() %>%
  dplyr::select(Patient, Timepoint, CellType, TotalReads, 
                InferredCells, UniqueSites, Gini, Shannon, Pielou, UC50, timePointDays) %>%
  dplyr::arrange(Patient, timePointDays) %>%
  dplyr::select(-timePointDays)

  
# Drop failed sample.
intSites <- intSites[-which(intSites$GTSP == 'GTSP1809'),]


# Flip to read relAbunds for older data.
# intSites <- bind_rows(lapply(split(intSites, paste(intSites$patient, intSites$cellType, intSites$timePoint)), function(x){
#  x$relAbundFrags <- x$relAbund
#  if(x$patient[1] %in% c('Patient 3', 'Patient 5', 'Patient 9')){
#    x$relAbund <- (x$reads / sum(x$reads))*100
#  }
# 
#  x
# }))


bind_rows(lapply(unique(intSites$cellType), function(x){
  tibble(cellType = x, maxRelAbund = sprintf("%.2f%%", max(subset(intSites, cellType == x & cellsPerSample >= 100)$relAbund)))
}))



# Create bad actor plot.  
#--------------------------------------------------------------------------------------------------
badActorPlot <- 
     dplyr::filter(intSites, ! timePoint %in% c('M1.5')) %>%
     dplyr::group_by(GTSP) %>%
     dplyr::mutate(abundPerGTSP = sum(estAbund)) %>%
     dplyr::ungroup() %>%
     dplyr::filter(abundPerGTSP >= 100, 
                   abs(nearestFeatureDist) <= 50000,
                   nearestFeature %in% c('SETBP1', 'PRDM16', 'MECOM')) %>%
     dplyr::mutate(cellType = paste0('  ', cellType, '    ')) %>%
     dplyr::select(patient, GTSP, cellType, timePoint, estAbund, relAbund, nearestFeature, posid) %>%
     dplyr::group_by(GTSP, nearestFeature) %>%
       dplyr::top_n(1, wt = relAbund) %>%
       dplyr::slice(1) %>%
     dplyr::ungroup() %>%
     #dplyr::mutate(patient = factor(patient, levels = c('B101', 'N102', 'N104', 'U103', 'AW', 'L102', 'L103'))) %>%
     ggplot(aes(x = timePoint, y = relAbund/100, fill = cellType)) +
     theme_bw() +
     geom_bar(stat = 'identity', position = position_dodge2(preserve = "single")) +
     scale_fill_manual(name = 'Cell type', values = brewer.pal(10, "Spectral")) +
     facet_grid(patient~nearestFeature) + 
     scale_y_continuous(labels = scales::percent_format(accuracy=.01)) +  
     labs(x = 'Time point', y = 'Relative abundance of most prominent clone') + 
     theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black"),
          strip.text = element_text(size = 12),
          axis.text.y =element_text(size = 12),
          axis.text.x =element_text(size = 12, angle = -90, hjust=0),
          axis.title = element_text(size = 16),
          legend.position = "bottom",
          legend.text=element_text(size=12),
          legend.title=element_blank()) + 
    guides(fill=guide_legend(ncol=6)) 

ggsave(badActorPlot, file = 'figures/SuppFig5.pdf', useDingbats = FALSE, height = 10, units = 'in')
ggsave(badActorPlot, file = 'figures/SuppFig5.svg', height = 10, units = 'in')



# Create plot data object.
plotData <-
  dplyr::group_by(intSites, patient, cellType, timePoint) %>%
  dplyr::summarise(nSites  = n_distinct(posid),
                   shannon = vegan::diversity(estAbund),
                   nSitesChao1 = round(vegan::estimateR(estAbund, index='chao')[2], 0),
                   days = timePointDays[1]) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(days) %>%
  dplyr::mutate(timePoint = factor(timePoint, levels = gtools::mixedsort(unique(timePoint)))) %>%
  dplyr::mutate(patient = factor(patient, levels = gtools::mixedsort(unique(patient)))) %>%
  dplyr::filter(timePoint != 'M1.5') 


# Fig3a
obsIntSites <- 
  ggplot(plotData, aes(timePoint, log10(nSites), group = cellType, color = cellType)) +
  theme_bw() +
  geom_point(size = 3) + 
  scale_color_manual(name = 'Cell type', values = rev(brewer.pal(10, "Paired"))) +
  geom_line() +
  facet_grid(patient~.) + 
  annotate('text', x = 10.7, y = 0.8, label = levels(plotData$patient), size = 4) +
  theme(strip.background = element_blank(), strip.text = element_blank(), legend.position="bottom") +
  labs(x = 'Time point', y = 'Unique Integrations (log10)') +
  guides(color=guide_legend(ncol=4)) +
  ylim(c(0,7)) +
  theme(text = element_text(size=18),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 5, l = 0)),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))


ggsave(obsIntSites, file = 'figures/3a.pdf', useDingbats = FALSE, height = 10, units = 'in')
ggsave(obsIntSites, file = 'figures/3a.svg', height = 10, units = 'in')



obsIntSitesChao1 <- 
  ggplot(plotData, aes(timePoint, log10(nSitesChao1), group = cellType, color = cellType)) +
  theme_bw() +
  geom_point(size = 3) + 
  #scale_color_manual(name = 'Cell type', values = rev(brewer.pal(10, "Paired"))) +
  scale_color_manual(name = 'Cell type', values = c('#A5CDE1','#1F78B3','#9958A2','#F48332','#B3D789','#349F48','#F59898','#E11F26','#FCBE6E')) +
  geom_line() +
  facet_grid(patient~.) + 
  annotate('text', x = 10.7, y = 0.8, label = levels(plotData$patient), size = 4) +
  theme(strip.background = element_blank(), strip.text = element_blank(), legend.position="bottom") +
  labs(x = 'Time point', y = 'Chao1 estimate of integrated cells (log10)') +
  guides(color=guide_legend(ncol=4)) +
  ylim(c(0,8)) +
  theme(text = element_text(size=18),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 5, l = 0)),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))

ggsave(obsIntSitesChao1, file = 'figures/3a_Chao1.pdf', useDingbats = FALSE, height = 10,  width = 5, units = 'in')
ggsave(obsIntSitesChao1, file = 'figures/3a_Chao1.svg', height = 10, units = 'in')



obsIntSites2 <- 
  ggplot(plotData, aes(cellType, nSites, fill = cellType)) +
  theme_bw() +
  geom_bar(stat='identity') +
  scale_fill_manual(name = 'Cell type', values = brewer.pal(10, "Spectral")) +
  labs(x = '', y = 'Integrations') +
  facet_grid(patient~timePoint, scales = 'free_y') +
  scale_y_continuous(labels = comma) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        strip.background = element_blank()) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


ChaoIntSites <- ggplot(plotData, aes(cellType, nSitesChao1, fill = cellType)) +
  theme_bw() +
  geom_bar(stat='identity') +
  scale_fill_manual(name = 'Cell type', values = brewer.pal(10, "Spectral")) +
  labs(x = '', y = 'Integrations (Chao1 estimate)') +
  facet_grid(patient~timePoint, scales = 'free_y') +
  scale_y_continuous(labels = comma) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) 


obsIntSites3 <- 
  ggplot(subset(plotData, timePoint != 'D0'), aes(timePoint, nSites, group = cellType, color = cellType)) +
  theme_bw() +
  geom_point(size = 3) + 
  scale_color_manual(name = 'Cell type', values = rev(brewer.pal(10, "Spectral"))) +
  geom_line() +
  scale_y_continuous(label=comma) +
  facet_grid(patient~.) + 
  annotate('text', x = 10, y = 20000, label = unique(plotData$patient), size = 4) +
  theme(strip.background = element_blank(), strip.text = element_blank(), legend.position="bottom") +
  labs(x = 'Time point', y = 'Integrations') +
  guides(color=guide_legend(ncol=4)) +
  theme(text = element_text(size=18),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 5, l = 0)),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))


obsIntSitesDiversity <- 
  ggplot(plotData, aes(timePoint, shannon, group = cellType, color = cellType)) +
  theme_bw() +
  geom_point(size = 3) + 
  scale_color_manual(name = 'Cell type', values = brewer.pal(10, "Spectral")) +
  geom_line() +
  facet_grid(patient~.) + 
  annotate('text', x = 11.2, y = 1.3, label = unique(plotData$patient), size = 4) +
  theme(strip.background = element_blank(), strip.text = element_blank(), legend.position="bottom") +
  labs(x = 'Time point', y = 'Shannon index') +
  guides(color=guide_legend(ncol=4)) +
  ylim(c(0,13)) +
  theme(text = element_text(size=18),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 5, l = 0)),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))


clonesToTrack <- 10

d <- dplyr::filter(intSites, 
                   cellType == 'Neutrophils', 
                   cellsPerSample >= 100, 
                   ! timePoint %in% 'M1.5') %>%
  dplyr::mutate(timePoint = factor(timePoint, levels = gtools::mixedsort(unique(timePoint)))) %>%
  dplyr::group_by(patient) %>%
  dplyr::arrange(desc(relAbund)) %>%
  dplyr::mutate(timePoints = list(gtools::mixedsort(unique(as.character(timePoint))))) %>%
  dplyr::filter(posid %in% unique(posid)[1:clonesToTrack]) %>%
  dplyr::ungroup()

max(d$relAbund)

d$labeledNearestFeature <- sub("BSCL2,HNRNPUL2-BSCL2 *", 'BSCL2 *', d$labeledNearestFeature)

timePoints <- gtools::mixedsort(unique(unlist(d$timePoints)))

Neutrophils_maps <- lapply(split(d, d$patient), function(x){
  x$timePoint <- factor(as.character(x$timePoint), levels = timePoints)
  x <- x[order(x$relAbund, decreasing = TRUE),]
  x$labeledNearestFeature <- as.character(x$labeledNearestFeature)
  x$labeledNearestFeature <- factor(x$labeledNearestFeature, levels = rev(unique(x$labeledNearestFeature)))
  
  ggplot(x, aes(timePoint, labeledNearestFeature, fill = relAbund/100)) +
    theme_bw() +
    scale_x_discrete(drop = FALSE) +
    geom_tile(color = 'white') +
    labs(x = 'Time post treatment', y = 'Nearest gene') +
    scale_fill_gradient(name = 'Relative abundance',
                        low = "#cce5ff", high = "#000fdb", 
                        # limits = c(0, 0.20), breaks = c(0.05, 0.10, 0.15, 0.20), label = scales::percent_format(accuracy = 0.01)) + 
                        limits = c(0, 0.025), breaks = c(0.01, 0.02, 0.025), label = scales::percent_format(accuracy = 0.1)) +
    labs(x = '', y = '') +
    theme(panel.border = element_blank(), 
          text = element_text(size=14),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"),
          legend.position="bottom") +
    guides(fill= guide_colorbar(barwidth = 15, title.position = 'top')) +
    ggtitle(x$patient[1])
})


gA <- ggplotGrob(Neutrophils_maps[[1]] + theme(legend.position = "none"))
gB <- ggplotGrob(Neutrophils_maps[[2]] + theme(legend.position = "none"))
gC <- ggplotGrob(Neutrophils_maps[[3]] + theme(legend.position = "none"))
gD <- ggplotGrob(Neutrophils_maps[[4]] + theme(legend.position = "none"))
gE <- ggplotGrob(Neutrophils_maps[[5]] + theme(legend.position = "none"))
gF <- ggplotGrob(Neutrophils_maps[[6]] + theme(legend.position = "none"))
gG <- ggplotGrob(Neutrophils_maps[[7]] + theme(legend.position = "none"))

legend <- get_legend(Neutrophils_maps[[4]])
legend <- as_ggplot(legend)

maxWidth = grid::unit.pmax(gA$widths[2:5], gB$widths[2:5], gC$widths[2:5], gD$widths[2:5], gE$widths[2:5], gF$widths[2:5], gG$widths[2:5])

gA$widths[2:5] <- as.list(maxWidth)
gB$widths[2:5] <- as.list(maxWidth)
gC$widths[2:5] <- as.list(maxWidth)
gD$widths[2:5] <- as.list(maxWidth)
gE$widths[2:5] <- as.list(maxWidth)
gF$widths[2:5] <- as.list(maxWidth)
gG$widths[2:5] <- as.list(maxWidth)


ggsave(grid.arrange(gA, gB, gC, gD, gE, gF, gG, legend, ncol=1), file = 'figures/3b.pdf', useDingbats = FALSE, height = 20, units = 'in')
ggsave(grid.arrange(gA, gB, gC, gD, gE, gF, gG, legend, ncol=1), file = 'figures/3b.svg', height = 20, units = 'in')




TopRelAbund_Neutrophils <- 
  dplyr::filter(intSites, 
                cellType == 'Neutrophils', 
                cellsPerSample >= 50, 
                ! timePoint %in% 'M1.5') %>%
  dplyr::mutate(timePoint = factor(timePoint, levels = gtools::mixedsort(unique(timePoint)))) %>%
  dplyr::group_by(patient) %>%
  dplyr::arrange(desc(relAbund), nchar(labeledNearestFeature)) %>%
  dplyr::mutate(posids = list(unique(posid)[1:clonesToTrack])) %>%
  dplyr::mutate(maxTimePoint = max(as.integer(timePoint)) + 0.55) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(patient) %>%
  dplyr::filter(posid %in% unlist(posids[1])) %>%
  dplyr::ungroup() %>%
  ggplot(aes(timePoint, labeledNearestFeature, fill = relAbund/100)) +
  theme_bw() +
  geom_tile() +
  geom_vline(aes(xintercept = maxTimePoint)) +
  labs(x = 'Time post treatment', y = 'Nearest gene') +
  facet_grid(patient~., scales = 'free') +
  scale_fill_gradient(name = 'Relative\nabundance', na.value = 'red',
                      low = "#cce5ff", high = "#000fdb",
                      labels = percent) +
  theme(strip.background = element_blank(),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))

save(list = ls(all.names = TRUE), file = 'data/report.RData', 
     envir = .GlobalEnv, compress = TRUE, compression_level = 9)


stop()

# Dev code below
#--------------------------------------------------------------------------------------------------
































# Gene enrichment between d0 and final patient time points.
#--------------------------------------------------------------------------------------------------

CPUs <- 30

intSites.stdChrom <- subset(intSites, seqnames %in% paste0('chr', c(1:23, 'X', 'Y')))
intSites.stdChrom <- mutate(intSites.stdChrom, nearestFeature = strsplit(intSites.stdChrom$nearestFeature, ','))
intSites.stdChrom <- unnest(intSites.stdChrom, nearestFeature)

genes <- tibble(gene = unique(intSites.stdChrom$nearestFeature))
genes$n <- ntile(1:nrow(genes), CPUs)

cluster <- makeCluster(CPUs)

nearestFeatureBoundaries <-
  bind_rows(parLapply(cluster, split(genes, genes$n), function(x){
    library(dplyr)
    library(GenomicRanges)
    bind_rows(lapply(x$gene, function(x2){
      o <- subset(gt23::hg38.refSeqGenesGRanges, toupper(name2) == toupper(x2))
      o <- o[! grepl('_', seqnames(o))]
      tibble(gene = x2, chr = as.character(seqnames(o))[1], start = min(start(o)), end = max(end(o)))
    }))
  }))

stopCluster(cluster)


patient_early_late_sites <- lapply(c('Patient 2',  'Patient 4', 'Patient 6',  'Patient 7'), function(x){
  o <- subset(intSites.stdChrom, patient == x)
  a <- subset(o, timePoint == 'D0')
  b <- subset(o, timePointDays == max(o$timePointDays) & cellType == 'PBMC')
  list(a, b)
})
names(patient_early_late_sites) <- c('Patient 2',  'Patient 4', 'Patient 6',  'Patient 7')


geneEnrichments.d0_vs_LastTimePoint  <-
  bind_rows(lapply(patient_early_late_sites, function(x){
    
    patientNearestFeatureBoundaries   <- subset(nearestFeatureBoundaries, gene %in% c(x[[1]]$nearestFeature, x[[2]]$nearestFeature))
    patientNearestFeatureBoundaries$n <- ntile(1:nrow(patientNearestFeatureBoundaries), CPUs)
    
    cluster <- makeCluster(CPUs)
    clusterExport(cluster, varlist = c('intSites.stdChrom', 'patientNearestFeatureBoundaries'), envir = environment())
    
    r <- bind_rows(parLapply(cluster, split(patientNearestFeatureBoundaries, patientNearestFeatureBoundaries$n), function(x2){
      library(dplyr)
      
      bind_rows(lapply(x2$gene, function(x3){
        g <- subset(patientNearestFeatureBoundaries, gene == x3)
        
        # Find patient sites within this gene.
        i1 <- subset(intSites.stdChrom, patient = x[[1]]$patient[1], seqnames == g$chr & start >= g$start & start <= g$end)
        i1.maxLateAbund <- max(subset(x[[2]], posid %in% i1$posid)$estAbund)
        
        m1 <- matrix(c(n_distinct(subset(x[[1]], posid %in% i1$posid)$posid),   n_distinct(subset(x[[2]], posid %in% i1$posid)$posid),
                       n_distinct(subset(x[[1]], ! posid %in% i1$posid)$posid), n_distinct(subset(x[[2]], ! posid %in% i1$posid)$posid)),
                     nrow = 2, byrow = TRUE, dimnames = list(c('In', 'Out'), c('Early', 'Late')))
        
        # Expand gene boundaries by 50K.
        i2 <- subset(intSites.stdChrom, patient = x[[1]]$patient[1], seqnames == g$chr & start >= (g$start - 50000) & start <= (g$end + 50000))
        i2.maxLateAbund <- max(subset(x[[2]], posid %in% i2$posid)$estAbund)
        
        m2 <- matrix(c(n_distinct(subset(x[[1]], posid %in% i2$posid)$posid), n_distinct(subset(x[[2]], posid %in% i2$posid)$posid),
                       n_distinct(subset(x[[1]], ! posid %in% i2$posid)$posid), n_distinct(subset(x[[2]], ! posid %in% i2$posid)$posid)),
                     nrow = 2, byrow = TRUE, dimnames = list(c('In', 'Out'), c('Early', 'Late')))
        
        matrix2string <- function(m) paste0('Early: ', m[1,1], ' in / ', m[2,1], ' out | Late: ', m[1,2], ' in / ', m[2,2], ' out')
        pChange <- function(m) ((m[1,2] / sum(m[,2])) - (m[1,1] / sum(m[,1])))*100
        
        tibble(pateint = x[[1]]$patient[1], gene = x3, mat1 = matrix2string(m1), p1 = fisher.test(m1)$p.val,
               p1.change = sprintf("%.2f%%", pChange(m1)), mat2 = matrix2string(m2), p2 = fisher.test(m2)$p.val,
               p2.change = sprintf("%.2f%%", pChange(m2)), maxLateAbund1 = i1.maxLateAbund, maxLateAbund2 = i2.maxLateAbund)
      }))
    }))
    
    stopCluster(cluster)
    gc()
    r
  }))

geneEnrichments.d0_vs_LastTimePoint$p1.adj <- p.adjust(geneEnrichments.d0_vs_LastTimePoint$p1)
geneEnrichments.d0_vs_LastTimePoint$p2.adj <- p.adjust(geneEnrichments.d0_vs_LastTimePoint$p2)
geneEnrichments.d0_vs_LastTimePoint$oncogene <- toupper(geneEnrichments.d0_vs_LastTimePoint$gene) %in% toupper(gt23::hg38.oncoGeneList)
geneEnrichments.d0_vs_LastTimePoint$lymphomaGene <- toupper(geneEnrichments.d0_vs_LastTimePoint$gene) %in% toupper(gt23::hg38.lymphomaGenesList)

geneEnrichments.d0_vs_LastTimePoint.p1 <- subset(geneEnrichments.d0_vs_LastTimePoint, p1.adj <= 0.05)
geneEnrichments.d0_vs_LastTimePoint.p2 <- subset(geneEnrichments.d0_vs_LastTimePoint, p2.adj <= 0.05)

wb <- createWorkbook()
sheet1 <- createSheet(wb, sheetName = 'in TU')
sheet2 <- createSheet(wb, sheetName = 'in TU +- 50KB')

o1 <- dplyr::select(geneEnrichments.d0_vs_LastTimePoint.p1, pateint, gene, p1.adj, oncogene, lymphomaGene, p1.change, maxLateAbund1, mat1) %>%
  dplyr::arrange(pateint, p1.adj) %>%
  dplyr::mutate(p1.adj = formatC(p1.adj, format = "e", digits = 2))
names(o1) <- c('Patient', 'Gene', 'Adjusted p-val', 'Oncogene', 'Lymphoma Gene', '% change (d0 vs final timepoint)', 'Max. abund. in final time point', 'Matrix')

o2 <- dplyr::select(geneEnrichments.d0_vs_LastTimePoint.p2, pateint, gene, p2.adj, oncogene, lymphomaGene, p2.change, maxLateAbund2, mat2) %>%
  dplyr::arrange(pateint, p2.adj) %>%
  dplyr::mutate(p2.adj = formatC(p2.adj, format = "e", digits = 2))
names(o2) <- c('Patient', 'Gene', 'Adjusted p-val', 'Oncogene', 'Lymphoma Gene', '% change (d0 vs final timepoint)', 'Max. abund. in final time point', 'Matrix')

addDataFrame(o1, sheet1)
addDataFrame(o2, sheet2)
saveWorkbook(wb, file = 'reports/d0vsLastTimePoint.xlsx')
rm(o1, o2, wb, sheet1, sheet2, intSites.stdChrom, genes)



# Shared sites
#--------------------------------------------------------------------------------------------------

d <- dplyr::group_by(intSites, patient, cellType, posid) %>%
     dplyr::mutate(count = n()) %>%
     dplyr::ungroup() %>%
     dplyr::filter(count > 1) %>%
     dplyr::select(patient, cellType, count, posid) %>%
     dplyr::group_by(patient, cellType) %>%
     dplyr::mutate(sites = list(unique(posid))) %>%
     dplyr::slice(1) %>%
     dplyr::ungroup()
  


invisible(lapply(split(d, d$patient), function(x){
  a <- lapply(1:nrow(x), function(x2) unlist(x[x2,]$sites))
  names(a) <- x$cellType
  
  png(file = paste0(x$patient[1], '_sharedSites.png'),  width = 5, height = 5, units = "in", res = 150)
  a2 <- a[names(a) %in% c("B cells", "Neutrophils", "NK cells", "PBMC", "T cells", "GM", "Whole Blood")]
  if(length(a2) > 5) a2 <- a[names(a) %in% c("B cells", "Neutrophils", "NK cells", "PBMC", "T cells", "GM")]
  if(length(a2) >= 2){
    upset(fromList(a2), mainbar.y.label='Integration sites', text.scale = 1.5)
    dev.off()
  }
}))



# Track clones with high relative abundances for each subject.
#--------------------------------------------------------------------------------------------------

clonesToTrack <- 10

# Define clones to track.
clones <- 
  dplyr::filter(intSites, cellsPerSample >= 50, 
                cellType %in% c('BM CD34', 'PBMC'),
                ! timePoint %in% 'M1.5') %>%
  dplyr::group_by(patient) %>%
  dplyr::arrange(desc(relAbund), nchar(labeledNearestFeature)) %>%
  dplyr::filter(! duplicated(posid)) %>%
  dplyr::select(patient, posid, relAbund, nearestFeature, labeledNearestFeature) %>%
  dplyr::slice(1:clonesToTrack) %>%
  dplyr::ungroup()

# Subset clones to match top clone selection criteria.
o <- dplyr::filter(intSites, cellsPerSample >= 50, 
                   cellType %in% c('BM CD34', 'PBMC'), 
                   ! timePoint %in% 'M1.5') %>% 
     dplyr::select(patient, timePoint, posid, relAbund, nearestFeature, labeledNearestFeature) 



# Helper function to create nearest gene indices to determine if nearest genes are repeated.
seenBefore <- function(x){
  s <- list()
  unname(sapply(x, function(x2){
    if(x2 %in% names(s)){
      s[[x2]] <<- s[[x2]]+1
    } else {
      s[[x2]] <<- 1
    }
    
    s[[x2]]
  }))
}

# Create plot data.
d <- bind_rows(lapply(split(o, paste(o$patient, o$timePoint)), function(x){
       a <- subset(clones, patient == x$patient[1])
       a$timePoint <- x$timePoint[1]
       a$relAbund <- 0
       
       x <- subset(x, posid %in% a$posid)
       if(nrow(x) == clonesToTrack){
         return(x)
       } else if(nrow(x) == 0){
         return(a)
       } else {
         return(bind_rows(a[-which(a$posid %in% x$posid),], x))
       }
    })) %>%
    dplyr::group_by(patient, timePoint) %>% 
    dplyr::arrange(desc(relAbund)) %>%
    dplyr::mutate(labeledNearestFeature = paste0(labeledNearestFeature , ' / ', seenBefore(labeledNearestFeature))) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(timePoint = factor(timePoint, levels = unique(gtools::mixedsort(timePoint))))


topClonalRelativeAbunds <- lapply(split(d, d$patient), function(x){
  pd <- ggplot2::position_dodge(0.15)
  ggplot(x, aes(timePoint, relAbund/100, color = labeledNearestFeature, group = labeledNearestFeature)) +
  theme_bw() +
  geom_point(position = pd, size = 3) + 
  geom_line(position = pd) +
  scale_x_discrete(drop=FALSE) +
  scale_color_manual(name = '', values = brewer.pal(12, "Paired")) +
  scale_y_continuous(labels = scales::percent_format(accuracy=0.02)) +
  labs(x = '', y = 'Relative Abundance') +
  ggtitle(x$patient[1]) +
  theme(text = element_text(size=14), 
        legend.text=element_text(size=9),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))
})
  



# Relabund plots
#--------------------------------------------------------------------------------------------------

numClones <- 10
intSites$posidLabel <- paste0(intSites$labeledNearestFeature, '\n', intSites$posid)
d <- subset(intSites, cellType == 'PBMC' & timePoint != 'M1.5')

abundantClones <- lapply(split(d, d$patient), function(x){

  o <- bind_rows(lapply(split(x, x$timePoint), function(x2){
         x2$totalCells <- ppNum(sum(x2$estAbund))
         
         x2 <- x2[order(x2$relAbund, decreasing = TRUE),][1:numClones,]
         x2$labeledNearestFeature <-  unlist(lapply(strsplit(x2$labeledNearestFeature, ','), '[[', 1))
         x2$timePoint2 <- paste0(c(x2$timePoint[1], '\n', x2$labeledNearestFeature), collapse = '\n')
         x2$nearestFeature2 <- paste('Top', as.character(1:numClones))
    
         lowAbundData <- dplyr::mutate(x2, nearestFeature2 = 'Low Abund',
                                       posidLabel = 'LowAbund',
                                       totalCells = totalCells[1],
                                       relAbund   = 100) %>% dplyr::slice(1)
    
         lowAbundData$relAbund <- 100 - sum(x2$relAbund)
         bind_rows(lowAbundData,  x2)
       }))
  
  o$nearestFeature2 <- factor(o$nearestFeature2, levels = c('Low Abund', paste('Top', as.character(1:10))))
  
  totalCellLabel <- unname(unlist(lapply(split(o, o$timePoint2), function(x) x$totalCells[1])))
  
  ggplot(o, aes(timePoint2, relAbund/100, fill=nearestFeature2)) +
    theme_bw() +
    scale_x_discrete(drop=FALSE) + 
    geom_bar(stat='identity') +
    ggtitle(paste0(o$patient[1], ' PBMC samples')) + 
    scale_fill_manual(name = '', values = c('#EEEEEE', rev(brewer.pal(10, "Spectral")))) +
    theme(legend.position='right') +
    scale_y_continuous(labels = scales::percent) + 
    guides(fill=guide_legend(ncol=1)) +
    labs(x = '', y = 'Relative clonal abundance') +
    annotate('text', x=1:length(totalCellLabel), y=1.04, label=totalCellLabel, size=2.7, angle=45, hjust=0.5) 
})



# Word clouds
#--------------------------------------------------------------------------------------------------

outputDir <- 'wordClouds'
o <- data.frame(subset(intSites, timePoint == 'Y1' & cellType == 'PBMC' & estAbund >= 2))
invisible(lapply(split(o, o$patient), function(x){
  x <- x[order(x$estAbund, decreasing = TRUE),]
  maxWords <- 200
  if(nrow(x) < maxWords) maxWords <- nrow(x)
  x <- x[1:maxWords,]
  w <- setNames(x$estAbund, x$labeledNearestFeature)
  
  #png(file = file.path(outputDir, paste0(x$patient[1], '.PBMC_Y1.', 
  #                                       x[maxWords,]$estAbund, '.', x[1,]$estAbund, '.wordCloud.png')), res = 150)
  
  svg(file = file.path(outputDir, paste0(x$patient[1], '.PBMC_Y1.', x[maxWords,]$estAbund, '.', x[1,]$estAbund, '.wordCloud.svg')))
  
  wordcloud(names(w), w, random.order=FALSE, colors=colorRampPalette(brewer.pal(12, "Paired"))(maxWords), rot.per=0, max.words = 200, scale = c(2, 0.2), min.freq = 2)
  
  dev.off()
}))


# All subjects, all cell types
o <- o[order(o$estAbund, decreasing = TRUE),]
maxWords <- 100
x <- o[1:maxWords,]
w <- setNames(x$estAbund, x$labeledNearestFeature)
cellType <- 'allCellTypes'
png(file = file.path(outputDir, paste0('allSubjects_allCellTypes.', 
                                       x[maxWords,]$estAbund, '.', x[1,]$estAbund, '.wordCloud.png')), res = 150)

wordcloud(names(w), w, random.order=FALSE, colors=colorRampPalette(brewer.pal(12, "Paired"))(maxWords), rot.per=0, max.words = 100, scale = c(2, 0.2))


dev.off()



# Heat maps
#--------------------------------------------------------------------------------------------------

#
# Create genomic heat maps
#

Rscript_path <- '/home/opt/R-3.4.0/bin/Rscript'

# Remove intSites in non-standard chromosomes since we already have to remove some sites
# in order not exceede the capacity of the calcualtion.

d <- intSites
d <- subset(d, seqnames %in% c(paste0('chr', 1:22), 'chrX', 'chrY'))

# Create a timepoint / cell type splitting factor.
d <- subset(d, timePoint == 'Y1' & cellType == 'PBMC')
d$s <- paste(d$timePoint, d$cellType, sep='.')

# Create the sample files for the genomic heatmap software.
# Rather than replicate samples, use the splitting factor as the sample otherwise
# the returned heatmap would have an unreasonable number of columns.


heatMapSamples <- bind_rows(lapply(split(d, d$patient), function(x){
  df <- do.call(rbind, lapply(split(x, x$s), function(x2){
    data.frame(sampleName=x2$GTSP[1],
               GTSP=x2$GTSP[1],
               patient=x2$patient[1],
               cellType=x2$cellType[1])
  }))
  
  df
}))


heatMapSites <- bind_rows(lapply(heatMapSamples$GTSP, function(x){
  o <- subset(d, GTSP == x)
  df <- data.frame(o$seqnames, o$strand, o$start, x, 'hg38') 
  names(df) <- c('seqnames', 'strand', 'position', 'sampleName', 'refGenome')
  df
}))

set.seed(1)
if(nrow(heatMapSites) > 20000) heatMapSites <- heatMapSites[sample(1:nrow(heatMapSites), 20000),]


write.table(heatMapSamples, quote = FALSE, row.names = FALSE, col.names = TRUE, file = 'heatMapSamples', sep=',')
write.table(heatMapSites, quote = FALSE, row.names = FALSE, col.names = TRUE, file = 'heatMapSites', sep=',')

comm <- paste(Rscript_path, './software/genomicHeatmapMaker-from_input/genomic_heatmap_from_file.R heatMapSamples', 
              '-c ./software/genomicHeatmapMaker-from_input/INSPIIRED.yml',
              '-o', file.path('./tmp', 'genomicHeatMap'),
              '-f heatMapSites',
              '-r hg38')



# Parse the heat map maker output.

genomicHeatmap <- within(
  list(), {
    heatmap_sample_info <- read.csv('heatMapSamples')
    gen_heatmap <- readRDS('tmp/genomicHeatMap/roc.res.rds')
    
    heatmap_scale <- seq(0.2, 0.8, 0.1)
    gen_heatmap_colors <- colorspace::diverge_hsv(
      length(heatmap_scale), h = c(240, 0), v = 1, power = 1)
    
    select_gen_features <- row.names(gen_heatmap$ROC)
    select_gen_features <- c(
      "boundary.dist", "start.dist", "general.width", "gene.width",
      "within_refSeq_gene", "refSeq_counts.10k", "refSeq_counts.100k",
      "refSeq_counts.1M", "GC.100", "GC.1k", "GC.10k", "GC.100k", "GC.1M",
      "CpG_counts.1k", "CpG_counts.10k", "CpG_density.10k", "CpG_density.100k",
      "CpG_density.1M", "DNaseI_count.1k", "DNaseI_count.10k",
      "DNaseI_count.100k", "DNaseI_count.1M")
    
    gen_heatmap$ROC <- gen_heatmap$ROC[select_gen_features,]
    
    #heatmap_sample_levels <- c("PCR_amplicon", "Huh7_d0", "Huh7_d3", "Huh7_w3", "Humanized_mouse_m1", "WAS_d0", "CART_d0")
    heatmap_sample_levels <- c("GTSP1574", "GTSP1831", "GTSP2367", "GTSP2361")
    #heatmap_sample_levels <- c("B101", "N102", "N104", "U103")
    
    heatmap_figure_labels <- heatmap_sample_levels
    
    stat_cuts <- c(0, 0.001, 0.01, 0.05, 1)
    gen_comp_stats <- structure(cut(
      gen_heatmap$pvalues$op[select_gen_features, 1],
      stat_cuts,
      labels = c("***", " **", " * ", "   "),
      include.lowest = TRUE),
      names = select_gen_features)
    gen_row_names <- paste0(names(gen_comp_stats), " - ", gen_comp_stats)
    
    plot_data <- gen_heatmap$ROC %>%
      reshape2::melt() %>%
      mutate(
        feat = Var1,
        comp.sym = gen_comp_stats[Var1],
        Var1 = paste0(Var1, " - ", comp.sym),
        Var1 = factor(Var1, levels = gen_row_names),
        Var2 = factor(Var2, levels = heatmap_sample_levels),
        grp = " ",
        sig = as.vector(gen_heatmap$pvalues$np[select_gen_features,]),
        sym = cut(
          sig, stat_cuts, labels = c("***", " **", " * ", "   "),
          include.lowest = TRUE))
    
    levels(plot_data$Var2) <- heatmap_figure_labels
    
    #hot fix
    plot_data$Var2 <- heatMapSamples[match(plot_data$Var2, heatMapSamples$GTSP),]$patient
    
    ## browser()
    
    gen_plot <- ggplot(plot_data, aes(x = Var2, y = Var1, fill = value)) +
      geom_tile(color = 'black') +
      geom_text(aes(label = sym), color = "black", size = 3, nudge_y = -0.15) +
      scale_x_discrete(position = "top") +
      scale_fill_gradient2(
        breaks = c(0.2, 0.4, 0.6, 0.8),
        low = gen_heatmap_colors[1],
        mid = gen_heatmap_colors[round(length(heatmap_scale)/2)],
        high = gen_heatmap_colors[length(heatmap_scale)],
        midpoint = 0.5) +
      guides(fill = guide_colorbar(
        title.position = "left", title.hjust = 0.5,
        direction = "horizontal")) +
      labs(x = NULL, y = NULL, fill = "ROC\nScore") +
      # custom_theme +
      theme(
        legend.position = "bottom",
        axis.text.y = element_text( angle = 0, hjust = 0, vjust = 0.5, size = 12),
        axis.text.x.top = element_text(
          angle = 90, hjust = 0, vjust = 0.5, size = 12),
        strip.placement = "outside",
        axis.line = element_blank(),
        axis.ticks.y = element_blank(),
        aspect.ratio = nrow(gen_heatmap$ROC)/ncol(gen_heatmap$ROC))
  })



comm <- paste('/home/opt/R-3.4.0/bin/Rscript ../software/EpigeneticHeatmapMaker-from_input/epi_heatmap_from_file.R heatMapSamples', 
              '-c  software/EpigeneticHeatmapMaker-from_input/INSPIIRED.yml',
              '-o ./tmp/epiHeatMap',
              '-t epiCellTypes',
              '-f heatMapSites',
              '-r hg38')

epiGenomicHeatmap <- within(
  list(), {
    heatmap_sample_info <- read.csv('heatMapSamples')
    gen_heatmap <- readRDS('tmp/epiHeatMap/roc.res.rds')
    
    heatmap_scale <- seq(0.2, 0.8, 0.1)
    gen_heatmap_colors <- colorspace::diverge_hsv(
      length(heatmap_scale), h = c(240, 115), v = 1, power = 1)
    
    browser() 
    select_gen_features <- row.names(gen_heatmap$ROC)
    # select_gen_features <- c(
    #   "boundary.dist", "start.dist", "general.width", "gene.width",
    #   "within_refSeq_gene", "refSeq_counts.10k", "refSeq_counts.100k",
    #   "refSeq_counts.1M", "GC.100", "GC.1k", "GC.10k", "GC.100k", "GC.1M",
    #   "CpG_counts.1k", "CpG_counts.10k", "CpG_density.10k", "CpG_density.100k",
    #   "CpG_density.1M", "DNaseI_count.1k", "DNaseI_count.10k",
    #   "DNaseI_count.100k", "DNaseI_count.1M")
    
    gen_heatmap$ROC <- gen_heatmap$ROC[select_gen_features,]
    
    #heatmap_sample_levels <- c("PCR_amplicon", "Huh7_d0", "Huh7_d3", "Huh7_w3", "Humanized_mouse_m1", "WAS_d0", "CART_d0")
    heatmap_sample_levels <- c("GTSP1574", "GTSP1831", "GTSP2367", "GTSP2361")
    #heatmap_sample_levels <- c("B101", "N102", "N104", "U103")
    
    heatmap_figure_labels <- heatmap_sample_levels
    
    stat_cuts <- c(0, 0.001, 0.01, 0.05, 1)
    gen_comp_stats <- structure(cut(
      gen_heatmap$pvalues$op[select_gen_features, 1],
      stat_cuts,
      labels = c("***", " **", " * ", "   "),
      include.lowest = TRUE),
      names = select_gen_features)
    gen_row_names <- paste0(names(gen_comp_stats), " - ", gen_comp_stats)
    
    plot_data <- gen_heatmap$ROC %>%
      reshape2::melt() %>%
      mutate(
        feat = Var1,
        comp.sym = gen_comp_stats[Var1],
        Var1 = paste0(Var1, " - ", comp.sym),
        Var1 = factor(Var1, levels = gen_row_names),
        Var2 = factor(Var2, levels = heatmap_sample_levels),
        grp = " ",
        sig = as.vector(gen_heatmap$pvalues$np[select_gen_features,]),
        sym = cut(
          sig, stat_cuts, labels = c("***", " **", " * ", "   "),
          include.lowest = TRUE))
    
    levels(plot_data$Var2) <- heatmap_figure_labels
    
    #hot fix
    plot_data$Var2 <- heatMapSamples[match(plot_data$Var2, heatMapSamples$GTSP),]$patient
    
    ## browser()
    
    gen_plot <- ggplot(plot_data, aes(x = Var2, y = Var1, fill = value)) +
      geom_tile(color = 'black') +
      geom_text(aes(label = sym), color = "black", size = 3, nudge_y = -0.15) +
      scale_x_discrete(position = "top") +
      scale_fill_gradient2(
        breaks = c(0.2, 0.4, 0.6, 0.8),
        low = gen_heatmap_colors[1],
        mid = gen_heatmap_colors[round(length(heatmap_scale)/2)],
        high = gen_heatmap_colors[length(heatmap_scale)],
        midpoint = 0.5) +
      guides(fill = guide_colorbar(
        title.position = "left", title.hjust = 0.5,
        direction = "horizontal")) +
      labs(x = NULL, y = NULL, fill = "ROC\nScore") +
      # custom_theme +
      theme(
        legend.position = "bottom",
        axis.text.y = element_text( angle = 0, hjust = 0, vjust = 0.5, size = 12),
        axis.text.x.top = element_text(
          angle = 90, hjust = 0, vjust = 0.5, size = 12),
        strip.placement = "outside",
        axis.line = element_blank(),
        axis.ticks.y = element_blank(),
        aspect.ratio = nrow(gen_heatmap$ROC)/ncol(gen_heatmap$ROC))
  })




sitesFileName <- file.path('../data/tmp', sub('samples', 'sites', file))
write.table(df, quote = FALSE, row.names = FALSE, col.names = TRUE, file = sitesFileName, sep=',')
dir.create(file.path('../data/tmp', sampleTable$patient[1]), showWarnings = FALSE)

write(paste('echo', sampleTable$patient[1]), file='../data/tmp/commandScript.sh', append = TRUE)
write(comm, file='../data/tmp/commandScript.sh', append = TRUE)



# For each of the generated sample files, create a script that will create the corresponding intSite data files.
# The genomic heatmap maker appears to be limted to 20,000 inSites.

#
# Create command script to create genomic heatmaps.
#

write('# Command script', file='../data/tmp/commandScript.sh', append = FALSE)

invisible(lapply(list.files(path='../data/tmp', pattern='\\.samples.csv'), function(file){
  sampleTable <- read.csv(file.path('../data/tmp', file))
  o <- subset(d, patient==sampleTable$patient[1] & s %in% sampleTable$GTSP)
  df <- data.frame(o$seqnames, o$strand, o$start, o$s, 'hg38') 
  names(df) <- c('seqnames', 'strand', 'position', 'sampleName', 'refGenome')
  
  if(nrow(df) > 20000) df <- df[sample(1:nrow(df), 20000),]
  
  sitesFileName <- file.path('../data/tmp', sub('samples', 'sites', file))
  write.table(df, quote = FALSE, row.names = FALSE, col.names = TRUE, file = sitesFileName, sep=',')
  dir.create(file.path('../data/tmp', sampleTable$patient[1]), showWarnings = FALSE)
  comm <- paste(Rscript_path, './software/genomicHeatmapMaker-from_input/genomic_heatmap_from_file.R', 
                file.path('../data/tmp', file),
                '-c ./software/genomicHeatmapMaker-from_input/INSPIIRED.yml',
                '-o', file.path('../data/tmp', sampleTable$patient[1]),
                '-f', sitesFileName,
                '-r hg38')
  write(paste('echo', sampleTable$patient[1]), file='../data/tmp/commandScript.sh', append = TRUE)
  write(comm, file='../data/tmp/commandScript.sh', append = TRUE)
}))



