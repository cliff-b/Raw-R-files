library(dplyr)
library(tidyr)
library(extrafont)
library(viridis)
library(cowplot)

setwd("/Users/Cliff/Documents/Kosuri Lab/NGS files/170224 R100 Hiseq-2/")
load("R100_full_DF2.Rda")
R100_full_DF <- mutate(R100_full_DF, "nCount" = count/readsinlib*1000000)
load("GFPBCs.Rda")



R100zero <- filter(R100_full_DF, grepl("0", R100_full_DF$condition))
R100zero <- separate(R100zero, condition, into = c("template", "condition"))
R100four <- filter(R100_full_DF, grepl("4", R100_full_DF$condition))
R100four <- separate(R100four, condition, into = c("template", "condition"))

R100_4 <- full_join(R100zero, R100four, by = c("barcode", "template"))
R100_4$count.x[is.na(R100_4$count.x)] <- 0
R100_4$count.y[is.na(R100_4$count.y)] <- 0
R100_4$readsinlib.x[is.na(R100_4$readsinlib.x)] <- 0
R100_4$readsinlib.y[is.na(R100_4$readsinlib.y)] <- 0
R100_4$nCount.x[is.na(R100_4$nCount.x)] <- 0
R100_4$nCount.y[is.na(R100_4$nCount.y)] <- 0
R100_4$condition.x[is.na(R100_4$condition.x)] <- 0
R100_4$condition.y[is.na(R100_4$condition.y)] <- 4

R100RNA <- filter(R100_4, grepl("RNA", R100_4$template)) %>%
  select(template, barcode, "rcount0" = nCount.x, "rcount4" = nCount.y)
g.RNA <- ggplot(sample_n(R100RNA, 100000, replace = FALSE), aes(rcount0, rcount4)) + geom_point(alpha = 0.1) + labs(x = "Normalized RNA counts at 0 hours", y = "Normalized RNA counts at 4 hours", title = "170410 R100 normalized RNA over time") + geom_abline(slope = 1)+
  scale_x_log10(limits = c(0.01,1200)) + scale_y_log10(limits = c(0.01, 1200))
g.RNA

R100DNA <- filter(R100_4, grepl("DNA", R100_4$template)) %>%
  select(template, barcode, "dcount0" = nCount.x, "dcount4" = nCount.y )
g.DNA <- ggplot(sample_n(R100D5, 300000, replace = FALSE), aes(dcount0, dcount4)) + geom_point(alpha = 0.1) + labs(x = "Normalized DNA counts at 0 hours", y = "Normalized DNA counts at 4 hours", title = "170410 R100 normalized DNA over time") + geom_abline(slope = 1)+
  scale_x_log10(limits = c(0.01,1200)) + scale_y_log10(limits = c(0.01, 1200))
g.DNA
R100D5 <-filter(R100DNA, dcount0 > 0.2, dcount4 > 0.2)


R100RD <- right_join(R100RNA,R100D5, by = c("barcode")) %>%
  mutate("RNADNA0" = rcount0/dcount0*100) %>%
  mutate("RNADNA4" = rcount4/dcount4*100) %>%
  mutate("RNADNA0" = ifelse(is.na(RNADNA0),0,RNADNA0)) %>%
  mutate("RNADNA4" =  ifelse(is.na(RNADNA4),0,RNADNA4))
g.RD <- ggplot(filter(R100RD, YN == "P"), aes(RNADNA0, RNADNA4, color = YN)) + geom_point(alpha = 0.1) + labs(x = "Normalized RNA/DNA counts at 0 hours", y = "Normalized RNA/DNA counts at 4 hours", title = "170410 R100 normalized RNA/DNA over time") + geom_abline(slope = 1)+# scale_x_continuous(limits=c(0,2000)) + scale_y_continuous(limits =c(0,2000))
  scale_x_log10(limits = c(1,2000)) + scale_y_log10(limits = c(1, 2000))
g.RD

increase <- filter(R100RD, RNADNA4 > 10* RNADNA0)
g.incre <- ggplot(sample_n(increase, 1000, replace = FALSE), aes(RNADNA0, RNADNA4, color = YN)) + geom_point(alpha = 0.1) +#  scale_x_continuous(limits=c(0,2000)) + scale_y_continuous(limits =c(0,2000))
  scale_x_log10(limits = c(1,2000)) + scale_y_log10(limits = c(1, 2000)) + scale_color_brewer(palette = "Set1")
g.incre

mapBCs <- read.csv("170123 R100 mapped_barcodes.csv", skip = 1)
colnames(mapBCs) <- c("X_peptide","Y_peptide","reads","num_bcs", "barcode")
mapBCs$barcode <- as.character(mapBCs$barcode)
mapBCs <- mutate(mapBCs, barcode=strsplit(barcode, ",")) %>%
  unnest(barcode)
mapBCs <- mapBCs[!(duplicated(mapBCs$barcode) | duplicated(mapBCs$barcode, fromLast = TRUE)), ]
mapBCs <-  separate(mapBCs, X_peptide, into=c("X_peptide", "X_codon"), sep = "_rev_") %>%
  separate(Y_peptide, into=c("Y_peptide", "Y_codon"), sep = -4) 
mapBCs <- select(mapBCs, X_peptide, Y_peptide, barcode)

GFPBCs <- select(GFPBCs, "X_peptide" = plasmid, "barcode" = RevCompBC)
GFPBCs$Y_peptide <- GFPBCs$X_peptide
mapBCs <- bind_rows(GFPBCs, mapBCs)

R100RD <- left_join(R100RD, mapBCs, by = "barcode")

R100RD$YN[!is.na(R100RD$X_peptide)] <- "Y"
R100RD$YN[grepl("p6", R100RD$X_peptide)] <- "P"

R100Gene <- matrix(nrow = length(R100RD$RNADNA0),ncol = 2)#, ncol = 3)
R100Gene[,1] <- floor(R100RD$RNADNA0)
R100Gene[,2] <- floor(R100RD$RNADNA4)
R100Gene[,3] <- R100RD$YN
colnames(R100Gene) <- c("RNADNA0","RNADNA4")#,"YN")
rownames(R100Gene) <- R100RD$barcode
spikes <- rownames(R100Gene)[grepl("P", R100Gene[,3])]
R100Gene <- R100Gene[,1:2]
R100Gene <- mapply(R100Gene, FUN=as.numeric)

nR100 <- newSeqExpressionSet(as.matrix(R100Gene), phenoData = data.frame(x = c(0, 4), row.names = colnames(R100Gene)))

nR100 <- RUVg(nR100,spikes, k=1, center = FALSE)

dR100 <- as.data.frame(normCounts(nR100))
dR100$YN <- R100RD[,12]

g.RD <- ggplot(filter(dR100, YN == "Y" | YN == "P"), aes(RNADNA0, RNADNA4, color=YN)) + geom_point(alpha = 0.1) + labs(x = "0 hours RNA/DNA", y = "4 hours RNA/DNA", title = "170414 R100 RUV normalized RNA/DNA over time") + geom_abline(slope=1)# +
  scale_x_log10(limits = c(1, 800)) + scale_y_log10(limits = c(1, 800))
g.RD

sumR100RD <- filter(R100RD, YN == "Y" | YN == "P")
sumR100RD <- group_by(sumR100RD, X_peptide, Y_peptide, YN) %>%
  summarise(mean_RNA0 = mean(rcount0), mean_RNA4 = mean(rcount4), mean_DNA0 = mean(dcount0), mean_DNA4 = mean(dcount4), mean_RD0 = mean(RNADNA0), mean_RD4 = mean(RNADNA4), sd_RD0 = sd(RNADNA0), sd_RD4 = sd(RNADNA4), med_RD0 = median(RNADNA0), med_RD4 = median(RNADNA4), bcnum = n_distinct(barcode))

g.sum <- ggplot(filter(sumR100RD, YN == "P" | bcnum >5), aes(med_RD0, med_RD4, color = YN)) + geom_point(alpha = 0.2)  + labs(x = "Mean RNA/DNA at 0 hours", y = "Mean RNA/DNA at 4 hours", title = "170411 R100-2 mean RNA/DNA for mapped BCs")
g.sum

g.sum2 <- ggplot(sumR100RD, aes((med_RD0 + med_RD4)/2, bcnum)) + geom_point(alpha = 0.4)
g.sum2
