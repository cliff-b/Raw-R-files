library(dplyr)
library(tidyr)
library(extrafont)
library(viridis)
library(cowplot)
library(RUVSeq)

setwd("/Users/Cliff/Documents/Kosuri Lab/NGS files/161010 R20 Synzip timecourse/")
load("R20bigDF.Rda")


GFPBCs <-read.csv("160719 constitutive GFP Barcodes.csv", header = FALSE)
GFPBCs <- data.frame(GFPBCs, rep("GAATTCTGGGACC", nrow(GFPBCs)))
colnames(GFPBCs) <- c("plasmid", "barcode","filler")
GFPBCs$RevCompBC <- paste(GFPBCs$barcode, GFPBCs$filler, sep = "")

revComp <- function(x)
  chartr("ATGC","TACG", sapply(lapply(strsplit(x, NULL), rev), paste, collapse=""))


bigDF <- data.frame(bigDF, "RevCompBC" = revComp(as.character(bigDF$barcode)))
bigDF <- left_join(bigDF, GFPBCs, by = "RevCompBC")
bigDF <- filter(bigDF, !is.na(bigDF$plasmid))
bigDF <- mutate(bigDF, ncount = count/readsinlib*1000000)

R20D5RNA <- filter(bigDF, grepl("RNA", bigDF$condition)) %>%
  select(plasmid, condition, RevCompBC, "rcount" = ncount ) %>%
  separate(condition, into = c("template", "condition"))
R20D5DNA <- filter(bigDF, grepl("DNA", bigDF$condition)) %>%
  select(plasmid, condition, RevCompBC, "dcount" = ncount ) %>%
  separate(condition, into = c("template", "condition"))
R20D5RD <- inner_join(R20D5RNA,R20D5DNA, by = c("RevCompBC","plasmid","condition")) %>%
  mutate("RNADNA" = rcount/dcount*100) %>%
  separate("plasmid", into = c("plasmid", "clone"), sep = "-") 

spGFP <- select(R20D5RD, c(plasmid,clone,condition,RNADNA)) %>%
  spread(key = condition, value = RNADNA) %>%
  unite(col = plasmid, plasmid, clone, sep = "-") %>%
  filter(plasmid != "p63-4")

row.names(spGFP) <- spGFP$plasmid
spGFP <- select(spGFP, -plasmid)

names <- row.names(spGFP)[1:30]
PCAset <- newSeqExpressionSet(t(as.matrix(spGFP)))
PCAset3 <- newSeqExpressionSet(as.matrix(spGFP))
PCAset3 <- newSeqExpressionSet(as.matrix(floor(spGFP)), phenoData =  data.frame(rep("gene", length(spGFP)), row.names = colnames(spGFP)))
PCAset2 <- RUVg(PCAset3, names, k=2, isLog = TRUE)

g.GFPBC <- ggplot(R20D5RD, aes(plasmid,RNADNA, fill=plasmid)) + geom_bar(stat="identity") + facet_wrap(~condition)
g.GFPBC

g.PCA <- plotPCA()
g.PCA

RUVdPCA <- data.frame(counts(PCAset2)-counts(PCAset3))

spGFP$plasmid <- row.names(spGFP)
spGFP <- separate(spGFP, plasmid, into = c("plasmid","junk"), sep = "-")
spGFP <- select(spGFP, -junk)
prcobj <- prcomp(spGFP[2:6], scale = TRUE, center = TRUE)
autoplot(prcobj, data = spGFP, colour = "plasmid")

RUVdPCA$plasmid <- row.names(RUVdPCA)
RUVdPCA <- separate(RUVdPCA, plasmid, into = c("plasmid","junk"), sep = "-")
RUVdPCA <- select(RUVdPCA, -junk)
prcobj2 <- prcomp(RUVdPCA[1:5], scale = TRUE, center = TRUE)
autoplot(prcobj2, data = RUVdPCA, colour = "plasmid")
