library(dplyr)
library(tidyr)
library(extrafont)
library(viridis)
library(cowplot)

setwd("/Users/Cliff/Documents/Kosuri Lab/NGS files/170224 R100 Hiseq-2/")

R100_RNA_0 <- read.table("counts_cDNA_0h.txt", header = FALSE)
R100_RNA_0 <- data.frame(R100_RNA_0, rep("RNA-0",nrow(R100_RNA_0)),rep(sum(R100_RNA_0$V1), nrow(R100_RNA_0)))
colnames(R100_RNA_0) <- c("count", "barcode", "condition","readsinlib")

R100_DNA_0 <- read.table("counts_pDNA_0h.txt", header = FALSE)
R100_DNA_0 <- data.frame(R100_DNA_0, rep("DNA-0",nrow(R100_DNA_0)),rep(sum(R100_DNA_0$V1), nrow(R100_DNA_0)))
colnames(R100_DNA_0) <- c("count", "barcode","condition","readsinlib")

R100_RNA_4 <- read.table("counts_cDNA_4h.txt", header = FALSE)
R100_RNA_4 <- data.frame(R100_RNA_4, rep("RNA-4",nrow(R100_RNA_4)),rep(sum(R100_RNA_4$V1), nrow(R100_RNA_4)))
colnames(R100_RNA_4) <- c("count", "barcode", "condition","readsinlib")

R100_DNA_4 <-read.table("counts_pDNA_4h.txt", header = FALSE)
R100_DNA_4 <- data.frame(R100_DNA_4, rep("DNA-4",nrow(R100_DNA_4)),rep(sum(R100_DNA_4$V1), nrow(R100_DNA_4)))
colnames(R100_DNA_4) <- c("count", "barcode","condition","readsinlib")


mapBCs <- read.csv("170123 R100 mapped_barcodes.csv", skip = 1)
colnames(mapBCs) <- c("X_peptide","Y_peptide","reads","num_bcs", "barcode")
mapBCs$barcode <- as.character(mapBCs$barcode)
mapBCs <- mutate(mapBCs, barcode=strsplit(barcode, ",")) %>%
  unnest(barcode)
mapBCs <- mapBCs[!(duplicated(mapBCs$barcode) | duplicated(mapBCs$barcode, fromLast = TRUE)), ]
mapBCs <-  separate(mapBCs, X_peptide, into=c("X_peptide", "X_codon"), sep = "_rev_") %>%
  separate(Y_peptide, into=c("Y_peptide", "Y_codon"), sep = -4) 
mapBCs <- select(mapBCs, X_peptide, Y_peptide, barcode)

R100_full_DF <- bind_rows(R100_RNA_0, R100_DNA_0, R100_RNA_4, R100_DNA_4)
remove(R100_RNA_0, R100_DNA_0, R100_RNA_4, R100_DNA_4)

load("R100_full_DF2.Rda")
load("/Users/Cliff/Documents/Kosuri Lab/NGS files/170208 R100 Hiseq-1/R100_full_DF.Rda")
load("GFPBCs.Rda")
GFPBCs <- select(GFPBCs, plasmid, "barcode" = RevCompBC)
mapGFP <- filter(left_join(R100_full_DF,GFPBCs, by = "barcode"),!is.na(plasmid)) %>%
  filter(plasmid != "p63-4") %>%
  separate(plasmid, into = c("plasmid", "clone"),sep = "-")
  
mapGFP <- mutate(mapGFP, "nCount" = count/readsinlib*1000000)

GFPRNA <- filter(mapGFP, grepl("RNA", mapGFP$condition)) %>%
  select(plasmid, clone, condition, barcode, "rcount" = nCount ) %>%
  separate(condition, into = c("template", "condition"))
GFPDNA <- filter(mapGFP, grepl("DNA", mapGFP$condition)) %>%
  select(plasmid, clone, condition, barcode, "dcount" = nCount ) %>%
  separate(condition, into = c("template", "condition"))
GFPRD <- inner_join(GFPRNA,GFPDNA, by = c("barcode","condition","plasmid","clone")) %>%
  mutate("RNAoDNA" = rcount/dcount*100) 
sumGFP <- group_by(GFPRD, plasmid, condition) %>%
  summarise(mean_rcount =  mean(rcount), mean_dcount = mean(dcount), mean_RNADNA = mean(RNAoDNA), sd_RNADNA=sd(RNAoDNA))
gfpERR <- aes(ymax=mean_RNADNA + sd_RNADNA, ymin=mean_RNADNA - sd_RNADNA)
g.sumGFP <- ggplot(sumGFP, aes(plasmid, mean_RNADNA, fill= plasmid))  +geom_bar(stat="identity")+ facet_wrap(~condition, labeller = labeller(condition=induction_labels), scales = "free") + geom_errorbar(gfpERR) +  #geom_point(mapping=aes(plasmid,RNAoDNA), data=GFPRD) +
  labs(title="170320 R100 GFP constructs Average RNA/DNA counts", x ="plasmid", y = "Mean RNA/DNA") + scale_y_continuous(limits=c(0,900)) + theme(text = element_text(family="Calibri"), strip.background = element_rect(fill="white"))
g.sumGFP



R100_full_DF <- mutate(R100_full_DF, "nCount" = count/readsinlib*1000000)
R100_full_DF <- filter(R100_full_DF, (count > 10) | grepl("RNA", R100_full_DF$condition))
R100_full_DF <- left_join(R100_full_DF, GFPBCs, by = "barcode")
R100_full_DF <- left_join(R100_full_DF, mapBCs, by = "barcode")
R100_full_DF3 <- filter(R100_full_DF, X_peptide != "<NA>", Y_peptide != "<NA>")
R100_full_DF3 <- filter(R100_full_DF3, condition != "RNA-8", condition != "DNA-8")

R100_4 <- full_join(R100_full_DF3, R100_full_DF2, by = c("barcode", "condition", "X_peptide", "Y_peptide"))
R100_4$count.x[is.na(R100_4$count.x)] <- 0
R100_4$count.y[is.na(R100_4$count.y)] <- 0
R100_4$readsinlib.x[is.na(R100_4$readsinlib.x)] <- 0
R100_4$readsinlib.y[is.na(R100_4$readsinlib.y)] <- 0
R100_4 <- mutate(R100_4, fcount = count.x + count.y, freads = readsinlib.x + readsinlib.y)
R100_4 <- mutate(R100_4, fnCount = fcount/freads*1000000)

R100RNA <- filter(R100_4, grepl("RNA", R100_4$condition)) %>%
  select(X_peptide, Y_peptide, condition, barcode, fcount, freads, "rcount" = fnCount ) %>%
  separate(condition, into = c("template", "condition"))
sumRNA <- group_by(R100RNA, X_peptide, Y_peptide, condition) %>%
  summarise(mean_RNA = mean(rcount), sd_R = sd(rcount), bcnum = n_distinct(barcode), tot_RNA = sum(fcount))
  
R100DNA <- filter(R100_4, grepl("DNA", R100_4$condition)) %>%
  select(X_peptide, Y_peptide, condition, barcode, fcount, "dcount" = fnCount ) %>%
  separate(condition, into = c("template", "condition"))
sumDNA <- group_by(R100DNA, X_peptide, Y_peptide, condition) %>%
  summarise(mean_DNA = mean(dcount), sd_D = sd(dcount), bcnum_D = n_distinct(barcode), tot_DNA = sum(fcount))
R100RD <- right_join(R100RNA,R100DNA, by = c("X_peptide","Y_peptide","barcode","condition")) %>%
  mutate("RNADNA" = rcount/dcount*100) %>%
  select(-template.x, -template.y) %>%
  mutate("RNADNA" = ifelse(is.na(RNADNA),0,RNADNA)) %>%
  mutate("rcount" =  ifelse(is.na(rcount),0,rcount))

######################################## Summary stats and graphs ######################################## 

sumR100RD <- group_by(R100RD, X_peptide, Y_peptide, condition) %>%
  summarise(mean_RNA = mean(rcount), sd_R = sd(rcount), mean_DNA = mean(dcount), sd_D = sd(dcount), mean_RD = mean(RNADNA), sd_RD = sd(RNADNA), med_RD = median(RNADNA), bcnum = n_distinct(barcode), SEM = sd(RNADNA)/sqrt(n()))
sumR100RD <- filter(sumR100RD, bcnum > 10)

induction_labels <- c("0" = "0 hours induction RNA/DNA", "4" = "4 hours induction RNA/DNA", "8" = "8 hour induction RNA/DNA")
induction_labels2 <- c("0" = "0 hours induction DNA", "4" = "4 hours induction DNA", "8" = "8 hour induction DNA")
induction_labels3 <- c("0" = "0 hours induction RNA", "4" = "4 hours induction RNA", "8" = "8 hour induction RNA")

g.sumR100RD <- ggplot(sumR100RD, aes(X_peptide, Y_peptide)) + geom_tile(aes(fill=SEM)) + facet_wrap(~condition, labeller = labeller(condition=induction_labels), scales = "free") +
  theme(text=element_text(family="Calibri"), strip.text=element_text(size = 10), strip.background=element_rect(fill = "White"), axis.text.x = element_text(angle=90,size =5),
        axis.text.y = element_text(size = 5), legend.title=element_text(size = 8), legend.title.align=0.5, legend.text=element_text(size =8)) +
  scale_fill_viridis(option="viridis", name = "Mean \nRNA to DNA  ", na.value="white", limits=c(0, 770)) + labs(title="Roman's 100 RNA/DNA barcode ratio", x="X peptide", y = "Y peptide")
g.sumR100RD

g.sumR100RDmed <- ggplot(sumR100RD, aes(X_peptide, Y_peptide)) + geom_tile(aes(fill=med_RD)) + facet_wrap(~condition, labeller = labeller(condition=induction_labels), scales = "free") +
  theme(text=element_text(family="Calibri"), strip.text=element_text(size = 10), strip.background=element_rect(fill = "White"), axis.text.x = element_text(angle=90, hjust=1,size =6),
        axis.text.y = element_text(size = 6), legend.title=element_text(size = 8), legend.title.align=0.5, legend.text=element_text(size =8)) +
  scale_fill_viridis(option="viridis", name = "Median \nRNA to DNA ", na.value="white") + labs(title="170214 Roman's 100 RNA/DNA barcode ratio (median)", x="X Peptide", y = "Y Peptide")
g.sumR100RDmed

g.sumR100D <- ggplot(sumR100RD, aes(X_peptide, Y_peptide)) + geom_tile(aes(fill=log(mean_DNA))) + facet_wrap(~condition, labeller = labeller(condition=induction_labels2), scales = "free") +
  theme(text=element_text(family="Calibri"), strip.text=element_text(size = 10), strip.background=element_rect(fill = "White"), axis.text.x = element_text(angle=90,size =6),
        axis.text.y = element_text(size = 6), legend.title=element_text(size = 8), legend.title.align=0.5, legend.text=element_text(size =8)) +
  scale_fill_viridis(option="viridis", name = "Log Mean\n DNA ", na.value="white") + labs(title="170214 Roman's 100 DNA counts", x="X Peptide", y = "Y Peptide")
g.sumR100D

g.sumR20RDsd <- ggplot(sumR100RD, aes(X_peptide, Y_peptide)) + geom_tile(aes(fill=sd_RD)) + facet_wrap(~condition, labeller = labeller(condition=induction_labels), scales = "free") +
  theme(text=element_text(family="Calibri"), strip.text=element_text(size = 10), strip.background=element_rect(fill = "White"), axis.text.x = element_text(angle=90,size =6),
        axis.text.y = element_text(size = 6), legend.title=element_text(size = 8), legend.title.align=0.5, legend.text=element_text(size =8)) +
  scale_fill_viridis(option="viridis", name = "Standard deviation \n of RNA to DNA ", na.value="white") + labs(title="170214 Roman's 100 RNA/DNA standard deviation of barcode ratio", x="X peptide", y = "Y peptide")
g.sumR20RDsd


g.sumR20RDbc <- ggplot(sumR100RD, aes(X_peptide, Y_peptide)) + geom_tile(aes(fill=bcnum)) + facet_wrap(~condition, labeller = labeller(condition=induction_labels), scales = "free") +
  theme(text=element_text(family="Calibri"), strip.text=element_text(size = 10), strip.background=element_rect(fill = "White"), axis.text.x = element_text(angle=90, size =8),
        axis.text.y = element_text(size = 6), legend.title=element_text(size = 8), legend.title.align=0.5, legend.text=element_text(size =8)) +
  scale_fill_viridis(option="viridis", name = "Number of\n barcodes ", na.value="white") + labs(title="170214 Roman's 100 number of barcodes per construct", x="X peptide", y = "Y peptide")
g.sumR20RDbc

R100RD <- mutate(R100RD, "combo" = paste(X_peptide,Y_peptide))
g.RDbyCon <- ggplot(filter(R100RD, X_peptide==">P4Q2950"), aes(Y_peptide, RNADNA, color=condition)) + geom_boxplot(outlier.shape=NA, color="red") + geom_jitter(alpha=0.4) +
  theme(axis.text.x = element_text(angle=90,size =6), strip.background=element_rect(fill = "White")) + facet_wrap(~condition, nrow=3, labeller = labeller(condition=induction_labels), scales = "free") +
  labs(title="170214 Roman's 100 >P4Q2950 RNA/DNA by barcode", x=">P4Q2950's partner", y = "RNA/DNA (AU)")
g.RDbyCon
