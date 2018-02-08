library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(viridis)
library(extrafont)

setwd("/Users/Cliff/Documents/Kosuri Lab/NGS files/160831 R20 Synzip run 1")

#read in a horrific amount of data and label it with library and number of reads from that library
R20_RNA_0 <- read.table("R20_RNA_0.txt", header = FALSE)
R20_RNA_0 <- data.frame(R20_RNA_0, rep("RNA-0",nrow(R20_RNA_0)),rep(sum(R20_RNA_0$V1), nrow(R20_RNA_0)))
colnames(R20_RNA_0) <- c("count", "barcode", "condition","readsinlib")

R20_DNA_0 <- read.table("R20_DNA_0.txt", header = FALSE)
R20_DNA_0 <- data.frame(R20_DNA_0, rep("DNA-0",nrow(R20_DNA_0)),rep(sum(R20_DNA_0$V1), nrow(R20_DNA_0)))
colnames(R20_DNA_0) <- c("count", "barcode","condition","readsinlib")

R20_RNA_a <- read.table("R20_RNA_A.txt", header = FALSE)
R20_RNA_a <- data.frame(R20_RNA_a, rep("RNA-A",nrow(R20_RNA_a)),rep(sum(R20_RNA_a$V1), nrow(R20_RNA_a)))
colnames(R20_RNA_a) <- c("count", "barcode", "condition","readsinlib")

R20_DNA_a <-read.table("R20_DNA_A.txt", header = FALSE)
R20_DNA_a <- data.frame(R20_DNA_a, rep("DNA-A",nrow(R20_DNA_a)),rep(sum(R20_DNA_a$V1), nrow(R20_DNA_a)))
colnames(R20_DNA_a) <- c("count", "barcode","condition","readsinlib")

R20_RNA_b <- read.table("R20_RNA_B.txt", header = FALSE)
R20_RNA_b <- data.frame(R20_RNA_b, rep("RNA-B",nrow(R20_RNA_b)),rep(sum(R20_RNA_b$V1), nrow(R20_RNA_b)))
colnames(R20_RNA_b) <- c("count", "barcode", "condition","readsinlib")

R20_DNA_b <- read.table("R20_DNA_B.txt", header = FALSE)
R20_DNA_b <- data.frame(R20_DNA_b, rep("DNA-B",nrow(R20_DNA_b)),rep(sum(R20_DNA_b$V1), nrow(R20_DNA_b)))
colnames(R20_DNA_b) <- c("count", "barcode","condition","readsinlib")

R20_RNA_c <- read.table("R20_RNA_C.txt", header = FALSE)
R20_RNA_c <- data.frame(R20_RNA_c, rep("RNA-C",nrow(R20_RNA_c)),rep(sum(R20_RNA_c$V1), nrow(R20_RNA_c)))
colnames(R20_RNA_c) <- c("count", "barcode", "condition","readsinlib")

R20_DNA_c <-read.table("R20_DNA_A.txt", header = FALSE)
R20_DNA_c <- data.frame(R20_DNA_c, rep("DNA-C",nrow(R20_DNA_c)),rep(sum(R20_DNA_c$V1), nrow(R20_DNA_c)))
colnames(R20_DNA_c) <- c("count", "barcode","condition","readsinlib")


GFPBCs <-read.csv("160719 constitutive GFP Barcodes.csv", header = FALSE)
GFPBCs <- data.frame(GFPBCs, rep("GAATTCTGGGACC", nrow(GFPBCs)))
colnames(GFPBCs) <- c("plasmid", "barcode","filler")
GFPBCs$RevCompBC <- paste(GFPBCs$barcode, GFPBCs$filler, sep = "")

################################get the mapped barcodes out into single BCs per column and remove duplicates##########
mapBCs <- read.csv("roman_counts.csv", skip = 1)
colnames(mapBCs) <- c("X_peptide","Y_peptide","reads","num_bcs", "RevCompBC")
mapBCs$RevCompBC <- as.character(mapBCs$RevCompBC)
mapBCs <- mutate(mapBCs, RevCompBC=strsplit(RevCompBC, ",")) %>%
  unnest(RevCompBC)
mapBCs <- mapBCs[!(duplicated(mapBCs$RevCompBC) | duplicated(mapBCs$RevCompBC, fromLast = TRUE)), ]


#make the big one and get rid of the duplicates
bigDF <- bind_rows(R20_RNA_0, R20_RNA_a, R20_RNA_b, R20_RNA_c, R20_DNA_0, R20_DNA_a, R20_DNA_b, R20_DNA_c)
remove(R20_RNA_0, R20_RNA_a, R20_RNA_b, R20_RNA_c, R20_DNA_0, R20_DNA_a, R20_DNA_b, R20_DNA_c)

#homemade ghettoass reverse complement function
revComp <- function(x)
  chartr("ATGC","TACG", sapply(lapply(strsplit(x, NULL), rev), paste, collapse=""))


bigDF <- data.frame(bigDF, "RevCompBC" = revComp(as.character(bigDF$barcode)))
bigDF <- mutate(bigDF, "nCount" = count/readsinlib*1000000)
DNA5DF <- filter(bigDF, (count > 10) | grepl("RNA", bigDF$condition))

D5Map <- filter(left_join(DNA5DF, mapBCs, by = "RevCompBC"), !is.na(reads))
  
R20D5RNA <- filter(D5Map, grepl("RNA", D5Map$condition)) %>%
  select(X_peptide, Y_peptide, condition, RevCompBC, "rcount" = nCount ) %>%
  separate(condition, into = c("template", "condition"))
R20D5DNA <- filter(D5Map, grepl("DNA", D5Map$condition)) %>%
  select(X_peptide, Y_peptide, condition, RevCompBC, "dcount" = nCount ) %>%
  separate(condition, into = c("template", "condition"))
R20D5RD <- inner_join(R20D5RNA,R20D5DNA, by = c("RevCompBC","X_peptide","Y_peptide","condition")) %>%
  mutate("RNADNA" = rcount/dcount*100) %>%
  separate("Y_peptide", into = c("Y_peptide", "Y_codon"), sep = "_") %>%
  separate ("X_peptide", into = c("X_peptide", "X_codon"), sep = "_")

R20D5RD$X_peptide <- as.factor(R20D5RD$X_peptide)
R20D5RD$Y_peptide <- as.factor(R20D5RD$Y_peptide)
R20D5RD$X_peptide <- factor(R20D5RD$X_peptide, levels = c(">GCN4", ">P1", ">P2",">P3",">P4",">P5",">P6",">P7",">P8",">P9",">P10",">P11",">P12",">P3mS",">P4mS",">P5mS",">P6mS",">P7mS",">P8mS",">P5SC1",">P5SC2",">P6SC1",">P6SC2"))
R20D5RD$Y_peptide <- factor(R20D5RD$Y_peptide, levels = rev(c(">GCN4", ">P1", ">P2",">P3",">P4",">P5",">P6",">P7",">P8",">P9",">P10",">P11",">P12",">P3mS",">P4mS",">P5mS",">P6mS",">P7mS",">P8mS",">P5SC1",">P5SC2",">P6SC1",">P6SC2")))


sumR20D5 <- group_by(R20D5RD, X_peptide, Y_peptide, condition) %>%
  summarise(mean_RNA = mean(rcount), sd_R = sd(rcount), mean_DNA = mean(dcount), sd_D = sd(dcount), mean_RD = mean(RNADNA), sd_RD = sd(RNADNA), med_RD = median(RNADNA), bcnum = n_distinct(RevCompBC))# %>%
  #to make graphs of subsets
  #filter(as.integer(X_peptide) >3 & as.integer(X_peptide) <14) %>%
  #filter(as.integer(Y_peptide) >10 & as.integer(Y_peptide) <22)
sumR20D5<-data.frame(sumR20D5,"recip" = rep(1.2, nrow(sumR20D5)))
for(i in 1:nrow(sumR20D5)){
  if(as.integer(sumR20D5$X_peptide[i]) <= 24 - as.integer(sumR20D5$Y_peptide[i])) {
  a<- filter(sumR20D5, X_peptide==sumR20D5[i,2],Y_peptide == sumR20D5[i,1], condition== sumR20D5[i,3])
  sumR20D5[i,12]<-a$med_RD
  }
}

R20D5RD <- mutate(R20D5RD, "combo" = paste(X_peptide,Y_peptide))

g.RDbyCon <- ggplot(filter(R20D5RD, X_peptide==">P4" & condition=="A"), aes(combo, RNADNA)) + geom_boxplot(outlier.shape=NA, color="red") + geom_jitter(alpha=0.4) +theme(axis.text.x = element_text(angle=90, hjust=1,size =6))
g.RDbyCon

g.RDbyCon2 <- ggplot(filter(R20D5RD), aes(combo, RNADNA)) + geom_boxplot(outlier.shape=NA, color="red") + geom_jitter(alpha=0.4) +theme(axis.text.x = element_text(angle=90, hjust=1,size =6))
g.RDbyCon2 +facet_grid(condition~.)

g.FCbyCon <- ggplot(filter(test2, X_peptide.x==">P4"), aes(combo.x, bRD)) + geom_jitter(alpha=0.4) +theme(axis.text.x = element_text(angle=90, hjust=1,size =6))
g.FCbyCon

g.RDbyCon2 <- ggplot(filter(R20D5RD), aes(combo, RNADNA)) + geom_boxplot(outlier.shape=NA, color="red") + geom_jitter(alpha=0.4) +theme(axis.text.x = element_text(angle=90, hjust=1,size =6))
g.RDbyCon2 +facet_grid(condition~.)



g.sumR20RDrecip <- ggplot(filter(sumR20D5, sumR20D5$recip != 1.2 & X_peptide != ">GCN4" & Y_peptide !=">GCN4"), aes(mean_RD, sd_RD)) + geom_point() +geom_smooth(method="lm", se=FALSE) +
  labs(title="160908 Roman's 20 RNA/DNA reciprocal counts", x="XY mean RNA/DNA counts", y = "XY sd RNA/DNA counts") + scale_color_manual()
g.sumR20RDrecip

g.sumR20RDrect <- ggplot(filter(sumR20D5, sumR20D5$recip != 1.2), aes(X_peptide, Y_peptide)) + geom_tile(aes(fill=mean_RD - recip))+ facet_wrap(~condition, labeller = labeller(condition=induction_labels2), scales = "free") +
  theme(text=element_text(family="Calibri"), strip.text=element_text(size = 10), strip.background=element_rect(fill = "White"), axis.text.x = element_text(angle=90, hjust=1,size =6),
        axis.text.y = element_text(size = 6), legend.title=element_text(size = 8), legend.title.align=0.5, legend.text=element_text(size =8)) +
  scale_fill_gradient2(low="Red", high = "Blue", mid ="gray69", name = "XY RNA/DNA - \nYX RNA/DNA ") + labs(title="160908 Roman's 20 XY RNA/DNA ratio minus YX RNA/DNA ratio", x="X peptide", y = "Y peptide")
g.sumR20RDrect

g.sumR20D5cor <- ggplot(filter(sumR20D5, mean_RD > 150 & recip > 150 | condition == "0"), aes(X_peptide, Y_peptide)) +geom_tile(aes(fill=mean_RD))+ facet_wrap(~condition, labeller = labeller(condition=induction_labels2), scales = "free",nrow = 2, ncol = 2) +
  theme(text=element_text(family="Calibri"), strip.text=element_text(size = 10), strip.background=element_rect(fill = "White"), axis.text.x = element_text(angle=90, hjust=1,size =6),
        axis.text.y = element_text(size = 6), legend.title=element_text(size = 8), legend.title.align=0.5, legend.text=element_text(size =8)) +
  scale_fill_viridis(option="viridis", name = "Both XY and\nYX >150 ", na.value="white") + labs(title="160908 Roman's 20 RNA/DNA of XY and YX >150 ", x="X peptide", y = "Y peptide")
g.sumR20D5cor

g.sumR20D5RD <- ggplot(sumR20D5, aes(X_peptide, Y_peptide)) + geom_tile(aes(fill=mean_RD)) + facet_wrap(~condition, labeller = labeller(condition=induction_labels2), scales = "free") +
  theme(text=element_text(family="Calibri"), strip.text=element_text(size = 10), strip.background=element_rect(fill = "White"), axis.text.x = element_text(angle=90, hjust=1,size =6),
        axis.text.y = element_text(size = 6), legend.title=element_text(size = 8), legend.title.align=0.5, legend.text=element_text(size =8)) +
  scale_fill_viridis(option="viridis", name = "Normalized \nRNA to DNA ", na.value="white") + labs(title="160908 Roman's 20 (DNA>10) RNA/DNA barcode ratio", x="X peptide", y = "Y peptide")
g.sumR20D5RD

g.sumR20D5sd <- ggplot(sumR20D5, aes(X_peptide, Y_peptide)) + geom_tile(aes(fill=sd_RD)) + facet_wrap(~condition, labeller = labeller(condition=induction_labels2), scales = "free") +
  theme(text=element_text(family="Calibri"), strip.text=element_text(size = 10), strip.background=element_rect(fill = "White"), axis.text.x = element_text(angle=90, hjust=1,size =6),
        axis.text.y = element_text(size = 6), legend.title=element_text(size = 8), legend.title.align=0.5, legend.text=element_text(size =8)) +
  scale_fill_viridis(option="inferno", name = "Standard \n Deviation of \nRNA to DNA", na.value="green") + labs(title="160906 Roman's 20 Standard Deviation (DNA >10) of RNA/DNA barcode ratio", x="X peptide", y = "Y peptide")
g.sumR20D5sd

g.sumR20D5D <- ggplot(sumR20D5, aes(X_peptide, Y_peptide)) + geom_tile(aes(fill=log10(mean_DNA))) + facet_wrap(~condition, labeller = labeller(condition=induction_labels3), scales = "free") +
  theme(text=element_text(family="Comic Sans MS"), strip.text=element_text(size = 10), strip.background=element_rect(fill = "White"), axis.text.x = element_text(angle=90, hjust=1,size =6),
        axis.text.y = element_text(size = 6), legend.title=element_text(size = 8), legend.title.align=0.5, legend.text=element_text(size =8), plot.title=element_text(family="Jokerman")) +
  scale_fill_viridis(option="inferno", name = "Normalized \n DNA", na.value="white") + labs(title="160906 Roman's 20 (DNA>5) DNA barcodes count per million", x="X peptide", y = "Y peptide")
g.sumR20D5D

g.sumR20D5D <- ggplot(sumR20D5, aes(X_peptide, Y_peptide)) + geom_tile(aes(fill=log10(sd_D))) + facet_wrap(~condition, labeller = labeller(condition=induction_labels3), scales = "free") +
  theme(text=element_text(family="Comic Sans MS"), strip.text=element_text(size = 10), strip.background=element_rect(fill = "White"), axis.text.x = element_text(angle=90, hjust=1,size =6),
        axis.text.y = element_text(size = 6), legend.title=element_text(size = 8), legend.title.align=0.5, legend.text=element_text(size =8), plot.title=element_text(family="Jokerman")) +
  scale_fill_viridis(option="inferno", name = "Normalized \n DNA \nLog(10)", na.value="white") + labs(title="160906 Roman's 20 (DNA>5) Standard Deviation of DNA barcodes count per million", x="X peptide", y = "Y peptide")
g.sumR20D5D

g.sumR20D5R <- ggplot(sumR20D5, aes(X_peptide, Y_peptide)) + geom_tile(aes(fill=log(mean_RNA))) + facet_wrap(~condition, labeller = labeller(condition=induction_labels3), scales = "free") +
  theme(text=element_text(family="Comic Sans MS"), strip.text=element_text(size = 10), strip.background=element_rect(fill = "White"), axis.text.x = element_text(angle=90, hjust=1,size =6),
        axis.text.y = element_text(size = 6), legend.title=element_text(size = 8), legend.title.align=0.5, legend.text=element_text(size =8), plot.title=element_text(family="Jokerman")) +
  scale_fill_viridis(option="inferno", name = "Normalized \n DNA", na.value="white") + labs(title="160906 Roman's 20 (Counts >5) RNA barcodes count per million", x="X peptide", y = "Y peptide")
g.sumR20D5R

########################GFP manipulations###################################
mapGFP <- filter(left_join(bigDF,GFPBCs, by = "RevCompBC"),!is.na(plasmid))
mapGFP <- separate(mapGFP,plasmid, into = c("plasmid", "clone"),sep = "-")

GFPRNA <- filter(mapGFP, grepl("RNA", mapGFP$condition)) %>%
  select(plasmid, clone, condition, RevCompBC, "rcount" = nCount ) %>%
  separate(condition, into = c("template", "condition"))
GFPDNA <- filter(mapGFP, grepl("DNA", mapGFP$condition)) %>%
  select(plasmid, clone, condition, RevCompBC, "dcount" = nCount ) %>%
  separate(condition, into = c("template", "condition"))
GFPRD <- inner_join(GFPRNA,GFPDNA, by = c("RevCompBC","condition","plasmid","clone")) %>%
  mutate("RNAoDNA" = rcount/dcount*100) %>%
  # get rid of p63-4 which seems to be mismatched
  filter(RevCompBC != "GTTACCAGAATTCTGGGACC")
sumGFP <- group_by(GFPRD, plasmid, clone) %>%
  summarise(mean_rcount =  mean(rcount), mean_dcount = mean(dcount), mean_RNADNA = mean(RNAoDNA), sd_RNADNA=sd(RNAoDNA))




R20RNA <- filter(fullMap, grepl("RNA", fullMap$condition)) %>%
  select(X_peptide, Y_peptide, condition, RevCompBC, "rcount" = nCount ) %>%
  separate(condition, into = c("template", "condition"))
R20DNA <- filter(fullMap, grepl("DNA", fullMap$condition)) %>%
  select(X_peptide, Y_peptide, condition, RevCompBC, "dcount" = nCount ) %>%
  separate(condition, into = c("template", "condition"))
R20RD <- inner_join(R20RNA,R20DNA, by = c("RevCompBC","X_peptide","Y_peptide","condition")) %>%
  mutate("RNADNA" = rcount/dcount*100) %>%
  separate("Y_peptide", into = c("Y_peptide", "Y_codon"), sep = "_") %>%
  separate ("X_peptide", into = c("X_peptide", "X_codon"), sep = "_")

#there's a handfull of #rev barcodes that appear in the Y peptide. This needs more investigation, but for now we'll drop it to make it look better.
R20RD <- filter(R20RD, !grepl("rev",R20RD$Y_codon))
R20RD$X_peptide <- as.factor(R20RD$X_peptide)
R20RD$Y_peptide <- as.factor(R20RD$Y_peptide)
R20RD$X_peptide <- factor(R20RD$X_peptide, levels = c(">GCN4", ">P1", ">P2",">P3",">P4",">P5",">P6",">P7",">P8",">P9",">P10",">P11",">P12",">P3mS",">P4mS",">P5mS",">P6mS",">P7mS",">P8mS",">P5SC1",">P5SC2",">P6SC1",">P6SC2"))
R20RD$Y_peptide <- factor(R20RD$Y_peptide, levels = rev(c(">GCN4", ">P1", ">P2",">P3",">P4",">P5",">P6",">P7",">P8",">P9",">P10",">P11",">P12",">P3mS",">P4mS",">P5mS",">P6mS",">P7mS",">P8mS",">P5SC1",">P5SC2",">P6SC1",">P6SC2")))



sumR20 <- group_by(R20RD, X_peptide, Y_peptide, condition) %>%
  summarise(mean_RNA = mean(rcount), mean_DNA = mean(dcount), mean_RD = mean(RNADNA), sd_RD = sd(RNADNA))

g.sumR20RD <- ggplot(sumR20, aes(X_peptide, Y_peptide)) + geom_tile(aes(fill=mean_RD)) + facet_wrap(~condition, labeller = labeller(condition=induction_labels2), scales = "free") +
  theme(text=element_text(family="Comic Sans MS"), strip.text=element_text(size = 10), strip.background=element_rect(fill = "White"), axis.text.x = element_text(angle=90, hjust=1,size =6),
        axis.text.y = element_text(size = 6), legend.title=element_text(size = 8), legend.title.align=0.5, legend.text=element_text(size =8), plot.title=element_text(family="Jokerman")) +
  scale_fill_viridis(option="inferno", name = "Normalized \nRNA to DNA \n(Log10)", na.value="white") + labs(title="160906 Roman's 20 RNA/DNA barcode ratio", x="X peptide", y = "Y peptide")
g.sumR20RD

R20RD$condition <- as.factor(R20RD$condition)



g.R20RD <- ggplot(R20RD, aes(X_peptide, Y_peptide)) + geom_tile(aes(fill=log10(RNADNA))) + facet_wrap(~condition)                     

g.R20RD
#bigDF <- mutate(bigDF, "ncount" = bigDF$count/)

induction_labels2 <- c("0" = "Uninduced RNA/DNA", "A" = "5uM DAPG 5ng/mL ATc RNA/DNA", "B" = "7.5uM DAPG 10ng/mL ATc RNA/DNA", "C" = "15uM DAPG 40ng/mL ATc RNA/DNA")
induction_labels3 <- c("0" = "Uninduced DNA", "A" = "5uM DAPG 5ng/mL ATc DNA", "B" = "7.5uM DAPG 10ng/mL ATc DNA", "C" = "15uM DAPG 40ng/mL ATc DNA")

g.GFPBCs <- ggplot(mapGFP, aes(count)) + geom_histogram(bins = 50, fill= "green3", color="black") + scale_x_log10() + labs(title="160901 R20 number of counts per GFP barcode", x ="Barcode counts", y = "# of Counts") +
  facet_wrap(~condition, labeller = labeller(condition=induction_labels)) + theme(text = element_text(family="Calibri"), strip.background = element_rect(fill="white"))
g.GFPBCs

g.GFPBCsbyp <- ggplot(mapGFP, aes(nCount, fill=plasmid)) + geom_histogram(bins = 50) + scale_x_continuous(limits=c(0,370)) + scale_y_continuous(limits = c(0,15))+ labs(title="160901 R20 number of normalized counts per GFP barcode", x ="Barcode counts", y = "# of Counts") +
  facet_wrap(~condition, labeller = labeller(condition=induction_labels), scales = "free") + theme(text = element_text(family="Calibri"), strip.background = element_rect(fill="white"))
g.GFPBCsbyp

g.GFPRD <- ggplot(GFPRD, aes(RNAoDNA, fill=plasmid)) + geom_histogram(bins = 50) + scale_x_continuous(limits=c(0,350)) + scale_y_continuous(limits = c(0,15))+ labs(title="160901 R20 normalized RNA/DNA GFP barcode counts", x ="Barcode counts", y = "# of Counts") +
  facet_wrap(~condition, labeller = labeller(condition=induction_labels2), scales = "free") + theme(text = element_text(family="Calibri"), strip.background = element_rect(fill="white"))
g.GFPRD

#, labeller = labeller(condition=induction_labels2), geom_errorbar(gfpERR, colour = "black") +
gfpERR <- aes(ymax=mean_RNADNA + sd_RNADNA, ymin=mean_RNADNA - sd_RNADNA)
g.sumGFP <- ggplot(GFPRD, aes(clone, DNA,color=condition)) + geom_point()+ facet_wrap(~plasmid, scales = "free") + 
   labs(title="160904 R20 Average RNA/DNA counts by clone", x ="clone", y = "Mean RNA/DNA") + scale_y_continuous(limits=c(0,350)) + theme(text = element_text(family="Calibri"), strip.background = element_rect(fill="white"))
g.sumGFP

g.sumGFP2 <- ggplot(sumGFP, aes(plasmid, mean_RNADNA, condition)) + geom_bar(stat="identity")#+ facet_wrap(~condition, labeller = labeller(condition=induction_labels2), scales = "free")
g.sumGFP2

#####################################################################################
###########################Below here be dragons#####################################
#####################################################################################

#currently a slow motherfucker
for (n in 1:nrow(mapBCs)) {
  while (mapBCs[n,4] > 1) {
    bc <-unlist(strsplit(mapBCs[n,5],","))
    mapBCs[n,4] <- mapBCs[n,4] - 1
    mapBCs[n,5] <- paste(bc[2:length(bc)], collapse = ",")
    mBC <- data.frame(mapBCs[n,1:3], "num_bcs" = 1, "RevCompBC" = bc[1])
    mapBCs <- bind_rows(mapBCs, mBC)
  }
}


R20_DNA_b <- data.frame(R20_DNA_b, "RevCompBC" = revComp(as.character(R20_DNA_b$barcode)))

R20_RNA_0 <- data.frame(R20_RNA_0, "RevCompBC" = revComp(as.character(R20_RNA_0$barcode)))

mapGFP <- filter(bigDF, bigDF$RevCompBC.1 %in% GFPBCs$RevCompBC)

mapBCs <- read.csv("roman_counts.csv", skip = 1)
colnames(mapBCs) <- c("X peptide","Y peptide","reads","# of bcs", "RevCompBC")
mapBCs$RevCompBC <- as.character(mapBCs$RevCompBC)



for (n in 1:nrow(GFPBCs)){
  mapDNA_b <- bind_rows(mapDNA_b, filter(R20_DNA_b, grepl(GFPBCs$RevCompBC[n], R20_DNA_b$RevCompBC)))
}

mapDNA_b <- filter(R20_RNA_0, R20_RNA_0$RevCompBC %in% GFPBCs$RevCompBC)
mapDNA_b <- mutate(mapDNA_b, "percReads" = mapDNA_b$count/7824366*100)

mapRNA_0 <- inner_join(GFPBCs, R20_RNA_0, by = "RevCompBC")

mapBCsRNA0 <- right_join(R20_RNA_0, GFPBCs, by = "RevCompBC")
fullDNAb <- right_join(R20_DNA_b, mapBCs, by = "RevCompBC")
write.table(sumdata, "/Users/Cliff/Documents/Kosuri Lab/NGS files/160831 R20 Synzip run 1/summary stats.txt")

#write.table(bigDF, "/Users/Cliff/Documents/Kosuri Lab/NGS files/160831 R20 Synzip run 1/160901 combined counts.txt")

#write.table(sumdata, "/Users/Cliff/Documents/Kosuri Lab/NGS files/160831 R20 Synzip run 1/160901 line and total counts.txt")

#g.R20R0dist <- ggplot(big_R20_R_0, aes(count)) + geom_histogram(bins = 250) + scale_x_log10() + labs(title="160831 R20 RNA uninduced", x ="Barcode counts", y = "# of Counts")
#g.R20R0dist

g.R20D0dist <- ggplot(big_R20_D_0, aes(count)) + geom_histogram(bins = 250) + scale_x_log10() + labs(title="160831 R20 DNA uninduced", x ="Barcode counts", y = "# of Counts")
g.R20D0dist

