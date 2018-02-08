library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(viridis)
library(extrafont)

#read in a horrific amount of data and label it with library and number of reads from that library
setwd("/Users/Cliff/Documents/Kosuri Lab/NGS files/160831 R20 Synzip run 1")
Syn_RNA_0 <- read.table("Syn_RNA_0.txt", header = FALSE)
Syn_RNA_0 <- data.frame(Syn_RNA_0, rep("RNA-0",nrow(Syn_RNA_0)),rep(sum(Syn_RNA_0$V1), nrow(Syn_RNA_0)))
colnames(Syn_RNA_0) <- c("count", "barcode", "condition","readsinlib")

Syn_DNA_0 <- read.table("Syn_DNA_0.txt", header = FALSE)
Syn_DNA_0 <- data.frame(Syn_DNA_0, rep("DNA-0",nrow(Syn_DNA_0)),rep(sum(Syn_DNA_0$V1), nrow(Syn_DNA_0)))
colnames(Syn_DNA_0) <- c("count", "barcode","condition","readsinlib")

Syn_RNA_a <- read.table("Syn_RNA_A.txt", header = FALSE)
Syn_RNA_a <- data.frame(Syn_RNA_a, rep("RNA-A",nrow(Syn_RNA_a)),rep(sum(Syn_RNA_a$V1), nrow(Syn_RNA_a)))
colnames(Syn_RNA_a) <- c("count", "barcode", "condition","readsinlib")

Syn_DNA_a <-read.table("Syn_DNA_A.txt", header = FALSE)
Syn_DNA_a <- data.frame(Syn_DNA_a, rep("DNA-A",nrow(Syn_DNA_a)),rep(sum(Syn_DNA_a$V1), nrow(Syn_DNA_a)))
colnames(Syn_DNA_a) <- c("count", "barcode","condition","readsinlib")

Syn_RNA_b <- read.table("Syn_RNA_B.txt", header = FALSE)
Syn_RNA_b <- data.frame(Syn_RNA_b, rep("RNA-B",nrow(Syn_RNA_b)),rep(sum(Syn_RNA_b$V1), nrow(Syn_RNA_b)))
colnames(Syn_RNA_b) <- c("count", "barcode", "condition","readsinlib")

Syn_DNA_b <- read.table("Syn_DNA_B.txt", header = FALSE)
Syn_DNA_b <- data.frame(Syn_DNA_b, rep("DNA-B",nrow(Syn_DNA_b)),rep(sum(Syn_DNA_b$V1), nrow(Syn_DNA_b)))
colnames(Syn_DNA_b) <- c("count", "barcode","condition","readsinlib")

Syn_RNA_c <- read.table("Syn_RNA_C.txt", header = FALSE)
Syn_RNA_c <- data.frame(Syn_RNA_c, rep("RNA-C",nrow(Syn_RNA_c)),rep(sum(Syn_RNA_c$V1), nrow(Syn_RNA_c)))
colnames(Syn_RNA_c) <- c("count", "barcode", "condition","readsinlib")

Syn_DNA_c <-read.table("Syn_DNA_A.txt", header = FALSE)
Syn_DNA_c <- data.frame(Syn_DNA_c, rep("DNA-C",nrow(Syn_DNA_c)),rep(sum(Syn_DNA_c$V1), nrow(Syn_DNA_c)))
colnames(Syn_DNA_c) <- c("count", "barcode","condition","readsinlib")

GFPBCs <-read.csv("160719 constitutive GFP Barcodes.csv", header = FALSE)
GFPBCs <- data.frame(GFPBCs, rep("GAATTCTGGGACC", nrow(GFPBCs)))
colnames(GFPBCs) <- c("plasmid", "barcode","filler")
GFPBCs$RevCompBC <- paste(GFPBCs$barcode, GFPBCs$filler, sep = "")

mapBCs <- read.csv("synzip_counts.csv", skip = 1)
colnames(mapBCs) <- c("X_peptide","Y_peptide","reads","num_bcs", "RevCompBC")
mapBCs$RevCompBC <- as.character(mapBCs$RevCompBC)
mapBCs <- mutate(mapBCs, RevCompBC=strsplit(RevCompBC, ",")) %>%
  unnest(RevCompBC)
mapBCs <- mapBCs[!(duplicated(mapBCs$RevCompBC) | duplicated(mapBCs$RevCompBC, fromLast = TRUE)), ]

###################Manipulate the data in to something useful############################

revComp <- function(x)
  chartr("ATGC","TACG", sapply(lapply(strsplit(x, NULL), rev), paste, collapse=""))

#make the big one
bigDF <- bind_rows(Syn_RNA_0, Syn_RNA_a, Syn_RNA_b, Syn_RNA_c, Syn_DNA_0, Syn_DNA_a, Syn_DNA_b, Syn_DNA_c)
remove(Syn_RNA_0, Syn_RNA_a, Syn_RNA_b, Syn_RNA_c, Syn_DNA_0, Syn_DNA_a, Syn_DNA_b, Syn_DNA_c)
bigDF <- data.frame(bigDF, "RevCompBC" = revComp(as.character(bigDF$barcode)))
bigDF <- mutate(bigDF, "nCount" = count/readsinlib*1000000)
bigDF <- filter(bigDF, count > 5 | grepl("RNA", bigDF$condition))




fullMap <- filter(left_join(bigDF, synBC, by = "RevCompBC"), !is.na(reads))
SynRNA <- filter(fullMap, grepl("RNA", fullMap$condition)) %>%
  select(X_peptide, Y_peptide, condition, RevCompBC, "rcount" = nCount ) %>%
  separate(condition, into = c("template", "condition"))
SynDNA <- filter(fullMap, grepl("DNA", fullMap$condition)) %>%
  select(X_peptide, Y_peptide, condition, RevCompBC, "dcount" = nCount ) %>%
  separate(condition, into = c("template", "condition"))
SynRD <- inner_join(SynRNA,SynDNA, by = c("RevCompBC","X_peptide","Y_peptide","condition")) %>%
  mutate("RNADNA" = rcount/dcount*100) %>%
  separate("Y_peptide", into = c("Y_peptide", "Y_codon"), sep = "_") %>%
  separate ("X_peptide", into = c("X_peptide", "X_codon"), sep = "_")


sumSyn <- group_by(SynRD, X_peptide, Y_peptide, condition) %>%
  summarise(mean_RNA = mean(rcount), mean_DNA = mean(dcount), mean_RD = mean(RNADNA), sd_RD = sd(RNADNA), bcnum= n_distinct(RevCompBC))


induction_labels2 <- c("0" = "Uninduced RNA/DNA", "A" = "5uM DAPG 5ng/mL ATc RNA/DNA", "B" = "7.5uM DAPG 10ng/mL ATc RNA/DNA", "C" = "15uM DAPG 40ng/mL ATc RNA/DNA")

g.sumSynRD <- ggplot(filter(sumSyn, bcnum>4), aes(X_peptide, Y_peptide)) + geom_tile(aes(fill=mean_RD)) + facet_wrap(~condition, labeller = labeller(condition=induction_labels2), scales = "free") +
  theme(text=element_text(family="Calibri"), strip.text=element_text(size = 10), strip.background=element_rect(fill = "White"), axis.text.x = element_text(angle=90, hjust=1,size =6),
        axis.text.y = element_text(size = 6), legend.title=element_text(size = 8), legend.title.align=0.5, legend.text=element_text(size =8)) +
  scale_fill_viridis(option="inferno", name = "Normalized \nRNA to DNA \n", na.value="white") + labs(title="160914 Human bzips in Synzips RNA/DNA barcode ratio", x="X peptide", y = "Y peptide")
g.sumSynRD

g.sumSynNum <- ggplot(sumSyn, aes(X_peptide, Y_peptide)) + geom_tile(aes(fill=bcnum)) + facet_wrap(~condition, labeller = labeller(condition=induction_labels2), scales = "free") +
  theme(text=element_text(family="Comic Sans MS"), strip.text=element_text(size = 10), strip.background=element_rect(fill = "White"), axis.text.x = element_text(angle=90, hjust=1,size =6),
        axis.text.y = element_text(size = 6), legend.title=element_text(size = 8), legend.title.align=0.5, legend.text=element_text(size =8), plot.title=element_text(family="Jokerman")) +
  scale_fill_viridis(option="inferno", name = "Normalized \nRNA to DNA \n(Log10)", na.value="white") + labs(title="160906 Synzip RNA/DNA barcode counts", x="X peptide", y = "Y peptide")
g.sumSynNum

g.sumSynD <- ggplot(sumSyn, aes(X_peptide, Y_peptide)) + geom_tile(aes(fill=mean_DNA)) + facet_wrap(~condition, labeller = labeller(condition=induction_labels2), scales = "free") +
  theme(text=element_text(family="Comic Sans MS"), strip.text=element_text(size = 10), strip.background=element_rect(fill = "White"), axis.text.x = element_text(angle=90, hjust=1,size =6),
        axis.text.y = element_text(size = 6), legend.title=element_text(size = 8), legend.title.align=0.5, legend.text=element_text(size =8), plot.title=element_text(family="Jokerman")) +
  scale_fill_viridis(option="inferno", name = "Normalized \nDNA \n(Log10)", na.value="white") + labs(title="160906 Synzip DNA barcode ratio", x="X peptide", y = "Y peptide")
g.sumSynD




#############################GFP manipulations####################################################
mapGFP <- filter(left_join(bigDF,GFPBCs, by = "RevCompBC"),!is.na(plasmid))
mapGFP <- separate(mapGFP,plasmid, into = c("plasmid", "clone"),sep = "-")

GFPRNA <- filter(mapGFP, grepl("RNA", mapGFP$condition)) %>%
  select(plasmid, clone, condition, RevCompBC, "rcount" = nCount ) %>%
  separate(condition, into = c("template", "condition"))
GFPDNA <- filter(mapGFP, grepl("DNA", mapGFP$condition)) %>%
  select(plasmid, clone, condition, RevCompBC, "dcount" = nCount ) %>%
  separate(condition, into = c("template", "condition"))
GFPRD <- inner_join(GFPRNA,GFPDNA, by = c("RevCompBC","condition","plasmid","clone")) %>%
  mutate("RNAoDNA" = rcount/dcount*100)%>%
  # get rid of p63-4 which seems to be mismatched
  filter(RevCompBC != "GTTACCAGAATTCTGGGACC")
sumGFP <- group_by(GFPRD, plasmid, condition) %>%
  summarise(mean_rcount =  mean(rcount), mean_dcount = mean(dcount), mean_RNADNA = mean(RNAoDNA), sd_RNADNA=sd(RNAoDNA))

induction_labels <- c("DNA-0" = "Uninduced DNA", "DNA-A" = "5uM DAPG 5ng/mL ATc DNA", "DNA-B" = "7.5uM DAPG 10ng/mL ATc DNA", "DNA-C" = "15uM DAPG 40ng/mL ATc DNA",
                      "RNA-0" = "Uninduced RNA", "RNA-A" = "5uM DAPG 5ng/mL ATc RNA", "RNA-B" = "7.5uM DAPG 10ng/mL ATc RNA", "RNA-C" = "15uM DAPG 40ng/mL ATc RNA")


g.GFPBCs <- ggplot(mapGFP, aes(count)) + geom_histogram(bins = 50, fill= "green3", color="black") + scale_x_log10() + labs(title="160913 Synzip number of counts per GFP barcode", x ="Barcode counts", y = "# of Counts") +
  facet_wrap(~condition, labeller = labeller(condition=induction_labels)) + theme(text = element_text(family="Calibri"), strip.background = element_rect(fill="white"))
g.GFPBCs

g.GFPBCsbyp <- ggplot(mapGFP, aes(nCount, fill=plasmid)) + geom_histogram(bins = 50) + scale_x_continuous(limits=c(0,370)) + scale_y_continuous(limits = c(0,15))+ labs(title="160913 synzip number of normalized counts per GFP barcode", x ="Barcode counts", y = "# of Counts") +
  facet_wrap(~condition, labeller = labeller(condition=induction_labels), scales = "free") + theme(text = element_text(family="Calibri"), strip.background = element_rect(fill="white"))
g.GFPBCsbyp

g.GFPRD <- ggplot(GFPRD, aes(RNAoDNA, fill=plasmid)) + geom_histogram(bins = 50) + scale_x_continuous(limits=c(0,350)) + scale_y_continuous(limits = c(0,15))+ labs(title="160913 Synzip normalized RNA/DNA GFP barcode counts", x ="Barcode counts", y = "# of Counts") +
  facet_wrap(~condition, labeller = labeller(condition=induction_labels2), scales = "free") + theme(text = element_text(family="Calibri"), strip.background = element_rect(fill="white"))
g.GFPRD

#, labeller = labeller(condition=induction_labels2), geom_errorbar(gfpERR, colour = "black") +
gfpERR <- aes(ymax=mean_RNADNA + sd_RNADNA, ymin=mean_RNADNA - sd_RNADNA)
g.sumGFP <- ggplot(sumGFP, aes(plasmid, mean_RNADNA,fill=plasmid)) + geom_bar(stat="identity")+ facet_wrap(~condition, scales = "free", labeller = labeller(condition=induction_labels2)) + geom_errorbar(gfpERR, colour = "black") +
  labs(title="160913 Synzip Average RNA/DNA counts by plasmid", x ="plasmid", y = "Mean RNA/DNA") + scale_y_continuous(limits=c(0,400)) + theme(text = element_text(family="Calibri"), strip.background = element_rect(fill="white"))
g.sumGFP

g.sumGFP2 <- ggplot(sumGFP, aes(plasmid, mean_RNADNA, condition)) + geom_bar(stat="identity")#+ facet_wrap(~condition, labeller = labeller(condition=induction_labels2), scales = "free")
g.sumGFP2








mapGFP <- filter(left_join(bigDF,GFPBCs, by = "RevCompBC"),!is.na(plasmid))
mapGFP <- separate(mapGFP,plasmid, into = c("plasmid", "clone"),sep = "-")