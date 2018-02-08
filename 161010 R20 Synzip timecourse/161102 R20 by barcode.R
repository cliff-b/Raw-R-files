library(dplyr)
library(tidyr)
library(extrafont)
library(viridis)
library(cowplot)

setwd("/Users/Cliff/Documents/Kosuri Lab/NGS files/161010 R20 Synzip timecourse/")

R20_RNA_0 <- read.table("counts_R20_RNA_0h.txt", header = FALSE)
R20_DNA_0 <- read.table("counts_R20_DNA_0h.txt", header = FALSE)
R20_RNA_a <- read.table("counts_R20_RNA_0.5h.txt", header = FALSE)
R20_DNA_a <- read.table("counts_R20_DNA_0.5h.txt", header = FALSE)
R20_RNA_b <- read.table("counts_R20_RNA_1h.txt", header = FALSE)
R20_DNA_b <- read.table("counts_R20_DNA_1h.txt", header = FALSE)
R20_RNA_c <- read.table("counts_R20_RNA_2h.txt", header = FALSE)
R20_DNA_c <-read.table("counts_R20_DNA_2h.txt", header = FALSE)
R20_RNA_d <- read.table("counts_R20_RNA_4h.txt", header = FALSE)
R20_DNA_d <-read.table("counts_R20_DNA_4h.txt", header = FALSE)


colnames(R20_RNA_0) <- c("RNA-0", "barcode")
colnames(R20_DNA_0) <- c("DNA-0", "barcode")
colnames(R20_RNA_a) <- c("RNA-A", "barcode")
colnames(R20_DNA_a) <- c("DNA-A", "barcode")
colnames(R20_RNA_b) <- c("RNA-B", "barcode")
colnames(R20_DNA_b) <- c("DNA-B", "barcode")
colnames(R20_RNA_c) <- c("RNA-C", "barcode")
colnames(R20_DNA_c) <- c("DNA-C", "barcode")
colnames(R20_RNA_d) <- c("RNA-D", "barcode")
colnames(R20_DNA_d) <- c("DNA-D", "barcode")


bigDF <- full_join(R20_RNA_0, R20_DNA_0, by = "barcode") %>%
  full_join(R20_RNA_a, by = "barcode") %>%
  full_join(R20_DNA_a, by = "barcode") %>%
  full_join(R20_RNA_b, by = "barcode") %>%
  full_join(R20_DNA_b, by = "barcode") %>%
  full_join(R20_RNA_c, by = "barcode") %>%
  full_join(R20_DNA_c, by = "barcode") %>%
  full_join(R20_RNA_d, by = "barcode") %>%
  full_join(R20_DNA_d, by = "barcode")

mapBCs <- read.csv("roman_counts.csv", skip = 1)
colnames(mapBCs) <- c("X_peptide","Y_peptide","reads","num_bcs", "RevCompBC")
mapBCs$RevCompBC <- as.character(mapBCs$RevCompBC)
mapBCs <- mutate(mapBCs, RevCompBC=strsplit(RevCompBC, ",")) %>%
  unnest(RevCompBC)
mapBCs <- mapBCs[!(duplicated(mapBCs$RevCompBC) | duplicated(mapBCs$RevCompBC, fromLast = TRUE)), ]


bigDF <- bigDF[,c(2,1,3,4,5,6,7,8,9,10,11)]
remove(R20_RNA_0, R20_RNA_a, R20_RNA_b, R20_RNA_c, R20_RNA_d, R20_DNA_0, R20_DNA_a, R20_DNA_b, R20_DNA_c, R20_DNA_d)
load(spreadR20Counts.Rda)

revComp <- function(x)
  chartr("ATGC","TACG", sapply(lapply(strsplit(x, NULL), rev), paste, collapse=""))

bigDF <- data.frame(bigDF, "RevCompBC" = revComp(as.character(bigDF$barcode)))
bigDF <- bigDF[,c(12,2,3,4,5,6,7,8,9,10,11)]
bigDF[is.na(bigDF)] <- 0

bigDF <- left_join(bigDF, mapBCs, by = "RevCompBC")

sumDF <- group_by(bigDF, X_peptide, Y_peptide) %>%
  summarise(sR0 = sum(RNA.0), sD0 = sum(DNA.0), sRa = sum(RNA.A), sDa = sum(DNA.A), sRb = sum(RNA.B), sDb = sum(DNA.B), sRc = sum(RNA.C), sDc = sum(DNA.C), sRd = sum(RNA.D), sDd = sum(DNA.D))

sumDF <- mutate(sumDF, RD0 = sR0/sD0, RDA = sRa/sRd, RDB = sRb/sDb, RDC = sRc/sDc, RDD = sRd/sDd) 
sumDF <- separate(sumDF, X_peptide, into = c("X_peptide", "rev", "X_codon"), sep = "_")
sumDF <- separate(sumDF, Y_peptide, into = c("Y_peptide", "Y_codon"), sep = "_")

g.R20RD <- ggplot(filter(sumDF, grepl(">P4", Y_peptide), !grepl("mS", Y_peptide)), aes(X_peptide, RDD)) + geom_boxplot(outlier.color = NA) + geom_jitter(alpha=0.3, aes(color=X_codon)) +
  theme(text=element_text(family="Calibri"), strip.text=element_text(size = 10), strip.background=element_rect(fill = "White"), axis.text.x = element_text(angle=90, hjust=1,size =6),
        axis.text.y = element_text(size = 6), legend.title=element_text(size = 8), legend.title.align=0.5, legend.text=element_text(size =8))
g.R20RD


g.sumR20RDmed <- ggplot(sumDF, aes(X_peptide, Y_peptide)) + geom_tile(aes(fill=RDD)) + #facet_wrap(~condition, labeller = labeller(condition=induction_labels2), scales = "free", ncol = 5) +
  theme(text=element_text(family="Calibri"), strip.text=element_text(size = 10), strip.background=element_rect(fill = "White"), axis.text.x = element_text(angle=90, hjust=1,size =6),
        axis.text.y = element_text(size = 6), legend.title=element_text(size = 8), legend.title.align=0.5, legend.text=element_text(size =8)) +
  scale_fill_viridis(option="viridis", name = "Median \nRNA to DNA ", na.value="white") + labs(title="161010 Roman's 20 RNA/DNA barcode ratio (median)", x="Protein X", y = "Protein Y")
g.sumR20RDmed

induction_labels2 <- c("0" = "0 hours induction RNA/DNA", "A" = "0.5 hours induction RNA/DNA", "B" = "1 hour induction RNA/DNA", "C" = "2 hour induction RNA/DNA", "D" = "4 hour induction RNA/DNA")


g.bccounts <- ggplot(bigDF, aes(RNA.C)) + geom_histogram(bins=300) + scale_x_continuous(limits = c(-1,302)) 
g.bccounts
