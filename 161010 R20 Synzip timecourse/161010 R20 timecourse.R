library(dplyr)
library(tidyr)
library(extrafont)
library(viridis)
library(cowplot)

setwd("/Users/Cliff/Documents/Kosuri Lab/NGS files/161010 R20 Synzip timecourse/")

R20_RNA_0 <- read.table("counts_R20_RNA_0h.txt", header = FALSE)
R20_RNA_0 <- data.frame(R20_RNA_0, rep("RNA-0",nrow(R20_RNA_0)),rep(sum(R20_RNA_0$V1), nrow(R20_RNA_0)))
colnames(R20_RNA_0) <- c("count", "barcode", "condition","readsinlib")

R20_DNA_0 <- read.table("counts_R20_DNA_0h.txt", header = FALSE)
R20_DNA_0 <- data.frame(R20_DNA_0, rep("DNA-0",nrow(R20_DNA_0)),rep(sum(R20_DNA_0$V1), nrow(R20_DNA_0)))
colnames(R20_DNA_0) <- c("count", "barcode","condition","readsinlib")

R20_RNA_a <- read.table("counts_R20_RNA_0.5h.txt", header = FALSE)
R20_RNA_a <- data.frame(R20_RNA_a, rep("RNA-A",nrow(R20_RNA_a)),rep(sum(R20_RNA_a$V1), nrow(R20_RNA_a)))
colnames(R20_RNA_a) <- c("count", "barcode", "condition","readsinlib")

R20_DNA_a <-read.table("counts_R20_DNA_0.5h.txt", header = FALSE)
R20_DNA_a <- data.frame(R20_DNA_a, rep("DNA-A",nrow(R20_DNA_a)),rep(sum(R20_DNA_a$V1), nrow(R20_DNA_a)))
colnames(R20_DNA_a) <- c("count", "barcode","condition","readsinlib")

R20_RNA_b <- read.table("counts_R20_RNA_1h.txt", header = FALSE)
R20_RNA_b <- data.frame(R20_RNA_b, rep("RNA-B",nrow(R20_RNA_b)),rep(sum(R20_RNA_b$V1), nrow(R20_RNA_b)))
colnames(R20_RNA_b) <- c("count", "barcode", "condition","readsinlib")

R20_DNA_b <- read.table("counts_R20_DNA_1h.txt", header = FALSE)
R20_DNA_b <- data.frame(R20_DNA_b, rep("DNA-B",nrow(R20_DNA_b)),rep(sum(R20_DNA_b$V1), nrow(R20_DNA_b)))
colnames(R20_DNA_b) <- c("count", "barcode","condition","readsinlib")

R20_RNA_c <- read.table("counts_R20_RNA_2h.txt", header = FALSE)
R20_RNA_c <- data.frame(R20_RNA_c, rep("RNA-C",nrow(R20_RNA_c)),rep(sum(R20_RNA_c$V1), nrow(R20_RNA_c)))
colnames(R20_RNA_c) <- c("count", "barcode", "condition","readsinlib")

R20_DNA_c <-read.table("counts_R20_DNA_2h.txt", header = FALSE)
R20_DNA_c <- data.frame(R20_DNA_c, rep("DNA-C",nrow(R20_DNA_c)),rep(sum(R20_DNA_c$V1), nrow(R20_DNA_c)))
colnames(R20_DNA_c) <- c("count", "barcode","condition","readsinlib")

R20_RNA_d <- read.table("counts_R20_RNA_4h.txt", header = FALSE)
R20_RNA_d <- data.frame(R20_RNA_d, rep("RNA-D",nrow(R20_RNA_d)),rep(sum(R20_RNA_d$V1), nrow(R20_RNA_d)))
colnames(R20_RNA_d) <- c("count", "barcode", "condition","readsinlib")

R20_DNA_d <-read.table("counts_R20_DNA_4h.txt", header = FALSE)
R20_DNA_d <- data.frame(R20_DNA_d, rep("DNA-D",nrow(R20_DNA_d)),rep(sum(R20_DNA_d$V1), nrow(R20_DNA_d)))
colnames(R20_DNA_d) <- c("count", "barcode","condition","readsinlib")

mapBCs <- read.csv("roman_counts.csv", skip = 1)
colnames(mapBCs) <- c("X_peptide","Y_peptide","reads","num_bcs", "RevCompBC")
mapBCs$RevCompBC <- as.character(mapBCs$RevCompBC)
mapBCs <- mutate(mapBCs, RevCompBC=strsplit(RevCompBC, ",")) %>%
  unnest(RevCompBC)
mapBCs <- mapBCs[!(duplicated(mapBCs$RevCompBC) | duplicated(mapBCs$RevCompBC, fromLast = TRUE)), ]


bigDF <- bind_rows(R20_RNA_0, R20_RNA_a, R20_RNA_b, R20_RNA_c, R20_RNA_d, R20_DNA_0, R20_DNA_a, R20_DNA_b, R20_DNA_c, R20_DNA_d)
remove(R20_RNA_0, R20_RNA_a, R20_RNA_b, R20_RNA_c, R20_RNA_d, R20_DNA_0, R20_DNA_a, R20_DNA_b, R20_DNA_c, R20_DNA_d)

load("bigDF.Rda")
bigDF <- spread(bigDF, condition, c(count, barcode))

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
  summarise(mean_RNA = mean(rcount), sd_R = sd(rcount), mean_DNA = mean(dcount), sd_D = sd(dcount), mean_RD = mean(RNADNA), sd_RD = sd(RNADNA), med_RD = median(RNADNA), bcnum = n_distinct(RevCompBC))

induction_labels2 <- c("0" = "0 hours induction RNA/DNA", "A" = "0.5 hours induction RNA/DNA", "B" = "1 hour induction RNA/DNA", "C" = "2 hour induction RNA/DNA", "D" = "4 hour induction RNA/DNA")

g.sumR20RD <- ggplot(filter(sumR20D5, !is.na(X_peptide)), aes(X_peptide, Y_peptide)) + geom_tile(aes(fill=mean_RD)) + facet_wrap(~condition, labeller = labeller(condition=induction_labels2), scales = "free") +
  theme(text=element_text(family="Calibri"), strip.text=element_text(size = 10), strip.background=element_rect(fill = "White"), axis.text.x = element_text(angle=90, hjust=1,size =6),
        axis.text.y = element_text(size = 6), legend.title=element_text(size = 8), legend.title.align=0.5, legend.text=element_text(size =8)) +
  scale_fill_viridis(option="viridis", name = "Normalized \nRNA to DNA ", na.value="white") + labs(title="161010 Roman's 20 RNA/DNA barcode ratio", x="X peptide", y = "Y peptide")
g.sumR20RD

g.sumR20RDmed <- ggplot(filter(sumR20D5, !is.na(X_peptide)), aes(X_peptide, Y_peptide)) + geom_tile(aes(fill=med_RD)) + facet_wrap(~condition, labeller = labeller(condition=induction_labels2), scales = "free") +
  theme(text=element_text(family="Calibri"), strip.text=element_text(size = 10), strip.background=element_rect(fill = "White"), axis.text.x = element_text(angle=90, hjust=1,size =6),
        axis.text.y = element_text(size = 6), legend.title=element_text(size = 8), legend.title.align=0.5, legend.text=element_text(size =8)) +
  scale_fill_viridis(option="viridis", name = "Median \nRNA to DNA ", na.value="white") + labs(title="161010 Roman's 20 RNA/DNA barcode ratio (median)", x="X peptide", y = "Y peptide")
g.sumR20RDmed

g.sumR20RDsd <- ggplot(filter(sumR20D5, !is.na(X_peptide)), aes(X_peptide, Y_peptide)) + geom_tile(aes(fill=sd_RD)) + facet_wrap(~condition, labeller = labeller(condition=induction_labels2), scales = "free") +
  theme(text=element_text(family="Calibri"), strip.text=element_text(size = 10), strip.background=element_rect(fill = "White"), axis.text.x = element_text(angle=90, hjust=1,size =6),
        axis.text.y = element_text(size = 6), legend.title=element_text(size = 8), legend.title.align=0.5, legend.text=element_text(size =8)) +
  scale_fill_viridis(option="viridis", name = "Standard deviation \n of RNA to DNA ", na.value="white") + labs(title="161010 Roman's 20 RNA/DNA standard deviation of barcode ratio", x="X peptide", y = "Y peptide")
g.sumR20RDsd


g.sumR20RDbc <- ggplot(filter(sumR20D5, !is.na(X_peptide)), aes(X_peptide, Y_peptide)) + geom_tile(aes(fill=bcnum)) + facet_wrap(~condition, labeller = labeller(condition=induction_labels2), scales = "free") +
  theme(text=element_text(family="Calibri"), strip.text=element_text(size = 10), strip.background=element_rect(fill = "White"), axis.text.x = element_text(angle=90, hjust=1,size =6),
        axis.text.y = element_text(size = 6), legend.title=element_text(size = 8), legend.title.align=0.5, legend.text=element_text(size =8)) +
  scale_fill_viridis(option="viridis", name = "Number of\n barcodes ", na.value="white") + labs(title="161010 Roman's 20 number of barcodes per construct", x="X peptide", y = "Y peptide")
g.sumR20RDbc

g.sdvmean <- ggplot(filter(sumR20D5, !is.na(X_peptide)), aes(mean_RD, sd_RD, color=condition)) +geom_point(alpha=0.4) + scale_y_continuous(limits = c(0,400))+ scale_x_continuous(limits=c(0,480))+ #, labels = c("Uninduced RNA/DNA","5uM DAPG 5ng/mL ATc RNA/DNA","7.5uM DAPG 10ng/mL ATc RNA/DNA","15uM DAPG 40ng/mL ATc RNA/DNA")) + 
  labs(title="161010 Roman's 20 Standard Deviation vs mean of RNA/DNA ratio", y="Standard Deviation", x = "Mean RNA/DNA ratio")
g.sdvmean
