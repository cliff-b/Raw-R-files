library(dplyr)
library(tidyr)
library(viridis)
library(cowplot)

setwd("/Users/Cliff/Documents/Kosuri Lab/NGS files/161010 R20 Synzip timecourse/")
load("SynzipDF.Rda")
load("Synzips_map.csv")
load("synzip_counts.csv")
Octsyns <- read.csv("OCTsynzip_counts.csv", skip = 2, header = FALSE)
colnames(Octsyns) <- c("X_peptide","Y_peptide","reads","num_bcs", "RevCompBC")
Octsyns$RevCompBC <- as.character(Octsyns$RevCompBC)
Octsyns <- mutate(Octsyns, RevCompBC=strsplit(RevCompBC, ",")) %>%
  unnest(RevCompBC)
Octsyns <- Octsyns[!(duplicated(Octsyns$RevCompBC) | duplicated(Octsyns$RevCompBC, fromLast = TRUE)), ]


maps <- read.csv("synzip_counts.csv", header = TRUE)
colnames(maps) <- c("X_peptide","Y_peptide","reads","num_bcs", "RevCompBC")
maps$RevCompBC <- as.character(maps$RevCompBC)
maps <- mutate(maps, RevCompBC=strsplit(RevCompBC, ",")) %>%
  unnest(RevCompBC)
maps <- maps[!(duplicated(maps$RevCompBC) | duplicated(maps$RevCompBC, fromLast = TRUE)), ]

revComp <- function(x)
  chartr("ATGC","TACG", sapply(lapply(strsplit(x, NULL), rev), paste, collapse=""))
maps <- data.frame(maps, "BC" = revComp(as.character(maps$RevCompBC)))

allmap <- full_join(maps, synbcs, by = c("BC"))
synbcs <- data.frame(synbcs, "RevCompBC" = revComp(as.character(synbcs$BC)))


bigDF <- mutate(bigDF, "nCount" = count/readsinlib*1000000)
bigDF <- left_join(bigDF, Octsyns, by = c("RevCompBC"))
bigDF <- select(bigDF, -X_peptide, -Y_peptide, -Reads, -BC)

setwd("/Users/Cliff/Documents/Kosuri Lab/NGS files/161025 synzip mapping and stuff")
syn_DNA_0 <- read.table("counts_Synzip_DNA_0h.txt", header = FALSE)
syn_DNA_0 <- data.frame(syn_DNA_0, rep("DNA-0",nrow(syn_DNA_0)),rep(sum(syn_DNA_0$V1), nrow(syn_DNA_0)))
colnames(syn_DNA_0) <- c("count", "barcode","condition","readsinlib")

syn_DNA_05 <- read.table("counts_Synzip_DNA_0.5h.txt", header = FALSE)
syn_DNA_05 <- data.frame(syn_DNA_05, rep("DNA-A",nrow(syn_DNA_05)),rep(sum(syn_DNA_05$V1), nrow(syn_DNA_05)))
colnames(syn_DNA_05) <- c("count", "barcode","condition","readsinlib")

syn_DNA_1 <- read.table("counts_Synzip_DNA_1h.txt", header = FALSE)
syn_DNA_1 <- data.frame(syn_DNA_1, rep("DNA-B",nrow(syn_DNA_1)),rep(sum(syn_DNA_1$V1), nrow(syn_DNA_1)))
colnames(syn_DNA_1) <- c("count", "barcode","condition","readsinlib")

syn_DNA_2 <- read.table("counts_Synzip_DNA_2h.txt", header = FALSE)
syn_DNA_2 <- data.frame(syn_DNA_2, rep("DNA-C",nrow(syn_DNA_2)),rep(sum(syn_DNA_2$V1), nrow(syn_DNA_2)))
colnames(syn_DNA_2) <- c("count", "barcode","condition","readsinlib")

syn_DNA_4 <- read.table("counts_Synzip_DNA_4h.txt", header = FALSE)
syn_DNA_4 <- data.frame(syn_DNA_4, rep("DNA-D",nrow(syn_DNA_4)),rep(sum(syn_DNA_4$V1), nrow(syn_DNA_4)))
colnames(syn_DNA_4) <- c("count", "barcode","condition","readsinlib")

synDNA <- bind_rows(syn_DNA_0,syn_DNA_05,syn_DNA_1,syn_DNA_2,syn_DNA_4)

R20_DNA_0 <- left_join(R20_DNA_0, synbcs, by = c("barcode" = "BC"))
R20bcs <- load("R20bigDF1.Rda")
