library(dplyr)
library(tidyr)
library(cowplot)
library(viridis)

setwd("/Users/Cliff/Documents/Kosuri Lab/NGS files/161025 synzip mapping and stuff")

synbcs1 <- read.table("SXYb-1.x-y.txt", header = TRUE)
synbcs2 <- read.table("SXYb-2.x-y.txt", header = TRUE)
synbcs3 <- read.table("SXYb-3.x-y.txt", header = TRUE)

synbcs <- bind_rows(synbcs1, synbcs2)
synbcs <- bind_rows(synbcs, synbcs3)

synbcs <- filter(synbcs, X != "<NA>", Y != "<NA>")
synbcs <- separate(synbcs, X, into = c("X_peptide", "X_codon"), sep = "_rev_") %>%
  separate(Y, into = c("Y_peptide", "Y_codon"), sep = -3)

synbcs$X_peptide[synbcs$X_peptide == "SYNZIP1[A]"] <- "SYNZIP01[A]"
synbcs$X_peptide[synbcs$X_peptide == "SYNZIP2[A]"] <- "SYNZIP02[A]"
synbcs$X_peptide[synbcs$X_peptide == "SYNZIP3[A]"] <- "SYNZIP03[A]"
synbcs$X_peptide[synbcs$X_peptide == "SYNZIP4[A]"] <- "SYNZIP04[A]"
synbcs$X_peptide[synbcs$X_peptide == "SYNZIP5[A]"] <- "SYNZIP05[A]"
synbcs$X_peptide[synbcs$X_peptide == "SYNZIP6[A]"] <- "SYNZIP06[A]"
synbcs$X_peptide[synbcs$X_peptide == "SYNZIP7[A]"] <- "SYNZIP07[A]"
synbcs$X_peptide[synbcs$X_peptide == "SYNZIP8[A]"] <- "SYNZIP08[A]"
synbcs$X_peptide[synbcs$X_peptide == "SYNZIP9[A]"] <- "SYNZIP09[A]"
synbcs$X_peptide[synbcs$X_peptide == "GCN4_1"] <- "GCN4"

synbcs$Y_peptide[synbcs$Y_peptide == "SYNZIP1[A]"] <- "SYNZIP01[A]"
synbcs$Y_peptide[synbcs$Y_peptide == "SYNZIP2[A]"] <- "SYNZIP02[A]"
synbcs$Y_peptide[synbcs$Y_peptide == "SYNZIP3[A]"] <- "SYNZIP03[A]"
synbcs$Y_peptide[synbcs$Y_peptide == "SYNZIP4[A]"] <- "SYNZIP04[A]"
synbcs$Y_peptide[synbcs$Y_peptide == "SYNZIP5[A]"] <- "SYNZIP05[A]"
synbcs$Y_peptide[synbcs$Y_peptide == "SYNZIP6[A]"] <- "SYNZIP06[A]"
synbcs$Y_peptide[synbcs$Y_peptide == "SYNZIP7[A]"] <- "SYNZIP07[A]"
synbcs$Y_peptide[synbcs$Y_peptide == "SYNZIP8[A]"] <- "SYNZIP08[A]"
synbcs$Y_peptide[synbcs$Y_peptide == "SYNZIP9[A]"] <- "SYNZIP09[A]"

synbcs <- separate(synbcs, X_peptide, into = c("X_peptide", "crap"), sep = "\\[")
synbcs <- separate(synbcs, Y_peptide, into = c("Y_peptide", "crap2"), sep = "\\[")
synbcs <- select(synbcs, -crap, -crap2)

synbcs <- distinct(synbcs, BC, X_peptide, Y_peptide, Reads)
otherbcs <- synbcs[(duplicated(synbcs$BC) | duplicated(synbcs$BC, fromLast = TRUE)),]
synbcs <- synbcs[!(duplicated(synbcs$BC) | duplicated(synbcs$BC, fromLast = TRUE)),]

g.bccounts <- ggplot(synbcs, aes(Reads)) + geom_histogram(bins = 51) + labs(title = "170522 number of reads per barcodes ", x = "number of reads per barcodes") + scale_x_continuous(limits = c(0,50))
g.bccounts

sumsynbcs <- group_by(synbcs, X_peptide, Y_peptide) %>%
  summarise(numbcs = n(), sumreads = sum(Reads))

g.tiles <- ggplot(sumsynbcs, aes(X_peptide, Y_peptide)) + geom_tile(aes(fill = numbcs)) + labs(title = "170519 Number of barcodes per X-Y pair", x = "X peptide", y = "Y peptide") +
  theme(text = element_text(family="Calibri"), axis.text.x = element_text(angle=90, size =6), axis.text.y = element_text(size=6)) + scale_fill_viridis(name = "Number of \nbarcodes")
g.tiles
