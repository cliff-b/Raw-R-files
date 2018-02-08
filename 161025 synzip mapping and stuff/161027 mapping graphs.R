library(dplyr)
library(tidyr)
library(cowplot)
library(viridis)
library(extrafont)

setwd("/Users/Cliff/Documents/Kosuri Lab/NGS files/161024 Synzip 3 codon mapping")

synbcs <- read.csv("synzip_counts.csv", skip = 2, header = FALSE)
colnames(synbcs) <- c("X_peptide", "Y_peptide","numcounts","numbcs","barcodes")
synbcs <- mutate(synbcs, barcodes=strsplit(as.character(barcodes),",")) %>%
  unnest(barcodes) %>%
  separate( X_peptide, into=c("X_peptide", "drop", "X_codon"), sep = "_") %>%
  separate(Y_peptide, into=c("Y_peptide", "Y_codon"), sep = "_") %>%
  select(-drop)

synbcs$X_peptide <- as.character(synbcs$X_peptide)
synbcs$Y_peptide <- as.character(synbcs$Y_peptide)

synbcs$X_peptide[synbcs$X_peptide == ">SYNZIP1[a]"] <- ">SYNZIP01[a]"
synbcs$X_peptide[synbcs$X_peptide == ">SYNZIP2[a]"] <- ">SYNZIP02[a]"
synbcs$X_peptide[synbcs$X_peptide == ">SYNZIP3[a]"] <- ">SYNZIP03[a]"
synbcs$X_peptide[synbcs$X_peptide == ">SYNZIP4[a]"] <- ">SYNZIP04[a]"
synbcs$X_peptide[synbcs$X_peptide == ">SYNZIP5[a]"] <- ">SYNZIP05[a]"
synbcs$X_peptide[synbcs$X_peptide == ">SYNZIP6[a]"] <- ">SYNZIP06[a]"
synbcs$X_peptide[synbcs$X_peptide == ">SYNZIP7[a]"] <- ">SYNZIP07[a]"
synbcs$X_peptide[synbcs$X_peptide == ">SYNZIP8[a]"] <- ">SYNZIP08[a]"
synbcs$X_peptide[synbcs$X_peptide == ">SYNZIP9[a]"] <- ">SYNZIP09[a]"

synbcs$Y_peptide[synbcs$Y_peptide == ">SYNZIP1[a]"] <- ">SYNZIP01[a]"
synbcs$Y_peptide[synbcs$Y_peptide == ">SYNZIP2[a]"] <- ">SYNZIP02[a]"
synbcs$Y_peptide[synbcs$Y_peptide == ">SYNZIP3[a]"] <- ">SYNZIP03[a]"
synbcs$Y_peptide[synbcs$Y_peptide == ">SYNZIP4[a]"] <- ">SYNZIP04[a]"
synbcs$Y_peptide[synbcs$Y_peptide == ">SYNZIP5[a]"] <- ">SYNZIP05[a]"
synbcs$Y_peptide[synbcs$Y_peptide == ">SYNZIP6[a]"] <- ">SYNZIP06[a]"
synbcs$Y_peptide[synbcs$Y_peptide == ">SYNZIP7[a]"] <- ">SYNZIP07[a]"
synbcs$Y_peptide[synbcs$Y_peptide == ">SYNZIP8[a]"] <- ">SYNZIP08[a]"
synbcs$Y_peptide[synbcs$Y_peptide == ">SYNZIP9[a]"] <- ">SYNZIP09[a]"



synbcs <- synbcs[!(duplicated(synbcs$barcodes) | duplicated(synbcs$barcodes, fromLast = TRUE)),]  

sumsynbcs <- group_by(synbcs, X_peptide, Y_peptide) %>%
  summarise(numbcs = n())

g.bccounts <- ggplot(sumsynbcs, aes(numbcs)) + geom_histogram(bins =100) + labs(title = "161024 number of barcodes per X-Y pair", x = "number of barcodes")
g.bccounts

g.tiles <- ggplot(sumsynbcs, aes(X_peptide, Y_peptide)) + geom_tile(aes(fill = numbcs)) + labs(title = "161024 Number of barcodes per X-Y pair", x = "X peptide", y = "Y peptide") +
  theme(text = element_text(family="Calibri"), axis.text.x = element_text(angle=90, size =6, hjust = .6), axis.text.y = element_text(size=6)) + scale_fill_viridis(name = "Number of \nbarcodes")
g.tiles

sumsynbcscodon <- group_by(synbcs, X_peptide, Y_peptide, X_codon, Y_codon) %>%
  summarise(numbcs = n())

g.tilesc1c2 <- ggplot(filter(sumsynbcscodon, X_codon == "c1" & Y_codon == "c2"), aes(X_peptide, Y_peptide)) + geom_tile(aes(fill = numbcs)) + labs(title = "161024 Number of barcodes per X c1-Y c2 pair", x = "X peptide", y = "Y peptide") +
  theme(text = element_text(family="Calibri"), axis.text.x = element_text(angle=90, size =6, hjust = .6), axis.text.y = element_text(size=6)) + scale_fill_viridis(name = "Number of \nbarcodes")
g.tilesc1c2

g.tilesc1c3 <- ggplot(filter(sumsynbcscodon, X_codon == "c1" & Y_codon == "c3"), aes(X_peptide, Y_peptide)) + geom_tile(aes(fill = numbcs)) + labs(title = "161024 Number of barcodes per X c1-Y c3 pair", x = "X peptide", y = "Y peptide") +
  theme(text = element_text(family="Calibri"), axis.text.x = element_text(angle=90, size =6, hjust = .6), axis.text.y = element_text(size=6)) + scale_fill_viridis(name = "Number of \nbarcodes")
g.tilesc1c3

g.tilesc1c3 <- ggplot(filter(sumsynbcscodon, X_codon == "c2" & Y_codon == "c1"), aes(X_peptide, Y_peptide)) + geom_tile(aes(fill = numbcs)) + labs(title = "161024 Number of barcodes per X c2-Y c1 pair", x = "X peptide", y = "Y peptide") +
  theme(text = element_text(family="Calibri"), axis.text.x = element_text(angle=90, size =6, hjust = .6), axis.text.y = element_text(size=6)) + scale_fill_viridis(name = "Number of \nbarcodes")
g.tilesc1c3

g.tilesc1c3 <- ggplot(filter(sumsynbcscodon, X_codon == "c2" & Y_codon == "c3"), aes(X_peptide, Y_peptide)) + geom_tile(aes(fill = numbcs)) + labs(title = "161024 Number of barcodes per X c2-Y c3 pair", x = "X peptide", y = "Y peptide") +
  theme(text = element_text(family="Calibri"), axis.text.x = element_text(angle=90, size =6, hjust = .6), axis.text.y = element_text(size=6)) + scale_fill_viridis(name = "Number of \nbarcodes")
g.tilesc1c3

g.tilesc3c1 <- ggplot(filter(sumsynbcscodon, X_codon == "c3" & Y_codon == "c1"), aes(X_peptide, Y_peptide)) + geom_tile(aes(fill = numbcs)) + labs(title = "161024 Number of barcodes per X c3-Y c1 pair", x = "X peptide", y = "Y peptide") +
  theme(text = element_text(family="Calibri"), axis.text.x = element_text(angle=90, size =6, hjust = .6), axis.text.y = element_text(size=6)) + scale_fill_viridis(name = "Number of \nbarcodes")
g.tilesc3c1

g.tilesc3c2 <- ggplot(filter(sumsynbcscodon, X_codon == "c2" & Y_codon == "c2"), aes(X_peptide, Y_peptide)) + geom_tile(aes(fill = numbcs)) + labs(title = "161024 Number of barcodes per X c2-Y c2 pair", x = "X peptide", y = "Y peptide") +
  theme(text = element_text(family="Calibri"), axis.text.x = element_text(angle=90, size =6, hjust = .6), axis.text.y = element_text(size=6)) + scale_fill_viridis(name = "Number of \nbarcodes")
g.tilesc3c2

g.box <- ggplot(sumsynbcscodon, aes(X_peptide, numbcs, color= X_codon))  + geom_boxplot(color="red", outlier.shape = NA) + geom_jitter(alpha = 0.2, size=1.5) + theme(text = element_text(family="Calibri"), axis.text.x = element_text(angle=90, size =6, hjust = .6)) +
  labs(title ="161024 Number of barcodes by X peptide", x = "X Peptide", y = "Number of BCs per Y Peptide") 
g.box

g.boxy <- ggplot(sumsynbcscodon, aes(Y_peptide, numbcs, color = Y_codon))  + geom_boxplot(color="red", outlier.shape = NA) + geom_jitter(alpha = 0.2, size=1.5) + theme(text = element_text(family="Calibri"), axis.text.x = element_text(angle=90, size =6, hjust = .6)) +
  labs(title ="161024 Number of barcodes by Y peptide", x = "Y Peptide", y = "Number of BCs per X Peptide") 
g.boxy
