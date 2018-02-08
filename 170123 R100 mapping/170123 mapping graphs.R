library(dplyr)
library(tidyr)
library(cowplot)
library(viridis)
library(extrafont)

setwd("/Users/Cliff/Documents/Kosuri Lab/NGS files/170123 R100 mapping")

R100bcs <- read.table("Codon-1.x-y.txt", sep = "", header = TRUE)
R100bcs <- bind_rows(R100bcs, read.table("Codon-2.x-y.txt", sep = "", header = TRUE))
R100bcs <- bind_rows(R100bcs, read.table("Codon-3.x-y.txt", sep = "", header = TRUE))

colnames(R100bcs) <- c("barcodes",  "X_peptide", "Y_peptide")
R100bcs <-  separate(R100bcs, X_peptide, into=c("X_peptide", "X_codon"), sep = "_rev_") %>%
  separate(Y_peptide, into=c("Y_peptide", "Y_codon"), sep = -3) 
 

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



R100bcsd <- R100bcs[!(duplicated(R100bcs$barcodes) | duplicated(R100bcs$barcodes, fromLast = TRUE)),]  
R100bcsmap <- filter(R100bcsd, X_codon != "<NA>", Y_peptide != "<NA>")
g.bcreads <- ggplot(R100bcsmap, aes(numcounts))  + labs(title = "170213 number of reads per BC", x = "# of reads per bc")
g.bcreads
+ geom_histogram(bins=59) + geom_point(data=psbcs, aes(nums,poisson), color="red", alpha = 0.6)

g.junk <-ggplot()

sumR100bcs <- group_by(R100bcsmap, X_peptide, Y_peptide) %>%
  summarise(numbcs = n())
sumR100bcs <- filter(sumR100bcs, X_peptide != Y_peptide)

g.bccounts <- ggplot(sumR100bcs, aes(numbcs)) + geom_histogram(bins =59) + labs(title = "170214 number of barcodes per X-Y pair", x = "number of barcodes") + scale_x_continuous(limits = c(0,300))
g.bccounts

g.tiles <- ggplot(sumR100bcs, aes(X_peptide, Y_peptide)) + geom_tile(aes(fill = numbcs)) + labs(title = "170123 Number of barcodes per X-Y pair", x = "X peptide", y = "Y peptide") +
  theme(text = element_text(family="Calibri"), axis.text.x = element_text(angle=90, size =6), axis.text.y = element_text(size=6)) + scale_fill_viridis(name = "Number of \nbarcodes")
g.tiles

sumsynbcscodon <- group_by(R100bcs, X_peptide, Y_peptide, X_codon, Y_codon) %>%
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
