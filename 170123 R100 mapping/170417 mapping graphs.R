library(dplyr)
library(tidyr)
library(cowplot)
library(viridis)
library(extrafont)

setwd("/Users/Cliff/Documents/Kosuri Lab/NGS files/170123 R100 mapping")

R100bcs1 <- read.table("counts_Codon-1_S1_L001_R2_001.txt", sep = "", header = FALSE)
colnames(R100bcs1) <- c("count", "barcode")
R100bcs2 <- read.table("counts_Codon-2_S2_L001_R2_001.txt", sep = "", header = FALSE)
colnames(R100bcs2) <- c("count", "barcode")
R100bcs3 <- read.table("counts_Codon-3_S3_L001_R2_001.txt", sep = "", header = FALSE)
colnames(R100bcs3) <- c("count", "barcode")

all_bcs <- full_join(R100bcs1,R100bcs2, by = "barcode")
all_bcs <- full_join(all_bcs, R100bcs3, by = "barcode")
good_bcs <- filter(all_bcs, (is.na(count.x) & is.na(count.y)) | (is.na(count.x) & is.na(count)) | (is.na(count.y) & is.na(count)))
good_bcs$count.x[is.na(good_bcs$count.x)] <- 0
good_bcs$count.y[is.na(good_bcs$count.y)] <- 0
good_bcs$count[is.na(good_bcs$count)] <- 0
good_bcs <- mutate(good_bcs, "full_count" = count.x + count.y + count)
good_bcs <- select(good_bcs, barcode, full_count)


bad_bcs <- setdiff(all_bcs, good_bcs)
bad_bcs$count.x[is.na(bad_bcs$count.x)] <- 0
bad_bcs$count.y[is.na(bad_bcs$count.y)] <- 0
bad_bcs$count[is.na(bad_bcs$count)] <- 0
bad_bcs <- mutate(bad_bcs, "full_count" = count.x + count.y + count)

x <- c(0:50)
zpt <- (.35666^x)/((exp(.35666)-1)*factorial(x))*11314308
dzpt <- data.frame(x, zpt)
dzpt <- mutate(dzpt, zpt = zpt + 1)

g.bad <- ggplot(bad_bcs, aes(full_count)) + geom_histogram(bins = 100) + scale_x_continuous(limits = c(0, 500)) 
g.bad
g.good <- ggplot(good_bcs, aes(full_count)) + geom_histogram(bins = 51) + geom_point(data = dzpt, aes(x,zpt), color = "red", alpha = 0.5) + scale_x_continuous(limits = c(0,50)) +
  labs(x = "# of reads", y = "counts", title ="170418 Number of reads per barcode in mapping with ZPT model")
g.good
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
