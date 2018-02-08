library(dplyr)
library(tidyr)
library(cowplot)
library(viridis)
library(extrafont)

setwd("/Users/Cliff/Documents/Kosuri Lab/NGS files/170509 R100-2 mapping")

R100bcs <- read.csv("mapped_barcodes.csv")

R100bcs <- read.table("R100-c1.x-y.txt", sep = "", header = TRUE)
R100bcs <- bind_rows(R100bcs, read.table("R100-c2.x-y.txt", sep = "", header = TRUE))
R100bcs <- bind_rows(R100bcs, read.table("R100-c3.x-y.txt", sep = "", header = TRUE))

R100counts <- read.csv("R100-c1.map.csv", header = FALSE)
R100counts <- bind_rows(R100counts, read.csv("R100-c2.map.csv", header = FALSE))
R100counts <- bind_rows(R100counts, read.csv("R100-c3.map.csv", header = FALSE))

colnames(R100bcs) <- c("barcodes",  "X_peptide", "Y_peptide")
colnames(R100counts) <- c("barcodes", "DNA_string", "counts")

bconly1 <- read.table("counts_R100-c1_R2.txt", sep = "", header = FALSE)
bconly2 <- read.table("counts_R100-c2_R2.txt", sep = "", header = FALSE)
bconly3 <- read.table("counts_R100-c3_R2.txt", sep = "", header = FALSE)
colnames(bconly1) <- c("count", "barcode")
colnames(bconly2) <- c("count", "barcode")
colnames(bconly3) <- c("count", "barcode")

all_bcs <- full_join(bconly1,bconly2, by = "barcode")
all_bcs <- full_join(all_bcs, bconly3, by = "barcode")
good_bcs <- filter(all_bcs, (is.na(count.x) & is.na(count.y)) | (is.na(count.x) & is.na(count)) | (is.na(count.y) & is.na(count)))
good_bcs$count.x[is.na(good_bcs$count.x)] <- 0
good_bcs$count.y[is.na(good_bcs$count.y)] <- 0
good_bcs$count[is.na(good_bcs$count)] <- 0
good_bcs <- mutate(good_bcs, "full_count" = count.x + count.y + count)
good_bcs <- select(good_bcs, barcode, full_count)

R100bcs <- left_join(R100bcs, R100counts, by = "barcodes") %>%
  select(-DNA_string)
R100bcs <- filter(R100bcs, !is.na(X_peptide), !is.na(Y_peptide))
R100bcs <-  separate(R100bcs, X_peptide, into=c("X_peptide", "X_codon"), sep = "_rev_") %>%
  separate(Y_peptide, into=c("Y_peptide", "Y_codon"), sep = -3) 

g.bcs <- ggplot(R100bcs, aes(num_reads)) + geom_histogram(bins = 51) + scale_x_continuous(limits=c(0,50)) +# geom_point(data = dzpt, aes(x,zpt), color = "red", alpha = 0.5) +
  labs(x = "# of reads per barcode", title = "170510 R100-2 mapping reads per barcode blind match")
g.bcs

x <- c(1:50)
zpt <- (5.29^x)/((exp(5.29)-1)*factorial(x))*2425100
dzpt <- data.frame(x, zpt)
dzpt <- mutate(dzpt,x = x -1, zpt = zpt + 1)


Mode <- function(x) { 
  ux <-unique(x)
  ux[which.max(tabulate(match(x,ux)))]
  }
R100bcs2 <- R100bcs[!(duplicated(R100bcs$barcodes) | duplicated(R100bcs$barcodes, fromLast = TRUE)), ]
kims2 <- kimsbcs[!(duplicated(kimsbcs$barcodes) | duplicated(kimsbcs$barcodes, fromLast = TRUE)), ]

bothmaps <- full_join(kims2, R100bcs2, by = "barcodes")
bothmaps$counts[is.na(bothmaps$counts)] <- 0
bothmaps$num_reads[is.na(bothmaps$num_reads)] <- 0
 
g.natevkim <- ggplot(bothmaps, aes(counts, num_reads)) + geom_point(alpha = 0.3) + labs(x = "Nates mapping counts", y = "Kims mapping counts", title = "170511 Comparison of barcode counts between mapping methods")
g.natevkim

sumR100bcs <- group_by(R100bcs, X_peptide, Y_peptide) %>%
  summarise(numbcs = n(), totcounts = sum(num_reads), sdcount = sd(num_reads))

g.tiles <- ggplot(filter(sumR100bcs, numbcs > 5), aes(X_peptide, Y_peptide)) + geom_tile(aes(fill = numbcs)) + labs(title = "170509 Kim's Number of barcodes per X-Y pair", x = "X peptide", y = "Y peptide") +
  theme(text = element_text(family="Calibri"), axis.text.x = element_text(angle=90, size =6), axis.text.y = element_text(size=6)) + scale_fill_viridis(name = "Number of \nbarcodes")
g.tiles


g.tcnts <- ggplot(sumR100bcs, aes(X_peptide, Y_peptide)) + geom_tile(aes(fill = totcounts)) + labs(title = "170509 Number of counts per X-Y pair", x = "X peptide", y = "Y peptide") +
  theme(text = element_text(family="Calibri"), axis.text.x = element_text(angle=90, size =6), axis.text.y = element_text(size=6)) + scale_fill_viridis(name = "Total number\n of counts")
g.tcnts


g.bccounts <- ggplot(sumR100bcs, aes(numbcs)) + geom_histogram(bins =151) + labs(title = "170510 number of barcodes per X-Y pair", x = "number of barcodes") + scale_x_continuous(limits = c(0,150))
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

g.box <- ggplot(sumR100bcs, aes(X_peptide, numbcs))  + geom_boxplot(color="red", outlier.shape = NA) + geom_jitter(alpha = 0.2, size=1.5) + theme(text = element_text(family="Calibri"), axis.text.x = element_text(angle=90, size =6, hjust = .6)) +
  labs(title ="170509 Number of barcodes by X peptide", x = "X Peptide", y = "Number of BCs per Y Peptide") 
g.box

g.boxy <- ggplot(sumR100bcs, aes(Y_peptide, numbcs))  + geom_boxplot(color="red", outlier.shape = NA) + geom_jitter(alpha = 0.2, size=1.5) + theme(text = element_text(family="Calibri"), axis.text.x = element_text(angle=90, size =6, hjust = .6)) +
  labs(title ="170509 Number of barcodes by Y peptide", x = "Y Peptide", y = "Number of BCs per X Peptide") 
g.boxy
