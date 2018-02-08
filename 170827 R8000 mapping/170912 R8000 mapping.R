library(dplyr)
library(tidyr)
library(cowplot)
library(viridis)

setwd("/Users/Cliff/Documents/Kosuri Lab/NGS Files/170827 R8000 mapping")

R8000 <- read.table("R8000-c1.x-y.txt", header = TRUE)
R8000 <- filter(R8000, !is.na(X), !is.na(Y))
R8000$X <- as.character(R8000$X)
R8000$Y <- as.character(R8000$Y)
R8000 <- filter(R8000, !(X != Y))


R8000 <- separate(R8000, X, into = c("X_peptide", "Y_peptide"), sep = "\\+")
R8000 <- separate(R8000, Y_peptide, into = c("Y_peptide", "Codon"), sep = "_")
R8000 <- separate(R8000, X_peptide, into = c("X_peptide", "X_group"), sep = "-")
R8000 <- separate(R8000, Y_peptide, into = c("Y_peptide", "Y_group"), sep = "-")
R8000$X_peptide <- as.factor(R8000$X_peptide)
R8000$Y_peptide <- as.factor(R8000$Y_peptide)

R8000$X_group[is.na(R8000$X_group) & grepl("mSN$",R8000$X_peptide)] <- "mSN"
R8000$X_group[is.na(R8000$X_group) & grepl("Hb$",R8000$X_peptide)] <- "Hb"
R8000$X_group[is.na(R8000$X_group) & grepl("Hq$",R8000$X_peptide)] <- "Hq"
R8000$X_group[is.na(R8000$X_group) & grepl("bH1$",R8000$X_peptide)] <- "bH1"
R8000$X_group[is.na(R8000$X_group) & grepl("mS$",R8000$X_peptide)] <- "mS"
R8000$X_group[is.na(R8000$X_group) & grepl("A$",R8000$X_peptide)] <- "A"
R8000$X_group[is.na(R8000$X_group) & grepl("SH$",R8000$X_peptide)] <- "SH"

sumR8000 <- group_by(R8000, X_peptide, Y_peptide, X_group) %>%
  summarise("totreads" = sum(Reads), "bcnum" = n())
gr1 <-c("4H1347","4H145","4H2619","4H2691","4H3070", "4H342","4H3632","4H496", "4H2745", "4H344", "4H363", "4H3907","4H3961","4H470","4H48","4H491")
gr2 <- c("4H1347","4H145","4H2619","4H2691","4H3070", "4H342","4H3632", "4H496")
gr3 <- c("4H1046","4H1072","4H2568","4H2619","4H2745","4H344","4H363","4H3907","4H3961","4H470")
gr4 <- c("4H1081","4H112","4H1158","4H118","4H170", "4H1831","4H2586", "4H2627", "4H2673", "4H2702", "4H3010", "4H3112", "4H3427", "4H346","4H3926","4H4073","4H457", "4H575", "4H683","4H754","4H780","4H797","4H855","4H918", "4H953", "4H968")

sumR8000 <- data.frame(sumR8000, "Hwrap")
g.bcnum <- ggplot(filter(sumR8000, X_peptide %in% gr4, Y_peptide %in% gr4), aes(X_peptide, Y_peptide)) + geom_tile(aes(fill = log10(bcnum))) + scale_fill_viridis() +
  theme(axis.text.x = element_text(angle = 90, size = 6), axis.text.y = element_text(size =6)) + facet_wrap(~X_group, scales = "free")
g.bcnum

g.bcnum <- ggplot(filter(sumR8000, X_group == "bA", grepl("4H",X_peptide)), aes(X_peptide, Y_peptide)) + geom_tile(aes(fill = log10(bcnum))) + scale_fill_viridis() +
  theme(axis.text.x = element_text(angle = 90, size =10), axis.text.y = element_text(size =10)) #+ facet_wrap(~X_group, scales = "free")
g.bcnum


g.bcnum.p <- ggplot(filter(sumR8000, grepl("P", X_peptide)), aes(X_peptide, Y_peptide)) + geom_tile(aes(fill = log10(bcnum))) + scale_fill_viridis(name = "Number of \n Barcodes \n log10") +
  theme(axis.text.x = element_text(angle = 90, size = 10), axis.text.y = element_text(size =10), strip.background = element_rect(fill = "white")) + facet_wrap(~X_group, scales = "free") +
  labs(x= "X Peptide", y = "Y Peptide", title = "170830 Number of Barcodes per P construct")
g.bcnum.p

o_peps2 <- filter(sumR8000, grepl("O", X_peptide))
o_peps$X_peptide <- as.factor(as.character(o_peps$X_peptide)) 
o_peps$Y_peptide <- as.factor(as.character(o_peps$Y_peptide)) 
o_peps$X_peptide <- factor(o_peps$X_peptide, c("O1"= "O1","O2" = "O2","O3" = "O3","O4" ="O4", "O5" = "O5", "O6" = "O6","O7" = "O7","O8" = "O8","O9" = "O9","O10" = "O10","O11" = "O11","O12" = "O12","O13"= "O13","O14"= "O14","O15" = "O15","O16" = "O16"))
o_peps$Y_peptide <- factor(o_peps$Y_peptide, rev(c("O1"= "O1","O2" = "O2","O3" = "O3","O4" ="O4", "O5" = "O5", "O6" = "O6","O7" = "O7","O8" = "O8","O9" = "O9","O10" = "O10","O11" = "O11","O12" = "O12","O13"= "O13","O14"= "O14","O15" = "O15","O16" = "O16")))

o_peps2 <- left_join(o_peps2, o_peps, by = c("X_peptide","Y_peptide", "X_group"))
g.check <- ggplot(o_peps2, aes(bcnum.x, bcnum.y)) + geom_point()
g.check

g.bcnum.o <- ggplot(o_peps, aes(X_peptide, Y_peptide)) + geom_tile(aes(fill = log10(bcnum))) + scale_fill_viridis(name = "Number of \n Barcodes \n log10") +
  theme(axis.text.x = element_text(angle = 90, size = 10), axis.text.y = element_text(size =10), strip.background = element_rect(fill = "white")) + facet_wrap(~X_group, scales = "free") +
  labs(x= "X Peptide", y = "Y Peptide", title = "170911 Number of Barcodes per O construct")
g.bcnum.o

g.bcnum.4h <- ggplot(filter(sumR8000, grepl("4H", X_peptide)), aes(X_peptide, Y_peptide)) + geom_tile(aes(fill = log10(bcnum))) + scale_fill_viridis(name = "Number of \n Barcodes \n log10") +
  theme(axis.text.x = element_text(angle = 90, size = 10), axis.text.y = element_text(size =10), strip.background = element_rect(fill = "white")) + facet_wrap(~X_group, scales = "free") +
  labs(x= "X Peptide", y = "Y Peptide", title = "170830 Number of Barcodes per 4H construct")
g.bcnum.4h

g.bcnum.g <- ggplot(filter(sumR8000, grepl("G", X_peptide)), aes(X_peptide, Y_peptide)) + geom_tile(aes(fill = log10(bcnum))) + scale_fill_viridis(name = "Number of \n Barcodes \n log10") +
  theme(axis.text.x = element_text(angle = 90, size = 10), axis.text.y = element_text(size =10), strip.background = element_rect(fill = "white")) + facet_wrap(~X_group, scales = "free") +
  labs(x= "X Peptide", y = "Y Peptide", title = "170830 Number of Barcodes per GCN construct")
g.bcnum.g

g.counts <- ggplot(R8000, aes(Reads)) + geom_histogram(bins = 126) + scale_x_continuous(limits = c(0,125)) + labs(x = "Reads per barcode", y = "histogram count", title = "170830 R8000 read depth per barcode" )
g.counts

g.bcper <- ggplot(sumR8000, aes(bcnum)) + geom_histogram(bins = 214) + scale_x_continuous(limits = c(0,213)) + labs(x = "Barcodes per construct pair", y = "histogram count", title = "170830 R8000 Barcodes per construct pair")
g.bcper

tempswit <- filter(R8000, (X != Y))

