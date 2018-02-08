library(dplyr)
library(tidyr)
library(cowplot)
library(viridis)
library(extrafont)

setwd("/Users/Cliff/Documents/Kosuri Lab/NGS files/161024 Synzip 3 codon mapping")

conReads = read.csv("161110 synzip maps.csv", header = FALSE)
colnames(conReads) <- c("barcode", "construct")

xref <- read.csv("synzip_x_ref.csv", header = FALSE)
colnames(xref) <- c("x_peptide","codon","sequence")
xbuff = "AGGAGAAGAGC"
xref$sequence <- toupper(xref$sequence)

xref$sequence <- substr(xref$sequence,16,length(xref$sequence))

test <-strsplit(as.character(xref$sequence),xbuff)
xref$sequence <- sapply(test, function(x) x[1])
xref <- separate(xref, x_peptide, into = c("X_peptide", "rev"), sep= "_")
xref <- select(xref, c(1,3))
xref <- distinct(xref)


