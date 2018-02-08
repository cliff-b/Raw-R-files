################################################# Standard libraries + viridis to make things pretty + extrafont for comic sans #################################################
library(dplyr)
library(tidyr)
library(extrafont)
library(viridis)
library(cowplot)

################################################# Read in the data and mash it together #################################################
setwd("/Users/Cliff/Documents/Kosuri Lab/NGS files/170524 R100-2-2 Hiseq/")

R100_RNA_0 <- read.table("sc-counts-R-0h.txt", header = FALSE)
R100_RNA_0 <- data.frame(R100_RNA_0, rep("RNA-0",nrow(R100_RNA_0)),rep(sum(R100_RNA_0$V2), nrow(R100_RNA_0)))
colnames(R100_RNA_0) <- c("barcode", "count","collapsedBC", "condition","readsinlib")

R100_DNA_0 <- read.table("sc-counts-D-0h.txt", header = FALSE)
R100_DNA_0 <- data.frame(R100_DNA_0, rep("DNA-0",nrow(R100_DNA_0)),rep(sum(R100_DNA_0$V2), nrow(R100_DNA_0)))
colnames(R100_DNA_0) <- c("barcode","count", "collapsedBC", "condition","readsinlib")

R100_RNA_4 <- read.table("sc-counts-R-4h.txt", header = FALSE)
R100_RNA_4 <- data.frame(R100_RNA_4, rep("RNA-4",nrow(R100_RNA_4)),rep(sum(R100_RNA_4$V2), nrow(R100_RNA_4)))
colnames(R100_RNA_4) <- c("barcode", "count", "collapsedBC", "condition","readsinlib")

R100_DNA_4 <-read.table("sc-counts-D-4h.txt", header = FALSE)
R100_DNA_4 <- data.frame(R100_DNA_4, rep("DNA-4",nrow(R100_DNA_4)),rep(sum(R100_DNA_4$V2), nrow(R100_DNA_4)))
colnames(R100_DNA_4) <- c("barcode", "count", "collapsedBC", "condition","readsinlib")

R100_full_DF <- bind_rows(R100_RNA_0, R100_DNA_0, R100_RNA_4, R100_DNA_4)

R100_old_map <- read.csv("/Users/Cliff/Documents/Kosuri Lab/NGS files/170123 R100 mapping/170123 R100 mapped_barcodes.csv", skip = 2)
colnames(R100_old_map) <- c("X_peptide", "Y_peptide", "reads", "numbcs", "barcodes")
R100_old_map$barcodes <- as.character(R100_old_map$barcodes)
R100_old_map <- mutate(R100_old_map, barcode = strsplit(barcodes, ",")) %>%
  unnest(barcode)
R100_old_map <- separate(R100_old_map, X_peptide, into = c("X_peptide", "X_codon"), sep = "_rev_")
R100_old_map <- separate(R100_old_map, Y_peptide, into = c("Y_peptide", "Y_codon"), sep = -3)
R100_old_map <- separate(R100_old_map, X_peptide, into = c("Junk","X_peptide"), sep = 1)
R100_old_map$X_peptide[R100_old_map$X_peptide == "GCN4_1"] <- "GCN4"
R100_old_map <- separate(R100_old_map, Y_peptide, into = c("Junk2","Y_peptide"), sep = 1)
R100_old_map <- separate(R100_old_map, Y_peptide, into = c("Y_peptide", "Junk3"), sep = -2)
R100_old_map <- select(R100_old_map, X_peptide, Y_peptide, reads, numbcs, barcode)

R100_map1 <- read.table("/Users/Cliff/Documents/Kosuri Lab/NGS files/170524 R100-2-2 Hiseq/mapping/R100-c1.x-y.txt", skip = 2)
colnames(R100_map1) <- c("barcode", "X_peptide", "Y_peptide", "reads")
R100_map2 <- read.table("/Users/Cliff/Documents/Kosuri Lab/NGS files/170524 R100-2-2 Hiseq/mapping/R100-c2.x-y.txt", skip = 2)
colnames(R100_map2) <- c("barcode", "X_peptide", "Y_peptide", "reads")
R100_map3 <- read.table("/Users/Cliff/Documents/Kosuri Lab/NGS files/170524 R100-2-2 Hiseq/mapping/R100-c3.x-y.txt", skip = 2)
colnames(R100_map3) <- c("barcode", "X_peptide", "Y_peptide", "reads")
R100_map4 <- read.table("/Users/Cliff/Documents/Kosuri Lab/NGS files/170524 R100-2-2 Hiseq/mapping/R100-c4.x-y.txt", skip = 2)
colnames(R100_map4) <- c("barcode", "X_peptide", "Y_peptide", "reads")
R100_map5 <- read.table("/Users/Cliff/Documents/Kosuri Lab/NGS files/170524 R100-2-2 Hiseq/mapping/R100-c5.x-y.txt", skip = 2)
colnames(R100_map5) <- c("barcode", "X_peptide", "Y_peptide", "reads")
R100_map6 <- read.table("/Users/Cliff/Documents/Kosuri Lab/NGS files/170524 R100-2-2 Hiseq/mapping/R100-c6.x-y.txt", skip = 2)
colnames(R100_map6) <- c("barcode", "X_peptide", "Y_peptide", "reads")

R100_map <- bind_rows(R100_map1, R100_map2, R100_map3, R100_map4, R100_map5, R100_map6)
R100_raw_map <- select(R100_raw_map, -collapsedBC)
R100_map <- distinct(R100_map, barcode, X_peptide, Y_peptide)
R100_map <- filter(R100_map, X_peptide != "<NA>", Y_peptide != "<NA>")
R100_map <- separate(R100_map, X_peptide, into = c("X_peptide", "X_codon"), sep = "_rev_")
R100_map$X_peptide[R100_map$X_peptide == "GCN4_1"] <- "GCN4"
R100_map <- separate(R100_map, Y_peptide, into = c("Y_peptide", "Y_codon"), sep = -3)

load("/Users/Cliff/Documents/Kosuri Lab/NGS files/170509 R100-2 mapping/R100bcs.Rda")
R100_map <- filter(R100bcs, X_peptide != "<NA>", Y_peptide != "<NA>")
R100map <- separate(R100map, X_peptide, into = c("X_peptide", "X_codon"), sep = "_rev_")
R100map$X_peptide[R100map$X_peptide == "GCN4_1"] <- "GCN4"
R100map <- separate(R100map, Y_peptide, into = c("Y_peptide", "Y_codon"), sep = -3)
R100map <- select(R100map, X_peptide, Y_peptide, barcodes, "mapread" = counts)


R100map <- distinct(R100map, barcodes, X_peptide, Y_peptide)
R100map <- R100map[!(duplicated(R100map$barcodes) | duplicated(R100map$barcodes, fromLast = TRUE)),]
otherbcs <- R100_map[(duplicated(R100_map$barcodes) | duplicated(R100_map$barcodes, fromLast = TRUE)),]

R100_full_DF <- mutate(R100_full_DF, "nCount" = count/readsinlib*1000000)
R100_full_DF <- inner_join(R100_full_DF, R100_map, by = c("barcode"))
R100_full_DF <- select(R100_full_DF, -collapsedBC)

revComp <- function(x)
  chartr("ATGC","TACG", sapply(lapply(strsplit(x, NULL), rev), paste, collapse=""))
R100_full_DF <- data.frame(R100_full_DF, "RevCompBC" = revComp(as.character(R100_full_DF$barcode)))


R100RNA <- filter(R100_full_DF, grepl("RNA", R100_full_DF$condition)) %>%
  select(X_peptide, Y_peptide, condition, barcode, "rcount" = nCount ) %>%
  separate(condition, into = c("template", "condition"))
sumRNA <- group_by(R100RNA, X_peptide, Y_peptide, condition) %>%
  summarise(mean_RNA = mean(rcount), sd_R = sd(rcount), bcnum = n_distinct(barcode))

R100DNA <- filter(R100_full_DF, grepl("DNA", R100_full_DF$condition)) %>%
  select(X_peptide, Y_peptide, condition, barcode, "dcount" = nCount ) %>%
  separate(condition, into = c("template", "condition"))
sumDNA <- group_by(R100DNA, X_peptide, Y_peptide, condition) %>%
  summarise(mean_DNA = mean(dcount), sd_D = sd(dcount), bcnum_D = n_distinct(barcode))
R100RD <- right_join(R100RNA,R100DNA, by = c("X_peptide","Y_peptide","barcode","condition")) %>%
  mutate("RNADNA" = rcount/dcount*100) %>%
  select(-template.x, -template.y) %>%
  mutate("RNADNA" = ifelse(is.na(RNADNA),0,RNADNA)) %>%
  mutate("rcount" =  ifelse(is.na(rcount),0,rcount))


bccor <- sumR100RD
bccor$X_peptide <- as.factor(bccor$X_peptide)
bccor$Y_peptide <- as.factor(bccor$Y_peptide)
bccor$recip <- mutate(filter(bccor, Y_peptide == X_peptide, X_peptide == Y_peptide), recip = med_RD)

bccor <- data.frame(bccor, recip=rep(-1, nrow(bccor)))
for(i in 1:nrow(bccor)){
  if(as.integer(bccor$X_peptide[i]) <= 113 - as.integer(bccor$Y_peptide[i])) {
  a <- filter(bccor, X_peptide==bccor[i,2],Y_peptide == bccor[i,1], condition == bccor[i,3])
  print(a)
    if(nrow(a)) {
      bccor[i,13] <- a$med_RD
    }
  }
}

bccor1 <- filter(bccor, recip != -1)
g.cor <- ggplot(bccor1, aes(med_RD, recip, color = condition)) + geom_point(alpha = 0.3) + scale_x_log10() + scale_y_log10() + labs(x = "median RNA/DNA", y = "reciprocal protein pair RNA/DNA", title = "170601 R100 median RNA/DNA against switch X and Y median RNA/DNA")
g.cor


sumR100RD <- group_by(R100RD, X_peptide, Y_peptide, condition) %>%
  summarise(mean_RNA = mean(rcount), sd_R = sd(rcount), mean_DNA = mean(dcount), sd_D = sd(dcount), mean_RD = mean(RNADNA), sd_RD = sd(RNADNA), med_RD = median(RNADNA), bcnum = n_distinct(barcode), SEM = sd(RNADNA)/sqrt(n()))



induction_labels <- c("0" = "0 hours induction RNA/DNA", "4" = "4 hours induction RNA/DNA")
induction_labels2 <- c("0" = "0 hours induction DNA", "4" = "4 hours induction DNA")
induction_labels3 <- c("0" = "0 hours induction RNA", "4" = "4 hours induction RNA")

g.sumR100RD <- ggplot(sumR100RD, aes(X_peptide, Y_peptide)) + geom_tile(aes(fill=mean_RD)) + facet_wrap(~condition, labeller = labeller(condition=induction_labels), scales = "free") +
  theme(text=element_text(family="Calibri"), strip.text=element_text(size = 10), strip.background=element_rect(fill = "White"), axis.text.x = element_text(angle=90,size =5),
        axis.text.y = element_text(size = 5), legend.title=element_text(size = 8), legend.title.align=0.5, legend.text=element_text(size =8)) +
  scale_fill_viridis(option="viridis", name = "Mean \nRNA to DNA  ", na.value="white", limits=c(0, 770)) + labs(title="Roman's 100 RNA/DNA barcode ratio", x="X peptide", y = "Y peptide")
g.sumR100RD

g.bccount <- ggplot(R100_full_DF, aes(count)) + geom_histogram() + labs(x = "Number of reads", y = "number of barcodes", title = "170525 R100-2-2 Hiseq Barcode counts") + scale_x_continuous(limits = c(0, 1000))
g.bccount

R100_raw_map <- full_join(R100_raw_map, R100_DNA_0, by ="barcode")
R100_raw_map$count.x[is.na(R100_raw_map$count.x)] <- 1
R100_raw_map$count.y[is.na(R100_raw_map$count.y)] <- 1

g.mapcount <- ggplot(R100_raw_map, aes(count.y, count.x)) + geom_point(alpha = 0.1) + scale_x_log10() + scale_y_log10() + labs(x = "Barcode-seq DNA barcode counts", y = "Mapping run barcode counts", title = "170526 comparison of barcode counts >=0")
g.mapcount
