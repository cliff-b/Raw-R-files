################################################# Standard libraries + viridis to make things pretty + extrafont for comic sans #################################################
library(dplyr)
library(tidyr)
library(extrafont)
library(viridis)
library(cowplot)
################################################# Read in the data and mash it together #################################################

setwd("/Users/Cliff/Google Drive/NGS files/171218 R8000-4/")

R8000_RNA_0 <- read.table("sc-counts-R8000-4_0h_RNA.txt", header = FALSE)
R8000_RNA_0 <- data.frame(R8000_RNA_0, rep("RNA-0",nrow(R8000_RNA_0)),rep(sum(R8000_RNA_0$V2), nrow(R8000_RNA_0)))
colnames(R8000_RNA_0) <- c("barcode", "count", "collapsed_bc", "condition","readsinlib")
R8000_RNA_0 <- mutate(R8000_RNA_0, count =  count/norm)

R8000_DNA_0 <- read.table("sc-counts-R8000-4_0h_DNA.txt", header = FALSE)
R8000_DNA_0 <- data.frame(R8000_DNA_0, rep("DNA-0",nrow(R8000_DNA_0)),rep(sum(R8000_DNA_0$V2), nrow(R8000_DNA_0)))
colnames(R8000_DNA_0) <- c("barcode", "count", "collapsed_bc", "condition","readsinlib")

R8000_RNA_6 <- read.table("sc-counts-R8000-4_6h_RNA.txt", header = FALSE)
R8000_RNA_6 <- data.frame(R8000_RNA_6, rep("RNA-6",nrow(R8000_RNA_6)),rep(sum(R8000_RNA_6$V2), nrow(R8000_RNA_6)))
colnames(R8000_RNA_6) <- c("barcode", "count", "collapsed_bc", "condition","readsinlib")

R8000_DNA_6 <-read.table("sc-counts-R8000-4_6h_DNA.txt", header = FALSE)
R8000_DNA_6 <- data.frame(R8000_DNA_6, rep("DNA-6",nrow(R8000_DNA_6)),rep(sum(R8000_DNA_6$V2), nrow(R8000_DNA_6)))
colnames(R8000_DNA_6) <- c("barcode", "count", "collapsed_bc", "condition","readsinlib")


R8000_full_DF <- bind_rows(R8000_RNA_0, R8000_DNA_0, R8000_RNA_6, R8000_DNA_6)


R8000_map <- read.table("R8000-c1.x-y.txt", skip = 1)
colnames(R8000_map) <- c("barcode", "X_peptide", "Y_peptide", "reads")
R8000_map <- filter(R8000_map, X_peptide != "<NA>", Y_peptide != "<NA>") %>%
  filter(as.character(X_peptide) == as.character(Y_peptide)) %>%
  select(-Y_peptide) %>%
  separate(X_peptide, into = c("X_peptide", "X_codon"), sep = "_") %>%
  separate(X_peptide, into = c("X_peptide","Y_peptide"), sep = "\\+") %>%
  separate(X_peptide, into = c("X_peptide", "X_group"), sep = "-") %>%
  separate(Y_peptide, into = c("Y_peptide", "Y_group"), sep = "-")

R8000_full_DF <- mutate(R8000_full_DF, "nCount" = count/readsinlib*1000000)
R8000_full_DF <- inner_join(R8000_full_DF, R8000_map, by = c("barcode"))


R8000RNA <- filter(R8000_full_DF, grepl("RNA", R8000_full_DF$condition)) %>%
  select(X_peptide, X_group, X_codon, Y_peptide, condition, barcode, "rcount" = nCount ) %>%
  separate(condition, into = c("template", "condition"))
sumRNA <- group_by(R100RNA, X_peptide, Y_peptide, condition) %>%
  summarise(mean_RNA = mean(rcount), sd_R = sd(rcount), bcnum = n_distinct(barcode))

R8000DNA <- filter(R8000_full_DF, grepl("DNA", R8000_full_DF$condition)) %>%
  select(X_peptide, X_group, X_codon, Y_peptide, condition, barcode, "dcount" = nCount ) %>%
  separate(condition, into = c("template", "condition"))
sumDNA <- group_by(R100DNA, X_peptide, Y_peptide, condition) %>%
  summarise(mean_DNA = mean(dcount), sd_D = sd(dcount), bcnum_D = n_distinct(barcode))
R8000RD <- right_join(R8000RNA,R8000DNA, by = c("X_peptide","Y_peptide","barcode","condition", "X_group", "X_codon")) %>%
  mutate("RNADNA" = rcount/dcount*100) %>%
  select(-template.x, -template.y) %>%
  mutate("RNADNA" = ifelse(is.na(RNADNA),0,RNADNA)) %>%
  mutate("rcount" =  ifelse(is.na(rcount),0,rcount))

sumR100RD <- group_by(R8000RD, X_peptide, Y_peptide, condition, X_group) %>%
  summarise(mean_RNA = mean(rcount), sd_R = sd(rcount), mean_DNA = mean(dcount), sd_D = sd(dcount), mean_RD = mean(RNADNA), sd_RD = sd(RNADNA), med_RD = median(RNADNA), sum_DNA = sum(dcount), sum_RNA = sum(rcount), sum_RD = sum(rcount)/(sum(dcount) +1), bcnum = n_distinct(barcode), SEM = sd(RNADNA)/sqrt(n()))

induction_labels <- c("0" = "0 hours induction RNA/DNA", "6" = "6 hours induction RNA/DNA")

g.sumOs <- ggplot(filter(sumR100RD, grepl("O",X_peptide)), aes(X_peptide, Y_peptide)) + geom_tile(aes(fill=med_RD)) + facet_grid(X_group~condition, labeller = labeller(condition=induction_labels), scales = "free") +
  theme(text=element_text(family="Calibri"), strip.text=element_text(size = 10), strip.background=element_rect(fill = "White"), axis.text.x = element_text(angle=90,size =5),
        axis.text.y = element_text(size = 5), legend.title=element_text(size = 8), legend.title.align=0.5, legend.text=element_text(size =8)) +
  scale_fill_viridis(option="viridis", name = "Med \nRNA to DNA  ", na.value="white") + labs(title="171218 Roman's 8000 RNA/DNA mean barcode ratio O constructs (bc >0) repeat", x="X peptide", y = "Y peptide")
g.sumOs

g.sumPs <- ggplot(filter(sumR100RD, grepl("P",X_peptide), !is.na(X_group)), aes(X_peptide, Y_peptide)) + geom_tile(aes(fill=med_RD)) + facet_grid(X_group~condition, labeller = labeller(condition=induction_labels)) +
  theme(text=element_text(family="Calibri"), strip.text=element_text(size = 10), strip.background=element_rect(fill = "White"), axis.text.x = element_text(angle=90,size =5),
        axis.text.y = element_text(size = 5), legend.title=element_text(size = 8), legend.title.align=0.5, legend.text=element_text(size =8)) +
  scale_fill_viridis(option="viridis", name = "Med \nRNA to DNA  ", na.value="white") + labs(title="171218 Roman's 8000 RNA/DNA med barcode ratio P constructs (bc >0) repeat", x="X peptide", y = "Y peptide")
g.sumPs



g.bccount <- ggplot(sumR100RD, aes(bcnum, fill = mean_RD)) + geom_histogram(bins = 101) + scale_x_continuous(limits = c(0,100)) + labs(title = "171218 R8000-4 Number of barcodes per construct", x = "number of barcodes")  +
  scale_fill_viridis()
g.bccount


setwd("/Users/Cliff/Documents/Kosuri Lab/NGS files/170926 R8000 barcode-1/")
load_obj <- function(f)
{
  load(f)
  get(ls()[ls() != "f"])
}
R8000_1_full <- load_obj("/Users/Cliff/Google Drive/NGS files/170926 R8000 barcode-1/R8000_full_DF.Rda")
R8000_1_full <- mutate(R8000_1_full, "nCount" = count/readsinlib*1000000)
R8000_1_full <- inner_join(R8000_1_full, R8000_map, by = "barcode")
R8000_1_RNA <- filter(R8000_1_full, grepl("RNA", R8000_1_full$condition)) %>%
  select(X_peptide, X_group, X_codon, Y_peptide, condition, barcode, "rcount" = nCount ) %>%
  separate(condition, into = c("template", "condition"))
R8000_1_DNA <- filter(R8000_1_full, grepl("DNA", R8000_1_full$condition)) %>%
  select(X_peptide, X_group, X_codon, Y_peptide, condition, barcode, "dcount" = nCount ) %>%
  separate(condition, into = c("template", "condition"))
R8000_1_RD <- right_join(R8000_1_RNA,R8000_1_DNA, by = c("X_peptide","Y_peptide","barcode","condition", "X_group", "X_codon")) %>%
  mutate("RNADNA" = rcount/dcount*100) %>%
  select(-template.x, -template.y) %>%
  mutate("RNADNA" = ifelse(is.na(RNADNA),0,RNADNA)) %>%
  mutate("rcount" =  ifelse(is.na(rcount),0,rcount))

insane_R8000 <- full_join(R8000RD, R8000_1_RD, by = c("condition", "barcode", "X_peptide","X_group","X_codon","Y_peptide"))
colnames(insane_R8000) <- c("X_peptide", "X_group", "X_codon", "Y_peptide", "condition", "barcode", "rcount-4","dcount-4","RNADNA-4","rcount-1","dcount-1","RNADNA-1")

R8000_3_RNA0 <- read.table("/Users/Cliff/Google Drive/NGS files/171113 R8000-3/sc-counts_R8000-3_0h_RNA.txt", header = FALSE)
R8000_3_RNA0 <- data.frame(R8000_3_RNA0, rep("RNA-0",nrow(R8000_3_RNA0)),rep(sum(R8000_3_RNA0$V2), nrow(R8000_3_RNA0)))
colnames(R8000_3_RNA0) <- c("barcode", "count", "collapsed_bc", "condition","readsinlib")

R8000_3_DNA0 <- read.table("/Users/Cliff/Google Drive/NGS files/171113 R8000-3/sc-counts_R8000-3_0h_DNA.txt", header = FALSE)
R8000_3_DNA0 <- data.frame(R8000_3_DNA0, rep("DNA-0",nrow(R8000_3_DNA0)),rep(sum(R8000_3_DNA0$V2), nrow(R8000_3_DNA0)))
colnames(R8000_3_DNA0) <- c("barcode", "count", "collapsed_bc", "condition","readsinlib")

R8000_3_RNA6 <- read.table("/Users/Cliff/Google Drive/NGS files/171113 R8000-3/sc-counts_R8000-3_6h_RNA.txt", header = FALSE)
R8000_3_RNA6 <- data.frame(R8000_3_RNA6, rep("RNA-6",nrow(R8000_3_RNA6)),rep(sum(R8000_3_RNA6$V2), nrow(R8000_3_RNA6)))
colnames(R8000_3_RNA6) <- c("barcode", "count", "collapsed_bc", "condition","readsinlib")

R8000_3_DNA6 <- read.table("/Users/Cliff/Google Drive/NGS files/171113 R8000-3/sc-counts_R8000-3_6h_DNA.txt", header = FALSE)
R8000_3_DNA6 <- data.frame(R8000_3_DNA6, rep("DNA-6",nrow(R8000_3_DNA6)),rep(sum(R8000_3_DNA6$V2), nrow(R8000_3_DNA6)))
colnames(R8000_3_DNA6) <- c("barcode", "count", "collapsed_bc", "condition","readsinlib")

R8000_3_full <- bind_rows(R8000_3_RNA0, R8000_3_DNA0, R8000_3_RNA6, R8000_3_DNA6)

R8000_3_full <- mutate(R8000_3_full, "nCount" = count/readsinlib*1000000)
R8000_3_full <- inner_join(R8000_3_full, R8000_map, by = "barcode")
R8000_3_RNA <- filter(R8000_3_full, grepl("RNA", R8000_3_full$condition)) %>%
  select(X_peptide, X_group, X_codon, Y_peptide, condition, barcode, "rcount" = nCount ) %>%
  separate(condition, into = c("template", "condition"))
R8000_3_DNA <- filter(R8000_3_full, grepl("DNA", R8000_3_full$condition)) %>%
  select(X_peptide, X_group, X_codon, Y_peptide, condition, barcode, "dcount" = nCount ) %>%
  separate(condition, into = c("template", "condition"))
R8000_3_RD <- right_join(R8000_3_RNA,R8000_3_DNA, by = c("X_peptide","Y_peptide","barcode","condition", "X_group", "X_codon")) %>%
  mutate("RNADNA" = rcount/dcount*100) %>%
  select(-template.x, -template.y) %>%
  mutate("RNADNA" = ifelse(is.na(RNADNA),0,RNADNA)) %>%
  mutate("rcount" =  ifelse(is.na(rcount),0,rcount))

insane_R8000 <- full_join(insane_R8000, R8000_3_RD, by = c("condition", "barcode", "X_peptide","X_group","X_codon","Y_peptide"))
colnames(insane_R8000) <- c("X_peptide", "X_group", "X_codon", "Y_peptide", "condition", "barcode", "rcount_4","dcount_4","RNADNA_4","rcount_1","dcount_1","RNADNA_1", "rcount_3","dcount_3","RNADNA_3")
insane_R8000$rcount_4[is.na(insane_R8000$rcount_4)] <- 0
insane_R8000$dcount_4[is.na(insane_R8000$dcount_4)] <- 0
insane_R8000$RNADNA_4[is.na(insane_R8000$RNADNA_4)] <- 0
insane_R8000$rcount_1[is.na(insane_R8000$rcount_1)] <- 0
insane_R8000$dcount_1[is.na(insane_R8000$dcount_1)] <- 0
insane_R8000$RNADNA_1[is.na(insane_R8000$RNADNA_1)] <- 0
insane_R8000$rcount_3[is.na(insane_R8000$rcount_3)] <- 0
insane_R8000$dcount_3[is.na(insane_R8000$dcount_3)] <- 0
insane_R8000$RNADNA_3[is.na(insane_R8000$RNADNA_3)] <- 0



sumInsane <- group_by(insane_R8000, X_peptide, Y_peptide, condition, X_group) %>%
  summarise(mean_RNA4 = mean(rcount_4), mean_DNA4 = mean(dcount_4), mean_RD4 = mean(RNADNA_4), sd_RD4 = sd(RNADNA_4), med_RD4 = median(RNADNA_4), sum_DNA4 = sum(dcount_4), sum_RNA4 = sum(rcount_4), sum_RD4 = sum(rcount_4)/(sum(dcount_4) +0.1), mean_RNA1 = mean(rcount_1), mean_DNA1 = mean(dcount_1), mean_RD1 = mean(RNADNA_1), sd_RD1 = sd(RNADNA_1), med_RD1 = median(RNADNA_1), sum_DNA1 = sum(dcount_1), sum_RNA1 = sum(rcount_1), sum_RD1 = sum(rcount_1)/(sum(dcount_1) +0.1), mean_RNA3 = mean(rcount_3), mean_DNA3 = mean(dcount_3), mean_RD3 = mean(RNADNA_3), sd_RD3 = sd(RNADNA_3), med_RD3 = median(RNADNA_3), sum_DNA3 = sum(dcount_3), sum_RNA3 = sum(rcount_3), sum_RD3 = sum(rcount_3)/(sum(dcount_3) +0.1), bcnum = n_distinct(barcode))



wrap <- gather(insane_R8000, library, counts, rcount_4, dcount_4, RNADNA_4, rcount_1, dcount_1, RNADNA_1, rcount_3, dcount_3, RNADNA_3)

#######DNA at time zero comparisons###############
corr_d14_0 <- data.frame(x = 5, y = 25, label =  paste("italic(r)==", round(cor(filter(insane_R8000, dcount_1 > 1, dcount_4 > 1, condition == 0)$dcount_1, filter(insane_R8000, dcount_1 > 1, dcount_4 > 1, condition == 0)$dcount_4, use = "complete.obs"), digits = 3)))
g.dna14_0 <- ggplot(filter(insane_R8000, dcount_4 > 1, dcount_1 >1, condition == 0), aes(dcount_1, dcount_4)) +geom_point(alpha=0.1) + scale_x_log10() + scale_y_log10() + labs(x = "DNA Replicate 1", y = "DNA Replicate 4") + geom_text( data = corr_d14_0, aes(x= x,y = y,label= label), parse = TRUE) #+ facet_wrap(condition~library, scales = "free")
g.dna14_0
corr_d13_0 <- data.frame(x = 5, y = 25, label =  paste("italic(r)==", round(cor(filter(insane_R8000, dcount_1 > 1, dcount_3 > 1, condition == 0)$dcount_1, filter(insane_R8000, dcount_1 > 1, dcount_3 > 1, condition == 0)$dcount_3, use = "complete.obs"), digits = 3)))
g.dna13_0 <- ggplot(filter(insane_R8000, dcount_3 > 1, dcount_1 >1, condition == 0), aes(dcount_1, dcount_3)) +geom_point(alpha=0.1) + scale_x_log10() + scale_y_log10() + labs(x = "DNA Replicate 1", y = "DNA Replicate 3") + geom_text( data = corr_d13_0, aes(x= x,y = y,label= label), parse = TRUE) #+ facet_wrap(condition~library, scales = "free")
g.dna13_0
corr_d34_0 <- data.frame(x = 5, y = 25, label =  paste("italic(r)==", round(cor(filter(insane_R8000, dcount_4 > 1, dcount_3 > 1, condition == 0)$dcount_4, filter(insane_R8000, dcount_4 > 1, dcount_3 > 1, condition == 0)$dcount_3, use = "complete.obs"), digits = 3)))
g.dna34_0 <- ggplot(filter(insane_R8000, dcount_3 > 1, dcount_4 >1, condition == 0), aes(dcount_3, dcount_4)) +geom_point(alpha=0.1) + scale_x_log10() + scale_y_log10() + labs(x = "DNA Replicate 3", y = "DNA Replicate 4") + geom_text( data = corr_d34_0, aes(x= x,y = y,label= label), parse = TRUE) #+ facet_wrap(condition~library, scales = "free")
g.dna34_0

#######RNA at time zero comparisons###############
corr_r14_0 <- data.frame(x = 2, y = 25, label =  paste("italic(r)==", round(cor(filter(insane_R8000, rcount_1 > 0.1, rcount_4 > 0.1, condition == 0)$rcount_1, filter(insane_R8000, rcount_1 > 0.1, rcount_4 > 0.1, condition == 0)$rcount_4, use = "complete.obs"), digits = 3)))
g.rna14_0 <- ggplot(filter(insane_R8000, rcount_4 > 0.1, rcount_1 >0.1, condition == 0), aes(rcount_1, rcount_4)) +geom_point(alpha=0.1) + scale_x_log10() + scale_y_log10() + labs(x = "RNA Replicate 1", y = "RNA Replicate 4") + geom_text( data = corr_r14_0, aes(x= x,y = y,label= label), parse = TRUE) #+ facet_wrap(condition~library, scales = "free")
g.rna14_0
corr_r13_0 <- data.frame(x = 1, y = 25, label =  paste("italic(r)==", round(cor(filter(insane_R8000, rcount_1 > 0.1, rcount_3 > 0.1, condition == 0)$rcount_1, filter(insane_R8000, rcount_1 > 0.1, rcount_3 > 0.1, condition == 0)$rcount_3, use = "complete.obs"), digits = 3)))
g.rna13_0 <- ggplot(filter(insane_R8000, rcount_3 > 0.1, rcount_1 > 0.1, condition == 0), aes(rcount_1, rcount_3)) +geom_point(alpha=0.1) + scale_x_log10() + scale_y_log10() + labs(x = "RNA Replicate 1", y = "RNA Replicate 3") + geom_text( data = corr_r13_0, aes(x= x,y = y,label= label), parse = TRUE) #+ facet_wrap(condition~library, scales = "free")
g.rna13_0
corr_r34_0 <- data.frame(x = 1, y = 25, label =  paste("italic(r)==", round(cor(filter(insane_R8000, rcount_4 > 0.1, rcount_3 > 0.1, condition == 0)$rcount_4, filter(insane_R8000, rcount_4 > 0.1, rcount_3 > 0.1, condition == 0)$rcount_3, use = "complete.obs"), digits = 3)))
g.rna34_0 <- ggplot(filter(insane_R8000, rcount_3 > 0.1, rcount_4 > 0.1, condition == 0), aes(rcount_3, rcount_4)) +geom_point(alpha=0.1) + scale_x_log10() + scale_y_log10() + labs(x = "RNA Replicate 3", y = "RNA Replicate 4") + geom_text( data = corr_r34_0, aes(x= x,y = y,label= label), parse = TRUE) #+ facet_wrap(condition~library, scales = "free")
g.rna34_0

bob <- plot_grid(NULL, g.dna14_0, g.dna34_0, g.rna14_0, NULL, g.dna13_0, g.rna34_0, g.rna13_0, NULL, ncol= 3)
ggdraw(bob) + draw_label("R8000 Biological replicates at time 0", x = 0.5, y =1, vjust = 1, hjust = 0.5, size =20)


#######DNA at time six comparisons###############
corr_d14_6 <- data.frame(x = 5, y = 25, label =  paste("italic(r)==", round(cor(filter(insane_R8000, dcount_1 > 1, dcount_4 > 1, condition == 6)$dcount_1, filter(insane_R8000, dcount_1 > 1, dcount_4 > 1, condition == 6)$dcount_4, use = "complete.obs"), digits = 3)))
g.dna14_6 <- ggplot(filter(insane_R8000, dcount_4 > 1, dcount_1 >1, condition == 6), aes(dcount_1, dcount_4)) +geom_point(alpha=0.1) + scale_x_log10() + scale_y_log10() + labs(x = "DNA Replicate 1", y = "DNA Replicate 4") + geom_text( data = corr_d14_6, aes(x= x,y = y,label= label), parse = TRUE) #+ facet_wrap(condition~library, scales = "free")
g.dna14_6
corr_d13_6 <- data.frame(x = 5, y = 25, label =  paste("italic(r)==", round(cor(filter(insane_R8000, dcount_1 > 1, dcount_3 > 1, condition == 6)$dcount_1, filter(insane_R8000, dcount_1 > 1, dcount_3 > 1, condition == 6)$dcount_3, use = "complete.obs"), digits = 3)))
g.dna13_6 <- ggplot(filter(insane_R8000, dcount_3 > 1, dcount_1 >1, condition == 6), aes(dcount_1, dcount_3)) +geom_point(alpha=0.1) + scale_x_log10() + scale_y_log10() + labs(x = "DNA Replicate 1", y = "DNA Replicate 3") + geom_text( data = corr_d13_6, aes(x= x,y = y,label= label), parse = TRUE) #+ facet_wrap(condition~library, scales = "free")
g.dna13_6
corr_d34_6 <- data.frame(x = 5, y = 25, label =  paste("italic(r)==", round(cor(filter(insane_R8000, dcount_4 > 1, dcount_3 > 1, condition == 6)$dcount_4, filter(insane_R8000, dcount_4 > 1, dcount_3 > 1, condition == 6)$dcount_3, use = "complete.obs"), digits = 3)))
g.dna34_6 <- ggplot(filter(insane_R8000, dcount_3 > 1, dcount_4 >1, condition == 6), aes(dcount_3, dcount_4)) +geom_point(alpha=0.1) + scale_x_log10() + scale_y_log10() + labs(x = "DNA Replicate 3", y = "DNA Replicate 4") + geom_text( data = corr_d34_6, aes(x= x,y = y,label= label), parse = TRUE) #+ facet_wrap(condition~library, scales = "free")
g.dna34_6

#######RNA at time six comparisons###############
corr_r14_6 <- data.frame(x = 1, y = 25, label =  paste("italic(r)==", round(cor(filter(insane_R8000, rcount_1 > 0.1, rcount_4 > 0.1, condition == 6)$rcount_1, filter(insane_R8000, rcount_1 > 0.1, rcount_4 > 0.1, condition == 6)$rcount_4, use = "complete.obs"), digits = 3)))
g.rna14_6 <- ggplot(filter(insane_R8000, rcount_4 > 0.1, rcount_1 >0.1, condition == 6), aes(rcount_1, rcount_4)) +geom_point(alpha=0.1) + scale_x_log10() + scale_y_log10() + labs(x = "RNA Replicate 1", y = "RNA Replicate 4") + geom_text( data = corr_r14_6, aes(x= x,y = y,label= label), parse = TRUE) #+ facet_wrap(condition~library, scales = "free")
g.rna14_6
corr_r13_6 <- data.frame(x = 1, y = 25, label =  paste("italic(r)==", round(cor(filter(insane_R8000, rcount_1 > 0.1, rcount_3 > 0.1, condition == 6)$rcount_1, filter(insane_R8000, rcount_1 > 0.1, rcount_3 > 0.1, condition == 6)$rcount_3, use = "complete.obs"), digits = 3)))
g.rna13_6 <- ggplot(filter(insane_R8000, rcount_3 > 0.1, rcount_1 > 0.1, condition == 6), aes(rcount_1, rcount_3)) +geom_point(alpha=0.1) + scale_x_log10() + scale_y_log10() + labs(x = "RNA Replicate 1", y = "RNA Replicate 3") + geom_text( data = corr_r13_6, aes(x= x,y = y,label= label), parse = TRUE) #+ facet_wrap(condition~library, scales = "free")
g.rna13_6
corr_r34_6 <- data.frame(x = 1, y = 25, label =  paste("italic(r)==", round(cor(filter(insane_R8000, rcount_4 > 0.1, rcount_3 > 0.1, condition == 6)$rcount_4, filter(insane_R8000, rcount_4 > 0.1, rcount_3 > 0.1, condition == 6)$rcount_3, use = "complete.obs"), digits = 3)))
g.rna34_6 <- ggplot(filter(insane_R8000, rcount_3 > 0.1, rcount_4 > 0.1, condition == 6), aes(rcount_3, rcount_4)) +geom_point(alpha=0.1) + scale_x_log10() + scale_y_log10() + labs(x = "RNA Replicate 3", y = "RNA Replicate 4") + geom_text( data = corr_r34_6, aes(x= x,y = y,label= label), parse = TRUE) #+ facet_wrap(condition~library, scales = "free")
g.rna34_6

tom <- plot_grid(NULL, g.dna14_6, g.dna34_6, g.rna14_6, NULL, g.dna13_6, g.rna34_6, g.rna13_6, NULL, ncol= 3)
ggdraw(tom) + draw_label("R8000 Biological replicates at time 6h", x = 0.5, y =1, vjust = 1, hjust = 0.5, size =20)

#######Summed DNA at time zero comparisons###############
corr_sd14_0 <- data.frame(x = 5, y = 25, label =  paste("italic(r)==", round(cor(filter(sumInsane, sum_DNA1 > 1, sum_DNA4 > 1, condition == 0)$sum_DNA1, filter(sumInsane, sum_DNA1 > 1, sum_DNA4 > 1, condition == 0)$sum_DNA4, use = "complete.obs"), digits = 3)))
g.sumdna14_0 <- ggplot(filter(sumInsane, sum_DNA4 > 1, sum_DNA1 >1, condition == 0), aes(sum_DNA1, sum_DNA4)) +geom_point(alpha=0.1) + scale_x_log10() + scale_y_log10() + labs(x = "summed DNA Replicate 1", y = "summed DNA Replicate 4") + geom_text( data = corr_sd14_0, aes(x= x,y = y,label= label), parse = TRUE) #+ facet_wrap(condition~library, scales = "free")
g.sumdna14_0
corr_sd13_0 <- data.frame(x = 5, y = 25, label =  paste("italic(r)==", round(cor(filter(sumInsane, sum_DNA1 > 1, sum_DNA3 > 1, condition == 0)$sum_DNA1, filter(sumInsane, sum_DNA1 > 1, sum_DNA3 > 1, condition == 0)$sum_DNA3, use = "complete.obs"), digits = 3)))
g.sumdna13_0 <- ggplot(filter(sumInsane, sum_DNA3 > 1, sum_DNA1 >1, condition == 0), aes(sum_DNA1, sum_DNA3)) +geom_point(alpha=0.1) + scale_x_log10() + scale_y_log10() + labs(x = "summed DNA Replicate 1", y = "summed DNA Replicate 3") + geom_text( data = corr_sd13_0, aes(x= x,y = y,label= label), parse = TRUE) #+ facet_wrap(condition~library, scales = "free")
g.sumdna13_0
corr_sd34_0 <- data.frame(x = 5, y = 25, label =  paste("italic(r)==", round(cor(filter(sumInsane, sum_DNA4 > 1, sum_DNA3 > 1, condition == 0)$sum_DNA4, filter(sumInsane, sum_DNA4 > 1, sum_DNA3 > 1, condition == 0)$sum_DNA3, use = "complete.obs"), digits = 3)))
g.sumdna34_0 <- ggplot(filter(sumInsane, sum_DNA3 > 1, sum_DNA4 >1, condition == 0), aes(sum_DNA3, sum_DNA4)) +geom_point(alpha=0.1) + scale_x_log10() + scale_y_log10() + labs(x = "summed DNA Replicate 3", y = "summed DNA Replicate 4") + geom_text( data = corr_sd34_0, aes(x= x,y = y,label= label), parse = TRUE) #+ facet_wrap(condition~library, scales = "free")
g.sumdna34_0

#######RNA at time zero comparisons###############
corr_sr14_0 <- data.frame(x = 2, y = 25, label =  paste("italic(r)==", round(cor(filter(sumInsane, sum_RNA1 > 0.1, sum_RNA4 > 0.1, condition == 0)$sum_RNA1, filter(sumInsane, sum_RNA1 > 0.1, sum_RNA4 > 0.1, condition == 0)$sum_RNA4, use = "complete.obs"), digits = 3)))
g.sumrna14_0 <- ggplot(filter(sumInsane, sum_RNA4 > 0.1, sum_RNA1 >0.1, condition == 0), aes(sum_RNA1, sum_RNA4)) +geom_point(alpha=0.1) + scale_x_log10() + scale_y_log10() + labs(x = "summed RNA Replicate 1", y = "summed RNA Replicate 4") + geom_text( data = corr_sr14_0, aes(x= x,y = y,label= label), parse = TRUE) #+ facet_wrap(condition~library, scales = "free")
g.sumrna14_0
corr_sr13_0 <- data.frame(x = 1, y = 25, label =  paste("italic(r)==", round(cor(filter(sumInsane, sum_RNA1 > 0.1, sum_RNA3 > 0.1, condition == 0)$sum_RNA1, filter(sumInsane, sum_RNA1 > 0.1, sum_RNA3 > 0.1, condition == 0)$sum_RNA3, use = "complete.obs"), digits = 3)))
g.sumrna13_0 <- ggplot(filter(sumInsane, sum_RNA3 > 0.1, sum_RNA1 > 0.1, condition == 0), aes(sum_RNA1, sum_RNA3)) +geom_point(alpha=0.1) + scale_x_log10() + scale_y_log10() + labs(x = "summed RNA Replicate 1", y = "summed RNA Replicate 3") + geom_text( data = corr_sr13_0, aes(x= x,y = y,label= label), parse = TRUE) #+ facet_wrap(condition~library, scales = "free")
g.sumrna13_0
corr_sr34_0 <- data.frame(x = 1, y = 25, label =  paste("italic(r)==", round(cor(filter(sumInsane, sum_RNA4 > 0.1, sum_RNA3 > 0.1, condition == 0)$sum_RNA4, filter(sumInsane, sum_RNA4 > 0.1, sum_RNA3 > 0.1, condition == 0)$sum_RNA3, use = "complete.obs"), digits = 3)))
g.sumrna34_0 <- ggplot(filter(sumInsane, sum_RNA3 > 0.1, sum_RNA4 > 0.1, condition == 0), aes(sum_RNA3, sum_RNA4)) +geom_point(alpha=0.1) + scale_x_log10() + scale_y_log10() + labs(x = "summed RNA Replicate 3", y = "summed RNA Replicate 4") + geom_text( data = corr_sr34_0, aes(x= x,y = y,label= label), parse = TRUE) #+ facet_wrap(condition~library, scales = "free")
g.sumrna34_0

sam <- plot_grid(NULL, g.sumdna14_0, g.sumdna34_0, g.sumrna14_0, NULL, g.sumdna13_0, g.sumrna34_0, g.sumrna13_0, NULL, ncol= 3)
ggdraw(sam) + draw_label("R8000 summed BC Biological replicates at time 0h", x = 0.5, y =1, vjust = 1, hjust = 0.5, size =20)

#######Summed DNA at time six comparisons###############
corr_sd14_6 <- data.frame(x = 5, y = 25, label =  paste("italic(r)==", round(cor(filter(sumInsane, sum_DNA1 > 1, sum_DNA4 > 1, condition == 6)$sum_DNA1, filter(sumInsane, sum_DNA1 > 1, sum_DNA4 > 1, condition == 6)$sum_DNA4, use = "complete.obs"), digits = 3)))
g.sumdna14_6 <- ggplot(filter(sumInsane, sum_DNA4 > 1, sum_DNA1 >1, condition == 6), aes(sum_DNA1, sum_DNA4)) +geom_point(alpha=0.1) + scale_x_log10() + scale_y_log10() + labs(x = "summed DNA Replicate 1", y = "summed DNA Replicate 4") + geom_text( data = corr_sd14_6, aes(x= x,y = y,label= label), parse = TRUE) #+ facet_wrap(condition~library, scales = "free")
g.sumdna14_6
corr_sd13_6 <- data.frame(x = 5, y = 25, label =  paste("italic(r)==", round(cor(filter(sumInsane, sum_DNA1 > 1, sum_DNA3 > 1, condition == 6)$sum_DNA1, filter(sumInsane, sum_DNA1 > 1, sum_DNA3 > 1, condition == 6)$sum_DNA3, use = "complete.obs"), digits = 3)))
g.sumdna13_6 <- ggplot(filter(sumInsane, sum_DNA3 > 1, sum_DNA1 >1, condition == 6), aes(sum_DNA1, sum_DNA3)) +geom_point(alpha=0.1) + scale_x_log10() + scale_y_log10() + labs(x = "summed DNA Replicate 1", y = "summed DNA Replicate 3") + geom_text( data = corr_sd13_6, aes(x= x,y = y,label= label), parse = TRUE) #+ facet_wrap(condition~library, scales = "free")
g.sumdna13_6
corr_sd34_6 <- data.frame(x = 5, y = 25, label =  paste("italic(r)==", round(cor(filter(sumInsane, sum_DNA4 > 1, sum_DNA3 > 1, condition == 6)$sum_DNA4, filter(sumInsane, sum_DNA4 > 1, sum_DNA3 > 1, condition == 6)$sum_DNA3, use = "complete.obs"), digits = 3)))
g.sumdna34_6 <- ggplot(filter(sumInsane, sum_DNA3 > 1, sum_DNA4 >1, condition == 6), aes(sum_DNA3, sum_DNA4)) +geom_point(alpha=0.1) + scale_x_log10() + scale_y_log10() + labs(x = "summed DNA Replicate 3", y = "summed DNA Replicate 4") + geom_text( data = corr_sd34_6, aes(x= x,y = y,label= label), parse = TRUE) #+ facet_wrap(condition~library, scales = "free")
g.sumdna34_6

#######RNA at time zero comparisons###############
corr_sr14_6 <- data.frame(x = 2, y = 25, label =  paste("italic(r)==", round(cor(filter(sumInsane, sum_RNA1 > 0.1, sum_RNA4 > 0.1, condition == 6)$sum_RNA1, filter(sumInsane, sum_RNA1 > 0.1, sum_RNA4 > 0.1, condition == 6)$sum_RNA4, use = "complete.obs"), digits = 3)))
g.sumrna14_6 <- ggplot(filter(sumInsane, sum_RNA4 > 0.1, sum_RNA1 >0.1, condition == 6), aes(sum_RNA1, sum_RNA4)) +geom_point(alpha=0.1) + scale_x_log10() + scale_y_log10() + labs(x = "summed RNA Replicate 1", y = "summed RNA Replicate 4") + geom_text( data = corr_sr14_6, aes(x= x,y = y,label= label), parse = TRUE) #+ facet_wrap(condition~library, scales = "free")
g.sumrna14_6
corr_sr13_6 <- data.frame(x = 1, y = 25, label =  paste("italic(r)==", round(cor(filter(sumInsane, sum_RNA1 > 0.1, sum_RNA3 > 0.1, condition == 6)$sum_RNA1, filter(sumInsane, sum_RNA1 > 0.1, sum_RNA3 > 0.1, condition == 6)$sum_RNA3, use = "complete.obs"), digits = 3)))
g.sumrna13_6 <- ggplot(filter(sumInsane, sum_RNA3 > 0.1, sum_RNA1 > 0.1, condition == 6), aes(sum_RNA1, sum_RNA3)) +geom_point(alpha=0.1) + scale_x_log10() + scale_y_log10() + labs(x = "summed RNA Replicate 1", y = "summed RNA Replicate 3") + geom_text( data = corr_sr13_6, aes(x= x,y = y,label= label), parse = TRUE) #+ facet_wrap(condition~library, scales = "free")
g.sumrna13_6
corr_sr34_6 <- data.frame(x = 1, y = 25, label =  paste("italic(r)==", round(cor(filter(sumInsane, sum_RNA4 > 0.1, sum_RNA3 > 0.1, condition == 6)$sum_RNA4, filter(sumInsane, sum_RNA4 > 0.1, sum_RNA3 > 0.1, condition == 6)$sum_RNA3, use = "complete.obs"), digits = 3)))
g.sumrna34_6 <- ggplot(filter(sumInsane, sum_RNA3 > 0.1, sum_RNA4 > 0.1, condition == 6), aes(sum_RNA3, sum_RNA4)) +geom_point(alpha=0.1) + scale_x_log10() + scale_y_log10() + labs(x = "summed RNA Replicate 3", y = "summed RNA Replicate 4") + geom_text( data = corr_sr34_6, aes(x= x,y = y,label= label), parse = TRUE) #+ facet_wrap(condition~library, scales = "free")
g.sumrna34_6


pim <- plot_grid(NULL, g.sumdna14_6, g.sumdna34_6, g.sumrna14_6, NULL, g.sumdna13_6, g.sumrna34_6, g.sumrna13_6, NULL, ncol= 3)
ggdraw(pim) + draw_label("R8000 summed BC Biological replicates at time 6h", x = 0.5, y =1, vjust = 1, hjust = 0.5, size =20)


#######RNA at time zero comparisons###############
corr_srd14_0 <- data.frame(x = 2, y = 25, label =  paste("italic(r)==", round(cor(filter(sumInsane, sum_RD1 > 0.1, sum_RD4 > 0.1, condition == 0)$sum_RD1, filter(sumInsane, sum_RD1 > 0.1, sum_RD4 > 0.1, condition == 0)$sum_RD4, use = "complete.obs"), digits = 3)))
g.sumrnadna14_0 <- ggplot(filter(sumInsane, sum_RD4 > 0.1, sum_DNA1 > 0.1, condition == 0), aes(sum_RD1, sum_RD4)) +geom_point(alpha=0.1) + scale_x_log10() + scale_y_log10() + labs(x = "summed RNA/DNA Replicate 1", y = "summed RNA/DNA Replicate 4") + geom_text( data = corr_srd14_0, aes(x= x,y = y,label= label), parse = TRUE) #+ facet_wrap(condition~library, scales = "free")
g.sumrnadna14_0
corr_srd13_0 <- data.frame(x = 1, y = 25, label =  paste("italic(r)==", round(cor(filter(sumInsane, sum_RD1 > 0.1, sum_RD3 > 0.1, condition == 0)$sum_RD1, filter(sumInsane, sum_RD1 > 0.1, sum_RD3 > 0.1, condition == 0)$sum_RD3, use = "complete.obs"), digits = 3)))
g.sumrnadna13_0 <- ggplot(filter(sumInsane, sum_RD3 > 0.1, sum_RD1 > 0.1, condition == 0), aes(sum_RD1, sum_RD3)) +geom_point(alpha=0.1) + scale_x_log10() + scale_y_log10() + labs(x = "summed RNA/DNA Replicate 1", y = "summed RNA/DNA Replicate 3") + geom_text( data = corr_srd13_0, aes(x= x,y = y,label= label), parse = TRUE) #+ facet_wrap(condition~library, scales = "free")
g.sumrnadna13_0
corr_srd34_0 <- data.frame(x = 1, y = 25, label =  paste("italic(r)==", round(cor(filter(sumInsane, sum_RD4 > 0.1, sum_RD3 > 0.1, condition == 0)$sum_RD4, filter(sumInsane, sum_RD4 > 0.1, sum_RD3 > 0.1, condition == 0)$sum_RD3, use = "complete.obs"), digits = 3)))
g.sumrnadna34_0 <- ggplot(filter(sumInsane, sum_RD3 > 0.1, sum_RD4 > 0.1, condition == 0), aes(sum_RD3, sum_RD4)) +geom_point(alpha=0.1) + scale_x_log10() + scale_y_log10() + labs(x = "summed RNA/DNA Replicate 3", y = "summed RNA/DNA Replicate 4") + geom_text( data = corr_srd34_0, aes(x= x,y = y,label= label), parse = TRUE) #+ facet_wrap(condition~library, scales = "free")
g.sumrnadna34_0
#######RNA at time zero comparisons###############
corr_srd14_6 <- data.frame(x = .4, y = 10, label =  paste("italic(r)==", round(cor(filter(sumInsane, sum_RD1 > 0.1, sum_RD4 > 0.1, condition == 6)$sum_RD1, filter(sumInsane, sum_RD1 > 0.1, sum_RD4 > 0.1, condition == 6)$sum_RD4, use = "complete.obs"), digits = 3)))
g.sumrnadna14_6 <- ggplot(filter(sumInsane, sum_RD4 > 0.1, sum_DNA1 > 0.1, condition == 6), aes(sum_RD1, sum_RD4)) +geom_point(alpha=0.1, fill = mean_DNA4) + scale_x_log10() + scale_y_log10() + labs(x = "summed RNA/DNA Replicate 1", y = "summed RNA/DNA Replicate 4", title = "171218 R8000 summed RNA/DNA counts between replicates") + geom_text( data = corr_srd14_6, aes(x= x,y = y,label= label), parse = TRUE) #+ facet_wrap(condition~library, scales = "free")
g.sumrnadna14_6
corr_srd13_6 <- data.frame(x = 1, y = 25, label =  paste("italic(r)==", round(cor(filter(sumInsane, sum_RD1 > 0.1, sum_RD3 > 0.1, condition == 6)$sum_RD1, filter(sumInsane, sum_RD1 > 0.1, sum_RD3 > 0.1, condition == 6)$sum_RD3, use = "complete.obs"), digits = 3)))
g.sumrnadna13_6 <- ggplot(filter(sumInsane, sum_RD3 > 0.1, sum_RD1 > 0.1, condition == 6), aes(sum_RD1, sum_RD3, color = filter(sumInsane, sum_RD3 > 0.1, sum_RD1 > 0.1, condition == 6)$mean_DNA4)) +geom_point(alpha=0.1) + scale_x_log10() + scale_y_log10() + labs(x = "summed RNA/DNA Replicate 1", y = "summed RNA/DNA Replicate 3") + geom_text( data = corr_srd13_6, aes(x= x,y = y,label= label), parse = TRUE) #+ facet_wrap(condition~library, scales = "free")
g.sumrnadna13_6
corr_srd34_6 <- data.frame(x = 1, y = 25, label =  paste("italic(r)==", round(cor(filter(sumInsane, sum_RD4 > 0.1, sum_RD3 > 0.1, condition == 6)$sum_RD4, filter(sumInsane, sum_RD4 > 0.1, sum_RD3 > 0.1, condition == 6)$sum_RD3, use = "complete.obs"), digits = 3)))
g.sumrnadna34_6 <- ggplot(filter(sumInsane, sum_RD3 > 0.1, sum_RD4 > 0.1, condition == 6), aes(sum_RD3, sum_RD4)) +geom_point(alpha=0.1) + scale_x_log10() + scale_y_log10() + labs(x = "summed RNA/DNA Replicate 3", y = "summed RNA/DNA Replicate 4") + geom_text( data = corr_srd34_6, aes(x= x,y = y,label= label), parse = TRUE) #+ facet_wrap(condition~library, scales = "free")
g.sumrnadna34_6

jim <- plot_grid(g.sumrnadna14_0, g.sumrnadna13_0, g.sumrnadna34_0, g.sumrnadna14_6, g.sumrnadna13_6, g.sumrnadna34_6, ncol= 3)
ggdraw(jim) + draw_label("R8000 summed BC Biological replicates at time 0h", x = 0.5, y =1, vjust = 1, hjust = 0.5, size =20) + draw_label("R8000 summed BC Biological replicates at time 6h", x = 0.5, y =0.5, vjust = 1, hjust = 0.5, size =20)

################################ GFP and normalization####################################
load(file = "GFPBCs.Rda")
colnames(GFPBCs) <- c("plasmid", "short", "fill", "barcode")
GFPBCs <- select(GFPBCs, plasmid, barcode)
GFPBCs <- separate(GFPBCs, plasmid, into = c("plasmid","number"), sep ="-")

R8000_GFP <- inner_join(GFPBCs, R8000_full_DF, by = "barcode")
R8000_GFP_R <- filter(R8000_GFP, grepl("RNA", condition)) %>%
  separate(condition, into = c("template","time"), sep = "-") %>%
  mutate(rCount = count/readsinlib *1000000) %>%
  select(-collapsed_bc)
R8000_GFP_D <- filter(R8000_GFP, grepl("DNA", condition)) %>%
  separate(condition, into = c("template","time"), sep = "-") %>%
  mutate(dCount = count/readsinlib *1000000) %>%
  select(-collapsed_bc)
R8000_GFP_RD <- left_join(R8000_GFP_D, R8000_GFP_R, by = c("plasmid","number","barcode","time")) %>%
  mutate(RNADNA = rCount/dCount*100) %>%
  select(-readsinlib.x, -readsinlib.y, -template.x, -template.y)

  
R8000_GFP_RD <- filter(R8000_GFP_RD, !(barcode == "GGAAACAGAATTCTGGGACC"), !(barcode =="AGAAACAGAATTCTGGGACC"))

sumGFP <- group_by(R8000_GFP_RD, plasmid, time) %>%
  summarise(medD = median(dCount), medR = median(rCount), medRD = median(RNADNA), meanRD = mean(RNADNA), sdRD = sd(RNADNA))

plasmid_labels <-  c("0" = "0h RNA/DNA", "6" = "6h RNA/DNA")
gfperr <- aes(ymax = meanRD + sdRD, ymin = meanRD - sdRD)
g.GFP <- ggplot(sumGFP, aes(plasmid,medRD, fill=plasmid)) + geom_bar(stat = "identity") + geom_errorbar(gfperr) + facet_wrap(~time, labeller = labeller(time = plasmid_labels)) + 
  labs(title = "171113 R8000 GFP constructs post-normalization", x = "Plasmid", y = "Median RNA/DNA") + geom_point(data = R8000_GFP_RD, aes(plasmid, RNADNA), alpha = 0.3, show.legend = FALSE) +
  theme(strip.background=element_rect(fill = "White")) + geom_text(data = R8000_GFP_RD, aes(plasmid, RNADNA, label = paste(plasmid,"-",number, sep = "")), nudge_x = 0.25, size = 2.5)
g.GFP

norm = sum(filter(sumGFP, time == 0)$medRD)/sum(filter(sumGFP, time == 6)$medRD)
