################################################# Standard libraries + viridis to make things pretty + extrafont for comic sans #################################################
library(dplyr)
library(tidyr)
library(extrafont)
library(viridis)
library(cowplot)

################################################# Read in the data and mash it together #################################################
setwd("/Users/Cliff/Documents/Kosuri Lab/NGS files/170926 R8000 barcode-1/")

R8000_RNA_0 <- read.table("sc-counts_R8000_0h_RNA.txt", header = FALSE)
R8000_RNA_0 <- data.frame(R8000_RNA_0, rep("RNA-0",nrow(R8000_RNA_0)),rep(sum(R8000_RNA_0$V2), nrow(R8000_RNA_0)))
colnames(R8000_RNA_0) <- c("barcode", "count", "collapsed_bc", "condition","readsinlib")

R8000_DNA_0 <- read.table("sc-counts_R8000_0h_DNA.txt", header = FALSE)
R8000_DNA_0 <- data.frame(R8000_DNA_0, rep("DNA-0",nrow(R8000_DNA_0)),rep(sum(R8000_DNA_0$V2), nrow(R8000_DNA_0)))
colnames(R8000_DNA_0) <- c("barcode", "count", "collapsed_bc", "condition","readsinlib")

R8000_RNA_6 <- read.table("sc-counts_R8000_6h_RNA.txt", header = FALSE)
R8000_RNA_6 <- data.frame(R8000_RNA_6, rep("RNA-6",nrow(R8000_RNA_6)),rep(sum(R8000_RNA_6$V2), nrow(R8000_RNA_6)))
colnames(R8000_RNA_6) <- c("barcode", "count", "collapsed_bc", "condition","readsinlib")

R8000_DNA_6 <-read.table("sc-counts_R8000_6h_DNA.txt", header = FALSE)
R8000_DNA_6 <- data.frame(R8000_DNA_6, rep("DNA-6",nrow(R8000_DNA_6)),rep(sum(R8000_DNA_6$V2), nrow(R8000_DNA_6)))
colnames(R8000_DNA_6) <- c("barcode", "count", "collapsed_bc", "condition","readsinlib")

R8000_full_DF <- bind_rows(R8000_RNA_0, R8000_DNA_0, R8000_RNA_6, R8000_DNA_6)

g.raw <- ggplot(R8000_full_DF, aes(condition, count, fill = condition)) + geom_violin()  + labs(title = "171004 R8000-1+2 counts distribution", y = "Number of reads per barcode", x = "Condition") + scale_fill_brewer(palette = "Set1") + 
  geom_boxplot(width = 0.025, fill = "Black") + stat_summary(fun.y = mean, geom = "point", color = "White", show.legend = FALSE)
g.raw


R8000_join <- full_join(R8000_RNA_0, R8000_RNA_6, by = "barcode")
R8000_join <- select(R8000_join, barcode, "RNA_0" = count.x, "RNA_6" = count.y)
R8000_join <- full_join(R8000_join, R8000_DNA_0, by = "barcode")
R8000_join <- full_join(R8000_join, R8000_DNA_6, by = "barcode")
R8000_join <- select(R8000_join, barcode, RNA_0, RNA_6, "DNA_0" = count.x, "DNA_6" = count.y)
R8000_join$RNA_0[is.na(R8000_join$RNA_0)] <-0
R8000_join$RNA_6[is.na(R8000_join$RNA_6)] <- 0
R8000_join$DNA_0[is.na(R8000_join$DNA_0)] <-0
R8000_join$DNA_6[is.na(R8000_join$DNA_6)] <- 0

g.R1vR6 <- ggplot(R8000_join, aes(RNA_0, RNA_6)) + geom_point(alpha = 0.1) + scale_x_log10() + scale_y_log10() + geom_density2d() + labs(title = "171004 R8000 RNA counts at 0h vs RNA counts at 6h", x = "RNA barcode counts at 0h", y = "RNA  barcode counts at 6h")
g.R1vR6

g.R1vD1 <- ggplot(R8000_join, aes(DNA_0, RNA_0)) + geom_point(alpha = 0.1) + scale_y_log10() + scale_x_log10() + labs(title = "171004 R8000 RNA counts at 0h vs DNA counts at 0h", x = "DNA barcode counts at 0h", y = "RNA  barcode counts at 0h")
g.R1vD1

g.R6vD6 <- ggplot(R8000_join, aes(DNA_6, RNA_6)) + geom_point(alpha = 0.1) + scale_y_log10() + scale_x_log10() + geom_density2d() + labs(title = "171004 R8000 RNA counts at 6h vs DNA counts at 6h", x = "DNA barcode counts at 6h", y = "RNA  barcode counts at 6h")
g.R6vD6

g.D1vD6 <- ggplot(R8000_join, aes(DNA_0, DNA_6)) + geom_point(alpha = 0.1) + scale_x_log10() + scale_y_log10() + labs(title = "171004 R8000 DNA counts at 0h vs DNA counts at 0h", x = "DNA barcode counts at 0h", y = "DNA  barcode counts at 6h")
g.D1vD6

R8000_full_DF_3 <- full_join(R8000_full_DF,R8000_full_DF_2, by = c("barcode", "condition"))
R8000_full_DF_3$count.x[is.na(R8000_full_DF_3$count.x)] <- 0
R8000_full_DF_3$count.y[is.na(R8000_full_DF_3$count.y)] <- 0
R8000_full_DF_3$readsinlib.y <- na.locf(R8000_full_DF_3$readsinlib.y)
R8000_full_DF_3$readsinlib.x <- na.locf(R8000_full_DF_3$readsinlib.x)
R8000_full_DF_3 <- mutate(R8000_full_DF_3, count = count.x + count.y, readsinlib = readsinlib.x + readsinlib.y)


R8000_map <- read.table("R8000-c1.x-y.txt", skip = 1)
colnames(R8000_map) <- c("barcode", "X_peptide", "Y_peptide", "reads")
R8000_map <- filter(R8000_map, X_peptide != "<NA>", Y_peptide != "<NA>") %>%
  filter(as.character(X_peptide) == as.character(Y_peptide)) %>%
  select(-Y_peptide) %>%
  separate(X_peptide, into = c("X_peptide", "X_codon"), sep = "_") %>%
  separate(X_peptide, into = c("X_peptide","Y_peptide"), sep = "\\+") %>%
  separate(X_peptide, into = c("X_peptide", "X_group"), sep = "-") %>%
  separate(Y_peptide, into = c("Y_peptide", "Y_group"), sep = "-")


R8000_full_DF_3 <- mutate(R8000_full_DF_3, "nCount" = count/readsinlib*1000000)
R8000_full_DF_3 <- inner_join(R8000_full_DF_3, R8000_map, by = c("barcode"))

all_DF <- full_join(R8000_full_DF,R8000_full_DF_2, by = c("barcode","X_peptide","Y_peptide", "X_group", "Y_group", "condition"))

g.all <- ggplot(all_DF, aes(nCount.x, nCount.y, color = condition)) + geom_point(alpha = 0.1) + labs (title = "170927 R8000 Biological replicates normalized counts", x = "Replicate 1 (RPM)", y = "Replicate 2 (RPM)") +
  geom_abline(slope = 1) + scale_x_continuous(limits = c(0,115)) + scale_y_continuous(limits = c(0,115)) + geom_text(data = corrlab, aes(x=x, y=y, label= label), color = "Black", parse = TRUE)
g.all

all_DF$nCount.x[is.na(all_DF$nCount.x)] <- 0
all_DF$nCount.y[is.na(all_DF$nCount.y)] <- 0
corr <- cor(all_DF$nCount.x, all_DF$nCount.y)
corrlab <- data.frame( x = 30, y = 90, label = paste("italic(r)==",round(corr, digits = 3)))

R8000RNA <- filter(R8000_full_DF_3, grepl("RNA", R8000_full_DF_3$condition)) %>%
  select(X_peptide, X_group, X_codon, Y_peptide, condition, barcode, "rcount" = nCount ) %>%
  separate(condition, into = c("template", "condition"))
sumRNA <- group_by(R100RNA, X_peptide, Y_peptide, condition) %>%
  summarise(mean_RNA = mean(rcount), sd_R = sd(rcount), bcnum = n_distinct(barcode))

R8000DNA <- filter(R8000_full_DF_3, grepl("DNA", R8000_full_DF_3$condition)) %>%
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
  summarise(mean_RNA = mean(rcount), sd_R = sd(rcount), mean_DNA = mean(dcount), sd_D = sd(dcount), mean_RD = mean(RNADNA), sd_RD = sd(RNADNA), med_RD = median(RNADNA), sum_DNA = sum(dcount), sum_RNA = sum(rcount), bcnum = n_distinct(barcode), SEM = sd(RNADNA)/sqrt(n()))

sumR8000C1 <- filter(sumR100RD, X_codon == "C1")
sumR8000C2 <- filter(sumR100RD, X_codon == "C2")
sumR8000C3 <- filter(sumR100RD, X_codon == "C3")

sumR8full <- full_join(sumR8000C1,sumR8000C2, by = c("X_peptide", "Y_peptide","condition","X_group"))
sumR8full <- full_join(sumR8full, sumR8000C3, by = c("X_peptide","Y_peptide","condition","X_group"))

sumR8full$mean_RD.x[is.na(sumR8full$mean_RD.x)] <- 0
sumR8full$mean_RD.y[is.na(sumR8full$mean_RD.y)] <- 0
sumR8full$mean_RD[is.na(sumR8full$mean_RD)] <- 0
g.RDC1C2 <- ggplot(filter(sumR8full, bcnum.x >5, bcnum.y >5), aes(mean_DNA.x, mean_DNA.y, color = log(bcnum.x + bcnum.y))) + geom_point(alpha = 0.2) +geom_text(data = corr1, aes(x=x, y=y, label= label), color = "Black", parse = TRUE) +
  labs(title = "170927 R8000 codon usage 1 vs codon usage 2 mean DNA (bc >5) ", x = "Mean DNA (Codon 1)", y = "Mean DNA (Codon 2)") +  scale_x_log10(limits = c(0.1,30)) + scale_y_log10(limits = c(0.1,30))  +scale_color_viridis(name = "Total number \nof barcodes" )
g.RDC1C2

corr1 <- cor(all_DF$nCount.x, all_DF$nCount.y)
corr1 <- data.frame( x = 30, y = 300, label = paste("italic(r)==",round(cor(filter(sumR8full, !is.na(mean_RD.x), !is.na(mean_RD.y))$mean_RD.x, filter(sumR8full, !is.na(mean_RD.x), !is.na(mean_RD.y))$mean_RD.y), digits = 3)))

g.RDC1C3 <- ggplot(filter(sumR8full, bcnum.x > 5, bcnum > 5), aes(mean_RD.x, mean_RD, color = condition)) + geom_point(alpha = 0.2) +geom_text(data = corr2, aes(x=x, y=y, label= label), color = "Black", parse = TRUE) +
  labs(title = "170927 R8000 codon usage 1 vs codon usage 3 mean RNA/DNA (bc>5)", x = "Mean RNA/DNA (Codon 1)", y = "Mean RNA/DNA (Codon 3)") + scale_x_log10() + scale_y_log10()
g.RDC1C3

corr2 <- data.frame( x = 30, y = 90, label = paste("italic(r)==",round(cor(filter(sumR8full, bcnum.x >5, bcnum >5)$mean_RD.x, filter(sumR8full, bcnum.x > 5, bcnum > 5)$mean_RD), digits = 3)))

g.RDC2C3 <- ggplot(filter(sumR8full, bcnum.y > 5, bcnum > 5), aes(mean_RD.y, mean_RD, color = condition)) + geom_point(alpha = 0.2) +geom_text(data = corr3, aes(x=x, y=y, label= label), color = "Black", parse = TRUE) +
  labs(title = "170927 R8000 codon usage 2 vs codon usage 3 mean RNA/DNA (bc>5)", x = "Mean RNA/DNA (Codon 2)", y = "Mean RNA/DNA (Codon 3)") + scale_x_log10() + scale_y_log10()
g.RDC2C3

corr3 <- data.frame( x = 30, y = 90, label = paste("italic(r)==",round(cor(filter(sumR8full, bcnum.y >5, bcnum >5)$mean_RD.y, filter(sumR8full, bcnum.y > 5, bcnum > 5)$mean_RD), digits = 3)))

sumR10comp <- filter(sumR100RD, grepl("O", X_peptide), is.na(X_group), condition == 6)
sumR10comp <- left_join(sumR10comp, mason2, by = c("X_peptide","Y_peptide"))

colrs <- c("c1" = "#FF0000","c2"="#000000")
g.mason <- ggplot(sumR10comp, aes(mean_RD, Tm)) + geom_point(aes(color="c2")) + labs(title = "171003 R8000 O and Drobnak constructs mean RNA/DNA vs reported Tm", x = "Mean RNA/DNA", y = "Reported Tm") + geom_smooth(formula = y ~ log(x), span = 9) + 
  geom_point(data = sumR10drob, aes(mean_RD, Tm, color = "c1")) + scale_color_manual(name = "Data set", breaks = c("c1","c2"), values = colrs, labels = c("Drobnak","Mason")) + scale_x_log10(limits = c(40, 1001))
g.mason

sumR10drob <- filter(right_join(sumR100RD, drobnak, by = c("X_peptide", "Y_peptide")), !is.na(X_peptide), condition == 6)

g.drob <- ggplot(sumR10drob, aes(mean_RD, Tm)) + geom_point()
g.drob

sumR100RD$X_group[is.na(sumR100RD$X_group) & grepl("mSN$",sumR100RD$X_peptide)] <- "mSN"
sumR100RD$X_group[is.na(sumR100RD$X_group) & grepl("Hb$",sumR100RD$X_peptide)] <- "Hb"
sumR100RD$X_group[is.na(sumR100RD$X_group) & grepl("Hq$",sumR100RD$X_peptide)] <- "Hq"
sumR100RD$X_group[is.na(sumR100RD$X_group) & grepl("bH1$",sumR100RD$X_peptide)] <- "bH1"
sumR100RD$X_group[is.na(sumR100RD$X_group) & grepl("mS$",sumR100RD$X_peptide)] <- "mS"
sumR100RD$X_group[is.na(sumR100RD$X_group) & grepl("A$",sumR100RD$X_peptide)] <- "A"
sumR100RD$X_group[is.na(sumR100RD$X_group) & grepl("SH$",sumR100RD$X_peptide)] <- "SH"

induction_labels <- c("0" = "0 hours induction RNA/DNA", "6" = "6 hours induction RNA/DNA")

os <-filter(sumR100RD, grepl("O",X_peptide))

g.sumR100RD <- ggplot(filter(sumR100RD, grepl("O",X_peptide), bcnum > 5), aes(X_peptide, Y_peptide)) + geom_tile(aes(fill=log(med_RD))) + facet_grid(X_group~condition, labeller = labeller(condition=induction_labels), scales = "free") +
  theme(text=element_text(family="Calibri"), strip.text=element_text(size = 10), strip.background=element_rect(fill = "White"), axis.text.x = element_text(angle=90,size =5),
        axis.text.y = element_text(size = 5), legend.title=element_text(size = 8), legend.title.align=0.5, legend.text=element_text(size =8)) +
  scale_fill_viridis(option="viridis", name = "Med \nRNA to DNA  ", na.value="white") + labs(title="Roman's 8000 RNA/DNA median barcode ratio O constructs (bc >5)", x="X peptide", y = "Y peptide")
g.sumR100RD

g.bccount <- ggplot(R100_full_DF, aes(count)) + geom_histogram() + labs(x = "Number of reads", y = "number of barcodes", title = "170525 R100-2-2 Hiseq Barcode counts") + scale_x_continuous(limits = c(0, 1000))
g.bccount

R100_raw_map <- full_join(R100_raw_map, R100_DNA_0, by ="barcode")
R100_raw_map$count.x[is.na(R100_raw_map$count.x)] <- 1
R100_raw_map$count.y[is.na(R100_raw_map$count.y)] <- 1

g.mapcount <- ggplot(R100_raw_map, aes(count.y, count.x)) + geom_point(alpha = 0.1) + scale_x_log10() + scale_y_log10() + labs(x = "Barcode-seq DNA barcode counts", y = "Mapping run barcode counts", title = "170526 comparison of barcode counts >=0")
g.mapcount





################################################### To compare against previously published data ###############################################
mason <- read.csv("mason-measured-Tm-both.csv", header = TRUE)
mason_names <- c(paste("O",1:16, sep = ""))
colnames(mason) <- mason_names
rownames(mason) <- mason_names
mason <- mason[2:17]
mason <- data.frame(rownames(mason),mason)
mason2 <- gather(mason, rownames.mason.)
colnames(mason2) <- c("X_peptide","Y_peptide","Tm") 

drobnak <- read.csv("Drobnak-JACS-Tm.csv", header = TRUE)
colnames(drobnak2) <- c("X_peptide","Y_peptide","Tm","bcipa_Tm")
drobnak2 <- select(drobnak, Y_peptide,X_peptide,Tm,bcipa_Tm)
drobnak <- rbind(drobnak,drobnak2)
