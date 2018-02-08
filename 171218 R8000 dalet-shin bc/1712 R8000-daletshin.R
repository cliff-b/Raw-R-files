################################################# Standard libraries + viridis to make things pretty + extrafont for comic sans #################################################
library(dplyr)
library(tidyr)
library(extrafont)
library(viridis)
library(cowplot)
library(nnet)
library(caret)
################################################# Read in the data and mash it together #################################################
setwd("/Users/Cliff/Google Drive/NGS files/171218 R8000 dalet-shin bc/")

R8000_RNA_0 <- read.table("sc-counts-R8000-9h-s-RNA.txt", header = FALSE)
R8000_RNA_0 <- data.frame(R8000_RNA_0, rep("RNA-9",nrow(R8000_RNA_0)),rep(sum(R8000_RNA_0$V2), nrow(R8000_RNA_0)))
colnames(R8000_RNA_0) <- c("barcode", "count", "collapsed_bc", "condition","readsinlib")
R8000_RNA_0 <- mutate(R8000_RNA_0, count =  count/norm)

R8000_DNA_0 <- read.table("sc-counts-R8000-9h-s-DNA.txt", header = FALSE)
R8000_DNA_0 <- data.frame(R8000_DNA_0, rep("DNA-9",nrow(R8000_DNA_0)),rep(sum(R8000_DNA_0$V2), nrow(R8000_DNA_0)))
colnames(R8000_DNA_0) <- c("barcode", "count", "collapsed_bc", "condition","readsinlib")

R8000_RNA_6 <- read.table("sc-counts-R8000-38h-d-RNA.txt", header = FALSE)
R8000_RNA_6 <- data.frame(R8000_RNA_6, rep("RNA-38",nrow(R8000_RNA_6)),rep(sum(R8000_RNA_6$V2), nrow(R8000_RNA_6)))
colnames(R8000_RNA_6) <- c("barcode", "count", "collapsed_bc", "condition","readsinlib")

R8000_DNA_6 <-read.table("sc-counts-R8000-38h-d-DNA.txt", header = FALSE)
R8000_DNA_6 <- data.frame(R8000_DNA_6, rep("DNA-38",nrow(R8000_DNA_6)),rep(sum(R8000_DNA_6$V2), nrow(R8000_DNA_6)))
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

g.bccount <- ggplot(sumR100RD, aes(bcnum, fill = mean_RD)) + geom_histogram(bins = 114) + scale_x_continuous(limits = c(0,113)) + labs(title = "171010 R8000 Number of barcodes per construct", x = "number of barcodes")  +
  scale_fill_viridis()
g.bccount

g.sumOs <- ggplot(filter(sumR100RD, grepl("O",X_peptide)), aes(X_peptide, Y_peptide)) + geom_tile(aes(fill=med_RD)) + facet_grid(X_group~condition, labeller = labeller(condition=induction_labels), scales = "free") +
  theme(text=element_text(family="Calibri"), strip.text=element_text(size = 10), strip.background=element_rect(fill = "White"), axis.text.x = element_text(angle=90,size =5),
        axis.text.y = element_text(size = 5), legend.title=element_text(size = 8), legend.title.align=0.5, legend.text=element_text(size =8)) +
  scale_fill_viridis(option="viridis", name = "Mean \nRNA to DNA  ", na.value="white") + labs(title="171006 Roman's 8000 RNA/DNA mean barcode ratio O constructs (bc >0)", x="X peptide", y = "Y peptide")
g.sumOs

g.sumPs <- ggplot(filter(sumR100RD, grepl("P",X_peptide), !is.na(X_group)), aes(X_peptide, Y_peptide)) + geom_tile(aes(fill=mean_RD)) + facet_grid(X_group~condition, labeller = labeller(condition=induction_labels)) +
  theme(text=element_text(family="Calibri"), strip.text=element_text(size = 10), strip.background=element_rect(fill = "White"), axis.text.x = element_text(angle=90,size =5),
        axis.text.y = element_text(size = 5), legend.title=element_text(size = 8), legend.title.align=0.5, legend.text=element_text(size =8)) +
  scale_fill_viridis(option="viridis", name = "Mean \nRNA to DNA  ", na.value="white") + labs(title="171006 Roman's 8000 RNA/DNA mean barcode ratio P constructs (bc >0)", x="X peptide", y = "Y peptide")
g.sumPs

g.sumPsx <- ggplot(filter(sumR100RD, grepl("P",X_peptide), X_group %in% grouper), aes(X_peptide, Y_peptide)) + geom_tile(aes(fill=sum_RD)) + facet_wrap(X_group~condition, labeller = labeller(condition=induction_labels), scales ="free") +
  theme(text=element_text(family="Calibri"), strip.text=element_text(size = 10), strip.background=element_rect(fill = "White"), axis.text.x = element_text(angle=90,size =5),
        axis.text.y = element_text(size = 5), legend.title=element_text(size = 8), legend.title.align=0.5, legend.text=element_text(size =8)) +
  scale_fill_viridis(option="viridis", name = "Mean \nRNA to DNA  ", na.value="white") + labs(title="171011 Roman's 8000 RNA/DNA summed RNA/DNA barcode ratio P constructs (bc >0)", x="X peptide", y = "Y peptide")
g.sumPsx



g.sum4Hs <- ggplot(filter(sumR100RD, grepl("4H",X_peptide), X_peptide %in% group2, Y_peptide %in% group2, bcnum >5, mean_DNA > 0.8), aes(X_peptide, Y_peptide)) + geom_tile(aes(fill=sum_RD)) + facet_wrap(condition~X_group, labeller = labeller(condition=induction_labels), scales = "free") +
  theme(text=element_text(family="Calibri"), strip.text=element_text(size = 10), strip.background=element_rect(fill = "White"), axis.text.x = element_text(angle=90,size =7),
        axis.text.y = element_text(size = 7), legend.title=element_text(size = 8), legend.title.align=0.5, legend.text=element_text(size =8)) +
  scale_fill_viridis(option="viridis", name = "Summed \nRNA to DNA  ", na.value="white") + labs(title="171011 Roman's 8000 RNA/DNA summed RNA/DNA barcode ratio 4H constructs (bc >0)", x="X peptide", y = "Y peptide")
g.sum4Hs

g.sumsvsmeans <- ggplot(sumR100RD, aes(mean_RD, sum_RD)) + geom_point(alpha = 0.1) + geom_text(data = coor5, aes(x, y, label = label), parse = TRUE) + scale_x_log10()+ scale_y_log10()
g.sumsvsmeans

coor5 <- data.frame(x = 300, y = 6, label =  paste("italic(r)==",round(cor(sumR100RD$mean_RD, sumR100RD$sum_RD), 3)))









sumR8000 <- select(sumR100RD, X_peptide, Y_peptide, condition, X_group, mean_RD, sd_RD, med_RD, sum_RD, bcnum)

sumR8000C1 <- filter(sumR100RD, X_codon == "C1")
sumR8000C2 <- filter(sumR100RD, X_codon == "C2")
sumR8000C3 <- filter(sumR100RD, X_codon == "C3")

sumR8full <- full_join(sumR8000C1,sumR8000C2, by = c("X_peptide", "Y_peptide","condition","X_group"))
sumR8full <- full_join(sumR8full, sumR8000C3, by = c("X_peptide","Y_peptide","condition","X_group"))

sumR8full$mean_RD.x[is.na(sumR8full$mean_RD.x)] <- 0
sumR8full$mean_RD.y[is.na(sumR8full$mean_RD.y)] <- 0
sumR8full$mean_RD[is.na(sumR8full$mean_RD)] <- 0

g.RDC1C2 <- ggplot(filter(sumR8full, bcnum.x >5, bcnum.y > 5), aes(sum_RD.x, sum_RD.y, color = log(bcnum.x + bcnum.y))) + geom_point(alpha = 0.2) +geom_text(data = corr1, aes(x=x, y=y, label= label), color = "Black", parse = TRUE) +
  labs(title = "171010 R8000 codon usage 1 vs codon usage 2 sum RNA/DNA (bc >5) ", x = "sum RNA/DNA (Codon 1)", y = "sum RNA/DNA (Codon 2)") +  scale_x_log10(limits = c(0.02,20)) + scale_y_log10(limits = c(0.02,20))  +scale_color_viridis(name = "Total number \nof barcodes \nlog2" )
g.RDC1C2

corr1 <- cor(all_DF$nCount.x, all_DF$nCount.y)
corr1 <- data.frame( x = .3, y = .3, label = paste("italic(r)==",round(cor(filter(sumR8full, bcnum.x >5 , bcnum.y > 5)$sum_RD.x, filter(sumR8full, bcnum.x>5, bcnum.y>5)$sum_RD.y), digits = 3)))

g.RDC1C3 <- ggplot(filter(sumR8full, bcnum.x > 5, bcnum > 5), aes(mean_RD.x, mean_RD, color = log10(bcnum.x + bcnum))) + geom_point(alpha = 0.2) +geom_text(data = corr2, aes(x=x, y=y, label= label), color = "Black", parse = TRUE) +
  labs(title = "171010 R8000 codon usage 1 vs codon usage 3 mean RNA/DNA (bc>5)", x = "Mean RNA/DNA (Codon 1)", y = "Mean RNA/DNA (Codon 3)") + scale_x_log10() + scale_y_log10() + scale_color_viridis(name = "Total number \nof barcodes \nlog10")
g.RDC1C3

corr2 <- data.frame( x = 30, y = 90, label = paste("italic(r)==",round(cor(filter(sumR8full, bcnum.x >5, bcnum >5)$mean_RD.x, filter(sumR8full, bcnum.x > 5, bcnum > 5)$mean_RD), digits = 3)))

g.RDC2C3 <- ggplot(filter(sumR8full, bcnum.y > 5, bcnum > 5), aes(mean_RD.y, mean_RD, color = log10(bcnum.y + bcnum))) + geom_point(alpha = 0.2) +geom_text(data = corr3, aes(x=x, y=y, label= label), color = "Black", parse = TRUE) +
  labs(title = "171010 R8000 codon usage 2 vs codon usage 3 mean RNA/DNA (bc>5)", x = "Mean RNA/DNA (Codon 2)", y = "Mean RNA/DNA (Codon 3)") + scale_x_log10() + scale_y_log10() + scale_color_viridis(name = "Total number \nof barcodes \nlog10")
g.RDC2C3

corr3 <- data.frame( x = 30, y = 90, label = paste("italic(r)==",round(cor(filter(sumR8full, bcnum.y >5, bcnum >5)$mean_RD.y, filter(sumR8full, bcnum.y > 5, bcnum > 5)$mean_RD), digits = 3)))

sumR10comp <- filter(sumR100RD, grepl("O", X_peptide), is.na(X_group), condition == 6)
sumR10comp <- left_join(sumR10comp, mason2, by = c("X_peptide","Y_peptide"))

colrs <- c("c1" = "#FF0000","c2"="#000000")
g.mason <- ggplot(sumR10comp, aes(mean_RD, Tm)) + geom_point(aes(color="c2")) + labs(title = "171009 R8000 O and Drobnak constructs mean RNA/DNA vs reported Tm", x = "Mean RNA/DNA", y = "Reported Tm") + geom_smooth(formula = y ~ log(x), span = 9) + 
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

grouper <- c("mSN","Hb","Hq","bH1","mS","A","SH")
group2 <- c("4H741","4H741","4H2143","4H2143","4H3277","4H3277","4H40","4H2915","4H121","4H2134","4H182","4H3057","4H319","4H2839","4H348","4H2588","4H477","4H2717","4H492","4H2692","4H2223","4H2503", "4H2438","4H2473")

induction_labels <- c("9" = "25C degree induction RNA/DNA", "38" = "17C degree induction RNA/DNA")

os <-filter(sumR100RD, grepl("O",X_peptide))

g.sumR100RD <- ggplot(filter(sumR100RD, grepl("O",X_peptide)), aes(X_peptide, Y_peptide)) + geom_tile(aes(fill=med_RD)) + facet_grid(X_group~condition, labeller = labeller(condition=induction_labels), scales = "free") +
  theme(text=element_text(family="Calibri"), strip.text=element_text(size = 10), strip.background=element_rect(fill = "White"), axis.text.x = element_text(angle=90,size =5),
        axis.text.y = element_text(size = 5), legend.title=element_text(size = 8), legend.title.align=0.5, legend.text=element_text(size =8)) +
  scale_fill_viridis(option="viridis", name = "Med \nRNA to DNA  ", na.value="white") + labs(title="171006 Roman's 8000 RNA/DNA median barcode ratio P constructs (bc >0)", x="X peptide", y = "Y peptide")
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


################################################Get GFPBCs ready to go#########################################################
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
  labs(title = "171010 R8000 GFP constructs post-normalization", x = "Plasmid", y = "Median RNA/DNA") + geom_point(data = R8000_GFP_RD, aes(plasmid, RNADNA), alpha = 0.3, show.legend = FALSE) +
  theme(strip.background=element_rect(fill = "White")) + geom_text(data = R8000_GFP_RD, aes(plasmid, RNADNA, label = paste(plasmid,"-",number, sep = "")), nudge_x = 0.25, size = 2.5)
g.GFP

norm = sum(filter(sumGFP, time == 0)$medRD)/sum(filter(sumGFP, time == 6)$medRD)

#############################weird neural net exploration#################################
seqR800 <- read.table("170307 all-8k R8000.fasta", header = FALSE)
peptide <- filter(seqR800, grepl(">", V1)) 
sequence <- filter(seqR800, !grepl(">",V1))
seqR8000 <- data.frame("peptide" = peptide,"sequence" = sequence)
colnames(seqR8000) <- c("peptide","sequence")
seqR8000 <- separate(seqR8000, peptide, into = c("peptide", "X_group"), sep = "-")
seq9000 <- seqR8000
position <- c("f0","g0","a1","b1","c1","d1","e1","f1","g1","a2","b2","c2","d2","e2","f2","g2","a3","b3","c3","d3","e3","f3","g3","a4","b4","c4","d4","e4","f4")

kmer = 3
for (i in seq(1, length(position)-kmer, by = 3)){
  seq9001 <- unite_(seq9000, position[i] ,i:i+2, sep = "")
}

seq9001 <- seq9000
seq9001 <- unite(seq9001, col = '9', 11:13, sep = "")

for (i in 1:length(position)){
  seq9000 <- separate(seq9000, sequence, into = c(position[i],"sequence"), sep = 1)
}
seq9000 <- separate(seq9000, peptide, into = c("junk","peptide"), sep = ">")
seq9000 <- select(seq9000,-junk, -sequence)
seq9000 <- filter(seq9000, peptide != "GCN4")

sumRD <- left_join(sumR100RD, seq9000, by = c("X_peptide" = "peptide", "X_group"))
sumRD <- left_join(sumRD, seq9000, by = c("Y_peptide" = "peptide", "X_group"))

set.seed(123)
sumRD$X_peptide <- as.factor(sumRD$X_peptide)
sumRD$Y_peptide <- as.factor(sumRD$Y_peptide)
sumRD2 <- sumRD
#sumRD <- mutate_each_(sumRD, funs(factor(.)), c(17:74))
sumRD[,17:74] <- lapply(sumRD[,17:74], factor)
sumRD <- ungroup(sumRD)
sumRD <- filter(sumRD, bcnum > 5, sum_DNA > 4, condition == 6)
sumRD <- mutate(sumRD, "log_sum_RD" = log(sum_RD))

train <- sample_frac(sumRD, size= 0.5, replace = FALSE)
test2 <- anti_join(sumRD, train)


### This seem to have too many variables to fit nicely? fit <- nnet(log_sum_RD ~  g0.x + a1.x + b1.x +c1.x + d1.x + e1.x + f1.x + g1.x + a2.x + b2.x + c2.x + d2.x + e2.x +f2.x +g2.x + a3.x + b3.x + c3.x + d3.x + e3.x + f3.x + g3.x + a4.x + b4.x + c4.x+ e4.x + f4.x + g0.y + a1.y + b1.y +c1.y + d1.y + e1.y + f1.y + g1.y + a2.y + b2.y + c2.y + d2.y + e2.y +f2.y +g2.y + a3.y + b3.y + c3.y + d3.y + e3.y + f3.y + g3.y + a4.y + b4.y + c4.y + e4.y + f4.y, data = train, size = 10, maxit = 500, linout=T, decay=0.01, MaxNWts = 15000)
fit <- nnet(sum_RD ~  g0.x + a1.x + b1.x +c1.x + e1.x + f1.x + g1.x + a2.x + b2.x + c2.x  + e2.x +f2.x +g2.x + a3.x + b3.x + c3.x + e3.x + f3.x + g3.x + a4.x + b4.x + c4.x+ e4.x + f4.x + g0.y + a1.y + b1.y +c1.y + e1.y + f1.y + g1.y + a2.y + b2.y + c2.y + e2.y +f2.y +g2.y + a3.y + b3.y + c3.y + e3.y + f3.y + g3.y + a4.y + b4.y + c4.y + e4.y + f4.y, data = train, size = 10, maxit = 500, linout=T, decay=0.01, MaxNWts = 15000)


fit <- nnet(sum_RD ~  g0.x + a1.x:a1.y + b1.x +c1.x + e1.x + f1.x + g1.x + a2.x + b2.x + c2.x  + e2.x +f2.x +g2.x + a3.x + b3.x + c3.x + e3.x + f3.x + g3.x + a4.x + b4.x + c4.x+ e4.x + f4.x + g0.y  + b1.y +c1.y + e1.y + f1.y + g1.y + a2.y + b2.y + c2.y + e2.y +f2.y +g2.y + a3.y + b3.y + c3.y + e3.y + f3.y + g3.y + a4.y + b4.y + c4.y + e4.y + f4.y, data = train, size = 10, maxit = 500, linout=T, decay=0.01, MaxNWts = 15000)


fit <- nnet(sum_RD ~  g0.x * a1.x * b1.x * c1.x * e1.x * f1.x * g1.x * a2.x * b2.x * c2.x  * e2.x * f2.x * g2.x * a3.x * b3.x * c3.x * e3.x * f3.x * g3.x * a4.x * b4.x * c4.x * e4.x * f4.x * g0.y * a1.y * b1.y * c1.y * e1.y * f1.y * g1.y * a2.y * b2.y * c2.y * e2.y * f2.y * g2.y * a3.y * b3.y * c3.y * e3.y * f3.y * g3.y * a4.y * b4.y * c4.y * e4.y * f4.y, data = train, size = 10, maxit = 500, linout=T, decay=0.01, MaxNWts = 15000)
fit <- nnet(log_sum_RD ~  g0.x + a1.x  + d1.x + e1.x + g1.x + a2.x + d2.x + e2.x +g2.x + a3.x + d3.x + e3.x + g3.x + a4.x + e4.x + g0.y + a1.y + d1.y + e1.y + g1.y + a2.y  + d2.y + e2.y +g2.y + a3.y  + d3.y + e3.y + g3.y + a4.y + e4.y, data = train, size = 10, maxit = 500, linout=T, decay=0.01, MaxNWts = 15000)
test2$predict <- predict(fit, test2, type = "raw")     
corr4 <- data.frame(x = 0, y = 1, label =  paste("italic(r)==", round(cor(test2$predict, test2$sum_RD, use = "complete.obs"), digits = 3)))

g.neural <- ggplot(test2, aes(sum_RD, predict)) + geom_point(alpha = 0.1) +geom_text( data = corr4, aes(x= x,y = y,label= label), parse = TRUE)

g.neural

fit <- nnet(log_sum_RD ~ . -..., data=..., ...)
