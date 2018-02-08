################################################# Standard libraries + viridis to make things pretty + extrafont for comic sans #################################################
library(dplyr)
library(tidyr)
library(extrafont)
library(viridis)
library(cowplot)
################################################# Read in the data and mash it together #################################################
setwd("/Users/Cliff/Documents/Kosuri Lab/NGS files/171113 R8000-3/")

R8000_RNA_0 <- read.table("sc-counts_R8000-3_0h_RNA.txt", header = FALSE)
R8000_RNA_0 <- data.frame(R8000_RNA_0, rep("RNA-0",nrow(R8000_RNA_0)),rep(sum(R8000_RNA_0$V2), nrow(R8000_RNA_0)))
colnames(R8000_RNA_0) <- c("barcode", "count", "collapsed_bc", "condition","readsinlib")
R8000_RNA_0 <- mutate(R8000_RNA_0, count =  count/norm)

R8000_DNA_0 <- read.table("sc-counts_R8000-3_0h_DNA.txt", header = FALSE)
R8000_DNA_0 <- data.frame(R8000_DNA_0, rep("DNA-0",nrow(R8000_DNA_0)),rep(sum(R8000_DNA_0$V2), nrow(R8000_DNA_0)))
colnames(R8000_DNA_0) <- c("barcode", "count", "collapsed_bc", "condition","readsinlib")

R8000_RNA_6 <- read.table("sc-counts_R8000-3_6h_RNA.txt", header = FALSE)
R8000_RNA_6 <- data.frame(R8000_RNA_6, rep("RNA-6",nrow(R8000_RNA_6)),rep(sum(R8000_RNA_6$V2), nrow(R8000_RNA_6)))
colnames(R8000_RNA_6) <- c("barcode", "count", "collapsed_bc", "condition","readsinlib")

R8000_DNA_6 <-read.table("sc-counts_R8000-3_6h_DNA.txt", header = FALSE)
R8000_DNA_6 <- data.frame(R8000_DNA_6, rep("DNA-6",nrow(R8000_DNA_6)),rep(sum(R8000_DNA_6$V2), nrow(R8000_DNA_6)))
colnames(R8000_DNA_6) <- c("barcode", "count", "collapsed_bc", "condition","readsinlib")

############################## get first run ################################################

setwd("/Users/Cliff/Documents/Kosuri Lab/NGS files/170926 R8000 barcode-1/")
load_obj <- function(f)
{
  load(f)
  get(ls()[ls() != "f"])
}
R8000_1_full <- load_obj("/Users/Cliff/Documents/Kosuri Lab/NGS files/170926 R8000 barcode-1/R8000_full_DF.Rda")
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

corr0 <- data.frame(x = 0.5, y = 300, label =  paste("italic(r)==", round(cor(filter(insane_R8000, dcount.x > 1, dcount.y > 1, condition == 0)$RNADNA.x, filter(insane_R8000, dcount.x > 1, dcount.y > 1, condition == 0)$RNADNA.y, use = "complete.obs"), digits = 3)))
corr0.5 <- data.frame(x = 5, y = 400, label =  paste("italic(r)==", round(cor(filter(insane_R8000, dcount.x > 1, dcount.y > 1, condition == 6)$RNADNA.x, filter(insane_R8000, dcount.x > 1, dcount.y > 1, condition == 6)$RNADNA.y, use = "complete.obs"), digits = 3)))

g.rawBios <- ggplot(filter(insane_R8000, dcount.x > 1, dcount.y > 1, condition == 0), aes(RNADNA.x, RNADNA.y)) + geom_point(alpha = 0.2) + scale_x_log10() + scale_y_log10() + labs(x =  "Replicate 2", y = "Replicate 1", title = "RNA/DNA 0h") +
  geom_text( data = corr0, aes(x= x,y = y,label= label), parse = TRUE)
g.rawBios

g.rawBios5 <- ggplot(filter(insane_R8000, dcount.x > 1, dcount.y > 1, condition == 6), aes(RNADNA.x, RNADNA.y)) + geom_point(alpha = 0.2) + scale_x_log10() + scale_y_log10() + labs(x =  "Replicate 2", y = "Replicate 1", title = "RNA/DNA 6h") +
  geom_text( data = corr0.5, aes(x= x,y = y,label= label), parse = TRUE)
g.rawBios5

plot_grid(g.rawBios, g.rawBios5)

corr1 <- data.frame(x = 5, y = 100, label =  paste("italic(r)==", round(cor(filter(insane_R8000, dcount.x > 1, dcount.y > 1, condition == 0)$dcount.x, filter(insane_R8000, dcount.x > 1, dcount.y > 1, condition == 0)$dcount.y, use = "complete.obs"), digits = 3)))
corr2 <- data.frame(x = 0.5, y = 100, label =  paste("italic(r)==", round(cor(filter(insane_R8000, dcount.x > 1, dcount.y > 1, condition == 0)$rcount.x, filter(insane_R8000, dcount.x > 1, dcount.y > 1, condition == 0)$rcount.y, use = "complete.obs"), digits = 3)))
corr3 <- data.frame(x = 5, y = 100, label =  paste("italic(r)==", round(cor(filter(insane_R8000, dcount.x > 1, dcount.y > 1, condition == 6)$dcount.x, filter(insane_R8000, dcount.x > 1, dcount.y > 1, condition == 6)$dcount.y, use = "complete.obs"), digits = 3)))
corr4 <- data.frame(x = 0.5, y = 100, label =  paste("italic(r)==", round(cor(filter(insane_R8000, dcount.x > 1, dcount.y > 1, condition == 6)$rcount.x, filter(insane_R8000, dcount.x > 1, dcount.y > 1, condition == 6)$rcount.y, use = "complete.obs"), digits = 3)))


g.rawBios1 <- ggplot(filter(insane_R8000, dcount.x > 1, dcount.y > 1, condition == 0), aes(dcount.x, dcount.y)) + geom_point(alpha = 0.2) + scale_x_log10() + scale_y_log10() + labs(x =  "Replicate 2", y = "Replicate 1", title = "DNA 0h") +
  geom_text( data = corr1, aes(x= x,y = y,label= label), parse = TRUE)
g.rawBios1
g.rawBios2 <- ggplot(filter(insane_R8000, rcount.x > 0.1, rcount.y > 0.1, condition == 0), aes(rcount.x, rcount.y)) + geom_point(alpha = 0.2) + scale_x_log10() + scale_y_log10() + labs(x =  "Replicate 2", y = "Replicate 1", title = "RNA 0h") +
  geom_text( data = corr2, aes(x= x,y = y,label= label), parse = TRUE)
g.rawBios2
g.rawBios3 <- ggplot(filter(insane_R8000, dcount.x > 1, dcount.y > 1, condition == 6), aes(dcount.x, dcount.y)) + geom_point(alpha = 0.2) + scale_x_log10() + scale_y_log10() + labs(x =  "Replicate 2", y = "Replicate 1", title = "DNA 6h") +
  geom_text( data = corr3, aes(x= x,y = y,label= label), parse = TRUE)
g.rawBios3
g.rawBios4 <- ggplot(filter(insane_R8000, rcount.x > 0.1, rcount.y > 0.1, condition == 6), aes(rcount.x, rcount.y)) + geom_point(alpha = 0.2) + scale_x_log10() + scale_y_log10() + labs(x =  "Replicate 2", y = "Replicate 1", title = "RNA 6h") +
  geom_text( data = corr4, aes(x= x,y = y,label= label), parse = TRUE)
g.rawBios4

plot_grid(g.rawBios1, g.rawBios2, g.rawBios3, g.rawBios4)

insane_sum <- group_by(insane_R8000, X_peptide, X_group, Y_peptide, condition) %>%
  summarise(mean_RDx = mean(RNADNA.x), mean_RDy = mean(RNADNA.y), med_RDx = median(RNADNA.x), med_RDy = median(RNADNA.y), sum_DNAx = sum(dcount.x), sum_DNAy = sum(dcount.y), sum_RNAx = sum(rcount.x), sum_RNAy = sum(rcount.y), sum_RDx = sum(rcount.x)/(sum(dcount.x) + 0.1), sum_RDy = sum(rcount.y)/(sum(dcount.y) + 0.1), bcnum = n_distinct(barcode), SEMx = sd(RNADNA.x)/sqrt(n()))

g.Biomean <- ggplot(filter(insane_sum, sum_DNAx > 10, sum_DNAy > 10, bcnum > 5), aes(mean_RDx, mean_RDy)) + geom_point(alpha = 0.2) + facet_wrap(~condition, scales = "free") + labs(x = "Mean RNA/DNA Biological Replicate 2", y = "Mean RNA/DNA Biological Replicate 1", title = "171115 R8000 biological replicates mean RNA/DNA") +
      scale_x_log10() + scale_y_log10()
g.Biomean

g.Biomean <- ggplot(filter(insane_sum, sum_DNAx > 10, sum_DNAy > 10, bcnum > 5), aes(sum_RDx, sum_RDy, color = log2(bcnum))) + geom_point(alpha = 0.2) + facet_wrap(~condition, scales = "free", labeller = labeller(condition=induction_labels) ) + labs(x = "Mean RNA/DNA Biological Replicate 2", y = "Mean RNA/DNA Biological Replicate 1", title = "171115 R8000 biological replicates summed RNA/DNA by construct") +
  scale_x_log10(limits = c(0.003, 5)) + scale_y_log10(limits = c(0.2, 12)) + scale_color_viridis(name = "Number of\nbarcodes \nlog_2") + theme(strip.background=element_rect(fill = "White")) 
g.Biomean

s
co1 <-  data.frame(x = 50, y = 100, label =  paste("italic(r)==", round(cor(filter(insane_sum, sum_DNAx > 10, sum_DNAy > 10, bcnum > 5, condition == 0)$sum_DNAx, filter(insane_sum, sum_DNAx > 10, sum_DNAy > 10, bcnum > 5, condition == 0)$sum_DNAy, use = "complete.obs"), digits = 3)))
co2 <-  data.frame(x = 10, y = 500, label =  paste("italic(r)==", round(cor(filter(insane_sum, sum_DNAx > 10, sum_DNAy > 10, bcnum > 5, condition == 0)$sum_RNAx, filter(insane_sum, sum_DNAx > 10, sum_DNAy > 10, bcnum > 5, condition == 0)$sum_RNAy, use = "complete.obs"), digits = 3)))
co3 <-  data.frame(x = 50, y = 150, label =  paste("italic(r)==", round(cor(filter(insane_sum, sum_DNAx > 10, sum_DNAy > 10, bcnum > 5, condition == 6)$sum_DNAx, filter(insane_sum, sum_DNAx > 10, sum_DNAy > 10, bcnum > 5, condition == 6)$sum_DNAy, use = "complete.obs"), digits = 3)))
co4 <-  data.frame(x = 10, y = 500, label =  paste("italic(r)==", round(cor(filter(insane_sum, sum_DNAx > 10, sum_DNAy > 10, bcnum > 5, condition == 6)$sum_RNAx, filter(insane_sum, sum_DNAx > 10, sum_DNAy > 10, bcnum > 5, condition == 6)$sum_RNAy, use = "complete.obs"), digits = 3)))




g.Biomean1 <- ggplot(filter(insane_sum, sum_DNAx > 10, sum_DNAy > 10, bcnum > 5, condition == 0), aes(sum_DNAx, sum_DNAy)) + geom_point(alpha = 0.2)  + labs(x = "Summed DNA Biological Replicate 2", y = "Summed DNA Biological Replicate 1", title = "DNA 0h") +
  scale_x_log10() + scale_y_log10() + geom_text( data = co1, aes(x= x,y = y,label= label), parse = TRUE)
g.Biomean1

g.Biomean2 <- ggplot(filter(insane_sum, sum_DNAx > 10, sum_DNAy > 10, bcnum > 5, condition == 0), aes(sum_RNAx, sum_RNAy)) + geom_point(alpha = 0.2)  + labs(x = "Summed RNA Biological Replicate 2", y = "Summed RNA Biological Replicate 1", title = "RNA 0h") +
  scale_x_log10() + scale_y_log10() + geom_text( data = co2, aes(x= x,y = y,label= label), parse = TRUE)
g.Biomean2

g.Biomean3 <- ggplot(filter(insane_sum, sum_DNAx > 10, sum_DNAy > 10, bcnum > 5, condition == 6), aes(sum_DNAx, sum_DNAy)) + geom_point(alpha = 0.2)  + labs(x = "Summed DNA Biological Replicate 2", y = "Summed DNA Biological Replicate 1", title = "DNA 6h") +
  scale_x_log10() + scale_y_log10() + geom_text( data = co3, aes(x= x,y = y,label= label), parse = TRUE)
g.Biomean3

g.Biomean4 <- ggplot(filter(insane_sum, sum_DNAx > 10, sum_DNAy > 10, bcnum > 5, condition == 6), aes(sum_RNAx, sum_RNAy)) + geom_point(alpha = 0.2)  + labs(x = "Summed RNA Biological Replicate 2", y = "Summed RNA Biological Replicate 1", title = "RNA 6h") +
  scale_x_log10() + scale_y_log10() + geom_text( data = co4, aes(x= x,y = y,label= label), parse = TRUE)
g.Biomean4


plot_grid(g.Biomean1, g.Biomean2, g.Biomean3, g.Biomean4)
################################################################################################################################################
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
R8000DNA <- filter(R8000_full_DF, grepl("DNA", R8000_full_DF$condition)) %>%
  select(X_peptide, X_group, X_codon, Y_peptide, condition, barcode, "dcount" = nCount ) %>%
  separate(condition, into = c("template", "condition"))
R8000RD <- right_join(R8000RNA,R8000DNA, by = c("X_peptide","Y_peptide","barcode","condition", "X_group", "X_codon")) %>%
  mutate("RNADNA" = rcount/dcount*100) %>%
  select(-template.x, -template.y) %>%
  mutate("RNADNA" = ifelse(is.na(RNADNA),0,RNADNA)) %>%
  mutate("rcount" =  ifelse(is.na(rcount),0,rcount))

sumR100RD <- group_by(R8000RD, X_peptide, Y_peptide, condition, X_group) %>%
  summarise(mean_RNA = mean(rcount), sd_R = sd(rcount), mean_DNA = mean(dcount), sd_D = sd(dcount), mean_RD = mean(RNADNA), sd_RD = sd(RNADNA), med_RD = median(RNADNA), sum_DNA = sum(dcount), sum_RNA = sum(rcount), sum_RD = sum(rcount)/(sum(dcount) +1), bcnum = n_distinct(barcode), SEM = sd(RNADNA)/sqrt(n()))

induction_labels <- c("0" = "0 hours induction RNA/DNA", "6" = "6 hours induction RNA/DNA")

g.bccount <- ggplot(sumR100RD, aes(bcnum, fill = mean_RD)) + geom_histogram(bins = 114) + scale_x_continuous(limits = c(0,113)) + labs(title = "171113 R8000-3 Number of barcodes per construct", x = "number of barcodes")  +
  scale_fill_viridis()
g.bccount

g.sumOs <- ggplot(filter(sumR100RD, grepl("O",X_peptide)), aes(X_peptide, Y_peptide)) + geom_tile(aes(fill=sum_RD)) + facet_grid(X_group~condition, labeller = labeller(condition=induction_labels), scales = "free") +
  theme(text=element_text(family="Calibri"), strip.text=element_text(size = 10), strip.background=element_rect(fill = "White"), axis.text.x = element_text(angle=90,size =5),
        axis.text.y = element_text(size = 5), legend.title=element_text(size = 8), legend.title.align=0.5, legend.text=element_text(size =8)) +
  scale_fill_viridis(option="viridis", name = "Mean \nRNA to DNA  ", na.value="white") + labs(title="171113 Roman's 8000-3 RNA/DNA mean barcode ratio O constructs (bc >0)", x="X peptide", y = "Y peptide")
g.sumOs





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
g.GFP <- ggplot(sumGFP, aes(plasmid,meanRD, fill=plasmid)) + geom_bar(stat = "identity") + geom_errorbar(gfperr) + facet_wrap(~time, labeller = labeller(time = plasmid_labels)) + 
  labs(title = "171113 R8000-3 GFP constructs pre-normalization", x = "Plasmid", y = "Median RNA/DNA") + geom_point(data = R8000_GFP_RD, aes(plasmid, RNADNA), alpha = 0.3, show.legend = FALSE) +
  theme(strip.background=element_rect(fill = "White")) #+ geom_text(data = R8000_GFP_RD, aes(plasmid, RNADNA, label = paste(plasmid,"-",number, sep = "")), nudge_x = 0.25, size = 2.5)
g.GFP

norm = sum(filter(sumGFP, time == 0)$medRD)/sum(filter(sumGFP, time == 6)$medRD)
