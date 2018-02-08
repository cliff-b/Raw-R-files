library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(viridis)
library(extrafont)
setwd("/Users/Cliff/Documents/Kosuri Lab/NGS files/160831 R20 Synzip run 1")

##############################Read in data from sequencing run and mapping files######################################
load("synzipbigDF.Rda")
load ("synzipmapBCs.Rda")
load("GFPBCs.Rda")

#####Remove DNA barcodes with less than 10 reads to lower noise 
#####Map with both the mapping run and known GFP barcodes
#####Get rid of non-mapped reads and rearrange for asthetics
bigDF <- filter(bigDF, (count > 5) | grepl("RNA", bigDF$condition))
bigDF <- left_join(bigDF, mapBCs, by = "RevCompBC")
bigDF <- left_join(bigDF, GFPBCs, by = "RevCompBC") %>%
  filter(!is.na(X_peptide) | !is.na(plasmid)) %>%
  select(X_peptide, Y_peptide, plasmid, condition, RevCompBC, count, nCount)

######Divide bigDF into RNA and DNA
######Merge RNA into DNA by construct and condition
######Replace NAs from having DNA but no RNA with 0s
R20RNA <- filter(bigDF, grepl("RNA", bigDF$condition)) %>%
  select(X_peptide, Y_peptide, plasmid, condition, RevCompBC, "rcount" = nCount ) %>%
  separate(condition, into = c("template", "condition"))
R20DNA <- filter(bigDF, grepl("DNA", bigDF$condition)) %>%
  select(X_peptide, Y_peptide, plasmid, condition, RevCompBC, "dcount" = nCount ) %>%
  separate(condition, into = c("template", "condition"))
R20RD <- right_join(R20RNA,R20DNA, by = c("RevCompBC","X_peptide","Y_peptide","plasmid","condition")) %>%
  mutate("RNADNA" = rcount/dcount*100) %>%
  separate("X_peptide", into = c("X_peptide", "X_codon"), sep = "\\[") %>%
  separate("X_peptide", into = c("X_peptide", "X_rev"), sep = "_") %>%
  separate("Y_peptide", into = c("Y_peptide", "Y_codon"), sep = "\\[") %>%
  select(-X_codon, -Y_codon, -X_rev, -template.x, -template.y)
R20RD$RNADNA[is.na(R20RD$RNADNA)] <- 0
R20RD$rcount[is.na(R20RD$rcount)] <- 0

################## helper functions to make pretty graphs #######################
R20RD$X_peptide <- as.factor(R20RD$X_peptide)
R20RD$Y_peptide <- as.factor(R20RD$Y_peptide)
R20RD$X_peptide <- factor(R20RD$X_peptide, levels = c(">GCN4", ">P1", ">P2",">P3",">P4",">P5",">P6",">P7",">P8",">P9",">P10",">P11",">P12",">P3mS",">P4mS",">P5mS",">P6mS",">P7mS",">P8mS",">P5SC1",">P5SC2",">P6SC1",">P6SC2"))
R20RD$Y_peptide <- factor(R20RD$Y_peptide, levels = rev(c(">GCN4", ">P1", ">P2",">P3",">P4",">P5",">P6",">P7",">P8",">P9",">P10",">P11",">P12",">P3mS",">P4mS",">P5mS",">P6mS",">P7mS",">P8mS",">P5SC1",">P5SC2",">P6SC1",">P6SC2")))
induction_labels2 <- c("0" = "Uninduced RNA/DNA", "A" = "5uM DAPG 5ng/mL ATc RNA/DNA", "B" = "7.5uM DAPG 10ng/mL ATc RNA/DNA", "C" = "15uM DAPG 40ng/mL ATc RNA/DNA")


sumR20D5 <- group_by(R20RD, X_peptide, Y_peptide, condition) %>%
  summarise(mean_RNA = mean(rcount), sd_R = sd(rcount), mean_DNA = mean(dcount), sd_D = sd(dcount), mean_RD = mean(RNADNA), sd_RD = sd(RNADNA), med_RD = median(RNADNA), bcnum = n_distinct(RevCompBC), SEM = sd(RNADNA)/sqrt(length(RNADNA)))


######################################################################################################################################
##############################################Create the graphs based off the summary data############################################
######################################################################################################################################

######### mean RNA/DNA ratio. first, filter to get rid of the GFP constitutive BCs ############
g.sumR20RD <- ggplot(filter(sumR20D5, !is.na(X_peptide)), aes(X_peptide, Y_peptide)) + geom_tile(aes(fill=mean_RD)) + facet_wrap(~condition, labeller = labeller(condition=induction_labels2), scales = "free") +
  theme(text=element_text(family="Calibri"), strip.text=element_text(size = 10), strip.background=element_rect(fill = "White"), axis.text.x = element_text(angle=90, hjust=1,size =6),
        axis.text.y = element_text(size = 6), legend.title=element_text(size = 8), legend.title.align=0.5, legend.text=element_text(size =8)) +
  scale_fill_viridis(option="viridis", name = "Normalized \nRNA to DNA ", na.value="white") + labs(title="160913 Synzip RNA/DNA barcode ratio", x="X peptide", y = "Y peptide")
g.sumR20RD

############ median RNA/DNA ratio ############
g.sumR20RDmed <- ggplot(filter(sumR20D5, !is.na(X_peptide) & SEM*6 < mean_RD ), aes(X_peptide, Y_peptide)) + geom_tile(aes(fill=med_RD)) + facet_wrap(~condition, labeller = labeller(condition=induction_labels2), scales = "free") +
  theme(text=element_text(family="Calibri"), strip.text=element_text(size = 10), strip.background=element_rect(fill = "White"), axis.text.x = element_text(angle=90, hjust=1,size =6),
        axis.text.y = element_text(size = 6), legend.title=element_text(size = 8), legend.title.align=0.5, legend.text=element_text(size =8)) +
  scale_fill_viridis(option="viridis", name = "Median \nRNA to DNA ", na.value="white") + labs(title="160913 Synzip RNA/DNA barcode ratio (median)", x="X peptide", y = "Y peptide")
g.sumR20RDmed

############## barcode counts per XY pair ############
g.sumR20RDbc <- ggplot(filter(sumR20D5, !is.na(X_peptide)), aes(X_peptide, Y_peptide)) + geom_tile(aes(fill=bcnum)) + facet_wrap(~condition, labeller = labeller(condition=induction_labels2), scales = "free") +
  theme(text=element_text(family="Calibri"), strip.text=element_text(size = 10), strip.background=element_rect(fill = "White"), axis.text.x = element_text(angle=90, hjust=1,size =6),
        axis.text.y = element_text(size = 6), legend.title=element_text(size = 8), legend.title.align=0.5, legend.text=element_text(size =8)) +
  scale_fill_viridis(option="viridis", name = "Number of\n barcodes ", na.value="white") + labs(title="160913 Synzip number of barcodes per construct", x="X peptide", y = "Y peptide")
g.sumR20RDbc

g.sdvmean <- ggplot(filter(sumR20D5, !is.na(X_peptide)), aes(mean_RD, sd_RD, color=condition)) +geom_point(alpha=0.4) + scale_y_continuous(limits = c(0,400))+ scale_x_continuous(limits=c(0,400))+ #, labels = c("Uninduced RNA/DNA","5uM DAPG 5ng/mL ATc RNA/DNA","7.5uM DAPG 10ng/mL ATc RNA/DNA","15uM DAPG 40ng/mL ATc RNA/DNA")) + 
  labs(title="160913 Synzip Standard Deviation vs mean of RNA/DNA ratio", y="Standard Deviation", x = "Mean RNA/DNA ratio")
g.sdvmean

g.medvmean <- ggplot(filter(sumR20D5, !is.na(X_peptide)), aes(mean_RD, med_RD, color=condition)) +geom_point(alpha=0.4) + #scale_y_continuous(limits = c(0,400))+ scale_x_continuous(limits=c(0,400))+ #, labels = c("Uninduced RNA/DNA","5uM DAPG 5ng/mL ATc RNA/DNA","7.5uM DAPG 10ng/mL ATc RNA/DNA","15uM DAPG 40ng/mL ATc RNA/DNA")) + 
  labs(title="160913 Synzip mean vs median of RNA/DNA ratio", y="Median", x = "Mean RNA/DNA ratio")
g.medvmean
