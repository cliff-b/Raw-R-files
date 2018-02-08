library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(viridis)
library(extrafont)
setwd("/Users/Cliff/Documents/Kosuri Lab/NGS files/160831 R20 Synzip run 1")

##############################Read in data from sequencing run and mapping files######################################
load("R20bigDF.Rda")
load ("R20mappedBCs.Rda")
load("GFPBCs.Rda")
  
#####Remove DNA barcodes with less than 10 reads to lower noise 
#####Map with both the mapping run and known GFP barcodes
#####Get rid of non-mapped reads and rearrange for asthetics
bigDF <- filter(bigDF, (count > 10) | grepl("RNA", bigDF$condition))
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
  separate("Y_peptide", into = c("Y_peptide", "Y_codon"), sep = "_") %>%
  separate ("X_peptide", into = c("X_peptide", "X_codon"), sep = "_") %>%
  select(-X_codon, -Y_codon, -template.x, -template.y)
R20RD$RNADNA[is.na(R20RD$RNADNA)] <- 0
R20RD$rcount[is.na(R20RD$rcount)] <- 0


######Create values for filtering
######Get rid of outliers
R20RD <- R20RD %>% 
  group_by(X_peptide, Y_peptide, condition) %>% 
  mutate("dstdev" = sd(dcount), "dmean" = mean(dcount), "rstdev" = sd(rcount), "rmean" = mean(rcount)) %>% 
  ungroup(R20RD)
R20RD <-  filter(R20RD, dcount < dmean + 3*dstdev & dcount > dmean - 3*dstdev & rcount < rmean +3*rstdev & rcount > rmean - 3*rstdev)

###################################Graphs for all barcodes##############################################
R20RD <- mutate(R20RD, "combo" = paste(X_peptide,Y_peptide))
g.RDbyCon <- ggplot(filter(R20RD, X_peptide==">P4" & condition=="A"), aes(combo, rcount)) + geom_boxplot(outlier.shape=NA, color="red") + geom_jitter(alpha=0.4) +theme(axis.text.x = element_text(angle=90, hjust=1,size =6))
g.RDbyCon

g.RbyCon <- ggplot(filter(R20RD, X_peptide==">P4" & condition !=0), aes(combo, rcount, color=condition)) + geom_boxplot(outlier.shape=NA, color="red") + geom_jitter(alpha=0.4) +theme(axis.text.x = element_text(angle=90, hjust=1,size =6)) +
  labs(title="160913 Roman's 20 RNA counts by construct", x="", y = "Reads (RNA)")
g.RbyCon

g.DbyCon <- ggplot(filter(R20RD, X_peptide==">P4" & condition !=0), aes(combo, dcount, color=condition)) + geom_boxplot(outlier.shape=NA, color="red") + geom_jitter(alpha=0.4) +theme(axis.text.x = element_text(angle=90, hjust=1,size =6)) +
  labs(title="160913 Roman's 20 DNA counts by construct", x="", y = "Reads (DNA)")
g.DbyCon

g.DbyCon <- ggplot(filter(R20RD, X_peptide==">P4" & condition !=0), aes(combo, RNADNA, color=condition)) + geom_boxplot(outlier.shape=NA, color="red") + geom_jitter(alpha=0.4) +theme(axis.text.x = element_text(angle=90, hjust=1,size =6)) +
  labs(title="160913 Roman's 20 RNA/DNA counts by construct", x="", y = "Reads (RNA/DNA)")
g.DbyCon

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
  scale_fill_viridis(option="viridis", name = "Normalized \nRNA to DNA ", na.value="white") + labs(title="160913 Roman's 20 RNA/DNA barcode ratio", x="X peptide", y = "Y peptide")
g.sumR20RD

############ median RNA/DNA ratio ############
g.sumR20RDmed <- ggplot(filter(sumR20D5, !is.na(X_peptide)), aes(X_peptide, Y_peptide)) + geom_tile(aes(fill=med_RD)) + facet_wrap(~condition, labeller = labeller(condition=induction_labels2), scales = "free") +
  theme(text=element_text(family="Calibri"), strip.text=element_text(size = 10), strip.background=element_rect(fill = "White"), axis.text.x = element_text(angle=90, hjust=1,size =6),
        axis.text.y = element_text(size = 6), legend.title=element_text(size = 8), legend.title.align=0.5, legend.text=element_text(size =8)) +
  scale_fill_viridis(option="viridis", name = "Median \nRNA to DNA ", na.value="white") + labs(title="160913 Roman's 20 RNA/DNA barcode ratio (median)", x="X peptide", y = "Y peptide")
g.sumR20RDmed

############ standard deviation RNA/DNA ratio ############
g.sumR20RDsd <- ggplot(filter(sumR20D5, !is.na(X_peptide)), aes(X_peptide, Y_peptide)) + geom_tile(aes(fill=sd_RD)) + facet_wrap(~condition, labeller = labeller(condition=induction_labels2), scales = "free") +
  theme(text=element_text(family="Calibri"), strip.text=element_text(size = 10), strip.background=element_rect(fill = "White"), axis.text.x = element_text(angle=90, hjust=1,size =6),
        axis.text.y = element_text(size = 6), legend.title=element_text(size = 8), legend.title.align=0.5, legend.text=element_text(size =8)) +
  scale_fill_viridis(option="viridis", name = "Standard deviation \n of RNA to DNA ", na.value="white") + labs(title="160913 Roman's 20 RNA/DNA standard deviation of barcode ratio", x="X peptide", y = "Y peptide")
g.sumR20RDsd

############ standard error RNA/DNA ratio ############
g.sumR20RDse <- ggplot(filter(sumR20D5, !is.na(X_peptide)), aes(X_peptide, Y_peptide)) + geom_tile(aes(fill=SEM)) + facet_wrap(~condition, labeller = labeller(condition=induction_labels2), scales = "free") +
  theme(text=element_text(family="Calibri"), strip.text=element_text(size = 10), strip.background=element_rect(fill = "White"), axis.text.x = element_text(angle=90, hjust=1,size =6),
        axis.text.y = element_text(size = 6), legend.title=element_text(size = 8), legend.title.align=0.5, legend.text=element_text(size =8)) +
  scale_fill_viridis(option="viridis", name = "Standard error of mean \n of RNA to DNA ", na.value="white") + labs(title="160913 Roman's 20 RNA/DNA standard error of the mean of barcode ratio", x="X peptide", y = "Y peptide")
g.sumR20RDse

############## barcode counts per XY pair ############
g.sumR20RDbc <- ggplot(filter(sumR20D5, !is.na(X_peptide)), aes(X_peptide, Y_peptide)) + geom_tile(aes(fill=bcnum)) + facet_wrap(~condition, labeller = labeller(condition=induction_labels2), scales = "free") +
  theme(text=element_text(family="Calibri"), strip.text=element_text(size = 10), strip.background=element_rect(fill = "White"), axis.text.x = element_text(angle=90, hjust=1,size =6),
        axis.text.y = element_text(size = 6), legend.title=element_text(size = 8), legend.title.align=0.5, legend.text=element_text(size =8)) +
  scale_fill_viridis(option="viridis", name = "Number of\n barcodes ", na.value="white") + labs(title="160913 Roman's 20 number of barcodes per construct", x="X peptide", y = "Y peptide")
g.sumR20RDbc

g.sdvmean <- ggplot(filter(sumR20D5, !is.na(X_peptide)), aes(mean_RD, sd_RD, color=condition)) +geom_point(alpha=0.4) + scale_y_continuous(limits = c(0,400))+ scale_x_continuous(limits=c(0,400))+ #, labels = c("Uninduced RNA/DNA","5uM DAPG 5ng/mL ATc RNA/DNA","7.5uM DAPG 10ng/mL ATc RNA/DNA","15uM DAPG 40ng/mL ATc RNA/DNA")) + 
  labs(title="160913 Roman's 20 Standard Deviation vs mean of RNA/DNA ratio", y="Standard Deviation", x = "Mean RNA/DNA ratio")
g.sdvmean

g.medvmean <- ggplot(filter(sumR20D5, !is.na(X_peptide) & condition == 0), aes(mean_RD, med_RD, color=condition)) +geom_point(alpha=0.4) + #scale_y_continuous(limits = c(0,400))+ scale_x_continuous(limits=c(0,400))+ #, labels = c("Uninduced RNA/DNA","5uM DAPG 5ng/mL ATc RNA/DNA","7.5uM DAPG 10ng/mL ATc RNA/DNA","15uM DAPG 40ng/mL ATc RNA/DNA")) + 
  labs(title="160913 Roman's 20 mean vs median of RNA/DNA ratio", y="Median", x = "Mean RNA/DNA ratio")
g.medvmean

mapGFP <- filter(left_join(bigDF,GFPBCs, by = "RevCompBC"),!is.na(plasmid))
mapGFP <- separate(mapGFP,plasmid, into = c("plasmid", "clone"),sep = "-")

GFPRNA <- filter(mapGFP, grepl("RNA", mapGFP$condition)) %>%
  select(plasmid, clone, condition, RevCompBC, "rcount" = nCount ) %>%
  separate(condition, into = c("template", "condition"))
GFPDNA <- filter(mapGFP, grepl("DNA", mapGFP$condition)) %>%
  select(plasmid, clone, condition, RevCompBC, "dcount" = nCount ) %>%
  separate(condition, into = c("template", "condition"))
GFPRD <- inner_join(GFPRNA,GFPDNA, by = c("RevCompBC","condition","plasmid","clone")) %>%
  mutate("RNAoDNA" = rcount/dcount*100) %>%
  # get rid of p63-4 which seems to be mismatched
  filter(RevCompBC != "GTTACCAGAATTCTGGGACC")
indStat <- group_by(GFPRD, condition) %>%
  summarise(meanRD = mean(RNAoDNA)) %>%
  mutate("norml" = max(meanRD)/meanRD)
nGFPRD <- full_join(GFPRD, indStat, by = "condition")
nGFPRD <- mutate(nGFPRD, "nRNAoDNA"=RNAoDNA*norml)

sumGFP <- group_by(nGFPRD, plasmid, clone) %>%
  summarise(mean_rcount =  mean(rcount), mean_dcount = mean(dcount), mean_RNADNA = mean(RNAoDNA), sd_RNADNA=sd(RNAoDNA), mean_nRD = mean(nRNAoDNA), sd_nRD = sd(nRNAoDNA))

sumGFP$clone <- factor(sumGFP$clone, levels = c("1"="1","2"="2","3"="3","4"="4","5"="5","6"="6","7"="7","8"="8","9"="9","10"="10","11"="11","12"="12","14"="14", "15"="15"))
gfpERR <- aes(ymax=mean_nRD + sd_nRD, ymin=mean_nRD - sd_nRD)
g.sumGFP <- ggplot(sumGFP, aes(clone, mean_nRD, fill= plasmid)) + geom_bar(stat="identity")+ facet_wrap(~plasmid, scales = "free") + geom_errorbar(gfpERR) +
  labs(title="160904 R20 GFP constructs Average RNA/DNA counts by clone after normalizing counts by condition", x ="clone", y = "Mean RNA/DNA") + scale_y_continuous(limits=c(0,390)) + theme(text = element_text(family="Calibri"), strip.background = element_rect(fill="white"))
g.sumGFP

g.sumGFPsd <- ggplot(sumGFP, aes(clone, sd_RNADNA-sd_nRD, fill= plasmid)) + geom_bar(stat="identity")+ facet_wrap(~plasmid, scales = "free") + 
  labs(title="160904 R20 GFP constructs by clone (across conditions) standard deviation - normalized standard deviation", x ="clone", y = "Mean RNA/DNA")  + theme(text = element_text(family="Calibri"), strip.background = element_rect(fill="white"))
g.sumGFPsd


