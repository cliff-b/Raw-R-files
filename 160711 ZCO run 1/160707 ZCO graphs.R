library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(viridis)
library(extrafont)
library(jpeg)
library(grid)

setwd("/Users/Cliff/Documents/Kosuri Lab/NGS files/160711 ZCO run 1")
df <- read.csv("160706 seq_results.csv",stringsAsFactors = FALSE)
img <- readJPEG("160721 background.jpg")
img <- matrix(rgb(img[,,1],img[,,2],img[,,3], 0.2), nrow =dim(img)[1])  
#Ghetto-ass labeller for ggplots
condition_names <- list("dna_0a0d"="0A,0D RNA/DNA", "dna_5d2_5a"="5D,2.5A RNA/DNA", "dna_5d5a"="5D,5A RNA/DNA", "dna_5d10a"="5D,10A RNA/DNA", "dna_5d25A"="5D,25A RNA/DNA", "dna_75d50a"="7.5D,50A RNA/DNA",
                        "cdna_0a0d"="0A,0D cDNA", "cdna_5d2_5a"="5D,2.5A cDNA", "cdna_5d5a"="5D,5A cDNA", "cdna_5d10a"="5D,10A cDNA", "cdna_5d25A"="5D,25A cDNA", "cdna_75d50a"="7.5D,50A cDNA",
                        "pdna_0a0d"="0A,0D pDNA", "pdna_5d2_5a"="5D,2.5A pDNA", "pdna_5d5a"="5D,5A pDNA", "pdna_5d10a"="5D,10A pDNA", "pdna_5d25A"="5D,25A pDNA", "pdna_75d50a"="7.5D,50A pDNA")
condition_labeller <- function(variable, value){
  variable <- as.character(variable)
  return(condition_names[value])
}

#split out the plasmid DNA from the dataframe
pDNAdf <- filter(df, grepl("pdna", condition), grepl("X", X) & grepl("Y",Y))

g.pDNA <- ggplot(pDNAdf, aes(X, Y)) + geom_tile(aes(fill=log10(count))) + facet_wrap(~condition, labeller=condition_labeller, ncol=3, nrow=2) +
  theme(text=element_text(family="Comic Sans MS"),  strip.background=element_rect(fill = "White"), axis.text.x = element_text(angle=90, hjust=1, size=6), axis.text.y = element_text(size=6), strip.text.y=element_text(size=8), legend.title=element_text(size = 8), legend.text=element_text(size =8)) +
  scale_fill_viridis(option = "inferno", name="Plasmid DNA \ncounts(Log10)") + xlab("X ZCO Construct") + ylab("Y ZCO Construct") + ggtitle("160711 DNA counts per barcode")
g.pDNA

#Work with the average of the four barcodes not just one. Should really have some stats on that. 
pDNAnorm1 <- arrange(pDNAdf, condition, X, Y) %>%
  separate(X, into=c("X", "BC1"), sep ="-") %>%
  separate(Y, into=c("Y", "BC2"), sep = "-")%>%
  group_by(condition, Y, X) %>%
  summarise(mean_x = mean(count), std_x = sd(count))
  
pDNAnorm1$condition = factor(pDNAnorm1$condition, levels=c("pdna_0a0d", "pdna_5d2_5a", "pdna_5d5a", "pdna_5d10a", "pdna_5d25A", "pdna_75d50a"))
g.pDNAnorm <- ggplot(pDNAnorm1, aes(X, Y)) + geom_tile(aes(fill=log10(mean_x))) + facet_wrap(~condition, labeller=condition_labeller, ncol=3, nrow=2) +
  theme(text=element_text(family="Comic Sans MS"),  strip.background=element_rect(fill = "White"), axis.text.x = element_text(angle=90, hjust=1, size=6), axis.text.y = element_text(size=6), strip.text.y=element_text(size=8), legend.title=element_text(size = 8), legend.text=element_text(size =8)) +
  scale_fill_viridis(option = "inferno", name="Plasmid DNA \ncount (log10)") + xlab("X ZCO Construct") + ylab("Y ZCO Construct")
g.pDNAnorm

g.pDNAstdev <- ggplot(pDNAnorm1, aes(X, Y)) + geom_tile(aes(fill=log10(std_x))) + facet_wrap(~condition, labeller=condition_labeller, ncol=3, nrow=2) +
  theme(text=element_text(family="Comic Sans MS"),  strip.background=element_rect(fill = "White"), axis.text.x = element_text(angle=90, hjust=1, size=6), axis.text.y = element_text(size=6), strip.text.y=element_text(size=8), legend.title=element_text(size = 8), legend.text=element_text(size =8), plot.title=element_text(color="red")) +
  scale_fill_viridis(option = "inferno", name="Stdev of\nPlasmid DNA \nBarcode \nCounts (log10)") + labs(title="Standard Deviation of Barcode Counts") + xlab("X ZCO Construct") + ylab("Y ZCO Construct")
g.pDNAstdev


#now do the same things with the cDNA
cDNAdf <- filter(df, grepl("cdna", condition), grepl("X",X) & grepl("Y",Y))

g.cDNA <- ggplot(cDNAdf, aes(X,Y)) + geom_tile(aes(fill=log10(count))) + facet_wrap(~condition, labeller = condition_labeller, ncol=3, nrow=2) +
  theme(text=element_text(family="Comic Sans MS"),  strip.background=element_rect(fill = "White"), axis.text.x = element_text(angle=90, hjust=1,size =6), axis.text.y = element_text(size = 6), strip.text.y=element_text(size = 8), legend.title=element_text(size = 8), legend.text=element_text(size =8)) +
  scale_fill_viridis(option = "inferno", name="RT'd RNA \ncounts (Log10)") + xlab("X ZCO Construct") + ylab("Y ZCO Construct")
g.cDNA

cDNAnorm1 <- arrange(cDNAdf, condition, X, Y) %>%
  separate(X, into=c("X", "BC1"), sep ="-") %>%
  separate(Y, into=c("Y", "BC2"), sep = "-")%>%
  group_by(condition, Y, X) %>%
  summarise(mean_x = mean(count), std_x = sd(count))

cDNAnorm1$condition = factor(cDNAnorm1$condition, levels=c("cdna_0a0d", "cdna_5d2_5a", "cdna_5d5a", "cdna_5d10a", "cdna_5d25A", "cdna_75d50a"))
g.cDNAnorm <- ggplot(cDNAnorm1, aes(X,Y)) + geom_tile(aes(fill=log10(mean_x))) + facet_wrap(~condition, labeller = condition_labeller, ncol=3, nrow=2) +
  theme(text=element_text(family="Comic Sans MS"),  strip.background=element_rect(fill = "White"), axis.text.x = element_text(angle=90, hjust=1,size =6), axis.text.y = element_text(size = 6), strip.text.y=element_text(size = 8), legend.title=element_text(size = 8), legend.text=element_text(size =8)) +
  scale_fill_viridis(option = "inferno", name="RT'd RNA \ncounts (log10)") + xlab("X ZCO Construct") + ylab("Y ZCO Construct") + ggtitle("ZCO cDNA averaged barcode counts")
g.cDNAnorm

g.cDNAstdev <- ggplot(cDNAnorm1, aes(X,Y)) + geom_tile(aes(fill=log10(std_x))) + facet_wrap(~condition, labeller = condition_labeller, ncol=3, nrow=2) +
  theme(text=element_text(family="Comic Sans MS"),  strip.background=element_rect(fill = "White"), axis.text.x = element_text(angle=90, hjust=1,size =6), axis.text.y = element_text(size = 6), strip.text.y=element_text(size = 8), legend.title=element_text(size = 8), legend.text=element_text(size =8)) +
  scale_fill_viridis(option = "inferno", name="Stdev of \nRT'd RNA \ncounts (log10)") + xlab("X ZCO Construct") + ylab("Y ZCO Construct") + ggtitle("ZCO cDNA standard deviation of barcode counts")
g.cDNAstdev

colnames(cDNAdf)[4] <- "count1"
cDNAdf$condition <- gsub("c","",cDNAdf$condition)
pDNAdf$condition <- gsub("p","",pDNAdf$condition)
bothDF <- left_join(pDNAdf,cDNAdf)
bothDF <- mutate(bothDF,normal = count1/count)

bothDFavg <- arrange(bothDF) %>%
  separate(X, into=c("X","BC1"), sep = "-") %>%
  separate(Y, into=c("Y","BC2"), sep = "-") %>%
  group_by(condition, X, Y) %>%
  summarise(mean_x=mean(normal, na.rm=TRUE), std_x=sd(normal, na.rm=TRUE))

bothDFavg$condition <- factor(bothDFavg$condition, levels = c("dna_0a0d","dna_5d2_5a","dna_5d5a","dna_5d10a","dna_5d25A","dna_75d50a"))

g.AVGnormal <- ggplot(bothDFavg, aes(X,Y)) + geom_tile(aes(fill=log10(mean_x))) + facet_wrap(~condition, labeller = condition_labeller, ncol=3, nrow=2, scales = "free") +
  theme(text=element_text(family="Comic Sans MS"), strip.text=element_text(size = 10), strip.background=element_rect(fill = "White"), axis.text.x = element_text(angle=90, hjust=1,size =6), axis.text.y = element_text(size = 6), legend.title=element_text(size = 8),
        legend.title.align=0.5, legend.text=element_text(size =8)) + 
  scale_fill_viridis(option="inferno", name = "Normalized \nRNA to DNA \n(Log10)", na.value="white") + xlab("X ZCO Construct") + ylab("Y ZCO Construct") + ggtitle("ZCO Averaged RNA/DNA Barcode Counts Ratio")
g.AVGnormal

g.STDnormal <- ggplot(bothDFavg, aes(X,Y)) + geom_tile(aes(fill=log10(std_x))) + facet_wrap(~condition, labeller = condition_labeller, ncol=3, nrow=2, scales = "free") +
  theme(text=element_text(family="Comic Sans MS"), strip.text=element_text(size = 10), strip.background=element_rect(fill = "White"), axis.text.x = element_text(angle=90, hjust=1,size =6), axis.text.y = element_text(size = 6), legend.title=element_text(size = 8),
        legend.title.align=0.5, legend.text=element_text(size =8)) + #annotation_custom(rasterGrob(background, width = unit(1,"npc"), height = unit(1,"npc")), -Inf, Inf, -Inf, Inf) +
  scale_fill_viridis(option="inferno", name = "Stdev of \nRNA to DNA \n(Log10)", na.value="white") + xlab("X ZCO Construct") + ylab("Y ZCO Construct") + ggtitle("ZCO Standard Deviation of RNA/DNA Barcode Counts Ratio")
g.STDnormal

g.ratio <- ggplot(bothDFavg, aes(mean_x, std_x)) + geom_point(alpha = 0.33) + labs(title = "160711 ZCO mean vs standard deviation", x= "mean RNA/DNA", y = "standard deviation of RNA/DNA") + 
  geom_smooth(method = "lm", se = FALSE) + annotation_custom(rasterGrob(img, width = unit(1,"npc"), height = unit(1,"npc")), -Inf, Inf, -Inf, Inf) + 
  theme(text=element_text(family="Comic Sans MS"), plot.title =  element_text(family = "Jokerman"))
g.ratio

g.normalize <- ggplot(bothDF, aes(X,Y)) + geom_tile(aes(fill=log10(normal)))  + facet_wrap(~condition, labeller = condition_labeller, ncol=3, nrow=2) +
  theme(text=element_text(family="Comic Sans MS"), strip.text=element_text(size = 10), strip.background=element_rect(fill = "White"), axis.text.x = element_text(angle=90, hjust=1,size =6), axis.text.y = element_text(size = 6), legend.title=element_text(size = 8),
        legend.title.align=0.5, legend.text=element_text(size =8)) + 
  scale_fill_viridis(option="inferno", name = "Normalized \nRNA to DNA \n(Log10)") + xlab("X ZCO Construct") + ylab("Y ZCO Construct")
g.normalize
