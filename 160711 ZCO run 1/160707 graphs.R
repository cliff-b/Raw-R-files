library(ggplot2)
library(dplyr)
library(cowplot)
options(stringsAsFactors=F)

setwd('/Users/Kimberly/Documents/projects/forCliff')

df <- read.csv('seq_results.csv', header=T)

condition_names <- list('pdna_0a0d'='0A,0D pDNA', 'pdna_5d10a'='5D,10A pDNA', 'pdna_5d25A'='5D,25A pDNA', 'pdna_5d2_5a'='5D,2.5A pDNA', 'pdna_5d5a'='5D,5A pDNA', 'pdna_75d50a'='7.5D,50A pDNA')
condition_labeller <- function(variable, value){
    return(condition_names[value])
}

a <- ggplot(filter(df, grepl('pdna', condition)), aes(X, Y)) + geom_tile(aes(fill=log10(count))) + 
    facet_wrap(~ condition, labeller=condition_labeller, ncol=3, nrow=2) +
    theme(axis.text.x = element_text(angle=90, hjust=1, size=6), axis.text.y = element_text(size=6), strip.text.y=element_text(size=8)) +
    viridis::scale_fill_viridis(name='Normalized\nread count, log10')
a

condition_names <- list('cdna_0a0d'='0A,0D cDNA', 'cdna_5d10a'='5D,10A cDNA', 'cdna_5d25A'='5D,25A cDNA', 'cdna_5d2_5a'='5D,2.5A cDNA', 'cdna_5d5a'='5D,5A cDNA', 'cdna_75d50a'='7.5D,50A cDNA')

b <- ggplot(filter(df, grepl('cdna', condition)), aes(X, Y)) + geom_tile(aes(fill=log10(count))) + 
    facet_wrap(~ condition, labeller=condition_labeller, nrow=2, ncol=3) +
    theme(axis.text.x = element_text(angle=90, hjust=1, size=6), axis.text.y = element_text(size=6), strip.text.y=element_text(size=8)) +
    viridis::scale_fill_viridis(name='Normalized\nread count, log10')
b
