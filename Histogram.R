library(gplots)
library(RColorBrewer)
library(genefilter)
library(ggplot2)
library(gplots)
library(BioCircos)

## Read in the data (one file has ~5000 differentially expressed genes, other file has 250)
mostDiff <- read.csv(file="./MostDifferentiallyExpressed.csv", header=TRUE, sep=",")
allDiff <- read.csv(file="./AllDifferentiallyExpressed.csv", header=TRUE, sep=",")

sortedchrome=order(nchar(as.character(allDiff$chromosome_name)),
                   allDiff$chromosome_name) 

allDiff$chromosome_name_f = factor(allDiff$chromosome_name,
                                   levels=c('1','2','3','4','5','6',
                                            '7','8','9','10','11','12',
                                            '13','14','15','16','17','18',
                                            '19','20','21','22','X','Y','MT'))
png(filename="histPlot1.png", width=2000, height=1000)

myColorPal = colorRampPalette(brewer.pal(11, "Spectral"))(length(chromosomes))
ggplot(allDiff) + geom_histogram(aes(x=log2FoldChange))+
  aes(fill= as.factor(chromosome_name_f))+
  facet_wrap(~chromosome_name_f) +
  theme_grey() +
  labs(fill = "Chromosome") +
  ggtitle("Differential Gene Expression By Chromosome") +
  theme(plot.title = element_text(hjust = 0.5, size=30))+
  scale_colour_brewer(palette = colorRampPalette(brewer.pal(3, "Spectral"))(5))

dev.off()
