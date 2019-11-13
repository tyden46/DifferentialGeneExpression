library(gplots)
library(RColorBrewer)
library(genefilter)
library(ggplot2)
library(gplots)
library(BioCircos)
library(MASS)

## Read in the data (one file has ~5000 differentially expressed genes, other file has 250)
mostDiff <- read.csv(file="./MostDifferentiallyExpressed.csv", header=TRUE, sep=",")
allDiff <- read.csv(file="./AllDifferentiallyExpressed.csv", header=TRUE, sep=",")

Chrom_1=which(allDiff$chromosome_name=='1')
Chrom_2=which(allDiff$chromosome_name=='2')
Chrom_3=which(allDiff$chromosome_name=='3')
Chrom_4=which(allDiff$chromosome_name=='4')
Chrom_5=which(allDiff$chromosome_name=='5')
Chrom_6=which(allDiff$chromosome_name=='6')
Chrom_7=which(allDiff$chromosome_name=='7')
Chrom_8=which(allDiff$chromosome_name=='8')
Chrom_9=which(allDiff$chromosome_name=='9')
Chrom_10=which(allDiff$chromosome_name=='10')
Chrom_11=which(allDiff$chromosome_name=='11')
Chrom_12=which(allDiff$chromosome_name=='12')
Chrom_13=which(allDiff$chromosome_name=='13')
Chrom_14=which(allDiff$chromosome_name=='14')
Chrom_15=which(allDiff$chromosome_name=='15')
Chrom_16=which(allDiff$chromosome_name=='16')
Chrom_17=which(allDiff$chromosome_name=='17')
Chrom_18=which(allDiff$chromosome_name=='18')
Chrom_19=which(allDiff$chromosome_name=='19')
Chrom_20=which(allDiff$chromosome_name=='20')
Chrom_21=which(allDiff$chromosome_name=='21')
Chrom_22=which(allDiff$chromosome_name=='22')
Chrom_X=which(allDiff$chromosome_name=='X')
Chrom_Y=which(allDiff$chromosome_name=='Y')
Chrom_MT=which(allDiff$chromosome_name=='MT')

totalGenomeLength=20412
Chrom_1_len=2058
Chrom_2_len=1309
Chrom_3_len=1078
Chrom_4_len=752
Chrom_5_len=876
Chrom_6_len=1048
Chrom_7_len=989
Chrom_8_len=677
Chrom_9_len=786	
Chrom_10_len=733
Chrom_11_len=1298
Chrom_12_len=1034
Chrom_13_len=327
Chrom_14_len=830
Chrom_15_len=613
Chrom_16_len=873
Chrom_17_len=1197
Chrom_18_len=270
Chrom_19_len=1472
Chrom_20_len=544
Chrom_21_len=234
Chrom_22_len=488
Chrom_X_len=842
Chrom_Y_len=71
Chrom_MT_len=13
chromosomeGeneNums= c(Chrom_1_len, Chrom_1_len, Chrom_2_len, Chrom_2_len,
                      Chrom_3_len, Chrom_3_len, Chrom_4_len, Chrom_4_len,
                      Chrom_5_len, Chrom_5_len, Chrom_6_len, Chrom_6_len,
                      Chrom_7_len, Chrom_7_len, Chrom_8_len, Chrom_8_len,
                      Chrom_9_len, Chrom_9_len, Chrom_10_len, Chrom_10_len,
                      Chrom_11_len, Chrom_11_len, Chrom_12_len, Chrom_12_len,
                      Chrom_13_len, Chrom_13_len, Chrom_14_len, Chrom_14_len, Chrom_15_len,
                      Chrom_15_len, Chrom_16_len, Chrom_16_len, Chrom_17_len,
                      Chrom_17_len, Chrom_18_len, Chrom_18_len, Chrom_19_len,
                      Chrom_19_len, Chrom_20_len, Chrom_20_len, Chrom_21_len,
                      Chrom_21_len, Chrom_22_len, Chrom_22_len, Chrom_X_len,
                      Chrom_X_len, Chrom_Y_len, Chrom_Y_len, Chrom_MT_len,
                      Chrom_MT_len)
mySum = (Chrom_1_len +
           Chrom_2_len +
           Chrom_3_len +
           Chrom_4_len +
           Chrom_5_len +
           Chrom_6_len +
           Chrom_7_len +
           Chrom_8_len +
           Chrom_9_len +
           Chrom_10_len +
           Chrom_11_len +
           Chrom_12_len +
           Chrom_13_len +
           Chrom_14_len +
           Chrom_15_len +
           Chrom_16_len +
           Chrom_17_len +
           Chrom_18_len +
           Chrom_19_len +
           Chrom_20_len +
           Chrom_21_len +
           Chrom_22_len +
           Chrom_X_len +
           Chrom_Y_len +
           Chrom_MT_len)
numGenes=length(rownames(allDiff))      

ratios=c(Chrom_1_len/mySum,
         Chrom_2_len/mySum,
         Chrom_3_len/mySum,
         Chrom_4_len/mySum,
         Chrom_5_len/mySum,
         Chrom_6_len/mySum,
         Chrom_7_len/mySum,
         Chrom_8_len/mySum,
         Chrom_9_len/mySum,
         Chrom_10_len/mySum,
         Chrom_11_len/mySum,
         Chrom_12_len/mySum,
         Chrom_13_len/mySum,
         Chrom_14_len/mySum,
         Chrom_15_len/mySum,
         Chrom_16_len/mySum,
         Chrom_17_len/mySum,
         Chrom_18_len/mySum,
         Chrom_19_len/mySum,
         Chrom_20_len/mySum,
         Chrom_21_len/mySum,
         Chrom_22_len/mySum,
         Chrom_X_len/mySum,
         Chrom_Y_len/mySum,
         Chrom_MT_len/mySum)
numDiffEx = c(length(Chrom_1),
              length(Chrom_2),
              length(Chrom_3),
              length(Chrom_4),
              length(Chrom_5),
              length(Chrom_6),
              length(Chrom_7),
              length(Chrom_8),
              length(Chrom_9),
              length(Chrom_10),
              length(Chrom_11),
              length(Chrom_12),
              length(Chrom_13),
              length(Chrom_14),
              length(Chrom_15),
              length(Chrom_16),
              length(Chrom_17),
              length(Chrom_18),
              length(Chrom_19),
              length(Chrom_20),
              length(Chrom_21),
              length(Chrom_22),
              length(Chrom_X),
              length(Chrom_Y),
              length(Chrom_MT))
expected=numGenes*ratios

index=1
chiSqVals=c()
pvals=c()
chiSq=0
pval=0
for (i in expected){
  chiSqu=(((numDiffEx[index]-expected[index])^2)/(expected[index]))
  pval=pchisq(chiSqu, df=1, lower.tail=FALSE)
  chiSqVals=c(chiSqVals, chiSqu)
  pvals=c(pvals, pval)
  index=index+1
  chiSq=0
  pval=0
}
listChrom=c(1:22, 'X', 'Y', 'MT')
myDF=data.frame(listChrom, expected, numDiffEx, chiSqVals, pvals)

negLog10P = -log10(pvals)
myDF=data.frame(myDF, negLog10P)
myDF$listChrom=as.character(myDF$listChrom)
myDF$expected=as.character(myDF$expected)
myDF$numDiffEx=as.character(myDF$numDiffEx)
myDF$chiSqVals=as.character(myDF$chiSqVals)
myDF$pvals=as.character(myDF$pvals)
myDF$listChrom = factor(myDF$listChrom,
                                   levels=c('1','2','3','4','5','6',
                                            '7','8','9','10','11','12',
                                            '13','14','15','16','17','18',
                                            '19','20','21','22','X','Y','MT'))
ggplot(myDF, aes(x=factor(listChrom), y=as.numeric(negLog10P)))+
  geom_col(fill="blue")+
  xlab("Chromosome")+
  ggtitle("Chi-Square Test of Locations of Differentially Expressed Genes")+
  ylab(expression('-Log'[10]*'P-value'))+
  scale_fill_brewer(palette="Spectral")+
  theme_grey()+
  theme(plot.title = element_text(hjust = 0.5, size=13))


counts <- data.frame(myDF$listChrom, as.numeric(myDF$expected), as.numeric(myDF$numDiffEx))
colnames(counts) <- c("chromosome", "expected", "observed")
thisChrom <- as.vector(myDF$listChrom)[1]
expRow <- c(thisChrom, "expected", counts[1, "expected"])
obsRow <- c(thisChrom, "observed", counts[1, "observed"])
chartTable=rbind(expRow, obsRow)
for (i in 2:25){
  thisChrom <- as.vector(myDF$listChrom)[i]
  expRow <- c(thisChrom, "expected", counts[i, "expected"])
  obsRow <- c(thisChrom, "observed", counts[i, "observed"])
  chartTable=rbind(chartTable, expRow)
  chartTable=rbind(chartTable, obsRow)
}
colnames(chartTable) <- c("chromosome", "condition", "value")
chartTable=as.data.frame(chartTable)
chartTable$chromosome=as.character(chartTable$chromosome)
chartTable$condition=as.character(chartTable$condition)
chartTable$value=as.numeric(as.character(chartTable$value))
ggplot(chartTable, aes(fill=condition, y=value, x=chromosome)) + 
  geom_bar(position="dodge", stat="identity")
ggplot(chartTable, aes(fill=condition, y=value/chromosomeGeneNums, x=chromosome)) + 
  geom_bar(position="dodge", stat="identity")

myTrack=BioCircosTracklist() # Initialize Empty BioCircosTracklist

## Add one track for each chromosome
chromosomes=c("MT")

index=1
i="MT"
## Define histogram/bars to be displayed
nbBars = which(allDiff$chromosome_name==i & allDiff$log2FoldChange<0)
barValues = allDiff[nbBars,]$log2FoldChange
barColor = colorRampPalette(brewer.pal(11, "Spectral"))(14)
## Add a track with bars on the i-th chromosome
startVec = allDiff[nbBars,]$start_position
endVec = allDiff[nbBars,]$end_position
myTrack = myTrack + BioCircosBarTrack(paste0("bars", i),
                                      chromosome=i,
                                      starts = startVec,
                                      ends = endVec,
                                      values = -1* allDiff[nbBars,]$log2FoldChange,
                                      color = "#00FF00", 
                                      minRadius = 0.0,
                                      maxRadius = 0.9)
index=index+1

myGenome = list("MT"=16569)
title=capture.output(cat('Location of Downregulated Genes on the Mitochondrial Chromosome'))
myTrack = myTrack + BioCircosTextTrack("testText", title, weight = "bolder", 
                                       x = -1.05, y = -1.3, size="12px")
BioCircos(myTrack, genome=myGenome, genomeFillColor = "darkblue",
          chrPad = 0.1, displayGenomeBorder = FALSE, yChr =  FALSE, genomeLabelDisplay=TRUE,
          genomeTicksScale=1657, genomeTicksLen = 5, genomeLabelTextColor = "#FFFFFF",
          genomeTicksDisplay = TRUE,  genomeLabelTextSize = 22, width='500px', height='500px')


myTrack=BioCircosTracklist() # Initialize Empty BioCircosTracklist

## Add one track for each chromosome
chromosomes=c("MT")

index=1
i="MT"
## Define histogram/bars to be displayed
nbBars = which(allDiff$chromosome_name==i & allDiff$log2FoldChange>0)
barValues = allDiff[nbBars,]$log2FoldChange
barColor = colorRampPalette(brewer.pal(11, "Spectral"))(14)
## Add a track with bars on the i-th chromosome
startVec = allDiff[nbBars,]$start_position
endVec = allDiff[nbBars,]$end_position
myTrack = myTrack + BioCircosBarTrack(paste0("bars", i),
                                      chromosome=i,
                                      starts = startVec,
                                      ends = endVec,
                                      values = allDiff[nbBars,]$log2FoldChange,
                                      color = "red", 
                                      minRadius = 0.0,
                                      maxRadius = 0.9)
index=index+1

myGenome = list("MT"=16569)
title=capture.output(cat('Location of Upregulated Genes on the Mitochondrial Chromosome'))
myTrack = myTrack + BioCircosTextTrack("testText", title, weight = "bolder", 
                                       x = -1.05, y = -1.3, size="12px")
BioCircos(myTrack, genome=myGenome, genomeFillColor = "darkblue",
          chrPad = 0.1, displayGenomeBorder = FALSE, yChr =  FALSE, genomeLabelDisplay=TRUE,
          genomeTicksScale=1657, genomeTicksLen = 5, genomeLabelTextColor = "#FFFFFF",
          genomeTicksDisplay = TRUE,  genomeLabelTextSize = 22, width='500px', height='500px')


