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

totalGenomeLength=3110748599
Chrom_1_len=248956422
Chrom_2_len=242193529
Chrom_3_len=198295559
Chrom_4_len=190214555
Chrom_5_len=181538259
Chrom_6_len=170805979
Chrom_7_len=159345973
Chrom_8_len=145138636
Chrom_9_len=138394717	
Chrom_10_len=133797422
Chrom_11_len=135086622
Chrom_12_len=133275309
Chrom_13_len=114364328
Chrom_14_len=107043718
Chrom_15_len=101991189
Chrom_16_len=90338345
Chrom_17_len=83257441
Chrom_18_len=80373285
Chrom_19_len=58617616
Chrom_20_len=64444167
Chrom_21_len=46709983
Chrom_22_len=50818468
Chrom_X_len=156040895
Chrom_Y_len=57227415
Chrom_MT_len=16569
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
noMT=data.frame(myDF[1:24, ])
noMT$listChrom=as.character(noMT$listChrom)
noMT$expected=as.character(noMT$expected)
noMT$numDiffEx=as.character(noMT$numDiffEx)
noMT$chiSqVals=as.character(noMT$chiSqVals)
noMT$pvals=as.character(noMT$pvals)
noMT$listChrom = factor(noMT$listChrom,
                                   levels=c('1','2','3','4','5','6',
                                            '7','8','9','10','11','12',
                                            '13','14','15','16','17','18',
                                            '19','20','21','22','X','Y','MT'))
ggplot(noMT, aes(x=factor(listChrom), y=as.numeric(negLog10P)))+
  geom_col(fill="blue")+
  xlab("Chromosome")+
  ggtitle("Chi-Square Test of Locations of Differentially Expressed Genes")+
  ylab(expression('-Log'[10]*'P-value'))+
  scale_fill_brewer(palette="Spectral")+
  theme_grey()+
  theme(plot.title = element_text(hjust = 0.5, size=13))
