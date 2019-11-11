library(gplots)
library(RColorBrewer)
library(genefilter)
library(ggplot2)
library(gplots)
library(BioCircos)

## Read in the data (one file has ~5000 differentially expressed genes, other file has 250)
mostDiff <- read.csv(file="./MostDifferentiallyExpressed.csv", header=TRUE, sep=",")
allDiff <- read.csv(file="./AllDifferentiallyExpressed.csv", header=TRUE, sep=",")


myTrack=BioCircosTracklist() # Initialize Empty BioCircosTracklist
SNPTrack = BioCircosSNPTrack("mostSigGenes", chromosomes=as.character(mostDiff$chromosome_name), 
                            positions=mostDiff$start_position,
                            values = mostDiff$log2FoldChange,
                            colors = as.character(mostDiff$myColors),
                            minRadius=1.12,
                            maxRadius = 1.5)
myTrack=myTrack+SNPTrack
## Add one track for each chromosome
chromosomes=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15",
              "16","17","18","19","20", "21", "22", "X")

index=1
for (i in chromosomes){
  ## Define histogram/bars to be displayed
  nbBars = which(allDiff$chromosome_name==i)
  barValues = allDiff[nbBars,]$log2FoldChange
  barColor = colorRampPalette(brewer.pal(11, "Spectral"))(length(chromosomes))[index]
  ## Add a track with bars on the i-th chromosome
  startVec = allDiff[nbBars,]$start_position
  endVec = allDiff[nbBars,]$end_position
  myTrack = myTrack + BioCircosBarTrack(paste0("bars", i),
                                        chromosome=i,
                                        starts = startVec-2000000, #Make bar slightly thicker so it is visible
                                        ends = endVec+2000000,
                                        values = allDiff[nbBars,]$log2FoldChange,
                                        color = barColor, 
                                        minRadius = 0.0,
                                        maxRadius = 0.9)
  index=index+1
}
myTrack = myTrack + BioCircosBackgroundTrack("myBackground",
                                             minRadius = 0.0,
                                             maxRadius = 1,
                                             borderColors = "#000000",
                                             borderSize = 0.6,
                                             fillColors = "#000000")
myTrack = myTrack + BioCircosBackgroundTrack("myBackground",
                                             minRadius = 1.1,
                                             maxRadius = 2,
                                             borderColors = "#000000",
                                             borderSize = 0.6,
                                             fillColors = "#000000")
BioCircos(myTrack, genomeFillColor = "Spectral",
          chrPad = 0.00, displayGenomeBorder = FALSE, yChr =  FALSE, genomeLabelTextColor = "#FFFFFF",
          genomeTicksDisplay = FALSE,  genomeLabelTextSize = 22, genomeLabelDy = -0.3,  width='1200px', height='1200px')

dev.off()