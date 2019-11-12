library(gplots)
library(RColorBrewer)
library(genefilter)
library(ggplot2)
library(gplots)
library(BioCircos)

## Read in the data
genes <- read.csv(file="./AllDifferentiallyExpressed.csv", header=TRUE, sep=",")

pdf(file="LabelledDifferentiallyExpressed.pdf",height = 12,width = 12)

sigUp <- which(genes$log2FoldChange>0 & genes$padj<4*(10^-14))
sigDown <- which(genes$log2FoldChange<0 & genes$padj<4*(10^-14))
nonSigUp <- which(genes$log2FoldChange>0 & genes$padj>4*(10^-14))
nonSigDown <- which(genes$log2FoldChange<0 & genes$padj>4*(10^-14))
genes$Significance <- "FDR<4e-14, Upregulated"
genes[ sigUp, "Significance"] <- "FDR<4e-14, Upregulated"
genes[ sigDown, "Significance"] <- "FDR<4e-14, Downregulated"
genes[ nonSigUp, "Significance"] <- "Less Significant, Upregulated"
genes[ nonSigDown, "Significance"] <- "Less Significant, Downregulated"

ggplot(genes, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(aes(color = Significance)) +
  scale_color_manual(values = c("green1", "red", "#b8d9b4",  "#d9b4b4")) +
  xlab(expression('Log'[2]*'Fold Change in Expression Value')) +
  ylab(expression('-Log'[10]*'P-Value')) +
  ggtitle("Differentially Expressed Genes in SLE") +
  theme(text = element_text(size = 20),
        legend.position = "bottom", plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10)) +
  geom_text_repel(
    data = subset(genes, padj < 4*(10^-14)),
    aes(label = symbol),
    size = 4,
    box.padding = unit(0.4, "lines"),
    point.padding = unit(0.3, "lines")
  )
dev.off()