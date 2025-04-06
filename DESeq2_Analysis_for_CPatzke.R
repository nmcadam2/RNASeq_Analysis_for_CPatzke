#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Written by Neil McAdams, Pfrender Lab, University of NotreDame - 2025
#Overviewed by Sheri Sanders, Department of Biological Sciences, University of NotreDame - 2025

#For general DESeq2 troubleshooting see the DESeq2 vignette: https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

#Libraries and Setup
#Libraries
library(limma)
library(edgeR)
library(readr)
library(variancePartition)
library(DESeq2)
library(pheatmap)
library(tidyverse)
library(dplyr)
library(RColorBrewer)

#Set Directories
res.dir <-"/Path/to/resultsDirectory"
dir <-"/Path/to/workingDirectory"
setwd(dir)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Load and Format Data
#Load in data
metaTable <- read.csv("designTable.csv"); head(metaTable)
countsTab <- read.csv("countsTable.csv", header = TRUE, row.names = 1); head(countsTab)

#Make counts integers
countsTab[] <- lapply(countsTab, function(x) if (is.numeric(x)) as.integer(x) else x); head(countsTab)

#Relevel 
metaTable$treatment<-relevel(as.factor(metaTable$treatment), "deltaCre")

#Remove Y chromosome genes - Experiment did not control for mouse sex so Y must be removed from analysis
#Literature reports few differences in mouse gene expression between chromosomal sexes at early developmental stages.
Ygenes <- readLines("YChromosome_genes.txt")
nrow(countsTab);countsTab <- countsTab %>% filter(!(rownames(countsTab) %in% Ygenes)); nrow(countsTab) # number should drop by ~340 in mouse samples

#Filter low counts - Using 15 for a slightly more stringent than the default reccomendation of DESeq2
drop = rowSums(countsTab < 15) > 6
workingCounts = countsTab[!drop, ]

#make metaData factors - skip if using continuous lib size covariate
for (i in 1:ncol(metaTable)) {metaTable[,i] <- as.factor(metaTable[,i])}

#write model - If you did the optional step of libsize as a covariate use ~libsize + group + treatment
design = ~group + treatment

#Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = workingCounts,
                              colData = metaTable,
                              design = design)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Run analysis on Object
#Collapse batches
ddsColl<-collapseReplicates(dds, dds$collapseBy);head(ddsColl@assays@data[["counts"]])

#Run DESeq2 model
dds = DESeq(ddsColl)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Exploratory visualization - technically optional but good practice to sanity check
#stablize variance
vsd = vst(dds, blind=FALSE)

#Dispersion
svg(paste(res.dir,"Dispersion.svg", sep = ""),width = 8, height = 7,bg = "white")
plotDispEsts(dds)
dev.off()

#Distance matrix
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
distanceDat<-log2(assay(vsd))
colnames(distanceDat)<-substring(vsd$sample,1,nchar(as.character(vsd$sample))-1)

sample_metadata <- data.frame(treatment = c("Cre", "Cre", "Cre","deltaCre", "deltaCre", "deltaCre"),
                              row.names = substring(vsd$sample,1,nchar(as.character(vsd$sample))-1))

ann_colors<-list(treatment = c(Cre="#FFD500", deltaCre="#005BBB"))

svg(paste(res.dir,"DistanceMap.svg", sep = ""),width = 8, height = 7, bg = "white") 
pheatmap(cor(distanceDat),
         clustering_distance_rows="euclidean",
         clustering_distance_cols="euclidean",
         cellwidth = 55,
         cellheight = 55,
         annotation_col = sample_metadata,
         annotation_row = sample_metadata,
         annotation_names_col = F,
         annotation_names_row = F,
         annotation_colors = ann_colors,
         labels_row = colnames(distanceDat),
         labels_col = c("","","","","",""),
         col = colors)
dev.off()


#PCA
pcaData <- plotPCA(vsd, intgroup=c("group", "treatment"), returnData=TRUE)
colnames(pcaData)<- c("PC1", "PC2","group.1","group","treatment","name")
percentVar <- round(100 * attr(pcaData, "percentVar"))

svg(paste(res.dir,"PCA.svg", sep = ""),width = 8, height = 7, bg = "white") 
ggplot(pcaData, aes(PC1, PC2, color=treatment, shape=group)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed()+
  theme_bw()+
  scale_color_manual(name = c("Cre", "deltaCre"), values=c("#FFD500","#005BBB"))
dev.off()

#MA - Two options here: DESeq2 Encourages the use shrinkage estimators like apeglm when plotting MA to avoid the use of arbitrary filtering thresholds
#For human / mouse co-cultured samples fold changes were so modest that shrinkage is inappropriate, as such only raw should be used with the co-cultured samples.
# source :https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#ma-plot
#Raw
svg(paste(res.dir,"MA_raw.svg", sep = ""),width = 8, height = 7,bg = "white") 
plotMA(dds, ylim=c(-5,5))
dev.off()

#Shrunken - Coef name shouldn't need changed but if you encounter an error check valid coef names with resultsNames(dds)
resLFC <- lfcShrink(dds, coef="treatment_Cre_vs_deltaCre", type="apeglm") # takes a moment!
svg(paste(res.dir,"MA_shrunken.svg", sep = ""),width = 8, height = 7,bg = "white") 
plotMA(resLFC, ylim=c(-2,2))
dev.off()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Results and result plotting
#Generate Results
results <- results(dds, contrast = c("treatment", "Cre", "deltaCre"))

#Filter Results - I opt for p-adj <= 0.05 and LFC of +/-0.5 per literature values of LFC in neuron cell culture
res.nona <- na.omit(results)
res.sig = res.nona[res.nona$padj <= 0.05,] ; nrow(res.sig) # for number of DE Genes at p <= 0.05
res.sig.lfc = res.sig[abs(res.sig$log2FoldChange) > 0.5,]; nrow(res.sig.lfc) # for number of DE Genes at p <= 0.05 & LFC =/- 0.5

write.csv(res.nona, paste(res.dir,"Results.csv", sep = ""), quote=FALSE)
write.csv(res.sig, paste(res.dir,"sigResults.csv", sep = ""), quote=FALSE)
write.csv(res.sig.lfc, paste(res.dir,"sig_fcResults.csv", sep = ""), quote=FALSE)

#Visualizations
#Volcano plot
source("old/fancy_volcano_plot.R") #Custom Volcano plot script prepared by the Patzke lab - I modified it slightly to allow for custom x axis limits and custom titles

svg(paste(res.dir,"Volcano.svg", sep = ""), width = 8, height = 7,bg = "white") 
fancy_volcano_plot(res.nona, 0.05, 0, 3, "Title")
dev.off()

#Heatmap
#data prep
sig_genes <- row.names(res.sig)
log_sig_counts <- log2(assay(vsd)[sig_genes, ])
colnames(log_sig_counts)<-substring(vsd$sample,1,nchar(as.character(vsd$sample))-1)

svg(paste(res.dir,"DE_Heatmap.svg", sep = ""), width = 8, height = 7, bg = "white") 
pheatmap(na.omit(log_sig_counts),
         cluster_rows=TRUE, cluster_cols=TRUE,   # Clustering
         cutree_rows = 2,
         annotation_col = sample_metadata,
         annotation_names_col = F,
         annotation_colors = ann_colors,
         scale="row",     # Scale by row (gene) to highlight relative expression differences
         show_rownames=F, # Optional: hide gene names for large datasets
         show_colnames = TRUE,
         main="Heatmap of Significant DE Genes")
dev.off()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%