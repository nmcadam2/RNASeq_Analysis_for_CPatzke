#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Written by Neil McAdams, Pfrender Lab, University of NotreDame - 2025
#Overviewed by Sheri Sanders, Department of Biological Sciences, University of NotreDame - 2025

#Libraries and Setup
#Libraries

library(AnnotationHub)
library(org.Mm.eg.db)
library(dplyr)
library(clusterProfiler)
library(data.table)
library(ggtangle)
library(enrichplot)

#Set Directories
dir <-"C:/Users/neilm/OneDrive/Desktop/Kabuki Syndrome RNASeq Analysis/"
setwd(dir)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#load data
data <-read.csv("Finished/Results/MouseHippocampus/DifferentialExpression/Results.csv", header = T)
gene.list<-data[,1]

#Convert gene symbol to entrez id
entrez.map<-mapIds(org.Mm.eg.db, keys = gene.list, column = "ENTREZID", keytype = "SYMBOL")
entrez.map

#Check: how many don't map
data$entrezid <- entrez.map[ data$X]
data <- data %>% filter(!is.na(entrezid))
sum(is.na(data$entrezid))

#Create gene list (lfc values and their respective entrez IDs sorted)
geneList = data[,3]
names(geneList) = as.character(data[,8])
geneList = sort(geneList, decreasing = T)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#gene set enrichment.
gseaDat <- gseGO(geneList     = geneList,
                 OrgDb        = org.Mm.eg.db,
                 ont          = "BP",
                 minGSSize    = 100,
                 maxGSSize    = 500,
                 pvalueCutoff = 0.05,
                 verbose      = FALSE,
                 nPermSimple = 10000)

#extract results and convert to gene symbols
resultTable<-setReadable(gseaDat, org.Mm.eg.db)
resultTable<-resultTable@result
write.csv(resultTable, "", quote = T)

#optional plotting of goHierarchy
#goplot(gseaDat)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Shortened Lists for Plotting BP Terms: These terms were selected following semantic similarity clustering to avoid plotting redundant terms.
#Select terms from hippocampus data 
mhPlotTerms<- c("immune effector process","neurotransmitter transport","regulation of immune effector process","adaptive immune response","cytokine-mediated signaling pathway","regulation of leukocyte activation","synaptic vesicle cycle","leukocyte cell-cell adhesion","embryonic epithelial tube formation","extrinsic apoptotic signaling pathway","exocytic process","neural tube development","stem cell differentiation","regulation of monoatomic ion transport","peptidyl-tyrosine modification","response to BMP","negative regulation of proteolysis","transforming growth factor beta receptor superfamily signaling pathway","cytoskeleton-dependent intracellular transport","epidermal cell differentiation","vesicle fusion","response to calcium ion","calcium ion homeostasis","regulation of innate immune response","negative regulation of immune system process","positive regulation of cytokine production","immune response-regulating signaling pathway","protein localization to synapse","regulation of gliogenesis","synapse assembly","collagen metabolic process","regulation of muscle contraction","autophagosome assembly","regulation of postsynaptic membrane neurotransmitter receptor levels","regulation of synapse structure or activity","homeostasis of number of cells","tumor necrosis factor superfamily cytokine production","tumor necrosis factor production","type I interferon production","interleukin-6 production","production of molecular mediator of immune response","cytokine production involved in immune response","regulation of postsynaptic membrane neurotransmitter receptor levels","regulation of synapse structure or activity","homeostasis of number of cells","regulation of membrane potential","regulation of postsynaptic membrane potential","action potential","regulation of type I interferon production","regulation of interleukin-6 production","regulation of production of molecular mediator of immune response","response to type II interferon","cellular response to type II interferon","immune response-regulating signaling pathway","canonical NF-kappaB signal transduction","cell surface receptor protein serine/threonine kinase signaling pathway","cell surface receptor signaling pathway via JAK-STAT","non-canonical NF-kappaB signal transduction","regulation of transmembrane receptor protein serine/threonine kinase signaling pathway","smoothened signaling pathway","immune effector process","signal release from synapse","neurotransmitter secretion","synaptic vesicle exocytosis","signal release","cellular response to lipid","regulation of canonical NF-kappaB signal transduction","regulation of synaptic plasticity","positive regulation of synaptic transmission","excitatory postsynaptic potential","chemical synaptic transmission, postsynaptic","synaptic transmission, glutamatergic","long-term synaptic potentiation","protein localization to synapse","protein localization to cell junction","receptor localization to synapse","protein-containing complex localization","regulation of gliogenesis","gliogenesis","neuron migration","neural tube development","synapse assembly","regulation of synapse organization","postsynapse organization","regulation of synapse assembly","regulation of postsynapse organization")

#Select terms from iN data
minPlotTerms<-c("microtubule bundle formation","cilium movement","cilium assembly","positive regulation of cell adhesion","positive regulation of cell motility","myeloid leukocyte differentiation","positive regulation of defense response","T cell activation","immune response-regulating signaling pathway","regulation of epithelial cell proliferation","angiogenesis","positive regulation of MAPK cascade","protein kinase B signaling","mesenchyme development","regulation of peptidyl-tyrosine phosphorylation","positive regulation of cytokine production","regulation of hormone levels","positive regulation of T cell activation","positive regulation of cytosolic calcium ion concentration","leukocyte proliferation","positive regulation of cell activation","regulation of metal ion transport","regulation of leukocyte proliferation","negative regulation of endopeptidase activity","tissue migration","immune effector process","regulation of fibroblast proliferation","regulation of body fluid levels","epithelial cell migration","positive regulation of secretion by cell","phagocytosis","adaptive immune response","regulation of cell-substrate adhesion","cellular response to peptide","regulation of reproductive process","lipid oxidation","fat cell differentiation","negative regulation of hydrolase activity","blood vessel diameter maintenance","regulation of DNA binding","leukocyte migration","determination of left/right symmetry","left/right pattern formation","leukocyte cell-cell adhesion","positive regulation of cell-cell adhesion","positive regulation of cell migration","T cell differentiation","inflammatory response","regulation of inflammatory response","positive regulation of immune response","transmembrane receptor protein tyrosine kinase signaling pathway","adenylate cyclase-modulating G protein-coupled receptor signaling pathway","skeletal system development","peptidyl-tyrosine phosphorylation","positive regulation of peptidyl-tyrosine phosphorylation","T cell proliferation","calcium ion transport","fibroblast proliferation","regulation of innate immune response","cell-matrix adhesion","extracellular matrix organization")

#plotting : treeplot
#Subset with terms of interest- Caution: These are specific to the terms used. The above terms can only be used for BP terms pulled from our Hippocampus data for example.
#Unfortunately there was an issue in GOsemsim for R4.3.3 and I had to update to R4.5.0 to perform the semantic similarity clustering used to cluster terms. Then downgrade again to finish the workflow in clusterProfiler. Otherwise i'd have used a single script.
selectDesc<- mhPlotTerms #mhPlotTerms or minPlotTerms
subsetGSEA <- filter(gseaDat, Description %in% selectDesc) #Guangchuang Yu's reccomended method for plotting select terms(

#For clustering of terms in the treepPlot
gseaDat2<-pairwise_termsim(subsetGSEA)

svg("",width = 16, height = 14,bg = "white")
treeplot(gseaDat2, hilight = F,cluster.params = list(color = c("black", "black", "black", "black", "black")))
dev.off()

#Plotting : cnet
gseaDatNet<-setReadable(gseaDat2, org.Mm.eg.db)
dataSig<- filter(data, padj < 0.05,abs(log2FoldChange) >= 0.5)
geneSig<- dataSig$X

#How many terms to include
showCat=5

#Build plot
p1<-enrichplot::cnetplot(gseaDatNet,
         shadowtext = 'gene',
         color.params = list(foldChange = geneList,color_category ='firebrick'),
         cex.params = list(gene_label = 0.6, category_label = 0.8),
         showCategory = showCat)

#Customize plot to only label genes that met our DE cutoffs
p1$layers
p1$layers[[4]]<-NULL
p2 <- p1 + ggraph::geom_node_text(data=p1$data[which( p1$data$name %in% c(p1$data$name[1:showCat],geneSig) ),],
                                  aes(x=x, y=y, label=name),
                                  color='black',
                                  bg.color="white",
                                  segment.size=0.5, repel=TRUE, size=5)

svg("./GSEA/CNETPlot_iNAll_MouseCC.svg",width = 16, height = 14,bg = "white")
p2
dev.off()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%