#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Written by Neil McAdams, Pfrender Lab, University of NotreDame - 2025
#Overviewed by Sheri Sanders, Department of Biological Sciences, University of NotreDame - 2025

#The freeCount shiny app created by Elizabeth Mae Brooks for GO enrichment is available at https://github.com/ElizabethBrooks/freeCount
# Load necessary library
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)

#Setwd#Sggplot2etwd
dir <-""
setwd(dir)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#read in data and fix truncated term strings
faOUT<-read.csv("freeCount_enrichmentResults.csv",header = T,sep = ",") # full results output from freeCount apps
panzOUT<-read.csv("Pannzer2Output.txt", header = T, sep = "\t", stringsAsFactors = FALSE) # Pannzer2 output file

#correct pannzer2 output to have the full id format
goMATCH<-panzOUT%>%
  mutate(GO.ID = sprintf("GO:%07d", as.integer(goid))) %>%
  select(GO.ID, desc) %>%
  distinct(GO.ID, .keep_all = T)

#correct truncated terms from the freeCount output
#if freeCount assigned a term ancestor not from the original pannzer2 file it will keep the truncated name to avoid NAs
faOUTcorrected<-faOUT%>%
  left_join(goMATCH, by = "GO.ID") %>%
  mutate(desc = if_else(is.na(desc), Term, desc)) %>%
  select(-Term)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Join with DE results
#Load in name2Go
n2g<-read.table("GeneName2GO.txt") #Cutsom input for reformatting
colnames(n2g)<-c("gene", "goterms")

#convert to long format
n2gLong <- n2g %>%
  separate_rows(goterms, sep = ",") %>%
  rename(GO.ID = goterms)

#join with enriched Go terms
Results <- n2gLong %>%
  inner_join(faOUTcorrected, by = "GO.ID")

#read in DE Results
DEG<-read.csv("DEresults.csv", header = T); head(DEG)
colnames(DEG)<-c("gene","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"); head(DEG)

#join by gene name
Final<-Results %>%
    inner_join(DEG, by ="gene")

write.csv(Final, "XX_GO_terms_Results.csv",row.names = F) #where XX is MF/BP/CC
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Filtering and Plotting
#filter for only DEG but all go terms regardless of significant enrichment
allterms.DEG<-Final[Final$padj<=0.05,]

#filter for only DEG and only significantly enriched GO terms
sigTerms.DEG<-allterms.DEG[allterms.DEG$weightFisher<=0.05,]
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#just sig genes
#reformat for plotting - this will count how many DE genes are within each go term
PlotDat_allterms.DEG<-allterms.DEG %>%
  group_by(GO.ID, desc, weightFisher) %>%
  summarise(gene_count = n_distinct(gene), .groups = 'drop') %>%
  arrange(desc(gene_count)) %>%
  slice_head(n=20)

#plot significantly DE genes but with no filter that they must be within significantly enriched GO terms
svg("XX_sigGenes.svg",width = 8, height = 7,bg = "white") 
ggplot(PlotDat_allterms.DEG, aes(x = gene_count, y = reorder(desc, gene_count), fill = weightFisher)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "powderblue", high = "darkblue", name = "weightFisher") +
  labs(
    title = "Top 20 XX GO Terms by Number of DE Genes",
    x = "Number of Differentially Expressed Genes",
    y = "GO Term"
  ) +
  theme_minimal()
dev.off()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#sig genes and sig enrichment
#reformat for plotting
PlotDat_sigterms.DEG<-sigTerms.DEG %>%
  group_by(GO.ID, desc, weightFisher) %>%
  summarise(gene_count = n_distinct(gene), .groups = 'drop') %>%
  arrange(desc(gene_count)) %>%
  slice_head(n=20)

#plot significantly enriched and significantly de genes
svg("XX_sigEnrichedsigGenes.svg",width = 8, height = 7,bg = "white") 
ggplot(PlotDat_sigterms.DEG, aes(x = gene_count, y = reorder(desc, gene_count), fill = weightFisher)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "powderblue", high = "darkblue", name = "weightFisher") +
  labs(
    title = "Top 20 Significantly Enriched XX GO Terms by Number of DE Genes",
    x = "Number of Differentially Expressed Genes",
    y = "GO Term"
  ) +
  theme_minimal()
dev.off()
