BiocManager::install("clusterProfiler")
BiocManager::install("pathview")
BiocManager::install("enrichplot")
library(clusterProfiler)
library(enrichplot)
# we use ggplot2 to add x axis labels (ex: ridgeplot)
library(ggplot2)

# SET THE DESIRED ORGANISM HERE
organism = "org.Mm.eg.db"
#BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)
############

library(WGCNA)
library(DESeq2)
library(GEOquery)
library(tidyverse)
library(CorLevelPlot)
library(gridExtra)
library(EnhancedVolcano)
library(IHW)
library(ggplot2)
allowWGCNAThreads() 

data <-  read.table(file = "STAR_gene_count.tsv", header = TRUE, row.names = 1)
phenoData <- read.table(file = "sample_sheet.txt", header = TRUE, row.names = 1)
phenoData$MucusGrowth <- as.factor(phenoData$MucusGrowth)
phenoData$Group <- as.factor(phenoData$Group)
# 2. QC - outlier detection ------------------------------------------------
# detect outlier genes

gsg <- goodSamplesGenes(t(data))
summary(gsg)
gsg$allOK

table(gsg$goodGenes)
table(gsg$goodSamples)

# remove genes that are detectd as outliers
data <- data[gsg$goodGenes == TRUE,]

# detect outlier samples - hierarchical clustering - method 1
htree <- hclust(dist(t(data)), method = "average")
plot(htree)

# pca - method 2

pca <- prcomp(t(data))
pca.dat <- pca$x

pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)

pca.dat <- as.data.frame(pca.dat)

ggplot(pca.dat, aes(PC1, PC2)) +
  geom_point() +
  geom_text(label = rownames(pca.dat)) +
  labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
       y = paste0('PC2: ', pca.var.percent[2], ' %'))

### NOTE: If there are batch effects observed, correct for them before moving ahead


# exclude outlier samples
#samples.to.be.excluded <- c('H01', 'E01')
samples.to.be.excluded <- c()
data.subset <- data[,!(colnames(data) %in% samples.to.be.excluded)]

# 3. Normalization ----------------------------------------------------------------------
# create a deseq2 dataset

# exclude outlier samples
colData <- phenoData %>% 
  filter(!row.names(.) %in% samples.to.be.excluded)


# fixing column names in colData
names(colData)
names(colData) <- gsub(':ch1', '', names(colData))
names(colData) <- gsub('\\s', '_', names(colData))
names(data.subset) <- gsub("^X", "", names(data.subset))
names(data.subset) <- gsub("\\.", "-", names(data.subset))
# making the rownames and column names identical
all(rownames(colData) %in% colnames(data.subset))
all(rownames(colData) == colnames(data.subset))

# create dds
dds <- DESeqDataSetFromMatrix(countData = data.subset,
                              colData = colData,
                              design = ~ Group) # not spcifying model

## remove all genes with counts < 15 in more than 75% of samples (31*0.75=23.25)
## suggested by WGCNA on RNAseq FAQ

dds75 <- dds[rowSums(counts(dds) >= 30) >= 4,]
#dds75 <- dds
nrow(dds75) # 15405 genes

# perform variance stabilization
dds_norm <- vst(dds75)

# get normalized counts
#norm.counts <- assay(dds_norm) %>% t()
dds_out <- DESeq(dds)
res <- results(dds_out, contrast=c("Group","BL_WSD","HF_Chow"))
#res <- results(dds_out)
#write.csv(res, file = "results_Group.csv", quote = FALSE)
write.csv(res, file = "results_Group.csv", quote = FALSE)

df = read.csv("lightgreen_results_Group_V1.csv", header=TRUE)
###################################
# Plot selected genes
###################################
library("ggbeeswarm")
geneCounts <- plotCounts(dds, gene = "ENSMUSG00000046259", intgroup = "Group",
                         returnData = TRUE)
ggplot(geneCounts, aes(x = Group, y = count, color = Group)) +
  scale_y_log10() +  geom_beeswarm(cex = 5)+ theme_bw()

####################################
#Volcano Plots
####################################
jpeg("VolcanoPlot_HF_Chow_WSD.jpg", height = 20, width = 20, units = 'in', res = 600)
EnhancedVolcano(df,
                lab = rownames(df),
                x = 'log2FoldChange',
                y = 'pvalue',
                xlab = bquote(~Log[2]~ 'fold change'),
                pointSize = 4.0,
                labSize = 6.0,
                colAlpha = 1,
                legendPosition = 'right',
                pCutoff = 0.05,
                legendLabSize = 12,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.75)
dev.off()

pdf(file="VolcanoPlot_BL_WSD.pdf", height = 20, width = 20)

EnhancedVolcano(df,
                lab = rownames(df),
                x = 'log2FoldChange',
                y = 'pvalue',
                xlab = bquote(~Log[2]~ 'fold change'),
                pointSize = 4.0,
                labSize = 6.0,
                colAlpha = 1,
                legendPosition = 'right',
                pCutoff = 0.05,
                legendLabSize = 12,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.75)
dev.off()

##########################################
#CLUSTER PROFILER
#########################################
df = read.csv("lightgreen_results_Group_V1.csv", header=TRUE)

# we want the log2 fold change 
original_gene_list <- df$log2FoldChange

# name the vector
names(original_gene_list) <- df$X

# omit any NA values 
gene_list<-na.omit(original_gene_list)


# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)
set.seed(34)  
gse <- gseGO(geneList=gene_list, 
             ont ="MF", 
             keyType = "ENSEMBL",
             minGSSize = 3, 
             maxGSSize = 50, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "BH")

require(DOSE)
pdf(file="LightCyan_GO_CC_dotplot.pdf", height = 20, width = 15)
dotplot(gse, showCategory=30, split=".sign") + facet_grid(.~.sign)
dev.off()

##################################################
#GSE-Plots
#################################################
values <- 2:length(gse$Description)-1# Your range of values
max_val <- length(gse$Description)-1
plot_list <- list()

for (i in values) {
  # Create the command as a string
  command <- paste("gseaplot(gse, by = 'all', title = gse$Description[", i, "], geneSetID =", i, ")")
  # Evaluate and execute the command
  p <- eval(parse(text = command))
  plot_list[[i]] <- p
}

pdf("GSE_GO_BP_plots.pdf")

# Print each plot to the PDF device
for (i in 1:max_val) {
  print(plot_list[[i]])
}

dev.off()
##################################################
# KEGG Analysis
##################################################
# Convert gene IDs for gseKEGG function
# We will lose some genes here because not all IDs will be converted
ids<-bitr(names(original_gene_list), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb=organism)
# remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
dedup_ids = ids[!duplicated(ids[c("ENSEMBL")]),]

# Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
df2 = df[df$X %in% dedup_ids$ENSEMBL,]

# Create a new column in df2 with the corresponding ENTREZ IDs
df2$Y = dedup_ids$ENTREZID

# Create a vector of the gene unuiverse
kegg_gene_list <- df2$log2FoldChange

# Name vector with ENTREZ ids
names(kegg_gene_list) <- df2$Y

# omit any NA values 
kegg_gene_list<-na.omit(kegg_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)
set.seed(34)
kegg_organism <- "mmu"
kk2 <- gseKEGG(geneList     = kegg_gene_list,
               organism     = kegg_organism,
               minGSSize    = 3,
               maxGSSize    = 750,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType       = "ncbi-geneid")

dotplot(kk2, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
# categorySize can be either 'pvalue' or 'geneNum'
cnetplot(kk2, categorySize="pvalue", foldChange=gene_list)

ridgeplot(kk2) + labs(x = "enrichment distribution")
