########################################################################################
#Over-Representation Analysis with ClusterProfiler
#link:https://learn.gencore.bio.nyu.edu/rna-seq-analysis/over-representation-analysis/
########################################################################################
BiocManager::install("clusterProfiler", version = "3.8")
BiocManager::install("pathview")
install.packages("wordcloud")
library(clusterProfiler)
library(wordcloud)

# SET THE DESIRED ORGANISM HERE
organism = "org.Mm.eg.db"
#BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

df = read.csv("lightgreen_results_Group_V1.csv", header=TRUE)
# we want the log2 fold change 
original_gene_list <- df$log2FoldChange

# name the vector
names(original_gene_list) <- df$X

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

# Exctract significant results (padj < 0.05)
sig_genes_df = subset(df, padj < 0.05)

# From significant results, we want to filter on log2fold change
genes <- sig_genes_df$log2FoldChange

# Name the vector
names(genes) <- sig_genes_df$X

# omit NA values
genes <- na.omit(genes)

# filter on min log2fold change (log2FoldChange > 2)
genes <- names(genes)[abs(genes) > 0.5]



go_enrich <- enrichGO(gene = genes,
                      universe = names(gene_list),
                      OrgDb = organism, 
                      keyType = 'ENSEMBL',
                      readable = T,
                      ont = "BP",
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.10)

BiocManager::install("enrichplot")
library(enrichplot)
upsetplot(go_enrich)

barplot(go_enrich, 
        drop = TRUE, 
        showCategory = 30, 
        title = "GO MF",
        font.size = 12)

# categorySize can be either 'pvalue' or 'geneNum'
cnetplot(go_enrich, categorySize="pvalue", foldChange=gene_list)
dotplot(go_enrich)

emapplot(go_enrich)
#----------------------------------------------------------------------------------
# 
#----------------------------------------------------------------------------------
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
write.csv(res, file = "results_Group_TRIAL.csv", quote = FALSE)


df = read.csv("results_Group_TRIAL.csv", header=TRUE)
# we want the log2 fold change 
original_gene_list <- df$log2FoldChange

# name the vector
names(original_gene_list) <- df$X

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

# Exctract significant results (padj < 0.05)
sig_genes_df = subset(df, padj < 0.05)

# From significant results, we want to filter on log2fold change
genes <- sig_genes_df$log2FoldChange

# Name the vector
names(genes) <- sig_genes_df$X

# omit NA values
genes <- na.omit(genes)

# filter on min log2fold change (log2FoldChange > 2)
#genes <- names(genes)[abs(genes) > 2]
genes <- names(genes)[genes < -2]


go_enrich <- enrichGO(gene = genes,
                      universe = names(gene_list),
                      OrgDb = organism, 
                      keyType = 'ENSEMBL',
                      readable = T,
                      ont = "BP",
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.10)

BiocManager::install("enrichplot")
library(enrichplot)
upsetplot(go_enrich)

barplot(go_enrich, 
        drop = TRUE, 
        showCategory = 30, 
        title = "GO Biological Pathways",
        font.size = 8)

# categorySize can be either 'pvalue' or 'geneNum'
cnetplot(go_enrich, categorySize="pvalue", foldChange=gene_list)

##############################################################
library("ggbeeswarm")
geneCounts <- plotCounts(dds75, gene = "ENSMUSG00000029371", intgroup = "Group",
                         returnData = TRUE)
ggplot(geneCounts, aes(x = Group, y = count, color = Group)) +
  scale_y_log10() +  geom_beeswarm(cex = 5)+ theme_bw()
