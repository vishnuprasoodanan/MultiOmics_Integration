# create PCoA using clr-transformed data
library(vegan)
library(ggplot2)
library(gridExtra)
library(data.table)
library(SummarizedExperiment)
library(easyGgplot2)
## In Rstudio, find the path to the directory where the current script is located.
current_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
destination = paste0(current_dir, "/my_plots.pdf")


genes <- read.table("allgenes_tr.txt", sep = "\t", row.names = 1, header = TRUE)

## load microbiome data (note, here we load microbiome data with sex covariate)
counts <- read.table("kraken2_species_readcount_transpose.txt", sep = "\t", row.names = 1, header = TRUE)

## Taxonomy table (to rowData)
tax     <- as.matrix(read.csv("taxonomy.txt", sep = "\t", row.names = 1))

## Metadata
samples <- as.matrix(read.table("metadata.txt", sep="\t", header = TRUE, row.names = 1))  

#----------------------------------------------------------------------------------------------------
counts <- as.matrix(apply(counts, 2, function(x) as.numeric(x)))
se <- SummarizedExperiment(assays = list(counts = counts),
                           colData = samples,
                           rowData = tax)
#tse <- as(se, "TreeSummarizedExperiment")
tse <- transformAssay(se, method = "clr", pseudocount = 0.05)
clr_assay <- assays(tse)$clr
clr_assay <- t(clr_assay)

# Calculates Euclidean distances between samples. Because taxa is in columns,
# it is used to compare different samples.
euclidean_dist <- vegan::vegdist(clr_assay, method = "euclidean")
euclidean_dist_mat <- as.matrix(euclidean_dist)

bc_distances <- bray_curtis_dist <- vegdist(as.matrix(t(genes)), method = "bray")
bc_distance_mat <- as.matrix(bc_distances)


pcoa.gene <- cmdscale(bc_distance_mat, k = nrow(bc_distance_mat)-1, eig = TRUE)
#ordiplot(pcoa.gene, type = "text", main = "PCoA for Species Data, Bray Distances")

pcoa.spe <- cmdscale(euclidean_dist_mat, k = nrow(euclidean_dist_mat)-1, eig = TRUE)
#ordiplot(pcoa.spe, type = "text", main = "PCoA for Species Data, Bray Distances")
#---------------------------------------------------------------------------------------------------
# Perform Procrustes analysis
pro_test <- protest(bc_distance_mat, euclidean_dist_mat, permutations = 999)
procrustes_result <- procrustes(X = bc_distance_mat, Y = euclidean_dist_mat, symmetric = TRUE)
plot(procrustes_result, kind = 1, type = "text")
plot(procrustes_result, kind = 2)


pro_test_pcoa <- protest(pcoa.gene, pcoa.spe, permutations = 999)
procrustes_result_pcoa <- procrustes(X = pcoa.gene, Y = pcoa.spe, symmetric = TRUE)
plot(procrustes_result_pcoa, kind = 1, type = "text")
plot(procrustes_result_pcoa, kind = 2)
#---------------------------------------------------------------------------------------------------

# Perform Mantel test
mantel_result <- mantel(bc_distance_mat, euclidean_dist_mat, method = "spearman", permutations = 999)
mantel_result_pcoa <- mantel(pcoa.gene$points, pcoa.spe$points, method = "spearman", permutations = 999)
# Display the Mantel test result
print(mantel_result)
print(mantel_result_pcoa)

# Run Mantel correlogram
mantel_result <- mantel.correlog(bc_distance_mat, euclidean_dist_mat, cutoff=TRUE, r.type="spearman", nperm=999, mult="holm", progressive=TRUE)
# Print the result
print(mantel_result)
# Plot Mantel correlogram
plot(mantel_result)
