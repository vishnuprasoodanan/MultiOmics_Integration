library(dplyr)
library(NMF)
library(RColorBrewer)
data <- read.table(file = "DISTANCE_MATRIX_v2.txt", header = TRUE, row.names = 1, sep = '\t')
result_matrix <- data.matrix(data)
aheatmap(result_matrix, color = "-BrBG:50", breaks = 0,
         main = "MASH-Pairwise genome distance",
         distfun = 'euclidean',
         hclustfun = "complete",
         fontsize=1,
         filename="FG_nmr_pathway_Association_eucl.pdf")
aheatmap(result_matrix, color = "-BrBG:50", breaks = 0, cellwidth = 10, cellheight =10, border_color = "white",
         main = "IR-Metabolite",
         distfun = 'euclidean',
         hclustfun = "complete",
         fontsize=10,
         filename="IR_metabolites_heatmap_eucl.pdf")
