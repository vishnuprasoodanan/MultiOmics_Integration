install.packages("msigdbr")
library(msigdbr)
library(clusterProfiler)
categories <- msigdbr_collections()
# all_gene_sets = msigdbr(species = "Mus musculus")
# head(all_gene_sets)
# msigdbr_species()
# 
# h_gene_sets = msigdbr(species = "mouse", category = "H")
# head(h_gene_sets)
# 
# cgp_gene_sets = msigdbr(species = "mouse", category = "C2", subcategory = "CGP")
# head(cgp_gene_sets)
# 
# all_gene_sets %>% dplyr::filter(gs_cat == "H") %>% head()
#--------------------------------------------------------------------------------------------------------------
####### function
load_gene_expr <- function(filename){
  genes <- data.frame(fread(filename,sep="\t",head=T), row.names = 1, check.names = F, stringsAsFactors = F)
  genes <- as.matrix(genes)
  
}
#--------------------------------------------------------------------------------------------------------------

## In Rstudio, find the path to the directory where the current script is located.
current_dir <- dirname(rstudioapi::getSourceEditorContext()$path)

#### load data
# Get the names of files ending with .txt in a directory
file_names <- list.files(path = current_dir, pattern = "\\.txt$", full.names = TRUE)

## load gene expression data
for (i in seq_along(file_names)) {
  genes <- load_gene_expr(file_names[i])
  gene_ids_vector <- row.names(genes)
  m_category <- character()
  m_category1 <- character()
  sub_category <- character()
  variables <- strsplit(basename(file_names[i]), '_')[[1]]
  # Loop through the rows of the data frame
  for (k in 1:nrow(categories)) {
    # Check if the second column is empty
    if (is.na(categories[k, 2]) || categories[k, 2] == "") {
      # Save the element in the first column to m_category
      m_category <- as.character(categories[k, 1])
      msigdbr_df = msigdbr(species = "mouse", category = m_category)
      msigdbr_t2g = msigdbr_df %>% dplyr::distinct(gs_name, ensembl_gene) %>% as.data.frame()
      genes_enrich <- enricher(gene = gene_ids_vector, TERM2GENE = msigdbr_t2g, pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.05)
      
      if (!is.null(genes_enrich)) {
        # Write the object to a text file
        output_file2 <- paste0(variables[1], m_category,"_gene_enrich.txt")
        write.table(genes_enrich@result, file = output_file2, sep = "\t", quote = FALSE)
        print(paste0("writing to ",output_file2))
      } else {
        cat("The object is NULL. Moving to the next line of code.\n")
      }

    } else {
      # Save elements in first and second columns to m_category1 and sub_category
      m_category1 <- as.character(categories[k, 1])
      sub_category <- as.character(categories[k, 2])
      msigdbr_df = msigdbr(species = "mouse", category = m_category1, subcategory = sub_category)
      msigdbr_t2g = msigdbr_df %>% dplyr::distinct(gs_name, ensembl_gene) %>% as.data.frame()
      genes_enrich <- enricher(gene = gene_ids_vector, TERM2GENE = msigdbr_t2g, pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.05)
      if (!is.null(genes_enrich)) {
        # Write the object to a text file
        output_file1 <- paste0(variables[1], m_category1, sub_category, "_gene_enrich.txt")
        write.table(genes_enrich@result, file = output_file1, sep = "\t", quote = FALSE)
        print(paste0("writing to ",output_file1))
      } else {
        cat("The object is NULL. Moving to the next line of code.\n")
      }
   }
  }
}
#---------------------------------------------------------------------------------------------------------
#Use the gene sets data frame for clusterProfiler with genes as gene symbols.
msigdbr_t2g = msigdbr_df %>% dplyr::distinct(gs_name, gene_symbol) %>% as.data.frame()
enricher(gene = gene_symbols_vector, TERM2GENE = msigdbr_t2g, ...)

#Use the gene sets data frame for fgsea.
msigdbr_list = split(x = msigdbr_df$gene_symbol, f = msigdbr_df$gs_name)
fgsea(pathways = msigdbr_list, ...)

#Use the gene sets data frame for GSVA.
msigdbr_list = split(x = msigdbr_df$gene_symbol, f = msigdbr_df$gs_name)
gsva(gset.idx.list = msigdbr_list, ...)

#----------------------------------------------------------------------------------------------
