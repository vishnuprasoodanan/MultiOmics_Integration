library(labdsv)
library(dplyr)
library(tidyr)
library(Boruta)
library(ANCOMBC)
library(phyloseq)
library(TreeSummarizedExperiment)
library(mia)
library(ggplot2)
library(gridExtra)
library(VennDiagram)

#set seed value
set.seed(123)
#---------Functions
Colors <- c("darkolivegreen4", "red", "darkgreen", "salmon4")
alpha_value <- 0.5
circle_alpha <- 0.8
##function to calculate indicator value using labdsv
library(labdsv)
indval_summary <- function(famdata){
  famdata$Status <- as.factor(famdata$Status)
  group_names <- sort(levels(as.factor(famdata$Status)))
  iva <- indval(famdata[,2:ncol(famdata)], famdata$Status)
  gr <- iva$maxcls[iva$pval <= 0.05]
  iv <- iva$indcls[iva$pval <= 0.05]
  pv <- iva$pval[iva$pval <= 0.05]
  fr <- apply(famdata[,2:ncol(famdata)] > 0, 2, sum)[iva$pval <= 0.05]
  indvalsummary <- data.frame(group = gr, indval = iv, pvalue = pv, freq = fr)
  indvalsummary <- indvalsummary[order(indvalsummary$group, -indvalsummary$indval),]
  subset_indvalsummary <- indvalsummary %>% filter(indval >= 0.6 & pvalue <= 0.05)
  return(subset_indvalsummary)
}
##function to run Boruta
boruta_function <- function(data, target, pValue = 0.05, mcAdj = TRUE, maxRuns = 100, doTrace = 0, getImp = getImpRfZ){
  boruta.train <- Boruta(data[, 2:ncol(data)], as.factor(data[[target]]), pValue = pValue, mcAdj = mcAdj, maxRuns = maxRuns, doTrace = doTrace, holdHistory = TRUE, getImp = getImp)
  boruta_decision <- as.data.frame(boruta.train$finalDecision)
  write.table(boruta_decision, file = "Boruta_results.txt", sep = "\t")
  
  colnames(boruta_decision)[1] <- "Decision"
  
  confirmed_df <- boruta_decision %>% filter(Decision == "Confirmed")
  
  return(confirmed_df)
}
# Split dataframe containing 'Status' in the first column as metadata, taxon, and count dataframes
split_dataframe <- function(df) {
  # Check if dataframe has at least one column
  if (ncol(df) < 1) {
    stop("Dataframe must have at least one column.")
  }
  
  # Part 1: Dataframe containing rownames and the first column
  sub_metadata <- data.frame(Status = df[, 1], row.names = rownames(df))
  
  # Part 2: Dataframe except the first column
  sub_count <- df[, -1, drop = FALSE]
  
  # Part 3: Dataframe with column names as rownames and same column names in the first column as 'Taxon'
  sub_taxa <- data.frame(Taxon = colnames(df), row.names = colnames(df))
  
  # Remove the first row of sub_taxa
  sub_taxa <- sub_taxa[-1, , drop = FALSE]
  
  return(list(sub_metadata = sub_metadata, sub_count = sub_count, sub_taxa = sub_taxa))
}
#function to run ANCOMBC
ancombc_func <- function(famdata){
  splitted_list <- split_dataframe(famdata)
  se <- SummarizedExperiment(assays = list(counts = t(splitted_list$sub_count)),
                             colData = splitted_list$sub_metadata,
                             rowData = splitted_list$sub_taxa)
  tse <- as(se, "TreeSummarizedExperiment")
  phy <- makePhyloseqFromTreeSE(tse)
  out = ancombc(
    phyloseq = phy, 
    formula = "Status", 
    p_adj_method = "fdr", 
    lib_cut = 0, 
    group = "Status", 
    struc_zero = TRUE, 
    neg_lb = TRUE, 
    tol = 1e-5, 
    max_iter = 100, 
    conserve = TRUE, 
    alpha = 0.01, 
    global = TRUE
  )
  res <- out$res
  ancombc_diff_abn <- out$res$diff_abn
  ancombc_diff_abn_sig <- ancombc_diff_abn[apply(ancombc_diff_abn, 1, function(x) any(x == TRUE)),]
  return(ancombc_diff_abn_sig)
}

#Function to create boxplots for each dataframe and save as PDF
create_boxplots <- function(df, filename) {
  if (nrow(df) == 0) {
    cat(paste("Dataframe", filename, "has zero elements. Skipping...\n"))
    return(NULL)
  }
  
  plots <- lapply(colnames(df)[-1], function(col) {
    p <- ggplot(df, aes(x = Status, y = !!sym(col), fill = Status)) +
      geom_boxplot(alpha = alpha_value, color = "black") +
      geom_jitter(shape = 16, size = 2, alpha = circle_alpha, aes(color = Status)) +
      scale_color_manual(values = Colors) +
      scale_fill_manual(values = Colors) +
      labs(title = col) +
      theme_minimal() +
      theme(legend.position = "none")
    return(p)
  })
  
  plots <- Filter(Negate(is.null), plots)
  if (length(plots) == 0) {
    cat(paste("No valid plots to create for", filename, ". Skipping...\n"))
    return(NULL)
  }
  
  pdf_file <- file.path(paste0(filename, ".pdf"))
  pdf(pdf_file, width = 12, height = 12)
  
  num_plots <- length(plots)
  num_pages <- ceiling(num_plots / 9)
  
  for (page in seq_len(num_pages)) {
    start_plot <- (page - 1) * 9 + 1
    end_plot <- min(page * 9, num_plots)
    
    if (start_plot <= end_plot) {
      plots_subset <- plots[start_plot:end_plot]
      
      if (length(plots_subset) < 9) {
        num_cols <- 3
        num_rows <- 3
        
        if (length(plots_subset) < 4) {
          num_cols <- 3
          num_rows <- 3
          
          if (length(plots_subset) < 2) {
            num_cols <- 3
            num_rows <- 3
          }
        }
        
        empty_plots <- (num_cols * num_rows) - length(plots_subset)
        for (i in 1:empty_plots) {
          plots_subset <- c(plots_subset, NULL)
        }
      } else {
        num_cols <- 3
        num_rows <- 3
      }
      
      plot_grid <- do.call(grid.arrange, c(plots_subset, ncol = num_cols, nrow = num_rows))
      print(plot_grid)
    }
  }
  
  dev.off()
  cat(paste("PDF file", pdf_file, "created successfully.\n"))
}

#function to merge rows of genome abundance dataframe with cluster_data and take average of columns
cluster_merge_avg <- function(famdata_common_4grp, cluster_data){
  famdata_common_4grp_t <- as.data.frame(t(famdata_common_4grp))
  famdata_common_4grp_t_merge <- merge(cluster_data, famdata_common_4grp_t, , by.x = "Bin", by.y = "row.names")
  famdata_common_4grp_t_merge <- famdata_common_4grp_t_merge[,-1]
  famdata_common_4grp_t_merge[,2:ncol(famdata_common_4grp_t_merge)] <- sapply(famdata_common_4grp_t_merge[,2:ncol(famdata_common_4grp_t_merge)], as.numeric)
  famdata_common_4grp_t_merge_avg <- famdata_common_4grp_t_merge %>%
    group_by(Cluster_info) %>%
    summarise_all(mean)
  famdata_common_4grp_t_merge_avg <- as.data.frame(famdata_common_4grp_t_merge_avg)
  rownames(famdata_common_4grp_t_merge_avg) <- famdata_common_4grp_t_merge_avg$Cluster_info
  famdata_common_4grp_t_merge_avg <- famdata_common_4grp_t_merge_avg[,-1]
  return(famdata_common_4grp_t_merge_avg)
}

merge_by_rownames <- function(df1, df2) {
  # Ensure rownames are set correctly
  if (is.null(rownames(df1)) || is.null(rownames(df2))) {
    stop("Both dataframes must have rownames.")
  }
  
  # Merge dataframes by rownames
  common_rownames <- intersect(rownames(df1), rownames(df2))
  
  # Subset dataframes by common rownames
  df1_subset <- df1[common_rownames, , drop = FALSE]
  df2_subset <- df2[common_rownames, , drop = FALSE]
  
  # Add the first column of df2 to df1 as the first column
  df1_with_df2_first_col <- cbind(df2_subset[, 1, drop = FALSE], df1_subset)
  
  return(df1_with_df2_first_col)
}

venn_plot <- function(common_4grp, common_fmt, common_md, filename){
  venn.diagram(
    x = list(common_4grp, common_fmt, common_md),
    category.names = c("4grp", "fmt", "md"),
    filename = filename,
    output = TRUE,
    imagetype = "tiff",
    height = 4800,
    width = 4800,
    resolution = 300,
    lwd = 2,
    fill = c("red", "green", "blue"),
    alpha = 0.5,
    cex = 2,
    fontfamily = "sans",
    fontface = "bold",
    cat.fontfamily = "sans",
    cat.fontface = "bold",
    cat.dist = c(0.1, 0.1, 0.1),
    cat.cex = 2,
    margin = 0.05
  )
}
################################################################################################################### 
#---------Main code
###################################################################################################################
#--------load data and metadata
bin_count <- read.csv("selected_bin_abundance.txt", sep = "\t", header = TRUE, row.names = 1)
metadata <- read.csv("Metadata.txt", sep = "\t", header = TRUE, row.names = 1)
cluster_data <- read.csv("99_clusters.txt", sep = "\t", header = TRUE)

#transpose bin_count and save as dataframe. compare and match row names of both dataframes and merge the transpose of bin_count with metadata. transpose of bin_count should be after metadata
bin_count <- as.data.frame(t(bin_count))
metadata_bin_count <- merge(metadata, bin_count, by = "row.names")
rownames(metadata_bin_count) <- metadata_bin_count$Row.names
famdata <- metadata_bin_count[,-1]

famdata_md <- famdata %>%
  separate(Status, into = c("Status0","Status1", "Status2"), sep = "_") %>%
  select(-Status1, -Status0)
colnames(famdata_md)[1] <- "Status"

famdata_fmt <- famdata %>%
  separate(Status, into = c("Status0","Status1", "Status2"), sep = "_") %>%
  select(-Status2, -Status0)
colnames(famdata_fmt)[1] <- "Status"

labdsv_4grp_indval <- indval_summary(famdata)
labdsv_fmt_indval <- indval_summary(famdata_fmt)
labdsv_md_indval <- indval_summary(famdata_md)

rownames_labdsv_4grp_indval <- rownames(labdsv_4grp_indval)
rownames_labdsv_fmt_indval <- rownames(labdsv_fmt_indval)
rownames_labdsv_md_indval <- rownames(labdsv_md_indval)

boruta_4grp_conf <- boruta_function(famdata, "Status")
boruta_fmt_conf <- boruta_function(famdata_fmt, "Status")
boruta_md_conf <- boruta_function(famdata_md, "Status")

rownames_boruta_4grp_conf <- rownames(boruta_4grp_conf)
rownames_boruta_fmt_conf <- rownames(boruta_fmt_conf)
rownames_boruta_md_conf <- rownames(boruta_md_conf)

abcombc_4grp_conf <- ancombc_func(famdata)
abcombc_fmt_conf <- ancombc_func(famdata_fmt)
abcombc_md_conf <- ancombc_func(famdata_md)

rownames_abcombc_4grp_conf <- abcombc_4grp_conf$taxon
rownames_abcombc_fmt_conf <- abcombc_fmt_conf$taxon
rownames_abcombc_md_conf <- abcombc_md_conf$taxon

#find the common elements in rownames_labdsv_4grp_indval and rownames_boruta_4grp_conf and rownames_abcombc_4grp_conf
common_4grp <- Reduce(intersect, list(rownames_labdsv_4grp_indval, rownames_boruta_4grp_conf, rownames_abcombc_4grp_conf))
common_fmt <- Reduce(intersect, list(rownames_labdsv_fmt_indval, rownames_boruta_fmt_conf, rownames_abcombc_fmt_conf))
common_md <- Reduce(intersect, list(rownames_labdsv_md_indval, rownames_boruta_md_conf, rownames_abcombc_md_conf))

#subset famdata by selecting the columns present in common_4grp
famdata_common_4grp <- famdata[,c("Status", common_4grp)]
famdata_common_fmt <- famdata_fmt[,c("Status", common_fmt)]
famdata_common_md <- famdata_md[,c("Status", common_md)]

#Create boxplots for each dataframe and save as PDF
create_boxplots(famdata_common_4grp, "famdata_common_4grp")
create_boxplots(famdata_common_fmt, "famdata_common_fmt")
create_boxplots(famdata_common_md, "famdata_common_md")

write.table(famdata_common_4grp, file = "famdata_common_4grp.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(famdata_common_fmt, file = "famdata_common_fmt.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(famdata_common_md, file = "famdata_common_md.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#plot a venn diagram to show common elements in common_4grp, common_fmt, and common_md. make the venn diagram more readable and better colors. 3 circles should be in 3 directions
venn_plot(common_4grp, common_fmt, common_md, "bin_venn_plot.tiff")

#combine the column names by taking the average of columns (of famdata_common_4grp) that has same values in cluster_data$Cluster_info
cluster_data$Bin <- gsub(".fa", "", cluster_data$Bin)
cluster_data$Cluster_info <- paste("Cluster99_", cluster_data$Cluster_info, sep = "")

famdata_common_4grp_clust_avg <- cluster_merge_avg(famdata_common_4grp, cluster_data)
famdata_common_4grp_clust_avg <- merge_by_rownames(as.data.frame(t(famdata_common_4grp_clust_avg)), metadata)
famdata_common_4grp_clust_names <- colnames(famdata_common_4grp_clust_avg)[-1]
create_boxplots(famdata_common_4grp_clust_avg, "famdata_common_4grp_cluster99")

famdata_common_fmt_clust_avg <- cluster_merge_avg(famdata_common_fmt, cluster_data)
famdata_common_fmt_clust_avg <- merge_by_rownames(as.data.frame(t(famdata_common_fmt_clust_avg)), metadata)
famdata_common_fmt_clust_names <- colnames(famdata_common_fmt_clust_avg)[-1]
create_boxplots(famdata_common_fmt_clust_avg, "famdata_common_fmt_cluster99")

famdata_common_md_clust_avg <- cluster_merge_avg(famdata_common_md, cluster_data)
famdata_common_md_clust_avg <- merge_by_rownames(as.data.frame(t(famdata_common_md_clust_avg)), metadata)
famdata_common_md_clust_names <- colnames(famdata_common_md_clust_avg)[-1]
create_boxplots(famdata_common_md_clust_avg, "famdata_common_md_cluster99")

venn_plot(famdata_common_4grp_clust_names, famdata_common_fmt_clust_names, famdata_common_md_clust_names,"bin_venn_plot_cluster99.tiff")
