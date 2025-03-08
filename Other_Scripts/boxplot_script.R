# Load required libraries
library(ggplot2)
library(dplyr)
library(ggpubr)  # For arranging multiple plots in one PDF
library(pgirmess)  # For kruskalmc
library(reshape2)  # For melting data
library(ggsignif)  # For adding significance annotations
library(tidyr)  # For data wrangling

# Read the input file
file_path <- readline(prompt = "Enter the path to your input file: ")
cazy_data <- read.table(file_path, header = TRUE, row.names = 1, stringsAsFactors = FALSE, sep = "\t")

# Convert all values to numeric except the 'Status' column
numeric_cols <- setdiff(names(cazy_data), "Status")
cazy_data[, numeric_cols] <- sapply(cazy_data[, numeric_cols], as.numeric)

# Initialize lists for storing results
kruskal_results_list <- list()
pairwise_results_list <- list()
boxplot_list <- list()  # To store plots for combining in a single PDF

# Loop through each numeric column for analysis
for (col in numeric_cols) {
  # Subset the dataframe
  subset_data <- cazy_data %>%
    dplyr::select(Status, all_of(col))
  
  # Perform Kruskal-Wallis test
  kruskal_result <- kruskal.test(subset_data[[col]] ~ subset_data$Status)
  kruskal_results_list[[col]] <- kruskal_result
  
  # Perform pairwise comparisons using kruskalmc for alpha = 0.01
  pairwise_result_01 <- tryCatch(
    {
      kruskalmc(subset_data[[col]], subset_data$Status, alpha = 0.01)
    },
    error = function(e) {
      cat("Error in pairwise comparisons for", col, ":", e$message, "\n")
      return(NULL)
    }
  )
  
  # Initialize sig_annotations
  sig_annotations <- NULL
  
  # Process results for alpha = 0.01
  if (!is.null(pairwise_result_01) && !is.null(pairwise_result_01$dif.com)) {
    pairwise_df <- as.data.frame(pairwise_result_01$dif.com)
    pairwise_df$comparison <- rownames(pairwise_df)
    
    # Split group names
    pairwise_df <- pairwise_df %>%
      separate(comparison, into = c("group1", "group2"), sep = "-")
    
    # Filter significant comparisons and add sign column
    sig_annotations <- pairwise_df %>%
      filter(stat.signif == TRUE) %>%
      select(group1, group2) %>%
      mutate(sign = 0.01)  # Assign alpha = 0.01
  }
  
  # Perform pairwise comparisons using kruskalmc for alpha = 0.05
  pairwise_result_05 <- tryCatch(
    {
      kruskalmc(subset_data[[col]], subset_data$Status, alpha = 0.05)
    },
    error = function(e) {
      cat("Error in pairwise comparisons for", col, ":", e$message, "\n")
      return(NULL)
    }
  )
  
  # Process results for alpha = 0.05
  if (!is.null(pairwise_result_05) && !is.null(pairwise_result_05$dif.com)) {
    pairwise_df_05 <- as.data.frame(pairwise_result_05$dif.com)
    pairwise_df_05$comparison <- rownames(pairwise_df_05)
    
    # Split group names
    pairwise_df_05 <- pairwise_df_05 %>%
      separate(comparison, into = c("group1", "group2"), sep = "-")
    
    # Filter significant comparisons
    new_annotations <- pairwise_df_05 %>%
      filter(stat.signif == TRUE) %>%
      select(group1, group2) %>%
      mutate(sign = 0.05)  # Assign alpha = 0.05
    
    # Check if any new significant comparisons appear at alpha = 0.05
    if (!is.null(sig_annotations)) {
      sig_annotations <- bind_rows(sig_annotations, anti_join(new_annotations, sig_annotations, by = c("group1", "group2")))
    } else {
      sig_annotations <- new_annotations
    }
  }
  
  # Ensure group1 and group2 are character vectors (avoid lists)
  if (!is.null(sig_annotations)) {
    sig_annotations <- sig_annotations %>%
      mutate(group1 = as.character(group1), group2 = as.character(group2))
  }
  
  # Create boxplot with statistical annotations
  p <- ggplot(subset_data, aes(x = Status, y = .data[[col]], fill = Status)) +
    geom_boxplot(alpha = 0.8, outlier.shape = NA) +  # Avoid duplicate outliers
    geom_jitter(shape = 16, position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), alpha = 0.6) +
    scale_fill_manual(values = c("darkolivegreen4", "red", "darkgreen", "salmon4")) +
    labs(title = paste("Boxplot of", col, "by Status"),
         x = "Status",
         y = col) +
    theme_minimal()
  
  # Add significance bars if significant comparisons exist
  if (!is.null(sig_annotations) && nrow(sig_annotations) > 0) {
    p <- p + geom_signif(
      comparisons = as.list(as.data.frame(t(sig_annotations[, c("group1", "group2")]))),
      annotations = ifelse(sig_annotations$sign == 0.01, "**", "*"),
      tip_length = 0.01,  # Keep the small tip length
      step_increase = 0.05,  # Adds vertical spacing between stacked bars
      textsize = 4 
    )
  }
  
  # Save individual boxplot as a separate PDF
  plot_file <- paste0("boxplot_", col, ".pdf")
  pdf(plot_file)
  print(p)
  dev.off()
  cat("Boxplot saved to", plot_file, "\n")
  
  # Store the plot in a list for later combination
  boxplot_list[[col]] <- p
}

# -------------------- 📌 NEW FEATURE: Save all boxplots in one PDF --------------------

if (length(boxplot_list) > 0) {
  cat("Creating all_boxplots.pdf with multiple plots on one page...\n")
  
  pdf("all_boxplots.pdf", width = 25, height = 30)
  print(ggarrange(plotlist = boxplot_list, ncol = 5, nrow = ceiling(length(boxplot_list) / 5)))
  dev.off()
  
  cat("All boxplots saved to all_boxplots.pdf\n")
}

# -------------------- VOLCANO PLOT (Violin + Boxplot) --------------------

# Prepare abundance data (Convert 'cazy_data' into long format)
cazy_data$Sample <- rownames(cazy_data)
long_data <- melt(cazy_data, id.vars = c("Sample", "Status"), variable.name = "Species", value.name = "Abundance")

# Define custom colors
custom_colors <- c("darkolivegreen4", "red", "darkgreen", "salmon4")

# Save violin-boxplot as a PDF
pdf("species_abundance_violin_boxplot.pdf", width = 50, height = 40)

p_violin <- ggplot(long_data, aes(x = Status, y = Abundance, fill = Status)) +
  geom_violin(alpha = 0.6, trim = FALSE) +  # Violin plot to show distribution
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.7, color = "black") +  # Boxplot inside violin plot
  geom_jitter(width = 0.2, size = 1, alpha = 0.7) +  # Jitter for individual data points
  facet_wrap(~Species, scales = "free_y") +  # Facet by species
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +  # Ensure y-axis starts at 0
  scale_fill_manual(values = custom_colors) +  # Apply custom colors
  labs(title = "Species Abundance Distribution",
       x = "Status", y = "Abundance") +
  theme_minimal()

print(p_violin)
dev.off()
cat("Violin-boxplot saved as species_abundance_violin_boxplot.pdf\n")
