# Load required libraries
library(dplyr)
library(tidyr)

diff_sp <- c("Alistipes shahii", "Alistipes senegalensis", "Bifidobacterium breve", "Gordonibacter pamelaeae", "Paraprevotella xylaniphila", "Barnesiella viscericola", "Ruminococcus gauvreauii", "Dysosmobacter sp. Marseille-Q4140", "Coprobacter secundus", "Selenomonas sputigena", "Qiania dongpingensis", "Bacteroides zhangwenhongii", "Bacteroides sp. DH3716P", "Bacteroides coprosuis", "Prevotella melaninogenica", "Bacteroides sp. BFG-257", "Eggerthella guodeyinii", "Bacteroides helcogenes", "Bacteroides heparinolyticus", "Prevotella intermedia", "Serratia marcescens", "Limosilactobacillus reuteri", "Bacteroides zoogleoformans", "Bacteroides caecimuris", "Bacteroides faecium", "Stenotrophomonas maltophilia", "Bacteroides faecis", "Flavonifractor plautii", "Desulfovibrio piger", "Treponema peruense", "Desulfovibrio sp. G11", "Desulfovibrio fairfieldensis", "Desulfovibrio desulfuricans", "Desulfovibrio vulgaris", "Lawsonia intracellularis", "Psychrodesulfovibrio sp. FT415", "Faecalibacterium sp. IP-1-18", "Megalodesulfovibrio gigas", "Desulfocurvibacter africanus", "Pseudodesulfovibrio portus", "Pseudodesulfovibrio tunisiensis", "Desulfovibrio subterraneus", "Pseudodesulfovibrio mercurii", "Desulfomicrobium sp. ZS1", "Desulfomicrobium orale", "Salidesulfovibrio onnuriiensis", "Desulfovibrio sulfodismutans", "Pseudodesulfovibrio indicus", "Solidesulfovibrio carbinoliphilus", "Solidesulfovibrio carbinolicus", "Lactococcus garvieae", "Pseudodesulfovibrio aespoeensis", "Solidesulfovibrio magneticus", "Oleidesulfovibrio alaskensis", "Candidatus Desulfovibrio trichonymphae", "Desulfovibrio sp. 86", "Ruthenibacterium lactatiformans", "Desulfomicrobium baculatum", "Oceanidesulfovibrio marinus", "Desulfovibrio ferrophilus", "Parabacteroides distasonis", "Pseudodesulfovibrio sp. zrk46", "Pseudodesulfovibrio sediminis", "Pseudodesulfovibrio cashew", "Lactobacillus helveticus", "Oscillibacter hominis", "Dysosmobacter welbionis", "Sorangium cellulosum", "Paradesulfovibrio bizertensis", "Christensenella minuta", "Streptococcus parauberis", "Variovorax paradoxus", "Bacteroides salyersiae", "Anaerostipes caccae", "Sellimonas intestinalis", "Bacteroides xylanisolvens", "Clostridium sp. C1", "Phocaeicola salanitronis", "Carnobacterium maltaromaticum", "Longicatena caecimuris", "Subdoligranulum variabile", "Faecalibaculum rodentium", "Parvimonas micra", "Alistipes megaguti", "Faecalibacterium sp. I2-3-92", "Megasphaera stantonii", "Vescimonas fastidiosa", "Rhodobacter capsulatus", "Pusillibacter faecalis", "Alistipes communis", "Alistipes onderdonkii", "Paenibacillus macerans", "Pseudomonas putida", "Desulfoluna sp. ASN36", "Listeria innocua", "Amedibacterium intestinale", "Blautia pseudococcoides", "Blautia producta", "Tannockella kyphosi", "Clostridium gelidum", "Gemella haemolysans", "Intestinibaculum porci", "Neobacillus sp. Marseille-Q6967", "Bifidobacterium longum", "Blautia obeum", "Ruminococcus bovis", "Bacteroides cellulosilyticus", "Enterocloster clostridioformis", "Pseudobutyrivibrio xylanivorans", "Ruminococcus bicirculans", "Bacteroides sp. M10", "Bacteroides ovatus", "Lachnoclostridium sp. YL32", "Bacteroides sp. HF-162", "Parabacteroides sp. CT06", "Bacteroides eggerthii", "Curtobacterium flaccumfaciens", "Akkermansia muciniphila", "Alcanivorax sp. N3-2A", "Pasteurella multocida", "[Clostridium] innocuum", "Desulfobulbus oralis")
mucus_corr_sp <- c("Acidilutibacter cellobiosedens", "Acutalibacter muris", "Amedibacterium intestinale", "Anaerocolumna chitinilytica", "Anaerocolumna sedimenticola", "Anaerostipes hadrus", "Anaerotignum propionicum", "Bacillus cereus", "Bacillus infantis", "Bacillus pumilus", "Bacteroides coprosuis", "Bacteroides eggerthii", "Berryella intestinalis", "Bifidobacterium breve", "Bifidobacterium longum", "Blattabacterium cuenoti", "Blautia obeum", "Blautia producta", "Blautia pseudococcoides", "Bulleidia sp. zg-1006", "Butyrivibrio hungatei", "Butyrivibrio proteoclasticus", "Carnobacterium divergens", "Carnobacterium maltaromaticum", "Catenibacterium mitsuokai", "Clostridioides difficile", "Clostridium manihotivorum", "Clostridium pasteurianum", "Clostridium sp. BNL1100", "Clostridium sp. MD294", "Clostridium tetani", "Desulfosporosinus acidiphilus", "Desulfovibrio sp. 86", "Dorea formicigenerans", "Enterocloster bolteae", "Erysipelothrix larvae", "Eubacterium ventriosum", "Finegoldia magna", "Flavonifractor plautii", "Gemella haemolysans", "Geosporobacter ferrireducens", "Herbinix luporum", "Intestinibaculum porci", "Lachnoanaerobaculum umeaense", "Lachnoclostridium phocaeense", "Lachnoclostridium sp. YL32", "Lachnospira eligens", "Lacrimispora xylanolytica", "Lactobacillus johnsonii", "Latilactobacillus sakei", "Ligilactobacillus murinus", "Longicatena caecimuris", "Marvinbryantia formatexigens", "Niallia circulans", "Oceanobacillus zhaokaii", "Paenibacillus donghaensis", "Paenibacillus polymyxa", "Paenibacillus sp. FSL H7-0357", "Parabacteroides distasonis", "Parabacteroides faecis", "Parabacteroides sp. CT06", "Pelosinus fermentans", "Peptacetobacter hiranonis", "Priestia megaterium", "Pseudobutyrivibrio xylanivorans", "Raoultibacter timonensis", "Ruminococcus bicirculans", "Ruminococcus bovis", "Ruminococcus gauvreauii", "Sedimentibacter sp. zth1", "Sellimonas intestinalis", "Simiaoa sunii", "Solibaculum mannosilyticum", "Solidesulfovibrio carbinolicus", "Solidesulfovibrio carbinoliphilus", "Sporofaciens musculi", "Sporomusa termitida", "Sporosarcina ureae", "Streptococcus gordonii", "Streptococcus pneumoniae", "Vallitalea guaymasensis", "Vallitalea pronyensis", "Vescimonas fastidiosa")
comm_diff_corr_sp <- c("Amedibacterium intestinale", "Bacteroides coprosuis", "Bacteroides eggerthii", "Bifidobacterium breve", "Bifidobacterium longum", "Blautia obeum", "Blautia producta", "Blautia pseudococcoides", "Carnobacterium maltaromaticum", "Desulfovibrio sp. 86", "Flavonifractor plautii", "Gemella haemolysans", "Intestinibaculum porci", "Lachnoclostridium sp. YL32", "Longicatena caecimuris", "Parabacteroides distasonis", "Parabacteroides sp. CT06", "Pseudobutyrivibrio xylanivorans", "Ruminococcus bicirculans", "Ruminococcus bovis", "Ruminococcus gauvreauii", "Sellimonas intestinalis", "Solidesulfovibrio carbinolicus", "Solidesulfovibrio carbinoliphilus", "Vescimonas fastidiosa")
# Plant, Animal, and Mucin vectors
Plant_Cell_Wall_Carbohydrates <- c("GH1", "GH2", "GH3", "GH4", "GH5", "GH8", "GH9", "GH11", "GH12", "GH15", "GH16", "GH17", "GH26", "GH27", "GH28", "GH29", "GH36", "GH39", "GH43", "GH44", "GH48", "GH51", "GH53", "GH55", "GH67", "GH74", "GH78", "GH93", "GH94", "GH95", "GH115", "GH117", "GH121", "PL1", "PL2", "PL6", "PL7", "PL9", "PL11", "PL15", "PL22")
Animal_Carbohydrates <- c("GH1", "GH2", "GH3", "GH4", "GH18", "GH19", "GH20", "GH29", "GH33", "GH38", "GH58", "GH79", "GH84", "GH85", "GH88", "GH89", "GH92", "GH95", "GH98", "GH99", "GH101", "GH105", "GH109", "GH110", "GH113", "PL6", "PL8", "PL12", "PL13", "PL21")
Sucrose_Fructans <- c("GH32", "GH68", "GH70", "GH91")
Mucin <- c("GH2", "GH20", "GH27", "GH29", "GH33", "GH35", "GH36", "GH42", "GH95", "GH84", "GH85", "GH89", "GH98","GH101", "GH110", "GH129")

# Function to combine plant, animal, sucrose, and mucins and return individual dataframes as well
split_and_summarize <- function(df, 
                                Plant_Cell_Wall_Carbohydrates, 
                                Animal_Carbohydrates, 
                                Sucrose_Fructans, 
                                Mucin) {
  
  # Initialize empty lists to hold column names for each category
  plant_cols <- c()
  animal_cols <- c()
  sucrose_cols <- c()
  mucin_cols <- c()
  
  # Iterate through column names of the input dataframe
  for (col_name in colnames(df)) {
    # Split the column name by underscore
    parts <- strsplit(col_name, "_")[[1]]
    # Use the first part if split, otherwise use the entire name
    first_part <- parts[1]
    
    # Check if the selected part is in any of the categories
    if (first_part %in% Plant_Cell_Wall_Carbohydrates) {
      plant_cols <- c(plant_cols, col_name)
    }
    if (first_part %in% Animal_Carbohydrates) {
      animal_cols <- c(animal_cols, col_name)
    }
    if (first_part %in% Sucrose_Fructans) {
      sucrose_cols <- c(sucrose_cols, col_name)
    }
    if (first_part %in% Mucin) {
      mucin_cols <- c(mucin_cols, col_name)
    }
  }
  
  # Create dataframes for each category
  Plant_Cell_Wall_Carbohydrates_df <- df %>% select(all_of(plant_cols))
  Animal_Carbohydrates_df <- df %>% select(all_of(animal_cols))
  Sucrose_Fructans_df <- df %>% select(all_of(sucrose_cols))
  Mucin_df <- df %>% select(all_of(mucin_cols))
  
  # Calculate the sum of values for each row across all columns
  Plant_Cell_Wall_Carbohydrates_sum <- rowSums(Plant_Cell_Wall_Carbohydrates_df)
  Animal_Carbohydrates_sum <- rowSums(Animal_Carbohydrates_df)
  Sucrose_Fructans_sum <- rowSums(Sucrose_Fructans_df)
  Mucin_sum <- rowSums(Mucin_df)
  
  # Combine the sums into a single dataframe
  combined_df <- data.frame(
    Plant_Cell_Wall_Carbohydrates = Plant_Cell_Wall_Carbohydrates_sum,
    Animal_Carbohydrates = Animal_Carbohydrates_sum,
    Sucrose_Fructans = Sucrose_Fructans_sum,
    Mucin = Mucin_sum
  )
  
  return(list(
    combined_df = combined_df,
    Plant_Cell_Wall_Carbohydrates_df = Plant_Cell_Wall_Carbohydrates_df,
    Animal_Carbohydrates_df = Animal_Carbohydrates_df,
    Sucrose_Fructans_df = Sucrose_Fructans_df,
    Mucin_df = Mucin_df
  ))
}

# Function to subset a dataframe if columns and rows are provided as vectors
subset_dataframe <- function(df, selected_columns, selected_rows) {
  # Check if selected_columns is empty
  if (length(selected_columns) == 1 && selected_columns == "") {
    selected_columns <- colnames(df)
  } else {
    # Exclude columns that are not in the dataframe
    selected_columns <- selected_columns[selected_columns %in% colnames(df)]
  }
  
  # Check if selected_rows is empty
  if (length(selected_rows) == 1 && selected_rows == "") {
    selected_rows <- 1:nrow(df)
  } else {
    # Exclude rows that are not in the dataframe
    selected_rows <- selected_rows[selected_rows %in% rownames(df)]
  }
  
  # Subset the dataframe
  df_subset <- df[selected_rows, selected_columns, drop = FALSE]
  
  return(df_subset)
}
# Read the correlation tables as two dataframes
all_corr_data <- read.csv("all_associations.txt", sep = "\t", header = TRUE)
sub_corr_data <- read.csv("sub_associations.txt", sep = "\t", header = TRUE)

create_heatmap <- function(df, output_file, h, w) {
  library(ggplot2)
  library(reshape2)
  library(scales)  # For rescale function
  
  # Convert row names to a column
  df <- cbind(Species = rownames(df), df)
  rownames(df) <- NULL
  
  # Melt the dataframe for easier plotting
  melted_df <- melt(df, id.vars = "Species", variable.name = "CAZy_Family", value.name = "Correlation")
  
  # Calculate the sizes based on the modulus of non-zero values
  max_abs_value <- max(abs(melted_df$Correlation), na.rm = TRUE)
  melted_df$Size <- ifelse(melted_df$Correlation != 0, abs(melted_df$Correlation) / max_abs_value, NA)
  
  # Create a color gradient for the correlation values
  melted_df$Color <- ifelse(melted_df$Correlation > 0, "positive", "negative")
  
  # Generate a color palette for positive and negative values
  color_palette <- scale_color_manual(values = c("positive" = "red4", "negative" = "midnightblue"))
  
  # Create the heatmap with ggplot2
  p <- ggplot(melted_df, aes(x = CAZy_Family, y = Species, size = Size, fill = Correlation)) +
    geom_point(shape = 21) +  # Use shape 21 to have fill without border
    scale_size_continuous(range = c(0, 5), limits = c(0, 1), breaks = seq(0, 1, 0.2), labels = scales::percent) +  # Adjust size range
    scale_fill_gradient2(low = "midnightblue", mid = "white", high = "red4", midpoint = 0, limits = c(-1, 1), oob = squish) +  # Set color gradient
    theme_minimal() +  # Use a minimal theme to remove grey background
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(title = "Correlation Heatmap", x = "CAZy Family", y = "Species", size = "Correlation Magnitude") +
    guides(color = "none", size = guide_legend(title = "Correlation Magnitude"))
  ggsave(filename = output_file, plot = p, device = "pdf", height = h, width = w, units = "in", limitsize = FALSE)
}

# Sanity check function
sanity_check <- function(df) {
  # Check column names
  expected_colnames <- c('X_features', 'Y_features', 'association', 'p.values', 'q.values')
  if (!all(colnames(df) == expected_colnames)) {
    stop("Column names are not in the expected order.")
  }
  
  # Check data types
  if (!all(sapply(df[, 1:2], is.character))) {
    stop("The first two columns are not character values.")
  }
  
  if (!all(sapply(df[, 3:5], is.numeric))) {
    stop("The 3rd, 4th, and 5th columns are not numeric values.")
  }
  
  return(TRUE)
}

# Perform sanity checks
sanity_check(all_corr_data)
sanity_check(sub_corr_data)

# Calculate how many times each unique value in the first two columns repeat and print it
count_repeats <- function(df) {
  x_counts <- count(df$X_features)
  y_counts <- count(df$Y_features)
  
  print(x_counts)
  print(y_counts)
}

count_repeats(all_corr_data)
count_repeats(sub_corr_data)

# Define the function to subset the dataframe
subset_nonzero <- function(df) {
  # Remove rows with all elements equal to zero
  df <- df[rowSums(df != 0) > 0, ]
  
  # Remove columns with all elements equal to zero
  df <- df[, colSums(df != 0) > 0]
  
  return(df)
}

# Function to subset the dataframe by selecting top 50 rows and columns with the largest number of non-zero values
subset_top_nonzero <- function(df) {
  # Count non-zero values in each row and column
  nonzero_rows <- rowSums(df != 0)
  nonzero_cols <- colSums(df != 0)
  
  # Select top 50 rows and columns with the largest number of non-zero values
  top_rows <- order(nonzero_rows, decreasing = TRUE)[1:min(50, length(nonzero_rows))]
  top_cols <- order(nonzero_cols, decreasing = TRUE)[1:min(50, length(nonzero_cols))]
  
  # Subset the dataframe
  subset_df <- df[top_rows, top_cols, drop = FALSE]
  
  return(subset_df)
}

# Function to report the number of positive, negative, and zero values in each row and column
report_values <- function(df) {
  if (nrow(df) == 0 || ncol(df) == 0) {
    stop("The dataframe has no rows or columns.")
  }
  
  # Ensure row names and column names are assigned
  if (is.null(rownames(df))) {
    rownames(df) <- paste0("Row", 1:nrow(df))
  }
  if (is.null(colnames(df))) {
    colnames(df) <- paste0("Col", 1:ncol(df))
  }
  
  report <- list()
  
  # Report for rows
  row_report <- data.frame(
    Row = rownames(df),
    Positive = rowSums(df > 0),
    Negative = rowSums(df < 0),
    Zero = rowSums(df == 0)
  )
  
  # Report for columns
  col_report <- data.frame(
    Column = colnames(df),
    Positive = colSums(df > 0),
    Negative = colSums(df < 0),
    Zero = colSums(df == 0)
  )
  
  report$rows <- row_report
  report$columns <- col_report
  
  return(report)
}

# Get unique features
species <- unique(all_corr_data$X_features)
cazy <- unique(all_corr_data$Y_features)
# Create an empty matrix
all_corr_data_pw_matrix <- matrix(NA, nrow = length(species), ncol = length(cazy))
rownames(all_corr_data_pw_matrix) <- species
colnames(all_corr_data_pw_matrix) <- cazy

# Fill the matrix with correlations
for (i in 1:nrow(all_corr_data)) {
  f1 <- all_corr_data$X_features[i]
  f2 <- all_corr_data$Y_features[i]
  corr <- all_corr_data$association[i]
  all_corr_data_pw_matrix[f1, f2] <- corr
}

all_qval_data_pw_matrix <- matrix(NA, nrow = length(species), ncol = length(cazy))
rownames(all_qval_data_pw_matrix) <- species
colnames(all_qval_data_pw_matrix) <- cazy

# Fill the matrix with correlations
for (i in 1:nrow(all_corr_data)) {
  f1 <- all_corr_data$X_features[i]
  f2 <- all_corr_data$Y_features[i]
  corr <- all_corr_data$q.values[i]
  all_qval_data_pw_matrix[f1, f2] <- corr
}

# Fill NA with zeros
all_corr_data_pw_matrix[is.na(all_corr_data_pw_matrix)] <- 0
all_qval_data_pw_matrix[is.na(all_qval_data_pw_matrix)] <- 0

# Step 4: Print the first 5 rows and columns of the matrices
print(all_corr_data_pw_matrix[1:5, 1:5])
print(all_qval_data_pw_matrix[1:5, 1:5])

# Step 5: Verify the order of row names and column names
if (identical(rownames(all_corr_data_pw_matrix), rownames(all_qval_data_pw_matrix)) &&
    identical(colnames(all_corr_data_pw_matrix), colnames(all_qval_data_pw_matrix))) {
  print("Order of row names and column names in all_corr_data_pw_matrix and all_qval_data_pw_matrix are in order")
}

# Step 6: Replace associations where q.values <= 0.05 with zero
all_corr_data_pw_matrix[all_qval_data_pw_matrix >= 0.05] <- 0

# Get unique features
species <- unique(sub_corr_data$X_features)
cazy <- unique(sub_corr_data$Y_features)
# Create an empty matrix
sub_corr_data_pw_matrix <- matrix(NA, nrow = length(species), ncol = length(cazy))
rownames(sub_corr_data_pw_matrix) <- species
colnames(sub_corr_data_pw_matrix) <- cazy

# Fill the matrix with correlations
for (i in 1:nrow(sub_corr_data)) {
  f1 <- sub_corr_data$X_features[i]
  f2 <- sub_corr_data$Y_features[i]
  corr <- sub_corr_data$association[i]
  sub_corr_data_pw_matrix[f1, f2] <- corr
}

sub_qval_data_pw_matrix <- matrix(NA, nrow = length(species), ncol = length(cazy))
rownames(sub_qval_data_pw_matrix) <- species
colnames(sub_qval_data_pw_matrix) <- cazy

# Fill the matrix with correlations
for (i in 1:nrow(sub_corr_data)) {
  f1 <- sub_corr_data$X_features[i]
  f2 <- sub_corr_data$Y_features[i]
  corr <- sub_corr_data$q.values[i]
  sub_qval_data_pw_matrix[f1, f2] <- corr
}

# Fill NA with zeros
sub_corr_data_pw_matrix[is.na(sub_corr_data_pw_matrix)] <- 0
sub_qval_data_pw_matrix[is.na(sub_qval_data_pw_matrix)] <- 0

# Step 4: Print the first 5 rows and columns of the matrices
print(sub_corr_data_pw_matrix[1:5, 1:5])
print(sub_qval_data_pw_matrix[1:5, 1:5])

# Step 5: Verify the order of row names and column names
if (identical(rownames(sub_corr_data_pw_matrix), rownames(sub_qval_data_pw_matrix)) &&
    identical(colnames(sub_corr_data_pw_matrix), colnames(sub_qval_data_pw_matrix))) {
  print("Order of row names and column names in sub_corr_data_pw_matrix and sub_qval_data_pw_matrix are in order")
}

# Step 6: Replace associations where q.values <= 0.05 with zero
sub_corr_data_pw_matrix[sub_qval_data_pw_matrix >= 0.05] <- 0
write.table(sub_corr_data_pw_matrix, "sub_corr_data_pw_matrix.txt", sep = "\t", col.names = NA)
write.table(all_corr_data_pw_matrix, "all_corr_data_pw_matrix.txt", sep = "\t", col.names = NA)
#----------------------------------------------------------------------------------------------------------------------------
# Step 2: Check the position of zeros in sub_corr_data_pw_matrix
zero_positions <- which(sub_corr_data_pw_matrix == 0, arr.ind = TRUE)

# Initialize a new matrix with the same dimensions and names as all_corr_data_pw_matrix filled with zeros
unique_in_all_corr <- matrix(0, nrow = nrow(all_corr_data_pw_matrix), ncol = ncol(all_corr_data_pw_matrix))
rownames(unique_in_all_corr) <- rownames(all_corr_data_pw_matrix)
colnames(unique_in_all_corr) <- colnames(all_corr_data_pw_matrix)

# Assign values from all_corr_data_pw_matrix to unique_in_all_corr at the positions of zeros in sub_corr_data_pw_matrix
for (pos in seq_len(nrow(zero_positions))) {
  row_idx <- zero_positions[pos, 1]
  col_idx <- zero_positions[pos, 2]
  unique_in_all_corr[row_idx, col_idx] <- all_corr_data_pw_matrix[row_idx, col_idx]
}

# Output the new matrix
print(unique_in_all_corr)
# Subset the dataframe using the function
subsetted_unique_in_all_corr <- subset_nonzero(unique_in_all_corr)
# Subset the dataframe
subsetted_unique_in_all_corr_sub <- subset_top_nonzero(subsetted_unique_in_all_corr)
# Report the values
subsetted_unique_in_all_corr_sub_report <- report_values(subsetted_unique_in_all_corr_sub)
result <- split_and_summarize(as.data.frame(subsetted_unique_in_all_corr), Plant_Cell_Wall_Carbohydrates, Animal_Carbohydrates, Sucrose_Fructans, Mucin)
Mucin_df_diff_sp <- subset_dataframe(result$Mucin_df, c(""), diff_sp)
Mucin_df_corr_sp <- subset_dataframe(result$Mucin_df, c(""), mucus_corr_sp)
comm_sp_df <- subset_dataframe(result$Mucin_df, c(""), comm_diff_corr_sp)
# Optionally: Write the new matrix to an output file
write.table(result$Mucin_df, "Mucin_unique_in_all_corr.txt", sep = "\t", col.names = NA)
write.table(result$Plant_Cell_Wall_Carbohydrates_df, "Plant_unique_in_all_corr.txt", sep = "\t", col.names = NA)
write.table(result$Animal_Carbohydrates_df, "Animal_unique_in_all_corr.txt", sep = "\t", col.names = NA)
write.table(result$Sucrose_Fructans_df, "Sucrose_unique_in_all_corr.txt", sep = "\t", col.names = NA)

write.table(subsetted_unique_in_all_corr, "unique_in_all_corr.txt", sep = "\t", col.names = NA)
write.table(subsetted_unique_in_all_corr_sub, "high_corr_unique_in_all_corr.txt", sep = "\t", col.names = NA)
write.table(subsetted_unique_in_all_corr_sub_report, "report_high_corr_unique_in_all_corr.txt", sep = "\t", col.names = NA)
create_heatmap(as.data.frame(subsetted_unique_in_all_corr_sub), "top_unique_in_all_heatmap.pdf", 15, 15)
create_heatmap(Mucin_df_diff_sp, "DiffSpecies_MucinCAZy_unique_in_all_heatmap.pdf", 20, 15)
create_heatmap(Mucin_df_corr_sp, "CorrSpecies_MucinCAZy_unique_in_all_heatmap.pdf", 20, 15)
create_heatmap(comm_sp_df, "Comm_diff_corr_Species_MucinCAZy_unique_in_all_heatmap.pdf", 15, 15)
#-----------------------------------------------------------------------------------------------------------------------------------------
# Step 2: Check the position of zeros in sub_corr_data_pw_matrix
zero_positions2 <- which(all_corr_data_pw_matrix == 0, arr.ind = TRUE)

# Initialize a new matrix with the same dimensions and names as all_corr_data_pw_matrix filled with zeros
unique_in_sub_corr <- matrix(0, nrow = nrow(sub_corr_data_pw_matrix), ncol = ncol(sub_corr_data_pw_matrix))
rownames(unique_in_sub_corr) <- rownames(sub_corr_data_pw_matrix)
colnames(unique_in_sub_corr) <- colnames(sub_corr_data_pw_matrix)

# Assign values from all_corr_data_pw_matrix to unique_in_sub_corr at the positions of zeros in sub_corr_data_pw_matrix
for (pos in seq_len(nrow(zero_positions2))) {
  row_idx <- zero_positions2[pos, 1]
  col_idx <- zero_positions2[pos, 2]
  unique_in_sub_corr[row_idx, col_idx] <- all_corr_data_pw_matrix[row_idx, col_idx]
}

# Output the new matrix
print(unique_in_sub_corr)
# Subset the dataframe using the function
subsetted_unique_in_sub_corr <- subset_nonzero(unique_in_sub_corr)

# Subset the dataframe
subsetted_unique_in_sub_corr_sub <- subset_top_nonzero(subsetted_unique_in_sub_corr)
# Report the values
subsetted_unique_in_sub_corr_sub_report <- report_values(subsetted_unique_in_sub_corr_sub)
# Optionally: Write the new matrix to an output file
write.table(subsetted_unique_in_sub_corr, "unique_in_sub_corr.txt", sep = "\t", col.names = NA)
write.table(subsetted_unique_in_sub_corr_sub, "high_corr_unique_in_sub_corr.txt", sep = "\t", col.names = NA)
write.table(subsetted_unique_in_sub_corr_sub_report, "report_high_corr_unique_in_sub_corr.txt", sep = "\t", col.names = NA)

#------------------------------------------------------------------------------------------------------------------------------------------
# Step 1: Identify positions in sub_corr_data_pw_matrix with values between -5 and 5
range_positions <- which(sub_corr_data_pw_matrix > -0.3 & sub_corr_data_pw_matrix < 0.3, arr.ind = TRUE)

# Check if range_positions is a matrix and has rows to process
if (is.matrix(range_positions) && nrow(range_positions) > 0) {
  # Initialize a new matrix with the same dimensions and names as all_corr_data_pw_matrix filled with zeros
  increase_in_all_corr <- matrix(0, nrow = nrow(all_corr_data_pw_matrix), ncol = ncol(all_corr_data_pw_matrix))
  rownames(increase_in_all_corr) <- rownames(all_corr_data_pw_matrix)
  colnames(increase_in_all_corr) <- colnames(all_corr_data_pw_matrix)
  
  # Assign values from all_corr_data_pw_matrix to increase_in_all_corr at the identified positions if they have changed to > 5 or < -5
  for (pos in seq_len(nrow(range_positions))) {
    row_idx <- range_positions[pos, 1]
    col_idx <- range_positions[pos, 2]
    all_corr_value <- all_corr_data_pw_matrix[row_idx, col_idx]
    
    if (all_corr_value > 0.3 || all_corr_value < -0.3) {
      increase_in_all_corr[row_idx, col_idx] <- all_corr_value
    }
  }
  
  # Output the new matrix
  print(increase_in_all_corr)
  # Subset the dataframe using the function
  subsetted_increase_in_all_corr <- subset_nonzero(increase_in_all_corr)
  # Subset the dataframe
  subsetted_increase_in_all_corr_sub <- subset_top_nonzero(subsetted_increase_in_all_corr)
  # Report the values
  subsetted_increase_in_all_corr_sub_report <- report_values(subsetted_increase_in_all_corr_sub)
  # Optionally: Write the new matrix to an output file
  write.table(subsetted_increase_in_all_corr, "Increased_in_all_corr.txt", sep = "\t", col.names = NA)
  write.table(subsetted_increase_in_all_corr_sub, "high_corr_Increased_in_all_corr.txt", sep = "\t", col.names = NA)
  write.table(subsetted_increase_in_all_corr_sub_report, "report_high_corr_Increased_in_all_corr.txt", sep = "\t", col.names = NA)
  
} else {
  message("No positions found in the specified range or invalid input data.")
}

#------------------------------------------------------------------------------------------------------------------------------------------
# Step 1: Identify positions in sub_corr_data_pw_matrix with values between -5 and 5
range_positions2 <- which(sub_corr_data_pw_matrix < -0.3 & sub_corr_data_pw_matrix > 0.3, arr.ind = TRUE)

# Check if range_positions is a matrix and has rows to process
if (is.matrix(range_positions2) && nrow(range_positions2) > 0) {
  # Initialize a new matrix with the same dimensions and names as all_corr_data_pw_matrix filled with zeros
  decrease_in_all_corr <- matrix(0, nrow = nrow(all_corr_data_pw_matrix), ncol = ncol(all_corr_data_pw_matrix))
  rownames(decrease_in_all_corr) <- rownames(all_corr_data_pw_matrix)
  colnames(decrease_in_all_corr) <- colnames(all_corr_data_pw_matrix)
  
  # Assign values from all_corr_data_pw_matrix to decrease_in_all_corr at the identified positions if they have changed to > 5 or < -5
  for (pos in seq_len(nrow(range_positions2))) {
    row_idx <- range_positions2[pos, 1]
    col_idx <- range_positions2[pos, 2]
    all_corr_value <- all_corr_data_pw_matrix[row_idx, col_idx]
    
    if (all_corr_value < 0.3 || all_corr_value > -0.3) {
      decrease_in_all_corr[row_idx, col_idx] <- all_corr_value
    }
  }
  
  # Output the new matrix
  print(decrease_in_all_corr)
  # Subset the dataframe using the function
  subsetted_decrease_in_all_corr <- subset_nonzero(decrease_in_all_corr)
  # Optionally: Write the new matrix to an output file
  write.table(subsetted_decrease_in_all_corr, "Decreased_in_all_corr.txt", sep = "\t", col.names = NA)
} else {
  message("No positions found in the specified range or invalid input data.")
}
