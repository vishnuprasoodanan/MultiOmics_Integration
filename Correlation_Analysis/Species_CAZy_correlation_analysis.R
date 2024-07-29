# Load required libraries
library(dplyr)
library(tidyr)

# Read the correlation tables as two dataframes
all_corr_data <- read.csv("all_associations.txt", sep = "\t", header = TRUE)
sub_corr_data <- read.csv("sub_associations.txt", sep = "\t", header = TRUE)

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
# Optionally: Write the new matrix to an output file
write.table(subsetted_unique_in_all_corr, "unique_in_all_corr.txt", sep = "\t", col.names = NA)

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
# Optionally: Write the new matrix to an output file
write.table(subsetted_unique_in_sub_corr, "unique_in_sub_corr.txt", sep = "\t", col.names = NA)

#------------------------------------------------------------------------------------------------------------------------------------------
# Step 1: Identify positions in sub_corr_data_pw_matrix with values between -5 and 5
range_positions <- which(sub_corr_data_pw_matrix > -0.5 & sub_corr_data_pw_matrix < 0.5, arr.ind = TRUE)

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
    
    if (all_corr_value > 0.5 || all_corr_value < -0.5) {
      increase_in_all_corr[row_idx, col_idx] <- all_corr_value
    }
  }
  
  # Output the new matrix
  print(increase_in_all_corr)
  # Subset the dataframe using the function
  subsetted_increase_in_all_corr <- subset_nonzero(increase_in_all_corr)
  # Optionally: Write the new matrix to an output file
  write.table(subsetted_increase_in_all_corr, "Increased_in_all_corr.txt", sep = "\t", col.names = NA)
} else {
  message("No positions found in the specified range or invalid input data.")
}

#------------------------------------------------------------------------------------------------------------------------------------------
# Step 1: Identify positions in sub_corr_data_pw_matrix with values between -5 and 5
range_positions2 <- which(sub_corr_data_pw_matrix < -0.5 & sub_corr_data_pw_matrix > 0.5, arr.ind = TRUE)

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
    
    if (all_corr_value < 0.5 || all_corr_value > -0.5) {
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


