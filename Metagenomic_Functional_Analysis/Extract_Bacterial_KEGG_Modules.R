# Install necessary packages if not already installed
if (!requireNamespace("KEGGREST", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("KEGGREST")
}

# List of objects to keep
keep_objects <- c("kegg_modules", "organisms", "bacterial_organisms")

# Remove all other objects
rm(list = setdiff(ls(), keep_objects))

library(KEGGREST)
library(dplyr)
library(stringr)

# Function to extract the desired module ID from a string

# Function to extract and validate module IDs
extract_and_validate_module_ids <- function(x) {
  # Split the string by '_'
  parts <- unlist(strsplit(x, "_"))
  
  # Find parts that start with 'M' and have exactly 5 digits
  module_ids <- grep("^M\\d{5}$", parts, value = TRUE)
  
  # Return unique and sorted module IDs
  return(sort(unique(module_ids)))
}

# Function to fetch module links for an organism
get_kegg_links <- function(organism) {
  tryCatch(
    {
      result <- keggLink("module", organism)
      return(result)
    },
    error = function(e) {
      return(paste("Error:", e$message))
    }
  )
}
get_kegg_links_with_retry <- function(organism) {
  repeat {
    result <- get_kegg_links(organism)
    if (all(result != "Error: Forbidden (HTTP 403).")) {
      return(result)
    }
    # Add a small delay to avoid spamming the server with requests
    Sys.sleep(1)
  }
}



# # Fetch the list of all KEGG modules
# kegg_modules <- keggList("module")
# write.table(kegg_modules, file = paste0("Keggrest_", "MODULE_list_", ".txt"), sep = "\t")
# # Fetch the list of all KEGG organisms
# organisms <- keggList("organism")
# write.table(organisms, file = paste0("Keggrest_", "ORGANISM_list", ".txt"), sep = "\t")
# # Identify bacterial organisms by filtering the list
# bacterial_organisms <- organisms[grep("Bacteria", organisms)]

organisms <- as.data.frame(organisms)
# Filter the 'organisms' data frame based on 'phylogeny' column matching 'bacterial_organisms'
bacterial_ids <- organisms %>%
  filter(phylogeny %in% bacterial_organisms) %>%
  select(organism)

# Print the resulting data frame with corresponding 'organism' values
print(bacterial_ids)

#bacterial_ids <- data.frame(organism = c('fbn', 'eco', 'bsu'))

# Initialize lists to store results and processed results
results <- list()
processed_results <- list()
# Loop through each organism in the bacterial_ids dataframe
for (organism in bacterial_ids$organism) {
  results[[organism]] <- get_kegg_links_with_retry(organism)
}

# Now process the results as required
for (organism in names(results)) {
  # Extract and validate module IDs, then store in results
  module_ids <- sapply(results[[organism]], extract_and_validate_module_ids)
  unique_sorted_ids <- unique(unlist(module_ids))
  results[[organism]] <- unique_sorted_ids
}

# Display the results
print(results)

# Combine all elements of 'results' into a single vector
combined_elements <- unlist(results)

# Uniquely sort the array
unique_sorted_elements <- sort(unique(combined_elements))

# Specify the file path
output_file <- "unique_sorted_elements.txt"

# Write the unique sorted elements to the file with each element on a new line
writeLines(unique_sorted_elements, con = output_file)

# Print message to indicate the task is complete
cat("Unique sorted elements have been written to", output_file, "\n")
