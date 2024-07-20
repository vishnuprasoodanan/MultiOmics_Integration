# Install necessary packages if not already installed
if (!requireNamespace("KEGGREST", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("KEGGREST")
}

library(KEGGREST)
library(dplyr)

# Function to fetch enzyme links for an organism
get_kegg_ec_links <- function(organism) {
  tryCatch(
    {
      result <- keggLink("enzyme", organism)
      return(result)
    },
    error = function(e) {
      return(paste("Error:", e$message))
    }
  )
}

# Function to retry fetching KEGG links until successful
get_kegg_links_with_retry <- function(organism) {
  repeat {
    result <- get_kegg_ec_links(organism)
    if (all(result != "Error: Forbidden (HTTP 403).")) {
      return(result)
    }
    # Add a small delay to avoid spamming the server with requests
    Sys.sleep(1)
  }
}

# Function to extract EC numbers from the KEGG links
extract_ec_numbers <- function(link_data) {
  ec_numbers <- unlist(lapply(link_data, function(x) {
    if (grepl("^ec:", x)) {
      return(sub("ec:", "", x))
    }
  }))
  return(ec_numbers[!is.na(ec_numbers)])
}
#------------------------------------------------------------------------

# Fetch the list of all KEGG organisms
organisms <- keggList("organism")
write.table(organisms, file = paste0("Keggrest_", "ORGANISM_list", ".txt"), sep = "\t")
# Identify bacterial organisms by filtering the list
bacterial_organisms <- organisms[grep("Bacteria", organisms)]

organisms <- as.data.frame(organisms)
# Filter the 'organisms' data frame based on 'phylogeny' column matching 'bacterial_organisms'
bacterial_ids <- organisms %>%
  filter(phylogeny %in% bacterial_organisms) %>%
  select(organism)

# Print the resulting data frame with corresponding 'organism' values
print(bacterial_ids)


# Example bacterial_ids dataframe
#bacterial_ids <- data.frame(organism = c('eco', 'bsu', 'hsa'))

# Initialize list to store results
results <- list()

# Loop through each organism in the bacterial_ids dataframe
for (organism in bacterial_ids$organism) {
  results[[organism]] <- get_kegg_links_with_retry(organism)
}

# Extract EC numbers and process results
ec_results <- list()
for (organism in names(results)) {
  ec_numbers <- extract_ec_numbers(results[[organism]])
  ec_results[[organism]] <- sort(unique(ec_numbers))
}

# Combine all EC numbers into a single vector
combined_ec_numbers <- unlist(ec_results)

# Uniquely sort the array
unique_sorted_ec_numbers <- sort(unique(combined_ec_numbers))

# Specify the file path
output_file <- "unique_sorted_ec_numbers.txt"

# Write the unique sorted EC numbers to the file with each element on a new line
writeLines(unique_sorted_ec_numbers, con = output_file)

# Print message to indicate the task is complete
cat("Unique sorted EC numbers have been written to", output_file, "\n")
