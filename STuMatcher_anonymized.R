# STuMatcher 1.0.0 https://github.com/4C554D43
# Desc: generates 5 closest cgmlst target matches for a given sample name
# 
# EXPERIMENTAL BUILD. USE AT YOUR OWN RISK
# Version 1.0.0,  27-11-2023: Initial version
#


# ----- CONFIG -----
#directory for input files. 
input_dir = r"(path)" 
filename_genotype =  "genotype.csv"
target = "target_name"
output_dir = r"(path)"
output_name = "association_report" 
create_logfile = TRUE
create_dir = TRUE #create new folder for each run?
# ----- END CONFIG -----

#include and init
if (!requireNamespace("devtools", quietly = TRUE)) {install.packages("devtools")}
if (!requireNamespace("readxl", quietly = TRUE)) {install.packages("readxl")}
if (!requireNamespace("dplyr", quietly = TRUE)) {install.packages("dplyr")}
library(readxl)
library(dplyr)

# ----- functions -----
timedate_string <- function()
{
  current_datetime <- Sys.time() 
  formatted_datetime <- format(current_datetime, "%Y-%m-%d_%H-%M-%S")
  return (as.character(formatted_datetime))
}#end timedate
print_timestamped <- function(first, second = NULL, third = NULL, fourth = NULL, fifth = NULL, sixth = NULL, se = NULL, ei = NULL, ne = NULL)
{
  message <- paste(first, second, third, fourth, fifth, sixth,se, ei, ne, collapse = " ")
  if ("\n" %in% message) {
    cat(timedate_string(), ": ", message, "\n")
  } else {
    cat(timedate_string(), ": ", message, "\n", sep = "")
  }
}#end print


# ----- main -----
#resolve I/O
if (create_dir == TRUE)
{
  output_dir = paste0(output_dir, timedate_string(), "_Report")
  dir.create(output_dir)
  if (!dir.exists(output_dir))
  {
    print_timestamped("[ERR] Could not create output directory")
    stop()
  }
}
sink(file.path(output_dir, r"(log.txt)"), split = TRUE, append = TRUE)
Sys.Date()
Sys.time()
R.Version()
path_genotype = file.path(input_dir, filename_genotype) #sample data matrix path
if (!file.exists((path_genotype)))
{
  print_timestamped("[ERR] An input file genotype file was not found!")
  stop()
}


genotype_matrix <- as.matrix(read.csv(path_genotype,header = TRUE, row.names = 1, sep=';'))

#cat("Target row: ", which(rownames(genotype_matrix) == target), "\n")

find_closest_matches <- function(genotype_matrix, target) {
  print_timestamped(paste("Assessing differences in", ncol(genotype_matrix), "columns."))
  print_timestamped("Target row: ", which(rownames(genotype_matrix) == target), "\n")
  # Calculate the number of differences for each row compared to the target row
  differences <- apply(genotype_matrix, 1, function(row) sum(row != genotype_matrix[target, ], na.rm = TRUE))
  
  # Exclude the target row itself
  differences[target] <- Inf
  
  # Find the 10 closest matches
  closest_matches_indices <- order(differences)[1:10]
  
  # Print the results
  print_timestamped(paste("Closest matches for", target, "(Position:", which(rownames(genotype_matrix) == target), "):"))
  
  for (match_number in seq_along(closest_matches_indices)) {
    index <- closest_matches_indices[match_number]
    row_name <- rownames(genotype_matrix)[index]
    difference_count <- differences[index]
    
    # Find different column names
    different_columns <- names(which(genotype_matrix[target, ] != genotype_matrix[row_name, ], arr.ind = TRUE))
    
    print_timestamped(paste("Match #", match_number, "- Row:", row_name, "- Differences:", difference_count))
    
    if (length(different_columns) > 0) {
      print_timestamped(paste("  Different columns:", paste(different_columns, collapse = ", ")), "\n")
    }
  }
}
remove_first_columns <- function(matrix, x) {
  if (x >= ncol(matrix)) {
    print("Error: The specified number of columns to remove is greater than or equal to the total number of columns.")
    return(NULL)
  }
  
  new_matrix <- matrix[, (x + 1):ncol(matrix), drop = FALSE]
  return(new_matrix)
}

genotype_matrix = remove_first_columns(genotype_matrix, 8) #exclude sequence type defining genes
find_closest_matches(genotype_matrix, target)


sink()



