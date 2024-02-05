#### Libraries ####
library(readr)
library(dplyr)
library(stringr)


#### Paths ####
in.path <- 'C:/Users/axi313/Documents/HIV_Infants_Reservoir_Prediction/Flow_datasets/'
out.path <- 'C:/Users/axi313/Documents/HIV_Infants_Reservoir_Prediction/Cleaned_datasets/'

#### Read in Files ####
t_b_data <- read_csv(paste0(in.path,"T-B dataset 12-4-23.csv"))
colnames(t_b_data)


#### Clean T_B_Cell Data ####

# Identify columns that do not contain 'MFI', 'Median', or ' (% of CD45)'
cols_to_keep <- !grepl("MFI|Median| \\(% of CD45\\)|Freq. ", colnames(t_b_data), ignore.case = TRUE)

# Subset the data frame to keep only the desired columns
t_b_data_cleaned <- t_b_data[, cols_to_keep]
colnames(t_b_data_cleaned) <- gsub("\\| ", "", colnames(t_b_data_cleaned))

# Check the structure of the cleaned data
str(t_b_data_cleaned)
colnames(t_b_data_cleaned)

#### Clean Column Names ####

# Current column names
current_names <- colnames(t_b_data_cleaned)
# Remove "QValue: " and " , " from the column names
current_names <- gsub("Q[0-9]+: ", "", current_names)  # Remove "QValue: "
current_names <- gsub(" , ", "", current_names)        # Remove " , "
current_names <- gsub(" \\(%\\)", "", current_names)        # Remove " , "
current_names <- gsub("\\(%\\)", "", current_names)        # Remove " , "
colnames(t_b_data_cleaned) <- current_names

# Add path to each column name to trace parent poppuations


names(t_b_data_cleaned) <- gsub("T cells", "Lymphocytes/Single Cells/Live cells/CD45/T cells", names(t_b_data_cleaned))
names(t_b_data_cleaned) <- gsub("CD4\\/", "Lymphocytes/Single Cells/Live cells/CD45/T cells/CD4/", names(t_b_data_cleaned))
names(t_b_data_cleaned) <- gsub("CD8\\/", "Lymphocytes/Single Cells/Live cells/CD45/T cells/CD8/", names(t_b_data_cleaned))

# Replace the column name 'CD4' with the full hierarchical name
names(t_b_data_cleaned)[names(t_b_data_cleaned) == "CD4"] <- "Lymphocytes/Single Cells/Live cells/CD45/T cells/CD4"
names(t_b_data_cleaned)[names(t_b_data_cleaned) == "CD8"] <- "Lymphocytes/Single Cells/Live cells/CD45/T cells/CD8"

colnames(t_b_data_cleaned)

#### Convert Percentages to Raw Counts

# Function to extract immediate parent name based on current column name
get_immediate_parent_name <- function(col_name) {
  parts <- strsplit(col_name, "/")[[1]]
  if (length(parts) > 1) {
    # Remove the last part to get the immediate parent
    return(paste(parts[-length(parts)], collapse = "/"))
  } else {
    # No parent exists (e.g., top-level category)
    return(NA)
  }
}

# Manually set Lymphocytes based on Total count first
t_b_data_cleaned$`Lymphocytes` <- (t_b_data_cleaned$`Lymphocytes` / 100) * t_b_data_cleaned$`Total count`

# Initialize a copy of the dataset for corrected counts
corrected_counts <- t_b_data_cleaned

# Loop through each column to adjust percentages based on immediate parent
for (col_name in names(t_b_data_cleaned)[-c(1:7)]) { # Exclude the first 7 columns that are not hierarchical
  immediate_parent_name <- get_immediate_parent_name(col_name)
  
  # Ensure we're working with percentage columns with a parent
  if (!is.na(immediate_parent_name) && immediate_parent_name %in% names(t_b_data_cleaned)) {
    # Calculate raw counts from percentages
    parent_counts <- corrected_counts[[immediate_parent_name]]
    corrected_counts[[col_name]] <- (t_b_data_cleaned[[col_name]] / 100) * parent_counts
  }
}

# Update original dataframe with corrected counts
t_b_data_cleaned <- corrected_counts

#### Normalize by converting to % of Lymphocytes/Single Cells/Live cells/CD45 column ####

base_column <- t_b_data_cleaned$`Lymphocytes/Single Cells/Live cells/CD45`

# Loop through each column from the 12th column onwards
for (i in 12:ncol(t_b_data_cleaned)) {
  t_b_data_cleaned[, i] <- (t_b_data_cleaned[, i] / base_column) * 100
}

# Dropping columns 1, 4, 5, 6, 7, 8, 9, 10 from t_b_data_cleaned
t_b_data_cleaned <- t_b_data_cleaned[, -c(1, 4, 6, 7, 8, 9, 10)]

#### Column Renaming for easy visualisations ####
colnames(t_b_data_cleaned)
names(t_b_data_cleaned)[names(t_b_data_cleaned) == "Lymphocytes/Single Cells/Live cells/CD45"] <- "T_B_panel_CD45_Raw_Counts"
# Remove "Lymphocytes/Single Cells/Live cells/CD45" from column names
names(t_b_data_cleaned) <- gsub("Lymphocytes/Single Cells/Live cells/CD45/", "", names(t_b_data_cleaned), fixed = TRUE)
names(t_b_data_cleaned) <- gsub("T cells/CD4", "CD4T", names(t_b_data_cleaned), fixed = TRUE)
names(t_b_data_cleaned) <- gsub("T cells/CD8", "CD8T", names(t_b_data_cleaned), fixed = TRUE)

names(t_b_data_cleaned) <- gsub("CD20/B cells", "B cells", names(t_b_data_cleaned), fixed = TRUE)
names(t_b_data_cleaned) <- gsub("CD20/CD19loCD20-/Plasmablasts", "Plasmablasts", names(t_b_data_cleaned), fixed = TRUE)
names(t_b_data_cleaned) <- gsub("B cells/Mature B cells", "Mature_B_Cells", names(t_b_data_cleaned), fixed = TRUE)
names(t_b_data_cleaned) <- gsub("B cells/Transitional", "Transitional_B_Cells", names(t_b_data_cleaned), fixed = TRUE)
names(t_b_data_cleaned)[names(t_b_data_cleaned) == "HIV status"] <- "Group"
colnames(t_b_data_cleaned)

output_file <- paste0(out.path, "T_B_Cell_data_cleaned_normalized_to_CD45.csv")

# Save the dataframe as a CSV file
write.csv(t_b_data_cleaned, output_file, row.names = FALSE)
