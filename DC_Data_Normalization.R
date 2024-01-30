#### Libraries ####
library(readxl)


#### Paths ####
in.path <- 'C:/Users/axi313/Documents/HIV_Infants_Reservoir_Prediction/Flow_datasets/'
out.path <- 'C:/Users/axi313/Documents/HIV_Infants_Reservoir_Prediction/Cleaned_datasets/'

#### Read in Files ####
DC_data <- read_excel(paste0(in.path,"Aurora DC Data 11-28-23.xlsx"))
colnames(DC_data)


#### Clean DC Data ####
# Identify columns that do not contain 'MFI', 'Median', or ' (% of CD45)
cols_to_keep <- !grepl("MFI|Median|Freq.|\\(% of CD45\\)|Count", colnames(DC_data), ignore.case = TRUE) | (seq_along(colnames(DC_data)) == 4)

# Subset the data frame to keep only the desired columns
DC_data_cleaned <- DC_data[, cols_to_keep]

# Check the structure of the cleaned data
str(DC_data_cleaned)
colnames(DC_data_cleaned)

#### Rename Columns as Parent Populations ####

# Current column names
current_names <- colnames(DC_data_cleaned)

# Define the patterns and their replacements for the new dataset
patterns_and_replacements <- list(
  "Total CD3-CD20-CD56-HLADR\\+" = "P7",
  "Total CD3-CD20-CD56-" = "P6",
  "Total CD3-CD20-" = "P5",
  "Total CD45\\+" = "P4",
  "Total LIVE\\+" = "P3",
  "Singlets" = "P2",
  "Total Count" = "P1",
  "Total DC" = "P8",
  "Total pDC" = "P9_1",
  "Total cDC" = "P9_2",
  "cDC1" = "P10_1",
  "cDC2" = "P10_2"
)

# Loop through each pattern and replace in the column names
for (pattern in names(patterns_and_replacements)) {
  replacement <- patterns_and_replacements[[pattern]]
  # Use gsub with a regular expression
  current_names <- gsub(paste0("^", pattern, "(?=[/ ]|$)"), replacement, current_names, perl = TRUE)
}

current_names
# Assign the new names to the data frame
colnames(DC_data_cleaned) <- current_names

# Check the new column names
colnames(DC_data_cleaned)

#### Further Cleaning ####

# Get current column names
current_names <- colnames(DC_data_cleaned)

# Remove "QValue: " and " , " from the column names
current_names <- gsub("Q[0-9]+: ", "", current_names)  # Remove "QValue: "
current_names <- gsub(" , ", "", current_names)        # Remove " , "

# Assign the new names to the data frame
colnames(DC_data_cleaned) <- current_names

# Check the new column names
colnames(DC_data_cleaned)

#### DC Total Cell Count Calculation ####

# Assuming DC_data_cleaned is your dataframe and 'P1' is the base count column

find_immediate_parent <- function(col, parent_cols) {
  if (col %in% c("P9_1 (%)", "P9_2 (%)")) {
    return("P8 (%)")
  } else if (col %in% c("P10_1 (%)", "P10_2 (%)")) {
    return("P9_2 (%)")
  } else {
    index <- match(col, parent_cols)
    if (!is.na(index) && index > 1) {
      return(parent_cols[index - 1])
    } else {
      return(NA)
    }
  }
}


# List of parent population columns in the order of hierarchy
parent_cols <- c("P1", "P2 (%)", "P3 (%)", "P4 (%)", "P5 (%)", "P6 (%)", "P7 (%)", "P8 (%)", "P9_1 (%)", "P9_2 (%)", "P10_1 (%)", "P10_2 (%)")

# Calculating Real Values for Parent Columns
for (col in parent_cols) {
  parent_col <- find_immediate_parent(col, parent_cols)
  if (!is.na(parent_col) && parent_col %in% colnames(DC_data_cleaned)) {
    DC_data_cleaned[[col]] <- (DC_data_cleaned[[col]] / 100) * DC_data_cleaned[[parent_col]]
  }
}


# Remove space (if present) and "(%)" from all column names
colnames(DC_data_cleaned) <- gsub(" \\(%\\)", "", colnames(DC_data_cleaned))
colnames(DC_data_cleaned)
# Update the parent_cols list accordingly
parent_cols <- gsub(" \\(%\\)|\\(%\\)", "", parent_cols)



# Process columns with one slash
for (col in colnames(DC_data_cleaned)) {
  if (grepl("/", col) && sum(charToRaw(col) == charToRaw("/")) == 1) {
    parent_col = gsub("/.*", "", col)
    DC_data_cleaned[[col]] <- (DC_data_cleaned[[col]] / 100) * DC_data_cleaned[[parent_col]]
  }
}

# Process columns with two slashes
for (col in colnames(DC_data_cleaned)) {
  if (grepl("/", col) && sum(charToRaw(col) == charToRaw("/")) == 2) {
    parent_col = gsub("/[^/]*$", "", col)
    DC_data_cleaned[[col]] <- (DC_data_cleaned[[col]] / 100) * DC_data_cleaned[[parent_col]]
  }
}

# Loop through each column starting from the 4th column
for (col in colnames(DC_data_cleaned)[4:ncol(DC_data_cleaned)]) {
  # Check if the column is numeric
  if (is.numeric(DC_data_cleaned[[col]])) {
    DC_data_cleaned[[col]] <- round(DC_data_cleaned[[col]])
  }
}

# Loop through each column starting from the 4th column
for (col in colnames(DC_data_cleaned)[4:ncol(DC_data_cleaned)]) {
  # Replace NA values with 0
  DC_data_cleaned[[col]] <- replace(DC_data_cleaned[[col]], is.na(DC_data_cleaned[[col]]), 0)
}

#### Convert cell counts to percentages of P4 (CD45) - Normalization ####

# Start from the 8th column (as the first 7 are PID, Age, Group, P1, P2, P3, and P4)
# Divide by the P4 column (7th column) and multiply by 100

# Find the index of the P4 column
p4_index <- which(colnames(DC_data_cleaned) == "P4")

# Convert cell counts to percentages for columns from P5 onwards
for (i in (p4_index + 1):ncol(DC_data_cleaned)) {
  DC_data_cleaned[, i] <- DC_data_cleaned[, i] / DC_data_cleaned[, p4_index] * 100
}

# Remove unneeded columns
columns_to_remove <- which(colnames(DC_data_cleaned) %in% c("P1", "P2", "P3","P5","P6","P7"))
DC_data_cleaned <- DC_data_cleaned[, -columns_to_remove]

# View the updated dataframe
colnames(DC_data_cleaned)

#### Renaming columns back to populations for further analysis ####

# Renaming columns with gsub
colnames(DC_data_cleaned)[colnames(DC_data_cleaned) == "P4"] <- "DCpanel_CD45_Raw_Counts"
colnames(DC_data_cleaned) <- gsub("^P8", "Total_DC", colnames(DC_data_cleaned))
colnames(DC_data_cleaned) <- gsub("^P9_1", "Total_pDC", colnames(DC_data_cleaned))
colnames(DC_data_cleaned) <- gsub("^P9_2", "Total_cDC", colnames(DC_data_cleaned))
colnames(DC_data_cleaned) <- gsub("^P10_1", "cDC1", colnames(DC_data_cleaned))
colnames(DC_data_cleaned) <- gsub("^P10_2", "cDC2", colnames(DC_data_cleaned))

# Checking the updated column names
colnames(DC_data_cleaned)

#### Save Output ####

output_file <- paste0(out.path, "DC_data_cleaned_normalized_to_CD45.csv")

# Save the dataframe as a CSV file
write.csv(DC_data_cleaned, output_file, row.names = FALSE)
