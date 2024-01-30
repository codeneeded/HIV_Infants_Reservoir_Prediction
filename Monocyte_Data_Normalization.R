#### Libraries ####
library(readxl)


#### Paths ####
in.path <- 'C:/Users/axi313/Documents/HIV_Infants_Reservoir_Prediction/Flow_datasets/'
out.path <- 'C:/Users/axi313/Documents/HIV_Infants_Reservoir_Prediction/Cleaned_datasets/'

#### Read in Files ####
m_data <- read_excel(paste0(in.path,"Aurora Monocyte Data 11-28-23.xlsx"))
colnames(m_data)


#### Clean Monocyte Data ####

# Identify columns that do not contain 'MFI', 'Median', or ' (% of CD45)'
cols_to_keep <- !grepl("MFI|Median| \\(% of CD45\\)|CD4-|Freq. of HLADR", colnames(m_data), ignore.case = TRUE)

# Subset the data frame to keep only the desired columns
m_data_cleaned <- m_data[, cols_to_keep]

# Check the structure of the cleaned data
str(m_data_cleaned)
colnames(m_data_cleaned)

#### Rename Columns as Parent Populations ####
# Current column names
current_names <- colnames(m_data_cleaned)

# Define the patterns and their replacements
patterns_and_replacements <- list(
  "Total Count" = "P1",
  "Singlets" = "P2",
  "Total LIVE\\+" = "P3",
  "Total CD45\\+" = "P4",
  "Total CD3-CD20-CD56-HLADR\\+" = "P7",
  "Total CD3-CD20-CD56-" = "P6",
  "Total CD3-CD20-" = "P5",
  "Total Monocytes" = "P8",
  "Classical Monocytes" = "P9_1",
  "Inflammatory Monocytes" = "P9_2",
  "Total NCM" = "P9_3",
  "CCR2- NCM" = "P10"
)

# Loop through each pattern and replace in the column names
for (pattern in names(patterns_and_replacements)) {
  replacement <- patterns_and_replacements[[pattern]]
  # Use gsub with a regular expression
  current_names <- gsub(paste0("^", pattern, "(?=[/ ]|$)"), replacement, current_names, perl = TRUE)
}
current_names

# Assign the new names to the data frame
colnames(m_data_cleaned) <- current_names

# Check the new column names
colnames(m_data_cleaned)

#### Monocyte Total Cell Count Calculation ####

# List of parent population columns in the order of hierarchy
parent_cols <- c("P1","P2 (%)", "P3 (%)", "P4 (%)", "P5 (%)", "P6 (%)", "P7 (%)", "P8 (%)", "P9_1 (%)", "P9_2 (%)", "P9_3 (%)", "P10 (%)")

find_immediate_parent <- function(col, parent_cols) {
  if (col %in% c("P9_1 (%)", "P9_2 (%)", "P9_3 (%)")) {
    return("P8 (%)")
  } else if (col == "P10 (%)") {
    return("P9_3 (%)")
  } else {
    index <- match(col, parent_cols)
    if (!is.na(index) && index > 1) {
      return(parent_cols[index - 1])
    } else {
      return(NA)
    }
  }
}


#Calculating Real Values for Parent Columns - We convert each P column on the basis of the previous Gate.
# Convert parent columns first
for (col in parent_cols) {
  parent_col <- find_immediate_parent(col, parent_cols)
  if (!is.na(parent_col) && parent_col %in% colnames(m_data_cleaned)) {
    m_data_cleaned[[col]] <- (m_data_cleaned[[col]] / 100) * m_data_cleaned[[parent_col]]
  }
}


# Remove space (if present) and "(%)" from all column names
colnames(m_data_cleaned) <- gsub(" \\(%\\)", "", colnames(m_data_cleaned))
colnames(m_data_cleaned)
# Update the parent_cols list accordingly
parent_cols <- gsub(" \\(%\\)|\\(%\\)", "", parent_cols)



# Process columns with one slash
for (col in colnames(m_data_cleaned)) {
  if (grepl("/", col) && sum(charToRaw(col) == charToRaw("/")) == 1) {
    parent_col = gsub("/.*", "", col)
    m_data_cleaned[[col]] <- (m_data_cleaned[[col]] / 100) * m_data_cleaned[[parent_col]]
  }
}

# Process columns with two slashes
for (col in colnames(m_data_cleaned)) {
  if (grepl("/", col) && sum(charToRaw(col) == charToRaw("/")) == 2) {
    parent_col = gsub("/[^/]*$", "", col)
    m_data_cleaned[[col]] <- (m_data_cleaned[[col]] / 100) * m_data_cleaned[[parent_col]]
  }
}

# Loop through each column starting from the 4th column
for (col in colnames(m_data_cleaned)[4:ncol(m_data_cleaned)]) {
  # Check if the column is numeric
  if (is.numeric(m_data_cleaned[[col]])) {
    m_data_cleaned[[col]] <- round(m_data_cleaned[[col]])
  }
}

# Loop through each column starting from the 4th column
for (col in colnames(m_data_cleaned)[4:ncol(m_data_cleaned)]) {
  # Replace NA values with 0
  m_data_cleaned[[col]] <- replace(m_data_cleaned[[col]], is.na(m_data_cleaned[[col]]), 0)
}

#### Convert cell counts to percentages of P4 (CD45) - Normalization ####

# Start from the 8th column (as the first 7 are PID, Age, Group, P1, P2, P3, and P4)
# Divide by the P4 column (7th column) and multiply by 100

# Find the index of the P4 column
p4_index <- which(colnames(m_data_cleaned) == "P4")

# Convert cell counts to percentages of P4 for columns from P5 onwards
for (i in (p4_index + 1):ncol(m_data_cleaned)) {
  m_data_cleaned[, i] <- m_data_cleaned[, i] / m_data_cleaned[, p4_index] * 100
}

# Remove unneeded columns
columns_to_remove <- which(colnames(m_data_cleaned) %in% c("P1", "P2", "P3","P5","P6","P7"))
m_data_cleaned <- m_data_cleaned[, -columns_to_remove]

# View the updated dataframe
colnames(m_data_cleaned)

#### Renaming columns back to populations for further analysis ####

# Renaming columns with gsub
colnames(m_data_cleaned)[colnames(m_data_cleaned) == "P4"] <- "Monopanel_CD45_Raw_Counts"
colnames(m_data_cleaned) <- gsub("^P8", "Total_Monocytes", colnames(m_data_cleaned))
colnames(m_data_cleaned) <- gsub("^P9_1", "Classical_Monocytes", colnames(m_data_cleaned))
colnames(m_data_cleaned) <- gsub("^P9_2", "Inflammatory_Monocytes", colnames(m_data_cleaned))
colnames(m_data_cleaned) <- gsub("^P9_3", "NC_Monocytes", colnames(m_data_cleaned))
colnames(m_data_cleaned) <- gsub("^P10", "NC_Monocytes_CCR2-", colnames(m_data_cleaned))



# Checking the updated column names
colnames(m_data_cleaned)

#### Save Output ####

output_file <- paste0(out.path, "Monocyte_data_cleaned_normalized_to_CD45.csv")

# Save the dataframe as a CSV file
write.csv(m_data_cleaned, output_file, row.names = FALSE)
