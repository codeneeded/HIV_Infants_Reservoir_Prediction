#### Libraries ####
library(readxl)


#### Paths ####
in.path <- 'C:/Users/axi313/Documents/HIV_Infants_Reservoir_Prediction/Flow_Datasets/'
out.path <- 'C:/Users/axi313/Documents/HIV_Infants_Reservoir_Prediction/Cleaned_Datasets/'

#### Read in Files ####

NK_data <- read_excel(paste0(in.path,"NK cell database 11-28-23.xlsx"))
colnames(NK_data)


#### NK Data Cleaning  ####

# Remove "| " from every column name
colnames(NK_data) <- gsub("\\| ", "", colnames(NK_data))

# Identify columns that do not contain 'MFI', 'Median', or ' (% of CD45)
cols_to_keep <- !grepl("MFI|Median|Freq.|\\(% of CD45\\)|Count", colnames(NK_data), ignore.case = TRUE) | (seq_along(colnames(NK_data)) == 4)

# Subset the data frame to keep only the desired columns
NK_data_cleaned <- NK_data[, cols_to_keep]

# Check the structure of the cleaned data
str(NK_data_cleaned)
colnames(NK_data_cleaned)

# Current column names
current_names <- colnames(NK_data_cleaned)


#### Rename Columns as Parent Poppulations ####

# Define the patterns and their replacements for the new dataset
patterns_and_replacements <- list(
  # Start with the most specific (longest) patterns
  "Singlets/LIVE\\+/CD45, SSC-A subset/CD3-CD20-/CD14-CD123-/HLADR-" = "P7",
  "Singlets/LIVE\\+/CD45, SSC-A subset/CD3-CD20-/CD14-CD123-" = "P6",
  "Singlets/LIVE\\+/CD45, SSC-A subset/CD3-CD20-" = "P5",
  "Singlets/LIVE\\+/CD45, SSC-A subset" = "P4",
  "Singlets/LIVE\\+" = "P3",
  # Then the shorter patterns
  "Singlets" = "P2",
  "Total Count" = "P1",
  "Total NK" = "P8"
  # Add more replacements if needed
)


# Loop through each pattern and replace in the column names
for (pattern in names(patterns_and_replacements)) {
  replacement <- patterns_and_replacements[[pattern]]
  # Use gsub with a regular expression
  current_names <- gsub(paste0("^", pattern, "(?=[/ ]|$)"), replacement, current_names, perl = TRUE)
}

current_names
# Assign the new names to the data frame
colnames(NK_data_cleaned) <- current_names

# Check the new column names
colnames(NK_data_cleaned)


# Define the patterns and their replacements for P9 subsets
p9_patterns_and_replacements <- list(
  "P8/CD56-CD16\\+ NK" = "P9_1",
  "P8/NKbright" = "P9_2",
  "P8/NKdimCD16\\+" = "P9_3",
  "P8/NKdimCD16-" = "P9_4"
)

# Get current column names
current_names <- colnames(NK_data_cleaned)

# Loop through each pattern and replace in the column names
for (pattern in names(p9_patterns_and_replacements)) {
  replacement <- p9_patterns_and_replacements[[pattern]]
  # Use gsub with a regular expression
  current_names <- gsub(paste0("^", pattern, "(?=[/ ]|$)"), replacement, current_names, perl = TRUE)
}

current_names
colnames(NK_data_cleaned)
# Assign the new names to the data frame
colnames(NK_data_cleaned) <- current_names

#### Further Cleaning ####

# Get current column names
current_names <- colnames(NK_data_cleaned)

# Remove "QValue: " and " , " from the column names
current_names <- gsub("Q[0-9]+: ", "", current_names)  # Remove "QValue: "
current_names <- gsub(" , ", "", current_names)        # Remove " , "

# Assign the new names to the data frame
colnames(NK_data_cleaned) <- current_names

# Check the new column names
colnames(NK_data_cleaned)

#### NK Total Cell Count Calculation ####

# List of parent population columns in the order of hierarchy
parent_cols <- c("P1","P2 (%)", "P3 (%)", "P4 (%)", "P5 (%)", "P6 (%)", "P7 (%)", "P8 (%)", "P9_1 (%)", "P9_2 (%)", "P9_3 (%)", "P9_4 (%)")

find_immediate_parent <- function(col, parent_cols) {
  if (col %in% c("P9_1 (%)", "P9_2 (%)", "P9_3 (%)","P9_4 (%)")) {
    return("P8 (%)")
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
  if (!is.na(parent_col) && parent_col %in% colnames(NK_data_cleaned)) {
    NK_data_cleaned[[col]] <- (NK_data_cleaned[[col]] / 100) * NK_data_cleaned[[parent_col]]
  }
}


# Remove space (if present) and "(%)" from all column names
colnames(NK_data_cleaned) <- gsub(" \\(%\\)", "", colnames(NK_data_cleaned))
colnames(NK_data_cleaned)
# Update the parent_cols list accordingly
parent_cols <- gsub(" \\(%\\)|\\(%\\)", "", parent_cols)



# Process columns with one slash
for (col in colnames(NK_data_cleaned)) {
  if (grepl("/", col) && sum(charToRaw(col) == charToRaw("/")) == 1) {
    parent_col = gsub("/.*", "", col) #Parent column is everything before the slash
    NK_data_cleaned[[col]] <- (NK_data_cleaned[[col]] / 100) * NK_data_cleaned[[parent_col]]
  }
}

# Loop through each column starting from the 4th column and round all the numbers to the nearest whole number
for (col in colnames(NK_data_cleaned)[4:ncol(NK_data_cleaned)]) {
  # Check if the column is numeric
  if (is.numeric(NK_data_cleaned[[col]])) {
    NK_data_cleaned[[col]] <- round(NK_data_cleaned[[col]])
  }
}

# Loop through each column starting from the 4th column and replace NA with 0
for (col in colnames(NK_data_cleaned)[4:ncol(NK_data_cleaned)]) {
  # Replace NA values with 0
  NK_data_cleaned[[col]] <- replace(NK_data_cleaned[[col]], is.na(NK_data_cleaned[[col]]), 0)
}

colnames(NK_data_cleaned)


#### Convert cell counts to percentages of P4 (CD45) - Normalization ####

# Start from the 8th column (as the first 7 are PID, Age, Group, P1, P2, P3, and P4)
# Divide by the P4 column (7th column) and multiply by 100

# Find the index of the P4 column
p4_index <- which(colnames(NK_data_cleaned) == "P4")

# Convert cell counts to percentages for columns from P5 onwards
for (i in (p4_index + 1):ncol(NK_data_cleaned)) {
  NK_data_cleaned[, i] <- NK_data_cleaned[, i] / NK_data_cleaned[, p4_index] * 100
}

# Remove uneeded columns
columns_to_remove <- which(colnames(NK_data_cleaned) %in% c("P1", "P2", "P3","P5","P6","P7"))
NK_data_cleaned <- NK_data_cleaned[, -columns_to_remove]

# View the updated dataframe
colnames(NK_data_cleaned)

#### Renaming columns back to populations for further analysis ####

# Renaming columns with gsub
colnames(NK_data_cleaned)[colnames(NK_data_cleaned) == "P4"] <- "CD45 Raw Counts"
colnames(NK_data_cleaned) <- gsub("^P8", "Total_NK", colnames(NK_data_cleaned))
colnames(NK_data_cleaned) <- gsub("^P9_1", "NK_CD56-CD16+", colnames(NK_data_cleaned))
colnames(NK_data_cleaned) <- gsub("^P9_2", "NKbright", colnames(NK_data_cleaned))
colnames(NK_data_cleaned) <- gsub("^P9_3", "NKdim_CD16+", colnames(NK_data_cleaned))
colnames(NK_data_cleaned) <- gsub("^P9_4", "NKdim_CD16-", colnames(NK_data_cleaned))

# Checking the updated column names
colnames(NK_data_cleaned)

#### Save Output ####

output_file <- paste0(out.path, "NK_data_cleaned_normalized_to_CD45.csv")

# Save the dataframe as a CSV file
write.csv(NK_data_cleaned, output_file, row.names = FALSE)
