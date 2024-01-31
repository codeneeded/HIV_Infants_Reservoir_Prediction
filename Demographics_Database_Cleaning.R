#### Libraries ####
library(readxl)
library(dplyr)


#### Paths ####
in.path <- 'C:/Users/axi313/Documents/HIV_Infants_Reservoir_Prediction/Clinical_Dataset/'
out.path <- 'C:/Users/axi313/Documents/HIV_Infants_Reservoir_Prediction/Cleaned_datasets/'

#### Read in Files ####
demographics <- read_excel(paste0(in.path,"R01_DB_complete 01-05-2024.xlsx"), sheet = 1)
tet_vaccine <- read_excel(paste0(in.path,"Tetanus and Measles Serology Database 6-2-23.xlsx"), sheet = 1)
m_vaccine <- read_excel(paste0(in.path,"Tetanus and Measles Serology Database 6-2-23.xlsx"), sheet = 2)

### Cleaning Demographics database ###

# Subset the first 9 columns as well as the HIV cp/ml column

# Select the first 9 columns
selected_columns <- 1:9

# Find the index of the "HIV cp/mL" column
hiv_col_index <- which(colnames(demographics) == "HIV cp/mL")
status_col_index <- which(colnames(demographics) == "PATIENT TYPE")

# Combine the indices of the first 9 columns with the index of the "HIV cp/mL" column
final_columns <- c(selected_columns, hiv_col_index,status_col_index)

# Create a new data frame with only the selected columns
demographics_selected <- demographics[, final_columns]

# Replace 'HIV cp/mL' values with 0 where 'PATIENT TYPE' is 'HIV Exposito'
demographics_selected <- demographics_selected %>%
  mutate(`HIV cp/mL` = if_else(`PATIENT TYPE` == "HIV Exposto", 0, `HIV cp/mL`))

#Subset the data based on non-NA values in the "HIV cp/mL" column
demographics_filtered <- demographics_selected[!is.na(demographics_selected[["HIV cp/mL"]]), ]

# Round AGE to the nearest Whole Number
demographics_filtered[["AGE IN MONTHS"]] <- round(demographics_filtered[["AGE IN MONTHS"]], digits = 0)

# Rename 'STUDY ID' to 'PID' and drop unwanted columns in demographics_filtered
demographics_filtered <- demographics_filtered %>%
  rename(PID = `STUDY ID`) %>%
  select(-`MOTHER ID`, -`VISIT N`, -`BIRTH DATE`, -`WEIGHT`, -`NEXT VISIT DATE`,-`PATIENT TYPE`)

# Rename values in the 'Sex' column
demographics_filtered <- demographics_filtered %>%
  mutate(SEX = case_when(
    SEX == "Masculino" ~ "Male",
    SEX == "Feminino"  ~ "Female",
    TRUE               ~ SEX  # Keeps original value if neither condition is met
  ))

demographics_filtered$`VISIT DATE` <- as.Date(demographics_filtered$`VISIT DATE`)

# Remove duplicate visits based on AGE IN MONTHS, keeping the earliest visit
demographics_filtered <- demographics_filtered %>%
  group_by(PID, `AGE IN MONTHS`) %>%
  arrange(PID, `AGE IN MONTHS`, `VISIT DATE`) %>%
  filter(row_number() == 1) %>%
  ungroup()

# Drop the 'VISIT DATE' column and rename other specified columns
demographics_filtered <- demographics_filtered %>%
  select(-`VISIT DATE`) %>%
  rename(Sex = SEX, Age = `AGE IN MONTHS`, `Viral Titre` = `HIV cp/mL`)

#### Save Output ####
output_file <- paste0(out.path, "Viral_Titres.csv")

# Save the dataframe as a CSV file
write.csv(demographics_filtered, output_file, row.names = FALSE)

#### Clean the vaccine databases ####
tet_vaccine <- tet_vaccine[1:4]
m_vaccine <- m_vaccine[1:4]

# Rename 'Age in mths' column to 'Age'
tet_vaccine <- tet_vaccine %>%
  rename(Age = `Age in mths`)
# Rename 'Age in mths' column to 'Age'
m_vaccine <- m_vaccine %>%
  rename(Age = `Age in mths`)

#### Save Output ####
tet_output_file <- paste0(out.path, "Tetanus_Serology_Clean.csv")
m_output_file <- paste0(out.path, "Measles_Serology_Clean.csv")

# Save the dataframe as a CSV file
write.csv(demographics_filtered, tet_output_file, row.names = FALSE)
write.csv(demographics_filtered, m_output_file, row.names = FALSE)
