#### Libraries ####
library(readxl)
library(readr)
library(limma)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(lmerTest)
library(broom.mixed)
library(pheatmap)

#### Paths ####

in.path <- 'C:/Users/axi313/Documents/HIV_Infants_Reservoir_Prediction/Cleaned_datasets/'
out.path <- 'C:/Users/axi313/Documents/HIV_Infants_Reservoir_Prediction/Total_Mixed_Model_Output/'

#### Read in Files ####
t_b_data <- read_csv(paste0(in.path,"T_B_Cell_data_cleaned_normalized_to_CD45.csv"))
m_data <- read_csv(paste0(in.path,"Monocyte_data_cleaned_normalized_to_CD45.csv"))
d_data <- read_csv(paste0(in.path,"DC_data_cleaned_normalized_to_CD45.csv"))
nk_data <- read_csv(paste0(in.path,"NK_data_cleaned_normalized_to_CD45.csv"))
demographics <- read_csv(paste0(in.path,"Viral_Titres.csv"))


#### Clean and combine datasets ####

# Get the names of all columns
column_names <- names(demographics)

# Swap the second and third column names
column_names[c(2, 3)] <- column_names[c(3, 2)]

# Reorder the dataframe based on the updated column names
demographics <- demographics[column_names]

# Inner join the demographics and all other dataframes based on PID and Age (and group when needed)
all_data <- inner_join(demographics, t_b_data, by = c("PID", "Age"))
all_data <- inner_join(all_data, m_data, by = c("PID", "Age","Group"))
all_data <- inner_join(all_data, d_data, by = c("PID", "Age","Group"))
all_data <- inner_join(all_data, nk_data, by = c("PID", "Age","Group"))

all_data <- as.data.frame(all_data)

# Assuming 'group_variable' is your factor variable for groups in 'data_frame'
all_data$Group <- as.factor(all_data$Group)
all_data$Group <- relevel(all_data$Group, ref = "HEU")
all_data$Sex <- as.factor(all_data$Sex)
all_data$PID <- as.factor(all_data$PID)

### Remove Raw Panel Counts
all_data <- select(all_data, 
                   -DCpanel_CD45_Raw_Counts, 
                   -T_B_panel_CD45_Raw_Counts,
                   -NKpanel_CD45_Raw_Counts,
                   -Monopanel_CD45_Raw_Counts)
str(all_data)

#### Linear Modeling with Limma ####

# Extracting relevant columns
covariate_columns <- c("PID", "Age", "Sex", "Viral Titre", "Group")
percentage_columns <- setdiff(names(all_data), covariate_columns)
covariates <- all_data[, covariate_columns]

# Combining PID and Age to create unique identifiers
covariates$unique_id <- paste(covariates$PID, covariates$Age, sep = "_")

# Extracting percentage data
percentage_columns <- setdiff(names(all_data), covariate_columns)
percentage_all_data <- all_data[, percentage_columns]


# Set row names of covariates to unique_id
rownames(covariates) <- covariates$unique_id
rownames(percentage_all_data) <- covariates$unique_id

# Transposing percentage_all_data
percentage_all_data_transposed <- t(percentage_all_data)

filtered_data <- percentage_all_data[complete.cases(percentage_all_data), ]

# Design matrix
design <- model.matrix(~ Group + `Viral Titre` + Age + Sex, data = covariates)
dim(design)

# Voom transformation and linear modeling
v <- voom(percentage_all_data_transposed, design, plot = TRUE)
corfit <- duplicateCorrelation(v, block = covariates$PID)
fit_voom <- lmFit(v, design, block = covariates$PID, correlation = corfit$consensus)
fit_voom <- eBayes(fit_voom)

# Extract results
results_with_voom <- topTable(fit_voom, number=100,coef = "GroupHEI")

# HEU vs HEI, +ve is higher in HEI, -ve is higher in HEU
#### Plot
results_df <- as.data.frame(results_with_voom)
results_df$feature <- rownames(results_df)


# Create the volcano plot
# Add a new column for absolute logFC to help with sorting
results_df$abs_logFC <- abs(results_df$logFC)

# Order the dataframe by absolute logFC and get the top 10
top_features <- results_df %>%
  arrange(desc(abs_logFC)) %>%
  slice(1:10) %>%
  .$feature  # Extract just the feature names

# Flag the top 10 features for labeling
results_df$label_flag <- ifelse(results_df$feature %in% top_features, TRUE, FALSE)

volcano_plot <- ggplot(results_df, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(aes(color = adj.P.Val < 0.05), alpha = 0.5) +
  geom_text_repel(data = filter(results_df, label_flag),  # Use filtered data for labeling
                  aes(label = feature), 
                  box.padding = 0.35, 
                  point.padding = 0.5,
                  segment.color = 'grey50') +
  scale_color_manual(values = c("red")) +
  theme_minimal() +
  labs(title = "Volcano Plot ALL Fold Change HEI Vs HEU", x = "Log2 Fold Change (HEI vs HEU)", y = "-Log10 P-value")


print(volcano_plot)
ggsave(paste0(out.path,"HEIvsHEU_all_data_Volcano.png"),bg='white',plot=volcano_plot)

#### Relationship with Viral Titre - Linear Mixed Model ####

all_data
covariates
percentage_columns

all_data$`Viral Titre` <- scale(all_data$`Viral Titre`)
str(all_data)

#Linear Mixed-Effects Models
# Check for NULL values in the relevant columns
sapply(all_data[c("Viral Titre", "Age", "Sex", "PID","Group")], is.null)

# Check for NA values
sapply(all_data[c("Viral Titre", "Age", "Sex", "PID","Group")], function(x) any(is.na(x)))


# Pre-define the part of the formula that remains constant
fixed_part_of_formula <- "~ `Viral Titre` + Age + Sex + Group + (1 | PID)"

# Initialize lists to store results
models_list <- list()

# Loop through each All Subset, ensuring column names are correctly handled
for (subset_name in percentage_columns) {
  # Escape subset_name if it contains special characters
  escaped_subset_name <- paste0("`", gsub("/", "\\/", subset_name), "`")
  
  # Construct the formula string with escaped column names
  formula_str <- paste(escaped_subset_name, fixed_part_of_formula)
  
  # Convert the string to a formula
  current_formula <- as.formula(formula_str)
  
  
  # Fit the model using the current formula
  model <- tryCatch({
    lmer(current_formula, data = all_data)
  }, error = function(e) {
    cat("Error in fitting model for subset:", subset_name, "\nError message:", e$message, "\n")
    return(NULL)  # Return NULL if there was an error fitting the model
  })
  
  # Store the model if successfully fitted
  if (!is.null(model)) {
    models_list[[subset_name]] <- model
  }
}

# Check Data Structure
summary(models_list[[3]])
summary(models_list[[3]])$coefficients

# After running this loop, models_list will contain the fitted models for each subset

#### Visualizations for Linear mixed Model ####


# Initialize an empty data frame to store the results
results_df <- data.frame(Subset = character(), Effect = character(), Estimate = numeric(), 
                         Std.Error = numeric(), P.Value = numeric(), stringsAsFactors = FALSE)

# Loop through each model to extract information
for (subset_name in names(models_list)) {
  model <- models_list[[subset_name]]
  summary_model <- summary(model)
  df <- as.data.frame(summary_model$coefficients)
  df$Subset <- subset_name
  df$Effect <- rownames(df)
  results_df <- rbind(results_df, df)
}

# Adjust column names if needed based on the structure of your summary
results_df <- results_df %>% 
  rename(Estimate = Estimate, Std.Error = `Std. Error`, P.Value = `Pr(>|t|)`) %>%
  select(Subset, Effect, Estimate, `Std.Error`, P.Value)





### Function for filtering top effects
filter_top_effects <- function(data_frame) {
  positive_effects <- data_frame %>% 
    arrange(desc(Estimate)) %>% 
    head(35)
  
  negative_effects <- data_frame %>% 
    arrange(Estimate) %>% 
    head(35)
  
  bind_rows(positive_effects, negative_effects)
}

# Filter for "Viral Titre" effects

viral_titre_effects <- results_df %>% 
  filter(Effect == "`Viral Titre`") %>%
  mutate(Significance = ifelse(P.Value < 0.05, "*", ""))

# Filter for top 10 effects
viral_titre_effects_filtered <- filter_top_effects(viral_titre_effects)

viral_effect <- ggplot(viral_titre_effects_filtered, aes(x = reorder(Subset, Estimate), y = Estimate)) +
  geom_col(fill = "skyblue") +
  geom_errorbar(aes(ymin = Estimate - `Std.Error`, ymax = Estimate + `Std.Error`), width = 0.4) +
  geom_text(aes(label = Significance), vjust = -0.5, colour = "red") +  # Color the stars red
  coord_flip() +
  labs(title = "Effect of Viral Titre on All Subsets", x = "All Subset", y = "Effect Size") +
  theme_minimal()

ggsave(paste0(out.path,"Viral_Effect_all_data_Coefficient_Plot.png"),bg='white',width=8, height=11.5,plot=viral_effect)

# Filter for "Age" effects (-ve means as age increases subset size decreases)

age_effects <- results_df %>% 
  filter(Effect == "Age") %>%
  mutate(Significance = ifelse(P.Value < 0.05, "*", ""))

# Filter for top 10 effects
age_effects_filtered <- filter_top_effects(age_effects)

# Plot
age_effect <- ggplot(age_effects_filtered, aes(x = reorder(Subset, Estimate), y = Estimate)) +
  geom_col(fill = "skyblue") +
  geom_errorbar(aes(ymin = Estimate - `Std.Error`, ymax = Estimate + `Std.Error`), width = 0.4) +
  geom_text(aes(label = Significance), vjust = -0.5, colour = "red") +
  coord_flip() +
  labs(title = "Effect of Age on All Subsets", x = "All Subset", y = "Effect Size") +
  theme_minimal()

ggsave(paste0(out.path,"Age_Effect_t_b_cell_Coefficient_Plot.png"), bg='white', width=8, height=11.5, plot=age_effect)


# Filter for "SexMale" effects if "SexMale" indicates the effect of being male vs. baseline (female)

sex_effects <- results_df %>% 
  filter(Effect == "SexMale") %>%
  mutate(Significance = ifelse(P.Value < 0.05, "*", ""))

# Filter for top 10 effects
sex_effects_filtered <- filter_top_effects(sex_effects)

# Plot
sex_effect <- ggplot(sex_effects_filtered, aes(x = reorder(Subset, Estimate), y = Estimate)) +
  geom_col(fill = ifelse(sex_effects_filtered$Estimate > 0, "lightblue", "pink")) +
  geom_errorbar(aes(ymin = Estimate - `Std.Error`, ymax = Estimate + `Std.Error`), width = 0.4) +
  geom_text(aes(label = Significance), vjust = -0.5, colour = "red") +
  coord_flip() +
  labs(title = "Effect of Being Male on All Subsets", x = "All Subset", y = "Effect Size Difference") +
  theme_minimal()

ggsave(paste0(out.path,"Sex_Effect_t_b_cell_Coefficient_Plot.png"), bg='white', width=8, height=11.5, plot=sex_effect)


# Filter for "GroupHEI" effects if "GroupHEI" indicates the effect of being in HEU vs. baseline
group_effects <- results_df %>% 
  filter(Effect == "GroupHEI") %>%
  mutate(Significance = ifelse(P.Value < 0.05, "*", ""))

# Filter for top 10 effects
group_effects_filtered <- filter_top_effects(group_effects)

# Plot
HEI_effect <- ggplot(group_effects_filtered, aes(x = reorder(Subset, Estimate), y = Estimate)) +
  geom_col(fill = ifelse(group_effects_filtered$Estimate > 0, "lightgreen", "salmon")) +
  geom_errorbar(aes(ymin = Estimate - `Std.Error`, ymax = Estimate + `Std.Error`), width = 0.4) +
  geom_text(aes(label = Significance), vjust = -0.5, colour = "red") +
  coord_flip() +
  labs(title = "Effect of HEI on All Subsets", x = "All Subset", y = "Effect Size Difference") +
  theme_minimal()

ggsave(paste0(out.path,"HIV_Effect_t_b_cell_Coefficient_Plot.png"), bg='white', width=8, height=11.5, plot=HEI_effect)

