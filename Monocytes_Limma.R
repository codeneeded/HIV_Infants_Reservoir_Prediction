#### Libraries ####
library(readxl)
library(readr)
library(limma)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(lme4)
library(broom.mixed)
library(pheatmap)
#### Paths ####

in.path <- 'C:/Users/axi313/Documents/HIV_Infants_Reservoir_Prediction/Cleaned_datasets/'
out.path <- 'C:/Users/axi313/Documents/HIV_Infants_Reservoir_Prediction/Monocyte_Linear_Model/'

#### Read in Files ####
m_data <- read_csv(paste0(in.path,"Monocyte_data_cleaned_normalized_to_CD45.csv"))
demographics <- read_csv(paste0(in.path,"Viral_Titres.csv"))


### Clean and combine datasets ###

# Get the names of all columns
column_names <- names(demographics)

# Swap the second and third column names
column_names[c(2, 3)] <- column_names[c(3, 2)]

# Reorder the dataframe based on the updated column names
demographics <- demographics[column_names]

# Inner join the demographics and m_data dataframes based on PID and Age
monocytes <- inner_join(demographics, m_data, by = c("PID", "Age"))

monocytes <- as.data.frame(monocytes)
str(monocytes)
#### Linear Modeling with Limma ####

# Extracting relevant columns
covariate_columns <- c("PID", "Age", "Sex", "Viral Titre", "Group")
percentage_columns <- setdiff(names(monocytes), covariate_columns)
covariates <- monocytes[, covariate_columns]

# Combining PID and Age to create unique identifiers
covariates$unique_id <- paste(covariates$PID, covariates$Age, sep = "_")

# Extracting percentage data
percentage_columns <- setdiff(names(monocytes), covariate_columns)
percentage_monocytes <- monocytes[, percentage_columns]


# Set row names of covariates to unique_id
rownames(covariates) <- covariates$unique_id
rownames(percentage_monocytes) <- covariates$unique_id

# Transposing percentage_monocytes
percentage_monocytes_transposed <- t(percentage_monocytes)


# Design matrix
design <- model.matrix(~ Group + `Viral Titre` + Age + Sex, data = covariates)
dim(design)

# Voom transformation and linear modeling
v <- voom(percentage_monocytes_transposed, design, plot = TRUE)
corfit <- duplicateCorrelation(v, block = covariates$PID)
fit_voom <- lmFit(v, design, block = covariates$PID, correlation = corfit$consensus)
fit_voom <- eBayes(fit_voom)


# Extract results
results_with_voom <- topTable(fit_voom, coef = "GroupHEU")

# HEU vs HEI, +ve is higher in HEU, -ve is higher in HEI
#### Plot
results_df <- as.data.frame(results_with_voom)
results_df$feature <- rownames(results_df)


# Create the volcano plot
volcano_plot <- ggplot(results_df, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(aes(color = adj.P.Val < 0.05), alpha = 0.5) +
  geom_text_repel(aes(label = feature), 
                  box.padding = 0.35, 
                  point.padding = 0.5,
                  segment.color = 'grey50') +
  scale_color_manual(values = c("black", "red")) +
  theme_minimal() +
  labs(title = "Volcano Plot", x = "Log2 Fold Change (HEU vs HEI)", y = "-Log10 P-value")

print(volcano_plot)

### Relationship with Viral Titre ###

#Subset Data to HEI Group

hei_data <- monocytes[monocytes$Group == "HEI", ]
hei_covariates <- hei_data[covariate_columns]
hei_percentage_monocytes <- hei_data[percentage_columns]
hei_data$`Viral Titre` <- scale(hei_data$`Viral Titre`)
str(hei_data)
#Linear Mixed-Effects Models
# Check for NULL values in the relevant columns
sapply(hei_data[c("Viral Titre", "Age", "Sex", "PID")], is.null)

# Check for NA values
sapply(hei_data[c("Viral Titre", "Age", "Sex", "PID")], function(x) any(is.na(x)))

# Convert 'Sex' to a factor if it's not already
hei_covariates$Sex <- as.factor(hei_covariates$Sex)
# Drop the 'Group' column from hei_covariates
hei_covariates <- hei_covariates[, !(names(hei_covariates) %in% "Group")]
percentage_columns

#######

# Ensure 'Sex' is correctly formatted as a factor, and 'PID' as a character or factor
hei_data$Sex <- as.factor(hei_data$Sex)
hei_data$PID <- as.factor(hei_data$PID)

# Pre-define the part of the formula that remains constant
fixed_part_of_formula <- "~ `Viral Titre` + Age + Sex + (1 | PID)"

# Initialize lists to store results
models_list <- list()
summaries_list <- list()

# Select just the first monocyte subset for this example
subset_names <- percentage_columns[1]  # Adjust this to test more subsets iteratively

# Loop through each monocyte subset, ensuring column names are correctly handled
for (subset_name in percentage_columns) {
  # Escape subset_name if it contains special characters
  escaped_subset_name <- paste0("`", gsub("/", "\\/", subset_name), "`")
  
  # Construct the formula string with escaped column names
  formula_str <- paste(escaped_subset_name, fixed_part_of_formula)
  
  # Convert the string to a formula
  current_formula <- as.formula(formula_str)
  
  # Fit the model using the current formula
  model <- tryCatch({
    lmer(current_formula, data = hei_data)
  }, error = function(e) {
    cat("Error in fitting model for subset:", subset_name, "\nError message:", e$message, "\n")
    return(NULL)  # Return NULL if there was an error fitting the model
  })
  
  # Store the model if successfully fitted
  if (!is.null(model)) {
    models_list[[subset_name]] <- model
  }
}

# After running this loop, models_list will contain the fitted models for each subset

#### Visualisations for Linear mIxed Model ####

# Initialize an empty data frame to store the results
results_df <- data.frame(
  Subset = character(),
  Estimate = numeric(),
  Std.Error = numeric(),
  CI.Low = numeric(),
  CI.High = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each model to extract information
for (subset_name in percentage_columns) {
  model <- models_list[[subset_name]]
  
  # Extract fixed effects estimates
  estimates <- summary(model)$coefficients
  
  # Check if `Viral Titre` is present in the model summary
  if ("`Viral Titre`" %in% rownames(estimates)) {
    # Extract estimate and standard error for `Viral Titre`
    estimate <- estimates["`Viral Titre`", "Estimate"]
    std_error <- estimates["`Viral Titre`", "Std. Error"]
    
    # Calculate confidence intervals
    conf_int <- confint(model, level = 0.95, method = "Wald")["`Viral Titre`", ]
    
    # Append the results to the results_df
    results_df <- rbind(results_df, data.frame(
      Subset = subset_name,
      Estimate = estimate,
      Std.Error = std_error,
      CI.Low = conf_int[1],
      CI.High = conf_int[2]
    ))
  }
}

### Coefficient Plot for Viral Load Effect Across Subsets

filtered_results_df <- results_df[results_df$Subset != "Monopanel_CD45_Raw_Counts", ]

# Plotting the filtered results
ggplot(filtered_results_df, aes(x = Estimate, y = reorder(Subset, Estimate))) +
  geom_point() +
  geom_errorbarh(aes(xmin = CI.Low, xmax = CI.High), height = 0.2) +
  labs(title = "Effect of Viral Load on Monocyte Subsets",
       x = "Effect Size (Estimate of Viral Load)", y = "Monocyte Subset") +
  theme_minimal()

### Heatmap


# Assuming models_list is your list of lmer model objects
p_values_list <- lapply(models_list, function(model) {
  # Extract the summary
  summary_model <- summary(model)
  
  # Find the row corresponding to 'Viral Titre' and extract the p-value
  # Note: Adjust if your summary structure differs
  if ("`Viral Titre`" %in% rownames(summary_model$coefficients)) {
    p_value <- summary_model$coefficients["`Viral Titre`", "Pr(>|t|)"]
  } else {
    p_value <- NA  # Assign NA if 'Viral Titre' not found
  }
  
  return(p_value)
})

# Add the p-values to your results dataframe
results_df$p_value <- unlist(p_values_list)

# Assuming 'results_df' is structured with 'Subset' and 'Estimate' columns
# Prepare the matrix for the heatmap
heatmap_matrix <- matrix(filtered_results_df$Estimate, nrow = nrow(filtered_results_df), ncol = 1,
                         dimnames = list(filtered_results_df$Subset, c("Viral Load Effect")))

# Plot the heatmap
pheatmap(heatmap_matrix,
         scale = "none",  # Or "row" depending on your preference
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         cluster_cols = FALSE,  # Disable column clustering
         color = colorRampPalette(c("blue", "white", "red"))(100),
         show_rownames = TRUE,
         show_colnames = TRUE,
         main = "Effect of Viral Load on Monocyte Subsets",
         legend_title = "Effect Size")









