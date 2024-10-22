# Map p-value to significance level symbols
# Example: "***" for p < 0.001, "**" for p < 0.01, etc.
map_p_value_to_pch <- function(p_value) {
  if (p_value < 0.001) {
    return("***")
  } else if (p_value < 0.01) {
    return("**")
  } else if (p_value <= 0.05) {
    return("*")
  } else {
    return(".")
  }
}

# Example usage:
# result <- map_p_value_to_pch(0.004)  # Should return "**"

# Confusion Matrix Function
# For classification results of the MUVR object, create a confusion matrix
confusionMatrix <- function(MVObj, model = 'max') {
  if (!any(class(MVObj) == 'Classification')) {
    stop('The MUVR object needs to be from a classification analysis')
  }
  nMod <- ifelse(model == 'min', 1, ifelse(model == 'mid', 2, 3))
  actual <- MVObj$inData$Y
  predicted <- MVObj$yClass[, nMod]
  table(actual = actual, predicted = predicted)
}

# Example usage:
# confusion_matrix <- confusionMatrix(MUVRObject)

# Histogram plotting for metabolomics data
metabo_histograms <- function(metabo_mat, metabo_names, plot_name, output_path = "results/") {
  pdf(file.path(output_path, paste0(plot_name, "_Histograms.pdf")), onefile = TRUE)
  for (name in 1:length(metabo_names)) {
    hist(metabo_mat[, name], main = paste("Histogram of", metabo_names[name]),
         xlab = paste("Intensity of", metabo_names[name]))
  }
  dev.off()
}

# Example usage:
# metabo_histograms(metabo_matrix, c("Metabolite1", "Metabolite2"), "my_analysis")

# Select specific columns from a DataFrame
select_columns <- function(df, cols) {
  df %>%
    select(all_of(cols))
}

# Example usage:
# selected_df <- select_columns(my_data, c("col1", "col2"))

# Rename columns in a DataFrame using base R
rename_columns_base <- function(df, new_names) {
  names(df)[names(df) %in% names(new_names)] <- new_names[names(df)[names(df) %in% names(new_names)]]
  df
}

# Example usage:
# renamed_df <- rename_columns_base(df, new_colnames)

# Flatten a correlation matrix
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor = cormat[ut],
    p = pmat[ut]
  )
}

# Example usage:
# flat_corr <- flattenCorrMatrix(cor_matrix, pval_matrix)

# Clean column names by removing semicolons
clean_colnames <- function(df) {
  colnames <- colnames(df)
  new_colnames <- sapply(colnames, function(colname) {
    if (grepl(";", colname)) {
      sub(";.*", "", colname)
    } else {
      colname
    }
  })
  colnames(df) <- new_colnames
  return(df)
}

# Example usage:
# cleaned_df <- clean_colnames(df)

# Preprocess data: normalize and optionally log-transform selected columns
preprocess_data <- function(data, start_col, end_col, transform_type = "normalized") {
  selected_data <- data[, start_col:end_col]
  process <- preProcess(selected_data, method = c("range"))
  norm_scale <- predict(process, selected_data)
  
  if (transform_type == "log") {
    norm_scale <- norm_scale %>%
      mutate(across(everything(), ~ log2(. + 1)))
  }
  
  data[, start_col:end_col] <- norm_scale
  return(data)
}

# Example usage:
# normalized_data <- preprocess_data(df, 2, 5, transform_type = "log")

# Find duplicated characters in a vector
find_duplicated_chars <- function(char_vector) {
  counts <- table(char_vector)
  duplicated_chars <- names(counts[counts > 1])
  return(duplicated_chars)
}

# Example usage:
# duplicated_characters <- find_duplicated_chars(c("a", "b", "b", "c"))

# Unadjusted logistic regression model analysis
unadjusted_model_analysis <- function(data, biomarkers, group_col) {
  results_table <- data.frame()
  n <- nrow(data)
  
  for (biomarker in biomarkers) {
    formula <- reformulate(biomarker, response = group_col)
    model <- glm(formula, family = binomial, data = data)
    model_summary <- summary(model)$coefficients
    
    conf_int <- confint(model)
    model_par <- as.data.frame(model_summary)
    colnames(model_par) <- c("Coefficient", "SE", "z", "p")
    
    model_par$CI_low <- conf_int[, 1]
    model_par$CI_upper <- conf_int[, 2]
    model_par$OR <- exp(model_par$Coefficient)
    model_par$CI_low_OR <- exp(model_par$CI_low)
    model_par$CI_upper_OR <- exp(model_par$CI_upper)
    
    results_table <- rbind(results_table, model_par)
  }
  
  return(results_table)
}

# Example usage:
# results <- unadjusted_model_analysis(data, biomarkers, "Group")

# Adjusted logistic regression model analysis
adjusted_model_analysis <- function(data, biomarkers, covariates, group_col) {
  results_table <- data.frame()
  n <- nrow(data)
  
  for (biomarker in biomarkers) {
    formula <- reformulate(c(biomarker, covariates), response = group_col)
    model <- glm(formula, family = binomial, data = data)
    model_summary <- summary(model)$coefficients
    
    conf_int <- confint(model)
    model_par <- as.data.frame(model_summary)
    colnames(model_par) <- c("Coefficient", "SE", "z", "p")
    
    model_par$CI_low <- conf_int[, 1]
    model_par$CI_upper <- conf_int[, 2]
    model_par$OR <- exp(model_par$Coefficient)
    model_par$CI_low_OR <- exp(model_par$CI_low)
    model_par$CI_upper_OR <- exp(model_par$CI_upper)
    
    results_table <- rbind(results_table, model_par)
  }
  
  return(results_table)
}

# Example usage:
# adjusted_results <- adjusted_model_analysis(data, biomarkers, covariates, "Group")

# Combined GLM, fold change, and ROC curve analysis
glm_model_analysis <- function(data, biomarkers, group_col, covariates = NULL, group1 = NULL, group2 = NULL) {
  results_table <- data.frame()
  roc_list <- list()
  n <- nrow(data)
  
  if (!is.null(group1) && !is.null(group2)) {
    log2fc_data <- calculate_fold_change(data, group_col, group1, group2)
  }
  
  for (biomarker in biomarkers) {
    formula <- if (is.null(covariates)) {
      reformulate(biomarker, response = group_col)
    } else {
      reformulate(c(biomarker, covariates), response = group_col)
    }
    
    model <- glm(formula, data = data, family = binomial)
    model_summary <- summary(model)$coefficients
    
    # Extract and process model results (omitted for brevity)
    # Append results to results_table
    
    # ROC calculation (omitted for brevity)
  }
  
  return(list(results_table = results_table, roc_list = roc_list))
}

# Example usage:
# combined_results <- glm_model_analysis(data, biomarkers, "Group", covariates = c("Age", "Gender"))

# Fold change calculation between two groups
calculate_fold_change <- function(data, group_col, group1, group2) {
  group1_data <- data %>% filter(.data[[group_col]] == group1)
  group2_data <- data %>% filter(.data[[group_col]] == group2)
  
  geo_mean <- function(x) exp(mean(log(x[x > 0]), na.rm = TRUE))
  
  fold_change <- apply(group1_data, 2, geo_mean) / apply(group2_data, 2, geo_mean)
  log2fc <- log2(fold_change)
  return(log2fc)
}

# Example usage:
# fold_change_result <- calculate_fold_change(data, "Group", "Group1", "Group2")
