############################# LOAD PACKAGES ###################################

library(tidyverse)  # collection of R packages designed for data science. It is 
# designed to make it easy to install and load multiple packages needed for data 
# analysis. The packages included are: ggplot2, dplyr, tidyr, readr, purrr,
# tibble, stringr and forcats
library(mice)       # Provides multiple imputation for missing data
library(caret)      # Provides a set of functions for training and evaluating 
#predictive models
library(themis)     # Provides functions for dealing with imbalanced data
library(randomForest)  # Provides functions for building random forests, an 
#ensemble learning method for classification, regression and other tasks
library(ggplot2)    # Provides a set of functions for data visualization
library(mltools)    # Provides a set of functions for machine learning tasks
library(ggtext)     # Provides a set of functions for using formatted text in 
#ggplot2 graphics
library(pROC)     # Provides functions for analyzing and visualizing ROC curves
library(plot3D)     # Provides functions for creating various types of 3D plots, 
# including surface plots, scatter plots, and wireframe plots. 


# Sets the seed for reproducibility
set.seed(123)


############################### DATA LOADING ###################################

# Set working directory and read in data from CSV file
setwd("C:/Users/Nera/Desktop/TFM_MUBBC_UAM")
df <- read.csv("ASC_CLASIFICACION.csv", sep = ";", header = TRUE)


############################ DATA PREPROCESSING ################################

# Convert columns 2-9 to numeric data type 
df <- df%>%mutate_at(vars(2:9), as.numeric)

# Recode factor variables for easier interpretation
df <- df %>%
  mutate(Exome = recode_factor(Exome, "0" = "Negative", "1" = "Positive")) %>%
  mutate(Gender = recode_factor(Gender, "1" = "Male", "2" = "Female")) %>%
  mutate(Epilepsy = recode_factor(Epilepsy, "0" = "No", "1" = "Probably", 
                                  "2" = "Yes")) %>%
  mutate(Regression = recode_factor(Regression, "0" = "No", "1" = "Yes")) %>%
  mutate(Delivery_problems = recode_factor(Delivery_problems, "0" = "No", 
                                           "1" = "Yes")) %>%
  mutate(Age_symptoms = recode_factor(Age_symptoms, "0" = "<1", "1" = "1-3", 
                                      "2" = ">3"))

# Remove rows with missing values in specific columns
df_no_ID <- df[,c(2:10)]
df_no_ID <- df_no_ID[complete.cases(df_no_ID[,
                                    c("Age_F_Birth", "Age_symptoms")]), ]

# Create a frequency table and plot a bar chart of the poportion of Exome values
table(df_no_ID$Exome)
my_colors <- c("#600000", "#004400")
freq_table <- table(df_no_ID$Exome)
barplot(freq_table, 
        xlab = "Exome", 
        ylab = "Frequency", 
        col = my_colors, 
        main = "",
        ylim = c(0, max(freq_table)*1.2))

title(main = c("Proportion positive /", "negative cases"), 
      line = 2, 
      cex.main = 1.3)

# Add the percentage values on top of each bar
for (i in 1:length(freq_table)) {
  text(x = i, 
       y = freq_table[i], 
       labels = paste0(round(prop.table(freq_table) * 100)[i], "%"),
       col = "black", 
       cex = 0.9, 
       pos = 3)
}

# Calculate rows with missing values in specific columns
mean(is.na(df_no_ID$Age_walking)) * 100
mean(is.na(df_no_ID$Regression)) * 100
mean(is.na(df_no_ID$Epilepsy)) * 100


############################## DATA IMPUTATION #################################

# Run multiple imputation on the data using the mice package
mice_object <- mice(df_no_ID[,-9], m = 5, maxit = 50, method = c("", "", "",
                                         "polyreg", "pmm", "logreg", "", ""))

# Retrieve the imputed data from the mice object
imputed_data <- complete(mice_object)

# Add the Exome column back to the imputed data
imputed_data <- cbind(imputed_data, df_no_ID$Exome)

# Rename the Exome column to match the original dataset
imputed_data <- dplyr::rename(imputed_data, "Exome" = "df_no_ID$Exome")

# Calculate density distribution of the original and imputed data
dens_orig <- density(df_no_ID$Age_walking, na.rm = TRUE)
dens_imp <- density(imputed_data$Age_walking)

# Create data frames for the original and imputed data
df_orig <- data.frame(x = dens_orig$x, y = dens_orig$y, type = "Original data")
df_imp <- data.frame(x = dens_imp$x, y = dens_imp$y, type = "Imputed data")

# Plot the density distributions using ggplot
ggplot(data = data.frame(x = c(dens_orig$x, dens_imp$x), y = c(dens_orig$y, 
       dens_imp$y), type = rep(c("Original data", "Imputed data"), 
       each = length(dens_orig$x))), aes(x = x, y = y, color = type)) +
  geom_line(linewidth = 0.5) +
  labs(x = "Age in months", y = "Density distribution", title = "Age walking") +
  scale_color_manual(values = c("blue", "red")) +
  theme(legend.position = "right") +
  theme_minimal()



############################## DATA BALANCING ##################################

# Create training and testing partitions to evaluate the balancing models
trainIndex_BALANCE <- createDataPartition(imputed_data$Exome, p = 0.7, 
                                          list = FALSE)
train_balance <- imputed_data[trainIndex_BALANCE, ]
test_balance <- imputed_data[-trainIndex_BALANCE, ]

# Set up the train control parameters
ctrl <- trainControl(method = "cv", number = 5, 
                     summaryFunction = twoClassSummary, classProbs = TRUE)

# Define the sampling methods to use
sampling_methods <- c("none", "under", "over", "smotenc")

# Create an empty list to store the results
results <- list()

# Loop through each sampling method and fit a random forest model
for (method in sampling_methods) {
  
  # Use the "none" method as the baseline (no method)
  if (method == "none") {
    model <- randomForest(Exome ~ ., data = train_balance, 
                          trControl = ctrl, importance = TRUE)
  } 
  # Use the undersampling method to balance the classes
  else if (method == "down") {
    train_balanced <- caret::downSample(x = train_balance[, -9], 
                                      y = train_balance$Exome, yname = "Exome")
    model <- randomForest(Exome ~ ., data = train_balanced, 
                          trControl = ctrl, importance = TRUE)
  } 
  # Use the oversampling method to balance the classes
  else if (method == "over") {
    train_balanced <- caret::upSample(x = train_balance[, -9], 
                                      y = train_balance$Exome, yname = "Exome")
    model <- randomForest(Exome ~ ., data = train_balanced, 
                          trControl = ctrl, importance = TRUE)
  } 
  # Use the SMOTE-NC method to balance the classes 
  else if (method == "smotenc") {
    train_balanced <- smotenc(train_balance, var = "Exome", k = 5, 
                              over_ratio = 1)
    model <- randomForest(Exome ~ ., data = train_balanced, 
                          trControl = ctrl, importance = TRUE)
  }
  
  # Store the results in the list
  results[[method]] <- model
}

# Create an empty list to store the evaluation results for each sampling method
evaluations <- list()

# For each sampling method, evaluate the model performance on the test data
for (method in sampling_methods) {
  
  # Get the model trained using the current sampling method
  model <- results[[method]]
  
  # Generate predicted probabilities for the test data
  pred_prob <- predict(model, newdata = test_balance, type = "prob")[, 2]
  
  # Convert predicted probabilities to predicted class levels (Positive or 
  # Negative)
  pred_class <- ifelse(pred_prob < 0.5, "Negative", "Positive")
  
  # Generate the confusion matrix using the predicted class levels and true 
  # classes
  cm <- confusionMatrix(data = as.factor(pred_class), 
                        reference = test_balance$Exome, mode = "prec_recall")
  
  # Store the evaluation metrics for the current sampling method in a dataframe
  evaluations[[method]] <- data.frame(
    Method = method,
    AUC = roc(test_balance$Exome, pred_prob)$auc,
    Precision = cm$byClass["Precision"],
    Sensitivity = cm$byClass["Sensitivity"],
    F1 = cm$byClass["F1"]
  )
  
}

# Combine the evaluation results for all sampling methods into a single data
# frame
evaluations_df <- do.call(rbind, evaluations)

# Identify the best sampling method based on the average performance across all 
# evaluation metrics
best_method <- evaluations_df[which.max(rowMeans(evaluations_df[, -1])),
                              "Method"]

# Reshape the evaluation data frame into a long formato for plotting
evaluations_long <- reshape2::melt(evaluations, id.vars = "Method")

# Create a grouped var plot of the evaluation metrics for each sampling method
ggplot(evaluations_long, aes(x = Method, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  xlab("Sampling Method") +
  ylab("Value") +
  ggtitle("Comparison of Sampling Methods") +
  scale_fill_discrete(name = "Evaluation Metric") 

# Select the sampling method with the best performance
if (best_method == "over") {
  # Upsample the minority class using caret's upSample function
  balanced_df <- caret::upSample(x = imputed_data[, -9], y = imputed_data$Exome, 
                                 yname = "Exome")
} else if (best_method == "under") {
  # Downsample the majority class using caret's downSample function
  balanced_df <- caret::downSample(x = imputed_data[, -9], 
                                   y = imputed_data$Exome)
} else if (best_method == "smotenc") {
  # Use SMOTENC algorithm to synthetically oversample the minority class
  balanced_df <- smotenc(imputed_data, var = "Exome", k = 5, over_ratio = 1)
} else if (best_method == "none") {
  # If the best method is "none", use the original imputed data without any 
  # sampling
  balanced_df <- imputed_data
}


############################### RF ALGORITHM ###################################

# Split the balanced data into training and test sets for the RF algorithm
trainIndex_RF <- createDataPartition(balanced_df$Exome, p = 0.7, list = FALSE)
train_RF <- balanced_df[trainIndex_RF, ]
test_RF <- balanced_df[-trainIndex_RF, ]

# Set the number of folds for cross-validation
k <- 5

# Define the values of maxnodes and ntree to test
maxnodes_values <- c(3, 4, 5, 6, 7)
ntree_values <- c(250, 500, 1000, 1500, 2000)
bestmtry <- c(1, 2, 3, 4, 5)

# Initialize matrix to store accuracies for each combination of hyperparameters
accuracy_matrix <- matrix(nrow=length(maxnodes_values)*length(ntree_values)*
                            length(bestmtry), ncol=4, dimnames=list(NULL, 
                            c("Maxnodes", "NTree", "Bestmtry", "Accuracy")))

# Loop over each combination of hyperparameters
index <- 1
for (maxnodes in maxnodes_values) {
  for (ntree in ntree_values) {
    for (bestmtry in bestmtry) {
      
      # Create k folds for cross-validation
      folds <- createFolds(train_RF$Exome, k = k)
      
      # Initialize variable to keep track of the average accuracy across folds
      avg_accuracy <- 0
      
      # Loop over each fold of cross-validation
      for (i in 1:k) {
        
        # Split the data into training and validation sets for this fold
        train_indices <- unlist(folds[-i])
        valid_indices <- folds[[i]]
        train_cv <- train_RF[train_indices, ]
        valid_cv <- train_RF[valid_indices, ]
        
        # Train random forest model on the training set
        rf_model <- randomForest(Exome ~ ., data = train_cv, importance = TRUE, 
                            ntree = ntree, maxnodes = maxnodes, mtry = bestmtry)
        
        # Make predictions on the validation set using the trained model
        rf_pred <- predict(rf_model, newdata = valid_cv)
        
        # Calculate the confusion matrix and accuracy score for this fold
        confusion_matrix <- table(valid_cv$Exome, rf_pred)
        accuracy <- sum(diag(confusion_matrix))/sum(confusion_matrix)
        
        # Add the accuracy to the running total for this combination of 
        # hyperparameters
        avg_accuracy <- avg_accuracy + accuracy
        
      }  # End of cross-validation loop
      
      # Calculate the average accuracy across all folds for this combination of 
      # hyperparameters
      avg_accuracy <- avg_accuracy/k
      
      # Update the accuracy matrix with the hyperparameters and accuracy
      accuracy_matrix[index, 1] <- maxnodes
      accuracy_matrix[index, 2] <- ntree
      accuracy_matrix[index, 3] <- bestmtry
      accuracy_matrix[index, 4] <- avg_accuracy
      
      # Increment the index for the next combination of hyperparameters
      index <- index + 1
      
    } 
  }
} 


# Find the row of the accuracy matrix with the highest accuracy
best_row <- which.max(accuracy_matrix[, 4])

# Extract the best hyperparameters from the accuracy matrix
best_maxnodes <- accuracy_matrix[best_row, 1]
best_ntree <- accuracy_matrix[best_row, 2]
best_mtry <- accuracy_matrix[best_row, 3]

# Train random forest model on the full training set using the best 
# hyperparameters
rf_model <- randomForest(Exome ~ ., data = train_RF, importance = TRUE, 
                         ntree = best_ntree, maxnodes = best_maxnodes,
                         mtry = best_mtry)

# Make predictions on the test set using the trained model
rf_pred <- predict(rf_model, newdata = test_RF)

# Calculate the confusion matrix and accuracy of the predictions
confusion_matrix <- table(test_RF$Exome, rf_pred)
accuracy <- sum(diag(confusion_matrix))/sum(confusion_matrix)

# Compute sensitivity and F1 score
sensitivity_val <- sensitivity(confusion_matrix)
f1_score_val <- F_meas(confusion_matrix)

# Compute ROC curve and ROC AUC score
roc_obj <- roc(test_RF$Exome, predict(rf_model, newdata = test_RF, 
                                      type = "prob")[,2])
roc_auc <- auc(roc_obj)

# Print the confusion matrix, accuracy, sensitivity, F1 score, and ROC AUC 
# score on the test set
cat("\nAccuracy on test set:", accuracy)
cat("\nSensitivity on test set:", sensitivity_val)
cat("\nF1 score on test set:", f1_score_val)
cat("\nROC AUC score on test set:", roc_auc)
cat("\n\nConfusion Matrix:")
print(confusion_matrix)

# Plot ROC curve
plot(roc_obj, main = "ROC curve")

# Create a data frame from the accuracy matrix and remove rows with NA values
accuracy_df <- as.data.frame(accuracy_matrix)
accuracy_df <- accuracy_df[complete.cases(accuracy_df),]

# Plot a 3D scatter plot withe Maxnodes and NTree as x and t axes
# Accuracy as z axis, and color representing Accuracy
scatter3D(x = accuracy_df$Maxnodes, y = accuracy_df$NTree, 
          z = accuracy_df$Accuracy,
          phi = 15, theta = 30, colvar = accuracy_df$Accuracy,
          pch = 16, cex = 2, type = "h", ticktype = "detailed", 
          xlab = "Maxnodes", ylab = "", zlab = "",
          main = "Accuracy vs Maxnodes and NTree")

# Add axis labels for Ntree and Accuracy
dims <- par("usr")
x1 <- dims[1]+ 0.75*diff(dims[1:2])
y1 <- dims[3]+ 0.17*diff(dims[3:4])
x2 <- dims[1]+ 0.49*diff(dims[1:2])
y2 <- dims[3]- 0.04*diff(dims[3:4])
text(x1,y1,expression(NTree),srt=60)
text(y2,x2,expression(Accuracy),srt=90)

# Calculate and print the variable importance of the model
var_imp <- importance(rf_model)
var_imp <- var_imp[order(var_imp[,4], decreasing = TRUE), ]
var_imp <- as.data.frame(var_imp)
print(var_imp)

# Plot the variable importance scores for each variable
varImpPlot(rf_model, main = "Variable importance of model", sort = T)

#Plot the tree for the model
reprtree::plot.getTree(rf_model)

# Define the function to calculate the percentage of positive cases
calc_percentage <- function(data) {
  # Split the patients to 7 groups based on the values for each variable on 
  # each node
  split1 <- ifelse(data$Age_walking <= 12,
                   ifelse(data$Age_walking <= 10.7,
                          ifelse(data$Age_F_Birth <= 35, "Group1", "Group2"),
                          ifelse(data$Age_M_Birth <= 29.5, "Group3", "Group4")),
                   ifelse(data$Age_F_Birth <= 39.7, 
                          ifelse(data$Epilepsy == "Yes", "Group5", "Group6"),
                          "Group7"))
  
  # Create a new data frame with the groups and Exome variable
  groups <- data.frame(Group = split1, Exome = data$Exome)
  
  # Calculate the percentage of positive cases for each group
  percentage_positive <- aggregate(Exome ~ Group, data = groups,
                                   FUN = function(x) 100 * sum(x == "Positive") 
                                   / length(x))
  
  # Use complete to ensure that all groups are included in each bootstrap sample
  percentage_positive <- complete(percentage_positive, Group = c("Group1", 
                  "Group2", "Group3", "Group4", "Group5", "Group6", "Group7"))
  
  return(percentage_positive$Exome)
}

# Check the patient distribution in each group
table(split1)

# Set the bootstraping value to 1000
n_boot <- 1000

# Perform the bootstrap resampling to estimate confidence intervals for 
# percentage of positive results
boot_results <- replicate(n_boot, {
  boot_data <- imputed_data[sample(nrow(imputed_data), replace = TRUE), ]
  calc_percentage(boot_data)
})

# Transpose the matrix so each row represents the percentage of positive results
# for a single group across all bootstrap samples
boot_results <- t(boot_results)

# Replace NA values in the secon column with 0, as there are missing values in
# some of the bootstrap samples
boot_results[, 2] <- replace(boot_results[, 2], is.na(boot_results[, 2]), 0)

# Calculate mean and confidence interval for each column
boot_mean <- apply(boot_results, 2, mean, na.rm = TRUE)
boot_ci <- t(apply(boot_results, 2, function(x) quantile(x, 
                                      probs = c(0.025, 0.975), na.rm = TRUE)))

# Combine mean and confidence interval into a data frame
boot_summary <- data.frame(Group = c("Group1", "Group2", "Group3", "Group4", 
                                     "Group5", "Group6", "Group7"),
                           Mean = boot_mean,
                           Lower_CI = boot_ci[,1],
                           Upper_CI = boot_ci[,2])
