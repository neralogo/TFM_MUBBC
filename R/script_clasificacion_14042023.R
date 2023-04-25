# Load necessary packages

library(dplyr)      # Provides a set of functions for data manipulation and 
#transformation
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
library(pROC)       # Provides functions for analyzing and visualizing ROC curves

# Sets the seed for reproducibility
set.seed(677)


############################### DATA LOADING ###################################

# Set working directory and read in data from CSV file
setwd("C:/Users/Nera/Desktop/TFM_MUBBC_UAM")
df <- read.csv("BBDD_ASC_CLASIFICACION.csv", sep = ";", header = TRUE)


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
df_no_ID <- df_no_ID[complete.cases(df_no_ID[, c("Age_F_Birth", "Age_symptoms")]), ]


# Create a frequency table and plot a bar chart of the poportion of Exome values
table(df_no_ID$Exome)
my_colors <- c("#600000", "#004400")
freq_table <- table(df_no_ID$Exome)
barplot(freq_table, main = "Proportion", xlab = "Exome", ylab = "Frequency",
        col = my_colors)

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

# Create a histogram of the "Age_walking" column in the original data
hist(df_no_ID$Age_walking, main = "Distribution before imputing", 
     xlab = "Age_walking")

#Create a histogram of the "Age_walking" column in the imputed data
hist(imputed_data$Age_walking, main = "Distribution after imputing", 
     xlab = "Age_walking")


############################## DATA BALANCING ##################################

# Create training and testing partitions to evaluate the balancing models
trainIndex_BALANCE <- createDataPartition(imputed_data$Exome, p = 0.7, list = FALSE)
train_BD <- imputed_data[trainIndex_BALANCE, ]
test_BD <- imputed_data[-trainIndex_BALANCE, ]

# Set up the train control parameters
ctrl <- trainControl(method = "cv", number = 5, summaryFunction = twoClassSummary, 
                     classProbs = TRUE)

# Define the sampling methods to use
sampling_methods <- c("none", "under", "over", "smotenc")

# Create an empty list to store the results
results <- list()

# Loop through each sampling method and fit a random forest model
for (method in sampling_methods) {
  
  # Use the "none" method as the baseline (no method)
  if (method == "none") {
    model <- caret::train(Exome ~ ., data = train_BD, method = "rf", 
                          trControl = ctrl, importance = TRUE)
  } 
  # Use the undersampling method to balance the classes
  else if (method == "down") {
    train_balanced <- caret::downSample(x = train_BD[, -9], y = train_BD$Exome, 
                                        yname = "Exome")
    model <- caret::train(Exome ~ ., data = train_balanced, method = "rf", 
                          trControl = ctrl, importance = TRUE)
  } 
  # Use the oversampling method to balance the classes
  else if (method == "over") {
    train_balanced <- caret::upSample(x = train_BD[, -9], y = train_BD$Exome, 
                                      yname = "Exome")
    model <- caret::train(Exome ~ ., data = train_balanced, method = "rf", 
                          trControl = ctrl, importance = TRUE)
  } 
  # Use the SMOTE-NC method to balance the classes 
  else if (method == "smotenc") {
    train_balanced <- smotenc(train_BD, var = "Exome", k = 5, over_ratio = 1)
    model <- caret::train(Exome ~ ., data = train_balanced, method = "rf", 
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
  
  # Generate predicted probabilites for the test data
  pred_prob <- predict(model, newdata = test_BD, type = "prob")[, 2]
  
  # Convertt predicted probabilites to predicted class levels (Positive or Negative)
  pred_class <- ifelse(pred_prob > 0.5, "Positive", "Negative")
  
  # Store the evaluation metrics for the current sampling method in a dataframe
  evaluations[[method]] <- data.frame(
    Method = method,
    ROC = roc(test_BD$Exome, pred_prob)$auc,
    Precision = caret::precision(as.factor(pred_class), test_BD$Exome, 
                                 positive = "Exome"),
    Sensitivity = caret::recall(as.factor(pred_class), test_BD$Exome, 
                                positive = "Exome"),
    F1 = caret::F_meas(as.factor(pred_class), test_BD$Exome, positive = "Exome"),
    MCC = mcc(as.factor(pred_class), test_BD$Exome)
  )
}

# Combine the evaluation results for all sampling methods into a single data frame
evaluations_df <- do.call(rbind, evaluations)

# Identify the best sampling method based on the average performance across all 
# evaluation metrics
best_method <- evaluations_df[which.max(rowMeans(evaluations_df[, -1])), "Method"]

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
  balanced_df <- caret::downSample(x = imputed_data[, -9], y = imputed_data$Exome)
} else if (best_method == "smotenc") {
  # Use SMOTENC algorithm to synthetically oversample the minority class
  balanced_df <- smotenc(imputed_data, var = "Exome", k = 5, over_ratio = 1)
} else if (best_method == "none") {
  # If the best method is "none", use the original imputed data without any sampling
  balanced_df <- imputed_data
}


############################### RF ALGORITHM ###################################

# Split the balanced data into training and test sets for the RF algorithm
trainIndex_RF <- createDataPartition(balanced_df$Exome, p = 0.7, list = FALSE)
train_RF <- balanced_df[trainIndex_RF, ]
test_RF <- balanced_df[-trainIndex_RF, ]

# Set the number of folds for cross-validation
k <- 10

# Initialize an array to store the accuracy scores for each fold of cross-validation
cv_accuracies <- rep(0, k)

# Define the values of maxnodes and ntree to test
maxnodes_values <- c(5, 7, 10)
ntree_values <- c(500, 1000, 1500)

# Loop over each combination of hyperparameters
for (maxnodes in maxnodes_values) {
  for (ntree in ntree_values) {
    
    # Create k folds for cross-validation
    folds <- createFolds(train_RF$Exome, k = k)
    
    # Loop over each fold of cross-validation
    for (i in 1:k) {
      
      # Split the data into training and validation sets for this fold
      train_indices <- unlist(folds[-i])
      valid_indices <- folds[[i]]
      train_cv <- train_RF[train_indices, ]
      valid_cv <- train_RF[valid_indices, ]
      
      # Train random forest model on the training set
      rf_model <- randomForest(Exome ~ ., data = train_cv, importance = TRUE, 
                               proximity = TRUE, ntree = ntree, maxnodes = maxnodes)
      
      # Make predictions on the validation set using the trained model
      rf_pred <- predict(rf_model, newdata = valid_cv)
      
      # Calculate the confusion matrix and accuracy score for this fold
      confusion_matrix <- table(valid_cv$Exome, rf_pred)
      accuracy <- sum(diag(confusion_matrix))/sum(confusion_matrix)
      cv_accuracies[i] <- accuracy
    }
    
    # Calculate the mean and standard deviation of the accuracy scores over all folds 
    # of cross-validation
    mean_accuracy <- mean(cv_accuracies)
    sd_accuracy <- sd(cv_accuracies)
    
    # Train a random forest model
    rf_model <- randomForest(Exome ~ ., data = train_RF, importance = TRUE, 
                             proximity = TRUE, ntree = ntree, maxnodes = maxnodes)
    
    # Make predictions on the test set using the trained model
    rf_pred <- predict(rf_model, newdata = test_RF)
    
    # Calculate the confusion matrix and accuracy of the predictions
    confusion_matrix <- table(test_RF$Exome, rf_pred)
    accuracy <- sum(diag(confusion_matrix))/sum(confusion_matrix)
    
    # Print the confusion matrix and accuracy on the test set, as well as the hyperparameters used
    cat("\nMaxnodes:", maxnodes, "NTree:", ntree)
    cat("\nConfusion Matrix:\n")
    print(confusion_matrix)
    cat("\nAccuracy on Test Set:", accuracy)
    
    # Print the mean accuracy and standard deviation over k folds
    cat("\nMean Accuracy over", k, "Folds:", mean_accuracy)
    cat("\nStandard Deviation of Accuracies over", k, "Folds:", sd_accuracy)
  }
}


# Train a random forest model
rf_model <- randomForest(Exome ~ ., data = train_RF, importance = TRUE, 
                         proximity = TRUE, ntree = 1000, mtry = 5, maxnodes = 7)

# Make predictions on the test set using the trained model
rf_pred <- predict(rf_model, newdata = test_RF)

# Calculate the confusion matrix and accuracy of the predictions
confusion_matrix <- table(test_RF$Exome, rf_pred)
accuracy <- sum(diag(confusion_matrix))/sum(confusion_matrix)

# Print the confusion matrix and accuracy on the test set
cat("Confusion Matrix:\n")
print(confusion_matrix)
cat("\nAccuracy on Test Set:", accuracy)

# Calculate and print the variable importance of the model
var_imp <- randomForest::importance(rf_model)
var_imp <- var_imp[order(var_imp[,4], decreasing = TRUE), ]
var_imp <- as.data.frame(var_imp)
print(var_imp)
