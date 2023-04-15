# Load required libraries
library(caret)
library(ROSE)
library(DMwR)
library(randomForest)
library(dplyr)
library(ggplot2)
library(pROC)
library(plyr)
library(ROCR)
library(smotefamily)
library(mltools)


##########################LOAD DDBB############################################
setwd("C:/Users/Nera/Desktop/TFM_MUBBC/6. PRUEBAS DECISION TREE")

df <- read.csv("data.csv", sep = ";", header = FALSE)
col_names <- c("Codigounico", "Edad", "Sexo", "REGRESION",
               "EDADINSINT", "SINTDIGESTIVOS", "EPILEPSIA",
               "ADOSss", "CI", "exoma")

colnames(df) <- col_names

df <- df%>%mutate_at(vars(2:9), as.numeric)
df$exoma <- as.factor(df$exoma)


#Comprobamos la división "equitativa" entre los pacientes con y sin resultados
#positivos y negativos del exoma
table(df$exoma)

#
df <- df[,c(2:10)]


##########################PREPROCESS DDBB#####################################

#df <- preProcess(df[,-9], method = "medianImpute", na.remove = FALSE, k = 5)



##########################BALANCE DDBB######################################
# Set seed for reproducibility
set.seed(677)

# Create a training-validation partition
trainIndex_BALANCE <- createDataPartition(df$exoma, p = 0.7, list = FALSE)
train_BD <- df[trainIndex_BALANCE, ]
test_BD <- df[-trainIndex_BALANCE, ]

# Define the control object for train function
ctrl <- trainControl(method = "cv", number = 5, summaryFunction = twoClassSummary, classProbs = TRUE)

# Define a list of sampling methods
sampling_methods <- c("none", "under", "over", "adasyn", "smote")

# Train random forest models using different balance methods
results <- list()
for (method in sampling_methods) {
  if (method == "none") {
    model <- caret::train(exoma ~ ., data = train_BD, method = "rf", trControl = ctrl, importance = TRUE)
  } else if (method == "down") {
    train_balanced <- caret::downSample(x = train_BD[, -9], y = train$exoma)
    train_balanced$exoma <- as.factor(train_balanced$Class)
    model <- caret::train(exoma ~ ., data = train_balanced, method = "rf", trControl = ctrl, importance = TRUE)
  } else if (method == "over") {
    train_balanced <- caret::upSample(x = train_BD[, -9], y = train_BD$exoma, yname = "exoma")
    model <- caret::train(exoma ~ ., data = train_balanced, method = "rf", trControl = ctrl, importance = TRUE)
  } else if (method == "adasyn") {
    train_balanced <- ADAS(train_BD[,-9], train_BD$exoma, K=5)
    model <- caret::train(class ~ ., data = train_balanced$data, method = "rf", trControl = ctrl, importance = TRUE)
  } else if (method == "smote") {
    train_balanced <- DMwR::SMOTE(exoma ~., train_BD, perc.over = 125, perc.under = 100)
    model <- caret::train(exoma ~ ., data = train_balanced, method = "rf", trControl = ctrl, importance = TRUE)
  }
  results[[method]] <- model
}

# Evaluate models on the validation set
evaluations <- list()
for (method in sampling_methods) {
  model <- results[[method]]
  pred_prob <- predict(model, newdata = test_BD, type = "prob")[, 2]
  pred_class <- ifelse(pred_prob > 0.5, "positive", "negative")
  evaluations[[method]] <- data.frame(
    Method = method,
    ROC = roc(test_BD$exoma, pred_prob)$auc,
    Precision = precision(as.factor(pred_class), test_BD$exoma, positive = "exoma"),
    Sensitivity = recall(as.factor(pred_class), test_BD$exoma, positive = "exoma"),
    F1 = caret::F_meas(as.factor(pred_class), test_BD$exoma, positive = "exoma"),
    MCC = mcc(as.factor(pred_class), test_BD$exoma)
  )
}


# Combine evaluations into a single data frame
evaluations_df <- do.call(rbind, evaluations)
best_method <- evaluations_df[which.max(rowMeans(evaluations_df[, -1])), "Method"]

# Print evaluation results
#print(evaluations_df)


#Plot the results
evaluations_long <- reshape2::melt(evaluations, id.vars = "Method")

ggplot(evaluations_long, aes(x = Method, y = value, fill = variable)) +
       geom_bar(stat = "identity", position = "dodge") +
       xlab("Sampling Method") +
       ylab("Value") +
       ggtitle("Comparison of Sampling Methods") +
       scale_fill_discrete(name = "Evaluation Metric")



# Create an empty plot
plot(NULL, xlim=c(0,1), ylim=c(0,1), xlab="False Positive Rate", ylab="True Positive Rate", main="Receiver Operating Characteristic (ROC) Curve")

# Add curves for each sampling method
for (method in sampling_methods) {
  model <- results[[method]]
  pred_prob <- predict(model, newdata = test_BD, type = "prob")[, 2]
  roc_obj <- roc(test_BD$exoma, pred_prob, plot=TRUE, col=ifelse(method=="none", "black", "red"))
  auc <- round(auc(roc_obj), 3)
  legend("bottomright", paste(method, " (AUC=", auc, ")", sep=""), col=ifelse(method=="none", "black", "red"), lty=1)
}


if (best_method == "adasyn") {
  # Apply ROSE to balance the data
  balanced_df <- ROSE(exoma ~ ., data = df, seed = 123)$data
} else if (best_method == "over") {
  # Apply undersampling to balance the data
  balanced_df <- caret::upSample(x = df[, -9], y = df$exoma, yname = "exoma")
} else if (best_method == "under") {
  # Apply oversampling to balance the data
  balanced_df <- caret::downSample(x = df[, -9], y = df$exoma)
} else if (best_method == "smote") {
  # Apply SMOTE to balance the data
  balanced_df <- SMOTE(exoma ~., df, perc.over = 125, perc.under = 100)
} else if (best_method == "none") {
  # No balancing method specified
  balanced_df <- df
}


##########################RANDOM FOREST APLICATION##############################

trainIndex_RF <- createDataPartition(balanced_df$exoma, p = 0.7, list = FALSE)
train_RF <- balanced_df[trainIndex_RF, ]
test_RF <- balanced_df[-trainIndex_RF, ]



# Set the number of folds for cross-validation
k <- 10

# Create an empty vector to store the accuracies for each fold
cv_accuracies <- rep(0, k)

# Create the folds using the createFolds function from the caret package
folds <- createFolds(train_RF$exoma, k = k)

# Perform k-fold cross-validation
for (i in 1:k) {
  
  # Split the data into train and validation sets for this fold
  train_indices <- unlist(folds[-i])
  valid_indices <- folds[[i]]
  train_cv <- train_RF[train_indices, ]
  valid_cv <- train_RF[valid_indices, ]
  
  # Train random forest model on the train_cv data
  rf_model <- randomForest(exoma ~ ., data = train_cv, importance = TRUE, proximity = TRUE, ntree = 500)
  
  # Make predictions on the validation set
  rf_pred <- predict(rf_model, newdata = valid_cv)
  
  # Calculate accuracy for this fold and store it in the cv_accuracies vector
  confusion_matrix <- table(valid_cv$exoma, rf_pred)
  accuracy <- sum(diag(confusion_matrix))/sum(confusion_matrix)
  cv_accuracies[i] <- accuracy
}

# Calculate the mean and standard deviation of the accuracies over all folds
mean_accuracy <- mean(cv_accuracies)
sd_accuracy <- sd(cv_accuracies)

# Train the final random forest model on the full train_RF data
rf_model <- randomForest(exoma ~ ., data = train_RF, importance = TRUE, 
                         proximity = TRUE, ntree = 1000, mtry = 5, maxnodes = 7)

# Make predictions on the test set
rf_pred <- predict(rf_model, newdata = test_RF)

# Print the confusion matrix, accuracy, mean accuracy, and standard deviation
confusion_matrix <- table(test_RF$exoma, rf_pred)
accuracy <- sum(diag(confusion_matrix))/sum(confusion_matrix)
cat("Confusion Matrix:\n")
print(confusion_matrix)
cat("\nAccuracy on Test Set:", accuracy)
cat("\nMean Accuracy over", k, "Folds:", mean_accuracy)
cat("\nStandard Deviation of Accuracies over", k, "Folds:", sd_accuracy)


plot.getTree(rf_model)



