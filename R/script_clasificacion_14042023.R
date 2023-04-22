# Load required libraries
library(caret)
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
setwd("C:/Users/Nera/Desktop/TFM_MUBBC_UAM")
df <- read.csv("BBDD_ASC_CLASIFICACION.csv", sep = ";", header = TRUE)


#Cambio a numerico para que limpie cualquier texto que hay en la BBDD
df <- df%>%mutate_at(vars(2:9), as.numeric)




df <- df %>%
  mutate(Exome = recode_factor(Exome, "0" = "Negative", "1" = "Positive")) %>%
  mutate(Gender = recode_factor(Gender, "1" = "Male", "2" = "Female")) %>%
  mutate(Epilepsy = recode_factor(Epilepsy, "0" = "No", "1" = "Probably", "2" = "Yes")) %>%
  mutate(Regression = recode_factor(Regression, "0" = "No", "1" = "Yes")) %>%
  mutate(Delivery_problems = recode_factor(Delivery_problems, "0" = "No", "1" = "Yes")) %>%
  mutate(Age_symptoms = recode_factor(Age_symptoms, "0" = "<1", "1" = "1-3", "2" = ">3"))


#Comprobamos la división "equitativa" entre los pacientes con y sin resultados
#positivos y negativos del exoma
table(df$Exome)

#We remove the identifiers 
df_no_ID <- df[,c(2:10)]


##########################PREPROCESS DDBB#####################################

#Eliminamos las filas de las variables que no vamos a imputar, como la edad
#de los padres
df_no_ID <- df_no_ID[complete.cases(df_no_ID[, c("Age_F_Birth", "Age_symptoms")]), ]


#Comprobamos el porcentaje de valores perdidos de las columnas que vamos a imputar:
mean(is.na(df_no_ID$Age_walking)) * 100
mean(is.na(df_no_ID$Regression)) * 100
mean(is.na(df_no_ID$Epilepsy)) * 100

#imputamos los valores perdidos (no tendremos en cuenta que una de las variables
#tiene más del 18% de NAs)
mice_object <- mice(df_no_ID[,-9], m = 5, maxit = 50, meth = c("logreg",
          "pmm", "pmm", "polyreg", "pmm", "logreg", "logreg", "polyreg"))

imputed_data <- complete(mice_object)

imputed_data <- cbind(imputed_data, df_no_ID$Exome)

imputed_data <- dplyr::rename(imputed_data, "Exome" = "df_no_ID$Exome")
##########################BALANCE DDBB######################################
# Set seed for reproducibility
set.seed(677)

# Create a training-validation partition
trainIndex_BALANCE <- createDataPartition(imputed_data$Exome, p = 0.7, list = FALSE)
train_BD <- imputed_data[trainIndex_BALANCE, ]
test_BD <- imputed_data[-trainIndex_BALANCE, ]

# Define the control object for train function
ctrl <- trainControl(method = "cv", number = 5, summaryFunction = twoClassSummary, classProbs = TRUE)

# Define a list of sampling methods
sampling_methods <- c("none", "under", "over", "smotenc")

# Train random forest models using different balance methods
results <- list()
for (method in sampling_methods) {
  if (method == "none") {
    model <- caret::train(Exome ~ ., data = train_BD, method = "rf", trControl = ctrl, importance = TRUE)
  } else if (method == "down") {
    train_balanced <- caret::downSample(x = train_BD[, -9], y = train_BD$Exome)
    train_balanced$Exome <- as.factor(train_balanced$Class)
    model <- caret::train(Exome ~ ., data = train_balanced, method = "rf", trControl = ctrl, importance = TRUE)
  } else if (method == "over") {
    train_balanced <- caret::upSample(x = train_BD[, -9], y = train_BD$Exome, yname = "Exome")
    model <- caret::train(Exome ~ ., data = train_balanced, method = "rf", trControl = ctrl, importance = TRUE)
  } else if (method == "smotenc") {
    train_balanced <- smotenc(train_BD, var = "Exome", k = 5, over_ratio = 1)
    model <- caret::train(Exome ~ ., data = train_balanced, method = "rf", trControl = ctrl, importance = TRUE)
  }
  results[[method]] <- model
}

# Evaluate models on the validation set
evaluations <- list()
for (method in sampling_methods) {
  model <- results[[method]]
  pred_prob <- predict(model, newdata = test_BD, type = "prob")[, 2]
  pred_class <- ifelse(pred_prob > 0.5, "Positive", "Negative")
  evaluations[[method]] <- data.frame(
    Method = method,
    ROC = roc(test_BD$Exome, pred_prob)$auc,
    Precision = precision(as.factor(pred_class), test_BD$Exome, positive = "Exome"),
    Sensitivity = recall(as.factor(pred_class), test_BD$Exome, positive = "Exome"),
    F1 = caret::F_meas(as.factor(pred_class), test_BD$Exome, positive = "Exome"),
    MCC = mcc(as.factor(pred_class), test_BD$Exome)
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


if (best_method == "over") {
  # Apply oversampling to balance the data
  balanced_df <- caret::upSample(x = imputed_data[, -9], y = imputed_data$Exome, yname = "Exome")
} else if (best_method == "under") {
  # Apply undersampling to balance the data
  balanced_df <- caret::downSample(x = imputed_data[, -9], y = imputed_data$Exome)
} else if (best_method == "smotenc") {
  # Apply SMOTENC to balance the data
  balanced_df <- smotenc(imputed_data, var = "Exome", k = 5, over_ratio = 1)
} else if (best_method == "none") {
  # No balancing method specified
  balanced_df <- imputed_data
}


##########################RANDOM FOREST APLICATION##############################

trainIndex_RF <- createDataPartition(balanced_df$Exome, p = 0.7, list = FALSE)
train_RF <- balanced_df[trainIndex_RF, ]
test_RF <- balanced_df[-trainIndex_RF, ]



# Set the number of folds for cross-validation
k <- 10

# Create an empty vector to store the accuracies for each fold
cv_accuracies <- rep(0, k)

# Create the folds using the createFolds function from the caret package
folds <- createFolds(train_RF$Exome, k = k)

# Perform k-fold cross-validation
for (i in 1:k) {
  
  # Split the data into train and validation sets for this fold
  train_indices <- unlist(folds[-i])
  valid_indices <- folds[[i]]
  train_cv <- train_RF[train_indices, ]
  valid_cv <- train_RF[valid_indices, ]
  
  # Train random forest model on the train_cv data
  rf_model <- randomForest(Exome ~ ., data = train_cv, importance = TRUE, proximity = TRUE, ntree = 500)
  
  # Make predictions on the validation set
  rf_pred <- predict(rf_model, newdata = valid_cv)
  
  # Calculate accuracy for this fold and store it in the cv_accuracies vector
  confusion_matrix <- table(valid_cv$Exome, rf_pred)
  accuracy <- sum(diag(confusion_matrix))/sum(confusion_matrix)
  cv_accuracies[i] <- accuracy
}

# Calculate the mean and standard deviation of the accuracies over all folds
mean_accuracy <- mean(cv_accuracies)
sd_accuracy <- sd(cv_accuracies)

# Train the final random forest model on the full train_RF data
rf_model <- randomForest(Exome ~ ., data = train_RF, importance = TRUE, 
                         proximity = TRUE, ntree = 1000, mtry = 5, maxnodes = 7)

# Make predictions on the test set
rf_pred <- predict(rf_model, newdata = test_RF)

# Print the confusion matrix, accuracy, mean accuracy, and standard deviation
confusion_matrix <- table(test_RF$Exome, rf_pred)
accuracy <- sum(diag(confusion_matrix))/sum(confusion_matrix)
cat("Confusion Matrix:\n")
print(confusion_matrix)
cat("\nAccuracy on Test Set:", accuracy)
cat("\nMean Accuracy over", k, "Folds:", mean_accuracy)
cat("\nStandard Deviation of Accuracies over", k, "Folds:", sd_accuracy)


# Get variable importance measures
var_imp <- randomForest::importance(rf_model)

# Sort the variables by importance (in descending order)
var_imp <- var_imp[order(var_imp[,4], decreasing = TRUE), ]

# Print the variable importance table
print(var_imp)



#########################Plot tree#############################################

plot.getTree(rf_model)



