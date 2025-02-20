library(ggplot2)
library(dplyr)
library(caret)

#load data
data <- read.table(("orf_localization/ecoli.data"))
head(data)
data <- as.data.frame(data)
colnames(data) <- c("Sequence ID","mcg", "gvh", "lip", "chg", "aac","alm1", "alm2", "class" )

#adjust data
data$class <- as.factor(data$class)
data$class <- as.numeric(data$class)
data$class <- as.factor(data$class)

#train test split
train_index <- createDataPartition(data$class, p = 0.7, list = FALSE)
train <- data[train_index, ]
test <- data[-train_index, ]

# Create a 10-fold cross validation
control <- trainControl(method = "cv", number = 10)

# Train Random Forest model using caret
rf_model <- train(class ~ mcg + gvh+lip+chg+aac+alm1+alm2,
                  data = train,
                  method = "rf",
                  trControl = control,
                  tuneLength = 5)

#print sumary
print(rf_model)

rf_pred <- predict(rf_model, test)

# confusion matrix
conf_matrix <- confusionMatrix(rf_pred, test$class)
print(conf_matrix)

#use helpere function
generate_conf_matrix(test$class,rf_pred)
