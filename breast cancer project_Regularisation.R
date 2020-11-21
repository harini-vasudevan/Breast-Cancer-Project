options(digits = 3)
library(matrixStats)
library(tidyverse)
library(caret)
library(dslabs)
data(brca)
# brca$y: a vector of sample classifications ("B" = benign or "M" = malignant)
# brca$x: a matrix of numeric features describing properties of the shape and size of cell nuclei extracted from biopsy microscope images
# Dimensions and properties
dim(brca[["x"]])  # 569 samples & 30 predictors
table(brca[["y"]])
mean(brca[["y"]]=="M")  # proportion of the samples that are malignant
which.max(colMeans(brca$x))  # column number that has the highest mean
which.min(colSds(brca$x))  # column number has the lowest standard deviation

#Scaling the matrix
X_centered <- sweep(brca$x, 2, colMeans(brca$x), FUN = "-")
X_scaled<- sweep(X_centered ,2,colSds(brca$x),FUN ="/")
sd(X_scaled[,1])
median(X_scaled[,1])

#Distance
d <- dist(X_scaled)
dist_BtoB <- as.matrix(d)[1, brca$y == "B"]  #the average distance between the first sample, which is benign, and other benign samples
mean(dist_BtoB[2:length(dist_BtoB)])
dist_BtoM <- as.matrix(d)[1, brca$y == "M"]  #the average distance between the first sample and malignant samples
mean(dist_BtoM)

#Heatmap of features
#heatmap of the relationship between features using the scaled matrix
d_features <- dist(t(X_scaled))
heatmap(as.matrix(d_features), labRow = NA,labCol = NA)

#Hierarchical clustering
h <- hclust(d_features)
groups <- cutree(h, k = 5)  #Cut the tree into 5 groups
split(names(groups), groups)
plot(groups)

#PCA: proportion of variance
#principal component analysis of the scaled matrix
pca <- prcomp(X_scaled)
summary(pca)

#PCA: plotting PCs
data.frame(pca$x[,1:2], type = brca$y) %>%
  ggplot(aes(PC1, PC2, color = type)) +
  geom_point()

#boxplot of the first 10 PCs grouped by tumor type
data.frame(type = brca$y, pca$x[,1:10]) %>%
  gather(key = "PC", value = "value", -type) %>%
  ggplot(aes(PC, value, fill = type)) +
  geom_boxplot()

#creating data set

set.seed(1, sample.kind = "Rounding")    # if using R 3.6 or later
test_index <- createDataPartition(brca$y, times = 1, p = 0.2, list = FALSE)
test_x <- X_scaled[test_index,]
test_y <- brca$y[test_index]
train_x <- X_scaled[-test_index,]
train_y <- brca$y[-test_index]

#Training and test sets
#Checking that the training and test sets have similar proportions of benign and malignant tumors.
mean(train_y== "B") # proportion of the training set is benign
mean(test_y== "B")  #proportion of the test set is benign
#K-means Clustering
#a matrix of observations x and a k-means object k - and assigns each row of x to a cluster from k
predict_kmeans <- function(x, k) {
  centers <- k$centers    # extract cluster centers
  # calculate distance to cluster centers
  distances <- sapply(1:nrow(x), function(i){
    apply(centers, 1, function(y) dist(rbind(x[i,], y)))
  })
  max.col(-t(distances))  # select cluster with min distance to center
}

set.seed(3, sample.kind = "Rounding")
k <- kmeans(train_x, centers=2)
kmeans_preds <- ifelse(predict_kmeans(test_x, k) == 1, "B", "M")
mean(kmeans_preds == test_y)

#proportion of benign tumors are correctly identified 
table(test_y,kmeans_preds)
sensitivity(factor(kmeans_preds), test_y, positive = "B")
sensitivity(factor(kmeans_preds), test_y, positive = "M")

#Logistic regression model

train_glm <- train(train_x, train_y,
                   method = "glm")
glm_preds <- predict(train_glm, test_x)
mean(glm_preds == test_y)


#LDA and QDA models
train_lda <- train(train_x, train_y,
                   method = "lda")
lda_preds <- predict(train_lda, test_x)
mean(lda_preds == test_y)

train_qda <- train(train_x, train_y,
                   method = "qda")
qda_preds <- predict(train_qda, test_x)
mean(qda_preds == test_y)

#Loess model

set.seed(5, sample.kind = "Rounding")
modelLookup("gamLoess")

#grid <- expand.grid(span = seq(0.15, 0.65, len = 10), degree = 1)

library(gam)
train_loess <- train(train_x, train_y, 
                     method = "gamLoess")
loess_preds <- predict(train_loess, test_x)
mean(loess_preds == test_y)

#K-nearest neighbors model
set.seed(7, sample.kind = "Rounding")
train_knn <- train(train_x, train_y, 
                   method = "knn",
                   tuneGrid = data.frame(k = seq(3, 21, 2)))
knn_preds <- predict(train_knn, test_x)
mean(knn_preds == test_y)
train_knn$results %>% 
  ggplot(aes(x = k, y = Accuracy)) +
  geom_line() +
  geom_point()
train_knn$bestTune

#Random forest model
set.seed(9, sample.kind = "Rounding")
train_rf <- train(train_x, train_y, 
                     method = "rf",tuneGrid = data.frame(mtry =c(3, 5, 7, 9)),importance = TRUE)
rf_preds <- predict(train_rf, test_x)
mean(rf_preds == test_y)
train_rf$bestTune
imp <- varImp(train_rf)

#Creating an ensemble
model <- c("kmeans_preds","glm_preds","lda_preds","qda_preds","rf_preds","knn_preds","loess_preds")
pred <- sapply(1:7, function(x){
  as.factor(get(model[x]))})
dim(pred)

pred <- as.data.frame(pred)
names(pred) <-c("kmeans_preds","glm_preds","lda_preds","qda_preds","rf_preds","knn_preds","loess_preds")
acc <- colMeans(as.matrix(pred)==test_y)
acc
mean(acc)

ensemble <- cbind(glm = glm_preds == "B", lda = lda_preds == "B", qda = qda_preds == "B", loess = loess_preds == "B", rf = rf_preds == "B", knn = knn_preds == "B", kmeans = kmeans_preds == "B")

ensemble_preds <- ifelse(rowMeans(ensemble) > 0.5, "B", "M")
mean(ensemble_preds == test_y)

models <- c("K means", "Logistic regression", "LDA", "QDA", "Loess", "K nearest neighbors", "Random forest", "Ensemble")
accuracy <- c(mean(kmeans_preds == test_y),
              mean(glm_preds == test_y),
              mean(lda_preds == test_y),
              mean(qda_preds == test_y),
              mean(loess_preds == test_y),
              mean(knn_preds == test_y),
              mean(rf_preds == test_y),
              mean(ensemble_preds == test_y))
data.frame(Model = models, Accuracy = accuracy)