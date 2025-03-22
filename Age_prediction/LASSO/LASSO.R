graphics.off()  
rm(list = ls()) 
setwd("pathway")
library(glmnet)
asv_data <- read.csv(file = "Binarized_data.csv", header = T, row.names = 1)
asv_data <- na.omit(asv_data)
asv_data <- t(asv_data)
age <- read.csv(file = "Sample_information.csv", header = T, row.names = 1)
age <-as.matrix(age[,2])
age <- as.numeric(age)
x <- asv_data
y <- age

#------------------------------------------------------------------
num_repeats <- 100
selected_features <- matrix(0, nrow(x), ncol(x))
selected_features <- rep(0, ncol(x)) 
for (i in 1:num_repeats) {
    set.seed(i) 
    lasso_fit <- cv.glmnet(x, y, family = "gaussian", alpha = 1, type.measure = "mse", nlambda = 100, nfolds=10) 
    selected <- coef(lasso_fit, s = "lambda.min")[-1] != 0  
    selected_features <- selected_features + selected  
}


threshold <- 0.5 
selected_frequencies <- selected_features / num_repeats
stable_features <- colnames(x)[selected_frequencies >= threshold]
print(stable_features)
write.csv(stable_features, file = "Result.csv", row.names = FALSE)
#-----------------------------------------------------------------------------------