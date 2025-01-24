set.seed(101)

# Read in Training Data
#args <- commandArgs(trailingOnly=TRUE)

load_Korrapati <- function() {
	args <- Sys.glob("Korrapati*.rds")

	Training_5th95th <- c()
	for (file in args) {
		predictors <- readRDS(file)
		predictors$sample <- rep(file, nrow(predictors))
		if (is.null(dim(Training_5th95th))) {
			Training_5th95th <- predictors
		} else {
			Training_5th95th <- rbind(Training_5th95th, predictors)
		}
	}
	Training_5th95th$is.cell <- as.numeric(Training_5th95th$is.cell)
	Training_5th95th$dataset <- "Korrapati"
	return(Training_5th95th)
}

load_MacPar <- function() {
	args <- Sys.glob("MacPar*.rds")

	Training_5th95th <- c()
	for (file in args) {
		predictors <- readRDS(file)
		predictors$sample <- rep(file, nrow(predictors))
		if (is.null(dim(Training_5th95th))) {
			Training_5th95th <- predictors
		} else {
			Training_5th95th <- rbind(Training_5th95th, predictors)
		}
	}
	Training_5th95th$is.cell <- as.numeric(Training_5th95th$is.cell)
	Training_5th95th$dataset <- "MacPar"
	return(Training_5th95th)
}

load_Denishenko <- function() {
	args <- Sys.glob("Denishenko*.rds")

	Training_5th95th <- c()
	for (file in args) {
		predictors <- readRDS(file)
		predictors$sample <- rep(file, nrow(predictors))
		if (is.null(dim(Training_5th95th))) {
			Training_5th95th <- predictors
		} else {
			Training_5th95th <- rbind(Training_5th95th, predictors)
		}
	}
	Training_5th95th$is.cell <- as.numeric(Training_5th95th$is.cell)
	Training_5th95th$dataset <- "Denishenko"
	return(Training_5th95th)
}

preprocessing <- function(data, probs=c(0.05,0.95)) {
	for (score in c("prop.intronic", "pct.mt", "pct.ribo", "gene.length", "n.exons")) {
		# rescale data based on 1st and 99th percentile
		thresholds <- quantile(data[,score], probs=probs)
		data[,score] <- (data[,score]-thresholds[1])/(thresholds[2]-thresholds[1])
	}
	return(data)
}


# Basic Model Set Up 
outcome = "is.cell" # 1 = cell, 0 = nucleus
predictors = c("pct.mt", "pct.ribo", "gene.length", "n.exons", "prop.intronic")
# pct.mt = % of expression from mt genes
# pct.ribo = % of expression from ribo genes
# gene.length = ratio of expr of 3rd quartile genes by length / expr of 1st quartile genes by length
# n.exons = ratio of expr of 3rd quartile genes by num of exons / expr of 1st quartile genes by num of exons

# Test Quantile Scaling
Training_5th95th <- rbind(preprocessing(load_Denishenko()), preprocessing(load_Korrapati()))
Training_1st99th <- rbind(preprocessing(load_Denishenko(), probs=c(0.01, 0.99)), preprocessing(load_Korrapati(), probs=c(0.01, 0.99)))

## Predictor Redundancy ##
Denis_data <- preprocessing(load_Denishenko())
Korra_data <- preprocessing(load_Korrapati())
MacPar_data <- preprocessing(load_MacPar())
pred_cors <- cor(rbind(Denis_data, Korra_data, MacPar_data)[,1:6])
require(pheatmap)
png("Redundancy_predictor_cors.png", width=6, height=6, units="in", res=300)
pheatmap(pred_cors)
dev.off()

### Logistic Regression ###
calc_accuracy <- function(scores=predict(model, newdata=MacPar_data, type="response"), truth=MacPar_data$is.cell, threshold=0.5) {
	if (threshold < 0.5) { # check for threshold errors
		threshold = 1-threshold
		print(paste("Warning: classification threshold must but > 0.5, has been changed to", threshold))
	}
	TP_c <- sum(scores > threshold & truth == max(truth))
	TP_n <- sum(scores <= 1-threshold & truth == min(truth))
	FP_c <- sum(scores > threshold & truth == min(truth))
	FP_n <- sum(scores <= 1-threshold & truth == max(truth))
	UNK_c <- sum(scores < threshold & scores > 1-threshold &  truth == max(truth))
	UNK_n <- sum(scores < threshold & scores > 1-threshold &  truth == min(truth))

	mann.whitney <- wilcox.test(scores[dataset$is.cell==max(dataset$is.cell)], scores[dataset$is.cell!=max(dataset$is.cell)])$statistic
	AUC <-  mann.whitney/prod(table(dataset$is.cell))

	return(list(accuracy=(TP_c+TP_n)/(TP_c+TP_n+FP_c+FP_n), 
		sc_precision=(TP_c)/(TP_c+FP_c), 
		sc_recall=TP_c/(TP_c+FP_n+UNK_c),
		sn_precision=(TP_n)/(TP_n+FP_n),
		sn_recall=TP_n/(TP_n+FP_c+UNK_n), 
		unclassified=(UNK_c+UNK_n)/length(scores),
		AUC=AUC))
}
# Compare Predictors / Combinations
compare_predictors <- function(training_data=Training_5th95th) {
	single_pred_out <- c()
	model <- glm(is.cell~pct.mt, family=binomial(link="logit"), data=training_data[,colnames(training_data) %in% c(outcome, predictors)])
	single_pred_out <- rbind(single_pred_out, c("pct.mt", round(unlist(calc_accuracy(model$fitted.values, truth=training_data$is.cell)), digits=4)))

	model <- glm(is.cell~pct.ribo, family=binomial(link="logit"), data=training_data[,colnames(training_data) %in% c(outcome, predictors)])
	single_pred_out <- rbind(single_pred_out, c("pct.ribo", round(unlist(calc_accuracy(model$fitted.values, truth=training_data$is.cell)), digits=4)))
	model <- glm(is.cell~prop.intronic, family=binomial(link="logit"), data=training_data[,colnames(training_data) %in% c(outcome, predictors)])
	single_pred_out <- rbind(single_pred_out, c("prop.intronic", round(unlist(calc_accuracy(model$fitted.values, truth=training_data$is.cell)), digits=4)))
	model <- glm(is.cell~gene.length, family=binomial(link="logit"), data=training_data[,colnames(training_data) %in% c(outcome, predictors)])
	single_pred_out <- rbind(single_pred_out, c("gene.length", round(unlist(calc_accuracy(model$fitted.values, truth=training_data$is.cell)), digits=4)))
	model <- glm(is.cell~n.exons, family=binomial(link="logit"), data=training_data[,colnames(training_data) %in% c(outcome, predictors)])
	single_pred_out <- rbind(single_pred_out, c("n.exons", round(unlist(calc_accuracy(model$fitted.values, truth=training_data$is.cell)), digits=4)))

	model <- glm(is.cell~pct.mt+pct.ribo+prop.intronic+n.exons+gene.length, family=binomial(link="logit"), data=training_data[,colnames(training_data) %in% c(outcome, predictors)])
	single_pred_out <- rbind(single_pred_out, c("all", round(unlist(calc_accuracy(model$fitted.values, truth=training_data$is.cell)), digits=4)))

	model <- glm(is.cell~pct.mt+pct.ribo+gene.length+n.exons, family=binomial(link="logit"), data=training_data[,colnames(training_data) %in% c(outcome, predictors)])
	single_pred_out <- rbind(single_pred_out, c("no_intronic", round(unlist(calc_accuracy(model$fitted.values, truth=training_data$is.cell)), digits=4)))

	model <- glm(is.cell~pct.mt+pct.ribo, family=binomial(link="logit"), data=training_data[,colnames(training_data) %in% c(outcome, predictors)])
	single_pred_out <- rbind(single_pred_out, c("mt+ribo", round(unlist(calc_accuracy(model$fitted.values, truth=training_data$is.cell)), digits=4)))
	return(single_pred_out)
}
predictors_5th95th_out <- compare_predictors(training_data=Training_5th95th)
predictors_1st99th_out <- compare_predictors(training_data=Training_1st99th)

write.table(predictors_5th95th_out, file="Scale_5th95ht_logistic_regression.out", sep="_")
write.table(predictors_1st99th_out, file="Scale_1st99th_logistic_regression.out", sep="_")

### Test with Liver Data ##
model_comparison_out <- c()

MacPar_data <- preprocessing(load_MacPar())

### Logistic Regression ###
Logistic_models <- list()
set.seed(3792)
# Combined
model <- glm(is.cell~pct.mt+pct.ribo+prop.intronic+n.exons+gene.length, family=binomial(link="logit"), data=Training_5th95th[,colnames(Training_5th95th) %in% c(outcome, predictors)])
summary(model)
model_comparison_out <- rbind(model_comparison_out, c("Logistic_all_train", round(unlist(calc_accuracy(
					scores=model$fitted.values, truth=Training_5th95th$is.cell), digits=4)))
model_comparison_out <- rbind(model_comparison_out, c("Logistic_all_test", round(unlist(calc_accuracy(
					scores=predict(model, newdata=MacPar_data, type="response")), truth=MacPar_data$is.cell), digits=4)))
Logistic_models[["all"]] <- model

# Simplified1
# This is for any well annotated genome
# But avoids the need to quantify intronic vs exonic expression
model <- glm(is.cell~pct.mt+pct.ribo+gene.length+n.exons, family=binomial(link="logit"), data=Training_5th95th[,colnames(Training_5th95th) %in% c(outcome, predictors)])
summary(model)
model_comparison_out <- rbind(model_comparison_out, c("Logistic_no-intronic_train", round(unlist(calc_accuracy(
					scores=model$fitted.values, truth=Training_5th95th$is.cell), digits=4)))
model_comparison_out <- rbind(model_comparison_out, c("Logistic_no-intronic_test", round(unlist(calc_accuracy(
					scores=predict(model, newdata=MacPar_data, type="response")), truth=MacPar_data$is.cell), digits=4)))
Logistic_models[["no-intronic"]] <- model

# Simplified2
# This is for genomes with unreliable gene models
model <- glm(is.cell~pct.mt+pct.ribo, family=binomial(link="logit"), data=Training_5th95th[,colnames(Training_5th95th) %in% c(outcome, predictors)])
summary(model)
model_comparison_out <- rbind(model_comparison_out, c("Logistic_mt-ribo_train", round(unlist(calc_accuracy(
					scores=model$fitted.values, truth=Training_5th95th$is.cell), digits=4)))
model_comparison_out <- rbind(model_comparison_out, c("Logistic_mt-ribo_test", round(unlist(calc_accuracy(
					scores=predict(model, newdata=MacPar_data, type="response")), truth=MacPar_data$is.cell), digits=4)))
Logistic_models[["mt+ribo"]] <- model
saveRDS(Logistic_models, "Logistic_models.rds")



### SVM ###
library(e1071)
set.seed(7391)
SVM_models <- list()
# Linear Kernel
model <- svm(as.factor(is.cell)~pct.mt+pct.ribo+n.exons+gene.length, data=Training_5th95th[,colnames(Training_5th95th) %in% c(outcome, predictors)], scale=FALSE, kernel="linear")
model_comparison_out <- rbind(model_comparison_out, c("SVM_linear_train", round(unlist(calc_accuracy(
					scores=model$fitted, truth=Training_5th95th$is.cell), digits=4)))
model_comparison_out <- rbind(model_comparison_out, c("SVM_linear_test", round(unlist(calc_accuracy(
					scores=predict(model, newdata=MacPar_data, type="response")), truth=MacPar_data$is.cell), digits=4)))
SVM_models[["linear"]] <- model

# Radial Kernel
model <- svm(as.factor(is.cell)~pct.mt+pct.ribo+n.exons+gene.length, data=Training_5th95th[,colnames(Training_5th95th) %in% c(outcome, predictors)], scale=FALSE, kernel="radial")
model_comparison_out <- rbind(model_comparison_out, c("SVM_radial_train", round(unlist(calc_accuracy(
					scores=model$fitted, truth=Training_5th95th$is.cell), digits=4)))
model_comparison_out <- rbind(model_comparison_out, c("SVM_radial_test", round(unlist(calc_accuracy(
					scores=predict(model, newdata=MacPar_data, type="response")), truth=MacPar_data$is.cell), digits=4)))
SVM_models[["radial"]] <- model


# Polynomial Kernel
model <- svm(as.factor(is.cell)~pct.mt+pct.ribo+n.exons+gene.length, data=Training_5th95th[,colnames(Training_5th95th) %in% c(outcome, predictors)], scale=FALSE, kernel="polynomial")
model_comparison_out <- rbind(model_comparison_out, c("SVM_polynomial_train", round(unlist(calc_accuracy(
					scores=model$fitted, truth=Training_5th95th$is.cell), digits=4)))
model_comparison_out <- rbind(model_comparison_out, c("SVM_polynomial_test", round(unlist(calc_accuracy(
					scores=predict(model, newdata=MacPar_data, type="response")), truth=MacPar_data$is.cell), digits=4)))
SVM_models[["polynomial"]] <- model
saveRDS(SVM_models, "SVM_models.rds")

### RandomForest ###
RF_models <- list()
# requires 32 GB of RAM???
require(randomForest)
set.seed(42)
model <- randomForest(as.factor(is.cell)~pct.mt+pct.ribo+gene.length+n.exons, data=Training_5th95th, proximity=FALSE, mtry=1, ntree=500, replace=TRUE, importance=TRUE)
model_comparison_out <- rbind(model_comparison_out, c("RF_train", round(unlist(calc_accuracy(
					scores=model$fitted, truth=Training_5th95th$is.cell), digits=4)))
model_comparison_out <- rbind(model_comparison_out, c("RF_test", round(unlist(calc_accuracy(
					scores=predict(model, newdata=MacPar_data, type="response")), truth=MacPar_data$is.cell), digits=4)))
RF_models[["RF"]] <- model

saveRDS(RF_models, "RF_model.rds")

# Add mclust???
model <- MclustDA(Training_5th95th[,predictors], Training_5th95th[,outcome])
model_comparison_out <- rbind(model_comparison_out, c("Mclust_train", round(unlist(calc_accuracy(
					scores=predict(model, Training_5th95th[,predictors])$z[,"1"], truth=Training_5th95th$is.cell)), digits=4)))
model_comparison_out <- rbind(model_comparison_out, c("Mclust_test", round(unlist(calc_accuracy(
					scores=predict(model, newdata=MacPar_data[,predictors], type="response")$z[,"1"], truth=MacPar_data$is.cell)), digits=4)))


write.table(model_comparison_out, "Model_comparison_table.out", col.names=TRUE, row.names=FALSE)

