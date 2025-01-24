# Test the performance of each model on varying proportions of SC:SN
set.seed(101)


load_MacPar <- function() {
        args <- Sys.glob("MacPar*.rds")

        alldata <- c()
        for (file in args) {
                predictors <- readRDS(file)
                predictors$sample <- rep(file, nrow(predictors))
                if (is.null(dim(alldata))) {
                        alldata <- predictors
                } else {
                        alldata <- rbind(alldata, predictors)
                }
        }
        alldata$is.cell <- as.numeric(alldata$is.cell)
        alldata$dataset <- "MacPar"
        return(alldata)
}


preprocessing <- function(data, probs=c(0.05,0.95)) {
        for (score in c("prop.intronic", "pct.mt", "pct.ribo", "gene.length", "n.exons")) {
                # rescale data based on 1st and 99th percentile
                thresholds <- quantile(data[,score], probs=probs)
                data[,score] <- (data[,score]-thresholds[1])/(thresholds[2]-thresholds[1])
        }
        return(data)
}

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

test_data_set <- load_MacPar()
Ns <- table(test_data_set$is.cell)
max_droplets <- min(Ns)

N_reps=10

model <- readRDS(".rds")

prefix="MacPar_Logistic_no-intronic"

sn_prop = seq(from=0.05, to=0.95, by=0.05)
out <- c()
set.seed(3781)
for (Pn in sn_prop) {
	N_sn = Pn*max_droplets
	N_sc = (1-Pn)*max_droplets
	for(rep in 1:N_reps) {
		cell_i = sample(which(test_data_set$is.cell == 1), N_sc, replace=FALSE)
		nuc_i = sample(which(test_data_set$is.cell == 0), N_sn, replace=FALSE)
		dataset <- preprocessing(test_data_set[c(cell_i, nuc_i),])
		Accuracy <- calc_accuracy(predict(model, newdata=dataset, type="response"), truth=dataset$is.cell, threshold=0.5)		
		out <- rbind(out, c(Pn, unlist(Accuracy)))
	}
	
}

write.table(out, paste(prefix, "RelativeFreq_output.txt", sep="_"))

averages = aggregate(out, by=list(out[,1]), mean)
std_errs = aggregate(out, by=list(out[,1]), sd)/sqrt(10)

png(paste(prefix, "RelativeFreq_results.png", sep="_"), width=4, height=4, units="in", res=150)
plot(averages[,1], averages[,3], type='b', xlab="nuclei proportion", ylab="Score")
arrows(averages[,1], averages[,3], averages[,1], averages[,3]+2*std_errs[,3], angle=90, length=0.1)
arrows(averages[,1], averages[,3], averages[,1], averages[,3]-2*std_errs[,3], angle=90, length=0.1)
lines(averages[,1], averages[,9], lty=2)
legend("topleft", c("accuracy", "AUC"), lty=c(1,2), bty="n")
dev.off()
