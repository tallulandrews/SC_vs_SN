require(Matrix)
Wu_dir = "/home/tandrew6/scratch/SC_vs_SN_Project/Wu_MouseKidney"
ReadWu <- function() {

	if (file.exists( paste(Wu_dir,"Wu_data.rds", sep="/") )) {
		return(readRDS(paste(Wu_dir,"Wu_data.rds", sep="/")))
	}

	a <- read.table(paste(Wu_dir,"/GSE119531_Healthy.combined.cell.annotation.txt.gz", sep=""), header=TRUE)
	mat <- read.delim(paste(Wu_dir,"/GSE119531_Healthy.combined.dge.txt.gz", sep=""))

	a[,1] <- sub("-", "\\.", a[,1])

	identical(colnames(mat), a[,1])

	a$sample <- sapply(strsplit( a[,1], "_" ), function(x){x[1]})
	a$is.cell <- rep(0, nrow(a));
	a[a$sample == "sCellDropseq","is.cell"] <- 1
	a$Phenotype <- "Healthy"

	a2 <- read.delim(paste(Wu_dir,"/GSE119531_UUO.cell.annotation.txt.gz", sep=""), header=TRUE)
	mat2 <- read.delim(paste(Wu_dir,"/GSE119531_UUO.dge.txt.gz", sep=""))

	a2[,1] <- sub("UUO", "UUO1", a2[,1])

	identical(colnames(mat2), a2[,1])

	a2$sample <- sapply(strsplit( a2[,1], "_" ), function(x){x[1]})
	a2$is.cell <- rep(0, nrow(a2));
	a2$Phenotype <- "UUO"

	all_genes <- c(rownames(mat), rownames(mat2))
	mat <- mat[match(all_genes, rownames(mat)),]
	mat[is.na(mat)] <- 0
	mat2 <- mat2[match(all_genes, rownames(mat2)),]
	mat2[is.na(mat2)] <- 0

	mat <- cbind(mat,mat2)
	a <- rbind(a, a2)

	require(Matrix)
	mat <- Matrix(as.matrix(mat))
	return(list(meta=a,counts=mat))
	saveRDS(list(meta=a,counts=mat), paste(Wu_dir,"Wu_data.rds", sep="/"))
}

data <- readRDS(paste(Wu_dir,"Wu_data.rds", sep="/"))


