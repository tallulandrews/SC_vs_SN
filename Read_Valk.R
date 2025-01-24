# From authors: Samples GSM6598984-87 are single-cell; samples GSM6598988-90 are single-nucleus.
require(Matrix)
Valk_dir="/home/tandrew6/scratch/SC_vs_SN_Project/Valk_HumanEar"
ReadValk <- function(){
	sc_samples <- c("GSM6598984_d75_WTC-SOX2", "GSM6598985_d100_WTC-SOX2", "GSM6598986_d75_WTC-GCaMP", "GSM6598987_d100_WTC-SOX2")
	sn_samples <- c("GSM6598988_d75_WA01_LUMC04i10", "GSM6598989_d110_WA01_LUMC04i10", "GSM6598990_d110_LUMC04i10")

	create_sparse_Mat <- function(Valk_dir, prefix="", geneid=c("ensg", "symbol")) {
	        require(Matrix)
	        pattern <- paste(Valk_dir, "/",prefix,"*", sep="")
	        files <- Sys.glob(pattern);
	        if (sum(grepl("mtx", files)) != 1) {
	                stop(paste(pattern, "does not include exactly one matrix."))
	        }
	        mat <- Matrix::readMM(files[grepl("mtx", files)])
	        genes <- read.table(files[grepl("features", files)])
	        barcodes <- read.table(files[grepl("barcodes", files)])
	        if (geneid[1] == "ensg") {
	                rownames(mat) <- genes[,1]
	        } else if (geneid[1] == "symbol") {
	                rownames(mat) <- genes[,2]
	        }
	        colnames(mat) <- barcodes[,1]
	        return(mat)
	}

	all_expr <- c()
	all_anno <- c()
	for (sample in c(sc_samples, sn_samples)) {
		mat <- create_sparse_Mat(paste(Valk_dir, sample, sep="/"), geneid="symbol")
	
		stuff <- strsplit(sample, "_")
		cellline <- sapply(stuff, function(x){unlist(x[3])})
		time <- sapply(stuff, function(x){unlist(x[2])})
		if (sample %in% sc_samples) {
			is.cell <- rep(1, ncol(mat))
		} else {
			mat <- mat[match(rownames(all_expr), rownames(mat)),]
			is.cell <- rep(0, ncol(mat))
		}
		all_expr <- cbind(all_expr, mat)
		this_anno <- data.frame(cell_line=as.character(cellline), 
					time=as.character(time), 
					is.cell=factor(is.cell, levels=c("0","1")), 
					stringsAsFactors = FALSE)
		all_anno <- rbind(all_anno, this_anno)
	}
	data = list(meta=all_anno, counts=all_expr)
	saveRDS(data, paste(Valk_dir, "Valk_data.rds", sep="/"))
	return(data)
}
data <- readRDS(paste(Valk_dir, "Valk_data.rds", sep="/"))
