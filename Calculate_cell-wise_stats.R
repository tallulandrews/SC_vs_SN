set.seed(101)

# Arguments: 
# 1. folder of my standardized intronic / exonic cellranger output.
# 2. Species
args <- commandArgs(trailingOnly=TRUE) 
data_dir <- args[1] # the folder with "cell_ranger_exon_only & cell_ranger_introns_and_exons subfolder OR this can be an R-script that loads the data from preprocessed files + metadata files for certain datasets where these folders aren't avaiable.
species <- args[2] # Mouse or Human
is.cells <- args[3] # Whether this dataset is cells or not (aka is nuclei) 
prefix <- args[4] # what should the output be called.

#data_dir <- "/home/tandrew6/scratch/SC_vs_SN_Project/Denishenko_MouseKidney/Cells1"
#species <- "Mmus"

match_species <- function(species) {
	if (species %in% c("Mmus", "Mus musculus", "mouse", "mice", "Mouse", "Mice")) {
		return("mouse")
	}
	if (species %in% c("Hsap", "Homo sapiens", "human", "Human")) {
		return("human")
	}
	stop(paste(species, "is not a recognized species."))
}

create_sparse_Mat <- function(dir, prefix="", geneid=c("ensg", "symbol")) {
	require(Matrix)
	pattern <- paste(dir, "/",prefix,"*", sep="")
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



### ---------------------  REad in Data & Calculate % intronic reads -------------------- ####
# 

intronic.prop <- NULL
is.cell <- NULL # HERE WE FIGURE OUT IS.CELL or not from metadata.
if (grepl("\\.R$", data_dir)) {
	# Using preprocessed files from author
	# rather than remapped files
	# so no intronic proportion
	source(data_dir) # This must load the data as: data = list(counts=X, meta=Y)
	intronic.prop <- rep(NA, ncol(data$counts))
	is.cell <- data$meta$is.cell
	intron_and_exon <- data$counts
} else {

	# Account for spelling error
	cellranger_both_dir <- paste(data_dir, "cellranger_introns_plus_exons", "filtered_feature_bc_matrix", sep="/")
	if (!file.exists(cellranger_both_dir)) {
		cellranger_both_dir <- sub("plus", "and", cellranger_both_dir)
	}

	# read in the data and subset them so we have the same cells in each dataset
	intron_and_exon <- create_sparse_Mat(cellranger_both_dir, geneid="symbol")
	exon_only <- create_sparse_Mat(paste(data_dir, "cellranger_exons_only", "filtered_feature_bc_matrix", sep="/"), geneid="symbol")

	olap <- colnames(intron_and_exon)[colnames(intron_and_exon) %in% colnames(exon_only)]
	N.olap <- sum(colnames(intron_and_exon) %in% colnames(exon_only))

	print(paste(round(N.olap/(ncol(intron_and_exon))*100), "% of cells (", N.olap, ") have exon-only info.", sep=""))
	print(paste(ncol(intron_and_exon)-N.olap, "cells will be removed."))

	intron_and_exon <- intron_and_exon[,olap]
	exon_only <- exon_only[,olap]

	identical(colnames(intron_and_exon), colnames(exon_only)) # Should be True

	# Calculate intronic proportion
	## % Intronic ##
	intronic.prop <- (colSums(intron_and_exon)-colSums(exon_only))/colSums(intron_and_exon);
	intronic.prop[intronic.prop<0] <- 0
}

### ----------------------------------------------------------------------------- ####

## % Mitochondrial ##
species <- match_species(species)
pct.mt <- -1;
if (species == "mouse") {
	pct.mt <- colSums(intron_and_exon[(grepl("^mt-", rownames(intron_and_exon))),])/colSums(intron_and_exon)*100
}
if (species == "human") {
	pct.mt <- colSums(intron_and_exon[(grepl("^MT-", rownames(intron_and_exon))),])/colSums(intron_and_exon)*100
}

## % Ribosomal ##
species <- match_species(species)
pct.ribo <- -1;
if (species == "mouse") {
	pct.ribo <- colSums(intron_and_exon[(grepl("^Rp[sl]", rownames(intron_and_exon))),])/colSums(intron_and_exon)*100
}
if (species == "human") {
	pct.ribo <- colSums(intron_and_exon[(grepl("^RP[LS]", rownames(intron_and_exon))),])/colSums(intron_and_exon)*100
}

#### Splicing ####
# Remember to exclude ribosomal & mitochondrial genes!
species <- match_species(species)
# Remove ribosomal & mitochondrial transcripts
trimmed_expr_mat <- c()
if (species == "mouse") {
	is.ribo <- grepl("^Rp[sl]", rownames(intron_and_exon))
	is.mt <- grepl("^mt-", rownames(intron_and_exon))
}
if (species == "human") {
	is.ribo <- grepl("^RP[SL]", rownames(intron_and_exon))
	is.mt <- grepl("^MT-", rownames(intron_and_exon))
}
trimmed_expr_mat <- intron_and_exon[!(is.ribo | is.mt),]

# Read gene data:
gene_stats <- read.table(paste("/home/tandrew6/scratch/ExternalData/Ensembl/Ensembl_gene_stats", species, "exon_length.tsv", sep="_"), header=TRUE)
gene_stats$intron_length <- gene_stats$genelength - gene_stats$totalExonlength
gene_stats <- gene_stats[match(rownames(trimmed_expr_mat), gene_stats$symbol),]

calc_quantile_expression <- function(expr_mat, scores, probs=c(0.25, 0.5, 0.75, 1)) {
	quantiles <- quantile(scores, probs=probs, na.rm=TRUE)
	quant_expr <- sapply(quantiles, FUN=function(threshold) {
				colSums(expr_mat[!is.na(scores) & scores <= threshold,], na.rm=TRUE)
				})
	quant_expr <- quant_expr/apply(quant_expr, 1, max)
	return(quant_expr)
}

# Intron Length
intron_length <- calc_quantile_expression(trimmed_expr_mat, gene_stats$intron_length)
# Number of Exons
num_exons <- calc_quantile_expression(trimmed_expr_mat, gene_stats$n_exons)
num_exons[num_exons == 0] <- min(num_exons[num_exons > 0]/2) #WHY? Does it actually become 0?
# Gene Length
gene_length <- calc_quantile_expression(trimmed_expr_mat, gene_stats$genelength)
gene_length[gene_length == 0] <- min(gene_length[gene_length > 0]/2) #WHY? Does it actually become 0?

#gene_stats <- read.table(paste("/home/tandrew6/scratch/ExternalData/Ensembl/Ensemble_gene_stats", species, "various.tsv", sep="_"), header=TRUE)

cor(intron_length, num_exons)
cor(intron_length, gene_length)
cor(num_exons, gene_length) # most different

# Ground Truth:
# Deal with getting is.cell from metadata for the preprocessed datasets
if (is.null(is.cell)) {
	is.cell = rep(is.cells, length(pct.mt))
}


##### ADD HERE ADDITIONAL PREDICTORS
# % of total expression for genes with Transmembrane 
# Same for protein families
# GC Content use same as gene_length above
# % expression for MALAT1
# % expression for NEAT1



predictors <- data.frame(prop.intronic=intronic.prop, pct.mt=pct.mt, pct.ribo=pct.ribo, gene.length=(1-gene_length[,3])/gene_length[,1], intron.length=intron_length,  n.exons=(1-num_exons[,3])/num_exons[,1], is.cell=is.cell)

saveRDS(predictors, paste(prefix, "predictors.rds", sep="_"))


