require(Matrix)
Bakken_dir = "/home/tandrew6/scratch/SC_vs_SN_Project/Bakken"
convert_to_Matrix<- function(df) {
	rownames(df) <- df[,1];
	df<- df[,-1]
	mat <- Matrix(as.matrix(df))
	return(mat)
}

ReadBakken <- function() {

	if (file.exists( paste(Bakken_dir,"Bakken_data.rds", sep="/") )) {
		return(readRDS(paste(Bakken_dir,"Bakken_data.rds", sep="/")))
	}

	a <- read.delim("/home/tandrew6/projects/def-tandrew6/tandrew6/SC_vs_SN/Bakken_MouseNeurons/GSE123454_GEO_seq_template_v2.1_Nuc_vs_Cell.csv", sep=",", header=TRUE)
	a$cellID <- gsub("-", ".", a[,1])
	cells_intron <- convert_to_Matrix(read.delim(paste(Bakken_dir, "/GSE123454_cells_intron_counts.csv.gz", sep=""), sep=","))
	cells_exon <- convert_to_Matrix(read.delim(paste(Bakken_dir, "/GSE123454_cells_exon_counts.csv.gz", sep=""), sep=","))
	
	cells_exon_only <- cells_exon
	cells_exon_and_intron <- cells_exon+cells_intron
	cells_meta <- a[match(colnames(cells_exon), a$cellID),]
	cells_meta$is.cell <- 1
	rownames(cells_meta) <- colnames(cells_exon)

	nuclei_intron <- convert_to_Matrix(read.delim(paste(Bakken_dir, "/GSE123454_nuclei_intron_counts.csv.gz", sep=""), sep=","))
	nuclei_exon <- convert_to_Matrix(read.delim(paste(Bakken_dir, "/GSE123454_nuclei_exon_counts.csv.gz", sep=""), sep=","))
	
	nuclei_exon_only <- nuclei_exon
	nuclei_exon_and_intron <- nuclei_exon+nuclei_intron
	nuclei_meta <- a[match(colnames(nuclei_exon), a$cellID),]
	nuclei_meta$is.cell <- 0
	rownames(nuclei_meta) <- colnames(nuclei_exon)

	out <- list(meta=rbind(cells_meta, nuclei_meta),intron_and_exon_counts=cbind(cells_exon_and_intron, nuclei_exon_and_intron), exon_only_counts=cbind(cells_exon_only, nuclei_exon_only))
	return(out)
	saveRDS(out, paste(Bakken_dir,"Bakken_data.rds", sep="/"))
}

data <- readRDS(paste(Bakken_dir,"Bakken_data.rds", sep="/"))


