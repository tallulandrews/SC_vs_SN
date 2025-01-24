require(Matrix)
Selewa_dir="/home/tandrew6/scratch/SC_vs_SN_Project/Selewa_HumanCardiomyocytes"
ReadSelewa <- function(){

	# sc_iPSCs
	sc_iPSCs <- read.delim(paste(Selewa_dir, "GSE129096_DROP_combined_allcells.tsv.gz", sep="/"))
	rownames(sc_iPSCs) <- sc_iPSCs[,1]
	sc_iPSCs <- sc_iPSCs[,-1]
	sc_iPSCs <- Matrix(as.matrix(sc_iPSCs))
	# Get Metadata
	is.cell <- rep(1, ncol(sc_iPSCs))
	assay <- rep("DropSeq", ncol(sc_iPSCs))
	info <- strsplit(colnames(sc_iPSCs), "_")
	time <- sapply(info, function(x){unlist(x)[2]})
	sample <- paste(time, sapply(info, function(x){unlist(x)[3]}), sep="_")
	metadata <- data.frame(assay=assay, is.cell=is.cell, day=time, sample=sample)
	rownames(metadata) <- colnames(sc_iPSCs)
	# Save data
	sc_iPSCs <- list(meta=metadata,counts=sc_iPSCs)
	saveRDS(sc_iPSCs, paste(Selewa_dir,"Selewa_HumanCardio_Cells.rds", sep="/"))

	# sn_iPSCs
	sn_iPSCs <- read.delim(paste(Selewa_dir,"GSE129096_DRONC_combined_allcells.tsv.gz", sep="/"))
	rownames(sn_iPSCs) <- sn_iPSCs[,1]
	sn_iPSCs <- sn_iPSCs[,-1]
	sn_iPSCs <- Matrix(as.matrix(sn_iPSCs))
	# Get Metadata
	is.cell <- rep(0, ncol(sn_iPSCs))
	assay <- rep("DroNcSeq", ncol(sn_iPSCs))
	info <- strsplit(colnames(sn_iPSCs), "_")
	time <- sapply(info, function(x){unlist(x)[2]})
	sample <- paste(time, sapply(info, function(x){unlist(x)[3]}), sep="_")
	metadata <- data.frame(assay=assay, is.cell=is.cell, day=time, sample=sample)
	rownames(metadata) <- colnames(sn_iPSCs)
	# Save data
	sn_iPSCs <- list(meta=metadata,counts=sn_iPSCs)
	saveRDS(sn_iPSCs, paste(Selewa_dir,"Selewa_HumanCardio_Nuclei.rds", sep="/"))

	sn_tissue <- read.delim(paste(Selewa_dir,"GSE129096_HumanHeartDroNc.tsv.gz", sep="/"))
	rownames(sn_tissue) <- sn_tissue[,1]
	sn_tissue <- sn_tissue[,-1]
	sn_tissue <- Matrix(as.matrix(sn_tissue))
	# Get Metadata
	is.cell <- rep(0, ncol(sn_tissue))
	assay <- rep("DroNcSeq", ncol(sn_tissue))
	info <- strsplit(colnames(sn_tissue), "_")
	time <- NA
	sample <- NA
	metadata <- data.frame(assay=assay, is.cell=is.cell, day=time, sample=sample)
	rownames(metadata) <- colnames(sn_tissue)
	# Save data
	sn_tissue <- list(meta=metadata,counts=sn_tissue)
	saveRDS(sn_tissue, paste(Selewa_dir,"Selewa_HumanCardio_TissueNuclei.rds",sep="/"))

        return(list(sc_iPSCs=sc_iPSCs, sn_iPSCs=sn_iPSCs, sn_tissue=sn_tissue))
}
sn_iPSCs<- readRDS(paste(Selewa_dir,"Selewa_HumanCardio_Nuclei.rds", sep="/"))
sc_iPSCs<- readRDS(paste(Selewa_dir,"Selewa_HumanCardio_Cells.rds", sep="/"))
common_genes <- intersect(rownames(sc_iPSCs$counts), rownames(sn_iPSCs$counts))
sn_iPSCs$counts <- sn_iPSCs$counts[match(common_genes, rownames(sn_iPSCs$counts)),]
sc_iPSCs$counts <- sc_iPSCs$counts[match(common_genes, rownames(sc_iPSCs$counts)),]
#sn_tissue<- readRDS(paste(Selewa_dir,"Selewa_HumanCardio_TissueNuclei.rds", sep="/")) # This one is ENSG IDs rather than symbols, it's also only 1 sample with no matching sc so just exclude it
data <- list(counts=cbind(sc_iPSCs$counts, sn_iPSCs$counts), 
		meta=rbind(sc_iPSCs$meta, sn_iPSCs$meta))
saveRDS(data, "Selewa_data.rds")


