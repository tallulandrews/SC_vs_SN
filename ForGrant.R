require(Seurat)
require(Matrix)
require(ggplot2)
color_scheme <- c("Neuronal progenitor"="#a6cee3",
			"Glutamatergic"="#1f78b4",
			"Interneuron"="#33a02c",
			"Radial glia"="#ff7f00",
			"Subplate"="grey50",
			"Microglia"="#b15928",
			"Endothelial"="#fb9a99",
			"Pericyte"="#cab2d6",
			"Endothelial-Pericyte"="#6a3d9a",
			"OPC"="goldenrod",
			"RBC"="firebrick",
			"VLMC"="#e31a1c",
			"Glial progenitor"="#fdbf6f")

# Cluster Annotation from supplementary table
SC_anno <- c("c0"="Glutamatergic-5","c1"="Interneuron-CGE", "c2"="Glutamatergic-1", "c3"="Interneuron-MGE", "c4"="Glutamatergic-4", "c5"="Glutamatergic-2", "c6"="Early Radial Glia", "c7"="Glutamatergic-7", "c8"="Cycling Progenitor", "c9"="Glutamatergic-3", "c10"="Late Radial Glia", "c11"="Multipotent glial progenitor", "c12"="Glutamatergic-6", "c13"="Subplate", "c14"="Neuronal progenitor", "c15"="Glutamatergic-8", "c16"="Microglia", "c17"="OPC", "c18"="Tuncated Radial Glia", "c19"="Pericyte", "c20"="Endothelial", "c21"="RBC", "c22"="Vascular-Leptomeningeal cell")
SC_anno_trunc <- c("c0"="Glutamatergic","c1"="Interneuron", "c2"="Glutamatergic", "c3"="Interneuron", "c4"="Glutamatergic", "c5"="Glutamatergic", "c6"="Radial glia", "c7"="Glutamatergic", "c8"="Cycling Progenitor", "c9"="Glutamatergic", "c10"="Radial glia", "c11"="Glial progenitor", "c12"="Glutamatergic", "c13"="Subplate", "c14"="Neuronal progenitor", "c15"="Glutamatergic", "c16"="Microglia", "c17"="OPC", "c18"="Radial glia", "c19"="Pericyte", "c20"="Endothelial", "c21"="RBC", "c22"="VLMC")
SN_anno <- c("c0"="Neuronal progenitor", "c1"="Glutamatergic-2", "c2"="Interneuron-1", "c3"="Glutamatergic-3", "c4"="Interneuron-2", "c5"="Radial glia", "c6"="Glutamatergic-4", "c7"="Cycling Progenitor", "c8"="OPC", "c9"="Interneuron-3", "c10"="Glutamatergic-5", "c11"="Subplate", "c12"="Microglia", "c13"="Endothelial-Pericyte")
SN_anno_trunc <- c("c0"="Neuronal progenitor", "c1"="Glutamatergic", "c2"="Interneuron", "c3"="Glutamatergic", "c4"="Interneuron", "c5"="Radial glia", "c6"="Glutamatergic", "c7"="Cycling Progenitor", "c8"="OPC", "c9"="Interneuron", "c10"="Glutamatergic", "c11"="Subplate", "c12"="Microglia", "c13"="Endothelial-Pericyte")



# Failed on my standard Interaction Session!
sc <- read.table("GSE162170_rna_counts.tsv.gz")
sc_meta <- read.table("GSE162170_rna_cell_metadata.txt.gz", header=T)

sc <- Matrix(as.matrix(sc))
obj <- CreateSeuratObject(sc, meta.data=sc_meta)
Idents(obj) <- obj@meta.data$seurat_clusters
obj <- RenameIdents(obj, SC_anno_trunc)
obj <- FindVariableFeatures(obj)
obj <- SCTransform(obj)
obj <- RunPCA(obj)
obj <- RunUMAP(obj, dims=1:20)
png("ForGrant_Cortex_SC.png", width=6, height=4, units="in", res=150)
DimPlot(obj)+scale_color_manual(values=color_scheme)
dev.off()



sn <- read.table("GSE162170_multiome_rna_counts.tsv.gz")
sn_meta <- read.table("GSE162170_multiome_cell_metadata.txt.gz", header=T)
rownames(sn_meta) <- sn_meta[,1]

sn <- Matrix(as.matrix(sn))
obj <- CreateSeuratObject(sn, meta.data=sn_meta)
Idents(obj) <- obj@meta.data$seurat_clusters
obj <- RenameIdents(obj, SN_anno_trunc)
obj <- FindVariableFeatures(obj)
obj <- SCTransform(obj)
obj <- RunPCA(obj)
obj <- RunUMAP(obj, dims=1:20)
png("ForGrant_Cortex_SN.png", width=6, height=4, units="in", res=150)
DimPlot(obj)+scale_color_manual(values=color_scheme)
dev.off()


