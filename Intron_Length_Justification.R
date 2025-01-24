
gene_stats_m <- read.table(paste("/home/tandrew6/scratch/ExternalData/Ensembl/Ensembl_gene_stats", "mouse", "exon_length.tsv", sep="_"), header=TRUE)
gene_stats_m$intron_length <- gene_stats_m$genelength - gene_stats_m$totalExonlength

gene_stats_h <- read.table(paste("/home/tandrew6/scratch/ExternalData/Ensembl/Ensembl_gene_stats", "human", "exon_length.tsv", sep="_"), header=TRUE)
gene_stats_h$intron_length <- gene_stats_h$genelength - gene_stats_h$totalExonlength

cors_h <- cor(gene_stats_h[,c("totalExonlength", "n_exons", "genelength", "intron_length")])
cors_m <- cor(gene_stats_m[,c("totalExonlength", "n_exons", "genelength", "intron_length")])

require(pheatmap)
png("human_gene_stats_cors.png", width=4, height=4, units="in", res=300)
pheatmap(cors_h)
dev.off()

png("mouse_gene_stats_cors.png", width=4, height=4, units="in", res=300)
pheatmap(cors_m)
dev.off()
