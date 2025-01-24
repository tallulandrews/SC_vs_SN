
args <- commandArgs(trailingOnly=TRUE)
species <- args[1] # "human" or "mouse"

match_species <- function(species) {
        if (species %in% c("Mmus", "Mus musculus", "mouse", "mice", "Mouse", "Mice")) {
                return("mouse")
        }
        if (species %in% c("Hsap", "Homo sapiens", "human", "Human")) {
                return("human")
        }
        stop(paste(species, "is not a recognized species."))
}

species <- match_species(species)

require("biomaRt")
ensembl = useMart("ensembl")

if (species == "human") {
	ensg <- useDataset("hsapiens_gene_ensembl", ensembl)
	attributes = listAttributes(ensg)

	ensg2symbol <- getBM(attributes=c('ensembl_gene_id', 
		"hgnc_symbol", 
		"chromosome_name",
		"start_position",
		"end_position",
		"strand",
		"transcription_start_site",
		"transcript_length",
		"percentage_gene_gc_content",
		"gene_biotype"
		)
		, mart = ensg)
}
if (species == "mouse") {
	ensg <- useDataset("mmusculus_gene_ensembl", ensembl)
	attributes = listAttributes(ensg)

	ensg2symbol <- getBM(attributes=c('ensembl_gene_id', 
		"hgnc_symbol", 
		"mgi_symbol",
		"chromosome_name",
		"start_position",
		"end_position",
		"strand",
		"transcription_start_site",
		"transcript_length",
		"percentage_gene_gc_content",
		"gene_biotype"
		)
		, mart = ensg)
	tmp <- ensg2symbol[,"hgnc_symbol"]
	ensg2symbol[,"hgnc_symbol"] <- ensg2symbol[,"mgi_symbol"]
	ensg2symbol[,"mgi_symbol"] <- tmp
	colnames(ensg2symbol) <- c('ensembl_gene_id', 
		"mgi_symbol",
                "hgnc_symbol",
                "chromosome_name",
                "start_position",
                "end_position",
                "strand",
                "transcription_start_site",
                "transcript_length",
                "percentage_gene_gc_content",
                "gene_biotype"
                )
}


write.table(ensg2symbol, paste("Ensemble_gene_stats", species, "various.tsv", sep="_"), sep="\t", row.names=F, col.names=T)

exon_info <- getBM(attributes=c('ensembl_gene_id', 'ensembl_exon_id', 'chromosome_name', 'exon_chrom_start', 'exon_chrom_end'), mart=ensg)
require(GenomicRanges)

granges <- GRanges(exon_info)

out <- c();
for (gene in unique(granges$ensembl_gene_id)) {
	exon_ranges <- granges[granges$ensembl_gene_id == gene,]
	exon_ranges <- reduce(exon_ranges)
	exon_length <- sum(exon_ranges@ranges@width)
	out <- rbind(out, c(gene, exon_length, length(exon_ranges)))
}

out <- out[match(ensg2symbol[,1],out[,1]),]
out <- cbind(out, ensg2symbol[,2]);
ensg2symbol_gene_length <- ensg2symbol[,"end_position"] - ensg2symbol[,"start_position"] + 1
out <- cbind(out, ensg2symbol_gene_length)
colnames(out) <- c("ensgid", "totalExonlength", "n_exons", "symbol", "genelength")
out <- unique(out)
write.table(out, paste("Ensembl_gene_stats", species, "exon_length.tsv", sep="_"), sep="\t", row.names=F, col.names=T);

exit();

### Below is DEFUNCT ###
#### Nope! ####
exon_seq <- getSequence(id=ensg2symbol$hgnc_symbol, type="hgnc_symbol", 
				seqType="gene_exon", mart=ensg)
saveRDS(exon_seq, "ensembl_exonic_sequences.rds")

whole_trans <- getSequence(id=ensg2symbol$hgnc_symbol, type="hgnc_symbol", 
					seqType="transcript_exon_intron", mart=ensg)

saveRDS(whole_trans, "ensembl_whole_transcript_sequences.rds")

### PROBLEM: genes don't match between whole trans and exon!

require(Biostrings)
gc_content  <- function(x) {
		freq <- alphabetFrequency(DNAString(x)); 
		gc <- sum(freq[c("G","C")])/sum(freq);
		return(gc)}

stats <- c();
for (i in 1:nrow(exon_seq)) {

	# Exon Length
	nt_exon <- nchar(exon_seq[i,1])

	# Intron Length
	nt_whole <- nchar(whole_trans[i,1])
	nt_intron <- nt_whole - nt_exon


	# Exonic sequnce first 400bp
	first_400 <- substr(exon_seq[i,1], 1, min(400, nt_exon))
	gc_exon_5pr <- gc_content(first_400)

	# Exonic sequence last 400bp
	last_400 <- substr(exon_seq[i,1], max(1, nt_exon-400+1), nt_exon)
	gc_exon_3pr <- gc_content(last_400)
	
	# Whole transcript first 400bp
	first_400 <- substr(whole_trans[i,1], 1, min(400, nt_whole))
	gc_whole_5pr <- gc_content(first_400)

	# Whole transcript last 400bp
	last_400 <- substr(whole_trans[i,1], max(1, nt_whole-400+1), nt_whole)
	gc_whole_3pr <- gc_content(last_400)

	stats <- rbind(stats, c(nt_exon, nt_intron, gc_exon_5pr, gc_exon_3pr, gc_whole_5pr, gc_whole_3pr))
}
rownames(stats) <- exon_seq[,2]
colnames(stats) <- c("length_exon", "length_intron", "exon_gc_5pr", "exon_gc_3pr", "trans_gc_5pr", "trans_gc_3pr")

write.table(stats, "Ensembl_transcript_stats.txt", sep="\t", row.names=T, col.names=T)
