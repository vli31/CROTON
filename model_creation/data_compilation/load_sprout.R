###############################
#### Create summary_df.tsv ####
###############################

# read in SPROUT Leenay files and output readable format
# slightly modified from `CrisprSet_Analysis.R` from Leenay et al.

# NOTE: just enough loading these libraries to parse RDS
library(S4Vectors)
library(Rsamtools)
library(CrispRVariants)  


get_alignment = function(crispr_set)
{
	alns = alns(crispr_set)
	nrep = length(alns)
	bam_names = c()
	cigars = c()
	seqs = c()
	for(i in 1:nrep) {
		aln = alns[i]
		n_total = length(aln@unlistData@cigar)
		bam_names = c(bam_names, rep(names(alns)[i], n_total) )
		seqs  = c(seqs, as.vector(aln@unlistData@elementMetadata@listData$seq) )
		cigars = c(cigars, as.vector(aln@unlistData@cigar) )
	}
	aln_df = data.frame(bam=bam_names, cigars=cigars, seqs=seqs)
	return(aln_df)
}


get_cigar_counts = function(crispr_set, min.count=10)
{
	#summcounts <- variantCounts(crispr_set, min.count = 10)
	summcounts = crispr_set$cigar_freqs
	row.avg = apply(summcounts, 1, mean)
	summcounts = summcounts[row.avg>min.count, ]
	summcounts
}


# Main entry to this script
main = function(wd, dd)
# wd : working directory
# dd : data directory
{
	# Load data
	bam_df.tbl = readRDS(file.path(dd, "Leenay_bam_df.rds"))
	bam_df = as.data.frame(bam_df.tbl)
	crispr_set_list = readRDS(file.path(dd, "Leenay_crispr_set_list.rds"))
	
	# Make folders
	count_dir = file.path(wd, "counts")
	aln_dir = file.path(wd, "aln")
	system(paste0("mkdir -p ", count_dir))
	system(paste0("mkdir -p ", aln_dir))
	
	# Placeholders
	genenames = c()
	refseq = c()
	bams = c()
	chrom = c()
	ranges = c()
	strands = c()

	# Match indices between bam_df and crispr_set
	crispr_set_gids = sapply(crispr_set_list, function(x) ifelse(class(x)=="CrisprSet", x$target$name, NA))
	bam_df_gids = bam_df$guide
	matched_crispr_set = match(bam_df_gids, crispr_set_gids)
	# if not matched everywhere (unless is NA), stop the script
	stopifnot(mean(matched_crispr_set == seq(1, length(matched_crispr_set)), na.rm=T) == 1)

	# Loop over each gRNA
	bam_colindex = seq(5, ncol(bam_df))	
	pb = utils::txtProgressBar(min = 0, max = nrow(bam_df), initial = 0, style=3)
	for(i in 1:nrow(bam_df)){
		utils::setTxtProgressBar(pb, i)
		crispr_set_index = matched_crispr_set[i]
		if(is.na(crispr_set_index)) next
		crispr_set = crispr_set_list[[crispr_set_index]]

		if(class(crispr_set)!="CrisprSet") next
		
		# information from bam_df
		genename =  as.character(bam_df$genename[i])
		genenames = c(genenames, genename)
		refseq = c(refseq, bam_df$reference[i])
		nonNA_bams = colnames(bam_df)[intersect(bam_colindex, which(!is.na(bam_df[i,])))]
		nonNA_bams = paste0(nonNA_bams, collapse=",")
		bams = c(bams, nonNA_bams)

		# information from crispr_set
		chrom = c(chrom, as.character(crispr_set$target@seqnames) )
		ranges = c(ranges, as.character(crispr_set$target@ranges) )
		strands = c(strands, as.character(crispr_set$target@strand) )
		
		# summary data by analyzing crispr_set
		cigar_counts = get_cigar_counts(crispr_set, min.count=20)
		aln_df = get_alignment(crispr_set)

		# write out
		write.table(cigar_counts, file=file.path(count_dir, paste0("counts-", genename, "-", i, ".txt")), sep="\t", quote=F )
		write.table(aln_df, file=file.path(aln_dir, paste0("aln-", genename, "-", i, ".txt")), sep="\t", quote=F )
	}

	sum_df = data.frame(index=seq(1, length(genenames)), 
			    genename=genenames,
			    refseq=refseq,
			    chrom=chrom, 
			    ranges=ranges, 
			    strand=strands, 
			    bams=bams,
			    stringsAsFactors=F
	)
	write.table(sum_df, file=file.path(wd, "summary_df.txt"), sep="\t", quote=F, row.names=F)
}

##############################
## Create Leenay_bam_df.csv ##
##############################

bamrds <- "Leenay_bam_df.rds" #downloaded from Leenay et. al. online repository
bam_df <- readRDS(bamrds)
write.csv(bam_df, file = "Leenay_bam_df.csv")