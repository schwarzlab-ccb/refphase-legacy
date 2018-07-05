
#Make a dataframe:

file.name.mapping.df <- data.frame(pc_id_region_file_name=c("Y_BTX004_SU_T1-R1_battenberg.allele.counts_chr8.txt",
															  "Y_BTX004_SU_T1-R2_battenberg.allele.counts_chr8.txt",
															  "Y_BTX005_SU_T1-R1_battenberg.allele.counts_chr14.txt",
															  "Y_BTX005_SU_T1-R2_battenberg.allele.counts_chr14.txt"), 
									chrom = c(8, 8, 14,14),
									battenberg_id_file_name=c("PD11461a_battenberg.allele.counts_chr8.txt",
															  "PD11461c_battenberg.allele.counts_chr8.txt",
															  "PD13596a_battenberg.allele.counts_chr14.txt",
															  "PD13596c_battenberg.allele.counts_chr14.txt"),
									stringsAsFactors = FALSE)

#Readin files from here:
pysam.allelecount.dir.path <- "/home/camp/watkint/working/watkint/data/yates/pysam_chr_allele_counts/"

#Write files out from here:
processed.allele.count.path <- "/home/camp/watkint/working/watkint/data/yates/processed_pysam_allele_counts/"
#Now what columns do we need to add to make the correct column
#Input
#Y_BTX004_SU_T1-R1_battenberg.allele.counts_chr8.txt
# 	A	C	G	T
# 10421	0	0	1	0
# 11879	0	16	0	0

#Output should look like :
#CHR	POS	Count_A	Count_C	Count_G	Count_T	Good_depth
# 11	70855	0	0	0	0	0
# 11	103253	0	0	0	0	0

battenberg.col.names.vec <- c("CHR", "POS", "Count_A", "Count_C", "Count_G", "Count_T", "Good_depth")
for (i in 1:nrow(file.name.mapping.df)) {


	#Read input
	input.file.path <- paste0(pysam.allelecount.dir.path, file.name.mapping.df[i, "pc_id_region_file_name"])
	input.file.df <- read.table(input.file.path, stringsAsFactors=FALSE, header=T, sep="\t")

	#Make Good_depth column
	input.file.df[,1] <- as.numeric(input.file.df[,1]) + 1
	good.depth.vec <- rowSums(input.file.df[, 2:5])
	output.df <- cbind(chrom = rep(file.name.mapping.df[i,]$chrom, length(good.depth.vec)), input.file.df, Good_depth=good.depth.vec)
	colnames(output.df) <- battenberg.col.names.vec

	print(head(output.df))

	output.file.path <- paste0(processed.allele.count.path, file.name.mapping.df[i,"battenberg_id_file_name"])
	write.table(output.df, file=output.file.path, sep="\t", quote=FALSE, row.names=FALSE)
	
}