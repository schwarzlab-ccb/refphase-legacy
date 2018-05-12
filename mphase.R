#!/usr/bin/env Rscript
# args = commandArgs(trailingOnly=TRUE)
# patient <- args[1]
# allow.dip.test.phasing.flag <- as.logical(args[2])
# output.directory.base.path <- args[3]
# data.format.to.use <- args[4]#Defulat for this is MEDICC, but can specify TRACERx_exome, TRACERx_snp_array or TRACERx_WGS
# snp.file.path <- args[5]#Points towards the snp files, of whatever format they are given that should be used based on data.format to use.
# seg.file.path <- args[6]#Pintts towards the seg files, of whatever format they should be determined byh the data.format.to.use flag.




#function taken from https://github.com/molgenis/molgenis-pipelines/wiki/How-to-source-another_file.R-from-within-your-R-script
LocationOfThisScript  <-  function() { # Function LocationOfThisScript returns the location of this .R script (may be needed to source other files in same dir)
	this.file = NULL
	# This file may be 'sourced'
	for (i in - (1:sys.nframe())) {
		if (identical(sys.function(i), base::source)) this.file = (normalizePath(sys.frame(i)$ofile))
	}#end for (i in - (1:sys.nframe())) {

	if (!is.null(this.file)) {

	 return(dirname(this.file))

	}#end if (!is.null(this.file)) {

	# But it may also be called from the command line
	cmd.args <-  commandArgs(trailingOnly = FALSE)
	cmd.args.trailing <- commandArgs(trailingOnly = TRUE)
	cmd.args <- cmd.args[seq.int(from = 1, length.out = length(cmd.args) - length(cmd.args.trailing))]
	res <- gsub("^(?:--file=(.*)|.*)$", "\\1", cmd.args)

	# If multiple --file arguments are given, R uses the last one
	res <- tail(res[res != ""], 1)
	if (0 < length(res)) {

		return(dirname(res))

	}#end if (0 < length(res)) {

	# Both are not the case. Maybe we are in an R GUI?
	return(NULL)
}#end LocationOfThisScript  <-  function() 

current.dir <-  LocationOfThisScript()

print(current.dir)

# #Debugging option parameters:
# patient <- "MET048"
# allow.dip.test.phasing.flag  <- FALSE
# output.directory.base.path <- "/home/camp/watkint/working/watkint/data/brastianos/release_v2/B_MET048/exome/CN_het"
# data.format.to.use <- "MEDICC"
# snp.file.path <- "/home/camp/watkint/working/watkint/data/brastianos/release_v2/B_MET048/exome/MSAI/MET048_joint_snps.txt"
# seg.file.path <- "/home/camp/watkint/working/watkint/data/brastianos/release_v2/B_MET048/exome/MSAI/MET048_joint_segs_all_ascat_info.txt"
# ##########
# mphase.script.path <- "/home/camp/watkint/working/watkint/scripts/multiregioncopynum/mphase.R"
# #MEDICC format for TRACERx example
# medicc.mphase.command <- paste("Rscript", mphase.script.path, patient, allow.dip.test.phasing.flag, data.format.to.use, snp.file.path, seg.file.path, sep=" ")

##Debugging option parameters:
patient <- "MET048"
allow.dip.test.phasing.flag = FALSE
output.directory.base.path <- "/home/camp/watkint/working/watkint/data/brastianos/release_v2/B_MET048/exome/CN_het_TRACERx_no_diptest/"
data.format.to.use <- "TRACERx_exome"
snp.file.path <- paste0(output.directory.base.path, "MET048_TRACERx_input_paths_df.txt")
seg.file.path <- "/home/camp/watkint/working/watkint/data/brastianos/release_v2/B_MET048/exome/MSAI/MET048_joint_segs_all_ascat_info.txt"
##########
mphase.script.path <- "/home/camp/watkint/working/watkint/scripts/multiregioncopynum/mphase.R"
#TRACERx format for TRACERx example
tracerx.mphase.command <- paste("Rscript", mphase.script.path, patient, allow.dip.test.phasing.flag, data.format.to.use, snp.file.path, seg.file.path, sep= " ")

# #Debugging option parameters:
# patient <- "MET048"
# allow.dip.test.phasing.flag  <- FALSE
# output.directory.base.path <- "/home/camp/watkint/working/watkint/data/brastianos/release_v2/B_MET048/exome/CN_het"
# data.format.to.use <- "MEDICC"
# snp.file.path <- "/home/camp/watkint/working/watkint/data/brastianos/release_v2/B_MET048/exome/MSAI/MET048_joint_snps.txt"
# seg.file.path <- "/home/camp/watkint/working/watkint/data/brastianos/release_v2/B_MET048/exome/MSAI/MET048_joint_segs_all_ascat_info.txt"
# ##########
# mphase.script.path <- "/home/camp/watkint/working/watkint/scripts/multiregioncopynum/mphase.R"
# #MEDICC format for TRACERx example
# medicc.mphase.command <- paste("Rscript", mphase.script.path, patient, allow.dip.test.phasing.flag, data.format.to.use, snp.file.path, seg.file.path, sep=" ")


print("patient")
print(patient)
print("allow.dip.test.phasing.flag")
print(allow.dip.test.phasing.flag)
print("output.directory.base.path")
print(output.directory.base.path)
print("data.format.to.use)")
print(data.format.to.use)
print("snp.file.path")
print(snp.file.path)
print("seg.file.path")
print(seg.file.path)



###Source relevant scripts, best to do this at the top to check for errors...
source(paste0(current.dir, "/data_input_and_QC_functions_v2.R"))
source(paste0(current.dir, "/minimum_consistent_segmentation_functions_v2.R"))
source(paste0(current.dir, "/allelic_imbalance_phasing_functions_v2.R"))
source(paste0(current.dir, "/pairwise_CN_heterogeneity_comparison_functions_v2.R"))
source(paste0(current.dir, "/cancer_cell_fraction_calculation_functions_v2.R"))
source(paste0(current.dir, "/final_output_segmentation_file_generation_functions_v2.R"))
source(paste0(current.dir, "/join_segments_functions_v2.R"))
source(paste0(current.dir, "/CN_heterogeneity_plotting_functions_v2.R"))



#So let's load in data etc....
#Todo refactor into function calls at a later date.
#Note that the afunctions to create the single dataframes that it accepts are in the data_input_and_QC_functions.R

if (data.format.to.use == "MEDICC") {

	snp.file.df <- read.table(snp.file.path, stringsAsFactors=FALSE, sep="\t", header=TRUE)
	seg.file.df <- read.table(seg.file.path, stringsAsFactors=FALSE, sep="\t", header=TRUE)

	#Designed to handle the MSAI nput files
	seg.file.df$sample <- gsub("\\.", "-", gsub("_LogR", "", seg.file.df$sample))
	seg.file.df$callState <- NULL
	seg.file.df$patient <- patient
	seg.file.df$sample_name <- paste(patient,sapply(seg.file.df$sample, function(x) strsplit(x, split="-")[[1]][2]), sep=":")

	processed.orig.cpn.segs.df <- errorCheckAndAddPloidyAdjRawCPNToSegDF(patient, seg.file.df)
	sample.vec <- unique(processed.orig.cpn.segs.df$sample)

	purity.est.vec <- processed.orig.cpn.segs.df[ ! duplicated(processed.orig.cpn.segs.df$sample), "ACF"]
	names(purity.est.vec) <- sample.vec

	ploidy.est.vec <- processed.orig.cpn.segs.df[ ! duplicated(processed.orig.cpn.segs.df$sample), "Ploidy"]
	names(ploidy.est.vec) <- sample.vec

	#Todo accept precomputed mininmum consistent segmentation
	min.cons.segs.df <- generateMinimumConsistentSegmentation(processed.orig.cpn.segs.df, chroms=1:22, col.to.use="ploidy_adj_tot_cpn")
	min.cons.segs.df <- cleanAndErrorCheckMinConsSegs(min.cons.segs.df)



	all.samp.SNP.list <- list()
	for ( samp in sample.vec ) {

		ASCAT.summary <- data.frame(chr = as.numeric(snp.file.df$Chr),
			                					pos = as.numeric(snp.file.df$Pos),
			   									BAF = as.numeric(snp.file.df[, paste0(gsub("-" ,"\\.",samp),"_BAF")]),
			   									BAFsegmented = rep(NA, length(snp.file.df$Chr)),
			   									logR = as.numeric(snp.file.df[, paste0(gsub("-", "\\.", samp), "_LogR")]),
			   									logRsegmented = rep(NA, length(snp.file.df$Chr)),
			   									stringsAsFactors = FALSE)

		seg <- processed.orig.cpn.segs.df[processed.orig.cpn.segs.df$sample == samp,]

		for ( j in 1:nrow(seg) ) {
			#print(seg[j,])
		  if(nrow(ASCAT.summary[ASCAT.summary$chr==seg[j,]$chr & ASCAT.summary$pos >= seg[j,]$startpos & ASCAT.summary$pos <= seg[j,]$endpos & ASCAT.summary$BAF != 2 ,]) > 0) {
		  	
		    if(seg[j,]$nMajor==seg[j,]$nMinor){

		      ASCAT.summary[ASCAT.summary$chr==seg[j,]$chr & ASCAT.summary$pos >= seg[j,]$startpos & ASCAT.summary$pos <= seg[j,]$endpos,]$BAFsegmented <- 0.5
		      
		    } else {
	
		      avg.seg.baf <- 0.5 - mean(abs(as.numeric(ASCAT.summary[ASCAT.summary$chr==seg[j,]$chr & ASCAT.summary$pos >= seg[j,]$startpos & ASCAT.summary$pos <= seg[j,]$endpos & ASCAT.summary$BAF != 2,]$BAF) - 0.5))
		      ASCAT.summary[ASCAT.summary$chr==seg[j,]$chr & ASCAT.summary$pos >= seg[j,]$startpos & ASCAT.summary$pos <= seg[j,]$endpos,]$BAFsegmented <- avg.seg.baf
		      
		    }
	 
		  }#end if(nrow(ASCAT.summary[ASCAT.summary$chrom==seg[j,]$chr & ASCAT.summary$pos...

		  if (nrow(ASCAT.summary[ASCAT.summary$chr==seg[j,]$chr & ASCAT.summary$pos >= seg[j,]$startpos & ASCAT.summary$pos <= seg[j,]$endpos ,]) > 0) {

		  	avg.seg.logr <- mean(ASCAT.summary[ASCAT.summary$chr==seg[j,]$chr & ASCAT.summary$pos >= seg[j,]$startpos & ASCAT.summary$pos <= seg[j,]$endpos, ]$logR)
		  	ASCAT.summary[ASCAT.summary$chr==seg[j,]$chr & ASCAT.summary$pos >= seg[j,]$startpos & ASCAT.summary$pos <= seg[j,]$endpos,]$logRsegmented <- avg.seg.logr

		  }

		}#end for ( j in 1:nrow(seg) ) {
		all.samp.SNP.list[[samp]] <- ASCAT.summary	                					
	}#end for ( samp in sample.vec ) {


	asas.SNPs.w.both.lr.baf.snp.lists <- cleanAndErrorCheckSNPsMPhase(all.samp.SNP.list, min.cons.segs.df, allow.logr.only.flag = FALSE)
	print(head(str(asas.SNPs.w.both.lr.baf.snp.lists[[1]], max = 1)))

	asas.SNPS.lr.and.or.baf.lists <- cleanAndErrorCheckSNPsMPhase(all.samp.SNP.list, min.cons.segs.df, allow.logr.only.flag = TRUE)



} else if (data.format.to.use == "TRACERx_exome" | data.format.to.use == "TRACERx_snp_array") {

	#The seg file path here will also be to the segments that were chosen to be provided to the MSAI scrpt as these are the best filtered as far as we know:
	#Go to the file that contains paths to all the seg file paths
	#Go to the file that contains paths to all the snp file paths

	#So currently looks like a bit of a mix betweeen loading in the MSAI segmentaton file
	#As well as going into the ASCAT direcories and actually pulling out the original ascat objects and the seg files.... we can do all of this in this loop

	f.joint.segs.df <- read.table(seg.file.path, stringsAsFactors=FALSE, sep="\t", header=TRUE)
	f.joint.segs.df$sample <- gsub("\\.", "-", gsub("_LogR", "", f.joint.segs.df$sample))
	f.joint.segs.df$callState <- NULL
	f.joint.segs.df$patient <- patient
	f.joint.segs.df$sample_name <- paste(patient,sapply(f.joint.segs.df$sample, function(x) strsplit(x, split="-")[[1]][2]), sep=":")
	array.type <- "exome"
	#Todo workout how much an issue this hospital identifier thing is, can we assume that for the purposes of this that the icentifier is going to be full? or will that break things?
	#For the moment keep as is and concentrate entirely on the brastianos dataset.
	full.pt <- paste("B_",patient, sep="")
	orig.cpn.segs.df <- f.joint.segs.df
	processed.orig.cpn.segs.df <- errorCheckAndAddPloidyAdjRawCPNToSegDF(patient, orig.cpn.segs.df)
	ploidy.est.vec <- processed.orig.cpn.segs.df[!duplicated(processed.orig.cpn.segs.df$sample),'Ploidy']
	names(ploidy.est.vec) <- unique(processed.orig.cpn.segs.df$sample)
	purity.est.vec <- processed.orig.cpn.segs.df[!duplicated(processed.orig.cpn.segs.df$sample),'ACF']
	names(purity.est.vec) <- unique(processed.orig.cpn.segs.df$sample)



	#Todo accept precomputed mininmum consistent segmentation
	min.cons.segs.df <- generateMinimumConsistentSegmentation(processed.orig.cpn.segs.df, chroms=1:22, col.to.use="ploidy_adj_tot_cpn")
	print(head(min.cons.segs.df))


	min.cons.segs.df <- cleanAndErrorCheckMinConsSegs(min.cons.segs.df)
	print(head(min.cons.segs.df))

	#This is wha we need to refactor, we just want a path to afile that contains all the data that these functions go and look for to try and grep the correct file paths
	#We can do this  here through refactoring, then can pull out that functionality and put in the existing CN heterogeneity TRACERx scripts that will actaully call 
	#This command line MPhase implemenation.

	#So we want to read in a file that looks like this:
	#patient_id sample_id ascat_bc_file_path 
	pt.snp.file.paths.df <- read.table(snp.file.path, stringsAsFactors=FALSE, sep="\t", header=TRUE)
	all.samp.SNP.list <- getSNPDataForTRACERxExomeSNPArray(pt.snp.file.paths.df)



	#all sample all seg logr baf lists is asas.lr.baf
	#All SNPs in these lists must have both a logr and baf value
	asas.SNPs.w.both.lr.baf.snp.lists <- cleanAndErrorCheckSNPs(all.samp.SNP.list, min.cons.segs.df, allow.logr.only.flag = FALSE)
	#Todo, is it necessary to keep both lists of SNPs like this?
	#SNPs that only have LogR are permitted in these lists, but it also contains all the SNPs with both BAF and LogR
	asas.SNPS.lr.and.or.baf.lists <- cleanAndErrorCheckSNPs(all.samp.SNP.list, min.cons.segs.df, allow.logr.only.flag = TRUE)


} else if (data.format.to.use == "TRACERx_WGS") {

	#As fas far as I can tell the whole segment stage of the for the WGS is identical to that for the exome and has been copied and pasted from the else if statement above:

	f.joint.segs.df <- read.table(seg.file.path, stringsAsFactors=FALSE, sep="\t", header=TRUE)
	f.joint.segs.df$sample <- gsub("\\.", "-", gsub("_LogR", "", f.joint.segs.df$sample))
	f.joint.segs.df$callState <- NULL
	f.joint.segs.df$patient <- patient
	f.joint.segs.df$sample_name <- paste(patient,sapply(f.joint.segs.df$sample, function(x) strsplit(x, split="-")[[1]][2]), sep=":")
	array.type <- "exome"
	#Todo workout how much an issue this hospital identifier thing is, can we assume that for the purposes of this that the icentifier is going to be full? or will that break things?
	#For the moment keep as is and concentrate entirely on the brastianos dataset.
	full.pt <- paste("B_",patient, sep="")
	orig.cpn.segs.df <- f.joint.segs.df
	processed.orig.cpn.segs.df <- errorCheckAndAddPloidyAdjRawCPNToSegDF(patient, orig.cpn.segs.df)
	ploidy.est.vec <- processed.orig.cpn.segs.df[!duplicated(processed.orig.cpn.segs.df$sample),'Ploidy']
	names(ploidy.est.vec) <- unique(processed.orig.cpn.segs.df$sample)
	purity.est.vec <- processed.orig.cpn.segs.df[!duplicated(processed.orig.cpn.segs.df$sample),'ACF']
	names(purity.est.vec) <- unique(processed.orig.cpn.segs.df$sample)


	#Todo accept precomputed mininmum consistent segmentation
	min.cons.segs.df <- generateMinimumConsistentSegmentation(processed.orig.cpn.segs.df, chroms=1:22, col.to.use="ploidy_adj_tot_cpn")
	print(head(min.cons.segs.df))


	min.cons.segs.df <- cleanAndErrorCheckMinConsSegs(min.cons.segs.df)
	print(head(min.cons.segs.df))

	#So we want to read in a file that looks like this:
	#patient_id sample_id tumor_baf_logr_file_path tumor_seg_data_path


	#Ok so let's make this file, what do we need:
	#For the TRACERx WGS format as it stands we don't care about QC directories as there aren't any...
	pt.snp.file.paths.df <- read.table(snp.file.path, stringsAsFactors=FALSE, sep="\t", header=TRUE)
	#So now we want to take the dataframe:
	all.samp.SNP.list <- getSNPDataForTRACERxWGS(pt.snp.file.paths.df)



	#all sample all seg logr baf lists is asas.lr.baf
	#All SNPs in these lists must have both a logr and baf value
	asas.SNPs.w.both.lr.baf.snp.lists <- cleanAndErrorCheckSNPs(all.samp.SNP.list, min.cons.segs.df, allow.logr.only.flag = FALSE)
	#Todo, is it necessary to keep both lists of SNPs like this?
	#SNPs that only have LogR are permitted in these lists, but it also contains all the SNPs with both BAF and LogR
	asas.SNPS.lr.and.or.baf.lists <- cleanAndErrorCheckSNPs(all.samp.SNP.list, min.cons.segs.df, allow.logr.only.flag = TRUE)




} else {

	print(data.format.to.use)
	stop("Incorrect data.format.to.use input to mphase")

}


#Todo add saving of all relevent files in correct positions....

#we can create a dataframe of the logr snps, we don't need to pass them all around in lists like the ccfs

#append dataframes


#####Get relevant libraries:
#Todo check this is the only required library.
library(diptest)
library(gplots)
library(gdata)




###############This is past the input stage and into the analysis:
#Set constants for all analyses here... eventually these can be command line set
MIN_NUM_HZ_SNP_PHASING <- 5
#PLOT_OUTPUT_DIR <- NULL

PLOT_OUTPUT_DIR <- paste(output.directory.base.path, "/", sep="")

if ( ! file.exists(PLOT_OUTPUT_DIR) ) {

	print(PLOT_OUTPUT_DIR)
	print("file doesn't exist")

}


GAMMA_VAL <- 1
MAX_CPN_FOR_SUBCLONALITY <- 10
MIN_CCF_TO_CONSIDER <- 0.15
MIN_CPN_TO_CONSIDER <- 0.15
MIN_NUM_CHRS <- 3
PVAL_THRESH_VS_MIN_CCF <- 0.001
PVAL_THRESH_VS_OTH_CHRS_SAME_CPN <- 0.01
GAIN_LOGR_THRESH <- log2(2.5 / 2)
LOSS_LOGR_THRESH <- log2(1.5 / 2)
AMP_LOGR_THRESH <- log2(4 / 2)
DEEP_LOSS_LOGR_THRESH <- log2(0.5 / 2)
PVAL_THRESH_VS_LOGR_FOR_CPN_EVENT <- 0.01
CHRS_TO_USE_VEC <- 1:22
CONF_INT_FOR_CPNX_BAF_VS_EXP <- 0.9999
PVAL_THRESH_TO_JOIN_PH_SEGS <- 0.001

MIN_PROP_SNP_SEP_THRESH <- 0.75
PVAL_COMB_WILCOX_PH_THRESH <- 0.01
PVAL_DIPTEST_AI_THRESH <- 0.01
ALLOW_DIPTEST_ALONE_PH_FLAG <- allow.dip.test.phasing.flag

if (is.null(output.directory.base.path)) {

	stop("No output path provided, aborting...")

}


#so now we need to read in the MSAI data:
if ( ! file.exists(snp.file.path)) {

	stop("No file at the snp.file.path provided in the second argument")

}

if ( ! file.exists(seg.file.path)) {

	stop("No file at the seg.file.path provided in the third argument")

}



#Todo, further refactor code so that everything is done at the level of a segment
#Todo, add options to phase using baf dip test not just ascat regions
#Todo, reduce the number of parameters per function
#Todo, Add function definitions documentaiton
#Todo, clarify the loops and next statements that are controlling what happens when there are less than 5 hz snps in a sample.
#Todo, check if the mixture phasing is ever being selected... in the example with B_MET002 no mixture phasing called.

asas.SNPs.w.both.lr.baf.snp.lists.path <- paste0(PLOT_OUTPUT_DIR, patient, "_asas_SNPs_w_both_lr_baf_snp_lists.RDS")
saveRDS(asas.SNPs.w.both.lr.baf.snp.lists, file = asas.SNPs.w.both.lr.baf.snp.lists.path)

min.cons.segs.df.path <- paste0(PLOT_OUTPUT_DIR, patient, "_min_cons_segs_df.RDS")
saveRDS(min.cons.segs.df, file = min.cons.segs.df.path)

ploidy.est.vec.path <- paste0(PLOT_OUTPUT_DIR, patient, "_ploidy_est_vec.RDS")
saveRDS(ploidy.est.vec, file = ploidy.est.vec.path)

purity.est.vec.path <- paste0(PLOT_OUTPUT_DIR, patient, "_purity_est_vec.RDS")
saveRDS(purity.est.vec, file = purity.est.vec.path)

all.seg.phasing.info.list <- getCombinedBAFDipTestPhasing(asas.lr.baf.lists = asas.SNPs.w.both.lr.baf.snp.lists,
                                            min.cons.segs.df = min.cons.segs.df,
                                            ploidy.est.vec = ploidy.est.vec,
                                            purity.est.vec = purity.est.vec,
                                            gamma.val = GAMMA_VAL,
                                            min.num.hz.snp.phasing = MIN_NUM_HZ_SNP_PHASING,
                                            min.prop.snp.sep.threshold = MIN_PROP_SNP_SEP_THRESH,
                                            dip.test.p.val.threshold = PVAL_DIPTEST_AI_THRESH,
                                            comb.wilcox.snp.sep.p.val.threshold = PVAL_COMB_WILCOX_PH_THRESH,
                                            allow.dip.test.without.ASCAT.to.phase = ALLOW_DIPTEST_ALONE_PH_FLAG
                                            )
all.seg.phasing.info.list.path <- paste0(PLOT_OUTPUT_DIR, patient, "_all_seg_phasing_info_list.RDS")
saveRDS(all.seg.phasing.info.list.path, file = all.seg.phasing.info.list.path)
#Todo unsure of what to do with the baf logr lists which have positions that only have logr here, at the moment exclude them but how
#does that affect later run
patient.phased.BAF.df.list <- createAllSegsPhasedBAFDataFrame(asas.lr.baf.lists = asas.SNPs.w.both.lr.baf.snp.lists,
                                                         sample.to.phase.seg.list = all.seg.phasing.info.list,
                                                         min.num.hz.snp.phasing = MIN_NUM_HZ_SNP_PHASING
                                                         )
patient.phased.BAF.df.list.path <- paste0(PLOT_OUTPUT_DIR, patient, "_phased_BAF_df_list.RDS")
saveRDS(patient.phased.BAF.df.list, file = patient.phased.BAF.df.list.path)

#Todo check wha happens to the missing rows, are these segments that contained less than 5 snps and thus didn't have something created?

sample.vec <- names(asas.SNPs.w.both.lr.baf.snp.lists)


#Get the logr dataframe list we'll use to add gain/loss calls:

#This will be used to work out a total raw and used for logr 
all.samp.logr.df.list <- list()

for ( samp in sample.vec ) {

	samp.logr.baf.df <- do.call(rbind, asas.SNPS.lr.and.or.baf.lists[[samp]])

	if ( nrow(samp.logr.baf.df[is.na(samp.logr.baf.df$BAF) ,]) > 0 ) {

		samp.logr.df <- samp.logr.baf.df[ ! is.na(samp.logr.baf.df$BAF) ,]

	} else if ( nrow(samp.logr.baf.df[samp.logr.baf.df$BAF == 2 ,]) > 0) {

		samp.logr.df <- samp.logr.baf.df[ samp.logr.baf.df$BAF == 2 ,]

	} else {

		print("No logr only SNPs, 0 row dataframe will be passed to functions.")
		samp.logr.df <- samp.logr.baf.df[ samp.logr.baf.df$BAF == 2 ,]

	}


	if (nrow(samp.logr.baf.df) > 0 ) {
		print("head(samp.logr.df)")
		print(head(samp.logr.df))
 		mock.cpnA.cpnB.logr.df <- calculateSNPsCopyNumberAB(logr.vec = samp.logr.df$logR,
	                                                        baf.vec = rep(0.5, nrow(samp.logr.df)),
	                                                        tumor.ploidy = ploidy.est.vec[samp],
	                                                        tumor.purity = purity.est.vec[samp],
	                                                        gamma.val = GAMMA_VAL)
 		samp.logr.df$tot_raw_un_ph_cpn <- unlist(rowSums(mock.cpnA.cpnB.logr.df))
 		print("head(samp.logr.df)")
 		print(head(samp.logr.df))
		#samp.logr.df$tot_raw_un_ph_cpn <- rowSums(mock.cpnA.cpnB.logr.df)

	}

	all.samp.logr.df.list[[samp]] <- samp.logr.df
	
}

## how to make total raw, even though no baf value





all.pw.cn.het.comp.lists <- list()
for (i in 1:(length(sample.vec) - 1)) {

  for(j in (i + 1):length(sample.vec)) {

    samp.one <- sample.vec[i]
    samp.two <- sample.vec[j]
    samp.one.two <- paste(samp.one, ":", samp.two, sep="")

    samp.one.two.comparison.list <- calculatePairwiseCNHeterogeneity(samp.one = samp.one,
                                                                  samp.two = samp.two,
                                                                  asas.both.lr.baf.lists = asas.SNPs.w.both.lr.baf.snp.lists,
                                                                  asas.lr.and.or.baf.lists = asas.SNPS.lr.and.or.baf.lists,
                                                                  min.cons.segs.df = min.cons.segs.df,
                                                                  ploidy.est.vec = ploidy.est.vec,
                                                                  purity.est.vec = purity.est.vec,
                                                                  all.seg.phasing.info.list = all.seg.phasing.info.list,
                                                                  gamma.val = GAMMA_VAL,
                                                                  min.num.hz.snp.phasing = MIN_NUM_HZ_SNP_PHASING,
                                                                  gain.threshold = GAIN_LOGR_THRESH,
                                                                  loss.threshold = LOSS_LOGR_THRESH,
                                                                  amp.threshold = AMP_LOGR_THRESH,
                                                                  deep.loss.threshold = DEEP_LOSS_LOGR_THRESH
                                                                  )
    #Todo, no need to keep reassigning here...
    for ( seg.name in names(samp.one.two.comparison.list) ) {

      all.pw.cn.het.comp.lists[[samp.one]][[seg.name]] <- samp.one.two.comparison.list[[seg.name]][[samp.one]]
      all.pw.cn.het.comp.lists[[samp.two]][[seg.name]] <- samp.one.two.comparison.list[[seg.name]][[samp.two]]
      all.pw.cn.het.comp.lists[[samp.one.two]][[seg.name]] <- samp.one.two.comparison.list[[seg.name]][[samp.one.two]]

    }

    all.pw.cn.het.comp.lists[[samp.one.two]]$pairwise_avg_cn_het_score <- samp.one.two.comparison.list[[seg.name]][[samp.one.two]]$pairwise_avg_cn_het_score

  }#end for (i in 1:(length(sample.vec) - 1)) {

}#end for(j in (i + 1):length(sample.vec)) {
all.pw.cn.het.comp.lists.path <- paste0(PLOT_OUTPUT_DIR, patient, "_all_pw_cn_het_comp_lists.RDS")
saveRDS(all.pw.cn.het.comp.lists, file=all.pw.cn.het.comp.lists.path)


#Let's get the CCF
seg.name.vec <- names(all.pw.cn.het.comp.lists[[1]])
#Check to make sure that all segments are present in all samples... might not be necessary todo, check.
for ( samp in sample.vec ) {

      seg.name.vec <- intersect(seg.name.vec, names(all.pw.cn.het.comp.lists[[samp]]))

}#end for ( samp in sample.vec ) {   


all.seg.CCF.list <- list()
all.seg.int.ob.xp.cpn.df.list <- list()
all.seg.CCF.diff.from.exp.list <- list()
all.seg.clonal.in.all.samps.list <- list()
all.seg.greater.than.two.diff.list <- list()

#Is it clonal on the same cpn or clonal on different segment..

#next assess for each chromosome and region, what are the the raw values underlying CN calls, these will be used to help test for subclonality.

#Check that has list structure so that not testing the the values from a segment against itself...
samps.int.cpn.raw.cpn.by.chr.list <- generateAllSamplesByChromIntCPNUnderlyingRawValuesLists(all.pw.cn.het.comp.lists, sample.vec, seg.name.vec)
samps.int.cpn.raw.cpn.by.chr.list.path <- paste0(PLOT_OUTPUT_DIR, patient, "_samps_int_cpn_raw_cpn_by_chr_list.RDS")
saveRDS(samps.int.cpn.raw.cpn.by.chr.list, file=samps.int.cpn.raw.cpn.by.chr.list.path)


for (cpnX in c("cpnA", "cpnB")) {

  #Here is where we're going to have iniatalise on all the lists
  all.seg.CCF.list[[cpnX]] <- list()
  all.seg.int.ob.xp.cpn.df.list[[cpnX]] <- list()
  all.seg.CCF.diff.from.exp.list[[cpnX]] <- list()
  all.seg.clonal.in.all.samps.list[[cpnX]] <- list()
  all.seg.greater.than.two.diff.list[[cpnX]] <- list()

  for ( seg.name in seg.name.vec ) {

    all.seg.CCF.list[[cpnX]][[seg.name]] <- list()
    all.seg.int.ob.xp.cpn.df.list[[cpnX]][[seg.name]] <- list()
    all.seg.CCF.diff.from.exp.list[[cpnX]][[seg.name]] <- list()
    all.seg.clonal.in.all.samps.list[[cpnX]][[seg.name]] <- list()
    all.seg.greater.than.two.diff.list[[cpnX]][[seg.name]] <- list()

  }

  for (seg.name in seg.name.vec) {

    if (is.na(all.pw.cn.het.comp.lists[[1]][[seg.name]]$ph_cpnA_vec[1])) {

      all.seg.CCF.list[[cpnX]][[seg.name]] <- NA
      all.seg.int.ob.xp.cpn.df.list[[cpnX]][[seg.name]] <- NA
      all.seg.CCF.diff.from.exp.list[[cpnX]][[seg.name]] <- NA
      all.seg.clonal.in.all.samps.list[[cpnX]][[seg.name]] <- NA
      all.seg.greater.than.two.diff.list[[cpnX]][[seg.name]] <- NA
      next 

    }

    #does both A and B, todo, refactor so only does one
    s.samps.int.ob.xp.cpn.df <- getSegObsExpCPNBySampleDF(all.pw.cn.het.comp.lists = all.pw.cn.het.comp.lists,
                                                          sample.vec = sample.vec,
                                                          seg.name = seg.name,
                                                          cpnX.vec = cpnX,
                                                          reduced.seg.flag = FALSE) 

    #does only cpnX
    seg.cpnX.diff.from.exp.list <- calculateWhetherSegCPNDiffFromIntegerExpectedAndAllMatchingSampCPN(all.pw.cn.het.comp.lists = all.pw.cn.het.comp.lists,
                                                                                                      s.samps.int.ob.xp.cpn.df = s.samps.int.ob.xp.cpn.df,
                                                                                                      seg.name = seg.name,
                                                                                                      samps.int.cpn.raw.cpn.by.chr.list = samps.int.cpn.raw.cpn.by.chr.list,
                                                                                                      cpnX = cpnX,
                                                                                                      chrs.can.use.for.ccf.test = CHRS_TO_USE_VEC,
                                                                                                      reduced.seg.flag = FALSE)
    #does only cpnX
    is.segment.clonal.in.all.samps.flag <- determineCCFClonalityOfSegment(s.samps.int.ob.xp.cpn.df = s.samps.int.ob.xp.cpn.df,
                                                                          seg.name = seg.name,
                                                                          all.samp.ph.BAF.df.list = patient.phased.BAF.df.list,
                                                                          purity.est.vec = purity.est.vec,
                                                                          cpnX = cpnX,
                                                                          min.CCF.to.consider.val = MIN_CCF_TO_CONSIDER,
                                                                          confidence.interval.threshold = CONF_INT_FOR_CPNX_BAF_VS_EXP)

    #does only cpnX
    is.there.cpnX.difference.greater.than.two.flag <- calculateTwoValueDiffFlag(s.samps.int.ob.xp.cpn.df = s.samps.int.ob.xp.cpn.df,
                                                                seg.cpnX.diff.from.exp.list = seg.cpnX.diff.from.exp.list,
                                                                cpnX = cpnX,
                                                                min.CCF.to.consider.val = MIN_CCF_TO_CONSIDER,
                                                                confidence.interval.threshold = CONF_INT_FOR_CPNX_BAF_VS_EXP,
                                                                min.num.other.chrs = MIN_NUM_CHRS)

    #does only cpnX
    seg.psc.low.high.cpn.and.use.preset.ccf.flag.vec <- getPossibleLowAndHighSubclonalCPNStates(seg.name = seg.name,
                                                                                                cpnX = cpnX,
                                                                                                s.samps.int.ob.xp.cpn.df = s.samps.int.ob.xp.cpn.df,
                                                                                                all.clonal.flag = is.segment.clonal.in.all.samps.flag,
                                                                                                two.value.diff.flag = is.there.cpnX.difference.greater.than.two.flag,
                                                                                                min.CCF.to.consider.val = MIN_CCF_TO_CONSIDER)
    seg.psc.cpnX.low  <- seg.psc.low.high.cpn.and.use.preset.ccf.flag.vec[1]
    seg.psc.cpnX.high  <- seg.psc.low.high.cpn.and.use.preset.ccf.flag.vec[2]
    seg.use.default.ccf <- as.logical(seg.psc.low.high.cpn.and.use.preset.ccf.flag.vec[3])

    #update segment for each sample 
    for ( samp  in sample.vec ) {

      all.seg.CCF.list[[cpnX]][[seg.name]][[samp]] <- list()
      ss.ph.cpnX.vals.vec <- all.pw.cn.het.comp.lists[[samp]][[seg.name]][[paste0("ph_", cpnX, "_vec")]]

      seg.samp.CCF.psc.cpn.df <- getSegCPNXMeanAndAllRawCCFsDF(ss.mean.cpnX.obs = s.samps.int.ob.xp.cpn.df[samp, paste0("mean_obs_", cpnX)],
                                                              ss.ph.cpnX.vals.vec = ss.ph.cpnX.vals.vec,
                                                              seg.psc.cpnX.low = seg.psc.cpnX.low,
                                                              seg.psc.cpnX.high = seg.psc.cpnX.high,
                                                              seg.use.default.ccf.flag = seg.use.default.ccf,
                                                              all.clonal.flag = is.segment.clonal.in.all.samps.flag,
                                                              max.cpn.val = MAX_CPN_FOR_SUBCLONALITY)

      all.seg.CCF.list[[cpnX]][[seg.name]][[samp]][[paste0("poss_sc_", cpnX, "_low")]] <- seg.samp.CCF.psc.cpn.df$poss_sc_cpnX_low
      all.seg.CCF.list[[cpnX]][[seg.name]][[samp]][[paste0("poss_sc_", cpnX, "_high")]] <- seg.samp.CCF.psc.cpn.df$poss_sc_cpnX_high
      all.seg.CCF.list[[cpnX]][[seg.name]][[samp]][[paste0("mean_frac_", cpnX, "_low")]] <- seg.samp.CCF.psc.cpn.df$mean_frac_cpnX_low
      all.seg.CCF.list[[cpnX]][[seg.name]][[samp]][[paste0("mean_frac_", cpnX, "_high")]] <- seg.samp.CCF.psc.cpn.df$mean_frac_cpnX_high
      all.seg.CCF.list[[cpnX]][[seg.name]][[samp]][[paste0(cpnX, "_orig_ph_obs_cpnA_vec")]]   <- all.pw.cn.het.comp.lists[[samp]][[seg.name]][[paste0("ph_cpnA_vec")]]
      all.seg.CCF.list[[cpnX]][[seg.name]][[samp]][[paste0(cpnX, "_orig_ph_obs_cpnB_vec")]]   <- all.pw.cn.het.comp.lists[[samp]][[seg.name]][[paste0("ph_cpnB_vec")]]
      all.seg.CCF.list[[cpnX]][[seg.name]][[samp]][[paste0("clean_ph_obs_", cpnX, "_vec")]]   <- seg.samp.CCF.psc.cpn.df$clean_ph_obs_cpnX
      all.seg.CCF.list[[cpnX]][[seg.name]][[samp]][[paste0("frac_", cpnX, "_low_vec")]]  <- seg.samp.CCF.psc.cpn.df$frac_cpnX_low 
      all.seg.CCF.list[[cpnX]][[seg.name]][[samp]][[paste0("frac_", cpnX, "_high_vec")]]  <- seg.samp.CCF.psc.cpn.df$frac_cpnX_high

    }#end  for ( samp  in sample.vec ) {

    all.seg.int.ob.xp.cpn.df.list[[cpnX]][[seg.name]] <- s.samps.int.ob.xp.cpn.df
    all.seg.clonal.in.all.samps.list[[cpnX]][[seg.name]] <- is.segment.clonal.in.all.samps.flag
    all.seg.greater.than.two.diff.list[[cpnX]][[seg.name]] <- is.there.cpnX.difference.greater.than.two.flag
    all.seg.CCF.diff.from.exp.list[[cpnX]][[seg.name]] <- seg.cpnX.diff.from.exp.list

  }#for (seg.name %in% seg.name.vec 

}#end for (cpnX in c("cpnA", "cpnB")) {

all.seg.CCF.list.path <- paste0(PLOT_OUTPUT_DIR, patient, "_all_seg_CCF_list.RDS")
saveRDS(all.seg.CCF.list, file=all.seg.CCF.list.path)

all.seg.int.ob.xp.cpn.df.list.path <- paste0(PLOT_OUTPUT_DIR, patient, "_all_seg_int_ob_xp_cpn_df_list.RDS")
saveRDS(all.seg.int.ob.xp.cpn.df.list, file=all.seg.int.ob.xp.cpn.df.list.path)

all.seg.CCF.diff.from.exp.list.path <- paste0(PLOT_OUTPUT_DIR, patient, "_all_seg_CCF_diff_from_exp_list.RDS")
saveRDS(all.seg.CCF.diff.from.exp.list, file=all.seg.CCF.diff.from.exp.list.path)

all.seg.clonal.in.all.samps.list.path <- paste0(PLOT_OUTPUT_DIR, patient, "_all_seg_clonal_in_all_samps_list.RDS")
saveRDS(all.seg.clonal.in.all.samps.list, file=all.seg.clonal.in.all.samps.list.path)

all.seg.greater.than.two.diff.list.path <- paste0(PLOT_OUTPUT_DIR, patient, "_all_seg_greater_than_two_diff_list.RDS")
saveRDS(all.seg.greater.than.two.diff.list, file=all.seg.greater.than.two.diff.list.path)

#Todo remove the dataframe initialisation from the function and modify it to work on a segment by segment basis.
#Todo rename clonal.pval to be more meaningful...
#Todo Check [1] "B_MET019_MT-M2:1:762496:53608176" to see if LOW/HIGH fraction working properly.

#Ok so our issue is we're always getting the low because the CCFs are not being recalculated correctly.
#We need to record whether a segment is clonal within a sample or not

#Let's initialise the first joint output:

#empty.unjoined.mcr.segs.df <- initialiseAllSegFullDataframe(seg.name.vec, sample.vec, make.row.names.seg.name.alone = FALSE)
#So no


empty.full.joint.segs.df <- initialiseFullOutputDataframe(seg.name.vec = names(all.pw.cn.het.comp.lists[[1]]),
														ploidy.est.vec = ploidy.est.vec,
														purity.est.vec = purity.est.vec,
														make.row.names.seg.name.alone = FALSE) 


#This adds the raw cpn for both segments that could be phased and those that couldn't
source("final_output_segmentation_file_generation_functions_v2.R")
new.cpn.full.joint.segs.df <- addNewRawAndIntCPNCallsToFinalSegmentation(seg.name.vec = names(all.pw.cn.het.comp.lists[[1]]),
																		sample.vec = names(purity.est.vec),
																		all.samp.logr.df.list = all.samp.logr.df.list,
																		all.samp.ph.BAF.df.list = patient.phased.BAF.df.list,
																		all.seg.CCF.vals.list = all.seg.CCF.list,
																		asas.lr.baf.lists = asas.SNPs.w.both.lr.baf.snp.lists,
																		all.seg.df = empty.full.joint.segs.df,
																		cpnX.vec = c("cpnA", "cpnB")
																		)


#We also need to add the original

#Add CCF information and within sample clonality information to the dataframe...
clonality.and.ccf.joint.segs.df <- addWithinSampleClonalityAndCCFToFinalSegmentation(seg.name.vec = names(all.pw.cn.het.comp.lists[[1]])
                                ,all.seg.CCF.vals.list = all.seg.CCF.list
                                ,all.seg.CCF.diff.from.exp.list = all.seg.CCF.diff.from.exp.list
                                ,all.seg.df = empty.full.joint.segs.df
                                ,sample.vec = names(purity.est.vec)
                                ,clonal.pval = PVAL_THRESH_VS_MIN_CCF
                                ,diff.from.oth.p.val = PVAL_THRESH_VS_OTH_CHRS_SAME_CPN
                                ,min.CCF.to.consider.val = MIN_CCF_TO_CONSIDER
                                ,min.cpn.to.consider.val = MIN_CPN_TO_CONSIDER
                                ,min.num.other.chrs = MIN_NUM_CHRS
                                ,chrs.to.use.vec = CHRS_TO_USE_VEC
                                ,cpnX.vec = c("cpnA", "cpnB")
                                ,max.cpn.val = MAX_CPN_FOR_SUBCLONALITY
                                )


#This will also add it to segments that were too small for 
full.joint.segs.df.w.cpn.vs.ploidy.df <- addLogrDefinedCPNEventsToOutputSegmentation(purity.vec = purity.vec,
																					ploidy.vec = ploidy.vec,	
                                                        							seg.name.vec = names(all.pw.cn.het.comp.lists[[1]]),
                                                        							all.seg.df = clonality.and.ccf.joint.segs.df,
                                                        							all.samp.logr.df.list = all.samp.logr.df.list,
                                                        							all.samp.ph.BAF.df.list = patient.phased.BAF.df.list,     
                                                       								gain.threshold = GAIN_LOGR_THRESH,
                                                        							loss.threshold = LOSS_LOGR_THRESH,
                                                        							amp.threshold = AMP_LOGR_THRESH,  
                                                        							deep.loss.threshold = DEEP_LOSS_LOGR_THRESH,
                                                        							tot.cpn.vs.ploidy.pval.thresh = PVAL_THRESH_VS_LOGR_FOR_CPN_EVENT)

#These are kept in separate functions as the original segmentaiton was refactored so it could run when segments are rejoined at the end of the script.
full.joint.segs.w.phasing.df <- addPhasingInformationToOutputSegmentation(sample.vec = names(purity.est.vec),
                                        seg.name.vec = names(all.pw.cn.het.comp.lists[[1]]), 
                                        sample.to.phase.seg.list = all.seg.phasing.info.list,
                                        all.pw.cn.het.comp.lists = all.pw.cn.het.comp.lists,
                                        all.segs.df = full.joint.segs.df.w.cpn.vs.ploidy.df)

full.joint.segs.w.phasing.and.raw.df <- addOriginalCPNValuesAndCallsVsPloidyToOutputSegmentation(sample.vec,
                                        seg.name.vec = names(all.pw.cn.het.comp.lists[[1]]), 
                                        processed.orig.cpn.segs.df,
                                        full.joint.segs.w.phasing.df)


#Todo plot how much of the original segmentation actually comes all the way through the code.
full.joint.segs.df <- full.joint.segs.w.phasing.and.raw.df

full.joint.segs.df.path <- paste0(PLOT_OUTPUT_DIR, patient, "_full_joint_segs_df.RDS")
saveRDS(full.joint.segs.df, file=full.joint.segs.df.path)


png(file = paste(PLOT_OUTPUT_DIR, patient, "_CCF_cpnA", ".png", sep = ""), width = 1200, height = 1400)

plotAllSampleCCFsAcrossGenome(all.pw.cn.het.comp.lists = all.pw.cn.het.comp.lists,
                              full.joint.segs.df = full.joint.segs.w.phasing.df,
                              sample.vec = sample.vec,
                              all.seg.CCF.vals.list = all.seg.CCF.list,
                              cpnX = "cpnA")

dev.off()

png(file = paste(PLOT_OUTPUT_DIR, patient, "_CCF_cpnB", ".png", sep = ""), width = 1200, height = 1400)

plotAllSampleCCFsAcrossGenome(all.pw.cn.het.comp.lists = all.pw.cn.het.comp.lists,
                              full.joint.segs.df = full.joint.segs.w.phasing.df,
                              sample.vec = sample.vec,
                              all.seg.CCF.vals.list = all.seg.CCF.list,
                              cpnX = "cpnB")
dev.off()

#Now let's plot baf
png(file = paste(PLOT_OUTPUT_DIR, patient, "_BAF_phased", ".png", sep=""), width = 1200, height = 1400)

plotAllSampleBAFsAcrossGenome(sample.vec, patient.phased.BAF.df.list, cpnA.allele.color = '#de2d2699', cpnB.allele.color= '#3182bd99')

dev.off()

#Now for these sample by sample comparisons the names of PNGs should contain the two samples being being compared really so we'll have to the plotting within the functions

plotOldNewComparisonOfAllCPN(sample.vec, full.joint.segs.df, plot.output.dir = PLOT_OUTPUT_DIR)


plotAllSamplePairwiseComparisonSegmentCPN(sample.vec, full.joint.segs.df, plot.output.dir = PLOT_OUTPUT_DIR)

testFinalOutputDataframe(full.joint.segs.df)

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        