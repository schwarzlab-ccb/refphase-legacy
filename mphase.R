#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
patient <- args[1]
allow.dip.test.phasing.flag <- as.logical(args[2])
output.directory.base.path <- args[3]
data.format.to.use <- args[4]#Default for this is MEDICC, but can specify TRACERx_exome, TRACERx_snp_array or TRACERx_WGS
snp.file.path <- args[5]#Points towards the snp files, of whatever format they are given that should be used based on data.format to use.
seg.file.path <- args[6]#Pintts towards the seg files, of whatever format they should be determined byh the data.format.to.use flag.
overwrite.flag  <- as.logical(args[7])
gamma.val <- args[8]


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

# #Debugging option parameters:
# patient <- "MET073"
# allow.dip.test.phasing.flag  <- FALSE
# output.directory.base.path <- paste("/home/camp/watkint/working/watkint/data/brastianos/release_v2/","B_", patient, "/exome/CN_het_TRACERx_no_diptest/", sep = "")
# data.format.to.use <- "TRACERx_exome"
# snp.file.path <- paste0(output.directory.base.path, patient, "_TRACERx_input_paths_df.txt")
# seg.file.path <- paste("/home/camp/watkint/working/watkint/data/brastianos/release_v2/", "B_", patient, "/exome/MSAI/", patient, "_joint_segs_all_ascat_info.txt", sep = "")
# overwrite.flag <- FALSE
# gamma.val <- 1
# ##########
#mphase.script.path <- "/home/camp/watkint/working/watkint/scripts/multiregioncopynum/mphase.R"
#TRACERx format for TRACERx example
#tracerx.mphase.command <- paste("Rscript", mphase.script.path, patient, allow.dip.test.phasing.flag, data.format.to.use, snp.file.path, seg.file.path, sep= " ")

# #Debugging option parameters:
# patient <- "ETX014"
# allow.dip.test.phasing.flag  <- FALSE
# # M_ETX014/exome/CN_het_MEDICC_no_diptest
# output.directory.base.path <- paste("/home/camp/watkint/working/watkint/data/murugaesu/release/","M_", patient, "/exome/CN_het_TRACERx_no_diptest/", sep = "")
# data.format.to.use <- "TRACERx_exome"
# snp.file.path <- paste0(output.directory.base.path, patient, "_TRACERx_input_paths_df.txt")
# seg.file.path <- paste("/home/camp/watkint/working/watkint/data/murugaesu/release/", "M_", patient, "/exome/MSAI/", patient, "_joint_segs_all_ascat_info.txt", sep = "")
# ##########
# mphase.script.path <- "/home/camp/watkint/working/watkint/scripts/multiregioncopynum/mphase.R"
# #TRACERx format for TRACERx example
# tracerx.mphase.command <- paste("Rscript", mphase.script.path, patient, allow.dip.test.phasing.flag, data.format.to.use, snp.file.path, seg.file.path, sep= " ")



# # #Debugging option parameters:
# # patient <- "MET048"
# allow.dip.test.phasing.flag  <- FALSE
# output.directory.base.path <- "/Users/watkint/Desktop/MET048_mPhase_data/CN_het_output/"
# data.format.to.use <- "MEDICC"
# snp.file.path <- "/Users/watkint/Desktop/MET048_mPhase_data/MET048_joint_snps.txt"
# seg.file.path <- "/Users/watkint/Desktop/MET048_mPhase_data/MET048_joint_segs_all_ascat_info.txt"
# ##########
# mphase.script.path <- "/Users/watkint/code/multiregioncopynum/mphase.R"
# #MEDICC format for TRACERx example
# medicc.mphase.command <- paste("Rscript", mphase.script.path, patient, allow.dip.test.phasing.flag, data.format.to.use, snp.file.path, seg.file.path, sep=" ")


# patient <- "PTX001"
# allow.dip.test.phasing.flag <- FALSE
# output.directory.base.path <- "/home/camp/watkint/working/watkint/data/gundem/release/G_PTX001/WGS/CN_het_MEDICC_no_diptest/"
# data.format.to.use <- "MEDICC"
# snp.file.path <- "/home/camp/watkint/working/watkint/data/gundem/release/G_PTX001/WGS/MSAI/PTX001_joint_snps.txt"
# seg.file.path <- "/home/camp/watkint/working/watkint/data/gundem/release/G_PTX001/WGS/MSAI/PTX001_joint_segs_all_ascat_info.txt"


# patient <- "BTX001"
# allow.dip.test.phasing.flag <- FALSE
# output.directory.base.path <- "/home/camp/watkint/working/watkint/data/yates/release/Y_BTX001/WGS//CN_het_TRACERx_WGS_no_diptest/"
# data.format.to.use <- "TRACERx_WGS"
# snp.file.path <- "/home/camp/watkint/working/watkint/data/yates/release/Y_BTX001/WGS//CN_het_TRACERx_WGS_no_diptest/Y_BTX001_input_snp_files.txt"
# seg.file.path <- "/home/camp/watkint/working/watkint/data/yates/release/Y_BTX001/WGS/MSAI/BTX001_joint_segs_all_ascat_info.txt"
#overwrite.flag <- FALSE


# patient <- "KTX167"
# allow.dip.test.phasing.flag <- FALSE
# output.directory.base.path <-  "/home/camp/watkint/working/watkint/data/turajlic_wgs/release/T_KTX167/WGS//CN_het_MEDICC_no_diptest/"
# data.format.to.use <-  "MEDICC"
# snp.file.path <- "/home/camp/watkint/working/watkint/data/turajlic_wgs/release/T_KTX167/WGS/MSAI/KTX167_joint_snps.txt"
# seg.file.path <- "/home/camp/watkint/working/watkint/data/turajlic_wgs/release/T_KTX167/WGS/MSAI/KTX167_joint_segs_all_ascat_info.txt"

# #"patient"
# patient <- "MTX005"
# allow.dip.test.phasing.flag <- FALSE
# output.directory.base.path <- paste0("/home/camp/watkint/working/watkint/data/martinez/release_v3/R_",patient,"/snp_array//CN_het_MEDICC_no_diptest/")
# data.format.to.use <- "MEDICC"
# snp.file.path <- paste0("/home/camp/watkint/working/watkint/data/martinez/release_v3/R_",patient,"/snp_array/MSAI/",patient,"_new_all_snps.txt")
# seg.file.path <- paste0("/home/camp/watkint/working/watkint/data/martinez/release_v3/R_",patient,"/snp_array/MSAI/",patient,"_joint_segs_all_ascat_info.txt")
# overwrite.flag <- TRUE
# gamma.val <- 0.5



# patient <- "PTX002"
# allow.dip.test.phasing.flag <- FALSE
# output.directory.base.path <- "/home/camp/watkint/working/watkint/data/gundem/release/G_PTX002/WGS/CN_het_TRACERx_WGS_no_diptest/"
# data.format.to.use <- "TRACERx_WGS"
# snp.file.path <- "/home/camp/watkint/working/watkint/data/gundem/release/G_PTX002/WGS/CN_het_TRACERx_WGS_no_diptest/G_PTX002_input_snp_files.txt"
# seg.file.path <- "/home/camp/watkint/working/watkint/data/gundem/release/G_PTX002/WGS/MSAI/PTX002_joint_segs_all_ascat_info.txt"
# overwrite.flag <- FALSE
# gamma.val <- 1


# patient  <-  "STX009"
# allow.dip.test.phasing.flag <- FALSE
# output.directory.base.path <- "/home/camp/watkint/working/watkint/data/sottoriva_snp/release_v3/C_STX009/snp_array//CN_het_MEDICC_no_diptest/"
# data.format.to.use <- "MEDICC"
# snp.file.path <- "/home/camp/watkint/working/watkint/data/sottoriva_snp/release_v3/C_STX009/snp_array/MSAI/STX009_final_all_snps.txt"
# seg.file.path <- "/home/camp/watkint/working/watkint/data/sottoriva_snp/release_v3/C_STX009/snp_array/MSAI/STX009_joint_segs_all_ascat_info.txt"
# overwrite.flag <- FALSE
# gamma.val <- 0.5


# patient  <-  "STX009"
# allow.dip.test.phasing.flag <- FALSE
# output.directory.base.path <- "/home/camp/watkint/working/watkint/data/sottoriva_snp/release_v3/C_STX009/snp_array//CN_het_MEDICC_no_diptest/"
# data.format.to.use <- "MEDICC"
# snp.file.path <- "/home/camp/watkint/working/watkint/data/sottoriva_snp/release_v3/C_STX009/snp_array/MSAI/STX009_final_all_snps.txt"
# seg.file.path <- "/home/camp/watkint/working/watkint/data/sottoriva_snp/release_v3/C_STX009/snp_array/MSAI/STX009_joint_segs_all_ascat_info.txt"
# overwrite.flag <- FALSE
# gamma.val <- 0.5

patient  <-  "LTX021"
allow.dip.test.phasing.flag <- FALSE
output.directory.base.path <- "/camp/project/proj-tracerx-lung/tctProjects/watkint/TRACERx100_lung_pan_can_analysis/release/A_LTX021/exome/CN_het_TRACERx_no_diptest/"
data.format.to.use <- "TRACERx_exome"
snp.file.path <- "/camp/project/proj-tracerx-lung/tctProjects/watkint/TRACERx100_lung_pan_can_analysis/release/A_LTX021/exome/CN_het_TRACERx_no_diptest/LTX021_TRACERx_input_paths_df.txt"
seg.file.path <- "/camp/project/proj-tracerx-lung/tctProjects/watkint/TRACERx100_lung_pan_can_analysis/release/A_LTX021/exome/MSAI/LTX021_joint_segs_all_ascat_info.txt"
overwrite.flag <- TRUE
gamma.val <- 1

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
print("overwrite.flag")
print(overwrite.flag)
print("gamma.val")
print(gamma.val)


#We need to be able to select the intersection of all regions or a full GISTIC compatible union with defined regions that will always be covered
#What proportion of the chosen genome, in this case hg19, arms are covered by the segments that have been inputted.
#We also need to add in the wgii calculations


###Source relevant scripts, best to do this at the top to check for errors...
source(paste0(current.dir, "/utility_MPhase_functions.R"))
source(paste0(current.dir, "/data_input_and_QC_functions_v2.R"))
source(paste0(current.dir, "/minimum_consistent_segmentation_functions.R"))
source(paste0(current.dir, "/consensus_phasing_from_multiple_samples_functions_v2.R"))
source(paste0(current.dir, "/generation_of_phased_and_total_cpn_functions_v2.R"))
source(paste0(current.dir, "/MSAI_detection_functions.R"))
source(paste0(current.dir, "/pairwise_cpn_comparison_and_heterogeneity_functions.R"))
source(paste0(current.dir, "/cancer_cell_fraction_calculation_functions_v2.R"))
source(paste0(current.dir, "/final_output_segmentation_file_generation_functions_v3.R"))
source(paste0(current.dir, "/segment_joining_functions.R"))
source(paste0(current.dir, "/cpn_heterogeneity_plotting_functions_v2.R"))

#####Get relevant libraries:
#Todo check this is the only required library.
library(diptest)
library(gplots)
library(gdata)
library(GenomicRanges)

###############This is past the input stage and into the analysis:
#Set constants for all analyses here... eventually these can be command line set
MIN_NUM_HZ_SNP_PHASING <- 5
#PLOT_OUTPUT_DIR <- NULL
OVERWRITE_FLAG <- as.logical(overwrite.flag)
PLOT_OUTPUT_DIR <- paste(output.directory.base.path, "/", sep="")
if ( ! file.exists(PLOT_OUTPUT_DIR) ) {

	print(PLOT_OUTPUT_DIR)
	print("directory doesn't exist")

}
GAMMA_VAL <- as.numeric(gamma.val)
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
PW_CPNA_CPNB_DIFF_AMPLITUDE_THRESH <- 0.001
PVAL_THRESH_TO_JOIN_PH_SEGS <- 0.001
MIN_PROP_SNP_SEP_THRESH <- 0.75
PVAL_COMB_WILCOX_PH_THRESH <- 0.01
PVAL_DIPTEST_AI_THRESH <- 0.01
ALLOW_DIPTEST_ALONE_PH_FLAG <- allow.dip.test.phasing.flag
MSAI_BAF_LTNPD_PVAL_THRESHOLD <- 0.001
MSAI_CPNX_PVAL_THRESHOLD <- 0.001
PVAL_THRESH_CONVERT_GAIN_TO_AMP  <- 0.01

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



###Load in the relevant genome file, #todo make this an option that can be set.
hg19.info.path <- paste0(current.dir, "/", "ucsc_hg19_chrominfo.txt")
hg19.info.df <- read.table(hg19.info.path, stringsAsFactors = FALSE, sep = "\t", header = TRUE)





#Paths for saving/loading objects depending on OVERWRITE_FLAG
#Reading in data object paths
asas.SNPS.lr.and.or.baf.lists.path <- paste0(PLOT_OUTPUT_DIR, patient, "_asas_SNPS_lr_and_or_baf_lists_lists.RDS")
asas.SNPs.w.both.lr.baf.snp.lists.path <- paste0(PLOT_OUTPUT_DIR, patient, "_asas_SNPs_w_both_lr_baf_snp_lists.RDS")
min.cons.segs.df.path <- paste0(PLOT_OUTPUT_DIR, patient, "_min_cons_segs_df.RDS")
processed.orig.cpn.segs.df.path <- paste0(PLOT_OUTPUT_DIR, patient, "_processed_orig_cpn_segs_df.RDS")

ploidy.est.vec.path <- paste0(PLOT_OUTPUT_DIR, patient, "_ploidy_est_vec.RDS")
purity.est.vec.path <- paste0(PLOT_OUTPUT_DIR, patient, "_purity_est_vec.RDS")
#Analysis object paths
all.seg.phasing.info.list.path <- paste0(PLOT_OUTPUT_DIR, patient, "_all_seg_phasing_info_list.RDS")
patient.phased.BAF.df.list.path <- paste0(PLOT_OUTPUT_DIR, patient, "_phased_BAF_df_list.RDS")
all.seg.MSAI.by.BAF.ltnpd.present.list.path <- paste0(PLOT_OUTPUT_DIR, patient, "_all_seg_MSAI_by_BAF_ltnpd_present_list.RDS")
all.samp.non.hz.logr.and.tot.cpn.df.list.path <- paste0(PLOT_OUTPUT_DIR, patient, "_all_samp_non_hz_logr_and_tot_cpn_df_list.RDS")
asas.ph.raw.cpn.lists.path <- paste0(PLOT_OUTPUT_DIR, patient, "_asas_ph_raw_cpn_lists.RDS")
samps.int.cpn.raw.cpn.by.chr.list.path <- paste0(PLOT_OUTPUT_DIR, patient, "_samps_int_cpn_raw_cpn_by_chr_list.RDS")
all.seg.CCF.list.path <- paste0(PLOT_OUTPUT_DIR, patient, "_all_seg_CCF_list.RDS")
all.seg.int.ob.xp.cpn.df.list.path <- paste0(PLOT_OUTPUT_DIR, patient, "_all_seg_int_ob_xp_cpn_df_list.RDS")
all.seg.CCF.diff.from.exp.list.path <- paste0(PLOT_OUTPUT_DIR, patient, "_all_seg_CCF_diff_from_exp_list.RDS")
all.seg.clonal.in.all.samps.list.path <- paste0(PLOT_OUTPUT_DIR, patient, "_all_seg_clonal_in_all_samps_list.RDS")
all.seg.greater.than.two.diff.list.path <- paste0(PLOT_OUTPUT_DIR, patient, "_all_seg_greater_than_two_diff_list.RDS")

#Load in the data: 
#Note that the afunctions to create the single dataframes that it accepts are in the data_input_and_QC_functions.R
if (OVERWRITE_FLAG == TRUE) {

	if (data.format.to.use == "MEDICC") {

		snp.file.df <- read.table(snp.file.path, stringsAsFactors = FALSE, sep = "\t", header = TRUE)
		seg.file.df <- read.table(seg.file.path, stringsAsFactors = FALSE, sep = "\t", header = TRUE)
		seg.file.df$sample <- gsub("\\.", "-", gsub("_LogR", "", seg.file.df$sample))
		seg.file.df$callState <- NULL
		seg.file.df$patient <- patient
		seg.file.df$sample_name <- paste(patient, sapply(seg.file.df$sample, function(x) strsplit(x, split = "-")[[1]][2]), sep = ":")
		processed.orig.cpn.segs.df <- errorCheckAndAddPloidyAdjRawCPNToSegDF(patient, seg.file.df)
		sample.vec <- unique(processed.orig.cpn.segs.df$sample)
		purity.est.vec <- processed.orig.cpn.segs.df[ ! duplicated(processed.orig.cpn.segs.df$sample), "ACF"]
		names(purity.est.vec) <- sample.vec
		ploidy.est.vec <- processed.orig.cpn.segs.df[ ! duplicated(processed.orig.cpn.segs.df$sample), "Ploidy"]
		names(ploidy.est.vec) <- sample.vec
		min.cons.segs.df <- generateMinimumConsistentSegmentation(processed.orig.cpn.segs.df, chroms = 1:22, col.to.use = "ploidy_adj_tot_cpn")
		min.cons.segs.df <- cleanAndErrorCheckMinConsSegs(min.cons.segs.df)
		colnames(min.cons.segs.df)[1:3] <- c("chr", "startpos", "endpos")
		all.samp.SNP.list <- list()

		for ( samp in sample.vec ) {

			print(samp)
			print(head(snp.file.df))

			ASCAT.summary <- data.frame(chr = as.numeric(snp.file.df$Chr),
				                		pos = as.numeric(snp.file.df$Pos),
				   						BAF = as.numeric(snp.file.df[, paste0(gsub("-" , "\\.", samp), "_BAF")]),
				   						BAFsegmented = rep(NA, length(snp.file.df$Chr)),
				   						logR = as.numeric(snp.file.df[, paste0(gsub("-", "\\.", samp), "_LogR")]),
				   						logRsegmented = rep(NA, length(snp.file.df$Chr)),
				   						stringsAsFactors = FALSE)

			print(head(ASCAT.summary))

			seg <- processed.orig.cpn.segs.df[processed.orig.cpn.segs.df$sample == samp,]

			for ( j in 1:nrow(seg) ) {
				#print(seg[j,])
			  if (nrow(ASCAT.summary[ASCAT.summary$chr == seg[j ,]$chr
			  						& ASCAT.summary$pos >= seg[j ,]$startpos 
			  						& ASCAT.summary$pos <= seg[j ,]$endpos 
			  						& ASCAT.summary$BAF != 2 ,]) > 0) {
			  	
			    if (seg[j ,]$nMajor == seg[j ,]$nMinor) {

			      ASCAT.summary[ASCAT.summary$chr == seg[j ,]$chr
			      			  & ASCAT.summary$pos >= seg[j ,]$startpos
			      			  & ASCAT.summary$pos <= seg[j ,]$endpos,]$BAFsegmented <- 0.5
			      
			    } else {
		
			      avg.seg.baf <- 0.5 - mean(abs(as.numeric(ASCAT.summary[ASCAT.summary$chr == seg[j ,]$chr
			      	 													& ASCAT.summary$pos >= seg[j ,]$startpos
			      	 													& ASCAT.summary$pos <= seg[j ,]$endpos
			      	 													& ASCAT.summary$BAF != 2 ,]$BAF)          - 0.5))

			      ASCAT.summary[ASCAT.summary$chr == seg[j ,]$chr 
			      			  & ASCAT.summary$pos >= seg[j ,]$startpos 
			      			  & ASCAT.summary$pos <= seg[j ,]$endpos ,]$BAFsegmented <- avg.seg.baf

			    }
		 
			  }#end if(nrow(ASCAT.summary[ASCAT.summary$chrom==seg[j,]$chr & ASCAT.summary$pos...

			  if (nrow(ASCAT.summary[ASCAT.summary$chr == seg[j ,]$chr
			                       & ASCAT.summary$pos >= seg[j ,]$startpos
			                       & ASCAT.summary$pos <= seg[j ,]$endpos ,]) > 0) {

			  	avg.seg.logr <- mean(ASCAT.summary[ASCAT.summary$chr == seg[j ,]$chr
			  					   & ASCAT.summary$pos >= seg[j ,]$startpos
			  					   & ASCAT.summary$pos <= seg[j ,]$endpos, ]$logR)

			  	ASCAT.summary[ASCAT.summary$chr == seg[j ,]$chr
			  			    & ASCAT.summary$pos >= seg[j ,]$startpos
			  			    & ASCAT.summary$pos <= seg[j ,]$endpos ,]$logRsegmented <- avg.seg.logr

			  }

			}#end for ( j in 1:nrow(seg) ) {

			all.samp.SNP.list[[samp]] <- ASCAT.summary

		}#end for ( samp in sample.vec ) {



	} else if (data.format.to.use == "TRACERx_exome" | data.format.to.use == "TRACERx_snp_array") {

		f.joint.segs.df <- read.table(seg.file.path, stringsAsFactors = FALSE, sep = "\t", header = TRUE)
		f.joint.segs.df$sample <- gsub("\\.", "-", gsub("_LogR", "", f.joint.segs.df$sample))
		f.joint.segs.df$callState <- NULL
		f.joint.segs.df$patient <- patient
		f.joint.segs.df$sample_name <- paste(patient, sapply(f.joint.segs.df$sample, function(x) strsplit(x, split = "-")[[1]][2]), sep = ":")

		# gr.f.joint.segs.df <- f.joint.segs.df[,c("chr", "startpos", "endpos")]
		# colnames(gr.f.joint.segs.df) <- c("Chromosome", "start", "end")

		# gr.f.joint.segs.gr <- makeGRangesFromDataFrame(gr.f.joint.segs.df)
		# gr.f.overlap.segs.gr <- countOverlaps(gr.f.joint.segs.gr)

		# dj.gr.joint.segs.gr <- disjoin(gr.f.joint.segs.gr)
		# joint.segs.df <- as.data.frame(dj.gr.joint.segs.gr)


		# #So now we want to match up the old values...

		# seg.test.df <- data.frame(Chromosome = chr.start.end.vec[1],
		#           start = as.numeric(chr.start.end.vec[2]),
		#           end = as.numeric(chr.start.end.vec[3]),
		#           stringsAsFactors = FALSE)
		# test.gen.range <- makeGRangesFromDataFrame(seg.test.df)

		if (data.format.to.use == "TRACERx_exome") {

			array.type <- "exome"

		} else if (data.format.to.use == "TRACERx_snp_array") {

			array.type <- "snp_array"
		}

		orig.cpn.segs.df <- f.joint.segs.df
		processed.orig.cpn.segs.df <- errorCheckAndAddPloidyAdjRawCPNToSegDF(patient, orig.cpn.segs.df)
		ploidy.est.vec <- processed.orig.cpn.segs.df[ ! duplicated(processed.orig.cpn.segs.df$sample), 'Ploidy']
		names(ploidy.est.vec) <- unique(processed.orig.cpn.segs.df$sample)
		purity.est.vec <- processed.orig.cpn.segs.df[ ! duplicated(processed.orig.cpn.segs.df$sample), 'ACF']
		names(purity.est.vec) <- unique(processed.orig.cpn.segs.df$sample)
		min.cons.segs.df <- generateMinimumConsistentSegmentation(processed.orig.cpn.segs.df, chroms = 1:22, col.to.use = "ploidy_adj_tot_cpn")
		min.cons.segs.df <- cleanAndErrorCheckMinConsSegs(min.cons.segs.df)
		#So we want to read in a file that looks like this:	#patient_id sample_id ascat_bc_file_path 
		pt.snp.file.paths.df <- read.table(snp.file.path, stringsAsFactors = FALSE, sep = "\t", header = TRUE)
		all.samp.SNP.list <- getSNPDataForTRACERxExomeSNPArray(pt.snp.file.paths.df)
		#all sample all seg logr baf lists is asas.lr.baf, All SNPs in these lists must have both a logr and baf value
		asas.SNPs.w.both.lr.baf.snp.lists <- cleanAndErrorCheckSNPs(all.samp.SNP.list, min.cons.segs.df, allow.logr.only.flag = FALSE)
		#SNPs that only have LogR are permitted in these lists, but it also contains all the SNPs with both BAF and LogR
		asas.SNPS.lr.and.or.baf.lists <- cleanAndErrorCheckSNPs(all.samp.SNP.list, min.cons.segs.df, allow.logr.only.flag = TRUE)


	} else if (data.format.to.use == "TRACERx_WGS") {

		f.joint.segs.df <- read.table(seg.file.path, stringsAsFactors = FALSE, sep = "\t", header = TRUE)
		f.joint.segs.df$sample <- gsub("\\.", "-", gsub("_LogR", "", f.joint.segs.df$sample))
		f.joint.segs.df$callState <- NULL
		f.joint.segs.df$patient <- patient
		f.joint.segs.df$sample_name <- paste(patient,sapply(f.joint.segs.df$sample, function(x) strsplit(x, split = "-")[[1]][2]), sep = ":")


		array.type <- "WGS"
		orig.cpn.segs.df <- f.joint.segs.df
		processed.orig.cpn.segs.df <- errorCheckAndAddPloidyAdjRawCPNToSegDF(patient, orig.cpn.segs.df)
		ploidy.est.vec <- processed.orig.cpn.segs.df[ ! duplicated(processed.orig.cpn.segs.df$sample), 'Ploidy']
		names(ploidy.est.vec) <- unique(processed.orig.cpn.segs.df$sample)
		purity.est.vec <- processed.orig.cpn.segs.df[ ! duplicated(processed.orig.cpn.segs.df$sample), 'ACF']
		names(purity.est.vec) <- unique(processed.orig.cpn.segs.df$sample)
		min.cons.segs.df <- generateMinimumConsistentSegmentation(processed.orig.cpn.segs.df, chroms = 1:22, col.to.use = "ploidy_adj_tot_cpn")
		min.cons.segs.df <- cleanAndErrorCheckMinConsSegs(min.cons.segs.df)
		all.samp.SNP.list <- getSNPDataForTRACERxWGS(snp.file.path)


	} else {

		print(data.format.to.use)
		stop("Incorrect data.format.to.use input to mphase")

	}


	chr22.j.segs.df <- f.joint.segs.df[f.joint.segs.df$chr==22 & f.joint.segs.df$sample%in%c("A_LTX021_SU_T1-R1", "A_LTX021_SU_T1-R2"),]
	f.joint.seg.grs <- makeGRangesFromDataFrame(chr22.j.segs.df,
				                         keep.extra.columns=FALSE,
				                         ignore.strand=FALSE,
				                         seqinfo=NULL,
				                         seqnames.field=c("seqnames", "seqname",
				                                          "chromosome", "chrom",
				                                          "chr", "chromosome_name",
				                                          "seqid"),
				                         start.field=c("Start_Position","start", "startpos"),
				                         end.field=c("End_Position", "endpos", "end", "stop"),
				                         strand.field="strand",
				                         starts.in.df.are.0based=FALSE)


getChromosomeArmGrange <- function(chr, chr.arm, chrom.info.df) {

	if (chr.arm == "p") {

		arm.df <- data.frame(
		chrom  = chr,
		start = chrom.info.df[chrom.info.df$chrom == chr,]$p_start,
		end = chrom.info.df[chrom.info.df$chrom == chr,]$p_end,
		stringsAsFactors = FALSE)

	} else if (chr.arm == "q") {

		arm.df <- data.frame(
		chrom  = chr,
		start = chrom.info.df[chrom.info.df$chrom ==chr,]$q_start,
		end = chrom.info.df[chrom.info.df$chrom == chr,]$q_end,
		stringsAsFactors = FALSE)

	}

	arm.gr <- makeGRangesFromDataFrame(arm.df,
		                        keep.extra.columns=FALSE,
		                        ignore.strand=FALSE,
		                        seqinfo=NULL,
		                        seqnames.field=c("seqnames", "seqname",
		                                         "chromosome", "chrom",
		                                         "chr", "chromosome_name",
		                                         "seqid"),
		                        start.field=c("Start_Position","start"),
		                        end.field=c("End_Position", "end", "stop"),
		                        strand.field="strand",
		                        starts.in.df.are.0based=FALSE)

	return(arm.gr)

}



##
##Params:
##
##Returns:
##
compareOverlapsBetweenSamples <- function(sample.vec, f.joint.segs.df, hg19.info.df, chrs.to.use.vec) {

	#For each sample how much of the chromosome arm
	chrs.to.use.vec <- 1:22
	#Subset the chromos
	chr.info.df <- hg19.info.df[hg19.info.df$chrom %in% chrs.to.use.vec,]
	p.name.vec <- chr.info.df$p_name
	q.name.vec <- chr.info.df$q_name
	arm.vec <- c(rbind(p.name.vec, q.name.vec)) 
	p.q.vec <- c("p", "q")

	#So now let's make an arm coverage data.frame
	arm.coverage.df <- data.frame(matrix(vector(), length(arm.vec), (length(sample.vec) + 2),
	                        dimnames = list(c(), c("arm", sample.vec, "overall_intersection"))),
	                        stringsAsFactors = FALSE,
	                        check.names = FALSE)
	arm.coverage.df$arm <- arm.vec
	#how 


	#Now what about intersection of all samples what proportion of the genome to they cover:
	samp.gr.list <- list()
	overall.intersection.gr <- NULL
	for ( samp in sample.vec ) {

		samp.gr <- makeGRangesFromDataFrame(f.joint.segs.df[f.joint.segs.df$sample==samp 
										  ,],
		                        keep.extra.columns=FALSE,
		                        ignore.strand=FALSE,
		                        seqinfo=NULL,
		                        seqnames.field=c("seqnames", "seqname",
		                                         "chromosome", "chrom",
		                                         "chr", "chromosome_name",
		                                         "seqid"),
		                        start.field=c("Start_Position","start","startpos"),
		                        end.field=c("End_Position", "end", "stop", "endpos"),
		                        strand.field="strand",
		                        starts.in.df.are.0based=FALSE)
		samp.gr.list[[samp]] <- samp.gr

		if ( is.null(overall.intersection.gr) ) {

			overall.intersection.gr <- samp.gr

		} else {

			overall.intersection.gr <- intersect(overall.intersection.gr, samp.gr)

		}#end ! is.null (overall.intersection.gr)

	}#en for ( samp in sample.vec ) {
	samp.gr.list[["overall_intersection"]] <- overall.intersection.gr
	#So now we want to take the intersection and compare the sizes of the

	for ( samp in names(samp.gr.list) ) {

		#Get sample granges:
		samp.gr <- samp.gr.list[[samp]]
		print(samp)
		print(samp.gr)
		for ( chr in chrs.to.use.vec ) {
			for ( arm in p.q.vec) {

				chr.arm.gr <- getChromosomeArmGrange(chr, arm, chr.info.df)
				chr.arm.overlap.gr <- intersect(samp.gr, chr.arm.gr)
				arm.coverage.df[arm.coverage.df$arm==paste0(chr, arm), samp] <- sum(as.numeric(width(chr.arm.overlap.gr))) / as.numeric(width(chr.arm.gr))

			}#end for ( arm in p.q.vec) {
		}# for ( chr in chrs.to.use.vec ) {

	}#end for ( samp in sample.vec ) {

	arm.coverage.df$diff_from_max <-  apply(arm.coverage.df[,2:ncol(arm.coverage.df)], 1, max) - arm.coverage.df$overall_intersection

	#Report max and union

	return(arm.coverage.df)
}

test.df <- data.frame(chr=c(1,1,1),
			start=c(1,2,5),
			end=c(3,6,8),
			stringsAsFactors=FALSE)

test.gr <- makeGRangesFromDataFrame(test.df,
                         keep.extra.columns=FALSE,
                         ignore.strand=FALSE,
                         seqinfo=NULL,
                         seqnames.field=c("seqnames", "seqname",
                                          "chromosome", "chrom",
                                          "chr", "chromosome_name",
                                          "seqid"),
                         start.field=c("Start_Position","start"),
                         end.field=c("End_Position", "end", "stop"),
                         strand.field="strand",
                         starts.in.df.are.0based=FALSE)


f.joint.segs.gr <- makeGRangesFromDataFrame(f.joint.segs.df,
		                         keep.extra.columns=FALSE,
		                         ignore.strand=FALSE,
		                         seqinfo=NULL,
		                         seqnames.field=c("seqnames", "seqname",
		                                          "chromosome", "chrom",
		                                          "chr", "chromosome_name",
		                                          "seqid"),
		                         start.field=c("Start_Position","start", "startpos"),
		                         end.field=c("End_Position", "end", "stop", "endpos"),
		                         strand.field="strand",
		                         starts.in.df.are.0based=FALSE)
disjoin(f.joint.segs.gr)

##
##Params:
##
##Returns:
##
getChromosomeArmsGRangeList <- function(current.chr, chrom.info.df) {

	arm.df.list <- list() 
	p.arm.df <- data.frame(
	chrom  = current.chr,
	start = chrom.info.df[chrom.info.df$chrom == current.chr,]$p_start,
	end = chrom.info.df[chrom.info.df$chrom == current.chr,]$p_end,
	stringsAsFactors = FALSE)

	count <- 1

	if(chrom.info.df[chrom.info.df$chrom == current.chr,]$acrocentric == FALSE) {

		arm.df.list[[count]] <- makeGRangesFromDataFrame(p.arm.df,
		                         keep.extra.columns=FALSE,
		                         ignore.strand=FALSE,
		                         seqinfo=NULL,
		                         seqnames.field=c("seqnames", "seqname",
		                                          "chromosome", "chrom",
		                                          "chr", "chromosome_name",
		                                          "seqid"),
		                         start.field=c("Start_Position","start"),
		                         end.field=c("End_Position", "end", "stop"),
		                         strand.field="strand",
		                         starts.in.df.are.0based=FALSE)

		count <- count + 1
	}



	q.arm.df <- data.frame(
	chrom  = current.chr,
	start = chrom.info.df[chrom.info.df$chrom == current.chr,]$q_start,
	end = chrom.info.df[chrom.info.df$chrom == current.chr,]$q_end,
	stringsAsFactors = FALSE)

	arm.df.list[[count]] <- makeGRangesFromDataFrame(q.arm.df,
	                         keep.extra.columns=FALSE,
	                         ignore.strand=FALSE,
	                         seqinfo=NULL,
	                         seqnames.field=c("seqnames", "seqname",
	                                          "chromosome", "chrom",
	                                          "chr", "chromosome_name",
	                                          "seqid"),
	                         start.field=c("Start_Position","start"),
	                         end.field=c("End_Position", "end", "stop"),
	                         strand.field="strand",
	                         starts.in.df.are.0based=FALSE)

	return(arm.df.list)
}

all.chrs.matched.grange.df <- data.frame(Sample=character(0),
										Chromosome=character(0),
										Start_Position=numeric(0), 
										End_Position=numeric(0),
										Num_Markers=numeric(0),
										Seg.CN=numeric(0),
										stringsAsFactors=FALSE)

#We also want to calculate what proportion of all arms is covered by the original segmentation for each sample:
for ( chr in CHRS_TO_USE_VEC ) {

	#Get the chromosome arms in non-acrocentric
	# for ( current.chr in chr.vec ) {

	curr.pt.chr.grange <- makeGRangesFromDataFrame(all.avg.samp.gistic.df[all.avg.samp.gistic.df$Sample == pt
											&  all.avg.samp.gistic.df$Chromosome == current.chr,],
	                         keep.extra.columns=TRUE,
	                         ignore.strand=FALSE,
	                         seqinfo=NULL,
	                         seqnames.field=c("seqnames", "seqname",
	                                          "chromosome", "chrom",
	                                          "chr", "chromosome_name",
	                                          "seqid"),
	                         start.field=c("Start_Position","start"),
	                         end.field=c("End_Position", "end", "stop"),
	                         strand.field="strand",
	                         starts.in.df.are.0based=FALSE)
	#This is done to stop intersect merging segments that start and end right after each other:
	#E.g.
	#1 200 300
	#1 300 350
	#Will no longer be merged 
	alternating.idx.vec <- (seq(1, length(curr.pt.chr.grange)) %% 2) + 1#plus one to match 1-based indexing
	plus.minus.vec <- c("+", "-")
	strand(curr.pt.chr.grange) <- plus.minus.vec[alternating.idx.vec]
	#Get the chromosome arm corresponding 
	chr.gr.list <- getChromosomeArmsGRangeList(current.chr, chrom.info.df)

	#So do this for the chromosome p and q arms (if q arm is preesent)
	for ( i in 1:length(chr.gr.list) ) {

		curr.chr.arm.gr <- chr.gr.list[[i]]
		curr.chr.arm.gr.plus <- curr.chr.arm.gr
		curr.chr.arm.gr.minus <- curr.chr.arm.gr
		strand(curr.chr.arm.gr.plus) <- "+"
		strand(curr.chr.arm.gr.minus) <- "-"
		curr.chr.arm.gr.strands <- c(curr.chr.arm.gr.plus, curr.chr.arm.gr.minus)
		print("curr.chr.arm.gr.strands")
		print(curr.chr.arm.gr.strands)
		#intersect taking into account strand to prevent merging of adjacent segments:
		overlapping.gr <- intersect(curr.pt.chr.grange, curr.chr.arm.gr.strands, ignore.strand = FALSE)
		print("overlapping.gr")
		print(overlapping.gr)
		#Reset starand as it isn't used for antyhing else
		strand(overlapping.gr) <- "*"
		#No need to take into account strand as there is no way two segments would be merged by the setdiff operation.
		non.overlapping.gr <- setdiff(curr.chr.arm.gr, curr.pt.chr.grange, ignore.strand = TRUE)
		print("non.overlapping.gr")
		print(non.overlapping.gr)
		# print("non.overlapping.gr")
		# print(non.overlapping.gr)
		#Ok let's do this as dataframes:
		#hits.obj <- findOverlaps(overlapping.gr, curr.pt.chr.grange, type = "any")
		#print("hits.obj")
		#print(hits.obj)
		#idx.vec <- subjectHits(hits.obj)
		#So now we don't have to add the copy number as 
		#chr.seg.arm.info.df <- as.data.frame(DataFrame(Num_Markers = curr.pt.chr.grange$Num_Markers[idx.vec], Seg.CN = curr.pt.chr.grange$Seg.CN[idx.vec]))
		overlapping.segs.df <- as.data.frame(overlapping.gr)[, c("seqnames", "start", "end")]
		overlapping.segs.info.df <- cbind(as.data.frame(rep(pt, nrow(overlapping.segs.df))), overlapping.segs.df, chr.seg.arm.info.df)
		colnames(overlapping.segs.info.df) <- c("Sample", "Chromosome", "Start_Position", "End_Position", "Num_Markers", "Seg.CN")
		print("overlapping.segs.info.df")
		print(overlapping.segs.info.df)
		#Now we need to sort out the nonoverlapping:
		non.overlapping.segs.df <- as.data.frame(non.overlapping.gr)[, c("seqnames", "start", "end")]
		non.overlapping.segs.info.df <- cbind(as.data.frame(rep(pt, nrow(non.overlapping.segs.df))),
											  non.overlapping.segs.df,
											  as.data.frame(rep(5, nrow(non.overlapping.segs.df))),
											  as.data.frame(rep(0.15, nrow(non.overlapping.segs.df))))
		colnames(non.overlapping.segs.info.df) <- c("Sample", "Chromosome", "Start_Position", "End_Position", "Num_Markers", "Seg.CN")

		print("non.overlapping.segs.info.df")
		print(non.overlapping.segs.info.df)

		full.chrom.arm.info.df <- rbind(overlapping.segs.info.df, non.overlapping.segs.info.df)
		full.chrom.arm.info.df <- full.chrom.arm.info.df[with(full.chrom.arm.info.df, order(Chromosome, Start_Position, End_Position)),]

		#Ok so we've got a whole arm filled in..., add it to the full dataset:
		all.chrs.matched.grange.df <- rbind(all.chrs.matched.grange.df, full.chrom.arm.info.df)
		print(all.chrs.matched.grange.df)
	}	


}


	#Now let's make the inetersection:
	#So we want an empty dataframe to hold all the rows that we're going to create for each chromosome arm:
	all.chrs.matched.grange.df <- data.frame(Sample=character(0),
											Chromosome=character(0),
											Start_Position=numeric(0), 
											End_Position=numeric(0),
											Num_Markers=numeric(0),
											Seg.CN=numeric(0),
											stringsAsFactors=FALSE)


	# chr.vec <- 1:22
	# for ( current.chr in chr.vec ) {

	# 	curr.pt.chr.grange <- makeGRangesFromDataFrame(all.avg.samp.gistic.df[all.avg.samp.gistic.df$Sample == pt
	# 											&  all.avg.samp.gistic.df$Chromosome == current.chr,],
	# 	                         keep.extra.columns=TRUE,
	# 	                         ignore.strand=FALSE,
	# 	                         seqinfo=NULL,
	# 	                         seqnames.field=c("seqnames", "seqname",
	# 	                                          "chromosome", "chrom",
	# 	                                          "chr", "chromosome_name",
	# 	                                          "seqid"),
	# 	                         start.field=c("Start_Position","start"),
	# 	                         end.field=c("End_Position", "end", "stop"),
	# 	                         strand.field="strand",
	# 	                         starts.in.df.are.0based=FALSE)
	# 	#This is done to stop intersect merging segments that start and end right after each other:
	# 	#E.g.
	# 	#1 200 300
	# 	#1 300 350
	# 	#Will no longer be merged 
	# 	alternating.idx.vec <- (seq(1, length(curr.pt.chr.grange)) %% 2) + 1#plus one to match 1-based indexing
	# 	plus.minus.vec <- c("+", "-")
	# 	strand(curr.pt.chr.grange) <- plus.minus.vec[alternating.idx.vec]
	# 	#Get the chromosome arm corresponding 
	# 	chr.gr.list <- getChromosomeArmsGRangeList(current.chr, chrom.info.df)

	# 	#So do this for the chromosome p and q arms (if q arm is preesent)
	# 	for ( i in 1:length(chr.gr.list) ) {

	# 		curr.chr.arm.gr <- chr.gr.list[[i]]
	# 		curr.chr.arm.gr.plus <- curr.chr.arm.gr
	# 		curr.chr.arm.gr.minus <- curr.chr.arm.gr
	# 		strand(curr.chr.arm.gr.plus) <- "+"
	# 		strand(curr.chr.arm.gr.minus) <- "-"
	# 		curr.chr.arm.gr.strands <- c(curr.chr.arm.gr.plus, curr.chr.arm.gr.minus)
	# 		print("curr.chr.arm.gr.strands")
	# 		print(curr.chr.arm.gr.strands)
	# 		#intersect taking into account strand to prevent merging of adjacent segments:
	# 		overlapping.gr <- intersect(curr.pt.chr.grange, curr.chr.arm.gr.strands, ignore.strand = FALSE)
	# 		print("overlapping.gr")
	# 		print(overlapping.gr)
	# 		#Reset starand as it isn't used for antyhing else
	# 		strand(overlapping.gr) <- "*"
	# 		#No need to take into account strand as there is no way two segments would be merged by the setdiff operation.
	# 		non.overlapping.gr <- setdiff(curr.chr.arm.gr, curr.pt.chr.grange, ignore.strand = TRUE)
	# 		print("non.overlapping.gr")
	# 		print(non.overlapping.gr)
	# 		# print("non.overlapping.gr")
	# 		# print(non.overlapping.gr)
	# 		#Ok let's do this as dataframes:
	# 		#hits.obj <- findOverlaps(overlapping.gr, curr.pt.chr.grange, type = "any")
	# 		#print("hits.obj")
	# 		#print(hits.obj)
	# 		#idx.vec <- subjectHits(hits.obj)
	# 		#So now we don't have to add the copy number as 
	# 		#chr.seg.arm.info.df <- as.data.frame(DataFrame(Num_Markers = curr.pt.chr.grange$Num_Markers[idx.vec], Seg.CN = curr.pt.chr.grange$Seg.CN[idx.vec]))
	# 		overlapping.segs.df <- as.data.frame(overlapping.gr)[, c("seqnames", "start", "end")]
	# 		overlapping.segs.info.df <- cbind(as.data.frame(rep(pt, nrow(overlapping.segs.df))), overlapping.segs.df, chr.seg.arm.info.df)
	# 		colnames(overlapping.segs.info.df) <- c("Sample", "Chromosome", "Start_Position", "End_Position", "Num_Markers", "Seg.CN")
	# 		print("overlapping.segs.info.df")
	# 		print(overlapping.segs.info.df)
	# 		#Now we need to sort out the nonoverlapping:
	# 		non.overlapping.segs.df <- as.data.frame(non.overlapping.gr)[, c("seqnames", "start", "end")]
	# 		non.overlapping.segs.info.df <- cbind(as.data.frame(rep(pt, nrow(non.overlapping.segs.df))),
	# 											  non.overlapping.segs.df,
	# 											  as.data.frame(rep(5, nrow(non.overlapping.segs.df))),
	# 											  as.data.frame(rep(0.15, nrow(non.overlapping.segs.df))))
	# 		colnames(non.overlapping.segs.info.df) <- c("Sample", "Chromosome", "Start_Position", "End_Position", "Num_Markers", "Seg.CN")

	# 		print("non.overlapping.segs.info.df")
	# 		print(non.overlapping.segs.info.df)

	# 		full.chrom.arm.info.df <- rbind(overlapping.segs.info.df, non.overlapping.segs.info.df)
	# 		full.chrom.arm.info.df <- full.chrom.arm.info.df[with(full.chrom.arm.info.df, order(Chromosome, Start_Position, End_Position)),]

	# 		#Ok so we've got a whole arm filled in..., add it to the full dataset:
	# 		all.chrs.matched.grange.df <- rbind(all.chrs.matched.grange.df, full.chrom.arm.info.df)
	# 		print(all.chrs.matched.grange.df)
	# 	}	

	# }










	#Todo work out the differences here:
	#All SNPs in these lists must have both a logr and baf value: all sample all seg logr baf lists is asas.lr.baf
	asas.SNPs.w.both.lr.baf.snp.lists <- cleanAndErrorCheckSNPs(all.samp.SNP.list, min.cons.segs.df, allow.logr.only.flag = FALSE)
	#SNPs that only have LogR are permitted in these lists, but it also contains all the SNPs with both BAF and LogR
	asas.SNPS.lr.and.or.baf.lists <- cleanAndErrorCheckSNPs(all.samp.SNP.list, min.cons.segs.df, allow.logr.only.flag = TRUE)

	asas.SNPs.w.both.lr.baf.snp.lists <- cleanAndErrorCheckSNPsMPhase(all.samp.SNP.list, min.cons.segs.df, allow.logr.only.flag = FALSE)
	asas.SNPS.lr.and.or.baf.lists <- cleanAndErrorCheckSNPsMPhase(all.samp.SNP.list, min.cons.segs.df, allow.logr.only.flag = TRUE)








	if ( length(asas.SNPs.w.both.lr.baf.snp.lists) <= 1 ) {

		print(str(asas.SNPs.w.both.lr.baf.snp.lists, max = 1))
		stop("Only one or fewer samples present after reading in data.")

	}

	#save:
	saveRDS(asas.SNPS.lr.and.or.baf.lists, file = asas.SNPS.lr.and.or.baf.lists.path)
	saveRDS(asas.SNPs.w.both.lr.baf.snp.lists, file = asas.SNPs.w.both.lr.baf.snp.lists.path)
	saveRDS(min.cons.segs.df, file = min.cons.segs.df.path)
	saveRDS(processed.orig.cpn.segs.df, file = processed.orig.cpn.segs.df.path) 
	saveRDS(ploidy.est.vec, file = ploidy.est.vec.path)
	saveRDS(purity.est.vec, file = purity.est.vec.path)

} else if (OVERWRITE_FLAG == FALSE) {#end if (OVERWRITE_FLAG == TRUE) {

	# if ( ! file.exists(processed.orig.cpn.segs.df.path) ) {

	# 	seg.file.df <- read.table(seg.file.path, stringsAsFactors = FALSE, sep = "\t", header = TRUE)
	# 	seg.file.df$sample <- gsub("\\.", "-", gsub("_LogR", "", seg.file.df$sample))
	# 	seg.file.df$callState <- NULL
	# 	seg.file.df$patient <- patient
	# 	seg.file.df$sample_name <- paste(patient, sapply(seg.file.df$sample, function(x) strsplit(x, split = "-")[[1]][2]), sep = ":")
	# 	processed.orig.cpn.segs.df <- errorCheckAndAddPloidyAdjRawCPNToSegDF(patient, seg.file.df)
	# 	saveRDS(processed.orig.cpn.segs.df, file = processed.orig.cpn.segs.df.path)

	# }

	#load original data:
	asas.SNPS.lr.and.or.baf.lists <- readRDS(file = asas.SNPS.lr.and.or.baf.lists.path)
	asas.SNPs.w.both.lr.baf.snp.lists <- readRDS(file = asas.SNPs.w.both.lr.baf.snp.lists.path)
	min.cons.segs.df <- readRDS(file = min.cons.segs.df.path)
	processed.orig.cpn.segs.df <- readRDS(file = processed.orig.cpn.segs.df.path)
	ploidy.est.vec <- readRDS(file = ploidy.est.vec.path)
	purity.est.vec <- readRDS(file = purity.est.vec.path)

}#end else if (OVERWRITE_FLAG == FALSE


#Now what about intersection of all samples what proportion of the genome to they cover:
samp.gr.list <- list()
overall.intersection.gr <- NULL
for ( samp in sample.vec ) {

	samp.gr <- makeGRangesFromDataFrame(f.joint.segs.df[f.joint.segs.df$sample==samp 
									  ,],
	                        keep.extra.columns=FALSE,
	                        ignore.strand=FALSE,
	                        seqinfo=NULL,
	                        seqnames.field=c("seqnames", "seqname",
	                                         "chromosome", "chrom",
	                                         "chr", "chromosome_name",
	                                         "seqid"),
	                        start.field=c("Start_Position","start","startpos"),
	                        end.field=c("End_Position", "end", "stop", "endpos"),
	                        strand.field="strand",
	                        starts.in.df.are.0based=FALSE)

	if ( is.null(overall.intersection.gr) ) {

		overall.intersection.gr <- samp.gr

	} else {

		overall.intersection.gr <- intersect(overall.intersection.gr, samp.gr)

	}#end ! is.null (overall.intersection.gr)

}#en for ( samp in sample.vec ) {

colnames(min.cons.segs.df)[1:3] <- c("chr", "start", "end")

min.cons.segs.gr <- makeGRangesFromDataFrame(min.cons.segs.df)

setdiff(min.cons.segs.gr,overall.intersection.gr)

overall.intersection.df <- as.data.frame(overall.intersection.gr)

overall.intersection.df$chr_pos <- paste(overall.intersection.df$seqnames, ":", overall.intersection.df$start, ":", overall.intersection.df$end, sep="")
chr.pos.vec <- paste0(min.cons.segs.df[,"chr"], ":", min.cons.segs.df[,"start"], ":", min.cons.segs.df[,"end"])
#if we're making a 
# Now we need to make sure that the asas.lr.baf.lists and sample.to.phase.seg.list 
# if (length(names(sample.to.phase.seg.list)) != length(asas.lr.baf.lists[[1]])) {

#   print("Missing or removed segments:")
#   print(names(sample.to.phase.seg.list)[ ! names(asas.lr.baf.lists[[1]]) %in% names(sample.to.phase.seg.list)])
#   print("Continuing creation of phased BAF dataframe")

# }

#Todo, Add function definitions documentaiton
#Todo, clarify the loops and next statements that are controlling what happens when there are less than 5 hz snps in a sample.
#Todo, check if the mixture phasing is ever being selected... in the example with B_MET002 no mixture phasing called.
seg.name.vec <- paste(min.cons.segs.df[, 1], min.cons.segs.df[, 2], min.cons.segs.df[, 3], sep = ":")
sample.vec <- names(asas.SNPs.w.both.lr.baf.snp.lists)

if (OVERWRITE_FLAG == TRUE) {

	print("Getting phasing to apply to all samples...")
	all.seg.phasing.info.list <- getCombinedBAFDipTestPhasing(asas.lr.baf.lists = asas.SNPs.w.both.lr.baf.snp.lists,
	                                            seg.name.vec = seg.name.vec,
	                                            ploidy.est.vec = ploidy.est.vec,
	                                            purity.est.vec = purity.est.vec,
	                                            gamma.val = GAMMA_VAL,
	                                            min.num.hz.snp.phasing = MIN_NUM_HZ_SNP_PHASING,
	                                            min.prop.snp.sep.threshold = MIN_PROP_SNP_SEP_THRESH,
	                                            dip.test.p.val.threshold = PVAL_DIPTEST_AI_THRESH,
	                                            comb.wilcox.snp.sep.p.val.threshold = PVAL_COMB_WILCOX_PH_THRESH,
	                                            allow.dip.test.without.ASCAT.to.phase = ALLOW_DIPTEST_ALONE_PH_FLAG
	                                            )
	
	saveRDS(all.seg.phasing.info.list, file = all.seg.phasing.info.list.path)


	#Todo unsure of what to do with the baf logr lists which have positions that only have logr here, at the moment exclude them but how
	#does that affect later run
	print("Generating phased BAF dataframes...")
	patient.phased.BAF.df.list <- createAllSegsPhasedBAFDataFrame(asas.lr.baf.lists = asas.SNPs.w.both.lr.baf.snp.lists,
                                                         sample.to.phase.seg.list = all.seg.phasing.info.list,
                                                         min.num.hz.snp.phasing = MIN_NUM_HZ_SNP_PHASING
                                                         )
	saveRDS(patient.phased.BAF.df.list, file = patient.phased.BAF.df.list.path)


	#Detect MSAI in the original BAF before it is processed alongside logr to get the Raw allele specifc cpn to be used in the rest of the analysis.
	print("Performing LTNPD MSAI detection on the BAF...")
	all.seg.MSAI.by.BAF.ltnpd.present.list <- getAllSegsMSAIDetectionFromBAFListLTNPD(seg.name.vec = seg.name.vec,
																				  sample.vec  = sample.vec,
																				  asas.SNPS.lr.and.or.baf.lists = asas.SNPS.lr.and.or.baf.lists,
																				  all.seg.phasing.info.list = all.seg.phasing.info.list,
																				  MSAI.BAF.LTPND.pval.thresh = MSAI_BAF_LTNPD_PVAL_THRESHOLD)
	saveRDS(all.seg.MSAI.by.BAF.ltnpd.present.list, file = all.seg.MSAI.by.BAF.ltnpd.present.list.path)

	#Todo add another method of detecting MSAI using only the BAF and a paired T test.
	# all.seg.MSAI.by.BAF.ltnpd.present.list <- getAllSegsMirroredSubclonalAllelicImbalanceDetectionFromBAFList(seg.name.vec = seg.name.vec,
	# 																										   sample.vec  = sample.vec,
	# 																										   asas.SNPS.lr.and.or.baf.lists = asas.SNPS.lr.and.or.baf.lists,
	# 																										   all.seg.phasing.info.list = all.seg.phasing.info.list,
	# 																										   MSAI.BAF.LTPND.pval.thresh = MSAI_BAF_LTNPD_PVAL_THRESHOLD)

	print("Generating datapoints with LogR total cpn dataframes...")
	#This is where the SNP filtering occurs. We can fix this by adding by passing an option of which data type is being analysed...
	all.samp.non.hz.logr.and.tot.cpn.df.list <- getAllSampLogrAlonePositionsLogrAndTotalCPNDFsList(sample.vec = sample.vec,
																				asas.SNPS.lr.and.or.baf.lists = asas.SNPS.lr.and.or.baf.lists,
																				ploidy.est.vec = ploidy.est.vec,
																				purity.est.vec = purity.est.vec,
																				gamma.val = GAMMA_VAL,
																				data.format.to.use = data.format.to.use)
	saveRDS(all.samp.non.hz.logr.and.tot.cpn.df.list, file = all.samp.non.hz.logr.and.tot.cpn.df.list.path)





	#Add this the unphaseable stuff too maybe? Not necessary.
	print("Applying phasing to all samples to get phased raw cpn...")
	asas.ph.raw.cpn.lists <- getAllSegAllSamplePhasedRawCPNList(sample.vec = sample.vec,
																ploidy.est.vec = ploidy.est.vec,
																purity.est.vec = purity.est.vec,
																min.cons.segs.df = min.cons.segs.df,
																all.seg.phasing.info.list = all.seg.phasing.info.list,
																asas.both.lr.baf.lists = asas.SNPs.w.both.lr.baf.snp.lists,
																gamma.val = GAMMA_VAL,
																min.num.hz.snp.phasing = MIN_NUM_HZ_SNP_PHASING)
	saveRDS(asas.ph.raw.cpn.lists, file = asas.ph.raw.cpn.lists.path)



	#Perform all pairwise comparisons of all samples, the use of much of this has been deprecated but we can revisit it,
	#We can add to this later if we want to redo the
	# print("Performing pairwise comparisons of phased raw cpn...")
	# all.pw.cn.het.comp.lists <- list()
	# for (i in 1:(length(sample.vec) - 1)) {

	#   for(j in (i + 1):length(sample.vec)) {

	#     samp.one <- sample.vec[i]
	#     samp.two <- sample.vec[j]
	#     samp.one.two <- paste(samp.one, ":", samp.two, sep="")
	#     s.one.two.comp.input.list <- list()
	#     s.one.two.comparison.list <- list()
	#     s.one.two.pw.het.list <- list()
	#     print(samp.one.two)
	#     for ( seg.name in seg.name.vec ) {
	#     	print(seg.name)
	#     	s.one.two.comp.input.list[[seg.name]] <- createAllSegInputForHetComparisons(seg.name = seg.name,
	#                                                   samp.one = samp.one,
	#                                                   samp.two = samp.two,
	#                                                   asas.ph.raw.cpn.lists = asas.ph.raw.cpn.lists,
	#                                                   asas.both.lr.baf.lists = asas.both.lr.baf.lists,
	#                                                   all.samp.non.hz.logr.and.tot.cpn.df.list = all.samp.non.hz.logr.and.tot.cpn.df.list,
	#                                                   min.num.hz.snp.phasing = MIN_NUM_HZ_SNP_PHASING)

	#     	#logr alone comparison: possible for all segments
	#     	#raw ph cpn comparison: possible only for segments thad had greater than MIN_NUM_HZ_SNP_PHASING in them
	#     	#raw tot cpn: possible for all segments
	#     	#processed call comparisons: currently possible only for segments with phased copy number
	#     	#Also there separate and combined logr and total cpn calls we can make that might come up as being different,
	#     	#We need to separate these out early on and then run the same functions twice over them.
	# 	  	s.one.two.comparison.list[[seg.name]] <- performSegPairwiseComparisons(seg.name = seg.name,
	# 	  	                                          					samp.one = samp.one,
	# 	  	                                          					samp.two = samp.two,
	# 	  	                                          					seg.samp.one.two.comp.data.list = s.one.two.comp.input.list[[seg.name]],
	# 	  	                                          					purity.est.vec = purity.est.vec,
	# 	  	                                          					ploidy.est.vec = ploidy.est.vec,
	# 	  	                                          					gain.threshold = GAIN_LOGR_THRESH,
	# 	  	                                          					amp.threshold = LOSS_LOGR_THRESH,
	# 	  	                                          					loss.threshold = AMP_LOGR_THRESH,
	# 	  	                                          					deep.loss.threshold = DEEP_LOSS_LOGR_THRESH)

	# 	  	s.one.two.pw.het.list[[seg.name]] <- calculateSegPairwiseCNHeterogeneityMeasures(seg.samp.comp.het.analysis.list = s.one.two.comparison.list[[seg.name]],
	# 	  																					 pval.for.diff.amplitude.thresh = PW_CPNA_CPNB_DIFF_AMPLITUDE_THRESH) 
	# 	}#end for ( seg.name in seg.name.vec ) {

	# 		all.pw.cn.het.comp.lists[[samp.one.two]] <- list()
	# 		all.pw.cn.het.comp.lists[[samp.one.two]]$comparison_input_list <- s.one.two.comp.input.list
	# 		all.pw.cn.het.comp.lists[[samp.one.two]]$comparison_output_list <- s.one.two.comparison.list
	# 		all.pw.cn.het.comp.lists[[samp.one.two]]$s_one_two_pw_het_list <- s.one.two.pw.het.list


	#   }#end for (i in 1:(length(sample.vec) - 1)) {

	# }#end for(j in (i + 1):length(sample.vec)) {
	# all.pw.cn.het.comp.lists.path <- paste0(PLOT_OUTPUT_DIR, patient, "_all_pw_cn_het_comp_lists.RDS")
	# saveRDS(all.pw.cn.het.comp.lists, file = all.pw.cn.het.comp.lists.path)



	#Let's get the CCF
	print("Checking all minimum consistent segments are present in each sample...")
	seg.name.vec <- names(asas.ph.raw.cpn.lists[[1]])
	#Check to make sure that all segments are present in all samples... might not be necessary todo, check.
	for ( samp in sample.vec ) {

	      seg.name.vec <- intersect(seg.name.vec, names(asas.ph.raw.cpn.lists[[samp]]))

	}#end for ( samp in sample.vec ) { 

	if (length(seg.name.vec) == 0) {

		stop("No shared segments for some reason, go back and check")

	}


	print("Calculating CCFs using phased raw cpn....")
	all.seg.CCF.list <- list()
	all.seg.int.ob.xp.cpn.df.list <- list()
	all.seg.CCF.diff.from.exp.list <- list()
	all.seg.clonal.in.all.samps.list <- list()
	all.seg.greater.than.two.diff.list <- list()

	samps.int.cpn.raw.cpn.by.chr.list <- generateAllSamplesByChromIntCPNUnderlyingRawValuesLists(asas.ph.raw.cpn.lists = asas.ph.raw.cpn.lists,
																								 sample.vec = sample.vec,
																								 seg.name.vec = seg.name.vec,
																								 cpnX.vec=c("cpnA", "cpnB"),
																								 reduced.seg.flag = FALSE)
	saveRDS(samps.int.cpn.raw.cpn.by.chr.list, file = samps.int.cpn.raw.cpn.by.chr.list.path)


	for (cpnX in c("cpnA", "cpnB")) {
	  print(cpnX)
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

	    if (is.na(asas.ph.raw.cpn.lists[[1]][[seg.name]]$ph_cpnA_vec[1])) {

	      all.seg.CCF.list[[cpnX]][[seg.name]] <- NA
	      all.seg.int.ob.xp.cpn.df.list[[cpnX]][[seg.name]] <- NA
	      all.seg.CCF.diff.from.exp.list[[cpnX]][[seg.name]] <- NA
	      all.seg.clonal.in.all.samps.list[[cpnX]][[seg.name]] <- NA
	      all.seg.greater.than.two.diff.list[[cpnX]][[seg.name]] <- NA
	      next 

	    }

	    #does both A and B, todo, refactor so only does one
	    s.samps.int.ob.xp.cpn.df <- getSegObsExpCPNBySampleDF(asas.ph.raw.cpn.lists = asas.ph.raw.cpn.lists,
	                                                          sample.vec = sample.vec,
	                                                          seg.name = seg.name,
	                                                          cpnX.vec = cpnX,
	                                                          reduced.seg.flag = FALSE) 

	    #does only cpnX
	    seg.cpnX.diff.from.exp.list <- calculateWhetherSegCPNDiffFromIntegerExpectedAndAllMatchingSampCPN(asas.ph.raw.cpn.lists = asas.ph.raw.cpn.lists,
	                                                                                                      s.samps.int.ob.xp.cpn.df = s.samps.int.ob.xp.cpn.df,
	                                                                                                      seg.name = seg.name,
	                                                                                                      samps.int.cpn.raw.cpn.by.chr.list = samps.int.cpn.raw.cpn.by.chr.list,
	                                                                                                      cpnX = cpnX,
	                                                                                                      chrs.can.use.for.ccf.test = CHRS_TO_USE_VEC,
	                                                                                                      reduced.seg.flag = FALSE)
	    #does only cpnX
	    is.segment.clonal.in.all.samps.flag <- determineCCFClonalityOfSegment(s.samps.int.ob.xp.cpn.df = s.samps.int.ob.xp.cpn.df,
	                                                                          seg.name = seg.name,
	                                                                          asas.ph.raw.cpn.lists = asas.ph.raw.cpn.lists,
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
	      #ss.ph.cpnX.vals.vec <- all.pw.cn.het.comp.lists[[samp]][[seg.name]][[paste0("ph_", cpnX, "_vec")]]
	      ss.ph.cpnX.vals.vec <- asas.ph.raw.cpn.lists[[samp]][[seg.name]][[paste0("ph_", cpnX, "_vec")]]

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
	      all.seg.CCF.list[[cpnX]][[seg.name]][[samp]][[paste0(cpnX, "_orig_ph_obs_cpnA_vec")]]   <- asas.ph.raw.cpn.lists[[samp]][[seg.name]][[paste0("ph_cpnA_vec")]]
	      all.seg.CCF.list[[cpnX]][[seg.name]][[samp]][[paste0(cpnX, "_orig_ph_obs_cpnB_vec")]]   <- asas.ph.raw.cpn.lists[[samp]][[seg.name]][[paste0("ph_cpnB_vec")]]
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

	#Save CCF data:
	saveRDS(all.seg.CCF.list, file = all.seg.CCF.list.path)
	saveRDS(all.seg.int.ob.xp.cpn.df.list, file = all.seg.int.ob.xp.cpn.df.list.path)
	saveRDS(all.seg.CCF.diff.from.exp.list, file = all.seg.CCF.diff.from.exp.list.path)
	saveRDS(all.seg.clonal.in.all.samps.list, file = all.seg.clonal.in.all.samps.list.path)
	saveRDS(all.seg.greater.than.two.diff.list, file = all.seg.greater.than.two.diff.list.path)


} else if (OVERWRITE_FLAG == FALSE) {

	# print("Getting phasing to apply to all samples...")
	# all.seg.phasing.info.list <- getCombinedBAFDipTestPhasing(asas.lr.baf.lists = asas.SNPs.w.both.lr.baf.snp.lists,
	#                                             seg.name.vec = seg.name.vec,
	#                                             ploidy.est.vec = ploidy.est.vec,
	#                                             purity.est.vec = purity.est.vec,
	#                                             gamma.val = GAMMA_VAL,
	#                                             min.num.hz.snp.phasing = MIN_NUM_HZ_SNP_PHASING,
	#                                             min.prop.snp.sep.threshold = MIN_PROP_SNP_SEP_THRESH,
	#                                             dip.test.p.val.threshold = PVAL_DIPTEST_AI_THRESH,
	#                                             comb.wilcox.snp.sep.p.val.threshold = PVAL_COMB_WILCOX_PH_THRESH,
	#                                             allow.dip.test.without.ASCAT.to.phase = ALLOW_DIPTEST_ALONE_PH_FLAG
	#                                             )
	# saveRDS(all.seg.phasing.info.list, file = all.seg.phasing.info.list.path)

	all.seg.phasing.info.list <- readRDS(file = all.seg.phasing.info.list.path)
	patient.phased.BAF.df.list <- readRDS(file = patient.phased.BAF.df.list.path)
	all.seg.MSAI.by.BAF.ltnpd.present.list <- readRDS(file = all.seg.MSAI.by.BAF.ltnpd.present.list.path)
	all.samp.non.hz.logr.and.tot.cpn.df.list <- readRDS(file = all.samp.non.hz.logr.and.tot.cpn.df.list.path)
	asas.ph.raw.cpn.lists <- readRDS(file = asas.ph.raw.cpn.lists.path)
	samps.int.cpn.raw.cpn.by.chr.list <- readRDS(file = samps.int.cpn.raw.cpn.by.chr.list.path)
	all.seg.CCF.list <- readRDS(file = all.seg.CCF.list.path)
	all.seg.int.ob.xp.cpn.df.list <- readRDS(file = all.seg.int.ob.xp.cpn.df.list.path)
	all.seg.CCF.diff.from.exp.list <- readRDS(file = all.seg.CCF.diff.from.exp.list.path)
	all.seg.clonal.in.all.samps.list <- readRDS(file = all.seg.clonal.in.all.samps.list.path)
	all.seg.greater.than.two.diff.list <- readRDS(file = all.seg.greater.than.two.diff.list.path)


}
#Todo remove the dataframe initialisation from the function and modify it to work on a segment by segment basis.

#Let's initialise the first joint output:
print("Initialising minimum consistent segmentation (MCS) output dataframe...")
empty.full.joint.segs.df <- initialiseFullOutputDataframe(seg.name.vec = names(asas.ph.raw.cpn.lists[[1]]),
														ploidy.est.vec = ploidy.est.vec,
														purity.est.vec = purity.est.vec,
														make.row.names.seg.name.alone = FALSE) 


#This adds the raw cpn for both segments that could be phased and those that couldn't
print("Adding new copy number calls to MCS output based on the phased raw copy number...")
new.cpn.full.joint.segs.df <- addNewRawAndIntCPNCallsToFinalSegmentation(seg.name.vec = names(asas.ph.raw.cpn.lists[[1]]),
																		sample.vec = names(purity.est.vec),
																		all.samp.logr.df.list = all.samp.non.hz.logr.and.tot.cpn.df.list,
																		all.samp.ph.BAF.df.list = patient.phased.BAF.df.list,
																		all.seg.CCF.vals.list = all.seg.CCF.list,
																		asas.lr.baf.lists = asas.SNPs.w.both.lr.baf.snp.lists,
																		all.seg.df = empty.full.joint.segs.df,
																		cpnX.vec = c("cpnA", "cpnB")
																		)


print("Mapping MCS segments back to original input and copy number states...")
full.joint.segs.w.old.cpn.vs.ploidy.df <- addOriginalSegmentInfoToFinalSegmentation(seg.name.vec = names(asas.ph.raw.cpn.lists[[1]]),
																					sample.vec = names(purity.est.vec),
																					processed.orig.cpn.segs.df = processed.orig.cpn.segs.df,
																					gain.threshold = GAIN_LOGR_THRESH,
																					loss.threshold = LOSS_LOGR_THRESH,
																					amp.threshold = AMP_LOGR_THRESH,  
																					deep.loss.threshold = DEEP_LOSS_LOGR_THRESH,
																					all.seg.df = new.cpn.full.joint.segs.df
																					)


print("Adding CCF information and subclonal cpn calls to MCS output...")
clonality.and.ccf.joint.segs.df <- addWithinSampleClonalityAndCCFToFinalSegmentation(seg.name.vec = names(asas.ph.raw.cpn.lists[[1]])
                                ,all.seg.CCF.vals.list = all.seg.CCF.list
                                ,all.seg.CCF.diff.from.exp.list = all.seg.CCF.diff.from.exp.list
                                ,all.seg.df = full.joint.segs.w.old.cpn.vs.ploidy.df
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

print("Adding total copy number estimates to MCS output for segments with only logr information...")
full.joint.segs.df.w.cpn.vs.ploidy.df <- addLogrDefinedCPNEventsToOutputSegmentation(sample.vec = names(purity.est.vec),
                                                        							seg.name.vec = names(asas.ph.raw.cpn.lists[[1]]),
                                                        							all.seg.df = clonality.and.ccf.joint.segs.df,
                                                        							all.samp.logr.df.list = all.samp.non.hz.logr.and.tot.cpn.df.list,
                                                        							all.samp.ph.BAF.df.list = patient.phased.BAF.df.list,     
                                                       								gain.threshold = GAIN_LOGR_THRESH,
                                                        							loss.threshold = LOSS_LOGR_THRESH,
                                                        							amp.threshold = AMP_LOGR_THRESH,  
                                                        							deep.loss.threshold = DEEP_LOSS_LOGR_THRESH,
                                                        							tot.cpn.vs.ploidy.pval.thresh = PVAL_THRESH_VS_LOGR_FOR_CPN_EVENT)

print(head(full.joint.segs.df.w.cpn.vs.ploidy.df[,"cpn_event_vs_ploidy"]))
print("Adding cpn event vs ploidy by logr measures to MCS output...")
full.joint.segs.logr.vs.ploidy.clonality.df <- determineAcrossSampleClonalityOfLogrEventsVsPloidy(sample.vec = names(purity.est.vec),
                                                        							seg.name.vec = names(asas.ph.raw.cpn.lists[[1]]),
                                                        							all.seg.df = full.joint.segs.df.w.cpn.vs.ploidy.df,
                                                        							all.samp.logr.df.list = all.samp.non.hz.logr.and.tot.cpn.df.list,
                                                        							all.samp.ph.BAF.df.list = patient.phased.BAF.df.list,     
                                                       								gain.threshold = GAIN_LOGR_THRESH,
                                                        							loss.threshold = LOSS_LOGR_THRESH,
                                                        							amp.threshold = AMP_LOGR_THRESH,  
                                                        							deep.loss.threshold = DEEP_LOSS_LOGR_THRESH,
                                                        							tot.cpn.vs.ploidy.pval.thresh = PVAL_THRESH_VS_LOGR_FOR_CPN_EVENT)


print("Testing whether any gains should be reclassified as amplifications...")
full.joint.check.amps.vs.ploidy.clonality.df <- retestSubclonalAmps(all.seg.df = full.joint.segs.logr.vs.ploidy.clonality.df,
																	all.samp.logr.df.list = all.samp.non.hz.logr.and.tot.cpn.df.list,
																	all.samp.ph.BAF.df.list = patient.phased.BAF.df.list,   
																	gain.threshold = GAIN_LOGR_THRESH,
																	amp.threshold = AMP_LOGR_THRESH,
																	p.val.for.converting.gain.to.amp = PVAL_THRESH_CONVERT_GAIN_TO_AMP)

print(head(full.joint.segs.logr.vs.ploidy.clonality.df[,c("across_samp_cl_cpn_vs_ploidy_present","across_samp_sc_cpn_vs_ploidy_present")]))


#These are kept in separate functions as the original segmentaiton was refactored so it could run when segments are rejoined at the end of the script.
print("Adding phasing information and MSAI detection by raw cpn calls to MCS output...")
full.joint.segs.w.phasing.df <- addPhasingInformationToOutputSegmentation(sample.vec = names(purity.est.vec),
                                        seg.name.vec = names(asas.ph.raw.cpn.lists[[1]]), 
                                        sample.to.phase.seg.list = all.seg.phasing.info.list,
                                        all.pw.cn.het.comp.lists = NULL,
                                        all.segs.df = full.joint.segs.logr.vs.ploidy.clonality.df)


#Now we want to add in the allelic imbalance clonal and subclonal inrormation to the final data frame.
print("Adding across sample clonal/subclonal summary to MCS output...")
full.joint.segs.w.phasing.and.raw.df <- addOriginalCPNValuesAndCallsVsPloidyToOutputSegmentation(sample.vec = sample.vec,
                                        seg.name.vec = names(asas.ph.raw.cpn.lists[[1]]), 
                                        processed.orig.cpn.segs.df = processed.orig.cpn.segs.df,
                                        all.segs.df = full.joint.segs.w.phasing.df)

print("Adding MSAI detected by LTNPD method to the MCS output...")
full.joint.segs.w.ltnpd.MSAI.df <- addMSAIDetectionByLTNPDToOutputSegmentation(sample.vec = sample.vec,
											seg.name.vec = names(asas.ph.raw.cpn.lists[[1]]),
											all.seg.MSAI.by.BAF.ltnpd.present.list = all.seg.MSAI.by.BAF.ltnpd.present.list,
											all.segs.df = full.joint.segs.w.phasing.and.raw.df)


print(head(full.joint.segs.w.ltnpd.MSAI.df[( ! is.na(full.joint.segs.w.ltnpd.MSAI.df$any_across_samp_MSAI)) &
											full.joint.segs.w.ltnpd.MSAI.df$any_across_samp_MSAI == TRUE,
											c("MSAI", "any_across_samp_MSAI")]))


print("Defining any parallel gain or loss events detected via MSAI and addingto MCS output...")
full.joint.segs.w.parallel.cpn.events.df <- addParallelEventsByMSAIToOutputSegmentation(all.segs.df = full.joint.segs.w.ltnpd.MSAI.df)


print(head(full.joint.segs.w.parallel.cpn.events.df[(! is.na(full.joint.segs.w.parallel.cpn.events.df$MSAI)),
													c("any_across_samp_MSAI", "MSAI",
														"cpn_event_vs_ploidy",
														"is_parallel_gain_and_or_amp",
														"is_parallel_loss",
														"is_parallel_LOH",
														"is_parallel_loss_and_or_LOH")]))

print("Adding comparison between original cpn calls and MPhase copy number calling to MCS output...")
full.joint.segs.within.and.changed.cpn.info.df <- addWithinSampleChangeToCpnAndAllelicImbalanceCalls(full.joint.segs.w.parallel.cpn.events.df)


print("Adding across sample estimates of allelic imbalance and LOH clonality to MCS output...")
full.joint.segs.w.across.samp.AI.stats.df <- addAcrossSampleAlleleSpecificClonalityCalls(all.seg.df = full.joint.segs.within.and.changed.cpn.info.df,
																						sample.vec = sample.vec,
																						seg.name.vec = seg.name.vec)

print(head(full.joint.segs.w.across.samp.AI.stats.df[!is.na(full.joint.segs.w.across.samp.AI.stats.df$within_samp_sc_cpnA_present), c("any_within_samp_sc_cpnA",
"any_within_samp_sc_cpnB",
"across_samp_cl_AI",
"across_samp_sc_AI",
"across_samp_cl_LOH",
"across_samp_sc_LOH")]))


full.joint.segs.w.across.plotting.cols.df <- addPanCanColumnDefinitionToSegDF(full.joint.segs.w.across.samp.AI.stats.df)

print(full.joint.segs.w.across.plotting.cols.df[, c("pt_seg_cl_loss_LOH_no_MSAI", "pt_seg_sc_loss_LOH_incl_MSAI", "pt_seg_cl_gain_amp_no_MSAI", "pt_seg_sc_gain_amp_incl_MSAI")])

full.joint.segs.df <- full.joint.segs.w.across.plotting.cols.df
#Remove segments without any datapoints.
#Todo, merge these int
#print("Removing any segments that do not have underlying datapoints, this is an artefact of the minimum consistent segmentation performed by MPhase...")
#full.joint.segs.df <- full.joint.segs.df[full.joint.segs.df$num_un_ph_hz_sites != 0 | full.joint.segs.df$num_ph_hz_sites != 0  | full.joint.segs.df$num_non_hz_sites != 0,]
#print("Removing segments that do not map back to any original segments")
#full.joint.segs.df <- full.joint.segs.df[full.joint.segs.df$orig_seg_name != "NA:NA:NA",]






full.joint.segs.df.path <- paste0(PLOT_OUTPUT_DIR, patient, "_full_joint_segs_df.RDS")
saveRDS(full.joint.segs.df, file=full.joint.segs.df.path)


####Plot all genome wide plots
print("Plotting CCFs for cpnA and cpnB across the genome using CCF definition of subclonality...")
png(file = paste(PLOT_OUTPUT_DIR, patient, "_CCF_w_CCF_thresh_cpnA", ".png", sep = ""), width = 1200, height = 1400)

plotAllSampleCCFsAcrossGenome(asas.ph.raw.cpn.lists = asas.ph.raw.cpn.lists,
                              full.joint.segs.df = full.joint.segs.w.phasing.df,
                              sample.vec = sample.vec,
                              all.seg.CCF.vals.list = all.seg.CCF.list,
                              cpnX = "cpnA",
                              use.cpn.or.ccf.thresh = "CCF")

dev.off()

png(file = paste(PLOT_OUTPUT_DIR, patient, "_CCF_w_CCF_thresh_cpnB", ".png", sep = ""), width = 1200, height = 1400)

plotAllSampleCCFsAcrossGenome(asas.ph.raw.cpn.lists = asas.ph.raw.cpn.lists,
                              full.joint.segs.df = full.joint.segs.w.phasing.df,
                              sample.vec = sample.vec,
                              all.seg.CCF.vals.list = all.seg.CCF.list,
                              cpnX = "cpnB",
                              use.cpn.or.ccf.thresh = "CCF")
dev.off()

####Plot all genome wide plots
print("Plotting CCFs for cpnA and cpnB across the genome using cpn definition of subclonality...")
png(file = paste(PLOT_OUTPUT_DIR, patient, "_CCF_w_cpn_thresh_cpnA", ".png", sep = ""), width = 1200, height = 1400)

plotAllSampleCCFsAcrossGenome(asas.ph.raw.cpn.lists = asas.ph.raw.cpn.lists,
                              full.joint.segs.df = full.joint.segs.w.phasing.df,
                              sample.vec = sample.vec,
                              all.seg.CCF.vals.list = all.seg.CCF.list,
                              cpnX = "cpnA",
                              use.cpn.or.ccf.thresh = "CPN")

dev.off()

png(file = paste(PLOT_OUTPUT_DIR, patient, "_CCF_w_cpn_thresh_cpnB", ".png", sep = ""), width = 1200, height = 1400)

plotAllSampleCCFsAcrossGenome(asas.ph.raw.cpn.lists = asas.ph.raw.cpn.lists,
                              full.joint.segs.df = full.joint.segs.w.phasing.df,
                              sample.vec = sample.vec,
                              all.seg.CCF.vals.list = all.seg.CCF.list,
                              cpnX = "cpnB",
                              use.cpn.or.ccf.thresh = "CPN")
dev.off()

png(file = paste(PLOT_OUTPUT_DIR, patient, "_BAF_phased", ".png", sep=""), width = 1200, height = 1400)

print("Plotting phased BAF across the genome...")
plotAllSampleBAFsAcrossGenome(sample.vec, patient.phased.BAF.df.list, cpnA.allele.color = '#de2d2699', cpnB.allele.color= '#3182bd99')

dev.off()

#Plot the region comparison across genome plots, will need their own directories:
old.vs.new.raw.cpn.dir <- paste(PLOT_OUTPUT_DIR, "old_raw_cpn_vs_MPhase_raw_cpn/", sep="")
if( ! file.exists(old.vs.new.raw.cpn.dir)) {

	dir.create(old.vs.new.raw.cpn.dir)

}#end if( ! file.exists(by.chrom.CCF.plot.dir)) {

print("Plotting comparison of old vs new MPhase raw cpn calls..")
plotOldNewComparisonOfAllCPN(sample.vec, full.joint.segs.df, plot.output.dir = old.vs.new.raw.cpn.dir)

pairwise.raw.cpn.dir <- paste(PLOT_OUTPUT_DIR, "pairwise_raw_cpn_comparisons/", sep = "")
if( ! file.exists(pairwise.raw.cpn.dir)) {

	dir.create(pairwise.raw.cpn.dir)

}#end if( ! file.exists(by.chrom.CCF.plot.dir)) {

print("Plotting pairwise sample comparison between MPhase raw cpn calls across the whole genome...")
plotAllSamplePairwiseComparisonSegmentCPN(sample.vec, full.joint.segs.df, plot.output.dir = pairwise.raw.cpn.dir)



#Plot by chromosome plots:

uniq.chrs <- unique(full.joint.segs.w.phasing.df$chr)
by.chrom.CCF.plot.dir <- paste(PLOT_OUTPUT_DIR, "CCF_by_chrom/", sep = "")
if( ! file.exists(by.chrom.CCF.plot.dir)) {

	dir.create(by.chrom.CCF.plot.dir)

}#end if( ! file.exists(by.chrom.CCF.plot.dir)) {
print("Plotting by chromosome CCF estimates for cpnA and cpnB using CCF subclonality thresholds...")
png(file = paste(by.chrom.CCF.plot.dir, "/", patient, "_CCF_w_CCF_thresh_cpnA_", "%02d", ".png", sep = ""), width = 1200, height = 1400)
for ( chrom in uniq.chrs ) {

	plotAllSampleCCFsAcrossGenome(asas.ph.raw.cpn.lists = asas.ph.raw.cpn.lists,
                              full.joint.segs.df = full.joint.segs.w.phasing.df[full.joint.segs.w.phasing.df$chr == chrom,] ,
                              sample.vec = sample.vec,
                              all.seg.CCF.vals.list = all.seg.CCF.list,
                              cpnX = "cpnA",
                              use.cpn.or.ccf.thresh = "CCF")


}
dev.off()

png(file = paste(by.chrom.CCF.plot.dir, "/", patient, "_CCF_w_CCF_thresh_cpnB_", "%02d", ".png", sep = ""), width = 1200, height = 1400)
for ( chrom in uniq.chrs ) {

	plotAllSampleCCFsAcrossGenome(asas.ph.raw.cpn.lists = asas.ph.raw.cpn.lists,
                              full.joint.segs.df = full.joint.segs.w.phasing.df[full.joint.segs.w.phasing.df$chr == chrom,] ,
                              sample.vec = sample.vec,
                              all.seg.CCF.vals.list = all.seg.CCF.list,
                              cpnX = "cpnB",
                              use.cpn.or.ccf.thresh = "CCF")


}
dev.off()

print("Plotting by chromosome CCF estimates for cpnA and cpnB using cpn subclonality thresholds...")
png(file = paste(by.chrom.CCF.plot.dir, "/", patient, "_CCF_w_cpn_thresh_cpnA_","%02d", ".png", sep = ""), width = 1200, height = 1400)
for ( chrom in uniq.chrs ) {

	plotAllSampleCCFsAcrossGenome(asas.ph.raw.cpn.lists = asas.ph.raw.cpn.lists,
                              full.joint.segs.df = full.joint.segs.w.phasing.df[full.joint.segs.w.phasing.df$chr == chrom,] ,
                              sample.vec = sample.vec,
                              all.seg.CCF.vals.list = all.seg.CCF.list,
                              cpnX = "cpnA",
                              use.cpn.or.ccf.thresh = "CPN")


}
dev.off()



png(file = paste(by.chrom.CCF.plot.dir, "/", patient, "_CCF_w_cpn_thresh_cpnB_", "%02d", ".png", sep = ""), width = 1200, height = 1400)
for ( chrom in uniq.chrs ) {

	plotAllSampleCCFsAcrossGenome(asas.ph.raw.cpn.lists = asas.ph.raw.cpn.lists,
                              full.joint.segs.df = full.joint.segs.w.phasing.df[full.joint.segs.w.phasing.df$chr == chrom,] ,
                              sample.vec = sample.vec,
                              all.seg.CCF.vals.list = all.seg.CCF.list,
                              cpnX = "cpnB",
                              use.cpn.or.ccf.thresh = "CPN")


}
dev.off()

by.chrom.BAF.plot.dir <- paste(PLOT_OUTPUT_DIR, "BAF_by_chrom/", sep="")
if( ! file.exists(by.chrom.BAF.plot.dir)) {

	dir.create(by.chrom.BAF.plot.dir)

}#end if( ! file.exists(by.chrom.CCF.plot.dir)) {

#Now let's plot baf
print("Plotting by chromosome phased BAF...")
png(file = paste(by.chrom.BAF.plot.dir, "/", patient, "_BAF_phased_", "%02d", ".png", sep=""), width = 1200, height = 1400)
for ( chrom in uniq.chrs ) {

	chr.pt.ph.BAF.df.list <- patient.phased.BAF.df.list

	for ( samp in sample.vec ) {

		chr.pt.ph.BAF.df.list[[samp]] <- chr.pt.ph.BAF.df.list[[samp]][chr.pt.ph.BAF.df.list[[samp]]$chrom == chrom,]

	}#end for ( samp in sample.vec ) {

	plotAllSampleBAFsAcrossGenome(sample.vec, chr.pt.ph.BAF.df.list, cpnA.allele.color = '#de2d2699', cpnB.allele.color= '#3182bd99')


}#end for ( chrom in uniq.chrs ) {
dev.off()

testFinalOutputDataframe(full.joint.segs.df)