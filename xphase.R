#!/usr/bin/env Rscript

library(optparse)
library(digest)

# option_list <- list(
#     make_option(c("--patient"), default="CB1002",
#                 help="Patient id"),
#     make_option(c("--dip-test"), default=FALSE,
#                 help="Use dip test (default=TRUE)"),
#     make_option(c("--output-dir"), default="../output/CB1002/mphase",
#                 help="Output directory"),
#     make_option(c("--medicc-output-dir"), default="../output/CB1002/mphase",
#                 help="Output directory"),
#     make_option(c("--input-format"), default="ASCAT",
#                 help="Input data format (default=MEDICC)"),
#     make_option(c("--snp-data"), default="../input/nb-schulte//CB1002_input_paths_df.txt",
#                 help="SNP data"),
#     make_option(c("--seg-data"), default="../input/nb-schulte//CB1002.segmentation",
#                 help="Segmentation data output from ASCAT"),
#     make_option(c("--overwrite"), default=TRUE,
#                 help="Overwrite output (default=TRUE)"),
#     make_option(c("--gamma"), default=1.0,
#                 help="Gamma (for seq=1, for arrays~=0.5 (default=1.0)"),
#     make_option(c("-j", "--cores"), default=1,
#                 help="CPU cores (default=1)"),
#     make_option(c("--segment-type"), default="intersection",
#                 help="Segment type, intersection or whole_arm (default=1)"),
#     make_option(c("--genome"), default="hg38",
#                 help="Human genome version (default=hg19, also supports hg38)"),
#     make_option(c("--seed"), default=2134234,
#                 help="Random seed")
#     )
# conf <- parse_args(OptionParser(option_list = option_list), convert_hyphens_to_underscores = TRUE)
# conf$id <- substr(digest(conf), 1, 8)
# set.seed(conf$seed)#Why is a seed required? I'm not sure anything non-deterministic should be occurring

# args = commandArgs(trailingOnly=TRUE)
# print("config")
# cat(paste(names(conf), conf, collapse="\n"))

# patient <- conf$patient
# #At the moment MPhase will only try to call MSAI and new copy number in areas of the genome identified as allelically imbalanced
# #Enable this flag to try and find more allelic imbalance but this will probably break.
# allow.dip.test.phasing.flag <- as.logical(conf$dip_test)
# output.directory.base.path <- conf$output_dir
# data.format.to.use <- conf$input_format #Default for this is MEDICC, but can specify TRACERx_exome, TRACERx_snp_array or TRACERx_WGS
# #all.samp.snp.list.path <- conf$cached_data #Points towards the snp files, of whatever format they are given that should be used based on data.format to use.
# snp.file.path <- conf$snp_data #Points towards the snp files, of whatever format they are given that should be used based on data.format to use.
# seg.file.path <- conf$seg_data #Pintts towards the seg files, of whatever format they should be determined byh the data.format.to.use flag.
# overwrite.flag  <- conf$overwrite
# gamma.val <- conf$gamma #for sequencing it's 1 and SNP array it's 0.4/0.5 but this is platform specific see ASCAT help pages.
# num.cores <- conf$cores #number of cores for the do parallel package functions
# intersection.or.whole.arm <- conf$segment_type #Can either be intersection or whole_arm, do we take the only intersection of segments from each sample
# genome <- conf$genome  # hg19 or hg38
# #or do attempt to fill in the gaps between every segment in each sample so that the entire genome is covered (so output from any sample can run in GISTIC
# #as all will cover the same genomic area)

#TBKW edit...
#full.patient <- args[1]
full.patient<- "A_LTX121"
seg.file.path <- "/camp/project/proj-tracerx-lung/tracerx/release_tx421_20190411/A_LTX121/exome/MPhase/LTX121_joint_segs_all_ascat_info.txt" #Points towards the seg files, of whatever format they should be determined byh the data.format.to.use flag.
snp.file.path <- "/camp/project/proj-tracerx-lung/tracerx/release_tx421_20190411/A_LTX121/exome/MPhase/LTX121_TRACERx_input_paths_df.txt"
patient <- substring(full.patient, 3, 8)
intersection.or.whole.arm <- "intersection"
genome <- "hg19"
allow.dip.test.phasing.flag <- FALSE
data.format.to.use <- "ASCAT"
overwrite.flag <- TRUE
gamma.val <- 1
num.cores <- 1
output.directory.base.path <- paste0("/camp/project/proj-tracerx-lung/tctProjects/watkint/tx421_copy_number_qc/xphase_runs/", full.patient, "/")
system(paste0("mkdir ", output.directory.base.path))
#Right we need a free run folder:
free.run.dir <- paste0(output.directory.base.path, "free_run/")
qc.choice.dir <- paste0(output.directory.base.path, "run_qc_choice/")
rescue.run.dir <- paste0(output.directory.base.path, "rescue_runs/")
system(paste0("mkdir ", free.run.dir))
system(paste0("mkdir ", qc.choice.dir))
system(paste0("mkdir ", rescue.run.dir))


###GLOBAL FLAG that alters the behaviour of some hacked functions so that they work with MPhase and XPhase
USING_XPHASE <- FALSE 



options(stringsAsFactors = FALSE)

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
print("num.cores")
print(num.cores)

#We need to be able to select the intersection of all regions or a full GISTIC compatible union with defined regions that will always be covered
#What proportion of the chosen genome, in this case hg19, arms are covered by the segments that have been inputted.
#We also need to add in the wgii calculations


###Source relevant scripts, best to do this at the top to check for errors...
if (is.null(current.dir)) {
	current.dir <- ""
} else {
	current.dir <- paste0(current.dir, "/")
}
source(paste0(current.dir, "utility_MPhase_functions.R"))
source(paste0(current.dir, "minimum_consistent_segmentation_functions.R"))
source(paste0(current.dir, "data_input_and_QC_functions_v2.R"))
source(paste0(current.dir, "consensus_phasing_from_multiple_samples_functions_v2.R"))
source(paste0(current.dir, "generation_of_phased_and_total_cpn_functions_v2.R"))
source(paste0(current.dir, "MSAI_detection_functions.R"))
source(paste0(current.dir, "cancer_cell_fraction_calculation_functions_v2.R"))
source(paste0(current.dir, "final_output_segmentation_file_generation_functions_v3.R"))
source(paste0(current.dir, "new_seg_joining_functions.R"))
source(paste0(current.dir, "cpn_heterogeneity_plotting_functions_v2.R"))
source(paste0(current.dir, "MEDICC_input_creation_functions.R"))
source(paste0(current.dir, "ascat.R"))


#####Get relevant libraries:
#Todo check this is the only required library.
suppressWarnings(suppressMessages(library(diptest)))
suppressWarnings(suppressMessages(library(gplots)))
suppressWarnings(suppressMessages(library(gdata)))
suppressWarnings(suppressMessages(library(GenomicRanges)))
suppressWarnings(suppressMessages(library(doParallel)))
options("stringsAsFactors" = FALSE)

###############This is past the input stage and into the analysis:
#Set constants for all analyses here... eventually these can be command line set
MIN_NUM_HZ_SNP_PHASING <- 5
#PLOT_OUTPUT_DIR <- NULL
OVERWRITE_FLAG <- as.logical(overwrite.flag)
PLOTTING_FLAG <- TRUE
#Pre joining overwrite flag
OVERWRITE_FLAG <- TRUE#Load in the pre-joining results and don't regenerate them if FALSE.
#POST JOINING FLAGS to run sections of code
JOINED_SEGS_ERROR_CHECKING_FLAG <- FALSE
OVERWRITE_INIT_MCS_ADD_INFO_FLAG <- TRUE#Using joined segs on second run throughof MPhase
OVERWRITE_ARM_INFO_FOR_JOIN_INIT_MCS_FLAG <- TRUE
OVERWRITE_JOIN_INIT_MCS_FLAG <- TRUE
OVERWRITE_JOINED_CCF_CALC_FLAG <- TRUE
OVERWRITE_JOINED_MCS_ADD_INFO_FLAG <- TRUE
CALCULATE_JOINED_ARM_FOCAL_FLAG <- TRUE
CALCULATE_FOCAL_EVENT_CLONALITY_FLAG <- TRUE
GENERATE_GISTIC_INPUT_FLAG <- FALSE

PLOT_OUTPUT_DIR <- paste(output.directory.base.path, "/", sep="")
if ( ! file.exists(PLOT_OUTPUT_DIR) ) {

	print(PLOT_OUTPUT_DIR)
	print("Output directory doesn't exist.")

}
NUM_CORES <- as.numeric(num.cores)

output_dir <- PLOT_OUTPUT_DIR
if (!dir.exists(output_dir))
  dir.create(output_dir, recursive = TRUE)

registerDoParallel(cores = NUM_CORES)
# print(foreach::getDoParWorkers())
# print(getDoParName())
# print(getDoParVersion())
#
# #Test the parallelisation as sometimes on some nodes parallel jobs file.
# print("Testing parallelisation....")
# test.parallel.vec <- sample(100:1200, NUM_CORES)
# itx <- iter(test.parallel.vec)
# test.parallel.list <- foreach(num = itx, .combine = c) %dopar% {
#
#
# 	test.sqrt <- sqrt(num)
# 	returnable.list <- list(test.sqrt)
# 	names(returnable.list) <- as.character(num)
# 	returnable.list
#
# }
# #So to test this, check the length of the list
# if (length(test.parallel.vec) != length(test.parallel.list) ) {
# 	stop("ERROR, length of returned results of parallel test are not the same as the input... PARALLELISATION NOT WORKING")
# }
# for (test.value in test.parallel.vec) {
# 	if ( ! as.character(test.value) %in% names(test.parallel.list)) {
# 		stop("Error, test parallel value absent from returned results of parallel test... PARALLELISATION NOT WORKING")
# 	}
# 	if ( sqrt(test.value) != as.numeric(test.parallel.list[as.character(test.value)])) {
# 		stop("Error, value returned by parallel test is not correct... PARALLELISATION NOT WORKING")
# 	}
# }
#
# print("non parallelised input and results")
# series.test.vec <- sqrt(test.parallel.vec)
# names(series.test.vec) <- test.parallel.vec
# print(series.test.vec)
# print("parallelised results:")
# print(unlist(test.parallel.list))

#TODO why is is this variable being reassigned here....
#intersection.or.whole.arm <- "intersection"



INTERSECTION_OR_ARM <- intersection.or.whole.arm
GAMMA_VAL <- as.numeric(gamma.val)
MAX_CPN_FOR_SUBCLONALITY <- 10
MIN_CCF_TO_CONSIDER <- 0.2  # This the median raw CCF different to the integer +/- this value.
MIN_CPN_TO_CONSIDER <- 0.2  # This the median raw cpn different to the integer +/- this value.
MIN_NUM_CHRS <- 2  # This is the minimum no of chromosomes used for comparisons to other other chromosomes when calling subclonal events
PVAL_THRESH_VS_MIN_CCF <- 0.001

PVAL_THRESH_VS_OTH_CHRS_SAME_CPN <- 0.01  # T test p value threshold for trying to call subclonal copy numnber...
GAIN_LOGR_THRESH <- log2(2.5 / 2)
LOSS_LOGR_THRESH <- log2(1.5 / 2)
AMP_LOGR_THRESH <- log2(4 / 2)
DEEP_LOSS_LOGR_THRESH <- log2(0.5 / 2)
PVAL_THRESH_VS_LOGR_FOR_CPN_EVENT <- 0.01
CHRS_TO_USE_VEC <- 1:22
CONF_INT_FOR_CPNX_BAF_VS_EXP <- 0.9999#used for CCF ca;culations
PW_CPNA_CPNB_DIFF_AMPLITUDE_THRESH <- 0.001#Legacy, used for pairwise comparisons between samples
PVAL_THRESH_TO_JOIN_PH_SEGS <- 0.001#P value thresholds for  joining segements on phased raw cpn values
PVAL_THRESH_TO_JOIN_TOT_SEGS <- 0.001#Some segments don't have phased copy number but we can try and join based on the total cpn
MIN_PROP_SNP_SEP_THRESH <- 0.75 #legacy minimum  proportion of snps to be separated by phasing, if fewer, try a more complicated phasing currently deactivated 11/09/2018
PVAL_COMB_WILCOX_PH_THRESH <- 0.01  #threshold for testing whether the BAFs are separated across all the samples when comparing phasings
PVAL_DIPTEST_AI_THRESH <- 0.01#Test to decide whether there is allelic imbalance, not being used as of  11/09/2018
ALLOW_DIPTEST_ALONE_PH_FLAG <- allow.dip.test.phasing.flag
MSAI_BAF_LTNPD_PVAL_THRESHOLD <- 0.001#P value threshold to call MSAI
MSAI_CPNX_PVAL_THRESHOLD <- 0.001#Legacy, can also detect MSAI from the phased copy number not ust the BAF
PVAL_THRESH_CONVERT_GAIN_TO_AMP  <- 0.01#Legacy, but might need to be re-implemented, is an amplification significantly different to a gain in a different region?

if (is.null(output.directory.base.path))
	stop("No output path provided, aborting...")

# so now we need to read in the MSAI data:
if (!file.exists(snp.file.path))
	stop("No file at the snp.file.path provided in the second argument")

if (!file.exists(seg.file.path))
	stop("No file at the seg.file.path provided in the third argument")

# Load in the relevant genome file
if (genome == "hg19") {
  info.path <- file.path(current.dir, "annot", "ucsc_hg19_chrominfo.txt")
  access.path <- file.path(current.dir, "annot", "access-10kb.hg19.bed")
} else if (genome == "hg38") {
  info.path <- file.path(current.dir, "annot", "ucsc_hg38_chrominfo.txt")
  access.path <- file.path(current.dir, "annot", "access-10kb.hg38.bed")
} else {
  stop("Genome version not recognized: ", genome)
}
chrom.info.df <- read.table(info.path, stringsAsFactors = FALSE, sep = "\t", header = TRUE)
# We don't need this until much later, don't try to load it now
#access.df <- read.table(access.path, sep="\t", header=F, stringsAsFactors=F)

#Paths for saving/loading objects depending on OVERWRITE_FLAG
#Reading in data object paths
input.snp.file.RDS.path <- paste0(PLOT_OUTPUT_DIR, patient, "_input_snp_file.RDS")
all.samp.SNP.list.path <- paste0(PLOT_OUTPUT_DIR, patient, "_all_samp_SNP_list.RDS")
asas.SNPS.lr.and.or.baf.lists.path <- paste0(PLOT_OUTPUT_DIR, patient, "_asas_SNPS_lr_and_or_baf_lists_lists.RDS")
asas.SNPs.w.both.lr.baf.snp.lists.path <- paste0(PLOT_OUTPUT_DIR, patient, "_asas_SNPs_w_both_lr_baf_snp_lists.RDS")
min.cons.segs.df.write.path <- paste0(PLOT_OUTPUT_DIR, patient, "_min_cons_segs_df.txt")
sorted.all.bp.all.samps.segs.anno.df.path <- paste0(PLOT_OUTPUT_DIR, patient, "_sorted_all_bp_all_samps_segs_anno.df.RDS")
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

print("Finished generating paths...")

#Load in the data:
#Note that the functions to create the single dataframes that it accepts are in the data_input_and_QC_functions.R
if (OVERWRITE_FLAG == TRUE) {

	init.f.joint.segs.df <- read.table(seg.file.path, stringsAsFactors = FALSE, sep = "\t", header = TRUE)

	#Legacy, but I wanted to add in the thresholds for gain/loss etc as they were called in the TRACER NEJM 2017 paper
	#Might be useful... it is aved as part of f.joint.segs.df
	f.joint.segs.df <- reformatInitialSegmentationInput(f.joint.segs.df = init.f.joint.segs.df)

	sample.vec <- unique(f.joint.segs.df$sample)
	ploidy.est.vec <- f.joint.segs.df[ ! duplicated(f.joint.segs.df$sample), 'Ploidy']
	names(ploidy.est.vec) <- unique(f.joint.segs.df$sample)
	purity.est.vec <- f.joint.segs.df[ ! duplicated(f.joint.segs.df$sample), 'ACF']
	names(purity.est.vec) <- unique(f.joint.segs.df$sample)



	if (data.format.to.use == "MEDICC") {
		cat("data format is MEDICC, reading in sites text file...\n")

		if (file.exists(input.snp.file.RDS.path)) {
			snp.file.df <- readRDS(input.snp.file.RDS.path)
		} else {
			snp.file.df <- read.table(snp.file.path, stringsAsFactors = FALSE, sep = "\t", header = TRUE)
			saveRDS(snp.file.df, file = input.snp.file.RDS.path)
		}
		all.samp.SNP.list <- readInMEDICCInputFormatParallel(snp.file.df = snp.file.df,
															processed.orig.cpn.segs.df = f.joint.segs.df,
															sample.vec = sample.vec)
	} else if (data.format.to.use == "ASCAT") {
    # Parse a bunch of rds files containing ascat.bc objects
    pt.snp.file.paths.df <- read.delim(snp.file.path,
                              stringsAsFactors = FALSE,
                              header = TRUE,
                              row.names = "sample")
    all.samp.SNP.list <- readInASCAT(sample_table = pt.snp.file.paths.df,
                                     segs = f.joint.segs.df,
                                     samples = sample.vec)
    #TBKW edit, I'm assuming that these ascat objects will also contain the germline values and am pulling 
    #them out for this "xphase" version that requires the germline BAFs to compare the allelic imbalance in each tumour
    #region segment vs the equivalent segment in the germline. This comparison is done to try and counter the effect
    #of noise and "reference bias" that might result in false allelic imbalance. The correction that is performed to avert this
    #is primitive as of 16/08/2019 and will require rewriting as well as testing. It should be noted the presence of 
    #bias observed by allelic imbalance is from separate work I've performed and eye-balling the data, it needs to be 
    #shown in the data...

    #Making the assumption here that all regions that were analysed had 
    gl.snp.df <- readInGermlineFromASCAT(ascat_bc_rds_path = pt.snp.file.paths.df[1,]$samp_ascat_bc_path)

    # # Temporarily remove everything except a few chromosomes
    # #test_chrs <- c("chr1", "chr2")
    # test_chrs <- c(1, 2)
    # all.samp.SNP.list <- lapply(all.samp.SNP.list, function(x) { str(x); return(x[x$chrom %in% test_chrs,]) })
    # f.joint.segs.df <- f.joint.segs.df[f.joint.segs.df$chr %in% test_chrs,]

    # cat("Subset of chromosomes 1 and 2 only.\n")
    # str(all.samp.SNP.list)
    # str(f.joint.segs.df)

	} else {

		print(data.format.to.use)
		stop("Incorrect data.format.to.use input to mphase")

	}

	if (length(unique(unlist(lapply(all.samp.SNP.list, nrow)))) != 1 | length(all.samp.SNP.list) < 2) {
		print("length(all.samp.SNP.list)")
		print(length(all.samp.SNP.list))
		print("length(unique(unlist(lapply(all.samp.SNP.list, nrow))))")
    print(length(unique(unlist(lapply(all.samp.SNP.list, nrow)))))
    print(lapply(all.samp.SNP.list, nrow))

		stop("Not all sites are present in all samples' sites information read in or too few samples are present at this stage, this needs to be fixed before running MPhase.")
	}

	#Save the loaded in data:
	saveRDS(all.samp.SNP.list, file = all.samp.SNP.list.path)
	#So we need the centromere postions for each arm:
	#Add in p_end and q_start

	#TBKW  20190818 edit so here is where the joint segmentation stuff has been dumped into a function to allow multiple
	#calls. Absolutely needs refactoring. 
	sorted.all.bp.all.samps.segs.anno.df <- createFullyAnnotatedJointSegmentation(f.joint.segs.df = f.joint.segs.df)
	#TBKW edit for xphase version
	anno.seg.df <- sorted.all.bp.all.samps.segs.anno.df
	# sorted.all.bp.all.samps.segs.anno.df.path <- paste0(PLOT_OUTPUT_DIR, patient, "_sorted_all_bp_all_samps_segs_anno_df.RDS")
	# saveRDS(sorted.all.bp.all.samps.segs.anno.df, file = sorted.all.bp.all.samps.segs.anno.df.path)
	#Now we want to save this annotated version:
	cpn.chr.start.end.df <- sorted.all.bp.all.samps.segs.anno.df[sorted.all.bp.all.samps.segs.anno.df$full_sample==sorted.all.bp.all.samps.segs.anno.df[1,]$full_sample,1:3]
	mock.cpn.matrix <- matrix(data=0, nrow=nrow(cpn.chr.start.end.df), ncol = length(sample.vec))
	colnames(mock.cpn.matrix) <- sample.vec
	min.cons.segs.df <- cbind(cpn.chr.start.end.df, mock.cpn.matrix)
	#end of 201908TBKW edit all this code 


	#Minimum consistent segmentation is finished! Let's get the SNPs and their LogR and BAF values amd assign to segments.
	#All SNPs in these lists must have both a logr and baf value: all sample all seg logr baf lists is asas.lr.baf
	ptm <- proc.time()
	print("cleanAndErrorCheckSNPsParallel first")
	#source(paste0(current.dir, "data_input_and_QC_functions_v2.R"))
	asas.SNPs.w.both.lr.baf.snp.lists <- cleanAndErrorCheckSNPsParallel(all.samp.SNP.list,
																		mcs.seg.to.use.df =min.cons.segs.df,
																		allow.logr.only.flag = FALSE,
																		data.type = data.format.to.use)
	print(str(asas.SNPs.w.both.lr.baf.snp.lists))
	if (length(unique(sapply(asas.SNPs.w.both.lr.baf.snp.lists, length))) != 1) {

		stop("Not all segments are present in all samples' snp assignations, likely something has gone wrong, is the parallelisation working correctly?")

	}
	print(asas.SNPs.w.both.lr.baf.snp.lists.path)
	saveRDS(asas.SNPs.w.both.lr.baf.snp.lists, file = asas.SNPs.w.both.lr.baf.snp.lists.path)
	time.taken <- proc.time() - ptm
	print("time.taken")
	print(time.taken)
	##
	#SNPs that only have LogR are permitted in these lists, but it also contains all the SNPs with both BAF and LogR
	ptm <- proc.time()
	print("cleanAndErrorCheckSNPsParallel second")

	asas.SNPS.lr.and.or.baf.lists <- cleanAndErrorCheckSNPsParallel(all.samp.SNP.list,
																	mcs.seg.to.use.df =min.cons.segs.df,
																	allow.logr.only.flag = TRUE,
																	data.type = data.format.to.use)

	if (length(unique(sapply(asas.SNPS.lr.and.or.baf.lists, length))) != 1) {

		stop("Not all segments are present in all samples' snp assignations, likely something has gone wrong, is the parallelisation working correctly?")

	}

	#print(asas.SNPS.lr.and.or.baf.lists.path)
	saveRDS(asas.SNPS.lr.and.or.baf.lists, file = asas.SNPS.lr.and.or.baf.lists.path)
	time.taken <- proc.time() - ptm
	print("time.taken")
	print(time.taken)

	if ( length(asas.SNPs.w.both.lr.baf.snp.lists) <= 1 ) {

		print(str(asas.SNPs.w.both.lr.baf.snp.lists, max = 1))
		stop("Only one or fewer samples present after reading in data.")

	}

	#######So here is where we should check if there are any segments that should be
	#save:
	print("Saving snps files....")
	print(asas.SNPS.lr.and.or.baf.lists.path)
	saveRDS(asas.SNPS.lr.and.or.baf.lists, file = asas.SNPS.lr.and.or.baf.lists.path)
	print(asas.SNPs.w.both.lr.baf.snp.lists.path)
	saveRDS(asas.SNPs.w.both.lr.baf.snp.lists, file = asas.SNPs.w.both.lr.baf.snp.lists.path)
	print(min.cons.segs.df.path)
	saveRDS(min.cons.segs.df, file = min.cons.segs.df.path)
	print(min.cons.segs.df.write.path)
	saveRDS(ploidy.est.vec, file = ploidy.est.vec.path)
	saveRDS(purity.est.vec, file = purity.est.vec.path)
	saveRDS(f.joint.segs.df, file = processed.orig.cpn.segs.df.path)

	saveRDS(sorted.all.bp.all.samps.segs.anno.df, file = sorted.all.bp.all.samps.segs.anno.df.path)

} else if (OVERWRITE_FLAG == FALSE) {#end if (OVERWRITE_FLAG == TRUE) {

	#load original data:
	all.samp.SNP.list <- readRDS(file = all.samp.SNP.list.path)
	asas.SNPS.lr.and.or.baf.lists <- readRDS(file = asas.SNPS.lr.and.or.baf.lists.path)
	asas.SNPs.w.both.lr.baf.snp.lists <- readRDS(file = asas.SNPs.w.both.lr.baf.snp.lists.path)
	sorted.all.bp.all.samps.segs.anno.df <- readRDS( file = sorted.all.bp.all.samps.segs.anno.df.path)
	min.cons.segs.df <- readRDS(file = min.cons.segs.df.path)
	f.joint.segs.df <- readRDS(file = processed.orig.cpn.segs.df.path)
	ploidy.est.vec <- readRDS(file = ploidy.est.vec.path)
	purity.est.vec <- readRDS(file = purity.est.vec.path)

}#end else if (OVERWRITE_FLAG == FALSE

#TODO move this up
#Let's come up with the original wgii score pre MPhase based on the segmentation:
#source(paste0(current.dir, "MEDICC_input_creation_functions.R"))

sample.vec <- names(asas.SNPs.w.both.lr.baf.snp.lists)
sample.pre.Mphase.wgii.list <- list()
for (samp in sample.vec) {

	samp.orig.seg.df <- f.joint.segs.df[f.joint.segs.df$sample %in% samp & ( ! f.joint.segs.df$chr %in% c(23,24)),]
	print(head(samp.orig.seg.df))
	samp.orig.formatted.seg.df <- samp.orig.seg.df[,c("sample","chr", "startpos", "endpos", "n.het", "cnTotal", "nMajor","nMinor", "Ploidy", "ACF", "nAraw", "nBraw")]
	samp.orig.formatted.seg.df <- samp.orig.formatted.seg.df[ ! is.na(samp.orig.formatted.seg.df$nMajor),]
	print(head(samp.orig.formatted.seg.df))
	pre.mphase.samp.wgii <- calculateWGII(seg.df = samp.orig.formatted.seg.df)
	print(samp)
	print("WGII:")
	print(pre.mphase.samp.wgii)
	sample.pre.Mphase.wgii.list[[samp]] <-  pre.mphase.samp.wgii

}#end for (samp in sample.vec) {

cat("All wGII\n")
print(str(sample.pre.Mphase.wgii.list))

#Based on the input and how many BAF values there are we need to set the subclonality thresholds
total.num.of.hz.sites.per.sample <- sum(sapply(asas.SNPs.w.both.lr.baf.snp.lists[[1]], nrow))
cat("total.num.of.hz.sites.per.sample\n")
print(total.num.of.hz.sites.per.sample)
print("Setting the pval thresholds for calling subclonality depending on the number of heterozygous sites assayed:")

if(total.num.of.hz.sites.per.sample < 100000) {

	PVAL_THRESH_VS_MIN_CCF <- 0.001
	PVAL_THRESH_VS_OTH_CHRS_SAME_CPN <- 0.01

} else if (total.num.of.hz.sites.per.sample > 300000 ) {

	PVAL_THRESH_VS_MIN_CCF <- 0.00000001
	PVAL_THRESH_VS_OTH_CHRS_SAME_CPN <- 0.0000000001

}

#So now we want

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
	source(paste0(current.dir, "consensus_phasing_from_multiple_samples_functions_v2.R"))
	#all.seg.phasing.info.list <- getCombinedBAFDipTestPhasing(asas.lr.baf.lists = asas.SNPs.w.both.lr.baf.snp.lists,
	#                                            seg.name.vec = seg.name.vec,
	#                                            ploidy.est.vec = ploidy.est.vec,
	#                                            purity.est.vec = purity.est.vec,
	#                                            gamma.val = GAMMA_VAL,
	#                                            min.num.hz.snp.phasing = MIN_NUM_HZ_SNP_PHASING,
	#                                            min.prop.snp.sep.threshold = MIN_PROP_SNP_SEP_THRESH,
	#                                            dip.test.p.val.threshold = PVAL_DIPTEST_AI_THRESH,
	#                                            comb.wilcox.snp.sep.p.val.threshold = PVAL_COMB_WILCOX_PH_THRESH,
	#                                            allow.dip.test.without.ASCAT.to.phase = ALLOW_DIPTEST_ALONE_PH_FLAG
	#                	                          )
	print("Getting phasing to apply to all samples using getCombinedBAFDipTestPhasingParallel...")
	all.seg.phasing.info.list <- getCombinedBAFDipTestPhasingParallel(asas.lr.baf.lists = asas.SNPs.w.both.lr.baf.snp.lists,
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
	if(length(patient.phased.BAF.df.list) <= 1) {

		print("Error in generating phased BAF dataframes")
	}
	saveRDS(patient.phased.BAF.df.list, file = patient.phased.BAF.df.list.path)


	#So this is where want to replace the rest of mphase more or less..
	#TBKW NEW WORK:
	print("Starting new TBKW edit... begun on 16/08/2019 to recall copy number using phasing.")
	source(paste0(current.dir, "ascat.R"))
	#So determine which SNPs are hz:
	snp.file.df <- gl.snp.df
	snp.file.df$snp_pos <- paste0(snp.file.df$chrom, "_", snp.file.df$pos)
	snp.file.df$snp_idx <- 1:nrow(snp.file.df)
	#get germline column
	gl.col.name <- "BAF"#grep("_GL_BAF", colnames(snp.file.df), value=TRUE)
	logr.only.snp.df <- snp.file.df[is.na(snp.file.df[,gl.col.name]),]
	actual.snp.file.df <- snp.file.df[! is.na(snp.file.df[, gl.col.name]),]
	print("nrow(gl.snp.df)")
	print(nrow(gl.snp.df))
	print("nrow(actual.snp.file.df)")
	print(nrow(actual.snp.file.df))
	gl.homo.snp.df <- actual.snp.file.df[! is.na(actual.snp.file.df[,gl.col.name]) & (as.numeric(actual.snp.file.df[,gl.col.name]) > 0.7 | as.numeric(actual.snp.file.df[,gl.col.name]) < 0.3) ,]
	gl.hz.snp.df <- actual.snp.file.df[! is.na(actual.snp.file.df[,gl.col.name]) & (as.numeric(actual.snp.file.df[,gl.col.name]) <= 0.7 &  as.numeric(actual.snp.file.df[,gl.col.name]) >= 0.3),]
	gl.hz.snp.df$hz_snp_idx <- 1:nrow(gl.hz.snp.df)
	anno.seg.df$start_pos <- paste0(anno.seg.df$Chromosome, "_", anno.seg.df$Start_Position)
	anno.seg.df$start_is_hz <- anno.seg.df$start_pos %in% gl.hz.snp.df$snp_pos
	anno.seg.df$start_is_homo <- anno.seg.df$start_pos %in% gl.homo.snp.df$snp_pos
	anno.seg.df$start_is_logr <- anno.seg.df$start_pos %in% logr.only.snp.df$snp_pos
	anno.seg.df$end_pos <- paste0(anno.seg.df$Chromosome, "_", anno.seg.df$End_Position)
	anno.seg.df$end_is_hz <- anno.seg.df$end_pos %in% gl.hz.snp.df$snp_pos
	anno.seg.df$end_is_homo <- anno.seg.df$end_pos %in% gl.homo.snp.df$snp_pos
	anno.seg.df$end_is_logr <- anno.seg.df$end_pos %in% logr.only.snp.df$snp_pos
	#make a quick summary column
	anno.seg.df$start_is_pos <- anno.seg.df$start_is_hz | anno.seg.df$start_is_homo | anno.seg.df$start_is_logr
	anno.seg.df$end_is_pos <- anno.seg.df$end_is_hz | anno.seg.df$end_is_homo | anno.seg.df$end_is_logr
	anno.seg.df$start_snp_idx <- NA
	anno.seg.df$end_snp_idx <- NA
	anno.seg.df[anno.seg.df$start_is_pos %in% TRUE ,]$start_snp_idx <- snp.file.df[match(anno.seg.df[anno.seg.df$start_is_pos %in% TRUE ,]$start_pos, snp.file.df$snp_pos),]$snp_idx
	anno.seg.df[anno.seg.df$end_is_pos %in% TRUE ,]$end_snp_idx <- snp.file.df[match(anno.seg.df[anno.seg.df$end_is_pos %in% TRUE ,]$end_pos, snp.file.df$snp_pos),]$snp_idx
	#potential column names of use: in_intersection in_disjoin_only in_arm_setdiff
	# #For no now let's only take those those segments that are in the intersection
	# int.anno.seg.df <- anno.seg.df[anno.seg.df$in_intersection%in%TRUE & anno.seg.df$full_sample %in% sample.vec[1],]

	#Let's try adding in the disjoin...
	int.anno.seg.df <- anno.seg.df[anno.seg.df$in_intersection%in%TRUE | anno.seg.df$in_disjoin_only%in%TRUE & anno.seg.df$full_sample %in% sample.vec[1],]
	#Let's just add the arm classificaitons to all the snp positions
	snp.file.df$arm <- NA
	for (ch in c(unique(int.anno.seg.df$chrom))) {
		#get start and end of p
		p.start <- chrom.info.df[chrom.info.df$chrom %in% ch,]$p_start
		p.end <- chrom.info.df[chrom.info.df$chrom %in% ch,]$p_end
		q.start <- chrom.info.df[chrom.info.df$chrom %in% ch,]$q_start
		q.end <- chrom.info.df[chrom.info.df$chrom %in% ch,]$q_end

		#head(snp.file.df)

		#gsample$Continent <- with(gsample, ifelse(MAKE=='HOLDEN', 'AUS', Continent))
		if(nrow(snp.file.df[snp.file.df$chrom %in% ch & snp.file.df$pos >= p.start & snp.file.df$pos <= p.end,]) > 0) {

			snp.file.df[snp.file.df$chrom %in% ch & snp.file.df$pos >= p.start & snp.file.df$pos <= p.end,]$arm<- paste0(ch, "p") 

		}

		if(nrow(snp.file.df[snp.file.df$chrom %in% ch & snp.file.df$pos >= q.start & snp.file.df$pos <= q.end,]) > 0) {

			snp.file.df[snp.file.df$chrom %in% ch & snp.file.df$pos >= q.start & snp.file.df$pos <= q.end,]$arm<- paste0(ch, "q")

		}

	}#end for (chr in unique(int.anno.seg.df$chrom)) {

	#So we need to write a function that finds the nearest SNP that is
	#1) in the same chromosome arm
	#2) not already the the start or end of another segment
	#3) Not within another segment...
	#Remove any of the SNPs that are already starting or ending segments...	
	all.snp.grs <- makeGRangesFromDataFrame(snp.file.df[ ! snp.file.df$snp_pos %in% int.anno.seg.df$start_pos & ! snp.file.df$snp_pos %in% int.anno.seg.df$end_pos,],
			                         keep.extra.columns=FALSE,
			                         ignore.strand=FALSE,
			                         seqinfo=NULL,
			                         seqnames.field=c("chrom"),
			                         start.field=c("pos"),
			                         end.field=c("pos"),
			                         strand.field="strand",
			                         starts.in.df.are.0based=FALSE)

	#Get indices of SNPs with missing start or end idxs
	int.anno.seg.df$strand <- "*"
	all.int.grs <- makeGRangesFromDataFrame(int.anno.seg.df,
			                         keep.extra.columns=FALSE,
			                         ignore.strand=FALSE,
			                         seqinfo=NULL,
			                         seqnames.field=c("Chromosome"),
			                         start.field=c("Start_Position"),
			                         end.field=c("End_Position"),
			                         strand.field="strand",
			                         starts.in.df.are.0based=FALSE)
	#
	lacking.snp.idx.vec <- which(int.anno.seg.df$start_is_pos %in% FALSE | int.anno.seg.df$end_is_pos %in% FALSE)
	#FIXME, potential improvement that could be written in would be to, once a SNP is being used
	#remove it from the pool of potential snps to be used as a potenial start and end. However, this might not
	#matter to much.
	int.anno.seg.df$nearest_preceding_snp <- NA
	int.anno.seg.df$nearest_following_snp <- NA 
	print("	for (i in 1:length(all.int.grs)) {")
	print("length(all.int.grs)")
	print(length(all.int.grs))
	for (i in 1:length(all.int.grs)) {
		#Check the seqnames part as SNPs really need to be in the same need to be in the same chromosome arm...?
		if( ! i %in% lacking.snp.idx.vec) {

			next

		}
		print("lacking snp idx:")
		print(i)

		curr.arm <- int.anno.seg.df[i,]$chr_arm

		if (! is.na(precede(all.int.grs[i], all.snp.grs)) && snp.file.df[precede(all.int.grs[i], all.snp.grs),]$arm %in% curr.arm) {

			int.anno.seg.df[i,]$nearest_preceding_snp <- start(all.snp.grs[precede(all.int.grs[i], all.snp.grs)])

		}

		if (! is.na(follow(all.int.grs[i], all.snp.grs)) && snp.file.df[follow(all.int.grs[i], all.snp.grs),]$arm %in% curr.arm) {

			int.anno.seg.df[i,]$nearest_following_snp <- start(all.snp.grs[follow(all.int.grs[i], all.snp.grs)])

		}

		#Any remaining examples are going to be the first or only 
		if(is.na(int.anno.seg.df[i,]$nearest_following_snp) & is.na(int.anno.seg.df[i,]$nearest_preceding_snp)) {

			print(int.anno.seg.df[i,c("Start_Position", "End_Position", "start_snp_idx", "end_snp_idx", "nearest_preceding_snp", "nearest_following_snp")])
			int.anno.seg.df[i,]$nearest_following_snp <- start(all.snp.grs[nearest(all.int.grs[i], all.snp.grs)])
			int.anno.seg.df[i,]$nearest_preceding_snp <- start(all.snp.grs[nearest(all.int.grs[i], all.snp.grs)])

			print(precede(all.int.grs[i], all.snp.grs))
			print(follow(all.int.grs[i], all.snp.grs))

			print(int.anno.seg.df[int.anno.seg.df$chr_arm %in% curr.arm, c("Start_Position", "End_Position", "start_snp_idx", "end_snp_idx", "nearest_preceding_snp", "nearest_following_snp")])


		}#end if(is.na(int.anno.seg.df[i,]$nearest_following_snp) & is.na(int.anno.seg.df[i,]$nearest_preceding_snp)) {


	}#end for (i in 1:length(all.int.grs)) {

	int.anno.seg.df[is.na(int.anno.seg.df$start_snp_idx) | is.na(int.anno.seg.df$end_snp_idx), c("Start_Position", "End_Position", "start_snp_idx", "end_snp_idx", "nearest_preceding_snp", "nearest_following_snp") ]
	int.anno.seg.df[int.anno.seg.df$start_is_pos %in% FALSE, ]$Start_Position <- int.anno.seg.df[int.anno.seg.df$start_is_pos %in% FALSE, ]$nearest_following_snp
	int.anno.seg.df[is.na(int.anno.seg.df$Start_Position), ]$Start_Position <- int.anno.seg.df[is.na(int.anno.seg.df$Start_Position), ]$nearest_preceding_snp
	int.anno.seg.df[int.anno.seg.df$end_is_pos %in% FALSE, ]$End_Position <- int.anno.seg.df[int.anno.seg.df$end_is_pos %in% FALSE, ]$nearest_preceding_snp
	int.anno.seg.df[is.na(int.anno.seg.df$End_Position), ]$End_Position <- int.anno.seg.df[is.na(int.anno.seg.df$End_Position), ]$nearest_following_snp

	f.int.anno.seg.df <- int.anno.seg.df
	print("for (i in 2:nrow(f.int.anno.seg.df)) {")
	for (i in 2:nrow(f.int.anno.seg.df)) {

		if(f.int.anno.seg.df[ (i - 1),]$chr_arm == f.int.anno.seg.df[i,]$chr_arm) {

			if(f.int.anno.seg.df[ (i - 1),]$End_Position > f.int.anno.seg.df[i, ]$Start_Position) {

				print(f.int.anno.seg.df[(i-1):i, 1:10])

				preceding.end.pos <- f.int.anno.seg.df[ (i - 1),]$End_Position
				f.int.anno.seg.df[ (i - 1),]$End_Position <- f.int.anno.seg.df[i, ]$Start_Position
				f.int.anno.seg.df[i, ]$Start_Position <- preceding.end.pos

			}#end if(f.int.anno.seg.df[ (i - 1),]$End_Position > f.int.anno.seg.df[i, ]$Start_Position) {

		}#end if(f.int.anno.seg.df[ (i - 1),]$chr_arm == f.int.anno.seg.df[i,]$chr_arm) {

	}#end for (i in 2:nrow(f.int.anno.seg.df)) {

	f.int.anno.seg.df$start_pos <- paste0(f.int.anno.seg.df$Chromosome, "_", f.int.anno.seg.df$Start_Position)
	f.int.anno.seg.df$start_is_hz <- f.int.anno.seg.df$start_pos %in% gl.hz.snp.df$snp_pos
	f.int.anno.seg.df$start_is_homo <- f.int.anno.seg.df$start_pos %in% gl.homo.snp.df$snp_pos
	f.int.anno.seg.df$start_is_logr <- f.int.anno.seg.df$start_pos %in% logr.only.snp.df$snp_pos
	f.int.anno.seg.df$end_pos <- paste0(f.int.anno.seg.df$Chromosome, "_", f.int.anno.seg.df$End_Position)
	f.int.anno.seg.df$end_is_hz <- f.int.anno.seg.df$end_pos %in% gl.hz.snp.df$snp_pos
	f.int.anno.seg.df$end_is_homo <- f.int.anno.seg.df$end_pos %in% gl.homo.snp.df$snp_pos
	f.int.anno.seg.df$end_is_logr <- f.int.anno.seg.df$end_pos %in% logr.only.snp.df$snp_pos

	#make a quick summary column
	f.int.anno.seg.df$start_is_pos <- f.int.anno.seg.df$start_is_hz | f.int.anno.seg.df$start_is_homo | f.int.anno.seg.df$start_is_logr
	f.int.anno.seg.df$end_is_pos <- f.int.anno.seg.df$end_is_hz | f.int.anno.seg.df$end_is_homo | f.int.anno.seg.df$end_is_logr
	#So now let's assign these new snps to be the starts and ends...
	f.int.anno.seg.df[f.int.anno.seg.df$start_is_pos %in% TRUE ,]$start_snp_idx <- snp.file.df[match(f.int.anno.seg.df[f.int.anno.seg.df$start_is_pos %in% TRUE ,]$start_pos, snp.file.df$snp_pos),]$snp_idx
	f.int.anno.seg.df[f.int.anno.seg.df$end_is_pos %in% TRUE ,]$end_snp_idx <- snp.file.df[match(f.int.anno.seg.df[f.int.anno.seg.df$end_is_pos %in% TRUE ,]$end_pos, snp.file.df$snp_pos),]$snp_idx

	#So we're going to need the indices of the BAF snps, only needs to be done once as will be the same for all samples
	phased.baf.info.df <- patient.phased.BAF.df.list[[1]]
	phased.baf.info.df$snp_pos <- paste0(phased.baf.info.df$chrom, "_", phased.baf.info.df $pos)
	#Now we need to know their idx in the original data frame
	phased.baf.info.df$snp_idx <- snp.file.df[match(phased.baf.info.df$snp_pos, snp.file.df$snp_pos),]$snp_idx
	#The variable was renamed here to allow debugging but is just continues.
	all.samp.j.seg.df <- f.int.anno.seg.df

	#Note there could be an issue here if not every segment contains a BAF point...
	#But, if we're using hte intersection this should be fine? Not sure on this point.
	#However it doesn't really matter as ASCAT makes sure every segment whether or not it had any baf point has an allele-specific value.
	print("	for (samp in sample.vec) { run ASCAT...")

	#So now there are several things we want to do here depending on the options passed to the tool
	#1) Using the phasing information obtained from all samples, allow a free run allowing ASCAT to call whatever purity and ploidy of all samples included
	#in the run it feels is most appropriate. This will potentially provide new and improved solutions that might be better than what was fed into mphase
	#particularly if the segmentation of the input did not make use of multi-region phasing data.
	#2)  
	#Initiate list for segmentation
	ph.ascat.output.list <- list()
	ph.ascat.segs.df.list <- list()
	for (samp in sample.vec) {
		print(samp)
		#load(pt.snp.file.paths.df[pt.snp.file.paths.df$sample %in% samp,]$samp_ascat_bc_path)
		#Here we load in the file path of the relevant samples ascat.bc object. For compatibility with Matt's code
		#This should accept the ascat.bc files in eithe RDS or RData format.
		load(pt.snp.file.paths.df[samp,]$samp_ascat_bc_path)

		if(grepl("RData", pt.snp.file.paths.df[samp,]$samp_ascat_bc_path)) {

		  print("ascat.bc object in RData file, loading...")
		  load(pt.snp.file.paths.df[samp,]$samp_ascat_bc_path)
		  #This will have loaded the ascat.bc object in

		} else {

		  ascat.bc <- readRDS(pt.snp.file.paths.df[samp,]$samp_ascat_bc_path)

		}

		#So let's get the phased BAF:
		samp.pha.baf.df <- patient.phased.BAF.df.list[[samp]]
		#So we're going to need the indices of the BAF snps...
		samp.pha.baf.df$snp_pos <- phased.baf.info.df$snp_pos
		samp.pha.baf.df$snp_idx <- phased.baf.info.df$snp_idx
		#So now we need the segments that refer to this sample (as these are joint segments
		#these should all be same for each sample regardless.)
		samp.seg.df <- all.samp.j.seg.df#[all.samp.j.seg.df$full_sample %in% samp,]#No need for everything to be reset as we have one set of joint segmentation.
		#Let's do the LogR
		tum.logr.vec <- ascat.bc[["Tumor_LogR"]][,1]
		tum.logr.segmented.vec <- ascat.bc[["Tumor_LogR_segmented"]][,1]

		for (i in 1:nrow(samp.seg.df)){
		  
		  tum.logr.segmented.vec[samp.seg.df[i,]$start_snp_idx:samp.seg.df[i,]$end_snp_idx] <- mean(tum.logr.vec[samp.seg.df[i,]$start_snp_idx:samp.seg.df[i,]$end_snp_idx])
		  
		}#end for (i in 1:nrow(rle.baf.pcf.df)){
		#Let's do the BAF recalculation
		tum.baf.vec <- ascat.bc[["Tumor_BAF"]][,1]
		tum.baf.segmented.vec <- ascat.bc[["Tumor_BAF_segmented"]][[1]][,1]
		#So now this is a little more complicated as the BAF segmented is only in the SNPs
		#Would be nice to make the compaison with the germline here as well. We have the data so could...
		for (i in 1:nrow(samp.seg.df)){ 

		  print(samp.seg.df[i,])
		  #So first get the indices of the BAF indices within the segment
		  seg.hz.snp.df <- gl.hz.snp.df[gl.hz.snp.df$snp_idx <= samp.seg.df[i,]$end_snp_idx & gl.hz.snp.df$snp_idx >= samp.seg.df[i,]$start_snp_idx,]
		  if( nrow(seg.hz.snp.df) == 0) {

		  	next

		  }#end if( nrow(seg.hz.snp.df) == 0) {
		  #Get indices of hz SNPs that are present ... so we need some form of mapping here to the smaller
		  #These data contain all samples from this patient:
		  seg.min.hz.snp.idx<- seg.hz.snp.df[1,]$hz_snp_idx
		  seg.max.hz.snp.idx<- seg.hz.snp.df[nrow(seg.hz.snp.df),]$hz_snp_idx
		  #Let's get the germline baf data specfically so that we can perform a crude correction
		  gl.col.name <- grep("BAF", colnames(seg.hz.snp.df))
		  gl.seg.hz.snp.df <- seg.hz.snp.df[, c(1, 2, gl.col.name)]
		  colnames(gl.seg.hz.snp.df)[3] <- "BAF"
		  #Sample specific!
		  #So let's get the SNPs that are within the phased baf df
		  seg.ph.hz.snp.df <- samp.pha.baf.df[samp.pha.baf.df$snp_idx <= samp.seg.df[i,]$end_snp_idx & samp.pha.baf.df$snp_idx >= samp.seg.df[i,]$start_snp_idx,]

		  if(nrow(seg.ph.hz.snp.df) != nrow(gl.seg.hz.snp.df)) {

		  	print("Number of phased SNPs is different to the actual number hz snps present in seg, perhaps to few to phase by MPhase?")
		  	print("nrow(seg.ph.hz.snp.df)")
		  	print(nrow(seg.ph.hz.snp.df) )

		  	print("nrow(gl.seg.hz.snp.df)")
		  	print(nrow(gl.seg.hz.snp.df))
		  	next
		  }#end if(nrow(seg.ph.hz.snp.df) != nrow(gl.seg.hz.snp.df)) {
		  #GL specific
		  gl.seg.hz.snp.df <- cbind(gl.seg.hz.snp.df, seg.ph.hz.snp.df[, c("maj_is_cpnA_logic_vec", "maj_is_cpnB_logic_vec", "snp_pos", "snp_idx")])
		  if(nrow(seg.hz.snp.df) == 0 | nrow(seg.ph.hz.snp.df) == 0){

		  	print("No heterozygous SNPs present within segment.")
		  	print(samp.seg.df[i,])

		  }
		  if (nrow(seg.hz.snp.df) != nrow(seg.ph.hz.snp.df)) {

		  	print("Differing numbers of hz snps in the original hz snp df and the phased. Not sure why?")
		  	print(samp.seg.df[i,])

		  }
		  #Ok so...
		  #Let's get the phasing and calculate difference for the sample and the germline
		  #We're assuming that these are in the same order at this point:
		  #germline segmentation value:
		  #Get groups of indices...
		  haplotype.A.snp.pos.vec <- seg.ph.hz.snp.df[seg.ph.hz.snp.df$maj_is_cpnA_logic_vec %in% 1,]$snp_pos
		  m.hap.A.BAF <- suppressWarnings(mean(seg.ph.hz.snp.df[seg.ph.hz.snp.df$snp_pos %in% haplotype.A.snp.pos.vec,]$BAF, na.rm=TRUE))
		  haplotype.B.snp.pos.vec <- seg.ph.hz.snp.df[seg.ph.hz.snp.df$maj_is_cpnA_logic_vec %in% 0,]$snp_pos
		  m.hap.B.BAF <- suppressWarnings(mean(seg.ph.hz.snp.df[seg.ph.hz.snp.df$snp_pos %in% haplotype.B.snp.pos.vec,]$BAF, na.rm=TRUE))

		  #Now let's do the same for the germline data.. this actually isn't necessary..
		  gl.m.hap.A.BAF <- suppressWarnings(mean(gl.seg.hz.snp.df[gl.seg.hz.snp.df$snp_pos %in% haplotype.A.snp.pos.vec,]$BAF, na.rm=TRUE)) - 0.5
		  gl.m.hap.B.BAF <- suppressWarnings(mean(gl.seg.hz.snp.df[gl.seg.hz.snp.df$snp_pos %in% haplotype.B.snp.pos.vec,]$BAF, na.rm=TRUE)) - 0.5

		  #So we'll actually work out the gl correction factor here?

		  #Let's get the equivalent germline snps 
		  #Default should be to make the BAF segmented 0.5 if we don't have any info w/regards to hz snps
		  seg.md.ph.BAF.vec<-c()
		  if(! is.na(m.hap.A.BAF) && m.hap.A.BAF > 0.5) {
			#mirror
			seg.md.ph.BAF.vec <- c(seg.md.ph.BAF.vec, abs(seg.ph.hz.snp.df[seg.ph.hz.snp.df$snp_pos %in% haplotype.A.snp.pos.vec,]$BAF - 1) - gl.m.hap.A.BAF)
			##We could attempt to add in a germline correction here...

		  } else if (! is.na(m.hap.A.BAF) && m.hap.A.BAF < 0.5) {

		  	seg.md.ph.BAF.vec <- c(seg.md.ph.BAF.vec, seg.ph.hz.snp.df[seg.ph.hz.snp.df$snp_pos %in% haplotype.A.snp.pos.vec,]$BAF - gl.m.hap.A.BAF)
		  	##We could attempt to add in a germline correction here...

		  }

		  if(! is.na(m.hap.B.BAF) && m.hap.B.BAF > 0.5) {
		  	#mirror
		  	seg.md.ph.BAF.vec <- c(seg.md.ph.BAF.vec, abs(seg.ph.hz.snp.df[seg.ph.hz.snp.df$snp_pos %in% haplotype.B.snp.pos.vec,]$BAF - 1) - gl.m.hap.B.BAF)
		  	##We could attempt to add in a germline correction here...

		  } else if (! is.na(m.hap.B.BAF) && m.hap.B.BAF < 0.5) {

		    seg.md.ph.BAF.vec <- c(seg.md.ph.BAF.vec, seg.ph.hz.snp.df[seg.ph.hz.snp.df$snp_pos %in% haplotype.B.snp.pos.vec,]$BAF - gl.m.hap.B.BAF)
		    ##We could attempt to add in a germline correction here...

		  }
		  #Format
		  if(length(seg.md.ph.BAF.vec) > 0) {

		  	seg.BAF.segmented.val <- rep(mean(seg.md.ph.BAF.vec, na.rm=TRUE), length(seg.min.hz.snp.idx:seg.max.hz.snp.idx))

		  } else {

		  	seg.BAF.segmented.val <- rep(0.5, length(seg.min.hz.snp.idx:seg.max.hz.snp.idx)) 

		  }

		  #Now let's update the segmented baf values with the correct indices in the vector
		  #So to do this we need to check the actual
		  tum.baf.segmented.vec[seg.min.hz.snp.idx:seg.max.hz.snp.idx] <- seg.BAF.segmented.val

		  #Now we could edit this to compare to the germline and work out whether we should actually call this as separate....

		}#end 	for (i in 1:nrow(samp.seg.df)){ 

		#First let's do the free runs:
		ascat.bc[["Tumor_LogR_segmented"]][,1] <- tum.logr.segmented.vec
		ascat.bc[["Tumor_BAF_segmented"]][[1]][,1] <- tum.baf.segmented.vec

		#End of section generating the spoofed ascat.bc object that utilises the phasing information. This should be made into a function.

		# #TODO there may be issues when no soution that passes ASCAT's goodness of fit thresholds.
		# #This will prevent any output in the runs that don't have a specified purity and ploidy. 
		# print("Doing the free run ASCAT steps....")

		# if (purity.est.vec[samp]  < 0.15) {

		# 	MINRHO <- 0.04#This is normally set within the ASCAT script but here I'm modifying based on the input...

		# } else {

		# 	MINRHO <- 0.10#These should probably be constants set at the to of the script.
		# }

		# setwd(free.run.dir)
		# ascat.plotSegmentedData(ascat.bc)
		# ascat.output = ascat.runAscat(ascat.bc,gamma=1)#, rho_manual = rho_solution, psi_manual = psi_solution)
		# save(ascat.bc, file=paste(samp, 'ascat.bc.RData',sep='.'))
		# save(ascat.output, file=paste(samp, 'ascat.output.RData',sep='.'))
		# #Calculate number of heterozygous positions per segment
		# tmp <- ascat.output$segments_raw
		# n.het <- 0#FIXME should really ad
		# acf <- ascat.output$aberrantcellfraction
		# ploidy <- ascat.output$ploidy
		
		# tmp <- cbind(tmp[,c(1:4)], n.het, tmp[,5] + tmp[,6], tmp[,5:6], rep(ploidy, nrow(tmp)), rep(acf, nrow(tmp)),tmp[,7],tmp[,8])
		# colnames(tmp)[c(6,9:12)] <- c('cnTotal', 'Ploidy', 'ACF', 'nAraw','nBraw')
		# chr <- sub('chr','', tmp[,2])
		# chr[chr %in% c('x','X')] <- 23
		# chr[chr %in% c('y', 'Y')] <- 24
		# tmp[,2] <- as.numeric(chr)
		# #add a column defining gain and losslength
		# callState <- rep('neutral', nrow(tmp))
		# callState[log2((tmp[,'nAraw']+tmp[,'nBraw'])/tmp[,'Ploidy']) >= log2(2.5/2)] <- 'Gain'
		# callState[log2((tmp[,'nAraw']+tmp[,'nBraw'])/tmp[,'Ploidy']) <= log2(1.5/2)] <- 'Loss'
		# callState[(tmp[,'nAraw']+tmp[,'nBraw']) >= ((2*tmp[,'Ploidy'])+1)] <- 'Amp'
		# callState[(tmp[,'nAraw']+tmp[,'nBraw']) <= 0.5] <- 'Del'

		# tmp <- cbind(tmp,callState)
		# seg <- tmp
		# save(seg, file=paste(samp,'ascat.seg.RData',sep='.'))

		# rm(ascat.output)
		# rm(seg)

		#So what we want to do here is to generate all of the solutions that ASCAT proposed using the phased data... 
		#TODO - pull out the code that saves the solutions
		#and get a loop that reruns ASCAT using all of them, each in its own directory.


		#OK so now let's do the fixed solution run... 
		setwd(qc.choice.dir)
		samp.ploidy <- ploidy.est.vec[samp]
		samp.purity <- purity.est.vec[samp]


		ascat.plotSegmentedData(ascat.bc)
		ascat.output = ascat.runAscat(ascat.bc,gamma=1, rho_manual = samp.purity, psi_manual = samp.ploidy)
		ph.ascat.output.list[[samp]] <- ascat.output
		save(ascat.bc, file=paste(samp, 'ascat.bc.RData',sep='.'))
		save(ascat.output, file=paste(samp, 'ascat.output.RData',sep='.'))
		#We need to work out where the MSAI is and flip it as well...
		#Calculate number of heterozygous positions per segment
		tmp <- ascat.output$segments_raw
		n.het <- 0#FIXME probably want to add these in.
		acf <- ascat.output$aberrantcellfraction
		ploidy <- ascat.output$ploidy
		
		tmp <- cbind(tmp[,c(1:4)], n.het, tmp[,5] + tmp[,6], tmp[,5:6], rep(ploidy, nrow(tmp)), rep(acf, nrow(tmp)),tmp[,7],tmp[,8])
		colnames(tmp)[c(6,9:12)] <- c('cnTotal', 'Ploidy', 'ACF', 'nAraw','nBraw')
		chr <- sub('chr','', tmp[,2])
		chr[chr %in% c('x','X')] <- 23
		chr[chr %in% c('y', 'Y')] <- 24
		tmp[,2] <- as.numeric(chr)

		#add a column defining gain and losslength
		callState <- rep('neutral', nrow(tmp))
		callState[log2((tmp[,'nAraw']+tmp[,'nBraw'])/tmp[,'Ploidy']) >= log2(2.5/2)] <- 'Gain'
		callState[log2((tmp[,'nAraw']+tmp[,'nBraw'])/tmp[,'Ploidy']) <= log2(1.5/2)] <- 'Loss'
		callState[(tmp[,'nAraw']+tmp[,'nBraw']) >= ((2*tmp[,'Ploidy'])+1)] <- 'Amp'
		callState[(tmp[,'nAraw']+tmp[,'nBraw']) <= 0.5] <- 'Del'

		samp.ph.ascat.segs.df <- cbind(tmp,callState)
		ph.ascat.segs.df.list[[samp]] <- samp.ph.ascat.segs.df
		#This ren
		seg <- samp.ph.ascat.segs.df
		save(seg, file=paste(samp,'ascat.seg.RData',sep='.'))
		#rm(ascat.output)
		#rm(seg)

		print("Next sample....")
	}#for (samp in sample.vec) {

}#end if (OVERWRITE_FLAG == TRUE) {


if (OVERWRITE_FLAG == TRUE) {
	
	#So now we're going to construct a list of dataframes that for every heterozygous SNP has
	#chrom
	#pos
	#BAF
	#LogR
	#ph_ascat_nAraw_w_minus
	#ph_ascat_nBraw_w_minus
	#ph_ascat_nAraw
	#ph_ascat_nBraw
	#for for each single sample separately:
	#This set of raw cpn will be queried by the phasing functions.
	#I'm now in tow minds as to whether ascat needs to output total copy number...
	#This list is being used as a global innappropriately at the moment FIXME!
	USING_XPHASE <- TRUE 
	PH_ALL_SAMP_ASCAT_RAW_CPN_DF_LIST <- list()
	for (samp in sample.vec) {

		load(pt.snp.file.paths.df[samp,]$samp_ascat_bc_path)
		if(grepl("RData", pt.snp.file.paths.df[samp,]$samp_ascat_bc_path)) {

		  print("ascat.bc object in RData file, loading...")
		  load(pt.snp.file.paths.df[samp,]$samp_ascat_bc_path)
		  #This will have loaded the ascat.bc object in

		} else {

		  ascat.bc <- readRDS(pt.snp.file.paths.df[samp,]$samp_ascat_bc_path)

		}
		normal_baf <- ascat.bc$Germline_BAF[[1]]
		is_homozygous <- normal_baf < 0.3 | normal_baf > 0.7
		is_heterozygous <- !(is.na(normal_baf) | is_homozygous)
		cat("# of heterozygous SNPs:", sum(is_heterozygous), "\n")

		# Get hz SNPs
		baf <- ascat.bc$Germline_BAF[[1]]
		baf[!is_heterozygous] <- NA

		chr.vec <- as.numeric(gsub("Y", "24", gsub("X", "23", gsub("chr", "", as.character(ascat.bc$SNPpos$chr)))))
		ph.ascat.nAraw.w.minus.vec <- ph.ascat.output.list[[samp]]$nArawfull$nArawfull
		ph.ascat.nBraw.w.minus.vec <- ph.ascat.output.list[[samp]]$nBrawfull$nBrawfull
		ph.ascat.nAraw <- ph.ascat.nAraw.w.minus.vec
		ph.ascat.nAraw[ph.ascat.nAraw < 0] <- 0
		ph.ascat.nBraw <- ph.ascat.nBraw.w.minus.vec
		ph.ascat.nBraw[ph.ascat.nBraw < 0] <- 0
		#Because the plotting requires flipped values make the phasing
		#and MSAI visible...
		new.ph.ascat.nAraw <- ph.ascat.nAraw
		new.ph.ascat.nBraw <- ph.ascat.nBraw
		new.ph.ascat.nAraw.w.minus.vec <- ph.ascat.nAraw.w.minus.vec
		new.ph.ascat.nBraw.w.minus.vec <- ph.ascat.nBraw.w.minus.vec
		# replace.idx.vec <- 1:length(ph.ascat.nAraw) %% 2 == 0
		# new.ph.ascat.nAraw[replace.idx.vec] <- ph.ascat.nBraw[replace.idx.vec]
		# new.ph.ascat.nBraw[replace.idx.vec] <- ph.ascat.nAraw[replace.idx.vec]
		# new.ph.ascat.nAraw.w.minus.vec[replace.idx.vec] <- ph.ascat.nBraw.w.minus.vec[replace.idx.vec]
		# new.ph.ascat.nBraw.w.minus.vec[replace.idx.vec] <- ph.ascat.nAraw.w.minus.vec[replace.idx.vec]

		PH_ALL_SAMP_ASCAT_RAW_CPN_DF_LIST[[samp]] <- data.frame(
			chr = chr.vec[is_heterozygous],
			pos = ascat.bc$SNPpos$pos[is_heterozygous],
			snp_pos = paste0(chr.vec[is_heterozygous], "_",ascat.bc$SNPpos$pos[is_heterozygous]),
			germline_BAF = ascat.bc$Germline_BAF$BAF[is_heterozygous],
			BAF = ascat.bc$Tumor_BAF$BAF[is_heterozygous],
			LogR = ascat.bc$Tumor_LogR[[1]][is_heterozygous],
			ph_ascat_nAraw_w_minus = new.ph.ascat.nAraw.w.minus.vec,
			ph_ascat_nBraw_w_minus = new.ph.ascat.nBraw.w.minus.vec,
			ph_ascat_nAraw = new.ph.ascat.nAraw,
			ph_ascat_nBraw = new.ph.ascat.nBraw,
			stringsAsFactors = FALSE)

	}#end for (samp in sample.vec) {




	#Paths for phased versions of data:
	ph.asas.SNPS.lr.and.or.baf.lists.path <- paste0(PLOT_OUTPUT_DIR, patient, "_ph_asas_SNPS_lr_and_or_baf_lists_lists.RDS")
	ph.asas.SNPs.w.both.lr.baf.snp.lists.path <- paste0(PLOT_OUTPUT_DIR, patient, "_ph_asas_SNPs_w_both_lr_baf_snp_lists.RDS")
	ph.min.cons.segs.df.write.path <- paste0(PLOT_OUTPUT_DIR, patient, "_ph_min_cons_segs_df.txt")
	ph.sorted.all.bp.all.samps.segs.anno.df.path <- paste0(PLOT_OUTPUT_DIR, patient, "_ph_sorted_all_bp_all_samps_segs_anno.df.RDS")
	ph.min.cons.segs.df.path <- paste0(PLOT_OUTPUT_DIR, patient, "_ph_min_cons_segs_df.RDS")
	ph.processed.orig.cpn.segs.df.path <- paste0(PLOT_OUTPUT_DIR, patient, "_ph_processed_orig_cpn_segs_df.RDS")
	ph.ploidy.est.vec.path <- paste0(PLOT_OUTPUT_DIR, patient, "_ph_ploidy_est_vec.RDS")
	ph.purity.est.vec.path <- paste0(PLOT_OUTPUT_DIR, patient, "_ph_purity_est_vec.RDS")
	#TODO add save these files when generated.


	#New df contain phased ASCAT segmentation for all samples. 
	ph.init.f.joint.segs.df <- do.call(rbind, ph.ascat.segs.df.list)
	ph.f.joint.segs.df <- reformatInitialSegmentationInput(f.joint.segs.df = ph.init.f.joint.segs.df)
	ph.sorted.all.bp.all.samps.segs.anno.df <- createFullyAnnotatedJointSegmentation(f.joint.segs.df = ph.f.joint.segs.df)

	#we need to update the purity and ploidy for each sample:
	ph.purity.est.vec <- c()
	ph.ploidy.est.vec <- c()
	for (samp in sample.vec) {

		ph.purity.est.vec <- c(ph.purity.est.vec, unique(ph.f.joint.segs.df[ph.f.joint.segs.df$sample %in%samp,]$ACF))
		ph.ploidy.est.vec <- c(ph.ploidy.est.vec, unique(ph.f.joint.segs.df[ph.f.joint.segs.df$sample %in%samp,]$Ploidy))

	}#end for (samp in sample.vec) {
	names(ph.purity.est.vec) <- sample.vec
	names(ph.ploidy.est.vec) <- sample.vec

	ph.cpn.chr.start.end.df <- ph.sorted.all.bp.all.samps.segs.anno.df[ph.sorted.all.bp.all.samps.segs.anno.df$full_sample==ph.sorted.all.bp.all.samps.segs.anno.df[1,]$full_sample,1:3]
	ph.mock.cpn.matrix <- matrix(data=0, nrow=nrow(ph.cpn.chr.start.end.df), ncol = length(sample.vec))
	colnames(ph.mock.cpn.matrix) <- sample.vec
	ph.min.cons.segs.df <- cbind(ph.cpn.chr.start.end.df, ph.mock.cpn.matrix)

	ph.seg.name.vec <- paste(ph.min.cons.segs.df[, 1], ph.min.cons.segs.df[, 2], ph.min.cons.segs.df[, 3], sep = ":")

	ptm <- proc.time()
	print("Phased cleanAndErrorCheckSNPsParallel first")
	#source(paste0(current.dir, "data_input_and_QC_functions_v2.R"))
	ph.asas.SNPs.w.both.lr.baf.snp.lists <- cleanAndErrorCheckSNPsParallel(all.samp.SNP.list,
																		mcs.seg.to.use.df =ph.min.cons.segs.df,
																		allow.logr.only.flag = FALSE,
																		data.type = data.format.to.use)
	print(str(ph.asas.SNPs.w.both.lr.baf.snp.lists))
	if (length(unique(sapply(ph.asas.SNPs.w.both.lr.baf.snp.lists, length))) != 1) {

		stop("Not all segments are present in all samples' snp assignations, likely something has gone wrong, is the parallelisation working correctly?")

	}
	print(ph.asas.SNPs.w.both.lr.baf.snp.lists.path)
	saveRDS(ph.asas.SNPs.w.both.lr.baf.snp.lists, file = ph.asas.SNPs.w.both.lr.baf.snp.lists.path)
	time.taken <- proc.time() - ptm
	print("time.taken")
	print(time.taken)
	##
	#SNPs that only have LogR are permitted in these lists, but it also contains all the SNPs with both BAF and LogR
	ptm <- proc.time()
	print("Phased cleanAndErrorCheckSNPsParallel second")

	ph.asas.SNPS.lr.and.or.baf.lists <- cleanAndErrorCheckSNPsParallel(all.samp.SNP.list,
																	mcs.seg.to.use.df =ph.min.cons.segs.df,
																	allow.logr.only.flag = TRUE,
																	data.type = data.format.to.use)

	if (length(unique(sapply(ph.asas.SNPS.lr.and.or.baf.lists, length))) != 1) {

		stop("Not all segments are present in all samples' snp assignations, likely something has gone wrong, is the parallelisation working correctly?")

	}

	print(ph.asas.SNPS.lr.and.or.baf.lists.path)
	saveRDS(ph.asas.SNPS.lr.and.or.baf.lists, file = ph.asas.SNPS.lr.and.or.baf.lists.path)
	time.taken <- proc.time() - ptm
	print("time.taken")
	print(time.taken)

	if ( length(ph.asas.SNPS.lr.and.or.baf.lists) <= 1 ) {

		print(str(ph.asas.SNPS.lr.and.or.baf.lists.path, max = 1))
		stop("Only one or fewer samples present after reading in data.")

	}



	print("Getting phasing to apply to all samples using getCombinedBAFDipTestPhasingParallel...")
	ph.all.seg.phasing.info.list <- getCombinedBAFDipTestPhasingParallel(asas.lr.baf.lists = ph.asas.SNPs.w.both.lr.baf.snp.lists,
	                                            seg.name.vec = ph.seg.name.vec,
	                                            ploidy.est.vec = ph.ploidy.est.vec,
	                                            purity.est.vec = ph.purity.est.vec,
	                                            gamma.val = GAMMA_VAL,
	                                            min.num.hz.snp.phasing = MIN_NUM_HZ_SNP_PHASING,
	                                            min.prop.snp.sep.threshold = MIN_PROP_SNP_SEP_THRESH,
	                                            dip.test.p.val.threshold = PVAL_DIPTEST_AI_THRESH,
	                                            comb.wilcox.snp.sep.p.val.threshold = PVAL_COMB_WILCOX_PH_THRESH,
	                                            allow.dip.test.without.ASCAT.to.phase = ALLOW_DIPTEST_ALONE_PH_FLAG
	                                            )

	#saveRDS(ph.all.seg.phasing.info.list, file = ph.all.seg.phasing.info.list.path)

	#Todo unsure of what to do with the baf logr lists which have positions that only have logr here, at the moment exclude them but how
	#does that affect later run
	print("Generating phased BAF dataframes...")
	ph.patient.phased.BAF.df.list <- createAllSegsPhasedBAFDataFrame(asas.lr.baf.lists = ph.asas.SNPs.w.both.lr.baf.snp.lists,
                                                         sample.to.phase.seg.list = ph.all.seg.phasing.info.list,
                                                         min.num.hz.snp.phasing = MIN_NUM_HZ_SNP_PHASING
                                                         )

	
}#end if (OVERWRITE_FLAG == TRUE) {

#So now we havve generated all input information based on the updated phased ASCAT run, let's reset the
#vairable names to allow the subclonality code to continue to work:
unph.f.joint.segs.df <- f.joint.segs.df
f.joint.segs.df <- ph.f.joint.segs.df

unph.sorted.all.bp.all.samps.segs.anno.df <- sorted.all.bp.all.samps.segs.anno.df
sorted.all.bp.all.samps.segs.anno.df <- ph.sorted.all.bp.all.samps.segs.anno.df

unph.purity.est.vec <- purity.est.vec
purity.est.vec <- ph.purity.est.vec

unph.ploidy.est.vec <- ploidy.est.vec
ploidy.est.vec <- ph.ploidy.est.vec

unph.min.cons.segs.df <- min.cons.segs.df
min.cons.segs.df <- ph.min.cons.segs.df

unph.seg.name.vec <- seg.name.vec
seg.name.vec <- ph.seg.name.vec

unph.asas.SNPs.w.both.lr.baf.snp.lists <- asas.SNPs.w.both.lr.baf.snp.lists
asas.SNPs.w.both.lr.baf.snp.lists <- ph.asas.SNPs.w.both.lr.baf.snp.lists

unph.asas.SNPS.lr.and.or.baf.lists <- asas.SNPS.lr.and.or.baf.lists
asas.SNPS.lr.and.or.baf.lists <- ph.asas.SNPS.lr.and.or.baf.lists

unph.ph.patient.phased.BAF.df.list <- patient.phased.BAF.df.list
patient.phased.BAF.df.list <- ph.patient.phased.BAF.df.list


if (OVERWRITE_FLAG == TRUE) {

	#Detect MSAI in the original BAF before it is processed alongside logr to get the Raw allele specifc cpn to be used in the rest of the analysis.
	print("Performing LTNPD MSAI detection on the BAF...")
	all.seg.MSAI.by.BAF.ltnpd.present.list <- getAllSegsMSAIDetectionFromBAFListLTNPD(seg.name.vec = seg.name.vec,
																				  sample.vec  = sample.vec,
																				  asas.SNPS.lr.and.or.baf.lists = asas.SNPS.lr.and.or.baf.lists,
																				  all.seg.phasing.info.list = all.seg.phasing.info.list,
																				  MSAI.BAF.LTPND.pval.thresh = MSAI_BAF_LTNPD_PVAL_THRESHOLD)
	#saveRDS(ph.all.seg.MSAI.by.BAF.ltnpd.present.list, file = ph.all.seg.MSAI.by.BAF.ltnpd.present.list.path)


	if (PLOTTING_FLAG == TRUE) {

		png(file = paste(PLOT_OUTPUT_DIR, patient, "_MSAI_detection_BAF_phased", ".png", sep=""), width = 1200, height = 1400)

		print("Plotting phased BAF across the genome...")
		plotAllSampleBAFsAndMSAIAcrossGenome(sample.vec = sample.vec,
											all.seg.MSAI.by.BAF.ltnpd.present.list = all.seg.MSAI.by.BAF.ltnpd.present.list,
											patient.phased.BAF.df.list = patient.phased.BAF.df.list,
											cpnA.allele.color = '#de2d2699',
											cpnB.allele.color= '#3182bd99')
		dev.off()

	} else if (PLOTTING_FLAG == FALSE) {

		print("Not plotting phased BAF as flag is FALSE")

	}

	print("Applying LTNPD MSAI detection to the original segments...")
	#So now we have the original ASCAT segmentation and we want to make the MEDICC input from it
	#So now we want to flip segments where MSAI has been detected and they were separated in the original data:
	#So we also want to flip the raw values:
	intersection.orig.mcs.df  <- sorted.all.bp.all.samps.segs.anno.df[sorted.all.bp.all.samps.segs.anno.df$in_intersection==TRUE,]
	intersection.seg.name.vec <- paste0(intersection.orig.mcs.df$Chromosome,
												":",
												intersection.orig.mcs.df$Start_Position,
												":",
												intersection.orig.mcs.df$End_Position)
	orig.intersection.mcs.w.MSAI.df <- intersection.orig.mcs.df
	unq.intersection.seg.name.vec <- unique(intersection.seg.name.vec)
	for (seg.name in unq.intersection.seg.name.vec) {

		#"5:24593514:36986185" MSAI in A_LTX121
		all.samp.matching.seg.vec <- paste0(sample.vec, "@", seg.name)

		#print(all.seg.MSAI.by.BAF.ltnpd.present.list[[seg.name]])

		if (all.seg.MSAI.by.BAF.ltnpd.present.list[[seg.name]]$MSAI_detection_feasible %in% TRUE) {

			#print("here")
			seg.MPhase.MSAI.df <- all.seg.MSAI.by.BAF.ltnpd.present.list[[seg.name]]$all_samps_BAF_and_ltnpd_df
			#print(seg.MPhase.MSAI.df)
			seg.MPhase.MSAI.ASCAT.df <- cbind(seg.MPhase.MSAI.df, intersection.orig.mcs.df[paste0(sample.vec, ":", seg.name),])
			seg.MPhase.MSAI.ASCAT.df$samp_seg_name <- paste0(sample.vec, ":", seg.name)
			seg.MPhase.MSAI.ASCAT.df$balanced  <- seg.MPhase.MSAI.ASCAT.df$nMinor == seg.MPhase.MSAI.ASCAT.df$nMajor
			AI.seg.MPhase.MSAI.ASCAT.df <- seg.MPhase.MSAI.ASCAT.df[seg.MPhase.MSAI.ASCAT.df$balanced == FALSE,]

			#Need to check MSAI is ocurring in the original segmentation allelic imbalance calls
			#Why checking both of these... FIXME...
			if ( FALSE %in% AI.seg.MPhase.MSAI.ASCAT.df$MSAI_pval_passed_thresh & TRUE %in% AI.seg.MPhase.MSAI.ASCAT.df$MSAI_pval_passed_thresh ) {

				switched.samp.seg.names.vec <- AI.seg.MPhase.MSAI.ASCAT.df[AI.seg.MPhase.MSAI.ASCAT.df$MSAI_pval_passed_thresh == TRUE & AI.seg.MPhase.MSAI.ASCAT.df$poss_dominant_allele == "cpnB" ,]$samp_seg_name
				old.nMinor.vec <- orig.intersection.mcs.w.MSAI.df[switched.samp.seg.names.vec,]$nMinor
				#old.nBraw.vec <- orig.intersection.mcs.w.MSAI.df[switched.samp.seg.names.vec,]$nBraw
				#orig.intersection.mcs.w.MSAI.df[switched.samp.seg.names.vec,]$nBraw <-  intersection.orig.mcs.df[switched.samp.seg.names.vec,]$nAraw
				orig.intersection.mcs.w.MSAI.df[switched.samp.seg.names.vec,]$nMinor <- intersection.orig.mcs.df[switched.samp.seg.names.vec,]$nMajor
				orig.intersection.mcs.w.MSAI.df[switched.samp.seg.names.vec,]$nMajor <- old.nMinor.vec

			}#end if ( FALSE %in% AI.seg.MPhase.MSAI.ASCAT.df$balanced & TRUE %in% AI.seg.MPhase.MSAI.ASCAT.df$balanced ) {

		}#end if (all.seg.MSAI.by.BAF.ltnpd.present.list[[seg.name]]$MSAI_detection_feasible %in% TRUE) {

	}#end for (seg.name in intersection.seg.name.vec) {
	print(orig.intersection.mcs.w.MSAI.df[orig.intersection.mcs.w.MSAI.df$nMinor > orig.intersection.mcs.w.MSAI.df$nMajor,])

	orig.intersection.mcs.w.MSAI.df.path <- paste0(PLOT_OUTPUT_DIR, patient, "_orig_intersection_mcs_w_MSAI_df.RDS")
	saveRDS(orig.intersection.mcs.w.MSAI.df, file = orig.intersection.mcs.w.MSAI.df.path)


	pre.MPhase.MEDICC.FASTA.file.list <- makeMEDICCInputListsFromSegmentationDF(segs.to.make.MEDICC.input.df = orig.intersection.mcs.w.MSAI.df,
											sample.vec = sample.vec,
											chrs.to.use.vec = CHRS_TO_USE_VEC,
											chr.col = "Chromosome",
											start.col = "Start_Position",
											end.col = "End_Position",
											cpnA.col = "nMajor",
											cpnB.col = "nMinor"
											)

	pre.MPhase.MEDICC.FASTA.file.list.path <- paste0(PLOT_OUTPUT_DIR, patient, "_pre_MPhase_MEDICC_FASTA_file_list.RDS")
	saveRDS(pre.MPhase.MEDICC.FASTA.file.list, file = pre.MPhase.MEDICC.FASTA.file.list.path)


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


	#TBKW edit, so this is where we want to update the raw copy number to reflect what
	#ASCAT has outputted
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

	print("Calculating CCFs using phased raw cpn...")
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

	  }#for (seg.name in seg.name.vec

	}#end for (cpnX in c("cpnA", "cpnB")) {

	#Save CCF data:
	saveRDS(all.seg.CCF.list, file = all.seg.CCF.list.path)
	saveRDS(all.seg.int.ob.xp.cpn.df.list, file = all.seg.int.ob.xp.cpn.df.list.path)
	saveRDS(all.seg.CCF.diff.from.exp.list, file = all.seg.CCF.diff.from.exp.list.path)
	saveRDS(all.seg.clonal.in.all.samps.list, file = all.seg.clonal.in.all.samps.list.path)
	saveRDS(all.seg.greater.than.two.diff.list, file = all.seg.greater.than.two.diff.list.path)


} else if (OVERWRITE_FLAG == FALSE) {

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



#Make a table that contains all the relevant information: subclonal calls... gain and loss vs ploidy etc etc etc
if (OVERWRITE_INIT_MCS_ADD_INFO_FLAG == TRUE) {

	#Let's initialise the first joint output:
	print("Initialising minimum consistent segmentation (MCS) output dataframe...")
	empty.full.joint.segs.df <- initialiseFullOutputDataframe(seg.name.vec = names(asas.ph.raw.cpn.lists[[1]]),
															ploidy.est.vec = ploidy.est.vec,
															purity.est.vec = purity.est.vec,
															make.row.names.seg.name.alone = FALSE)


	#This adds the raw cn for both segments that could be phased and those that couldn't
	print("Adding new copy number calls to MCS output based on the phased raw copy number...")
	# new.cpn.full.joint.segs.df <- addNewRawAndIntCPNCallsToFinalSegmentation(seg.name.vec = names(asas.ph.raw.cpn.lists[[1]]),
	# 																		sample.vec = names(purity.est.vec),
	# 																		all.samp.logr.df.list = all.samp.non.hz.logr.and.tot.cpn.df.list,
	# 																		all.samp.ph.BAF.df.list = patient.phased.BAF.df.list,
	# 																		all.seg.CCF.vals.list = all.seg.CCF.list,
	# 																		asas.lr.baf.lists = asas.SNPs.w.both.lr.baf.snp.lists,
	# 																		all.seg.df = empty.full.joint.segs.df,
	# 																		cpnX.vec = c("cpnA", "cpnB")
	# 																		)

	# source(paste0(current.dir, "/final_output_segmentation_file_generation_functions_v3.R"))
	# source(paste0(current.dir, "/utility_MPhase_functions.R"))
	new.cpn.full.joint.segs.df <- addNewRawAndIntCPNCallsToFinalSegmentationParallel(seg.name.vec = names(asas.ph.raw.cpn.lists[[1]]),
																			sample.vec = names(purity.est.vec),
																			all.samp.logr.df.list = all.samp.non.hz.logr.and.tot.cpn.df.list,
																			all.samp.ph.BAF.df.list = patient.phased.BAF.df.list,
																			all.seg.CCF.vals.list = all.seg.CCF.list,
																			asas.lr.baf.lists = asas.SNPs.w.both.lr.baf.snp.lists,
																			all.seg.df = empty.full.joint.segs.df,
																			cpnX.vec = c("cpnA", "cpnB")
																			)




	# print("Mapping MCS segments back to original input and copy number states...")
	full.joint.segs.w.old.cpn.vs.ploidy.df <- addOriginalSegmentInfoToFinalSegmentationNew(new.cpn.full.joint.segs.df, sorted.all.bp.all.samps.segs.anno.df)





	print("Adding CCF information and subclonal cpn calls to MCS output...")
	clonality.and.ccf.joint.segs.df <- addWithinSampleClonalityAndCCFToFinalSegmentation(seg.name.vec = names(asas.ph.raw.cpn.lists[[1]])
	                                ,all.seg.CCF.vals.list = all.seg.CCF.list
	                                ,all.seg.CCF.diff.from.exp.list = all.seg.CCF.diff.from.exp.list
	                                ,all.seg.df = full.joint.segs.w.old.cpn.vs.ploidy.df#new.cpn.full.joint.segs.df
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
	#source(paste0(current.dir, "/final_output_segmentation_file_generation_functions_v3.R"))
	#source(paste0(current.dir, "/utility_MPhase_functions.R"))
	full.joint.segs.df.w.cpn.vs.ploidy.df <- addLogrDefinedCPNEventsToOutputSegmentationParallel(sample.vec = names(purity.est.vec),
	                                                        							seg.name.vec = names(asas.ph.raw.cpn.lists[[1]]),
	                                                        							all.seg.df = clonality.and.ccf.joint.segs.df,
	                                                        							ploidy.est.vec = ploidy.est.vec,
	                                                        							purity.est.vec = purity.est.vec,
	                                                        							all.samp.logr.df.list = all.samp.non.hz.logr.and.tot.cpn.df.list,
	                                                        							all.samp.ph.BAF.df.list = patient.phased.BAF.df.list,
	                                                       								gain.threshold = GAIN_LOGR_THRESH,
	                                                        							loss.threshold = LOSS_LOGR_THRESH,
	                                                        							amp.threshold = AMP_LOGR_THRESH,
	                                                        							deep.loss.threshold = DEEP_LOSS_LOGR_THRESH,
	                                                        							tot.cpn.vs.ploidy.pval.thresh = PVAL_THRESH_VS_LOGR_FOR_CPN_EVENT)

	print(head(full.joint.segs.df.w.cpn.vs.ploidy.df[,"cpn_event_vs_ploidy"]))
	print("Adding cpn event vs ploidy by logr measures to MCS output...")
	#source(paste0(current.dir, "/final_output_segmentation_file_generation_functions_v3.R"))
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
	# full.joint.check.amps.vs.ploidy.clonality.df <- retestSubclonalAmps(all.seg.df = full.joint.segs.logr.vs.ploidy.clonality.df,
	# 																	all.samp.logr.df.list = all.samp.non.hz.logr.and.tot.cpn.df.list,
	# 																	all.samp.ph.BAF.df.list = patient.phased.BAF.df.list,
	# 																	gain.threshold = GAIN_LOGR_THRESH,
	# 																	amp.threshold = AMP_LOGR_THRESH,
	# 																	p.val.for.converting.gain.to.amp = PVAL_THRESH_CONVERT_GAIN_TO_AMP)

	print(head(full.joint.segs.logr.vs.ploidy.clonality.df[,c("across_samp_cl_cpn_vs_ploidy_present","across_samp_sc_cpn_vs_ploidy_present")]))


	#These are kept in separate functions as the original segmentaiton was refactored so it could run when segments are rejoined at the end of the script.
	print("Adding phasing information and MSAI detection by raw cpn calls to MCS output...")
	full.joint.segs.w.phasing.df <- addPhasingInformationToOutputSegmentation(sample.vec = names(purity.est.vec),
	                                        seg.name.vec = names(asas.ph.raw.cpn.lists[[1]]),
	                                        sample.to.phase.seg.list = all.seg.phasing.info.list,
	                                        all.pw.cn.het.comp.lists = NULL,
	                                        all.segs.df = full.joint.segs.logr.vs.ploidy.clonality.df)


	# #Now we want to add in the allelic imbalance clonal and subclonal inrormation to the final data frame.
	# print("Adding across sample clonal/subclonal summary to MCS output...")
	# full.joint.segs.w.phasing.and.raw.df <- addOriginalCPNValuesAndCallsVsPloidyToOutputSegmentation(sample.vec = sample.vec,
	#                                         seg.name.vec = names(asas.ph.raw.cpn.lists[[1]]),
	#                                         processed.orig.cpn.segs.df = processed.orig.cpn.segs.df,
	#                                         all.segs.df = full.joint.segs.w.phasing.df)

	print("Adding MSAI detected by LTNPD method to the MCS output...")
	full.joint.segs.w.ltnpd.MSAI.df <- addMSAIDetectionByLTNPDToOutputSegmentation(sample.vec = sample.vec,
												seg.name.vec = names(asas.ph.raw.cpn.lists[[1]]),
												all.seg.MSAI.by.BAF.ltnpd.present.list = all.seg.MSAI.by.BAF.ltnpd.present.list,
												all.segs.df = full.joint.segs.w.phasing.df)#full.joint.segs.w.phasing.and.raw.df)


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
	full.joint.segs.within.and.changed.cpn.info.df <- addWithinSampleChangeToCpnAndAllelicImbalanceCalls(all.seg.df = full.joint.segs.w.parallel.cpn.events.df,
																										sample.vec = sample.vec,
																										seg.name.vec = unique(seg.name.vec))



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

	#Legacy mostly, trying to make definitions that Emilia was doing and put them ion here to ensure they were correct
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



} else if (OVERWRITE_INIT_MCS_ADD_INFO_FLAG == FALSE) {

	full.joint.segs.df.path <- paste0(PLOT_OUTPUT_DIR, patient, "_full_joint_segs_df.RDS")
	full.joint.segs.df <- readRDS(file = full.joint.segs.df.path)

}#end  else if (OVERWRITE_INIT_MCS_ADD_INFO_FLAG == FALSE) {

if (PLOTTING_FLAG == TRUE) {

	####Plot all genome wide plots
	print("Plotting CCFs for cpnA and cpnB across the genome using cpn definition of subclonality...")
	png(file = paste(PLOT_OUTPUT_DIR, patient, "_CCF_w_cpn_thresh_cpnA", ".png", sep = ""), width = 1200, height = 1400)

	plotAllSampleCCFsAcrossGenome(asas.ph.raw.cpn.lists = asas.ph.raw.cpn.lists,
	                              full.joint.segs.df = full.joint.segs.df,
	                              sample.vec = sample.vec,
	                              seg.name.vec = seg.name.vec,
	                              all.seg.CCF.vals.list = all.seg.CCF.list,
	                              cpnX = "cpnA",
	                              use.cpn.or.ccf.thresh = "CPN")

	dev.off()

	png(file = paste(PLOT_OUTPUT_DIR, patient, "_CCF_w_cpn_thresh_cpnB", ".png", sep = ""), width = 1200, height = 1400)

	plotAllSampleCCFsAcrossGenome(asas.ph.raw.cpn.lists = asas.ph.raw.cpn.lists,
	                              full.joint.segs.df = full.joint.segs.df,
	                              sample.vec = sample.vec,
	                              seg.name.vec = seg.name.vec,
	                              all.seg.CCF.vals.list = all.seg.CCF.list,
	                              cpnX = "cpnB",
	                              use.cpn.or.ccf.thresh = "CPN")
	dev.off()

	################################################

} else if (PLOTTING_FLAG == FALSE) {


	print("Not plotting CCFs for cpnA and cpnB across the genome using cpn definition of subclonality as plotting flag is FALSE.")


}

# We need a .tsv version of this data for interoperability with python
write.table(full.joint.segs.df, file=file.path(PLOT_OUTPUT_DIR, "full_joint_segs.tsv"), row.names=FALSE, sep="\t")
cat("MPhase comparison with snpphase: Finished for now. Got the <patient>_full_joint_segs_df.RDS object.\n")


# Hack! Hope this works...
full.joint.segs.nona <- full.joint.segs.df[!is.na(full.joint.segs.df$ph_int_cpnA),]
# post.MPhase.MCS.MEDICC.FASTA.files.list <- makeMEDICCInputListsFromSegmentationDF(segs.to.make.MEDICC.input.df = pj.fjs.df[ ( ! is.na(pj.fjs.df$n.het)) & pj.fjs.df$n.het >= 5 ,],
#post.MPhase.MCS.MEDICC.FASTA.files.list <- makeMEDICCInputListsFromSegmentationDF(segs.to.make.MEDICC.input.df = all.samp.mcs.joined.cpnA.cpnB.df,
post.MPhase.MCS.MEDICC.FASTA.files.list <- makeMEDICCInputListsFromSegmentationDF(segs.to.make.MEDICC.input.df = full.joint.segs.nona,
										sample.vec = sample.vec,
										chrs.to.use.vec = CHRS_TO_USE_VEC,
										chr.col = "chrom",
										start.col = "startpos",
										end.col = "endpos",
										cpnA.col = "ph_int_cpnA",
										cpnB.col = "ph_int_cpnB"
										)
post.MPhase.MCS.MEDICC.FASTA.files.list.path <- paste0(PLOT_OUTPUT_DIR, patient, "_post_MPhase_MCS_MEDICC_FASTA_file_list.RDS")
#post.MPhase.MEDICC.FASTA.files.list
saveRDS(post.MPhase.MCS.MEDICC.FASTA.files.list, file = post.MPhase.MCS.MEDICC.FASTA.files.list.path)

createMEDICCInputFiles(post.MPhase.MCS.MEDICC.FASTA.files.list, sample.vec, conf$medicc_output_dir)


#quit(save = "no", status = 0, runLast = FALSE)


#So there are various things that can be improved.

# ###Additional calling that needs adding:
# #wgii
# #ploidy
# #fraction LOH
#
# #focal and arm level event definition
#
#
#
#
#
#
# if (OVERWRITE_ARM_INFO_FOR_JOIN_INIT_MCS_FLAG == TRUE) {
#
#
# 	chr.arms.present.vec <- unique(full.joint.segs.df$chr_arm)
#
#
# 	print("Make data frame of raw data behind each chromosome arm in which we will join the segments...")
# 	#Make a data onject containing info needed to join segments together
# 	arm.iterable.df <- expand.grid(sample.vec, chr.arms.present.vec, stringsAsFactors=FALSE)
# 	colnames(arm.iterable.df) <- c("sample", "arm")
# 	arm.iterable.df <- arm.iterable.df[with(arm.iterable.df, order(sample, arm)), ]
# 	samp.arm.itx <- iter(arm.iterable.df, by="row")
#
# 	samp.arms <- foreach(samp.arm.df = samp.arm.itx, .combine = c) %dopar% {
#
# 	# samp.arms <- list()
# 	# for (samp in sample.vec) {
# 	# 	print(samp)
# 	# 	for (chr.arm in chr.arms.present.vec) {
#
# 			samp <- samp.arm.df$sample
# 			chr.arm <- samp.arm.df$arm
# 			#print(chr.arm)
# 			#So now we are within a chromosome arm we want to make a seg.name.vec that matches only the chromosome
# 			samp.chr.arm.segs.df <- full.joint.segs.df[full.joint.segs.df$chr_arm %in% chr.arm
# 							                         & full.joint.segs.df$sample %in% samp, ]
# 			samp.chr.arm.seg.name.vec <- samp.chr.arm.segs.df$seg_name
#
# 			joinable.df <- 	data.frame(chrom=numeric(0),
# 						pos=numeric(0),
# 						BAF=numeric(0),
# 						BAFsegmented=numeric(0),
# 						logR=numeric(0),
# 						logRsegmented=numeric(0),
# 						sample = character(0),
# 						orig_seg_name=character(0),
# 						chr_pos = character(0),
# 						tot_raw_un_ph_cpn = character(0),
# 						ph_cpnA  = numeric(0),
# 						ph_cpnB = numeric(0),
# 						stringsAsFactors=FALSE)
# 			#Note that if there are no SNPs within a segment then it won't be recorded in this dataframe.
# 			for (seg.name in samp.chr.arm.seg.name.vec) {
# 				#print(seg.name)
# 				#Add the logr values
# 				orig.seg.lr.and.or.baf.df <- asas.SNPS.lr.and.or.baf.lists[[samp]][[seg.name]]
# 				if(nrow(orig.seg.lr.and.or.baf.df) > 0) {
# 					#print("Adding original logr and baf...")
# 					#print(nrow(orig.seg.lr.and.or.baf.df))
# 					orig.seg.lr.and.or.baf.df$sample <- samp
# 					orig.seg.lr.and.or.baf.df$orig_seg_name <- seg.name
# 					orig.seg.lr.and.or.baf.df$chr_pos <- paste0(orig.seg.lr.and.or.baf.df$chrom, "_", orig.seg.lr.and.or.baf.df$pos)
# 					orig.seg.lr.and.or.baf.df$tot_raw_un_ph_cpn <- NA
# 					orig.seg.lr.and.or.baf.df$ph_cpnA <- NA
# 					orig.seg.lr.and.or.baf.df$ph_cpnB <- NA
#
# 					#Now let's get the cpn values
# 					#Let's add the total copy number
# 					#print("Adding total copy number...")
# 					#print(nrow(orig.seg.lr.and.or.baf.df))
# 					samp.logr.tot.cpn.df <- all.samp.non.hz.logr.and.tot.cpn.df.list[[samp]]
# 					samp.logr.tot.cpn.df <- samp.logr.tot.cpn.df[ ! is.na(samp.logr.tot.cpn.df$chrom) , , drop = FALSE ]
#
# 					if(nrow(samp.logr.tot.cpn.df) > 0) {
#
# 						samp.logr.tot.cpn.df$chr_pos <- paste0(samp.logr.tot.cpn.df$chrom, "_", samp.logr.tot.cpn.df$pos)
# 						orig.seg.lr.and.or.baf.df$tot_raw_un_ph_cpn <- samp.logr.tot.cpn.df[match(orig.seg.lr.and.or.baf.df$chr_pos, samp.logr.tot.cpn.df$chr_pos),]$tot_raw_un_ph_cpn
#
# 					}#end if(nrow(samp.logr.tot.cpn.df) > 0) {
#
# 					#Let's add phased copy number
# 					seg.ph.cpn.list <- asas.ph.raw.cpn.lists[[samp]][[seg.name]]
# 					if (  ! anyNA(seg.ph.cpn.list[1]$chr_vec) ) {
# 						#print("Adding phased copy number...")
# 						#print(nrow(orig.seg.lr.and.or.baf.df))
# 						#make a dataframe:
# 						samp.seg.ph.cpn.df <- data.frame(chrom = seg.ph.cpn.list$chr_vec,
# 									pos = seg.ph.cpn.list$pos_vec,
# 									ph_cpnA = seg.ph.cpn.list$ph_cpnA_vec,
# 									ph_cpnB = seg.ph.cpn.list$ph_cpnB_vec,
# 									chr_pos = paste0(seg.ph.cpn.list$chr_vec, "_", seg.ph.cpn.list$pos_vec),
# 									stringsAsFactors = FALSE)
# 						matching.vec <- match(samp.seg.ph.cpn.df$chr_pos, orig.seg.lr.and.or.baf.df$chr_pos)
# 						matching.vec <- matching.vec[ ! is.na(matching.vec) ]
# 						orig.seg.lr.and.or.baf.df[matching.vec, c("ph_cpnA", "ph_cpnB")] <- samp.seg.ph.cpn.df[, c("ph_cpnA", "ph_cpnB")]
#
# 					}#end if ( ! is.na(seg.ph.cpn.list$chr_vec) ) {
# 					#print(nrow(orig.seg.lr.and.or.baf.df))
# 					orig.seg.lr.and.or.baf.df <- orig.seg.lr.and.or.baf.df[! is.na(orig.seg.lr.and.or.baf.df$chrom), , drop = FALSE]
# 					#print(nrow(orig.seg.lr.and.or.baf.df))
#
# 					#Now let's make sure that all the phased copy number has total cpn values
# 					#Cases with phased copy number
# 					# if(nrow(orig.seg.lr.and.or.baf.df) > 0) {
#
# 					# 	STOP
# 					# }
#
#
# 					#So let's calculate missing total copy number where this phased copy number:
# 					if(nrow(orig.seg.lr.and.or.baf.df[is.na(orig.seg.lr.and.or.baf.df$tot_raw_un_ph_cpn)
# 											& (! is.na(orig.seg.lr.and.or.baf.df$ph_cpnA))
# 											& (! is.na(orig.seg.lr.and.or.baf.df$ph_cpnB)), c('ph_cpnA', 'ph_cpnB')]) > 0){
#
# 						orig.seg.lr.and.or.baf.df[is.na(orig.seg.lr.and.or.baf.df$tot_raw_un_ph_cpn)
# 											& (! is.na(orig.seg.lr.and.or.baf.df$ph_cpnA))
# 											& (! is.na(orig.seg.lr.and.or.baf.df$ph_cpnB)), ]$tot_raw_un_ph_cpn <- rowSums(orig.seg.lr.and.or.baf.df[is.na(orig.seg.lr.and.or.baf.df$tot_raw_un_ph_cpn)
# 																																				& (! is.na(orig.seg.lr.and.or.baf.df$ph_cpnA))
# 																																				& (! is.na(orig.seg.lr.and.or.baf.df$ph_cpnB)), c('ph_cpnA', 'ph_cpnB')])
#
# 					}#end if(nrow(orig.seg.lr.and.or.baf.df[is.na(orig.seg.lr.and.or.baf.df$tot_raw_un_ph_cpn) ...
#
# 					#Now let's calculate total copy number if it's missing and there's no phased copy number:
# 					#For this we'll have to call the
# 					if(nrow(orig.seg.lr.and.or.baf.df[is.na(orig.seg.lr.and.or.baf.df$tot_raw_un_ph_cpn)
# 											&  is.na(orig.seg.lr.and.or.baf.df$ph_cpnA)
# 											&  is.na(orig.seg.lr.and.or.baf.df$ph_cpnB),]) > 0) {
#
# 						missing.tot.cpn.logr.vec <- orig.seg.lr.and.or.baf.df[is.na(orig.seg.lr.and.or.baf.df$tot_raw_un_ph_cpn)
# 											& (is.na(orig.seg.lr.and.or.baf.df$ph_cpnA))
# 											& (is.na(orig.seg.lr.and.or.baf.df$ph_cpnB)),]$logR
#
# 						orig.seg.lr.and.or.baf.df[is.na(orig.seg.lr.and.or.baf.df$tot_raw_un_ph_cpn)
# 											& (is.na(orig.seg.lr.and.or.baf.df$ph_cpnA))
# 											& (is.na(orig.seg.lr.and.or.baf.df$ph_cpnB)), ]$tot_raw_un_ph_cpn <- rowSums(calculateSNPsRawCopyNumberAB(logr.vec = missing.tot.cpn.logr.vec,
# 																																			  baf.vec = rep(0.5, length(missing.tot.cpn.logr.vec)),
# 																																			  tumor.ploidy = ploidy.est.vec[samp],
# 																																			  tumor.purity = purity.est.vec[samp],
# 																																			  gamma.val = GAMMA_VAL))
#
# 					}#end if(nrow(orig.seg.lr.and.or.baf.df[is.na(orig.seg.lr.and.or.baf.df$tot_raw_un_ph_cpn) ...
#
# 					joinable.df <- rbind(joinable.df, orig.seg.lr.and.or.baf.df, stringsAsFactors = FALSE)
#
# 				}#end if(nrow(orig.seg.lr.and.or.baf.df) > 0) {
#
# 			}#end for (seg.name in samp.chr.arm.seg.name.vec) {
# 				print(paste0(samp, "_", chr.arm))
#
#
# 				#samp.arms[[samp]][[chr.arm]] <- joinable.df
#
# 				returnable.list <- list(joinable.df)
# 				names(returnable.list) <- paste0(samp, "_", chr.arm)
# 				returnable.list
#
# 	}
# 	# 	}#end for (chr.arm in chr.arms.present.vec) {
#
# 	# }#for (samp in sample.vec) {
#
#
# 	#Quickly convert the results:
# 	parallel.samp.arms <- samp.arms
#
#
#
# 	samp.arms <- list()
# 	for(samp in sample.vec) {
#
# 		samp.arms[[samp]] <- list()
#
# 		for(chr.arm in chr.arms.present.vec) {
#
#
# 			samp.arms[[samp]][[chr.arm]] <- parallel.samp.arms[[paste0(samp, "_", chr.arm)]]
#
# 		}
#
# 	}
#
#
# 	all.sample.arms.underlying.raw.and.phased.data.path <- paste0(PLOT_OUTPUT_DIR, patient, "_all_sample_arm_underlying_raw_and_phased_data.RDS")
# 	saveRDS(samp.arms, file = all.sample.arms.underlying.raw.and.phased.data.path )
#
#
# } else if (OVERWRITE_ARM_INFO_FOR_JOIN_INIT_MCS_FLAG == FALSE) {
#
# 	all.sample.arms.underlying.raw.and.phased.data.path <- paste0(PLOT_OUTPUT_DIR, patient, "_all_sample_arm_underlying_raw_and_phased_data.RDS")
# 	samp.arms <- readRDS(file = all.sample.arms.underlying.raw.and.phased.data.path )
#
# }
#
# ########################################################################
# # Join segments
# ########################################################################
#
# if (OVERWRITE_JOIN_INIT_MCS_FLAG == TRUE) {
#
# 	MIN_PH_CPN <- 5
# 	NUM_TOT_CPN <- 5
# 	CPN_AVG_THRESH <- 0.75
# 	cpnX.vec <- c("ph_cpnA", "ph_cpnB")
# 	chr.arms.present.vec <- unique(full.joint.segs.df$chr_arm)
#
# 	#Debugging values:
# 	#full_sample <- sapply(rownames(all.bp.all.samps.segs.df), function(x) strsplit(x, split = ":")[[1]][1])
# 	# cpnX.string <- "ph_cpnA"
# 	# chr.arm <- "3p"
#
# 	print("Join segments within each sample's arms if they meet significance criteria....")
# 	#This needs to be parallelised... as it takes a while.
# 	#Let's try and debug it a bit more first, run it on the tx100.
#
# 	JoinedArmMultiResultClass <- function(samp_chr_arm = NULL, joined_samp_arms=NULL, joined_samp_arms_inferred_ph_cpn=NULL)
# 	{
# 	  me <- list(
# 	  	samp_chr_arm = samp_chr_arm,
# 	  	joined_samp_arms = joined_samp_arms,
# 	    joined_samp_arms_inferred_ph_cpn = joined_samp_arms_inferred_ph_cpn
# 	  )
# 	  # Set the name for the class
# 	  class(me) <- append(class(me),"MultiResultClass")
# 	  return(me)
# 	}
#
# 	#source("/home/camp/watkint/working/watkint/scripts/multiregioncopynum/new_seg_joining_functions.R")
# 	source(paste0(current.dir, "new_seg_joining_functions.R"))
#
# 	joined.samp.arm.iterable.df <- expand.grid(sample=sample.vec, arm=chr.arms.present.vec, stringsAsFactors=FALSE)
# 	joined.samp.arm.iterable.df <- joined.samp.arm.iterable.df[with(joined.samp.arm.iterable.df, order(sample, arm)), ]
# 	joined.samp.arm.itx <- iter(joined.samp.arm.iterable.df, by="row")
#   #cat("Savign object and starting foreach: joined.samp.arms.list\n")
#   #saveRDS(joined.samp.arm.iterable.df, file = "joined_samp_arm_iterable_df.rds")
#   joined.samp.arms.list <- list()
#
#   for (i in seq_len(nrow(joined.samp.arm.iterable.df))) {
#     joined.samp.arm.df <- joined.samp.arm.iterable.df[i,]
#     #print(joined.samp.arm.df)
# 	#joined.samp.arms.list <- foreach(joined.samp.arm.df = joined.samp.arm.itx, .export=ls(envir=globalenv()), .combine = c) %dopar% {
#
#       require("GenomicRanges")
#       new.sample.vec <- sample.vec#[length(sample.vec)]
# 	#new.chr.arms.present.vec <- chr.arms.present.vec[34:34]
# 	# joined.samp.arms <- list()
# 	# joined.samp.arms.inferred.ph.cpn <- list()
# 	# for( samp in sample.vec) {
#
# 	# 	joined.samp.arms[[samp]] <- list()
# 	# 	joined.samp.arms.inferred.ph.cpn[[samp]] <- list()
#
# 	# 	for ( chr.arm in chr.arms.present.vec) {
# 			#print("here")
# 			#Parallel version:
# 			samp <- joined.samp.arm.df$sample
# 			chr.arm <- joined.samp.arm.df$arm
# 			#End parallel vesion
# 			print(paste0(samp, "_", chr.arm))
# 			chr.arm.joined.list <- list()
#
# 			for (cpnX.string in c("ph_cpnA", "ph_cpnB")) {
# 				# print(cpnX.string)
# 				# print(samp)
# 				# print(chr.arm)
# 				samp.arm.sites.df <- samp.arms[[samp]][[chr.arm]]
# 				samp.chr.arm.segs.df <- full.joint.segs.df[full.joint.segs.df$chr_arm %in% chr.arm
# 								                         & full.joint.segs.df$sample %in% samp, ]
# 				samp.chr.arm.seg.name.vec <- samp.chr.arm.segs.df$seg_name
#
# 				#samp.arm.sites.df <- joinable.df
#
# 				chr.vec <- as.numeric(sapply(samp.chr.arm.seg.name.vec, function(x) strsplit(x, split = ":")[[1]][1]))
# 				start.vec <- as.numeric(sapply(samp.chr.arm.seg.name.vec, function(x) strsplit(x, split = ":")[[1]][2]))
# 				end.vec <- as.numeric(sapply(samp.chr.arm.seg.name.vec, function(x) strsplit(x, split = ":")[[1]][3]))
# 				arm.grs <- GRanges(seqnames=chr.vec,
# 								ranges=IRanges(start.vec, end.vec),
# 								strand=rep("*", length(chr.vec)))
#
# 				if(length(getOverlapGRanges(arm.grs, arm.grs)) != length(arm.grs)) {
#
# 					print(getOverlapGRanges(arm.grs, arm.grs))
# 					stop("Error, overlapping granges within an arm!")
#
# 				}
# 				current.arm.grs <- arm.grs
# 				orig.joining.vals.mat.and.arm.grs.list <- creatingArmJoiningInfoMatrixAndGranges(current.arm.grs, samp.arm.sites.df, cpnX.string)
# 				orig.joining.vals.mat <- orig.joining.vals.mat.and.arm.grs.list[[1]]
# 				orig.joining.vals.mat <- orig.joining.vals.mat[order(orig.joining.vals.mat[,"startpos"]),]
# 				orig.current.arm.grs <- orig.joining.vals.mat.and.arm.grs.list[[2]]
#
# 				old.length.current.arm.grs <- 0
# 				while (length(current.arm.grs) != old.length.current.arm.grs & length(current.arm.grs) != 1) {
#
# 					joining.vals.mat.and.arm.grs.list <- creatingArmJoiningInfoMatrixAndGranges(current.arm.grs, samp.arm.sites.df, cpnX.string)
# 					joining.vals.mat <- joining.vals.mat.and.arm.grs.list[[1]]
# 					current.arm.grs <- joining.vals.mat.and.arm.grs.list[[2]]
# 					# print("############")
# 					# print(current.arm.grs)
# 					# print(joining.vals.mat)
# 					# print("%%%%%%%%%%%%%")
# 					# readline(prompt="Press [enter] to continue")
# 					old.length.current.arm.grs <- length(current.arm.grs)
# 					current.arm.grs <- attemptToJoinOneSegmentWithinChromArm(current.arm.grs, joining.vals.mat, samp.arm.sites.df, sorted.arm.grs, cpnX.string, MIN_PH_CPN, NUM_TOT_CPN, CPN_AVG_THRESH)
#
# 					if(length(getOverlapGRanges(current.arm.grs, current.arm.grs)) != length(current.arm.grs)) {
# 						print(getOverlapGRanges(current.arm.grs, current.arm.grs))
# 						stop("Error, overlapping granges within an arm!")
# 					}
# 					# print(length(current.arm.grs))
# 					# print(old.length.current.arm.grs)
# 				} #end while (length(current.arm.grs) == old.length.current.arm.grs ) {
#
# 				final.joining.vals.mat.and.arm.grs.list <- creatingArmJoiningInfoMatrixAndGranges(current.arm.grs, samp.arm.sites.df, cpnX.string)
# 				final.joining.vals.mat <- final.joining.vals.mat.and.arm.grs.list[[1]]
# 				final.joining.vals.mat <- final.joining.vals.mat[order(final.joining.vals.mat[,"startpos"]), , drop = FALSE]
# 				final.current.arm.grs <- final.joining.vals.mat.and.arm.grs.list[[2]]
#
# 				sorted.joining.vals.mat <- joining.vals.mat[order(joining.vals.mat[,"startpos"]), ]
# 				##Serial version. commented out:
# 				#joined.samp.arms[[samp]][[chr.arm]][[cpnX.string]] <- as.data.frame(final.joining.vals.mat, stringsAsFactors=FALSE)
# 				###End Serial version.
# 				chr.arm.joined.list[[cpnX.string]] <- as.data.frame(final.joining.vals.mat, stringsAsFactors=FALSE)
#
# 			} # end for (cpnX.string in cpnX.vec) {
# 			#As a first pass
# 			chr.arm.joined.w.int.cpn.list <- list()
# 			for (cpnX.string in c("ph_cpnA", "ph_cpnB")) {
#
# 				cpnX.chr.joined.df <- chr.arm.joined.list[[cpnX.string]]
# 				cpnX.chr.joined.df[, paste0("int_",cpnX.string)] <- -1
# 				#So now for each row calculate the total and phased copy number
# 				for (k in 1:nrow(cpnX.chr.joined.df)) {
# 					#print(cpnX.chr.joined.df[k,])
# 					if (cpnX.chr.joined.df[k,"num_ph_cpn"] > MIN_PH_CPN) {
# 						#Let's get the relevant copy number values from the samp.sites.df
# 						cpnX.chr.joined.df[k, paste0("int_",cpnX.string)] <- ifelse(cpnX.chr.joined.df[k, ]$median_ph_cpn < 0,
# 																					0,
# 																					round(cpnX.chr.joined.df[k, ]$median_ph_cpn) )
#
# 					}
#
# 				}
# 				chr.arm.joined.w.int.cpn.list[[cpnX.string]] <- cpnX.chr.joined.df
#
# 			}
# 			#So now we want to decide which segment we should use to fill all the any blank row.
# 			#The segment will have failed to have been joined in both cpnA and cpnB dfs as it can only joined by tot cpn
# 			#However the segments around it may have differet total copy numbers...
#
# 			#So actually the function call starts up here:
#
# 			cpnA.joined.segs.df <- chr.arm.joined.w.int.cpn.list[[cpnX.vec[1]]]
# 			cpnB.joined.segs.df <- chr.arm.joined.w.int.cpn.list[[cpnX.vec[2]]]
# 			#Will dataframe to return
# 			cpnA.no.blank.df <- cpnA.joined.segs.df
# 			cpnB.no.blank.df <- cpnB.joined.segs.df
#
#       #write.csv(head(cpnA.joined.segs.df), file="head-cpnA.joined.segs.df.csv")
#       #print(head(cpnA.joined.segs.df))
# 			cpnA.joined.segs.grs <- makeGRangesFromDataFrame(cpnA.joined.segs.df,
# 					                         		keep.extra.columns=TRUE,
# 					                         		ignore.strand=TRUE,
# 					                         		seqinfo=NULL,
# 					                         		seqnames.field="chrom",
# 					                         		start.field="startpos",
# 					                         		end.field="endpos",
# 					                         		starts.in.df.are.0based=FALSE)
# 			#Checking for -1 as some of segments that have >= MIN_PH_CPN still lack ph cpn estimates, need to check this further. TODO
# 			only.ph.cpnA.joined.segs.grs <- cpnA.joined.segs.grs[cpnA.joined.segs.grs$num_ph_cpn >= MIN_PH_CPN & cpnA.joined.segs.grs$int_ph_cpnA != -1 ]
#       #write.csv(head(cpnB.joined.segs.df), file="head-cpnB.joined.segs.df.csv")
#       #print(head(cpnB.joined.segs.df))
# 			cpnB.joined.segs.grs <- makeGRangesFromDataFrame(cpnB.joined.segs.df,
# 					                         		keep.extra.columns=TRUE,
# 					                         		ignore.strand=TRUE,
# 					                         		seqinfo=NULL,
# 					                         		seqnames.field="chrom",
# 					                         		start.field="startpos",
# 					                         		end.field="endpos",
# 					                         		starts.in.df.are.0based=FALSE)
# 			#Checking for -1 as some of segments that have >= MIN_PH_CPN still lack ph cpn estimates, need to check this further. TODO
# 			only.ph.cpnB.joined.segs.grs <- cpnB.joined.segs.grs[cpnB.joined.segs.grs$num_ph_cpn >= MIN_PH_CPN & cpnB.joined.segs.grs$int_ph_cpnB != -1 ]
# 			cpnX.joined.segs.df.list <- list(cpnA.joined.segs.df, cpnB.joined.segs.df)
# 			cpnX.no.blank.df.list <- list(cpnA.no.blank.df, cpnB.no.blank.df)
#
#
# 			#So this is really a function call:
#
# 			for (j in 1:length(cpnX.joined.segs.df.list)) {
#
# 				cpnX.joined.segs.df <- cpnX.joined.segs.df.list[[j]]
#
# 				for (k in 1:nrow(cpnX.joined.segs.df)) {
#
# 					#print(k)
#
# 					if (cpnX.joined.segs.df[k, "num_ph_cpn"] < MIN_PH_CPN | cpnX.joined.segs.df[k, ncol(cpnX.joined.segs.df)] == -1 ) {
# 						#Let's get the relevant copy number values from the samp.sites.df
# 						#So now we want to find out which of the neighboring segments it's closest to in both dfs.
# 						#Make Query Grange:
#             #write.csv(head(cpnX.joined.segs.df[k,]), file="head-cpnX.joined.segs.df.csv")
#             #print(head(cpnX.joined.segs.df[k,]))
#
# 						query.seg.gr <- makeGRangesFromDataFrame(cpnX.joined.segs.df[k,],
# 						                         		keep.extra.columns=TRUE,
# 						                         		ignore.strand=TRUE,
# 						                         		seqinfo=NULL,
# 						                         		seqnames.field="chrom",
# 						                         		start.field="startpos",
# 						                         		end.field="endpos",
# 						                         		starts.in.df.are.0based=FALSE)
# 						#Note that there is the case that the segment is actaully between two others that don't have phased copy number
# 						#Therefore we should also include the nearest segment with phaseable cpn...
# 						#Or we can just subset the data.frame to remove all other segs that don't have phased copy number?
# 						cpnA.precede.and.follow.gr.idx.vec <- c(follow(query.seg.gr, only.ph.cpnA.joined.segs.grs), precede(query.seg.gr, only.ph.cpnA.joined.segs.grs))
# 						cpnA.sort.non.na.gr.pre.and.fol.idx.vec <- sort(cpnA.precede.and.follow.gr.idx.vec[ ! is.na(cpnA.precede.and.follow.gr.idx.vec)])
# 						#edge case where all segments are joined into a single segment so there are no preceding or following... get the nearest:
# 						if (length(cpnA.sort.non.na.gr.pre.and.fol.idx.vec) == 0){
#
# 							cpnA.sort.non.na.gr.pre.and.fol.idx.vec <- nearest(query.seg.gr, only.ph.cpnA.joined.segs.grs)
# 						}
# 						adj.ph.cpnA.grs <- only.ph.cpnA.joined.segs.grs[cpnA.sort.non.na.gr.pre.and.fol.idx.vec]
# 						cpnA.tot.cpn.diff.vec <- abs(adj.ph.cpnA.grs$median_tot_cpn - query.seg.gr$median_tot_cpn)
# 						cpnA.val <- adj.ph.cpnA.grs[which(cpnA.tot.cpn.diff.vec == min(cpnA.tot.cpn.diff.vec))]$int_ph_cpnA
# 						cpnB.precede.and.follow.gr.idx.vec <- c(follow(query.seg.gr, only.ph.cpnB.joined.segs.grs), precede(query.seg.gr, only.ph.cpnB.joined.segs.grs))
# 						cpnB.sort.non.na.gr.pre.and.fol.idx.vec <- sort(cpnB.precede.and.follow.gr.idx.vec[ ! is.na(cpnB.precede.and.follow.gr.idx.vec)])
# 						#edge case where all segments are joined into a single segment so there are no preceding or following... get the nearest:
# 						if (length(cpnB.sort.non.na.gr.pre.and.fol.idx.vec) == 0){
#
# 							cpnB.sort.non.na.gr.pre.and.fol.idx.vec <- nearest(query.seg.gr, only.ph.cpnB.joined.segs.grs)
#
# 						}
# 						adj.ph.cpnB.grs <- only.ph.cpnB.joined.segs.grs[cpnB.sort.non.na.gr.pre.and.fol.idx.vec]
# 						cpnB.tot.cpn.diff.vec <- abs(adj.ph.cpnB.grs$median_tot_cpn - query.seg.gr$median_tot_cpn)
# 						cpnB.val <- adj.ph.cpnB.grs[which(cpnB.tot.cpn.diff.vec == min(cpnB.tot.cpn.diff.vec))]$int_ph_cpnB
# 						#now we want the total rounded to
#             cat("i=", i, "\tj=", j, "\tk=", k, "\tcpnA.val=", cpnA.val, "\tcpnB.val=", cpnB.val, "\n")
# 						query.seg.int.tot.cpn <- round(query.seg.gr$median_tot_cpn)
# 						if (query.seg.int.tot.cpn < 0) {
#
# 							query.seg.int.tot.cpn <- 0
#
# 						}
# 						if(cpnA.val == 0 & cpnB.val == 0) {
#
# 							final.cpnA.val <- 0
# 							final.cpnB.val <- 0
#
# 						} else {
#
# 							top.ratio.val <- max(cpnA.val, cpnB.val) / sum(cpnA.val, cpnB.val)
# 							min.ratio.val <- min(cpnA.val, cpnB.val) / sum(cpnA.val, cpnB.val)
# 							if (cpnA.val >= cpnB.val ) {
#
# 								final.cpnA.val <- round(top.ratio.val * query.seg.int.tot.cpn)
# 								final.cpnB.val <- round(min.ratio.val * query.seg.int.tot.cpn)
#
# 							} else if (cpnA.val < cpnB.val) {
#
# 								final.cpnA.val <- round(min.ratio.val * query.seg.int.tot.cpn)
# 								final.cpnB.val <- round(top.ratio.val * query.seg.int.tot.cpn)
# 							}
#
# 						}
# 						for (x in 1:length(cpnX.no.blank.df.list)) {
# 						#So now we need to know which one is larger:
#
# 							if( nrow(cpnX.no.blank.df.list[[x]][cpnX.no.blank.df.list[[x]]$startpos == start(query.seg.gr) &
# 											 cpnX.no.blank.df.list[[x]]$endpos == end(query.seg.gr) ,]) > 0){
#
# 								if(x == 1) {
#
# 									final.int.ph.cpn.val <- final.cpnA.val
# 									int.ph.cpn.string <- "int_ph_cpnA"
#
# 								} else {
#
# 									final.int.ph.cpn.val <- final.cpnB.val
# 									int.ph.cpn.string <- "int_ph_cpnB"
#
# 								}
#
# 								cpnX.no.blank.df.list[[x]][cpnX.no.blank.df.list[[x]]$startpos == start(query.seg.gr) &
# 												 cpnX.no.blank.df.list[[x]]$endpos == end(query.seg.gr), int.ph.cpn.string] <- final.int.ph.cpn.val
# 												 #print("!!!!!!!!@@@@@@@@@@@@!!!!!!!!!!!@@@@@@@@@@@@@@@@@@@@")
#
# 							}#end if( nrow(cpnX.no.blank.df.list[[x]][cpnX.no.blank.df.list[[x]]$startpos == start(query.seg.gr) & ....
#
# 						}#end for (x in 1:length(cpnX.no.blank.df.list)) {
#
# 					} #end if (cpnA.joined.segs.df[k, "num_ph_cpn"] < MIN_PH_CPN) {
#
# 				}#end for (k in 1:nrow(cpnA.joined.segs.df)) {
#
# 			}#end for (j in 1:length(cpnX.joined.segs.df.list)) {
#
# 			if (any(cpnX.no.blank.df.list[[1]] %in% -1)) {
#
# 				stop("cpnX.no.blank.df.list[[1]] == -1")
#
#
# 			}
# 			if (any(cpnX.no.blank.df.list[[2]] %in% -1)) {
#
# 				stop("cpnX.no.blank.df.list[[2]] == -1")
#
#
# 			}
#
# 			#Non-parallel version:
# 			# joined.samp.arms.inferred.ph.cpn[[samp]][[chr.arm]][["ph_cpnA"]] <- cpnX.no.blank.df.list[[1]]
# 			# joined.samp.arms.inferred.ph.cpn[[samp]][[chr.arm]][["ph_cpnB"]] <- cpnX.no.blank.df.list[[2]]
#
# 	# 	#The two braces below this are for non parallel version
# 	# 	}#end for ( chr.arm in chr.arms.present.vec) {
#
# 	# } # end for( samp in sample.vec) {
#
# 		#Parallel workversion
# 		joined.samp.arms.ret.mrc.inferred.ph.cpn <- list()
# 		joined.samp.arms.ret.mrc.inferred.ph.cpn[[chr.arm]] <- list()
# 		joined.samp.arms.ret.mrc.inferred.ph.cpn[[chr.arm]][["ph_cpnA"]] <- cpnX.no.blank.df.list[[1]]
# 		joined.samp.arms.ret.mrc.inferred.ph.cpn[[chr.arm]][["ph_cpnB"]] <- cpnX.no.blank.df.list[[2]]
# 		joined.arm.mrc <- JoinedArmMultiResultClass(samp_chr_arm = paste0(samp, "_", chr.arm),
# 													joined_samp_arms=chr.arm.joined.w.int.cpn.list,
# 													joined_samp_arms_inferred_ph_cpn=joined.samp.arms.ret.mrc.inferred.ph.cpn)
# 		#returnable.list <- list(joined.arm.mrc)
# 		#names(returnable.list) <- paste0(samp, "_", chr.arm)
# 		#return(returnable.list)
# 		#Endparallel version
#
#     joined.samp.arms.list[[paste0(samp, "_", chr.arm)]] <- joined.arm.mrc
#
# 	}#Parallel version closing brace.
#
# 	#So now we need a function to unscramble the parallelised results stuff:
#
# 	joined.samp.arms <- list()
# 	joined.samp.arms.inferred.ph.cpn <- list()
#
# 	for(samp in sample.vec) {
#
# 		joined.samp.arms[[samp]] <- list()
# 		joined.samp.arms.inferred.ph.cpn[[samp]] <- list()
#
# 		for(chr.arm in chr.arms.present.vec) {
#
# 			joined.samp.arms.inferred.ph.cpn[[samp]][[chr.arm]] <- joined.samp.arms.list[[paste0(samp, "_", chr.arm)]]$joined_samp_arms_inferred_ph_cpn[[chr.arm]]
# 			joined.samp.arms[[samp]][[chr.arm]] <- joined.samp.arms.list[[paste0(samp, "_", chr.arm)]]$joined_samp_arms_inferred_ph_cpn
#
# 		}
#
# 	}
#
# 	str(samp.arms, max=1)
# 	str(joined.samp.arms.inferred.ph.cpn, max=1)
#
# 	all.sample.joined.cpnA.cnpB.cpnTot.arm.segs.path <- paste0(PLOT_OUTPUT_DIR, patient, "_all_sample_joined_separate_cpnA_cnpB_cpnTot_arm_segs.RDS")
# 	saveRDS(joined.samp.arms, file = all.sample.joined.cpnA.cnpB.cpnTot.arm.segs.path )
#
# 	all.sample.joined.cpnA.cnpB.arm.segs.path <- paste0(PLOT_OUTPUT_DIR, patient, "_all_sample_joined_separate_cpnA_cnpB_arm_segs.RDS")
# 	saveRDS(joined.samp.arms.inferred.ph.cpn, file = all.sample.joined.cpnA.cnpB.arm.segs.path )
#
# 	all.samp.arm.cpnA.cpnB.list <- list()
# 	#Remake separate segments...? Needed for plotting?
# 	print("Merging separate joined cpnA and cpnB tracks into single set of segments for each sample.")
# 	for ( samp in sample.vec ) {
#
# 		#initialise list sample output:
# 		all.samp.arm.cpnA.cpnB.list[[samp]] <- list()
#
# 		samp.joined.arms.list <- joined.samp.arms.inferred.ph.cpn[[samp]]
# 		arms.present.vec <- names(samp.joined.arms.list)
#
# 		for (arm in arms.present.vec) {
#
# 			#print(str(samp.joined.arms.list[[arm]], max=1))
# 			samp.cpnA.arm.df <- samp.joined.arms.list[[arm]][["ph_cpnA"]]
# 			samp.cpnA.arm.grs <- makeGRangesFromDataFrame(samp.cpnA.arm.df,
# 					                         		keep.extra.columns=TRUE,
# 					                         		seqnames.field="chrom",
# 					                         		start.field="startpos",
# 					                         		end.field="endpos",
# 					                         		ignore.strand = TRUE,
# 					                         		starts.in.df.are.0based=FALSE)
# 			samp.cpnB.arm.df <- samp.joined.arms.list[[arm]][["ph_cpnB"]]
# 			samp.cpnB.arm.grs <- makeGRangesFromDataFrame(samp.cpnB.arm.df,
# 					                         		keep.extra.columns=TRUE,
# 					                         		ignore.strand=TRUE,
# 					                         		seqinfo=NULL,
# 					                         		seqnames.field="chrom",
# 					                         		start.field="startpos",
# 					                         		end.field="endpos",
# 					                         		starts.in.df.are.0based=FALSE)
# 			#print(samp.cpnA.arm.df)
# 			#print(samp.cpnB.arm.df)
# 			no.mcols.samp.cpnA.arm.grs <- samp.cpnA.arm.grs
# 			mcols(no.mcols.samp.cpnA.arm.grs) <- NULL
# 			no.mcols.samp.cpnB.arm.grs <- samp.cpnB.arm.grs
# 			mcols(no.mcols.samp.cpnB.arm.grs) <- NULL
# 			samp.arm.dj.grs <- disjoin(append(no.mcols.samp.cpnA.arm.grs, no.mcols.samp.cpnB.arm.grs), ignore.strand = TRUE)
# 			#print(samp.arm.dj.grs)
# 			dj.arm.cpnA.grs <- getOverlapGRanges(samp.arm.dj.grs, samp.cpnA.arm.grs)
# 			dj.arm.cpnB.grs <- getOverlapGRanges(samp.arm.dj.grs, samp.cpnB.arm.grs)
#
# 			#So now we've made the
# 			# print("Check making num_tot_cpn")
# 			# print(dj.arm.cpnA.grs$num_tot_cpn)
# 			# print(dj.arm.cpnB.grs$num_tot_cpn)
# 			# print(cbind(dj.arm.cpnA.grs$num_tot_cpn, dj.arm.cpnB.grs$num_tot_cpn))
# 			# print(apply(cbind(dj.arm.cpnA.grs$num_tot_cpn, dj.arm.cpnB.grs$num_tot_cpn), 1, min))
#
# 			arm.mcs.cpn.df <- data.frame(sample = rep(samp, length(samp.arm.dj.grs)),
# 						chrom = as.numeric(as.character(seqnames(samp.arm.dj.grs))),
# 						arm = rep(arm,  length(samp.arm.dj.grs)),
# 						startpos = start(samp.arm.dj.grs),
# 						endpos = end(samp.arm.dj.grs),
# 						num_hz = apply(cbind(dj.arm.cpnA.grs$num_ph_cpn, dj.arm.cpnB.grs$num_ph_cpn), 1, min),#Todo reconsider how we do this, may actually want to map it back to the original dataframe
# 						num_non_hz = apply(cbind(dj.arm.cpnA.grs$num_tot_cpn, dj.arm.cpnB.grs$num_tot_cpn), 1, min),#Todo reconsider how we do this, may actually want to map it back to the original dataframe
# 						mean_ph_cpnA = dj.arm.cpnA.grs$mean_ph_cpn,
# 						median_ph_cpnA = dj.arm.cpnA.grs$median_ph_cpn,
# 						int_ph_cpnA = dj.arm.cpnA.grs$int_ph_cpnA,
# 						mean_ph_cpnB = dj.arm.cpnB.grs$mean_ph_cpn,
# 						median_ph_cpnB = dj.arm.cpnB.grs$median_ph_cpn,
# 						int_ph_cpnB = dj.arm.cpnB.grs$int_ph_cpnB,
# 						int_tot_ph_cpn = dj.arm.cpnB.grs$int_ph_cpnB + dj.arm.cpnA.grs$int_ph_cpnA,
# 						mean_tot_cpn = dj.arm.cpnA.grs$mean_tot_cpn,
# 						median_tot_cpn = dj.arm.cpnA.grs$median_tot_cpn,
# 						stringsAsFactors = FALSE)
# 			#print(arm.mcs.cpn.df)
# 			all.samp.arm.cpnA.cpnB.list[[samp]][[arm]] <- arm.mcs.cpn.df
#
# 		}#end for (arm in arms.present.vec) {
#
# 	}#end for ( samp in sample.vec ) {
#
# 	all.samp.arm.cpnA.cpnB.list.path <- paste0(PLOT_OUTPUT_DIR, patient, "_all_samp_arm_cpnA_cpnB_list.RDS")
# 	saveRDS(all.samp.arm.cpnA.cpnB.list, file = all.samp.arm.cpnA.cpnB.list.path )
#
# 	#Make a per sample joined cpnA and cpnB data.frame in the same format as ASCAT (ie segments wherever either cpnA or cpnB change.)
# 	samp.joined.cpnA.cpnB.segs.df.list <- list()
# 	for(samp in sample.vec) {
#
# 		samp.joined.cpnA.cpnB.segs.df.list[[samp]] <- do.call(rbind, all.samp.arm.cpnA.cpnB.list[[samp]])
#
# 	}# for(samp in sample.vec) {
# 	all.samp.joined.cpnA.cpnB.df <- do.call(rbind, samp.joined.cpnA.cpnB.segs.df.list)
#
# 	all.samp.joined.cpnA.cpnB.df.path <- paste0(PLOT_OUTPUT_DIR, patient, "_all_samp_joined_cpnA_cpnB_df.RDS")
# 	saveRDS(all.samp.joined.cpnA.cpnB.df, file = all.samp.joined.cpnA.cpnB.df.path )
#
# 	print("Making new MCS segmentation based on joined segments... this will be used to generate MPhase MEDICC input...")
# 	#Now we want to make a per sample joined cpnA and cpnB dataframe where all breakpoints are shared between samples. Much like full.joint.segs.df above.
# 	all.samp.mcs.joined.cpnA.cpnB.df.list  <- list()
# 	#First thing we need to do is to come up with a disjoint set of Granges that incorporates all the breakpoints from each sample.
# 	all.samp.joined.cpnA.cpnB.grs <- makeGRangesFromDataFrame(all.samp.joined.cpnA.cpnB.df,
# 					                         	keep.extra.columns=TRUE,
# 					                         	ignore.strand=TRUE,
# 					                         	seqinfo=NULL,
# 					                         	seqnames.field="chrom",
# 					                         	start.field="startpos",
# 					                         	end.field="endpos",
# 					                         	starts.in.df.are.0based=FALSE)
# 	all.samp.joined.cpnA.cpnB.for.dj.grs <- all.samp.joined.cpnA.cpnB.grs
# 	mcols(all.samp.joined.cpnA.cpnB.for.dj.grs) <- NULL
# 	all.samp.joined.cpnA.cpnB.dj.grs <- disjoin(all.samp.joined.cpnA.cpnB.for.dj.grs)
# 	for (samp in sample.vec) {
# 		all.samp.joined.cpnA.cpnB.for.dj.grs
#
# 		samp.joined.cpnA.cpnB.grs <- all.samp.joined.cpnA.cpnB.grs[ all.samp.joined.cpnA.cpnB.grs$sample == samp ]
# 		samp.with.all.samp.dj.grs <- getOverlapGRanges(all.samp.joined.cpnA.cpnB.dj.grs, samp.joined.cpnA.cpnB.grs)
#
# 		samp.mcs.over.all.samp.cpn.df <- data.frame(sample = rep(samp, length(samp.with.all.samp.dj.grs)),
# 					chrom = as.numeric(as.character(seqnames(samp.with.all.samp.dj.grs))),
# 					arm = samp.with.all.samp.dj.grs$arm,
# 					startpos = start(samp.with.all.samp.dj.grs),
# 					endpos = end(samp.with.all.samp.dj.grs),
# 					num_hz = samp.with.all.samp.dj.grs$num_hz,#Todo reconsider how we do this, may actually want to map it back to the original dataframe
# 					num_non_hz = samp.with.all.samp.dj.grs$num_non_hz,#Todo reconsider how we do this, may actually want to map it back to the original dataframe
# 					mean_ph_cpnA = samp.with.all.samp.dj.grs$mean_ph_cpnA,
# 					median_ph_cpnA = samp.with.all.samp.dj.grs$median_ph_cpnA,
# 					int_ph_cpnA = samp.with.all.samp.dj.grs$int_ph_cpnA,
# 					mean_ph_cpnB = samp.with.all.samp.dj.grs$mean_ph_cpnB,
# 					median_ph_cpnB = samp.with.all.samp.dj.grs$median_ph_cpnB,
# 					int_ph_cpnB = samp.with.all.samp.dj.grs$int_ph_cpnB,
# 					int_tot_ph_cpn = samp.with.all.samp.dj.grs$int_ph_cpnB + samp.with.all.samp.dj.grs$int_ph_cpnA,
# 					mean_tot_cpn = samp.with.all.samp.dj.grs$mean_tot_cpn,
# 					median_tot_cpn = samp.with.all.samp.dj.grs$median_tot_cpn,
# 					stringsAsFactors = FALSE)
#
# 		all.samp.mcs.joined.cpnA.cpnB.df.list[[samp]] <-  samp.mcs.over.all.samp.cpn.df
# 	}#end for (samp in sample.vec) {
# 	all.samp.mcs.joined.cpnA.cpnB.df <- do.call(rbind, all.samp.mcs.joined.cpnA.cpnB.df.list)
# 	#So now we want to make the MEDICC input RDS based on this dataset's joined MPhase MCS results.
# 	all.samp.mcs.joined.cpnA.cpnB.df.path <- paste0(PLOT_OUTPUT_DIR, patient, "_all_samp_mcs_joined_cpnA_cpnB_df.RDS")
# 	saveRDS(all.samp.mcs.joined.cpnA.cpnB.df, file = all.samp.mcs.joined.cpnA.cpnB.df.path )
#
#
#
# 	print("Making post-MPhase MEDICC Input...")
# 	#So now we want to make the MEDICC input RDS based on this dataset's joined MPhase MCS results.
# 	post.MPhase.MEDICC.FASTA.files.list <- makeMEDICCInputListsFromSegmentationDF(segs.to.make.MEDICC.input.df = all.samp.mcs.joined.cpnA.cpnB.df,
# 											sample.vec = sample.vec,
# 											chrs.to.use.vec = CHRS_TO_USE_VEC,
# 											chr.col = "chrom",
# 											start.col = "startpos",
# 											end.col = "endpos",
# 											cpnA.col = "int_ph_cpnA",
# 											cpnB.col = "int_ph_cpnB"
#
# 											)
#
# 	post.MPhase.MEDICC.FASTA.files.list.path <- paste0(PLOT_OUTPUT_DIR, patient, "_post_MPhase_MEDICC_FASTA_file_list.RDS")
# 	saveRDS(post.MPhase.MEDICC.FASTA.files.list, file = post.MPhase.MEDICC.FASTA.files.list.path )
#
# } else if (OVERWRITE_JOIN_INIT_MCS_FLAG == FALSE) {
#
# 	all.sample.joined.cpnA.cnpB.cpnTot.arm.segs.path <- paste0(PLOT_OUTPUT_DIR, patient, "_all_sample_joined_separate_cpnA_cnpB_cpnTot_arm_segs.RDS")
# 	joined.samp.arms <- readRDS(file = all.sample.joined.cpnA.cnpB.cpnTot.arm.segs.path )
#
# 	all.sample.joined.cpnA.cnpB.arm.segs.path <- paste0(PLOT_OUTPUT_DIR, patient, "_all_sample_joined_separate_cpnA_cnpB_cpnTot_arm_segs.RDS")
# 	joined.samp.arms.inferred.ph.cpn <- readRDS(file = all.sample.joined.cpnA.cnpB.arm.segs.path)
#
# 	all.samp.arm.cpnA.cpnB.list.path <- paste0(PLOT_OUTPUT_DIR, patient, "_all_samp_arm_cpnA_cpnB_list.RDS")
# 	all.samp.arm.cpnA.cpnB.list <- readRDS(file = all.samp.arm.cpnA.cpnB.list.path )
#
# 	all.samp.joined.cpnA.cpnB.df.path <- paste0(PLOT_OUTPUT_DIR, patient, "_all_samp_joined_cpnA_cpnB_df.RDS")
# 	all.samp.joined.cpnA.cpnB.df <- readRDS(file = all.samp.joined.cpnA.cpnB.df.path )
#
# 	all.samp.mcs.joined.cpnA.cpnB.df.path <- paste0(PLOT_OUTPUT_DIR, patient, "_all_samp_mcs_joined_cpnA_cpnB_df.RDS")
# 	all.samp.mcs.joined.cpnA.cpnB.df <- readRDS(file = all.samp.mcs.joined.cpnA.cpnB.df.path )
#
# 	post.MPhase.MEDICC.FASTA.files.list.path <- paste0(PLOT_OUTPUT_DIR, patient, "_post_MPhase_MEDICC_FASTA_file_list.RDS")
# 	post.MPhase.MEDICC.FASTA.files.list <- readRDS(file = post.MPhase.MEDICC.FASTA.files.list.path )
#
#
# }
#
# ########################################################################
# # Generate some copy number plots
# #
# # <patient>_logr_BAF_raw_cpn_and_segments.png
# ########################################################################
#
# if (PLOTTING_FLAG == TRUE) {
#
# 	#We want to give each of relevant columns the same names in these dataframes and the plot them:
# 	#Chosem column names are:
# 	plottable.col.name.vec <- c("chr", "startpos", "endpos", "int_ph_cpnA", "int_ph_cpnB", "tot_int_cpn", "sample")
#
# 	#Original segments
# 	plottable.orig.joint.segs.df <- f.joint.segs.df[, c("chr","startpos", "endpos", "nMajor", "nMinor", "cnTotal", "sample")]
# 	plottable.orig.joint.segs.df <- plottable.orig.joint.segs.df[plottable.orig.joint.segs.df$chr %in% 1:22,]
# 	colnames(plottable.orig.joint.segs.df) <- plottable.col.name.vec
# 	#Full.joint.df segments
# 	plottable.mphase.mcs.segs.df <- full.joint.segs.df[, c("chr","startpos", "endpos", "ph_int_cpnA", "ph_int_cpnB", "total_int_cpn", "sample")]
# 	plottable.mphase.mcs.segs.df[is.na(plottable.mphase.mcs.segs.df)] <- -1
# 	colnames(plottable.mphase.mcs.segs.df) <- plottable.col.name.vec
#
# 	#Mphase joined segments
# 	plottable.mphase.samp.cpnA.cpnB.df <- all.samp.joined.cpnA.cpnB.df[, c("chrom","startpos", "endpos", "int_ph_cpnA", "int_ph_cpnB", "int_tot_ph_cpn", "sample")]
# 	colnames(plottable.mphase.samp.cpnA.cpnB.df) <- plottable.col.name.vec
#
# 	print("Plotting before and after joining copy number profiles for each sample....")
#
# 	png(file = paste(PLOT_OUTPUT_DIR, patient, "_mphase_mcs_cpn_across_genome", ".png", sep=""), width = 1200, height = 1400)
# 	plotCPNSegDFAccrossGenome(segs.df = plottable.mphase.mcs.segs.df, sample.vec = sample.vec, cpn.cols =c("int_ph_cpnA","int_ph_cpnB"), color.vec = c("green", "black"))
# 	dev.off()
#
# 	png(file = paste(PLOT_OUTPUT_DIR, patient, "_orig_ASCAT_across_genome", ".png", sep=""), width = 1200, height = 1400)
# 	plotCPNSegDFAccrossGenome(segs.df = plottable.orig.joint.segs.df, sample.vec = sample.vec, cpn.cols =c("int_ph_cpnA","int_ph_cpnB"), color.vec = c("green", "black"))
# 	dev.off()
#
# 	png(file = paste(PLOT_OUTPUT_DIR, patient, "_mphase_joined_across_genome", ".png", sep=""), width = 1200, height = 1400)
# 	plotCPNSegDFAccrossGenome(segs.df = plottable.mphase.samp.cpnA.cpnB.df, sample.vec = sample.vec, cpn.cols =c("int_ph_cpnA","int_ph_cpnB"), color.vec = c("green", "black"))
# 	dev.off()
#
# 	print("Plotting sample by sample logr, BAF, phased cpn, and cpn segmentation pre and post joining...")
# 	for (samp in sample.vec) {
#
# 		png(file = paste(PLOT_OUTPUT_DIR, samp, "_logr_BAF_raw_cpn_and_segments", ".png", sep=""), width = 1200, height = 1400)
#
# 			plotSampleLogRBAFCPNABsAndMSAIAcrossGenome(samp = samp,
# 												full.joint.segs.df = full.joint.segs.df,
# 												orig.segs.df = plottable.orig.joint.segs.df,
# 												mcs.segs.df = plottable.mphase.mcs.segs.df,
# 												new.joined.segs.df = plottable.mphase.samp.cpnA.cpnB.df,
# 												asas.ph.raw.cpn.lists = asas.ph.raw.cpn.lists,
# 												all.seg.CCF.vals.list = all.seg.CCF.list,
# 												all.seg.MSAI.by.BAF.ltnpd.present.list = all.seg.MSAI.by.BAF.ltnpd.present.list,
# 												patient.phased.BAF.df.list = patient.phased.BAF.df.list,
# 												asas.SNPs.w.both.lr.baf.snp.lists = asas.SNPs.w.both.lr.baf.snp.lists,
# 												asas.SNPS.lr.and.or.baf.lists = asas.SNPS.lr.and.or.baf.lists,
# 												cpnA.allele.color = '#de2d2699',
# 												cpnB.allele.color= '#3182bd99')
#
# 		dev.off()
#
# 	}
#
#
# } else if (PLOTTING_FLAG == FALSE) {
#
# 	print("Not plotting before and after joining copy number profiles for each sample as PLOTTING_FLAG is FALSE")
#
#
# }
#
# #We want to add error checking here:
# #Does ploidy change drastically from the original?
# #Is there no LOH or too much LOH?
# #Are there biologically infeasible amounts of homozygous deletion? If so we need to flag this.
#
# if ( JOINED_SEGS_ERROR_CHECKING_FLAG == TRUE) {
#
# 	#TODO refactor all this so that it works with the MPhase script itself
#
#
# 	fjs.df$orig_int_ploidy_recalculated <- NA
# 	fjs.df$mphase_int_ploidy <- NA
# 	fjs.df$orig_raw_ploidy_recalculated <- NA
# 	fjs.df$mphase_raw_ploidy <- NA
# 	fjs.df$mphase_int_ploidy <- NA
# 	fjs.df$num_na_mphase_tot_cpn_seg_remove <- NA
# 	fjs.df$num_na_orig_tot_cpn_seg_remove <- NA
# 	fjs.df$num_na_orig_or_mphase_tot_cpn_segs_to_remove <- NA
# 	fjs.df$orig_raw_ploidy_recalculated <- NA
# 	fjs.df$mphase_raw_ploidy <- NA
# 	fjs.df$orig_frac_LOH <- NA
# 	fjs.df$MPhase_frac_LOH <- NA
# 	fjs.df$ploidy_difference <- NA
# 	fjs.df$LOH_difference <- NA
# 	fjs.df$ploidy_issue_flag <- NA
# 	fjs.df$frac_LOH_issue_flag <- NA
#
#
# 	intersection.col <- grep("intersection", colnames(fjs.df), value=TRUE)
#
#
# 	int.fjs.df <- fjs.df[fjs.df[, "in_intersection_w_orig_segs"] %in% TRUE,]
#
#
# 	#get unique sample names
# 	sample.vec <- unique(fjs.df$sample)
# 	#frac LOH
# 	for (samp in sample.vec) {
#
# 	  samp.int.fjs.df <- int.fjs.df[int.fjs.df$sample == samp,]
#
#
# 	  num.na.mphase.tot.cpn.segs <- nrow(samp.int.fjs.df[ is.na(samp.int.fjs.df$total_int_cpn), ])
# 	  num.na.orig.tot.cpn.segs <- nrow(samp.int.fjs.df[ is.na(samp.int.fjs.df$orig_cnTotal),  ])
# 	  tot.na.orig.and.mphase.tot.cpn.segs <- nrow(samp.int.fjs.df[ is.na(samp.int.fjs.df$total_int_cpn) | is.na(samp.int.fjs.df$orig_cnTotal), ])
# 	  samp.int.fjs.df[samp.int.fjs.df$sample == samp,]$num_na_mphase_tot_cpn_seg_remove <- num.na.mphase.tot.cpn.segs
# 	  samp.int.fjs.df[samp.int.fjs.df$sample == samp,]$num_na_orig_tot_cpn_seg_remove <- num.na.orig.tot.cpn.segs
# 	  samp.int.fjs.df[samp.int.fjs.df$sample == samp,]$num_na_orig_or_mphase_tot_cpn_segs_to_remove <- tot.na.orig.and.mphase.tot.cpn.segs
#
# 	  samp.int.fjs.df <- samp.int.fjs.df[ (! is.na(samp.int.fjs.df$total_int_cpn)) & ( ! is.na(samp.int.fjs.df$orig_cnTotal)), ]
#
# 	  samp.int.fjs.grs <- makeGRangesFromDataFrame(samp.int.fjs.df[samp.int.fjs.df$sample == samp,],
# 	                                 keep.extra.columns=TRUE,
# 	                                 ignore.strand=FALSE,
# 	                                 seqinfo=NULL,
# 	                                 seqnames.field=c("chrom"),
# 	                                 start.field=c("Start_Position","start", "startpos"),
# 	                                 end.field=c("End_Position", "end", "stop", "endpos"),
# 	                                 strand.field="strand",
# 	                                 starts.in.df.are.0based=FALSE)
# 	  samp.int.fjs.grs$seg_width <- width(samp.int.fjs.grs)
# 	  print(head(samp.int.fjs.grs))
#
# 	  #Do we need to remove any segments that wil otherwise contribute to an NA
# 	  tot.genome.size <- sum(as.numeric(samp.int.fjs.grs$seg_width))
#
#
# 	  print("tot.genome.size")
# 	  print(tot.genome.size)
# 	  samp.int.fjs.grs$prop_of_tot_genome_size <- as.numeric(samp.int.fjs.grs$seg_width) / tot.genome.size
#
#
# 	  #Get the original raw ploidy...
# 	  samp.int.fjs.grs$tot_mphase_raw_contribute_to_ploidy <- as.numeric(samp.int.fjs.grs$total_raw_cpn) * as.numeric(samp.int.fjs.grs$prop_of_tot_genome_size)
# 	  samp.int.fjs.grs$tot_orig_raw_contribute_to_ploidy <- (as.numeric(samp.int.fjs.grs$orig_cpnA_raw) + as.numeric(samp.int.fjs.grs$orig_cpnB_raw)) * as.numeric(samp.int.fjs.grs$prop_of_tot_genome_size)
#
#
#
# 	  samp.int.fjs.grs$tot_mphase_int_contribute_to_ploidy <- as.numeric(samp.int.fjs.grs$total_int_cpn) * as.numeric(samp.int.fjs.grs$prop_of_tot_genome_size)
# 	  samp.int.fjs.grs$tot_orig_int_contribute_to_ploidy <- as.numeric(samp.int.fjs.grs$orig_cnTotal) * as.numeric(samp.int.fjs.grs$prop_of_tot_genome_size)
#
#
#
#
# 	  mphase.int.ploidy <- sum(samp.int.fjs.grs$tot_mphase_int_contribute_to_ploidy)
# 	  orig.int.ploidy <- sum(samp.int.fjs.grs$tot_orig_int_contribute_to_ploidy)
# 	  print("orig and MPhase int ploidies:")
# 	  print(orig.int.ploidy)
# 	  print(mphase.int.ploidy)
#
# 	  #Now get ploidy.
# 	  mphase.raw.ploidy <- sum(samp.int.fjs.grs$tot_mphase_raw_contribute_to_ploidy)
# 	  orig.raw.ploidy <- sum(samp.int.fjs.grs$tot_orig_raw_contribute_to_ploidy)
#
# 	  print("orig and MPhase raw ploidies:")
# 	  print(orig.raw.ploidy)
# 	  print(mphase.raw.ploidy)
# 	  #Get the new raw ploidy...
#
# 	  fjs.df[fjs.df$sample == samp,]$orig_int_ploidy_recalculated <- orig.int.ploidy
# 	  fjs.df[fjs.df$sample == samp,]$mphase_int_ploidy <- mphase.int.ploidy
# 	  fjs.df[fjs.df$sample == samp,]$orig_raw_ploidy_recalculated <- orig.raw.ploidy
# 	  fjs.df[fjs.df$sample == samp,]$mphase_raw_ploidy <- mphase.raw.ploidy
#
# 	  orig.frac.LOH  <- sum(samp.int.fjs.grs[ samp.int.fjs.grs$min_int_cpnX %in% 0]$prop_of_tot_genome_size)
# 	  mphase.frac.LOH <- sum(samp.int.fjs.grs[ samp.int.fjs.grs$orig_nMinor %in% 0 | samp.int.fjs.grs$orig_nMajor %in% 0]$prop_of_tot_genome_size)
# 	  print("Orig and MPhase fraction LOH:")
# 	  print(orig.frac.LOH)
# 	  print(mphase.frac.LOH)
#
# 	  fjs.df[fjs.df$sample == samp,]$orig_frac_LOH <- orig.frac.LOH
# 	  fjs.df[fjs.df$sample == samp,]$MPhase_frac_LOH <- mphase.frac.LOH
# 	  fjs.df[fjs.df$sample == samp,]$ploidy_difference <- mphase.int.ploidy - orig.int.ploidy
# 	  fjs.df[fjs.df$sample == samp,]$LOH_difference <- mphase.frac.LOH - orig.frac.LOH
#
# 	  #Problem sample flags:
# 	  fjs.df[fjs.df$sample == samp,]$ploidy_issue_flag <- abs(mphase.int.ploidy - orig.int.ploidy) >= 0.25
# 	  fjs.df[fjs.df$sample == samp,]$frac_LOH_issue_flag <- mphase.frac.LOH < 0.03
#
#
#
# 	}#end for (samp in sample.vec) {
#
#
#
#
#
# } else if (JOINED_SEGS_ERROR_CHECKING_FLAG == FALSE) {
#
# 	print("Not checking joined segs for errors as JOINED_SEGS_ERROR_CHECKING_FLAG is set to FALSE")
#
#
# }
#
#
#
# all.samp.joined.cpnA.cpnB.df <- readRDS(file = all.samp.joined.cpnA.cpnB.df.path )
#
#
# if (CALCULATE_JOINED_ARM_FOCAL_FLAG == TRUE) {
#
# 	getLogRValsForSampleSeg <- function(all.samp.SNP.list, samp ,chrom, start, end) {
#
# 		samp.df <- all.samp.SNP.list[[samp]]
# 		logr.vec <- samp.df[samp.df$chrom == as.numeric(chrom) &
# 				samp.df$pos >= as.numeric(start) &
# 				samp.df$pos <= as.numeric(end) ,]$logR
#
# 		return(logr.vec)
# 	}
#
# 	threshold.df <- data.frame(threshold_name = c("gain", "amp", "loss", "deep_loss"),
# 								threshold = c(GAIN_LOGR_THRESH, AMP_LOGR_THRESH,LOSS_LOGR_THRESH, DEEP_LOSS_LOGR_THRESH),
# 								cpn_type = c("gain", "gain", "loss", "loss"),
# 								stringsAsFactors = FALSE)
#
#
# 	all.samp.joined.cpnA.cpnB.df$callState <- apply(all.samp.joined.cpnA.cpnB.df, 1, function(x) {
#
# 		samp.seg.logr.vec <- getLogRValsForSampleSeg(all.samp.SNP.list = all.samp.SNP.list,
# 								   			samp = x["sample"],
# 								   			chrom = x["chrom"],
# 								   			start = x["startpos"],
# 								   			end = x["endpos"])
#
# 		cpn.event.to.return <- "neutral"
# 		for (i in 1:nrow(threshold.df)) {
#
# 			cpn.event.df <- checkSegGainOrLossSignificance(logr.vec = samp.seg.logr.vec,
# 										samp.to.check = x["sample"],
# 										ploidy.est.vec = ploidy.est.vec,
# 										purity.est.vec = purity.est.vec,
# 										threshold = threshold.df[i, ]$threshold,
# 										cpn.change.type = threshold.df[i, ]$cpn_type,
# 										min.num.for.logr.test = 3)
#
# 			if (is.na(cpn.event.df[, 7])) {
#
# 				return("neutral")
#
# 			} else if (cpn.event.df[, 7] < PVAL_THRESH_VS_LOGR_FOR_CPN_EVENT) {
#
# 				cpn.event.to.return <- threshold.df[i, ]$threshold_name
# 			}
#
# 		}#end for (i in 1:nrow(threshold.df)) {
#
# 		return(cpn.event.to.return)
# 		})
#
# 	all.samp.joined.cpnA.cpnB.df$MPhase_joined_seg_width <- as.numeric(all.samp.joined.cpnA.cpnB.df$endpos) - as.numeric(all.samp.joined.cpnA.cpnB.df$startpos)
# 	all.samp.joined.cpnA.cpnB.df$LOH_present <- "none"
# 	# all.samp.joined.cpnA.cpnB.df[all.samp.joined.cpnA.cpnB.df$int_ph_cpnB == 0
# 	# 							| all.samp.joined.cpnA.cpnB.df$int_ph_cpnA == 0,]$LOH_present  <- "LOH"
#
# 	all.samp.joined.cpnA.cpnB.df <- transform(all.samp.joined.cpnA.cpnB.df, LOH_present = ifelse(int_ph_cpnB == 0 & int_ph_cpnA == 0, "LOH", LOH_present))
#
# 	arm.size.by.samp.df <- aggregate(MPhase_joined_seg_width~sample+arm, all.samp.joined.cpnA.cpnB.df, sum)
# 	colnames(arm.size.by.samp.df)[colnames(arm.size.by.samp.df) == "MPhase_joined_seg_width" ] <- "total_samp_arm_size"
# 	all.samp.joined.cpnA.cpnB.df$prop_of_samp_arm <- 0
# 	all.samp.joined.cpnA.cpnB.w.arm.df <- merge(all.samp.joined.cpnA.cpnB.df, arm.size.by.samp.df, by=c("arm", "sample"))
# 	all.samp.joined.cpnA.cpnB.w.arm.df$prop_of_samp_arm <- all.samp.joined.cpnA.cpnB.w.arm.df$MPhase_joined_seg_width / all.samp.joined.cpnA.cpnB.w.arm.df$total_samp_arm_size
#
#
# 	cpn.size.by.samp.arm.df <- aggregate(MPhase_joined_seg_width~sample+arm+callState, all.samp.joined.cpnA.cpnB.df, sum)
# 	LOH.by.samp.arm.df <- aggregate(MPhase_joined_seg_width~sample+arm+LOH_present, all.samp.joined.cpnA.cpnB.df, sum)
#
# 	#Remove ttoal arm size as it has already been added:
# 	arm.size.by.samp.df$arm_gain_percent <- 0
# 	arm.size.by.samp.df$arm_loss_percent <- 0
# 	arm.size.by.samp.df$arm_LOH_percent <- 0
# 	for( i in 1:nrow(arm.size.by.samp.df)) {
#
# 		arm.size.by.samp.df[i,]$arm_gain_percent <-  sum(cpn.size.by.samp.arm.df[cpn.size.by.samp.arm.df$sample == arm.size.by.samp.df[i,]$sample
# 																		& cpn.size.by.samp.arm.df$arm == arm.size.by.samp.df[i,]$arm
# 																		& cpn.size.by.samp.arm.df$callState %in% c("amp", "gain"),]$MPhase_joined_seg_width) / arm.size.by.samp.df[i,]$total_samp_arm_size
#
# 		arm.size.by.samp.df[i,]$arm_loss_percent <-   sum(cpn.size.by.samp.arm.df[cpn.size.by.samp.arm.df$sample == arm.size.by.samp.df[i,]$sample
# 																		& cpn.size.by.samp.arm.df$arm == arm.size.by.samp.df[i,]$arm
# 																		& cpn.size.by.samp.arm.df$callState %in% c("deep_loss", "loss"),]$MPhase_joined_seg_width) / arm.size.by.samp.df[i,]$total_samp_arm_size
#
# 		arm.size.by.samp.df[i,]$arm_LOH_percent <-   sum(LOH.by.samp.arm.df[LOH.by.samp.arm.df$sample == arm.size.by.samp.df[i,]$sample
# 																		& LOH.by.samp.arm.df$arm == arm.size.by.samp.df[i,]$arm
# 																		& LOH.by.samp.arm.df$LOH_present %in% "LOH",]$MPhase_joined_seg_width) / arm.size.by.samp.df[i,]$total_samp_arm_size
#
# 	}
#
# 	INITIAL_ARM_EVENT_THRESHOLD <- 0.75
# 	LOWER_ARM_EVENT_THRESHOLD <- 0.7
# 	arm.size.by.samp.df$simple_MPhase_arm_event_vs_ploidy <- "none"
# 	arm.size.by.samp.df$simple_MPhase_arm_event_vs_ploidy_lower_thresh <- "none"
# 	arm.size.by.samp.df$simple_MPhase_LOH_arm_event <- "none"
# 	arm.size.by.samp.df$simple_MPhase_LOH_arm_event_lower_thresh <- "none"
#
# 	arm.size.by.samp.df <- transform(arm.size.by.samp.df,
# 										simple_MPhase_arm_event_vs_ploidy = ifelse(arm_gain_percent >= INITIAL_ARM_EVENT_THRESHOLD, "gain", simple_MPhase_arm_event_vs_ploidy))
# 	arm.size.by.samp.df <- transform(arm.size.by.samp.df,
# 										simple_MPhase_arm_event_vs_ploidy = ifelse(arm_loss_percent >= INITIAL_ARM_EVENT_THRESHOLD, "loss", simple_MPhase_arm_event_vs_ploidy))
# 	arm.size.by.samp.df <- transform(arm.size.by.samp.df,
# 										simple_MPhase_arm_event_vs_ploidy = ifelse((arm_loss_percent +  arm_gain_percent) < INITIAL_ARM_EVENT_THRESHOLD, "neutral", simple_MPhase_arm_event_vs_ploidy))
# 	arm.size.by.samp.df <- transform(arm.size.by.samp.df,
# 										simple_MPhase_LOH_arm_event = ifelse(arm_LOH_percent >= INITIAL_ARM_EVENT_THRESHOLD, "LOH", simple_MPhase_LOH_arm_event))
# 	#Get arms from
# 	gain.arm.vec <- arm.size.by.samp.df[arm.size.by.samp.df$simple_MPhase_arm_event_vs_ploidy == "gain",]$arm
# 	loss.arm.vec <- arm.size.by.samp.df[arm.size.by.samp.df$simple_MPhase_arm_event_vs_ploidy == "loss",]$arm
# 	LOH.arm.vec <- arm.size.by.samp.df[arm.size.by.samp.df$simple_MPhase_LOH_arm_event == "LOH",]$arm
#
#
# 	arm.size.by.samp.df <- transform(arm.size.by.samp.df,
# 									simple_MPhase_arm_event_vs_ploidy_lower_thresh = ifelse(arm %in% gain.arm.vec & arm_gain_percent >= LOWER_ARM_EVENT_THRESHOLD, "gain", simple_MPhase_arm_event_vs_ploidy))
# 	arm.size.by.samp.df <- transform(arm.size.by.samp.df,
# 									simple_MPhase_arm_event_vs_ploidy_lower_thresh = ifelse(arm %in% loss.arm.vec & arm_loss_percent >= LOWER_ARM_EVENT_THRESHOLD, "loss", simple_MPhase_arm_event_vs_ploidy))
# 	arm.size.by.samp.df <- transform(arm.size.by.samp.df,
# 									simple_MPhase_LOH_arm_event_lower_thresh = ifelse(arm %in% LOH.arm.vec & arm_LOH_percent >= LOWER_ARM_EVENT_THRESHOLD, "LOH", simple_MPhase_LOH_arm_event_lower_thresh))
#
# 	#Ok so now add these definitions to the original segmentation:
#
# 	arm.size.by.samp.df <- arm.size.by.samp.df[, grep("total_samp_arm_size", colnames(arm.size.by.samp.df), invert=TRUE)]
# 	all.samp.joined.cpnA.cpnB.arm.calls.df <- merge(all.samp.joined.cpnA.cpnB.w.arm.df, arm.size.by.samp.df, by=c("arm", "sample"))
#
# 	#We need to add a unique joined segment name that will allow us to map back once merge these segments with mcs segments:
# 	all.samp.joined.cpnA.cpnB.arm.calls.df$uniq_joined_seg_name <- paste0("joined/",all.samp.joined.cpnA.cpnB.arm.calls.df$sample,
# 																		"@",
# 																		all.samp.joined.cpnA.cpnB.arm.calls.df$chrom,
# 																		":",
# 																		all.samp.joined.cpnA.cpnB.arm.calls.df$startpos,
# 																		":",
# 																		all.samp.joined.cpnA.cpnB.arm.calls.df$endpos)
#
# 	#So let's define arm events as less than a proportion of the total arm size, let's say 15%?
# 	FOCAL_EVENT_PROP_OF_ARM <- 0.25
# 	FOCAL_EVENT_SIZE_THRESH <- 25000000
#
# 	all.samp.j.focal.arm.df <- transform(all.samp.joined.cpnA.cpnB.arm.calls.df, prop_arm_or_focal_in_joined = ifelse(prop_of_samp_arm <= FOCAL_EVENT_PROP_OF_ARM, "focal", "arm"))
# 	all.samp.j.focal.arm.df <- transform(all.samp.j.focal.arm.df, size_arm_or_focal_in_joined = ifelse(MPhase_joined_seg_width <= FOCAL_EVENT_SIZE_THRESH, "focal", "arm"))
#
# 	#So now define the ACROSS sample clonality of the events:
# 	#We need a measure of size always to allow us to step down the
# 	#So now we need to take these and look at the MCS version:
#
# 	all.samp.j.focal.arm.df.path <- paste0(PLOT_OUTPUT_DIR, patient, "_all_samp_j_focal_arm_df_path.RDS")
# 	saveRDS(all.samp.j.focal.arm.df, file = all.samp.j.focal.arm.df.path)
#
#
# } else if (CALCULATE_JOINED_ARM_FOCAL_FLAG == FALSE) {
#
# 	all.samp.j.focal.arm.df.path <- paste0(PLOT_OUTPUT_DIR, patient, "_all_samp_j_focal_arm_df_path.RDS")
# 	all.samp.j.focal.arm.df <- readRDS(file = all.samp.j.focal.arm.df.path)
#
# }
#
#
# #Now we need to free up some memory, so let's rm and then garbage collect some large variables...
# print("Removing large variables that are no longer required and performing garbage collection.")
# if(exists("snp.file.df")) { rm(snp.file.df) }
# if(exists("asas.SNPS.lr.and.or.baf.lists")) { rm(asas.SNPS.lr.and.or.baf.lists) }
# if(exists("asas.SNPs.w.both.lr.baf.snp.lists")) { rm(asas.SNPs.w.both.lr.baf.snp.lists) }
# if(exists("all.seg.CCF.list")) { rm(all.seg.CCF.list) }
# if(exists("all.seg.int.ob.xp.cpn.df.list")) { rm(all.seg.int.ob.xp.cpn.df.list) }
# if(exists("all.seg.CCF.diff.from.exp.list")) { rm(all.seg.CCF.diff.from.exp.list) }
# if(exists("all.seg.clonal.in.all.samps.list")) { rm(all.seg.clonal.in.all.samps.list) }
# if(exists("all.seg.greater.than.two.diff.list")) { rm(all.seg.greater.than.two.diff.list) }
# if(exists("asas.ph.raw.cpn.lists")) { rm(asas.ph.raw.cpn.lists) }
# if(exists("all.samp.non.hz.logr.and.tot.cpn.df.list")) { rm(all.samp.non.hz.logr.and.tot.cpn.df.list) }
# if(exists("all.seg.phasing.info.list")) { rm(all.seg.phasing.info.list) }
# if(exists("patient.phased.BAF.df.list")) { rm(patient.phased.BAF.df.list) }
# if(exists("all.seg.MSAI.by.BAF.ltnpd.present.list")) { rm(all.seg.MSAI.by.BAF.ltnpd.present.list) }
# if(exists("post.MPhase.MEDICC.FASTA.files.list")) { rm(post.MPhase.MEDICC.FASTA.files.list) }
# if(exists("joined.samp.arms")) { rm(joined.samp.arms) }
# if(exists("joined.samp.arms.inferred.ph.cpn")) { rm(joined.samp.arms.inferred.ph.cpn) }
# if(exists("all.samp.arm.cpnA.cpnB.list")) { rm(all.samp.arm.cpnA.cpnB.list) }
# if(exists("samp.arms")) { rm(samp.arms) }
# gc()
#
#
#
# ########################################################################
# # From Tom:
# # there are output files that contain the letters “pj”
# # these are from essentially a second run through of the algorithm
# # that runs it after the original segments produced by ascat are phased and then rejoined using the phasing info and run through the whole algorithm again
# # so there are actually two files in the output that are named in exactly the same way
# # for now just use the LTX121_full_joint_segs_df.RDS
# ########################################################################
#
# #So now we need to calculate whether any of these are gains or losses, we can do this in the same way as in the original MPhase or just apply the function as we did:
# print("Generating new paths:")
# #Reading in data object paths
# asas.SNPS.lr.and.or.baf.lists.pj.path <- paste0(PLOT_OUTPUT_DIR, patient, "_pj_asas_SNPS_lr_and_or_baf_lists_lists.RDS")
# asas.SNPs.w.both.lr.baf.snp.lists.pj.path <- paste0(PLOT_OUTPUT_DIR, patient, "_pj_asas_SNPs_w_both_lr_baf_snp_lists.RDS")
# min.cons.segs.df.write.pj.path <- paste0(PLOT_OUTPUT_DIR, patient, "_pj_min_cons_segs_df.txt")
# sorted.all.bp.all.samps.segs.anno.df.pj.path <- paste0(PLOT_OUTPUT_DIR, patient, "_pj_sorted_all_bp_all_samps_segs_anno.df.txt")
# min.cons.segs.df.pj.path <- paste0(PLOT_OUTPUT_DIR, patient, "_pj_min_cons_segs_df.RDS")
# processed.orig.cpn.segs.df.pj.path <- paste0(PLOT_OUTPUT_DIR, patient, "_pj_processed_orig_cpn_segs_df.RDS")
# #Analysis object paths
# all.seg.phasing.info.list.pj.path <- paste0(PLOT_OUTPUT_DIR, patient, "_pj_all_seg_phasing_info_list.RDS")
# patient.phased.BAF.df.list.pj.path <- paste0(PLOT_OUTPUT_DIR, patient, "_pj_phased_BAF_df_list.RDS")
# all.seg.MSAI.by.BAF.ltnpd.present.list.pj.path <- paste0(PLOT_OUTPUT_DIR, patient, "_pj_all_seg_MSAI_by_BAF_ltnpd_present_list.RDS")
# all.samp.non.hz.logr.and.tot.cpn.df.list.pj.path <- paste0(PLOT_OUTPUT_DIR, patient, "_pj_all_samp_non_hz_logr_and_tot_cpn_df_list.RDS")
# asas.ph.raw.cpn.lists.pj.path <- paste0(PLOT_OUTPUT_DIR, patient, "_pj_asas_ph_raw_cpn_lists.RDS")
# samps.int.cpn.raw.cpn.by.chr.list.pj.path <- paste0(PLOT_OUTPUT_DIR, patient, "_pj_samps_int_cpn_raw_cpn_by_chr_list.RDS")
# all.seg.CCF.list.pj.path <- paste0(PLOT_OUTPUT_DIR, patient, "_pj_all_seg_CCF_list.RDS")
# all.seg.int.ob.xp.cpn.df.list.pj.path <- paste0(PLOT_OUTPUT_DIR, patient, "_pj_all_seg_int_ob_xp_cpn_df_list.RDS")
# all.seg.CCF.diff.from.exp.list.pj.path <- paste0(PLOT_OUTPUT_DIR, patient, "_pj_all_seg_CCF_diff_from_exp_list.RDS")
# all.seg.clonal.in.all.samps.list.pj.path <- paste0(PLOT_OUTPUT_DIR, patient, "_pj_all_seg_clonal_in_all_samps_list.RDS")
# all.seg.greater.than.two.diff.list.pj.path <- paste0(PLOT_OUTPUT_DIR, patient, "_pj_all_seg_greater_than_two_diff_list.RDS")
#
#
#
#
# #All SNPs in these lists must have both a logr and baf value: all sample all seg logr baf lists is asas.lr.baf
#
# if (OVERWRITE_JOINED_CCF_CALC_FLAG == TRUE) {
#
# 	#Do we need this loading in here?
# 	all.samp.mcs.joined.cpnA.cpnB.df <- readRDS(file = all.samp.mcs.joined.cpnA.cpnB.df.path)
#
# 	#cpn.chr.start.end.df <- sorted.all.bp.all.samps.segs.anno.df[sorted.all.bp.all.samps.segs.anno.df$full_sample==sorted.all.bp.all.samps.segs.anno.df[1,]$full_sample,1:3]
# 	cpn.chr.start.end.df.pj <- all.samp.mcs.joined.cpnA.cpnB.df[, c("chrom", "startpos", "endpos")]
# 	mock.cpn.matrix.pj <- matrix(data=0, nrow=nrow(cpn.chr.start.end.df.pj), ncol = length(sample.vec))
# 	colnames(mock.cpn.matrix.pj) <- sample.vec
# 	min.cons.segs.df.pj <- cbind(cpn.chr.start.end.df.pj, mock.cpn.matrix.pj)
#
# 	seg.name.vec.pj <- paste(min.cons.segs.df.pj[, 1], min.cons.segs.df.pj[, 2], min.cons.segs.df.pj[, 3], sep = ":")
#
# 	print("cleanAndErrorCheckSNPsParallel for hz sites alone on the joined segments:")
# 	source(paste0(current.dir, "data_input_and_QC_functions_v2.R"))
# 	asas.SNPs.w.both.lr.baf.snp.lists.pj <- cleanAndErrorCheckSNPsParallel(all.samp.SNP.list = all.samp.SNP.list,
# 																		mcs.seg.to.use.df = min.cons.segs.df.pj,
# 																		allow.logr.only.flag = FALSE,
# 																		data.type = data.format.to.use)
# 	if (length(unique(sapply(asas.SNPs.w.both.lr.baf.snp.lists.pj, length))) != 1) {
#
# 		stop("Not all segments are present in all samples' snp assignations, likely something has gone wrong, is the parallelisation working correctly?")
#
# 	}
# 	saveRDS(asas.SNPs.w.both.lr.baf.snp.lists.pj, file = asas.SNPs.w.both.lr.baf.snp.lists.pj.path)
#
# 	print("cleanAndErrorCheckSNPsParallel for both non-hz and hz sites on the joined segments:")
#
# 	asas.SNPS.lr.and.or.baf.lists.pj <- cleanAndErrorCheckSNPsParallel(all.samp.SNP.list = all.samp.SNP.list,
# 																	mcs.seg.to.use.df = min.cons.segs.df.pj,
# 																	allow.logr.only.flag = TRUE,
# 																	data.type = data.format.to.use)
#
# 	print("Removing all.samp.SNP.list as no longer required and performing garbage collection.")
# 	if(exists("all.samp.SNP.list")) { rm(all.samp.SNP.list) }
# 	gc()
#
#
# 	if (length(unique(sapply(asas.SNPS.lr.and.or.baf.lists.pj, length))) != 1) {
#
# 		stop("Not all segments are present in all samples' snp assignations, likely something has gone wrong, is the parallelisation working correctly?")
#
# 	}
# 	saveRDS(asas.SNPS.lr.and.or.baf.lists.pj, file = asas.SNPS.lr.and.or.baf.lists.pj.path)
#
# 	print("Getting phasing to apply to all samples using getCombinedBAFDipTestPhasingParallel...")
# 	source(paste0(current.dir, "consensus_phasing_from_multiple_samples_functions_v2.R"))
# 	all.seg.phasing.info.list.pj <- getCombinedBAFDipTestPhasingParallel(asas.lr.baf.lists = asas.SNPs.w.both.lr.baf.snp.lists.pj,
# 	                                            seg.name.vec = seg.name.vec.pj,
# 	                                            ploidy.est.vec = ploidy.est.vec,
# 	                                            purity.est.vec = purity.est.vec,
# 	                                            gamma.val = GAMMA_VAL,
# 	                                            min.num.hz.snp.phasing = MIN_NUM_HZ_SNP_PHASING,
# 	                                            min.prop.snp.sep.threshold = MIN_PROP_SNP_SEP_THRESH,
# 	                                            dip.test.p.val.threshold = PVAL_DIPTEST_AI_THRESH,
# 	                                            comb.wilcox.snp.sep.p.val.threshold = PVAL_COMB_WILCOX_PH_THRESH,
# 	                                            allow.dip.test.without.ASCAT.to.phase = ALLOW_DIPTEST_ALONE_PH_FLAG
# 	                                            )
#
# 	saveRDS(all.seg.phasing.info.list.pj, file = all.seg.phasing.info.list.pj.path)
#
# 	print("Generating phased BAF dataframes...")
# 	patient.phased.BAF.df.list.pj <- createAllSegsPhasedBAFDataFrame(asas.lr.baf.lists = asas.SNPs.w.both.lr.baf.snp.lists.pj,
# 	                                                     sample.to.phase.seg.list = all.seg.phasing.info.list.pj,
# 	                                                     min.num.hz.snp.phasing = MIN_NUM_HZ_SNP_PHASING
# 	                                                     )
# 	if(length(patient.phased.BAF.df.list.pj) <= 1) {
#
# 		print("Error in generating phased BAF dataframes")
#
#
# 	}
# 	saveRDS(patient.phased.BAF.df.list.pj, file = patient.phased.BAF.df.list.pj.path)
#
#
# 	#Detect MSAI in the original BAF before it is processed alongside logr to get the Raw allele specifc cpn to be used in the rest of the analysis.
# 	print("Performing LTNPD MSAI detection on the BAF...")
# 	all.seg.MSAI.by.BAF.ltnpd.present.list.pj <- getAllSegsMSAIDetectionFromBAFListLTNPD(seg.name.vec = seg.name.vec.pj,
# 																				  sample.vec  = sample.vec,
# 																				  asas.SNPS.lr.and.or.baf.lists = asas.SNPS.lr.and.or.baf.lists.pj,
# 																				  all.seg.phasing.info.list = all.seg.phasing.info.list.pj,
# 																				  MSAI.BAF.LTPND.pval.thresh = MSAI_BAF_LTNPD_PVAL_THRESHOLD)
# 	saveRDS(all.seg.MSAI.by.BAF.ltnpd.present.list.pj, file = all.seg.MSAI.by.BAF.ltnpd.present.list.pj.path)
#
#
# 	print("Generating datapoints with LogR total cpn dataframes...")
# 	#This is where the SNP filtering occurs. We can fix this by adding by passing an option of which data type is being analysed...
# 	all.samp.non.hz.logr.and.tot.cpn.df.list.pj <- getAllSampLogrAlonePositionsLogrAndTotalCPNDFsList(sample.vec = sample.vec,
# 																				asas.SNPS.lr.and.or.baf.lists = asas.SNPS.lr.and.or.baf.lists.pj,
# 																				ploidy.est.vec = ploidy.est.vec,
# 																				purity.est.vec = purity.est.vec,
# 																				gamma.val = GAMMA_VAL,
# 																				data.format.to.use = data.format.to.use)
# 	saveRDS(all.samp.non.hz.logr.and.tot.cpn.df.list.pj, file = all.samp.non.hz.logr.and.tot.cpn.df.list.pj.path)
#
#
#
#
#
# 	if (PLOTTING_FLAG == TRUE) {
#
# 		png(file = paste(PLOT_OUTPUT_DIR, patient, "_MSAI_detection_BAF_phased_pj", ".png", sep=""), width = 1200, height = 1400)
#
# 		print("Plotting phased BAF across the genome...")
# 		plotAllSampleBAFsAndMSAIAcrossGenome(sample.vec, all.seg.MSAI.by.BAF.ltnpd.present.list.pj, patient.phased.BAF.df.list.pj, cpnA.allele.color = '#de2d2699', cpnB.allele.color= '#3182bd99')
#
# 		dev.off()
#
# 	} else if (PLOTTING_FLAG == FALSE) {
#
# 		print("Not plotting phased BAF as flag is FALSE")
#
# 	}
#
#
# 	#Add this the unphaseable stuff too maybe? Not necessary.
# 	print("Applying phasing to all samples to get phased raw cpn...")
# 	asas.ph.raw.cpn.lists.pj <- getAllSegAllSamplePhasedRawCPNList(sample.vec = sample.vec,
# 																ploidy.est.vec = ploidy.est.vec,
# 																purity.est.vec = purity.est.vec,
# 																min.cons.segs.df = min.cons.segs.df.pj,
# 																all.seg.phasing.info.list = all.seg.phasing.info.list.pj,
# 																asas.both.lr.baf.lists = asas.SNPs.w.both.lr.baf.snp.lists.pj,
# 																gamma.val = GAMMA_VAL,
# 																min.num.hz.snp.phasing = MIN_NUM_HZ_SNP_PHASING)
# 	saveRDS(asas.ph.raw.cpn.lists.pj, file = asas.ph.raw.cpn.lists.pj.path)
#
# 	print("Calculating CCFs using phased raw cpn...")
# 	all.seg.CCF.list.pj <- list()
# 	all.seg.int.ob.xp.cpn.df.list.pj <- list()
# 	all.seg.CCF.diff.from.exp.list.pj <- list()
# 	all.seg.clonal.in.all.samps.list.pj <- list()
# 	all.seg.greater.than.two.diff.list.pj <- list()
#
# 	samps.int.cpn.raw.cpn.by.chr.list.pj <- generateAllSamplesByChromIntCPNUnderlyingRawValuesLists(asas.ph.raw.cpn.lists = asas.ph.raw.cpn.lists.pj,
# 																								 sample.vec = sample.vec,
# 																								 seg.name.vec = seg.name.vec.pj,
# 																								 cpnX.vec=c("cpnA", "cpnB"),
# 																								 reduced.seg.flag = FALSE)
# 	saveRDS(samps.int.cpn.raw.cpn.by.chr.list.pj, file = samps.int.cpn.raw.cpn.by.chr.list.pj.path)
#
# 	for (cpnX in c("cpnA", "cpnB")) {
# 	  print(cpnX)
# 	  #Here is where we're going to have iniatalise on all the lists
# 	  all.seg.CCF.list.pj[[cpnX]] <- list()
# 	  all.seg.int.ob.xp.cpn.df.list.pj[[cpnX]] <- list()
# 	  all.seg.CCF.diff.from.exp.list.pj[[cpnX]] <- list()
# 	  all.seg.clonal.in.all.samps.list.pj[[cpnX]] <- list()
# 	  all.seg.greater.than.two.diff.list.pj[[cpnX]] <- list()
#
# 	  for ( seg.name in seg.name.vec.pj ) {
#
# 	    all.seg.CCF.list.pj[[cpnX]][[seg.name]] <- list()
# 	    all.seg.int.ob.xp.cpn.df.list.pj[[cpnX]][[seg.name]] <- list()
# 	    all.seg.CCF.diff.from.exp.list.pj[[cpnX]][[seg.name]] <- list()
# 	    all.seg.clonal.in.all.samps.list.pj[[cpnX]][[seg.name]] <- list()
# 	    all.seg.greater.than.two.diff.list.pj[[cpnX]][[seg.name]] <- list()
#
# 	  }
#
# 	  for (seg.name in seg.name.vec.pj) {
#
# 	    if (is.na(asas.ph.raw.cpn.lists.pj[[1]][[seg.name]]$ph_cpnA_vec[1])) {
#
# 	      all.seg.CCF.list.pj[[cpnX]][[seg.name]] <- NA
# 	      all.seg.int.ob.xp.cpn.df.list.pj[[cpnX]][[seg.name]] <- NA
# 	      all.seg.CCF.diff.from.exp.list.pj[[cpnX]][[seg.name]] <- NA
# 	      all.seg.clonal.in.all.samps.list.pj[[cpnX]][[seg.name]] <- NA
# 	      all.seg.greater.than.two.diff.list.pj[[cpnX]][[seg.name]] <- NA
# 	      next
#
# 	    }
#
# 	    #does both A and B, todo, refactor so only does one
# 	    s.samps.int.ob.xp.cpn.df.pj <- getSegObsExpCPNBySampleDF(asas.ph.raw.cpn.lists = asas.ph.raw.cpn.lists.pj,
# 	                                                          sample.vec = sample.vec,
# 	                                                          seg.name = seg.name,
# 	                                                          cpnX.vec = cpnX,
# 	                                                          reduced.seg.flag = FALSE)
#
# 	    #does only cpnX
# 	    seg.cpnX.diff.from.exp.list.pj <- calculateWhetherSegCPNDiffFromIntegerExpectedAndAllMatchingSampCPN(asas.ph.raw.cpn.lists = asas.ph.raw.cpn.lists.pj,
# 	                                                                                                      s.samps.int.ob.xp.cpn.df = s.samps.int.ob.xp.cpn.df.pj,
# 	                                                                                                      seg.name = seg.name,
# 	                                                                                                      samps.int.cpn.raw.cpn.by.chr.list = samps.int.cpn.raw.cpn.by.chr.list.pj,
# 	                                                                                                      cpnX = cpnX,
# 	                                                                                                      chrs.can.use.for.ccf.test = CHRS_TO_USE_VEC,
# 	                                                                                                      reduced.seg.flag = FALSE)
# 	    #does only cpnX
# 	    is.segment.clonal.in.all.samps.flag.pj <- determineCCFClonalityOfSegment(s.samps.int.ob.xp.cpn.df = s.samps.int.ob.xp.cpn.df.pj,
# 	                                                                          seg.name = seg.name,
# 	                                                                          asas.ph.raw.cpn.lists = asas.ph.raw.cpn.lists.pj,
# 	                                                                          all.samp.ph.BAF.df.list = patient.phased.BAF.df.list.pj,
# 	                                                                          purity.est.vec = purity.est.vec,
# 	                                                                          cpnX = cpnX,
# 	                                                                          min.CCF.to.consider.val = MIN_CCF_TO_CONSIDER,
# 	                                                                          confidence.interval.threshold = CONF_INT_FOR_CPNX_BAF_VS_EXP)
#
# 	    #does only cpnX
# 	    is.there.cpnX.difference.greater.than.two.flag.pj <- calculateTwoValueDiffFlag(s.samps.int.ob.xp.cpn.df = s.samps.int.ob.xp.cpn.df.pj,
# 	                                                                seg.cpnX.diff.from.exp.list = seg.cpnX.diff.from.exp.list.pj,
# 	                                                                cpnX = cpnX,
# 	                                                                min.CCF.to.consider.val = MIN_CCF_TO_CONSIDER,
# 	                                                                confidence.interval.threshold = CONF_INT_FOR_CPNX_BAF_VS_EXP,
# 	                                                                min.num.other.chrs = MIN_NUM_CHRS)
#
# 	    #does only cpnX
# 	    seg.psc.low.high.cpn.and.use.preset.ccf.flag.vec.pj <- getPossibleLowAndHighSubclonalCPNStates(seg.name = seg.name,
# 	                                                                                                cpnX = cpnX,
# 	                                                                                                s.samps.int.ob.xp.cpn.df = s.samps.int.ob.xp.cpn.df.pj,
# 	                                                                                                all.clonal.flag = is.segment.clonal.in.all.samps.flag.pj,
# 	                                                                                                two.value.diff.flag = is.there.cpnX.difference.greater.than.two.flag.pj,
# 	                                                                                                min.CCF.to.consider.val = MIN_CCF_TO_CONSIDER)
#
#
# 	    seg.psc.cpnX.low.pj  <- seg.psc.low.high.cpn.and.use.preset.ccf.flag.vec.pj[1]
# 	    seg.psc.cpnX.high.pj  <- seg.psc.low.high.cpn.and.use.preset.ccf.flag.vec.pj[2]
# 	    seg.use.default.ccf.pj <- as.logical(seg.psc.low.high.cpn.and.use.preset.ccf.flag.vec.pj[3])
#
# 	    #update segment for each sample
# 	    for ( samp  in sample.vec ) {
#
# 	      all.seg.CCF.list.pj[[cpnX]][[seg.name]][[samp]] <- list()
# 	      #ss.ph.cpnX.vals.vec <- all.pw.cn.het.comp.lists[[samp]][[seg.name]][[paste0("ph_", cpnX, "_vec")]]
# 	      ss.ph.cpnX.vals.vec.pj <- asas.ph.raw.cpn.lists.pj[[samp]][[seg.name]][[paste0("ph_", cpnX, "_vec")]]
#
# 	      seg.samp.CCF.psc.cpn.df.pj <- getSegCPNXMeanAndAllRawCCFsDF(ss.mean.cpnX.obs = s.samps.int.ob.xp.cpn.df.pj[samp, paste0("mean_obs_", cpnX)],
# 	                                                              ss.ph.cpnX.vals.vec = ss.ph.cpnX.vals.vec.pj,
# 	                                                              seg.psc.cpnX.low = seg.psc.cpnX.low.pj,
# 	                                                              seg.psc.cpnX.high = seg.psc.cpnX.high.pj,
# 	                                                              seg.use.default.ccf.flag = seg.use.default.ccf.pj,
# 	                                                              all.clonal.flag = is.segment.clonal.in.all.samps.flag.pj,
# 	                                                              max.cpn.val = MAX_CPN_FOR_SUBCLONALITY)
#
# 	      all.seg.CCF.list.pj[[cpnX]][[seg.name]][[samp]][[paste0("poss_sc_", cpnX, "_low")]] <- seg.samp.CCF.psc.cpn.df.pj$poss_sc_cpnX_low
# 	      all.seg.CCF.list.pj[[cpnX]][[seg.name]][[samp]][[paste0("poss_sc_", cpnX, "_high")]] <- seg.samp.CCF.psc.cpn.df.pj$poss_sc_cpnX_high
# 	      all.seg.CCF.list.pj[[cpnX]][[seg.name]][[samp]][[paste0("mean_frac_", cpnX, "_low")]] <- seg.samp.CCF.psc.cpn.df.pj$mean_frac_cpnX_low
# 	      all.seg.CCF.list.pj[[cpnX]][[seg.name]][[samp]][[paste0("mean_frac_", cpnX, "_high")]] <- seg.samp.CCF.psc.cpn.df.pj$mean_frac_cpnX_high
# 	      all.seg.CCF.list.pj[[cpnX]][[seg.name]][[samp]][[paste0(cpnX, "_orig_ph_obs_cpnA_vec")]]   <- asas.ph.raw.cpn.lists.pj[[samp]][[seg.name]][[paste0("ph_cpnA_vec")]]
# 	      all.seg.CCF.list.pj[[cpnX]][[seg.name]][[samp]][[paste0(cpnX, "_orig_ph_obs_cpnB_vec")]]   <- asas.ph.raw.cpn.lists.pj[[samp]][[seg.name]][[paste0("ph_cpnB_vec")]]
# 	      all.seg.CCF.list.pj[[cpnX]][[seg.name]][[samp]][[paste0("clean_ph_obs_", cpnX, "_vec")]]   <- seg.samp.CCF.psc.cpn.df.pj$clean_ph_obs_cpnX
# 	      all.seg.CCF.list.pj[[cpnX]][[seg.name]][[samp]][[paste0("frac_", cpnX, "_low_vec")]]  <- seg.samp.CCF.psc.cpn.df.pj$frac_cpnX_low
# 	      all.seg.CCF.list.pj[[cpnX]][[seg.name]][[samp]][[paste0("frac_", cpnX, "_high_vec")]]  <- seg.samp.CCF.psc.cpn.df.pj$frac_cpnX_high
#
# 	    }#end  for ( samp  in sample.vec ) {
#
# 	    all.seg.int.ob.xp.cpn.df.list.pj[[cpnX]][[seg.name]] <- s.samps.int.ob.xp.cpn.df.pj
# 	    all.seg.clonal.in.all.samps.list.pj[[cpnX]][[seg.name]] <- is.segment.clonal.in.all.samps.flag.pj
# 	    all.seg.greater.than.two.diff.list.pj[[cpnX]][[seg.name]] <- is.there.cpnX.difference.greater.than.two.flag.pj
# 	    all.seg.CCF.diff.from.exp.list.pj[[cpnX]][[seg.name]] <- seg.cpnX.diff.from.exp.list.pj
#
# 	  }#for (seg.name %in% seg.name.vec
#
# 	}#end for (cpnX in c("cpnA", "cpnB")) {
#
# 	#Save CCF data:
# 	saveRDS(all.seg.CCF.list.pj, file = all.seg.CCF.list.pj.path)
# 	saveRDS(all.seg.int.ob.xp.cpn.df.list.pj, file = all.seg.int.ob.xp.cpn.df.list.pj.path)
# 	saveRDS(all.seg.CCF.diff.from.exp.list.pj, file = all.seg.CCF.diff.from.exp.list.pj.path)
# 	saveRDS(all.seg.clonal.in.all.samps.list.pj, file = all.seg.clonal.in.all.samps.list.pj.path)
# 	saveRDS(all.seg.greater.than.two.diff.list.pj, file = all.seg.greater.than.two.diff.list.pj.path)
#
#
# } else if (OVERWRITE_JOINED_CCF_CALC_FLAG == FALSE) {
#
# 	all.seg.phasing.info.list.pj <- readRDS(file = all.seg.phasing.info.list.pj.path)
# 	patient.phased.BAF.df.list.pj <- readRDS(file = patient.phased.BAF.df.list.pj.path)
# 	all.seg.MSAI.by.BAF.ltnpd.present.list.pj <- readRDS(file = all.seg.MSAI.by.BAF.ltnpd.present.list.pj.path)
# 	all.samp.non.hz.logr.and.tot.cpn.df.list.pj <- readRDS(file = all.samp.non.hz.logr.and.tot.cpn.df.list.pj.path)
# 	asas.ph.raw.cpn.lists.pj <- readRDS(file = asas.ph.raw.cpn.lists.pj.path)
# 	samps.int.cpn.raw.cpn.by.chr.list.pj <- readRDS(file = samps.int.cpn.raw.cpn.by.chr.list.pj.path)
# 	all.seg.CCF.list.pj <- readRDS(file = all.seg.CCF.list.pj.path)
# 	all.seg.int.ob.xp.cpn.df.list.pj <- readRDS(file = all.seg.int.ob.xp.cpn.df.list.pj.path)
# 	all.seg.CCF.diff.from.exp.list.pj <- readRDS(file = all.seg.CCF.diff.from.exp.list.pj.path)
# 	all.seg.clonal.in.all.samps.list.pj <- readRDS(file = all.seg.clonal.in.all.samps.list.pj.path)
# 	all.seg.greater.than.two.diff.list.pj <- readRDS(file = all.seg.greater.than.two.diff.list.pj.path)
#
# }
# # #Todo remove the dataframe initialisation from the function and modify it to work on a segment by segment basis.
#
#
# if (OVERWRITE_JOINED_MCS_ADD_INFO_FLAG == TRUE) {
#
# 	#Let's initialise the first joint output:
# 	print("Initialising minimum consistent segmentation (MCS) output dataframe...")
# 	empty.full.joint.segs.df.pj <- initialiseFullOutputDataframe(seg.name.vec = names(asas.ph.raw.cpn.lists.pj[[1]]),
# 															ploidy.est.vec = ploidy.est.vec,
# 															purity.est.vec = purity.est.vec,
# 															make.row.names.seg.name.alone = FALSE)
#
#
# 	#This adds the raw cn for both segments that could be phased and those that couldn't
# 	print("Adding new copy number calls to MCS output based on the phased raw copy number...")
# 	# new.cpn.full.joint.segs.df <- addNewRawAndIntCPNCallsToFinalSegmentation(seg.name.vec = names(asas.ph.raw.cpn.lists[[1]]),
# 	# 																		sample.vec = names(purity.est.vec),
# 	# 																		all.samp.logr.df.list = all.samp.non.hz.logr.and.tot.cpn.df.list,
# 	# 																		all.samp.ph.BAF.df.list = patient.phased.BAF.df.list,
# 	# 																		all.seg.CCF.vals.list = all.seg.CCF.list,
# 	# 																		asas.lr.baf.lists = asas.SNPs.w.both.lr.baf.snp.lists,
# 	# 																		all.seg.df = empty.full.joint.segs.df,
# 	# 																		cpnX.vec = c("cpnA", "cpnB")
# 	# 																		)
#
# 	# source(paste0(current.dir, "/final_output_segmentation_file_generation_functions_v3.R"))
# 	# source(paste0(current.dir, "/utility_MPhase_functions.R"))
# 	new.cpn.full.joint.segs.df.pj <- addNewRawAndIntCPNCallsToFinalSegmentationParallel(seg.name.vec = names(asas.ph.raw.cpn.lists.pj[[1]]),
# 																			sample.vec = names(purity.est.vec),
# 																			all.samp.logr.df.list = all.samp.non.hz.logr.and.tot.cpn.df.list.pj,
# 																			all.samp.ph.BAF.df.list = patient.phased.BAF.df.list.pj,
# 																			all.seg.CCF.vals.list = all.seg.CCF.list.pj,
# 																			asas.lr.baf.lists = asas.SNPs.w.both.lr.baf.snp.lists.pj,
# 																			all.seg.df = empty.full.joint.segs.df.pj,
# 																			cpnX.vec = c("cpnA", "cpnB")
# 																			)
#
#
# 	# print("Mapping MCS segments back to original input and copy number states...")
# 	# full.joint.segs.w.old.cpn.vs.ploidy.df.pj <- addOriginalSegmentInfoToFinalSegmentationNew(new.cpn.full.joint.segs.df.pj, sorted.all.bp.all.samps.segs.anno.df.pj)
#
# 	#sorted.all.bp.all.samps.segs.anno.df.pj
#
# 	print("Adding CCF information and subclonal cpn calls to MCS output...")
# 	clonality.and.ccf.joint.segs.df.pj <- addWithinSampleClonalityAndCCFToFinalSegmentation(seg.name.vec = names(asas.ph.raw.cpn.lists.pj[[1]])
# 	                                ,all.seg.CCF.vals.list = all.seg.CCF.list.pj
# 	                                ,all.seg.CCF.diff.from.exp.list = all.seg.CCF.diff.from.exp.list.pj
# 	                                ,all.seg.df = new.cpn.full.joint.segs.df.pj#full.joint.segs.w.old.cpn.vs.ploidy.df.pj#new.cpn.full.joint.segs.df
# 	                                ,sample.vec = names(purity.est.vec)
# 	                                ,clonal.pval = PVAL_THRESH_VS_MIN_CCF
# 	                                ,diff.from.oth.p.val = PVAL_THRESH_VS_OTH_CHRS_SAME_CPN
# 	                                ,min.CCF.to.consider.val = MIN_CCF_TO_CONSIDER
# 	                                ,min.cpn.to.consider.val = MIN_CPN_TO_CONSIDER
# 	                                ,min.num.other.chrs = MIN_NUM_CHRS
# 	                                ,chrs.to.use.vec = CHRS_TO_USE_VEC
# 	                                ,cpnX.vec = c("cpnA", "cpnB")
# 	                                ,max.cpn.val = MAX_CPN_FOR_SUBCLONALITY
# 	                                )
#
# 	print("Adding total copy number estimates to MCS output for segments with only logr information...")
# 	#source(paste0(current.dir, "/final_output_segmentation_file_generation_functions_v3.R"))
# 	#source(paste0(current.dir, "/utility_MPhase_functions.R"))
# 	full.joint.segs.df.w.cpn.vs.ploidy.df.pj <- addLogrDefinedCPNEventsToOutputSegmentationParallel(sample.vec = names(purity.est.vec),
# 	                                                        							seg.name.vec = names(asas.ph.raw.cpn.lists.pj[[1]]),
# 	                                                        							all.seg.df = clonality.and.ccf.joint.segs.df.pj,
# 	                                                        							ploidy.est.vec = ploidy.est.vec,
# 	                                                        							purity.est.vec = purity.est.vec,
# 	                                                        							all.samp.logr.df.list = all.samp.non.hz.logr.and.tot.cpn.df.list.pj,
# 	                                                        							all.samp.ph.BAF.df.list = patient.phased.BAF.df.list.pj,
# 	                                                       								gain.threshold = GAIN_LOGR_THRESH,
# 	                                                        							loss.threshold = LOSS_LOGR_THRESH,
# 	                                                        							amp.threshold = AMP_LOGR_THRESH,
# 	                                                        							deep.loss.threshold = DEEP_LOSS_LOGR_THRESH,
# 	                                                        							tot.cpn.vs.ploidy.pval.thresh = PVAL_THRESH_VS_LOGR_FOR_CPN_EVENT)
#
# 	print(head(full.joint.segs.df.w.cpn.vs.ploidy.df.pj[,"cpn_event_vs_ploidy"]))
# 	print("Adding cpn event vs ploidy by logr measures to MCS output...")
# 	#source(paste0(current.dir, "/final_output_segmentation_file_generation_functions_v3.R"))
# 	full.joint.segs.logr.vs.ploidy.clonality.df.pj <- determineAcrossSampleClonalityOfLogrEventsVsPloidy(sample.vec = names(purity.est.vec),
# 	                                                        							seg.name.vec = names(asas.ph.raw.cpn.lists.pj[[1]]),
# 	                                                        							all.seg.df = full.joint.segs.df.w.cpn.vs.ploidy.df.pj,
# 	                                                        							all.samp.logr.df.list = all.samp.non.hz.logr.and.tot.cpn.df.list.pj,
# 	                                                        							all.samp.ph.BAF.df.list = patient.phased.BAF.df.list.pj,
# 	                                                       								gain.threshold = GAIN_LOGR_THRESH,
# 	                                                        							loss.threshold = LOSS_LOGR_THRESH,
# 	                                                        							amp.threshold = AMP_LOGR_THRESH,
# 	                                                        							deep.loss.threshold = DEEP_LOSS_LOGR_THRESH,
# 	                                                        							tot.cpn.vs.ploidy.pval.thresh = PVAL_THRESH_VS_LOGR_FOR_CPN_EVENT)
#
#
# 	print("Testing whether any gains should be reclassified as amplifications...")
# 	# full.joint.check.amps.vs.ploidy.clonality.df <- retestSubclonalAmps(all.seg.df = full.joint.segs.logr.vs.ploidy.clonality.df,
# 	# 																	all.samp.logr.df.list = all.samp.non.hz.logr.and.tot.cpn.df.list,
# 	# 																	all.samp.ph.BAF.df.list = patient.phased.BAF.df.list,
# 	# 																	gain.threshold = GAIN_LOGR_THRESH,
# 	# 																	amp.threshold = AMP_LOGR_THRESH,
# 	# 																	p.val.for.converting.gain.to.amp = PVAL_THRESH_CONVERT_GAIN_TO_AMP)
#
# 	print(head(full.joint.segs.logr.vs.ploidy.clonality.df.pj[,c("across_samp_cl_cpn_vs_ploidy_present","across_samp_sc_cpn_vs_ploidy_present")]))
#
#
# 	#These are kept in separate functions as the original segmentaiton was refactored so it could run when segments are rejoined at the end of the script.
# 	print("Adding phasing information and MSAI detection by raw cpn calls to MCS output...")
# 	full.joint.segs.w.phasing.df.pj <- addPhasingInformationToOutputSegmentation(sample.vec = names(purity.est.vec),
# 	                                        seg.name.vec = names(asas.ph.raw.cpn.lists.pj[[1]]),
# 	                                        sample.to.phase.seg.list = all.seg.phasing.info.list.pj,
# 	                                        all.pw.cn.het.comp.lists = NULL,
# 	                                        all.segs.df = full.joint.segs.logr.vs.ploidy.clonality.df.pj)
#
#
# 	# #Now we want to add in the allelic imbalance clonal and subclonal inrormation to the final data frame.
# 	# print("Adding across sample clonal/subclonal summary to MCS output...")
# 	# full.joint.segs.w.phasing.and.raw.df <- addOriginalCPNValuesAndCallsVsPloidyToOutputSegmentation(sample.vec = sample.vec,
# 	#                                         seg.name.vec = names(asas.ph.raw.cpn.lists[[1]]),
# 	#                                         processed.orig.cpn.segs.df = processed.orig.cpn.segs.df,
# 	#                                         all.segs.df = full.joint.segs.w.phasing.df)
#
# 	print("Adding MSAI detected by LTNPD method to the MCS output...")
# 	full.joint.segs.w.ltnpd.MSAI.df.pj <- addMSAIDetectionByLTNPDToOutputSegmentation(sample.vec = sample.vec,
# 												seg.name.vec = names(asas.ph.raw.cpn.lists.pj[[1]]),
# 												all.seg.MSAI.by.BAF.ltnpd.present.list = all.seg.MSAI.by.BAF.ltnpd.present.list.pj,
# 												all.segs.df = full.joint.segs.w.phasing.df.pj)#full.joint.segs.w.phasing.and.raw.df)
#
#
# 	print(head(full.joint.segs.w.ltnpd.MSAI.df.pj[( ! is.na(full.joint.segs.w.ltnpd.MSAI.df.pj$any_across_samp_MSAI)) &
# 												full.joint.segs.w.ltnpd.MSAI.df.pj$any_across_samp_MSAI == TRUE,
# 												c("MSAI", "any_across_samp_MSAI")]))
#
#
# 	print("Defining any parallel gain or loss events detected via MSAI and addingto MCS output...")
# 	full.joint.segs.w.parallel.cpn.events.df.pj <- addParallelEventsByMSAIToOutputSegmentation(all.segs.df = full.joint.segs.w.ltnpd.MSAI.df.pj)
#
#
# 	print(head(full.joint.segs.w.parallel.cpn.events.df.pj[(! is.na(full.joint.segs.w.parallel.cpn.events.df.pj$MSAI)),
# 														c("any_across_samp_MSAI", "MSAI",
# 															"cpn_event_vs_ploidy",
# 															"is_parallel_gain_and_or_amp",
# 															"is_parallel_loss",
# 															"is_parallel_LOH",
# 															"is_parallel_loss_and_or_LOH")]))
#
# 	print("Adding comparison between original cpn calls and MPhase copy number calling to MCS output...")
# 	full.joint.segs.within.and.changed.cpn.info.df.pj <- addWithinSampleChangeToCpnAndAllelicImbalanceCalls(full.joint.segs.w.parallel.cpn.events.df.pj,
# 																										sample.vec = sample.vec,
# 																										seg.name.vec = seg.name.vec.pj)
#
#
#
#
#
# 	print("Adding across sample estimates of allelic imbalance and LOH clonality to MCS output...")
# 	full.joint.segs.w.across.samp.AI.stats.df.pj <- addAcrossSampleAlleleSpecificClonalityCalls(all.seg.df = full.joint.segs.within.and.changed.cpn.info.df.pj,
# 																							sample.vec = sample.vec,
# 																							seg.name.vec = seg.name.vec.pj)
#
# 	print(head(full.joint.segs.w.across.samp.AI.stats.df.pj[!is.na(full.joint.segs.w.across.samp.AI.stats.df.pj$within_samp_sc_cpnA_present), c("any_within_samp_sc_cpnA",
# 	"any_within_samp_sc_cpnB",
# 	"across_samp_cl_AI",
# 	"across_samp_sc_AI",
# 	"across_samp_cl_LOH",
# 	"across_samp_sc_LOH")]))
#
#
# 	full.joint.segs.w.across.plotting.cols.df.pj <- addPanCanColumnDefinitionToSegDF(full.joint.segs.w.across.samp.AI.stats.df.pj)
#
# 	print(full.joint.segs.w.across.plotting.cols.df.pj[, c("pt_seg_cl_loss_LOH_no_MSAI", "pt_seg_sc_loss_LOH_incl_MSAI", "pt_seg_cl_gain_amp_no_MSAI", "pt_seg_sc_gain_amp_incl_MSAI")])
#
# 	full.joint.segs.df.pj <- full.joint.segs.w.across.plotting.cols.df.pj
# 	#Remove segments without any datapoints.
# 	#Todo, merge these int
# 	#print("Removing any segments that do not have underlying datapoints, this is an artefact of the minimum consistent segmentation performed by MPhase...")
# 	#full.joint.segs.df <- full.joint.segs.df[full.joint.segs.df$num_un_ph_hz_sites != 0 | full.joint.segs.df$num_ph_hz_sites != 0  | full.joint.segs.df$num_non_hz_sites != 0,]
# 	#print("Removing segments that do not map back to any original segments")
# 	#full.joint.segs.df <- full.joint.segs.df[full.joint.segs.df$orig_seg_name != "NA:NA:NA",]
#
#
# 	full.joint.segs.df.pj.path <- paste0(PLOT_OUTPUT_DIR, patient, "_full_joint_segs_df_pj.RDS")
# 	saveRDS(full.joint.segs.df.pj, file=full.joint.segs.df.pj.path)
#
#
# 	##
# } else if (OVERWRITE_JOINED_MCS_ADD_INFO_FLAG == FALSE) {
#
#
# 	full.joint.segs.df.pj.path <- paste0(PLOT_OUTPUT_DIR, patient, "_full_joint_segs_df_pj.RDS")
# 	full.joint.segs.df.pj <- readRDS(file = full.joint.segs.df.pj.path )
#
#
#
# }
#
#
# if (PLOTTING_FLAG == TRUE) {
#
#
# 	####Plot all genome wide plots
# 	print("Plotting CCFs for cpnA and cpnB across the genome using cpn definition of subclonality...")
# 	png(file = paste(PLOT_OUTPUT_DIR, patient, "_CCF_w_cpn_thresh_cpnA_pj", ".png", sep = ""), width = 1200, height = 1400)
#
# 	plotAllSampleCCFsAcrossGenome(asas.ph.raw.cpn.lists = asas.ph.raw.cpn.lists.pj,
# 	                              full.joint.segs.df = full.joint.segs.df.pj,
# 	                              sample.vec = sample.vec,
# 	                              all.seg.CCF.vals.list = all.seg.CCF.list.pj,
# 	                              cpnX = "cpnA",
# 	                              use.cpn.or.ccf.thresh = "CPN")
#
# 	dev.off()
#
# 	png(file = paste(PLOT_OUTPUT_DIR, patient, "_CCF_w_cpn_thresh_cpnB_pj", ".png", sep = ""), width = 1200, height = 1400)
#
# 	plotAllSampleCCFsAcrossGenome(asas.ph.raw.cpn.lists = asas.ph.raw.cpn.lists.pj,
# 	                              full.joint.segs.df = full.joint.segs.df.pj,
# 	                              sample.vec = sample.vec,
# 	                              all.seg.CCF.vals.list = all.seg.CCF.list.pj,
# 	                              cpnX = "cpnB",
# 	                              use.cpn.or.ccf.thresh = "CPN")
# 	dev.off()
#
# 	################################################
#
# } else if (PLOTTING_FLAG == FALSE) {
#
#
# 	print("Not plotting CCFs for cpnA and cpnB across the genome using cpn definition of subclonality as plotting flag is FALSE.")
#
#
# }
#
#
# #### So #TODO FIX THIS:
# #We need a unique name from the joined segs that we want to add to
# if (CALCULATE_FOCAL_EVENT_CLONALITY_FLAG == TRUE) {
#
# 	all.pj.segs.merged.w.mcs.df.list <- list()
# 	for ( samp in sample.vec ) {
#
# 		#Make Granges of the joined segments so that we can get the overlap of joint segments with them:
# 		samp.pj.seg.gr <- makeGRangesFromDataFrame(all.samp.j.focal.arm.df[all.samp.j.focal.arm.df$sample == samp,],
# 	                               keep.extra.columns=TRUE,
# 	                               ignore.strand=FALSE,
# 	                               seqinfo=NULL,
# 	                               seqnames.field=c("chrom"),
# 	                               start.field=c("Start_Position","start", "startpos"),
# 	                               end.field=c("End_Position", "end", "stop", "endpos"),
# 	                               strand.field="strand",
# 	                               starts.in.df.are.0based=FALSE)
#
#
# 		samp.full.joint.segs.gr.pj <- makeGRangesFromDataFrame(full.joint.segs.df.pj[full.joint.segs.df.pj$sample == samp,],
# 	                               keep.extra.columns=TRUE,
# 	                               ignore.strand=FALSE,
# 	                               seqinfo=NULL,
# 	                               seqnames.field=c("chrom"),
# 	                               start.field=c("Start_Position","start", "startpos"),
# 	                               end.field=c("End_Position", "end", "stop", "endpos"),
# 	                               strand.field="strand",
# 	                               starts.in.df.are.0based=FALSE)
#
# 		samp.merged.pj.segs.w.mcs.gr <- getOverlapGRanges(samp.full.joint.segs.gr.pj, samp.pj.seg.gr)
# 		all.pj.segs.merged.w.mcs.df.list[[samp]] <- as.data.frame(samp.merged.pj.segs.w.mcs.gr)
# 		#Convert back to dataframe:
#
# 	}
# 	#Now bind these back together
#
# 	all.pj.segs.merged.w.mcs.df <- do.call(rbind, all.pj.segs.merged.w.mcs.df.list)
# 	rownames(all.pj.segs.merged.w.mcs.df) <- rownames(full.joint.segs.df.pj)
# 	#Tidy data
# 	#Remove extra columns:
# 	extra.cols.to.remove.vec <- c("sample.1", "strand", "chrom", "query_hit_idx", "subject_hit_idx", "query_start", "query_end", "subject_start", "subject_end", "overlap_start")
# 	all.pj.segs.merged.w.mcs.df <- all.pj.segs.merged.w.mcs.df[, colnames(all.pj.segs.merged.w.mcs.df)[ ! colnames(all.pj.segs.merged.w.mcs.df) %in% extra.cols.to.remove.vec]]
# 	#rename columns:
# 	colnames(all.pj.segs.merged.w.mcs.df)[colnames(all.pj.segs.merged.w.mcs.df) == "seqnames"] <- "chrom"
# 	colnames(all.pj.segs.merged.w.mcs.df)[colnames(all.pj.segs.merged.w.mcs.df) == "start"] <- "startpos"
# 	colnames(all.pj.segs.merged.w.mcs.df)[colnames(all.pj.segs.merged.w.mcs.df) == "end"] <- "endpos"
# 	colnames(all.pj.segs.merged.w.mcs.df)[colnames(all.pj.segs.merged.w.mcs.df) == "callState"] <- "MPhase_joined_seg_vs_ploidy"
#
# 	#TODO fix the nameing conventions for phcpnA/B in pre and post segmentation... Not ideal at the moment
# 	all.pj.segs.merged.w.mcs.df <- transform(all.pj.segs.merged.w.mcs.df, nMinor = ifelse(is.na(nMinor), min(int_ph_cpnB, int_ph_cpnA), nMinor))
# 	all.pj.segs.merged.w.mcs.df <- transform(all.pj.segs.merged.w.mcs.df, nMajor = ifelse(is.na(nMinor), max(int_ph_cpnB, int_ph_cpnA), nMajor))
# 	all.pj.segs.merged.w.mcs.df <- transform(all.pj.segs.merged.w.mcs.df, nAraw = ifelse(is.na(nAraw), max(median_ph_cpnB, median_ph_cpnA), nAraw))
# 	all.pj.segs.merged.w.mcs.df <- transform(all.pj.segs.merged.w.mcs.df, nBraw = ifelse(is.na(nBraw), min(median_ph_cpnB, median_ph_cpnA), nBraw))
# 	all.pj.segs.merged.w.mcs.df <- transform(all.pj.segs.merged.w.mcs.df, ph_raw_cpnA = ifelse(is.na(ph_raw_cpnA), median_ph_cpnA, ph_raw_cpnA))
# 	all.pj.segs.merged.w.mcs.df <- transform(all.pj.segs.merged.w.mcs.df, ph_raw_cpnB = ifelse(is.na(ph_raw_cpnB), median_ph_cpnB, ph_raw_cpnB))
# 	all.pj.segs.merged.w.mcs.df <- transform(all.pj.segs.merged.w.mcs.df, min_int_cpnX = ifelse(is.na(min_int_cpnX), min(int_ph_cpnB, int_ph_cpnA), min_int_cpnX))
# 	all.pj.segs.merged.w.mcs.df <- transform(all.pj.segs.merged.w.mcs.df, max_int_cpnX = ifelse(is.na(max_int_cpnX), max(int_ph_cpnB, int_ph_cpnA), max_int_cpnX))
# 	all.pj.segs.merged.w.mcs.df <- transform(all.pj.segs.merged.w.mcs.df, cpn_event_vs_ploidy = ifelse(is.na(cpn_event_vs_ploidy), MPhase_joined_seg_vs_ploidy, cpn_event_vs_ploidy))
# 	all.pj.segs.merged.w.mcs.df <- transform(all.pj.segs.merged.w.mcs.df, total_raw_cpn = ifelse(is.na(total_raw_cpn), ph_raw_cpnA + ph_raw_cpnB, total_raw_cpn))
# 	#Todo, decide on median or mean:
# 	all.pj.segs.merged.w.mcs.df <- transform(all.pj.segs.merged.w.mcs.df, ph_raw_cpnA = ifelse(is.na(ph_raw_cpnA), mean_ph_cpnA, ph_raw_cpnA))
# 	all.pj.segs.merged.w.mcs.df <- transform(all.pj.segs.merged.w.mcs.df, ph_raw_cpnB = ifelse(is.na(ph_raw_cpnB), mean_ph_cpnB, ph_raw_cpnB))
# 	all.pj.segs.merged.w.mcs.df <- transform(all.pj.segs.merged.w.mcs.df, ph_raw_cpnB = ifelse(ph_raw_cpnB < 0, 0, ph_raw_cpnB))
# 	all.pj.segs.merged.w.mcs.df <- transform(all.pj.segs.merged.w.mcs.df, ph_raw_cpnA = ifelse(ph_raw_cpnA < 0, 0, ph_raw_cpnA))
#
# 	#TODO should move this potentially, perhaps into the final segmentation bit of the code/
# 	#Calculating the raw ploidy and contribution of each segment to that ploidy.
# 	all.pj.segs.merged.w.mcs.df$MPhase_raw_ploidy <- 0
# 	for( samp in sample.vec) {
#
# 		samp.fjs.df.pj <- all.pj.segs.merged.w.mcs.df[all.pj.segs.merged.w.mcs.df$sample == samp,]
# 		samp.fjs.grs.pj <- makeGRangesFromDataFrame(samp.fjs.df.pj,
# 		                               keep.extra.columns=TRUE,
# 		                               ignore.strand=FALSE,
# 		                               seqinfo=NULL,
# 		                               seqnames.field=c("chrom"),
# 		                               start.field=c("Start_Position","start", "startpos"),
# 		                               end.field=c("End_Position", "end", "stop", "endpos"),
# 		                               strand.field="strand",
# 		                               starts.in.df.are.0based=FALSE)
# 		samp.fjs.grs.pj$seg_width <- width(samp.fjs.grs.pj)
# 		samp.pj.tot.genome.size <- sum(as.numeric(samp.fjs.grs.pj$seg_width))
# 		samp.fjs.grs.pj$prop_of_tot_genome_size <- as.numeric(samp.fjs.grs.pj$seg_width) / samp.pj.tot.genome.size
# 		samp.fjs.grs.pj$tot_mphase_raw_contribute_to_ploidy <- as.numeric(samp.fjs.grs.pj$total_raw_cpn) * as.numeric(samp.fjs.grs.pj$prop_of_tot_genome_size)
# 		samp.mphase.raw.ploidy <- sum(samp.fjs.grs.pj$tot_mphase_raw_contribute_to_ploidy)
# 		all.pj.segs.merged.w.mcs.df[all.pj.segs.merged.w.mcs.df$sample == samp,]$MPhase_raw_ploidy <- samp.mphase.raw.ploidy
#
# 	}
#
#
# 	pj.mcs.df <- all.pj.segs.merged.w.mcs.df
# 	#MCS sample and segment name required:
# 	pj.mcs.df$samp_seg_name <- paste0("MCS/",
# 										pj.mcs.df$sample,
# 										"@",
# 										pj.mcs.df$chrom,
# 										":",
# 										pj.mcs.df$startpos,
# 										":",
# 										pj.mcs.df$endpos)
#
# 	#Make focal in joined id column:
# 	pj.mcs.df <- transform(pj.mcs.df, pj_focal_by_size_and_prop = ifelse(
# 							prop_arm_or_focal_in_joined %in% "focal" & size_arm_or_focal_in_joined %in% "focal",
# 							"focal",
# 							"arm")
# 							)
#
#
# 	#potential.cl.foc.mcs.segs.vec all the samp_seg_names that match the condition of the subclonal focal we're looking for
# 	#all chr/arm post joining minimum consistent segmentation segmentation dataframe
# 	#focal column that will be updated by this function
# 	#What is the type of the event we're looking to capture and the string will be placed in the focal id name
# 	markClonalFocalEvents <- function(potential.cl.foc.mcs.segs.vec, pj.mcs.df, col.to.update, event.type) {
#
# 		arm.vec <- unique(pj.mcs.df$arm)
# 		num.samples  <- length(unique(pj.mcs.df$sample))
#
# 		for (arm in arm.vec ) {
#
# 			arm.pj.mcs.df <- pj.mcs.df[pj.mcs.df$arm %in% arm,]
# 			arm.pj.mcs.df$old_focal_gain_definition <- arm.pj.mcs.df[, col.to.update]
# 			#So now we need to create the column that tell sus whether we actually can look for runs:
# 			arm.pj.mcs.df$potential_new_focal <- arm.pj.mcs.df$samp_seg_name %in% potential.cl.foc.mcs.segs.vec
# 			arm.unclassified.focal.mat <- matrix(FALSE, nrow=length(unique(arm.pj.mcs.df$sample)), ncol=length(unique(arm.pj.mcs.df$seg_name)))
# 			colnames(arm.unclassified.focal.mat) <- unique(arm.pj.mcs.df$seg_name)
# 			rownames(arm.unclassified.focal.mat) <- unique(arm.pj.mcs.df$sample)
# 			#arm.unclassified.focal.mat.allele <- matrix(FALSE, nrow=length(unique(arm.pj.mcs.df$sample)), ncol=length(unique(arm.pj.mcs.df$seg_name)))
# 			#colnames(arm.unclassified.focal.mat.allele)  <- colnames(arm.unclassified.focal.mat) <- unique(arm.pj.mcs.df$seg_name)
# 			#rownames(arm.unclassified.focal.mat.allele)  <- rownames(arm.unclassified.focal.mat) <- unique(arm.pj.mcs.df$sample)
# 			for(samp in unique(arm.pj.mcs.df$sample)) {
# 				for(seg in unique(arm.pj.mcs.df$seg_name)){
# 					arm.unclassified.focal.mat[samp, seg] <- arm.pj.mcs.df[arm.pj.mcs.df$sample == samp & arm.pj.mcs.df$seg_name == seg,]$potential_new_focal
# 					#arm.unclassified.focal.mat.allele[samp, seg] <- arm.pj.mcs.df[arm.pj.mcs.df$sample == samp & arm.pj.mcs.df$seg_name == seg,]$dominant_allele
# 				}#end for(seg in unique(arm.pj.mcs.df$seg_name)){
# 			}#end for(samp in unique(arm.pj.mcs.df$sample)) {
#
# 			f.acr.samp.vec <- colSums(arm.unclassified.focal.mat) == num.samples
# 			unclassified.cl.focals.remain <- any(f.acr.samp.vec)
#
# 			while (unclassified.cl.focals.remain) {
# 				print(unclassified.cl.focals.remain)
#
# 				cl.foc.rle.df <- with(rle(f.acr.samp.vec),
# 									data.frame(number = values,
# 												start = cumsum(lengths) - lengths + 1,
# 												end = cumsum(lengths))[order(values),])
#
# 				curr.cl.fcl.rle.idx <- min(which(cl.foc.rle.df$number == TRUE))
# 				curr.start.cl.fcl.idx <- cl.foc.rle.df[curr.cl.fcl.rle.idx,]$start
# 				curr.end.cl.fcl.idx <- cl.foc.rle.df[curr.cl.fcl.rle.idx,]$end
# 				#So now we want the middle
# 				mid.cl.fcl.idx <- curr.start.cl.fcl.idx + ceiling((curr.end.cl.fcl.idx - curr.start.cl.fcl.idx) / 2)
#
# 				#So now we need to know whether there are contiguous segments
# 				centre.seg.name <- colnames(arm.unclassified.focal.mat)[mid.cl.fcl.idx]
# 				samps.w.centre.seg <- rownames(arm.unclassified.focal.mat[arm.unclassified.focal.mat[,mid.cl.fcl.idx]==TRUE,mid.cl.fcl.idx, drop = FALSE])
# 				print(centre.seg.name)
# 				joined.seg.size.vec <- arm.pj.mcs.df[arm.pj.mcs.df$sample %in% samps.w.centre.seg & arm.pj.mcs.df$seg_name %in% centre.seg.name,]$MPhase_joined_seg_width
# 				print(joined.seg.size.vec)
# 				names(joined.seg.size.vec) <- arm.pj.mcs.df[arm.pj.mcs.df$sample %in% samps.w.centre.seg & arm.pj.mcs.df$seg_name %in% centre.seg.name,]$uniq_joined_seg_name
# 				joined.focal.seg.matching.back.to <- names(which.max(joined.seg.size.vec))
# 				print(joined.focal.seg.matching.back.to)
# 				#So now we have the segment we're matching back to, we need its start and end so we can give these segments this id
# 				chr.start.end.vec <- as.numeric(strsplit(strsplit(joined.focal.seg.matching.back.to, split="@")[[1]][2], split = ":")[[1]])
# 				#So now we can make the new ID
# 				new.focal.id <- paste0("focal_", event.type,"_id","@", chr.start.end.vec[1],":",chr.start.end.vec[2],":",chr.start.end.vec[3])
# 				#We can also rejig the matrix to eliminate this from a potential new focal.
# 				#The colons and $ are to help only look for startpositions or end positions as appropriate, otherwise short numbers e.g. 1 or 22 produce errors.
# 				start.col.idx <- grep(paste0(":",chr.start.end.vec[2],":"), colnames(arm.unclassified.focal.mat))
# 				end.col.idx <- grep(paste0(":",chr.start.end.vec[3], "$"), colnames(arm.unclassified.focal.mat))
# 				print(start.col.idx)
# 				print(end.col.idx)
# 				print(colnames(arm.unclassified.focal.mat))
# 				seg.names.to.update <- colnames(arm.unclassified.focal.mat)[start.col.idx:end.col.idx]
# 				arm.unclassified.focal.mat[,start.col.idx:end.col.idx] <- FALSE
# 				#Update the original matrix with the new focal ids
# 				arm.pj.mcs.df[arm.pj.mcs.df$sample%in%samps.w.centre.seg & arm.pj.mcs.df$seg_name%in%seg.names.to.update, col.to.update] <- new.focal.id
# 				f.acr.samp.vec <- colSums(arm.unclassified.focal.mat) == length(unique(arm.pj.mcs.df$sample))
# 				unclassified.cl.focals.remain <- any(f.acr.samp.vec)
#
# 			}#end while (unclassified.focals.remain) {
#
# 				#Update dataframe with this arms results
# 				pj.mcs.df[pj.mcs.df$arm %in% arm, col.to.update] <- arm.pj.mcs.df[, col.to.update]
#
# 		}#end for (arm in arm.vec ) {
#
# 		#Return fully updated dataframe
# 		return(pj.mcs.df)
#
# 	}#end markClonalFocalEvents <- function(potential.cl.foc.mcs.segs.vec, pj.mcs.df, col.to.update, event.type) {
#
# 	#Check if event focal in all regions:
# 	pj.mcs.df$cl_focal_mcs_seg_cpn_vs_ploidy <- "none"
# 	pj.mcs.df$cl_focal_mcs_seg_LOH <- "none"
# 	pj.mcs.df$cl_focal_mcs_seg_LOH_or_loss <- "none"
#
# 	#markClonalFocalEvents <- function(potential.cl.foc.mcs.segs.vec, pj.mcs.df, col.to.update, event.type) {
#
# 	#Focal clonal gains/amp vs ploidy:
# 	#So now we need to come up with any potential events fitting the correct description:
# 	pj.mcs.gain.amp.samp.seg.name.vec <- pj.mcs.df[pj.mcs.df$pj_focal_by_size_and_prop %in% c("focal")
# 									& pj.mcs.df$cpn_event_vs_ploidy %in% c("gain", "amp")
# 									,]$samp_seg_name
# 	pj.mcs.df <- markClonalFocalEvents(potential.cl.foc.mcs.segs.vec = pj.mcs.gain.amp.samp.seg.name.vec,
# 							pj.mcs.df = pj.mcs.df,
# 							col.to.update = "cl_focal_mcs_seg_cpn_vs_ploidy",
# 							event.type = "gain")
# 	#View
# 	pj.mcs.df[pj.mcs.df$cl_focal_mcs_seg_cpn_vs_ploidy !="none",c(colnames(pj.mcs.df)[1:4],"cl_focal_mcs_seg_cpn_vs_ploidy")]
#
# 	#Note we're updating the same column here but it should be fine as doing clonals...
# 	#Should this be done on the minor allele alone really?
# 	#Focal clonal loss/deep_loss vs ploidy:
# 	pj.mcs.loss.deep.loss.samp.seg.name.vec <- pj.mcs.df[pj.mcs.df$pj_focal_by_size_and_prop %in% c("focal")
# 									& pj.mcs.df$cpn_event_vs_ploidy %in% c("loss", "deep_loss")
# 									,]$samp_seg_name
# 	pj.mcs.df <- markClonalFocalEvents(potential.cl.foc.mcs.segs.vec = pj.mcs.loss.deep.loss.samp.seg.name.vec,
# 							pj.mcs.df = pj.mcs.df,
# 							col.to.update = "cl_focal_mcs_seg_cpn_vs_ploidy",
# 							event.type = "loss")
# 	pj.mcs.df[pj.mcs.df$cl_focal_mcs_seg_cpn_vs_ploidy !="none",c(colnames(pj.mcs.df)[1:4],"cl_focal_mcs_seg_cpn_vs_ploidy")]
#
#
# 	#Focal LOH of any kind
# 	#Relevant question here of whether we should actually have a totally different joining thing for LOH?
# 	#Should this be done on the minor allele alone really?
# 	pj.mcs.focal.LOH.samp.seg.name.vec <- pj.mcs.df[pj.mcs.df$pj_focal_by_size_and_prop %in% c("focal")
# 												& pj.mcs.df$within_samp_cl_LOH_present %in% TRUE
# 												& pj.mcs.df$across_samp_cl_LOH %in% TRUE, ]$samp_seg_name
# 	pj.mcs.df <- markClonalFocalEvents(potential.cl.foc.mcs.segs.vec = pj.mcs.focal.LOH.samp.seg.name.vec,
# 							pj.mcs.df = pj.mcs.df,
# 							col.to.update = "cl_focal_mcs_seg_LOH",
# 							event.type = "LOH")
# 	pj.mcs.df[pj.mcs.df$cl_focal_mcs_seg_LOH !="none",c(colnames(pj.mcs.df)[1:4],"cl_focal_mcs_seg_LOH")]
#
# 	#Focal clonal loss/deep_loss and or LOH vs ploidy vs ploidy:
# 	#In order to keep a consistent definition with everything else this is simply either clonal LOH OR clonal loss, not combinations of loss/gainLOH to get to clonal etc.
# 	pj.mcs.df$cl_focal_mcs_seg_LOH_or_loss <- pj.mcs.df$cl_focal_mcs_seg_LOH
# 	pj.mcs.df[pj.mcs.df$cl_focal_mcs_seg_LOH_or_loss == "none"
# 			 & grepl("loss", pj.mcs.df$cl_focal_mcs_seg_cpn_vs_ploidy),]$cl_focal_mcs_seg_LOH_or_loss <- pj.mcs.df[pj.mcs.df$cl_focal_mcs_seg_LOH_or_loss == "none"
# 																										& grepl("loss", pj.mcs.df$cl_focal_mcs_seg_cpn_vs_ploidy),]$cl_focal_mcs_seg_cpn_vs_ploidy
# 	#First we need to change the IDs of the events so they match each other:
# 	pj.mcs.df$cl_focal_mcs_seg_LOH_or_loss <- gsub("focal_loss_id|focal_LOH_id", "focal_loss_or_LOH_id", pj.mcs.df$cl_focal_mcs_seg_LOH_or_loss)
# 	pj.mcs.df[pj.mcs.df$cl_focal_mcs_seg_LOH_or_loss !="none",c(colnames(pj.mcs.df)[1:4],"cl_focal_mcs_seg_LOH_or_loss")]
#
# 	#Now, we need to recalculate the sizes of these events as they may have got larger if next to one another?
# 	#For now let's just leave them... they probably have a separate origin anyway.
#
# 	#Now let's do the subclonal focals:
# 	#potential.sc.foc.mcs.segs.vec all the samp_seg_names that match the condition of the subclonal focal we're looking for
# 	#all chr/arm post joining minimum consistent segmentation segmentation dataframe
# 	#focal column that will be updated by this function
# 	#What is the type of the event we're looking to capture and the string will be placed in the focal id name
# 	markSubclonalFocalEvents <- function(potential.sc.foc.mcs.segs.vec, pj.mcs.df, col.to.update, event.type) {
#
# 		arm.vec <- unique(pj.mcs.df$arm)
# 		num.samples <- length(unique(pj.mcs.df$sample))
#
# 		for (arm in arm.vec ) {
#
# 			arm.pj.mcs.df <- pj.mcs.df[pj.mcs.df$arm %in% arm,]
# 			arm.pj.mcs.df$old_focal_gain_definition <- arm.pj.mcs.df[, col.to.update]
# 			#So now we need to create the column that tell sus whether we actually can look for runs:
# 			arm.pj.mcs.df$potential_new_focal <- arm.pj.mcs.df$samp_seg_name %in% potential.sc.foc.mcs.segs.vec
#
# 			arm.unclassified.focal.mat <- matrix(FALSE, nrow=num.samples, ncol=length(unique(arm.pj.mcs.df$seg_name)))
# 			colnames(arm.unclassified.focal.mat) <- unique(arm.pj.mcs.df$seg_name)
# 			rownames(arm.unclassified.focal.mat) <- unique(arm.pj.mcs.df$sample)
# 			#arm.unclassified.focal.mat.allele <- matrix(FALSE, nrow=length(unique(arm.pj.mcs.df$sample)), ncol=length(unique(arm.pj.mcs.df$seg_name)))
# 			#colnames(arm.unclassified.focal.mat.allele)  <- colnames(arm.unclassified.focal.mat) <- unique(arm.pj.mcs.df$seg_name)
# 			#rownames(arm.unclassified.focal.mat.allele)  <- rownames(arm.unclassified.focal.mat) <- unique(arm.pj.mcs.df$sample)
# 			for(samp in unique(arm.pj.mcs.df$sample)) {
# 				for(seg in unique(arm.pj.mcs.df$seg_name)){
# 					arm.unclassified.focal.mat[samp, seg] <- arm.pj.mcs.df[arm.pj.mcs.df$sample == samp & arm.pj.mcs.df$seg_name == seg,]$potential_new_focal
# 					#arm.unclassified.focal.mat.allele[samp, seg] <- arm.pj.mcs.df[arm.pj.mcs.df$sample == samp & arm.pj.mcs.df$seg_name == seg,]$dominant_allele
# 				}#end for(seg in unique(arm.pj.mcs.df$seg_name)){
# 			}#end for(samp in unique(arm.pj.mcs.df$sample)) {
#
# 			unclassified.focals.remain <- any(arm.unclassified.focal.mat==TRUE)
# 			while (unclassified.focals.remain) {
# 				print(arm.unclassified.focal.mat)
# 				max.seg.name <- names(which.max(colSums(arm.unclassified.focal.mat)))
# 				samps.w.max.seg <-rownames(arm.unclassified.focal.mat[arm.unclassified.focal.mat[,max.seg.name]==TRUE,max.seg.name, drop = FALSE])
# 				joined.seg.size.vec <- arm.pj.mcs.df[arm.pj.mcs.df$sample %in% samps.w.max.seg & arm.pj.mcs.df$seg_name %in% max.seg.name,]$MPhase_joined_seg_width
# 				print("joined.seg.size.vec")
# 				print(joined.seg.size.vec)
# 				names(joined.seg.size.vec) <- arm.pj.mcs.df[arm.pj.mcs.df$sample %in% samps.w.max.seg & arm.pj.mcs.df$seg_name %in% max.seg.name,]$uniq_joined_seg_name
# 				joined.focal.seg.matching.back.to <- names(which.max(joined.seg.size.vec))
# 				print(joined.seg.size.vec)
# 				#So now we have the segment we're matching back to, we need its start and end so we can give these segments this id
# 				chr.start.end.vec <- as.numeric(strsplit(strsplit(joined.focal.seg.matching.back.to, split="@")[[1]][2], split = ":")[[1]])
# 				#So now we can make the new ID
# 				new.focal.id <- paste0("focal_", event.type,"_id","@", chr.start.end.vec[1],":",chr.start.end.vec[2],":",chr.start.end.vec[3])
# 				#We can also rejig the matrix to eliminate this from a potential new focal.
# 				start.col.idx <- grep(paste0(":",chr.start.end.vec[2],":"), colnames(arm.unclassified.focal.mat))
# 				end.col.idx <- grep(paste0(":",chr.start.end.vec[3], "$"), colnames(arm.unclassified.focal.mat))
# 				seg.names.to.update <- colnames(arm.unclassified.focal.mat)[start.col.idx:end.col.idx]
# 				arm.unclassified.focal.mat[,start.col.idx:end.col.idx] <- FALSE
# 				#Update the original matrix with the new focal ids
# 				arm.pj.mcs.df[arm.pj.mcs.df$sample%in%samps.w.max.seg & arm.pj.mcs.df$seg_name%in%seg.names.to.update, col.to.update] <- new.focal.id
# 				unclassified.focals.remain <- any(arm.unclassified.focal.mat==TRUE)
# 				print(new.focal.id)
# 				print("##########3")
#
# 			}#end while (unclassified.focals.remain) {
#
# 			#Update the dataframe with this arm's results
# 			pj.mcs.df[pj.mcs.df$arm %in% arm, col.to.update] <- arm.pj.mcs.df[, col.to.update]
#
# 		}#end for (arm in arm.vec ) {
#
# 		#return updated dataframe
# 		return(pj.mcs.df)
#
# 	}#end markClonalFocalEvents <- function() {
#
#
#
# 	pj.mcs.df$sc_focal_mcs_seg_cpn_vs_ploidy <- "none"
# 	pj.mcs.df$sc_focal_mcs_seg_LOH <- "none"
# 	pj.mcs.df$sc_focal_mcs_seg_LOH_or_loss <- "none"
#
#
#
# 	###SUBCLONAL FOCAL GAINS:
# 	pj.mcs.sc.gain.amp.samp.seg.name.vec <- pj.mcs.df[pj.mcs.df$pj_focal_by_size_and_prop %in% c("focal")
# 							& pj.mcs.df$cpn_event_vs_ploidy %in% c("gain", "amp")
# 							& pj.mcs.df$cl_focal_mcs_seg_cpn_vs_ploidy %in% c("none")
# 							,]$samp_seg_name
# 	pj.mcs.df <- markSubclonalFocalEvents(potential.sc.foc.mcs.segs.vec = pj.mcs.sc.gain.amp.samp.seg.name.vec,
# 										   pj.mcs.df = pj.mcs.df,
# 										   col.to.update = "sc_focal_mcs_seg_cpn_vs_ploidy",
# 										   event.type = "gain")
#
# 	###SUBCLONAL FOCAL LOSSES:
# 	pj.mcs.sc.loss.deep.loss.samp.seg.name.vec <- pj.mcs.df[pj.mcs.df$pj_focal_by_size_and_prop %in% c("focal")
# 							& pj.mcs.df$cpn_event_vs_ploidy %in% c("loss", "deep_loss")
# 							& pj.mcs.df$cl_focal_mcs_seg_cpn_vs_ploidy %in% c("none")
# 							,]$samp_seg_name
#
# 	pj.mcs.df <- markSubclonalFocalEvents(potential.sc.foc.mcs.segs.vec = pj.mcs.sc.loss.deep.loss.samp.seg.name.vec,
# 										   pj.mcs.df = pj.mcs.df,
# 										   col.to.update = "sc_focal_mcs_seg_cpn_vs_ploidy",
# 										   event.type = "loss")
#
# 	pj.mcs.df[pj.mcs.df$cl_focal_mcs_seg_LOH_or_loss !="none",c(colnames(pj.mcs.df)[1:4],"sc_focal_mcs_seg_cpn_vs_ploidy")]
#
#
# 	#SUBCLONAL FOCAL LOH:
# 	#Note these might need to be defined on the minor allele alone, at the moment there are segments that are introduced which are
# 	#are based on the change in copy number...?
# 	####FOCAL LOH
# 	pj.mcs.sc.LOH.samp.seg.name.vec <- pj.mcs.df[pj.mcs.df$pj_focal_by_size_and_prop %in% c("focal")
# 									& pj.mcs.df$within_samp_cl_LOH_present %in% TRUE
# 									& pj.mcs.df$cl_focal_mcs_seg_LOH %in% "none"
# 									& pj.mcs.df$across_samp_cl_LOH %in% FALSE,]$samp_seg_name
# 	pj.mcs.df <- markSubclonalFocalEvents(potential.sc.foc.mcs.segs.vec = pj.mcs.sc.LOH.samp.seg.name.vec,
# 										   pj.mcs.df = pj.mcs.df,
# 										   col.to.update = "sc_focal_mcs_seg_LOH",
# 										   event.type = "LOH")
# 	pj.mcs.df[pj.mcs.df$sc_focal_mcs_seg_LOH !="none",c(colnames(pj.mcs.df)[1:4],"sc_focal_mcs_seg_LOH")]
#
# 	###SUBCLONALFOCAL LOH OR LOSS
# 	#FOCAL LOSS AND/OR LOH
# 	pj.mcs.sc.loss.or.LOH.samp.seg.name.vec <- pj.mcs.df[pj.mcs.df$pj_focal_by_size_and_prop %in% c("focal")
# 			  						 & (pj.mcs.df$within_samp_cl_LOH_present %in% TRUE & pj.mcs.df$across_samp_cl_LOH %in% FALSE
# 			  						 & pj.mcs.df$cl_focal_mcs_seg_LOH_or_loss %in% "none")
# 			  						 | pj.mcs.df$cpn_event_vs_ploidy %in% c("loss", "deep_loss")
# 			  						,]$samp_seg_name
# 	pj.mcs.df <- markSubclonalFocalEvents(potential.sc.foc.mcs.segs.vec = pj.mcs.sc.LOH.samp.seg.name.vec,
# 										   pj.mcs.df = pj.mcs.df,
# 										   col.to.update = "sc_focal_mcs_seg_LOH_or_loss",
# 										   event.type = "loss_or_LOH")
# 	pj.mcs.df[pj.mcs.df$sc_focal_mcs_seg_LOH_or_loss !="none",c(colnames(pj.mcs.df)[1:4],"sc_focal_mcs_seg_LOH_or_loss")]
#
#
# 	#Now save the
# 	pj.full.joint.segs.w.cl.sc.focal.df <- pj.mcs.df
# 	pj.full.joint.segs.w.cl.sc.focal.df.path <- paste(PLOT_OUTPUT_DIR, patient, "_pj_full_joint_segs_df_w_cl_sc_focal_counts.RDS", sep="")
# 	saveRDS(pj.full.joint.segs.w.cl.sc.focal.df, file = pj.full.joint.segs.w.cl.sc.focal.df.path)
#
# } else if (CALCULATE_FOCAL_EVENT_CLONALITY_FLAG == FALSE ) {
#
#
# 	#Let's fix the current bugs:
# 	pj.full.joint.segs.w.cl.sc.focal.df.path <- paste(PLOT_OUTPUT_DIR, patient, "_pj_full_joint_segs_df_w_cl_sc_focal_counts.RDS", sep="")
# 	pj.full.joint.segs.w.cl.sc.focal.df <- readRDS(file = pj.full.joint.segs.w.cl.sc.focal.df.path)
#
#
# }
#
#
#
# #Create the GISTIC input for clonal and subclonal events that we want:
# if (GENERATE_GISTIC_INPUT_FLAG == TRUE) {
#
# 	source(paste0(current.dir, "GISTIC_input_creation_functions.R"))
# 	#MSAI events will be split out (by double counting) into their own tracks.
# 	#So in order to create meaningful GISTIC input for MPhase results we are going to create four ouputs for each patient
# 	#In order to incorporate MSAI events into the analysis these MSAI events will then be added to patients
# 	#who previously only had neutral cpn at this position.
#
# 	#Gains
# 	#Gain MSAI
#
# 	#Loss
# 	#Loss MSAI
#
# 	#As we wish both losses and LOH to be reflected in these data  rather than simply taking the minimum copy numnber for clonal and subclonal values
# 	#As we will do with the gains we will take the minor allele in call cases and use this for the creation of the GISTIC input, this will pick up both
# 	#losses vs ploidy (as we are always taking the minor allele, not the phased cpnB) and LOH event
#
#
# 	#So these are MCS segments within a patient but, due to the work on adding the segments that were missing they should all have
# 	#Let's add the a column that is going to hold the cpn for GISTIC:
# 	#Let's make copies of pt.summ.df that we can use as the ouput:
# 	#Process the data and get it ready for gnerating the GISTIC input:
# 	pj.full.joint.segs.w.cl.sc.focal.df.path <- paste(PLOT_OUTPUT_DIR, patient, "_pj_full_joint_segs_df_w_cl_sc_focal_counts.RDS", sep="")
# 	pj.full.joint.segs.w.cl.sc.focal.df <- readRDS(file = pj.full.joint.segs.w.cl.sc.focal.df.path)
# 	pj.fjs.df <- pj.full.joint.segs.w.cl.sc.focal.df
# 	#Make sure everything has a raw value. TODO fix this earlier in the script.
# 	pj.fjs.df <- transform(pj.fjs.df, ph_raw_cpnA = ifelse(is.na(ph_raw_cpnA), mean_ph_cpnA, ph_raw_cpnA))
# 	pj.fjs.df <- transform(pj.fjs.df, ph_raw_cpnB = ifelse(is.na(ph_raw_cpnB), mean_ph_cpnB, ph_raw_cpnB))
# 	pj.fjs.df <- transform(pj.fjs.df, ph_raw_cpnB = ifelse(ph_raw_cpnB < 0, 0, ph_raw_cpnB))
# 	pj.fjs.df <- transform(pj.fjs.df, ph_raw_cpnB = ifelse(ph_raw_cpnA < 0, 0, ph_raw_cpnA))
#
# 	pj.fjs.df$seg_cn_for_gistic <- 0.7
#
# 	pj.fjs.df$min_raw_cpn_for_GISTIC_loss <- 1
# 	pj.seg.name.vec <- unique(pj.fjs.df$seg_name)
# 	pj.fjs.df$cl_loss_or_LOH_for_GISTIC <- FALSE
#
# 	for ( pj.seg.name in pj.seg.name.vec ) {
# 		print(pj.seg.name)
# 		seg.df <- pj.fjs.df[pj.fjs.df$seg_name %in% pj.seg.name,]
#
# 		if (all(seg.df$min_int_cpnX %in% 0) | all(seg.df$cpn_event_vs_ploidy %in% c("loss", "deep_loss"))) {
#
# 			pj.fjs.df[pj.fjs.df$seg_name %in% pj.seg.name ,]$cl_loss_or_LOH_for_GISTIC <- TRUE
#
# 		}#end if (all(seg.df$min_int_cpnX %in% 0) | all(seg.df$cpn_event_vs_ploidy %in% c("loss", "deep_loss"))) {
#
# 		pj.fjs.df[pj.fjs.df$seg_name %in% pj.seg.name, ]$min_raw_cpn_for_GISTIC_loss <- min(c(seg.df[seg.df$seg_name %in% pj.seg.name]$ph_raw_cpnA,
# 																						  seg.df[seg.df$seg_name %in% pj.seg.name]$ph_raw_cpnB))
#
# 	}#end for ( pj.seg.name in pj.seg.name.vec ) {
# 	#Fixed values and vectors that we'll use for all the different GISITC event criteria:
# 	neutral.seg.name.vec <- unique(pj.fjs.df[ (! is.na(pj.fjs.df$across_samp_cl_cpn_vs_ploidy_present))
# 									& pj.fjs.df$across_samp_cl_cpn_vs_ploidy_present %in% "neutral",]$seg_name)
# 	neutral.min.raw.cpn.df <- collapseSegsByPatientAndMinMaxMean(all.segs.df = pj.fjs.df[pj.fjs.df$seg_name %in% neutral.seg.name.vec,],
# 										 markers.to.use.as.default = "all_sites",
# 										 cpn.to.use = "min_raw_cpn_for_GISTIC_loss",
# 										 min.max.mean.to.use = "mean")
# 	mean.neutral.minor.cpn.for.GISTIC.loss.or.LOH.val <-  mean(neutral.min.raw.cpn.df$Seg.CN)
#
#
# 	all.seg.name.vec <- unique(pj.fjs.df$seg_name)
# 	all.samp.seg.name.vec <- unique(pj.fjs.df$samp_seg_name)
#
# 	neutral.seg.name.vec <- unique(pj.fjs.df[ (! is.na(pj.fjs.df$across_samp_cl_cpn_vs_ploidy_present))
# 									& pj.fjs.df$across_samp_cl_cpn_vs_ploidy_present %in% "neutral",]$seg_name)
# 	neutral.cpn.df <- collapseSegsByPatientAndMinMaxMean(all.segs.df = pj.fjs.df[pj.fjs.df$seg_name %in% neutral.seg.name.vec,],
# 										 markers.to.use.as.default = "all_sites",
# 										 cpn.to.use = "total_raw_cpn",
# 										 min.max.mean.to.use = "mean")
# 	#Use this value to set the segments we're not interested in for gains:
# 	mean.neutral.cpn.for.gains.val <- mean(neutral.cpn.df$Seg.CN)
#
# 	num.markers.by.seg.vec <- collapseSegsByPatientAndMinMaxMean(all.segs.df = pj.fjs.df,
# 										 markers.to.use.as.default = "all_sites",
# 										 cpn.to.use = "total_raw_cpn",
# 										 min.max.mean.to.use = "mean")$Num_Markers[1:length(pj.seg.name.vec)]
#
# 	#So now we have All the values to fill in each of the following dataframes that will be our output for GISTIC input.
#
#
# 	#CLONAL EVENTS
# 	################################################################3s
# 	#CLONAL GAINS
# 	pt.cl.gains.cn.gistic.summ.df <- pj.fjs.df
# 	pt.cl.gains.cn.gistic.summ.df$is_cl_GISTIC_gain <- FALSE
#
# 	cl.gain.seg.name.vec <- unique(pj.fjs.df[ (! is.na(pj.fjs.df$across_samp_cl_cpn_vs_ploidy_present))
# 									& pj.fjs.df$across_samp_cl_cpn_vs_ploidy_present != "subclonal_events_present"
# 									& (! pj.fjs.df$any_across_samp_MSAI %in% TRUE)
# 									& pj.fjs.df$across_samp_cl_cpn_vs_ploidy_present != "neutral"
# 									& pj.fjs.df$across_samp_cl_cpn_vs_ploidy_present%in% c("gain", "amp"),]$seg_name)#is not equal to subclonal_events_present or neutral or NA
#
# 	#Now update the clonal output segments:
# 	if (length(cl.gain.seg.name.vec) > 0) {
#
# 		max.cpn.cl.gains.df <- collapseSegsByPatientAndMinMaxMean(all.segs.df = pj.fjs.df[pj.fjs.df$seg_name %in% cl.gain.seg.name.vec,],
# 												 markers.to.use.as.default = "all_sites",
# 												 cpn.to.use = "total_raw_cpn",
# 												 min.max.mean.to.use = "max")
#
# 		pt.cl.gains.cn.gistic.summ.df[pt.cl.gains.cn.gistic.summ.df$seg_name %in% cl.gain.seg.name.vec,]$seg_cn_for_gistic <- max.cpn.cl.gains.df$Seg.CN
# 		pt.cl.gains.cn.gistic.summ.df[pt.cl.gains.cn.gistic.summ.df$seg_name %in% cl.gain.seg.name.vec,]$is_cl_GISTIC_gain <- TRUE
#
# 	}#end if (length(cl.gain.seg.name.vec) > 0) {
# 	all.non.cl.gain.seg.name.vec <- all.seg.name.vec[ ! all.seg.name.vec %in% cl.gain.seg.name.vec]
# 	pt.cl.gains.cn.gistic.summ.df[pt.cl.gains.cn.gistic.summ.df$seg_name %in% all.non.cl.gain.seg.name.vec,]$seg_cn_for_gistic <- mean.neutral.cpn.for.gains.val
#
# 	pt.cl.gains.cn.gistic.red.summ.df <- pt.cl.gains.cn.gistic.summ.df[match(pj.seg.name.vec, pt.cl.gains.cn.gistic.summ.df$seg_name),]
# 	pt.cl.gains.gistic.input.df <- data.frame(Sample = rep(patient, nrow(pt.cl.gains.cn.gistic.red.summ.df)),
# 			Chromosome = as.numeric(pt.cl.gains.cn.gistic.red.summ.df$chrom),
# 			Start_Position = pt.cl.gains.cn.gistic.red.summ.df$startpos,
# 			End_Position = pt.cl.gains.cn.gistic.red.summ.df$endpos,
# 			 Num_Markers = rep(0, nrow(pt.cl.gains.cn.gistic.red.summ.df)),#This will be filled in using Granges when all the results are collated to make into the input
# 			 #As we want to make sure that all starts and ends are counted as "Markers" so as not to confuse GISTIC
# 			 Seg.CN = pt.cl.gains.cn.gistic.red.summ.df$seg_cn_for_gistic,
# 			 is_cl_GISTIC_gain = pt.cl.gains.cn.gistic.red.summ.df$is_cl_GISTIC_gain,
# 			 stringsAsFactors = FALSE)
# 	pt.cl.gains.gistic.ouput.df.path <- paste(PLOT_OUTPUT_DIR, patient, "_pj_GISTIC_cl_gain_input.RDS", sep="")
# 	saveRDS(pt.cl.gains.gistic.input.df, file = pt.cl.gains.gistic.ouput.df.path)
# 	#So now we need to convert this to being the input for GISTIC
#
# 	################################################################3s
# 	#CLONAL LOSSES:
# 	#Now let's do the clonal losses/LOH
# 	pt.cl.losses.or.loh.cn.gistic.summ.df <- pj.fjs.df
# 	pt.cl.losses.or.loh.cn.gistic.summ.df$is_cl_GISTIC_loss_or_LOH <- FALSE
#
# 	cl.loss.or.loh.seg.name.vec <- unique(pj.fjs.df[pj.fjs.df$cl_loss_or_LOH_for_GISTIC %in% TRUE,]$seg_name)
# 	if (length(cl.loss.or.loh.seg.name.vec) > 0) {
#
# 		min.cpn.cl.losses.df <- collapseSegsByPatientAndMinMaxMean(all.segs.df = pj.fjs.df[pj.fjs.df$seg_name %in% cl.loss.or.loh.seg.name.vec,],
# 												 markers.to.use.as.default = "all_sites",
# 												 cpn.to.use = "min_raw_cpn_for_GISTIC_loss",
# 												 min.max.mean.to.use = "min")
#
# 		pt.cl.losses.or.loh.cn.gistic.summ.df[pt.cl.losses.or.loh.cn.gistic.summ.df$seg_name %in% cl.loss.or.loh.seg.name.vec,]$seg_cn_for_gistic <- min.cpn.cl.losses.df$Seg.CN
# 		pt.cl.losses.or.loh.cn.gistic.summ.df[pt.cl.losses.or.loh.cn.gistic.summ.df$seg_name %in% cl.loss.or.loh.seg.name.vec,]$is_cl_GISTIC_loss_or_LOH <- TRUE
# 	}
# 	all.non.cl.loss.or.loh.seg.name.vec <- all.seg.name.vec[ ! all.seg.name.vec %in% cl.loss.or.loh.seg.name.vec]
# 	pt.cl.losses.or.loh.cn.gistic.summ.df[pt.cl.losses.or.loh.cn.gistic.summ.df$seg_name %in% all.non.cl.loss.or.loh.seg.name.vec,]$seg_cn_for_gistic <- mean.neutral.minor.cpn.for.GISTIC.loss.or.LOH.val
# 	#So now we need to convert this to being the input for GISTIC
# 	pt.cl.losses.or.loh.cn.gistic.red.summ.df <- pt.cl.losses.or.loh.cn.gistic.summ.df[match(pj.seg.name.vec, pt.cl.losses.or.loh.cn.gistic.summ.df$seg_name),]
# 	pt.cl.losses.or.LOH.gistic.input.df <- data.frame(Sample = rep(patient, nrow(pt.cl.losses.or.loh.cn.gistic.red.summ.df)),
# 			Chromosome = as.numeric(pt.cl.losses.or.loh.cn.gistic.red.summ.df$chrom),
# 			Start_Position = pt.cl.losses.or.loh.cn.gistic.red.summ.df$startpos,
# 			End_Position = pt.cl.losses.or.loh.cn.gistic.red.summ.df$endpos,
# 			 Num_Markers = rep(0, nrow(pt.cl.losses.or.loh.cn.gistic.red.summ.df)),#This will be filled in using Granges when all the results are collated to make into the input
# 			 #As we want to make sure that all starts and ends are counted as "Markers" so as not to confuse GISTIC
# 			 Seg.CN = pt.cl.losses.or.loh.cn.gistic.red.summ.df$seg_cn_for_gistic,
# 			 is_cl_GISTIC_gain = pt.cl.losses.or.loh.cn.gistic.red.summ.df$is_cl_GISTIC_loss_or_LOH,
# 			 stringsAsFactors = FALSE)
# 	pt.cl.losses.or.LOH.gistic.input.df.path <- paste(PLOT_OUTPUT_DIR, patient, "_pj_GISTIC_cl_loss_or_LOH_input.RDS", sep="")
# 	saveRDS(pt.cl.losses.or.LOH.gistic.input.df, file = pt.cl.losses.or.LOH.gistic.input.df.path)
#
# 	#SUBCLONAL EVENTS:
# 	################################################################3s
# 	#SUBCLONAL GAINS
# 	pt.sc.gains.cn.gistic.summ.df <- pj.fjs.df
# 	pt.sc.gains.cn.gistic.summ.df$is_sc_GISTIC_gain <- FALSE
# 	sc.gain.seg.name.vec <-  unique(pj.fjs.df[ (! is.na(pj.fjs.df$across_samp_cl_cpn_vs_ploidy_present))
# 									& (pj.fjs.df$across_samp_cl_cpn_vs_ploidy_present == "subclonal_events_present" | pj.fjs.df$any_across_samp_MSAI %in% TRUE)
# 									& pj.fjs.df$cpn_event_vs_ploidy != "neutral"
# 									& pj.fjs.df$cpn_event_vs_ploidy %in% c("gain", "amp"),]$seg_name)
#
# 	#So now if there are parallel gains we need to take gains coming from one allele out of the sc.gain.seg.name.vec to avoid selecting a the same max cpn
# 	#In both the MSAI output and the sc.gain output for GISTIC.
# 	if (length(sc.gain.seg.name.vec) > 0) {
#
# 		max.cpn.sc.gains.df <- collapseSegsByPatientAndMinMaxMean(all.segs.df = pj.fjs.df[pj.fjs.df$seg_name %in% sc.gain.seg.name.vec,],
# 												 markers.to.use.as.default = "all_sites",
# 												 cpn.to.use = "total_raw_cpn",
# 												 min.max.mean.to.use = "max")
#
# 		max.cpn.sc.gains.seg.name.vec <- paste0(max.cpn.sc.gains.df$Chromosome, ":", max.cpn.sc.gains.df$Start_Position, ":", max.cpn.sc.gains.df$End_Position)
# 		pt.sc.gains.cn.gistic.summ.df[pt.sc.gains.cn.gistic.summ.df$seg_name %in% max.cpn.sc.gains.seg.name.vec,]$seg_cn_for_gistic <- max.cpn.sc.gains.df$Seg.CN
# 		pt.sc.gains.cn.gistic.summ.df[pt.sc.gains.cn.gistic.summ.df$seg_name %in% max.cpn.sc.gains.seg.name.vec,]$is_sc_GISTIC_gain <- TRUE
#
# 	}
# 	all.non.sc.gain.seg.name.vec <- all.seg.name.vec[ ! all.seg.name.vec %in% sc.gain.seg.name.vec]
# 	pt.sc.gains.cn.gistic.summ.df[pt.sc.gains.cn.gistic.summ.df$seg_name %in% all.non.sc.gain.seg.name.vec,]$seg_cn_for_gistic <- mean.neutral.cpn.for.gains.val
#
# 	pt.sc.gains.gistic.input.df <- data.frame(Sample = rep(patient, nrow(pt.cl.losses.or.loh.cn.gistic.red.summ.df)),
# 			Chromosome = as.numeric(pt.cl.losses.or.loh.cn.gistic.red.summ.df$chrom),
# 			Start_Position = pt.cl.losses.or.loh.cn.gistic.red.summ.df$startpos,
# 			End_Position = pt.cl.losses.or.loh.cn.gistic.red.summ.df$endpos,
# 			Num_Markers = rep(0, nrow(pt.cl.losses.or.loh.cn.gistic.red.summ.df)),#This will be filled in using Granges when all the results are collated to make into the input
# 			#As we want to make sure that all starts and ends are counted as "Markers" so as not to confuse GISTIC
# 			Seg.CN = pt.cl.losses.or.loh.cn.gistic.red.summ.df$seg_cn_for_gistic,
# 			is_cl_GISTIC_gain = pt.cl.losses.or.loh.cn.gistic.red.summ.df$is_cl_GISTIC_loss_or_LOH,
# 			stringsAsFactors = FALSE)
# 	pt.sc.gains.gistic.input.df.path <- paste(PLOT_OUTPUT_DIR, patient, "_pj_GISTIC_sc_gain_input.RDS", sep="")
# 	saveRDS(pt.sc.gains.gistic.input.df, file = pt.sc.gains.gistic.input.df.path)
#
#
# 	################################################################3s
# 	#MSAI GAINS
# 	pt.msai.gain.cn.gistic.summ.df <- pj.fjs.df
# 	pt.msai.gain.cn.gistic.summ.df$is_MSAI_GISTIC_gain <- FALSE
# 	msai.gain.samp.seg.name.vec <-  unique(pj.fjs.df[ (! is.na(pj.fjs.df$across_samp_cl_cpn_vs_ploidy_present))
# 									& pj.fjs.df$is_parallel_gain_and_or_amp %in% TRUE
# 									& pj.fjs.df$dominant_allele %in% "cpnB",]$samp_seg_name)#is not equal to subclonal_events_present or neutral or NA
#
# 	#So now let's get number values for the MSAI selected segments:
# 	max.cpn.sc.gain.msai.seg.name.vec <- c()
# 	if ( length(msai.gain.samp.seg.name.vec) > 0) {
#
# 		max.cpn.msai.gains.df <- collapseSegsByPatientAndMinMaxMean(all.segs.df = pj.fjs.df[pj.fjs.df$samp_seg_name %in% msai.gain.samp.seg.name.vec,],
# 												 markers.to.use.as.default = "all_sites",
# 												 cpn.to.use = "total_raw_cpn",
# 												 min.max.mean.to.use = "max")
# 		max.cpn.sc.gain.msai.seg.name.vec <- paste0(max.cpn.msai.gains.df$Chromosome, ":", max.cpn.msai.gains.df$Start_Position, ":", max.cpn.msai.gains.df$End_Position)
# 		pt.msai.gain.cn.gistic.summ.df[pt.msai.gain.cn.gistic.summ.df$seg_name %in% max.cpn.sc.gain.msai.seg.name.vec,]$seg_cn_for_gistic <- max.cpn.msai.gains.df$Seg.CN
# 		pt.msai.gain.cn.gistic.summ.df[pt.msai.gain.cn.gistic.summ.df$seg_name %in% max.cpn.sc.gains.seg.name.vec,]$is_MSAI_GISTIC_gain <- TRUE
#
# 	}
# 	all.non.msai.gain.seg.name.vec <- all.seg.name.vec[ ! all.seg.name.vec %in% max.cpn.sc.gain.msai.seg.name.vec]
# 	pt.msai.gain.cn.gistic.summ.df[pt.msai.gain.cn.gistic.summ.df$seg_name %in% all.non.msai.gain.seg.name.vec,]$seg_cn_for_gistic <- mean.neutral.cpn.for.gains.val
#
# 	pt.msai.gains.gistic.input.df <- data.frame(Sample = rep(patient, nrow(pt.msai.gain.cn.gistic.summ.df)),
# 			Chromosome = as.numeric(pt.msai.gain.cn.gistic.summ.df$chrom),
# 			Start_Position = pt.msai.gain.cn.gistic.summ.df$startpos,
# 			End_Position = pt.msai.gain.cn.gistic.summ.df$endpos,
# 			Num_Markers = rep(0, nrow(pt.msai.gain.cn.gistic.summ.df)),#This will be filled in using Granges when all the results are collated to make into the input
# 			#As we want to make sure that all starts and ends are counted as "Markers" so as not to confuse GISTIC
# 			Seg.CN = pt.msai.gain.cn.gistic.summ.df$seg_cn_for_gistic,
# 			is_cl_GISTIC_gain = pt.msai.gain.cn.gistic.summ.df$is_MSAI_GISTIC_gain,
# 			stringsAsFactors = FALSE)
# 	pt.msai.gains.gistic.input.df.path <- paste(PLOT_OUTPUT_DIR, patient, "_pj_GISTIC_msai_gain_input.RDS", sep="")
# 	saveRDS(pt.msai.gains.gistic.input.df, file = pt.msai.gains.gistic.input.df.path)
#
#
# 	################################################################3s
# 	#SUBCLONAL LOSSES and LOH
# 	#Subclonal loss or LOH
# 	pt.sc.losses.or.loh.cn.gistic.summ.df <- pj.fjs.df
# 	pt.sc.losses.or.loh.cn.gistic.summ.df$is_sc_GISTIC_loss_or_LOH <- FALSE
#
# 	sc.loss.or.loh.seg.name.vec <- unique(pj.fjs.df[ (! is.na(pj.fjs.df$across_samp_cl_cpn_vs_ploidy_present))
# 									& ((pj.fjs.df$across_samp_cl_cpn_vs_ploidy_present == "subclonal_events_present" | pj.fjs.df$any_across_samp_MSAI %in% TRUE)
# 									& pj.fjs.df$cpn_event_vs_ploidy %in% c("loss", "deep_loss", "neutral")) | pj.fjs.df$across_samp_sc_LOH %in% TRUE,]$seg_name)
#
# 	#So now how to pull out the MSAI...? Find segments that have both cpnA and cpnB as the dominant allele?
# 	min.cpn.sc.losses.seg.name.vec <- c()
# 	if( length(sc.loss.or.loh.seg.name.vec) > 0 ) {
#
# 		min.cpn.sc.losses.df <- collapseSegsByPatientAndMinMaxMean(all.segs.df = pj.fjs.df[pj.fjs.df$seg_name %in% sc.loss.or.loh.seg.name.vec,],
# 												 markers.to.use.as.default = "all_sites",
# 												 cpn.to.use = "min_raw_cpn_for_GISTIC_loss",
# 												 min.max.mean.to.use = "min")
#
# 		min.cpn.sc.losses.seg.name.vec <- paste0(min.cpn.sc.losses.df$Chromosome,":", min.cpn.sc.losses.df$Start_Position, ":", min.cpn.sc.losses.df$End_Position)
#
# 		pt.sc.losses.or.loh.cn.gistic.summ.df[pt.sc.losses.or.loh.cn.gistic.summ.df$seg_name %in% min.cpn.sc.losses.seg.name.vec,]$seg_cn_for_gistic <- min.cpn.sc.losses.df$Seg.CN
# 		pt.sc.losses.or.loh.cn.gistic.summ.df[pt.sc.losses.or.loh.cn.gistic.summ.df$seg_name %in% min.cpn.sc.losses.seg.name.vec,]$is_sc_GISTIC_loss_or_LOH <- TRUE
#
#
# 	}
# 	all.non.sc.loss.or.loh.seg.name.vec <- all.seg.name.vec[ ! all.seg.name.vec %in% sc.loss.or.loh.seg.name.vec]
# 	pt.sc.losses.or.loh.cn.gistic.summ.df[pt.sc.losses.or.loh.cn.gistic.summ.df$seg_name %in% all.non.sc.loss.or.loh.seg.name.vec,]$seg_cn_for_gistic <- mean.neutral.minor.cpn.for.GISTIC.loss.or.LOH.val
# 	#So now we need to convert this to being the input for GISTIC
# 	pt.sc.losses.or.loh.cn.gistic.red.summ.df <- pt.sc.losses.or.loh.cn.gistic.summ.df[match(pj.seg.name.vec, pt.sc.losses.or.loh.cn.gistic.summ.df$seg_name),]
# 	pt.sc.losses.or.LOH.gistic.input.df <- data.frame(Sample = rep(patient, nrow(pt.cl.losses.or.loh.cn.gistic.red.summ.df)),
# 			Chromosome = as.numeric(pt.sc.losses.or.loh.cn.gistic.red.summ.df$chrom),
# 			Start_Position = pt.sc.losses.or.loh.cn.gistic.red.summ.df$startpos,
# 			End_Position = pt.sc.losses.or.loh.cn.gistic.red.summ.df$endpos,
# 			 Num_Markers = rep(0, nrow(pt.sc.losses.or.loh.cn.gistic.red.summ.df)),#This will be filled in using Granges when all the results are collated to make into the input
# 			 #As we want to make sure that all starts and ends are counted as "Markers" so as not to confuse GISTIC
# 			 Seg.CN = pt.sc.losses.or.loh.cn.gistic.red.summ.df$seg_cn_for_gistic,
# 			 is_cl_GISTIC_gain = pt.sc.losses.or.loh.cn.gistic.red.summ.df$is_sc_GISTIC_loss_or_LOH,
# 			 stringsAsFactors = FALSE)
# 	pt.sc.losses.or.LOH.gistic.input.df.path <- paste(PLOT_OUTPUT_DIR, patient, "_pj_GISTIC_sc_loss_or_LOH_input.RDS", sep="")
# 	saveRDS(pt.sc.losses.or.LOH.gistic.input.df, file = pt.sc.losses.or.LOH.gistic.input.df.path)
#
#
# 	################################################################3s
# 	#MSAI LOSSES and LOH
# 	#MSAI losses or LOH (will be a repeat of a subset of the sc loss and loh)
# 	pt.MSAI.losses.or.LOH.gistic.summ.df <- pj.fjs.df
# 	pt.MSAI.losses.or.LOH.gistic.summ.df$is_MSAI_loss_or_LOH  <- FALSE
#
# 	msai.loss.or.loh.samp.seg.name.vec <-  unique(pj.fjs.df[ (! is.na(pj.fjs.df$across_samp_cl_cpn_vs_ploidy_present))
# 									& pj.fjs.df$is_parallel_loss_and_or_LOH %in% TRUE
# 									& pj.fjs.df$dominant_allele %in% "cpnB",]$samp_seg_name)#is not equal to subclonal_events_present or neutral or NA
#
#
# 	min.cpn.sc.losses.or.LOH.msai.seg.name.vec <- c()
# 	if (length(msai.loss.or.loh.samp.seg.name.vec) > 0 ) {
#
# 		min.cpn.msai.losses.or.loh.df <- collapseSegsByPatientAndMinMaxMean(all.segs.df = pj.fjs.df[pj.fjs.df$samp_seg_name %in% msai.loss.or.loh.samp.seg.name.vec,],
# 												 markers.to.use.as.default = "all_sites",
# 												 cpn.to.use = "min_raw_cpn_for_GISTIC_loss",
# 												 min.max.mean.to.use = "min")
#
# 		#Get the subclonal and MSAI gains and losses actual segments:
# 		min.cpn.sc.losses.or.LOH.msai.seg.name.vec <- paste0(min.cpn.msai.losses.or.loh.df$Chromosome,":", min.cpn.msai.losses.or.loh.df$Start_Position, ":", min.cpn.msai.losses.or.loh.df$End_Position)
#
# 		#Now update th  output segments:
# 		pt.MSAI.losses.or.LOH.gistic.summ.df[pt.MSAI.losses.or.LOH.gistic.summ.df$seg_name %in% min.cpn.sc.losses.or.LOH.msai.seg.name.vec,]$seg_cn_for_gistic <- min.cpn.msai.losses.or.loh.df$Seg.CN
# 		pt.MSAI.losses.or.LOH.gistic.summ.df[pt.MSAI.losses.or.LOH.gistic.summ.df$seg_name %in% min.cpn.sc.losses.or.LOH.msai.seg.name.vec,]$is_MSAI_loss_or_LOH <- TRUE
#
#
# 	}
# 	all.non.msai.losses.or.loh.seg.name.vec <- all.seg.name.vec[ ! all.seg.name.vec %in% min.cpn.sc.losses.or.LOH.msai.seg.name.vec]
# 	pt.MSAI.losses.or.LOH.gistic.summ.df[pt.MSAI.losses.or.LOH.gistic.summ.df$seg_name %in% all.non.msai.losses.or.loh.seg.name.vec,]$seg_cn_for_gistic <- mean.neutral.minor.cpn.for.GISTIC.loss.or.LOH.val
#
# 	pt.msai.losses.or.loh.gistic.input.df <- data.frame(Sample = rep(patient, nrow(pt.MSAI.losses.or.LOH.gistic.summ.df)),
# 			Chromosome = as.numeric(pt.MSAI.losses.or.LOH.gistic.summ.df$chrom),
# 			Start_Position = pt.MSAI.losses.or.LOH.gistic.summ.df$startpos,
# 			End_Position = pt.MSAI.losses.or.LOH.gistic.summ.df$endpos,
# 			Num_Markers = rep(0, nrow(pt.MSAI.losses.or.LOH.gistic.summ.df)),#This will be filled in using Granges when all the results are collated to make into the input
# 			#As we want to make sure that all starts and ends are counted as "Markers" so as not to confuse GISTIC
# 			Seg.CN = pt.MSAI.losses.or.LOH.gistic.summ.df$seg_cn_for_gistic,
# 			is_cl_GISTIC_gain = pt.MSAI.losses.or.LOH.gistic.summ.df$is_MSAI_loss_or_LOH,
# 			stringsAsFactors = FALSE)
# 	pt.msai.losses.or.loh.gistic.input.df.path <- paste(PLOT_OUTPUT_DIR, patient, "_pj_GISTIC_msai_losses_or_LOH_input.RDS", sep="")
# 	saveRDS(pt.msai.losses.or.loh.gistic.input.df, file = pt.msai.losses.or.loh.gistic.input.df.path)
#
#
# } else if (GENERATE_GISTIC_INPUT_FLAG == FALSE) {
#
#
#
#
# }
#
#
#
# #Do we need this loading in here?
# all.samp.mcs.joined.cpnA.cpnB.df.path <- paste0(PLOT_OUTPUT_DIR, patient, "_all_samp_mcs_joined_cpnA_cpnB_df.RDS")
# all.samp.mcs.joined.cpnA.cpnB.df <- readRDS(file = all.samp.mcs.joined.cpnA.cpnB.df.path)
#
# #cpn.chr.start.end.df <- sorted.all.bp.all.samps.segs.anno.df[sorted.all.bp.all.samps.segs.anno.df$full_sample==sorted.all.bp.all.samps.segs.anno.df[1,]$full_sample,1:3]
# cpn.chr.start.end.df.pj <- all.samp.mcs.joined.cpnA.cpnB.df[, c("chrom", "startpos", "endpos")]
# mock.cpn.matrix.pj <- matrix(data=0, nrow=nrow(cpn.chr.start.end.df.pj), ncol = length(sample.vec))
# colnames(mock.cpn.matrix.pj) <- sample.vec
# min.cons.segs.df.pj <- cbind(cpn.chr.start.end.df.pj, mock.cpn.matrix.pj)
#
# seg.name.vec.pj <- unique(paste(min.cons.segs.df.pj[, 1], min.cons.segs.df.pj[, 2], min.cons.segs.df.pj[, 3], sep = ":"))
#
# rownames(pj.full.joint.segs.w.cl.sc.focal.df) <- rownames(full.joint.segs.df.pj)
#
#
# pj.full.joint.segs.w.cl.sc.focal.df <- addWithinSampleChangeToCpnAndAllelicImbalanceCalls(all.seg.df = pj.full.joint.segs.w.cl.sc.focal.df,
# 																						 sample.vec = sample.vec,
# 																						 seg.name.vec = seg.name.vec.pj)
#
#
# # pj.full.joint.segs.w.cl.sc.focal.wd.sc.df.path <- paste(PLOT_OUTPUT_DIR, patient, "_pj_full_joint_segs_df_w_cl_sc_focal_counts_ws_sc.RDS", sep="")
# # pj.full.joint.segs.w.cl.sc.focal.df <- readRDS(file = pj.full.joint.segs.w.cl.sc.focal.wd.sc.df.path)
#
# pj.full.joint.segs.w.cl.sc.focal.df$within_samp_cpnA_passed <- FALSE
# pj.full.joint.segs.w.cl.sc.focal.df$within_samp_cpnB_passed <- FALSE
#
# #both the low and high possible subclonal cpn states are either both <= integer ploidy or >= integer ploidy
# #get a vector for each:
# less.than.floor.raw.ploidy.vec <- ifelse((pj.full.joint.segs.w.cl.sc.focal.df$psc_cpnB_low <= floor(pj.full.joint.segs.w.cl.sc.focal.df$MPhase_raw_ploidy)
# 										& pj.full.joint.segs.w.cl.sc.focal.df$psc_cpnB_high <= floor(pj.full.joint.segs.w.cl.sc.focal.df$MPhase_raw_ploidy))
# 									   | (pj.full.joint.segs.w.cl.sc.focal.df$psc_cpnB_low <= 2
# 									     & pj.full.joint.segs.w.cl.sc.focal.df$psc_cpnB_high <= 2)  ,TRUE, FALSE)
# less.than.floor.raw.ploidy.vec[less.than.floor.raw.ploidy.vec %in% NA] <- FALSE
#
# greater.than.floor.raw.ploidy.vec <- ifelse(pj.full.joint.segs.w.cl.sc.focal.df$psc_cpnA_low >= floor(pj.full.joint.segs.w.cl.sc.focal.df$MPhase_raw_ploidy)
# 					& pj.full.joint.segs.w.cl.sc.focal.df$psc_cpnA_high >= round(pj.full.joint.segs.w.cl.sc.focal.df$MPhase_raw_ploidy), TRUE, FALSE)
# greater.than.floor.raw.ploidy.vec[greater.than.floor.raw.ploidy.vec %in% NA] <- FALSE
#
#
#
# #So now do sc_cpnA:
# init.sc.cpnA.vec <- pj.full.joint.segs.w.cl.sc.focal.df$within_samp_sc_cpnA_present_using_min_cpn_diff
# init.sc.cpnA.vec[init.sc.cpnA.vec %in% NA] <- FALSE
# init.sc.cpnB.vec <- pj.full.joint.segs.w.cl.sc.focal.df$within_samp_sc_cpnB_present_using_min_cpn_diff
# init.sc.cpnB.vec[init.sc.cpnB.vec %in% NA] <- FALSE
#
# #Let's make a
# less.than.floor.raw.ploidy.vec <- floor(pj.full.joint.segs.w.cl.sc.focal.df$MPhase_raw_ploidy)
# less.than.floor.raw.ploidy.vec[less.than.floor.raw.ploidy.vec <= 2] <- 2.5
# less.than.floor.raw.ploidy.vec <- less.than.floor.raw.ploidy.vec * 2
#
# #we’re also going to restrict the highest possible subclonal cn state we’ll examine to 2x integer ploidy
# less.than.thresh.state.vec <- ifelse(pj.full.joint.segs.w.cl.sc.focal.df$total_int_cpn <= less.than.floor.raw.ploidy.vec, TRUE, FALSE)
# less.than.thresh.state.vec[less.than.thresh.state.vec %in% NA] <- FALSE
#
#
# cpn.event.is.gain.or.amp.vec <- pj.full.joint.segs.w.cl.sc.focal.df$cpn_event_vs_ploidy %in% c("gain","amp")
# cpn.event.is.loss.vec <- pj.full.joint.segs.w.cl.sc.focal.df$cpn_event_vs_ploidy %in% c("loss")
#
#
# pj.full.joint.segs.w.cl.sc.focal.df$within_samp_cpnA_passed <- init.sc.cpnA.vec & greater.than.floor.raw.ploidy.vec & less.than.thresh.state.vec & cpn.event.is.gain.or.amp.vec
# pj.full.joint.segs.w.cl.sc.focal.df$within_samp_cpnB_passed <- init.sc.cpnB.vec & less.than.floor.raw.ploidy.vec & less.than.thresh.state.vec & cpn.event.is.loss.vec
#
#
# #Now let's sort out the possible cytoband that shouldn't be kept in the results:
# #This is he accessible genome from CNVkit, including to try and remove the seemingly erroneous LOH present in at the end of some chromosomes that seems to
# #be related to potentially repetitive regions close to the cetnromeres.
# #TODO this should actually go at the top of the script and run all the way through with this...
#
# access.df <- read.table(access.path, sep="\t", header=F, stringsAsFactors=F)
# access.df$V1 <- gsub("chr","", access.df$V1)
# access.df$V1 <- gsub("X","23", access.df$V1)
# access.df$V1 <- gsub("Y","24", access.df$V1)
#
# #Make a reduced name df to add these to:
# pj.fjs.df <- pj.full.joint.segs.w.cl.sc.focal.df
#
# for (chrom in CHRS_TO_USE_VEC) {
# 	#So we need the arm starts and ends that we're currently using in the segmentation
# 	chr.uniq.arms.vec <- unique(pj.fjs.df[pj.fjs.df$chrom %in% chrom,]$arm)
#
# 	#So we need to the accessible regions of the genome that match this
# 	chr.access.df <- access.df[access.df$V1 %in% chrom,]
#
# 	for (chr.arm in chr.uniq.arms.vec) {
#
# 		#So now we need the min and max positions for each arm:
#
#
#
# 		print(chr.arm)
# 		curr.arm.start <- min(pj.fjs.df[pj.fjs.df$arm %in% chr.arm,]$startpos)
# 		curr.arm.end <- max(pj.fjs.df[pj.fjs.df$arm %in% chr.arm,]$endpos)
# 		print(curr.arm.start)
# 		print(curr.arm.end)
#
# 		#Get the relevant positions from the acessible parts of the genome df .
#
# 		#So now we gotta get the split in
# 		arm.chr.access.df <- chr.access.df[chr.access.df$V2 >= curr.arm.start & chr.access.df$V3 <= curr.arm.end, ]
#
# 		if (nrow(arm.chr.access.df) == 0) {
#
# 			arm.chr.access.df <- unique(rbind(chr.access.df[which.min(abs(chr.access.df$V2 - curr.arm.start)),],
# 					chr.access.df[which.min(abs(chr.access.df$V3 - curr.arm.end)),]))
#
#
# 		}
#
# 		red.arm.start <- min(arm.chr.access.df[arm.chr.access.df$V3 <= curr.arm.end,]$V2)
# 		red.arm.end <- max(arm.chr.access.df[arm.chr.access.df$V3 <= curr.arm.end,]$V3)
# 		print(red.arm.start)
# 		print(red.arm.end)
#
# 		#So now update these values in the original dataframe:
# 		pj.fjs.df[pj.fjs.df$arm %in% chr.arm & pj.fjs.df$endpos == curr.arm.end,]$endpos <- red.arm.end
# 		pj.fjs.df[pj.fjs.df$arm %in% chr.arm & pj.fjs.df$startpos == curr.arm.start,]$startpos <- red.arm.start
#
# 	}#end 	for(chr.arm in chr.uniq.arms.vec) {
#
#
# }#end for (chrom in CHRS_TO_USE_VEC) {
#
# # pj.full.joint.segs.w.cl.sc.focal.df[pj.full.joint.segs.w.cl.sc.focal.df$chrom %in% 1, c("startpos", "endpos", "width")]
# # pj.fjs.df[pj.fjs.df$chrom %in% 1, c("startpos", "endpos", "width")]
# pj.fjs.df$seg_width <- as.numeric(pj.fjs.df$endpos) - as.numeric(pj.fjs.df$startpos)
#
# #Now we should be able to drop some rows that meet the following conditions:
# #If the end is less than the start:
# pj.full.joint.segs.w.cl.sc.focal.red.df <- pj.fjs.df[as.numeric(pj.fjs.df$endpos) > as.numeric(pj.fjs.df$startpos),]
# pj.full.joint.segs.w.cl.sc.focal.wd.sc.red.df.path <- paste(PLOT_OUTPUT_DIR, patient, "_pj_full_joint_segs_df_w_cl_sc_focal_counts_ws_sc_red.RDS", sep="")
# saveRDS(pj.full.joint.segs.w.cl.sc.focal.red.df, file = pj.full.joint.segs.w.cl.sc.focal.wd.sc.red.df.path)
#
# #source(paste0(current.dir, "MEDICC_input_creation_functions.R"))
#
# #For each sample let's calculat the samples WGII
# post.mphase.wgii.list <- list()
# for (samp in sample.vec) {
#
# 	#Let's subset the dataframe so we get the correct input for the dataframe
#
# 	samp.pj.fjs.df <- pj.full.joint.segs.w.cl.sc.focal.red.df[pj.full.joint.segs.w.cl.sc.focal.red.df$sample %in% samp,]
# 	samp.formatted.seg.df <- samp.pj.fjs.df[,c("sample","chrom", "startpos", "endpos", "n.het", "cnTotal", "nMajor","nMinor", "Ploidy", "ACF", "nAraw", "nBraw")]
# 	samp.formatted.seg.df <- samp.formatted.seg.df[ ! is.na(samp.formatted.seg.df$nMajor),]
# 	print(samp.formatted.seg.df)
# 	samp.wgii <- calculateWGII(seg.df = samp.formatted.seg.df,)
#
# 	post.mphase.wgii.list[[samp]] <- samp.wgii
#
# 	#pj.full.joint.segs.w.cl.sc.focal.red.df[pj.full.joint.segs.w.cl.sc.focal.red.df$sample %in% samp ,]$mphase_WGII <- samp.wgii
#
# }
#
# #Ok so now we need a version of the wgii function that does the exact same thing that only keeps the gains and losses versus ploidy that we have
#
# #aggregate the size of segments by the copy number states on a per chromosome basis
#
# #Reformat the dataframe so that gains and amps can be handled by aggregate
# pt.df.for.cpn.vs.ploidy.wgii.df <- pj.full.joint.segs.w.cl.sc.focal.red.df
# pt.df.for.cpn.vs.ploidy.wgii.df$cpn_event_vs_ploidy <- gsub("amp", "gain", pt.df.for.cpn.vs.ploidy.wgii.df$cpn_event_vs_ploidy)
# pt.df.for.cpn.vs.ploidy.wgii.df$cpn_event_vs_ploidy <- gsub("deep_loss", "loss", pt.df.for.cpn.vs.ploidy.wgii.df$cpn_event_vs_ploidy)
# tot.chrom.size.df <- aggregate(seg_width ~ chrom + sample, pt.df.for.cpn.vs.ploidy.wgii.df, sum)
# rownames(tot.chrom.size.df) <- paste0(tot.chrom.size.df$sample, "_", tot.chrom.size.df$chrom)
# empty.all.pt.df <- expand.grid(chrom=1:22, sample=sample.vec)
#
#
#
# getWeightedProportionOfGenome <- function(segs.subset.to.interest.df, empty.all.pt.summ.df, all.pt.tot.chrom.size.df) {
#
#
# 	summed.segs.by.interest.df <- aggregate(seg_width ~ chrom + sample, segs.subset.to.interest.df, sum)
# 	all.pt.summed.segs.by.interest.by.chrom.df <- merge(empty.all.pt.summ.df, summed.segs.by.interest.df, by=c("chrom", "sample"), all=TRUE)
# 	rownames(all.pt.summed.segs.by.interest.by.chrom.df) <- paste0(all.pt.summed.segs.by.interest.by.chrom.df$sample, "_", all.pt.summed.segs.by.interest.by.chrom.df$chrom)
# 	all.pt.summed.segs.by.interest.by.chrom.df <- all.pt.summed.segs.by.interest.by.chrom.df[rownames(all.pt.tot.chrom.size.df),]
# 	all.pt.summed.segs.by.interest.by.chrom.df[is.na(all.pt.summed.segs.by.interest.by.chrom.df$seg_width),]$seg_width <- 0
# 	all.pt.summed.segs.by.interest.by.chrom.df$prop_event_width <- 0
#
# 	w.seg.df <- all.pt.summed.segs.by.interest.by.chrom.df
# 	sample.vec <- unique(w.seg.df$sample)
#
# 	for (samp in sample.vec) {
#
# 		w.seg.df[w.seg.df$sample == samp,]$prop_event_width <- as.numeric(w.seg.df[w.seg.df$sample == samp,]$seg_width)/(all.pt.tot.chrom.size.df[all.pt.tot.chrom.size.df$sample == samp,]$seg_width) * (as.numeric(all.pt.tot.chrom.size.df[all.pt.tot.chrom.size.df$sample == samp,]$seg_width)/
# 																												sum(as.numeric(all.pt.tot.chrom.size.df[all.pt.tot.chrom.size.df$sample == samp,]$seg_width)))
#
# 	}#end for (samp in sample.vec) {
#
# 	pt.mphase.wgii.df <- aggregate(prop_event_width ~ sample, w.seg.df, sum)
#
# 	return(pt.mphase.wgii.df)
#
# }
#
# #Weighted MPhase GII
# cpn.event.by.samp.chrom.df <- aggregate(seg_width ~ chrom + cpn_event_vs_ploidy + sample, pt.df.for.cpn.vs.ploidy.wgii.df, sum)
# #So now we want to get only the sectinos of the genome that have events
# mphase.wgii.cpn.event.sep.loss.gain.df <- cpn.event.by.samp.chrom.df[cpn.event.by.samp.chrom.df$cpn_event_vs_ploidy %in% c("gain", "loss"),]
# all.pt.mphase.w.GII.df <- getWeightedProportionOfGenome(segs.subset.to.interest.df = mphase.wgii.cpn.event.sep.loss.gain.df,
# 							empty.all.pt.summ.df = empty.all.pt.df,
# 							all.pt.tot.chrom.size.df = tot.chrom.size.df)
#
# #Weighted MPhase LOH
# LOH.event.by.samp.chrom.df <- aggregate(seg_width ~ chrom + within_samp_cl_LOH_present + sample, pt.df.for.cpn.vs.ploidy.wgii.df, sum)
# only.LOH.event.by.samp.chrom.df <- LOH.event.by.samp.chrom.df[LOH.event.by.samp.chrom.df$within_samp_cl_LOH_present %in% TRUE,]
# all.pt.mphase.w.LOH.df <- getWeightedProportionOfGenome(segs.subset.to.interest.df = only.LOH.event.by.samp.chrom.df,
# 							empty.all.pt.summ.df = empty.all.pt.df,
# 							all.pt.tot.chrom.size.df = tot.chrom.size.df)
#
# #Weighted Loss alone
# ag.loss.event.by.samp.chrom.df <- aggregate(seg_width ~ chrom + cpn_event_vs_ploidy + sample, pt.df.for.cpn.vs.ploidy.wgii.df, sum)
# loss.event.by.samp.chrom.df <- ag.loss.event.by.samp.chrom.df[ag.loss.event.by.samp.chrom.df$cpn_event_vs_ploidy %in% "loss",]
# all.pt.mphase.w.loss.df <- getWeightedProportionOfGenome(segs.subset.to.interest.df = loss.event.by.samp.chrom.df,
# 																empty.all.pt.summ.df = empty.all.pt.df,
# 																all.pt.tot.chrom.size.df = tot.chrom.size.df)
#
# #Weighted MPhase Loss or LOH
# ag.LOH.and.loss.event.by.samp.chrom.df <- aggregate(seg_width ~ chrom + within_samp_cl_LOH_present + cpn_event_vs_ploidy + sample, pt.df.for.cpn.vs.ploidy.wgii.df, sum)
# LOH.and.loss.event.by.samp.chrom.df <- ag.LOH.and.loss.event.by.samp.chrom.df[ag.LOH.and.loss.event.by.samp.chrom.df$within_samp_cl_LOH_present %in% TRUE |
# 																			ag.LOH.and.loss.event.by.samp.chrom.df$cpn_event_vs_ploidy %in% "loss",]
# all.pt.mphase.w.loss.or.LOH.df <- getWeightedProportionOfGenome(segs.subset.to.interest.df = LOH.and.loss.event.by.samp.chrom.df,
# 							empty.all.pt.summ.df = empty.all.pt.df,
# 							all.pt.tot.chrom.size.df = tot.chrom.size.df)
#
# #Weighted MPhase Gain
# ag.gain.event.by.samp.chrom.df <- aggregate(seg_width ~ chrom + cpn_event_vs_ploidy + sample, pt.df.for.cpn.vs.ploidy.wgii.df, sum)
# gain.event.by.samp.chrom.df <- ag.gain.event.by.samp.chrom.df[ag.gain.event.by.samp.chrom.df$cpn_event_vs_ploidy %in% "gain",]
# all.pt.mphase.w.gain.df <- getWeightedProportionOfGenome(segs.subset.to.interest.df = gain.event.by.samp.chrom.df,
# 							empty.all.pt.summ.df = empty.all.pt.df,
# 							all.pt.tot.chrom.size.df = tot.chrom.size.df)
#
# #So we're going to come up with a frew more eve
# #Weighted MPhase GII
# cpn.and.loh.event.by.samp.chrom.df <- aggregate(seg_width ~ chrom + cpn_event_vs_ploidy + sample + within_samp_cl_LOH_present, pt.df.for.cpn.vs.ploidy.wgii.df, sum)
# #So now we want to get only the sectinos of the genome that have events
# mphase.wgii.cpn.event.loss.gain.LOH.df <- cpn.and.loh.event.by.samp.chrom.df[cpn.and.loh.event.by.samp.chrom.df$cpn_event_vs_ploidy %in% c("gain", "loss") |
# 																	cpn.and.loh.event.by.samp.chrom.df$within_samp_cl_LOH_present %in% TRUE,]
# all.pt.mphase.w.LOH.gain.loss.df <- getWeightedProportionOfGenome(segs.subset.to.interest.df = mphase.wgii.cpn.event.loss.gain.LOH.df,
# 							empty.all.pt.summ.df = empty.all.pt.df,
# 							all.pt.tot.chrom.size.df = tot.chrom.size.df)
#
#
#
# pj.full.joint.segs.w.cl.sc.focal.red.df$orig_w_GII <- 0
# pj.full.joint.segs.w.cl.sc.focal.red.df$orig_w_FLOH <- 0
# pj.full.joint.segs.w.cl.sc.focal.red.df$mphase_samp_integer_w_GII <- 0
# pj.full.joint.segs.w.cl.sc.focal.red.df$mphase_samp_integer_w_FLOH <- 0
# pj.full.joint.segs.w.cl.sc.focal.red.df$MPhase_samp_cpn_vs_ploidy_w_GII <- 0
# pj.full.joint.segs.w.cl.sc.focal.red.df$MPhase_samp_w_LOH <- 0
# pj.full.joint.segs.w.cl.sc.focal.red.df$MPhase_samp_w_gain <- 0
# pj.full.joint.segs.w.cl.sc.focal.red.df$MPhase_samp_w_loss <- 0
# pj.full.joint.segs.w.cl.sc.focal.red.df$MPhase_samp_w_loss_or_LOH <- 0
# pj.full.joint.segs.w.cl.sc.focal.red.df$MPhase_w_comb_LOH_vs_ploidy <- 0
#
# for (samp in sample.vec) {
#
# 	#These are calculated using the original WGII calculaton function, written by Nicolai and modified by me.
# 	pj.full.joint.segs.w.cl.sc.focal.red.df[pj.full.joint.segs.w.cl.sc.focal.red.df$sample == samp,]$orig_w_GII <- sample.pre.Mphase.wgii.list[[samp]][,"wGII"]
# 	pj.full.joint.segs.w.cl.sc.focal.red.df[pj.full.joint.segs.w.cl.sc.focal.red.df$sample == samp,]$orig_w_FLOH <- sample.pre.Mphase.wgii.list[[samp]][,"wFLOH"]
# 	pj.full.joint.segs.w.cl.sc.focal.red.df[pj.full.joint.segs.w.cl.sc.focal.red.df$sample == samp,]$mphase_samp_integer_w_GII <- post.mphase.wgii.list[[samp]][,"wGII"]
# 	pj.full.joint.segs.w.cl.sc.focal.red.df[pj.full.joint.segs.w.cl.sc.focal.red.df$sample == samp,]$mphase_samp_integer_w_FLOH <- post.mphase.wgii.list[[samp]][,"wFLOH"]
#
# 	#So now the new weighted wgii measures:
# 	pj.full.joint.segs.w.cl.sc.focal.red.df[pj.full.joint.segs.w.cl.sc.focal.red.df$sample == samp,]$MPhase_samp_cpn_vs_ploidy_w_GII <- all.pt.mphase.w.GII.df[all.pt.mphase.w.GII.df$sample == samp ,]$prop_event_width
# 	pj.full.joint.segs.w.cl.sc.focal.red.df[pj.full.joint.segs.w.cl.sc.focal.red.df$sample == samp,]$MPhase_samp_w_LOH <- all.pt.mphase.w.LOH.df[all.pt.mphase.w.LOH.df$sample == samp ,]$prop_event_width
# 	pj.full.joint.segs.w.cl.sc.focal.red.df[pj.full.joint.segs.w.cl.sc.focal.red.df$sample == samp,]$MPhase_samp_w_gain <-  all.pt.mphase.w.gain.df[all.pt.mphase.w.gain.df$sample == samp ,]$prop_event_width
# 	pj.full.joint.segs.w.cl.sc.focal.red.df[pj.full.joint.segs.w.cl.sc.focal.red.df$sample == samp,]$MPhase_samp_w_loss <-  all.pt.mphase.w.loss.df[all.pt.mphase.w.loss.df$sample == samp ,]$prop_event_width
# 	pj.full.joint.segs.w.cl.sc.focal.red.df[pj.full.joint.segs.w.cl.sc.focal.red.df$sample == samp,]$MPhase_samp_w_loss_or_LOH <-  all.pt.mphase.w.loss.or.LOH.df[all.pt.mphase.w.loss.or.LOH.df$sample == samp ,]$prop_event_width
# 	pj.full.joint.segs.w.cl.sc.focal.red.df[pj.full.joint.segs.w.cl.sc.focal.red.df$sample == samp,]$MPhase_w_comb_LOH_vs_ploidy <-  all.pt.mphase.w.LOH.gain.loss.df[all.pt.mphase.w.LOH.gain.loss.df$sample == samp ,]$prop_event_width
# }
#
# pj.full.joint.segs.w.cl.sc.focal.wd.sc.red.df.wgii.path <- paste(PLOT_OUTPUT_DIR, patient, "_pj_full_joint_segs_df_w_cl_sc_focal_counts_ws_sc_red_wgii.RDS", sep="")
# saveRDS(pj.full.joint.segs.w.cl.sc.focal.red.df, file = pj.full.joint.segs.w.cl.sc.focal.wd.sc.red.df.wgii.path)
#
#
# #So do we want
# print("Making post-MPhase MCS MEDICC Input...")
# post.MPhase.MCS.MEDICC.FASTA.files.list <- makeMEDICCInputListsFromSegmentationDF(segs.to.make.MEDICC.input.df = pj.fjs.df[ ( ! is.na(pj.fjs.df$n.het)) & pj.fjs.df$n.het >= 5 ,],
# 										sample.vec = sample.vec,
# 										chrs.to.use.vec = CHRS_TO_USE_VEC,
# 										chr.col = "chrom",
# 										start.col = "startpos",
# 										end.col = "endpos",
# 										cpnA.col = "ph_int_cpnA",
# 										cpnB.col = "ph_int_cpnB"
# 										)
# post.MPhase.MCS.MEDICC.FASTA.files.list.path <- paste0(PLOT_OUTPUT_DIR, patient, "_post_MPhase_MCS_MEDICC_FASTA_file_list.RDS")
# saveRDS(post.MPhase.MCS.MEDICC.FASTA.files.list, file = post.MPhase.MCS.MEDICC.FASTA.files.list.path)
#
# ########################################################################
# # Only figures are output from here. No saved data files (RDS).
# ########################################################################
#
# #Plot the BAFs again for the final thing
# png(file = paste(PLOT_OUTPUT_DIR, patient, "_BAF_phased", ".png", sep=""), width = 1200, height = 1400)
#
# print("Plotting phased BAF across the genome...")
# plotAllSampleBAFsAcrossGenome(sample.vec, patient.phased.BAF.df.list.pj, cpnA.allele.color = '#de2d2699', cpnB.allele.color= '#3182bd99')
#
# dev.off()
#
#
# #Plot the CCFs for the final thing:
#
# #########OLD FIGURES.........................
# source(paste0(current.dir, "cpn_heterogeneity_plotting_functions_v2.R"))
# ####Plot all genome wide plots
# print("Plotting CCFs for cpnA and cpnB across the genome using CCF definition of subclonality...")
# png(file = paste(PLOT_OUTPUT_DIR, patient, "_CCF_w_cpn_thresh_cpnA_pj", ".png", sep = ""), width = 1200, height = 1400)
#
# plotAllSampleCCFsAcrossGenome(asas.ph.raw.cpn.lists = asas.ph.raw.cpn.lists.pj,
#                               full.joint.segs.df = pj.full.joint.segs.w.cl.sc.focal.df,
#                               sample.vec = sample.vec,
#                               seg.name.vec = seg.name.vec.pj,
#                               all.seg.CCF.vals.list = all.seg.CCF.list.pj,
#                               cpnX = "cpnA",
#                               use.cpn.or.ccf.thresh = "CPN")
#
# dev.off()
#
# png(file = paste(PLOT_OUTPUT_DIR, patient, "_CCF_w_cpn_thresh_cpnB_pj", ".png", sep = ""), width = 1200, height = 1400)
#
# plotAllSampleCCFsAcrossGenome(asas.ph.raw.cpn.lists = asas.ph.raw.cpn.lists.pj,
#                               full.joint.segs.df = pj.full.joint.segs.w.cl.sc.focal.df,
#                               sample.vec = sample.vec,
#                               seg.name.vec = seg.name.vec.pj,
#                               all.seg.CCF.vals.list = all.seg.CCF.list.pj,
#                               cpnX = "cpnB",
#                               use.cpn.or.ccf.thresh = "CPN")
# dev.off()
#
#
# #So now we also want to plot these but only with final threshold held in the passing columns:
#
#
# png(file = paste(PLOT_OUTPUT_DIR, patient, "_CCF_w_cpn_thresh_cpnB_pj", ".png", sep = ""), width = 1200, height = 1400)
#
# plotAllSampleCCFsAcrossGenome(asas.ph.raw.cpn.lists = asas.ph.raw.cpn.lists.pj,
#                               full.joint.segs.df = pj.full.joint.segs.w.cl.sc.focal.df,
#                               sample.vec = sample.vec,
#                               seg.name.vec = seg.name.vec.pj,
#                               all.seg.CCF.vals.list = all.seg.CCF.list.pj,
#                               cpnX = "cpnB",
#                               use.cpn.or.ccf.thresh = "CPN")
# dev.off()
#
# #Plot all the data region by region again;
# #Chosem column names are:
# plottable.col.name.vec <- c("chr", "startpos", "endpos", "int_ph_cpnA", "int_ph_cpnB", "tot_int_cpn", "sample")
#
# #Original segments
# plottable.orig.joint.segs.df <- f.joint.segs.df[, c("chr","startpos", "endpos", "nMajor", "nMinor", "cnTotal", "sample")]
# plottable.orig.joint.segs.df <- plottable.orig.joint.segs.df[plottable.orig.joint.segs.df$chr %in% 1:22,]
# colnames(plottable.orig.joint.segs.df) <- plottable.col.name.vec
# #Full.joint.df segments
# plottable.mphase.mcs.segs.df <- full.joint.segs.df[, c("chr","startpos", "endpos", "ph_int_cpnA", "ph_int_cpnB", "total_int_cpn", "sample")]
# plottable.mphase.mcs.segs.df[is.na(plottable.mphase.mcs.segs.df)] <- -1
# colnames(plottable.mphase.mcs.segs.df) <- plottable.col.name.vec
#
# #Mphase joined segments
# plottable.mphase.samp.cpnA.cpnB.df <- all.samp.joined.cpnA.cpnB.df[, c("chrom","startpos", "endpos", "int_ph_cpnA", "int_ph_cpnB", "int_tot_ph_cpn", "sample")]
# colnames(plottable.mphase.samp.cpnA.cpnB.df) <- plottable.col.name.vec
# print("Plotting sample by sample logr, BAF, phased cpn, and cpn segmentation pre and post joining...")
#
# #MPhase minimum consistent segments:
# plottable.pj.mphase.mcs.df <- pj.full.joint.segs.w.cl.sc.focal.df[,  c("chrom","startpos", "endpos", "ph_raw_cpnA", "ph_raw_cpnB", "int_tot_ph_cpn", "sample")]
# plottable.pj.mphase.mcs.df <- transform(plottable.pj.mphase.mcs.df, ph_raw_cpnB = ifelse(ph_raw_cpnB < 0, 0, ph_raw_cpnB))
# plottable.pj.mphase.mcs.df <- transform(plottable.pj.mphase.mcs.df, ph_raw_cpnA = ifelse(ph_raw_cpnA < 0, 0, ph_raw_cpnA))
# colnames(plottable.pj.mphase.mcs.df) <- c("chr","startpos", "endpos", "ph_raw_cpnA", "ph_raw_cpnB", "int_tot_ph_cpn", "sample")
#
#
# source(paste0(current.dir, "cpn_heterogeneity_plotting_functions_v2.R"))
#
# #We're missing some stuff here:
#
# asas.SNPs.w.both.lr.baf.snp.lists.pj <- readRDS(file = asas.SNPs.w.both.lr.baf.snp.lists.pj.path)
# asas.SNPS.lr.and.or.baf.lists.pj <- readRDS(file = asas.SNPS.lr.and.or.baf.lists.pj.path)
# all.seg.MSAI.by.BAF.ltnpd.present.list.pj <- readRDS(file = all.seg.MSAI.by.BAF.ltnpd.present.list.pj.path)
# patient.phased.BAF.df.list.pj <- readRDS(file = patient.phased.BAF.df.list.pj.path)
# all.seg.CCF.list.pj <- readRDS(file = all.seg.CCF.list.pj.path)
# asas.ph.raw.cpn.lists.pj <- readRDS(file = asas.ph.raw.cpn.lists.pj.path)
#
# for (samp in sample.vec) {
#
# 	png(file = paste(PLOT_OUTPUT_DIR, samp, "_logr_BAF_raw_cpn_and_segments_pj", ".png", sep=""), width = 1200, height = 1400)
#
# 		plotSampleLogRBAFCPNABsAndMSAIAcrossGenome(samp = samp,
# 											full.joint.segs.df = full.joint.segs.df.pj,
# 											orig.segs.df = plottable.orig.joint.segs.df,
# 											mcs.segs.df = plottable.mphase.mcs.segs.df,
# 											new.joined.segs.df = plottable.mphase.samp.cpnA.cpnB.df,
# 											new.mcs.segs.df = plottable.pj.mphase.mcs.df,
# 											asas.ph.raw.cpn.lists = asas.ph.raw.cpn.lists.pj,
# 											all.seg.CCF.vals.list = all.seg.CCF.list.pj,
# 											all.seg.MSAI.by.BAF.ltnpd.present.list = all.seg.MSAI.by.BAF.ltnpd.present.list.pj,
# 											patient.phased.BAF.df.list = patient.phased.BAF.df.list.pj,
# 											asas.SNPs.w.both.lr.baf.snp.lists = asas.SNPs.w.both.lr.baf.snp.lists.pj,
# 											asas.SNPS.lr.and.or.baf.lists = asas.SNPS.lr.and.or.baf.lists.pj,
# 											cpnA.allele.color = '#de2d2699',
# 											cpnB.allele.color= '#3182bd99')
#
# 	dev.off()
#
# }
#
# cat("Finished running mphase.\n")
#
