#!/usr/bin/env Rscript

library(optparse)
library(digest)



if (!interactive()) {

  option_list <- list( make_option(c("--patient"), default="CB1002",
                  help="Patient id"),
      make_option(c("--dip-test"), default=FALSE,
                  help="Use dip test (default=TRUE)"),
      make_option(c("--output-dir"), default="../output/CB1002/mphase",
                  help="Output directory"),
      make_option(c("--medicc-output-dir"), default="../output/CB1002/mphase",
                  help="Output directory"),
      make_option(c("--input-format"), default="ASCAT",
                  help="Input data format (default=ASCAT)"),
      make_option(c("--snp-data"), default="../input/nb-schulte//CB1002_input_paths_df.txt",
                  help="SNP data"),
      make_option(c("--snp-data-prefix"), default="",
                  help="Directory prefix for the paths to the SNP data in the file above"),
      make_option(c("--seg-data"), default="../input/nb-schulte//CB1002.segmentation",
                  help="Segmentation data output from ASCAT"),
      make_option(c("--overwrite"), default=TRUE,
                  help="Overwrite output (default=TRUE)"),
      make_option(c("--gamma"), default=1.0,
                  help="Gamma (for seq=1, for arrays~=0.5 (default=1.0)"),
      make_option(c("-j", "--cores"), default=1,
                  help="CPU cores (default=1)"),
      make_option(c("--segment-type"), default="intersection",
                  help="Segment type, intersection or whole_arm (default=1)"),
      make_option(c("--genome"), default="hg38",
                  help="Human genome version (default=hg19, also supports hg38)"),
      make_option(c("--validate-input"), default=TRUE,
                  help="Validate input files and abort if input is not valid"),
      make_option(c("--seed"), default=2134234,
                  help="Random seed")
      )

  #conf <- parse_args(OptionParser(option_list = option_list, convert_hyphens_to_underscores = TRUE))
  #Hack to get it working on TRACERx 
  if(packageVersion("optparse") != "1.3.2") {

    conf <- parse_args(OptionParser(option_list = option_list), convert_hyphens_to_underscores = TRUE)

  } else {

    conf <- parse_args(OptionParser(option_list = option_list))
    names(conf) <- gsub("\\-", "_", names(conf))
  }
  conf$id <- substr(digest(conf), 1, 8)
  args = commandArgs(trailingOnly=TRUE)

} else {

  source("debugging_conf.R")

}

message("Arguments:")
cat(paste(names(conf), conf, collapse="\n"))

set.seed(conf$seed)#Why is a seed required? I'm not sure anything non-deterministic should be occurring


patient <- conf$patient
#At the moment MPhase will only try to call MSAI and new copy number in areas of the genome identified as allelically imbalanced
#Enable this flag to try and find more allelic imbalance but this will probably break.
allow.dip.test.phasing.flag <- as.logical(conf$dip_test)
output.directory.base.path <- conf$output_dir
data.format.to.use <- conf$input_format #Default for this is MEDICC, but can specify TRACERx_exome, TRACERx_snp_array or TRACERx_WGS

#all.samp.snp.list.path <- conf$cached_data #Points towards the snp files, of whatever format they are given that should be used based on data.format to use.
snp.file.path <- conf$snp_data #Points towards the snp files, of whatever format they are given that should be used based on data.format to use.
seg.file.path <- conf$seg_data #Pintts towards the seg files, of whatever format they should be determined byh the data.format.to.use flag.
overwrite.flag  <- conf$overwrite
gamma.val <- conf$gamma #for sequencing it's 1 and SNP array it's 0.4/0.5 but this is platform specific see ASCAT help pages.
num.cores <- conf$cores #number of cores for the do parallel package functions
intersection.or.whole.arm <- conf$segment_type #Can either be intersection or whole_arm, do we take the only intersection of segments from each sample
genome <- conf$genome  # hg19 or hg38
#or do attempt to fill in the gaps between every segment in each sample so that the entire genome is covered (so output from any sample can run in GISTIC
#as all will cover the same genomic area)

#TBKW edit...
#full.patient <- args[1]
# full.patient<- "A_LTX121"
# seg.file.path <- "/camp/project/proj-tracerx-lung/tracerx/release_tx421_20190411/A_LTX121/exome/MPhase/LTX121_joint_segs_all_ascat_info.txt" #Points towards the seg files, of whatever format they should be determined byh the data.format.to.use flag.
# snp.file.path <- "/camp/project/proj-tracerx-lung/tracerx/release_tx421_20190411/A_LTX121/exome/MPhase/LTX121_TRACERx_input_paths_df.txt"
# patient <- substring(full.patient, 3, 8)
# intersection.or.whole.arm <- "intersection"
# genome <- "hg19"
# allow.dip.test.phasing.flag <- FALSE
# data.format.to.use <- "ASCAT"
# overwrite.flag <- TRUE
# gamma.val <- 1
# num.cores <- 1
# output.directory.base.path <- paste0("/camp/project/proj-tracerx-lung/tctProjects/watkint/tx421_copy_number_qc/xphase_runs/", full.patient, "/")
# system(paste0("mkdir ", output.directory.base.path))
#Right we need a free run folder:
free.run.dir <- file.path(output.directory.base.path, "free_run/")
qc.choice.dir <- file.path(output.directory.base.path, "run_qc_choice/")
rescue.run.dir <- file.path(output.directory.base.path, "rescue_runs/")

dir.create(free.run.dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc.choice.dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rescue.run.dir, recursive = TRUE, showWarnings = FALSE)

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
source(paste0(current.dir, "validate_input.R"))


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
#OJO change back!
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

#Load in the data:
#Note that the functions to create the single dataframes that it accepts are in the data_input_and_QC_functions.R
if (OVERWRITE_FLAG == TRUE) {

	init.f.joint.segs.df <- read.table(seg.file.path, stringsAsFactors = FALSE, sep = "\t", header = TRUE)
  
  # This function will stop the script if the input is not valid
  # FIXME: not working yet
  if (conf$validate_input)
    validate_seg_all_samples(init.f.joint.segs.df)

	# > head(init.f.joint.segs.df) CB001
	#             sample chr  startpos    endpos n.het cnTotal nMajor nMinor Ploidy  ACF nAraw nBraw callState
	# 1 CB1001-11-D01-X1   1    138781  16490694   418       1      1      0    1.9 0.88     1     0        NA

	#                   sample chr  startpos    endpos n.het cnTotal nMajor nMinor   Ploidy ACF    nAraw    nBraw callState
	# 1 A_LTX121_SU_T1.R1_LogR   1    762273  26136282 17836       2      2      0 3.174119 0.5 2.044263 0.000000      Loss


	#Legacy, but I wanted to add in the thresholds for gain/loss etc as they were called in the TRACER NEJM 2017 paper
	#Might be useful... it is aved as part of f.joint.segs.df
	f.joint.segs.df <- reformatInitialSegmentationInput(f.joint.segs.df = init.f.joint.segs.df)
	# > head(f.joint.segs.df) CB001
	#             sample chr  startpos    endpos n.het cnTotal nMajor nMinor Ploidy  ACF nAraw nBraw patient sample_name cpn_call arm_or_focal_event arm_event_present within_chr_arm                                uniq_id
	# 1 CB1001-11-D01-X1   1    138781  16490694   418       1      1      0    1.9 0.88     1     0  CB1001   CB1001:11     loss               none              none     not_in_arm     CB1001-11-D01-X1@1:138781:16490694

	# > head(f.joint.segs.df)
	#              sample chr  startpos    endpos n.het cnTotal nMajor nMinor   Ploidy ACF    nAraw    nBraw patient sample_name cpn_call arm_or_focal_event arm_event_present within_chr_arm                                 uniq_id
	# 1 A_LTX121_SU_T1-R1   1    762273  26136282 17836       2      2      0 3.174119 0.5 2.044263 0.000000  LTX121   LTX121:R1     loss               none              none     not_in_arm     A_LTX121_SU_T1-R1@1:762273:26136282
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
    #                  patient_id                                                        samp_ascat_bc_path 
    # CB1001-11-D01-X1     CB1001 ../input/modified-nb-schulte/CB1001-11-D01-X1/ascat_bc_nochr_modified.rds
    #                   patient_id                                                                       samp_ascat_bc_path
    # A_LTX121_SU_T1-R1     LTX121 ../input/20190617_MPhase_example_data/A_LTX121_SU_T1-R1/A_LTX121_SU_T1-R1.ascat.bc.RData

    if (nchar(conf$snp_data_prefix) > 0) {
       pt.snp.file.paths.df$samp_ascat_bc_path <- file.path(conf$snp_data_prefix, pt.snp.file.paths.df$samp_ascat_bc_path)           
    } 
    #                  patient_id                                                           samp_ascat_bc_path
    # CB1001-11-D01-X1     CB1001 ../../input/modified-nb-schulte/CB1001-11-D01-X1/ascat_bc_nochr_modified.rds

    #                   patient_id                                                                          samp_ascat_bc_path
    # A_LTX121_SU_T1-R1     LTX121 ../../input/20190617_MPhase_example_data/A_LTX121_SU_T1-R1/A_LTX121_SU_T1-R1.ascat.bc.RData

    all.samp.SNP.list <- readInASCAT(sample_table = pt.snp.file.paths.df,
                                     segs = f.joint.segs.df,
                                     samples = sample.vec)
    # List of 5
    #  $ CB1001-11-D01-X1:'data.frame':	675964 obs. of  6 variables:
    #   ..$ chrom        : num [1:675964] 1 1 1 1 1 1 1 1 1 1 ...
    #   ..$ pos          : int [1:675964] 13116 13118 13273 13284 13380 138781 138817 138829 139060 139294 ...
    #   ..$ BAF          : num [1:675964] NA NA NA NA NA NA NA NA NA NA ...
    #   ..$ logR         : num [1:675964] 0.573 0.602 0.622 0.77 0.505 ...
    #   ..$ BAFsegmented : num [1:675964] NA NA NA NA NA ...
    #   ..$ logRsegmented: num [1:675964] -0.838 -0.838 -0.838 -0.838 -0.838 ...


    # List of 2
    #  $ A_LTX121_SU_T1-R1:'data.frame':	938946 obs. of  6 variables:
    #   ..$ chrom        : num [1:938946] 1 1 1 1 1 1 1 1 1 1 ...
    #   ..$ pos          : int [1:938946] 762273 762496 762497 777237 777318 777335 777361 777456 783071 783147 ...
    #   ..$ BAF          : num [1:938946] NA NA NA NA NA NA NA NA NA NA ...
    #   ..$ logR         : num [1:938946] -0.0542 -0.2103 -0.2103 -0.4974 -0.3978 ...
    #   ..$ BAFsegmented : num [1:938946] 0.231 0.231 0.231 0.231 0.231 ...
    #   ..$ logRsegmented: num [1:938946] -0.335 -0.335 -0.335 -0.335 -0.335 ...



    #TBKW edit, I'm assuming that these ascat objects will also contain the germline values and am pulling 
    #them out for this "xphase" version that requires the germline BAFs to compare the allelic imbalance in each tumour
    #region segment vs the equivalent segment in the germline. This comparison is done to try and counter the effect
    #of noise and "reference bias" that might result in false allelic imbalance. The correction that is performed to avert this
    #is primitive as of 16/08/2019 and will require rewriting as well as testing. It should be noted the presence of 
    #bias observed by allelic imbalance is from separate work I've performed and eye-balling the data, it needs to be 
    #shown in the data...

    #Making the assumption here that all regions that were analysed had 
    gl.snp.df <- readInGermlineFromASCAT(ascat_bc_rds_path = pt.snp.file.paths.df[1,]$samp_ascat_bc_path,
    									snp_pos_vec = paste0(all.samp.SNP.list[[1]]$chrom
    														,"_"
    														,all.samp.SNP.list[[1]]$pos))

    # > head(gl.snp.df) CB
    #   chrom    pos BAF logR BAFsegmented logRsegmented  snp_pos
    # 1     1  13116  NA    0           NA    -0.8384372  1_13116
    #ALTX
    #   chrom    pos BAF logR BAFsegmented logRsegmented  snp_pos
    # 1     1 762273  NA    0           NA    -0.3346204 1_762273

    #Strange difference here in that the dat coming from tracerx has mostly NAs in the logr column with occasional zeros, mostly acompanied
    #by a BAF values but not always:
    # > min(gl.snp.df$BAF,na.rm=TRUE)
    # [1] 0.300244
    # > max(gl.snp.df$BAF,na.rm=TRUE)
    # [1] 0.7
    #So it seems that SNPs whose BAF falls outside the 0.3 to 0.7 threshold are given a zero value and have the BAF set to zero but those
    #with it wihin have zero and their ogiinal BAF 
    # > nrow(gl.snp.df)
    # [1] 938946

    #Similary the min and max in the CB is 0.3 and 0.7 but every datapoint has a zero set as logr
    # > nrow(gl.snp.df)
    # [1] 675964

	} else {

		print(data.format.to.use)
		stop("Incorrect data.format.to.use input to mphase")

	}

  # Immediately validate input SNP information (BAF, logR and more)
  #validate_snp_input(all.samp.SNP.list, gl.snp.df)

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
	#                                      Chromosome Start_Position End_Position CB1001-11-D01-X1 CB1001-16-D01-X1 CB1001-15-D01-X1 CB1001-17-D01-X1 CB1001-18-D01-X1
	# CB1001-11-D01-X1:1:1:138780                   1              1       138780                0                0                0                0                0
	# CB1001-11-D01-X1:1:138781:16421937            1         138781     16421937                0                0                0                0                0

	# > head(min.cons.segs.df)
	#                                         Chromosome Start_Position End_Position A_LTX121_SU_T1-R1 A_LTX121_SU_T1-R2
	# A_LTX121_SU_T1-R1:1:1:762272                     1              1       762272                 0                 0

	#end of 201908TBKW edit all this code 


	#Minimum consistent segmentation is finished! Let's get the SNPs and their LogR and BAF values amd assign to segments.
	#All SNPs in these lists must have both a logr and baf value: all sample all seg logr baf lists is asas.lr.baf
	ptm <- proc.time()
	message("cleanAndErrorCheckSNPsParallel first")
	#source(paste0(current.dir, "data_input_and_QC_functions_v2.R"))
	asas.SNPs.w.both.lr.baf.snp.lists <- cleanAndErrorCheckSNPsParallel(all.samp.SNP.list,
																		mcs.seg.to.use.df =min.cons.segs.df,
																		allow.logr.only.flag = FALSE,
																		data.type = data.format.to.use)

	# List of 5
	#  $ CB1001-11-D01-X1:List of 186
	# > str(asas.SNPs.w.both.lr.baf.snp.lists[[1]][[2]], max=2)
	# 'data.frame':	294 obs. of  6 variables:
	#  $ chrom        : num  1 1 1 1 1 1 1 1 1 1 ...
	#  $ pos          : int  970892 973673 974175 981169 1049886 1087765 1091327 1522693 1541864 1704310 ...
	#  $ BAF          : num  1 0 0.1 0.9 0.867 ...
	#  $ logR         : num  -0.526 -1.879 -0.363 -1.287 -1.773 ...
	#  $ BAFsegmented : num  0.107 0.107 0.107 0.107 0.107 ...
	#  $ logRsegmented: num  -0.838 -0.838 -0.838 -0.838 -0.838 ...

	#  List of 2
	#   $ A_LTX121_SU_T1-R1:List of 418

	# >  str(asas.SNPs.w.both.lr.baf.snp.lists[[1]][[2]], max=2)
	# 'data.frame':	755 obs. of  6 variables:
	#  $ chrom        : num  1 1 1 1 1 1 1 1 1 1 ...
	#  $ pos          : int  871020 879317 880390 881627 892460 894573 896064 897730 902108 906272 ...
	#  $ BAF          : num  0.257 0.248 0.2 0.866 0.143 ...
	#  $ logR         : num  -0.271 -0.36 -0.385 -0.208 -0.571 ...
	#  $ BAFsegmented : num  0.231 0.231 0.231 0.231 0.231 ...
	#  $ logRsegmented: num  -0.335 -0.335 -0.335 -0.335 -0.335 ...
	# >
	if (length(unique(sapply(asas.SNPs.w.both.lr.baf.snp.lists, length))) != 1) {

		stop("Not all segments are present in all samples' snp assignations, likely something has gone wrong, is the parallelisation working correctly?")

	}
	saveRDS(asas.SNPs.w.both.lr.baf.snp.lists, file = asas.SNPs.w.both.lr.baf.snp.lists.path)
	time.taken <- proc.time() - ptm
	message("Runtime of cleanAndErrorCheckSNPsParallel:")
	print(time.taken)
	##
	#SNPs that only have LogR are permitted in these lists, but it also contains all the SNPs with both BAF and LogR
	ptm <- proc.time()
	message("cleanAndErrorCheckSNPsParallel second")

	asas.SNPS.lr.and.or.baf.lists <- cleanAndErrorCheckSNPsParallel(all.samp.SNP.list,
																	mcs.seg.to.use.df =min.cons.segs.df,
																	allow.logr.only.flag = TRUE,
																	data.type = data.format.to.use)
	#CB
	# > str(asas.SNPS.lr.and.or.baf.lists[[1]][[2]], max=1)
	# 'data.frame':	7923 obs. of  6 variables:
	#  $ chrom        : num  1 1 1 1 1 1 1 1 1 1 ...
	#  $ pos          : int  138781 138817 138829 139060 139294 732966 732994 733201 733202 733225 ...
	#  $ BAF          : num  NA NA NA NA NA NA NA NA NA NA ...
	#  $ logR         : num  -0.3305 -0.1493 0.073 -0.3676 -0.0592 ...
	#  $ BAFsegmented : num  0.107 0.107 0.107 0.107 0.107 ...
	#  $ logRsegmented: num  -0.838 -0.838 -0.838 -0.838 -0.838 ...
	#ALTX121
	# > str(asas.SNPS.lr.and.or.baf.lists[[1]][[2]], max=1)
	# 'data.frame':	7923 obs. of  6 variables:
	#  $ chrom        : num  1 1 1 1 1 1 1 1 1 1 ...
	#  $ pos          : int  138781 138817 138829 139060 139294 732966 732994 733201 733202 733225 ...
	#  $ BAF          : num  NA NA NA NA NA NA NA NA NA NA ...
	#  $ logR         : num  -0.3305 -0.1493 0.073 -0.3676 -0.0592 ...
	#  $ BAFsegmented : num  0.107 0.107 0.107 0.107 0.107 ...
	#  $ logRsegmented: num  -0.838 -0.838 -0.838 -0.838 -0.838 ...

	if (length(unique(sapply(asas.SNPS.lr.and.or.baf.lists, length))) != 1) {

		stop("Not all segments are present in all samples' snp assignations, likely something has gone wrong, is the parallelisation working correctly?")

	}

	#print(asas.SNPS.lr.and.or.baf.lists.path)
	saveRDS(asas.SNPS.lr.and.or.baf.lists, file = asas.SNPS.lr.and.or.baf.lists.path)
	time.taken <- proc.time() - ptm
	message("Runtime of cleanAndErrorCheckSNPsParallel:")
	print(time.taken)

	if ( length(asas.SNPs.w.both.lr.baf.snp.lists) <= 1 ) {

		print(str(asas.SNPs.w.both.lr.baf.snp.lists, max = 1))
		stop("Only one or fewer samples present after reading in data.")

	}

	#######So here is where we should check if there are any segments that should be
	#save:
	saveRDS(asas.SNPS.lr.and.or.baf.lists, file = asas.SNPS.lr.and.or.baf.lists.path)
	saveRDS(asas.SNPs.w.both.lr.baf.snp.lists, file = asas.SNPs.w.both.lr.baf.snp.lists.path)
	saveRDS(min.cons.segs.df, file = min.cons.segs.df.path)
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
	samp.orig.formatted.seg.df <- samp.orig.seg.df[,c("sample","chr", "startpos", "endpos", "n.het", "cnTotal", "nMajor","nMinor", "Ploidy", "ACF", "nAraw", "nBraw")]
	samp.orig.formatted.seg.df <- samp.orig.formatted.seg.df[ ! is.na(samp.orig.formatted.seg.df$nMajor),]
	pre.mphase.samp.wgii <- calculateWGII(seg.df = samp.orig.formatted.seg.df)
	sample.pre.Mphase.wgii.list[[samp]] <-  pre.mphase.samp.wgii

}#end for (samp in sample.vec) {

#Based on the input and how many BAF values there are we need to set the subclonality thresholds
total.num.of.hz.sites.per.sample <- sum(sapply(asas.SNPs.w.both.lr.baf.snp.lists[[1]], nrow))
message("Setting the pval thresholds for calling subclonality depending on the number of heterozygous sites assayed:")

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

	#source(paste0(current.dir, "consensus_phasing_from_multiple_samples_functions_v2.R"))
	# all.seg.phasing.info.list <- getCombinedBAFDipTestPhasing(asas.lr.baf.lists = asas.SNPs.w.both.lr.baf.snp.lists,
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
	message("Getting phasing to apply to all samples using getCombinedBAFDipTestPhasingParallel...")
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

	# #SO note both CB and LTX have the following warnings:
	# 29: In regularize.values(x, y, ties, missing(ties)) :
 #  	collapsing to unique 'x' values
 #  	#However, LTX has many more of them with fewer regions? Why?

 	#Also both seem to be broken as there is never anything returned by samp_most_significant_by_diptest in either.

	saveRDS(all.seg.phasing.info.list, file = all.seg.phasing.info.list.path)


	#Todo unsure of what to do with the baf logr lists which have positions that only have logr here, at the moment exclude them but how
	#does that affect later run
	message("Generating phased BAF dataframes...")
	patient.phased.BAF.df.list <- createAllSegsPhasedBAFDataFrame(asas.lr.baf.lists = asas.SNPs.w.both.lr.baf.snp.lists,
                                                         sample.to.phase.seg.list = all.seg.phasing.info.list,
                                                         min.num.hz.snp.phasing = MIN_NUM_HZ_SNP_PHASING
                                                         )
	if(length(patient.phased.BAF.df.list) <= 1) {

		stop("Error in generating phased BAF dataframes")
	}
	saveRDS(patient.phased.BAF.df.list, file = patient.phased.BAF.df.list.path)

	#Now we have the phasing, we want to check whether, using it, we think that there is sufficient separation
	#of the bands of allelic balance to consider a segment "allelically imbalanced"
	#If not set the bafsegmented to 0.5 
	#So we may want to reintroduce the log truncated normal stuff here...
	


	#So this is where want to create a phased ascat.bc object:
	#TBKW NEW WORK:
	message("Starting new TBKW edit... begun on 16/08/2019 to recall copy number using phasing.")
	source(paste0(current.dir, "ascat.R"))
	#So determine which SNPs are hz:
	snp.file.df <- gl.snp.df



	snp.file.df$snp_pos <- paste0(snp.file.df$chrom, "_", snp.file.df$pos)
	snp.file.df$snp_idx <- 1:nrow(snp.file.df)

	# #CB
	# > head(snp.file.df)
	#   chrom    pos BAF logR BAFsegmented logRsegmented  snp_pos snp_idx
	# 1     1  13116  NA    0           NA    -0.8384372  1_13116       1
	# 2     1  13118  NA    0           NA    -0.8384372  1_13118       2
	# 3     1  13273  NA    0           NA    -0.8384372  1_13273       3
	# 4     1  13284  NA    0           NA    -0.8384372  1_13284       4
	# 5     1  13380  NA    0           NA    -0.8384372  1_13380       5
	# 6     1 138781  NA    0           NA    -0.8384372 1_138781       6

	# #LTX121
	# > head(snp.file.df)
	#   chrom    pos BAF logR BAFsegmented logRsegmented  snp_pos
	# 1     1 762273  NA    0           NA    -0.3346204 1_762273
	# 2     1 762496  NA   NA           NA    -0.3346204 1_762496
	# 3     1 762497  NA   NA           NA    -0.3346204 1_762497
	# 4     1 777237  NA   NA           NA    -0.3346204 1_777237
	# 5     1 777318  NA   NA           NA    -0.3346204 1_777318
	# 6     1 777335  NA   NA           NA    -0.3346204 1_777335
	#Differences can be seen here and are noted when the gl.snp.df has been read in

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
	#Note anno seg is sorted.all.bp.all.samps.segs.anno.df
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

	# > head(int.anno.seg.df[,ncol(int.anno.seg.df):(ncol(int.anno.seg.df)-12)])
	#                                      end_snp_idx start_snp_idx end_is_pos start_is_pos end_is_logr end_is_homo end_is_hz    end_pos start_is_logr start_is_homo start_is_hz  start_pos chr_arm
	# CB1001-11-D01-X1:1:138781:16421937          7928             6       TRUE         TRUE       FALSE       FALSE      TRUE 1_16421937          TRUE         FALSE       FALSE   1_138781      1p
	# CB1001-11-D01-X1:1:16421938:16490693          NA            NA      FALSE        FALSE       FALSE       FALSE     FALSE 1_16490693         FALSE         FALSE       FALSE 1_16421938      1p
	# CB1001-11-D01-X1:1:16490694:16535960          NA          7973      FALSE         TRUE       FALSE       FALSE     FALSE 1_16535960         FALSE         FALSE        TRUE 1_16490694      1p
	# CB1001-11-D01-X1:1:16535961:16948870        8287          7979       TRUE         TRUE        TRUE       FALSE     FALSE 1_16948870          TRUE         FALSE       FALSE 1_16535961      1p
	# CB1001-11-D01-X1:1:16948871:17232700          NA            NA      FALSE        FALSE       FALSE       FALSE     FALSE 1_17232700         FALSE         FALSE       FALSE 1_16948871      1p
	# CB1001-11-D01-X1:1:17232701:17243971        8487          8477       TRUE         TRUE       FALSE       FALSE      TRUE 1_17243971         FALSE         FALSE        TRUE 1_17232701      1p

	# > head(int.anno.seg.df[,ncol(int.anno.seg.df):(ncol(int.anno.seg.df)-12)])
	#                                         end_snp_idx start_snp_idx end_is_pos start_is_pos end_is_logr end_is_homo end_is_hz     end_pos start_is_logr start_is_homo start_is_hz   start_pos chr_arm
	# A_LTX121_SU_T1-R1:1:762273:26136282           17836             1       TRUE         TRUE        TRUE       FALSE     FALSE  1_26136282          TRUE         FALSE       FALSE    1_762273      1p
	# A_LTX121_SU_T1-R1:1:26136283:26136297            NA            NA      FALSE        FALSE       FALSE       FALSE     FALSE  1_26136297         FALSE         FALSE       FALSE  1_26136283      1p
	# A_LTX121_SU_T1-R1:1:26136298:121310061        50594         17837       TRUE         TRUE        TRUE       FALSE     FALSE 1_121310061          TRUE         FALSE       FALSE  1_26136298      1p
	# A_LTX121_SU_T1-R1:1:144052752:145096660       50645         50595       TRUE         TRUE        TRUE       FALSE     FALSE 1_145096660          TRUE         FALSE       FALSE 1_144052752      1q
	# A_LTX121_SU_T1-R1:1:145096661:145109582          NA            NA      FALSE        FALSE       FALSE       FALSE     FALSE 1_145109582         FALSE         FALSE       FALSE 1_145096661      1q
	# A_LTX121_SU_T1-R1:1:145109583:145477446       50734         50646       TRUE         TRUE        TRUE       FALSE     FALSE 1_145477446         FALSE         FALSE        TRUE 1_145109583      1q
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

	# #IN CB there are some SNPs that do not have an arm assiged:
	# >  length(unique(snp.file.df$arm))
	# [1] 41
	# > unique(snp.file.df$arm
	# + )
	#  [1] "1p"  "1q"  "2p"  "2q"  "3p"  "3q"  "4p"  "4q"  "5p"  "5q"  
	#  "6p"  "6q"  "7p"  "7q"  "8p"  "8q"  "9p"  "9q"  "10p" "10q" "11p" 
	#  "11q" "12p" "12q" "13q" "14q" "15q" "16p" "16q" "17p" "17q" "18p" 
	#  "18q" "19p" "19q" "20p" "20q" "21p" "21q" "22q" NA

	#vs LTX121

	# > length(unique(snp.file.df$arm))
	# [1] 40
	# > unique(snp.file.df$arm)
	#  [1] "1p"  "1q"  "2p"  "2q"  "3p"  "3q"  "4p"  "4q"  "5p"  "5q"
	#    "6p"  "6q"  "7p"  "7q"  "8p"  "8q"  "9p"  "9q"  "10p" "10q" "11p"
	#     "11q" "12p" "12q" "13q" "14q" "15q" "16p" "16q" "17p" "17q" "18p"
	#      "18q" "19p" "19q" "20p" "20q" "21p" "21q" "22q"

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
	#CB
	# > head(all.snp.grs)
	# GRanges object with 6 ranges and 0 metadata columns:
	#     seqnames    ranges strand
	#        <Rle> <IRanges>  <Rle>
	#   1        1     13116      *
	#   2        1     13118      *
	#   3        1     13273      *
	#   4        1     13284      *
	#   5        1     13380      *
	#   7        1    138817      *
	#   -------
	#LTX121
	# > head(all.snp.grs)
	# GRanges object with 6 ranges and 0 metadata columns:
	#     seqnames    ranges strand
	#        <Rle> <IRanges>  <Rle>
	#   2        1    762496      *
	#   3        1    762497      *
	#   4        1    777237      *
	#   5        1    777318      *
	#   6        1    777335      *
	#   7        1    777361      *

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

	# > head(all.int.grs)
	# GRanges object with 6 ranges and 0 metadata columns:
	#                                        seqnames            ranges strand
	#                                           <Rle>         <IRanges>  <Rle>
	#     CB1001-11-D01-X1:1:138781:16421937        1   138781-16421937      *
	#   CB1001-11-D01-X1:1:16421938:16490693        1 16421938-16490693      *
	#   CB1001-11-D01-X1:1:16490694:16535960        1 16490694-16535960      *
	#   CB1001-11-D01-X1:1:16535961:16948870        1 16535961-16948870      *
	#   CB1001-11-D01-X1:1:16948871:17232700        1 16948871-17232700      *
	#   CB1001-11-D01-X1:1:17232701:17243971        1 17232701-17243971      *

	# > head(all.int.grs)
	# GRanges object with 6 ranges and 0 metadata columns:
	#                                           seqnames              ranges strand
	#                                              <Rle>           <IRanges>  <Rle>
	#       A_LTX121_SU_T1-R1:1:762273:26136282        1     762273-26136282      *
	#     A_LTX121_SU_T1-R1:1:26136283:26136297        1   26136283-26136297      *
	#    A_LTX121_SU_T1-R1:1:26136298:121310061        1  26136298-121310061      *
	#   A_LTX121_SU_T1-R1:1:144052752:145096660        1 144052752-145096660      *
	#   A_LTX121_SU_T1-R1:1:145096661:145109582        1 145096661-145109582      *
	#   A_LTX121_SU_T1-R1:1:145109583:145477446        1 145109583-145477446      *
	#
	lacking.snp.idx.vec <- which(int.anno.seg.df$start_is_pos %in% FALSE | int.anno.seg.df$end_is_pos %in% FALSE)
	#FIXME, potential improvement that could be written in would be to, once a SNP is being used
	#remove it from the pool of potential snps to be used as a potenial start and end. However, this might not
	#matter to much.
	int.anno.seg.df$nearest_preceding_snp <- NA
	int.anno.seg.df$nearest_following_snp <- NA 
	message("Lacking snp idx:")
	for (i in 1:length(all.int.grs)) {
		#Check the seqnames part as SNPs really need to be in the same need to be in the same chromosome arm...?
		if( ! i %in% lacking.snp.idx.vec) {
			next
		}
    message(" ", i, appendLF = FALSE)

		curr.arm <- int.anno.seg.df[i,]$chr_arm
    i_grs <- all.int.grs[i]
    #Precde is wrong here follow is correct, check the granges documentation.
    #before_i <- precede(i_grs, all.snp.grs)
    before_i <- follow(i_grs, all.snp.grs)
		if (! is.na(before_i) ) {
			if(snp.file.df[before_i,]$arm %in% curr.arm) {
				int.anno.seg.df[i,]$nearest_preceding_snp <- start(all.snp.grs[before_i])

			}
		}

    #follow is wrong here precede is correct, check the granges documentation.
    after_i <- precede(i_grs, all.snp.grs)
		if (! is.na(after_i) ) {
			if( snp.file.df[after_i,]$arm %in% curr.arm) {
				int.anno.seg.df[i,]$nearest_following_snp <- start(all.snp.grs[after_i])
			}
		}

		#Any remaining examples are going to be the first or only 
		if(is.na(int.anno.seg.df[i,]$nearest_following_snp) & is.na(int.anno.seg.df[i,]$nearest_preceding_snp)) {
			int.anno.seg.df[i,]$nearest_following_snp <- start(all.snp.grs[nearest(i_grs, all.snp.grs)])
			int.anno.seg.df[i,]$nearest_preceding_snp <- start(all.snp.grs[nearest(i_grs, all.snp.grs)])
		}
	}

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

	test <- f.int.anno.seg.df[(! is.na(f.int.anno.seg.df$nearest_preceding_snp)) & ( ! is.na(f.int.anno.seg.df$nearest_following_snp)),]
	nrow(test)
	nrow(test[test$nearest_following_snp > test$nearest_preceding_snp,])
	nrow(test[test$nearest_following_snp == test$nearest_preceding_snp,])
	nrow(test[test$nearest_following_snp < test$nearest_preceding_snp,])#THis is a bug!
	#Note there could be an issue here if not every segment contains a BAF point...
	#But, if we're using hte intersection this should be fine? Not sure on this point.
	#However it doesn't really matter as ASCAT makes sure every segment whether or not it had any baf point has an allele-specific value.
	message("Running ASCAT on all samples...")

	#So now there are several things we want to do here depending on the options passed to the tool
	#1) Using the phasing information obtained from all samples, allow a free run allowing ASCAT to call whatever purity and ploidy of all samples included
	#in the run it feels is most appropriate. This will potentially provide new and improved solutions that might be better than what was fed into mphase
	#particularly if the segmentation of the input did not make use of multi-region phasing data.
	#2)  
	#Initiate list for segmentation
	ph.ascat.output.list <- list()
	ph.ascat.segs.df.list <- list()
	#/fast/work/groups/ag_schwarz/Projects/project-refphase/tom/tools/multiregioncopynum
	base.dir <- getwd()
	for (samp in sample.vec) {
		message("Loading sample: ", samp)
		#Here we load in the file path of the relevant samples ascat.bc object. For compatibility with Matt's code
		#This should accept the ascat.bc files in eithe RDS or RData format.

    if (grepl("RData", pt.snp.file.paths.df[samp,]$samp_ascat_bc_path)) {
		#OJO: fix this warning
		# [1] "CB1001-11-D01-X1"
		# Error in gzfile(file, "rb") :
		#   (converted from warning) cannot open compressed file '../../input/modified-nb-schulte/CB1001-11-D01-X1/ascat_bc_nochr_modified.rds', probable reason 'No such file or directory'
		# > options(warn=1)
		# > traceback()
		# 6: doWithOneRestart(return(expr), restart)
		# 5: withOneRestart(expr, restarts[[1L]])
		# 4: withRestarts({
		#        .Internal(.signalCondition(simpleWarning(msg, call), msg,
		#            call))
		#        .Internal(.dfltWarn(msg, call))
		#    }, muffleWarning = function() NULL)
		# 3: .signalSimpleWarning("cannot open compressed file '../../input/modified-nb-schulte/CB1001-11-D01-X1/ascat_bc_nochr_modified.rds', probable reason 'No such file or directory'",
		#        base::quote(gzfile(file, "rb")))
		# 2: gzfile(file, "rb")
		# 1: readRDS(pt.snp.file.paths.df[samp, ]$samp_ascat_bc_path)

		  message("ascat.bc object in RData file, loading...")
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

		if (any(is.na(tum.logr.vec))) {
			message("There are NA values in the ascat.bc Tumor_LogR field which need to be removed. Revisit generation of this input. For changing the NAs to 0.")
			na.idx.vec <- which(is.na(tum.logr.vec))
			tum.logr.vec[na.idx.vec] <- 0
		}

		tum.logr.segmented.vec <- ascat.bc[["Tumor_LogR_segmented"]][,1]

		if(any(is.na(tum.logr.segmented.vec))) {
			message("There are NA values in the ascat.bc Tumor_LogR_segmented field which need to be removed. Revisit generation of this input. For changing the NAs to 0.")
			na.segmented.idx.vec <- which(is.na(tum.logr.segmented.vec))
			tum.logr.segmented.vec[na.segmented.idx.vec] <- 0
		}


		for (i in 1:nrow(samp.seg.df)){
		  tum.logr.segmented.vec[samp.seg.df[i,]$start_snp_idx:samp.seg.df[i,]$end_snp_idx] <- mean(tum.logr.vec[samp.seg.df[i,]$start_snp_idx:samp.seg.df[i,]$end_snp_idx])
		}#end for (i in 1:nrow(rle.baf.pcf.df)){
		#Let's do the BAF recalculation
		tum.baf.vec <- ascat.bc[["Tumor_BAF"]][,1]
		tum.baf.segmented.vec <- ascat.bc[["Tumor_BAF_segmented"]][[1]][,1]

		#So now this is a little more complicated as the BAF segmented is only in the SNPs
		#Would be nice to make the compaison with the germline here as well. We have the data so could...
		for (i in 1:nrow(samp.seg.df)){ 

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

		  	message("Number of phased SNPs is different to the actual number hz snps present in seg, perhaps to few to phase by MPhase?")
		  	message("nrow(seg.ph.hz.snp.df) = ", nrow(seg.ph.hz.snp.df), " nrow(gl.seg.hz.snp.df) = ", nrow(gl.seg.hz.snp.df))
		  	next
		  }#end if(nrow(seg.ph.hz.snp.df) != nrow(gl.seg.hz.snp.df)) {
		  #GL specific
		  gl.seg.hz.snp.df <- cbind(gl.seg.hz.snp.df, seg.ph.hz.snp.df[, c("maj_is_cpnA_logic_vec", "maj_is_cpnB_logic_vec", "snp_pos", "snp_idx")])
		  if(nrow(seg.hz.snp.df) == 0 | nrow(seg.ph.hz.snp.df) == 0){

		  	message("No heterozygous SNPs present within segment.")
		  	#print(samp.seg.df[i,])

		  }
		  if (nrow(seg.hz.snp.df) != nrow(seg.ph.hz.snp.df)) {

		  	message("Differing numbers of hz snps in the original hz snp df and the phased. Not sure why?")
		  	#print(samp.seg.df[i,])

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

		  #So we actually
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


		if(any(is.na(tum.logr.segmented.vec))) {

			message("There are NA values in the tum.logr.segmented.vec after applying phasing which need to be removed. Revisit generation of this input. For changing the NAs to 0.")
			na.segmented.idx.vec <- which(is.na(tum.logr.segmented.vec))
			tum.logr.segmented.vec[na.segmented.idx.vec] <- 0

		}

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
		samp.ploidy <- ploidy.est.vec[samp]
		samp.purity <- purity.est.vec[samp]

		setwd(qc.choice.dir)
		ascat.plotSegmentedData(ascat.bc)
		#So now somehow there's a missing value resulting from running ascat on the ascat.bc object.
		# List of 14
		#  $ Tumor_LogR          :'data.frame':	675964 obs. of  1 variable:
		#  $ Tumor_BAF           :'data.frame':	675964 obs. of  1 variable:
		#  $ Tumor_LogR_segmented: num [1:675964, 1] -0.838 -0.838 -0.838 -0.838 -0.838 ...
		#   ..- attr(*, "dimnames")=List of 2
		#  $ Tumor_BAF_segmented :List of 1
		#  $ Germline_LogR       :'data.frame':	675964 obs. of  1 variable:
		#  $ Germline_BAF        :'data.frame':	675964 obs. of  1 variable:
		#  $ SNPpos              :'data.frame':	675964 obs. of  2 variables:
		#  $ ch                  :List of 22
		#  $ chr                 :List of 28
		#  $ chrs                : chr [1:22] "1" "2" "3" "4" ...
		#  $ samples             : chr "CB1001-11-D01-X1"
		#  $ gender              : chr "XX"
		#  $ sexchromosomes      : chr [1:2] "X" "Y"
		#  $ failedarrays        : NULL
		#>
		# ascat.bc$BAFsegmented
		# any(is.na(ascat.bc$Tumor_LogR[,1]))
		# which(is.na(ascat.bc$Tumor_LogR[,1]))

		# any(is.na(ascat.bc$Tumor_LogR_segmented[,1]))
		missing.logr.segmented.idx.vec <- which(is.na(ascat.bc$Tumor_LogR_segmented[,1]))
		checking.cont.range <- missing.logr.segmented.idx.vec - c(missing.logr.segmented.idx.vec[2:length(missing.logr.segmented.idx.vec)], missing.logr.segmented.idx.vec[length(missing.logr.segmented.idx.vec)] + 1)  

    #any(is.na(ascat.bc$Tumor_LogR[,1]))

		ascat.output = ascat.runAscat(ascat.bc,gamma=1, rho_manual = samp.purity, psi_manual = samp.ploidy)
		setwd(base.dir)
		#setwd("/fast/work/groups/ag_schwarz/Projects/project-refphase/tom/tools/multiregioncopynum/")
		ph.ascat.output.list[[samp]] <- ascat.output
		save(ascat.bc, file=file.path(qc.choice.dir, paste(samp, 'ascat.bc.RData',sep='.')))
		save(ascat.output, file=file.path(qc.choice.dir, paste(samp, 'ascat.output.RData',sep='.')))
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
		save(seg, file=file.path(qc.choice.dir, paste(samp,'ascat.seg.RData',sep='.')))
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
		if(grepl("RData", pt.snp.file.paths.df[samp,]$samp_ascat_bc_path)) {
		  print("ascat.bc object in RData file, loading...")
		  load(pt.snp.file.paths.df[samp,]$samp_ascat_bc_path)
		  #This will have loaded the ascat.bc object in

		} else {
		  ascat.bc <- readRDS(pt.snp.file.paths.df[samp,]$samp_ascat_bc_path)
		}

    # Validate the input immediately
    if (conf$validate_input)
      validate_ascat_bc_input(ascat.bc)

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

	#OJO.... FIX THIS
	#So we're going to drop the singe bp segments that we get here... with a bit of effort they are rescuable
	#https://support.bioconductor.org/p/119589/
	#Or could be merged into a larger segment as they produce the following error:
	# Error in validObject(.Object) : invalid class "GRanges" object:
	#     'x@strand' is not parallel to 'x'
	# > traceback()
	# 15: stop(msg, ": ", errors, domain = NA)
	# 14: validObject(.Object)
	# 13: initialize(value, ...)
	# 12: initialize(value, ...)
	# 11: new(Class, seqnames = seqnames, ranges = ranges, strand = strand,
	#         elementMetadata = mcols, seqinfo = seqinfo)
	# 10: new_GRanges("GRanges", seqnames = seqnames, ranges = ranges,
	#         strand = strand, mcols = mcols, seqlengths = seqlengths,
	#         seqinfo = seqinfo)
	# 9: GRanges(ans_seqnames, ans_ranges, strand = ans_strand, ans_mcols,
	#        seqinfo = ans_seqinfo)
	# 8: makeGRangesFromDataFrame(f.joint.segs.df[f.joint.segs.df$sample ==
	#        samp, ], keep.extra.columns = FALSE, ignore.strand = FALSE,
	#        seqinfo = NULL, seqnames.field = c("seqnames", "seqname",
	#            "chromosome", "chrom", "chr", "chromosome_name", "seqid"),
	#        start.field = c("Start_Position", "start", "startpos"), end.field = c("End_Position",
	#            "end", "stop", "endpos"), strand.field = "strand", starts.in.df.are.0based = FALSE) at data_input_and_QC_functions_v2.R#856
	# 7: getAllPatientSamplesAndIntersectionGRangesList(sample.vec, f.joint.segs.df) at data_input_and_QC_functions_v2.R#904
	# 6: compareOverlapsBetweenSamples(sample.vec, f.joint.segs.df, chrom.info.df,
	#        chrs.to.use.vec = CHRS_TO_USE_VEC) at minimum_consistent_segmentation_functions.R#306
	# 5: createFullyAnnotatedJointSegmentation(f.joint.segs.df = ph.f.joint.segs.df) at xphase.R#1436
	# 4: eval(ei, envir)
	# 3: eval(ei, envir)
	# 2: withVisible(eval(ei, envir))
	ph.f.joint.segs.df <- ph.f.joint.segs.df[ph.f.joint.segs.df$startpos != ph.f.joint.segs.df$endpos,]
	rownames(ph.f.joint.segs.df) <- 1:nrow(ph.f.joint.segs.df)





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
	#print(str(ph.asas.SNPs.w.both.lr.baf.snp.lists))
	if (length(unique(sapply(ph.asas.SNPs.w.both.lr.baf.snp.lists, length))) != 1) {

		stop("Not all segments are present in all samples' snp assignations, likely something has gone wrong, is the parallelisation working correctly?")

	}
	#print(ph.asas.SNPs.w.both.lr.baf.snp.lists.path)
	saveRDS(ph.asas.SNPs.w.both.lr.baf.snp.lists, file = ph.asas.SNPs.w.both.lr.baf.snp.lists.path)
	time.taken <- proc.time() - ptm
	print("time.taken")
	print(time.taken)
	##
	#SNPs that only have LogR are permitted in these lists, but it also contains all the SNPs with both BAF and LogR
	ptm <- proc.time()
	message("Phased cleanAndErrorCheckSNPsParallel second")

	ph.asas.SNPS.lr.and.or.baf.lists <- cleanAndErrorCheckSNPsParallel(all.samp.SNP.list,
																	mcs.seg.to.use.df =ph.min.cons.segs.df,
																	allow.logr.only.flag = TRUE,
																	data.type = data.format.to.use)

	if (length(unique(sapply(ph.asas.SNPS.lr.and.or.baf.lists, length))) != 1) {

		stop("Not all segments are present in all samples' snp assignations, likely something has gone wrong, is the parallelisation working correctly?")

	}

	#print(ph.asas.SNPS.lr.and.or.baf.lists.path)
	saveRDS(ph.asas.SNPS.lr.and.or.baf.lists, file = ph.asas.SNPS.lr.and.or.baf.lists.path)
	time.taken <- proc.time() - ptm
	print("time.taken")
	print(time.taken)

	if ( length(ph.asas.SNPS.lr.and.or.baf.lists) <= 1 ) {

		#print(str(ph.asas.SNPS.lr.and.or.baf.lists.path, max = 1))
		stop("Only one or fewer samples present after reading in data.")

	}



	message("Getting phasing to apply to all samples using getCombinedBAFDipTestPhasingParallel...")
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
	message("Generating phased BAF dataframes...")
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
	message("Performing LTNPD MSAI detection on the BAF...")
	all.seg.MSAI.by.BAF.ltnpd.present.list <- getAllSegsMSAIDetectionFromBAFListLTNPD(seg.name.vec = seg.name.vec,
																				  sample.vec  = sample.vec,
																				  asas.SNPS.lr.and.or.baf.lists = asas.SNPS.lr.and.or.baf.lists,
																				  all.seg.phasing.info.list = all.seg.phasing.info.list,
																				  MSAI.BAF.LTPND.pval.thresh = MSAI_BAF_LTNPD_PVAL_THRESHOLD)
	#saveRDS(ph.all.seg.MSAI.by.BAF.ltnpd.present.list, file = ph.all.seg.MSAI.by.BAF.ltnpd.present.list.path)


	if (PLOTTING_FLAG == TRUE) {

		png(file = paste(PLOT_OUTPUT_DIR, patient, "_MSAI_detection_BAF_phased", ".png", sep=""), width = 1200, height = 1400)

		message("Plotting phased BAF across the genome...")
		plotAllSampleBAFsAndMSAIAcrossGenome(sample.vec = sample.vec,
											all.seg.MSAI.by.BAF.ltnpd.present.list = all.seg.MSAI.by.BAF.ltnpd.present.list,
											patient.phased.BAF.df.list = patient.phased.BAF.df.list,
											cpnA.allele.color = '#de2d2699',
											cpnB.allele.color= '#3182bd99')
		dev.off()

	} else if (PLOTTING_FLAG == FALSE) {

		message("Not plotting phased BAF as flag is FALSE")

	}

	message("Applying LTNPD MSAI detection to the original segments...")
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
	#print(orig.intersection.mcs.w.MSAI.df[orig.intersection.mcs.w.MSAI.df$nMinor > orig.intersection.mcs.w.MSAI.df$nMajor,])

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

	message("Generating datapoints with LogR total cpn dataframes...")
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
	message("Applying phasing to all samples to get phased raw cpn...")
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
	message("Checking all minimum consistent segments are present in each sample...")
	seg.name.vec <- names(asas.ph.raw.cpn.lists[[1]])
	#Check to make sure that all segments are present in all samples... might not be necessary todo, check.
	for ( samp in sample.vec ) {

	      seg.name.vec <- intersect(seg.name.vec, names(asas.ph.raw.cpn.lists[[samp]]))

	}#end for ( samp in sample.vec ) {

	if (length(seg.name.vec) == 0) {

		stop("No shared segments for some reason, go back and check")

	}

	message("Calculating CCFs using phased raw cpn...")
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
	  #print(cpnX)
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
	message("Initialising minimum consistent segmentation (MCS) output dataframe...")
	empty.full.joint.segs.df <- initialiseFullOutputDataframe(seg.name.vec = names(asas.ph.raw.cpn.lists[[1]]),
															ploidy.est.vec = ploidy.est.vec,
															purity.est.vec = purity.est.vec,
															make.row.names.seg.name.alone = FALSE)


	#This adds the raw cn for both segments that could be phased and those that couldn't
	message("Adding new copy number calls to MCS output based on the phased raw copy number...")
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





	message("Adding CCF information and subclonal cpn calls to MCS output...")
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

	message("Adding total copy number estimates to MCS output for segments with only logr information...")
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

	#print(head(full.joint.segs.df.w.cpn.vs.ploidy.df[,"cpn_event_vs_ploidy"]))
	message("Adding cpn event vs ploidy by logr measures to MCS output...")
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


	message("Testing whether any gains should be reclassified as amplifications...")
	# full.joint.check.amps.vs.ploidy.clonality.df <- retestSubclonalAmps(all.seg.df = full.joint.segs.logr.vs.ploidy.clonality.df,
	# 																	all.samp.logr.df.list = all.samp.non.hz.logr.and.tot.cpn.df.list,
	# 																	all.samp.ph.BAF.df.list = patient.phased.BAF.df.list,
	# 																	gain.threshold = GAIN_LOGR_THRESH,
	# 																	amp.threshold = AMP_LOGR_THRESH,
	# 																	p.val.for.converting.gain.to.amp = PVAL_THRESH_CONVERT_GAIN_TO_AMP)

	#print(head(full.joint.segs.logr.vs.ploidy.clonality.df[,c("across_samp_cl_cpn_vs_ploidy_present","across_samp_sc_cpn_vs_ploidy_present")]))


	#These are kept in separate functions as the original segmentaiton was refactored so it could run when segments are rejoined at the end of the script.
	message("Adding phasing information and MSAI detection by raw cpn calls to MCS output...")
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

	message("Adding MSAI detected by LTNPD method to the MCS output...")
	full.joint.segs.w.ltnpd.MSAI.df <- addMSAIDetectionByLTNPDToOutputSegmentation(sample.vec = sample.vec,
												seg.name.vec = names(asas.ph.raw.cpn.lists[[1]]),
												all.seg.MSAI.by.BAF.ltnpd.present.list = all.seg.MSAI.by.BAF.ltnpd.present.list,
												all.segs.df = full.joint.segs.w.phasing.df)#full.joint.segs.w.phasing.and.raw.df)


	#print(head(full.joint.segs.w.ltnpd.MSAI.df[( ! is.na(full.joint.segs.w.ltnpd.MSAI.df$any_across_samp_MSAI)) &
	#											full.joint.segs.w.ltnpd.MSAI.df$any_across_samp_MSAI == TRUE,
	#											c("MSAI", "any_across_samp_MSAI")]))


	message("Defining any parallel gain or loss events detected via MSAI and addingto MCS output...")
	full.joint.segs.w.parallel.cpn.events.df <- addParallelEventsByMSAIToOutputSegmentation(all.segs.df = full.joint.segs.w.ltnpd.MSAI.df)


	#print(head(full.joint.segs.w.parallel.cpn.events.df[(! is.na(full.joint.segs.w.parallel.cpn.events.df$MSAI)),
	#													c("any_across_samp_MSAI", "MSAI",
	#														"cpn_event_vs_ploidy",
	#														"is_parallel_gain_and_or_amp",
	#														"is_parallel_loss",
	#														"is_parallel_LOH",
	#														"is_parallel_loss_and_or_LOH")]))

	message("Adding comparison between original cpn calls and MPhase copy number calling to MCS output...")
	full.joint.segs.within.and.changed.cpn.info.df <- addWithinSampleChangeToCpnAndAllelicImbalanceCalls(all.seg.df = full.joint.segs.w.parallel.cpn.events.df,
																										sample.vec = sample.vec,
																										seg.name.vec = unique(seg.name.vec))



	message("Adding across sample estimates of allelic imbalance and LOH clonality to MCS output...")
	full.joint.segs.w.across.samp.AI.stats.df <- addAcrossSampleAlleleSpecificClonalityCalls(all.seg.df = full.joint.segs.within.and.changed.cpn.info.df,
																							sample.vec = sample.vec,
																							seg.name.vec = seg.name.vec)

	#print(head(full.joint.segs.w.across.samp.AI.stats.df[!is.na(full.joint.segs.w.across.samp.AI.stats.df$within_samp_sc_cpnA_present), c("any_within_samp_sc_cpnA",
	#"any_within_samp_sc_cpnB",
	#"across_samp_cl_AI",
	#"across_samp_sc_AI",
	#"across_samp_cl_LOH",
	#"across_samp_sc_LOH")]))

	#Legacy mostly, trying to make definitions that Emilia was doing and put them ion here to ensure they were correct
	full.joint.segs.w.across.plotting.cols.df <- addPanCanColumnDefinitionToSegDF(full.joint.segs.w.across.samp.AI.stats.df)

	#print(full.joint.segs.w.across.plotting.cols.df[, c("pt_seg_cl_loss_LOH_no_MSAI", "pt_seg_sc_loss_LOH_incl_MSAI", "pt_seg_cl_gain_amp_no_MSAI", "pt_seg_sc_gain_amp_incl_MSAI")])

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
	message("Plotting CCFs for cpnA and cpnB across the genome using cpn definition of subclonality...")
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
	message("Not plotting CCFs for cpnA and cpnB across the genome using cpn definition of subclonality as plotting flag is FALSE.")
}

# We need a .tsv version of this data for interoperability with python
write.table(full.joint.segs.df, file=file.path(PLOT_OUTPUT_DIR, "full_joint_segs.tsv"), row.names=FALSE, sep="\t")
message("MPhase comparison with snpphase: Finished for now. Got the <patient>_full_joint_segs_df.RDS object.")

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
