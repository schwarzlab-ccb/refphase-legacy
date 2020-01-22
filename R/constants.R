# Define some global constants

MIN_NUM_HZ_SNP_PHASING <- 5
OVERWRITE_FLAG <- conf$overwrite
PLOTTING_FLAG <- FALSE
# POST JOINING FLAGS to run sections of code
JOINED_SEGS_ERROR_CHECKING_FLAG <- FALSE
# Using joined segs on second run throughof MPhase
OVERWRITE_INIT_MCS_ADD_INFO_FLAG <- TRUE
OVERWRITE_ARM_INFO_FOR_JOIN_INIT_MCS_FLAG <- TRUE
OVERWRITE_JOIN_INIT_MCS_FLAG <- TRUE
OVERWRITE_JOINED_CCF_CALC_FLAG <- TRUE
OVERWRITE_JOINED_MCS_ADD_INFO_FLAG <- TRUE
CALCULATE_JOINED_ARM_FOCAL_FLAG <- TRUE
CALCULATE_FOCAL_EVENT_CLONALITY_FLAG <- TRUE
GENERATE_GISTIC_INPUT_FLAG <- FALSE
PLOT_OUTPUT_DIR <- paste(output.directory.base.path, "/", sep = "")
INTERSECTION_OR_ARM <- intersection.or.whole.arm
GAMMA_VAL <- as.numeric(gamma.val)
MAX_CPN_FOR_SUBCLONALITY <- 10
# This the median raw CCF different to the integer +/- this value.
MIN_CCF_TO_CONSIDER <- 0.2
# This the median raw cpn different to the integer +/- this value.
MIN_CPN_TO_CONSIDER <- 0.2
# This is the minimum no of chromosomes used for comparisons to other other
# chromosomes when calling subclonal events
MIN_NUM_CHRS <- 2
GAIN_LOGR_THRESH <- log2(2.5 / 2)
LOSS_LOGR_THRESH <- log2(1.5 / 2)
AMP_LOGR_THRESH <- log2(4 / 2)
DEEP_LOSS_LOGR_THRESH <- log2(0.5 / 2)
PVAL_THRESH_VS_LOGR_FOR_CPN_EVENT <- 0.01
CHRS_TO_USE_VEC <- 1:22
# used for CCF ca;culations
CONF_INT_FOR_CPNX_BAF_VS_EXP <- 0.9999
# Legacy, used for pairwise comparisons between samples
PW_CPNA_CPNB_DIFF_AMPLITUDE_THRESH <- 0.001
# P value thresholds for  joining segements on phased raw cpn values
PVAL_THRESH_TO_JOIN_PH_SEGS <- 0.001
# Some segments don't have phased copy number but we can try and join based on
# the total cpn
PVAL_THRESH_TO_JOIN_TOT_SEGS <- 0.001
# legacy minimum  proportion of snps to be separated by phasing, if fewer, try
# a more complicated phasing currently deactivated 11/09/2018
MIN_PROP_SNP_SEP_THRESH <- 0.75
# threshold for testing whether the BAFs are separated across all the samples
# when comparing phasings
PVAL_COMB_WILCOX_PH_THRESH <- 0.01
# Test to decide whether there is allelic imbalance, not being used as of
# 11/09/2018
PVAL_DIPTEST_AI_THRESH <- 0.01
ALLOW_DIPTEST_ALONE_PH_FLAG <- allow.dip.test.phasing.flag
# P value threshold to call MSAI
MSAI_BAF_LTNPD_PVAL_THRESHOLD <- 0.001
# Legacy, can also detect MSAI from the phased copy number not ust the BAF
MSAI_CPNX_PVAL_THRESHOLD <- 0.001
# Legacy, but might need to be re-implemented, is an amplification
# significantly different to a gain in a different region?
PVAL_THRESH_CONVERT_GAIN_TO_AMP <- 0.01
# GLOBAL FLAG that alters the behaviour of some hacked functions so that they
# work with MPhase and XPhase
USING_XPHASE <- FALSE

PVAL_THRESH_VS_MIN_CCF <- 0.001
# T test p value threshold for trying to call subclonal copy numnber...
PVAL_THRESH_VS_OTH_CHRS_SAME_CPN <- 0.01
