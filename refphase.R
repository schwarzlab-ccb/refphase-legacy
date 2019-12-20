#!/usr/bin/env Rscript

library(optparse)
library(digest)

options(stringsAsFactors = FALSE)

# Returns the location of this .R script (may be needed to source other files in same dir)
# function taken from https://github.com/molgenis/molgenis-pipelines/wiki/How-to-source-another_file.R-from-within-your-R-script
LocationOfThisScript <- function() {
  this.file <- NULL
  # This file may be 'sourced'
  for (i in -(1:sys.nframe())) {
    if (identical(sys.function(i), base::source)) {
      this.file <- (normalizePath(sys.frame(i)$ofile))
    }
  }

  if (!is.null(this.file)) {
    return(dirname(this.file))
  }

  # But it may also be called from the command line
  cmd.args <- commandArgs(trailingOnly = FALSE)
  cmd.args.trailing <- commandArgs(trailingOnly = TRUE)
  cmd.args <- cmd.args[seq.int(from = 1, length.out = length(cmd.args) - length(cmd.args.trailing))]
  res <- gsub("^(?:--file=(.*)|.*)$", "\\1", cmd.args)

  # If multiple --file arguments are given, R uses the last one
  res <- tail(res[res != ""], 1)
  if (0 < length(res)) {
    return(dirname(res))
  }

  # Both are not the case. Maybe we are in an R GUI?
  return(NULL)
}


if (!interactive()) {
  option_list <- list(
    make_option(c("--patient"),
      default = "CB1002",
      help = "Patient id"
    ),
    make_option(c("--dip-test"),
      default = FALSE, action = "store_true",
      help = "Use dip test [default = %default]"
    ),
    make_option(c("--output-dir"),
      default = "../output/CB1002/mphase",
      help = "Output directory"
    ),
    make_option(c("--medicc-output-dir"),
      default = NA, type = "character",
      help = "Output directory for MEDICC input files"
    ),
    make_option(c("--input-format"),
      default = "ASCAT",
      help = "Input data format [default = %default]"
    ),
    make_option(c("--snp-data"),
      default = "../input/nb-schulte//CB1002_input_paths_df.txt",
      help = "SNP data"
    ),
    make_option(c("--snp-data-prefix"),
      default = "",
      help = "Directory prefix for the paths to the SNP data in the file above"
    ),
    make_option(c("--seg-data"),
      default = "../input/nb-schulte//CB1002.segmentation",
      help = "Segmentation data output from ASCAT"
    ),
    make_option(c("--gamma"),
      default = 1.0, type = "numeric",
      help = "Gamma (for seq=1, for arrays~=0.5 [default = %default]"
    ),
    make_option(c("-j", "--cores"),
      default = 1, type = "integer",
      help = "CPU cores [default = %default]"
    ),
    make_option(c("--segment-type"),
      default = "intersection",
      help = "Segment type, intersection or whole_arm [default = %default]"
    ),
    make_option(c("--genome"),
      default = "hg38",
      help = "Human genome version [default = %default] (supports hg19 and hg38)"
    ),
    make_option(c("--overwrite"),
      default = TRUE, action = "store_true",
      help = "Overwrite output [default = %default]"
    ),
    make_option(c("--no-overwrite"),
      action = "store_false", dest = "overwrite",
      help = "Do not overwrite previous output"
    ),
    make_option("--cached-dir",
      default = NA,
      help = "Set this to load intermediate data from the given directory"
    ),
    make_option("--debug-dir",
      default = NA, type="character",
      help = "Set this to save intermediate data to the given directory for later debugging"
    ),
    make_option(c("--validate-input"),
      default = TRUE, action = "store_true", dest = "validate_input",
      help = "Validate input files and abort if input is not valid"
    ),
    make_option(c("--no-validate-input"),
      action = "store_false", dest = "validate_input",
      help = "Do not validate input"
    ),
    make_option(c("--seed"),
      default = 2134234,
      help = "Random seed"
    )
  )

  # Hack to get it working on TRACERx
  if (packageVersion("optparse") != "1.3.2") {
    conf <- parse_args(OptionParser(option_list = option_list), convert_hyphens_to_underscores = TRUE)
  } else {
    conf <- parse_args(OptionParser(option_list = option_list))
    names(conf) <- gsub("\\-", "_", names(conf))
  }
  conf$id <- substr(digest(conf), 1, 8)
  args <- commandArgs(trailingOnly = TRUE)

  dput(conf)
} else {
  current.dir <- LocationOfThisScript()
  source(file.path(current.dir, "debugging_conf.R"))
}

message("Arguments:")
for (key in names(conf))
  message(key, " = ", conf[[key]], " (class = ", class(conf[[key]]), ")")

set.seed(conf$seed) # Why is a seed required? I'm not sure anything non-deterministic should be occurring


patient <- conf$patient
# At the moment MPhase will only try to call MSAI and new copy number in areas of the genome identified as allelically imbalanced
# Enable this flag to try and find more allelic imbalance but this will probably break.
allow.dip.test.phasing.flag <- conf$dip_test
output.directory.base.path <- conf$output_dir
data.format.to.use <- conf$input_format

# all.samp.snp.list.path <- conf$cached_data #Points towards the snp files, of whatever format they are given that should be used based on data.format to use.
snp.file.path <- conf$snp_data # Points towards the snp files, of whatever format they are given that should be used based on data.format to use.
seg.file.path <- conf$seg_data # Pintts towards the seg files, of whatever format they should be determined byh the data.format.to.use flag.
gamma.val <- conf$gamma # for sequencing it's 1 and SNP array it's 0.4/0.5 but this is platform specific see ASCAT help pages.
intersection.or.whole.arm <- conf$segment_type # Can either be intersection or whole_arm, do we take the only intersection of segments from each sample
genome <- conf$genome # hg19 or hg38
# or do attempt to fill in the gaps between every segment in each sample so that the entire genome is covered (so output from any sample can run in GISTIC
# as all will cover the same genomic area)

# Right we need a free run folder:
qc.choice.dir <- file.path(output.directory.base.path, "run_qc_choice")
dir.create(qc.choice.dir, recursive = TRUE, showWarnings = FALSE)


# We need to be able to select the intersection of all regions or a full GISTIC compatible union with defined regions that will always be covered
# What proportion of the chosen genome, in this case hg19, arms are covered by the segments that have been inputted.
# We also need to add in the wgii calculations

### Source relevant scripts, best to do this at the top to check for errors...
current.dir <- LocationOfThisScript()
if (is.null(current.dir)) {
  current.dir <- ""
} else {
  current.dir <- paste0(current.dir, "/")
}
source(file.path(current.dir, "utils.R"))
source(file.path(current.dir, "minimum_consistent_segmentation.R"))
source(file.path(current.dir, "input_and_qc.R"))
source(file.path(current.dir, "consensus_phasing.R"))
source(file.path(current.dir, "copy_number.R"))
source(file.path(current.dir, "msai.R"))
source(file.path(current.dir, "cancer_cell_fraction.R"))
source(file.path(current.dir, "prepare_output_segmentation.R"))
source(file.path(current.dir, "new_seg_joining_functions.R"))
source(file.path(current.dir, "plot_cpn_heterogeneity.R"))
source(file.path(current.dir, "output_medicc.R"))
source(file.path(current.dir, "ascat.R"))

# New files to reduce the length of this script
source(file.path(current.dir, "annotate_joint_segs.R"))
source(file.path(current.dir, "input.R"))
source(file.path(current.dir, "validate_input.R"))
source(file.path(current.dir, "phase_ascat.R"))
source(file.path(current.dir, "input_ascat.R"))
source(file.path(current.dir, "input_sequenza.R"))
source(file.path(current.dir, "input_medicc.R")) # Not used currently
source(file.path(current.dir, "ascat_utils.R"))
source(file.path(current.dir, "ccf_from_raw_cpn.R"))
source(file.path(current.dir, "add_info_to_output_table.R"))

suppressWarnings(suppressMessages(library("diptest")))
suppressWarnings(suppressMessages(library("doParallel")))
suppressWarnings(suppressMessages(library("gdata")))
suppressWarnings(suppressMessages(library("GenomicRanges")))
suppressWarnings(suppressMessages(library("gplots")))

# When debugging interactively, we need to start the debugger after sourcing
# the code above, so that we can set breakpoints on functions of interest
if (interactive()) {
  browser()
}

# -------------------------------------------------------------------------
# Set constants for the analysis
# -------------------------------------------------------------------------
source(file.path(current.dir, "constants.R"))

# -------------------------------------------------------------------------
# Create output dirs and check that input files exist
# -------------------------------------------------------------------------
if (!file.exists(PLOT_OUTPUT_DIR)) {
  message("Output directory doesn't exist:", PLOT_OUTPUT_DIR)
}

output_dir <- PLOT_OUTPUT_DIR
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

if (is.null(output.directory.base.path)) {
  stop("No output path provided, aborting...")
}

# Make sure the input files exist before going too far
if (!file.exists(snp.file.path)) {
  stop("No file at the snp.file.path provided")
}

if (!file.exists(seg.file.path)) {
  stop("No file at the seg.file.path provided")
}


registerDoParallel(cores = conf$cores)

# Load in the relevant genome file
# We don't need this until much later, don't try to load it now
if (genome == "hg19") {
  info.path <- file.path(current.dir, "annot", "ucsc_hg19_chrominfo.txt")
} else if (genome == "hg38") {
  info.path <- file.path(current.dir, "annot", "ucsc_hg38_chrominfo.txt")
} else {
  stop("Genome version not recognized: ", genome)
}
chrom.info.df <- read.table(info.path, stringsAsFactors = FALSE, sep = "\t", header = TRUE)

# Paths for saving/loading objects
if (!is.na(conf$debug_dir)) {
  dir.create(conf$debug_dir, recursive = TRUE, showWarnings = FALSE)

  input.snp.file.rds.path <- file.path(conf$debug_dir, "input_snp_file.rds")
  all.samp.SNP.list.path <- file.path(conf$debug_dir, "all_samp_SNP_list.rds")
  asas.SNPS.lr.and.or.baf.lists.path <- file.path(conf$debug_dir, "asas_SNPS_lr_and_or_baf_lists_lists.rds")
  asas.SNPs.w.both.lr.baf.snp.lists.path <- file.path(conf$debug_dir, "asas_SNPs_w_both_lr_baf_snp_lists.rds")
  min.cons.segs.df.write.path <- file.path(conf$debug_dir, "min_cons_segs_df.txt")
  sorted.all.bp.all.samps.segs.anno.df.path <- file.path(conf$debug_dir, "sorted_all_bp_all_samps_segs_anno.df.rds")
  min.cons.segs.df.path <- file.path(conf$debug_dir, "min_cons_segs_df.rds")
  processed.orig.cpn.segs.df.path <- file.path(conf$debug_dir, "processed_orig_cpn_segs_df.rds")

  ploidy.est.vec.path <- file.path(conf$debug_dir, "ploidy_est_vec.rds")
  purity.est.vec.path <- file.path(conf$debug_dir, "purity_est_vec.rds")
  # Analysis object paths
  all.seg.phasing.info.list.path <- file.path(conf$debug_dir, "all_seg_phasing_info_list.rds")
  anno.seg.df.path <- file.path(conf$debug_dir, "anno_seg_df_path.rds")
  patient.phased.BAF.df.list.path <- file.path(conf$debug_dir, "phased_BAF_df_list.rds")
  all.seg.MSAI.by.BAF.ltnpd.present.list.path <- file.path(conf$debug_dir, "all_seg_MSAI_by_BAF_ltnpd_present_list.rds")
  all.samp.non.hz.logr.and.tot.cpn.df.list.path <- file.path(conf$debug_dir, "all_samp_non_hz_logr_and_tot_cpn_df_list.rds")
  asas.ph.raw.cpn.lists.path <- file.path(conf$debug_dir, "asas_ph_raw_cpn_lists.rds")
  samps.int.cpn.raw.cpn.by.chr.list.path <- file.path(conf$debug_dir, "samps_int_cpn_raw_cpn_by_chr_list.rds")
  all.seg.CCF.list.path <- file.path(conf$debug_dir, "all_seg_CCF_list.rds")
  all.seg.int.ob.xp.cpn.df.list.path <- file.path(conf$debug_dir, "all_seg_int_ob_xp_cpn_df_list.rds")
  all.seg.CCF.diff.from.exp.list.path <- file.path(conf$debug_dir, "all_seg_CCF_diff_from_exp_list.rds")
  all.seg.clonal.in.all.samps.list.path <- file.path(conf$debug_dir, "all_seg_clonal_in_all_samps_list.rds")
  all.seg.greater.than.two.diff.list.path <- file.path(conf$debug_dir, "all_seg_greater_than_two_diff_list.rds")
  pt.snp.file.paths.path <- file.path(conf$debug_dir, "pt_snp_file_paths_df.rds")
}

# -------------------------------------------------------------------------
# Main code starts here
# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
# Load in the data:
# -------------------------------------------------------------------------

use_cached_input <- FALSE
if (!use_cached_input) {
  input_data <- read_input(seg.file.path, snp.file.path, conf$input_format)
  saveRDS(input_data, file = file.path(PLOT_OUTPUT_DIR, "input_data.rds"))
} else {
  message("Loading cached input (may take a minute)")
  input_data <- readRDS(file.path(PLOT_OUTPUT_DIR, "input_data.rds"))
  message("done!")
}

# Unwrap the list returned from the read_input() function
all.samp.SNP.list <- input_data$all.samp.SNP.list
anno.seg.df <- input_data$anno.seg.df
asas.SNPS.lr.and.or.baf.lists <- input_data$asas.SNPS.lr.and.or.baf.lists
asas.SNPs.w.both.lr.baf.snp.lists <- input_data$asas.SNPs.w.both.lr.baf.snp.lists
f.joint.segs.df <- input_data$f.joint.segs.df
min.cons.segs.df <- input_data$min.cons.segs.df
ploidy.est.vec <- input_data$ploidy.est.vec
pt.snp.file.paths.df <- input_data$pt.snp.file.paths.df
purity.est.vec <- input_data$purity.est.vec
sorted.all.bp.all.samps.segs.anno.df <- input_data$sorted.all.bp.all.samps.segs.anno.df

if (!is.na(conf$debug_dir)) {
  # Save input data for later debugging
  message("Save output for later")
  saveRDS(all.samp.SNP.list, file = all.samp.SNP.list.path)
  saveRDS(anno.seg.df, file = anno.seg.df.path)
  saveRDS(asas.SNPS.lr.and.or.baf.lists, file = asas.SNPS.lr.and.or.baf.lists.path)
  saveRDS(asas.SNPs.w.both.lr.baf.snp.lists, file = asas.SNPs.w.both.lr.baf.snp.lists.path)
  saveRDS(f.joint.segs.df, file = processed.orig.cpn.segs.df.path)
  saveRDS(min.cons.segs.df, file = min.cons.segs.df.path)
  saveRDS(ploidy.est.vec, file = ploidy.est.vec.path)
  saveRDS(pt.snp.file.paths.df, file = pt.snp.file.paths.path)
  saveRDS(purity.est.vec, file = purity.est.vec.path)
  saveRDS(sorted.all.bp.all.samps.segs.anno.df, file = sorted.all.bp.all.samps.segs.anno.df.path)
}

#message("Number of rows in each sample (asas.SNPs.w.both.lr.baf.snp.lists)")
snps_per_segment <- sapply(asas.SNPs.w.both.lr.baf.snp.lists, function(x) {
  sapply(x, nrow)
})
message("Segments with enough hz positions: >= ", MIN_NUM_HZ_SNP_PHASING)
print(snps_per_segment[apply(snps_per_segment >= MIN_NUM_HZ_SNP_PHASING, 1, all), ])
message(sum(apply(snps_per_segment >= MIN_NUM_HZ_SNP_PHASING, 1, all)), "/", nrow(snps_per_segment))

message("Setting the pval thresholds for calling subclonality depending on the number of heterozygous sites assayed.")
total.num.of.hz.sites.per.sample <- sum(sapply(asas.SNPs.w.both.lr.baf.snp.lists[[1]], nrow))
if (total.num.of.hz.sites.per.sample > 300000) {
  PVAL_THRESH_VS_MIN_CCF <- 0.00000001
  PVAL_THRESH_VS_OTH_CHRS_SAME_CPN <- 0.0000000001
}

seg.name.vec <- paste(min.cons.segs.df[, "Chromosome"], min.cons.segs.df[, "Start_Position"], min.cons.segs.df[, "End_Position"], sep = ":")
sample.vec <- names(asas.SNPs.w.both.lr.baf.snp.lists)

message("Getting phasing to apply to all samples using getCombinedBAFDipTestPhasingParallel")
all.seg.phasing.info.list <- getCombinedBAFDipTestPhasingParallel(
  asas.lr.baf.lists = asas.SNPs.w.both.lr.baf.snp.lists,
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

# Also both seem to be broken as there is never anything returned by samp_most_significant_by_diptest in either.

if (!is.na(conf$debug_dir)) {
  saveRDS(all.seg.phasing.info.list, file = all.seg.phasing.info.list.path)
}

# Todo unsure of what to do with the baf logr lists which have positions that only have logr here, at the moment exclude them but how
# does that affect later run
message("Generating phased BAF dataframes...")
patient.phased.BAF.df.list <- createAllSegsPhasedBAFDataFrame(
  asas.lr.baf.lists = asas.SNPs.w.both.lr.baf.snp.lists,
  sample.to.phase.seg.list = all.seg.phasing.info.list,
  min.num.hz.snp.phasing = MIN_NUM_HZ_SNP_PHASING
)

if (length(patient.phased.BAF.df.list) <= 1) {
  stop("Error in generating phased BAF dataframes")
}

if (!is.na(conf$debug_dir)) {
  saveRDS(patient.phased.BAF.df.list, file = patient.phased.BAF.df.list.path)
}

# Now we have the phasing, we want to check whether, using it, we think that
# there is sufficient separation of the bands of allelic balance to consider a
# segment "allelically imbalanced". If not set the bafsegmented to 0.5
# So we may want to reintroduce the log truncated normal stuff here...

# So this is where want to create a phased ascat.bc object:
# TBKW NEW WORK:
#message("Starting new TBKW edit... begun on 16/08/2019 to recall copy number using phasing.")

ret <- annotate_joint_segs(anno.seg.df, all.samp.SNP.list[[1]])
all.samp.j.seg.df <- ret$all.samp.j.seg.df
phased.baf.info.df <- ret$phased.baf.info.df

# Note there could be an issue here if not every segment contains a BAF point...
# But, if we're using hte intersection this should be fine? Not sure on this point.
# However it doesn't really matter as ASCAT makes sure every segment whether or
# not it had any baf point has an allele-specific value.
message("************************************************************************")
message("Running ASCAT on all samples")
message("************************************************************************")

# So now there are several things we want to do here depending on the options
# passed to the tool:
# 1) Using the phasing information obtained from all samples, allow a free run
#    allowing ASCAT to call whatever purity and ploidy of all samples included
#    in the run it feels is most appropriate. This will potentially provide new
#    and improved solutions that might be better than what was fed into mphase
#    particularly if the segmentation of the input did not make use of
#    multi-region phasing data.
# 2) Initiate list for segmentation
ph.ascat.input.list <- list()
ph.ascat.output.list <- list()
ph.ascat.segs.df.list <- list()
base.dir <- getwd()

USING_XPHASE <- TRUE
PH_ALL_SAMP_ASCAT_RAW_CPN_DF_LIST <- list()

for (samp in sample.vec) {
  message("Preparing sample for ASCAT: ", samp)

  if (conf$input_format == "ASCAT") {
    ascat.bc <- read_ascat_bc(pt.snp.file.paths.df[samp, ]$samp_ascat_bc_path)
  } else {
    # Create an ascat.bc object from the existing snp data, rather than
    # reading an existing ascat.bc object. This also removes NAs.
    ascat.bc <- snps_to_ascat_bc(all.samp.SNP.list[[samp]], samp)
  }

  ascat.bc.phased <- phase_ascat(
    ascat_bc = ascat.bc,
    phased_bafs = patient.phased.BAF.df.list[[samp]],
    segs = all.samp.j.seg.df,
    phased_baf_info = phased.baf.info.df
  )

  # OK so now let's do the fixed solution run...
  samp.ploidy <- ploidy.est.vec[samp]
  samp.purity <- purity.est.vec[samp]

  if (conf$validate_input)
    validate_ascat_bc_input(ascat.bc.phased)

  ph.ascat.input.list[[samp]] <- ascat.bc.phased

  setwd(qc.choice.dir)
  ascat.plotSegmentedData(ascat.bc.phased)
  ascat.output <- ascat.runAscat(ascat.bc.phased, gamma = 1, rho_manual = samp.purity, psi_manual = samp.ploidy)
  setwd(base.dir)

  ph.ascat.output.list[[samp]] <- ascat.output

  ascat.bc <- ascat.bc.phased
  saveRDS(ascat.bc, file = file.path(qc.choice.dir, paste(samp, "ascat.bc.rds", sep = ".")))
  saveRDS(ascat.output, file = file.path(qc.choice.dir, paste(samp, "ascat.output.rds", sep = ".")))

  # We need to work out where the MSAI is and flip it as well...
  # Calculate number of heterozygous positions per segment
  segs_raw <- ascat.output$segments_raw
  message("Number of segments in ASCAT rerun: ", nrow(segs_raw))
  n.het <- 0 # FIXME probably want to add these in.
  acf <- ascat.output$aberrantcellfraction
  ploidy <- ascat.output$ploidy

  segs_raw <- cbind(segs_raw[, c(1:4)], n.het, segs_raw[, 5] + segs_raw[, 6], segs_raw[, 5:6], rep(ploidy, nrow(segs_raw)), rep(acf, nrow(segs_raw)), segs_raw[, 7], segs_raw[, 8])
  colnames(segs_raw)[c(6, 9:12)] <- c("cnTotal", "Ploidy", "ACF", "nAraw", "nBraw")

  segs_raw[, "chr"] <- chrom2integer(segs_raw[, "chr"])

  # add a column defining gain and losslength
  callState <- rep("neutral", nrow(segs_raw))
  callState[log2((segs_raw[, "nAraw"] + segs_raw[, "nBraw"]) / segs_raw[, "Ploidy"]) >= log2(2.5 / 2)] <- "Gain"
  callState[log2((segs_raw[, "nAraw"] + segs_raw[, "nBraw"]) / segs_raw[, "Ploidy"]) <= log2(1.5 / 2)] <- "Loss"
  callState[(segs_raw[, "nAraw"] + segs_raw[, "nBraw"]) >= ((2 * segs_raw[, "Ploidy"]) + 1)] <- "Amp"
  callState[(segs_raw[, "nAraw"] + segs_raw[, "nBraw"]) <= 0.5] <- "Del"

  samp.ph.ascat.segs.df <- cbind(segs_raw, callState)
  ph.ascat.segs.df.list[[samp]] <- samp.ph.ascat.segs.df

  seg <- samp.ph.ascat.segs.df
  saveRDS(seg, file = file.path(qc.choice.dir, paste(samp, "ascat.seg.rds", sep = ".")))

  ph.ascat.nAraw.w.minus.vec <- ph.ascat.output.list[[samp]]$nArawfull$nArawfull
  ph.ascat.nBraw.w.minus.vec <- ph.ascat.output.list[[samp]]$nBrawfull$nBrawfull
  ph.ascat.nAraw <- ph.ascat.nAraw.w.minus.vec
  ph.ascat.nAraw[ph.ascat.nAraw < 0] <- 0
  ph.ascat.nBraw <- ph.ascat.nBraw.w.minus.vec
  ph.ascat.nBraw[ph.ascat.nBraw < 0] <- 0

  # Because the plotting requires flipped values make the phasing
  # and MSAI visible...
  new.ph.ascat.nAraw <- ph.ascat.nAraw
  new.ph.ascat.nBraw <- ph.ascat.nBraw
  new.ph.ascat.nAraw.w.minus.vec <- ph.ascat.nAraw.w.minus.vec
  new.ph.ascat.nBraw.w.minus.vec <- ph.ascat.nBraw.w.minus.vec
  # replace.idx.vec <- 1:length(ph.ascat.nAraw) %% 2 == 0
  # new.ph.ascat.nAraw[replace.idx.vec] <- ph.ascat.nBraw[replace.idx.vec]
  # new.ph.ascat.nBraw[replace.idx.vec] <- ph.ascat.nAraw[replace.idx.vec]
  # new.ph.ascat.nAraw.w.minus.vec[replace.idx.vec] <- ph.ascat.nBraw.w.minus.vec[replace.idx.vec]
  # new.ph.ascat.nBraw.w.minus.vec[replace.idx.vec] <- ph.ascat.nAraw.w.minus.vec[replace.idx.vec]

  ascat.bc <- ph.ascat.input.list[[samp]]
  keep_idx <- as.numeric(rownames(ascat.bc$Tumor_BAF_segmented[[1]]))
  stopifnot(length(keep_idx) == length(new.ph.ascat.nAraw.w.minus.vec))

  chr.vec <- chrom2integer(as.character(ascat.bc$SNPpos$chr))

  PH_ALL_SAMP_ASCAT_RAW_CPN_DF_LIST[[samp]] <- data.frame(
    chr = chr.vec[keep_idx],
    pos = ascat.bc$SNPpos$pos[keep_idx],
    snp_pos = paste0(chr.vec[keep_idx], "_", ascat.bc$SNPpos$pos[keep_idx]),
    BAF = ascat.bc$Tumor_BAF$BAF[keep_idx],
    LogR = ascat.bc$Tumor_LogR[[1]][keep_idx],
    ph_ascat_nAraw_w_minus = new.ph.ascat.nAraw.w.minus.vec,
    ph_ascat_nBraw_w_minus = new.ph.ascat.nBraw.w.minus.vec,
    ph_ascat_nAraw = new.ph.ascat.nAraw,
    ph_ascat_nBraw = new.ph.ascat.nBraw,
    stringsAsFactors = FALSE
  )
}

# Paths for phased versions of data:
ph.min.cons.segs.df.write.path <- paste0(PLOT_OUTPUT_DIR, patient, "_ph_min_cons_segs_df.txt")
ph.sorted.all.bp.all.samps.segs.anno.df.path <- paste0(PLOT_OUTPUT_DIR, patient, "_ph_sorted_all_bp_all_samps_segs_anno.df.RDS")
ph.min.cons.segs.df.path <- paste0(PLOT_OUTPUT_DIR, patient, "_ph_min_cons_segs_df.RDS")
ph.processed.orig.cpn.segs.df.path <- paste0(PLOT_OUTPUT_DIR, patient, "_ph_processed_orig_cpn_segs_df.RDS")
ph.ploidy.est.vec.path <- paste0(PLOT_OUTPUT_DIR, patient, "_ph_ploidy_est_vec.RDS")
ph.purity.est.vec.path <- paste0(PLOT_OUTPUT_DIR, patient, "_ph_purity_est_vec.RDS")

# New df contain phased ASCAT segmentation for all samples.
ph.init.f.joint.segs.df <- do.call(rbind, ph.ascat.segs.df.list)
ph.f.joint.segs.df <- reformatInitialSegmentationInput(f.joint.segs.df = ph.init.f.joint.segs.df)

# So we're going to drop the singe bp segments that we get here... with a bit of effort they are rescuable
ph.f.joint.segs.df <- ph.f.joint.segs.df[ph.f.joint.segs.df$startpos < ph.f.joint.segs.df$endpos, ]
rownames(ph.f.joint.segs.df) <- 1:nrow(ph.f.joint.segs.df)

ph.sorted.all.bp.all.samps.segs.anno.df <- createFullyAnnotatedJointSegmentation(f.joint.segs.df = ph.f.joint.segs.df, samples = sample.vec)

# we need to update the purity and ploidy for each sample:
ph.purity.est.vec <- c()
ph.ploidy.est.vec <- c()

for (samp in sample.vec) {
  ph.purity.est.vec <- c(ph.purity.est.vec, unique(ph.f.joint.segs.df[ph.f.joint.segs.df$sample == samp, ]$ACF))
  ph.ploidy.est.vec <- c(ph.ploidy.est.vec, unique(ph.f.joint.segs.df[ph.f.joint.segs.df$sample == samp, ]$Ploidy))
}

names(ph.purity.est.vec) <- sample.vec
names(ph.ploidy.est.vec) <- sample.vec

ph.cpn.chr.start.end.df <- ph.sorted.all.bp.all.samps.segs.anno.df[ph.sorted.all.bp.all.samps.segs.anno.df$full_sample == ph.sorted.all.bp.all.samps.segs.anno.df[1, ]$full_sample, c("Chromosome", "Start_Position", "End_Position")]
ph.mock.cpn.matrix <- matrix(data = 0, nrow = nrow(ph.cpn.chr.start.end.df), ncol = length(sample.vec))
colnames(ph.mock.cpn.matrix) <- sample.vec
ph.min.cons.segs.df <- cbind(ph.cpn.chr.start.end.df, ph.mock.cpn.matrix)

ph.seg.name.vec <- paste(ph.min.cons.segs.df[, "Chromosome"], ph.min.cons.segs.df[, "Start_Position"], ph.min.cons.segs.df[, "End_Position"], sep = ":")

ph.asas.SNPs.w.both.lr.baf.snp.lists <- cleanAndErrorCheckSNPsParallel(
  all.samp.SNP.list,
  mcs.seg.to.use.df = ph.min.cons.segs.df,
  samples = sample.vec,
  allow.logr.only.flag = FALSE,
  data.type = data.format.to.use
)

if (length(unique(sapply(ph.asas.SNPs.w.both.lr.baf.snp.lists, length))) != 1) {
  stop("Not all segments are present in all samples' snp assignations, likely something has gone wrong, is the parallelisation working correctly?")
}

# SNPs that only have LogR are permitted in these lists, but it also contains all the SNPs with both BAF and LogR
ph.asas.SNPS.lr.and.or.baf.lists <- cleanAndErrorCheckSNPsParallel(all.samp.SNP.list,
  mcs.seg.to.use.df = ph.min.cons.segs.df,
  samples = sample.vec,
  allow.logr.only.flag = TRUE,
  data.type = data.format.to.use
)

if (!is.na(conf$debug_dir)) {
  saveRDS(ph.asas.SNPs.w.both.lr.baf.snp.lists, file = file.path(conf$debug_dir, "ph_asas_SNPs_w_both_lr_baf_snp_lists.RDS"))
  saveRDS(ph.asas.SNPS.lr.and.or.baf.lists, file = file.path(conf$debug_dir, "ph_asas_SNPS_lr_and_or_baf_lists_lists.RDS"))
}

if (length(unique(sapply(ph.asas.SNPS.lr.and.or.baf.lists, length))) != 1) {
  stop("Not all segments are present in all samples' snp assignations, likely something has gone wrong, is the parallelisation working correctly?")
}

if (length(ph.asas.SNPS.lr.and.or.baf.lists) <= 1) {
  stop("Only one or fewer samples present after reading in data.")
}

message("Getting phasing to apply to all samples using getCombinedBAFDipTestPhasingParallel...")
ph.all.seg.phasing.info.list <- getCombinedBAFDipTestPhasingParallel(
  asas.lr.baf.lists = ph.asas.SNPs.w.both.lr.baf.snp.lists,
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

# Todo unsure of what to do with the baf logr lists which have positions that only have logr here, at the moment exclude them but how
# does that affect later run
message("Generating phased BAF dataframes...")
ph.patient.phased.BAF.df.list <- createAllSegsPhasedBAFDataFrame(
  asas.lr.baf.lists = ph.asas.SNPs.w.both.lr.baf.snp.lists,
  sample.to.phase.seg.list = ph.all.seg.phasing.info.list,
  min.num.hz.snp.phasing = MIN_NUM_HZ_SNP_PHASING
)


# # These unphased copies of the variables are never used
# unph.f.joint.segs.df <- f.joint.segs.df
# unph.sorted.all.bp.all.samps.segs.anno.df <- sorted.all.bp.all.samps.segs.anno.df
# unph.purity.est.vec <- purity.est.vec
# unph.ploidy.est.vec <- ploidy.est.vec
# unph.min.cons.segs.df <- min.cons.segs.df
# unph.seg.name.vec <- seg.name.vec
# unph.asas.SNPs.w.both.lr.baf.snp.lists <- asas.SNPs.w.both.lr.baf.snp.lists
# unph.asas.SNPS.lr.and.or.baf.lists <- asas.SNPS.lr.and.or.baf.lists
# unph.ph.patient.phased.BAF.df.list <- patient.phased.BAF.df.list

# So now we havve generated all input information based on the updated phased
# ASCAT run, let's reset the vairable names to allow the subclonality code to
# continue to work:
f.joint.segs.df <- ph.f.joint.segs.df
sorted.all.bp.all.samps.segs.anno.df <- ph.sorted.all.bp.all.samps.segs.anno.df
purity.est.vec <- ph.purity.est.vec
ploidy.est.vec <- ph.ploidy.est.vec
min.cons.segs.df <- ph.min.cons.segs.df
seg.name.vec <- ph.seg.name.vec
asas.SNPs.w.both.lr.baf.snp.lists <- ph.asas.SNPs.w.both.lr.baf.snp.lists
asas.SNPS.lr.and.or.baf.lists <- ph.asas.SNPS.lr.and.or.baf.lists
patient.phased.BAF.df.list <- ph.patient.phased.BAF.df.list


# Detect MSAI in the original BAF before it is processed alongside logr to get the Raw allele specifc cpn to be used in the rest of the analysis.
message("************************************************************************")
message("Performing LTNPD MSAI detection on the BAF...")
message("************************************************************************")
all.seg.MSAI.by.BAF.ltnpd.present.list <- getAllSegsMSAIDetectionFromBAFListLTNPD(
  seg.name.vec = seg.name.vec,
  sample.vec = sample.vec,
  asas.SNPS.lr.and.or.baf.lists = asas.SNPS.lr.and.or.baf.lists,
  all.seg.phasing.info.list = all.seg.phasing.info.list,
  MSAI.BAF.LTPND.pval.thresh = MSAI_BAF_LTNPD_PVAL_THRESHOLD
)
# saveRDS(ph.all.seg.MSAI.by.BAF.ltnpd.present.list, file = ph.all.seg.MSAI.by.BAF.ltnpd.present.list.path)


if (isTRUE(PLOTTING_FLAG)) {
  png(file = file.path(PLOT_OUTPUT_DIR, "MSAI_detection_BAF_phased.png"), width = 1200, height = 1400)

  message("Plotting phased BAF across the genome...")
  plotAllSampleBAFsAndMSAIAcrossGenome(
    sample.vec = sample.vec,
    all.seg.MSAI.by.BAF.ltnpd.present.list = all.seg.MSAI.by.BAF.ltnpd.present.list,
    patient.phased.BAF.df.list = patient.phased.BAF.df.list,
    cpnA.allele.color = "#de2d2699",
    cpnB.allele.color = "#3182bd99",
    use.xphase = USING_XPHASE
  )
  dev.off()
} else {
  message("Not plotting phased BAF as flag is FALSE")
}

message("Applying LTNPD MSAI detection to the original segments...")
# So now we have the original ASCAT segmentation and we want to make the MEDICC input from it
# So now we want to flip segments where MSAI has been detected and they were separated in the original data:
# So we also want to flip the raw values:
intersection.orig.mcs.df <- sorted.all.bp.all.samps.segs.anno.df[sorted.all.bp.all.samps.segs.anno.df$in_intersection == TRUE, ]
intersection.seg.name.vec <- paste(intersection.orig.mcs.df$Chromosome, intersection.orig.mcs.df$Start_Position, intersection.orig.mcs.df$End_Position, sep = ":")
orig.intersection.mcs.w.MSAI.df <- intersection.orig.mcs.df
unq.intersection.seg.name.vec <- unique(intersection.seg.name.vec)
for (seg.name in unq.intersection.seg.name.vec) {
  all.samp.matching.seg.vec <- paste0(sample.vec, "@", seg.name)

  if (isTRUE(all.seg.MSAI.by.BAF.ltnpd.present.list[[seg.name]]$MSAI_detection_feasible)) {
    seg.MPhase.MSAI.df <- all.seg.MSAI.by.BAF.ltnpd.present.list[[seg.name]]$all_samps_BAF_and_ltnpd_df
    seg.MPhase.MSAI.ASCAT.df <- cbind(seg.MPhase.MSAI.df, intersection.orig.mcs.df[paste0(sample.vec, ":", seg.name), ])
    seg.MPhase.MSAI.ASCAT.df$samp_seg_name <- paste0(sample.vec, ":", seg.name)
    seg.MPhase.MSAI.ASCAT.df$balanced <- seg.MPhase.MSAI.ASCAT.df$nMinor == seg.MPhase.MSAI.ASCAT.df$nMajor
    AI.seg.MPhase.MSAI.ASCAT.df <- seg.MPhase.MSAI.ASCAT.df[seg.MPhase.MSAI.ASCAT.df$balanced == FALSE, ]

    # Need to check MSAI is ocurring in the original segmentation allelic imbalance calls
    # Why checking both of these... FIXME...
    if (FALSE %in% AI.seg.MPhase.MSAI.ASCAT.df$MSAI_pval_passed_thresh & TRUE %in% AI.seg.MPhase.MSAI.ASCAT.df$MSAI_pval_passed_thresh) {
      switched.samp.seg.names.vec <- AI.seg.MPhase.MSAI.ASCAT.df[AI.seg.MPhase.MSAI.ASCAT.df$MSAI_pval_passed_thresh == TRUE & AI.seg.MPhase.MSAI.ASCAT.df$poss_dominant_allele == "cpnB", ]$samp_seg_name
      old.nMinor.vec <- orig.intersection.mcs.w.MSAI.df[switched.samp.seg.names.vec, ]$nMinor
      # old.nBraw.vec <- orig.intersection.mcs.w.MSAI.df[switched.samp.seg.names.vec,]$nBraw
      # orig.intersection.mcs.w.MSAI.df[switched.samp.seg.names.vec,]$nBraw <-  intersection.orig.mcs.df[switched.samp.seg.names.vec,]$nAraw
      orig.intersection.mcs.w.MSAI.df[switched.samp.seg.names.vec, ]$nMinor <- intersection.orig.mcs.df[switched.samp.seg.names.vec, ]$nMajor
      orig.intersection.mcs.w.MSAI.df[switched.samp.seg.names.vec, ]$nMajor <- old.nMinor.vec
    }
  }
}

if (!is.na(conf$debug_dir)) {
  saveRDS(orig.intersection.mcs.w.MSAI.df, file = file.path(conf$debug_dir, "orig_intersection_mcs_w_MSAI_df.rds"))
}

pre.MPhase.MEDICC.FASTA.file.list <- makeMEDICCInputListsFromSegmentationDF(
  segs.to.make.MEDICC.input.df = orig.intersection.mcs.w.MSAI.df,
  sample.vec = sample.vec,
  chrs.to.use.vec = CHRS_TO_USE_VEC,
  chr.col = "Chromosome",
  start.col = "Start_Position",
  end.col = "End_Position",
  cpnA.col = "nMajor",
  cpnB.col = "nMinor"
)

if (!is.na(conf$debug_dir)) {
  saveRDS(pre.MPhase.MEDICC.FASTA.file.list, file = file.path(conf$debug_dir, "pre_MPhase_MEDICC_FASTA_file_list.rds"))
}


# Todo add another method of detecting MSAI using only the BAF and a paired T test.
# all.seg.MSAI.by.BAF.ltnpd.present.list <- getAllSegsMirroredSubclonalAllelicImbalanceDetectionFromBAFList(seg.name.vec = seg.name.vec,
# 																										   sample.vec  = sample.vec,
# 																										   asas.SNPS.lr.and.or.baf.lists = asas.SNPS.lr.and.or.baf.lists,
# 																										   all.seg.phasing.info.list = all.seg.phasing.info.list,
# 																										   MSAI.BAF.LTPND.pval.thresh = MSAI_BAF_LTNPD_PVAL_THRESHOLD)

message("Generating datapoints with LogR total cpn dataframes...")
# This is where the SNP filtering occurs. We can fix this by adding by passing an option of which data type is being analysed...
all.samp.non.hz.logr.and.tot.cpn.df.list <- getAllSampLogrAlonePositionsLogrAndTotalCPNDFsList(
  sample.vec = sample.vec,
  asas.SNPS.lr.and.or.baf.lists = asas.SNPS.lr.and.or.baf.lists,
  ploidy.est.vec = ploidy.est.vec,
  purity.est.vec = purity.est.vec,
  gamma.val = GAMMA_VAL,
  data.format.to.use = data.format.to.use
)
if (!is.na(conf$debug_dir)) {
  saveRDS(all.samp.non.hz.logr.and.tot.cpn.df.list, file = all.samp.non.hz.logr.and.tot.cpn.df.list.path)
}


# TBKW edit, so this is where we want to update the raw copy number to reflect what
# ASCAT has outputted
# Add this the unphaseable stuff too maybe? Not necessary.
message("Applying phasing to all samples to get phased raw cpn...")
asas.ph.raw.cpn.lists <- getAllSegAllSamplePhasedRawCPNList(
  sample.vec = sample.vec,
  ploidy.est.vec = ploidy.est.vec,
  purity.est.vec = purity.est.vec,
  min.cons.segs.df = min.cons.segs.df,
  all.seg.phasing.info.list = all.seg.phasing.info.list,
  asas.both.lr.baf.lists = asas.SNPs.w.both.lr.baf.snp.lists,
  gamma.val = GAMMA_VAL,
  min.num.hz.snp.phasing = MIN_NUM_HZ_SNP_PHASING,
  use.xphase = USING_XPHASE,
  ascat.raw.cpn = PH_ALL_SAMP_ASCAT_RAW_CPN_DF_LIST
)
if (!is.na(conf$debug_dir)) {
  saveRDS(asas.ph.raw.cpn.lists, file = asas.ph.raw.cpn.lists.path)
}

# Let's get the CCF
message("Checking all minimum consistent segments are present in each sample...")
seg.name.vec <- names(asas.ph.raw.cpn.lists[[1]])

# Check to make sure that all segments are present in all samples... might not be necessary todo, check.
for (samp in sample.vec) {
  seg.name.vec <- intersect(seg.name.vec, names(asas.ph.raw.cpn.lists[[samp]]))
} # end for ( samp in sample.vec ) {

if (length(seg.name.vec) == 0) {
  stop("No shared segments for some reason, go back and check")
}

ccf_data <- ccf_from_raw_cpn(asas.ph.raw.cpn.lists, sample.vec, seg.name.vec)
all.seg.CCF.list <- ccf_data$all.seg.CCF.list
all.seg.CCF.diff.from.exp.list <- ccf_data$all.seg.CCF.diff.from.exp.list

# Make a table that contains all the relevant information: subclonal calls... gain and loss vs ploidy etc etc etc
full.joint.segs.df <- add_info_to_output_table(
  asas.ph.raw.cpn.lists,
  ploidy.est.vec, purity.est.vec,
  all.samp.non.hz.logr.and.tot.cpn.df.list,
  patient.phased.BAF.df.list,
  asas.SNPs.w.both.lr.baf.snp.lists,
  all.seg.CCF.list,
  all.seg.CCF.diff.from.exp.list,
  all.seg.phasing.info.list,
  all.seg.MSAI.by.BAF.ltnpd.present.list
)

# Remove segments without any datapoints.
# Todo, merge these int
# print("Removing any segments that do not have underlying datapoints, this is an artefact of the minimum consistent segmentation performed by MPhase...")
# full.joint.segs.df <- full.joint.segs.df[full.joint.segs.df$num_un_ph_hz_sites != 0 | full.joint.segs.df$num_ph_hz_sites != 0  | full.joint.segs.df$num_non_hz_sites != 0,]
# print("Removing segments that do not map back to any original segments")
# full.joint.segs.df <- full.joint.segs.df[full.joint.segs.df$orig_seg_name != "NA:NA:NA",]
saveRDS(full.joint.segs.df, file = file.path(PLOT_OUTPUT_DIR, "full_joint_segs.rds"))

if (isTRUE(PLOTTING_FLAG)) {
  message("Plotting CCFs for cpnA and cpnB across the genome using cpn definition of subclonality...")

  png(file = file.path(PLOT_OUTPUT_DIR, "CCF_w_cpn_thresh_cpnA.png"), width = 1200, height = 1400)
  plotAllSampleCCFsAcrossGenome(
    asas.ph.raw.cpn.lists = asas.ph.raw.cpn.lists,
    full.joint.segs.df = full.joint.segs.df,
    sample.vec = sample.vec,
    seg.name.vec = seg.name.vec,
    all.seg.CCF.vals.list = all.seg.CCF.list,
    cpnX = "cpnA",
    use.cpn.or.ccf.thresh = "CPN"
  )
  dev.off()

  png(file = file.path(PLOT_OUTPUT_DIR, "CCF_w_cpn_thresh_cpnB.png"), width = 1200, height = 1400)
  plotAllSampleCCFsAcrossGenome(
    asas.ph.raw.cpn.lists = asas.ph.raw.cpn.lists,
    full.joint.segs.df = full.joint.segs.df,
    sample.vec = sample.vec,
    seg.name.vec = seg.name.vec,
    all.seg.CCF.vals.list = all.seg.CCF.list,
    cpnX = "cpnB",
    use.cpn.or.ccf.thresh = "CPN"
  )
  dev.off()

} else {
  message("Not plotting CCFs for cpnA and cpnB across the genome using cpn definition of subclonality as plotting flag is FALSE.")
}

# We need a .tsv version of this data for interoperability with python
write.table(full.joint.segs.df, file = file.path(PLOT_OUTPUT_DIR, "full_joint_segs.tsv"), row.names = FALSE, sep = "\t")

# Write out a simplified version of the segmentation
# From discussion with Tom:
#   Also, just to be 100% certain, for mphase you use nMajor, nMinor and
#   cnTotal, and for xphase you use orig_nMajor, orig_nMinor, and orig_cnTotal?
#   And then flip the numbers if dominant_allele == cpnB
#simple_segs <- full.joint.segs.df[!is.na(full.joint.segs.df$orig_nMajor) & full.joint.segs.df$samp_used_to_phase != "too_few_SNPs",
#simple_segs <- full.joint.segs.df[!is.na(full.joint.segs.df$n.het),
simple_segs <- full.joint.segs.df[!is.na(full.joint.segs.df$cnTotal) & !is.na(full.joint.segs.df$orig_cnTotal),
  c("sample", "chr", "startpos", "endpos", "orig_cnTotal", "orig_nMajor", "orig_nMinor", "dominant_allele", "samp_used_to_phase", "Ploidy", "ACF")]
# Because the dominant allele is only set when there is allelic imbalance, we
# force it to a value (= "no_ai") here so that we don't have to deal with NAs
# in this column
simple_segs$dominant_allele[is.na(simple_segs$dominant_allele)] <- "no_ai"
simple_segs$cpn_total <- simple_segs$orig_cnTotal
simple_segs$cpn_a <- ifelse(simple_segs$dominant_allele == "cpnA", simple_segs$orig_nMajor, simple_segs$orig_nMinor)
simple_segs$cpn_b <- ifelse(simple_segs$dominant_allele == "cpnA", simple_segs$orig_nMinor, simple_segs$orig_nMajor)
simple_segs$is_msai <- simple_segs$dominant_allele == "cpnB"
simple_segs$is_reference <- simple_segs$sample == simple_segs$samp_used_to_phase
simple_segs$ploidy <- simple_segs$Ploidy
simple_segs$purity <- simple_segs$ACF
# Remove some columns we don't need anymore
simple_segs[, c("orig_cnTotal", "orig_nMajor", "orig_nMinor", "dominant_allele", "samp_used_to_phase", "Ploidy", "ACF")] <- NULL
write.table(simple_segs, file = file.path(PLOT_OUTPUT_DIR, "segmentation.tsv"), row.names = FALSE, sep = "\t", quote = FALSE)

# Write an extra table with sample-wide information like purity and ploidy
sample_info <- unique(full.joint.segs.df[, c("sample", "Ploidy", "ACF")])
colnames(sample_info)[colnames(sample_info) == "Ploidy"] <- "ploidy"
colnames(sample_info)[colnames(sample_info) == "ACF"] <- "purity"
write.table(sample_info, file = file.path(PLOT_OUTPUT_DIR, "sample_info.tsv"), row.names = FALSE, sep = "\t", quote = FALSE)

if (!is.na(conf$medicc_output_dir)) {
  # Hack! Hope this works...
  full.joint.segs.nona <- full.joint.segs.df[!is.na(full.joint.segs.df$ph_int_cpnA), ]
  post.MPhase.MCS.MEDICC.FASTA.files.list <- makeMEDICCInputListsFromSegmentationDF(
    segs.to.make.MEDICC.input.df = full.joint.segs.nona,
    sample.vec = sample.vec,
    chrs.to.use.vec = CHRS_TO_USE_VEC,
    chr.col = "chrom",
    start.col = "startpos",
    end.col = "endpos",
    cpnA.col = "ph_int_cpnA",
    cpnB.col = "ph_int_cpnB"
  )

  if (!is.na(conf$debug_dir)) {
    saveRDS(post.MPhase.MCS.MEDICC.FASTA.files.list, file = file.path(conf$debug_dir, "post_MPhase_MCS_MEDICC_FASTA_file_list.rds"))
  }

  createMEDICCInputFiles(post.MPhase.MCS.MEDICC.FASTA.files.list, sample.vec, conf$medicc_output_dir)
}
