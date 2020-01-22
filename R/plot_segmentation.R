#!/usr/bin/env Rscript
#
# Plot the original ASCAT segmentation as well as the post-mphase segmentation

library(optparse)
library(digest)

option_list <- list(
  make_option(c("-p", "--patient"),
    default = "CB1002",
    help = "Patient identifier"
  ),
  make_option(c("-i", "--input-dir"),
    default = "../output/CB1002/mphase",
    help = "Input and output dir"
  ),
  make_option(c("-d", "--seed"),
    default = 531,
    help = "Set the random seed at the beginning of the script"
  )
)
args <- parse_args(OptionParser(option_list = option_list), convert_hyphens_to_underscores = TRUE)
conf <- args
conf$id <- substr(digest(args), 1, 8)
set.seed(args$seed)


LocationOfThisScript <- function() { # Function LocationOfThisScript returns the location of this .R script (may be needed to source other files in same dir)
  this.file <- NULL
  # This file may be 'sourced'
  for (i in -(1:sys.nframe())) {
    if (identical(sys.function(i), base::source)) this.file <- (normalizePath(sys.frame(i)$ofile))
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

current.dir <- LocationOfThisScript()

source(file.path(current.dir, "cpn_heterogeneity_plotting_functions_v2.R"))

patient <- conf$patient
PLOT_OUTPUT_DIR <- paste0(conf$input_dir, "/")

processed.orig.cpn.segs.df.path <- paste0(PLOT_OUTPUT_DIR, patient, "_processed_orig_cpn_segs_df.RDS")
f.joint.segs.df <- readRDS(file = processed.orig.cpn.segs.df.path)

full.joint.segs.df.path <- paste0(PLOT_OUTPUT_DIR, patient, "_full_joint_segs_df.RDS")
full.joint.segs.df <- readRDS(file = full.joint.segs.df.path)

plotCN <- function(seg, samples, plot_cols = c("int_ph_cpnA", "int_ph_cpnB")) {
  seg <- seg[seg$chr %in% 1:22, ]
  plotCPNSegDFAccrossGenome(segs.df = seg, sample.vec = samples, cpn.cols.to.plot = plot_cols, color.vec = c("green", "black"))
}

sample.vec <- unique(f.joint.segs.df$sample)

png(file = paste0(PLOT_OUTPUT_DIR, patient, "_orig_ASCAT_across_genome.png"), width = 1200, height = 1400)
plotCN(f.joint.segs.df, sample.vec, plot_cols = c("nMajor", "nMinor"))
dev.off()

png(file = paste0(PLOT_OUTPUT_DIR, patient, "_mphase_mcs_cpn_across_genome.png"), width = 1200, height = 1400)
plotCN(full.joint.segs.df, sample.vec, plot_cols = c("ph_int_cpnA", "ph_int_cpnB"))
dev.off()
