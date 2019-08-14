#!/usr/bin/env Rscript

library(optparse)
library(digest)

option_list <- list(
    make_option(c("-i", "--input-dir"), default="",
                help="Input directory containing the output of Mphase"),
    make_option(c("-p", "--patient-id"), default="",
                help="Patient ID"),
    make_option(c("-o", "--output-dir"), default="",
                help="Input directory containing the output of Mphase"),
    make_option(c("-f", "--plot-format"), default="png",
                help="Input directory containing the output of Mphase"),
    make_option(c("-d", "--seed"), default=531,
                help="Set the random seed at the beginning of the script")
    )
args <- parse_args(OptionParser(option_list = option_list), convert_hyphens_to_underscores = TRUE)
conf <- args
conf$id <- substr(digest(args), 1, 8)
set.seed(args$seed)

source("install/src/multiregioncopynum/snpphase_plot_utils.R")

if (!dir.exists(conf$output_dir))
  dir.create(conf$output_dir, recursive = TRUE, showWarnings = FALSE)

# Read in all the data we need to plot
baf_logr <- readRDS(file.path(conf$input_dir, paste0(conf$patient_id, "_phased_BAF_df_list.RDS")))
segs <- readRDS(file.path(conf$input_dir, paste0(conf$patient_id, "_full_joint_segs_df.RDS")))
msai <- readRDS(file.path(conf$input_dir, paste0(conf$patient_id, "_all_seg_MSAI_by_BAF_ltnpd_present_list.RDS")))
str(baf_logr)
str(segs)
str(msai)
segs$seg_center <- segs$startpos + (segs$endpos - segs$startpos) / 2

# Convert chromosomal positions to absolute positions for plotting the whole genome
chrom_ends <- cumsum(by(segs$endpos, segs$chr, max))
chrom_starts <- c(0, chrom_ends[-length(chrom_ends)])
names(chrom_starts) <- names(chrom_ends)
chrom_center <- chrom_starts + (chrom_ends - chrom_starts) / 2
segs$chrom_startpos <- chrom_starts[segs$chr]

# Chromosome start and end positions in genome-wide absolute positions
chrom_info <- data.frame(start = chrom_starts,
                         end = chrom_ends,
                         center = chrom_center,
                         label = names(chrom_starts))

num_subfigures <- 2  # BAF and logR subfigures only
#chrom_to_plot <- 1
#chrom_to_plot <- NA  # NA == whole genome

# Split the output figure into one row per sample
#split.screen(c(length(samples) * num_subfigures, 1))
samples <- names(baf_logr)
#dir.create(conf$output_dir, "by-chromosome")

output_plot <- file.path(conf$output_dir, paste0("all_chromosomes.pdf"))
pdf(output_plot, width = 10, height = 3 * length(samples))
for (chrom_to_plot in names(chrom_starts)) {
  par(cex = 0.8)
  # Header: Leave space for the chromosome name
  # coordinates are: left, right, bottom, and top
  #heading_height_fraction <- 1 / (length(samples) + 0.5)
  heading_height_fraction <- 0.05
  #split_coords <- matrix(c(0, 1, 0, 1 - heading_height_fraction, 0, 1, 0, heading_height_fraction), ncol = 4, byrow = TRUE)
  split_coords <- matrix(c(0, 1, 1 - heading_height_fraction, 1, 0, 1, 0, 1 - heading_height_fraction), ncol = 4, byrow = TRUE)
  outer_screens <- split.screen(split_coords)
  screen(outer_screens[[1]])
  # Plot header
  plot_chromosome_header(conf$patient_id, paste0("chr", chrom_to_plot))
  screen(outer_screens[[2]])

  # Split horizontally into figures for each sample
  sample_screens <- split.screen(c(length(samples), 1))
  for (sample_idx in seq_along(samples)) {
    sample_name <- samples[[sample_idx]]
    data_sample <- baf_logr[[sample_idx]]
    segs_sample <- segs[segs$sample == sample_name, ]

    # Split the plotting area vertically (sample name + info on left, plots on right)
    screen(sample_screens[[sample_idx]])
    split_coords <- matrix(c(0, 0.1, 0, 1, 0.1, 1, 0, 1), ncol = 4, byrow = TRUE)
    screen_idxs <- split.screen(split_coords)

    screen(screen_idxs[[1]])
    plot_sample_info(sample_name, segs_sample$Ploidy[[1]], segs_sample$ACF[[1]])

    screen(screen_idxs[[2]])
    subfigure_idx <- split.screen(c(num_subfigures, 1))

    # Remove the many NA regions
    segs_sample <- segs_sample[!is.na(segs_sample$cnTotal), ]

    # Keep data for this sample
    data_sample <- data_sample[data_sample$chrom == chrom_to_plot, ]
    segs_sample <- segs_sample[segs_sample$chr == chrom_to_plot, ]

    # logR track
    screen(subfigure_idx[[1]])
    plot_track_logr(data_sample, segs_sample)

    # BAF track
    screen(subfigure_idx[[2]])
    plot_track_baf(data_sample, segs_sample)

  }

  close.screen(all = TRUE)
}
dev.off()

# Whole genome plot
output_plot <- file.path(conf$output_dir, "whole_genome.pdf")
pdf(output_plot, width = 10, height = 3 * length(samples))
par(cex = 0.8)
split.screen(c(length(samples), 1))

for (sample_idx in seq_along(samples)) {
  sample_name <- samples[[sample_idx]]
  data_sample <- baf_logr[[sample_idx]]
  segs_sample <- segs[segs$sample == sample_name, ]

  # Split the plotting area vertically (sample name + info on left, plots on right)
  screen(sample_idx)
  split_coords <- matrix(c(0, 0.1, 0, 1, 0.1, 1, 0, 1), ncol = 4, byrow = TRUE)
  screen_idxs <- split.screen(split_coords)

  screen(screen_idxs[[1]])
  plot_sample_info(sample_name, segs_sample$Ploidy[[1]], segs_sample$ACF[[1]])

  screen(screen_idxs[[2]])
  subfigure_idx <- split.screen(c(num_subfigures, 1))

  # Remove the many NA regions
  segs_sample <- segs_sample[!is.na(segs_sample$cnTotal), ]

  # Adjust all positions to genome-wide absolute positions
  data_sample$pos <- data_sample$pos + chrom_starts[data_sample$chrom]
  segs_sample$startpos <- segs_sample$startpos + segs_sample$chrom_startpos
  segs_sample$endpos <- segs_sample$endpos + segs_sample$chrom_startpos
  segs_sample$seg_center <- segs_sample$seg_center + segs_sample$chrom_startpos

  # logR track
  screen(subfigure_idx[[1]])
  plot_track_logr_whole_genome(data_sample, segs_sample, chrom_info)

  # BAF track
  screen(subfigure_idx[[2]])
  plot_track_baf_whole_genome(data_sample, segs_sample, chrom_info)

}

close.screen(all = TRUE)
dev.off()
