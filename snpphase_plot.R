#!/usr/bin/env Rscript

library(digest)
library(GenomicRanges)
library(optparse)

option_list <- list(
    make_option(c("-s", "--original-segmentation"), default="",
                help="Pre-refphase segmentation (from sequenza)"),
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
orig_segs <- read.table(conf$original_segmentation, stringsAsFactors = FALSE, sep = "\t", header = TRUE)
orig_segs$seg_center <- orig_segs$startpos + (orig_segs$endpos - orig_segs$startpos) / 2

baf_logr <- readRDS(file.path(conf$input_dir, paste0(conf$patient_id, "_phased_BAF_df_list.RDS")))
segs <- readRDS(file.path(conf$input_dir, paste0(conf$patient_id, "_full_joint_segs_df.RDS")))
segs$seg_center <- segs$startpos + (segs$endpos - segs$startpos) / 2
segs$seg_idx <- seq(nrow(segs))
segs$is_reference <- segs$samp_used_to_phase == segs$sample

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

num_subfigures <- 5  # BAF and logR subfigures only + copy numbers

# Split the output figure into one row per sample
samples <- names(baf_logr)

output_plot <- file.path(conf$output_dir, paste0("all_chromosomes.pdf"))
pdf(output_plot, width = 10, height = 1 + 4 * length(samples))
for (chrom_to_plot in names(chrom_starts)) {
  par(cex = 0.8)
  # Header: Leave space for the chromosome name
  # coordinates are: left, right, bottom, and top
  #heading_height_fraction <- 1 / (length(samples) + 0.5)
  heading_height_fraction <- 1 / (1 + 5 * length(samples))
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
    orig_segs_sample <- orig_segs[orig_segs$sample == sample_name, ]

    # Split the plotting area vertically (sample name + info on left, plots on right)
    screen(sample_screens[[sample_idx]])
    split_coords <- matrix(c(0, 0.1, 0, 1, 0.1, 1, 0, 1), ncol = 4, byrow = TRUE)
    screen_idxs <- split.screen(split_coords)

    screen(screen_idxs[[1]])
    plot_sample_info(sample_name, segs_sample$Ploidy[[1]], segs_sample$ACF[[1]])

    screen(screen_idxs[[2]])
    sub_split_coords <- matrix(c(0, 1, 0.65, 1,
                                 0, 1, 0.3, 0.65,
                                 0, 1, 0.2, 0.3,
                                 0, 1, 0.1, 0.2,
                                 0, 1, 0, 0.1), ncol = 4, byrow = TRUE)
    subfigure_idx <- split.screen(sub_split_coords)
    #subfigure_idx <- split.screen(c(num_subfigures, 1))

    # Remove the many NA regions
    segs_sample <- segs_sample[!is.na(segs_sample$cnTotal), ]

    # Keep data for this sample
    data_sample <- data_sample[data_sample$chrom == chrom_to_plot, ]
    segs_sample <- segs_sample[segs_sample$chr == chrom_to_plot, ]
    orig_segs_sample <- orig_segs_sample[orig_segs_sample$chr == chrom_to_plot, ]

    # ------------------------------------------------------------------------
    # Calculate which segments have been changed by refphase
    gr_orig <- makeGRangesFromDataFrame(orig_segs_sample,
                                        keep.extra.columns = TRUE,
                                        ignore.strand = FALSE,
                                        seqnames.field = "chr",
                                        start.field = "startpos",
                                        end.field = "endpos")

    stopifnot(length(gr_orig) == nrow(orig_segs_sample))

    gr_new <- makeGRangesFromDataFrame(segs_sample,
                                       keep.extra.columns = TRUE,
                                       ignore.strand = FALSE,
                                       seqnames.field = "chr",
                                       start.field = "startpos",
                                       end.field = "endpos")

    stopifnot(length(gr_new) == nrow(segs_sample))

    ov <- findOverlaps(gr_new, gr_orig)

    new_cns <- mcols(gr_new)[from(ov), c("orig_cnTotal", "orig_nMinor", "orig_nMajor")]
    colnames(new_cns) <- c("cnTotal", "nMinor", "nMajor")
    orig_cns <- mcols(gr_orig)[to(ov), c("cnTotal", "nMinor", "nMajor")]

    str(orig_cns)
    str(new_cns)

    na_row <- apply(is.na(orig_cns), 1, any) | apply(is.na(new_cns), 1, any)

    orig_cns <- orig_cns[!na_row,]
    new_cns <- new_cns[!na_row,]

    str(orig_cns)
    str(new_cns)

    write.csv(cbind(ov[!na_row, ], new_cns, orig_cns), file="all_cns.csv")

    seg_changed <- apply(as.matrix(orig_cns) != as.matrix(new_cns), 1, any)
    message("# of segments changed=", sum(seg_changed))
    message("# of segments not changed=", sum(!seg_changed))

    new_changed <- by(seg_changed, from(ov)[!na_row], any)
    str(new_changed)
    message("# of segments changed=", sum(unlist(new_changed)))

    changed_idx <- mcols(gr_new)[from(ov)[!na_row],"seg_idx"][new_changed]
    print(changed_idx)
    segs_sample$is_changed <- FALSE
    segs_sample$is_changed[segs_sample$seg_idx %in% changed_idx] <- TRUE

    str(segs_sample)
    # ------------------------------------------------------------------------

    # logR track
    screen(subfigure_idx[[1]])
    plot_track_logr(data_sample, segs_sample)

    # BAF track
    screen(subfigure_idx[[2]])
    plot_track_baf(data_sample, segs_sample)

    # Orig CNs
    screen(subfigure_idx[[3]])
    plot_track_cn(orig_segs_sample, "orig")

    # Refphase CNs
    screen(subfigure_idx[[4]])
    plot_track_cn(data.frame(startpos = segs_sample$startpos,
                             endpos = segs_sample$endpos,
                             seg_center = segs_sample$seg_center,
                             cnTotal = segs_sample$orig_cnTotal,
                             nMajor = segs_sample$orig_nMajor,
                             nMinor = segs_sample$orig_nMinor,
                             is_changed = segs_sample$is_changed),
                  track_name = "new")

    # Add an x-axis with genomic co-ordinates
    screen(subfigure_idx[[5]])
    plot_track_axis(min(segs_sample$startpos), max(segs_sample$endpos))

  }

  close.screen(all = TRUE)
}
dev.off()

# Whole genome plot
output_plot <- file.path(conf$output_dir, "whole_genome.pdf")
pdf(output_plot, width = 10, height = 3 * length(samples))
par(cex = 0.8)
split.screen(c(length(samples), 1))
num_subfigures <- 2

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
