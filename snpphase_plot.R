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
    make_option(c("-d", "--seed"), default=531,
                help="Set the random seed at the beginning of the script")
    )
args <- parse_args(OptionParser(option_list = option_list), convert_hyphens_to_underscores = TRUE)
conf <- args
conf$id <- substr(digest(args), 1, 8)
set.seed(args$seed)

plot_track_logr <- function(data_sample, segs_sample) {
  par(mar = c(1, 5, 1, 0.1), xaxs = "i")
  plot(data_sample$pos, data_sample$logR, ylab = "logR", xaxt = "n", col = adjustcolor("grey", 0.6), pch = 20)
  #lines(x = data_sample$pos, y = data_sample$logRsegmented, col=ifelse(is.na(data_sample$logRsegmented), "white", "black"))
  segments(x0 = segs_sample$startpos, y0 = segs_sample$median_non_hz_logr, x1 = segs_sample$endpos, col = "black")
  abline(v = segs_sample$startpos, lty = 2, col = "black")
  text(x = segs_sample$seg_center, y = 0, labels = segs_sample$cnTotal)
}

plot_track_baf <- function(data_sample, segs_sample) {
  par(mar = c(2, 5, 0.1, 0.1), xaxs = "i")
  phasing_col <- ifelse(data_sample$maj_is_cpnB_logic_vec == 1, "blue", "orange")
  #phasing_col <- ifelse(data_sample$, "blue", "yellow")
  phasing_col[data_sample$sample_phased_by_this_sample == "none_no_allelic_imbalance_called"] <- "grey"
  phasing_col <- adjustcolor(phasing_col, 0.6)
  plot(data_sample$pos, data_sample$BAF, ylab = "BAF", col = phasing_col, ylim = c(0, 1), pch = 20)
  #lines(x = data_sample$pos, y = data_sample$BAFsegmented, col="black")

  # Mean segments
  mean_col <- ifelse(segs_sample$within_samp_cl_AI_present, "black", NA)
  segments(x0 = segs_sample$startpos, y0 = segs_sample$mean_BAF + 0.5, x1 = segs_sample$endpos, col = mean_col)
  segments(x0 = segs_sample$startpos, y0 = -1 * segs_sample$mean_BAF + 0.5, x1 = segs_sample$endpos, col = mean_col)

  # Segment boundaries
  abline(v = segs_sample$startpos, lty = 2, col = "black")

  # Original copy numbers
  text(x = segs_sample$seg_center, y = 0.9, labels = segs_sample$nMajor)
  text(x = segs_sample$seg_center, y = 0.1, labels = segs_sample$nMinor)

  # Changed copy numbers
  changed_cn <- segs_sample$nMinor_change | segs_sample$nMajor_change
  changed_cn[is.na(changed_cn)] <- FALSE
  if (sum(changed_cn, na.rm = TRUE) > 0) {
    message("Found some adjusted copy numbers!")
    text(x = segs_sample$seg_center[changed_cn], y = 0.6, labels = segs_sample$orig_nMajor[changed_cn])
    text(x = segs_sample$seg_center[changed_cn], y = 0.4, labels = segs_sample$orig_nMinor[changed_cn])
  }
}

plot_sample_info <- function(sample, ploidy, cellularity) {
  par(mar = c(0,0,0,0))
  plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')

  text(x = 0.5, y = 0.5, adj = 0.5,
       labels = paste0(sample, "\n", "cellularity:", round(cellularity, 2), "  ploidy:", round(ploidy, 2)),
       cex = 1.3, col = "black", srt = 90)
}

# Read in all the data we need to plot
baf_logr <- readRDS(file.path(conf$input_dir, paste0(conf$patient_id, "_phased_BAF_df_list.RDS")))
segs <- readRDS(file.path(conf$input_dir, paste0(conf$patient_id, "_full_joint_segs_df.RDS")))
msai <- readRDS(file.path(conf$input_dir, paste0(conf$patient_id, "_all_seg_MSAI_by_BAF_ltnpd_present_list.RDS")))
str(baf_logr)
str(segs)
segs$seg_center <- segs$startpos + (segs$endpos - segs$startpos) / 2

samples <- names(baf_logr)
output_plot <- file.path(conf$output_dir, "whole_genome.pdf")
pdf(output_plot, width=10, height=15)
par(cex = 0.7)

num_subfigures <- 2  # BAF and logR subfigures only
chrom_to_plot <- 9

# Split the output figure into one row per sample
#split.screen(c(length(samples) * num_subfigures, 1))
split.screen(c(length(samples), 1))

for (sample_idx in seq_along(samples)) {
  sample <- samples[[sample_idx]]
  data_sample <- baf_logr[[sample_idx]]
  segs_sample <- segs[segs$sample == sample,]

  # Split the plotting area vertically (sample name + info on left, plots on right)
  screen(sample_idx)
  #screen_idxs <- split.screen(matrix(c(0, 0.1, 1, 1, 0.1, 1, 1, 1), ncol = 4, byrow = TRUE))
  split_coords <- matrix(c(0, 0.1, 0, 1, 0.1, 1, 0, 1), ncol = 4, byrow = TRUE)
  screen_idxs <- split.screen(split_coords)

  screen(screen_idxs[[1]])
  plot_sample_info(sample, segs_sample$Ploidy[[1]], segs_sample$ACF[[1]])

  screen(screen_idxs[[2]])
  subfigure_idx <- split.screen(c(num_subfigures, 1))

  # Remove the many NA regions
  segs_sample <- segs_sample[!is.na(segs_sample$cnTotal),]

  # Just plot one chromosome for now
  data_sample <- data_sample[data_sample$chrom == chrom_to_plot, ]
  segs_sample <- segs_sample[segs_sample$chr == chrom_to_plot, ]

  # Further split this plot into separate plots for logR, BAF and llratio(not necessary with mphase?)
  print(table(data_sample$sample_phased_by_this_sample))

  # logR track
  screen(subfigure_idx[[1]])
  plot_track_logr(data_sample, segs_sample)

  # BAF track
  screen(subfigure_idx[[2]])
  plot_track_baf(data_sample, segs_sample)

}

close.screen(all = TRUE)
dev.off()

