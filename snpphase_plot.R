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
  par(mar = c(1, 5, 0, 0))
  plot(data_sample$pos, data_sample$logR, ylab = "logR", xaxt = "n", col = adjustcolor("grey", 0.6), pch = 20)
  #lines(x = data_sample$pos, y = data_sample$logRsegmented, col=ifelse(is.na(data_sample$logRsegmented), "white", "black"))
  segments(x0 = segs_sample$startpos, y0 = segs_sample$median_non_hz_logr, x1 = segs_sample$endpos, col="black")
  abline(v = segs_sample$startpos, lty = 2, col = "black")
  text(x = segs_sample$seg_center, y = 0, labels = segs_sample$cnTotal)
}

plot_track_baf <- function(data_sample, segs_sample) {
  par(mar = c(1, 5, 0, 0))
  phasing_col <- ifelse(data_sample$maj_is_cpnB_logic_vec == 1, "blue", "yellow")
  #phasing_col <- ifelse(data_sample$, "blue", "yellow")
  phasing_col[data_sample$sample_phased_by_this_sample == "none_no_allelic_imbalance_called"] <- "grey"
  phasing_col <- adjustcolor(phasing_col, 0.6)
  plot(data_sample$pos, data_sample$BAF, ylab = "BAF", col = phasing_col, ylim = c(0, 1), xaxt = "n", pch = 20)
  #lines(x = data_sample$pos, y = data_sample$BAFsegmented, col="black")

  # Mean segments
  segments(x0 = segs_sample$startpos, y0 = segs_sample$mean_BAF + 0.5, x1 = segs_sample$endpos, col="black")
  segments(x0 = segs_sample$startpos, y0 = -1 * segs_sample$mean_BAF + 0.5, x1 = segs_sample$endpos, col="black")

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

# Read in all the data we need to plot
baf_logr <- readRDS(file.path(conf$input_dir, paste0(conf$patient_id, "_phased_BAF_df_list.RDS")))
segs <- readRDS(file.path(conf$input_dir, paste0(conf$patient_id, "_full_joint_segs_df.RDS")))
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
split.screen(c(length(samples) * num_subfigures, 1))

for (sample_idx in seq_along(samples)) {
  # Further split this plot into separate plots for logR, BAF and llratio(not necessary with mphase?)
  screen_idx <- 1 + (sample_idx - 1) * num_subfigures
  sample <- samples[[sample_idx]]
  data_sample <- baf_logr[[sample_idx]]
  segs_sample <- segs[segs$sample == sample,]
  # Remove the many NA regions
  segs_sample <- segs_sample[!is.na(segs_sample$cnTotal),]

  # Just plot one chromosome for now
  data_sample <- data_sample[data_sample$chrom == chrom_to_plot, ]
  segs_sample <- segs_sample[segs_sample$chr == chrom_to_plot, ]

  print(table(data_sample$sample_phased_by_this_sample))

  # logR track
  screen(screen_idx)
  plot_track_logr(data_sample, segs_sample)

  # BAF track
  screen(screen_idx + 1)
  plot_track_baf(data_sample, segs_sample)

}

close.screen(all = TRUE)
dev.off()

