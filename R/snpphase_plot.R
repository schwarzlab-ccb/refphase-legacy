#!/usr/bin/env Rscript

library(digest)
library(GenomicRanges)
library(optparse)

options(error = traceback)

option_list <- list(
  make_option(c("-s", "--original-segmentation"),
    default = "",
    help = "Pre-refphase segmentation (from sequenza)"
  ),
  make_option(c("-i", "--input-dir"),
    default = "",
    help = "Input directory containing the output of Xphase"
  ),
  make_option(c("--input-dir-old"),
    default = "",
    help = "Input directory containing the output of Mphase"
  ),
  make_option(c("-p", "--patient-id"),
    default = "",
    help = "Patient ID"
  ),
  make_option(c("-r", "--raw"),
    default = FALSE,
    help = "Plot raw CNs instead of integer CNs"
  ),
  make_option("--scale-segs",
    default = FALSE,
    help = "Rescale segment lengths to all be equal (for debugging)"
  ),
  make_option(c("-o", "--output-dir"),
    default = "",
    help = "Input directory containing the output of Mphase"
  ),
  make_option(c("-f", "--plot-format"),
    default = "png",
    help = "Input directory containing the output of Mphase"
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

# ------------------------------------------------------------------------
# Utility functions for plotting
# ------------------------------------------------------------------------
cn_text <- function(x, y, labels, cex = 1.5, col = "black") {
  text(x = x, y = y, labels = labels, cex = cex * 1.1, col = "white")
  text(x = x, y = y, labels = labels, cex = cex, col = col)
}

plot_track_axis <- function(xlim) {
  par(mar = c(2, 5, 0, 1))
  plot(y = c(0, 0), x = xlim, xlim = xlim, type = "n", ylab = "", xlab = "", axes = FALSE)
  axis(1, at = NULL)
}

# Plot CNs only in this track
plot_track_cn <- function(segs, track_name, xlim) {
  if (missing(xlim)) {
    xlim <- c(min(segs$startpos), max(segs$endpos))
  }

  par(mar = c(0.1, 5, 0.1, 1))
  num_segs <- nrow(segs)
  plot(1, ylim = c(-1, 1), xlim = xlim, col = "white", ylab = paste0(track_name, "\n(", num_segs, ")"), xaxt = "n", yaxt = "n")
  abline(v = unique(c(segs$startpos, segs$endpos)), lty = 2, col = "black")
  total_cn <- as.character(segs$cnTotal)
  major <- as.character(segs$nMajor)
  minor <- as.character(segs$nMinor)
  full_label <- paste(total_cn, major, minor, sep = "\n")

  col <- rep("grey", nrow(segs))
  if ("is_changed" %in% colnames(segs)) {
    col[segs$is_changed] <- "black"
  }
  cn_text(x = segs$seg_center, y = 0, labels = full_label, cex = 0.7, col = col)
}

plot_track_logr <- function(data_sample, segs_sample) {
  par(mar = c(0.1, 5, 0.1, 1))
  plot(data_sample$pos, data_sample$logR, ylab = "logR", xaxt = "n", col = adjustcolor("darkgrey", 0.2), pch = 20, ylim = c(-2, 2))
  segments(x0 = segs_sample$startpos, y0 = segs_sample$median_non_hz_logr, x1 = segs_sample$endpos, col = "black")
  abline(v = unique(c(segs_sample$startpos, segs_sample$endpos)), lty = 2, col = "black")

  # MPhase copy numbers
  total_cn <- as.character(segs_sample$orig_cnTotal)

  cn_text(x = segs_sample$seg_center, y = 1.5, labels = total_cn)
}

plot_track_baf <- function(data_sample, segs_sample, msai) {
  par(mar = c(0.1, 5, 0.1, 1))
  phasing_col <- ifelse(data_sample$maj_is_cpnB_logic_vec == 1, "blue", "orange")
  phasing_col[data_sample$sample_phased_by_this_sample == "none_no_allelic_imbalance_called"] <- "grey"
  phasing_col <- adjustcolor(phasing_col, 0.6)
  plot(data_sample$pos, data_sample$BAF, ylab = "BAF", col = phasing_col, ylim = c(0, 1), pch = 20, xaxt = "n", yaxt = "n")
  axis(2, seq(0, 1, by = 0.25))

  # Mean segments
  mean_col <- ifelse(segs_sample$orig_nMajor != segs_sample$orig_nMinor, "black", NA)
  segments(x0 = segs_sample$startpos, y0 = segs_sample$mean_BAF + 0.5, x1 = segs_sample$endpos, col = mean_col)
  segments(x0 = segs_sample$startpos, y0 = -1 * segs_sample$mean_BAF + 0.5, x1 = segs_sample$endpos, col = mean_col)

  # Segment boundaries
  abline(v = unique(c(segs_sample$startpos, segs_sample$endpos)), lty = 2, col = "black")

  # MPhase copy numbers
  major_cn <- as.character(segs_sample$orig_nMajor)
  minor_cn <- as.character(segs_sample$orig_nMinor)

  # major_cn <- as.character(segs_sample$nMajor)
  # minor_cn <- as.character(segs_sample$nMinor)

  cn_text(x = segs_sample$seg_center, y = 0.9, labels = major_cn)
  cn_text(x = segs_sample$seg_center, y = 0.1, labels = minor_cn)

  # Indicate which sample is the reference
  segs_sample$dominant_allele[is.na(segs_sample$dominant_allele)] <- "na"
  if (any(segs_sample$dominant_allele == "cpnB")) {
    segments(
      x0 = segs_sample$startpos[segs_sample$dominant_allele == "cpnB"], y0 = 1.02,
      x1 = segs_sample$endpos[segs_sample$dominant_allele == "cpnB"], col = "darkgreen", lwd = 4, lend = 1
    )
  }
  cn_text(x = segs_sample$seg_center, y = 0.5, labels = ifelse(segs_sample$dominant_allele == "cpnB", "msai", ""), col = "darkgreen")

  # Indicate which sample is the reference
  if (any(segs_sample$is_reference)) {
    segments(
      x0 = segs_sample$startpos[segs_sample$is_reference], y0 = -0.02,
      x1 = segs_sample$endpos[segs_sample$is_reference], col = "red", lwd = 4, lend = 1
    )
  }
}

# --------------------------------------------------------------------------
# Whole genome versions of the tracks

plot_track_logr_whole_genome <- function(data_sample, segs_sample, chrom_info) {
  par(mar = c(1, 5, 1, 0.1), xaxs = "i")
  plot(data_sample$pos, data_sample$logR, ylab = "logR", xaxt = "n", col = adjustcolor("darkgrey", 0.2), pch = ".", ylim = c(-2, 2))
  # lines(x = data_sample$pos, y = data_sample$logRsegmented, col=ifelse(is.na(data_sample$logRsegmented), "white", "black"))
  segments(x0 = segs_sample$startpos, y0 = segs_sample$median_non_hz_logr, x1 = segs_sample$endpos, col = "black")
  abline(v = chrom_info$start, lty = 2, col = "black")
  cn_text(x = chrom_info$center, y = 1.5, labels = chrom_info$label, cex = 1.0)
}

plot_track_baf_whole_genome <- function(data_sample, segs_sample, chrom_info) {
  par(mar = c(2, 5, 0.1, 0.1), xaxs = "i")
  phasing_col <- ifelse(data_sample$maj_is_cpnB_logic_vec == 1, "blue", "orange")
  phasing_col[data_sample$sample_phased_by_this_sample == "none_no_allelic_imbalance_called"] <- "grey"
  phasing_col <- adjustcolor(phasing_col, 0.6)
  plot(data_sample$pos, data_sample$BAF, ylab = "BAF", col = phasing_col, ylim = c(0, 1), pch = ".")

  # Mean segments
  mean_col <- ifelse(segs_sample$within_samp_cl_AI_present, "black", NA)
  segments(x0 = segs_sample$startpos, y0 = segs_sample$mean_BAF + 0.5, x1 = segs_sample$endpos, col = mean_col)
  segments(x0 = segs_sample$startpos, y0 = -1 * segs_sample$mean_BAF + 0.5, x1 = segs_sample$endpos, col = mean_col)

  # Chromosome boundaries
  abline(v = chrom_info$start, lty = 2, col = "black")
}

# --------------------------------------------------------------------------

plot_sample_info <- function(sample, ploidy, cellularity) {
  par(mar = c(0, 0, 0, 0))
  plot(c(0, 1), c(0, 1), ann = F, bty = "n", type = "n", xaxt = "n", yaxt = "n")

  text(
    x = 0.5, y = 0.5, adj = 0.5,
    labels = paste0(sample, "\n", "cellularity:", round(cellularity, 2), "  ploidy:", round(ploidy, 2)),
    cex = 1.3, col = "black", srt = 90
  )
}

plot_chromosome_header <- function(sample, chromosome, cex = 1.5) {
  par(mar = c(0, 0, 0, 0))
  plot(c(0, 1), c(0, 1), ann = F, bty = "n", type = "n", xaxt = "n", yaxt = "n")

  text(
    x = 0.5, y = 0.5, adj = 0.5, labels = paste0(sample, ": ", chromosome),
    cex = cex, col = "black"
  )
}

# ------------------------------------------------------------------------
# Calculate which segments have been changed by refphase
get_seg_changes <- function(new_seg, old_seg) {

  # cat("------ new:")
  # str(new_seg)
  # cat("------ old:")
  # str(old_seg)

  new_seg$seg_idx <- seq(nrow(new_seg))
  gr_new <- makeGRangesFromDataFrame(new_seg,
    keep.extra.columns = TRUE,
    ignore.strand = FALSE,
    seqnames.field = "chr",
    start.field = "startpos",
    end.field = "endpos"
  )
  stopifnot(length(gr_new) == nrow(new_seg))

  gr_orig <- makeGRangesFromDataFrame(old_seg,
    keep.extra.columns = TRUE,
    ignore.strand = FALSE,
    seqnames.field = "chr",
    start.field = "startpos",
    end.field = "endpos"
  )
  stopifnot(length(gr_orig) == nrow(old_seg))

  ov <- findOverlaps(gr_new, gr_orig)
  new_cns <- mcols(gr_new)[from(ov), c("cnTotal", "nMinor", "nMajor")]
  orig_cns <- mcols(gr_orig)[to(ov), c("cnTotal", "nMinor", "nMajor")]

  na_row <- apply(is.na(orig_cns), 1, any) | apply(is.na(new_cns), 1, any)

  orig_cns <- orig_cns[!na_row, ]
  new_cns <- new_cns[!na_row, ]

  # write.csv(cbind(ov[!na_row, ], new_cns, orig_cns), file="all_cns.csv")

  message("# of segments not changed=", nrow(new_cns))
  seg_changed <- apply(as.matrix(orig_cns) != as.matrix(new_cns), 1, any)
  # message("# of segments changed=", sum(seg_changed))
  message("# of segments not changed=", sum(!seg_changed))

  new_changed <- by(seg_changed, from(ov)[!na_row], any)
  message("# of segments changed=", sum(unlist(new_changed)))

  changed_idx <- mcols(gr_new)[from(ov)[!na_row], "seg_idx"][new_changed]
  # print(changed_idx)
  is_changed <- FALSE
  is_changed[new_seg$seg_idx %in% changed_idx] <- TRUE

  return(is_changed)
}

get_seg_new_ai <- function(new_seg, old_seg) {
  new_seg$seg_idx <- seq(nrow(new_seg))
  gr_new <- makeGRangesFromDataFrame(new_seg,
    keep.extra.columns = TRUE,
    ignore.strand = FALSE,
    seqnames.field = "chr",
    start.field = "startpos",
    end.field = "endpos"
  )
  stopifnot(length(gr_new) == nrow(new_seg))

  gr_orig <- makeGRangesFromDataFrame(old_seg,
    keep.extra.columns = TRUE,
    ignore.strand = FALSE,
    seqnames.field = "chr",
    start.field = "startpos",
    end.field = "endpos"
  )
  stopifnot(length(gr_orig) == nrow(old_seg))

  ov <- findOverlaps(gr_new, gr_orig)
  new_cns <- mcols(gr_new)[from(ov), c("cnTotal", "nMinor", "nMajor")]
  orig_cns <- mcols(gr_orig)[to(ov), c("cnTotal", "nMinor", "nMajor")]

  na_row <- apply(is.na(orig_cns), 1, any) | apply(is.na(new_cns), 1, any)

  orig_cns <- orig_cns[!na_row, ]
  new_cns <- new_cns[!na_row, ]

  # write.csv(cbind(ov[!na_row, ], new_cns, orig_cns), file="all_cns.csv")

  message("# of segments not changed=", nrow(new_cns))
  orig_cns <- as.matrix(orig_cns)
  new_cns <- as.matrix(new_cns)

  message("CN for checking new_ai")
  print(orig_cns)
  print(new_cns)

  orig_balanced <- orig_cns[, "nMinor"] == orig_cns[, "nMajor"]
  new_not_balanced <- new_cns[, "nMinor"] != new_cns[, "nMajor"]
  seg_changed <- orig_balanced & new_not_balanced
  print(seg_changed)
  print("--------------------------")
  # message("# of segments changed=", sum(seg_changed))
  message("# of segments not changed=", sum(!seg_changed))

  new_changed <- by(seg_changed, from(ov)[!na_row], any)
  message("# of segments changed=", sum(unlist(new_changed)))

  changed_idx <- mcols(gr_new)[from(ov)[!na_row], "seg_idx"][new_changed]
  # print(changed_idx)
  is_changed <- FALSE
  is_changed[new_seg$seg_idx %in% changed_idx] <- TRUE

  print(is_changed)
  message("Num changed:", sum(is_changed))
  print("--------------------------")

  return(is_changed)
}

if (!dir.exists(conf$output_dir)) {
  dir.create(conf$output_dir, recursive = TRUE, showWarnings = FALSE)
}

# Read in all the data we need to plot
orig_segs <- read.table(conf$original_segmentation, stringsAsFactors = FALSE, sep = "\t", header = TRUE)
orig_segs$seg_center <- orig_segs$startpos + (orig_segs$endpos - orig_segs$startpos) / 2
orig_segs$sample <- gsub("_LogR$", "", orig_segs$sample)
orig_segs$sample <- gsub(".", "-", orig_segs$sample, fixed = TRUE)

cat("----------- orig segs after reading\n")
# str(orig_segs)

# Out from mphase (the version before Tom's most recent changes that add ASCAT
# runs)
old_segs <- read.table(file.path(conf$input_dir_old, "full_joint_segs.tsv"), stringsAsFactors = FALSE, sep = "\t", header = TRUE)
old_segs$seg_center <- old_segs$startpos + (old_segs$endpos - old_segs$startpos) / 2
old_segs$dominant_allele[is.na(old_segs$dominant_allele)] <- "na"

baf_logr <- readRDS(file.path(conf$input_dir, paste0(conf$patient_id, "_phased_BAF_df_list.RDS")))
segs <- readRDS(file.path(conf$input_dir, paste0(conf$patient_id, "_full_joint_segs_df.RDS")))
segs$seg_center <- segs$startpos + (segs$endpos - segs$startpos) / 2
segs$is_reference <- segs$samp_used_to_phase == segs$sample
segs$dominant_allele[is.na(segs$dominant_allele)] <- "na"

print(unique(sort(orig_segs$sample)))
print(unique(sort(segs$sample)))
print(unique(sort(old_segs$sample)))
stopifnot(all(old_segs$sample %in% orig_segs$sample))
stopifnot(all(segs$sample %in% orig_segs$sample))

# Size distribution of segments
pdf(file.path(conf$output_dir, "segment-size-distribution-orig.pdf"))
hist(orig_segs$endpos - orig_segs$startpos, xlab = "Segment width")
dev.off()

pdf(file.path(conf$output_dir, "segment-size-distribution-oldrefphase.pdf"))
hist(old_segs$endpos - old_segs$startpos, xlab = "Segment width")
dev.off()

pdf(file.path(conf$output_dir, "segment-size-distribution-newrefphase.pdf"))
hist(segs$endpos - segs$startpos, xlab = "Segment width")
dev.off()

keep_orig <- !is.na(orig_segs$cnTotal)
keep_old <- !is.na(old_segs$nMinor)
keep_new <- !is.na(segs$orig_nMinor)

# Dump some basic stats
stats_file <- file.path(conf$output_dir, "misc-stats.tsv")
stats <- data.frame(
  cn_diff_total = c(NA, 0, 0),
  cn_diff_changed = c(NA, 0, 0),
  cn_diff_new_ai = c(NA, 0, 0),
  num_segments = c(
    sum(!is.na(orig_segs$cnTotal)),
    sum(!is.na(old_segs$nMinor)),
    sum(!is.na(segs$orig_nMinor))
  ),
  num_segments_with_na = c(
    nrow(orig_segs),
    nrow(old_segs),
    nrow(segs)
  ),
  num_segments_under_1kb = c(
    sum(orig_segs[keep_orig, ]$endpos - orig_segs[keep_orig, ]$startpos < 1000),
    sum(old_segs[keep_old, ]$endpos - old_segs[keep_old, ]$startpos < 1000),
    sum(segs[keep_new, ]$endpos - segs[keep_new, ]$startpos < 1000)
  ),
  num_homozygous_deletion = c(
    sum(orig_segs[keep_orig, ]$cnTotal == 0, na.rm = TRUE),
    sum(old_segs[keep_old, ]$cnTotal == 0, na.rm = TRUE),
    sum(segs[keep_new, ]$orig_cnTotal == 0, na.rm = TRUE)
  ),
  row.names = c("original", "oldrefphase", "newrefphase")
)

stats$fraction_segments_under_1kb <- stats$num_segments_under_1kb / stats$num_segments
stats$fraction_homozygous_deletion <- stats$num_homozygous_deletion / stats$num_segments

# Convert chromosomal positions to absolute positions for plotting the whole genome
chrom_ends <- cumsum(by(segs$endpos, segs$chr, max))
chrom_starts <- c(0, chrom_ends[-length(chrom_ends)])
names(chrom_starts) <- names(chrom_ends)
chrom_center <- chrom_starts + (chrom_ends - chrom_starts) / 2
segs$chrom_startpos <- chrom_starts[segs$chr]

# Chromosome start and end positions in genome-wide absolute positions
chrom_info <- data.frame(
  start = chrom_starts,
  end = chrom_ends,
  center = chrom_center,
  label = names(chrom_starts)
)

num_subfigures <- 6 # BAF and logR subfigures only + copy numbers

# Split the output figure into one row per sample
samples <- names(baf_logr)

output_plot <- file.path(conf$output_dir, "all_chromosomes.pdf")
pdf(output_plot, width = 10, height = 1 + 4 * length(samples))
for (chrom_to_plot in names(chrom_starts)) {
  par(cex = 0.8)
  # Header: Leave space for the chromosome name
  # coordinates are: left, right, bottom, and top
  # heading_height_fraction <- 1 / (length(samples) + 0.5)
  heading_height_fraction <- 1 / (1 + 4 * length(samples))
  # split_coords <- matrix(c(0, 1, 0, 1 - heading_height_fraction, 0, 1, 0, heading_height_fraction), ncol = 4, byrow = TRUE)
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
    message("Sample: ", sample_name, " Chromosome: ", chrom_to_plot)
    data_sample <- baf_logr[[sample_idx]]
    segs_sample <- segs[segs$sample == sample_name, ]
    orig_segs_sample <- orig_segs[orig_segs$sample == sample_name, ]
    old_segs_sample <- old_segs[old_segs$sample == sample_name, ]

    # message("-------------------- segs_sample")
    # str(segs_sample)
    # message("-------------------- orig_segs_sample")
    # str(orig_segs_sample)
    # message("-------------------- old_segs_sample")
    # str(old_segs_sample)

    # Split the plotting area vertically (sample name + info on left, plots on right)
    screen(sample_screens[[sample_idx]])
    split_coords <- matrix(c(0, 0.1, 0, 1, 0.1, 1, 0, 1), ncol = 4, byrow = TRUE)
    screen_idxs <- split.screen(split_coords)

    screen(screen_idxs[[1]])
    plot_sample_info(sample_name, segs_sample$Ploidy[[1]], segs_sample$ACF[[1]])

    screen(screen_idxs[[2]])
    sub_split_coords <- matrix(c(
      0, 1, 0.7, 1,
      0, 1, 0.4, 0.7,
      0, 1, 0.3, 0.4,
      0, 1, 0.2, 0.3,
      0, 1, 0.1, 0.2,
      0, 1, 0, 0.1
    ), ncol = 4, byrow = TRUE)
    subfigure_idx <- split.screen(sub_split_coords)
    # subfigure_idx <- split.screen(c(num_subfigures, 1))

    # Remove the many NA regions
    segs_sample <- segs_sample[!is.na(segs_sample$nMajor), ]
    old_segs_sample <- old_segs_sample[!is.na(old_segs_sample$nMajor), ]

    # Keep data for this sample
    data_sample <- data_sample[data_sample$chrom == chrom_to_plot, ]
    segs_sample <- segs_sample[segs_sample$chr == chrom_to_plot, ]
    orig_segs_sample <- orig_segs_sample[orig_segs_sample$chr == chrom_to_plot, ]
    old_segs_sample <- old_segs_sample[old_segs_sample$chr == chrom_to_plot, ]

    # message("-------------------- segs_sample")
    # str(segs_sample)
    # message("-------------------- orig_segs_sample")
    # str(orig_segs_sample)
    # message("-------------------- old_segs_sample")
    # str(old_segs_sample)


    # Adjust segment widths
    if (conf$scale_segs) {
      # segs_sample <- segs[segs$sample == sample_name, ]
      # orig_segs_sample <- orig_segs[orig_segs$sample == sample_name, ]
      # old_segs_sample <- old_segs[old_segs$sample == sample_name, ]

      data_chrom <- data[data$chrom == chrom_to_plot, ]
      segs_chrom <- segs[segs$chr == chrom_to_plot, ]
      orig_segs_chrom <- orig_segs[orig_segs$chr == chrom_to_plot, ]
      old_segs_chrom <- old_segs[old_segs$chr == chrom_to_plot, ]

      allbreaks <- sort(unique(c(
        segs_chrom$startpos,
        segs_chrom$endpos,
        old_segs_chrom$startpos,
        old_segs_chrom$endpos,
        orig_segs_chrom$startpos,
        orig_segs_chrom$endpos
      )))

      # allbreaks <- sort(unique(c(segs_sample$startpos,
      #                           segs_sample$endpos,
      #                           old_segs_sample$startpos,
      #                           old_segs_sample$endpos,
      #                           orig_segs_sample$startpos,
      #                           orig_segs_sample$endpos)))

      all_segs <- list()
      all_segs$startpos <- allbreaks[-length(allbreaks)]
      all_segs$endpos <- allbreaks[-1]
      all_segs$width <- all_segs$endpos - all_segs$startpos
      all_segs$scaling_factor <- 1 / all_segs$width
      all_segs$startidx <- seq_len(length(all_segs$startpos)) - 1
      all_segs$endidx <- all_segs$startidx + 1

      plot_range <- c(min(all_segs$startidx), max(all_segs$endidx))

      # message("----------------- all_segs")
      # print(all_segs)
      # Correct logr and baf data
      for (segidx in seq_len(length(all_segs$startpos))) {
        # which_seg <- segs_sample$startpos[[segidx]] <= data_sample$pos & segs_sample$endpos[[segidx]] >= data_sample$pos
        which_pos <- data_sample$pos >= all_segs$startpos[[segidx]] & data_sample$pos <= all_segs$endpos[[segidx]]
        if (sum(which_pos) == 0) {
          next
        }
        data_sample$pos[which_pos] <- all_segs$startidx[[segidx]] + ((data_sample$pos[which_pos] - all_segs$startpos[[segidx]]) * all_segs$scaling_factor[[segidx]])
      }

      for (segidx in seq_len(length(all_segs$startpos))) {
        which_pos <- segs_sample$startpos == all_segs$startpos[[segidx]]
        segs_sample$startpos[which_pos] <- all_segs$startidx[[segidx]]

        which_pos <- segs_sample$endpos == all_segs$endpos[[segidx]]
        segs_sample$endpos[which_pos] <- all_segs$endidx[[segidx]]

        segs_sample$seg_center <- segs_sample$startpos + ((segs_sample$endpos - segs_sample$startpos) / 2)
      }

      for (segidx in seq_len(length(all_segs$startpos))) {
        which_pos <- orig_segs_sample$startpos == all_segs$startpos[[segidx]]
        orig_segs_sample$startpos[which_pos] <- all_segs$startidx[[segidx]]

        which_pos <- orig_segs_sample$endpos == all_segs$endpos[[segidx]]
        orig_segs_sample$endpos[which_pos] <- all_segs$endidx[[segidx]]

        orig_segs_sample$seg_center <- orig_segs_sample$startpos + ((orig_segs_sample$endpos - orig_segs_sample$startpos) / 2)
      }

      for (segidx in seq_len(length(all_segs$startpos))) {
        which_pos <- old_segs_sample$startpos == all_segs$startpos[[segidx]]
        old_segs_sample$startpos[which_pos] <- all_segs$startidx[[segidx]]

        which_pos <- old_segs_sample$endpos == all_segs$endpos[[segidx]]
        old_segs_sample$endpos[which_pos] <- all_segs$endidx[[segidx]]

        old_segs_sample$seg_center <- old_segs_sample$startpos + ((old_segs_sample$endpos - old_segs_sample$startpos) / 2)
      }
    } else {
      plot_range <- c(min(data_sample$pos), max(data_sample$pos))
    }

    message("Plot range:", plot_range[[1]], "-", plot_range[[2]])

    # ------------------------------------------------------------------------
    # Swap CNs when we have MSAI
    if (conf$raw) {
      segs_sample$nMajor <- round(ifelse(segs_sample$dominant_allele == "cpnB", segs_sample$orig_cpnB_raw, segs_sample$orig_cpnA_raw), 2)
      segs_sample$nMinor <- round(ifelse(segs_sample$dominant_allele == "cpnB", segs_sample$orig_cpnA_raw, segs_sample$orig_cpnB_raw), 2)
      segs_sample$cnTotal <- round(segs_sample$nMajor + segs_sample$nMinor, 2)

      old_segs_sample$nMajor <- round(old_segs_sample$ph_raw_cpnA, 2)
      old_segs_sample$nMinor <- round(old_segs_sample$ph_raw_cpnB, 2)
      old_segs_sample$cnTotal <- round(old_segs_sample$nMajor + old_segs_sample$nMinor, 2)
    } else {
      segs_sample$nMajor <- ifelse(segs_sample$dominant_allele == "cpnB", segs_sample$orig_nMinor, segs_sample$orig_nMajor)
      segs_sample$nMinor <- ifelse(segs_sample$dominant_allele == "cpnB", segs_sample$orig_nMajor, segs_sample$orig_nMinor)
      segs_sample$cnTotal <- segs_sample$orig_cnTotal

      nMajor <- ifelse(old_segs_sample$dominant_allele == "cpnB", old_segs_sample$nMinor, old_segs_sample$nMajor)
      nMinor <- ifelse(old_segs_sample$dominant_allele == "cpnB", old_segs_sample$nMajor, old_segs_sample$nMinor)
      old_segs_sample$nMajor <- nMajor
      old_segs_sample$nMinor <- nMinor
    }

    # ------------------------------------------------------------------------
    # Indicate which segments have changed
    segs_sample$is_changed <- get_seg_changes(segs_sample, orig_segs_sample)
    old_segs_sample$is_changed <- get_seg_changes(old_segs_sample, orig_segs_sample)

    # Count the number of segments that are different to the original
    stats["newrefphase", "cn_diff_total"] <- stats["newrefphase", "cn_diff_total"] + length(segs_sample$is_changed)
    stats["oldrefphase", "cn_diff_total"] <- stats["oldrefphase", "cn_diff_total"] + length(old_segs_sample$is_changed)

    stats["newrefphase", "cn_diff_changed"] <- stats["newrefphase", "cn_diff_changed"] + sum(segs_sample$is_changed, na.rm = TRUE)
    stats["oldrefphase", "cn_diff_changed"] <- stats["oldrefphase", "cn_diff_changed"] + sum(old_segs_sample$is_changed, na.rm = TRUE)

    segs_sample$is_new_ai <- get_seg_new_ai(segs_sample, orig_segs_sample)
    old_segs_sample$is_new_ai <- get_seg_new_ai(old_segs_sample, orig_segs_sample)

    stats["newrefphase", "cn_diff_new_ai"] <- stats["newrefphase", "cn_diff_new_ai"] + sum(segs_sample$is_new_ai, na.rm = TRUE)
    stats["oldrefphase", "cn_diff_new_ai"] <- stats["oldrefphase", "cn_diff_new_ai"] + sum(old_segs_sample$is_new_ai, na.rm = TRUE)

    # logR track
    screen(subfigure_idx[[1]])
    plot_track_logr(data_sample, segs_sample)

    # BAF track
    screen(subfigure_idx[[2]])
    plot_track_baf(data_sample, segs_sample)

    # Orig CNs
    screen(subfigure_idx[[3]])
    plot_track_cn(orig_segs_sample, "orig", xlim = plot_range)

    # Refphase CNs
    screen(subfigure_idx[[4]])
    plot_track_cn(segs_sample, track_name = "new", xlim = plot_range)

    # pre-ASCAT mphase CNs
    screen(subfigure_idx[[5]])
    plot_track_cn(old_segs_sample, "noascat", xlim = plot_range)

    # Add an x-axis with genomic co-ordinates
    screen(subfigure_idx[[6]])
    plot_track_axis(plot_range)
  }

  close.screen(all = TRUE)
}
dev.off()

# Finally write the stats
stats$patient <- conf$patient_id
write.table(stats, file = stats_file, sep = "\t", col.names = NA)

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
