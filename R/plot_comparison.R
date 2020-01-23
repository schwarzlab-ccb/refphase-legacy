##!/usr/bin/env Rscript
#
#library("digest")
#library("optparse")
#
#options(error = traceback)
#
#option_list <- list(
#  make_option("--snps",
#    default = "",
#    help = "SNP data from refphase (rds)"
#  ),
#  make_option(c("-s", "--segmentations"),
#    default = "",
#    help = "Comma-separated segmentations in refphase segmentations.tsv format"
#  ),
#  make_option(c("-n", "--seg-names"),
#    default = "",
#    help = "Comma-separated list of segmentation names that match the files in --segmentations"
#  ),
#  make_option("--sample-info",
#    default = "",
#    help = "Per-sample information including purity and ploidy"
#  ),
#  make_option(c("-p", "--patient-id"),
#    default = "",
#    help = "Patient ID"
#  ),
#  make_option(c("-r", "--raw"),
#    default = FALSE,
#    help = "Plot raw CNs instead of integer CNs"
#  ),
#  make_option("--scale-segs",
#    default = FALSE,
#    help = "Rescale segment lengths to all be equal (for debugging)"
#  ),
#  make_option(c("-o", "--output-dir"),
#    default = "",
#    help = "Input directory containing the output of Mphase"
#  ),
#  make_option(c("-d", "--seed"),
#    default = 531,
#    help = "Set the random seed at the beginning of the script"
#  )
#)
#args <- parse_args(OptionParser(option_list = option_list), convert_hyphens_to_underscores = TRUE)
#conf <- args
#conf$id <- substr(digest(args), 1, 8)
#set.seed(args$seed)
#
#library("GenomicRanges")
#library("gtools")

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
  # We don't have much space for track names so let's just display the first
  # few characters
  short_track_name <- substring(track_name, 1, 5)
  plot(1, ylim = c(-1, 1), xlim = xlim, col = "white", ylab = paste0(short_track_name, "\n(", num_segs, ")"), xaxt = "n", yaxt = "n")

  # Color the background of each segment to make it clear where there are
  # missing segments
  polygon(x = unlist(mapply(c, segs$startpos, segs$endpos, segs$endpos, segs$startpos, NA, SIMPLIFY = FALSE)), y = rep(c(rep(c(-1, 1), each = 2), NA), nrow(segs)), col = "grey95", border = NA)

  abline(v = unique(c(segs$startpos, segs$endpos)), lty = 2, col = "black")
  total_cn <- as.character(segs$cpn_total)
  major <- as.character(segs$cpn_a)
  minor <- as.character(segs$cpn_b)
  full_label <- paste(total_cn, major, minor, sep = "\n")

  cn_text(x = segs$seg_center, y = 0, labels = full_label, cex = 0.7, col = "black")

  # Indicate if the CN has changed
  if (any(segs$is_changed)) {
    segments(
      x0 = segs$startpos[segs$is_changed], y0 = -0.98,
      x1 = segs$endpos[segs$is_changed], col = "yellow", lwd = 4, lend = 1
    )
  }
}

plot_track_logr <- function(data_sample, segs_sample, max_values = 6000) {
  par(mar = c(0.1, 5, 0.1, 1))
  # Downsample logR values to keep file sizes reasonable
  data_sample <- data_sample[sample.int(nrow(data_sample), min(max_values, nrow(data_sample))),]
  plot(data_sample$pos, data_sample$logR, ylab = "logR", xaxt = "n", col = adjustcolor("darkgrey", 0.2), pch = 16, cex = 0.75, ylim = c(min(c(-2, data_sample$logR)), max(c(2, data_sample$logR))))
  abline(v = unique(c(segs_sample$startpos, segs_sample$endpos)), lty = 2, col = "black")

  # Fit a loess regression line
  fit <- loess(data_sample$logR ~ data_sample$pos, span = 0.05)
  points(x = data_sample$pos, y = predict(fit, data_sample$pos), col = "black", pch = ".")
}

plot_track_baf <- function(data_sample, segs_sample, msai) {
  par(mar = c(0.1, 5, 0.1, 1))

  phasing_col <- data_sample$phasing
  phasing_col[phasing_col == "a"] <- "orange"
  phasing_col[phasing_col == "b"] <- "blue"
  phasing_col[phasing_col == "none"] <- "grey"
  phasing_col <- adjustcolor(phasing_col, 0.6)
  plot(data_sample$pos, data_sample$BAF, ylab = "BAF", col = phasing_col, ylim = c(0, 1), pch = 16, cex = 0.75, yaxt = "n", xaxt = "n")
  axis(2, seq(0, 1, by = 0.25))

#  # Mean segments
  #phasing_col <- ifelse(data_sample$maj_is_cpnB_logic_vec == 1, "blue", "orange")
  #phasing_col[data_sample$sample_phased_by_this_sample == "none_no_allelic_imbalance_called"] <- "grey"
  #phasing_col <- adjustcolor(phasing_col, 0.6)
#  mean_col <- ifelse(segs_sample$orig_nMajor != segs_sample$orig_nMinor, "black", NA)
#  segments(x0 = segs_sample$startpos, y0 = segs_sample$mean_BAF + 0.5, x1 = segs_sample$endpos, col = mean_col)
#  segments(x0 = segs_sample$startpos, y0 = -1 * segs_sample$mean_BAF + 0.5, x1 = segs_sample$endpos, col = mean_col)
#
  # Segment boundaries
  abline(v = unique(c(segs_sample$startpos, segs_sample$endpos)), lty = 2, col = "black")

  # Indicate which sample is the reference
  if (any(segs_sample$is_msai)) {
    segments(
      x0 = segs_sample$startpos[segs_sample$is_msai], y0 = 1.02,
      x1 = segs_sample$endpos[segs_sample$is_msai], col = "darkgreen", lwd = 4, lend = 1
    )
  }
  cn_text(x = segs_sample$seg_center, y = 0.5, labels = ifelse(segs_sample$is_msai, "msai", ""), col = "darkgreen")

  # Indicate which sample is the reference
  if (any(segs_sample$is_reference)) {
    segments(
      x0 = segs_sample$startpos[segs_sample$is_reference], y0 = -0.02,
      x1 = segs_sample$endpos[segs_sample$is_reference], col = "red", lwd = 4, lend = 1
    )
  }
}

# --------------------------------------------------------------------------

plot_sample_info <- function(sample, ploidy, cellularity) {
  par(mar = c(0, 0, 0, 0))
  plot(c(0, 1), c(0, 1), ann = FALSE, bty = "n", type = "n", xaxt = "n", yaxt = "n")

  text(
    x = 0.5, y = 0.5, adj = 0.5,
    labels = paste0(sample, "\n", "cellularity:", round(cellularity, 2), "  ploidy:", round(ploidy, 2)),
    cex = 1.3, col = "black", srt = 90
  )
}

plot_chromosome_header <- function(sample, chromosome, cex = 1.5) {
  par(mar = c(0, 0, 0, 0))
  plot(c(0, 1), c(0, 1), ann = FALSE, bty = "n", type = "n", xaxt = "n", yaxt = "n")

  text(
    x = 0.5, y = 0.5, adj = 0.5, labels = paste0(sample, ": ", chromosome),
    cex = cex, col = "black"
  )
}

# ------------------------------------------------------------------------
# Calculate which segments have been changed by refphase
#
# Returns a logical vector the length of new_seg, indicating which segments in
# the new segmentation do not have the same copy number as the corresponding
# segment region in the old segmentation
get_seg_changes <- function(new_seg, old_seg) {

  print("new_seg")
  str(new_seg)
  print("old_seg")
  str(old_seg)

  new_seg$seg_idx <- seq_len(nrow(new_seg))
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
  new_cns <- mcols(gr_new)[from(ov), c("cpn_total", "cpn_b", "cpn_a")]
  orig_cns <- mcols(gr_orig)[to(ov), c("cpn_total", "cpn_b", "cpn_a")]

  print(cbind(new_cns, orig_cns))

  na_row <- apply(is.na(orig_cns), 1, any) | apply(is.na(new_cns), 1, any)

  orig_cns <- orig_cns[!na_row, ]
  new_cns <- new_cns[!na_row, ]

  # Default to saying that the region is changed, and then identify the regions
  # that have NOT been changed. Should deal with regions that have no
  # overlapping segment in the original.
  is_changed <- rep(TRUE, nrow(new_seg))

  seg_not_changed <- apply(as.matrix(orig_cns) == as.matrix(new_cns), 1, all)
  new_not_changed <- by(seg_not_changed, from(ov)[!na_row], all)

  is_changed[] <- !unlist(new_not_changed)

  print(cbind(new_seg, is_changed=is_changed))

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
  new_cns <- mcols(gr_new)[from(ov), c("cpn_total", "cpn_b", "cpn_a")]
  orig_cns <- mcols(gr_orig)[to(ov), c("cpn_total", "cpn_b", "cpn_a")]

  na_row <- apply(is.na(orig_cns), 1, any) | apply(is.na(new_cns), 1, any)

  orig_cns <- orig_cns[!na_row, ]
  new_cns <- new_cns[!na_row, ]

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
  message("# of segments not changed=", sum(!seg_changed))

  new_changed <- by(seg_changed, from(ov)[!na_row], any)
  message("# of segments changed=", sum(unlist(new_changed)))

  changed_idx <- mcols(gr_new)[from(ov)[!na_row], "seg_idx"][new_changed]
  is_changed <- FALSE
  is_changed[new_seg$seg_idx %in% changed_idx] <- TRUE

  print(is_changed)
  message("Num changed:", sum(is_changed))
  print("--------------------------")

  return(is_changed)
}

# Do reference phasing according to the reference samples that are indicated in
# the input segmentation seg
# baf_threshold could possibly be changed to account for reference bias
phase_regions <- function(snps, seg, baf_threshold = 0.5) {

  for (samp in names(snps))
    snps[[samp]]$phasing <- "none"

  phased_regions <- seg[seg$is_reference, c("sample", "chr", "startpos", "endpos")]
  for (region_idx in seq_len(nrow(phased_regions))) {
    region <- phased_regions[region_idx, ]

    snp_idx_region <- which(snps[[region$sample]]$chr == region$chr & snps[[region$sample]]$pos >= region$startpos & snps[[region$sample]]$pos <= region$endpos)
    snps[[region$sample]][snp_idx_region, "phasing"] <- ifelse(snps[[region$sample]][snp_idx_region, "BAF"] > baf_threshold, "a", "b")

    # Get the phasing from our reference
    for (samp in names(snps)) {
      if (samp == region$sample)
        next

      snps[[samp]][snp_idx_region, "phasing"] <- snps[[region$sample]][snp_idx_region, "phasing"]
    }
  }

  snps
}

# --------------------------------------------------------------------------
# Whole genome versions of the tracks

plot_track_logr_whole_genome <- function(data_sample, segs_sample, chrom_info) {

  par(mar = c(1, 5, 1, 0.1), xaxs = "i")
  plot(data_sample$pos, data_sample$logR, ylab = "logR", xaxt = "n", col = adjustcolor("darkgrey", 0.2), pch = ".", ylim = c(-2, 2))

  fit <- loess(data_sample$logR ~ data_sample$pos, span = 0.001)
  points(x = data_sample$pos, y = predict(fit, data_sample$pos), col = "black", pch = ".")
  # We don't have the meana right now
  #segments(x0 = segs_sample$startpos, y0 = segs_sample$median_non_hz_logr, x1 = segs_sample$endpos, col = "black")
  abline(v = chrom_info$start, lty = 2, col = "black")
  cn_text(x = chrom_info$center, y = 1.5, labels = chrom_info$label, cex = 1.0)
}

plot_track_baf_whole_genome <- function(data_sample, segs_sample, chrom_info) {
  par(mar = c(2, 5, 0.1, 0.1), xaxs = "i")

  phasing_col <- data_sample$phasing
  phasing_col[phasing_col == "a"] <- "orange"
  phasing_col[phasing_col == "b"] <- "blue"
  phasing_col[phasing_col == "none"] <- "grey"
  phasing_col <- adjustcolor(phasing_col, 0.6)
  plot(data_sample$pos, data_sample$BAF, ylab = "BAF", col = phasing_col, ylim = c(0, 1), pch = ".", xaxt = "n")

#  # Mean segments
#  mean_col <- ifelse(segs_sample$within_samp_cl_AI_present, "black", NA)
#  segments(x0 = segs_sample$startpos, y0 = segs_sample$mean_BAF + 0.5, x1 = segs_sample$endpos, col = mean_col)
#  segments(x0 = segs_sample$startpos, y0 = -1 * segs_sample$mean_BAF + 0.5, x1 = segs_sample$endpos, col = mean_col)

  # Chromosome boundaries
  abline(v = chrom_info$start, lty = 2, col = "black")
}

# Plot CNs only in this track
plot_track_cn_whole_genome <- function(segs, chrom_info, xlim = c(NA, NA), max_cn = 10) {
  par(mar = c(2, 5, 0.1, 0.1), xaxs = "i")
  num_segs <- nrow(segs)

  # Set a max CN for the sake of plotting
  segs$cpn_a <- pmin(segs$cpn_a, max_cn)
  segs$cpn_b <- pmin(segs$cpn_b, max_cn)

  # We don't have much space for track names so let's just display the first
  # few characters
  plot(1, ylim = c(min(c(0, segs$cpn_a, segs$cpn_b)),
                   max(c(segs$cpn_a, segs$cpn_b))),
          xlim = c(min(c(segs$startpos, xlim[[1]]), na.rm = TRUE),
                   max(c(segs$endpos, xlim[[2]]), na.rm = TRUE)),
          ylab = "CN", type = "n")

  for (i in seq_len(nrow(segs))) {
    lines(x = c(segs$startpos[[i]], segs$endpos[[i]]), y = rep(segs$cpn_a[[i]], 2), col = adjustcolor("green", 0.4), lwd = 2, lend = 1)
    lines(x = c(segs$startpos[[i]], segs$endpos[[i]]), y = rep(segs$cpn_b[[i]], 2), col = adjustcolor("blue", 0.4), lwd = 2, lend = 1)
  }

  abline(v = chrom_info$start, lty = 2, col = "black")
}

# ------------------------------------------------------------------------
# Start the main script
# ------------------------------------------------------------------------

#' @export
plot_comparison <- function(conf) {

  if (!dir.exists(conf$output_dir)) {
    dir.create(conf$output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  # The BAF/LogR data along the genome
  snps <- readRDS(conf$snps)$all.samp.SNP.list

  #str(snps)

  # A list of copy number segmentations to plot as individual tracks
  seg_input_files <- strsplit(conf$segmentations, ",")[[1]]
  seg_names <- strsplit(conf$seg_names, ",")[[1]]

  stopifnot(length(seg_input_files) == length(seg_names))

  segs <- list()
  for (seg_idx in seq_along(seg_input_files)) {
    seg <- read.delim(seg_input_files[[seg_idx]], sep = "\t", stringsAsFactors = FALSE)

    # First we deal with the seg format that is used for refphase
    if ("nMajor" %in% colnames(seg)) {
      colnames(seg)[colnames(seg) == "nMajor"] <- "cpn_a"
      colnames(seg)[colnames(seg) == "nMinor"] <- "cpn_b"
      seg <- cbind(seg, cpn_total = seg$cpn_a + seg$cpn_b)
    }
    if ("cnTotal" %in% colnames(seg)) {
      colnames(seg)[colnames(seg) == "cnTotal"] <- "cpn_total"
    }
    seg$seg_center <- seg$startpos + ((seg$endpos - seg$startpos) / 2)
    # Make sure we have the data that we need to plot
    if (seg_idx == 1) {
      # We need some extra information from the first set of segmentatsions (msai + reference)
      stopifnot(c("sample", "chr", "startpos", "endpos", "cpn_a", "cpn_b", "is_msai", "is_reference") %in% colnames(seg))
    } else {
      stopifnot(c("sample", "chr", "startpos", "endpos", "cpn_a", "cpn_b") %in% colnames(seg))
    }

    segs[[seg_idx]] <- seg
  }

  sample_info <- read.delim(conf$sample_inf, sep = "\t", stringsAsFactors = FALSE, row.names = "sample")

  # We consider the first segmentation that is passed in as the "main"
  # segmentation, the one that matches the BAF/logR values and the one that is
  # marked on the BAF/logR plots (reference + msai events)
  main_seg <- segs[[1]]
  #chroms <- mixedsort(unique(sapply(segs, "[[", "chr")))

  # Assign phasing of regions with MSAI. I should directly carry over the phasing
  # from refphase instead of recalculating this here.
  snps <- phase_regions(snps, main_seg)

  chroms <- 1:22
  sample_names <- sort(names(snps))

  output_plot <- file.path(conf$output_dir, "all_chromosomes.pdf")
  pdf(output_plot, width = 10, height = 1 + 4 * length(snps))

  for (chrom_to_plot in chroms) {
    par(cex = 0.8)
    # Header: Leave space for the chromosome name
    # coordinates are: left, right, bottom, and top
    # heading_height_fraction <- 1 / (length(samples) + 0.5)
    heading_height_fraction <- 1 / (1 + 4 * length(sample_names))
    # split_coords <- matrix(c(0, 1, 0, 1 - heading_height_fraction, 0, 1, 0, heading_height_fraction), ncol = 4, byrow = TRUE)
    split_coords <- matrix(c(0, 1, 1 - heading_height_fraction, 1, 0, 1, 0, 1 - heading_height_fraction), ncol = 4, byrow = TRUE)
    outer_screens <- split.screen(split_coords)
    screen(outer_screens[[1]])
    # Plot header
    plot_chromosome_header(conf$patient_id, paste0("chr", chrom_to_plot))
    screen(outer_screens[[2]])

    # Split horizontally into figures for each sample
    sample_screens <- split.screen(c(length(sample_names), 1))
    for (sample_idx in seq_along(sample_names)) {
      sample_name <- sample_names[[sample_idx]]
      message("Sample: ", sample_name, " Chromosome: ", chrom_to_plot)
      data_sample <- snps[[sample_name]]
      data_sample <- data_sample[data_sample$chr == chrom_to_plot, ]

      # Split the plotting area vertically (sample name + info on left, plots on right)
      screen(sample_screens[[sample_idx]])
      split_coords <- matrix(c(0, 0.1, 0, 1, 0.1, 1, 0, 1), ncol = 4, byrow = TRUE)
      screen_idxs <- split.screen(split_coords)

      screen(screen_idxs[[1]])
      plot_sample_info(sample_name, sample_info[sample_name, "ploidy"], sample_info[sample_name, "purity"])

      screen(screen_idxs[[2]])

      # Heights are in arbitrary units, and will be used to calculate fractional
      # heights of each track. Only the relative difference between heights
      # matters.
      track_heights <- c(30,  # logR
                         30,  # baf
                         rep(10, length(segs)),  # segmentations
                         10)  # axis

      track_height_fractions <- rev(cumsum(rev(track_heights) / sum(track_heights)))

      sub_split_coords <- matrix(c(
        rep(0, length(track_height_fractions)),
        rep(1, length(track_height_fractions)),
        c(track_height_fractions[-1], 0),
        track_height_fractions
      ), ncol = 4, byrow = FALSE)

      subfigure_idx <- split.screen(sub_split_coords)

      plot_range <- c(min(data_sample$pos), max(data_sample$pos))

      main_seg_sample <- main_seg[main_seg$sample == sample_name & main_seg$chr == chrom_to_plot,]

      # logR track
      screen(subfigure_idx[[1]])
      plot_track_logr(data_sample, main_seg_sample)

      # BAF track
      screen(subfigure_idx[[2]])
      plot_track_baf(data_sample, main_seg_sample)

      seg_orig <- segs[[2]]
      for (seg_idx in seq_along(segs)) {
        seg <- segs[[seg_idx]]
        seg_chrom <- seg[seg$sample == sample_name & seg$chr == chrom_to_plot, ]
        seg_orig_chrom <- seg_orig[seg_orig$sample == sample_name & seg_orig$chr == chrom_to_plot, ]

        seg_chrom$is_changed <- rep(FALSE, nrow(seg_chrom))
        if (seg_idx == 1) {
          seg_chrom$is_changed <- get_seg_changes(seg_chrom, seg_orig_chrom)
        }

        # Orig CNs
        screen(subfigure_idx[[2 + seg_idx]])
        plot_track_cn(seg_chrom, seg_names[[seg_idx]], xlim = plot_range)
      }

      # Add an x-axis with genomic co-ordinates
      screen(subfigure_idx[[3 + length(segs)]])
      plot_track_axis(plot_range)
    }

    close.screen(all = TRUE)
  }
  dev.off()

  # ------------------------------------------------------------------------
  # Convert chromosomal positions to absolute positions for plotting the whole genome
  segs <- main_seg
  chrom_ends <- cumsum(as.numeric(by(segs$endpos, segs$chr, max)))
  chrom_starts <- c(0, chrom_ends[-length(chrom_ends)])
  names(chrom_starts) <- names(chrom_ends)
  chrom_center <- chrom_starts + (chrom_ends - chrom_starts) / 2
  segs$chrom_startpos <- chrom_starts[segs$chr]

  # Chromosome start and end positions in genome-wide absolute positions
  chrom_info <- data.frame(
    start = chrom_starts,
    end = chrom_ends,
    center = chrom_center,
    label = as.character(unique(segs$chr))
  )

  # Whole genome plot
  output_plot <- file.path(conf$output_dir, "whole_genome.pdf")
  pdf(output_plot, width = 10, height = 4 * length(snps))
  par(cex = 0.8)
  split.screen(c(length(snps), 1))
  num_subfigures <- 3

  for (sample_idx in seq_along(sample_names)) {
    sample_name <- sample_names[[sample_idx]]
    data_sample <- snps[[sample_name]]
    segs_sample <- segs[segs$sample == sample_name, ]

    # Split the plotting area vertically (sample name + info on left, plots on right)
    screen(sample_idx)
    split_coords <- matrix(c(0, 0.1, 0, 1, 0.1, 1, 0, 1), ncol = 4, byrow = TRUE)
    screen_idxs <- split.screen(split_coords)

    screen(screen_idxs[[1]])
    plot_sample_info(sample_name, sample_info[sample_name, "ploidy"], sample_info[sample_name, "purity"])

    screen(screen_idxs[[2]])
    split_coords <- matrix(c(0, 1, 0.7, 1,
                             0, 1, 0.4, 0.7,
                             0, 1, 0.0, 0.4), ncol = 4, byrow = TRUE)
    subfigure_idx <- split.screen(split_coords)

    # Remove the many NA regions
    segs_sample <- segs_sample[!is.na(segs_sample$cpn_total), ]

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

    # Copy number track
    screen(subfigure_idx[[3]])
    plot_track_cn_whole_genome(segs_sample, chrom_info, xlim = c(min(data_sample$pos), max(data_sample$pos)))
  }
  close.screen(all = TRUE)
  dev.off()

  print(warnings())
}
