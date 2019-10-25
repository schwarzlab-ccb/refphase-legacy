
cn_text <- function(x, y, labels, cex = 1.5, col = "black") {
  text(x = x, y = y, labels = labels, cex = cex * 1.1, col = "white")
  text(x = x, y = y, labels = labels, cex = cex, col = col)
}

plot_track_axis <- function(xmin, xmax) {
  par(mar = c(2, 5, 0, 1))
  plot(x=c(0, 0), y=c(xmin, xmax), xlim=c(xmin, xmax), type = "n", ylab = "", xlab = "", axes = FALSE)
  axis(1, at = NULL)
}

# Plot CNs only in this track
plot_track_cn <- function(segs, track_name, show_axis = FALSE) {
  message(track_name)
  par(mar = c(ifelse(show_axis, 2, 0.1), 5, 0.1, 1))
  plot(1, ylim=c(-1, 1), xlim=c(min(segs$startpos), max(segs$endpos)), col = "white", ylab = track_name, xaxt = ifelse(show_axis, "s", "n"), yaxt = "n")
  abline(v = unique(c(segs$startpos, segs$endpos)), lty = 2, col = "black")
  total_cn <- as.character(segs$cnTotal)
  major <- as.character(segs$nMajor)
  minor <- as.character(segs$nMinor)
  full_label <- paste(major, minor, total_cn, sep = "/")

  col <- rep("grey", nrow(segs))
  if ("is_changed" %in% colnames(segs))
    col[segs$is_changed] <- "black"
  cn_text(x = segs$seg_center, y = 0, labels = full_label, cex = 1.0, col = col)
}

plot_track_logr <- function(data_sample, segs_sample) {
  par(mar = c(0.1, 5, 2, 1))
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
  #phasing_col <- ifelse(data_sample$, "blue", "yellow")
  phasing_col[data_sample$sample_phased_by_this_sample == "none_no_allelic_imbalance_called"] <- "grey"
  phasing_col <- adjustcolor(phasing_col, 0.6)
  plot(data_sample$pos, data_sample$BAF, ylab = "BAF", col = phasing_col, ylim = c(0, 1), pch = 20, xaxt = "n", yaxt = "n")
  axis(2, seq(0, 1, by=0.25))
  #lines(x = data_sample$pos, y = data_sample$BAFsegmented, col="black")

  # Mean segments
  mean_col <- ifelse(segs_sample$orig_nMajor != segs_sample$orig_nMinor, "black", NA)
  segments(x0 = segs_sample$startpos, y0 = segs_sample$mean_BAF + 0.5, x1 = segs_sample$endpos, col = mean_col)
  segments(x0 = segs_sample$startpos, y0 = -1 * segs_sample$mean_BAF + 0.5, x1 = segs_sample$endpos, col = mean_col)

  # Segment boundaries
  abline(v = unique(c(segs_sample$startpos, segs_sample$endpos)), lty = 2, col = "black")

  # MPhase copy numbers
  major_cn <- as.character(segs_sample$orig_nMajor)
  minor_cn <- as.character(segs_sample$orig_nMinor)
  
  #major_cn <- as.character(segs_sample$nMajor)
  #minor_cn <- as.character(segs_sample$nMinor)

  cn_text(x = segs_sample$seg_center, y = 0.9, labels = major_cn)
  cn_text(x = segs_sample$seg_center, y = 0.1, labels = minor_cn)

  # Indicate which sample is the reference
  #print(head(segs_sample))
  #print(segs_sample$dominant_allele == "cpnB")
  segs_sample$dominant_allele[is.na(segs_sample$dominant_allele)] <- "na"
  if (any(segs_sample$dominant_allele == "cpnB")) {
    segments(x0 = segs_sample$startpos[segs_sample$dominant_allele == "cpnB"], y0 = 1.02,
             x1 = segs_sample$endpos[segs_sample$dominant_allele == "cpnB"], col = "darkgreen", lwd = 4, lend = 1)
  }
  cn_text(x = segs_sample$seg_center, y = 0.5, labels = ifelse(segs_sample$dominant_allele == "cpnB", "msai", ""), col = "darkgreen")

  # Indicate which sample is the reference
  if (any(segs_sample$is_reference))
    segments(x0 = segs_sample$startpos[segs_sample$is_reference], y0 = -0.02,
             x1 = segs_sample$endpos[segs_sample$is_reference], col = "red", lwd = 4, lend = 1)

}

# --------------------------------------------------------------------------
# Whole genome versions of the tracks

plot_track_logr_whole_genome <- function(data_sample, segs_sample, chrom_info) {
  par(mar = c(1, 5, 1, 0.1), xaxs = "i")
  plot(data_sample$pos, data_sample$logR, ylab = "logR", xaxt = "n", col = adjustcolor("darkgrey", 0.2), pch = ".", ylim = c(-2, 2))
  #lines(x = data_sample$pos, y = data_sample$logRsegmented, col=ifelse(is.na(data_sample$logRsegmented), "white", "black"))
  segments(x0 = segs_sample$startpos, y0 = segs_sample$median_non_hz_logr, x1 = segs_sample$endpos, col = "black")
  abline(v = chrom_info$start, lty = 2, col = "black")
  cn_text(x = chrom_info$center, y = 1, labels = chrom_info$label, cex=1.0)
}

plot_track_baf_whole_genome <- function(data_sample, segs_sample, chrom_info) {
  par(mar = c(2, 5, 0.1, 0.1), xaxs = "i")
  phasing_col <- ifelse(data_sample$maj_is_cpnB_logic_vec == 1, "blue", "orange")
  phasing_col[data_sample$sample_phased_by_this_sample == "none_no_allelic_imbalance_called"] <- "grey"
  phasing_col <- adjustcolor(phasing_col, 0.6)
  plot(data_sample$pos, data_sample$BAF, ylab = "BAF", col = phasing_col, ylim = c(0, 1), pch = ".")
  #lines(x = data_sample$pos, y = data_sample$BAFsegmented, col="black")

  # Mean segments
  mean_col <- ifelse(segs_sample$within_samp_cl_AI_present, "black", NA)
  segments(x0 = segs_sample$startpos, y0 = segs_sample$mean_BAF + 0.5, x1 = segs_sample$endpos, col = mean_col)
  segments(x0 = segs_sample$startpos, y0 = -1 * segs_sample$mean_BAF + 0.5, x1 = segs_sample$endpos, col = mean_col)

  # Chromosome boundaries
  abline(v = chrom_info$start, lty = 2, col = "black")
}

# --------------------------------------------------------------------------

plot_sample_info <- function(sample, ploidy, cellularity) {
  par(mar = c(0,0,0,0))
  plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')

  text(x = 0.5, y = 0.5, adj = 0.5,
       labels = paste0(sample, "\n", "cellularity:", round(cellularity, 2), "  ploidy:", round(ploidy, 2)),
       cex = 1.3, col = "black", srt = 90)
}

plot_chromosome_header <- function(sample, chromosome, cex = 1.5) {
  par(mar = c(0,0,0,0))
  plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')

  text(x = 0.5, y = 0.5, adj = 0.5, labels = paste0(sample, ": ", chromosome),
       cex = cex, col = "black")
}
