plotAllSampleBAFsAndMSAIAcrossGenomeForPaperFigures <- function(full.joint.segs.df, sample.vec, patient.phased.BAF.df.list, cpnA.allele.color = "#de2d2699", cpnB.allele.color = "#3182bd99") {
  samp.BAF.df <- patient.phased.BAF.df.list[[1]]
  chr.start.end.df <- samp.BAF.df[, c("chrom", "pos", "pos")]
  chr.start.end.df <- apply(chr.start.end.df, 2, as.numeric)

  chr.vec <- unique(samp.BAF.df[, "chrom"])

  chroms.length.vec <- getChromLengthsVec(chr.start.end.df)
  names(chroms.length.vec) <- chr.vec
  chrom.segs.to.plot.df <- getCumulativeChromLengthsForPlottingDF(chr.start.end.df, chroms.length.vec, chr.vec)
  min.x <- 0
  min.y <- -0.5
  max.y <- 8.5
  max.x <- as.numeric(max(chrom.segs.to.plot.df[, 3]))

  num.samples <- length(sample.vec)
  first.samp <- sample.vec[1]
  last.samp <- sample.vec[length(sample.vec)]
  # sampJoinedSeg <- funJoinAdjacentSegments(copyITHRDS$sampComparisonList) is the equivalent

  par(mfrow = c(num.samples, 1))

  for (samp in sample.vec) {

    # Set out border size for sub plot
    if (samp == first.samp) {
      par(mar = c(0, 0, 0, 0))
    } else if (samp == last.samp) {
      par(mar = c(0, 0, 0, 0))
    } else {
      par(mar = c(0, 0, 0, 0))
    }

    samp.BAF.df <- patient.phased.BAF.df.list[[samp]]
    samp.BAF.df$color_to_plot <- ifelse(samp.BAF.df$maj_is_cpnA_logic_vec == 1, cpnA.allele.color, cpnB.allele.color)
    samp.BAF.df$mphase_segmented_BAF <- 0.5

    if (nrow(samp.BAF.df[grepl("none", samp.BAF.df$sample_phased_by_this_sample), ]) > 0) {
      samp.BAF.df[grepl("none", samp.BAF.df$sample_phased_by_this_sample), ]$ph_median_segmented_BAF <- 0.5
    }
    rownames(samp.BAF.df) <- 1:nrow(samp.BAF.df)
    plot(chrom.segs.to.plot.df[, 2], samp.BAF.df$BAF,
      col = samp.BAF.df$color_to_plot,
      pch = 16,
      ylim = c(0, 1),
      cex = 1,
      main = samp,
      ylab = "BAF",
      xlab = "",
      las = 1,
      xaxt = "n"
    )

    mtext(
      side = 3,
      at = getChromMidpointsVec(chrom.segs.to.plot.df, unique(chr.start.end.df[, 1])),
      text = sort(unique(chrom.segs.to.plot.df[, 1])),
      cex = seq(0.6, 0.4, length.out = length(unique(chrom.segs.to.plot.df[, 1]))),
      line = -1
    )

    for (chr in sort(unique(as.numeric(chrom.segs.to.plot.df[, 1])))) {
      abline(v = max(as.numeric(chrom.segs.to.plot.df[chrom.segs.to.plot.df[, 1] == chr, 3])), lwd = 0.5)
    }

    # points(chrom.segs.to.plot.df[, 2], samp.BAF.df$BAFsegmented, col = 'black', ylim = c(0 , 1), pch = 15, cex = 0.25)
    # points(chrom.segs.to.plot.df[, 2], samp.BAF.df$ph_median_segmented_BAF, col = '#78c679', ylim = c(0 , 1), pch = 15, cex = 0.5)

    # points(chrom.segs.to.plot.df[, 2], 1 - samp.BAF.df$BAFsegmented, col = 'black', ylim = c(0 , 1), pch = 15, cex = 0.25)
    # oints(chrom.segs.to.plot.df[, 2], 1 - samp.BAF.df$ph_median_segmented_BAF, col = '#78c679', ylim = c(0 , 1), pch = 15, cex = 0.5)

    box()
  }
}


plotAllSampleBAFsAcrossGenome <- function(sample.vec, patient.phased.BAF.df.list, cpnA.allele.color = "#de2d2699", cpnB.allele.color = "#3182bd99") {
  samp.BAF.df <- patient.phased.BAF.df.list[[1]]
  chr.start.end.df <- samp.BAF.df[, c("chrom", "pos", "pos")]
  chr.start.end.df <- apply(chr.start.end.df, 2, as.numeric)

  chr.vec <- unique(samp.BAF.df[, "chrom"])

  chroms.length.vec <- getChromLengthsVec(chr.start.end.df)
  names(chroms.length.vec) <- chr.vec
  chrom.segs.to.plot.df <- getCumulativeChromLengthsForPlottingDF(chr.start.end.df, chroms.length.vec, chr.vec)
  min.x <- 0
  min.y <- -0.5
  max.y <- 8.5
  max.x <- as.numeric(max(chrom.segs.to.plot.df[, 3]))

  num.samples <- length(sample.vec)
  first.samp <- sample.vec[1]
  last.samp <- sample.vec[length(sample.vec)]
  # sampJoinedSeg <- funJoinAdjacentSegments(copyITHRDS$sampComparisonList) is the equivalent

  par(mfrow = c(num.samples, 1))

  for (samp in sample.vec) {

    # Set out border size for sub plot
    if (samp == first.samp) {
      par(mar = c(2, 3, 4, 3))
    } else if (samp == last.samp) {
      par(mar = c(4, 3, 2, 3))
    } else {
      par(mar = c(3, 3, 3, 3))
    }

    samp.BAF.df <- patient.phased.BAF.df.list[[samp]]
    samp.BAF.df$color_to_plot <- ifelse(samp.BAF.df$maj_is_cpnA_logic_vec == 1, cpnA.allele.color, cpnB.allele.color)
    samp.BAF.df$mphase_segmented_BAF <- 0.5

    if (nrow(samp.BAF.df[grepl("none", samp.BAF.df$sample_phased_by_this_sample), ]) > 0) {
      samp.BAF.df[grepl("none", samp.BAF.df$sample_phased_by_this_sample), ]$ph_median_segmented_BAF <- 0.5
    }

    rownames(samp.BAF.df) <- 1:nrow(samp.BAF.df)
    plot(chrom.segs.to.plot.df[, 2], samp.BAF.df$BAF,
      col = samp.BAF.df$color_to_plot,
      pch = 16,
      ylim = c(0, 1),
      cex = 1,
      main = samp,
      ylab = "BAF",
      xlab = "",
      las = 1,
      xaxt = "n"
    )

    mtext(
      side = 3,
      at = getChromMidpointsVec(chrom.segs.to.plot.df, unique(chr.start.end.df[, 1])),
      text = sort(unique(chrom.segs.to.plot.df[, 1])),
      cex = seq(0.6, 0.4, length.out = length(unique(chrom.segs.to.plot.df[, 1]))),
      line = -1
    )

    for (chr in sort(unique(as.numeric(chrom.segs.to.plot.df[, 1])))) {
      abline(v = max(as.numeric(chrom.segs.to.plot.df[chrom.segs.to.plot.df[, 1] == chr, 3])), lwd = 0.5)
    }

    points(chrom.segs.to.plot.df[, 2], samp.BAF.df$BAFsegmented, col = "black", ylim = c(0, 1), pch = 15, cex = 0.25)
    points(chrom.segs.to.plot.df[, 2], samp.BAF.df$ph_median_segmented_BAF, col = "#78c679", ylim = c(0, 1), pch = 15, cex = 0.5)

    points(chrom.segs.to.plot.df[, 2], 1 - samp.BAF.df$BAFsegmented, col = "black", ylim = c(0, 1), pch = 15, cex = 0.25)
    points(chrom.segs.to.plot.df[, 2], 1 - samp.BAF.df$ph_median_segmented_BAF, col = "#78c679", ylim = c(0, 1), pch = 15, cex = 0.5)

    box()
  }
}


plotOldNewComparisonOfAllCPN <- function(sample.vec, full.joint.segs.df, plot.output.dir) {
  for (samp in sample.vec) {
    png(file = paste(plot.output.dir, samp, "_orig_vs_mphase_raw_cpn_comparison", ".png", sep = ""), width = 1200, height = 1400)
    par(mfrow = c(1, 2))

    samp.proc.segs.df <- full.joint.segs.df[full.joint.segs.df$sample %in% samp, , drop = FALSE]
    # defineYmax and Xmax
    max.cpnA.val <- min(c(5, ceiling(max(c(as.numeric(samp.proc.segs.df$nAraw), as.numeric(samp.proc.segs.df$orig_cpnA_raw)), na.rm = TRUE))))
    max.cpnB.val <- min(c(5, ceiling(max(c(as.numeric(samp.proc.segs.df$nBraw), as.numeric(samp.proc.segs.df$orig_cpnB_raw)), na.rm = TRUE))))
    plot(ifelse(as.numeric(samp.proc.segs.df$nAraw) >= max.cpnA.val, max.cpnA.val, as.numeric(samp.proc.segs.df$nAraw)),
      ifelse(as.numeric(samp.proc.segs.df$orig_cpnA_raw) >= max.cpnA.val, max.cpnA.val, as.numeric(samp.proc.segs.df$orig_cpnA_raw)),
      pch = 16, ,
      cex = log(as.numeric(samp.proc.segs.df$n.het), 5), ,
      col = ifelse(samp.proc.segs.df$sc_cpnA_present == "TRUE", "#de2d2699", "#3182bd99"),
      las = 1,
      xlab = "new cpnA raw",
      ylab = "original cpnA raw",
      xlim = c(0, max.cpnA.val),
      ylim = c(0, max.cpnA.val),
      main = samp
    )

    abline(0, 1)
    abline(h = seq(0, 20, by = 1), lty = "dashed")
    abline(v = seq(0, 20, by = 1), lty = "dashed")
    par()
    plot(ifelse(as.numeric(samp.proc.segs.df$nBraw) >= max.cpnB.val, max.cpnB.val, as.numeric(samp.proc.segs.df$nBraw)),
      ifelse(as.numeric(samp.proc.segs.df$orig_cpnB_raw) >= max.cpnB.val, max.cpnB.val, as.numeric(samp.proc.segs.df$orig_cpnB_raw)),
      pch = 16, cex = log(as.numeric(samp.proc.segs.df$n.het), 10),
      col = ifelse(samp.proc.segs.df$sc_cpnB_present == "TRUE", "#de2d2699", "#3182bd99"),
      las = 1,
      xlab = "new cpnB raw",
      ylab = "original cpnB raw",
      xlim = c(0, max.cpnB.val),
      ylim = c(0, max.cpnB.val),
      main = samp
    )

    abline(h = seq(0, 20, by = 1), lty = "dashed")
    abline(v = seq(0, 20, by = 1), lty = "dashed")
    abline(0, 1)

    dev.off()
  }
}


plotAllSamplePairwiseComparisonSegmentCPN <- function(sample.vec, full.joint.segs.df, plot.output.dir) {
  # next, compare each and every samp to each other
  par(mfrow = c(2, 2))
  for (x in 1:(length(sample.vec) - 1)) {
    for (y in (x + 1):length(sample.vec)) {
      x.axis.samp <- sample.vec[x]
      y.axis.samp <- sample.vec[y]

      png(file = paste(plot.output.dir, y.axis.samp, "_vs_", x.axis.samp, "_cpnA_cpnB_comparison", ".png", sep = ""), width = 1200, height = 1400)
      par(mfrow = c(2, 2))

      x.axis.samp.proc.segs.df <- full.joint.segs.df[full.joint.segs.df$sample %in% x.axis.samp, , drop = FALSE]
      y.axis.samp.proc.segs.df <- full.joint.segs.df[full.joint.segs.df$sample %in% y.axis.samp, , drop = FALSE]
      max.cpnA.val <- min(c(5, ceiling(max(c(as.numeric(x.axis.samp.proc.segs.df$nAraw), as.numeric(y.axis.samp.proc.segs.df$nAraw)), na.rm = TRUE))))
      max.cpnB.val <- min(c(5, ceiling(max(c(as.numeric(x.axis.samp.proc.segs.df$nBraw), as.numeric(y.axis.samp.proc.segs.df$nBraw)), na.rm = TRUE))))
      x.axis.samp.proc.segs.df$nAraw <- ifelse(as.numeric(x.axis.samp.proc.segs.df$nAraw) > max.cpnA.val, max.cpnA.val, as.numeric(x.axis.samp.proc.segs.df$nAraw))
      y.axis.samp.proc.segs.df$nAraw <- ifelse(as.numeric(y.axis.samp.proc.segs.df$nAraw) > max.cpnA.val, max.cpnA.val, as.numeric(y.axis.samp.proc.segs.df$nAraw))
      x.axis.samp.proc.segs.df$nBraw <- ifelse(as.numeric(x.axis.samp.proc.segs.df$nBraw) > max.cpnB.val, max.cpnB.val, as.numeric(x.axis.samp.proc.segs.df$nBraw))
      y.axis.samp.proc.segs.df$nBraw <- ifelse(as.numeric(y.axis.samp.proc.segs.df$nBraw) > max.cpnB.val, max.cpnB.val, as.numeric(y.axis.samp.proc.segs.df$nBraw))


      plot(x.axis.samp.proc.segs.df$ph_raw_cpnA,
        y.axis.samp.proc.segs.df$ph_raw_cpnA,
        pch = 16,
        cex = log(as.numeric(x.axis.samp.proc.segs.df$n.het), 10),
        col = ifelse(x.axis.samp.proc.segs.df$sc_cpnA_present == "TRUE" | y.axis.samp.proc.segs.df$sc_cpnA_present == "TRUE", "#de2d2699", "#3182bd99"),
        las = 1,
        xlab = paste(x.axis.samp, "Araw"),
        ylab = paste(y.axis.samp, "Araw"),
        main = paste(y.axis.samp, x.axis.samp),
        xlim = c(0, max.cpnA.val),
        ylim = c(0, max.cpnA.val)
      )
      abline(0, 1)
      abline(h = seq(0, 20, by = 1), lty = "dashed")
      abline(v = seq(0, 20, by = 1), lty = "dashed")

      plot(x.axis.samp.proc.segs.df$ph_raw_cpnB,
        y.axis.samp.proc.segs.df$ph_raw_cpnB,
        pch = 16,
        cex = log(as.numeric(x.axis.samp.proc.segs.df$n.het), 10),
        col = ifelse(x.axis.samp.proc.segs.df$sc_cpnB_present == "TRUE" | y.axis.samp.proc.segs.df$sc_cpnB_present == "TRUE", "#de2d2699", "#3182bd99"),
        las = 1,
        xlab = paste(x.axis.samp, "Braw"),
        ylab = paste(y.axis.samp, "Braw"),
        main = paste(y.axis.samp, x.axis.samp),
        xlim = c(0, max.cpnB.val),
        ylim = c(0, max.cpnB.val)
      )
      abline(0, 1)
      abline(h = seq(0, 20, by = 1), lty = "dashed")
      abline(v = seq(0, 20, by = 1), lty = "dashed")

      plot(x.axis.samp.proc.segs.df$frac_psc_cpnA_high,
        y.axis.samp.proc.segs.df$frac_psc_cpnA_high,
        pch = 16,
        cex = log(as.numeric(x.axis.samp.proc.segs.df$n.het), 10),
        col = ifelse(x.axis.samp.proc.segs.df$sc_cpnA_present == "TRUE" | y.axis.samp.proc.segs.df$sc_cpnA_present == "TRUE", "#de2d2699", "#3182bd99"),
        las = 1,
        xlab = paste(x.axis.samp, "frac_psc_cpnA_high"),
        ylab = paste(y.axis.samp, "frac_psc_cpnA_high"),
        main = paste(y.axis.samp, x.axis.samp),
        xlim = c(0, 1.5),
        ylim = c(0, 1.5)
      )
      abline(0, 1)
      abline(h = seq(0, 20, by = 1), lty = "dashed")
      abline(v = seq(0, 20, by = 1), lty = "dashed")

      plot(x.axis.samp.proc.segs.df$frac_psc_cpnB_high,
        y.axis.samp.proc.segs.df$frac_psc_cpnB_high,
        pch = 16,
        cex = log(as.numeric(x.axis.samp.proc.segs.df$n.het), 10),
        col = ifelse(x.axis.samp.proc.segs.df$sc_cpnB_present == "TRUE" | y.axis.samp.proc.segs.df$sc_cpnB_present == "TRUE", "#de2d2699", "#3182bd99"),
        las = 1,
        xlab = paste(x.axis.samp, "frac_psc_cpnB_high"),
        ylab = paste(y.axis.samp, "frac_psc_cpnB_high"),
        main = paste(y.axis.samp, x.axis.samp),
        xlim = c(0, 1.5),
        ylim = c(0, 1.5)
      )
      abline(0, 1)
      abline(h = seq(0, 20, by = 1), lty = "dashed")
      abline(v = seq(0, 20, by = 1), lty = "dashed")

      dev.off()
    }
  }
}


# These small functions should bve quick to rewrite:
# was fun.chrom.length
getChromLengthsVec <- function(seg) {
  # This function returns the length of chromosome 1,...,22
  # seg: first three columns of the minimum conistent samp matrix
  chrom.length <- c()

  for (i in 1:22) {
    sub.chrom <- subset(seg, seg[, 1] == i)
    chrom.length <- c(chrom.length, sub.chrom[nrow(sub.chrom), 3])
  }

  return(chrom.length)
}


# Was fun.add.chrom
getCumulativeChromLengthsForPlottingDF <- function(seg, chrom.length, chr.vec) {
  # This function adds the length of chromosome 1,...,(i - 1) to chromsomes 2,...,22
  # seg: first three columns of the minimum conistent samp matrix
  # chrom.length: output of fun.chrom.length(...)
  # print("called getCumulativeChromLengthsForPlottingDF")
  # print(str(seg, max=1))

  if (is.null(nrow(seg))) {
    seg <- matrix(seg, ncol = 3)
  }

  chr.vec <- as.numeric(chr.vec)
  start.chr <- chr.vec[1]

  if (start.chr %in% seg[, 1]) {
    seg.tmp <- subset(seg, seg[, 1] == start.chr)

    if (length(chr.vec) == 1) {
      return(seg.tmp)
    }
  } else {
    seg.tmp <- matrix(NA, nr = 0, nc = ncol(seg))
  }

  for (i in 2:length(chr.vec)) {
    if (!(chr.vec[i] %in% seg[, 1])) {
      next
    }

    sub.chrom <- subset(seg, seg[, 1] == chr.vec[i])
    sub.chrom[, 2] <- sub.chrom[, 2] + sum(chrom.length[1:(chr.vec[i] - 1)])
    sub.chrom[, 3] <- sub.chrom[, 3] + sum(chrom.length[1:(chr.vec[i] - 1)])
    seg.tmp <- rbind(seg.tmp, sub.chrom)
  }

  return(seg.tmp)
}


# was fun.chrom.mean
getChromMidpointsVec <- function(seg, chr.vec) {
  # This function finds the middle of the positions for each chromosome.
  # seg: Matrix with three columns, containing the positions.
  chrom.mean <- c()

  for (i in chr.vec) {
    sub.chrom <- subset(seg, seg[, 1] == i)
    chrom.mean <- c(chrom.mean, (min(sub.chrom[, 2]) + max(sub.chrom[, 3])) / 2)
  }

  return(chrom.mean)
}


# Used for plottin gthe copy nmber states alone.
plotCPNSegDFAccrossGenome <- function(segs.df,
                                      sample.vec,
                                      cpn.cols.to.plot = c("int_ph_cpnA"),
                                      color.vec = c("black")) {
  stopifnot(all(cpn.cols.to.plot %in% colnames(segs.df)))

  samp <- sample.vec[1]
  first.segs.df <- segs.df[segs.df$sample == samp, ]
  first.chr.start.end.df <- first.segs.df[, c("chr", "startpos", "endpos")]
  first.chr.start.end.df <- apply(first.chr.start.end.df, 2, as.numeric)
  chr.vec <- unique(first.segs.df[, "chr"])
  first.chroms.length.vec <- getChromLengthsVec(first.chr.start.end.df)
  names(first.chroms.length.vec) <- unique(first.chr.start.end.df[, 1])

  first.chrom.segs.to.plot.df <- getCumulativeChromLengthsForPlottingDF(first.chr.start.end.df, first.chroms.length.vec, chr.vec)
  min.x <- 0
  min.y <- -1.05
  max.y <- 8.5
  max.x <- as.numeric(max(first.chrom.segs.to.plot.df[, 3]))

  num.samples <- length(sample.vec)
  first.samp <- sample.vec[1]
  last.samp <- sample.vec[length(sample.vec)]

  par(mfrow = c(num.samples, 1))

  # Set out border size for sub plot
  for (samp in sample.vec) {
    if (samp == first.samp) {
      par(mar = c(2, 3, 4, 3))
    } else if (samp == last.samp) {
      par(mar = c(4, 3, 2, 3))
    } else {
      par(mar = c(3, 3, 3, 3))
    }

    plot(1,
      xlim = c(min.x, max.x),
      ylim = c(min.y, max.y),
      type = "n",
      xaxt = "n",
      yaxs = "i",
      yaxt = "n",
      xaxs = "i",
      ylab = "",
      xlab = "",
      lwd = 2,
      main = paste0(samp, " phased allele specific copy number"),
      bty = "n"
    )
    box()

    samp.seg.df <- segs.df[segs.df$sample == samp, ]
    samp.chr.start.end.df <- samp.seg.df[, c("chr", "startpos", "endpos")]
    samp.chr.start.end.df <- apply(samp.chr.start.end.df, 2, as.numeric)
    samp.chroms.length.vec <- getChromLengthsVec(samp.chr.start.end.df)
    names(samp.chroms.length.vec) <- unique(samp.chr.start.end.df[, 1])
    samp.chrom.segs.to.plot.df <- getCumulativeChromLengthsForPlottingDF(samp.chr.start.end.df, samp.chroms.length.vec, chr.vec)

    for (k in 1:nrow(samp.seg.df)) {
      for (j in 1:length(cpn.cols.to.plot)) {
        segments(
          x0 = samp.chrom.segs.to.plot.df[k, "startpos"],
          y0 = samp.seg.df[k, cpn.cols.to.plot[j]] - ((j - 1) * 0.05),
          x1 = samp.chrom.segs.to.plot.df[k, "endpos"],
          y1 = samp.seg.df[k, cpn.cols.to.plot[j]] - ((j - 1) * 0.05),
          col = color.vec[j],
          lwd = 3,
          lend = 2
        )
      }
    }

    axis(side = 2, las = 1)
    abline(h = seq((min.y + 0.5), (max.y - 0.5), by = 1), lty = "dashed", lwd = 0.75, col = "gray")

    mtext(
      side = 3,
      at = getChromMidpointsVec(first.chrom.segs.to.plot.df[, 1:3], unique(first.chr.start.end.df[, 1])),
      text = sort(unique(first.chrom.segs.to.plot.df[, 1])),
      cex = seq(0.6, 0.4, length.out = length(unique(first.chrom.segs.to.plot.df[, 1]))),
      line = -1
      # ,lwd=0.5
    )

    for (i in sort(unique(as.numeric(first.chrom.segs.to.plot.df[, 1])))) {
      abline(
        v = max(as.numeric(first.chrom.segs.to.plot.df[first.chrom.segs.to.plot.df[, 1] == i, 3])),
        lwd = 0.5
      )
    }
  }
}


plotSampleLogRBAFCPNABsAndMSAIAcrossGenome <- function(samp,
                                                       full.joint.segs.df,
                                                       orig.segs.df,
                                                       mcs.segs.df,
                                                       new.joined.segs.df,
                                                       new.mcs.segs.df = NULL,
                                                       asas.ph.raw.cpn.lists,
                                                       all.seg.CCF.vals.list,
                                                       all.seg.MSAI.by.BAF.ltnpd.present.list,
                                                       patient.phased.BAF.df.list,
                                                       asas.SNPs.w.both.lr.baf.snp.lists,
                                                       asas.SNPS.lr.and.or.baf.lists.path,
                                                       cpnA.allele.color = "#de2d2699",
                                                       cpnB.allele.color = "#3182bd99") {

  print("Called plotSampleLogRBAFCPNABsAndMSAIAcrossGenome....")
  print("params: samp")
  print(samp)
  seg.name.vec <- names(all.seg.MSAI.by.BAF.ltnpd.present.list)
  chr.vec <- as.numeric(sapply(seg.name.vec, function(x) strsplit(x, split = ":")[[1]][1]))
  start.vec <- as.numeric(sapply(seg.name.vec, function(x) strsplit(x, split = ":")[[1]][2]))
  end.vec <- as.numeric(sapply(seg.name.vec, function(x) strsplit(x, split = ":")[[1]][3]))

  chr.start.end.df <- data.frame(chrom = as.numeric(chr.vec), startpos = as.numeric(start.vec), endpos = as.numeric(end.vec), stringsAsFactors = FALSE)

  # samp.BAF.df <- patient.phased.BAF.df.list[[1]]
  # chr.start.end.df <- samp.BAF.df[, c("chrom", "pos", "pos")]
  # chr.start.end.df <- apply(chr.start.end.df, 2, as.numeric)

  chr.vec <- unique(chr.start.end.df[, "chrom"])

  chroms.length.vec <- getChromLengthsVec(chr.start.end.df)
  names(chroms.length.vec) <- chr.vec
  chrom.segs.to.plot.df <- getCumulativeChromLengthsForPlottingDF(chr.start.end.df, chroms.length.vec, chr.vec)
  rownames(chrom.segs.to.plot.df) <- seg.name.vec

  min.x <- 0
  min.y <- -0.5
  max.y <- 8.5
  max.x <- as.numeric(max(chrom.segs.to.plot.df[, 3]))

  num.samples <- length(sample.vec)
  first.samp <- sample.vec[1]
  last.samp <- sample.vec[length(sample.vec)]
  # sampJoinedSeg <- funJoinAdjacentSegments(copyITHRDS$sampComparisonList) is the equivalent
  cpnA.allele.color <- "#de2d2699"
  cpnB.allele.color <- "#3182bd99"
  # png(file = paste(PLOT_OUTPUT_DIR, patient, "_BAF_phased", ".png", sep=""), width = 1200, height = 1400)


  track.name.vec <- c("logr only sites", "logr hz sites", "BAF", "cpnA", "cpnB", "mcs allele specific cpn", "final allele specific cpn", "final mcs allele specific cpn raw", "original allele specific cpn")
  num.tracks <- length(track.name.vec)
  first.track <- track.name.vec[1]
  last.track <- track.name.vec[length(track.name.vec)]

  par(mfrow = c(num.tracks, 1))

  for (track.name in track.name.vec) {
    print(track.name)
    # Set out border size for sub plot
    if (track.name == first.track) {
      par(mar = c(2, 3, 4, 3))
    } else if (track.name == last.track) {
      par(mar = c(4, 3, 2, 3))
    } else {
      par(mar = c(3, 3, 3, 3))
    }

    if (track.name == "logr only sites") {
      print("plotting logr only sites")

      min.x <- 0
      min.y <- -2
      max.y <- 3
      max.x <- as.numeric(max(chrom.segs.to.plot.df[, 3]))
      plot(1,
        xlim = c(min.x, max.x),
        ylim = c(min.y, max.y),
        type = "n",
        xaxt = "n",
        yaxs = "i",
        yaxt = "n",
        xaxs = "i",
        ylab = "",
        xlab = "",
        lwd = 2,
        main = paste0(samp, " logR from non-hz sites"),
        bty = "n"
      )
      box()

      seg.name.vec <- names(asas.SNPs.w.both.lr.baf.snp.lists[[samp]]) # names(all.seg.CCF.vals.list)

      # print(seg.name.vec)
      for (seg.name in seg.name.vec) {
        # print(seg.name)
        # Will be NA if too few SNPs and can't be plotted, length of NA == 1, do this check to avoid warnings about checking only first element of list.


        chr.start.end.vec <- unlist(strsplit(seg.name, split = ":"))

        # print("print(str(seg.ccf.cpnX.frac.vec))")
        # print(str(seg.ccf.cpnX.frac.vec))
        seg.to.plot.orig.df <- asas.SNPS.lr.and.or.baf.lists.path[[samp]][[seg.name]]

        seg.to.plot.orig.df <- seg.to.plot.orig.df[ !seg.to.plot.orig.df$BAF %in% c(NA, 2), ]
        seg.to.plot.df <- seg.to.plot.orig.df[, c("chrom", "pos", "pos", "logR", "logR", "logR", "logRsegmented")]

        # print(head(seg.to.plot.df))
        if (nrow(seg.to.plot.df) == 0) {
          next
        } # end if(is.na(seg.to.plot.df[1,1])) {

        # toPlot <- cbind(logRBAFSample1$chrom,logRBAFSample1$pos,logRBAFSample1$pos,ABSample1,ABSample2)
        seg.to.plot.df[, 1:3] <- getCumulativeChromLengthsForPlottingDF(apply(seg.to.plot.df[, 1:3, drop = FALSE], 2, as.numeric), chroms.length.vec, chr.vec)

        seg.to.plot.df <- data.frame(seg.to.plot.df, stringsAsFactors = FALSE)
        colnames(seg.to.plot.df) <- c("chrom", "pos", "pos", "logR", "logR", "logR", "logRsegmented")


        points(
          x = as.numeric(seg.to.plot.df$pos), y = as.numeric(seg.to.plot.df$logR),
          col = "grey",
          cex = .5,
          pch = 16
        ) # ifelse(MSAIstatus==1,'#238b4550',ifelse(PValsPaired<0.001&totalPval<0.001,'#3182bd50','#fc927250')),cex=0.75,pch=16)
        par(lend = 1)

        segments(
          x0 = min(as.numeric(seg.to.plot.df$pos)),
          y0 = seg.to.plot.df[, "logRsegmented"],
          x1 = max(as.numeric(seg.to.plot.df$pos)),
          y1 = seg.to.plot.df[, "logRsegmented"],
          col = "#de2d26",
          lwd = 3
        )
      } # end   for ( seg.name in seg.name.vec ) {

      axis(side = 2, las = 1)
      abline(h = seq((min.y + 0.5), (max.y - 0.5), by = 1), lty = "dashed", lwd = 0.75, col = "gray")

      mtext(
        side = 3,
        at = getChromMidpointsVec(chrom.segs.to.plot.df[, 1:3], unique(chr.start.end.df[, 1])),
        text = sort(unique(chrom.segs.to.plot.df[, 1])),
        cex = seq(0.6, 0.4, length.out = length(unique(chrom.segs.to.plot.df[, 1]))),
        line = -1
        # ,lwd=0.5
      )

      for (i in sort(unique(as.numeric(chrom.segs.to.plot.df[, 1]))))
      {
        abline(
          v = max(as.numeric(chrom.segs.to.plot.df[chrom.segs.to.plot.df[, 1] == i, 3])),
          lwd = 0.5
        )
        # ,lty='dashed')
      } # en

      # box()
    } else if (track.name == "logr hz sites") {
      print("plotting hz logr sites")

      min.x <- 0
      min.y <- -2
      max.y <- 3
      max.x <- as.numeric(max(chrom.segs.to.plot.df[, 3]))
      plot(1,
        xlim = c(min.x, max.x),
        ylim = c(min.y, max.y),
        type = "n",
        xaxt = "n",
        yaxs = "i",
        yaxt = "n",
        xaxs = "i",
        ylab = "",
        xlab = "",
        lwd = 2,
        main = paste0(samp, " logr from hz sites"),
        bty = "n"
      )
      box()

      seg.name.vec <- names(asas.SNPs.w.both.lr.baf.snp.lists[[samp]]) # names(all.seg.CCF.vals.list)

      # print(seg.name.vec)
      for (seg.name in seg.name.vec) {
        # print(seg.name)
        # Will be NA if too few SNPs and can't be plotted, length of NA == 1, do this check to avoid warnings about checking only first element of list.


        chr.start.end.vec <- unlist(strsplit(seg.name, split = ":"))

        # print("print(str(seg.ccf.cpnX.frac.vec))")
        # print(str(seg.ccf.cpnX.frac.vec))
        seg.to.plot.orig.df <- asas.SNPs.w.both.lr.baf.snp.lists[[samp]][[seg.name]]

        seg.to.plot.df <- seg.to.plot.orig.df[, c("chrom", "pos", "pos", "logR", "logR", "logR", "logRsegmented")]

        # print(head(seg.to.plot.df))
        if (nrow(seg.to.plot.df) == 0) {
          next
        } # end if(is.na(seg.to.plot.df[1,1])) {

        # toPlot <- cbind(logRBAFSample1$chrom,logRBAFSample1$pos,logRBAFSample1$pos,ABSample1,ABSample2)
        seg.to.plot.df[, 1:3] <- getCumulativeChromLengthsForPlottingDF(apply(seg.to.plot.df[, 1:3, drop = FALSE], 2, as.numeric), chroms.length.vec, chr.vec)

        seg.to.plot.df <- data.frame(seg.to.plot.df, stringsAsFactors = FALSE)
        colnames(seg.to.plot.df) <- c("chrom", "pos", "pos", "logR", "logR", "logR", "logRsegmented")


        points(
          x = as.numeric(seg.to.plot.df$pos), y = as.numeric(seg.to.plot.df$logR),
          col = "grey",
          cex = .5,
          pch = 16
        ) # ifelse(MSAIstatus==1,'#238b4550',ifelse(PValsPaired<0.001&totalPval<0.001,'#3182bd50','#fc927250')),cex=0.75,pch=16)
        par(lend = 1)

        segments(
          x0 = min(as.numeric(seg.to.plot.df$pos)),
          y0 = seg.to.plot.df[, "logRsegmented"],
          x1 = max(as.numeric(seg.to.plot.df$pos)),
          y1 = seg.to.plot.df[, "logRsegmented"],
          col = "#de2d26",
          lwd = 3
        )
      } # end   for ( seg.name in seg.name.vec ) {

      axis(side = 2, las = 1)
      abline(h = seq((min.y + 0.5), (max.y - 0.5), by = 1), lty = "dashed", lwd = 0.75, col = "gray")

      mtext(
        side = 3,
        at = getChromMidpointsVec(chrom.segs.to.plot.df[, 1:3], unique(chr.start.end.df[, 1])),
        text = sort(unique(chrom.segs.to.plot.df[, 1])),
        cex = seq(0.6, 0.4, length.out = length(unique(chrom.segs.to.plot.df[, 1]))),
        line = -1
        # ,lwd=0.5
      )

      for (i in sort(unique(as.numeric(chrom.segs.to.plot.df[, 1]))))
      {
        abline(
          v = max(as.numeric(chrom.segs.to.plot.df[chrom.segs.to.plot.df[, 1] == i, 3])),
          lwd = 0.5
        )
        # ,lty='dashed')
      } # en

      # box()
    } else if (track.name == "BAF") {
      print("plotting BAFs...")

      samp.BAF.df <- patient.phased.BAF.df.list[[samp]]
      samp.BAF.df <- samp.BAF.df[, c("chrom", "pos", "pos", "BAF", "BAFsegmented", "logR", "logRsegmented", "samp", "sample_phased_by_this_sample", "maj_is_cpnA_logic_vec", "maj_is_cpnB_logic_vec", "ph_median_segmented_BAF")]
      samp.BAF.df.to.plot <- getCumulativeChromLengthsForPlottingDF(samp.BAF.df, chroms.length.vec, chr.vec)
      samp.BAF.df <- samp.BAF.df.to.plot

      samp.BAF.df$color_to_plot <- ifelse(samp.BAF.df$maj_is_cpnA_logic_vec == 1, cpnA.allele.color, cpnB.allele.color)
      samp.BAF.df$mphase_segmented_BAF <- 0.5

      if (nrow(samp.BAF.df[grepl("none", samp.BAF.df$sample_phased_by_this_sample), ]) > 0) {
        samp.BAF.df[grepl("none", samp.BAF.df$sample_phased_by_this_sample), ]$ph_median_segmented_BAF <- 0.5
      }

      rownames(samp.BAF.df) <- 1:nrow(samp.BAF.df)


      plot(samp.BAF.df[, 2], samp.BAF.df$BAF,
        col = samp.BAF.df$color_to_plot,
        xlim = c(min.x, max.x),
        pch = 16,
        ylim = c(0, 1),
        cex = 0.5, ,
        xaxt = "n",
        yaxs = "i",
        yaxt = "n",
        xaxs = "i",
        main = samp,
        ylab = "BAF",
        xlab = "",
        las = 1,
        xaxt = "n",
        bty = "n"
      )

      mtext(
        side = 3,
        at = getChromMidpointsVec(chrom.segs.to.plot.df, unique(chr.start.end.df[, 1])),
        text = sort(unique(chrom.segs.to.plot.df[, 1])),
        cex = seq(0.6, 0.4, length.out = length(unique(chrom.segs.to.plot.df[, 1]))),
        line = -1
      )

      msai.present.seg.color <- "red"
      non.msai.seg.color <- "yellow"
      problem.seg.color <- "black"

      reference.samp.color <- "green"
      non.reference.samp.color <- "blue"
      no.detection.possible.color <- "grey"
      mixture.phasing.color <- "pink"


      for (chr in sort(unique(as.numeric(chrom.segs.to.plot.df[, 1])))) {

        # plot the segments of MSAI detection:
        chrom.seg.name.vec <- grep(paste0("^", chr, ":"), names(all.seg.MSAI.by.BAF.ltnpd.present.list), value = TRUE)
        # print(chrom.seg.name.vec)

        for (chrom.seg.name in chrom.seg.name.vec) {

          # MSAI seg colors
          if (length(all.seg.MSAI.by.BAF.ltnpd.present.list[[chrom.seg.name]]$MSAI_present_any_sample_by_thresh) > 0) {
            # Color for the MSAI
            # print(all.seg.MSAI.by.BAF.ltnpd.present.list[[chrom.seg.name]]$MSAI_present_any_sample_by_thresh %in% TRUE)
            msai.seg.color <- ifelse(all.seg.MSAI.by.BAF.ltnpd.present.list[[chrom.seg.name]]$MSAI_present_any_sample_by_thresh %in% TRUE, msai.present.seg.color, non.msai.seg.color)
          } else { # end if ( length(all.seg.MSAI.by.BAF.ltnpd.present.list[[chrom.seg.name]]$MSAI_present_any_sample_by_thresh) > 0 ) {

            msai.seg.color <- problem.seg.color
          } # end else

          # Plot whether MSAI detected
          segments(
            x0 = chrom.segs.to.plot.df[chrom.seg.name, 2],
            y0 = 0,
            x1 = chrom.segs.to.plot.df[chrom.seg.name, 3],
            y1 = 0,
            col = msai.seg.color,
            lwd = 6
          )

          # Plot if reference
          if (length(all.seg.MSAI.by.BAF.ltnpd.present.list[[chrom.seg.name]]$MSAI_detection_feasible) > 0) {
            if (all.seg.MSAI.by.BAF.ltnpd.present.list[[chrom.seg.name]]$MSAI_detection_feasible %in% TRUE) {
              seg.msai.df <- all.seg.MSAI.by.BAF.ltnpd.present.list[[chrom.seg.name]]$MSAI_by_sample_df

              if (grepl("mixture", all.seg.MSAI.by.BAF.ltnpd.present.list[[chrom.seg.name]]$MSAI_detection_info)) {
                ref.seg.color <- mixture.phasing.color
              } else if (seg.msai.df[seg.msai.df$sample_name == samp, ]$sample_is_ref == TRUE) {
                ref.seg.color <- reference.samp.color
              } else {
                ref.seg.color <- non.reference.samp.color
              }
            } else {
              ref.seg.color <- no.detection.possible.color
            }
          } else {
            ref.seg.color <- problem.seg.color
          }
          segments(
            x0 = chrom.segs.to.plot.df[chrom.seg.name, 2],
            y0 = 0.05,
            x1 = chrom.segs.to.plot.df[chrom.seg.name, 3],
            y1 = 0.05,
            col = ref.seg.color,
            lwd = 6
          )
        }

        abline(v = max(as.numeric(chrom.segs.to.plot.df[chrom.segs.to.plot.df[, 1] == chr, 3])), lwd = 0.5)

        points(samp.BAF.df[, 2], samp.BAF.df$BAFsegmented, col = "black", ylim = c(0, 1), pch = 15, cex = 0.25)
        points(samp.BAF.df[, 2], samp.BAF.df$ph_median_segmented_BAF, col = "#78c679", ylim = c(0, 1), pch = 15, cex = 0.5)

        points(samp.BAF.df[, 2], 1 - samp.BAF.df$BAFsegmented, col = "black", ylim = c(0, 1), pch = 15, cex = 0.25)
        points(samp.BAF.df[, 2], 1 - samp.BAF.df$ph_median_segmented_BAF, col = "#78c679", ylim = c(0, 1), pch = 15, cex = 0.5)
      } # end for ( chr in sort(unique(as.numeric(chrom.segs.to.plot.df[, 1])))) {
      box()
    } else if (track.name == "cpnA" | track.name == "cpnB") {
      cpnX <- track.name
      print("Plotting ph cpn")
      # "B_MET004_MT-M1"
      Colours <- colorpanel(100, "steelblue4", "#deebf7", "steelblue4")

      # filter out the segments that are unphased
      phased.joint.segs.df <- full.joint.segs.df[full.joint.segs.df$samp_used_to_phase != "too_few_SNPs", ]

      min.x <- 0
      min.y <- -0.5
      max.y <- 8.5
      max.x <- as.numeric(max(chrom.segs.to.plot.df[, 3]))

      plot(1,
        xlim = c(min.x, max.x),
        ylim = c(min.y, max.y),
        type = "n",
        xaxt = "n",
        yaxs = "i",
        yaxt = "n",
        xaxs = "i",
        ylab = "",
        xlab = "",
        lwd = 2,
        main = paste0(samp, " ", cpnX, " allele"),
        bty = "n"
      )
      box()

      samp.vals.list <- asas.ph.raw.cpn.lists[[samp]]
      samp.one.vals.list <- asas.ph.raw.cpn.lists[[first.samp]]
      seg.name.vec <- phased.joint.segs.df$seg_name # names(all.seg.CCF.vals.list)
      # print(seg.name.vec)
      for (seg.name in seg.name.vec) {
        # print(seg.name)
        # Will be NA if too few SNPs and can't be plotted, length of NA == 1, do this check to avoid warnings about checking only first element of list.
        if (length(all.seg.CCF.vals.list[[cpnX]][[seg.name]]) == 1) {
          next
        }

        chr.start.end.vec <- unlist(strsplit(seg.name, split = ":"))
        seg.samp.cpnX.vals.list <- samp.vals.list[[seg.name]]
        # print("str(seg.samp.cpnX.vals.list, max=1)")
        # print(str(seg.samp.cpnX.vals.list, max=1))
        seg.samp.one.cpnX.vals.list <- samp.one.vals.list[[seg.name]]
        # print("str(seg.samp.one.cpnX.vals.list, max=1)")
        # print(str(seg.samp.one.cpnX.vals.list, max=1))
        seg.ccf.cpnX.frac.vec <- all.seg.CCF.vals.list[[cpnX]][[seg.name]][[samp]][[paste0("frac_", cpnX, "_low_vec")]]

        # print("print(str(seg.ccf.cpnX.frac.vec))")
        # print(str(seg.ccf.cpnX.frac.vec))

        seg.to.plot.df <- cbind(
          seg.samp.one.cpnX.vals.list$chr_vec,
          seg.samp.one.cpnX.vals.list$pos_vec,
          seg.samp.one.cpnX.vals.list$pos_vec,
          seg.samp.cpnX.vals.list[[paste0("ph_", cpnX, "_vec")]],
          seg.samp.cpnX.vals.list[[paste0("ph_", cpnX, "_vec")]],
          seg.samp.cpnX.vals.list[[paste0("ph_", cpnX, "_vec")]],
          seg.ccf.cpnX.frac.vec
        )

        if (is.na(seg.to.plot.df[1, 1])) {
          next
        } # end if(is.na(seg.to.plot.df[1,1])) {

        # toPlot <- cbind(logRBAFSample1$chrom,logRBAFSample1$pos,logRBAFSample1$pos,ABSample1,ABSample2)
        seg.to.plot.df[, 1:3] <- getCumulativeChromLengthsForPlottingDF(apply(seg.to.plot.df[, 1:3], 2, as.numeric), chroms.length.vec, chr.vec)
        seg.to.plot.df <- data.frame(seg.to.plot.df, stringsAsFactors = FALSE)
        colnames(seg.to.plot.df) <- c("chr", "pos", "pos", "Sample1nA", "Sample1nA", "Sample1nA", "CCF")

        seg.to.plot.df$Sample1nA <- ifelse(as.numeric(seg.to.plot.df$Sample1nA) > max.y, (max.y - 0.25), as.numeric(seg.to.plot.df$Sample1nA))

        colourCode <- ifelse(round(as.numeric(seg.to.plot.df$CCF) * 100) > 100, 100, round(as.numeric(seg.to.plot.df$CCF) * 100))
        colourCode <- ifelse(colourCode == 0, 1, colourCode)


        points(
          x = as.numeric(seg.to.plot.df$pos), y = as.numeric(seg.to.plot.df$Sample1nA),
          col = Colours[colourCode],
          cex = .5,
          pch = 16
        ) # ifelse(MSAIstatus==1,'#238b4550',ifelse(PValsPaired<0.001&totalPval<0.001,'#3182bd50','#fc927250')),cex=0.75,pch=16)
        par(lend = 1)
        segColour <- as.numeric(all.seg.CCF.vals.list[[cpnX]][[seg.name]][[samp]][[paste0("mean_frac_", cpnX, "_low")]])
        segColour <- Colours[(ifelse(segColour > 1, 1, segColour) * 100)]

        # So now very quickly we can check wehter that segment has been called significantly subclonal or not:
        # Here we want to implement the two options.. either use the CCF threshold or the CPN threshold



        sc.present <- full.joint.segs.df[paste0(samp, ":", seg.name), paste0("within_samp_sc_", cpnX, "_present_using_min_cpn_diff")]


        if (sc.present == TRUE) {
          ccf.bar.color <- "yellow"
        } else {
          ccf.bar.color <- "#525252"
        }


        segments(
          x0 = min(as.numeric(seg.to.plot.df$pos)),
          y0 = all.seg.CCF.vals.list[[cpnX]][[seg.name]][[samp]][[paste0("poss_sc_", cpnX, "_low")]],
          x1 = max(as.numeric(seg.to.plot.df$pos)),
          y1 = all.seg.CCF.vals.list[[cpnX]][[seg.name]][[samp]][[paste0("poss_sc_", cpnX, "_low")]],
          col = "#de2d26",
          lwd = 3
        )

        segments(
          x0 = min(as.numeric(seg.to.plot.df$pos)),
          y0 = all.seg.CCF.vals.list[[cpnX]][[seg.name]][[samp]][[paste0("poss_sc_", cpnX, "_high")]],
          x1 = max(as.numeric(seg.to.plot.df$pos)),
          y1 = all.seg.CCF.vals.list[[cpnX]][[seg.name]][[samp]][[paste0("poss_sc_", cpnX, "_high")]],
          col = "#31a354",
          lwd = 3
        )

        segments(
          x0 = min(as.numeric(seg.to.plot.df$pos)),
          y0 = median(c(as.numeric(seg.to.plot.df$Sample1nA))),
          x1 = max(as.numeric(seg.to.plot.df$pos)), ,
          y1 = median(c(as.numeric(seg.to.plot.df$Sample1nA))),
          col = ccf.bar.color,
          lwd = 3
        )
      } # end   for ( seg.name in seg.name.vec ) {

      axis(side = 2, las = 1)
      abline(h = seq((min.y + 0.5), (max.y - 0.5), by = 1), lty = "dashed", lwd = 0.75, col = "gray")

      mtext(
        side = 3,
        at = getChromMidpointsVec(chrom.segs.to.plot.df[, 1:3], unique(chr.start.end.df[, 1])),
        text = sort(unique(chrom.segs.to.plot.df[, 1])),
        cex = seq(0.6, 0.4, length.out = length(unique(chrom.segs.to.plot.df[, 1]))),
        line = -1
        # ,lwd=0.5
      )

      for (i in sort(unique(as.numeric(chrom.segs.to.plot.df[, 1]))))
      {
        abline(
          v = max(as.numeric(chrom.segs.to.plot.df[chrom.segs.to.plot.df[, 1] == i, 3])),
          lwd = 0.5
        )
        # ,lty='dashed')
      } # en

      # box()
    } else if (track.name %in% c("mcs allele specific cpn", "final allele specific cpn", "final mcs allele specific cpn raw", "original allele specific cpn")) {
      if (track.name == "final allele specific cpn") {
        print("Plotting final allele specific segs...")
        segs.df <- new.joined.segs.df
        cpn.cols.to.plot <- c("int_ph_cpnA", "int_ph_cpnB")
      } else if (track.name == "original allele specific cpn") {
        print("Plotting original allele specific cpn...")
        segs.df <- orig.segs.df
        cpn.cols.to.plot <- c("int_ph_cpnA", "int_ph_cpnB")
      } else if (track.name == "mcs allele specific cpn") {
        print("Plotting MCS segs...")
        segs.df <- mcs.segs.df
        cpn.cols.to.plot <- c("int_ph_cpnA", "int_ph_cpnB")
      } else if (track.name == "final mcs allele specific cpn raw") {
        if (!is.null(new.mcs.segs.df)) {
          print("Plotting post-joining MCS segs...")
          segs.df <- new.mcs.segs.df
          cpn.cols.to.plot <- c("ph_raw_cpnA", "ph_raw_cpnB")
        } else {
          print("No post-joining MCS segs passed to function, skipping track")
          next
        }
      }


      color.vec <- c("green", "black")

      first.segs.df <- segs.df[segs.df$sample == samp, ]
      first.chr.start.end.df <- first.segs.df[, c("chr", "startpos", "endpos")]
      first.chr.start.end.df <- apply(first.chr.start.end.df, 2, as.numeric)
      # print("head(chr.start.end.df)")
      chr.vec <- unique(first.segs.df[, "chr"])
      # print(head(chr.start.end.df))
      first.chroms.length.vec <- getChromLengthsVec(first.chr.start.end.df)
      names(first.chroms.length.vec) <- unique(first.chr.start.end.df[, 1])

      # print("chroms.length.vec")
      # print(chroms.length.vec)
      first.chrom.segs.to.plot.df <- getCumulativeChromLengthsForPlottingDF(first.chr.start.end.df, first.chroms.length.vec, chr.vec)

      plot(1,
        xlim = c(min.x, max.x),
        ylim = c(min.y, max.y),
        type = "n",
        xaxt = "n",
        yaxs = "i",
        yaxt = "n",
        xaxs = "i",
        ylab = "",
        xlab = "",
        lwd = 2,
        main = paste0(samp, " ", track.name),
        bty = "n"
      )
      box()

      samp.seg.df <- segs.df[segs.df$sample == samp, ]
      samp.chr.start.end.df <- samp.seg.df[, c("chr", "startpos", "endpos")]
      samp.chr.start.end.df <- apply(samp.chr.start.end.df, 2, as.numeric)
      samp.chroms.length.vec <- getChromLengthsVec(samp.chr.start.end.df)
      names(samp.chroms.length.vec) <- unique(samp.chr.start.end.df[, 1])
      samp.chrom.segs.to.plot.df <- getCumulativeChromLengthsForPlottingDF(samp.chr.start.end.df, samp.chroms.length.vec, chr.vec)

      # print(head(samp.chrom.segs.to.plot.df))

      for (k in 1:nrow(samp.seg.df)) {
        for (j in 1:length(cpn.cols.to.plot)) {
          segments(
            x0 = samp.chrom.segs.to.plot.df[k, "startpos"],
            y0 = samp.seg.df[k, cpn.cols.to.plot[j]] - ((j - 1) * 0.05),
            x1 = samp.chrom.segs.to.plot.df[k, "endpos"],
            y1 = samp.seg.df[k, cpn.cols.to.plot[j]] - ((j - 1) * 0.05),
            col = color.vec[j],
            lwd = 3,
            lend = 2
          )
        } # end for(cpn.col in cpn.cols.to.plot) {
      } # end   for ( seg.name in seg.name.vec ) {

      axis(side = 2, las = 1)
      abline(h = seq((min.y + 0.5), (max.y - 0.5), by = 1), lty = "dashed", lwd = 0.75, col = "gray")

      mtext(
        side = 3,
        at = getChromMidpointsVec(first.chrom.segs.to.plot.df[, 1:3], unique(first.chr.start.end.df[, 1])),
        text = sort(unique(first.chrom.segs.to.plot.df[, 1])),
        cex = seq(0.6, 0.4, length.out = length(unique(first.chrom.segs.to.plot.df[, 1]))),
        line = -1
        # ,lwd=0.5
      )

      for (i in sort(unique(as.numeric(first.chrom.segs.to.plot.df[, 1]))))
      {
        abline(
          v = max(as.numeric(first.chrom.segs.to.plot.df[first.chrom.segs.to.plot.df[, 1] == i, 3])),
          lwd = 0.5
        )
        # ,lty='dashed')
      } # en
    }
  } # end tracname vec
}


plotAllSampleCCFsAcrossGenome <- function(asas.ph.raw.cpn.lists,
                                          full.joint.segs.df,
                                          sample.vec,
                                          seg.name.vec,
                                          all.seg.CCF.vals.list,
                                          cpnX,
                                          use.cpn.or.ccf.thresh = "CCF",
                                          use.pass.column.to.plot = FALSE) {
  Colours <- colorpanel(100, "steelblue4", "#deebf7", "steelblue4")
  # First get the consensus segmentation, probably want to filter out the smallest segments:
  # filter out the segments that are unphased
  phased.joint.segs.df <- full.joint.segs.df[full.joint.segs.df$samp_used_to_phase != "too_few_SNPs", ]
  # phased.joint.segs.df <- full.joint.segs.df
  chr.start.end.df <- phased.joint.segs.df[, c("chr", "startpos", "endpos")]
  chr.start.end.df <- apply(chr.start.end.df, 2, as.numeric)
  # print("head(chr.start.end.df)")

  chr.vec <- unique(full.joint.segs.df[, "chr"])
  # print(head(chr.start.end.df))
  chroms.length.vec <- getChromLengthsVec(chr.start.end.df)

  names(chroms.length.vec) <- unique(chr.start.end.df[, 1])

  # print("chroms.length.vec")
  # print(chroms.length.vec)
  chrom.segs.to.plot.df <- getCumulativeChromLengthsForPlottingDF(chr.start.end.df, chroms.length.vec, chr.vec)

  # print("head(chrom.segs.to.plot.df)")
  # print(head(chrom.segs.to.plot.df))
  min.x <- 0
  min.y <- -0.5
  max.y <- 8.5
  max.x <- as.numeric(max(chrom.segs.to.plot.df[, 3]))


  num.samples <- length(sample.vec)
  first.samp <- sample.vec[1]
  last.samp <- sample.vec[length(sample.vec)]
  # sampJoinedSeg <- funJoinAdjacentSegments(copyITHRDS$sampComparisonList) is the equivalent

  par(mfrow = c(num.samples, 1))

  for (samp in sample.vec) {

    # Set out border size for sub plot
    if (samp == first.samp) {
      par(mar = c(2, 3, 4, 3))
    } else if (samp == last.samp) {
      par(mar = c(4, 3, 2, 3))
    } else {
      par(mar = c(3, 3, 3, 3))
    }
    plot(1,
      xlim = c(min.x, max.x),
      ylim = c(min.y, max.y),
      type = "n",
      xaxt = "n",
      yaxs = "i",
      yaxt = "n",
      xaxs = "i",
      ylab = "",
      xlab = "",
      lwd = 2,
      main = paste0(samp, " ", cpnX, " allele"),
      bty = "n"
    )
    box()

    # "B_MET004_MT-M1"
    samp.vals.list <- asas.ph.raw.cpn.lists[[samp]]
    samp.one.vals.list <- asas.ph.raw.cpn.lists[[first.samp]]
    seg.name.vec <- phased.joint.segs.df$seg_name # names(all.seg.CCF.vals.list)
    for (seg.name in seg.name.vec) {

      # Will be NA if too few SNPs and can't be plotted, length of NA == 1, do this check to avoid warnings about checking only first element of list.
      if (length(all.seg.CCF.vals.list[[cpnX]][[seg.name]]) == 1) {
        next
      }

      chr.start.end.vec <- unlist(strsplit(seg.name, split = ":"))
      seg.samp.cpnX.vals.list <- samp.vals.list[[seg.name]]
      # print("str(seg.samp.cpnX.vals.list, max=1)")
      # print(str(seg.samp.cpnX.vals.list, max=1))
      seg.samp.one.cpnX.vals.list <- samp.one.vals.list[[seg.name]]
      # print("str(seg.samp.one.cpnX.vals.list, max=1)")
      # print(str(seg.samp.one.cpnX.vals.list, max=1))
      seg.ccf.cpnX.frac.vec <- all.seg.CCF.vals.list[[cpnX]][[seg.name]][[samp]][[paste0("frac_", cpnX, "_low_vec")]]

      # print("print(str(seg.ccf.cpnX.frac.vec))")
      # print(str(seg.ccf.cpnX.frac.vec))

      seg.to.plot.df <- cbind(
        seg.samp.one.cpnX.vals.list$chr_vec,
        seg.samp.one.cpnX.vals.list$pos_vec,
        seg.samp.one.cpnX.vals.list$pos_vec,
        seg.samp.cpnX.vals.list[[paste0("ph_", cpnX, "_vec")]],
        seg.samp.cpnX.vals.list[[paste0("ph_", cpnX, "_vec")]],
        seg.samp.cpnX.vals.list[[paste0("ph_", cpnX, "_vec")]],
        seg.ccf.cpnX.frac.vec
      )

      if (is.na(seg.to.plot.df[1, 1])) {
        next
      } # end if(is.na(seg.to.plot.df[1,1])) {

      # toPlot <- cbind(logRBAFSample1$chrom,logRBAFSample1$pos,logRBAFSample1$pos,ABSample1,ABSample2)
      seg.to.plot.df[, 1:3] <- getCumulativeChromLengthsForPlottingDF(apply(seg.to.plot.df[, 1:3], 2, as.numeric), chroms.length.vec, chr.vec)
      seg.to.plot.df <- data.frame(seg.to.plot.df, stringsAsFactors = FALSE)
      colnames(seg.to.plot.df) <- c("chr", "pos", "pos", "Sample1nA", "Sample1nA", "Sample1nA", "CCF")

      seg.to.plot.df$Sample1nA <- ifelse(as.numeric(seg.to.plot.df$Sample1nA) > max.y, (max.y - 0.25), as.numeric(seg.to.plot.df$Sample1nA))

      colourCode <- ifelse(round(as.numeric(seg.to.plot.df$CCF) * 100) > 100, 100, round(as.numeric(seg.to.plot.df$CCF) * 100))
      colourCode <- ifelse(colourCode == 0, 1, colourCode)


      points(
        x = as.numeric(seg.to.plot.df$pos), y = as.numeric(seg.to.plot.df$Sample1nA),
        col = Colours[colourCode],
        cex = .5,
        pch = 16
      ) # ifelse(MSAIstatus==1,'#238b4550',ifelse(PValsPaired<0.001&totalPval<0.001,'#3182bd50','#fc927250')),cex=0.75,pch=16)
      par(lend = 1)
      segColour <- as.numeric(all.seg.CCF.vals.list[[cpnX]][[seg.name]][[samp]][[paste0("mean_frac_", cpnX, "_low")]])
      segColour <- Colours[(ifelse(segColour > 1, 1, segColour) * 100)]

      # So now very quickly we can check wehter that segment has been called significantly subclonal or not:
      # Here we want to implement the two options.. either use the CCF threshold or the CPN threshold

      if (use.cpn.or.ccf.thresh == "CCF") {
        sc.present <- full.joint.segs.df[paste0(samp, ":", seg.name), paste0("within_samp_sc_", cpnX, "_present")]
      } else if (use.cpn.or.ccf.thresh == "CPN") {
        sc.present <- full.joint.segs.df[paste0(samp, ":", seg.name), paste0("within_samp_sc_", cpnX, "_present_using_min_cpn_diff")]

        if (sc.present & use.pass.column.to.plot == TRUE) {
          sc.present <- full.joint.segs.df[paste0(samp, ":", seg.name), paste0("within_samp_", cpnX, "_passed")]
        }
      }

      if (sc.present == TRUE) {
        ccf.bar.color <- "yellow"
      } else {
        ccf.bar.color <- "#525252"
      }

      segments(
        x0 = min(as.numeric(seg.to.plot.df$pos)),
        y0 = all.seg.CCF.vals.list[[cpnX]][[seg.name]][[samp]][[paste0("poss_sc_", cpnX, "_low")]],
        x1 = max(as.numeric(seg.to.plot.df$pos)),
        y1 = all.seg.CCF.vals.list[[cpnX]][[seg.name]][[samp]][[paste0("poss_sc_", cpnX, "_low")]],
        col = "#de2d26",
        lwd = 3,
        lend = 2
      )

      segments(
        x0 = min(as.numeric(seg.to.plot.df$pos)),
        y0 = all.seg.CCF.vals.list[[cpnX]][[seg.name]][[samp]][[paste0("poss_sc_", cpnX, "_high")]],
        x1 = max(as.numeric(seg.to.plot.df$pos)),
        y1 = all.seg.CCF.vals.list[[cpnX]][[seg.name]][[samp]][[paste0("poss_sc_", cpnX, "_high")]],
        col = "#31a354",
        lwd = 3,
        lend = 2
      )

      segments(
        x0 = min(as.numeric(seg.to.plot.df$pos)),
        y0 = median(c(as.numeric(seg.to.plot.df$Sample1nA))),
        x1 = max(as.numeric(seg.to.plot.df$pos)), y1 = median(c(as.numeric(seg.to.plot.df$Sample1nA))),
        col = ccf.bar.color,
        lwd = 3,
        lend = 2
      )
    } # end   for ( seg.name in seg.name.vec ) {

    axis(side = 2, las = 1)
    abline(h = seq((min.y + 0.5), (max.y - 0.5), by = 1), lty = "dashed", lwd = 0.75, col = "gray")

    mtext(
      side = 3,
      at = getChromMidpointsVec(chrom.segs.to.plot.df[, 1:3], unique(chr.start.end.df[, 1])),
      text = sort(unique(chrom.segs.to.plot.df[, 1])),
      cex = seq(0.6, 0.4, length.out = length(unique(chrom.segs.to.plot.df[, 1]))),
      line = -1
      # ,lwd=0.5
    )

    for (i in sort(unique(as.numeric(chrom.segs.to.plot.df[, 1]))))
    {
      abline(
        v = max(as.numeric(chrom.segs.to.plot.df[chrom.segs.to.plot.df[, 1] == i, 3])),
        lwd = 0.5
      )
      # ,lty='dashed')
    } # en
  } # end for ( samp in sample.vec ) {
}


plotAllSampleBAFsAndMSAIAcrossGenome <- function(sample.vec,
                                                 all.seg.MSAI.by.BAF.ltnpd.present.list,
                                                 patient.phased.BAF.df.list,
                                                 cpnA.allele.color = "#de2d2699",
                                                 cpnB.allele.color = "#3182bd99",
                                                 use.xphase = FALSE) {
  seg.name.vec <- names(all.seg.MSAI.by.BAF.ltnpd.present.list)
  chr.vec <- as.numeric(sapply(seg.name.vec, function(x) strsplit(x, split = ":")[[1]][1]))
  start.vec <- as.numeric(sapply(seg.name.vec, function(x) strsplit(x, split = ":")[[1]][2]))
  end.vec <- as.numeric(sapply(seg.name.vec, function(x) strsplit(x, split = ":")[[1]][3]))

  chr.start.end.df <- data.frame(chrom = as.numeric(chr.vec), startpos = as.numeric(start.vec), endpos = as.numeric(end.vec), stringsAsFactors = FALSE)

  chr.vec <- unique(chr.start.end.df[, "chrom"])

  chroms.length.vec <- getChromLengthsVec(chr.start.end.df)
  names(chroms.length.vec) <- chr.vec
  chrom.segs.to.plot.df <- getCumulativeChromLengthsForPlottingDF(chr.start.end.df, chroms.length.vec, chr.vec)
  rownames(chrom.segs.to.plot.df) <- seg.name.vec

  min.x <- 0
  min.y <- -0.5
  max.y <- 8.5
  max.x <- as.numeric(max(chrom.segs.to.plot.df[, 3]))

  num.samples <- length(sample.vec)
  first.samp <- sample.vec[1]
  last.samp <- sample.vec[length(sample.vec)]
  # sampJoinedSeg <- funJoinAdjacentSegments(copyITHRDS$sampComparisonList) is the equivalent
  cpnA.allele.color <- "#de2d2699"
  cpnB.allele.color <- "#3182bd99"

  par(mfrow = c(num.samples, 1))

  for (samp in sample.vec) {

    # Set out border size for sub plot
    message("Processing: ", samp, " (First sample: ", first.samp, ")")
    if (samp == first.samp) {
      par(mar = c(2, 3, 4, 3))
    } else if (samp == last.samp) {
      par(mar = c(4, 3, 2, 3))
    } else {
      par(mar = c(3, 3, 3, 3))
    }

    samp.BAF.df <- patient.phased.BAF.df.list[[samp]]
    samp.BAF.df <- samp.BAF.df[, c("chrom", "pos", "pos", "BAF", "BAFsegmented", "logR", "logRsegmented", "samp", "sample_phased_by_this_sample", "maj_is_cpnA_logic_vec", "maj_is_cpnB_logic_vec", "ph_median_segmented_BAF", "BAF_upper_in_ref")]
    samp.BAF.df.to.plot <- getCumulativeChromLengthsForPlottingDF(samp.BAF.df, chroms.length.vec, chr.vec)
    samp.BAF.df <- samp.BAF.df.to.plot

    if (!use.xphase) {
      samp.BAF.df$color_to_plot <- ifelse(samp.BAF.df$maj_is_cpnA_logic_vec == 1, cpnA.allele.color, cpnB.allele.color)
    } else {
      samp.BAF.df$color_to_plot <- ifelse(samp.BAF.df$BAF_upper_in_ref, cpnA.allele.color, cpnB.allele.color)
    }

    samp.BAF.df$mphase_segmented_BAF <- 0.5

    if (nrow(samp.BAF.df[grepl("none", samp.BAF.df$sample_phased_by_this_sample), ]) > 0) {
      samp.BAF.df[grepl("none", samp.BAF.df$sample_phased_by_this_sample), ]$ph_median_segmented_BAF <- 0.5
    }

    rownames(samp.BAF.df) <- 1:nrow(samp.BAF.df)
    plot(samp.BAF.df[, 2],
      samp.BAF.df$BAF,
      col = samp.BAF.df$color_to_plot,
      pch = 16,
      ylim = c(0, 1),
      cex = 1,
      main = samp,
      ylab = "BAF",
      xlab = "",
      las = 1,
      xaxt = "n"
    )

    mtext(
      side = 3,
      at = getChromMidpointsVec(chrom.segs.to.plot.df, unique(chr.start.end.df[, 1])),
      text = sort(unique(chrom.segs.to.plot.df[, 1])),
      cex = seq(0.6, 0.4, length.out = length(unique(chrom.segs.to.plot.df[, 1]))),
      line = -1
    )

    msai.present.seg.color <- "red"
    non.msai.seg.color <- "yellow"
    problem.seg.color <- "black"

    reference.samp.color <- "green"
    non.reference.samp.color <- "blue"
    no.detection.possible.color <- "grey"
    mixture.phasing.color <- "pink"


    for (chr in sort(unique(as.numeric(chrom.segs.to.plot.df[, 1])))) {

      # plot the segments of MSAI detection:
      chrom.seg.name.vec <- grep(paste0("^", chr, ":"), names(all.seg.MSAI.by.BAF.ltnpd.present.list), value = TRUE)
      # print(chrom.seg.name.vec)

      for (chrom.seg.name in chrom.seg.name.vec) {

        # MSAI seg colors
        if (length(all.seg.MSAI.by.BAF.ltnpd.present.list[[chrom.seg.name]]$MSAI_present_any_sample_by_thresh) > 0) {
          # Color for the MSAI
          # print(all.seg.MSAI.by.BAF.ltnpd.present.list[[chrom.seg.name]]$MSAI_present_any_sample_by_thresh %in% TRUE)
          msai.seg.color <- ifelse(all.seg.MSAI.by.BAF.ltnpd.present.list[[chrom.seg.name]]$MSAI_present_any_sample_by_thresh %in% TRUE, msai.present.seg.color, non.msai.seg.color)
        } else { # end if ( length(all.seg.MSAI.by.BAF.ltnpd.present.list[[chrom.seg.name]]$MSAI_present_any_sample_by_thresh) > 0 ) {

          msai.seg.color <- problem.seg.color
        } # end else

        # Plot whether MSAI detected
        segments(
          x0 = chrom.segs.to.plot.df[chrom.seg.name, 2], y0 = 0,
          x1 = chrom.segs.to.plot.df[chrom.seg.name, 3], y1 = 0,
          col = msai.seg.color, lwd = 6
        )

        # Plot if reference
        if (length(all.seg.MSAI.by.BAF.ltnpd.present.list[[chrom.seg.name]]$MSAI_detection_feasible) > 0) {
          if (all.seg.MSAI.by.BAF.ltnpd.present.list[[chrom.seg.name]]$MSAI_detection_feasible %in% TRUE) {
            seg.msai.df <- all.seg.MSAI.by.BAF.ltnpd.present.list[[chrom.seg.name]]$MSAI_by_sample_df

            if (grepl("mixture", all.seg.MSAI.by.BAF.ltnpd.present.list[[chrom.seg.name]]$MSAI_detection_info)) {
              ref.seg.color <- mixture.phasing.color
            } else if (seg.msai.df[seg.msai.df$sample_name == samp, ]$sample_is_ref == TRUE) {
              ref.seg.color <- reference.samp.color
            } else {
              ref.seg.color <- non.reference.samp.color
            }
          } else {
            ref.seg.color <- no.detection.possible.color
          }
        } else {
          ref.seg.color <- problem.seg.color
        }
        segments(
          x0 = chrom.segs.to.plot.df[chrom.seg.name, 2], y0 = 0.05,
          x1 = chrom.segs.to.plot.df[chrom.seg.name, 3], y1 = 0.05,
          col = ref.seg.color, lwd = 6
        )
      }

      abline(v = max(as.numeric(chrom.segs.to.plot.df[chrom.segs.to.plot.df[, 1] == chr, 3])), lwd = 0.5)

      points(samp.BAF.df[, 2], samp.BAF.df$BAFsegmented, col = "black", ylim = c(0, 1), pch = 15, cex = 0.25)
      points(samp.BAF.df[, 2], samp.BAF.df$ph_median_segmented_BAF, col = "#78c679", ylim = c(0, 1), pch = 15, cex = 0.5)

      points(samp.BAF.df[, 2], 1 - samp.BAF.df$BAFsegmented, col = "black", ylim = c(0, 1), pch = 15, cex = 0.25)
      points(samp.BAF.df[, 2], 1 - samp.BAF.df$ph_median_segmented_BAF, col = "#78c679", ylim = c(0, 1), pch = 15, cex = 0.5)
    } # end for ( chr in sort(unique(as.numeric(chrom.segs.to.plot.df[, 1])))) {
    box()
  } # end for ( samp in sample.vec ) {
}


plotSampleLogRBAFCPNABsAndMSAIAcrossGenomeFigure <- function(samp,
                                                             full.joint.segs.df,
                                                             orig.segs.df,
                                                             mcs.segs.df,
                                                             new.joined.segs.df,
                                                             new.mcs.segs.df = NULL,
                                                             asas.ph.raw.cpn.lists,
                                                             all.seg.CCF.vals.list,
                                                             all.seg.MSAI.by.BAF.ltnpd.present.list,
                                                             patient.phased.BAF.df.list,
                                                             asas.SNPs.w.both.lr.baf.snp.lists,
                                                             asas.SNPS.lr.and.or.baf.lists.path,
                                                             cpnA.allele.color = "#de2d2699",
                                                             cpnB.allele.color = "#3182bd99") {
  print("Called plotSampleLogRBAFCPNABsAndMSAIAcrossGenome....")
  print("params: samp")
  print(samp)
  seg.name.vec <- names(all.seg.MSAI.by.BAF.ltnpd.present.list)
  chr.vec <- as.numeric(sapply(seg.name.vec, function(x) strsplit(x, split = ":")[[1]][1]))
  start.vec <- as.numeric(sapply(seg.name.vec, function(x) strsplit(x, split = ":")[[1]][2]))
  end.vec <- as.numeric(sapply(seg.name.vec, function(x) strsplit(x, split = ":")[[1]][3]))

  chr.start.end.df <- data.frame(chrom = as.numeric(chr.vec), startpos = as.numeric(start.vec), endpos = as.numeric(end.vec), stringsAsFactors = FALSE)

  chr.vec <- unique(chr.start.end.df[, "chrom"])

  chroms.length.vec <- getChromLengthsVec(chr.start.end.df)
  names(chroms.length.vec) <- chr.vec
  chrom.segs.to.plot.df <- getCumulativeChromLengthsForPlottingDF(chr.start.end.df, chroms.length.vec, chr.vec)
  rownames(chrom.segs.to.plot.df) <- seg.name.vec

  min.x <- 0
  min.y <- -0.5
  max.y <- 6.5
  max.x <- as.numeric(max(chrom.segs.to.plot.df[, 3]))

  num.samples <- length(sample.vec)
  first.samp <- sample.vec[1]
  last.samp <- sample.vec[length(sample.vec)]
  # sampJoinedSeg <- funJoinAdjacentSegments(copyITHRDS$sampComparisonList) is the equivalent
  cpnA.allele.color <- "#de2d2699"
  cpnB.allele.color <- "#3182bd99"

  track.name.vec <- c("logr only sites", "logr hz sites", "BAF", "cpnA", "cpnB", "mcs allele specific cpn", "final allele specific cpn", "final mcs allele specific cpn raw", "original allele specific cpn")
  num.tracks <- length(track.name.vec)
  first.track <- track.name.vec[1]
  last.track <- track.name.vec[length(track.name.vec)]

  par(mfrow = c(num.tracks, 1))

  for (track.name in track.name.vec) {
    print(track.name)
    # Set out border size for sub plot
    if (track.name == first.track) {
      par(mar = c(0, 3, 0, 1))
    } else if (track.name == last.track) {
      par(mar = c(0, 3, 0, 1))
    } else {
      par(mar = c(0, 3, 0, 1))
    }

    if (track.name == "logr only sites") {
      print("plotting logr only sites")

      min.x <- 0
      min.y <- -2
      max.y <- 3
      max.x <- as.numeric(max(chrom.segs.to.plot.df[, 3]))
      plot(1,
        xlim = c(min.x, max.x),
        ylim = c(min.y, max.y),
        type = "n",
        xaxt = "n",
        yaxs = "i",
        yaxt = "n",
        xaxs = "i",
        ylab = "",
        xlab = "",
        lwd = 2,
        main = NULL # paste0(samp, " logR from non-hz sites")
        , bty = "n"
      )
      box()

      seg.name.vec <- names(asas.SNPs.w.both.lr.baf.snp.lists[[samp]]) # names(all.seg.CCF.vals.list)

      # print(seg.name.vec)
      for (seg.name in seg.name.vec) {
        # print(seg.name)
        # Will be NA if too few SNPs and can't be plotted, length of NA == 1, do this check to avoid warnings about checking only first element of list.


        chr.start.end.vec <- unlist(strsplit(seg.name, split = ":"))

        # print("print(str(seg.ccf.cpnX.frac.vec))")
        # print(str(seg.ccf.cpnX.frac.vec))
        seg.to.plot.orig.df <- asas.SNPS.lr.and.or.baf.lists.path[[samp]][[seg.name]]

        seg.to.plot.orig.df <- seg.to.plot.orig.df[ !seg.to.plot.orig.df$BAF %in% c(NA, 2), ]
        seg.to.plot.df <- seg.to.plot.orig.df[, c("chrom", "pos", "pos", "logR", "logR", "logR", "logRsegmented")]

        # print(head(seg.to.plot.df))
        if (nrow(seg.to.plot.df) == 0) {
          next
        } # end if(is.na(seg.to.plot.df[1,1])) {

        # toPlot <- cbind(logRBAFSample1$chrom,logRBAFSample1$pos,logRBAFSample1$pos,ABSample1,ABSample2)
        seg.to.plot.df[, 1:3] <- getCumulativeChromLengthsForPlottingDF(apply(seg.to.plot.df[, 1:3, drop = FALSE], 2, as.numeric), chroms.length.vec, chr.vec)

        seg.to.plot.df <- data.frame(seg.to.plot.df, stringsAsFactors = FALSE)
        colnames(seg.to.plot.df) <- c("chrom", "pos", "pos", "logR", "logR", "logR", "logRsegmented")


        points(
          x = as.numeric(seg.to.plot.df$pos), y = as.numeric(seg.to.plot.df$logR),
          col = "grey",
          cex = .5,
          pch = 16
        ) # ifelse(MSAIstatus==1,'#238b4550',ifelse(PValsPaired<0.001&totalPval<0.001,'#3182bd50','#fc927250')),cex=0.75,pch=16)
        par(lend = 1)

        # segments(x0 = min(as.numeric(seg.to.plot.df$pos))
        #          ,y0 = seg.to.plot.df[, "logRsegmented"]
        #          ,x1 = max(as.numeric(seg.to.plot.df$pos))
        #          ,y1 = seg.to.plot.df[, "logRsegmented"]
        #          ,col = "#de2d26"
        #          ,lwd = 3)
      } # end   for ( seg.name in seg.name.vec ) {

      # print(seg.name.vec)
      for (seg.name in seg.name.vec) {
        # print(seg.name)
        # Will be NA if too few SNPs and can't be plotted, length of NA == 1, do this check to avoid warnings about checking only first element of list.


        chr.start.end.vec <- unlist(strsplit(seg.name, split = ":"))

        # print("print(str(seg.ccf.cpnX.frac.vec))")
        # print(str(seg.ccf.cpnX.frac.vec))
        seg.to.plot.orig.df <- asas.SNPs.w.both.lr.baf.snp.lists[[samp]][[seg.name]]

        seg.to.plot.df <- seg.to.plot.orig.df[, c("chrom", "pos", "pos", "logR", "logR", "logR", "logRsegmented")]

        # print(head(seg.to.plot.df))
        if (nrow(seg.to.plot.df) == 0) {
          next
        } # end if(is.na(seg.to.plot.df[1,1])) {

        # toPlot <- cbind(logRBAFSample1$chrom,logRBAFSample1$pos,logRBAFSample1$pos,ABSample1,ABSample2)
        seg.to.plot.df[, 1:3] <- getCumulativeChromLengthsForPlottingDF(apply(seg.to.plot.df[, 1:3, drop = FALSE], 2, as.numeric), chroms.length.vec, chr.vec)

        seg.to.plot.df <- data.frame(seg.to.plot.df, stringsAsFactors = FALSE)
        colnames(seg.to.plot.df) <- c("chrom", "pos", "pos", "logR", "logR", "logR", "logRsegmented")

        segments(
          x0 = min(as.numeric(seg.to.plot.df$pos)),
          y0 = seg.to.plot.df[, "logRsegmented"],
          x1 = max(as.numeric(seg.to.plot.df$pos)),
          y1 = seg.to.plot.df[, "logRsegmented"],
          col = "#de2d26",
          lwd = 3
        )
      } # end   for ( seg.name in seg.name.vec ) {

      # axis(side = 2, las = 1)
      abline(h = seq((min.y + 0.5), (max.y - 0.5), by = 1), lty = "dashed", lwd = 0.75, col = "gray")

      # mtext(side = 3
      #       ,at = getChromMidpointsVec(chrom.segs.to.plot.df[,1:3], unique(chr.start.end.df[,1]))
      #       ,text = sort(unique(chrom.segs.to.plot.df[,1]))
      #       ,cex = seq(0.6, 0.4, length.out = length(unique(chrom.segs.to.plot.df[,1])))
      #       ,line = -1
      #       #,lwd=0.5
      # )

      for (i in sort(unique(as.numeric(chrom.segs.to.plot.df[, 1]))))
      {
        abline(
          v = max(as.numeric(chrom.segs.to.plot.df[chrom.segs.to.plot.df[, 1] == i, 3])),
          lwd = 0.5
        )
        # ,lty='dashed')
      } # en

      # box()
    } else if (track.name == "logr hz sites") {
      print("plotting hz logr sites")

      min.x <- 0
      min.y <- -2
      max.y <- 3
      max.x <- as.numeric(max(chrom.segs.to.plot.df[, 3]))
      plot(1,
        xlim = c(min.x, max.x),
        ylim = c(min.y, max.y),
        type = "n",
        xaxt = "n",
        yaxs = "i",
        yaxt = "n",
        xaxs = "i",
        ylab = "",
        xlab = "",
        lwd = 2,
        main = NULL # paste0(samp, " logr from hz sites")
        , bty = "n"
      )
      box()

      seg.name.vec <- names(asas.SNPs.w.both.lr.baf.snp.lists[[samp]]) # names(all.seg.CCF.vals.list)

      # print(seg.name.vec)
      for (seg.name in seg.name.vec) {
        # print(seg.name)
        # Will be NA if too few SNPs and can't be plotted, length of NA == 1, do this check to avoid warnings about checking only first element of list.


        chr.start.end.vec <- unlist(strsplit(seg.name, split = ":"))

        # print("print(str(seg.ccf.cpnX.frac.vec))")
        # print(str(seg.ccf.cpnX.frac.vec))
        seg.to.plot.orig.df <- asas.SNPs.w.both.lr.baf.snp.lists[[samp]][[seg.name]]

        seg.to.plot.df <- seg.to.plot.orig.df[, c("chrom", "pos", "pos", "logR", "logR", "logR", "logRsegmented")]

        # print(head(seg.to.plot.df))
        if (nrow(seg.to.plot.df) == 0) {
          next
        } # end if(is.na(seg.to.plot.df[1,1])) {

        # toPlot <- cbind(logRBAFSample1$chrom,logRBAFSample1$pos,logRBAFSample1$pos,ABSample1,ABSample2)
        seg.to.plot.df[, 1:3] <- getCumulativeChromLengthsForPlottingDF(apply(seg.to.plot.df[, 1:3, drop = FALSE], 2, as.numeric), chroms.length.vec, chr.vec)

        seg.to.plot.df <- data.frame(seg.to.plot.df, stringsAsFactors = FALSE)
        colnames(seg.to.plot.df) <- c("chrom", "pos", "pos", "logR", "logR", "logR", "logRsegmented")


        points(
          x = as.numeric(seg.to.plot.df$pos), y = as.numeric(seg.to.plot.df$logR),
          col = "grey",
          cex = .5,
          pch = 16
        ) # ifelse(MSAIstatus==1,'#238b4550',ifelse(PValsPaired<0.001&totalPval<0.001,'#3182bd50','#fc927250')),cex=0.75,pch=16)
        par(lend = 1)

        # segments(x0 = min(as.numeric(seg.to.plot.df$pos))
        #          ,y0 = seg.to.plot.df[, "logRsegmented"]
        #          ,x1 = max(as.numeric(seg.to.plot.df$pos))
        #          ,y1 = seg.to.plot.df[, "logRsegmented"]
        #          ,col = "#de2d26"
        #          ,lwd = 3)
      } # end   for ( seg.name in seg.name.vec ) {

      # print(seg.name.vec)
      for (seg.name in seg.name.vec) {
        # print(seg.name)
        # Will be NA if too few SNPs and can't be plotted, length of NA == 1, do this check to avoid warnings about checking only first element of list.


        chr.start.end.vec <- unlist(strsplit(seg.name, split = ":"))

        # print("print(str(seg.ccf.cpnX.frac.vec))")
        # print(str(seg.ccf.cpnX.frac.vec))
        seg.to.plot.orig.df <- asas.SNPs.w.both.lr.baf.snp.lists[[samp]][[seg.name]]

        seg.to.plot.df <- seg.to.plot.orig.df[, c("chrom", "pos", "pos", "logR", "logR", "logR", "logRsegmented")]

        # print(head(seg.to.plot.df))
        if (nrow(seg.to.plot.df) == 0) {
          next
        } # end if(is.na(seg.to.plot.df[1,1])) {

        # toPlot <- cbind(logRBAFSample1$chrom,logRBAFSample1$pos,logRBAFSample1$pos,ABSample1,ABSample2)
        seg.to.plot.df[, 1:3] <- getCumulativeChromLengthsForPlottingDF(apply(seg.to.plot.df[, 1:3, drop = FALSE], 2, as.numeric), chroms.length.vec, chr.vec)

        seg.to.plot.df <- data.frame(seg.to.plot.df, stringsAsFactors = FALSE)
        colnames(seg.to.plot.df) <- c("chrom", "pos", "pos", "logR", "logR", "logR", "logRsegmented")

        segments(
          x0 = min(as.numeric(seg.to.plot.df$pos)),
          y0 = seg.to.plot.df[, "logRsegmented"],
          x1 = max(as.numeric(seg.to.plot.df$pos)),
          y1 = seg.to.plot.df[, "logRsegmented"],
          col = "#de2d26",
          lwd = 3
        )
      } # end   for ( seg.name in seg.name.vec ) {

      # axis(side = 2, las = 1)
      abline(h = seq((min.y + 0.5), (max.y - 0.5), by = 1), lty = "dashed", lwd = 0.75, col = "gray")

      # mtext(side = 3
      #       ,at = getChromMidpointsVec(chrom.segs.to.plot.df[,1:3], unique(chr.start.end.df[,1]))
      #       ,text = sort(unique(chrom.segs.to.plot.df[,1]))
      #       ,cex = seq(0.6, 0.4, length.out = length(unique(chrom.segs.to.plot.df[,1])))
      #       ,line = -1
      #       #,lwd=0.5
      # )

      for (i in sort(unique(as.numeric(chrom.segs.to.plot.df[, 1]))))
      {
        abline(
          v = max(as.numeric(chrom.segs.to.plot.df[chrom.segs.to.plot.df[, 1] == i, 3])),
          lwd = 0.5
        )
        # ,lty='dashed')
      } # en

      # box()
    } else if (track.name == "BAF") {
      print("plotting BAFs...")

      samp.BAF.df <- patient.phased.BAF.df.list[[samp]]
      samp.BAF.df <- samp.BAF.df[, c("chrom", "pos", "pos", "BAF", "BAFsegmented", "logR", "logRsegmented", "samp", "sample_phased_by_this_sample", "maj_is_cpnA_logic_vec", "maj_is_cpnB_logic_vec", "ph_median_segmented_BAF")]
      samp.BAF.df.to.plot <- getCumulativeChromLengthsForPlottingDF(samp.BAF.df, chroms.length.vec, chr.vec)
      samp.BAF.df <- samp.BAF.df.to.plot

      samp.BAF.df$color_to_plot <- ifelse(samp.BAF.df$maj_is_cpnA_logic_vec == 1, cpnA.allele.color, cpnB.allele.color)
      samp.BAF.df$mphase_segmented_BAF <- 0.5

      if (nrow(samp.BAF.df[grepl("none", samp.BAF.df$sample_phased_by_this_sample), ]) > 0) {
        samp.BAF.df[grepl("none", samp.BAF.df$sample_phased_by_this_sample), ]$ph_median_segmented_BAF <- 0.5
      }

      rownames(samp.BAF.df) <- 1:nrow(samp.BAF.df)


      plot(samp.BAF.df[, 2], samp.BAF.df$BAF,
        col = samp.BAF.df$color_to_plot,
        xlim = c(min.x, max.x),
        pch = 16,
        ylim = c(0, 1),
        cex = 0.5, ,
        xaxt = "n",
        yaxs = "i",
        yaxt = "n",
        xaxs = "i",
        main = NULL, # samp,
        ylab = "BAF",
        xlab = "",
        las = 1,
        xaxt = "n",
        bty = "n"
      )

      # mtext(side = 3,
      #       at = getChromMidpointsVec(chrom.segs.to.plot.df, unique(chr.start.end.df[,1])),
      #       text = sort(unique(chrom.segs.to.plot.df[, 1])),
      #       cex = seq(0.6, 0.4, length.out = length(unique(chrom.segs.to.plot.df[, 1]))),
      #       line = -1)

      msai.present.seg.color <- "red"
      non.msai.seg.color <- "yellow"
      problem.seg.color <- "black"

      reference.samp.color <- "green"
      non.reference.samp.color <- "blue"
      no.detection.possible.color <- "grey"
      mixture.phasing.color <- "pink"


      for (chr in sort(unique(as.numeric(chrom.segs.to.plot.df[, 1])))) {

        # plot the segments of MSAI detection:
        chrom.seg.name.vec <- grep(paste0("^", chr, ":"), names(all.seg.MSAI.by.BAF.ltnpd.present.list), value = TRUE)
        # print(chrom.seg.name.vec)

        for (chrom.seg.name in chrom.seg.name.vec) {

          # MSAI seg colors
          if (length(all.seg.MSAI.by.BAF.ltnpd.present.list[[chrom.seg.name]]$MSAI_present_any_sample_by_thresh) > 0) {
            # Color for the MSAI
            # print(all.seg.MSAI.by.BAF.ltnpd.present.list[[chrom.seg.name]]$MSAI_present_any_sample_by_thresh %in% TRUE)
            msai.seg.color <- ifelse(all.seg.MSAI.by.BAF.ltnpd.present.list[[chrom.seg.name]]$MSAI_present_any_sample_by_thresh %in% TRUE, msai.present.seg.color, non.msai.seg.color)
          } else { # end if ( length(all.seg.MSAI.by.BAF.ltnpd.present.list[[chrom.seg.name]]$MSAI_present_any_sample_by_thresh) > 0 ) {

            msai.seg.color <- problem.seg.color
          } # end else

          # Plot whether MSAI detected
          segments(
            x0 = chrom.segs.to.plot.df[chrom.seg.name, 2],
            y0 = 0.02,
            x1 = chrom.segs.to.plot.df[chrom.seg.name, 3],
            y1 = 0.02,
            col = msai.seg.color,
            lwd = 3
          )

          # Plot if reference
          if (length(all.seg.MSAI.by.BAF.ltnpd.present.list[[chrom.seg.name]]$MSAI_detection_feasible) > 0) {
            if (all.seg.MSAI.by.BAF.ltnpd.present.list[[chrom.seg.name]]$MSAI_detection_feasible %in% TRUE) {
              seg.msai.df <- all.seg.MSAI.by.BAF.ltnpd.present.list[[chrom.seg.name]]$MSAI_by_sample_df

              if (grepl("mixture", all.seg.MSAI.by.BAF.ltnpd.present.list[[chrom.seg.name]]$MSAI_detection_info)) {
                ref.seg.color <- mixture.phasing.color
              } else if (seg.msai.df[seg.msai.df$sample_name == samp, ]$sample_is_ref == TRUE) {
                ref.seg.color <- reference.samp.color
              } else {
                ref.seg.color <- non.reference.samp.color
              }
            } else {
              ref.seg.color <- no.detection.possible.color
            }
          } else {
            ref.seg.color <- problem.seg.color
          }
          segments(
            x0 = chrom.segs.to.plot.df[chrom.seg.name, 2],
            y0 = 0.05,
            x1 = chrom.segs.to.plot.df[chrom.seg.name, 3],
            y1 = 0.05,
            col = ref.seg.color,
            lwd = 3
          )
        }

        abline(v = max(as.numeric(chrom.segs.to.plot.df[chrom.segs.to.plot.df[, 1] == chr, 3])), lwd = 0.5)
        abline(h = seq((min.y + 0.5), (max.y - 0.5), by = 1), lty = "dashed", lwd = 0.75, col = "gray")


        if (!(samp == "U_LTX144_SU_T1-R2" & grepl("^2:", chrom.seg.name))) {
          points(samp.BAF.df[, 2], samp.BAF.df$BAFsegmented, col = "black", ylim = c(0, 1), pch = 15, cex = 0.25)

          points(samp.BAF.df[, 2], 1 - samp.BAF.df$BAFsegmented, col = "black", ylim = c(0, 1), pch = 15, cex = 0.25)
        }
        points(samp.BAF.df[, 2], samp.BAF.df$ph_median_segmented_BAF, col = "#78c679", ylim = c(0, 1), pch = 15, cex = 0.5)
        points(samp.BAF.df[, 2], 1 - samp.BAF.df$ph_median_segmented_BAF, col = "#78c679", ylim = c(0, 1), pch = 15, cex = 0.5)
      } # end for ( chr in sort(unique(as.numeric(chrom.segs.to.plot.df[, 1])))) {
      box()
    } else if (track.name == "cpnA" | track.name == "cpnB") {
      cpnX <- track.name
      print("Plotting ph cpn")
      # "B_MET004_MT-M1"
      Colours <- colorpanel(100, "steelblue4", "#deebf7", "steelblue4")

      # filter out the segments that are unphased
      phased.joint.segs.df <- full.joint.segs.df[full.joint.segs.df$samp_used_to_phase != "too_few_SNPs", ]

      min.x <- 0
      min.y <- -0.5
      max.y <- 6.5
      max.x <- as.numeric(max(chrom.segs.to.plot.df[, 3]))

      plot(1,
        xlim = c(min.x, max.x),
        ylim = c(min.y, max.y),
        type = "n",
        xaxt = "n",
        yaxs = "i",
        yaxt = "n",
        xaxs = "i",
        ylab = "",
        xlab = "",
        lwd = 2,
        main = NULL # paste0(samp, " ", cpnX, " allele")
        , bty = "n"
      )
      box()

      samp.vals.list <- asas.ph.raw.cpn.lists[[samp]]
      samp.one.vals.list <- asas.ph.raw.cpn.lists[[first.samp]]
      seg.name.vec <- phased.joint.segs.df$seg_name # names(all.seg.CCF.vals.list)
      # print(seg.name.vec)

      # axis(side = 2, las = 1)
      abline(h = seq((min.y + 0.5), (max.y - 0.5), by = 1), lty = "dashed", lwd = 0.75, col = "gray")

      # mtext(side = 3
      #       ,at = getChromMidpointsVec(chrom.segs.to.plot.df[,1:3], unique(chr.start.end.df[,1]))
      #       ,text = sort(unique(chrom.segs.to.plot.df[,1]))
      #       ,cex = seq(0.6, 0.4, length.out = length(unique(chrom.segs.to.plot.df[,1])))
      #       ,line = -1
      #       #,lwd=0.5
      # )

      for (i in sort(unique(as.numeric(chrom.segs.to.plot.df[, 1]))))
      {
        abline(
          v = max(as.numeric(chrom.segs.to.plot.df[chrom.segs.to.plot.df[, 1] == i, 3])),
          lwd = 0.5
        )
        # ,lty='dashed')
      } # en
      for (seg.name in seg.name.vec) {
        # print(seg.name)
        # Will be NA if too few SNPs and can't be plotted, length of NA == 1, do this check to avoid warnings about checking only first element of list.
        if (length(all.seg.CCF.vals.list[[cpnX]][[seg.name]]) == 1) {
          next
        }

        chr.start.end.vec <- unlist(strsplit(seg.name, split = ":"))
        seg.samp.cpnX.vals.list <- samp.vals.list[[seg.name]]
        # print("str(seg.samp.cpnX.vals.list, max=1)")
        # print(str(seg.samp.cpnX.vals.list, max=1))
        seg.samp.one.cpnX.vals.list <- samp.one.vals.list[[seg.name]]
        # print("str(seg.samp.one.cpnX.vals.list, max=1)")
        # print(str(seg.samp.one.cpnX.vals.list, max=1))
        seg.ccf.cpnX.frac.vec <- all.seg.CCF.vals.list[[cpnX]][[seg.name]][[samp]][[paste0("frac_", cpnX, "_low_vec")]]

        # print("print(str(seg.ccf.cpnX.frac.vec))")
        # print(str(seg.ccf.cpnX.frac.vec))

        seg.to.plot.df <- cbind(
          seg.samp.one.cpnX.vals.list$chr_vec,
          seg.samp.one.cpnX.vals.list$pos_vec,
          seg.samp.one.cpnX.vals.list$pos_vec,
          seg.samp.cpnX.vals.list[[paste0("ph_", cpnX, "_vec")]],
          seg.samp.cpnX.vals.list[[paste0("ph_", cpnX, "_vec")]],
          seg.samp.cpnX.vals.list[[paste0("ph_", cpnX, "_vec")]],
          seg.ccf.cpnX.frac.vec
        )

        if (is.na(seg.to.plot.df[1, 1])) {
          next
        } # end if(is.na(seg.to.plot.df[1,1])) {

        # toPlot <- cbind(logRBAFSample1$chrom,logRBAFSample1$pos,logRBAFSample1$pos,ABSample1,ABSample2)
        seg.to.plot.df[, 1:3] <- getCumulativeChromLengthsForPlottingDF(apply(seg.to.plot.df[, 1:3], 2, as.numeric), chroms.length.vec, chr.vec)
        seg.to.plot.df <- data.frame(seg.to.plot.df, stringsAsFactors = FALSE)
        colnames(seg.to.plot.df) <- c("chr", "pos", "pos", "Sample1nA", "Sample1nA", "Sample1nA", "CCF")

        seg.to.plot.df$Sample1nA <- ifelse(as.numeric(seg.to.plot.df$Sample1nA) > max.y, (max.y - 0.25), as.numeric(seg.to.plot.df$Sample1nA))

        # colourCode <- ifelse(round(as.numeric(seg.to.plot.df$CCF) * 100) > 100, 100, round(as.numeric(seg.to.plot.df$CCF) * 100))
        colourCode <- ifelse(round(as.numeric(seg.to.plot.df$CCF) * 100) > 100, 100, 100)

        colourCode <- ifelse(colourCode == 0, 1, colourCode)


        points(
          x = as.numeric(seg.to.plot.df$pos), y = as.numeric(seg.to.plot.df$Sample1nA),
          col = Colours[colourCode],
          cex = .5,
          pch = 16
        ) # ifelse(MSAIstatus==1,'#238b4550',ifelse(PValsPaired<0.001&totalPval<0.001,'#3182bd50','#fc927250')),cex=0.75,pch=16)
        par(lend = 1)
        segColour <- as.numeric(all.seg.CCF.vals.list[[cpnX]][[seg.name]][[samp]][[paste0("mean_frac_", cpnX, "_low")]])
        segColour <- Colours[(ifelse(segColour > 1, 1, segColour) * 100)]

        # So now very quickly we can check wehter that segment has been called significantly subclonal or not:
        # Here we want to implement the two options.. either use the CCF threshold or the CPN threshold



        sc.present <- full.joint.segs.df[paste0(samp, ":", seg.name), paste0("within_samp_sc_", cpnX, "_present_using_min_cpn_diff")]


        if (sc.present == TRUE) {

          # ccf.bar.color <- "yellow"
        } else {
          ccf.bar.color <- "#525252"
        }


        # segments(x0 = min(as.numeric(seg.to.plot.df$pos))
        #          ,y0 = all.seg.CCF.vals.list[[cpnX]][[seg.name]][[samp]][[paste0("poss_sc_", cpnX, "_low")]]
        #          ,x1 = max(as.numeric(seg.to.plot.df$pos))
        #          ,y1 = all.seg.CCF.vals.list[[cpnX]][[seg.name]][[samp]][[paste0("poss_sc_", cpnX, "_low")]]
        #          ,col = "#de2d26"
        #          ,lwd = 3)

        # segments(x0 = min(as.numeric(seg.to.plot.df$pos))
        #          ,y0 = all.seg.CCF.vals.list[[cpnX]][[seg.name]][[samp]][[paste0("poss_sc_", cpnX, "_high")]]
        #          ,x1 = max(as.numeric(seg.to.plot.df$pos))
        #          ,y1 = all.seg.CCF.vals.list[[cpnX]][[seg.name]][[samp]][[paste0("poss_sc_", cpnX, "_high")]]
        #          ,col = "#31a354"
        #          ,lwd = 3)

        if (track.name == "cpnA") {
          ccf.bar.color <- "tomato1"
        } else if (track.name == "cpnB") {
          ccf.bar.color <- "steelblue2"
        }

        segments(
          x0 = min(as.numeric(seg.to.plot.df$pos)),
          y0 = median(c(as.numeric(seg.to.plot.df$Sample1nA))),
          x1 = max(as.numeric(seg.to.plot.df$pos)), ,
          y1 = median(c(as.numeric(seg.to.plot.df$Sample1nA))),
          col = ccf.bar.color,
          lwd = 3
        )
      } # end   for ( seg.name in seg.name.vec ) {

      for (seg.name in seg.name.vec) {
        # print(seg.name)
        # Will be NA if too few SNPs and can't be plotted, length of NA == 1, do this check to avoid warnings about checking only first element of list.
        if (length(all.seg.CCF.vals.list[[cpnX]][[seg.name]]) == 1) {
          next
        }

        chr.start.end.vec <- unlist(strsplit(seg.name, split = ":"))
        seg.samp.cpnX.vals.list <- samp.vals.list[[seg.name]]
        # print("str(seg.samp.cpnX.vals.list, max=1)")
        # print(str(seg.samp.cpnX.vals.list, max=1))
        seg.samp.one.cpnX.vals.list <- samp.one.vals.list[[seg.name]]
        # print("str(seg.samp.one.cpnX.vals.list, max=1)")
        # print(str(seg.samp.one.cpnX.vals.list, max=1))
        seg.ccf.cpnX.frac.vec <- all.seg.CCF.vals.list[[cpnX]][[seg.name]][[samp]][[paste0("frac_", cpnX, "_low_vec")]]

        # print("print(str(seg.ccf.cpnX.frac.vec))")
        # print(str(seg.ccf.cpnX.frac.vec))

        seg.to.plot.df <- cbind(
          seg.samp.one.cpnX.vals.list$chr_vec,
          seg.samp.one.cpnX.vals.list$pos_vec,
          seg.samp.one.cpnX.vals.list$pos_vec,
          seg.samp.cpnX.vals.list[[paste0("ph_", cpnX, "_vec")]],
          seg.samp.cpnX.vals.list[[paste0("ph_", cpnX, "_vec")]],
          seg.samp.cpnX.vals.list[[paste0("ph_", cpnX, "_vec")]],
          seg.ccf.cpnX.frac.vec
        )

        if (is.na(seg.to.plot.df[1, 1])) {
          next
        } # end if(is.na(seg.to.plot.df[1,1])) {

        # toPlot <- cbind(logRBAFSample1$chrom,logRBAFSample1$pos,logRBAFSample1$pos,ABSample1,ABSample2)
        seg.to.plot.df[, 1:3] <- getCumulativeChromLengthsForPlottingDF(apply(seg.to.plot.df[, 1:3], 2, as.numeric), chroms.length.vec, chr.vec)
        seg.to.plot.df <- data.frame(seg.to.plot.df, stringsAsFactors = FALSE)
        colnames(seg.to.plot.df) <- c("chr", "pos", "pos", "Sample1nA", "Sample1nA", "Sample1nA", "CCF")
        # So now very quickly we can check wehter that segment has been called significantly subclonal or not:
        # Here we want to implement the two options.. either use the CCF threshold or the CPN threshold
        sc.present <- full.joint.segs.df[paste0(samp, ":", seg.name), paste0("within_samp_sc_", cpnX, "_present_using_min_cpn_diff")]


        if (sc.present == TRUE) {

          # ccf.bar.color <- "yellow"
        } else {
          ccf.bar.color <- "#525252"
        }


        # segments(x0 = min(as.numeric(seg.to.plot.df$pos))
        #          ,y0 = all.seg.CCF.vals.list[[cpnX]][[seg.name]][[samp]][[paste0("poss_sc_", cpnX, "_low")]]
        #          ,x1 = max(as.numeric(seg.to.plot.df$pos))
        #          ,y1 = all.seg.CCF.vals.list[[cpnX]][[seg.name]][[samp]][[paste0("poss_sc_", cpnX, "_low")]]
        #          ,col = "#de2d26"
        #          ,lwd = 3)

        # segments(x0 = min(as.numeric(seg.to.plot.df$pos))
        #          ,y0 = all.seg.CCF.vals.list[[cpnX]][[seg.name]][[samp]][[paste0("poss_sc_", cpnX, "_high")]]
        #          ,x1 = max(as.numeric(seg.to.plot.df$pos))
        #          ,y1 = all.seg.CCF.vals.list[[cpnX]][[seg.name]][[samp]][[paste0("poss_sc_", cpnX, "_high")]]
        #          ,col = "#31a354"
        #          ,lwd = 3)

        if (track.name == "cpnA") {
          ccf.bar.color <- "tomato1"
        } else if (track.name == "cpnB") {
          ccf.bar.color <- "steelblue2"
        }

        segments(
          x0 = min(as.numeric(seg.to.plot.df$pos)),
          y0 = median(c(as.numeric(seg.to.plot.df$Sample1nA))),
          x1 = max(as.numeric(seg.to.plot.df$pos)), ,
          y1 = median(c(as.numeric(seg.to.plot.df$Sample1nA))),
          col = ccf.bar.color,
          lwd = 3
        )
      } # end   for ( seg.name in seg.name.vec ) {



      # box()
    } else if (track.name %in% c("mcs allele specific cpn", "final allele specific cpn", "final mcs allele specific cpn raw", "original allele specific cpn")) {
      if (track.name == "final allele specific cpn") {
        print("Plotting final allele specific segs...")
        segs.df <- new.joined.segs.df
        cpn.cols.to.plot <- c("int_ph_cpnA", "int_ph_cpnB")
      } else if (track.name == "original allele specific cpn") {
        print("Plotting original allele specific cpn...")
        segs.df <- orig.segs.df
        cpn.cols.to.plot <- c("int_ph_cpnA", "int_ph_cpnB")
      } else if (track.name == "mcs allele specific cpn") {
        print("Plotting MCS segs...")
        segs.df <- mcs.segs.df
        cpn.cols.to.plot <- c("int_ph_cpnA", "int_ph_cpnB")
      } else if (track.name == "final mcs allele specific cpn raw") {
        if (!is.null(new.mcs.segs.df)) {
          print("Plotting post-joining MCS segs...")
          segs.df <- new.mcs.segs.df
          cpn.cols.to.plot <- c("ph_raw_cpnA", "ph_raw_cpnB")
        } else {
          print("No post-joining MCS segs passed to function, skipping track")
          next
        }
      }


      color.vec <- c("tomato1", "steelblue2")

      first.segs.df <- segs.df[segs.df$sample == samp, ]
      first.chr.start.end.df <- first.segs.df[, c("chr", "startpos", "endpos")]
      first.chr.start.end.df <- apply(first.chr.start.end.df, 2, as.numeric)
      # print("head(chr.start.end.df)")
      chr.vec <- unique(first.segs.df[, "chr"])
      # print(head(chr.start.end.df))
      first.chroms.length.vec <- getChromLengthsVec(first.chr.start.end.df)
      names(first.chroms.length.vec) <- unique(first.chr.start.end.df[, 1])

      # print("chroms.length.vec")
      # print(chroms.length.vec)
      first.chrom.segs.to.plot.df <- getCumulativeChromLengthsForPlottingDF(first.chr.start.end.df, first.chroms.length.vec, chr.vec)

      plot(1,
        xlim = c(min.x, max.x),
        ylim = c(min.y, max.y),
        type = "n",
        xaxt = "n",
        yaxs = "i",
        yaxt = "n",
        xaxs = "i",
        ylab = "",
        xlab = "",
        lwd = 2,
        main = NULL # paste0(samp, " ", track.name)
        , bty = "n"
      )
      box()

      samp.seg.df <- segs.df[segs.df$sample == samp, ]
      samp.chr.start.end.df <- samp.seg.df[, c("chr", "startpos", "endpos")]
      samp.chr.start.end.df <- apply(samp.chr.start.end.df, 2, as.numeric)
      samp.chroms.length.vec <- getChromLengthsVec(samp.chr.start.end.df)
      names(samp.chroms.length.vec) <- unique(samp.chr.start.end.df[, 1])
      samp.chrom.segs.to.plot.df <- getCumulativeChromLengthsForPlottingDF(samp.chr.start.end.df, samp.chroms.length.vec, chr.vec)

      # print(head(samp.chrom.segs.to.plot.df))

      # axis(side = 2, las = 1)
      abline(h = seq((min.y + 0.5), (max.y - 0.5), by = 1), lty = "dashed", lwd = 0.75, col = "gray")

      # mtext(side = 3
      #       ,at = getChromMidpointsVec(first.chrom.segs.to.plot.df[,1:3], unique(first.chr.start.end.df[,1]))
      #       ,text = sort(unique(first.chrom.segs.to.plot.df[,1]))
      #       ,cex = seq(0.6, 0.4, length.out = length(unique(first.chrom.segs.to.plot.df[,1])))
      #       ,line = -1
      #       #,lwd=0.5
      #       )

      for (i in sort(unique(as.numeric(first.chrom.segs.to.plot.df[, 1]))))
      {
        abline(
          v = max(as.numeric(first.chrom.segs.to.plot.df[first.chrom.segs.to.plot.df[, 1] == i, 3])),
          lwd = 0.5
        )
        # ,lty='dashed')
      } # en


      for (k in 1:nrow(samp.seg.df)) {
        for (j in 1:length(cpn.cols.to.plot)) {
          segments(
            x0 = samp.chrom.segs.to.plot.df[k, "startpos"],
            y0 = samp.seg.df[k, cpn.cols.to.plot[j]] - ((j - 1) * 0.1),
            x1 = samp.chrom.segs.to.plot.df[k, "endpos"],
            y1 = samp.seg.df[k, cpn.cols.to.plot[j]] - ((j - 1) * 0.1),
            col = color.vec[j],
            lwd = 2.5,
            lend = 2
          )
        } # end for(cpn.col in cpn.cols.to.plot) {
      } # end   for ( seg.name in seg.name.vec ) {
    }
  } # end tracname vec
}
