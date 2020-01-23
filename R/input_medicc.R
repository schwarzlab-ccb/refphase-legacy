readInMEDICCInputFormatParallel <- function(snp.file.df, processed.orig.cpn.segs.df, sample.vec) {
  # TODO add check that every single all.samp.SNP.list has the same SNPs..
  # Note error checking of removing NAs from the min.cons.seg.df should already have been performed
  # registerDoParallel(cores = num.cores)
  # print(foreach::getDoParWorkers())
  # print(getDoParName())
  # print(getDoParVersion())
  message("Calling readInMEDICCInputFormatParallel")
  message("Reading in sites data for:")
  print(sample.vec)
  message("Starting parallel work:")

  num.pre.parallel.samps <- length(unique(sample.vec))
  itx <- iterators::iter(sample.vec)
  all.samp.SNP.df <- foreach(samp = itx, .combine = rbind) %dopar% {
    ASCAT.summary <- data.frame(
      chr = as.numeric(snp.file.df$Chr),
      pos = as.numeric(snp.file.df$Pos),
      BAF = as.numeric(snp.file.df[, paste0(gsub("-", "\\.", samp), "_BAF")]),
      BAFsegmented = rep(NA, length(snp.file.df$Chr)),
      logR = as.numeric(snp.file.df[, paste0(gsub("-", "\\.", samp), "_LogR")]),
      logRsegmented = rep(NA, length(snp.file.df$Chr)),
      stringsAsFactors = FALSE
    )
    ASCAT.summary$sample <- samp

    if (nrow(processed.orig.cpn.segs.df) == 0) {
      stop("processed.orig.cpn.segs.df has no rows...")
    }
    seg <- processed.orig.cpn.segs.df[processed.orig.cpn.segs.df$sample == samp, ]
    if (nrow(processed.orig.cpn.segs.df) == 0) {
      stop("processed.orig.cpn.segs.df subsetted by sample has no rows...")
    }

    for (j in 1:nrow(seg)) {
      # Are there any BAF values that aren't set at 2 (and thus shouldn't be
      # used for allele specific cpn) present in this segment?
      if (nrow(ASCAT.summary[ASCAT.summary$chr == seg[j, ]$chr
      & ASCAT.summary$pos >= seg[j, ]$startpos
      & ASCAT.summary$pos <= seg[j, ]$endpos
      & ASCAT.summary$BAF != 2, ]) > 0) {

        # Is this segment balanced, if so set the BAFsegmented to 0.5
        if (seg[j, , ]$nMajor == seg[j, ]$nMinor) {
          ASCAT.summary[ASCAT.summary$chr == seg[j, ]$chr
          & ASCAT.summary$pos >= seg[j, ]$startpos
          & ASCAT.summary$pos <= seg[j, ]$endpos, ]$BAFsegmented <- 0.5
        } else {
          avg.seg.baf <- 0.5 - mean(abs(as.numeric(ASCAT.summary[ASCAT.summary$chr == seg[j, ]$chr
          & ASCAT.summary$pos >= seg[j, ]$startpos
          & ASCAT.summary$pos <= seg[j, ]$endpos
          & ASCAT.summary$BAF != 2, ]$BAF) - 0.5))

          ASCAT.summary[ASCAT.summary$chr == seg[j, ]$chr
          & ASCAT.summary$pos >= seg[j, ]$startpos
          & ASCAT.summary$pos <= seg[j, ]$endpos, ]$BAFsegmented <- avg.seg.baf
        }
      }

      # We now need to calculate the average logr of any snps that are present
      if (nrow(ASCAT.summary[ASCAT.summary$chr == seg[j, ]$chr
      & ASCAT.summary$pos >= seg[j, ]$startpos
      & ASCAT.summary$pos <= seg[j, ]$endpos, ]) > 0) {
        avg.seg.logr <- mean(ASCAT.summary[ASCAT.summary$chr == seg[j, ]$chr
        & ASCAT.summary$pos >= seg[j, ]$startpos
        & ASCAT.summary$pos <= seg[j, ]$endpos, ]$logR)

        ASCAT.summary[ASCAT.summary$chr == seg[j, ]$chr
        & ASCAT.summary$pos >= seg[j, ]$startpos
        & ASCAT.summary$pos <= seg[j, ]$endpos, ]$logRsegmented <- avg.seg.logr
      }
    }
    ASCAT.summary
  }

  num.post.parallel.samps <- length(unique(all.samp.SNP.df$sample))

  if (num.post.parallel.samps != num.pre.parallel.samps) {
    message("parallelisation is not working properly, some parallel workers are not returning jobs correctly....")
    message("Either investigate why some cores are not returning results, change to a different node or run in series.")
    message("Running in series now...")
  }

  message("Assigning to list in expected format:")
  all.samp.SNP.list <- list()
  for (samp in sample.vec) {
    message(samp)
    samp.w.segmented.info.SNP.df <- all.samp.SNP.df[all.samp.SNP.df$sample == samp, ]
    print(head(samp.w.segmented.info.SNP.df))
    samp.w.segmented.info.SNP.df <- samp.w.segmented.info.SNP.df[, c("chr", "pos", "BAF", "BAFsegmented", "logR", "logRsegmented")]
    all.samp.SNP.list[[samp]] <- samp.w.segmented.info.SNP.df
  }

  return(all.samp.SNP.list)
}
