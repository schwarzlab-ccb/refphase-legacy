ascat_set_bafsegmented <- function(ascat, phased_snps, segs) {
  warning("FIXME: Using a simple re-implementation of code that calculates BAF segmented values for ascat.bc")

  # Reset all BAF_segmented values to 0.5 so we know which ones aren't updated
  ascat$Tumor_BAF_segmented[[1]][, 1] <- 0.5

  # Chromosome + position of each element in Tumor_BAF_segmented
  ascat_bafseg_pos <- ascat$SNPpos[as.numeric(rownames(ascat$Tumor_BAF_segmented[[1]])), ]

  for (i in seq_len(nrow(segs))) {
    seg <- segs[i, ]

    idx_phased <- phased_snps$chrom == seg$Chromosome & phased_snps$pos >= seg$Start_Position & phased_snps$pos <= seg$End_Position
    phased_snps_seg <- phased_snps[idx_phased, ]
    pre_phased_bafs <- unique(phased_snps_seg$ph_median_segmented_BAF)
    if (length(pre_phased_bafs) == 0) {
      pre_phased_bafs <- 0.5
    }

    idx_ascat <- ascat_bafseg_pos$chrs == seg$Chromosome & ascat_bafseg_pos$pos >= seg$Start_Position & ascat_bafseg_pos$pos <= seg$End_Position

    if (any(idx_ascat)) {
      ascat$Tumor_BAF_segmented[[1]][idx_ascat, 1] <- pre_phased_bafs
    }
  }

  if (any(is.na(ascat$Tumor_BAF_segmented[[1]][, 1]))) {
    logdebug(paste0("Number of NA bafsegmented values: ", sum(is.na(ascat$Tumor_BAF_segmented[[1]][, 1])), " (they will be removed)"))
    # Remove NA elements
    bafseg <- ascat$Tumor_BAF_segmented[[1]][, 1, drop = FALSE]
    ascat$Tumor_BAF_segmented[[1]] <- ascat$Tumor_BAF_segmented[[1]][!is.na(bafseg), 1, drop = FALSE]
  }

  ascat
}


# Apply phasing to an existing ascat.bc object and recalculate all necessary
# values within the object
phase_ascat <- function(ascat_bc, phased_bafs, segs, phased_baf_info) {
  # So let's get the phased BAF:
  samp.pha.baf.df <- phased_bafs

  # So we're going to need the indices of the BAF snps...
  samp.pha.baf.df$snp_pos <- phased_baf_info$snp_pos
  samp.pha.baf.df$snp_idx <- phased_baf_info$snp_idx

  # So now we need the segments that refer to this sample (as these are joint segments
  # these should all be same for each sample regardless.)
  samp.seg.df <- segs

  tum.logr.vec <- ascat_bc[["Tumor_LogR"]][, 1]
  if (any(is.na(tum.logr.vec))) {
    logwarn("There are NA values in the ascat.bc Tumor_LogR field which need to be removed. Revisit generation of this input. For changing the NAs to 0.")
    na.idx.vec <- which(is.na(tum.logr.vec))
    tum.logr.vec[na.idx.vec] <- 0
  }

  tum.logr.segmented.vec <- ascat_bc[["Tumor_LogR_segmented"]][, 1]
  if (any(is.na(tum.logr.segmented.vec))) {
    logwarn("There are NA values in the ascat.bc Tumor_LogR_segmented field which need to be removed. Revisit generation of this input. For changing the NAs to 0.")
    na.segmented.idx.vec <- which(is.na(tum.logr.segmented.vec))
    tum.logr.segmented.vec[na.segmented.idx.vec] <- 0
  }
  for (i in 1:nrow(samp.seg.df)) {
    tum.logr.segmented.vec[samp.seg.df[i, ]$start_snp_idx:samp.seg.df[i, ]$end_snp_idx] <- mean(tum.logr.vec[samp.seg.df[i, ]$start_snp_idx:samp.seg.df[i, ]$end_snp_idx])
  }

  ascat_bc <- ascat_set_bafsegmented(ascat_bc, samp.pha.baf.df, samp.seg.df)

  if (any(is.na(tum.logr.segmented.vec))) {
    logwarn("There are NA values in the tum.logr.segmented.vec after applying phasing which need to be removed. Revisit generation of this input. For changing the NAs to 0.")
    na.segmented.idx.vec <- which(is.na(tum.logr.segmented.vec))
    tum.logr.segmented.vec[na.segmented.idx.vec] <- 0
  }

  ascat_bc[["Tumor_LogR_segmented"]][, 1] <- tum.logr.segmented.vec

  ascat_bc
}
