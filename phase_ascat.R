# Apply phasing to an existing ascat.bc object and recalculate all necessary
# values within the object
phase_ascat <- function(ascat_bc, phased_bafs, segs, phased_baf_info) {

  message("*** phased_bafs ***")
  str(phased_bafs)

  message("*** phased_baf_info ***")
  str(phased_baf_info)

  #So let's get the phased BAF:
  samp.pha.baf.df <- phased_bafs

  #So we're going to need the indices of the BAF snps...
  samp.pha.baf.df$snp_pos <- phased_baf_info$snp_pos
  samp.pha.baf.df$snp_idx <- phased_baf_info$snp_idx

  #So now we need the segments that refer to this sample (as these are joint segments
  #these should all be same for each sample regardless.)
  samp.seg.df <- segs

  #Let's do the LogR
  tum.logr.vec <- ascat_bc[["Tumor_LogR"]][,1]
  if (any(is.na(tum.logr.vec))) {
    message("There are NA values in the ascat.bc Tumor_LogR field which need to be removed. Revisit generation of this input. For changing the NAs to 0.")
    na.idx.vec <- which(is.na(tum.logr.vec))
    tum.logr.vec[na.idx.vec] <- 0
  }

  tum.logr.segmented.vec <- ascat_bc[["Tumor_LogR_segmented"]][,1]
  if(any(is.na(tum.logr.segmented.vec))) {
    message("There are NA values in the ascat.bc Tumor_LogR_segmented field which need to be removed. Revisit generation of this input. For changing the NAs to 0.")
    na.segmented.idx.vec <- which(is.na(tum.logr.segmented.vec))
    tum.logr.segmented.vec[na.segmented.idx.vec] <- 0
  }

  for (i in 1:nrow(samp.seg.df)) {
    tum.logr.segmented.vec[samp.seg.df[i,]$start_snp_idx:samp.seg.df[i,]$end_snp_idx] <- mean(tum.logr.vec[samp.seg.df[i,]$start_snp_idx:samp.seg.df[i,]$end_snp_idx])
  }

  # Let's do the BAF recalculation
  tum.baf.segmented.vec <- ascat_bc[["Tumor_BAF_segmented"]][[1]][,1]
  names(tum.baf.segmented.vec) <- rownames(ascat_bc[["Tumor_BAF_segmented"]][[1]])

  message("********************************* Disable BAF segmentation update ***********")
  #baf_segmented <- ascat_set_bafsegmented(tum.baf.segmented.vec, samp.seg.df)
  #ascat_set_bafsegmented(ascat_bc, samp.pha.baf.df, samp.seg.df)
  baf_segmented <- tum.baf.segmented.vec

  if(any(is.na(tum.logr.segmented.vec))) {
    message("There are NA values in the tum.logr.segmented.vec after applying phasing which need to be removed. Revisit generation of this input. For changing the NAs to 0.")
    na.segmented.idx.vec <- which(is.na(tum.logr.segmented.vec))
    tum.logr.segmented.vec[na.segmented.idx.vec] <- 0
  }

  #First let's do the free runs:
  message("Setting LogR segmented")
  ascat_bc[["Tumor_LogR_segmented"]][,1] <- tum.logr.segmented.vec
  #message("tum.baf.segmented.vec length =", length(tum.baf.segmented.vec))
  message("ascat.bc[[Tumor_BAF_segmented]][[1]] length =", nrow(ascat_bc[["Tumor_BAF_segmented"]][[1]]))
  #ascat_bc[["Tumor_BAF_segmented"]][[1]]
  #ascat_bc[["Tumor_BAF_segmented"]][[1]] <- matrix(baf_segmented, nrow=1, dimnames(list(as.character(seq_along(baf_segmented)))))

  ascat_phased <- ascat_bc

  #End of section generating the spoofed ascat.bc object that utilises the phasing information. This should be made into a function.

  return(ascat_phased)

}

ascat_set_bafsegmented <- function(ascat, phased_snps, segs) {
  str(ascat)
  str(phased_snps)
  str(segs)

  for (i in seq_len(nrow(segs))) { 
  }

  return(ascat)
}


#  message("Number of NAs: ", sum(is.na(ascat$Tumor_BAF$BAF)))
#  message("Number of non-NAs: ", sum(!is.na(ascat$Tumor_BAF$BAF)))
#  str(ascat)
#
#  ascat_baf <- ascat[["Tumor_BAF_segmented"]][[1]]
#  ascat_baf_idx <- as.numeric(rownames(ascat_baf))
#  #print(ascat_baf_idx)
#
#  baf_segmented <- numeric(nrow(snps))
#
#  baf_idx <- as.integer(baf_segs)
#
#  snps$phased_baf <- ((snps$BAF - 0.5) * (-1 * snps$maj_is_cpnB_logic_vec)) + 0.5
#
#  str(snps)
#  str(segs)
#
#  for (i in seq_len(nrow(segs))) { 
#
#    seg <- segs[i,]
#    snps_idx_for_seg <- snps$snp_idx >= seg$start_snp_idx & snps$snp_idx <= seg$end_snp_idx
#    snps_seg <- snps[snps_idx_for_seg,]
#
#
#    mean_baf_seg <- mean(snps_seg$phased_baf, na.rm = TRUE)
#    mean_baf_seg2 <- mean(snps_seg$BAF, na.rm = TRUE)
#
#    baf_segmented[snps_idx_for_seg] <- mean_baf_seg
#
#    message("segidx: ", i, "\tmean = ", mean_baf_seg)
#    message("segidx: ", i, "\tmean notfphased = ", mean_baf_seg2)

    #baf_seg <- 

    #baf_seg <- mean(

    #snp_idx <- as.integer(names(bafs))

    #baf_seg <- baf[snp_idx >= seg$start_snp_idx & snp_idx <= seg$end_snp_idx]


#    #So first get the indices of the BAF indices within the segment
#    #seg.hz.snp.df <- gl.hz.snp.df[gl.hz.snp.df$snp_idx <= seg$end_snp_idx & gl.hz.snp.df$snp_idx >= seg$start_snp_idx,]
#    #if (nrow(seg.hz.snp.df) == 0) {
#    #  next
#    #}
#
#    #Get indices of hz SNPs that are present ... so we need some form of mapping here to the smaller
#    #These data contain all samples from this patient:
#    seg.min.hz.snp.idx <- seg.hz.snp.df[1,]$hz_snp_idx
#    seg.max.hz.snp.idx <- seg.hz.snp.df[nrow(seg.hz.snp.df),]$hz_snp_idx
#
##    #Let's get the germline baf data specfically so that we can perform a crude correction
##    gl.col.name <- grep("BAF", colnames(seg.hz.snp.df))
##    gl.seg.hz.snp.df <- seg.hz.snp.df[, c(1, 2, gl.col.name)]
##    colnames(gl.seg.hz.snp.df)[3] <- "BAF"
##
##    #Sample specific!
##    #So let's get the SNPs that are within the phased baf df
##    seg.ph.hz.snp.df <- samp.pha.baf.df[samp.pha.baf.df$snp_idx <= seg$end_snp_idx & samp.pha.baf.df$snp_idx >= seg$start_snp_idx,]
##
##    if(nrow(seg.ph.hz.snp.df) != nrow(gl.seg.hz.snp.df)) {
##      message("Number of phased SNPs is different to the actual number hz snps present in seg, perhaps to few to phase by MPhase?")
##      message("nrow(seg.ph.hz.snp.df) = ", nrow(seg.ph.hz.snp.df), " nrow(gl.seg.hz.snp.df) = ", nrow(gl.seg.hz.snp.df))
##      next
##    }
##
##    #GL specific
##    gl.seg.hz.snp.df <- cbind(gl.seg.hz.snp.df, seg.ph.hz.snp.df[, c("maj_is_cpnA_logic_vec", "maj_is_cpnB_logic_vec", "snp_pos", "snp_idx")])
##
##    if (nrow(seg.hz.snp.df) == 0 | nrow(seg.ph.hz.snp.df) == 0) {
##      message("No heterozygous SNPs present within segment.")
##    }
##
##    if (nrow(seg.hz.snp.df) != nrow(seg.ph.hz.snp.df)) {
##      message("Differing numbers of hz snps in the original hz snp df and the phased. Not sure why.")
##    }
#
#    #Ok so...
#    #Let's get the phasing and calculate difference for the sample and the germline
#    #We're assuming that these are in the same order at this point:
#    #germline segmentation value:
#    #Get groups of indices...
#    haplotype.A.snp.pos.vec <- seg.ph.hz.snp.df[seg.ph.hz.snp.df$maj_is_cpnA_logic_vec == 1,]$snp_pos
#    m.hap.A.BAF <- suppressWarnings(mean(seg.ph.hz.snp.df[seg.ph.hz.snp.df$snp_pos %in% haplotype.A.snp.pos.vec,]$BAF, na.rm=TRUE))
#    haplotype.B.snp.pos.vec <- seg.ph.hz.snp.df[seg.ph.hz.snp.df$maj_is_cpnA_logic_vec == 0,]$snp_pos
#    m.hap.B.BAF <- suppressWarnings(mean(seg.ph.hz.snp.df[seg.ph.hz.snp.df$snp_pos %in% haplotype.B.snp.pos.vec,]$BAF, na.rm=TRUE))
#
##    #Now let's do the same for the germline data.. this actually isn't necessary..
##    gl.m.hap.A.BAF <- suppressWarnings(mean(gl.seg.hz.snp.df[gl.seg.hz.snp.df$snp_pos %in% haplotype.A.snp.pos.vec,]$BAF, na.rm=TRUE)) - 0.5
##    gl.m.hap.B.BAF <- suppressWarnings(mean(gl.seg.hz.snp.df[gl.seg.hz.snp.df$snp_pos %in% haplotype.B.snp.pos.vec,]$BAF, na.rm=TRUE)) - 0.5
#
#    #So we actually
#    #So we'll actually work out the gl correction factor here?
#
#    #Let's get the equivalent germline snps 
#    #Default should be to make the BAF segmented 0.5 if we don't have any info w/regards to hz snps
#    seg.md.ph.BAF.vec <- c()
#    # We could attempt to add in a germline correction here...
#    if(! is.na(m.hap.A.BAF) && m.hap.A.BAF > 0.5) {
#      #seg.md.ph.BAF.vec <- c(seg.md.ph.BAF.vec, abs(seg.ph.hz.snp.df[seg.ph.hz.snp.df$snp_pos %in% haplotype.A.snp.pos.vec,]$BAF - 1) - gl.m.hap.A.BAF)
#      seg.md.ph.BAF.vec <- c(seg.md.ph.BAF.vec, abs(seg.ph.hz.snp.df[seg.ph.hz.snp.df$snp_pos %in% haplotype.A.snp.pos.vec,]$BAF - 1))
#    } else if (! is.na(m.hap.A.BAF) && m.hap.A.BAF < 0.5) {
#      #seg.md.ph.BAF.vec <- c(seg.md.ph.BAF.vec, seg.ph.hz.snp.df[seg.ph.hz.snp.df$snp_pos %in% haplotype.A.snp.pos.vec,]$BAF - gl.m.hap.A.BAF)
#      seg.md.ph.BAF.vec <- c(seg.md.ph.BAF.vec, seg.ph.hz.snp.df[seg.ph.hz.snp.df$snp_pos %in% haplotype.A.snp.pos.vec,]$BAF)
#    }
#
#    #We could attempt to add in a germline correction here...
#    if(! is.na(m.hap.B.BAF) && m.hap.B.BAF > 0.5) {
#      #seg.md.ph.BAF.vec <- c(seg.md.ph.BAF.vec, abs(seg.ph.hz.snp.df[seg.ph.hz.snp.df$snp_pos %in% haplotype.B.snp.pos.vec,]$BAF - 1) - gl.m.hap.B.BAF)
#      seg.md.ph.BAF.vec <- c(seg.md.ph.BAF.vec, abs(seg.ph.hz.snp.df[seg.ph.hz.snp.df$snp_pos %in% haplotype.B.snp.pos.vec,]$BAF - 1))
#    } else if (! is.na(m.hap.B.BAF) && m.hap.B.BAF < 0.5) {
#      #seg.md.ph.BAF.vec <- c(seg.md.ph.BAF.vec, seg.ph.hz.snp.df[seg.ph.hz.snp.df$snp_pos %in% haplotype.B.snp.pos.vec,]$BAF - gl.m.hap.B.BAF)
#      seg.md.ph.BAF.vec <- c(seg.md.ph.BAF.vec, seg.ph.hz.snp.df[seg.ph.hz.snp.df$snp_pos %in% haplotype.B.snp.pos.vec,]$BAF)
#    }
#
#    #Format
#    if (length(seg.md.ph.BAF.vec) > 0) {
#      seg.BAF.segmented.val <- rep(mean(seg.md.ph.BAF.vec, na.rm=TRUE), length(seg.min.hz.snp.idx:seg.max.hz.snp.idx))
#    } else {
#      seg.BAF.segmented.val <- rep(0.5, length(seg.min.hz.snp.idx:seg.max.hz.snp.idx)) 
#    }
#
#    #Now let's update the segmented baf values with the correct indices in the vector
#    #So to do this we need to check the actual
#    bafs[seg.min.hz.snp.idx:seg.max.hz.snp.idx] <- seg.BAF.segmented.val
#    #rownames(tum.baf.segmented.vec)
#
#    #Now we could edit this to compare to the germline and work out whether we should actually call this as separate....
#
#  }
#
#  str(baf_segmented)
#
#  baf_segmented
#
#}
