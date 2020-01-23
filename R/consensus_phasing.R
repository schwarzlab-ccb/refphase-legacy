getCombinedBAFDipTestPhasing <- function(asas.lr.baf.lists,
                                         seg.name.vec,
                                         ploidy.est.vec,
                                         purity.est.vec,
                                         gamma.val = 1,
                                         min.num.hz.snp.phasing = 5,
                                         min.prop.snp.sep.threshold = 0.75,
                                         dip.test.p.val.threshold = 0.01,
                                         comb.wilcox.snp.sep.p.val.threshold = 0.01,
                                         allow.dip.test.without.ASCAT.to.phase = FALSE) {
  sample.vec <- names(asas.lr.baf.lists)
  sample.to.phase.seg.list <- list()

  for (seg.name in seg.name.vec) {
    first.ss.lr.baf.df <- asas.lr.baf.lists[[sample.vec[1]]][[seg.name]]

    if (nrow(first.ss.lr.baf.df) < min.num.hz.snp.phasing) {
      next
    }

    # For each sample how good is it at phasing all samples for this segment.
    all.samp.seg.phasing.ability.list <- getAssessmentOfEachSamplesAbilityToPhaseAllOthersForSeg(
      asas.lr.baf.lists = asas.lr.baf.lists,
      sample.vec = sample.vec,
      seg.name = seg.name,
      ploidy.est.vec = ploidy.est.vec,
      purity.est.vec = purity.est.vec,
      gamma.val = gamma.val,
      min.num.hz.snp.phasing = min.num.hz.snp.phasing
    )
    sample.to.phase.seg.list[[seg.name]]$all_samp_phasing_comparison_list <- all.samp.seg.phasing.ability.list
    a.s.sum.phas.df <- data.frame(
      diptest_for_AI_pval = sapply(all.samp.seg.phasing.ability.list, function(x) x$diptest_for_AI_pval),
      ASCAT_called_as_AI = sapply(all.samp.seg.phasing.ability.list, function(x) x$ASCAT_called_as_AI),
      stringsAsFactors = FALSE
    )
    rownames(a.s.sum.phas.df) <- sample.vec
    # Reorder by which has best separation according to th eoriginal diptest, then subset to choose only samples that pass the threshold for having allelic imbalance.
    a.s.sum.phas.df <- a.s.sum.phas.df[with(a.s.sum.phas.df, order(diptest_for_AI_pval)), ]
    # Now select the samples that we will use to call for phasing.
    if (allow.dip.test.without.ASCAT.to.phase == TRUE) {
      poss.phas.s.df <- a.s.sum.phas.df[a.s.sum.phas.df$ASCAT_called_as_AI == TRUE | a.s.sum.phas.df$diptest_for_AI_pval < dip.test.p.val.threshold, , drop = FALSE]
    } else {
      poss.phas.s.df <- a.s.sum.phas.df[a.s.sum.phas.df$ASCAT_called_as_AI == TRUE, , drop = FALSE]
    }
    # Now assess the samples that can actualy be phased...
    # initialise the seg mixture phasing to NULL so that, as a default, it won't be used:
    seg.mixture.phasing.list <- NULL
    if (nrow(poss.phas.s.df) > 0) {
      print(poss.phas.s.df)
      if (nrow(poss.phas.s.df) == 1) {
        # If only one sample passes then choose this sample as the reference
        sample.to.phase.seg <- rownames(poss.phas.s.df)[1]
      } else if (nrow(poss.phas.s.df) > 1) {
        # If there are multiple alance instances called choose the top one and check mixture phasing
        sample.to.phase.seg <- rownames(poss.phas.s.df)[1]
        poss.phas.sample.vec <- rownames(poss.phas.s.df)
        print(poss.phas.sample.vec)
        all.poss.phas.samp.seg.phasing.ability.list <- getAssessmentOfEachSamplesAbilityToPhaseAllOthersForSeg(
          asas.lr.baf.lists = asas.lr.baf.lists,
          sample.vec = poss.phas.sample.vec,
          seg.name = seg.name,
          ploidy.est.vec = ploidy.est.vec,
          purity.est.vec = purity.est.vec,
          gamma.val = gamma.val,
          min.num.hz.snp.phasing = min.num.hz.snp.phasing
        )

        phas.s.df <- data.frame(
          diptest_for_AI_pval = sapply(all.poss.phas.samp.seg.phasing.ability.list, function(x) x$diptest_for_AI_pval),
          ASCAT_called_as_AI = sapply(all.poss.phas.samp.seg.phasing.ability.list, function(x) x$ASCAT_called_as_AI),
          mean_prop_SNPs_sep_across_samples = sapply(all.poss.phas.samp.seg.phasing.ability.list, function(x) x$mean_prop_SNPs_sep_across_samples),
          comb_pval_SNP_sep_across_samps = sapply(all.poss.phas.samp.seg.phasing.ability.list, function(x) x$comb_wilcox_sep_pval_across_samps),
          stringsAsFactors = FALSE
        )

        seg.mixture.phasing.list <- getMixturePhasingResultsForSegment(
          sample.to.phase.seg = sample.to.phase.seg,
          sample.vec = poss.phas.sample.vec,
          poss.phas.s.df = phas.s.df,
          asas.lr.baf.lists = asas.lr.baf.lists,
          all.samp.seg.phasing.ability.list = all.poss.phas.samp.seg.phasing.ability.list,
          seg.name = seg.name,
          ploidy.est.vec = ploidy.est.vec,
          purity.est.vec = purity.est.vec,
          gamma.val = gamma.val,
          comb.wilcox.snp.sep.p.val.threshold = comb.wilcox.snp.sep.p.val.threshold,
          min.prop.snp.sep.threshold = min.prop.snp.sep.threshold
        )
      } else {
        stop("There's an error in phasing the phasing logic")
      }
    } else {
      sample.to.phase.seg <- "none_no_allelic_imbalance_called"
    }
    # Update output with the choice of phasing:
    if (is.null(seg.mixture.phasing.list)) { # mixture phasing was not required

      if (grepl("none", sample.to.phase.seg) == TRUE) {
        # As no phasing was determined to be possible we will simply take
        # alternate indices to be cpnA or cpnB as being in the upper band in
        # the in the reference.
        all.phaseable.positions.idx.vec <- c()
        maj.is.cpnA.idx.vec <- all.samp.seg.phasing.ability.list[[1]]$maj_is_cpnA_idx_vec
        maj.is.cpnB.idx.vec <- all.samp.seg.phasing.ability.list[[1]]$maj_is_cpnB_idx_vec
        if (length(maj.is.cpnA.idx.vec[!is.na(maj.is.cpnA.idx.vec)]) > 0) {
          all.phaseable.positions.idx.vec <- c(all.phaseable.positions.idx.vec, maj.is.cpnA.idx.vec)
        }
        if (length(maj.is.cpnB.idx.vec[!is.na(maj.is.cpnB.idx.vec)]) > 0) {
          all.phaseable.positions.idx.vec <- c(all.phaseable.positions.idx.vec, maj.is.cpnB.idx.vec)
        }
        sorted.all.phas.idx.vec <- sort(all.phaseable.positions.idx.vec)
        sample.to.phase.seg.list[[seg.name]]$samp_chosen_to_phase_seg <- sample.to.phase.seg

        if (length(sorted.all.phas.idx.vec) > 0) {
          sample.to.phase.seg.list[[seg.name]]$maj_is_cpnA_idx_vec <- sorted.all.phas.idx.vec[seq(1, length(all.phaseable.positions.idx.vec), 2)]
        } else {
          sample.to.phase.seg.list[[seg.name]]$maj_is_cpnA_idx_vec <- NA
        }
        if (length(sorted.all.phas.idx.vec) > 1) {
          sample.to.phase.seg.list[[seg.name]]$maj_is_cpnB_idx_vec <- sorted.all.phas.idx.vec[seq(2, length(all.phaseable.positions.idx.vec), 2)]
        } else {
          sample.to.phase.seg.list[[seg.name]]$maj_is_cpnB_idx_vec <- NA
        }
      } else {
        # As there is a phasing possible Take the appropriate sample selected
        sample.to.phase.seg.list[[seg.name]]$samp_chosen_to_phase_seg <- sample.to.phase.seg
        sample.to.phase.seg.list[[seg.name]]$maj_is_cpnA_idx_vec <- all.samp.seg.phasing.ability.list[[sample.to.phase.seg]]$maj_is_cpnA_idx_vec
        sample.to.phase.seg.list[[seg.name]]$maj_is_cpnB_idx_vec <- all.samp.seg.phasing.ability.list[[sample.to.phase.seg]]$maj_is_cpnB_idx_vec
      }
    } else {
      # mixture phasing has worked/been performed and is required, take the
      # appropriate results from the mixture phasing.
      sample.to.phase.seg.list[[seg.name]]$samp_chosen_to_phase_seg <- "mixture"
      sample.to.phase.seg.list[[seg.name]]$samp_most_significant_by_diptest <- seg.mixture.phasing.list$samp_most_significant_by_diptest
      sample.to.phase.seg.list[[seg.name]]$maj_is_cpnA_idx_vec <- seg.mixture.phasing.list$maj_is_cpnA_idx_vec
      sample.to.phase.seg.list[[seg.name]]$maj_is_cpnB_idx_vec <- seg.mixture.phasing.list$maj_is_cpnB_idx_vec
      sample.to.phase.seg.list[[seg.name]]$all_samp_mixture_phasing_comparison_list <- seg.mixture.phasing.list$all_samp_mixture_phasing_comparison_list
    }
  }

  return(sample.to.phase.seg.list)
}

#' @importFrom foreach foreach %do% %dopar%
getCombinedBAFDipTestPhasingParallel <- function(asas.lr.baf.lists,
                                                 seg.name.vec,
                                                 ploidy.est.vec,
                                                 purity.est.vec,
                                                 gamma.val = 1,
                                                 min.num.hz.snp.phasing = 5,
                                                 min.prop.snp.sep.threshold = 0.75,
                                                 dip.test.p.val.threshold = 0.01,
                                                 comb.wilcox.snp.sep.p.val.threshold = 0.01,
                                                 allow.dip.test.without.ASCAT.to.phase = FALSE) {
  sample.vec <- names(asas.lr.baf.lists)
  sample.to.phase.seg.list <- list()

  MultiResultClass <- function() {
    me <- list(
      samp_chosen_to_phase_seg = NA,
      samp_most_significant_by_diptest = NA,
      maj_is_cpnA_idx_vec = NA,
      maj_is_cpnB_idx_vec = NA,
      all_samp_mixture_phasing_comparison_list = NA
    )
    # Set the name for the class
    class(me) <- append(class(me), "MultiResultClass")
    return(me)
  }

  itx <- iterators::iter(seg.name.vec)
  # sample.to.phase.seg.mrc.list <- foreach(seg.name = itx, .combine = c) %dopar% {
  sample.to.phase.seg.mrc.list <- foreach(seg.name = itx, .combine = c) %do% {
    ret.mrc <- MultiResultClass()
    first.ss.lr.baf.df <- asas.lr.baf.lists[[sample.vec[1]]][[seg.name]]
    if (nrow(first.ss.lr.baf.df) < min.num.hz.snp.phasing) {

      # next;#If phasng not possible for this segment then move onto the next
      ret.mrc <- c()
    } else { # end if( nrow(ref.ss.lr.baf.df) < min.hz.snp.phasing)
      # For each sample how good is it at phasing all samples for this segment.
      all.samp.seg.phasing.ability.list <- getAssessmentOfEachSamplesAbilityToPhaseAllOthersForSeg(
        asas.lr.baf.lists = asas.lr.baf.lists,
        sample.vec = sample.vec,
        seg.name = seg.name,
        ploidy.est.vec = ploidy.est.vec,
        purity.est.vec = purity.est.vec,
        gamma.val = gamma.val,
        min.num.hz.snp.phasing = min.num.hz.snp.phasing
      )
      sample.to.phase.seg.list[[seg.name]]$all_samp_phasing_comparison_list <- all.samp.seg.phasing.ability.list
      a.s.sum.phas.df <- data.frame(
        diptest_for_AI_pval = sapply(all.samp.seg.phasing.ability.list, function(x) x$diptest_for_AI_pval),
        ASCAT_called_as_AI = sapply(all.samp.seg.phasing.ability.list, function(x) x$ASCAT_called_as_AI),
        AI_MAD = sapply(all.samp.seg.phasing.ability.list, function(x) x$AI_MAD),
        stringsAsFactors = FALSE
      )
      rownames(a.s.sum.phas.df) <- sample.vec
      # Reorder by which has best separation according to th eoriginal diptest,
      # then subset to choose only samples that pass the threshold for having
      # allelic imbalance.
      a.s.sum.phas.df <- a.s.sum.phas.df[with(a.s.sum.phas.df, order(diptest_for_AI_pval, -AI_MAD)), ]
      # Now select the samples that we will use to call for phasing.
      if (allow.dip.test.without.ASCAT.to.phase == TRUE) {
        poss.phas.s.df <- a.s.sum.phas.df[a.s.sum.phas.df$ASCAT_called_as_AI == TRUE | a.s.sum.phas.df$diptest_for_AI_pval < dip.test.p.val.threshold, , drop = FALSE]
      } else {
        poss.phas.s.df <- a.s.sum.phas.df[a.s.sum.phas.df$ASCAT_called_as_AI == TRUE, , drop = FALSE]
      }
      # Now assess the samples that can actualy be phased...
      # initialise the seg mixture phasing to NULL so that, as a default, it
      # won't be used:
      seg.mixture.phasing.list <- NULL
      if (nrow(poss.phas.s.df) > 0) {
        if (nrow(poss.phas.s.df) == 1) {
          # If only one sample passes then choose this sample as the reference
          sample.to.phase.seg <- rownames(poss.phas.s.df)[1]
        } else if (nrow(poss.phas.s.df) > 1) {
          # If there are multiple alance instances called choose the top one
          # and check mixture phasing
          sample.to.phase.seg <- rownames(poss.phas.s.df)[1]
          poss.phas.sample.vec <- rownames(poss.phas.s.df)
          all.poss.phas.samp.seg.phasing.ability.list <- getAssessmentOfEachSamplesAbilityToPhaseAllOthersForSeg(
            asas.lr.baf.lists = asas.lr.baf.lists,
            sample.vec = poss.phas.sample.vec,
            seg.name = seg.name,
            ploidy.est.vec = ploidy.est.vec,
            purity.est.vec = purity.est.vec,
            gamma.val = gamma.val,
            min.num.hz.snp.phasing = min.num.hz.snp.phasing
          )

          phas.s.df <- data.frame(
            diptest_for_AI_pval = sapply(all.poss.phas.samp.seg.phasing.ability.list, function(x) x$diptest_for_AI_pval),
            ASCAT_called_as_AI = sapply(all.poss.phas.samp.seg.phasing.ability.list, function(x) x$ASCAT_called_as_AI),
            mean_prop_SNPs_sep_across_samples = sapply(all.poss.phas.samp.seg.phasing.ability.list, function(x) x$mean_prop_SNPs_sep_across_samples),
            comb_pval_SNP_sep_across_samps = sapply(all.poss.phas.samp.seg.phasing.ability.list, function(x) x$comb_wilcox_sep_pval_across_samps),
            stringsAsFactors = FALSE
          )

          seg.mixture.phasing.list <- getMixturePhasingResultsForSegment(
            sample.to.phase.seg = sample.to.phase.seg,
            sample.vec = poss.phas.sample.vec,
            poss.phas.s.df = phas.s.df,
            asas.lr.baf.lists = asas.lr.baf.lists,
            all.samp.seg.phasing.ability.list = all.poss.phas.samp.seg.phasing.ability.list,
            seg.name = seg.name,
            ploidy.est.vec = ploidy.est.vec,
            purity.est.vec = purity.est.vec,
            gamma.val = gamma.val,
            comb.wilcox.snp.sep.p.val.threshold = comb.wilcox.snp.sep.p.val.threshold,
            min.prop.snp.sep.threshold = min.prop.snp.sep.threshold
          )
        } else {
          stop("There's an error in phasing the phasing logic")
        }
      } else {
        sample.to.phase.seg <- "none_no_allelic_imbalance_called"
      }
      # Update output with the choice of phasing:
      if (is.null(seg.mixture.phasing.list)) {
        # mixture phasing was not required

        if (isTRUE(grepl("none", sample.to.phase.seg))) {
          # As no phasing was determined to be possible we will simply take
          # alternate indices to be cpnA or cpnB as being in the upper band in
          # the in the reference.
          all.phaseable.positions.idx.vec <- c()
          maj.is.cpnA.idx.vec <- all.samp.seg.phasing.ability.list[[1]]$maj_is_cpnA_idx_vec
          maj.is.cpnB.idx.vec <- all.samp.seg.phasing.ability.list[[1]]$maj_is_cpnB_idx_vec
          if (length(maj.is.cpnA.idx.vec[!is.na(maj.is.cpnA.idx.vec)]) > 0) {
            all.phaseable.positions.idx.vec <- c(all.phaseable.positions.idx.vec, maj.is.cpnA.idx.vec)
          }
          if (length(maj.is.cpnB.idx.vec[!is.na(maj.is.cpnB.idx.vec)]) > 0) {
            all.phaseable.positions.idx.vec <- c(all.phaseable.positions.idx.vec, maj.is.cpnB.idx.vec)
          }
          sorted.all.phas.idx.vec <- sort(all.phaseable.positions.idx.vec)
          ret.mrc$samp_chosen_to_phase_seg <- sample.to.phase.seg

          if (length(sorted.all.phas.idx.vec) > 0) {
            ret.mrc$maj_is_cpnA_idx_vec <- sorted.all.phas.idx.vec[seq(1, length(all.phaseable.positions.idx.vec), 2)]
          } else {
            ret.mrc$maj_is_cpnA_idx_vec <- NA
          }
          if (length(sorted.all.phas.idx.vec) > 1) {
            ret.mrc$maj_is_cpnB_idx_vec <- sorted.all.phas.idx.vec[seq(2, length(all.phaseable.positions.idx.vec), 2)]
          } else {
            ret.mrc$maj_is_cpnB_idx_vec <- NA
          }
        } else {
          # As there is a phasing possible Take the appropriate sample selected
          ret.mrc$samp_chosen_to_phase_seg <- sample.to.phase.seg
          ret.mrc$maj_is_cpnA_idx_vec <- all.samp.seg.phasing.ability.list[[sample.to.phase.seg]]$maj_is_cpnA_idx_vec
          ret.mrc$maj_is_cpnB_idx_vec <- all.samp.seg.phasing.ability.list[[sample.to.phase.seg]]$maj_is_cpnB_idx_vec
        }
      } else {
        # mixture phasing has worked/been performed and is required, take the appropriate results from the mixture phasing.
        ret.mrc$samp_chosen_to_phase_seg <- "mixture"
        ret.mrc$samp_most_significant_by_diptest <- seg.mixture.phasing.list$samp_most_significant_by_diptest
        ret.mrc$maj_is_cpnA_idx_vec <- seg.mixture.phasing.list$maj_is_cpnA_idx_vec
        ret.mrc$maj_is_cpnB_idx_vec <- seg.mixture.phasing.list$maj_is_cpnB_idx_vec
        ret.mrc$all_samp_mixture_phasing_comparison_list <- seg.mixture.phasing.list$all_samp_mixture_phasing_comparison_list
      }
    }

    returnable.list <- list(ret.mrc)
    names(returnable.list) <- seg.name
    returnable.list
  }

  logdebug("Finished paralell portion of getCombinedBAFDipTestPhasingParallel, combining results into correct format...")
  sample.to.phase.seg.list <- list()
  for (seg.name in names(sample.to.phase.seg.mrc.list)) {
    if (length(sample.to.phase.seg.mrc.list[[seg.name]]) > 0) {
      sample.to.phase.seg.list[[seg.name]]$samp_chosen_to_phase_seg <- sample.to.phase.seg.mrc.list[[seg.name]]$samp_chosen_to_phase_seg
      sample.to.phase.seg.list[[seg.name]]$samp_most_significant_by_diptest <- sample.to.phase.seg.mrc.list[[seg.name]]$samp_most_significant_by_diptest
      sample.to.phase.seg.list[[seg.name]]$maj_is_cpnA_idx_vec <- sample.to.phase.seg.mrc.list[[seg.name]]$maj_is_cpnA_idx_vec
      sample.to.phase.seg.list[[seg.name]]$maj_is_cpnB_idx_vec <- sample.to.phase.seg.mrc.list[[seg.name]]$maj_is_cpnB_idx_vec
      sample.to.phase.seg.list[[seg.name]]$all_samp_mixture_phasing_comparison_list <- sample.to.phase.seg.mrc.list[[seg.name]]$all_samp_mixture_phasing_comparison_list
    }
  }

  return(sample.to.phase.seg.list)
}


assessPhasingAcrossSelectedSamples <- function(seg.name,
                                               sample.comparison.list,
                                               samp,
                                               asas.lr.baf.lists,
                                               samples.to.test.vec,
                                               ref.ss.maj.is.cpnA.idx.vec,
                                               ref.ss.maj.is.cpnB.idx.vec,
                                               ploidy.est.vec,
                                               purity.est.vec,
                                               gamma.val,
                                               use.xphase = FALSE,
                                               ascat.raw.cpn) {
  other.ss.phasing.df <- data.frame(
    separation_wilcox_pval = rep(1, length(samples.to.test.vec)),
    prop_separated = numeric(length(samples.to.test.vec)),
    concordant_w_ref = numeric(length(samples.to.test.vec)),
    stringsAsFactors = FALSE
  )
  rownames(other.ss.phasing.df) <- samples.to.test.vec

  for (other.samp in samples.to.test.vec) {
    oth.ss.lr.baf.df <- asas.lr.baf.lists[[other.samp]][[seg.name]]

    if (!use.xphase) {
      oth.ss.cpnA.cpnB.df <- calculateSNPsRawCopyNumberAB(
        logr.vec = oth.ss.lr.baf.df$logR,
        baf.vec = oth.ss.lr.baf.df$BAF,
        tumor.ploidy = ploidy.est.vec[other.samp],
        tumor.purity = purity.est.vec[other.samp],
        gamma.val = gamma.val
      )
    } else {
      oth.ss.cpnA.cpnB.df <- getPhasedASCATRawCopyNumberAB(
        ref.ss.lr.baf.df = oth.ss.lr.baf.df,
        sample = other.samp,
        sample.ph.hz.ascat.df.list = ascat.raw.cpn
      )
    }

    # old code had some checks here to make sure that these have the same
    # number of rows, but if everything holds they should remain the same
    # anyway...
    other.samp.wilcox.test <- wilcoxTestWithErrorHandling(c(oth.ss.cpnA.cpnB.df[ref.ss.maj.is.cpnA.idx.vec, 1], oth.ss.cpnA.cpnB.df[ref.ss.maj.is.cpnB.idx.vec, 2]),
      c(oth.ss.cpnA.cpnB.df[ref.ss.maj.is.cpnA.idx.vec, 2], oth.ss.cpnA.cpnB.df[ref.ss.maj.is.cpnB.idx.vec, 1]),
      paired = TRUE
    )

    oth.ss.snps.same.as.ref <- length(which(c(oth.ss.cpnA.cpnB.df[ref.ss.maj.is.cpnA.idx.vec, 1], oth.ss.cpnA.cpnB.df[ref.ss.maj.is.cpnB.idx.vec, 2]) -
      c(oth.ss.cpnA.cpnB.df[ref.ss.maj.is.cpnA.idx.vec, 2], oth.ss.cpnA.cpnB.df[ref.ss.maj.is.cpnB.idx.vec, 1]) > 0)) / nrow(oth.ss.cpnA.cpnB.df)
    oth.ss.snps.diff.to.ref <- length(which(c(oth.ss.cpnA.cpnB.df[ref.ss.maj.is.cpnA.idx.vec, 1], oth.ss.cpnA.cpnB.df[ref.ss.maj.is.cpnB.idx.vec, 2]) -
      c(oth.ss.cpnA.cpnB.df[ref.ss.maj.is.cpnA.idx.vec, 2], oth.ss.cpnA.cpnB.df[ref.ss.maj.is.cpnB.idx.vec, 1]) < 0)) / nrow(oth.ss.cpnA.cpnB.df)
    # Updating the other sample segment phasing df
    other.ss.phasing.df[other.samp, ]$separation_wilcox_pval <- other.samp.wilcox.test$p.val
    other.ss.phasing.df[other.samp, ]$prop_separated <- max(oth.ss.snps.same.as.ref, oth.ss.snps.diff.to.ref)
    other.ss.phasing.df[other.samp, ]$concordant_w_ref <- ifelse(oth.ss.snps.same.as.ref >= oth.ss.snps.diff.to.ref, 1, -1)
  }

  sample.comparison.list[[samp]]$comb_wilcox_sep_pval_across_samps <- pchisq((sum(log(other.ss.phasing.df$separation_wilcox_pval)) * -2),
    df = length(other.ss.phasing.df$separation_wilcox_pval) * 2,
    lower.tail = FALSE
  )
  sample.comparison.list[[samp]]$oth_samp_wilcox_sep_pvals_vec <- other.ss.phasing.df$separation_wilcox_pval
  sample.comparison.list[[samp]]$oth_samp_SNPs_sep_prop_vec <- other.ss.phasing.df$prop_separated
  sample.comparison.list[[samp]]$maj_is_cpnA_idx_vec <- ref.ss.maj.is.cpnA.idx.vec
  sample.comparison.list[[samp]]$maj_is_cpnB_idx_vec <- ref.ss.maj.is.cpnB.idx.vec
  sample.comparison.list[[samp]]$mean_prop_SNPs_sep_across_samples <- mean(other.ss.phasing.df$prop_separated)
  sample.comparison.list[[samp]]$oth_samp_cpnA_cpnB_match_ref_vec <- other.ss.phasing.df$concordant_w_ref

  return(sample.comparison.list)
}


# For this segment go through and calculate how good each sample is at phasing
# all the other samples
# This will give us a list that for each sample contains another list comparing
# its ability to phase all others via a dip test on the major cpn vs all other
# samples and the wilcox tests on all others using
# the "phasing" given by the reference sample.
getAssessmentOfEachSamplesAbilityToPhaseAllOthersForSeg <- function(asas.lr.baf.lists,
                                                                    sample.vec,
                                                                    seg.name,
                                                                    ploidy.est.vec,
                                                                    purity.est.vec,
                                                                    gamma.val = 1,
                                                                    min.num.hz.snp.phasing = 5,
                                                                    use.xphase = FALSE,
                                                                    ascat.raw.cpn) {
  # other sample segment phasing df, used when working out the optimum sample
  # in all regions, for each segment the same dataframe is reused to save time
  # allocating memory...

  num.other.samps <- length(sample.vec) - 1
  other.ss.phasing.df <- data.frame(
    separation_wilcox_pval = rep(1, num.other.samps),
    prop_separated = numeric(num.other.samps),
    concordant_w_ref = numeric(num.other.samps),
    stringsAsFactors = FALSE
  )
  sample.comparison.list <- list()

  for (samp in sample.vec) {
    ref.ss.lr.baf.df <- asas.lr.baf.lists[[samp]][[seg.name]]

    if (nrow(ref.ss.lr.baf.df) < min.num.hz.snp.phasing) {
      # Might actally want this to be a break as they should all have too few
      # snps as it is joint segmentation with the same SNPs...
      next
    }

    # renaming the dataframe rows for the other samples, now we are trying a
    # new reference sample
    rownames(other.ss.phasing.df) <- sample.vec[ !sample.vec %in% samp ]

    if (!use.xphase) {
      ref.ss.cpnA.cpnB.df <- calculateSNPsRawCopyNumberAB(
        logr.vec = ref.ss.lr.baf.df$logR,
        baf.vec = ref.ss.lr.baf.df$BAF,
        tumor.ploidy = ploidy.est.vec[samp],
        tumor.purity = purity.est.vec[samp],
        gamma.val = gamma.val
      )
    } else {
      ref.ss.cpnA.cpnB.df <- getPhasedASCATRawCopyNumberAB(
        ref.ss.lr.baf.df = ref.ss.lr.baf.df,
        sample = samp,
        sample.ph.hz.ascat.df.list = ascat.raw.cpn
      )
    }

    ref.ss.max.cpn.vec <- apply(ref.ss.cpnA.cpnB.df, 1, max)
    ref.ss.is.max.cpn.cpnA.cpnB.logical.df <- ifelse(ref.ss.cpnA.cpnB.df >= ref.ss.max.cpn.vec, 1, 0)
    ref.ss.maj.is.cpnA.idx.vec <- which(ref.ss.is.max.cpn.cpnA.cpnB.logical.df[, 1] == 1)
    ref.ss.maj.is.cpnB.idx.vec <- which(ref.ss.is.max.cpn.cpnA.cpnB.logical.df[, 2] == 1)

    sample.comparison.list[[samp]]$diptest_for_AI_pval <- as.numeric(diptest::dip.test(ref.ss.lr.baf.df$BAF)$p.val)
    sample.comparison.list[[samp]]$AI_MAD <- mean(abs(ref.ss.lr.baf.df$BAF - 0.5))

    # Changed to 2 due to idxing but... as only uysed for plotting ignoring.
    sample.comparison.list[[samp]]$ASCAT_called_as_AI <- ifelse(median(ref.ss.lr.baf.df$BAFsegmented, na.rm = TRUE) %in% 0.5, FALSE, TRUE)
    samples.to.test.vec <- sample.vec[ !sample.vec %in% samp]
    sample.comparison.list <- assessPhasingAcrossSelectedSamples(
      seg.name = seg.name,
      sample.comparison.list = sample.comparison.list,
      samp = samp,
      asas.lr.baf.lists = asas.lr.baf.lists,
      samples.to.test.vec = samples.to.test.vec,
      ref.ss.maj.is.cpnA.idx.vec = ref.ss.maj.is.cpnA.idx.vec,
      ref.ss.maj.is.cpnB.idx.vec = ref.ss.maj.is.cpnB.idx.vec,
      ploidy.est.vec = ploidy.est.vec,
      purity.est.vec = purity.est.vec,
      gamma.val = gamma.val
    )
  }

  return(sample.comparison.list)
}


getMixturePhasingResultsForSegment <- function(sample.to.phase.seg,
                                               sample.vec,
                                               poss.phas.s.df,
                                               asas.lr.baf.lists,
                                               all.samp.seg.phasing.ability.list,
                                               seg.name,
                                               ploidy.est.vec,
                                               purity.est.vec,
                                               gamma.val,
                                               comb.wilcox.snp.sep.p.val.threshold,
                                               min.prop.snp.sep.threshold,
                                               use.xphase = FALSE,
                                               ascat.raw.cpn) {

  # This is currently dead code... as of 20190818
  return(NULL)
}
