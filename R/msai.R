# This is a new version
getAllSegsAIorMSAIDetectionFromBAFListLTNPD <- function(seg.name.vec, sample.vec, asas.SNPS.lr.and.or.baf.lists, all.seg.phasing.info.list, MSAI.BAF.LTPND.pval.thresh) {
  all.seg.MSAI.by.BAF.ltnpd.present.list <- list()
  for (seg.name in seg.name.vec) {
    all.seg.MSAI.by.BAF.ltnpd.present.list[[seg.name]] <- list()
    all.seg.MSAI.by.BAF.ltnpd.present.list[[seg.name]]$seg_name <- seg.name
    #message(seg.name)
    if (seg.name %in% names(all.seg.phasing.info.list)) {
      ref.samp <- all.seg.phasing.info.list[[seg.name]]$samp_chosen_to_phase_seg
      # Determine whether we can detect MSAI and fill in the info fields, resent ref.samp if mixture phasing was used
      if (grepl("none", ref.samp) == TRUE) {
        all.seg.MSAI.by.BAF.ltnpd.present.list[[seg.name]]$MSAI_detection_feasible <- FALSE
        all.seg.MSAI.by.BAF.ltnpd.present.list[[seg.name]]$MSAI_detection_info <- "no_reference_although_phasing_attempted"
        next
      } else if (ref.samp == "mixture") {
        all.seg.MSAI.by.BAF.ltnpd.present.list[[seg.name]]$MSAI_detection_feasible <- TRUE
        all.seg.MSAI.by.BAF.ltnpd.present.list[[seg.name]]$MSAI_detection_info <- "mixture_phasing_using_top_ranked_samp_by_diptest" # Note might be an issue here, is it only
        # samples that made it through
        ref.samp <- all.seg.phasing.info.list[[seg.name]]$samp_most_significant_by_diptest
      } else {
        all.seg.MSAI.by.BAF.ltnpd.present.list[[seg.name]]$MSAI_detection_feasible <- TRUE
        all.seg.MSAI.by.BAF.ltnpd.present.list[[seg.name]]$MSAI_detection_info <- "phasing_from_single_sample_used" # Note might be an issue here, is it only
      }
      nr <- length(asas.SNPS.lr.and.or.baf.lists[[1]][[seg.name]]$BAF)
      lower.band.col.names <- paste0(sample.vec, "_ltnpd_in_lower")
      upper.band.col.names <- paste0(sample.vec, "_ltnpd_in_upper")
      interleaved.col.names <- c(rbind(sample.vec, lower.band.col.names, upper.band.col.names))
      all.samp.seg.baf.df <- data.frame(matrix(vector(), nr, length(interleaved.col.names),
        dimnames = list(c(), interleaved.col.names)
      ),
      stringsAsFactors = FALSE,
      check.names = FALSE
      )

      for (samp in sample.vec) {
        all.samp.seg.baf.df[[samp]] <- asas.SNPS.lr.and.or.baf.lists[[samp]][[seg.name]]$BAF
        all.samp.seg.baf.df[[paste0(samp, "_ltnpd_in_lower")]] <- getLogTruncNormProbabilityDensitiesBAFvsLowerOrUpper(lower.or.upper.mu = 0.33, baf.vector = asas.SNPS.lr.and.or.baf.lists[[samp]][[seg.name]]$BAF)
        all.samp.seg.baf.df[[paste0(samp, "_ltnpd_in_upper")]] <- getLogTruncNormProbabilityDensitiesBAFvsLowerOrUpper(lower.or.upper.mu = 0.66, baf.vector = asas.SNPS.lr.and.or.baf.lists[[samp]][[seg.name]]$BAF)
      }

      seg.all.samp.MSAI.results.df <- data.frame(
        sample_name = sample.vec,
        sample_is_ref = logical(length(sample.vec)),
        MSAI_pval = numeric(length(sample.vec)),
        MSAI_pval_passed_thresh = logical(length(sample.vec)),
        poss_dominant_allele = rep("cpnA", length(sample.vec)), # Need to check whether we end up calling alleic imbalance or not..
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
      rownames(seg.all.samp.MSAI.results.df) <- sample.vec
      seg.all.samp.MSAI.results.df[seg.all.samp.MSAI.results.df$sample_name == ref.samp, ]$sample_is_ref <- TRUE
      seg.all.samp.MSAI.results.df[seg.all.samp.MSAI.results.df$sample_name == ref.samp, ]$MSAI_pval <- 1

      for (test.samp in sample.vec[ !sample.vec %in% ref.samp]) {
        msai.ltnpd.vec <- log(exp(all.samp.seg.baf.df[, paste0(test.samp, "_ltnpd_in_lower")] + all.samp.seg.baf.df[, paste0(ref.samp, "_ltnpd_in_upper")]) +
          exp(all.samp.seg.baf.df[, paste0(test.samp, "_ltnpd_in_upper")] + all.samp.seg.baf.df[, paste0(ref.samp, "_ltnpd_in_lower")]))
        no.msai.ltnpd.vec <- log(exp(all.samp.seg.baf.df[, paste0(test.samp, "_ltnpd_in_upper")] + all.samp.seg.baf.df[, paste0(ref.samp, "_ltnpd_in_upper")]) +
          exp(all.samp.seg.baf.df[, paste0(test.samp, "_ltnpd_in_lower")] + all.samp.seg.baf.df[, paste0(ref.samp, "_ltnpd_in_lower")]))

        ratio.msai.vs.no.msai.ltnpd.vec <- msai.ltnpd.vec - no.msai.ltnpd.vec
        seg.all.samp.MSAI.results.df[test.samp, ]$MSAI_pval <- tTestWithErrorHandling(ratio.msai.vs.no.msai.ltnpd.vec, alternative = "greater")$p.value
        # Now we can introduce another test to look for allelic imbalance here, if the tTest alternative is sent to

        if (seg.all.samp.MSAI.results.df[test.samp, ]$MSAI_pval <= MSAI.BAF.LTPND.pval.thresh) {
          seg.all.samp.MSAI.results.df[test.samp, ]$MSAI_pval_passed_thresh <- TRUE
          seg.all.samp.MSAI.results.df[test.samp, ]$poss_dominant_allele <- "cpnB"
        }
      }

      all.seg.MSAI.by.BAF.ltnpd.present.list[[seg.name]]$MSAI_present_any_sample_by_thresh <- any(seg.all.samp.MSAI.results.df$MSAI_pval_passed_thres)
      all.seg.MSAI.by.BAF.ltnpd.present.list[[seg.name]]$MSAI_by_sample_df <- seg.all.samp.MSAI.results.df
      all.seg.MSAI.by.BAF.ltnpd.present.list[[seg.name]]$all_samps_BAF_and_ltnpd_df <- seg.all.samp.MSAI.results.df
    } else {
      all.seg.MSAI.by.BAF.ltnpd.present.list[[seg.name]]$MSAI_detection_feasible <- FALSE
      all.seg.MSAI.by.BAF.ltnpd.present.list[[seg.name]]$MSAI_detection_info <- "no_phasing_was_attempted"
      all.seg.MSAI.by.BAF.ltnpd.present.list[[seg.name]]$MSAI_present_any_sample_by_thresh <- NA
    }
  }

  return(all.seg.MSAI.by.BAF.ltnpd.present.list)
}


getAllSegsMSAIDetectionFromBAFListLTNPD <- function(seg.name.vec, sample.vec, asas.SNPS.lr.and.or.baf.lists, all.seg.phasing.info.list, MSAI.BAF.LTPND.pval.thresh) {
  all.seg.MSAI.by.BAF.ltnpd.present.list <- list()
  for (seg.name in seg.name.vec) {
    all.seg.MSAI.by.BAF.ltnpd.present.list[[seg.name]] <- list()
    all.seg.MSAI.by.BAF.ltnpd.present.list[[seg.name]]$seg_name <- seg.name
    if (seg.name %in% names(all.seg.phasing.info.list)) {
      ref.samp <- all.seg.phasing.info.list[[seg.name]]$samp_chosen_to_phase_seg
      # Determine whether we can detect MSAI and fill in the info fields, resent ref.samp if mixture phasing was used

      if (grepl("none", ref.samp) == TRUE) {
        all.seg.MSAI.by.BAF.ltnpd.present.list[[seg.name]]$MSAI_detection_feasible <- FALSE
        all.seg.MSAI.by.BAF.ltnpd.present.list[[seg.name]]$MSAI_detection_info <- "no_reference_although_phasing_attempted"
        next
      } else if (ref.samp == "mixture") {
        all.seg.MSAI.by.BAF.ltnpd.present.list[[seg.name]]$MSAI_detection_feasible <- TRUE
        all.seg.MSAI.by.BAF.ltnpd.present.list[[seg.name]]$MSAI_detection_info <- "mixture_phasing_using_top_ranked_samp_by_diptest" # Note might be an issue here, is it only
        # samples that made it throigh to
        ref.samp <- all.seg.phasing.info.list[[seg.name]]$samp_most_significant_by_diptest
      } else {
        all.seg.MSAI.by.BAF.ltnpd.present.list[[seg.name]]$MSAI_detection_feasible <- TRUE
        all.seg.MSAI.by.BAF.ltnpd.present.list[[seg.name]]$MSAI_detection_info <- "phasing_from_single_sample_used" # Note might be an issue here, is it only
      }

      nr <- length(asas.SNPS.lr.and.or.baf.lists[[1]][[seg.name]]$BAF)
      lower.band.col.names <- paste0(sample.vec, "_ltnpd_in_lower")
      upper.band.col.names <- paste0(sample.vec, "_ltnpd_in_upper")
      interleaved.col.names <- c(rbind(sample.vec, lower.band.col.names, upper.band.col.names))
      all.samp.seg.baf.df <- data.frame(matrix(vector(), nr, length(interleaved.col.names),
        dimnames = list(c(), interleaved.col.names)
      ),
      stringsAsFactors = FALSE,
      check.names = FALSE
      )

      for (samp in sample.vec) {
        all.samp.seg.baf.df[[samp]] <- asas.SNPS.lr.and.or.baf.lists[[samp]][[seg.name]]$BAF
        all.samp.seg.baf.df[[paste0(samp, "_ltnpd_in_lower")]] <- getLogTruncNormProbabilityDensitiesBAFvsLowerOrUpper(lower.or.upper.mu = 0.33, baf.vector = asas.SNPS.lr.and.or.baf.lists[[samp]][[seg.name]]$BAF)
        all.samp.seg.baf.df[[paste0(samp, "_ltnpd_in_upper")]] <- getLogTruncNormProbabilityDensitiesBAFvsLowerOrUpper(lower.or.upper.mu = 0.66, baf.vector = asas.SNPS.lr.and.or.baf.lists[[samp]][[seg.name]]$BAF)
      }

      seg.all.samp.MSAI.results.df <- data.frame(
        sample_name = sample.vec,
        sample_is_ref = logical(length(sample.vec)),
        MSAI_pval = numeric(length(sample.vec)),
        MSAI_pval_passed_thresh = logical(length(sample.vec)),
        poss_dominant_allele = rep("cpnA", length(sample.vec)), # Need to check whether we end up calling alleic imbalance or not..
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
      rownames(seg.all.samp.MSAI.results.df) <- sample.vec
      seg.all.samp.MSAI.results.df[seg.all.samp.MSAI.results.df$sample_name == ref.samp, ]$sample_is_ref <- TRUE
      seg.all.samp.MSAI.results.df[seg.all.samp.MSAI.results.df$sample_name == ref.samp, ]$MSAI_pval <- 1

      for (test.samp in sample.vec[ !sample.vec %in% ref.samp]) {
        msai.ltnpd.vec <- log(exp(all.samp.seg.baf.df[, paste0(test.samp, "_ltnpd_in_lower")] + all.samp.seg.baf.df[, paste0(ref.samp, "_ltnpd_in_upper")]) +
          exp(all.samp.seg.baf.df[, paste0(test.samp, "_ltnpd_in_upper")] + all.samp.seg.baf.df[, paste0(ref.samp, "_ltnpd_in_lower")]))
        no.msai.ltnpd.vec <- log(exp(all.samp.seg.baf.df[, paste0(test.samp, "_ltnpd_in_upper")] + all.samp.seg.baf.df[, paste0(ref.samp, "_ltnpd_in_upper")]) +
          exp(all.samp.seg.baf.df[, paste0(test.samp, "_ltnpd_in_lower")] + all.samp.seg.baf.df[, paste0(ref.samp, "_ltnpd_in_lower")]))

        ratio.msai.vs.no.msai.ltnpd.vec <- msai.ltnpd.vec - no.msai.ltnpd.vec
        seg.all.samp.MSAI.results.df[test.samp, ]$MSAI_pval <- tTestWithErrorHandling(ratio.msai.vs.no.msai.ltnpd.vec, alternative = "greater")$p.value
        # Now we can introduce another test to look for allelic imbalance here, if the tTest alternative is sent to

        if (seg.all.samp.MSAI.results.df[test.samp, ]$MSAI_pval <= MSAI.BAF.LTPND.pval.thresh) {
          seg.all.samp.MSAI.results.df[test.samp, ]$MSAI_pval_passed_thresh <- TRUE
          seg.all.samp.MSAI.results.df[test.samp, ]$poss_dominant_allele <- "cpnB"
        }
      }

      all.seg.MSAI.by.BAF.ltnpd.present.list[[seg.name]]$MSAI_present_any_sample_by_thresh <- any(seg.all.samp.MSAI.results.df$MSAI_pval_passed_thres)
      all.seg.MSAI.by.BAF.ltnpd.present.list[[seg.name]]$MSAI_by_sample_df <- seg.all.samp.MSAI.results.df
      all.seg.MSAI.by.BAF.ltnpd.present.list[[seg.name]]$all_samps_BAF_and_ltnpd_df <- seg.all.samp.MSAI.results.df
    } else {
      all.seg.MSAI.by.BAF.ltnpd.present.list[[seg.name]]$MSAI_detection_feasible <- FALSE
      all.seg.MSAI.by.BAF.ltnpd.present.list[[seg.name]]$MSAI_detection_info <- "no_phasing_was_attempted"
      all.seg.MSAI.by.BAF.ltnpd.present.list[[seg.name]]$MSAI_present_any_sample_by_thresh <- NA
    }
  }

  return(all.seg.MSAI.by.BAF.ltnpd.present.list)
}


# Using the same sigma as was used in the NEJM
getLogTruncNormProbabilityDensitiesBAFvsLowerOrUpper <- function(lower.or.upper.mu, sigma.to.use = 0.07280109889280519, baf.vector) {
  nr <- length(baf.vector)
  samp.ltnpd.df <- data.frame(
    hz_snp_name = 1:nr,
    BAF = baf.vector,
    mu = rep(lower.or.upper.mu, nr),
    sigma = rep(sigma.to.use, nr),
    leftmost_truncation = rep(0, nr),
    rightmost_truncation = rep(1, nr),
    alpha = rep(0, nr),
    beta = rep(0, nr),
    baf_val_as_sds_from_mu = rep(0, nr),
    alpha_cdf = rep(0, nr),
    beta_cdf = rep(0, nr),
    baf_lpdf = rep(0, nr),
    ltn_prob_density = rep(0, nr),
    stringsAsFactors = FALSE
  )
  # Re-expressing values in terms of standard deviations from mean (mu) of normal distribution
  samp.ltnpd.df$alpha <- (samp.ltnpd.df$leftmost_truncation - samp.ltnpd.df$mu) / samp.ltnpd.df$sigma
  samp.ltnpd.df$beta <- (samp.ltnpd.df$rightmost_truncation - samp.ltnpd.df$mu) / samp.ltnpd.df$sigma
  samp.ltnpd.df$baf_val_as_sds_from_mu <- (samp.ltnpd.df$BAF - samp.ltnpd.df$mu) / samp.ltnpd.df$sigma

  samp.ltnpd.df$alpha_cdf <- pnorm(samp.ltnpd.df$alpha)
  samp.ltnpd.df$beta_cdf <- pnorm(samp.ltnpd.df$beta)
  samp.ltnpd.df$baf_lpdf <- log(dnorm(samp.ltnpd.df$baf_val_as_sds_from_mu) - log(samp.ltnpd.df$sigma)) # change height of distribution according to sigma
  samp.ltnpd.df$ltn_prob_density <- samp.ltnpd.df$baf_lpdf - log(samp.ltnpd.df$beta_cdf - samp.ltnpd.df$alpha_cdf) # change height of prob dist according to truncation

  return(samp.ltnpd.df$ltn_prob_density)
}


testPairOfSamplesSegForMSAIWithCPN <- function(seg.name, samp.one, samp.two, asas.ph.raw.cpn.lists) {
  phased.s.one.s.cpnA.cpnB.list <- asas.ph.raw.cpn.lists[[samp.one]][[seg.name]]
  phased.s.one.s.cpnA.cpnB.df <- as.matrix(cbind(phased.s.one.s.cpnA.cpnB.list$ph_cpnA_vec, phased.s.one.s.cpnA.cpnB.list$ph_cpnB_vec))
  phased.s.two.s.cpnA.cpnB.list <- asas.ph.raw.cpn.lists[[samp.one]][[seg.name]]
  phased.s.two.s.cpnA.cpnB.df <- as.matrix(cbind(phased.s.two.s.cpnA.cpnB.list$ph_cpnA_vec, phased.s.two.s.cpnA.cpnB.list$ph_cpnB_vec))

  if (tTestWithErrorHandling(phased.s.one.s.cpnA.cpnB.df[, 1], phased.s.one.s.cpnA.cpnB.df[, 2], paired = TRUE)$p.val < 0.001 &
    tTestWithErrorHandling(phased.s.two.s.cpnA.cpnB.df[, 1], phased.s.two.s.cpnA.cpnB.df[, 2], paired = TRUE)$p.val < 0.001) {
    if ((mean(phased.s.one.s.cpnA.cpnB.df[, 1] - phased.s.one.s.cpnA.cpnB.df[, 2]) > 0 &
      mean(phased.s.two.s.cpnA.cpnB.df[, 1] - phased.s.two.s.cpnA.cpnB.df[, 2]) < 0) |
      (mean(phased.s.one.s.cpnA.cpnB.df[, 2] - phased.s.one.s.cpnA.cpnB.df[, 1]) < 0 &
        mean(phased.s.two.s.cpnA.cpnB.df[, 2] - phased.s.two.s.cpnA.cpnB.df[, 1]) > 0)) {

      # s.one.two.comparison[[seg.name]][[samp.one.two.string]]$MSAI <- 1
      return(TRUE)
    }
  }

  return(FALSE)
}
