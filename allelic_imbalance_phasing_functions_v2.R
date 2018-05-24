##
##Params:
##
##Returns:
##
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

  for ( seg.name in seg.name.vec ) {

    first.ss.lr.baf.df <- asas.lr.baf.lists[[sample.vec[1]]][[seg.name]]
    if (nrow(first.ss.lr.baf.df) < min.num.hz.snp.phasing) {
      next;#If phasng not possible for this segment then move onto the next
    }#end if( nrow(ref.ss.lr.baf.df) < min.hz.snp.phasing)
    #For each sample how good is it at phasing all samples for this segment.
    all.samp.seg.phasing.ability.list <- getAssessmentOfEachSamplesAbilityToPhaseAllOthersForSeg(asas.lr.baf.lists = asas.lr.baf.lists,
                                                                                                sample.vec = sample.vec,
                                                                                                seg.name = seg.name,
                                                                                                ploidy.est.vec = ploidy.est.vec,
                                                                                                purity.est.vec = purity.est.vec,
                                                                                                gamma.val = gamma.val,
                                                                                                min.num.hz.snp.phasing = min.num.hz.snp.phasing)
    sample.to.phase.seg.list[[seg.name]]$all_samp_phasing_comparison_list <- all.samp.seg.phasing.ability.list
    a.s.sum.phas.df <- data.frame(mean_prop_SNPs_sep_across_samples = sapply(all.samp.seg.phasing.ability.list, function(x) x$mean_prop_SNPs_sep_across_samples),
                                  diptest_for_AI_pval = sapply(all.samp.seg.phasing.ability.list, function(x) x$diptest_for_AI_pval),
                                  ASCAT_called_as_AI = sapply(all.samp.seg.phasing.ability.list, function(x) x$ASCAT_called_as_AI),
                                  comb_pval_SNP_sep_across_samps = sapply(all.samp.seg.phasing.ability.list, function(x) x$comb_wilcox_sep_pval_across_samps),
                                  stringsAsFactors = FALSE)
    rownames(a.s.sum.phas.df) <- sample.vec
    #Reorder by which has best separation according to th eoriginal diptest, then subset to choose only samples that pass the threshold for having allelic imbalance.
    a.s.sum.phas.df <- a.s.sum.phas.df[with(a.s.sum.phas.df, order(diptest_for_AI_pval)), ]

    #Now select the samples that we will use to call for phasing.
    if (allow.dip.test.without.ASCAT.to.phase == TRUE) {

      poss.phas.s.df <- a.s.sum.phas.df[a.s.sum.phas.df$ASCAT_called_as_AI == TRUE | a.s.sum.phas.df$diptest_for_AI_pval < dip.test.p.val.threshold, ,drop = FALSE]

    } else {

      poss.phas.s.df <- a.s.sum.phas.df[a.s.sum.phas.df$ASCAT_called_as_AI == TRUE, ,drop = FALSE] 

    }

    #If s more than one sample that passes either ASCAT or diptest
    seg.mixture.phasing.list <- NULL
    if (nrow(poss.phas.s.df) > 0) {
      #initialise sthe seg mixture phasing to NULL so that, as a default, it won't be used:

      if (nrow(poss.phas.s.df) == 1) {
        #If only one sample passes then choose this sample as the reference
        sample.to.phase.seg <- rownames(poss.phas.s.df)[1]

      } else if (nrow(poss.phas.s.df) > 1) {
        #If there are multiple alance instances called choose the top one and check mixture phasing
        sample.to.phase.seg <- rownames(poss.phas.s.df)[1]
        seg.mixture.phasing.list <- getMixturePhasingResultsForSegment(sample.to.phase.seg, sample.vec, poss.phas.s.df, asas.lr.baf.lists, all.samp.seg.phasing.ability.list, seg.name, ploidy.est.vec, purity.est.vec, gamma.val, comb.wilcox.snp.sep.p.val.threshold, min.prop.snp.sep.threshold)

      }  else {

        stop("There's an error in phasing the phasing logic")

      }

    } else {#end if (nrow(poss.phas.s.df) > 0)

      sample.to.phase.seg <- "none_no_allelic_imbalance_called"

    }
    #Update output with the choice of phasing:
    if ( is.null(seg.mixture.phasing.list) ) {#mixture phasing was not required

      if( grepl("none", sample.to.phase.seg) == TRUE) {#so 
        #As no phasing was determined to be possible we will simply take alternate indices to be cpnA or cpnB as being in the upper band in the in the reference.
        all.phaseable.positions.idx.vec <- c()
        maj.is.cpnA.idx.vec <- all.samp.seg.phasing.ability.list[[1]]$maj_is_cpnA_idx_vec
        maj.is.cpnB.idx.vec <- all.samp.seg.phasing.ability.list[[1]]$maj_is_cpnB_idx_vec
        if(length(maj.is.cpnA.idx.vec[!is.na(maj.is.cpnA.idx.vec)]) > 0 ){
          all.phaseable.positions.idx.vec <- c(all.phaseable.positions.idx.vec, maj.is.cpnA.idx.vec)
        }
        if(length(maj.is.cpnB.idx.vec[!is.na(maj.is.cpnB.idx.vec)]) > 0 ){
          all.phaseable.positions.idx.vec <- c(all.phaseable.positions.idx.vec, maj.is.cpnB.idx.vec)
        }
        sorted.all.phas.idx.vec <- sort(all.phaseable.positions.idx.vec)
        sample.to.phase.seg.list[[seg.name]]$samp_chosen_to_phase_seg <- sample.to.phase.seg

        if(length(sorted.all.phas.idx.vec) > 0) {

          sample.to.phase.seg.list[[seg.name]]$maj_is_cpnA_idx_vec <- sorted.all.phas.idx.vec[seq(1, length(all.phaseable.positions.idx.vec),  2)]

        } else {

          sample.to.phase.seg.list[[seg.name]]$maj_is_cpnA_idx_vec <- NA
          
        } 
        if(length(sorted.all.phas.idx.vec) > 1) {

          sample.to.phase.seg.list[[seg.name]]$maj_is_cpnB_idx_vec <- sorted.all.phas.idx.vec[seq(2, length(all.phaseable.positions.idx.vec),  2)]

        } else {

          sample.to.phase.seg.list[[seg.name]]$maj_is_cpnB_idx_vec <- NA

        }

      } else {#end if( grepl("none", sample.to.phase.seg) == TRUE)
        #As there is a phasing possible Take the appropriate sample selected
        sample.to.phase.seg.list[[seg.name]]$samp_chosen_to_phase_seg <- sample.to.phase.seg
        sample.to.phase.seg.list[[seg.name]]$maj_is_cpnA_idx_vec <- all.samp.seg.phasing.ability.list[[sample.to.phase.seg]]$maj_is_cpnA_idx_vec
        sample.to.phase.seg.list[[seg.name]]$maj_is_cpnB_idx_vec <- all.samp.seg.phasing.ability.list[[sample.to.phase.seg]]$maj_is_cpnB_idx_vec

      }#end if( grepl("none", sample.to.phase.seg) == FALSE)

    } else {#end if ( is.null(seg.mixture.phasing.list) ) 
      #mixture phasing has worked/been performed and is required, take the appropriate results from the mixture phasing.
      sample.to.phase.seg.list[[seg.name]]$samp_chosen_to_phase_seg <- "mixture"
      sample.to.phase.seg.list[[seg.name]]$samp_most_significant_by_diptest <- seg.mixture.phasing.list$samp_most_significant_by_diptest
      sample.to.phase.seg.list[[seg.name]]$maj_is_cpnA_idx_vec <- seg.mixture.phasing.list$maj_is_cpnA_idx_vec
      sample.to.phase.seg.list[[seg.name]]$maj_is_cpnB_idx_vec <- seg.mixture.phasing.list$maj_is_cpnB_idx_vec
      sample.to.phase.seg.list[[seg.name]]$all_samp_mixture_phasing_comparison_list <- seg.mixture.phasing.list$all_samp_mixture_phasing_comparison_list
    }#end if( ! is.null(seg.mixture.phasing.list) )

  }#end for ( seg.name in seg.name.vec ) {

  return(sample.to.phase.seg.list)

}#end getCombinedBAFDipTestPhasing





##
##Params:
##
##Returns:
##
getAllSegsMSAIDetectionFromBAFListLTNPD <- function(seg.name.vec, sample.vec, asas.SNPS.lr.and.or.baf.lists, all.seg.phasing.info.list, MSAI.BAF.LTPND.pval.thresh) {

  all.seg.MSAI.by.BAF.ltnpd.present.list <- list()
  for ( seg.name in seg.name.vec ) {

    all.seg.MSAI.by.BAF.ltnpd.present.list[[seg.name]] <- list()  
    all.seg.MSAI.by.BAF.ltnpd.present.list[[seg.name]]$seg_name <- seg.name
    print(seg.name)
    if (seg.name %in% names(all.seg.phasing.info.list)) {

      ref.samp <- all.seg.phasing.info.list[[seg.name]]$samp_chosen_to_phase_seg
      #Determine whether we can detect MSAI and fill in the info fields, resent ref.samp if mixture phasing was used
      if ( grepl("none", ref.samp) == TRUE ) {

        all.seg.MSAI.by.BAF.ltnpd.present.list[[seg.name]]$MSAI_detection_feasible <- FALSE
        all.seg.MSAI.by.BAF.ltnpd.present.list[[seg.name]]$MSAI_detection_info <- "no_reference_although_phasing_attempted"
        next

      } else if (ref.samp == "mixture") {# end if ( grepl("none", ref.samp) == TRUE ) { 

        all.seg.MSAI.by.BAF.ltnpd.present.list[[seg.name]]$MSAI_detection_feasible <- TRUE
        all.seg.MSAI.by.BAF.ltnpd.present.list[[seg.name]]$MSAI_detection_info <- "mixture_phasing_using_top_ranked_samp_by_diptest"#Note might be an issue here, is it only
        #samples that made it throigh to 
        ref.samp <- all.seg.phasing.info.list[[seg.name]]$samp_most_significant_by_diptest

      } else {# end if ( grepl("none", ref.samp) == TRUE ) { 

        all.seg.MSAI.by.BAF.ltnpd.present.list[[seg.name]]$MSAI_detection_feasible <- TRUE
        all.seg.MSAI.by.BAF.ltnpd.present.list[[seg.name]]$MSAI_detection_info <- "phasing_from_single_sample_used"#Note might be an issue here, is it only

      }# end else if (ref.samp == "mixture") {
      nr <- length(asas.SNPs.w.both.lr.baf.snp.lists[[1]][[seg.name]]$BAF)
      lower.band.col.names <- paste0(sample.vec, "_ltnpd_in_lower")
      upper.band.col.names <- paste0(sample.vec, "_ltnpd_in_upper") 
      interleaved.col.names <- c(rbind(sample.vec, lower.band.col.names, upper.band.col.names))
      all.samp.seg.baf.df <- data.frame(matrix(vector(), nr, length(interleaved.col.names),
                              dimnames = list(c(), interleaved.col.names)),
                              stringsAsFactors = FALSE,
                              check.names = FALSE)

      for ( samp in sample.vec ) {

        all.samp.seg.baf.df[[samp]] <- asas.SNPs.w.both.lr.baf.snp.lists[[samp]][[seg.name]]$BAF
        all.samp.seg.baf.df[[paste0(samp, "_ltnpd_in_lower")]] <- getLogTruncNormProbabilityDensitiesBAFvsLowerOrUpper(lower.or.upper.mu = 0.33, baf.vector = asas.SNPs.w.both.lr.baf.snp.lists[[samp]][[seg.name]]$BAF)
        all.samp.seg.baf.df[[paste0(samp, "_ltnpd_in_upper")]] <- getLogTruncNormProbabilityDensitiesBAFvsLowerOrUpper(lower.or.upper.mu = 0.66, baf.vector = asas.SNPs.w.both.lr.baf.snp.lists[[samp]][[seg.name]]$BAF)

      }#end for ( samp in sample.vec ) {

      seg.all.samp.MSAI.results.df <- data.frame(sample_name = sample.vec,
                            sample_is_ref = logical(length(sample.vec)),
                            MSAI_pval = numeric(length(sample.vec)),
                            MSAI_pval_passed_thresh = logical(length(sample.vec)),
                            poss_dominant_allele = rep("cpnA", length(sample.vec)),#Need to check whether we end up calling alleic imbalance or not..
                            stringsAsFactors = FALSE,
                            check.names = FALSE)
      rownames(seg.all.samp.MSAI.results.df) <- sample.vec
      seg.all.samp.MSAI.results.df[seg.all.samp.MSAI.results.df$sample_name == ref.samp,]$sample_is_ref <- TRUE
      seg.all.samp.MSAI.results.df[seg.all.samp.MSAI.results.df$sample_name == ref.samp,]$MSAI_pval <- 1 

      for (test.samp in sample.vec[ ! sample.vec %in% ref.samp]) { 

        msai.ltnpd.vec <- log( exp(all.samp.seg.baf.df[, paste0(test.samp, "_ltnpd_in_lower")] + all.samp.seg.baf.df[, paste0(ref.samp, "_ltnpd_in_upper")]) +
                         exp(all.samp.seg.baf.df[, paste0(test.samp, "_ltnpd_in_upper")] + all.samp.seg.baf.df[, paste0(ref.samp, "_ltnpd_in_lower")]) )
        no.msai.ltnpd.vec <- log( exp(all.samp.seg.baf.df[, paste0(test.samp, "_ltnpd_in_upper")] + all.samp.seg.baf.df[, paste0(ref.samp, "_ltnpd_in_upper")]) +
                          exp(all.samp.seg.baf.df[, paste0(test.samp, "_ltnpd_in_lower")] + all.samp.seg.baf.df[, paste0(ref.samp, "_ltnpd_in_lower")]) )
        
        ratio.msai.vs.no.msai.ltnpd.vec <- msai.ltnpd.vec - no.msai.ltnpd.vec
        seg.all.samp.MSAI.results.df[test.samp,]$MSAI_pval <- tTestWithErrorHandling(ratio.msai.vs.no.msai.ltnpd.vec, alternative = "greater")$p.value

        if (seg.all.samp.MSAI.results.df[test.samp,]$MSAI_pval <= MSAI.BAF.LTPND.pval.thresh) {

          seg.all.samp.MSAI.results.df[test.samp,]$MSAI_pval_passed_thresh <- TRUE
          seg.all.samp.MSAI.results.df[test.samp,]$poss_dominant_allele <- "cpnB"

        }
        
      }#end for (test.samp in sample.vec[ ! sample.vec %in% ref.samp]) { 

      all.seg.MSAI.by.BAF.ltnpd.present.list[[seg.name]]$MSAI_present_any_sample_by_thresh <- any(seg.all.samp.MSAI.results.df[test.samp,]$MSAI_pval_passed_thres)
      all.seg.MSAI.by.BAF.ltnpd.present.list[[seg.name]]$MSAI_by_sample_df <- seg.all.samp.MSAI.results.df
      all.seg.MSAI.by.BAF.ltnpd.present.list[[seg.name]]$all_samps_BAF_and_ltnpd_df <- all.samp.seg.baf.df

    } else { #end if (seg.name %in% names(all.seg.phasing.info.list)) {

      all.seg.MSAI.by.BAF.ltnpd.present.list[[seg.name]]$MSAI_detection_feasible <- FALSE
      all.seg.MSAI.by.BAF.ltnpd.present.list[[seg.name]]$MSAI_detection_info <- "no_phasing_was_attempted"
      all.seg.MSAI.by.BAF.ltnpd.present.list[[seg.name]]$MSAI_present_any_sample_by_thresh <- NA
    }#end if ( ! seg.name %in% names(all.seg.phasing.info.list)) {

  }#end for ( seg.name in seg.name.vec ) {

  return(all.seg.MSAI.by.BAF.ltnpd.present.list)

}#end getAllSegsMirroredSubclonalAllelicImbalanceDetectionFromBAFList (...)


##
##Params:
##
##Returns:
##
getAllSegsMSAIDetectionFromBAFListLTNPD <- function(seg.name.vec, sample.vec, asas.SNPS.lr.and.or.baf.lists, all.seg.phasing.info.list, MSAI.BAF.LTPND.pval.thresh) {

  all.seg.MSAI.by.BAF.ltnpd.present.list <- list()
  for ( seg.name in seg.name.vec ) {

    all.seg.MSAI.by.BAF.ltnpd.present.list[[seg.name]] <- list()  
    all.seg.MSAI.by.BAF.ltnpd.present.list[[seg.name]]$seg_name <- seg.name
    print(seg.name)
    if (seg.name %in% names(all.seg.phasing.info.list)) {

      ref.samp <- all.seg.phasing.info.list[[seg.name]]$samp_chosen_to_phase_seg
      #Determine whether we can detect MSAI and fill in the info fields, resent ref.samp if mixture phasing was used
      if ( grepl("none", ref.samp) == TRUE ) {

        all.seg.MSAI.by.BAF.ltnpd.present.list[[seg.name]]$MSAI_detection_feasible <- FALSE
        all.seg.MSAI.by.BAF.ltnpd.present.list[[seg.name]]$MSAI_detection_info <- "no_reference_although_phasing_attempted"
        next

      } else if (ref.samp == "mixture") {# end if ( grepl("none", ref.samp) == TRUE ) { 

        all.seg.MSAI.by.BAF.ltnpd.present.list[[seg.name]]$MSAI_detection_feasible <- TRUE
        all.seg.MSAI.by.BAF.ltnpd.present.list[[seg.name]]$MSAI_detection_info <- "mixture_phasing_using_top_ranked_samp_by_diptest"#Note might be an issue here, is it only
        #samples that made it throigh to 
        ref.samp <- all.seg.phasing.info.list[[seg.name]]$samp_most_significant_by_diptest

      } else {# end if ( grepl("none", ref.samp) == TRUE ) { 

        all.seg.MSAI.by.BAF.ltnpd.present.list[[seg.name]]$MSAI_detection_feasible <- TRUE
        all.seg.MSAI.by.BAF.ltnpd.present.list[[seg.name]]$MSAI_detection_info <- "phasing_from_single_sample_used"#Note might be an issue here, is it only

      }# end else if (ref.samp == "mixture") {
      nr <- length(asas.SNPs.w.both.lr.baf.snp.lists[[1]][[seg.name]]$BAF)
      lower.band.col.names <- paste0(sample.vec, "_ltnpd_in_lower")
      upper.band.col.names <- paste0(sample.vec, "_ltnpd_in_upper") 
      interleaved.col.names <- c(rbind(sample.vec, lower.band.col.names, upper.band.col.names))
      all.samp.seg.baf.df <- data.frame(matrix(vector(), nr, length(interleaved.col.names),
                              dimnames = list(c(), interleaved.col.names)),
                              stringsAsFactors = FALSE,
                              check.names = FALSE)

      for ( samp in sample.vec ) {

        all.samp.seg.baf.df[[samp]] <- asas.SNPs.w.both.lr.baf.snp.lists[[samp]][[seg.name]]$BAF
        all.samp.seg.baf.df[[paste0(samp, "_ltnpd_in_lower")]] <- getLogTruncNormProbabilityDensitiesBAFvsLowerOrUpper(lower.or.upper.mu = 0.33, baf.vector = asas.SNPs.w.both.lr.baf.snp.lists[[samp]][[seg.name]]$BAF)
        all.samp.seg.baf.df[[paste0(samp, "_ltnpd_in_upper")]] <- getLogTruncNormProbabilityDensitiesBAFvsLowerOrUpper(lower.or.upper.mu = 0.66, baf.vector = asas.SNPs.w.both.lr.baf.snp.lists[[samp]][[seg.name]]$BAF)

      }#end for ( samp in sample.vec ) {

      seg.all.samp.MSAI.results.df <- data.frame(sample_name = sample.vec,
                            sample_is_ref = logical(length(sample.vec)),
                            MSAI_pval = numeric(length(sample.vec)),
                            MSAI_pval_passed_thresh = logical(length(sample.vec)),
                            poss_dominant_allele = rep("cpnA", length(sample.vec)),#Need to check whether we end up calling alleic imbalance or not..
                            stringsAsFactors = FALSE,
                            check.names = FALSE)
      rownames(seg.all.samp.MSAI.results.df) <- sample.vec
      seg.all.samp.MSAI.results.df[seg.all.samp.MSAI.results.df$sample_name == ref.samp,]$sample_is_ref <- TRUE
      seg.all.samp.MSAI.results.df[seg.all.samp.MSAI.results.df$sample_name == ref.samp,]$MSAI_pval <- 1 

      for (test.samp in sample.vec[ ! sample.vec %in% ref.samp]) { 

        msai.ltnpd.vec <- log( exp(all.samp.seg.baf.df[, paste0(test.samp, "_ltnpd_in_lower")] + all.samp.seg.baf.df[, paste0(ref.samp, "_ltnpd_in_upper")]) +
                         exp(all.samp.seg.baf.df[, paste0(test.samp, "_ltnpd_in_upper")] + all.samp.seg.baf.df[, paste0(ref.samp, "_ltnpd_in_lower")]) )
        no.msai.ltnpd.vec <- log( exp(all.samp.seg.baf.df[, paste0(test.samp, "_ltnpd_in_upper")] + all.samp.seg.baf.df[, paste0(ref.samp, "_ltnpd_in_upper")]) +
                          exp(all.samp.seg.baf.df[, paste0(test.samp, "_ltnpd_in_lower")] + all.samp.seg.baf.df[, paste0(ref.samp, "_ltnpd_in_lower")]) )
        
        ratio.msai.vs.no.msai.ltnpd.vec <- msai.ltnpd.vec - no.msai.ltnpd.vec
        seg.all.samp.MSAI.results.df[test.samp,]$MSAI_pval <- tTestWithErrorHandling(ratio.msai.vs.no.msai.ltnpd.vec, alternative = "greater")$p.value

        if (seg.all.samp.MSAI.results.df[test.samp,]$MSAI_pval <= MSAI.BAF.LTPND.pval.thresh) {

          seg.all.samp.MSAI.results.df[test.samp,]$MSAI_pval_passed_thresh <- TRUE
          seg.all.samp.MSAI.results.df[test.samp,]$poss_dominant_allele <- "cpnB"

        }
        
      }#end for (test.samp in sample.vec[ ! sample.vec %in% ref.samp]) { 

      all.seg.MSAI.by.BAF.ltnpd.present.list[[seg.name]]$MSAI_present_any_sample_by_thresh <- any(seg.all.samp.MSAI.results.df[test.samp,]$MSAI_pval_passed_thres)
      all.seg.MSAI.by.BAF.ltnpd.present.list[[seg.name]]$MSAI_by_sample_df <- seg.all.samp.MSAI.results.df
      all.seg.MSAI.by.BAF.ltnpd.present.list[[seg.name]]$all_samps_BAF_and_ltnpd_df <- seg.all.samp.MSAI.results.df

    } else { #end if (seg.name %in% names(all.seg.phasing.info.list)) {

      all.seg.MSAI.by.BAF.ltnpd.present.list[[seg.name]]$MSAI_detection_feasible <- FALSE
      all.seg.MSAI.by.BAF.ltnpd.present.list[[seg.name]]$MSAI_detection_info <- "no_phasing_was_attempted"
      all.seg.MSAI.by.BAF.ltnpd.present.list[[seg.name]]$MSAI_present_any_sample_by_thresh <- NA
    }#end if ( ! seg.name %in% names(all.seg.phasing.info.list)) {

  }#end for ( seg.name in seg.name.vec ) {

  return(all.seg.MSAI.by.BAF.ltnpd.present.list)

}#end getAllSegsMirroredSubclonalAllelicImbalanceDetectionFromBAFList (...)


getAllSampLogrAlonePositionsLogrAndTotalCPNDFsList <- function(sample.vec, asas.SNPS.lr.and.or.baf.lists, ploidy.est.vec, purity.est.vec, gamma.val) {

  all.samp.logr.df.list <- list()

  for ( samp in sample.vec ) {

    samp.logr.baf.df <- do.call(rbind, asas.SNPS.lr.and.or.baf.lists[[samp]])

    if ( nrow(samp.logr.baf.df[is.na(samp.logr.baf.df$BAF) ,]) > 0 ) {

      samp.logr.df <- samp.logr.baf.df[ ! is.na(samp.logr.baf.df$BAF) ,]

    } else if ( nrow(samp.logr.baf.df[samp.logr.baf.df$BAF == 2 ,]) > 0) {

      samp.logr.df <- samp.logr.baf.df[ samp.logr.baf.df$BAF == 2 ,]

    } else {

      print("No logr only SNPs, 0 row dataframe will be passed to functions.")
      samp.logr.df <- samp.logr.baf.df[ samp.logr.baf.df$BAF == 2 ,]

    }


    if (nrow(samp.logr.baf.df) > 0 ) {
      print("head(samp.logr.df)")
      print(head(samp.logr.df))
      mock.cpnA.cpnB.logr.df <- calculateSNPsCopyNumberAB(logr.vec = samp.logr.df$logR,
                                                            baf.vec = rep(0.5, nrow(samp.logr.df)),
                                                            tumor.ploidy = ploidy.est.vec[samp],
                                                            tumor.purity = purity.est.vec[samp],
                                                            gamma.val = gamma.val)
      samp.logr.df$tot_raw_un_ph_cpn <- unlist(rowSums(mock.cpnA.cpnB.logr.df))
      print("head(samp.logr.df)")
      print(head(samp.logr.df))
      #samp.logr.df$tot_raw_un_ph_cpn <- rowSums(mock.cpnA.cpnB.logr.df)

    }

    all.samp.logr.df.list[[samp]] <- samp.logr.df
    
  }


  return(all.samp.logr.df.list)

}#end getAllSampLogrAlonePositionsLogrAndTotalCPNDFsList





#Currently not used, but could be an option to allow sigmas to be set baded on read depth.. todo, test this
calculateBinonmialNoiseBasedOnReadDepth <- function(read.depth, mu = 0.5, pessimism.factor = 2) {
  #read.depth is the number of trials
  binomial.noise <- sqrt(read.depth * mu * (1 - mu))
  corrected.binomial.noise <- binomial.noise * pessimism.factor
  
  return(corrected.binomial.noise)

}#end calculateBinonmialNoiseBasedOnReadDepth


##
##Params:
##
##Returns:
##
getLogTruncNormProbabilityDensitiesBAFvsLowerOrUpper <- function(lower.or.upper.mu, sigma.to.use = 0.07280109889280519 , baf.vector) {#Using the same sigma as was used in the NEJM
  nr <- length(baf.vector)
  samp.ltnpd.df <- data.frame(hz_snp_name = 1:nr,
                         BAF = baf.vector,
                         mu = rep(lower.or.upper.mu, nr),
                         sigma = rep(sigma.to.use, nr),
                         leftmost_truncation = rep(0, nr),
                         rightmost_truncation = rep(1, nr),
                         alpha = rep(0,  nr),
                         beta = rep(0, nr),
                         baf_val_as_sds_from_mu = rep(0, nr),
                         alpha_cdf = rep(0, nr),
                         beta_cdf = rep(0, nr),
                         baf_lpdf = rep(0, nr),
                         ltn_prob_density = rep(0, nr),
                         stringsAsFactors = FALSE)
  #Re-expressing values in terms of standard deviations from mean (mu) of normal distribution
  samp.ltnpd.df$alpha <- (samp.ltnpd.df$leftmost_truncation - samp.ltnpd.df$mu) / samp.ltnpd.df$sigma
  samp.ltnpd.df$beta <- (samp.ltnpd.df$rightmost_truncation -  samp.ltnpd.df$mu) / samp.ltnpd.df$sigma
  samp.ltnpd.df$baf_val_as_sds_from_mu <- (samp.ltnpd.df$BAF - samp.ltnpd.df$mu) / samp.ltnpd.df$sigma

  samp.ltnpd.df$alpha_cdf <- pnorm(samp.ltnpd.df$alpha)
  samp.ltnpd.df$beta_cdf <- pnorm(samp.ltnpd.df$beta)
  samp.ltnpd.df$baf_lpdf <- log(dnorm(samp.ltnpd.df$baf_val_as_sds_from_mu) - log(samp.ltnpd.df$sigma))#change height of distribution according to sigma
  samp.ltnpd.df$ltn_prob_density <- samp.ltnpd.df$baf_lpdf - log(samp.ltnpd.df$beta_cdf - samp.ltnpd.df$alpha_cdf)#change height of prob dist according to truncation
  
  return(samp.ltnpd.df$ltn_prob_density)
}




##
##Params:
##
##Returns:
##
calculateSNPsRawCopyNumberAB <- function(logr.vec, baf.vec, tumor.ploidy, tumor.purity, gamma.val = 1) {

  cpnA <- (tumor.purity - 1 + baf.vec * 2 ^ (logr.vec / gamma.val) * ((1 - tumor.purity) * 2 + tumor.purity * tumor.ploidy)) / tumor.purity
  cpnB <- (tumor.purity - 1 - (baf.vec - 1) * 2 ^ (logr.vec / gamma.val) * ((1 - tumor.purity) * 2 + tumor.purity * tumor.ploidy)) / tumor.purity

  return(cbind(cpnA, cpnB))

}#end calculateSNPCopyNumberAB <- function(


##
##Params:
##
##Returns:
##
assessPhasingAcrossSelectedSamples <- function(seg.name,
                                              sample.comparison.list,
                                              samp,
                                              asas.lr.baf.lists,
                                              samples.to.test.vec,
                                              ref.ss.maj.is.cpnA.idx.vec,
                                              ref.ss.maj.is.cpnB.idx.vec,
                                              ploidy.est.vec,
                                              purity.est.vec,
                                              gamma.val) {

  other.ss.phasing.df <- data.frame(separation_wilcox_pval = rep(1, length(samples.to.test.vec)),
                                    prop_separated = numeric(length(samples.to.test.vec)),
                                    concordant_w_ref = numeric(length(samples.to.test.vec)),
                                    stringsAsFactors = FALSE)
  rownames(other.ss.phasing.df) <- samples.to.test.vec
  
  for ( other.samp in samples.to.test.vec ) {

    oth.ss.lr.baf.df <- asas.lr.baf.lists[[other.samp]][[seg.name]]

    oth.ss.cpnA.cpnB.df <- calculateSNPsCopyNumberAB(logr.vec=oth.ss.lr.baf.df$logR, baf.vec=oth.ss.lr.baf.df$BAF, tumor.ploidy=ploidy.est.vec[other.samp], tumor.purity=purity.est.vec[other.samp], gamma.val=gamma.val)

    #old code had some checks here to make sure that these have the same number of rows, but if everything holds they should remain the same anyway...
    other.samp.wilcox.test <-  wilcox.test( c(oth.ss.cpnA.cpnB.df[ref.ss.maj.is.cpnA.idx.vec,  1], oth.ss.cpnA.cpnB.df[ref.ss.maj.is.cpnB.idx.vec, 2]), 
                                            c(oth.ss.cpnA.cpnB.df[ref.ss.maj.is.cpnA.idx.vec,  2], oth.ss.cpnA.cpnB.df[ref.ss.maj.is.cpnB.idx.vec, 1]), paired=TRUE)

    oth.ss.snps.same.as.ref <- length(which(c(oth.ss.cpnA.cpnB.df[ref.ss.maj.is.cpnA.idx.vec, 1],oth.ss.cpnA.cpnB.df[ref.ss.maj.is.cpnB.idx.vec, 2]) - 
                                            c(oth.ss.cpnA.cpnB.df[ref.ss.maj.is.cpnA.idx.vec, 2],oth.ss.cpnA.cpnB.df[ref.ss.maj.is.cpnB.idx.vec, 1]) > 0)) / nrow(oth.ss.cpnA.cpnB.df)
    oth.ss.snps.diff.to.ref <- length(which(c(oth.ss.cpnA.cpnB.df[ref.ss.maj.is.cpnA.idx.vec, 1],oth.ss.cpnA.cpnB.df[ref.ss.maj.is.cpnB.idx.vec, 2]) - 
                                            c(oth.ss.cpnA.cpnB.df[ref.ss.maj.is.cpnA.idx.vec, 2],oth.ss.cpnA.cpnB.df[ref.ss.maj.is.cpnB.idx.vec, 1]) < 0)) / nrow(oth.ss.cpnA.cpnB.df)
    #Updating the other sample segment phasing df
    other.ss.phasing.df[other.samp,] <- c(other.samp.wilcox.test$p.val, max(oth.ss.snps.same.as.ref, oth.ss.snps.diff.to.ref), ifelse(oth.ss.snps.same.as.ref >= oth.ss.snps.diff.to.ref,1,-1))

  }#end for ( other.samp in sample.vec[ ! sample%in%sample.vec] )

  sample.comparison.list[[samp]]$comb_wilcox_sep_pval_across_samps <- pchisq((sum(log(other.ss.phasing.df$separation_wilcox_pval))*-2), df=length(other.ss.phasing.df$separation_wilcox_pval)*2, lower.tail=F)
  sample.comparison.list[[samp]]$oth_samp_wilcox_sep_pvals_vec <- other.ss.phasing.df$separation_wilcox_pval
  sample.comparison.list[[samp]]$oth_samp_SNPs_sep_prop_vec <- other.ss.phasing.df$prop_separated
  sample.comparison.list[[samp]]$maj_is_cpnA_idx_vec <- ref.ss.maj.is.cpnA.idx.vec
  sample.comparison.list[[samp]]$maj_is_cpnB_idx_vec <- ref.ss.maj.is.cpnB.idx.vec
  sample.comparison.list[[samp]]$mean_prop_SNPs_sep_across_samples <- mean(other.ss.phasing.df$prop_separated)
  sample.comparison.list[[samp]]$oth_samp_cpnA_cpnB_match_ref_vec <- other.ss.phasing.df$concordant_w_ref


  return(sample.comparison.list)
}


#So here is our first subfunction:
#For this segment go through and calculate how good each sample is at phasing all the other samples
#This will give us a list that for each sample contains another list comparing its ability to phase all others via a dip test on the major cpn vs all other samples and the wilcox tests on all others using
#the "phasing" given by the reference sample.
##
##Params:
##
##Returns:
##
getAssessmentOfEachSamplesAbilityToPhaseAllOthersForSeg <- function(asas.lr.baf.lists,
                                                                    sample.vec,
                                                                    seg.name,
                                                                    ploidy.est.vec,
                                                                    purity.est.vec,
                                                                    gamma.val=1,
                                                                    min.num.hz.snp.phasing=5){
  #other sample segment phasing df, used when working out the optimum sample in all regions, for each segment the same dataframe is reused to save time allocating memory...
  num.other.samps <- length(sample.vec) - 1
  other.ss.phasing.df <- data.frame(separation_wilcox_pval = rep(1, num.other.samps),
                                    prop_separated = numeric(num.other.samps),
                                    concordant_w_ref = numeric(num.other.samps),
                                    stringsAsFactors = FALSE)
  sample.comparison.list <- list()

  for ( samp in sample.vec ) {

      #ref sample seg logr baf df
      ref.ss.lr.baf.df <- asas.lr.baf.lists[[samp]][[seg.name]]
      if ( nrow(ref.ss.lr.baf.df) < min.num.hz.snp.phasing) {

        next;#Might actally want this to be a break as they should all have too few snps as it is joint segmentation with the same SNPs...

      }#end if( nrow(ref.ss.lr.baf.df) < min.hz.snp.phasing)
      #renaming the dataframe rows for the other samples, now we are trying a new reference sample
      rownames(other.ss.phasing.df) <- sample.vec[ ! sample.vec %in% samp ]

      #Get phasing, reference sample segment cpn A cpn B
      ref.ss.cpnA.cpnB.df <- calculateSNPsRawCopyNumberAB(logr.vec=ref.ss.lr.baf.df$logR,
                                                          baf.vec=ref.ss.lr.baf.df$BAF,
                                                          tumor.ploidy=ploidy.est.vec[samp],
                                                          tumor.purity=purity.est.vec[samp],
                                                          gamma.val=gamma.val)
      ref.ss.max.cpn.vec <- apply(ref.ss.cpnA.cpnB.df, 1, max)
      ref.ss.is.max.cpn.cpnA.cpnB.logical.df <- ifelse(ref.ss.cpnA.cpnB.df>=ref.ss.max.cpn.vec, 1, 0)
      ref.ss.maj.is.cpnA.idx.vec <- which(ref.ss.is.max.cpn.cpnA.cpnB.logical.df[,1]==1)
      ref.ss.maj.is.cpnB.idx.vec <- which(ref.ss.is.max.cpn.cpnA.cpnB.logical.df[,2]==1)

      sample.comparison.list[[samp]]$diptest_for_AI_pval <- as.numeric(dip.test(ref.ss.lr.baf.df$BAF)$p.val)
      sample.comparison.list[[samp]]$ASCAT_called_as_AI <- ifelse(ref.ss.lr.baf.df$BAFsegmented[1]!=0.5,TRUE,FALSE)
      samples.to.test.vec <- sample.vec[ ! sample.vec %in%  samp]
      sample.comparison.list <- assessPhasingAcrossSelectedSamples(seg.name = seg.name,
                                                                  sample.comparison.list = sample.comparison.list,
                                                                  samp = samp,
                                                                  asas.lr.baf.lists = asas.lr.baf.lists,
                                                                  samples.to.test.vec = samples.to.test.vec,
                                                                  ref.ss.maj.is.cpnA.idx.vec = ref.ss.maj.is.cpnA.idx.vec,
                                                                  ref.ss.maj.is.cpnB.idx.vec = ref.ss.maj.is.cpnB.idx.vec,
                                                                  ploidy.est.vec = ploidy.est.vec,
                                                                  purity.est.vec = purity.est.vec,
                                                                  gamma.val = gamma.val)

  }#end for( samp in sample.vec ) 

  return(sample.comparison.list)

}#end of getAssessmentOfEachSamplesAbilityToPhaseAllOthersForSeg <- function(asas.lr.baf.lists, sample.vec, seg.idx)




##
##Params:
##
##Returns:
##
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
                                              min.prop.snp.sep.threshold) {
  #Check if the mixture phsing might actually improve things, if not continue on
  if (as.numeric(poss.phas.s.df[sample.to.phase.seg, ,drop = FALSE]$comb_pval_SNP_sep_across_samps) < as.numeric(comb.wilcox.snp.sep.p.val.threshold)
     & as.numeric(poss.phas.s.df[sample.to.phase.seg, ,drop = FALSE]$mean_prop_SNPs_sep_across_samples) > as.numeric(min.prop.snp.sep.threshold)) {
    
    return(NULL)
    
  }#end check of wthether mixture phasing should be applied.
  phas.ss.lr.baf.df <- asas.lr.baf.lists[[sample.to.phase.seg]][[seg.name]]
  phas.ss.cpnA.cpnB.df <- calculateSNPsCopyNumberAB(logr.vec=phas.ss.lr.baf.df$logR, baf.vec=phas.ss.lr.baf.df$BAF, tumor.ploidy=ploidy.est.vec[sample.to.phase.seg], tumor.purity=purity.est.vec[sample.to.phase.seg], gamma.val=gamma.val)
  other.samps.vec <- sample.vec[ ! sample.vec %in% sample.to.phase.seg ]
  oth.s.cnpA.cpnB.match.ref.vec <- all.samp.seg.phasing.ability.list[[sample.to.phase.seg]]$oth_samp_cpnA_cpnB_match_ref_vec
  other.passing.samples.vec <- other.samps.vec[other.samps.vec %in% rownames(poss.phas.s.df)]
  #Adding in a 1 to represent the sample to phase seg value and that its orrientation matches its own orientation
  oth.s.cnpA.cpnB.match.ref.vec <- c(1, oth.s.cnpA.cpnB.match.ref.vec)
  all.passing.samples.vec <- c(sample.to.phase.seg, other.passing.samples.vec)
  names(oth.s.cnpA.cpnB.match.ref.vec) <- all.passing.samples.vec

  new.maj.is.cpnA.idx.vec  <- c()
  new.maj.is.cpnB.idx.vec  <- c()
  change.from.ref <- c()
  #Now, for each SNP work out whether it should be major or minor across all samples...?
  for ( SNP.idx in 1:nrow(phas.ss.cpnA.cpnB.df) ) {

    if (SNP.idx %in% all.samp.seg.phasing.ability.list[[sample.to.phase.seg]]$maj_is_cpnA_idx_vec) {
      current.maj.cpn.col.idx <- 1
      current.min.cpn.col.idx <- 2
    }#end if (SNP.idx%in%all.samp.seg.phasing.ability.list[[sample.to.phase.seg]]$maj_is_cpnA_idx_vec)

    if (SNP.idx %in% all.samp.seg.phasing.ability.list[[sample.to.phase.seg]]$maj_is_cpnB_idx_vec) {
      current.maj.cpn.col.idx <- 2
      current.min.cpn.col.idx <- 1
    }# end if (SNP.idx%in%all.samp.seg.phasing.ability.list[[sample.to.phase.seg]]$maj_is_cpnB_idx_vec)
    
    all.s.SNP.sep.scores <- c()

    for ( samp in all.passing.samples.vec ) {

      ss.lr.baf.df <- asas.lr.baf.lists[[samp]][[seg.name]]
      ss.cpnA.cpnB.df <- calculateSNPsCopyNumberAB(logr.vec=ss.lr.baf.df$logR, baf.vec=ss.lr.baf.df$BAF, tumor.ploidy=ploidy.est.vec[samp], tumor.purity=purity.est.vec[samp], gamma.val=gamma.val)
      ss.SNP.cpnA <- ss.cpnA.cpnB.df[SNP.idx, current.maj.cpn.col.idx]
      ss.SNP.cpnB <- ss.cpnA.cpnB.df[SNP.idx, current.min.cpn.col.idx]

      s.cpnAB.match.ref.factor <- oth.s.cnpA.cpnB.match.ref.vec[samp]
      #sample snp sep score
      s.SNP.sep.score <- (ss.SNP.cpnA - ss.SNP.cpnB) * s.cpnAB.match.ref.factor
      all.s.SNP.sep.scores <- c(all.s.SNP.sep.scores, s.SNP.sep.score)

    }#end for ( samp in all.passing.samples.vec ) {

    #
    if (mean(all.s.SNP.sep.scores) >= 0) {

      change.from.ref <- c(change.from.ref, 0)

      if (current.maj.cpn.col.idx == 1) {

        new.maj.is.cpnA.idx.vec <- c(new.maj.is.cpnA.idx.vec, SNP.idx)

      } else if (current.maj.cpn.col.idx == 2) {

        new.maj.is.cpnB.idx.vec <- c(new.maj.is.cpnB.idx.vec, SNP.idx)

      }

    }#end if (mean(all.SNP.sep.score) > 0) 

    if (mean(all.s.SNP.sep.scores) < 0) {

      change.from.ref <- c(change.from.ref, 1)

      if (current.maj.cpn.col.idx == 1) {

        new.maj.is.cpnB.idx.vec <- c(new.maj.is.cpnB.idx.vec, SNP.idx)

      } else if (current.maj.cpn.col.idx == 2) {

        new.maj.is.cpnA.idx.vec <- c(new.maj.is.cpnA.idx.vec, SNP.idx)

      }

    }#end if (mean(all.SNP.sep.score) < 0) 

  }#end for ( SNP.idx in 1:nrow(phas.ss.cpnA.cpnB.df) )
  #So have now come up with whether whether cpnA or cpnB are the Major cpn through all samples allowing for MSAI
  output.list <- list()
  output.list$sample_to_phase_seg <- "mixture"
  output.list$samp_most_significant_by_diptest <- rownames(poss.phas.s.df)[1]
  output.list$maj_is_cpnA_idx_vec <- new.maj.is.cpnA.idx.vec
  output.list$maj_is_cpnB_idx_vec <- new.maj.is.cpnB.idx.vec
  #Now check the ability of the combined mixture snp major cpn A or B analysis accross all samples
  output.list$all_samp_mixture_phasing_comparison_list <- assessPhasingAcrossSelectedSamples(seg.name=seg.name,sample.comparison.list=list(),
                                                                                                              samp = sample.to.phase.seg,
                                                                                                              asas.lr.baf.lists,
                                                                                                              sample.vec,
                                                                                                              ref.ss.maj.is.cpnA.idx.vec = new.maj.is.cpnA.idx.vec,
                                                                                                              ref.ss.maj.is.cpnB.idx.vec = new.maj.is.cpnB.idx.vec,
                                                                                                              ploidy.est.vec,
                                                                                                              purity.est.vec,
                                                                                                              gamma.val = gamma.val)

  return(output.list)
}




getAllSegAllSamplePhasedRawCPNList <- function(sample.vec, ploidy.est.vec, purity.est.vec, min.cons.segs.df, all.seg.phasing.info.list, asas.both.lr.baf.lists, gamma.val = 1, min.num.hz.snp.phasing = 5) {

    seg.name.vec <- paste(min.cons.segs.df[ , 1], min.cons.segs.df[ , 2], min.cons.segs.df[ , 3], sep=":")
    #Initialise all lists 
    all.samples.phased.cpn.list <- list()
    for ( samp in sample.vec ) {
      all.samples.phased.cpn.list[[samp]] <- list()
      for ( seg.name in seg.name.vec ) {
        all.samples.phased.cpn.list[[samp]][[seg.name]] <- list()
        all.samples.phased.cpn.list[[samp]][[seg.name]]$chr_vec <- NA 
        all.samples.phased.cpn.list[[samp]][[seg.name]]$pos_vec <- NA
        all.samples.phased.cpn.list[[samp]][[seg.name]]$logr_vec <- NA
        all.samples.phased.cpn.list[[samp]][[seg.name]]$ph_cpnA_vec <- NA
        all.samples.phased.cpn.list[[samp]][[seg.name]]$ph_cpnB_vec <- NA
        all.samples.phased.cpn.list[[samp]][[seg.name]]$mean_ph_cpnA <- NA
        all.samples.phased.cpn.list[[samp]][[seg.name]]$mean_ph_cpnB <- NA
        all.samples.phased.cpn.list[[samp]][[seg.name]]$un_ph_cpn_total
        all.samples.phased.cpn.list[[samp]][[seg.name]]$un_ph_cpn_total_mean <- NA
    }#end for ( seg.name in seg.name.vec ) {
  }#end for ( samp in sample.vec ) {
  #Now add all phased copy number usin gthe all.se.phasing.info.list
    for ( samp in sample.vec) {

      for ( seg.name in seg.name.vec ) {
        #We may have the problem here that some segments don't have any snps in them at all, check for this perhaps? see if fails
          samp.b.lr.baf.df <- asas.both.lr.baf.lists[[samp]][[seg.name]]
          if (nrow(samp.b.lr.baf.df) == 0) {
            next
          }#end   if (nrow(samp.b.lr.baf.df) == 0) {

          samp.seg.cpnA.cpnB.df <- calculateSNPsCopyNumberAB(logr.vec = samp.b.lr.baf.df$logR,
                                                            baf.vec = samp.b.lr.baf.df$BAF,
                                                            tumor.ploidy = ploidy.est.vec[samp],
                                                            tumor.purity = purity.est.vec[samp],
                                                            gamma.val = gamma.val)
   
        if (nrow(samp.b.lr.baf.df) < min.num.hz.snp.phasing) {

          all.samples.phased.cpn.list[[samp]][[seg.name]]$chr_vec <- samp.b.lr.baf.df$chrom
          all.samples.phased.cpn.list[[samp]][[seg.name]]$pos_vec <- samp.b.lr.baf.df$pos
          all.samples.phased.cpn.list[[samp]][[seg.name]]$logr_vec <- samp.b.lr.baf.df$logR
          all.samples.phased.cpn.list[[samp]][[seg.name]]$un_ph_cpn_total <- samp.seg.cpnA.cpnB.df[, 1] + samp.seg.cpnA.cpnB.df[, 2]
          all.samples.phased.cpn.list[[samp]][[seg.name]]$un_ph_cpn_total_mean <- mean(samp.seg.cpnA.cpnB.df[, 1] + samp.seg.cpnA.cpnB.df[, 2])


        } else {#end if (nrow(samp.b.lr.baf.df) > min.num.hz.snp.phasing) { now applying phasing to those that were phased.

          samp.seg.upper.band.cpn.vals.vec <- apply(samp.seg.cpnA.cpnB.df, 1, max)
          samp.seg.median.upper.band.cpn <- median(samp.seg.upper.band.cpn.vals.vec)
          samp.seg.cpnA.cpnB.is.upper.vec <- ifelse(samp.seg.cpnA.cpnB.df == samp.seg.upper.band.cpn.vals.vec, 1, 0)
          phased.samp.seg.cpnA.cpnB.df <- samp.seg.cpnA.cpnB.df
          #Get indices to phase cpnA
          if (length(all.seg.phasing.info.list[[seg.name]]$maj_is_cpnA_idx_vec) > 0) {

               is.upper.and.cpnA.in.ref.idx.vec <- all.seg.phasing.info.list[[seg.name]]$maj_is_cpnA_idx_vec
               phased.samp.seg.cpnA.cpnB.df[is.upper.and.cpnA.in.ref.idx.vec, 1] <- phased.samp.seg.cpnA.cpnB.df[is.upper.and.cpnA.in.ref.idx.vec, 1]
               phased.samp.seg.cpnA.cpnB.df[is.upper.and.cpnA.in.ref.idx.vec, 2] <- phased.samp.seg.cpnA.cpnB.df[is.upper.and.cpnA.in.ref.idx.vec, 2]

          }
          #Get indices to phase cpnB
          if (length(all.seg.phasing.info.list[[seg.name]]$maj_is_cpnB_idx_vec) > 0) {

               is.lower.and.cpnB.in.ref.idx.vec <- all.seg.phasing.info.list[[seg.name]]$maj_is_cpnB_idx_vec
               phased.samp.seg.cpnA.cpnB.df[is.lower.and.cpnB.in.ref.idx.vec, 1] <- phased.samp.seg.cpnA.cpnB.df[is.lower.and.cpnB.in.ref.idx.vec, 2]
               phased.samp.seg.cpnA.cpnB.df[is.lower.and.cpnB.in.ref.idx.vec, 2] <- phased.samp.seg.cpnA.cpnB.df[is.lower.and.cpnB.in.ref.idx.vec, 1]
               
          } 
          #Store this information in the list:
          all.samples.phased.cpn.list[[samp]][[seg.name]]$chr_vec <- samp.b.lr.baf.df$chrom
          all.samples.phased.cpn.list[[samp]][[seg.name]]$pos_vec <- samp.b.lr.baf.df$pos
          all.samples.phased.cpn.list[[samp]][[seg.name]]$logr_vec <- samp.b.lr.baf.df$logR
          all.samples.phased.cpn.list[[samp]][[seg.name]]$ph_cpnA_vec <- as.numeric(phased.samp.seg.cpnA.cpnB.df[,1])
          all.samples.phased.cpn.list[[samp]][[seg.name]]$ph_cpnB_vec <- as.numeric(phased.samp.seg.cpnA.cpnB.df[,2])
          all.samples.phased.cpn.list[[samp]][[seg.name]]$mean_ph_cpnA <- mean(as.numeric(phased.samp.seg.cpnA.cpnB.df[,1]))
          all.samples.phased.cpn.list[[samp]][[seg.name]]$mean_ph_cpnB <- mean(as.numeric(phased.samp.seg.cpnA.cpnB.df[,2]))

        }#end if(nrow(samp.b.lr.baf.df) > min.num.hz.snp.phasing)
      }#end for ( seg.name in chr.start.end.string.vec ) {
    }#end for   for ( samp in sample.vec) {
    return(all.samples.phased.cpn.list)
}#end getAllSegAllSamplePhasedRawCPNList








#This is used to generate CCF values so needs to be added to the sample.to.phase.seg.list
createAllSegsPhasedBAFDataFrame <- function(asas.lr.baf.lists, sample.to.phase.seg.list, min.num.hz.snp.phasing = 5) {

  sample.vec <- names(asas.lr.baf.lists)
 
  #Now we need to make sure that the asas.lr.baf.lists and sample.to.phase.seg.list 
  if (length(names(sample.to.phase.seg.list)) != length(asas.lr.baf.lists[[1]])) {

    print("Missing or removed segments:")
    print(names(sample.to.phase.seg.list)[ ! names(asas.lr.baf.lists[[1]]) %in% names(sample.to.phase.seg.list)])
    print("Continuing creation of phased BAF dataframe")

  }

  seg.name.vec <- names(sample.to.phase.seg.list)
  all.samp.ph.BAF.df.list <- list()

  for ( seg.name in seg.name.vec ) {

    for ( samp in sample.vec ) {
      #Todo, make sure all segments are referred to by name, define this earlier on. Makes no sense to sometimes be using
      #indices and other times be using the name format, causes confusion and bugs.
      samp.BAF.df <- asas.lr.baf.lists[[samp]][[seg.name]]

      if (nrow(samp.BAF.df) < min.num.hz.snp.phasing) {

        print(seg.name)
        print("Has fewer than 5 heterozygous SNPs, not updating corresponding BAFs in table")
        next;

      }#end if(nrow(samp.BAF.df)<5) {
      
      #initialise new columns
      samp.BAF.df$samp <- samp
      samp.BAF.df$sample_phased_by_this_sample <-  sample.to.phase.seg.list[[seg.name]]$samp_chosen_to_phase_seg
      samp.BAF.df$maj_is_cpnA_logic_vec <- 0
      samp.BAF.df$maj_is_cpnB_logic_vec <- 0

      if (length(sample.to.phase.seg.list[[seg.name]]$maj_is_cpnA_idx_vec) > 1) {

        samp.BAF.df[sample.to.phase.seg.list[[seg.name]]$maj_is_cpnA_idx_vec,]$maj_is_cpnA_logic_vec <- 1 

      }#end if (length(sample.to.phase.seg.list[[seg.name]]$maj_is_cpnA_idx_vec) > 1) {

      if (length(sample.to.phase.seg.list[[seg.name]]$maj_is_cpnB_idx_vec) > 1) {

        samp.BAF.df[sample.to.phase.seg.list[[seg.name]]$maj_is_cpnB_idx_vec,]$maj_is_cpnB_logic_vec <- 1 

      }#end if (length(sample.to.phase.seg.list[[seg.name]]$maj_is_cpnB_idx_vec) > 1) {

      if (nrow(samp.BAF.df[is.na(samp.BAF.df$maj_is_cpnA_logic_vec), ,drop=FALSE]) > 0) {

        print("samp.BAF.df")
        print(samp.BAF.df)
        print("nrow(samp.BAF.df)")
        print(nrow(samp.BAF.df))
        stop("Check error in sample BAF df generation")
        
      }#end if (nrow(samp.BAF.df[is.na(samp.BAF.df$maj_is_cpnA_idx_vec), ,drop=FALSE])>0)){

      samp.BAF.df$ph_median_segmented_BAF <- min(median(samp.BAF.df[samp.BAF.df$maj_is_cpnB_logic_vec == 1 ,]$BAF),
                                           median(samp.BAF.df[samp.BAF.df$maj_is_cpnB_logic_vec == 0 ,]$BAF))

      all.samp.ph.BAF.df.list[[samp]] <- rbind(all.samp.ph.BAF.df.list[[samp]], samp.BAF.df)

    }#end for ( samp in sample.vec ) {

  }#end for ( i in 1:length(chr.start.end.string.vec) ) {

  return(all.samp.ph.BAF.df.list) 

}#end createSegPhasedBAFTable <- function(....)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          