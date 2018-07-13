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
  # print("sample.vec <- names(asas.lr.baf.lists)")
  # print(sample.vec)

  sample.to.phase.seg.list <- list()

  for ( seg.name in seg.name.vec ) {
        print(seg.name)
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
    # print("all.samp.seg.phasing.ability.list")
    # print(str(all.samp.seg.phasing.ability.list, max=1))
    a.s.sum.phas.df <- data.frame(diptest_for_AI_pval = sapply(all.samp.seg.phasing.ability.list, function(x) x$diptest_for_AI_pval),
                                  ASCAT_called_as_AI = sapply(all.samp.seg.phasing.ability.list, function(x) x$ASCAT_called_as_AI),
                                  stringsAsFactors = FALSE)
    rownames(a.s.sum.phas.df) <- sample.vec
    #Reorder by which has best separation according to th eoriginal diptest, then subset to choose only samples that pass the threshold for having allelic imbalance.
    a.s.sum.phas.df <- a.s.sum.phas.df[with(a.s.sum.phas.df, order(diptest_for_AI_pval)), ]
    # print("a.s.sum.phas.df ")
    # print(a.s.sum.phas.df)
    #Now select the samples that we will use to call for phasing.
    if (allow.dip.test.without.ASCAT.to.phase == TRUE) {

      poss.phas.s.df <- a.s.sum.phas.df[a.s.sum.phas.df$ASCAT_called_as_AI == TRUE | a.s.sum.phas.df$diptest_for_AI_pval < dip.test.p.val.threshold, ,drop = FALSE]

    } else {

      poss.phas.s.df <- a.s.sum.phas.df[a.s.sum.phas.df$ASCAT_called_as_AI == TRUE, ,drop = FALSE] 

    }
    # print("poss.phas.s.df after")
    # print(poss.phas.s.df)
    #Now assess the samples that can actualy be phased...
    #initialise the seg mixture phasing to NULL so that, as a default, it won't be used:
    seg.mixture.phasing.list <- NULL
    if (nrow(poss.phas.s.df) > 0) {

      if (nrow(poss.phas.s.df) == 1) {
        #If only one sample passes then choose this sample as the reference
        sample.to.phase.seg <- rownames(poss.phas.s.df)[1]

      } else if (nrow(poss.phas.s.df) > 1) {
        #If there are multiple alance instances called choose the top one and check mixture phasing
        sample.to.phase.seg <- rownames(poss.phas.s.df)[1]
        poss.phas.sample.vec <- rownames(poss.phas.s.df)
        print(poss.phas.sample.vec)
        all.poss.phas.samp.seg.phasing.ability.list <- getAssessmentOfEachSamplesAbilityToPhaseAllOthersForSeg(asas.lr.baf.lists = asas.lr.baf.lists,
                                                                                                    sample.vec = poss.phas.sample.vec,
                                                                                                    seg.name = seg.name,
                                                                                                    ploidy.est.vec = ploidy.est.vec,
                                                                                                    purity.est.vec = purity.est.vec,
                                                                                                    gamma.val = gamma.val,
                                                                                                    min.num.hz.snp.phasing = min.num.hz.snp.phasing)

        phas.s.df <- data.frame(diptest_for_AI_pval = sapply(all.poss.phas.samp.seg.phasing.ability.list, function(x) x$diptest_for_AI_pval),
                                      ASCAT_called_as_AI = sapply(all.poss.phas.samp.seg.phasing.ability.list, function(x) x$ASCAT_called_as_AI),
                                      mean_prop_SNPs_sep_across_samples = sapply(all.poss.phas.samp.seg.phasing.ability.list, function(x) x$mean_prop_SNPs_sep_across_samples),
                                      comb_pval_SNP_sep_across_samps = sapply(all.poss.phas.samp.seg.phasing.ability.list, function(x) x$comb_wilcox_sep_pval_across_samps),
                                      stringsAsFactors = FALSE)

        seg.mixture.phasing.list <- getMixturePhasingResultsForSegment(sample.to.phase.seg = sample.to.phase.seg,
                                                                      sample.vec = poss.phas.sample.vec,
                                                                      poss.phas.s.df = phas.s.df,
                                                                      asas.lr.baf.lists = asas.lr.baf.lists,
                                                                      all.samp.seg.phasing.ability.list = all.poss.phas.samp.seg.phasing.ability.list,
                                                                      seg.name = seg.name,
                                                                      ploidy.est.vec = ploidy.est.vec,
                                                                      purity.est.vec = purity.est.vec,
                                                                      gamma.val = gamma.val,
                                                                      comb.wilcox.snp.sep.p.val.threshold = comb.wilcox.snp.sep.p.val.threshold,
                                                                      min.prop.snp.sep.threshold = min.prop.snp.sep.threshold)

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

  #DEBUG
  # asas.lr.baf.lists = asas.SNPs.w.both.lr.baf.snp.lists
  # seg.name.vec = seg.name.vec
  # ploidy.est.vec = ploidy.est.vec
  # purity.est.vec = purity.est.vec
  # gamma.val = GAMMA_VAL
  # min.num.hz.snp.phasing = MIN_NUM_HZ_SNP_PHASING
  # min.prop.snp.sep.threshold = MIN_PROP_SNP_SEP_THRESH
  # dip.test.p.val.threshold = PVAL_DIPTEST_AI_THRESH
  # comb.wilcox.snp.sep.p.val.threshold = PVAL_COMB_WILCOX_PH_THRESH
  # allow.dip.test.without.ASCAT.to.phase = ALLOW_DIPTEST_ALONE_PH_FLAG

  source("/home/camp/watkint/working/watkint/scripts/multiregioncopynum/consensus_phasing_from_multiple_samples_functions_v2.R")
  source("/home/camp/watkint/working/watkint/scripts/multiregioncopynum/utility_MPhase_functions.R")

  sample.vec <- names(asas.lr.baf.lists)
  sample.to.phase.seg.list <- list()

  MultiResultClass <- function()
  {
    me <- list(
      samp_chosen_to_phase_seg=NA,
      samp_most_significant_by_diptest=NA,
      maj_is_cpnA_idx_vec=NA,
      maj_is_cpnB_idx_vec=NA,
      all_samp_mixture_phasing_comparison_list=NA
    )
    # Set the name for the class
    class(me) <- append(class(me),"MultiResultClass")
    return(me)
  }



  itx <- iter(seg.name.vec)
  sample.to.phase.seg.mrc.list <- foreach(seg.name = itx, .combine = c) %dopar% {
  #sample.to.phase.seg.list <- list()
  #for ( seg.name in seg.name.vec ) {
    ret.mrc <- MultiResultClass()
    print(seg.name)
    first.ss.lr.baf.df <- asas.lr.baf.lists[[sample.vec[1]]][[seg.name]]
    if (nrow(first.ss.lr.baf.df) < min.num.hz.snp.phasing) {

      #next;#If phasng not possible for this segment then move onto the next
      ret.mrc <- c()

    } else {#end if( nrow(ref.ss.lr.baf.df) < min.hz.snp.phasing)
    #For each sample how good is it at phasing all samples for this segment.
      all.samp.seg.phasing.ability.list <- getAssessmentOfEachSamplesAbilityToPhaseAllOthersForSeg(asas.lr.baf.lists = asas.lr.baf.lists,
                                                                                                  sample.vec = sample.vec,
                                                                                                  seg.name = seg.name,
                                                                                                  ploidy.est.vec = ploidy.est.vec,
                                                                                                  purity.est.vec = purity.est.vec,
                                                                                                  gamma.val = gamma.val,
                                                                                                  min.num.hz.snp.phasing = min.num.hz.snp.phasing)
      sample.to.phase.seg.list[[seg.name]]$all_samp_phasing_comparison_list <- all.samp.seg.phasing.ability.list
      a.s.sum.phas.df <- data.frame(diptest_for_AI_pval = sapply(all.samp.seg.phasing.ability.list, function(x) x$diptest_for_AI_pval),
                                    ASCAT_called_as_AI = sapply(all.samp.seg.phasing.ability.list, function(x) x$ASCAT_called_as_AI),
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
      #Now assess the samples that can actualy be phased...
      #initialise the seg mixture phasing to NULL so that, as a default, it won't be used:
      seg.mixture.phasing.list <- NULL
      if (nrow(poss.phas.s.df) > 0) {

        if (nrow(poss.phas.s.df) == 1) {
          #If only one sample passes then choose this sample as the reference
          sample.to.phase.seg <- rownames(poss.phas.s.df)[1]

        } else if (nrow(poss.phas.s.df) > 1) {
          #If there are multiple alance instances called choose the top one and check mixture phasing
          sample.to.phase.seg <- rownames(poss.phas.s.df)[1]
          poss.phas.sample.vec <- rownames(poss.phas.s.df)
          all.poss.phas.samp.seg.phasing.ability.list <- getAssessmentOfEachSamplesAbilityToPhaseAllOthersForSeg(asas.lr.baf.lists = asas.lr.baf.lists,
                                                                                                      sample.vec = poss.phas.sample.vec,
                                                                                                      seg.name = seg.name,
                                                                                                      ploidy.est.vec = ploidy.est.vec,
                                                                                                      purity.est.vec = purity.est.vec,
                                                                                                      gamma.val = gamma.val,
                                                                                                      min.num.hz.snp.phasing = min.num.hz.snp.phasing)

          phas.s.df <- data.frame(diptest_for_AI_pval = sapply(all.poss.phas.samp.seg.phasing.ability.list, function(x) x$diptest_for_AI_pval),
                                        ASCAT_called_as_AI = sapply(all.poss.phas.samp.seg.phasing.ability.list, function(x) x$ASCAT_called_as_AI),
                                        mean_prop_SNPs_sep_across_samples = sapply(all.poss.phas.samp.seg.phasing.ability.list, function(x) x$mean_prop_SNPs_sep_across_samples),
                                        comb_pval_SNP_sep_across_samps = sapply(all.poss.phas.samp.seg.phasing.ability.list, function(x) x$comb_wilcox_sep_pval_across_samps),
                                        stringsAsFactors = FALSE)

          seg.mixture.phasing.list <- getMixturePhasingResultsForSegment(sample.to.phase.seg = sample.to.phase.seg,
                                                                        sample.vec = poss.phas.sample.vec,
                                                                        poss.phas.s.df = phas.s.df,
                                                                        asas.lr.baf.lists = asas.lr.baf.lists,
                                                                        all.samp.seg.phasing.ability.list = all.poss.phas.samp.seg.phasing.ability.list,
                                                                        seg.name = seg.name,
                                                                        ploidy.est.vec = ploidy.est.vec,
                                                                        purity.est.vec = purity.est.vec,
                                                                        gamma.val = gamma.val,
                                                                        comb.wilcox.snp.sep.p.val.threshold = comb.wilcox.snp.sep.p.val.threshold,
                                                                        min.prop.snp.sep.threshold = min.prop.snp.sep.threshold)

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
          ret.mrc$samp_chosen_to_phase_seg <- sample.to.phase.seg

          if(length(sorted.all.phas.idx.vec) > 0) {

           ret.mrc$maj_is_cpnA_idx_vec <- sorted.all.phas.idx.vec[seq(1, length(all.phaseable.positions.idx.vec),  2)]

          } else {

            ret.mrc$maj_is_cpnA_idx_vec <- NA
            
          } 
          if(length(sorted.all.phas.idx.vec) > 1) {

            ret.mrc$maj_is_cpnB_idx_vec <- sorted.all.phas.idx.vec[seq(2, length(all.phaseable.positions.idx.vec),  2)]

          } else {

            ret.mrc$maj_is_cpnB_idx_vec <- NA

          }

        } else {#end if( grepl("none", sample.to.phase.seg) == TRUE)
          #As there is a phasing possible Take the appropriate sample selected
          ret.mrc$samp_chosen_to_phase_seg <- sample.to.phase.seg
          ret.mrc$maj_is_cpnA_idx_vec <- all.samp.seg.phasing.ability.list[[sample.to.phase.seg]]$maj_is_cpnA_idx_vec
          ret.mrc$maj_is_cpnB_idx_vec <- all.samp.seg.phasing.ability.list[[sample.to.phase.seg]]$maj_is_cpnB_idx_vec

        }#end if( grepl("none", sample.to.phase.seg) == FALSE)

      } else {#end if ( is.null(seg.mixture.phasing.list) ) 
        #mixture phasing has worked/been performed and is required, take the appropriate results from the mixture phasing.
        ret.mrc$samp_chosen_to_phase_seg <- "mixture"
        ret.mrc$samp_most_significant_by_diptest <- seg.mixture.phasing.list$samp_most_significant_by_diptest
        ret.mrc$maj_is_cpnA_idx_vec <- seg.mixture.phasing.list$maj_is_cpnA_idx_vec
        ret.mrc$maj_is_cpnB_idx_vec <- seg.mixture.phasing.list$maj_is_cpnB_idx_vec
        ret.mrc$all_samp_mixture_phasing_comparison_list <- seg.mixture.phasing.list$all_samp_mixture_phasing_comparison_list
      }#end if( ! is.null(seg.mixture.phasing.list) )

    }

      returnable.list <- list(ret.mrc)
      names(returnable.list) <- seg.name
      returnable.list


  }#end for ( seg.name in seg.name.vec ) {

  print("Finished paralell portion of getCombinedBAFDipTestPhasingParallel, combining results into correct format...")
  sample.to.phase.seg.list <- list()
  for ( seg.name in names(sample.to.phase.seg.mrc.list)) {

    if( length(sample.to.phase.seg.mrc.list[[seg.name]]) > 0) {

      sample.to.phase.seg.list[[seg.name]]$samp_chosen_to_phase_seg <- sample.to.phase.seg.mrc.list[[seg.name]]$samp_chosen_to_phase_seg
      sample.to.phase.seg.list[[seg.name]]$samp_most_significant_by_diptest <- sample.to.phase.seg.mrc.list[[seg.name]]$samp_most_significant_by_diptest
      sample.to.phase.seg.list[[seg.name]]$maj_is_cpnA_idx_vec <- sample.to.phase.seg.mrc.list[[seg.name]]$maj_is_cpnA_idx_vec
      sample.to.phase.seg.list[[seg.name]]$maj_is_cpnB_idx_vec <- sample.to.phase.seg.mrc.list[[seg.name]]$maj_is_cpnB_idx_vec
      sample.to.phase.seg.list[[seg.name]]$all_samp_mixture_phasing_comparison_list <- sample.to.phase.seg.mrc.list[[seg.name]]$all_samp_mixture_phasing_comparison_list

    }

  }

  return(sample.to.phase.seg.list)

}#end getCombinedBAFDipTestPhasing


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

    oth.ss.cpnA.cpnB.df <- calculateSNPsRawCopyNumberAB(logr.vec = oth.ss.lr.baf.df$logR,
                                                        baf.vec = oth.ss.lr.baf.df$BAF,
                                                        tumor.ploidy = ploidy.est.vec[other.samp],
                                                        tumor.purity = purity.est.vec[other.samp],
                                                        gamma.val = gamma.val)

    #old code had some checks here to make sure that these have the same number of rows, but if everything holds they should remain the same anyway...
    other.samp.wilcox.test <-  wilcoxTestWithErrorHandling( c(oth.ss.cpnA.cpnB.df[ref.ss.maj.is.cpnA.idx.vec,  1], oth.ss.cpnA.cpnB.df[ref.ss.maj.is.cpnB.idx.vec, 2]), 
                                                            c(oth.ss.cpnA.cpnB.df[ref.ss.maj.is.cpnA.idx.vec,  2], oth.ss.cpnA.cpnB.df[ref.ss.maj.is.cpnB.idx.vec, 1]), paired = TRUE)

    oth.ss.snps.same.as.ref <- length(which(c(oth.ss.cpnA.cpnB.df[ref.ss.maj.is.cpnA.idx.vec, 1], oth.ss.cpnA.cpnB.df[ref.ss.maj.is.cpnB.idx.vec, 2]) - 
                                            c(oth.ss.cpnA.cpnB.df[ref.ss.maj.is.cpnA.idx.vec, 2], oth.ss.cpnA.cpnB.df[ref.ss.maj.is.cpnB.idx.vec, 1]) > 0)) / nrow(oth.ss.cpnA.cpnB.df)
    oth.ss.snps.diff.to.ref <- length(which(c(oth.ss.cpnA.cpnB.df[ref.ss.maj.is.cpnA.idx.vec, 1], oth.ss.cpnA.cpnB.df[ref.ss.maj.is.cpnB.idx.vec, 2]) - 
                                            c(oth.ss.cpnA.cpnB.df[ref.ss.maj.is.cpnA.idx.vec, 2], oth.ss.cpnA.cpnB.df[ref.ss.maj.is.cpnB.idx.vec, 1]) < 0)) / nrow(oth.ss.cpnA.cpnB.df)
    #Updating the other sample segment phasing df
    other.ss.phasing.df[other.samp,]$separation_wilcox_pval  <- other.samp.wilcox.test$p.val
    other.ss.phasing.df[other.samp,]$prop_separated <- max(oth.ss.snps.same.as.ref, oth.ss.snps.diff.to.ref)
    other.ss.phasing.df[other.samp,]$concordant_w_ref <- ifelse(oth.ss.snps.same.as.ref >= oth.ss.snps.diff.to.ref, 1, -1)

  }#end for ( other.samp in sample.vec[ ! sample%in%sample.vec] )

  sample.comparison.list[[samp]]$comb_wilcox_sep_pval_across_samps <- pchisq((sum(log(other.ss.phasing.df$separation_wilcox_pval)) * -2),
                                                                            df = length(other.ss.phasing.df$separation_wilcox_pval) * 2,
                                                                            lower.tail = FALSE)
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
      # print(sample.vec)
      # print(samp)
      # print(seg.name)
      # print(str(asas.lr.baf.lists[[samp]][[seg.name]], max=1))
      ref.ss.lr.baf.df <- asas.lr.baf.lists[[samp]][[seg.name]]
      # print("ref.ss.lr.baf.df$BAFsegmented")
      # print(head(ref.ss.lr.baf.df$BAFsegmented))
      # print("ref.ss.lr.baf.df$BAF")
      # print(head(ref.ss.lr.baf.df$BAF))
      # print(nrow(ref.ss.lr.baf.df))
      # print(head(ref.ss.lr.baf.df$logr))
      if ( nrow(ref.ss.lr.baf.df) < min.num.hz.snp.phasing) {

        next#Might actally want this to be a break as they should all have too few snps as it is joint segmentation with the same SNPs...

      }#end if( nrow(ref.ss.lr.baf.df) < min.hz.snp.phasing)
      #renaming the dataframe rows for the other samples, now we are trying a new reference sample
      rownames(other.ss.phasing.df) <- sample.vec[ ! sample.vec %in% samp ]

      #Get phasing, reference sample segment cpn A cpn B
      ref.ss.cpnA.cpnB.df <- calculateSNPsRawCopyNumberAB(logr.vec = ref.ss.lr.baf.df$logR,
                                                          baf.vec = ref.ss.lr.baf.df$BAF,
                                                          tumor.ploidy = ploidy.est.vec[samp],
                                                          tumor.purity = purity.est.vec[samp],
                                                          gamma.val = gamma.val)
      ref.ss.max.cpn.vec <- apply(ref.ss.cpnA.cpnB.df, 1, max)
      ref.ss.is.max.cpn.cpnA.cpnB.logical.df <- ifelse(ref.ss.cpnA.cpnB.df >= ref.ss.max.cpn.vec, 1, 0)
      ref.ss.maj.is.cpnA.idx.vec <- which(ref.ss.is.max.cpn.cpnA.cpnB.logical.df[, 1] == 1)
      ref.ss.maj.is.cpnB.idx.vec <- which(ref.ss.is.max.cpn.cpnA.cpnB.logical.df[, 2] == 1)

      sample.comparison.list[[samp]]$diptest_for_AI_pval <- as.numeric(dip.test(ref.ss.lr.baf.df$BAF)$p.val)
      # print(head(ref.ss.lr.baf.df))
      # print("ifelse(ref.ss.lr.baf.df$BAFsegmented[1] != 0.5, TRUE, FALSE)")
      # print(ifelse(ref.ss.lr.baf.df$BAFsegmented[1] != 0.5, TRUE, FALSE))
      sample.comparison.list[[samp]]$ASCAT_called_as_AI <- ifelse(median(ref.ss.lr.baf.df$BAFsegmented, na.rm=TRUE) %in% 0.5, FALSE, TRUE)#Changed to 2 due to idxing but... as only uysed for plotting ignoring.
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

  return(NULL)
  print("Applying mixture phasing...")
  phas.ss.lr.baf.df <- asas.lr.baf.lists[[sample.to.phase.seg]][[seg.name]]
  phas.ss.cpnA.cpnB.df <- calculateSNPsRawCopyNumberAB(logr.vec=phas.ss.lr.baf.df$logR, baf.vec=phas.ss.lr.baf.df$BAF, tumor.ploidy=ploidy.est.vec[sample.to.phase.seg], tumor.purity=purity.est.vec[sample.to.phase.seg], gamma.val=gamma.val)
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
      ss.cpnA.cpnB.df <- calculateSNPsRawCopyNumberAB(logr.vec = ss.lr.baf.df$logR,
                                                      baf.vec = ss.lr.baf.df$BAF,
                                                      tumor.ploidy = ploidy.est.vec[samp],
                                                      tumor.purity = purity.est.vec[samp],
                                                      gamma.val = gamma.val)
      ss.SNP.cpnA <- ss.cpnA.cpnB.df[SNP.idx, current.maj.cpn.col.idx]
      ss.SNP.cpnB <- ss.cpnA.cpnB.df[SNP.idx, current.min.cpn.col.idx]

      s.cpnAB.match.ref.factor <- oth.s.cnpA.cpnB.match.ref.vec[samp]
      #sample snp sep score
      s.SNP.sep.score <- (ss.SNP.cpnA - ss.SNP.cpnB) * s.cpnAB.match.ref.factor
      all.s.SNP.sep.scores <- c(all.s.SNP.sep.scores, s.SNP.sep.score)

    }#end for ( samp in all.passing.samples.vec ) {

    if (mean(all.s.SNP.sep.scores) > 0) {

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
                                                                                                              asas.lr.baf.lists = asas.lr.baf.lists,
                                                                                                              samples.to.test.vec = sample.vec,
                                                                                                              ref.ss.maj.is.cpnA.idx.vec = new.maj.is.cpnA.idx.vec,
                                                                                                              ref.ss.maj.is.cpnB.idx.vec = new.maj.is.cpnB.idx.vec,
                                                                                                              ploidy.est.vec = ploidy.est.vec,
                                                                                                              purity.est.vec = purity.est.vec,
                                                                                                              gamma.val = gamma.val)
  print("returning from mixture phasing.")
  return(output.list)
}


