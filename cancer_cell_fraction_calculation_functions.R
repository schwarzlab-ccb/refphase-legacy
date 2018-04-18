#Note needs tTestWithErrorHandling which is in CN_heterogeneity_utility_functions.R


generateAllSegsRoundedCPNSameOrDifferentAcrossSamplesLists <- function(all.pw.cn.het.comp.lists, sample.vec, seg.name.vec) {
      #todo change these to preallocated dataframes if possible.. not ideally suited as theses are vectors of different lengths.
      seg.cpn.across.samps.list$same_cpnA_cpnB_seg_name_vec <- c()     
      seg.cpn.across.samps.list$diff_cpnA_cpnB_seg_name_vec <- c()
      seg.cpn.across.samps.list$same_cpnA_seg_name_vec <- c()
      seg.cpn.across.samps.list$diff_cnpA_seg_name_vec <- c()
      seg.cpn.across.samps.list$same_cpnB_seg_name_vec <- c()
      seg.cpn.across.samps.list$diff_cpnB_seg_name_vec <- c()

      for ( seg.name in seg.name.vec ) {
            #samp
            ss.mean.int.cpnA <- c()
            ss.mean.int.cpnB <- c()

            for ( samp in samp.vec ) {

                  ss.mean.int.cpnA <- c(ss.mean.int.cpnA , round(all.pw.cn.het.comp.lists$sample_comparison_list[[samp]][[seg.name]]$mean_ph_cpnA))
                  ss.mean.int.cpnB <- c(ss.mean.int.cpnB , round(all.pw.cn.het.comp.lists$sample_comparison_list[[samp]][[seg.name]]$mean_ph_cpnB))

            }#end for ( samp in samp.vec ) {    

            num.uniq.cpnA <- length(unique(ss.int.mean.cpnA))
            num.uniq.cpnB <- length(unique(ss.int.mean.cpnB))

            #Check if total cpn the same across all samples for seg
            if (num.uniq.cpnA == 1 & num.uniq.cpnB == 1) {

                  seg.cpn.across.samps.list$same_cpnA_cpnB_seg_name_vec <- c(seg.cpn.across.samps.list$same_cpnA_cpnB_seg_name_vec, seg.name)

            } else {

                  seg.cpn.across.samps.list$diff_cpnA_cpnB_seg_name_vec <- c(seg.cpn.across.samps.list$diff_cpnA_cpnB_seg_name_vec, seg.name)

            }
            #Check if cpnA the same across alls samples for seg
            if (num.uniq.cpnA == 1) {

                  seg.cpn.across.samps.list$same_cpnA_seg_name_vec <- c(seg.cpn.across.samps.list$same_cpnA_seg_name_vec, seg.name)

            } else {

                  seg.cpn.across.samps.list$diff_cnpA_seg_name_vec <- c(seg.cpn.across.samps.list$diff_cnpA_seg_name_vec, seg.name)

            }
            #Check if cpnB the same across alls samples for seg
            if (num.uniq.cpnB == 1 ) {

                  seg.cpn.across.samps.list$same_cpnB_seg_name_vec <- c(seg.cpn.across.samps.list$same_cpnB_seg_name_vec, seg.name)

            } else {

                  seg.cpn.across.samps.list$diff_cpnB_seg_name_vec <- c(seg.cpn.across.samps.list$diff_cpnB_seg_name_vec, seg.name)

            }
           
      }#end for ( seg.name in seg.names ) {

      return(seg.cpn.across.samps.list)

}#end generateAllSegsVectorListRoundedSameOrDifferentCPN <- function(...)



##
##Params:
##
##Returns:
##
generateAllSamplesByChromIntCPNUnderlyingRawValuesLists <- function(all.pw.cn.het.comp.lists, sample.vec, seg.names) {

      samps.chr.int.cpn.raw.vals.list <- list()
      #s.samps.obs.exp.cpn.df <- data.frame()
      #But really the datastructure organisation overall needs to change if we speed.
      for ( samp in sample.vec ) {

            for ( seg.name in seg.name.vec ) {

                  samp.seg.int.cpnA.val <- round(all.pw.cn.het.comp.lists[[samp]][[seg.name]]$mean_ph_cpnA)
                  samp.seg.int.cpnB.val <- round(all.pw.cn.het.comp.lists[[samp]][[seg.name]]$mean_ph_cpnB)
                  chr <- unlist(strsplit(seg.name, split=":"))[1]

                  sample.chr.int.cpn.raw.vals.list[[samp]][[as.character(samp.seg.int.cpnA.val)]][[chr]]$all_seg_cpnA_vec <- c(
                        sample.chr.int.cpn.raw.vals.list[[samp]][[as.character(samp.seg.int.cpnA.val)]][[chr]]$all_seg_cpnA_vec, 
                        all.pw.cn.het.comp.lists$sample_comparison_list[[samp]][[seg.name]]$ph_cpnA_vec)

                  sample.chr.int.cpn.raw.vals.list[[samp]][[as.character(samp.seg.int.cpnB.val)]][[chr]]$all_seg_cpnB_vec <- c(
                        sample.chr.int.cpn.raw.vals.list[[samp]][[as.character(samp.seg.int.cpnB.val)]][[chr]]$all_seg_cpnB_vec, 
                        all.pw.cn.het.comp.lists$sample_comparison_list[[samp]][[seg.name]]$ph_cpnB_vec)
            }#end for ( seg.name in seg.name.vec ) {

      }#end  for ( samp in sample.vec ) {

      return(samps.chr.int.cpn.raw.vals.list)
}#end generateAllSamplesChromsIntCpnAvgRawValuesList <- function()





##
##Params:
##
##Returns:
##
#Todo no sure this is well written...
calculateWhetherSegCPNDiffFromIntegerExpectedAndAllMatchingSampCPN <- function(samp.pvals.diff.from.exp.list,
                                                                              all.pw.cn.het.comp.lists,
                                                                              s.samps.int.ob.xp.cpn.df,
                                                                              seg.name,
                                                                              chrs.for.ccf.comp,
                                                                              samps.int.cpn.raw.cpn.by.chr.list,
                                                                              cpnA.or.cpnB){

      sample.vec <- rownames(s.samps.obs.exp.cpn.df)
      all.match.cpn.vals.vec <- c()
      chrs.w.same.cpn.vec <- c()
      #Todo, perhapsrefactor into a returnable dataframe...
      for ( samp in sample.vec ) {
            
            for ( chr in chrs.for.ccf.comp ) {

                  s.chr.cpn.vec <- samps.int.cpn.raw.cpn.by.chr.list[[samp]][[chr]][[s.samps.int.ob.xp.cpn.df[samp, paste0("round_int_", cpnA.or.cpnB)]]][[paste0("all_seg_", cpnA.or.cpnB, "_vec")]]
                  if ( length(s.chr.cpn.vec) != 0 ) {

                        chrs.w.same.cpn.vec <- c(chrs.w.same.cpn.vec, as.character(chr))
                        all.match.cpn.vals.vec <- c(all.match.cpn.vals.vec, s.chr.cpn.vec)

                  }#end if ( length(s.chr.cpn.vec) != 0 )

            }#endfor ( chr in chrs.for.ccf.comp ) {

            seg.raw.cpn.vals <- all.pw.cn.het.comp.lists$sample_comparison_list[[samp]][[seg.name]][[paste0("ph_",cpnA.or.cpnB,"_vec"]]

            t.test.raw.seg.cpn.vs.all.matching.cpn.raw <- tTestWithErrorHandling(seg.raw.cpn.vals, all.match.cpn.vals.vec)
            t.test.raw.seg.vs.cpn.int <- tTestWithErrorHandling(seg.raw.cpn.vals, s.samps.int.ob.xp.cpn.df[samp, paste0("round_int_", cpnA.or.cpnB)])
            median.seg.raw.cpn.greater.than.matched.cpn <- median(seg.raw.cpn.vals) > median(all.match.cpn.vals.vec)
            median.seg.raw.cpn.greater.seg.int.cpn   <- median(seg.raw.cpn.vals) > s.samps.int.ob.xp.cpn.df[samp, paste0("round_int_", cpnA.or.cpnB)]

            samp.pvals.diff.from.exp.list[[samp]][[paste0("round_int_", cpnA.or.cpnB)]] <- s.samps.int.ob.xp.cpn.df[samp, paste0("round_int_", cpnA.or.cpnB)]
            samp.pvals.diff.from.exp.list[[samp]][[paste0(cpnA.or.cpnB, "_pval_vs_oth_cpn_match_chrs_raw_cpn")]] <- t.test.raw.seg.cpn.vs.all.matching.cpn.raw$p.value
            samp.pvals.diff.from.exp.list[[samp]][[paste0(cpnA.or.cpnB, "_pval_vs_round_int_cpn")]] <- t.test.raw.seg.vs.cpn.int$p.value
            samp.pvals.diff.from.exp.list[[samp]][[paste0(cpnA.or.cpnB, "_num_chrs_w_match_int_cpn")]]  <- chrs.w.same.cpn.vec
            samp.pvals.diff.from.exp.list[[samp]][[paste0("median_", cpnA.or.cpnB, "_greater_than_match_raw_cpn")]] <- median.seg.raw.cpn.greater.than.matched.cpn
            samp.pvals.diff.from.exp.list[[samp]][[paste0("median_", cpnA.or.cpnB, "_greater_than_match_int_cpn")]] <- median.seg.raw.cpn.greater.seg.int.cpn

      }#end for ( samp in sample.vec ) {

      return(samp.pvals.diff.from.exp.list)

}





##
##Params:
##
##Returns:
##
getSegObsExpCPNBySampleDF <- function(all.pw.cn.het.comp.lists, sample.vec, seg.name){

      num.samples <- length(sample.vec)
      s.samps.obs.exp.cpn.df <- data.frame(round_int_cpnA = numeric(num.samples),
                                          round_int_cpnB = numeric(num.samples),
                                          mean_obs_cpnA = numeric(num.samples),
                                          mean_obs_cpnB = numeric(num.samples),
                                          exp_cpnA = numeric(num.samples),
                                          exp_cpnB = numeric(num.samples),
                                          stringsAsFactors = FALSE)
      rownames(seg.samps.obs.exp.cpn.df) <- sample.vec    

      for ( samp in sample.vec ) {

            s.samps.obs.exp.cpn.df[samp, "round_int_cpnA"] <- round(all.pw.cn.het.comp.lists$sample_comparison_list[[samp]][[seg.name]]$mean_ph_cpnA)
            s.samps.obs.exp.cpn.df[samp, "round_int_cpnB"] <- round(all.pw.cn.het.comp.lists$sample_comparison_list[[samp]][[seg.name]]$mean_ph_cpnB)
            #Calculate observed
            if ( all.pw.cn.het.comp.lists$sample_comparison_list[[samp]][[seg.name]]$mean_ph_cpnA < 0 ) {

                  s.samps.obs.exp.cpn.df[samp, "obs_mean_cpnA"]  <- 0

            } else {

                  s.samps.obs.exp.cpn.df[samp, "obs_mean_cpnA"]  <- all.pw.cn.het.comp.lists$sample_comparison_list[[samp]][[seg.name]]$mean_ph_cpnA

            }
            if (  all.pw.cn.het.comp.lists$sample_comparison_list[[samp]][[seg.name]]$mean_ph_cpnB < 0 ) {

                  s.samps.obs.exp.cpn.df[samp, "obs_mean_cpnB"] <-  0

            } else {

                  s.samps.obs.exp.cpn.df[samp, "obs_mean_cpnB"] <- all.pw.cn.het.comp.lists$sample_comparison_list[[samp]][[seg.name]]$mean_ph_cpnB)

            }
            #Calculate expected
            s.samps.obs.exp.cpn.df[samp, "exp_cpnA"] <- max(0, s.samps.obs.exp.cpn.df[samp, "round_int_cpnA"])

            s.samps.obs.exp.cpn.df[samp, "exp_cpnB"] <- max(0, s.samps.obs.exp.cpn.df[samp, "round_int_cpnB"])


      }#End for ( samp in sample.vec ) {

      return(s.samps.obs.exp.cpn.df)

}#end getSegObsExpCPNBySampleDF <- function(...



getCCFractionNeeded <- function(obs.val, low.exp.val, high.exp.val) {

  if (obs.val > high.exp.val) {

    if (obs.val / high.exp.val != Inf) {

      return(c(0, obs.val / high.exp.val))

    }#end if (obs.val / high.exp.val != Inf) {
    
    if(obs.val / high.exp.val == Inf) {

      return(c(1, 0))

    }#end if(obs.val / high.exp.val == Inf) {

  }#end if (obs.val > high.exp.val)
  
  if (obs.val < low.exp.val) {

    return(c(1, 0))

  }#if (obs.val < low.exp.val) {

  # because my brain isn't working, i'm going to do a brute force 2digit CCF calc
  CCF.to.use <- 0.01

  est.val <- ((CCF.to.use * low.exp.val)) + ((1 - CCF.to.use) * high.exp.val) 
  est.val <- abs(obs.val - est.val)
  
  for (CCF.val in seq(0.0, 1.5, by=.01)) {

    c.est.val <- ((CCF.val * low.exp.val)) + ((1 - CCF.val) * high.exp.val) 
    c.est.val <- abs(obs.val - c.est.val)

    if (c.est.val < est.val) {

      est.val <- c.est.val
      CCF.to.use <- CCF.val

    }#end if(c.est.val < est.val) {

  }#end for (CCF.val in seq(0.0, 1.5, by=.01)) {
  
  return(c(CCF.to.use, (1 - CCF.to.use)))

}#end getCCFractionNeeded <- function(...)



calcExpectedLogR <- function(exp_intA, exp_intB, tumor.ploidy, tumor.purity) {

  expectedLogR <- log2((2 * (1 - tumor.purity) + tumor.purity * (exp_intA + exp_intB)) / (2 * (1 - tumor.purity) + (tumor.purity * tumor.ploidy)))

  return(expectedLogR)
}

calcExpectedBAF <- function(exp_intA, exp_intB, tumor.purity) {

  expectedBAF <- (1 - tumor.purity + (tumor.purity * exp_intB)) / (2 - 2 * tumor.purity + tumor.purity * (exp_intA + exp_intB))

  return(expectedBAF)
}


#Can be refactored further...
determineCCFClonalityOfSegment <- function(s.samps.int.ob.xp.cpn.df,
                    samp.vec,
                    seg.name,
                    all.samp.ph.BAF.df.list,
                    ploidy.est.vec,
                    purity.est.vec,
                    min.CCF.to.consider.val) {
  #all the same so can just take the first value... not sure it needs to be like this, could just be unique
  all.clonal.flag <- TRUE
  ss.round.mean.obs.cpn <- round(s.samps.int.ob.xp.cpn.df[, paste0("mean_obs_", cpnX.string)])[1]
  non.clonal.CPN.samps.vec <- c()
  # first check whether all regions are clonal
  for ( samp in samp.vec ) {
        #sample seg phased cpnX vector of values...
        ss.ph.cpnX.vec <- all.pw.cn.het.comp.lists$sample_comparison_list[[samp]][[seg.name]][[paste0("ph_", cpnX.string, "_vec")]]
        samp.ph.cpnX.vs.samp.round.mean.obs.cpnX.test <- tTestWithErrorHandling(ss.ph.cpnX.vec, mu = ss.round.mean.obs.cpn, conf.level = confidence.interval.theshold)

        if (samp.ph.cpnX.vs.samp.round.mean.obs.cpnX.test$p.val < 1 - confidence.interval.theshold) {

              all.clonal.flag <- FALSE
              non.clonal.samps.vec <- c(non.clonal.samps.vec, samp)

        }#end if (samp.ph.cpnX.vs.samp.round.mean.obs.cpnX.test$p.val < 1 - confidence.interval.theshold) {

  }#for ( samp in samp.vec ) {
  all.samp.clonal.t.test <- TRUE
  non.clonal.BAF.samps.vec <- c()
  seg.chr.start.end.vec <- unlist(strsplit(seg.name, split=":"))
  chr <- as.numeric(seg.chr.start.end.vec[1])
  start <- as.numeric(seg.chr.start.end.vec[2])
  end <- as.numeric(seg.chr.start.end.vec[3])
  #Todo check with Nicky that we want to use this...
  for ( samp in sample.vec ) {

        #Get BAF values from the dataframe
        samp.ph.baf.df <- all.samp.ph.BAF.df.list[[samp]]
        s.samp.ph.baf.df <- samp.ph.baf.df[samp.ph.baf.df$chrom == chr & samp.ph.baf.df$pos >= start & samp.ph.baf.df$pos <= end, ,drop=FALSE]

        samp.maj.is.cpnA.cpnB.BAF.vs.exp.BAF.test <- tTestWithErrorHandling(c(s.samp.ph.baf.df[s.samp.ph.baf.df$maj_is_cpnA_logic_vec==1,"BAF"],
                                              s.samp.ph.baf.df[s.samp.ph.baf.df$maj_is_cpnB_logic_vec==1,"BAF"]]),
                                            mu = calcExpectedBAF(s.samps.int.ob.xp.cpn.df[samp, "exp_cpnA"],
                                                                 s.samps.int.ob.xp.cpn.df[samp, "exp_cpnB"],
                                                                 ploidy.est.vec[samp],
                                                                 purity.est.vec[samp]))

        if (samp.maj.is.cpnA.cpnB.BAF.vs.exp.BAF.test$p.val < (1 - confidence.interval.theshold)) {

              clonal.samp.t.test <- FALSE
              non.clonal.BAF.samps.vec <- c(non.clonal.BAF.samps.vec, samp)

        }#if (samp.maj.is.cpnA.cpnB.BAF.vs.exp.BAF.test$p.val < (1 - confidence.interval.theshold)) {

  }#end for ( samp in sample.vec ) {
  #Finally check what the CCF is of each region if we assume subclonality
  max.CCF <- 0
  poss.sc.cpnX.low <- min(unique(floor(s.samps.int.ob.xp.cpn.df[, paste0("mean_obs_", cpnX.string)])))
  poss.sc.cpnX.high <- max(unique(ceiling(s.samps.int.ob.xp.cpn.df[, paste0("mean_obs_", cpnX.string)])))

  for (samp in sample.vec) {

        ss.mean.osbs.ph.cpnX  <- s.samps.int.ob.xp.cpn.df[samp, paste0("mean_obs_", cpnX.string)]
        ss.fracs.to.use.vec <- getCCFractionNeeded(ss.mean.osbs.ph.cpnX, poss.sc.cpnX.low, poss.sc.cpnX.high)
        ss.min.CCF <- min(ss.fracs.to.use.vec)

        if (ss.min.CCF > max.CCF) {

              max.CCF <- ss.min.CCF

        }#end if (seg.min.CCF > maxCCF) {

  }#end or (samp in sample.vec) {
  #Now determine whether to explore subclonality:
  if (all.clonal.flag == TRUE | all.samp.clonal.t.test == TRUE | max.CCF < min.CCF.to.consider.val) {

        all.clonal.flag <- TRUE

  }#end if (all.clonal.flag == TRUE | all.samp.clonal.t.test == TRUE | max.CCF min.CCF.to.consider.val){

  return(all.clonal.flag)

}#determineCCFClonalityOfSegment <- function(...)





##
##Params:
##
##Returns:
##
updateSegCCFValList <- function(seg.CCF.val.list,
                sample.vec,
                s.samps.obs.exp.cpn.df,
                all.pw.cn.het.comp.lists,
                cpnX.string,
                ss.poss.sc.cpnX.low,
                ss.poss.sc.cpnX.high,
                use.preset.mean.frac.cpnX.low.high
                ) {

  seg.CCF.val.list[[seg.name]][[paste0("poss_sc_", cpnX.string, "_low"]] <- ss.poss.sc.cpnX.low
  seg.CCF.val.list[[seg.name]][[paste0("poss_sc_", cpnX.string, "_high"]] <- ss.poss.sc.cpnX.high

  #Get and save mean fractions based on obs

  for ( samp in sample.vec ) {

    #Generate mean values:
    ss.mean.cpnX.obs <- s.samps.obs.exp.cpn.df[samp, paste0("mean_obs_", cpnX.string)]
    #Clean the mean...
    ss.mean.cpnX.obs <- ifelse(ss.mean.cpnX.obs < 0, 0, ss.mean.cpnX.obs)

    if (use.preset.mean.frac.cpnX.low.high == TRUE) {

      ss.mean.frac.cpnX.low <- 0
      ss.mean.frac.cpnX.high <- 1

    } else {

      ss.fracs.to.use.vec <- getCCFractionNeeded(ss.mean.cpnX.obs, ss.poss.sc.cpnX.low, ss.poss.sc.cpnX.high)
      ss.mean.frac.cpnX.low <- ss.fracs.to.use.vec[1]
      ss.mean.frac.cpnX.low <- ss.fracs.to.use.vec[2]

    }#end if (use.preset.mean.frac.cpnX.low.high == TRUE) {

    #Save mean values:
    seg.CCF.val.list[[seg.name]][[samp]][[paste0("mean_obs_", cpnX.string)]] <- ss.mean.cpnX.obs
    seg.CCF.val.list[[seg.name]][[samp]][[paste0("mean_frac_", cpnX.string, "_low")]] <- ss.mean.frac.cpnX.low
    seg.CCF.val.list[[seg.name]][[samp]][[paste0("mean_frac_", cpnX.string, "_high")]] <- ss.mean.frac.cpnX.high

    #Generate raw values:
    ss.ph.cpnX.vals.vec <- all.pw.cn.het.comp.lists$sample_comparison_list[[samp]][[seg.name]][[paste0("ph_", cpnX.string, "_vec")]]
    ss.ph.cpnX.vals.vec <- ifelse(ss.ph.cpnX.vals.vec < 0, 0, ss.ph.cpnX.vals.vec)
    ss.ph.cpnX.vals.vec <- ifelse(ss.ph.cpnX.vals.vec >= max.cpn.val, max.cpn.val, ss.ph.cpnX.vals.vec.no.neg)
    ss.fracs.to.use.vec <- sapply(ss.ph.cpnX.vals.vec, getCCFractionNeeded, ss.mean.frac.cpnX.low,  ss.mean.frac.cpnX.high) 
    ss.fracs.to.use.vec <-  ss.fracs.to.use.vec[ss.fracs.to.use.vec == "Inf"] <- 1
    ss.frac.poss.sc.cpnX.low.vec <- apply(ss.fracs.to.use.no.inf.vec, 2, min)
    ss.frac.poss.sc.cpnX.high.vec <- apply(ss.fracs.to.use.no.inf.vec, 2, max)

    #Save raw values:
    seg.CCF.val.list[[seg.name]][[samp]][[paste0("clean_obs_", cpnX.string, "_vec")]]   <- ss.ph.cpnX.vals.vec
    seg.CCF.val.list[[seg.name]][[samp]][[paste0("frac_", cpnX.string, "_low_vec")]]  <- ss.frac.poss.sc.cpnX.low.vec
    seg.CCF.val.list[[seg.name]][[samp]][[paste0("frac_", cpnX.string, "_high_vec")]]  <- ss.frac.poss.sc.cpnX.high.vec

  }#end for ( samp in sample.vec ) {

  return(seg.CCF.val.list)

}#end calculateWhetherSegCPNDiffFromIntegerExpectedAndAllMatchingSampCPN <- function(




##
##Params:
##
##Returns:
##
calculateTwoValueDiffFlag <- function(s.samps.int.ob.xp.cpn.df, seg.cpnA.cnpB.diff.exp.list, cpnX.string, min.CCF.to.consider.val, confidence.interval.threshold) {

  two.value.diff.flag <- TRUE
  mean.cpnX.vec <- s.samps.int.ob.xp.cpn.df[, paste0("mean_obs_", cpnX.string)]
  low.int.mean.cpnX <- min(floor(mean.cpnX.vec))
  top.int.mean.cpnX <- max(ceiling(mean.cpnX.vec))

  if (abs(min(top.int.mean.cpnX - mean.cpnX.vec)) >= abs(max(low.int.mean.cpnX - mean.cpnX.vec))) {

      ss.poss.sc.cpnX.low <- low.int.mean.cpnX
      ss.poss.sc.cpnX.high <- top.int.mean.cpnX

      #sis this choice of ss.poss.sc.cpnX.high, too conservative?
      samps.poss.more.than.sc.cnpX.high.vec  <- rownames(s.samps.int.ob.xp.cpn.df[s.samps.int.ob.xp.cpn.df[, paste0("mean_obs_", cpnX.string)] >= ss.poss.sc.cpnX.high, ])

      for ( samp in samps.poss.more.than.sc.cnpX.high.vec ) {

        samp.max.CCF <- abs(ss.poss.sc.cpnX.high - samps.int.ob.xp.cpn.df[samp, paste0("mean_obs_", cpnX.string)])
        raw.diff.from.int.cpnX.pval <- seg.cpnA.cnpB.diff.exp.list[[samp]][[paste0(cpnX.string, "_pval_vs_round_int_cpn")]]
        diff.from.int.flag <- raw.diff.from.int.cpnX.pval < (1 - confidence.interval.threshold)

        diff.from.other.chrs.flag <- TRUE
        if (length(seg.cpnA.cnpB.diff.exp.list[[samp]][[paste0(cpnX.string, "_chrs_w_match_int_cpn_vec")]]) > min.num.other.chrs) {

          diff.from.other.chrs.flag <- seg.cpnA.cnpB.diff.exp.list[[samp]][[paste0(cpnX.string, "_pval_vs_oth_cpn_match_chrs_raw_cpn")]] < (1 - confidence.interval.threshold)

        }#end if (length(seg.cpnA.cnpB.diff.exp.list[[samp]][[paste0(cpnX.string, "_chrs_w_match_int_cpn_vec")]]) > min.num.other.chrs) {

        if (diff.from.int.flag == TRUE & diff.from.other.chrs.flag == TRUE & samp.max.CCF >= min.CCF.to.consider.val) {

          two.value.diff.flag <- FALSE

        }

      }#end for ( samp in sample.vec ) {

  }#end if(abs(min(top.int.mean.cpnX - mean.cpnX.vec)) >= abs(max(low.int.mean.cpnX - mean.cpnX.vec))) {

  if (abs(min(top.int.mean.cpnX - mean.cpnX.vec)) < abs(max(low.int.mean.cpnX - mean.cpnX.vec))) {

    ss.poss.sc.cpnX.low  <- round(mean.cpnX.vec)[1]
    ss.poss.sc.cpnX.high  <- max(unique(ceiling(s.samps.int.ob.xp.cpn.df[, paste0("mean_obs_", cpnX.string)])))

    # is this choice of nA1, too conservative
    samps.poss.less.than.sc.cnpX.low.vec <- rownames(s.samps.int.ob.xp.cpn.df[s.samps.int.ob.xp.cpn.df[, paste0("mean_obs_", cpnX.string)] <= ss.poss.sc.cpnX.low, ])
    for (samp in samps.poss.less.than.sc.cnpX.low.vec) {
      
      samp.max.CCF      <-  abs(ss.poss.sc.cpnX.low - samps.int.ob.xp.cpn.df[samp, paste0("mean_obs_", cpnX.string)])
      raw.diff.from.int.cpnX.pval <- seg.cpnA.cnpB.diff.exp.list[[samp]][[paste0(cpnX.string, "_pval_vs_round_int_cpn")]]
      diff.from.other.chrs.flag <- TRUE
      if (length(seg.cpnA.cnpB.diff.exp.list[[samp]][[paste0(cpnX.string, "_chrs_w_match_int_cpn_vec")]]) > min.num.other.chrs) { 

        diff.from.other.chrs.flag <- seg.cpnA.cnpB.diff.exp.list[[samp]][[paste0(cpnX.string, "_pval_vs_oth_cpn_match_chrs_raw_cpn")]] < (1 - confidence.interval.threshold)

      }#end if (length(seg.cpnA.cnpB.diff.exp.list[[samp]][[paste0(cpnX.string, "_chrs_w_match_int_cpn_vec")]]) > min.num.other.chrs) { 
      
      if (diff.from.int.flag == TRUE & diff.from.other.chrs.flag == TRUE & samp.max.CCF >= min.CCF.to.consider.val) {

        two.value.diff.flag <- FALSE

      }#end if (diff.from.int.flag == TRUE & diff.from.other.chrs.flag == TRUE & samp.max.CCF >= min.CCF.to.consider.val) {

    }#end for (samp in samps.poss.less.than.sc.cnpX.low.vec) {

  }#end if (abs(min(top.int.mean.cpnX - mean.cpnX.vec)) < abs(max(low.int.mean.cpnX - mean.cpnX.vec))) {

  return(two.value.diff.flag)

}









getIndividualCCFestimates_withPvals <- function(seg.name,
                                                sample.vec,
                                                all.pw.cn.het.comp.lists,
                                                seg.r.cpn.same.diff.across.samps.list,
                                                all.samp.ph.BAF.df.list,
                                                max.cpn.val = 5,
                                                min.CCF.to.consider.val = 0.15,
                                                chrs.can.use.for.ccf.comp = c(1:22),
                                                confidence.interval.threshold = 0.9999,
                                                min.num.other.chrs = 3) {

      if (is.na(all.pw.cn.het.comp.lists$sample_comparison_list[[samp]][[seg.name]]$ph_cpnA_vec)) {

            return(NA)
      }
      chr <- strsplit(seg.name, split=":")[1]
      s.samps.int.ob.xp.cpn.df <- getSegObsExpCPNBySampleDF(all.pw.cn.het.comp.lists, sample.vec, seg.name)
      chrs.for.ccf.comp.vec <- chrs.can.use.for.ccf.comp[ ! chrs.can.use.for.ccf.comp%in% chr] #Going to be used to assess CCFs across the each sample

      seg.CCF.val.list <- list()

      for ( cpnX.string in c("cpnA", "cpnB") ) {

            seg.cpnA.cnpB.diff.exp.list <- calculateWhetherSegCPNDiffFromIntegerExpectedAndAllMatchingSampCPN(samp.pvals.diff.from.exp.list=seg.cpnA.cnpB.diff.exp.list,
                                                                                                              all.pw.cn.het.comp.lists,
                                                                                                              s.samps.int.ob.xp.cpn.df,
                                                                                                              seg.name,
                                                                                                              chrs.for.ccf.comp.vec,
                                                                                                              samps.int.cpn.raw.cpn.by.chr.list,
                                                                                                              cpnA.or.cpnB=cpnX.string)
            #setup what we're going to use for all the potential conditions here and tests here so we can pass it to the 
            s.num.unq.mean.obs.cpn <- length(unique(s.samps.int.ob.xp.cpn.df[, paste0("mean_obs_", cpnX.string)]
            s.floor.obs.vec <- s.samps.int.ob.xp.cpn.df[, paste0("mean_obs_", cpnX.string)])
            s.ceil.obs.vec <- s.samps.int.ob.xp.cpn.df[, paste0("mean_obs_", cpnX.string)])

            #Not sure that this ever called.... check
            if (length(unique(s.samps.int.ob.xp.cpn.df[, paste0("mean_obs_", cpnX.string)])) == 1) {

                  #add these things to the segCCFvalues
                  seg.CCF.val.list <- updateSegCCFValList(seg.CCF.val.list = seg.CCF.val.list, sample.vec = sample.vec, s.samps.obs.exp.cpn.df = s.samps.int.ob.xp.cpn.df, 
                                                          all.pw.cn.het.comp.lists = all.pw.cn.het.comp.lists, cpnX.string = cpnX.string,
                                                          ###These parameters vary:
                                                          ss.poss.sc.cpnX.low = 1,#Should really be ploidy n this case.
                                                          ss.poss.sc.cpnX.high = unique(s.samps.obs.exp.cpn.df[, paste0("mean_obs_", cpnX.string)]),
                                                          use.preset.mean.frac.cpnX.low.high = TRUE) 

            }#end if (length(unique(s.samps.obs.exp.cpn.df[, paste0("mean_obs_",cpnX)])) == 1) {

            # let's check whether the ceiling and floor cpn lead to only two options:
            if (length(unique(c(floor(s.samps.int.ob.xp.cpn.df[, paste0("mean_obs_", cpnX.string)]),
                              ceiling(s.samps.int.ob.xp.cpn.df[, paste0("mean_obs_", cpnX.string)])))) == 2) {

                  # check whether only one state when rounding
                  if (length(unique(round(s.samps.int.ob.xp.cpn.df[, paste0("mean_obs_", cpnX.string)]))) == 1) {

                        all.clonal.flag <- determineCCFClonalityOfSegment(s.samps.int.ob.xp.cpn.df = s.samps.int.ob.xp.cpn.df,
                                            sample.vec = sample.vec,
                                            seg.name = seg.name,
                                            all.samp.ph.BAF.df.list= all.samp.ph.BAF.df.list,
                                            ploidy.est.vec = ploidy.est.vec,
                                            purity.est.vec = purity.est.vec,
                                            min.CCF.to.consider.val = min.CCF.to.consider.val)

                        if (all.clonal.flag == TRUE) {

                              #add these things to the segCCFvalues
                              seg.CCF.val.list <- updateSegCCFValList(seg.CCF.val.list = seg.CCF.val.list, sample.vec = sample.vec, s.samps.obs.exp.cpn.df = s.samps.int.ob.xp.cpn.df, 
                                                                  all.pw.cn.het.comp.lists = all.pw.cn.het.comp.lists, cpnX.string = cpnX.string,
                                                                  ###These parameters vary:
                                                                  ss.poss.sc.cpnX.low =  round(s.samps.obs.exp.cpn.df[ 1 , paste0("mean_obs_", cpnX.string)]),
                                                                  ss.poss.sc.cpnX.high = round(s.samps.obs.exp.cpn.df[ 1 , paste0("mean_obs_", cpnX.string)]),
                                                                  use.preset.mean.frac.cpnX.low.high = FALSE) 

                        } else {#if not all clonal: 

                              #add these things to the segCCFvalues
                              seg.CCF.val.list <- updateSegCCFValList(seg.CCF.val.list = seg.CCF.val.list, sample.vec = sample.vec, s.samps.obs.exp.cpn.df = s.samps.int.ob.xp.cpn.df, 
                                                                  all.pw.cn.het.comp.lists = all.pw.cn.het.comp.lists, cpnX.string = cpnX.string,
                                                                  ###These parameters vary:
                                                                  ss.poss.sc.cpnX.low =  min(unique(floor(s.samps.obs.exp.cpn.df[ , paste0("mean_obs_", cpnX.string)]))),
                                                                  ss.poss.sc.cpnX.high = max(unique(ceiling(s.samps.obs.exp.cpn.df[ , paste0("mean_obs_", cpnX.string)]))),
                                                                  use.preset.mean.frac.cpnX.low.high = FALSE) 

                        }#end if (all.clonal.flag == TRUE) {

                  }#end if (length(unique(round(s.samps.int.ob.xp.cpn.df[, paste0("mean_obs_", cpnX.string)]))) == 1) {

                  if (length(unique(round(s.samps.int.ob.xp.cpn.df[, paste0("mean_obs_", cpnX.string)]))) != 1) {


                        #add these things to the segCCFvalues
                        seg.CCF.val.list <- updateSegCCFValList(seg.CCF.val.list = seg.CCF.val.list, sample.vec = sample.vec, s.samps.obs.exp.cpn.df = s.samps.int.ob.xp.cpn.df, 
                                                            all.pw.cn.het.comp.lists = all.pw.cn.het.comp.lists, cpnX.string = cpnX.string,
                                                            ###These parameters vary:
                                                            ss.poss.sc.cpnX.low =  min(unique(floor(s.samps.obs.exp.cpn.df[ , paste0("mean_obs_", cpnX.string)]))),
                                                            ss.poss.sc.cpnX.high = max(unique(ceiling(s.samps.obs.exp.cpn.df[ , paste0("mean_obs_", cpnX.string)]))),
                                                            use.preset.mean.frac.cpnX.low.high = FALSE) 


                  }#end if (length(unique(round(s.samps.int.ob.xp.cpn.df[, paste0("mean_obs_", cpnX.string)]))) != 1) {

            }#end if (length(unique(c(floor(s.samps.int.ob.xp.cpn.df[, paste0("mean_obs_", cpnX.string))]),
                  #ceiling(s.samps.int.ob.xp.cpn.df[, paste0("mean_obs_", cpnX.string)]))) == 2) {

            # let's check whether the ceiling and floor cpn lead to only two options:
            if (length(unique(c(floor(s.samps.int.ob.xp.cpn.df[, paste0("mean_obs_", cpnX.string)]),
                              ceiling(s.samps.int.ob.xp.cpn.df[, paste0("mean_obs_", cpnX.string)])))) > 2) {


                  #check whether there is ony one state when rounding
                  if (length(unique(round(s.samps.int.ob.xp.cpn.df[, paste0("mean_obs_", cpnX.string)]))) == 1) {

                        all.clonal.flag <- determineCCFClonalityOfSegment(s.samps.int.ob.xp.cpn.df = s.samps.int.ob.xp.cpn.df,
                                            sample.vec = sample.vec,
                                            seg.name = seg.name,
                                            all.samp.ph.BAF.df.list= all.samp.ph.BAF.df.list,
                                            ploidy.est.vec = ploidy.est.vec,
                                            purity.est.vec = purity.est.vec,
                                            min.CCF.to.consider.val = min.CCF.to.consider.val)

                        if (all.clonal.flag == TRUE) {

                              seg.CCF.val.list <- updateSegCCFValList(seg.CCF.val.list = seg.CCF.val.list, sample.vec = sample.vec, s.samps.obs.exp.cpn.df = s.samps.int.ob.xp.cpn.df, 
                                                                  all.pw.cn.het.comp.lists = all.pw.cn.het.comp.lists, cpnX.string = cpnX.string,
                                                                  ###These parameters vary:
                                                                  ss.poss.sc.cpnX.low =  round(s.samps.obs.exp.cpn.df[ 1 , paste0("mean_obs_", cpnX.string)]),
                                                                  ss.poss.sc.cpnX.high = round(s.samps.obs.exp.cpn.df[ 1 , paste0("mean_obs_", cpnX.string)]),
                                                                  use.preset.mean.frac.cpnX.low.high = FALSE) 

                        } else {#end if (all.clonal.flag == TRUE) {

                              two.value.diff.flag <- calculateTwoValueDiffFlag(s.samps.int.ob.xp.cpn.df,
                                                                              seg.cpnA.cnpB.diff.exp.list,
                                                                              cpnX.string,
                                                                              min.CCF.to.consider.val,
                                                                              confidence.interval.threshold)


                              if (two.value.diff.flag == TRUE) {

                                    mean.cpnX.vec <- s.samps.int.ob.xp.cpn.df[, paste0("mean_obs_", cpnX.string)]
                                    #Do we twist up or do we twist down
                                    if (abs(min(max(ceiling(mean.cpnX.vec - mean.cpnX.vec)) < abs(max(min(floor(mean.cpnX.vec)) - mean.cpnX.vec))){
                                          #take the ceiling
                                          seg.CCF.val.list <- updateSegCCFValList(seg.CCF.val.list = seg.CCF.val.list, sample.vec = sample.vec, s.samps.obs.exp.cpn.df = s.samps.int.ob.xp.cpn.df, 
                                                                             all.pw.cn.het.comp.lists = all.pw.cn.het.comp.lists, cpnX.string = cpnX.string,
                                                                             ###These parameters vary:
                                                                             ss.poss.sc.cpnX.low =  round(s.samps.obs.exp.cpn.df[ 1 , paste0("mean_obs_", cpnX.string)]),
                                                                             ss.poss.sc.cpnX.high = max(unique(ceiling(s.samps.obs.exp.cpn.df[ , paste0("mean_obs_", cpnX.string)]))),
                                                                             use.preset.mean.frac.cpnX.low.high = FALSE)  


                                    }#end if (abs(min(max(ceiling(mean.cpnX.vec - mean.cpnX.vec)) < abs(max(min(floor(mean.cpnX.vec)) - mean.cpnX.vec))){

                                    if (abs(min(max(ceiling(mean.cpnX.vec - mean.cpnX.vec)) >= abs(max(min(floor(mean.cpnX.vec)) - mean.cpnX.vec))) {      

                                          # take the floor
                                          seg.CCF.val.list <- updateSegCCFValList(seg.CCF.val.list = seg.CCF.val.list, sample.vec = sample.vec, s.samps.obs.exp.cpn.df = s.samps.int.ob.xp.cpn.df, 
                                                                             all.pw.cn.het.comp.lists = all.pw.cn.het.comp.lists, cpnX.string = cpnX.string,
                                                                             ###These parameters vary:
                                                                             ss.poss.sc.cpnX.low =  min(unique(floor(s.samps.obs.exp.cpn.df[ 1 , paste0("mean_obs_", cpnX.string)]))),
                                                                             ss.poss.sc.cpnX.high = max(unique(round(s.samps.obs.exp.cpn.df[ , paste0("mean_obs_", cpnX.string)]))),
                                                                             use.preset.mean.frac.cpnX.low.high = FALSE)  



                                    }#end if(abs(min(max(ceiling(mean.cpnX.vec - mean.cpnX.vec)) >= abs(max(min(floor(mean.cpnX.vec)) - mean.cpnX.vec))) { 


                              } else {#end: if (two.value.diff.flag == TRUE)# start of: two.value.diff.flag==FALSE

                                    mean.cpnX.vec <- s.samps.int.ob.xp.cpn.df[, paste0("mean_obs_", cpnX.string)]
                                    # do we twist up or down
                                    if ((abs(min(max(ceiling(mean.cpnX.vec - mean.cpnX.vec)) < abs(max(min(floor(mean.cpnX.vec)) - mean.cpnX.vec)))) {
                                          #take the ceiling
                                          seg.CCF.val.list <- updateSegCCFValList(seg.CCF.val.list = seg.CCF.val.list, sample.vec = sample.vec, s.samps.obs.exp.cpn.df = s.samps.int.ob.xp.cpn.df, 
                                                                             all.pw.cn.het.comp.lists = all.pw.cn.het.comp.lists, cpnX.string = cpnX.string,
                                                                             ###These parameters vary:
                                                                             ss.poss.sc.cpnX.low =  min(unique(floor(s.samps.obs.exp.cpn.df[ , paste0("mean_obs_", cpnX.string)]))),
                                                                             ss.poss.sc.cpnX.high = max(unique(ceiling(s.samps.obs.exp.cpn.df[ , paste0("mean_obs_", cpnX.string)]))),
                                                                             use.preset.mean.frac.cpnX.low.high = FALSE)  

                                    }
                                    if (abs(min(max(ceiling(mean.cpnX.vec)) - mean.cpnX.vec)) >= abs(max(min(floor(mean.cpnX.vec)) - mean.cpnX.vec))) {
                                          #take the floor #TODO, check with Nicky but these values are the same as the < than, BUG?
                                          seg.CCF.val.list <- updateSegCCFValList(seg.CCF.val.list = seg.CCF.val.list, sample.vec = sample.vec, s.samps.obs.exp.cpn.df = s.samps.int.ob.xp.cpn.df, 
                                                                             all.pw.cn.het.comp.lists = all.pw.cn.het.comp.lists, cpnX.string = cpnX.string,
                                                                             ###These parameters vary:
                                                                             ss.poss.sc.cpnX.low =  min(unique(floor(s.samps.obs.exp.cpn.df[ , paste0("mean_obs_", cpnX.string)]))),
                                                                             ss.poss.sc.cpnX.high = max(unique(ceiling(s.samps.obs.exp.cpn.df[ , paste0("mean_obs_", cpnX.string)]))),
                                                                             use.preset.mean.frac.cpnX.low.high = FALSE)  


                                    }#end if (abs(min(max(ceiling(mean.cpnX.vec)) - mean.cpnX.vec)) >= abs(max(min(floor(mean.cpnX.vec)) - mean.cpnX.vec))) {

                              }#end  two.value.diff.flag==FALSE

                        }#end  else {#end if (all.clonal.flag == FALSE) {

                  }#end if (length(unique(round(s.samps.int.ob.xp.cpn.df[, paste0("mean_obs_", cpnX.string)]))) == 1){

                  #check whether there are only two states when rounding
                  if (length(unique(round(s.samps.int.ob.xp.cpn.df[, paste0("mean_obs_", cpnX.string)]))) == 2) {

                        seg.CCF.val.list <- updateSegCCFValList(seg.CCF.val.list = seg.CCF.val.list, sample.vec = sample.vec, s.samps.obs.exp.cpn.df = s.samps.int.ob.xp.cpn.df, 
                                                           all.pw.cn.het.comp.lists = all.pw.cn.het.comp.lists, cpnX.string = cpnX.string,
                                                           ###These parameters vary:
                                                           ss.poss.sc.cpnX.low =  min(round(s.samps.obs.exp.cpn.df[, paste0("mean_obs_", cpnX.string)])),
                                                           ss.poss.sc.cpnX.high = max(round(s.samps.obs.exp.cpn.df[ , paste0("mean_obs_", cpnX.string)])),
                                                           use.preset.mean.frac.cpnX.low.high = FALSE)      


                  }#end if (length(unique(round(s.samps.int.ob.xp.cpn.df[, paste0("mean_obs_", cpnX.string)]))) == 2){

                  #check whether there are more than two states when rounding
                  if (length(unique(round(s.samps.int.ob.xp.cpn.df[, paste0("mean_obs_", cpnX.string)]))) > 2) {

                        seg.CCF.val.list <- updateSegCCFValList(seg.CCF.val.list = seg.CCF.val.list, sample.vec = sample.vec, s.samps.obs.exp.cpn.df = s.samps.int.ob.xp.cpn.df, 
                                                           all.pw.cn.het.comp.lists = all.pw.cn.het.comp.lists, cpnX.string = cpnX.string,
                                                           ###These parameters vary:
                                                           ss.poss.sc.cpnX.low =  min(unique(floor(s.samps.obs.exp.cpn.df[, paste0("mean_obs_", cpnX.string)] + min.CCF.to.consider))),
                                                           ss.poss.sc.cpnX.high = max(unique(ceiling(s.samps.obs.exp.cpn.df[ , paste0("mean_obs_", cpnX.string)] + min.CCF.to.consider))),
                                                           use.preset.mean.frac.cpnX.low.high = FALSE)      

                  }#end if (length(unique(round(s.samps.int.ob.xp.cpn.df[, paste0("mean_obs_", cpnX.string)]))) == 2){

            }#end if (length(unique(c(floor(s.samps.int.ob.xp.cpn.df[, paste0("mean_obs_", cpnX.string))]),
                  #ceiling(s.samps.int.ob.xp.cpn.df[, paste0("mean_obs_", cpnX.string)]))) > 2) {
    
      }#end for ( cpnX.string in c("cpnA", "cpnB") )

}#end getIndividualCCFestimates_withPvals <- function(...)

#Todo break this up into multiple functions.
getSegNewCCFEstimate <- function( all.pw.cn.het.comp.lists
                                 ,ploidy.est.vec
                                 ,purity.est.vec
                                 ,highest.cpn.trusted = 10
                                 ,pval.threshold = 0.01
                                 ,pval.to.separate = 0.001
                                 ,gain.pval.to.use = 'gainPval'
                                 ,loss.pval.to.use = 'lossPval'
                                 ,min.CCF.to.consider = 0.15
                                 ,confidence.interval = 0.99999
) {   

      # need to obtain an ancestral state for each copy number 
      # (initially, simply do this by voting (will fail for subclonal GD events)).

      sample.vec <- names(ploidy.est.vec)
      #old code got them from all.pw.cn.het.comp.lists, shouldn't bve nece
      seg.name.vec <- names(all.pw.cn.het.comp.lists$sample_comparison_list[[sample.vec[1]]])

      #Check to make sure that all segments are present in all samples... why is this necessary
      #Todo, check why this bug is ocurring
      for ( samp in sample.vec ) {

            seg.name.vec <- intersect(seg.name.vec, names((all.pw.cn.het.comp.lists$sample_comparison_list[[samp]]))

      }#end for ( samp in sample.vec ) {     
      #todo change these to preallocated dataframes if possible      
      seg.r.cpn.same.diff.across.samps.list <- generateAllSegsRoundedCPNSameOrDifferentAcrossSamplesLists(all.pw.cn.het.comp.lists, sample.vec, seg.name.vec)
      #next assess for each chromosome and region, what are the different raw-values
      samps.int.cpn.raw.cpn.by.chr.list <- generateAllSamplesByChromIntCPNUnderlyingRawValuesLists(all.pw.cn.het.comp.lists, sample.vec, seg.name.vec)

      max.CCF.vals.list <- list()
      max.cpn <- highest.cpn.trusted







}#end getSegNewCCFEstimate           <- function(