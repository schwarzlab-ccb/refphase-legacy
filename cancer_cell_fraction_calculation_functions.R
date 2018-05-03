#Note needs tTestWithErrorHandling which is in CN_heterogeneity_utility_functions.R


generateAllSegsRoundedCPNSameOrDifferentAcrossSamplesLists <- function(all.pw.cn.het.comp.lists, sample.vec, seg.name.vec) {
      #todo change these to preallocated dataframes if possible.. not ideally suited as theses are vectors of different lengths.
      seg.cpn.across.samps.list <- list()
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

            for ( samp in sample.vec ) {

                  ss.mean.int.cpnA <- c(ss.mean.int.cpnA , round(all.pw.cn.het.comp.lists[[samp]][[seg.name]]$mean_ph_cpnA))
                  ss.mean.int.cpnB <- c(ss.mean.int.cpnB , round(all.pw.cn.het.comp.lists[[samp]][[seg.name]]$mean_ph_cpnB))

            }#end for ( samp in samp.vec ) {    

            num.uniq.cpnA <- length(unique(ss.mean.int.cpnA))
            num.uniq.cpnB <- length(unique(ss.mean.int.cpnB))

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
            if (num.uniq.cpnB == 1) {

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
generateAllSamplesByChromIntCPNUnderlyingRawValuesLists <- function(all.pw.cn.het.comp.lists, sample.vec, seg.name.vec, cpnX.vec=c("cpnA", "cpnB"), reduced.seg.flag=FALSE) {
      sample.chr.int.cpn.raw.vals.list <- list()
      #s.samps.obs.exp.cpn.df <- data.frame()
      #But really the datastructure organisation overall needs to change if we speed.
      for ( samp in sample.vec ) {

            for ( seg.name in seg.name.vec ) {

                  if (reduced.seg.flag == FALSE) {
                    if (is.na(all.pw.cn.het.comp.lists[[1]][[seg.name]]$ph_cpnA_vec[1]) | is.na(all.pw.cn.het.comp.lists[[1]][[seg.name]]$ph_cpnB_vec[1])) {

                          next
                    }#end if (is.na(all.pw.cn.het.comp.lists[[1]][[seg.name]]$ph_cpnA_vec[1]) | is.na(all.pw.cn.het.comp.lists[[1]][[seg.name]]$ph_cpnB_vec[1])) {
                  }

                  for ( cpnX in cpnX.vec ) {

                    if(reduced.seg.flag == FALSE) {
                      samp.seg.int.cpnX.val <- round(all.pw.cn.het.comp.lists[[samp]][[seg.name]][[paste0("mean_ph_",cpnX)]])

                      chr <- unlist(strsplit(seg.name, split=":"))[1]
                      sample.chr.int.cpn.raw.vals.list[[samp]][[as.character(samp.seg.int.cpnX.val)]][[chr]][[paste0("all_seg_",cpnX,"_vec")]] <- c(
                            sample.chr.int.cpn.raw.vals.list[[samp]][[as.character(samp.seg.int.cpnX.val)]][[chr]][[paste0("all_seg_",cpnX,"_vec")]], 
                            all.pw.cn.het.comp.lists[[samp]][[seg.name]][[paste0("ph_",cpnX,"_vec")]])

                    } else {

                      #print(str(all.pw.cn.het.comp.lists, max=1))
                      #print(str(all.pw.cn.het.comp.lists[[1]], max=1))
                      samp.seg.int.cpnX.val <- round(mean(all.pw.cn.het.comp.lists[[seg.name]][[paste0(cpnX,"_orig_ph_obs_",cpnX,"_vec")]]))
                      #print("print(samp.seg.int.cpnX.val)")
                      #print(samp.seg.int.cpnX.val)
                      chr <- unlist(strsplit(seg.name, split=":"))[1]
                      sample.chr.int.cpn.raw.vals.list[[samp]][[as.character(samp.seg.int.cpnX.val)]][[chr]][[paste0("all_seg_",cpnX,"_vec")]] <- c(
                            sample.chr.int.cpn.raw.vals.list[[samp]][[as.character(samp.seg.int.cpnX.val)]][[chr]][[paste0("all_seg_",cpnX,"_vec")]], 
                            all.pw.cn.het.comp.lists[[seg.name]][[paste0(cpnX,"_orig_ph_obs_",cpnX,"_vec")]])


                    }
                  }#end  for ( cpnX in cpnX.vec ) {

            }#end for ( seg.name in seg.name.vec ) {

      }#end  for ( samp in sample.vec ) {

      return(sample.chr.int.cpn.raw.vals.list)

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
                                                                              cpnA.or.cpnB,
                                                                              reduced.seg.flag=FALSE){

      sample.vec <- rownames(s.samps.int.ob.xp.cpn.df)
      #Todo, perhapsrefactor into a returnable dataframe...

      #Todo refacor so doesn't need the reduced seg flag.
      if (reduced.seg.flag == FALSE) {
        if (is.na(all.pw.cn.het.comp.lists[[1]][[seg.name]]$ph_cpnA_vec[1]) | is.na(all.pw.cn.het.comp.lists[[1]][[seg.name]]$ph_cpnB_vec[1])) {

            next
        }
      }

      for ( samp in sample.vec ) {

            chrs.w.same.cpn.vec <- c()
            all.match.cpn.vals.vec <- c()


            for ( chr in chrs.for.ccf.comp ) {

                  s.chr.cpn.vec <- samps.int.cpn.raw.cpn.by.chr.list[[samp]][[as.character(s.samps.int.ob.xp.cpn.df[samp, paste0("round_int_", cpnA.or.cpnB)])]][[as.character(chr)]][[paste0("all_seg_", cpnA.or.cpnB, "_vec")]]

                  if ( length(s.chr.cpn.vec) != 0 ) {

                        chrs.w.same.cpn.vec <- c(chrs.w.same.cpn.vec, as.character(chr))
                        all.match.cpn.vals.vec <- c(all.match.cpn.vals.vec, s.chr.cpn.vec)

                  }#end if ( length(s.chr.cpn.vec) != 0 )

            }#endfor ( chr in chrs.for.ccf.comp ) {

            if  (reduced.seg.flag == FALSE) {

              seg.raw.cpn.vals <- all.pw.cn.het.comp.lists[[samp]][[seg.name]][[paste0("ph_",cpnA.or.cpnB,"_vec")]]

            } else {

              seg.raw.cpn.vals <- all.pw.cn.het.comp.lists[[seg.name]][[paste0(cpnA.or.cpnB,"_orig_ph_obs_",cpnA.or.cpnB,"_vec")]]

            }

            t.test.raw.seg.cpn.vs.all.matching.cpn.raw <- tTestWithErrorHandling(seg.raw.cpn.vals, all.match.cpn.vals.vec)
            t.test.raw.seg.vs.cpn.int <- tTestWithErrorHandling(seg.raw.cpn.vals, mu = s.samps.int.ob.xp.cpn.df[samp, paste0("round_int_", cpnA.or.cpnB)])
            median.seg.raw.cpn.greater.than.matched.cpn <- median(seg.raw.cpn.vals) > median(all.match.cpn.vals.vec)
            median.seg.raw.cpn.greater.seg.int.cpn   <- median(seg.raw.cpn.vals) > s.samps.int.ob.xp.cpn.df[samp, paste0("round_int_", cpnA.or.cpnB)]

            samp.pvals.diff.from.exp.list[[samp]][[paste0("round_int_", cpnA.or.cpnB)]] <- s.samps.int.ob.xp.cpn.df[samp, paste0("round_int_", cpnA.or.cpnB)]
            samp.pvals.diff.from.exp.list[[samp]][[paste0(cpnA.or.cpnB, "_pval_vs_oth_cpn_match_chrs_raw_cpn")]] <- t.test.raw.seg.cpn.vs.all.matching.cpn.raw$p.value
            samp.pvals.diff.from.exp.list[[samp]][[paste0(cpnA.or.cpnB, "_pval_vs_round_int_cpn")]] <- t.test.raw.seg.vs.cpn.int$p.value
            #print("chrs.w.same.cpn.vec")
            #print(chrs.w.same.cpn.vec)
            #print(str(samp.pvals.diff.from.exp.list, max=1))
            samp.pvals.diff.from.exp.list[[samp]][[paste0(cpnA.or.cpnB, "_num_chrs_w_match_int_cpn_vec")]]  <- chrs.w.same.cpn.vec
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
getSegObsExpCPNBySampleDF <- function(all.pw.cn.het.comp.lists, sample.vec, seg.name, cpnX.vec=c("cpnA", "cpnB"), reduced.seg.flag=FALSE){

      num.samples <- length(sample.vec)
      s.samps.obs.exp.cpn.df <- data.frame(round_int_cpnA = numeric(num.samples),
                                          round_int_cpnB = numeric(num.samples),
                                          mean_obs_cpnA = numeric(num.samples),
                                          mean_obs_cpnB = numeric(num.samples),
                                          exp_cpnA = numeric(num.samples),
                                          exp_cpnB = numeric(num.samples),
                                          stringsAsFactors = FALSE)
      rownames(s.samps.obs.exp.cpn.df) <- sample.vec    
      for ( cpnX in cpnX.vec ) {

        for ( samp in sample.vec ) {
              
              #print(str(all.pw.cn.het.comp.lists[[seg.name]][[paste0(cpnX,"_orig_ph_obs_",cpnX,"_vec")]], max=1))

              if (reduced.seg.flag == FALSE) {

                mean.cpnX <- all.pw.cn.het.comp.lists[[samp]][[seg.name]][[paste0("mean_ph_",cpnX)]]
                s.samps.obs.exp.cpn.df[samp, paste0("round_int_",cpnX)] <- round(mean.cpnX)


              } else {
                mean.cpnX <- mean(all.pw.cn.het.comp.lists[[seg.name]][[paste0(cpnX,"_orig_ph_obs_",cpnX,"_vec")]])
                s.samps.obs.exp.cpn.df[samp, paste0("round_int_",cpnX)] <- round(mean.cpnX)
              }
              #print("mean.cpnX")
              #print(mean.cpnX)
              #Calculate observed
              if ( mean.cpnX < 0 ) {

                    s.samps.obs.exp.cpn.df[samp, paste0("mean_obs_",cpnX)]  <- 0

              } else {

                    s.samps.obs.exp.cpn.df[samp, paste0("mean_obs_",cpnX)]  <- mean.cpnX

              }
              #Calculate expected
              s.samps.obs.exp.cpn.df[samp, paste0("exp_",cpnX)] <- max(0, s.samps.obs.exp.cpn.df[samp, paste0("round_int_",cpnX)])

        }#End for ( samp in sample.vec ) {
          
      }#end for ( cpnX in cpnX.vec ) {
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



calcExpectedBAF <- function(exp_intA, exp_intB, tumor.purity) {

  expectedBAF <- (1 - tumor.purity + (tumor.purity * exp_intB)) / (2 - 2 * tumor.purity + tumor.purity * (exp_intA + exp_intB))

  return(expectedBAF)
  
}#end calcExpectedBAF


#Can be refactored further...
determineCCFClonalityOfSegment <- function(s.samps.int.ob.xp.cpn.df,
                    seg.name,
                    all.samp.ph.BAF.df.list,
                    purity.est.vec,
                    cpnX.string,
                    min.CCF.to.consider.val,
                    confidence.interval.threshold) {
  #all the same so can just take the first value... not sure it needs to be like this, could just be unique
  sample.vec <- names(purity.est.vec)
  all.clonal.flag <- TRUE
  ss.round.mean.obs.cpn <- round(s.samps.int.ob.xp.cpn.df[, paste0("mean_obs_", cpnX.string)])[1]
  #Todo, remove non.clonal.CPN.samps.vec as not used..
  non.clonal.CPN.samps.vec <- c()
  # first check whether all regions are clonal
  for ( samp in sample.vec ) {
        #sample seg phased cpnX vector of values...
        ss.ph.cpnX.vec <- all.pw.cn.het.comp.lists[[samp]][[seg.name]][[paste0("ph_", cpnX.string, "_vec")]]
        samp.ph.cpnX.vs.samp.round.mean.obs.cpnX.test <- tTestWithErrorHandling(ss.ph.cpnX.vec, mu = ss.round.mean.obs.cpn, conf.level = confidence.interval.threshold)

        if (samp.ph.cpnX.vs.samp.round.mean.obs.cpnX.test$p.val < 1 - confidence.interval.threshold) {

              all.clonal.flag <- FALSE
              non.clonal.CPN.samps.vec <- c(non.clonal.CPN.samps.vec, samp)

        }#end if (samp.ph.cpnX.vs.samp.round.mean.obs.cpnX.test$p.val < 1 - confidence.interval.threshold) {

  }#for ( samp in sample.vec ) {
  all.samp.clonal.t.test <- TRUE
  #Todo, remove non.clonal.BAF.samps.vec as not used..
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
                                              s.samp.ph.baf.df[s.samp.ph.baf.df$maj_is_cpnB_logic_vec==1,"BAF",]),
                                            mu = calcExpectedBAF(s.samps.int.ob.xp.cpn.df[samp, "exp_cpnA"],
                                                                 s.samps.int.ob.xp.cpn.df[samp, "exp_cpnB"],
                                                                 purity.est.vec[samp]))

        if (samp.maj.is.cpnA.cpnB.BAF.vs.exp.BAF.test$p.val < (1 - confidence.interval.threshold)) {

              clonal.samp.t.test <- FALSE
              non.clonal.BAF.samps.vec <- c(non.clonal.BAF.samps.vec, samp)

        }#if (samp.maj.is.cpnA.cpnB.BAF.vs.exp.BAF.test$p.val < (1 - confidence.interval.threshold)) {

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
updateSegCCFValList <- function(seg.name,
                seg.CCF.val.list,
                sample.vec,
                s.samps.int.ob.xp.cpn.df,
                all.pw.cn.het.comp.lists,
                cpnX.string,
                max.cpn.val,
                ss.poss.sc.cpnX.low,
                ss.poss.sc.cpnX.high,
                use.preset.mean.frac.cpnX.low.high=FALSE
                ) {

  seg.CCF.val.list[[seg.name]][[samp]][[paste0("poss_sc_", cpnX.string, "_low")]] <- ss.poss.sc.cpnX.low
  seg.CCF.val.list[[seg.name]][[samp]][[paste0("poss_sc_", cpnX.string, "_high")]] <- ss.poss.sc.cpnX.high


  #Get and save mean fractions based on obs

  for ( samp in sample.vec ) {
    #Generate mean values:
    ss.mean.cpnX.obs <- s.samps.int.ob.xp.cpn.df[samp, paste0("mean_obs_", cpnX.string)]

    #Clean the mean...
    ss.mean.cpnX.obs <- ifelse(ss.mean.cpnX.obs < 0, 0, ss.mean.cpnX.obs)

    if (use.preset.mean.frac.cpnX.low.high == TRUE) {

      ss.mean.frac.cpnX.low <- 0
      ss.mean.frac.cpnX.high <- 1

    } else {

      ss.fracs.to.use.vec <- getCCFractionNeeded(ss.mean.cpnX.obs, ss.poss.sc.cpnX.low, ss.poss.sc.cpnX.high)
      ss.mean.frac.cpnX.low <- ss.fracs.to.use.vec[1]
      ss.mean.frac.cpnX.high <- ss.fracs.to.use.vec[2]

    }#end if (use.preset.mean.frac.cpnX.low.high == TRUE) {
    #Save mean values:
    #Generate raw values:
    ss.ph.cpnX.vals.vec <- all.pw.cn.het.comp.lists[[samp]][[seg.name]][[paste0("ph_", cpnX.string, "_vec")]]
    ss.ph.cpnX.vals.vec <- ifelse(ss.ph.cpnX.vals.vec < 0, 0, ss.ph.cpnX.vals.vec)
    ss.ph.cpnX.vals.vec <- ifelse(ss.ph.cpnX.vals.vec >= max.cpn.val, max.cpn.val, ss.ph.cpnX.vals.vec)
    ss.fracs.to.use.vec <- sapply(ss.ph.cpnX.vals.vec, getCCFractionNeeded, ss.poss.sc.cpnX.low, ss.poss.sc.cpnX.high)
    ss.fracs.to.use.vec[ss.fracs.to.use.vec == "Inf"] <- 1
    ss.frac.poss.sc.cpnX.low.vec <- apply(ss.fracs.to.use.vec, 2, min)
    ss.frac.poss.sc.cpnX.high.vec <- apply(ss.fracs.to.use.vec, 2, max)
    num.vals <- length(ss.ph.cpnX.vals.vec)
    #Save raw values:
    seg.CCF.val.list[[seg.name]][[samp]][[paste0("mean_obs_", cpnX.string)]] <- ss.mean.cpnX.obs
    #Note we only need one of these....
    seg.CCF.val.list[[seg.name]][[samp]][[paste0("mean_frac_", cpnX.string, "_low")]] <- rep(ss.mean.frac.cpnX.low, length(num.vals))
    seg.CCF.val.list[[seg.name]][[samp]][[paste0("mean_frac_", cpnX.string, "_high")]] <- rep(ss.mean.frac.cpnX.high, length(num.vals))

    seg.CCF.val.list[[seg.name]][[samp]][[paste0(cpnX.string,"_orig_ph_obs_cpnA_vec")]]   <- all.pw.cn.het.comp.lists[[samp]][[seg.name]][[paste0("ph_cpnA_vec")]]
    seg.CCF.val.list[[seg.name]][[samp]][[paste0(cpnX.string,"_orig_ph_obs_cpnB_vec")]]   <- all.pw.cn.het.comp.lists[[samp]][[seg.name]][[paste0("ph_cpnB_vec")]]
    seg.CCF.val.list[[seg.name]][[samp]][[paste0("clean_ph_obs_", cpnX.string, "_vec")]]   <- ss.ph.cpnX.vals.vec
    seg.CCF.val.list[[seg.name]][[samp]][[paste0("frac_", cpnX.string, "_low_vec")]]  <- ss.frac.poss.sc.cpnX.low.vec
    seg.CCF.val.list[[seg.name]][[samp]][[paste0("frac_", cpnX.string, "_high_vec")]]  <- ss.frac.poss.sc.cpnX.high.vec
    seg.CCF.val.list[[seg.name]][[samp]][[paste0("poss_sc_", cpnX.string, "_low")]] <- rep(ss.poss.sc.cpnX.low, length(num.vals))
    seg.CCF.val.list[[seg.name]][[samp]][[paste0("poss_sc_", cpnX.string, "_high")]] <- rep(ss.poss.sc.cpnX.high, length(num.vals))

  }#end for ( samp in sample.vec ) {

  return(seg.CCF.val.list)

}#end updateSegCCFValList <- function(




##
##Params:
##
##Returns:
##
calculateTwoValueDiffFlag <- function(s.samps.int.ob.xp.cpn.df, seg.cpnA.cnpB.diff.exp.list, cpnX.string, min.CCF.to.consider.val, confidence.interval.threshold) {

  two.value.diff.flag <- TRUE
  samps.mean.cpnX.vec <- s.samps.int.ob.xp.cpn.df[, paste0("mean_obs_", cpnX.string)]
  low.int.mean.cpnX <- min(floor(mean.cpnX.vec))
  top.int.mean.cpnX <- max(ceiling(mean.cpnX.vec))

  if (abs(min(top.int.mean.cpnX - samps.mean.cpnX.vec)) >= abs(max(low.int.mean.cpnX - samps.mean.cpnX.vec))) {

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
  #Todo, make this return a single value (split into two functions ideally)
  return(two.value.diff.flag)

}


chooseLowHighSubclonalCPNWithOnlyTwoPossibilities <- function(s.samps.int.ob.xp.cpn.df,
                                                              seg.name,
                                                              all.samp.ph.BAF.df.list,
                                                              purity.est.vec,
                                                              cpnX.string,
                                                              min.CCF.to.consider.val,
                                                              confidence.interval.threshold)
                                                              {
  samps.mean.cpnX.vec <- s.samps.int.ob.xp.cpn.df[, paste0("mean_obs_", cpnX.string)]
  # check whether only one state when rounding
  if (length(unique(round(samps.mean.cpnX.vec))) == 1) {

        all.clonal.flag <- determineCCFClonalityOfSegment(s.samps.int.ob.xp.cpn.df = s.samps.int.ob.xp.cpn.df,
                            seg.name = seg.name,
                            all.samp.ph.BAF.df.list= all.samp.ph.BAF.df.list,
                            purity.est.vec = purity.est.vec,
                            cpnX.string = cpnX.string,
                            min.CCF.to.consider.val = min.CCF.to.consider.val,
                            confidence.interval.threshold = confidence.interval.threshold)

        if (all.clonal.flag == TRUE) {

          ss.poss.sc.cpnX.low =  round(samps.mean.cpnX.vec)
          ss.poss.sc.cpnX.high = round(samps.mean.cpnX.vec)

        } else {#if not all clonal: 

          ss.poss.sc.cpnX.low =  min(unique(floor(samps.mean.cpnX.vec)))
          ss.poss.sc.cpnX.high = max(unique(ceiling(samps.mean.cpnX.vec)))

        }#end if (all.clonal.flag == TRUE) {

  }#end if (length(unique(round(s.samps.int.ob.xp.cpn.df[, paste0("mean_obs_", cpnX.string)]))) == 1) {

  if (length(unique(round(samps.mean.cpnX.vec))) != 1) {

      ss.poss.sc.cpnX.low =  min(unique(floor(samps.mean.cpnX.vec)))
      ss.poss.sc.cpnX.high = max(unique(ceiling(samps.mean.cpnX.vec)))

  }#end if (length(unique(round(s.samps.int.ob.xp.cpn.df[, paste0("mean_obs_", cpnX.string)]))) != 1) {


  return(c(ss.poss.sc.cpnX.low, ss.poss.sc.cpnX.high))

}



chooseLowHighSubclonalCPNWithMoreThanTwoPossibilities <- function(s.samps.int.ob.xp.cpn.df,
                                                                  seg.cpnA.cnpB.diff.exp.list,
                                                                  seg.name,
                                                                  all.samp.ph.BAF.df.list,
                                                                  purity.est.vec,
                                                                  cpnX.string,
                                                                  min.CCF.to.consider.val,
                                                                  confidence.interval.threshold) {

  samps.mean.cpnX.vec <- s.samps.int.ob.xp.cpn.df[, paste0("mean_obs_", cpnX.string)]
  #check whether there is ony one state when rounding, can this ever be called?
  if (length(unique(round(samps.mean.cpnX.vec))) == 1) {

        all.clonal.flag <- determineCCFClonalityOfSegment(s.samps.int.ob.xp.cpn.df = s.samps.int.ob.xp.cpn.df,
                            seg.name = seg.name,
                            all.samp.ph.BAF.df.list= all.samp.ph.BAF.df.list,
                            purity.est.vec = purity.est.vec,
                            cpnX.string = cpnX.string,
                            min.CCF.to.consider.val = min.CCF.to.consider.val,
                            confidence.interval.threshold = confidence.interval.threshold)

        if (all.clonal.flag == TRUE) {

          ss.poss.sc.cpnX.low =  round(samps.mean.cpnX.vec)
          ss.poss.sc.cpnX.high = round(samps.mean.cpnX.vec)

        } else {#end if (all.clonal.flag == TRUE) {
              #Not all clonal but all cpn rounds to the same cpn, question is:
              two.value.diff.flag <- calculateTwoValueDiffFlag(s.samps.int.ob.xp.cpn.df = s.samps.int.ob.xp.cpn.df,
                                                              seg.cpnA.cnpB.diff.exp.list = seg.cpnA.cnpB.diff.exp.list,
                                                              cpnX.string = cpnX.string,
                                                              min.CCF.to.consider.val = min.CCF.to.consider.val,
                                                              confidence.interval.threshold = confidence.interval.threshold)

              max.ceiling.cpnX <- max(ceiling(samps.mean.cpnX.vec))
              min.floor.cpnX <- min(floor(samps.mean.cpnX.vec))
              #minimum/maximum absolute difference from maximum possible subclonal copy number
              min.diff.from.max.psc.cpn <- abs(min(max.ceiling.cpnX - samps.mean.cpnX.vec))
              max.diff.from.min.psc.cpn <- abs(max(min.floor.cpnX - samps.mean.cpnX.vec))

              if (two.value.diff.flag == TRUE) {
                    #Do we twist up or do we twist down
                    if (min.diff.from.max.psc.cpn < max.diff.from.min.psc.cpn) {
                        #take the ceiling
                        ss.poss.sc.cpnX.low = round(samps.mean.cpnX.vec)
                        ss.poss.sc.cpnX.high = max(unique(ceiling(samps.mean.cpnX.vec)))

                    } else {#(min.diff.from.max.psc.cpn >= max.diff.from.min.psc.cpn) {      
                        # take the floor
                        ss.poss.sc.cpnX.low =  min(unique(floor(samps.mean.cpnX.vec)))
                        ss.poss.sc.cpnX.high = max(unique(round(samps.mean.cpnX.vec)))

                    }#end if (min.diff.from.max.psc.cpn >= max.diff.from.min.psc.cpn) {    

              } else {#end if (two.value.diff.flag == TRUE) {

                    # do we twist up or down
                    if (min.diff.from.max.psc.cpn < max.diff.from.min.psc.cpn) {
                        #take the ceiling
                        ss.poss.sc.cpnX.low =  min(unique(floor(samps.mean.cpnX.vec)))
                        ss.poss.sc.cpnX.high = max(unique(ceiling(samps.mean.cpnX.vec)))

                    } else {#(min.diff.from.max.psc.cpn >= max.diff.from.min.psc.cpn) {  
                        #take the floor #TODO, check with Nicky but these values are the same as the < than, BUG?
                        #Maybe don't treat as one continuous factor... we're treating a loss and a gain as one event...
                        #Flag if gain and losses either here or in the segmentation...
                        #If you have three samples one wth cpn 1 and cpn4 and cpn8, cpn4 would be completely subclonal
                        ss.poss.sc.cpnX.low =  min(unique(floor(samps.mean.cpnX.vec)))
                        ss.poss.sc.cpnX.high = max(unique(ceiling(samps.mean.cpnX.vec)))

                    }#end if (min.diff.from.max.psc.cpn >= max.diff.from.min.psc.cpn) {    

              }#end  two.value.diff.flag==FALSE

        }#end  else {#end if (all.clonal.flag == FALSE) {

  }#end if (length(unique(round(samps.mean.cpnX.vec))) == 1) {
  #check whether there are only two states when rounding
  if (length(unique(round(samps.mean.cpnX.vec))) == 2) {

    ss.poss.sc.cpnX.low =  min(round(samps.mean.cpnX.vec))
    ss.poss.sc.cpnX.high = max(round(samps.mean.cpnX.vec))

  }#end if (length(unique(round(samps.mean.cpnX.vec))) == 2) {
  #check whether there are more than two states when rounding
  if (length(unique(round(samps.mean.cpnX.vec))) > 2) {

    ss.poss.sc.cpnX.low =  min(unique(floor(samps.mean.cpnX.vec + min.CCF.to.consider.val)))
    ss.poss.sc.cpnX.high = max(unique(ceiling(samps.mean.cpnX.vec + min.CCF.to.consider.val)))

  }#end if (length(unique(round(samps.mean.cpnX.vec))) > 2) {

    return(c(ss.poss.sc.cpnX.low, ss.poss.sc.cpnX.high))

}


getIndividualCCFestimates_withPvals <- function(seg.name,
                                                purity.est.vec,
                                                all.pw.cn.het.comp.lists,
                                                all.samp.ph.BAF.df.list,
                                                samps.int.cpn.raw.cpn.by.chr.list,
                                                max.cpn.val = 5,
                                                min.CCF.to.consider.val = 0.15,
                                                chrs.can.use.for.ccf.comp = c(1:22),
                                                confidence.interval.threshold = 0.9999,
                                                min.num.other.chrs = 3) {

      if (is.na(all.pw.cn.het.comp.lists[[1]][[seg.name]]$ph_cpnA_vec[1])) {

            return(NA)
      }
      sample.vec <- names(purity.est.vec)
      chr <- strsplit(seg.name, split=":")[1]
      s.samps.int.ob.xp.cpn.df <- getSegObsExpCPNBySampleDF(all.pw.cn.het.comp.lists, sample.vec, seg.name)
      chrs.for.ccf.comp.vec <- chrs.can.use.for.ccf.comp[ ! chrs.can.use.for.ccf.comp%in% chr] #Going to be used to assess CCFs across the each sample
      #Need to initialise all the lists outside the loop now we are doing cpnA cpnB in a loop and cpnB entries would wipe cpnA ones...
      seg.CCF.val.list <- list()
      seg.CCF.val.list[[seg.name]] <- list()
      seg.cpnA.cnpB.diff.exp.list <- list()
      for ( samp in sample.vec ) {
        seg.CCF.val.list[[seg.name]][[samp]] <- list()
        seg.cpnA.cnpB.diff.exp.list[[samp]] <- list()
      }
      for ( cpnX.string in c("cpnA", "cpnB") ) {

            seg.cpnA.cnpB.diff.exp.list <- calculateWhetherSegCPNDiffFromIntegerExpectedAndAllMatchingSampCPN(samp.pvals.diff.from.exp.list=seg.cpnA.cnpB.diff.exp.list,
                                                                                                              all.pw.cn.het.comp.lists,
                                                                                                              s.samps.int.ob.xp.cpn.df,
                                                                                                              seg.name,
                                                                                                              chrs.for.ccf.comp.vec,
                                                                                                              samps.int.cpn.raw.cpn.by.chr.list,
                                                                                                              cpnA.or.cpnB=cpnX.string)
            
            use.preset.mean.frac.cpnX.low.high <- FALSE #Todo, decide whether this should be default parameter
            samps.mean.cpnX.vec <- s.samps.int.ob.xp.cpn.df[, paste0("mean_obs_", cpnX.string)]
            #Catch the edge case that all the seg's mean cpnX are equal
            if (length(unique(samps.mean.cpnX.vec)) == 1) {

              ss.poss.sc.cpnX.low = 1#Should really be ploidy n this case.
              ss.poss.sc.cpnX.high = unique(samps.mean.cpnX.vec)
              low.high.ss.poss.cpnX.vec <- c(ss.poss.sc.cpnX.low, ss.poss.sc.cpnX.low)
              use.preset.mean.frac.cpnX.low.high = TRUE 

            }#end if (length(unique(round(samps.mean.cpnX.vec))) == 1) {
            if (length(unique(c(floor(samps.mean.cpnX.vec), ceiling(samps.mean.cpnX.vec)))) == 2) {

               low.high.ss.poss.cpnX.vec <- chooseLowHighSubclonalCPNWithOnlyTwoPossibilities(s.samps.int.ob.xp.cpn.df,
                                                                             seg.name,
                                                                             all.samp.ph.BAF.df.list,
                                                                             purity.est.vec,
                                                                             cpnX.string,
                                                                             min.CCF.to.consider.val,
                                                                             confidence.interval.threshold)

            }#end if (length(unique(c(floor(samps.mean.cpnX.vec), ceiling(samps.mean.cpnX.vec)))) == 2) {
            if (length(unique(c(floor(samps.mean.cpnX.vec), ceiling(samps.mean.cpnX.vec)))) > 2) {

              low.high.ss.poss.cpnX.vec <- chooseLowHighSubclonalCPNWithMoreThanTwoPossibilities(s.samps.int.ob.xp.cpn.df,
                                                                    seg.cpnA.cnpB.diff.exp.list,
                                                                    seg.name,
                                                                    all.samp.ph.BAF.df.list,
                                                                    purity.est.vec,
                                                                    cpnX.string,
                                                                    min.CCF.to.consider.val,
                                                                    confidence.interval.threshold)

            }#end (length(unique(c(floor(samps.mean.cpnX.vec), ceiling(samps.mean.cpnX.vec)))) > 2) {
                  #ceiling(s.samps.int.ob.xp.cpn.df[, paste0("mean_obs_", cpnX.string)]))) > 2) {
        seg.CCF.val.list <- updateSegCCFValList(seg.name = seg.name, seg.CCF.val.list = seg.CCF.val.list, sample.vec = sample.vec, s.samps.int.ob.xp.cpn.df = s.samps.int.ob.xp.cpn.df, 
                                           all.pw.cn.het.comp.lists = all.pw.cn.het.comp.lists, cpnX.string = cpnX.string, max.cpn.val = max.cpn.val,
                                           ###These parameters vary:
                                           ss.poss.sc.cpnX.low = low.high.ss.poss.cpnX.vec[1],
                                           ss.poss.sc.cpnX.high = low.high.ss.poss.cpnX.vec[2],
                                           use.preset.mean.frac.cpnX.low.high = use.preset.mean.frac.cpnX.low.high)   
      }#end for ( cpnX.string in c("cpnA", "cpnB") )
      #Todo refactor this seg_diff_form expected so it's generated by another function call...
      seg.CCF.val.list[[seg.name]]$seg_diff_from_expected_list <- seg.cpnA.cnpB.diff.exp.list
      return(seg.CCF.val.list)

}#end getIndividualCCFestimates_withPvals <- function(...)










                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  