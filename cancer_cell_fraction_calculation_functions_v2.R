##
##Params:
##
##Returns:
##
generateAllSamplesByChromIntCPNUnderlyingRawValuesLists <- function(asas.ph.raw.cpn.lists, sample.vec, seg.name.vec, cpnX.vec=c("cpnA", "cpnB"), reduced.seg.flag=FALSE) {
      sample.chr.int.cpn.raw.vals.list <- list()
      #s.samps.obs.exp.cpn.df <- data.frame()
      #But really the datastructure organisation overall needs to change if we speed.
      for ( samp in sample.vec ) {

            for ( seg.name in seg.name.vec ) {

                  if (reduced.seg.flag == FALSE) {

                    if (is.na(asas.ph.raw.cpn.lists[[1]][[seg.name]]$ph_cpnA_vec[1]) | is.na(asas.ph.raw.cpn.lists[[1]][[seg.name]]$ph_cpnB_vec[1])) {

                          next
                    }#end if (is.na(asas.ph.raw.cpn.lists[[1]][[seg.name]]$ph_cpnA_vec[1]) | is.na(asas.ph.raw.cpn.lists[[1]][[seg.name]]$ph_cpnB_vec[1])) {
                  }#if (reduced.seg.flag == FALSE) {

                  for ( cpnX in cpnX.vec ) {

                    if(reduced.seg.flag == FALSE) {
                      samp.seg.int.cpnX.val <- round(asas.ph.raw.cpn.lists[[samp]][[seg.name]][[paste0("mean_ph_",cpnX)]])

                      chr <- unlist(strsplit(seg.name, split=":"))[1]
                      sample.chr.int.cpn.raw.vals.list[[samp]][[as.character(samp.seg.int.cpnX.val)]][[chr]][[paste0("all_seg_",cpnX,"_vec")]] <- c(
                            sample.chr.int.cpn.raw.vals.list[[samp]][[as.character(samp.seg.int.cpnX.val)]][[chr]][[paste0("all_seg_",cpnX,"_vec")]], 
                            asas.ph.raw.cpn.lists[[samp]][[seg.name]][[paste0("ph_",cpnX,"_vec")]])

                    } else {

                      samp.seg.int.cpnX.val <- round(mean(asas.ph.raw.cpn.lists[[seg.name]][[paste0(cpnX,"_orig_ph_obs_",cpnX,"_vec")]]))

                      chr <- unlist(strsplit(seg.name, split=":"))[1]
                      sample.chr.int.cpn.raw.vals.list[[samp]][[as.character(samp.seg.int.cpnX.val)]][[chr]][[paste0("all_seg_",cpnX,"_vec")]] <- c(
                            sample.chr.int.cpn.raw.vals.list[[samp]][[as.character(samp.seg.int.cpnX.val)]][[chr]][[paste0("all_seg_",cpnX,"_vec")]], 
                            asas.ph.raw.cpn.lists[[seg.name]][[paste0(cpnX,"_orig_ph_obs_",cpnX,"_vec")]])

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
calculateWhetherSegCPNDiffFromIntegerExpectedAndAllMatchingSampCPN <- function(asas.ph.raw.cpn.lists,
                                                                              s.samps.int.ob.xp.cpn.df,
                                                                              seg.name,
                                                                              samps.int.cpn.raw.cpn.by.chr.list,
                                                                              cpnX,
                                                                              chrs.can.use.for.ccf.test,
                                                                              reduced.seg.flag = FALSE) {

      samp.pvals.diff.from.exp.list <- list()

      sample.vec <- rownames(s.samps.int.ob.xp.cpn.df)

      for ( samp in sample.vec ) {

        samp.pvals.diff.from.exp.list[[samp]] <- list()

      }

      chr <- strsplit(seg.name, split=":")[1]
      chrs.for.ccf.comp.vec <- chrs.can.use.for.ccf.test[ ! chrs.can.use.for.ccf.test %in% chr] #Going to be used to assess CCFs across the each sample
      
      #Todo refacor so doesn't need the reduced seg flag.
      if (reduced.seg.flag == FALSE) {
        if (is.na(asas.ph.raw.cpn.lists[[1]][[seg.name]]$ph_cpnA_vec[1]) | is.na(asas.ph.raw.cpn.lists[[1]][[seg.name]]$ph_cpnB_vec[1])) {

            next
        }
      }

      for ( samp in sample.vec ) {

            chrs.w.same.cpn.vec <- c()
            all.match.cpn.vals.vec <- c()


            for ( chr in chrs.for.ccf.comp.vec ) {

                  s.chr.cpn.vec <- samps.int.cpn.raw.cpn.by.chr.list[[samp]][[as.character(s.samps.int.ob.xp.cpn.df[samp, paste0("round_int_", cpnX)])]][[as.character(chr)]][[paste0("all_seg_", cpnX, "_vec")]]

                  if ( length(s.chr.cpn.vec) != 0 ) {

                        chrs.w.same.cpn.vec <- c(chrs.w.same.cpn.vec, as.character(chr))
                        all.match.cpn.vals.vec <- c(all.match.cpn.vals.vec, s.chr.cpn.vec)

                  }#end if ( length(s.chr.cpn.vec) != 0 )

            }#endfor ( chr in chrs.for.ccf.comp.vec ) {

            if (reduced.seg.flag == FALSE) {

              seg.raw.cpn.vals <- asas.ph.raw.cpn.lists[[samp]][[seg.name]][[paste0("ph_",cpnX,"_vec")]]

            } else {

              seg.raw.cpn.vals <- asas.ph.raw.cpn.lists[[seg.name]][[paste0(cpnX,"_orig_ph_obs_",cpnX,"_vec")]]

            }

            t.test.raw.seg.cpn.vs.all.matching.cpn.raw <- tTestWithErrorHandling(seg.raw.cpn.vals, all.match.cpn.vals.vec)
            t.test.raw.seg.vs.cpn.int <- tTestWithErrorHandling(seg.raw.cpn.vals, mu = s.samps.int.ob.xp.cpn.df[samp, paste0("round_int_", cpnX)])
            median.seg.raw.cpn.greater.than.matched.cpn <- median(seg.raw.cpn.vals) > median(all.match.cpn.vals.vec)
            median.seg.raw.cpn.greater.seg.int.cpn   <- median(seg.raw.cpn.vals) > s.samps.int.ob.xp.cpn.df[samp, paste0("round_int_", cpnX)]

            samp.pvals.diff.from.exp.list[[samp]][[paste0("round_int_", cpnX)]] <- s.samps.int.ob.xp.cpn.df[samp, paste0("round_int_", cpnX)]
            samp.pvals.diff.from.exp.list[[samp]][[paste0(cpnX, "_pval_vs_oth_cpn_match_chrs_raw_cpn")]] <- t.test.raw.seg.cpn.vs.all.matching.cpn.raw$p.value
            samp.pvals.diff.from.exp.list[[samp]][[paste0(cpnX, "_pval_vs_round_int_cpn")]] <- t.test.raw.seg.vs.cpn.int$p.value
            #print("chrs.w.same.cpn.vec")
            #print(chrs.w.same.cpn.vec)
            #print(str(samp.pvals.diff.from.exp.list, max=1))
            samp.pvals.diff.from.exp.list[[samp]][[paste0(cpnX, "_num_chrs_w_match_int_cpn_vec")]]  <- chrs.w.same.cpn.vec
            samp.pvals.diff.from.exp.list[[samp]][[paste0("median_", cpnX, "_greater_than_match_raw_cpn")]] <- median.seg.raw.cpn.greater.than.matched.cpn
            samp.pvals.diff.from.exp.list[[samp]][[paste0("median_", cpnX, "_greater_than_match_int_cpn")]] <- median.seg.raw.cpn.greater.seg.int.cpn

      }#end for ( samp in sample.vec ) {

      return(samp.pvals.diff.from.exp.list)

}



##
##Params:
##
##Returns:
##
getSegObsExpCPNBySampleDF <- function(asas.ph.raw.cpn.lists, sample.vec, seg.name, cpnX.vec = c("cpnA", "cpnB"), reduced.seg.flag = FALSE) {

      num.samples <- length(sample.vec)
      s.samps.obs.exp.cpn.df <- data.frame(round_int_cpnA = numeric(num.samples),#rounded
                                          round_int_cpnB = numeric(num.samples),
                                          mean_obs_cpnA = numeric(num.samples),
                                          mean_obs_cpnB = numeric(num.samples),
                                          exp_cpnA = numeric(num.samples),#this expected is just to remove 0 values whi
                                          exp_cpnB = numeric(num.samples),
                                          stringsAsFactors = FALSE)
      rownames(s.samps.obs.exp.cpn.df) <- sample.vec    
      for ( cpnX in cpnX.vec ) {

        for ( samp in sample.vec ) {
              
              #print(str(asas.ph.raw.cpn.lists[[seg.name]][[paste0(cpnX,"_orig_ph_obs_",cpnX,"_vec")]], max=1))

              if (reduced.seg.flag == FALSE) {

                mean.cpnX <- asas.ph.raw.cpn.lists[[samp]][[seg.name]][[paste0("mean_ph_",cpnX)]]
                s.samps.obs.exp.cpn.df[samp, paste0("round_int_",cpnX)] <- round(mean.cpnX)


              } else {
                mean.cpnX <- mean(asas.ph.raw.cpn.lists[[seg.name]][[paste0(cpnX,"_orig_ph_obs_",cpnX,"_vec")]])
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






getPossibleLowAndHighSubclonalCPNStates <- function(seg.name,
                                                cpnX,
                                                s.samps.int.ob.xp.cpn.df,
                                                all.clonal.flag,
                                                two.value.diff.flag,
                                                min.CCF.to.consider.val) {

      #Need to initialise all the lists outside the loop now we are doing cpnA cpnB in a loop and cpnB entries would wipe cpnA ones...
      use.preset.mean.frac.cpnX.low.high <- FALSE #Todo, decide whether this should be default parameter
      samps.mean.cpnX.vec <- s.samps.int.ob.xp.cpn.df[, paste0("mean_obs_", cpnX)]

      #Catch the edge case that all the seg's mean cpnX are equal
      if (length(unique(samps.mean.cpnX.vec)) == 1) {

        ss.poss.sc.cpnX.low <- 1#Should really be ploidy n this case.
        ss.poss.sc.cpnX.high <- unique(samps.mean.cpnX.vec)
        low.high.ss.poss.cpnX.vec <- c(ss.poss.sc.cpnX.low, ss.poss.sc.cpnX.low)
        #So what does this mean exactly?s
        use.preset.mean.frac.cpnX.low.high  <-  TRUE 

      }#end if (length(unique(round(samps.mean.cpnX.vec))) == 1) {
      #Is it a simple case there are only two cpn states found by floor and ceiling on the phased cpnX for all samples.  
      if (length(unique(c(floor(samps.mean.cpnX.vec), ceiling(samps.mean.cpnX.vec)))) == 2) {

         low.high.ss.poss.cpnX.vec <- chooseLowHighSubclonalCPNWithOnlyTwoPossibilities(s.samps.int.ob.xp.cpn.df,
                                                                                        cpnX,
                                                                                        all.clonal.flag,
                                                                                        two.value.diff.flag)

      }#end if (length(unique(c(floor(samps.mean.cpnX.vec), ceiling(samps.mean.cpnX.vec)))) == 2) {
      #Is it a more complex case where there are more than two cpn states on the phased cpnX over all the samples.  
      if (length(unique(c(floor(samps.mean.cpnX.vec), ceiling(samps.mean.cpnX.vec)))) > 2) {

        low.high.ss.poss.cpnX.vec <- chooseLowHighSubclonalCPNWithMoreThanTwoPossibilities(s.samps.int.ob.xp.cpn.df,
                                                                                          cpnX,
                                                                                          all.clonal.flag,
                                                                                          two.value.diff.flag,
                                                                                          min.CCF.to.consider.val)

      }#end (length(unique(c(floor(samps.mean.cpnX.vec), ceiling(samps.mean.cpnX.vec)))) > 2) {
            #ceiling(s.samps.int.ob.xp.cpn.df[, paste0("mean_obs_", cpnX)]))) > 2) {

      return(c(low.high.ss.poss.cpnX.vec, use.preset.mean.frac.cpnX.low.high))

}#end getIndividualCCFestimates_withPvals <- function(...)




#Can be refactored further...
determineCCFClonalityOfSegment <- function(s.samps.int.ob.xp.cpn.df,
                    asas.ph.raw.cpn.lists,
                    seg.name,
                    all.samp.ph.BAF.df.list,
                    purity.est.vec,
                    cpnX,
                    min.CCF.to.consider.val,
                    confidence.interval.threshold) {
  sample.vec <- names(purity.est.vec)
  ss.round.mean.obs.cpn <- round(s.samps.int.ob.xp.cpn.df[, paste0("mean_obs_", cpnX)])[1]
  # first check whether all regions are clonal
  all.samp.cpnX.clonal.t.test <- TRUE
  for ( samp in sample.vec ) {
        #sample seg phased cpnX vector of values...
        ss.ph.cpnX.vec <- asas.ph.raw.cpn.lists[[samp]][[seg.name]][[paste0("ph_", cpnX, "_vec")]]
        samp.ph.cpnX.vs.samp.round.mean.obs.cpnX.test <- tTestWithErrorHandling(ss.ph.cpnX.vec, mu = ss.round.mean.obs.cpn, conf.level = confidence.interval.threshold)

        if (samp.ph.cpnX.vs.samp.round.mean.obs.cpnX.test$p.val < 1 - confidence.interval.threshold) {
              #print(" (samp.ph.cpnX.vs.samp.round.mean.obs.cpnX.test$p.val < 1 - confidence.interval.threshold) { so all clonal FALSE")
              all.samp.cpnX.clonal.t.test <- FALSE

        }#end if (samp.ph.cpnX.vs.samp.round.mean.obs.cpnX.test$p.val < 1 - confidence.interval.threshold) {

  }#for ( samp in sample.vec ) {
  seg.chr.start.end.vec <- unlist(strsplit(seg.name, split=":"))
  chr <- as.numeric(seg.chr.start.end.vec[1])
  start <- as.numeric(seg.chr.start.end.vec[2])
  end <- as.numeric(seg.chr.start.end.vec[3])
  all.samp.BAF.clonal.t.test <- TRUE
  for ( samp in sample.vec ) {
        #Get BAF values from the dataframe
        samp.ph.baf.df <- all.samp.ph.BAF.df.list[[samp]]
        s.samp.ph.baf.df <- samp.ph.baf.df[samp.ph.baf.df$chrom == chr & samp.ph.baf.df$pos >= start & samp.ph.baf.df$pos <= end, ,drop=FALSE]
        samp.maj.is.cpnA.cpnB.BAF.vs.exp.BAF.test <- tTestWithErrorHandling(c(s.samp.ph.baf.df[s.samp.ph.baf.df$maj_is_cpnA_logic_vec == 1, "BAF"],#tTestWithErrorHandling(c
                                                                          1 - s.samp.ph.baf.df[s.samp.ph.baf.df$maj_is_cpnB_logic_vec == 1, "BAF"]),
                                                                          mu = calcExpectedBAF(s.samps.int.ob.xp.cpn.df[samp, "exp_cpnB"],
                                                                                               s.samps.int.ob.xp.cpn.df[samp, "exp_cpnA"],
                                                                                                purity.est.vec[samp]))

        if (samp.maj.is.cpnA.cpnB.BAF.vs.exp.BAF.test$p.val < (1 - confidence.interval.threshold)) {

              all.samp.BAF.clonal.t.test <- FALSE

        }#if (samp.maj.is.cpnA.cpnB.BAF.vs.exp.BAF.test$p.val < (1 - confidence.interval.threshold)) {

  }#end for ( samp in sample.vec ) {
  #Finally check what the CCF is of each region if we assume subclonality
  max.CCF <- 0
  poss.sc.cpnX.low <- min(unique(floor(s.samps.int.ob.xp.cpn.df[, paste0("mean_obs_", cpnX)])))
  poss.sc.cpnX.high <- max(unique(ceiling(s.samps.int.ob.xp.cpn.df[, paste0("mean_obs_", cpnX)])))
  for (samp in sample.vec) {
        print(s.samps.int.ob.xp.cpn.df)
        ss.mean.obs.ph.cpnX  <- s.samps.int.ob.xp.cpn.df[samp, paste0("mean_obs_", cpnX)]
        ss.fracs.to.use.vec <- getCCFractionNeeded(ss.mean.obs.ph.cpnX, poss.sc.cpnX.low, poss.sc.cpnX.high)
        ss.min.CCF <- min(ss.fracs.to.use.vec)

        if (ss.min.CCF > max.CCF) {

              max.CCF <- ss.min.CCF

        }#end if (seg.min.CCF > maxCCF) {

  }#end or (samp in sample.vec) {
  #Now determine whether to explore subclonality:
  if (all.samp.cpnX.clonal.t.test == TRUE | all.samp.BAF.clonal.t.test == TRUE | max.CCF < min.CCF.to.consider.val) {

    final.all.clonal.flag <- TRUE

  } else {#end if (all.clonal.flag == TRUE | all.samp.clonal.t.test == TRUE | max.CCF min.CCF.to.consider.val){

    final.all.clonal.flag <- FALSE

  }

  return(final.all.clonal.flag)

}#determineCCFClonalityOfSegment <- function(...)


##
##Params:
##
##Returns:
##
getSegCPNXMeanAndAllRawCCFsDF <- function(ss.mean.cpnX.obs, ss.ph.cpnX.vals.vec, seg.psc.cpnX.low, seg.psc.cpnX.high, seg.use.default.ccf.flag, all.clonal.flag, max.cpn.val) {


  ss.mean.cpnX.obs <- ifelse(ss.mean.cpnX.obs < 0, 0, ss.mean.cpnX.obs)

  if (seg.use.default.ccf == TRUE) {

    ss.mean.frac.cpnX.low <- 0
    ss.mean.frac.cpnX.high <- 1

  } else {

    ss.fracs.to.use.vec <- getCCFractionNeeded(ss.mean.cpnX.obs, seg.psc.cpnX.low, seg.psc.cpnX.high)
    ss.mean.frac.cpnX.low <- ss.fracs.to.use.vec[1]
    ss.mean.frac.cpnX.high <- ss.fracs.to.use.vec[2]

  }#end if (use.preset.mean.frac.cpnX.low.high == TRUE) {

  ss.ph.cpnX.vals.vec <- ifelse(ss.ph.cpnX.vals.vec < 0, 0, ss.ph.cpnX.vals.vec)
  ss.ph.cpnX.vals.vec <- ifelse(ss.ph.cpnX.vals.vec >= max.cpn.val, max.cpn.val, ss.ph.cpnX.vals.vec)
  ss.fracs.to.use.vec <- sapply(ss.ph.cpnX.vals.vec, getCCFractionNeeded, seg.psc.cpnX.low, seg.psc.cpnX.high)
  ss.fracs.to.use.vec[ss.fracs.to.use.vec == "Inf"] <- 1

  #Todo, refactor this to take into account whether clonal or not, if clonal:
  if (all.clonal.flag == TRUE) { 

    ss.frac.poss.sc.cpnX.low.vec <- apply(ss.fracs.to.use.vec, 2, min)
    ss.frac.poss.sc.cpnX.high.vec <- apply(ss.fracs.to.use.vec, 2, max)

  } else {  #if non clonal:

    ss.frac.poss.sc.cpnX.low.vec <- ss.fracs.to.use.vec[1 ,]
    ss.frac.poss.sc.cpnX.high.vec <- ss.fracs.to.use.vec[2 ,]

  }

  num.rows <- length(ss.ph.cpnX.vals.vec)

  seg.cpnX.psc.cpn.CCF.df <- data.frame(clean_mean_cpnX = rep(ss.mean.cpnX.obs, num.rows),
                                       clean_ph_obs_cpnX = ss.ph.cpnX.vals.vec,
                                        poss_sc_cpnX_low = rep(seg.psc.cpnX.low, num.rows),
                                        poss_sc_cpnX_high = rep(seg.psc.cpnX.high, num.rows),
                                        mean_frac_cpnX_low = rep(ss.mean.frac.cpnX.low, num.rows),
                                        mean_frac_cpnX_high = rep(ss.mean.frac.cpnX.high, num.rows),
                                        frac_cpnX_low = ss.frac.poss.sc.cpnX.low.vec,
                                        frac_cpnX_high = ss.frac.poss.sc.cpnX.high.vec,
                                        stringsAsFactors = FALSE
                                        )

  return(seg.cpnX.psc.cpn.CCF.df)
}



##
##Params:
##
##Returns:
##
calculateTwoValueDiffFlag <- function(s.samps.int.ob.xp.cpn.df, seg.cpnX.diff.from.exp.list, cpnX, min.CCF.to.consider.val, confidence.interval.threshold, min.num.other.chrs) {

  two.value.diff.flag <- TRUE
  samps.mean.cpnX.vec <- s.samps.int.ob.xp.cpn.df[, paste0("mean_obs_", cpnX)]
  low.int.mean.cpnX <- min(floor(samps.mean.cpnX.vec))
  top.int.mean.cpnX <- max(ceiling(samps.mean.cpnX.vec))

  if (abs(min(top.int.mean.cpnX - samps.mean.cpnX.vec)) >= abs(max(low.int.mean.cpnX - samps.mean.cpnX.vec))) {

      ss.poss.sc.cpnX.low <- low.int.mean.cpnX
      ss.poss.sc.cpnX.high <- top.int.mean.cpnX

      #sis this choice of ss.poss.sc.cpnX.high, too conservative?
      samps.poss.more.than.sc.cnpX.high.vec  <- rownames(s.samps.int.ob.xp.cpn.df[s.samps.int.ob.xp.cpn.df[, paste0("mean_obs_", cpnX)] >= ss.poss.sc.cpnX.high, ])

      for ( samp in samps.poss.more.than.sc.cnpX.high.vec ) {
        samp.max.CCF <- abs(ss.poss.sc.cpnX.high - s.samps.int.ob.xp.cpn.df[samp, paste0("mean_obs_", cpnX)])
        raw.diff.from.int.cpnX.pval <- seg.cpnX.diff.from.exp.list[[samp]][[paste0(cpnX, "_pval_vs_round_int_cpn")]]
        diff.from.int.flag <- raw.diff.from.int.cpnX.pval < (1 - confidence.interval.threshold)

        diff.from.other.chrs.flag <- TRUE
        if (length(seg.cpnX.diff.from.exp.list[[samp]][[paste0(cpnX, "_chrs_w_match_int_cpn_vec")]]) > min.num.other.chrs) {

          diff.from.other.chrs.flag <- seg.cpnX.diff.from.exp.list[[samp]][[paste0(cpnX, "_pval_vs_oth_cpn_match_chrs_raw_cpn")]] < (1 - confidence.interval.threshold)

        }#end if (length(seg.cpnX.diff.from.exp.list[[samp]][[paste0(cpnX, "_chrs_w_match_int_cpn_vec")]]) > min.num.other.chrs) {

        if (diff.from.int.flag == TRUE & diff.from.other.chrs.flag == TRUE & samp.max.CCF >= min.CCF.to.consider.val) {

          two.value.diff.flag <- FALSE

        }

      }#end for ( samp in sample.vec ) {

  }#end if(abs(min(top.int.mean.cpnX - mean.cpnX.vec)) >= abs(max(low.int.mean.cpnX - mean.cpnX.vec))) {

  if (abs(min(top.int.mean.cpnX - samps.mean.cpnX.vec)) < abs(max(low.int.mean.cpnX - samps.mean.cpnX.vec))) {

    ss.poss.sc.cpnX.low  <- round(samps.mean.cpnX.vec)[1]
    ss.poss.sc.cpnX.high  <- max(unique(ceiling(s.samps.int.ob.xp.cpn.df[, paste0("mean_obs_", cpnX)])))

    # is this choice of ss.poss.sc.cpnX.low, too conservative
    samps.poss.less.than.sc.cnpX.low.vec <- rownames(s.samps.int.ob.xp.cpn.df[s.samps.int.ob.xp.cpn.df[, paste0("mean_obs_", cpnX)] <= ss.poss.sc.cpnX.low, ])
    for (samp in samps.poss.less.than.sc.cnpX.low.vec) {
      
      samp.max.CCF      <-  abs(ss.poss.sc.cpnX.low - s.samps.int.ob.xp.cpn.df[samp, paste0("mean_obs_", cpnX)])
      raw.diff.from.int.cpnX.pval <- seg.cpnX.diff.from.exp.list[[samp]][[paste0(cpnX, "_pval_vs_round_int_cpn")]]
      diff.from.int.flag <- raw.diff.from.int.cpnX.pval < (1 - confidence.interval.threshold)

      diff.from.other.chrs.flag <- TRUE
      if (length(seg.cpnX.diff.from.exp.list[[samp]][[paste0(cpnX, "_chrs_w_match_int_cpn_vec")]]) > min.num.other.chrs) { 

        diff.from.other.chrs.flag <- seg.cpnX.diff.from.exp.list[[samp]][[paste0(cpnX, "_pval_vs_oth_cpn_match_chrs_raw_cpn")]] < (1 - confidence.interval.threshold)

      }#end if (length(seg.cpnX.diff.from.exp.list[[samp]][[paste0(cpnX, "_chrs_w_match_int_cpn_vec")]]) > min.num.other.chrs) { 
      
      if (diff.from.int.flag == TRUE & diff.from.other.chrs.flag == TRUE & samp.max.CCF >= min.CCF.to.consider.val) {

        two.value.diff.flag <- FALSE

      }#end if (diff.from.int.flag == TRUE & diff.from.other.chrs.flag == TRUE & samp.max.CCF >= min.CCF.to.consider.val) {

    }#end for (samp in samps.poss.less.than.sc.cnpX.low.vec) {

  }#end if (abs(min(top.int.mean.cpnX - mean.cpnX.vec)) < abs(max(low.int.mean.cpnX - mean.cpnX.vec))) {
  #Todo, make this return a single value (split into two functions ideally)
  return(two.value.diff.flag)

}


chooseLowHighSubclonalCPNWithOnlyTwoPossibilities <- function(s.samps.int.ob.xp.cpn.df,
                                                              cpnX,
                                                              all.clonal.flag,
                                                              two.value.diff.flag) {
  samps.mean.cpnX.vec <- s.samps.int.ob.xp.cpn.df[, paste0("mean_obs_", cpnX)]
  # check whether only one state when rounding
  if (length(unique(round(samps.mean.cpnX.vec))) == 1) {

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
                                                                  cpnX,
                                                                  all.clonal.flag,
                                                                  two.value.diff.flag,
                                                                  min.CCF.to.consider.val) {

  samps.mean.cpnX.vec <- s.samps.int.ob.xp.cpn.df[, paste0("mean_obs_", cpnX)]
  #check whether there is ony one state when rounding, can this ever be called?
  if (length(unique(round(samps.mean.cpnX.vec))) == 1) {

        if (all.clonal.flag == TRUE) {

          ss.poss.sc.cpnX.low =  round(samps.mean.cpnX.vec)
          ss.poss.sc.cpnX.high = round(samps.mean.cpnX.vec)

        } else {#end if (all.clonal.flag == TRUE) {

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
                        #If you have three samples one wth cpn 1 and cpn4 and cpn8, cpn4 would be completely subclonal, not right.
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
    #Todo, check his with Nicky, should be added or subtracted to floor?
    ss.poss.sc.cpnX.low =  min(unique(floor(samps.mean.cpnX.vec + min.CCF.to.consider.val)))
    ss.poss.sc.cpnX.high = max(unique(ceiling(samps.mean.cpnX.vec + min.CCF.to.consider.val)))

  }#end if (length(unique(round(samps.mean.cpnX.vec))) > 2) {

    return(c(ss.poss.sc.cpnX.low, ss.poss.sc.cpnX.high))

}