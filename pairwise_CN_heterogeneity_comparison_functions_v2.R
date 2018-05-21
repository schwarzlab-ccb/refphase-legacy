
calculateSNPsCopyNumberAB <- function(logr.vec, baf.vec, tumor.ploidy, tumor.purity, gamma.val) {

  cpnA <- (tumor.purity - 1 + baf.vec * 2 ^ (logr.vec / gamma.val) * ((1 - tumor.purity) * 2 + tumor.purity * tumor.ploidy)) / tumor.purity
  cpnB <- (tumor.purity - 1 - (baf.vec - 1) * 2 ^ (logr.vec / gamma.val) * ((1 - tumor.purity) * 2 + tumor.purity * tumor.ploidy)) / tumor.purity
  #print("Returning from funCalcN_withBAF")

  return(cbind(cpnA , cpnB))
  
}

#
tTestWithErrorHandling <- function(x, y = NULL, alternative = c("two.sided", "less", "greater"),
                                   mu = 0, paired = FALSE, var.equal = FALSE, conf.level = 0.95, errors.not.permitted=""){
  
  test.results <- tryCatch({
    return(t.test(x=x, y=y, alternative=alternative, mu=mu, paired=paired, var.equal=var.equal, conf.level=conf.level))},
    error = function(err) {
      message(err)
      write(paste(err$message, "\n"), stderr())
      write("\n", stderr())
      if(grepl("NULL", err$message) & "NULL"%in% errors.not.permitted){
        stop("t.test ERROR with a NULL")
      } else if (grepl("enough", err$message) & "enough"%in% errors.not.permitted){
        stop("t.test ERROR with not enough")
      } else if (grepl("constant", err$message) & "constant" %in% errors.not.permitted){
        stop("t.test ERROR with constant data")
      } else if(grepl("enough|constant|NULL", err$message) & "other" %in% errors.not.permitted){
        stop(err$message)
      }
      # cat("\n")
      # paste("Input for t test")
      # print("x")
      # print(x)
      # print("y")
      # print(y)
      # print("mu")
      # print(mu)
      # print("paired")
      # print(paired)
      # print("conf.level")
      # print(conf.level)
      test.results <- list()
      test.results[["stat"]] <- 1
      test.results[["conf.int"]] <- c(NA, NA)
      test.results[["p.value"]] <- 1
      return(test.results)
    })
  return(test.results)
}

initialiseSampOneTwoSegPairwiseComparisonList <- function(samp.one, samp.two, samp.one.two.string){

     output.list <- list()

     output.list[[samp.one]]$chr_vec <- NA
     output.list[[samp.one]]$pos_vec <- NA
     output.list[[samp.two]]$chr_vec <- NA
     output.list[[samp.two]]$pos_vec <- NA

     output.list[[samp.one]]$ph_cpnA_vec <- NA
     output.list[[samp.one]]$ph_cpnB_vec <- NA
     output.list[[samp.one]]$mean_ph_cpnA <- NA
     output.list[[samp.one]]$mean_ph_cpnB <- NA
     output.list[[samp.one]]$un_ph_cpn_total <- NA
     output.list[[samp.one]]$un_ph_cpn_total_mean <- NA

     output.list[[samp.two]]$ph_cpnA_vec <- NA
     output.list[[samp.two]]$ph_cpnB_vec <- NA
     output.list[[samp.two]]$mean_ph_cpnA <- NA
     output.list[[samp.two]]$mean_ph_cpnB <- NA
     output.list[[samp.two]]$un_ph_cpn_total <- NA
     output.list[[samp.two]]$un_ph_cpn_total_mean <- NA


     output.list[[samp.one.two.string]]$un_ph_cpnA_cpnB_het_comb_paired_pval    <- NA
     output.list[[samp.one.two.string]]$un_ph_cpnA_cpnB_het_comb_unpaired_pval  <- NA
     output.list[[samp.one.two.string]]$un_ph_combinedUnPairedPval      <- NA
     output.list[[samp.one.two.string]]$ph_cpnA_het_paired_pval        <- NA
     output.list[[samp.one.two.string]]$ph_cpnB_het_paired_pval        <- NA
     output.list[[samp.one.two.string]]$ph_cpnA_het_unpaired_pval      <- NA
     output.list[[samp.one.two.string]]$ph_cpnB_het_unpaired_pval      <- NA
     output.list[[samp.one.two.string]]$ph_cpnA_cpnB_het_comb_paired_pval   <- NA
     output.list[[samp.one.two.string]]$ph_cpnA_cpnB_het_comb_unpaired_pval <- NA
     output.list[[samp.one.two.string]]$ph_tot_cpn_paired_pval <- NA 
     output.list[[samp.one.two.string]]$ph_cpnA_diff_by_size_het          <- NA
     output.list[[samp.one.two.string]]$ph_cpnB_diff_by_size_het          <- NA
     output.list[[samp.one.two.string]]$comb_cpnA_cpnB_diff_by_size_het     <- NA
     output.list[[samp.one.two.string]]$MSAI     <- NA

  
     return(output.list)
}

calcExpectedLogR <- function(exp.intA, exp.intB, tumor.ploidy, tumor.purity) {

  exp.logr <- log2((2 * (1 - tumor.purity) + tumor.purity * (exp.intA + exp.intB)) / (2 * (1 - tumor.purity) + (tumor.purity * tumor.ploidy)))

  return(exp.logr)

}#end calcExpectedLogR <- function(exp.intA, exp.intB, tumor.ploidy, tumor.purity) {


calcExpectedBAF <- function(exp.intA, exp.intB, tumor.purity) {

  exp.BAF <- (1 - tumor.purity + (tumor.purity * exp.intB)) / (2 - 2 * tumor.purity + tumor.purity * (exp.intA + exp.intB))

  return(exp.BAF)

}#end calcExpectedBAF <- function(exp.intA, exp.intB, tumor.purity)



checkSegGainOrLossSignificance <- function(ss.baf.lr.baf.df, samp.to.check, ploidy.est.vec, purity.est.vec, threshold, cpn.change.type="gain", min.num.for.logr.test = 3) {

  min.thresh.raw.exp.samp       <-  ploidy.est.vec[samp.to.check] * (2 ^ (threshold))
  min.thresh.logr.exp.for.samp   <-  calcExpectedLogR(min.thresh.raw.exp.samp, 0, ploidy.est.vec[samp.to.check], purity.est.vec[samp.to.check])
  median.seg.logr                  <- median(ss.baf.lr.baf.df$logR)
  
  if (nrow(ss.baf.lr.baf.df) >= min.num.for.logr.test) {

     gen.t.test.results        <- tTestWithErrorHandling(ss.baf.lr.baf.df$logR, mu = min.thresh.logr.exp.for.samp)

    if (cpn.change.type=="gain") {

          specific.t.test.alt <- "greater"

    } else if (cpn.change.type=="loss") {

          specific.t.test.alt <- "less"

    }

    specific.t.test.result <- tTestWithErrorHandling(ss.baf.lr.baf.df$logR, mu = min.thresh.logr.exp.for.samp, alternative = specific.t.test.alt)

    low.obs.logr.conf.int  <- specific.t.test.result$conf.int[1]
    high.obs.logr.conf.int <- specific.t.test.result$conf.int[2]
    gen.t.test.pval     <- gen.t.test.results$p.value
    alt.t.test.pval     <- specific.t.test.result$p.value
    
  }#end if (nrow(ss.baf.lr.baf.df) >= 3) {

  #Todo, refactor this magic number...
  if (nrow(ss.baf.lr.baf.df) < 3) {

    low.obs.logr.conf.int  <- NA
    high.obs.logr.conf.int <- NA
    gen.t.test.pval         <- NA
    alt.t.test.pval     <- NA

  }#end f (nrow(ss.baf.lr.baf.df) < 3) {
  
  out <- cbind(min.thresh.raw.exp.samp,
               min.thresh.logr.exp.for.samp,
               median.seg.logr,
               low.obs.logr.conf.int,
               high.obs.logr.conf.int,
               gen.t.test.pval,
               alt.t.test.pval)  
  return(out)
  
}#end checkSegGainOrLossSignificance <- function(ss.baf.lr.baf.df, samp.to.check, ploidy.est.vec, purity.est.vec, threshold) {






calculatePairwiseLogrCNHeterogeneity <- function(samp.one, samp.two, asas.both.lr.baf.lists, asas.lr.and.or.baf.lists, min.cons.segs.df, ploidy.est.vec, purity.est.vec) {





}



  ##Todo, check fif the final segmentation functions can also use this instead of calling the check GainSegloss repeatedly
##Params:
##
##Returns:
##
getSampLogRSegCPNEventVec <- function(samp, samp.seg.ao.lr.baf.df, purity.est.vec, ploidy.est.vec, gain.threshold, amp.threshold, loss.threshold, deep.loss.threshold) {

  cpn.event.names.vec <- c('median_seg_logr','low_ci_logr_vs_exp_logr','high_ci_logr_vs_exp_logr'
                                     ,'min_exp_logr_gain','gain_pval'
                                     ,'min_exp_logr_amp','amp_pval'
                                     ,'min_exp_logr_loss','loss_pval'
                                     ,'min_exp_logr_deep_loss','deep_loss_pval')

  gain.samp.df <- checkSegGainOrLossSignificance(samp.seg.ao.lr.baf.df, samp, ploidy.est.vec, purity.est.vec, gain.threshold, "gain")
  amp.samp.df <- checkSegGainOrLossSignificance(samp.seg.ao.lr.baf.df, samp, ploidy.est.vec, purity.est.vec, amp.threshold, "gain")
  loss.samp.df <- checkSegGainOrLossSignificance(samp.seg.ao.lr.baf.df, samp, ploidy.est.vec, purity.est.vec, loss.threshold, "loss")
  deep.loss.samp.df <- checkSegGainOrLossSignificance(samp.seg.ao.lr.baf.df, samp, ploidy.est.vec, purity.est.vec, deep.loss.threshold, "loss")
  samp.one.cpn.event.vec <- c(gain.samp.df[, 3:5, drop = FALSE]
                           ,gain.samp.df[, 2, drop = FALSE], gain.samp.df[, 7, drop = FALSE]
                           ,amp.samp.df[, 2,drop = FALSE], amp.samp.df[, 7, drop = FALSE]
                           ,loss.samp.df[, 2, drop = FALSE], loss.samp.df[, 7, drop = FALSE]
                           ,deep.loss.samp.df[, 2, drop = FALSE], deep.loss.samp.df[, 7, drop = FALSE])
  names(samp.cpn.event.vec) <- cpn.event.names.vec

  return(samp.cpn.event.vec)

}#end getSampLogRSegCPNEventVec (..)



##
##Params:
##
##Returns:
##
testPairOfSamplesSegForMSAIWithCPN <- function(seg.name, samp.one, samp.two, all.samp.seg.phased.raw.cpn.list) {

  phased.s.one.s.cpnA.cpnB.list <- all.samp.seg.phased.raw.cpn.list[[seg.name]][[samp.one]]
  phased.s.one.s.cpnA.cpnB.df <- as.matrix(cbind(phased.s.one.s.cpnA.cpnB.list$ph_cpnA_vec, phased.s.one.s.cpnA.cpnB.list$ph_cpnB_vec))
  phased.s.two.s.cpnA.cpnB.list <- all.samp.seg.phased.raw.cpn.list[[seg.name]][[samp.two]]
  phased.s.two.s.cpnA.cpnB.df <- as.matrix(cbind(phased.s.two.s.cpnA.cpnB.list$ph_cpnA_vec, phased.s.two.s.cpnA.cpnB.list$ph_cpnB_vec))

  if (tTestWithErrorHandling(phased.s.one.s.cpnA.cpnB.df[, 1], phased.s.one.s.cpnA.cpnB.df[, 2], paired=TRUE)$p.val < 0.001 &
       tTestWithErrorHandling(phased.s.two.s.cpnA.cpnB.df[, 1], phased.s.two.s.cpnA.cpnB.df[, 2], paired=TRUE)$p.val < 0.001) {

       if( (mean(phased.s.one.s.cpnA.cpnB.df[, 1] - phased.s.one.s.cpnA.cpnB.df[, 2]) > 0 & 
            mean(phased.s.two.s.cpnA.cpnB.df[, 1] - phased.s.two.s.cpnA.cpnB.df[, 2]) < 0) | 
           (mean(phased.s.one.s.cpnA.cpnB.df[, 2] - phased.s.one.s.cpnA.cpnB.df[, 1]) < 0 & 
            mean(phased.s.two.s.cpnA.cpnB.df[, 2] - phased.s.two.s.cpnA.cpnB.df[, 1]) > 0) ){

            #s.one.two.comparison[[seg.name]][[samp.one.two.string]]$MSAI <- 1
            return(TRUE)

       }

  }#end if (tTestWithErrorHandling(phased.s.one.s.cpnA.cpnB.df[, 1], phased.s.one.s.cpnA.cpnB.df[, 2], paired=TRUE)$p.val < 0.001 &
       #tTestWithErrorHandling(phased.s.two.s.cpnA.cpnB.df[, 1], phased.s.two.s.cpnA.cpnB.df[,2], paired=TRUE)$p.val < 0.001) {

  return(FALSE)

}#testPairOfSamplesSegForMSAIWithCPN function(....)




##Note for now we're missing heterozygous SNPs that couldn't be phased...
##Params:
##
##Returns:
##
createAllSegInputForHetComparisons <- function(seg.name,
                                              samp.one,
                                              samp.two,
                                              all.samp.seg.phased.raw.cpn.list,
                                              asas.both.lr.baf.lists,
                                              all.samp.non.hz.logr.and.tot.cpn.df.list,
                                              min.num.hz.snp.phasing) {

  seg.samp.one.two.comp.data.list <- list()
  seg.samp.one.two.comp.data.list[[samp.one]] <- list()
  seg.samp.one.two.comp.data.list[[samp.two]] <- list()
  data.present.vec <- c(FALSE, FALSE, FALSE)
  names(data.present.vec) <- ("hz_phased", "non_hz", "hz_unphased")
  #Any unphased cpnA cpnB if it's present, ie there were hz positions but fewer than than min.num.hz.snp.phasing so they weren't phased
  # s.one.s.b.lr.baf.df <- asas.both.lr.baf.lists[[samp.one]][[seg.name]]
  # s.two.s.b.lr.baf.df <- asas.both.lr.baf.lists[[samp.two]][[seg.name]]

  seg.name.split.vec <- unlist(strsplit(seg.name, split=":"))
  samp.one.non.hz.df <- all.samp.non.hz.logr.and.tot.cpn.df.list[[samp.one]]
  samp.seg.one.non.hz.df <- samp.one.non.hz.df[samp.one.non.hz.df$chrom == as.numeric(seg.name.split.vec[1]) &
                                               samp.one.non.hz.df$pos >= as.numeric(seg.name.split.vec[2]) &
                                               samp.one.non.hz.df$pos <= as.numeric(seg.name.split.vec[3])
                                               ,]
  samp.two.non.hz.df <- all.samp.non.hz.logr.and.tot.cpn.df.list[[samp.two]]
  samp.seg.two.non.hz.df <- samp.two.non.hz.df[samp.two.non.hz.df$chrom == as.numeric(seg.name.split.vec[1]) &
                                               samp.two.non.hz.df$pos >= as.numeric(seg.name.split.vec[2]) &
                                               samp.two.non.hz.df$pos <= as.numeric(seg.name.split.vec[3])
                                               ,]

  print("head(samp.seg.one.non.hz.df)")
  print(head(samp.seg.one.non.hz.df))
  if (nrow(samp.seg.one.non.hz.df) > 0) {

    data.present.vec["non_hz"] <- TRUE

    seg.samp.one.two.comp.data.list[[samp.one]]$non_hz_chr_vec <- samp.one.non.hz.df$chrom
    seg.samp.one.two.comp.data.list[[samp.one]]$non_hz_pos_vec <- samp.one.non.hz.df$pos

    seg.samp.one.two.comp.data.list[[samp.one]]$non_hz_logr_vec <- samp.one.non.hz.df$logR
    seg.samp.one.two.comp.data.list[[samp.one]]$non_hz_tot_raw_cpn_vec <- samp.one.non.hz.df$tot_raw_un_ph_cpn

    seg.samp.one.two.comp.data.list[[samp.two]]$non_hz_chr_vec <- samp.two.non.hz.df$chrom
    seg.samp.one.two.comp.data.list[[samp.two]]$non_hz_pos_vec <- samp.two.non.hz.df$pos

    seg.samp.one.two.comp.data.list[[samp.two]]$non_hz_logr_vec <- samp.two.non.hz.df$logR
    seg.samp.one.two.comp.data.list[[samp.two]]$non_hz_tot_raw_cpn_vec <- samp.two.non.hz.df$tot_raw_un_ph_cpn


  }#end   

  phased.s.one.s.cpnA.cpnB.list <- all.samp.seg.phased.raw.cpn.list[[seg.name]][[samp.one]]
  print("str(phased.s.one.s.cpnA.cpnB.list)")
  print(str(phased.s.one.s.cpnA.cpnB.list, max =1))
  num.non.na.ph <- length(phased.s.one.s.cpnA.cpnB.list$ph_cpnA_vec[ ! is.na(phased.s.one.s.cpnA.cpnB.list$ph_cpnA_vec)])
  print("num.non.na.ph")
  print(num.non.na.ph)
  if ( num.non.na.ph >= min.num.hz.snp.phasing) {

    data.present.vec["hz_phased"] <- TRUE

    phased.s.one.s.cpnA.cpnB.df <- as.matrix(cbind(phased.s.one.s.cpnA.cpnB.list$ph_cpnA_vec, phased.s.one.s.cpnA.cpnB.list$ph_cpnB_vec))
    phased.s.two.s.cpnA.cpnB.list <- all.samp.seg.phased.raw.cpn.list[[seg.name]][[samp.two]]
    phased.s.two.s.cpnA.cpnB.df <- as.matrix(cbind(phased.s.two.s.cpnA.cpnB.list$ph_cpnA_vec, phased.s.two.s.cpnA.cpnB.list$ph_cpnB_vec))

    seg.samp.one.two.comp.data.list[[samp.one]]$ph_logr_vec <- phased.s.one.s.cpnA.cpnB.list$logr_vec
    seg.samp.one.two.comp.data.list[[samp.one]]$ph_cpnA_vec <- as.numeric(phased.s.one.s.cpnA.cpnB.df[, 1])
    seg.samp.one.two.comp.data.list[[samp.one]]$ph_cpnB_vec <- as.numeric(phased.s.one.s.cpnA.cpnB.df[, 2])
    seg.samp.one.two.comp.data.list[[samp.one]]$mean_ph_cpnA <- mean(as.numeric(phased.s.one.s.cpnA.cpnB.df[, 1]))
    seg.samp.one.two.comp.data.list[[samp.one]]$mean_ph_cpnB <- mean(as.numeric(phased.s.one.s.cpnA.cpnB.df[, 2]))

    seg.samp.one.two.comp.data.list[[samp.one]]$ph_logr_vec <- phased.s.two.s.cpnA.cpnB.list$logr_vec
    seg.samp.one.two.comp.data.list[[samp.two]]$ph_cpnA_vec <- as.numeric(phased.s.two.s.cpnA.cpnB.df[, 1])
    seg.samp.one.two.comp.data.list[[samp.two]]$ph_cpnB_vec <- as.numeric(phased.s.two.s.cpnA.cpnB.df[, 2])
    seg.samp.one.two.comp.data.list[[samp.two]]$mean_ph_cpnA <- mean(as.numeric(phased.s.two.s.cpnA.cpnB.df[, 1]))
    seg.samp.one.two.comp.data.list[[samp.two]]$mean_ph_cpnB <- mean(as.numeric(phased.s.two.s.cpnA.cpnB.df[, 2]))

  } else if (num.non.na.ph < min.num.hz.snp.phasing & num.non.na.ph > 0 ) {

    data.present.vec["hz_unphased"] <- TRUE
    #not sure about this add them to
    seg.samp.one.two.comp.data.list[[samp.one]]$un_ph_hz_logr_vec <- phased.s.one.s.cpnA.cpnB.list$logr_vec
    seg.samp.one.two.comp.data.list[[samp.one]]$un_ph_hz_tot_raw_cpn_vec <- phased.s.one.s.cpnA.cpnB.list$tot_raw_un_ph_cpn

    seg.samp.one.two.comp.data.list[[samp.two]]$un_ph_hz_logr_vec <- phased.s.one.s.cpnA.cpnB.list$logr_vec
    seg.samp.one.two.comp.data.list[[samp.two]]$un_ph_hz_tot_raw_cpn_vec <- phased.s.one.s.cpnA.cpnB.list$tot_raw_un_ph_cpn


  } 
  
  seg.samp.one.two.comp.data.list$data_present_vec <- data.present.vec

  return(seg.samp.one.two.comp.data.list)


}#end createAllSegInputForHetComparisons



##
##Params:
##
##Returns:
##
getSubsetSegInfoForHetComparison <- function(analysis.name, samp.one, samp.two, seg.samp.one.two.comp.data.list){


  data.present.vec <- seg.samp.one.two.comp.data.list$data_present_vec
  data.for.analysis.list <- list()

  if (analysis.name == "hz_phased" & data.present.vec["hz_phased"] == TRUE) {
    #Analysis
    samp.one.logr.vec <- seg.samp.one.two.comp.data.list[[samp.one]]$ph_logr_vec
    samp.two.logr.vec <- seg.samp.one.two.comp.data.list[[samp.two]]$ph_logr_vec
    samp.one.tot.raw.cpn.vec <- seg.samp.one.two.comp.data.list[[samp.one]]$non_hz_tot_raw_cpn_vec
    samp.two.tot.raw.cpn.vec <- seg.samp.one.two.comp.data.list[[samp.two]]$non_hz_tot_raw_cpn_vec
    samp.one.ph.cpnA.vec <- seg.samp.one.two.comp.data.list[[samp.one]]$ph_cpnA_vec
    samp.one.ph.cpnB.vec <- seg.samp.one.two.comp.data.list[[samp.one]]$ph_cpnB_vec
    samp.two.ph.cpnA.vec <- seg.samp.one.two.comp.data.list[[samp.two]]$ph_cpnA_vec
    samp.two.ph.cpnB.vec <- seg.samp.one.two.comp.data.list[[samp.two]]$ph_cpnB_vec


  } else if (analysis.name == "hz_phased" & data.present.vec["hz_phased"] == FALSE) {

    return(NULL)

  }

  if (analysis.name == "non_hz" & data.present.vec["non_hz"] == TRUE) {

    samp.one.logr.vec <- seg.samp.one.two.comp.data.list[[samp.one]]$non_hz_logr_vec
    samp.two.logr.vec <- seg.samp.one.two.comp.data.list[[samp.two]]$non_hz_logr_vec
    samp.one.tot.raw.cpn.vec <- seg.samp.one.two.comp.data.list[[samp.one]]$non_hz_tot_raw_cpn_vec
    samp.two.tot.raw.cpn.vec <- seg.samp.one.two.comp.data.list[[samp.two]]$non_hz_tot_raw_cpn_vec
    samp.one.ph.cpnA.vec <- NA
    samp.one.ph.cpnB.vec <- NA
    samp.two.ph.cpnA.vec <- NA
    samp.two.ph.cpnB.vec <- NA

  } else if (analysis.name == "non_hz" & data.present.vec["non_hz"] == FALSE) {

    return(NULL)

  }

  if (analysis.name == "combined" & data.present.vec["non_hz"] == TRUE & data.present.vec["hz_phased"] == TRUE) {
    #Make sure we add in the unphased heterozygous snps here if they are present

    samp.one.logr.vec <- c(seg.samp.one.two.comp.data.list[[samp.one]]$ph_logr_vec, seg.samp.one.two.comp.data.list[[samp.one]]$non_hz_logr_vec)
    samp.two.logr.vec <- c(seg.samp.one.two.comp.data.list[[samp.two]]$ph_logr_vec, seg.samp.one.two.comp.data.list[[samp.two]]$non_hz_logr_vec)
    samp.one.tot.raw.cpn.vec <- c(seg.samp.one.two.comp.data.list[[samp.one]]$non_hz_tot_raw_cpn_vec, seg.samp.one.two.comp.data.list[[samp.one]]$non_hz_tot_raw_cpn_vec) 
    samp.two.tot.raw.cpn.vec <- c(seg.samp.one.two.comp.data.list[[samp.two]]$non_hz_tot_raw_cpn_vec, seg.samp.one.two.comp.data.list[[samp.one]]$non_hz_tot_raw_cpn_vec
    samp.one.ph.cpnA.vec <- seg.samp.one.two.comp.data.list[[samp.one]]$ph_cpnA_vec
    samp.one.ph.cpnB.vec <- seg.samp.one.two.comp.data.list[[samp.one]]$ph_cpnB_vec
    samp.two.ph.cpnA.vec <- seg.samp.one.two.comp.data.list[[samp.two]]$ph_cpnA_vec
    samp.two.ph.cpnB.vec <- seg.samp.one.two.comp.data.list[[samp.two]]$ph_cpnB_vec

    if (data.present.vec["hz_unphased"] == TRUE) {

      samp.one.logr.vec <- c(samp.one.logr.vec, seg.samp.one.two.comp.data.list[[samp.one]]$un_ph_hz_logr_vec)
      samp.two.logr.vec <- c(samp.two.logr.vec, seg.samp.one.two.comp.data.list[[samp.two]]$un_ph_hz_logr_vec)
      samp.one.tot.raw.cpn.vec <- c(samp.one.tot.raw.cpn.vec, seg.samp.one.two.comp.data.list[[samp.one]]$un_ph_hz_tot_raw_cpn_vec)
      samp.two.tot.raw.cpn.vec <- c(samp.two.tot.raw.cpn.vec, seg.samp.one.two.comp.data.list[[samp.two]]$un_ph_hz_tot_raw_cpn_vec)

    }#end if (data.present.vec["hz_unphased"] == TRUE) {

  } else {

    return(NULL)

  }

  data.for.analysis.list$samp_one_logr_vec <- samp.one.logr.vec
  data.for.analysis.list$samp_two_logr_vec <- samp.two.logr.vec
  data.for.analysis.list$samp_one_tot_raw_cpn_vec <- samp.one.tot.raw.cpn.vec
  data.for.analysis.list$samp_two_tot_raw_cpn_vec <- samp.two.tot.raw.cpn.vec
  data.for.analysis.list$samp_one_ph_cpnA_vec <- samp.one.ph.cpnA.vec
  data.for.analysis.list$samp_one_ph_cpnB_vec <- samp.one.ph.cpnB.vec
  data.for.analysis.list$samp_two_ph_cpnA_vec <- samp.two.ph.cpnA.vec
  data.for.analysis.list$samp_two_ph_cpnB_vec <- samp.two.ph.cpnB.vec

  return(data.for.analysis.list)
}


##
##Params:
##
##Returns:
##
performSegPairwiseComparisons <- function(seg.name, samp.one, samp.two, seg.samp.one.two.comp.data.list) {

  #So we need to decide how we're going to perform all these comparisons, there are three types:
  #logr alone comparison: possible for all segments
  #raw ph cpn comparison: possible only for segments thad had greater than MIN_NUM_HZ_SNP_PHASING in them
  #raw tot cpn: possible for all segments
  #processed call comparisons: currently possible only for segments with phased copy number
  #Also there separate and combined logr and total cpn calls we can make that might come up as being different,
  #We need to separate these out early on and then run the same functions twice over them.

  #First do the hz alone if present
  analysis.type.vec <- c("hz_phased", "non_hz", "combined")

  data.present.vec <- seg.samp.one.two.comp.data.list$data_present_vec

  for ( analysis.name %in% analysis.type.vec ) {



  }


  #If we can perform total copy number comparisons, do so here:

  #logr comparison data



  return(object)

}









# seg.name.split.vec <- strsplit(chr.start.end.string.vec[i], split=":")
# seg.size <- as.numeric(unlist(seg.name.split.vec)[3]) - as.numeric(unlist(seg.name.split.vec)[2])

# pw.CN.het.score.vec <- c()
# seg.size.vec <- c()

# seg.cpnA.het.size.by.diff <- 0            
# seg.cpnB.het.size.by.diff <- 0


# #Todo make the 0.001 a parameter..
# if (seg.paired.p.vals.vec[1] < 0.001) {

#      cpnA.diff.amplitude <- abs(mean(as.numeric(phased.s.one.s.cpnA.cpnB.df[, 1])) - mean(as.numeric(phased.s.two.s.cpnA.cpnB.df[, 1])))
#      seg.cpnA.het.size.by.diff <- seg.size * cpnA.diff.amplitude

# }#end if (seg.paired.p.vals.vec[1] < 0.001) {

# if (seg.paired.p.vals.vec[2] < 0.001) {

#      cpnB.diff.amplitude <- abs(mean(as.numeric(phased.s.one.s.cpnA.cpnB.df[, 2])) - mean(as.numeric(phased.s.two.s.cpnA.cpnB.df[, 2])))
#      seg.cpnB.het.size.by.diff <- seg.size * cpnB.diff.amplitude

# }#end if (seg.paired.p.vals.vec[2] < 0.001) { 


#   seg.tot.cpnA.cpnB.het <- seg.cpnA.het.size.by.diff + seg.cpnB.het.size.by.diff
#   pw.CN.het.score.vec <- c(pw.CN.het.score.vec, seg.tot.cpnA.cpnB.het)
#   seg.size.vec <- c(seg.size.vec, seg.size)


#   s.one.two.comparison[[seg.name]][[samp.one.two.string]]$ph_cpnA_diff_by_size_het <- seg.cpnA.het.size.by.diff
#   s.one.two.comparison[[seg.name]][[samp.one.two.string]]$ph_cpnB_diff_by_size_het <- seg.cpnB.het.size.by.diff
#   s.one.two.comparison[[seg.name]][[samp.one.two.string]]$comb_cpnA_cpnB_diff_by_size_het <- seg.cpnA.het.size.by.diff + seg.cpnB.het.size.by.diff




# #So now do the p value tests:

# seg.paired.p.vals.vec <- c(tTestWithErrorHandling(s.one.s.cpnA.cpnB.df[, 1], s.two.s.cpnA.cpnB.df[, 1], paired=TRUE)$p.val,
#                           tTestWithErrorHandling(s.one.s.cpnA.cpnB.df[, 2], s.two.s.cpnA.cpnB.df[, 2], paired=TRUE)$p.val)
# combined.cpnA.cpnB.paired.pval <- seg.paired.p.vals.vec[1]
# seg.unpaired.p.vals.vec <- c(tTestWithErrorHandling(s.one.s.cpnA.cpnB.df[, 1], s.two.s.cpnA.cpnB.df[,1])$p.val,
#                               tTestWithErrorHandling(s.one.s.cpnA.cpnB.df[, 2], s.two.s.cpnA.cpnB.df[,2])$p.val)
# combined.cpnA.cpnB.unpaired.pval <- seg.unpaired.p.vals.vec[1]#pchisq((sum(log(PValsUnPaired))*-2), df=length(PValsUnPaired)*2, lower.tail=F)

# # also get a p-value that relates to both samples combined
# total.cpn.paired.pval <- tTestWithErrorHandling(c(phased.s.one.s.cpnA.cpnB.df[, 1] + phased.s.one.s.cpnA.cpnB.df[, 2]),
#                                                 c(phased.s.two.s.cpnA.cpnB.df[, 1] + phased.s.two.s.cpnA.cpnB.df[, 2]), paired = TRUE)$p.val

# #Froms the phased copy number

# seg.paired.p.vals.vec <- c(tTestWithErrorHandling(phased.s.one.s.cpnA.cpnB.df[, 1],
#                                                   phased.s.two.s.cpnA.cpnB.df[, 1], paired = TRUE)$p.val,
#                            tTestWithErrorHandling(phased.s.one.s.cpnA.cpnB.df[, 2],
#                                                   phased.s.two.s.cpnA.cpnB.df[, 2], paired = TRUE)$p.val)

# seg.comb.paired.p.val <- pchisq((sum(log(seg.paired.p.vals.vec)) * - 2), df = length(seg.paired.p.vals.vec) * 2, lower.tail = FALSE)

# seg.unpaired.p.vals.vec <- c(tTestWithErrorHandling(phased.s.one.s.cpnA.cpnB.df[, 1],
#                                                     phased.s.two.s.cpnA.cpnB.df[, 1])$p.val,
#                              tTestWithErrorHandling(phased.s.one.s.cpnA.cpnB.df[, 2],
#                                                     phased.s.two.s.cpnA.cpnB.df[, 2])$p.val)

# seg.comb.unpaired.p.val <- pchisq((sum(log(seg.unpaired.p.vals.vec)) * - 2), df = length(seg.unpaired.p.vals.vec) * 2, lower.tail = FALSE)

# # also get a p-value that relates to both samples combined
# total.cpn.paired.pval <- tTestWithErrorHandling(c(phased.s.one.s.cpnA.cpnB.df[, 1] + phased.s.one.s.cpnA.cpnB.df[, 2]),
#                                                 c(phased.s.two.s.cpnA.cpnB.df[, 1] + phased.s.two.s.cpnA.cpnB.df[, 2]), paired = TRUE)$p.val


# s.one.two.comparison[[seg.name]][[samp.one.two.string]]$un_ph_cpnA_cpnB_het_comb_paired_pval    <- combined.cpnA.cpnB.paired.pval
# s.one.two.comparison[[seg.name]][[samp.one.two.string]]$un_ph_cpnA_cpnB_het_comb_unpaired_pval  <- combined.cpnA.cpnB.unpaired.pval
# s.one.two.comparison[[seg.name]][[samp.one.two.string]]$un_ph_combinedUnPairedPval      <- total.cpn.paired.pval

# s.one.two.comparison[[seg.name]][[samp.one.two.string]]$un_ph_cpnA_cpnB_het_comb_paired_pval    <- NA
# s.one.two.comparison[[seg.name]][[samp.one.two.string]]$un_ph_cpnA_cpnB_het_comb_unpaired_pval  <- NA

# s.one.two.comparison[[seg.name]][[samp.one.two.string]]$ph_cpnA_het_paired_pval <- seg.paired.p.vals.vec[1]
# s.one.two.comparison[[seg.name]][[samp.one.two.string]]$ph_cpnB_het_paired_pval <- seg.paired.p.vals.vec[2]
# s.one.two.comparison[[seg.name]][[samp.one.two.string]]$ph_cpnA_het_unpaired_pval <- seg.unpaired.p.vals.vec[1]
# s.one.two.comparison[[seg.name]][[samp.one.two.string]]$ph_cpnB_het_unpaired_pval <- seg.unpaired.p.vals.vec[2]
# s.one.two.comparison[[seg.name]][[samp.one.two.string]]$ph_cpnA_cpnB_het_comb_paired_pval <- seg.comb.paired.p.val
# s.one.two.comparison[[seg.name]][[samp.one.two.string]]$ph_cpnA_cpnB_het_comb_unpaired_pval <- seg.comb.unpaired.p.val
# s.one.two.comparison[[seg.name]][[samp.one.two.string]]$ph_tot_cpn_paired_pval <- total.cpn.paired.pval









# #Todo harmonise > or >= for minimum SNP number throughout code.
# #Todo add in 0.001 threshold as a parameter that can be varied
# calculatePairwiseCNHeterogeneity <- function(samp.one,
#                                              samp.two,
#                                              asas.both.lr.baf.lists,
#                                              asas.lr.and.or.baf.lists,
#                                              seg.name.vec,
#                                              ploidy.est.vec,
#                                              purity.est.vec,
#                                              all.seg.phasing.info.list,
#                                              all.samp.seg.phased.raw.cpn.list,
#                                              gamma.val,
#                                              gain.threshold = log2(2.5 / 2),
#                                              loss.threshold = log2(1.5 / 2),
#                                              amp.threshold = log2(4 / 2),
#                                              deep.loss.threshold = log2(0.5 / 2),
#                                              min.num.hz.snp.phasing = 5
#                                              ) {


#   #debug parameters 
#   # samp.one <- "B_MET048_MT-M1"
#   # samp.two <- "B_MET048_T1-R1"


#   #Todo, why are we using seg.name or the index of the seg, choose one throughout code and stick to it 
# 	for ( seg.name in seg.name.vec ) {
#     #This is the sample comparison ppart 
#           samp.one.two.string <- paste(samp.one, samp.two, sep=":")
#           s.one.s.b.lr.baf.df <- asas.both.lr.baf.lists[[samp.one]][[seg.name]]
#           s.two.s.b.lr.baf.df <- asas.both.lr.baf.lists[[samp.two]][[seg.name]]

#           s.one.s.ao.lr.baf.df <- asas.lr.and.or.baf.lists[[samp.one]][[seg.name]]
#           s.two.s.ao.lr.baf.df <- asas.lr.and.or.baf.lists[[samp.two]][[seg.name]]

#           s.one.two.comparison[[seg.name]] <- initialiseSampOneTwoSegPairwiseComparisonList(samp.one, samp.two, samp.one.two.string)

#           print(nrow(s.one.s.b.lr.baf.df))

          
#           if (nrow(s.one.s.b.lr.baf.df) < min.num.hz.snp.phasing) {


#                s.one.s.cpnA.cpnB.df <- calculateSNPsCopyNumberAB(logr.vec = s.one.s.ao.lr.baf.df$logR,
#                                                                  baf.vec = s.one.s.ao.lr.baf.df$BAF,
#                                                                  tumor.ploidy = ploidy.est.vec[samp.one],
#                                                                  tumor.purity = purity.est.vec[samp.one],
#                                                                  gamma.val = gamma.val)
#                s.two.s.cpnA.cpnB.df <- calculateSNPsCopyNumberAB(logr.vec = s.two.s.ao.lr.baf.df$logR,
#                                                                  baf.vec = s.two.s.ao.lr.baf.df$BAF,
#                                                                  tumor.ploidy = ploidy.est.vec[samp.two],
#                                                                  tumor.purity = purity.est.vec[samp.two],
#                                                                  gamma.val = gamma.val)

#                seg.paired.p.vals.vec <- c(tTestWithErrorHandling(s.one.s.cpnA.cpnB.df[, 1], s.two.s.cpnA.cpnB.df[, 1], paired=TRUE)$p.val,
#                                          tTestWithErrorHandling(s.one.s.cpnA.cpnB.df[, 2], s.two.s.cpnA.cpnB.df[, 2], paired=TRUE)$p.val)
#                combined.cpnA.cpnB.paired.pval <- seg.paired.p.vals.vec[1]
#                seg.unpaired.p.vals.vec <- c(tTestWithErrorHandling(s.one.s.cpnA.cpnB.df[, 1], s.two.s.cpnA.cpnB.df[,1])$p.val,
#                                              tTestWithErrorHandling(s.one.s.cpnA.cpnB.df[, 2], s.two.s.cpnA.cpnB.df[,2])$p.val)
#                combined.cpnA.cpnB.unpaired.pval <- seg.unpaired.p.vals.vec[1]#pchisq((sum(log(PValsUnPaired))*-2), df=length(PValsUnPaired)*2, lower.tail=F)
               
#                # also get a p-value that relates to both samples combined
#                total.cpn.paired.pval <- tTestWithErrorHandling(c(phased.s.one.s.cpnA.cpnB.df[, 1] + phased.s.one.s.cpnA.cpnB.df[, 2]),
#                                                                c(phased.s.two.s.cpnA.cpnB.df[, 1] + phased.s.two.s.cpnA.cpnB.df[, 2]), paired = TRUE)$p.val


#                seg.name.split.vec <- strsplit(chr.start.end.string.vec[i], split=":")
#                seg.size <- as.numeric(unlist(seg.name.split.vec)[3]) - as.numeric(unlist(seg.name.split.vec)[2])

#                s.one.two.comparison[[seg.name]][[samp.one]]$chr_vec <- s.one.s.b.lr.baf.df$chrom
#                s.one.two.comparison[[seg.name]][[samp.one]]$pos_vec <- s.one.s.b.lr.baf.df$pos
#                s.one.two.comparison[[seg.name]][[samp.two]]$chr_vec <- s.one.s.b.lr.baf.df$chrom
#                s.one.two.comparison[[seg.name]][[samp.two]]$pos_vec <- s.one.s.b.lr.baf.df$pos

#                s.one.two.comparison[[seg.name]][[samp.one]]$un_ph_cpn_total <- s.one.s.cpnA.cpnB.df[, 1] + s.one.s.cpnA.cpnB.df[, 2]
#                s.one.two.comparison[[seg.name]][[samp.one]]$un_ph_cpn_total_mean <- mean(s.one.s.cpnA.cpnB.df[, 1] + s.one.s.cpnA.cpnB.df[, 2])

#                s.one.two.comparison[[seg.name]][[samp.two]]$un_ph_cpn_total <- s.two.s.cpnA.cpnB.df[, 1] + s.two.s.cpnA.cpnB.df[, 2]
#                s.one.two.comparison[[seg.name]][[samp.two]]$un_ph_cpn_total_mean <- mean(s.two.s.cpnA.cpnB.df[, 1] + s.two.s.cpnA.cpnB.df[, 2])

#                s.one.two.comparison[[seg.name]][[samp.one.two.string]]$un_ph_cpnA_cpnB_het_comb_paired_pval    <- combined.cpnA.cpnB.paired.pval
#                s.one.two.comparison[[seg.name]][[samp.one.two.string]]$un_ph_cpnA_cpnB_het_comb_unpaired_pval  <- combined.cpnA.cpnB.unpaired.pval
#                s.one.two.comparison[[seg.name]][[samp.one.two.string]]$un_ph_combinedUnPairedPval      <- total.cpn.paired.pval


#                samp.one.cpn.event.vec <- getSampLogRSegCPNEventVec(samp.one, s.one.s.ao.lr.baf.df, purity.est.vec, ploidy.est.vec,
#                                                                   gain.threshold, amp.threshold, loss.threshold, deep.loss.threshold)
#                s.one.two.comparison[[seg.name]][[samp.one]]$seg_cpn_event_status_vec <- samp.one.cpn.event.vec

#                samp.two.cpn.event.vec <- getSampLogRSegCPNEventVec(samp.two, s.two.s.ao.lr.baf.df, purity.est.vec, ploidy.est.vec,
#                                                                   gain.threshold, amp.threshold, loss.threshold, deep.loss.threshold)
#                s.one.two.comparison[[seg.name]][[samp.two]]$seg_cpn_event_status_vec <- samp.two.cpn.event.vec

#                next;

#           }#end if (nrow(s.one.s.b.lr.baf.df) < min.num.hz.snp.phasing) {

#           if (nrow(s.one.s.b.lr.baf.df) > min.num.hz.snp.phasing) {

#               phased.s.one.s.cpnA.cpnB.list <- all.samp.seg.phased.raw.cpn.list[[seg.name]][[samp.one]]
#               phased.s.one.s.cpnA.cpnB.df <- as.matrix(cbind(phased.s.one.s.cpnA.cpnB.list$ph_cpnA_vec, phased.s.one.s.cpnA.cpnB.list$ph_cpnB_vec))
#               phased.s.two.s.cpnA.cpnB.list <- all.samp.seg.phased.raw.cpn.list[[seg.name]][[samp.two]]
#               phased.s.two.s.cpnA.cpnB.df <- as.matrix(cbind(phased.s.two.s.cpnA.cpnB.list$ph_cpnA_vec, phased.s.two.s.cpnA.cpnB.list$ph_cpnB_vec))

#                seg.paired.p.vals.vec <- c(tTestWithErrorHandling(phased.s.one.s.cpnA.cpnB.df[, 1],
#                                                                  phased.s.two.s.cpnA.cpnB.df[, 1], paired = TRUE)$p.val,
#                                           tTestWithErrorHandling(phased.s.one.s.cpnA.cpnB.df[, 2],
#                                                                  phased.s.two.s.cpnA.cpnB.df[, 2], paired = TRUE)$p.val)

#                seg.comb.paired.p.val <- pchisq((sum(log(seg.paired.p.vals.vec)) * - 2), df = length(seg.paired.p.vals.vec) * 2, lower.tail = FALSE)

#                seg.unpaired.p.vals.vec <- c(tTestWithErrorHandling(phased.s.one.s.cpnA.cpnB.df[, 1],
#                                                                    phased.s.two.s.cpnA.cpnB.df[, 1])$p.val,
#                                             tTestWithErrorHandling(phased.s.one.s.cpnA.cpnB.df[, 2],
#                                                                    phased.s.two.s.cpnA.cpnB.df[, 2])$p.val)

#                seg.comb.unpaired.p.val <- pchisq((sum(log(seg.unpaired.p.vals.vec)) * - 2), df = length(seg.unpaired.p.vals.vec) * 2, lower.tail = FALSE)

#                # also get a p-value that relates to both samples combined
#                total.cpn.paired.pval <- tTestWithErrorHandling(c(phased.s.one.s.cpnA.cpnB.df[, 1] + phased.s.one.s.cpnA.cpnB.df[, 2]),
#                                                                c(phased.s.two.s.cpnA.cpnB.df[, 1] + phased.s.two.s.cpnA.cpnB.df[, 2]), paired = TRUE)$p.val

#                seg.cpnA.het.size.by.diff <- 0            
#                seg.cpnB.het.size.by.diff <- 0

#                seg.name.split.vec <- strsplit(chr.start.end.string.vec[i], split=":")
#                seg.size <- as.numeric(unlist(seg.name.split.vec)[3]) - as.numeric(unlist(seg.name.split.vec)[2])

#                #Todo make the 0.001 a parameter..
#                if (seg.paired.p.vals.vec[1] < 0.001) {

#                     cpnA.diff.amplitude <- abs(mean(as.numeric(phased.s.one.s.cpnA.cpnB.df[, 1])) - mean(as.numeric(phased.s.two.s.cpnA.cpnB.df[, 1])))
#                     seg.cpnA.het.size.by.diff <- seg.size * cpnA.diff.amplitude

#                }#end if (seg.paired.p.vals.vec[1] < 0.001) {

#                if (seg.paired.p.vals.vec[2] < 0.001) {

#                     cpnB.diff.amplitude <- abs(mean(as.numeric(phased.s.one.s.cpnA.cpnB.df[, 2])) - mean(as.numeric(phased.s.two.s.cpnA.cpnB.df[, 2])))
#                     seg.cpnB.het.size.by.diff <- seg.size * cpnB.diff.amplitude

#                }#end if (seg.paired.p.vals.vec[2] < 0.001) { 

#                seg.tot.cpnA.cpnB.het <- seg.cpnA.het.size.by.diff + seg.cpnB.het.size.by.diff
#                pw.CN.het.score.vec <- c(pw.CN.het.score.vec, seg.tot.cpnA.cpnB.het)
#                seg.size.vec <- c(seg.size.vec, seg.size)

#                s.one.two.comparison[[seg.name]][[samp.one]]$chr_vec <- s.one.s.b.lr.baf.df$chrom
#                s.one.two.comparison[[seg.name]][[samp.one]]$pos_vec <- s.one.s.b.lr.baf.df$pos

#                s.one.two.comparison[[seg.name]][[samp.one]]$ph_cpnA_vec <- as.numeric(phased.s.one.s.cpnA.cpnB.df[, 1])
#                s.one.two.comparison[[seg.name]][[samp.one]]$ph_cpnB_vec <- as.numeric(phased.s.one.s.cpnA.cpnB.df[, 2])
#                s.one.two.comparison[[seg.name]][[samp.one]]$mean_ph_cpnA <- mean(as.numeric(phased.s.one.s.cpnA.cpnB.df[, 1]))
#                s.one.two.comparison[[seg.name]][[samp.one]]$mean_ph_cpnB <- mean(as.numeric(phased.s.one.s.cpnA.cpnB.df[, 2]))
#                s.one.two.comparison[[seg.name]][[samp.one]]$un_ph_cpn_total <- NA
#                s.one.two.comparison[[seg.name]][[samp.one]]$un_ph_cpn_total_mean <- NA

#                s.one.two.comparison[[seg.name]][[samp.two]]$ph_cpnA_vec <- as.numeric(phased.s.two.s.cpnA.cpnB.df[, 1])
#                s.one.two.comparison[[seg.name]][[samp.two]]$ph_cpnB_vec <- as.numeric(phased.s.two.s.cpnA.cpnB.df[, 2])
#                s.one.two.comparison[[seg.name]][[samp.two]]$mean_ph_cpnA <- mean(as.numeric(phased.s.two.s.cpnA.cpnB.df[, 1]))
#                s.one.two.comparison[[seg.name]][[samp.two]]$mean_ph_cpnB <- mean(as.numeric(phased.s.two.s.cpnA.cpnB.df[, 2]))
#                s.one.two.comparison[[seg.name]][[samp.two]]$un_ph_cpn_total <- NA
#                s.one.two.comparison[[seg.name]][[samp.two]]$un_ph_cpn_total_mean <- NA

#                s.one.two.comparison[[seg.name]][[samp.one.two.string]]$un_ph_cpnA_cpnB_het_comb_paired_pval    <- NA
#                s.one.two.comparison[[seg.name]][[samp.one.two.string]]$un_ph_cpnA_cpnB_het_comb_unpaired_pval  <- NA

#                s.one.two.comparison[[seg.name]][[samp.one.two.string]]$ph_cpnA_het_paired_pval <- seg.paired.p.vals.vec[1]
#                s.one.two.comparison[[seg.name]][[samp.one.two.string]]$ph_cpnB_het_paired_pval <- seg.paired.p.vals.vec[2]
#                s.one.two.comparison[[seg.name]][[samp.one.two.string]]$ph_cpnA_het_unpaired_pval <- seg.unpaired.p.vals.vec[1]
#                s.one.two.comparison[[seg.name]][[samp.one.two.string]]$ph_cpnB_het_unpaired_pval <- seg.unpaired.p.vals.vec[2]
#                s.one.two.comparison[[seg.name]][[samp.one.two.string]]$ph_cpnA_cpnB_het_comb_paired_pval <- seg.comb.paired.p.val
#                s.one.two.comparison[[seg.name]][[samp.one.two.string]]$ph_cpnA_cpnB_het_comb_unpaired_pval <- seg.comb.unpaired.p.val
#                s.one.two.comparison[[seg.name]][[samp.one.two.string]]$ph_tot_cpn_paired_pval <- total.cpn.paired.pval
#                s.one.two.comparison[[seg.name]][[samp.one.two.string]]$ph_cpnA_diff_by_size_het <- seg.cpnA.het.size.by.diff
#                s.one.two.comparison[[seg.name]][[samp.one.two.string]]$ph_cpnB_diff_by_size_het <- seg.cpnB.het.size.by.diff
#                s.one.two.comparison[[seg.name]][[samp.one.two.string]]$comb_cpnA_cpnB_diff_by_size_het <- seg.cpnA.het.size.by.diff + seg.cpnB.het.size.by.diff


#                s.one.two.comparison[[seg.name]][[samp.one.two.string]]$MSAI <- testPairOfSamplesSegForMSAIWithCPN(seg.name, samp.one, samp.two, all.samp.seg.phased.raw.cpn.list) 

#                samp.one.cpn.event.vec <- getSampLogRSegCPNEventVec(samp.one, s.one.s.ao.lr.baf.df, purity.est.vec, ploidy.est.vec,
#                                                                   gain.threshold, amp.threshold, loss.threshold, deep.loss.threshold)
#                s.one.two.comparison[[seg.name]][[samp.one]]$seg_cpn_event_status_vec <- samp.one.cpn.event.vec

#                samp.two.cpn.event.vec <- getSampLogRSegCPNEventVec(samp.two, s.two.s.ao.lr.baf.df, purity.est.vec, ploidy.est.vec,
#                                                                   gain.threshold, amp.threshold, loss.threshold, deep.loss.threshold)
#                s.one.two.comparison[[seg.name]][[samp.two]]$seg_cpn_event_status_vec <- samp.two.cpn.event.vec


#           }#end if if (nrow(logRBAFSample1) > min.num.hz.snp.phasing)

# 	}#end for ( i in 1:length(chr.start.end.string.vec) )

#      average.pairwise.het.over.all.segs <- sum(pw.CN.het.score.vec) / sum(seg.size.vec)

#      s.one.two.comparison[[samp.one.two.string]]$pairwise_avg_cn_het_score <- average.pairwise.het.over.all.segs

#      return(s.one.two.comparison)

# }#End calculatePairwiseCNHeterogeneity <- function(samp






# #Todo harmonise > or >= for minimum SNP number throughout code.
# #Todo add in 0.001 threshold as a parameter that can be varied
# calculatePairwiseCNHeterogeneity <- function(samp.one,
#                                              samp.two,
#                                              asas.both.lr.baf.lists,
#                                              asas.lr.and.or.baf.lists,
#                                              seg.name.vec,
#                                              ploidy.est.vec,
#                                              purity.est.vec,
#                                              all.seg.phasing.info.list,
#                                              all.samp.seg.phased.raw.cpn.list,
#                                              gamma.val,
#                                              gain.threshold = log2(2.5 / 2),
#                                              loss.threshold = log2(1.5 / 2),
#                                              amp.threshold = log2(4 / 2),
#                                              deep.loss.threshold = log2(0.5 / 2),
#                                              min.num.hz.snp.phasing = 5
#                                              ) {

#   print("HERE")
#   #debug parameters 
#   # samp.one <- "B_MET048_MT-M1"
#   # samp.two <- "B_MET048_T1-R1"

#   s.one.two.comparison.input.list <- list()
#   #Todo, why are we using seg.name or the index of the seg, choose one throughout code and stick to it 
#   for ( seg.name in seg.name.vec ) {
#     #This is the sample comparison ppart 

#     seg.name <- seg.name.vec[35]
#         print(seg.name)
 
#     s.one.two.comparison.input.list[[seg.name]] <- createAllSegInputForHetComparisons(seg.name,
#                                                   samp.one,
#                                                   samp.two,
#                                                   all.samp.seg.phased.raw.cpn.list,
#                                                   asas.both.lr.baf.lists,
#                                                   all.samp.non.hz.logr.and.tot.cpn.df.list,
#                                                   min.num.hz.snp.phasing)
#     print("Returned:")
#     print(str(s.one.two.comparison.input.list[[seg.name]], max=1))
#     print(str(s.one.two.comparison.input.list[[seg.name]][[1]], max=1))
#     STOP

#   }#end for ( seg.name in seg.name.vec ) {



#   return(s.one.two.comparison.input.list)
# }# end calculatePairwiseCNHeterogeneity


#one issue we need to do do comparisons with and iwthout the extra logr and then combined like w've done in the final dataframe,
#ideally these should be done so that the same funciton that does the comparison for the final joint segmentation should work for
#
 #we can do allele specific tests and total copy bnumber tests,,,]


#let's pull out the total copy number stuff from both sections and compare


#We want the pawixs classification of copy number wteen all regions but do we also want the actual measure based of them to be computed at the same time?

#Not really as we might have many different measures we woish to apply to the same data... let's try and keep it separate...





                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              