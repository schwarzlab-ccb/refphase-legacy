
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



checkSegGainOrLossSignificance <- function(logr.vec, samp.to.check, ploidy.est.vec, purity.est.vec, threshold, cpn.change.type = "gain", min.num.for.logr.test = 3) {

  min.thresh.raw.exp.samp       <-  ploidy.est.vec[samp.to.check] * (2 ^ (threshold))
  min.thresh.logr.exp.for.samp   <-  calcExpectedLogR(min.thresh.raw.exp.samp, 0, ploidy.est.vec[samp.to.check], purity.est.vec[samp.to.check])
  median.seg.logr                  <- median(logr.vec)
  
  if (length(logr.vec) >= min.num.for.logr.test) {

     gen.t.test.results        <- tTestWithErrorHandling(logr.vec, mu = min.thresh.logr.exp.for.samp)

    if (cpn.change.type=="gain") {

          specific.t.test.alt <- "greater"

    } else if (cpn.change.type=="loss") {

          specific.t.test.alt <- "less"

    }

    specific.t.test.result <- tTestWithErrorHandling(logr.vec, mu = min.thresh.logr.exp.for.samp, alternative = specific.t.test.alt)

    low.obs.logr.conf.int  <- specific.t.test.result$conf.int[1]
    high.obs.logr.conf.int <- specific.t.test.result$conf.int[2]
    gen.t.test.pval     <- gen.t.test.results$p.value
    alt.t.test.pval     <- specific.t.test.result$p.value
    
  }#end if (nrow(ss.baf.lr.baf.df) >= 3) {

  #Todo, refactor this magic number...
  if (length(logr.vec) < 3) {

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


##Todo, check fif the final segmentation functions can also use this instead of calling the check GainSegloss repeatedly
##Params:
##
##Returns:
##
getSampLogRSegCPNEventVec <- function(samp, samp.logr.vec, purity.est.vec, ploidy.est.vec, gain.threshold, amp.threshold, loss.threshold, deep.loss.threshold) {

  cpn.event.names.vec <- c('median_seg_logr','low_ci_logr_vs_exp_logr','high_ci_logr_vs_exp_logr'
                                     ,'min_exp_logr_gain','gain_pval'
                                     ,'min_exp_logr_amp','amp_pval'
                                     ,'min_exp_logr_loss','loss_pval'
                                     ,'min_exp_logr_deep_loss','deep_loss_pval')

  gain.samp.df <- checkSegGainOrLossSignificance(samp.logr.vec, samp, ploidy.est.vec, purity.est.vec, gain.threshold, "gain")
  amp.samp.df <- checkSegGainOrLossSignificance(samp.logr.vec, samp, ploidy.est.vec, purity.est.vec, amp.threshold, "gain")
  loss.samp.df <- checkSegGainOrLossSignificance(samp.logr.vec, samp, ploidy.est.vec, purity.est.vec, loss.threshold, "loss")
  deep.loss.samp.df <- checkSegGainOrLossSignificance(samp.logr.vec, samp, ploidy.est.vec, purity.est.vec, deep.loss.threshold, "loss")
  samp.cpn.event.vec <- c(gain.samp.df[, 3:5, drop = FALSE]
                           ,gain.samp.df[, 2, drop = FALSE], gain.samp.df[, 7, drop = FALSE]
                           ,amp.samp.df[, 2,drop = FALSE], amp.samp.df[, 7, drop = FALSE]
                           ,loss.samp.df[, 2, drop = FALSE], loss.samp.df[, 7, drop = FALSE]
                           ,deep.loss.samp.df[, 2, drop = FALSE], deep.loss.samp.df[, 7, drop = FALSE])
  names(samp.cpn.event.vec) <- cpn.event.names.vec

  return(samp.cpn.event.vec)

}#end getSampLogRSegCPNEventVec (..)


##Todo, make
##Params:
##
##Returns:
##
testPairOfSamplesSegForMSAIWithCPN <- function(seg.name, samp.one, samp.two, asas.ph.raw.cpn.lists) {

  phased.s.one.s.cpnA.cpnB.list <- asas.ph.raw.cpn.lists[[samp.one]][[seg.name]]
  phased.s.one.s.cpnA.cpnB.df <- as.matrix(cbind(phased.s.one.s.cpnA.cpnB.list$ph_cpnA_vec, phased.s.one.s.cpnA.cpnB.list$ph_cpnB_vec))
  phased.s.two.s.cpnA.cpnB.list <- asas.ph.raw.cpn.lists[[samp.one]][[seg.name]]
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
                                              asas.ph.raw.cpn.lists,
                                              asas.both.lr.baf.lists,
                                              all.samp.non.hz.logr.and.tot.cpn.df.list,
                                              min.num.hz.snp.phasing) {

  seg.samp.one.two.comp.data.list <- list()
  seg.samp.one.two.comp.data.list[[samp.one]] <- list()
  seg.samp.one.two.comp.data.list[[samp.two]] <- list()
  data.present.vec <- c(FALSE, FALSE, FALSE)
  names(data.present.vec) <- c("hz_phased", "non_hz", "hz_unphased")
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

  print(seg.name.split.vec)
  print("seg.name")
  print(seg.name)
  print(samp.one)
  print(" asas.ph.raw.cpn.lists[[seg.name]][[samp.one]]")
  print(str(asas.ph.raw.cpn.lists[[samp.one]][[seg.name]], max = 1))
  phased.s.one.s.cpnA.cpnB.list <- asas.ph.raw.cpn.lists[[samp.one]][[seg.name]]

  num.non.na.ph <- length(phased.s.one.s.cpnA.cpnB.list$ph_cpnA_vec[ ! is.na(phased.s.one.s.cpnA.cpnB.list$ph_cpnA_vec)])

  if ( num.non.na.ph >= min.num.hz.snp.phasing) {

    data.present.vec["hz_phased"] <- TRUE

    phased.s.one.s.cpnA.cpnB.df <- as.matrix(cbind(phased.s.one.s.cpnA.cpnB.list$ph_cpnA_vec, phased.s.one.s.cpnA.cpnB.list$ph_cpnB_vec))
    phased.s.two.s.cpnA.cpnB.list <- asas.ph.raw.cpn.lists[[samp.two]][[seg.name]]
    phased.s.two.s.cpnA.cpnB.df <- as.matrix(cbind(phased.s.two.s.cpnA.cpnB.list$ph_cpnA_vec, phased.s.two.s.cpnA.cpnB.list$ph_cpnB_vec))

    seg.samp.one.two.comp.data.list[[samp.one]]$ph_logr_vec <- phased.s.one.s.cpnA.cpnB.list$logr_vec
    seg.samp.one.two.comp.data.list[[samp.one]]$ph_cpnA_vec <- as.numeric(phased.s.one.s.cpnA.cpnB.df[, 1])
    seg.samp.one.two.comp.data.list[[samp.one]]$ph_cpnB_vec <- as.numeric(phased.s.one.s.cpnA.cpnB.df[, 2])
    seg.samp.one.two.comp.data.list[[samp.one]]$mean_ph_cpnA <- mean(as.numeric(phased.s.one.s.cpnA.cpnB.df[, 1]))
    seg.samp.one.two.comp.data.list[[samp.one]]$mean_ph_cpnB <- mean(as.numeric(phased.s.one.s.cpnA.cpnB.df[, 2]))

    seg.samp.one.two.comp.data.list[[samp.two]]$ph_logr_vec <- phased.s.two.s.cpnA.cpnB.list$logr_vec
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
getSubsetSegInfoForHetComparison <- function(analysis.type, samp.one, samp.two, seg.samp.one.two.comp.data.list) {

  data.present.vec <- seg.samp.one.two.comp.data.list$data_present_vec
  data.for.analysis.list <- list()

  if(all(data.present.vec == FALSE)) {

    return(NULL)
    
  }

  if ((analysis.type == "hz_phased" & data.present.vec["hz_phased"] == TRUE) | (analysis.type == "combined" & data.present.vec["non_hz"] == FALSE & data.present.vec["hz_phased"] == TRUE)) {
    #Analysis
    print("hzphased and is true")
    samp.one.logr.vec <- seg.samp.one.two.comp.data.list[[samp.one]]$ph_logr_vec
    samp.two.logr.vec <- seg.samp.one.two.comp.data.list[[samp.two]]$ph_logr_vec
    samp.one.tot.raw.cpn.vec <- seg.samp.one.two.comp.data.list[[samp.one]]$non_hz_tot_raw_cpn_vec
    samp.two.tot.raw.cpn.vec <- seg.samp.one.two.comp.data.list[[samp.two]]$non_hz_tot_raw_cpn_vec
    samp.one.ph.cpnA.vec <- seg.samp.one.two.comp.data.list[[samp.one]]$ph_cpnA_vec
    samp.one.ph.cpnB.vec <- seg.samp.one.two.comp.data.list[[samp.one]]$ph_cpnB_vec
    samp.two.ph.cpnA.vec <- seg.samp.one.two.comp.data.list[[samp.two]]$ph_cpnA_vec
    samp.two.ph.cpnB.vec <- seg.samp.one.two.comp.data.list[[samp.two]]$ph_cpnB_vec


  } else if (analysis.type == "hz_phased" & data.present.vec["hz_phased"] == FALSE) {

    return(NULL)

  }

  if ((analysis.type == "non_hz" & data.present.vec["non_hz"] == TRUE) | (analysis.type == "combined" & data.present.vec["non_hz"] == TRUE & data.present.vec["hz_phased"] == FALSE)){

    samp.one.logr.vec <- seg.samp.one.two.comp.data.list[[samp.one]]$non_hz_logr_vec
    samp.two.logr.vec <- seg.samp.one.two.comp.data.list[[samp.two]]$non_hz_logr_vec
    samp.one.tot.raw.cpn.vec <- seg.samp.one.two.comp.data.list[[samp.one]]$non_hz_tot_raw_cpn_vec
    samp.two.tot.raw.cpn.vec <- seg.samp.one.two.comp.data.list[[samp.two]]$non_hz_tot_raw_cpn_vec
    samp.one.ph.cpnA.vec <- NA
    samp.one.ph.cpnB.vec <- NA
    samp.two.ph.cpnA.vec <- NA
    samp.two.ph.cpnB.vec <- NA

  } else if (analysis.type == "non_hz" & data.present.vec["non_hz"] == FALSE) {

    return(NULL)

  }

  if (analysis.type == "combined" & data.present.vec["non_hz"] == TRUE & data.present.vec["hz_phased"] == TRUE) {
    #Make sure we add in the unphased heterozygous snps here if they are present

    samp.one.logr.vec <- c(seg.samp.one.two.comp.data.list[[samp.one]]$ph_logr_vec, seg.samp.one.two.comp.data.list[[samp.one]]$non_hz_logr_vec)
    samp.two.logr.vec <- c(seg.samp.one.two.comp.data.list[[samp.two]]$ph_logr_vec, seg.samp.one.two.comp.data.list[[samp.two]]$non_hz_logr_vec)
    samp.one.tot.raw.cpn.vec <- c(seg.samp.one.two.comp.data.list[[samp.one]]$non_hz_tot_raw_cpn_vec, seg.samp.one.two.comp.data.list[[samp.one]]$non_hz_tot_raw_cpn_vec) 
    samp.two.tot.raw.cpn.vec <- c(seg.samp.one.two.comp.data.list[[samp.two]]$non_hz_tot_raw_cpn_vec, seg.samp.one.two.comp.data.list[[samp.one]]$non_hz_tot_raw_cpn_vec)
    samp.one.ph.cpnA.vec <- seg.samp.one.two.comp.data.list[[samp.one]]$ph_cpnA_vec
    samp.one.ph.cpnB.vec <- seg.samp.one.two.comp.data.list[[samp.one]]$ph_cpnB_vec
    samp.two.ph.cpnA.vec <- seg.samp.one.two.comp.data.list[[samp.two]]$ph_cpnA_vec
    samp.two.ph.cpnB.vec <- seg.samp.one.two.comp.data.list[[samp.two]]$ph_cpnB_vec

  }

  if (analysis.type == "combined" & data.present.vec["hz_unphased"] == TRUE) {

    samp.one.logr.vec <- c(samp.one.logr.vec, seg.samp.one.two.comp.data.list[[samp.one]]$un_ph_hz_logr_vec)
    samp.two.logr.vec <- c(samp.two.logr.vec, seg.samp.one.two.comp.data.list[[samp.two]]$un_ph_hz_logr_vec)
    samp.one.tot.raw.cpn.vec <- c(samp.one.tot.raw.cpn.vec, seg.samp.one.two.comp.data.list[[samp.one]]$un_ph_hz_tot_raw_cpn_vec)
    samp.two.tot.raw.cpn.vec <- c(samp.two.tot.raw.cpn.vec, seg.samp.one.two.comp.data.list[[samp.two]]$un_ph_hz_tot_raw_cpn_vec)

  }#end if (data.present.vec["hz_unphased"] == TRUE) {

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
performSegPairwiseComparisons <- function(seg.name,
                                          samp.one,
                                          samp.two,
                                          seg.samp.one.two.comp.data.list,
                                          purity.est.vec,
                                          ploidy.est.vec,
                                          gain.threshold,
                                          amp.threshold,
                                          loss.threshold,
                                          deep.loss.threshold) {

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
  print("data.present.vec")
  print(data.present.vec)

  seg.het.data.list <- list() 
  seg.het.data.list$data_present_vec <- data.present.vec


  for ( analysis.type in analysis.type.vec ) {

    s.data.list <- getSubsetSegInfoForHetComparison(analysis.type, samp.one, samp.two, seg.samp.one.two.comp.data.list)

    if (is.null(s.data.list)) {
      print("Analysis not performed:")
      print(analysis.type)
      print(seg.name)
      next
    }

    seg.name.split.vec <- strsplit(seg.name, split=":")
    seg.het.data.list$seg_size <- as.numeric(unlist(seg.name.split.vec)[3]) - as.numeric(unlist(seg.name.split.vec)[2])

    if (analysis.type == "hz_phased" & data.present.vec[analysis.type] == TRUE) {

      seg.het.data.list$cpnA_paired_samp_one_two_pval <- tTestWithErrorHandling(s.data.list$samp_one_ph_cpnA_vec, s.data.list$samp_two_ph_cpnA_vec, paired = TRUE)$p.val                  
      seg.het.data.list$cpnB_paired_samp_one_two_pval <- tTestWithErrorHandling(s.data.list$samp_one_ph_cpnB_vec, s.data.list$samp_two_ph_cpnB_vec, paired = TRUE)$p.val
      seg.het.data.list$seg_comb_cpnAB_paired_pval <- pchisq((sum(log(c(seg.het.data.list$cpnA_paired_samp_one_two_pval, seg.het.data.list$cpnB_paired_samp_one_two_pval))) * - 2),
                                                              df = 2 * 2,
                                                              lower.tail = FALSE)

      seg.het.data.list$cpnA_unpaired_samp_one_two_pval <- tTestWithErrorHandling(s.data.list$samp_one_ph_cpnA_vec, s.data.list$samp_two_ph_cpnA_vec)$p.val
      seg.het.data.list$cpnB_unpaired_samp_one_two_pval <- tTestWithErrorHandling(s.data.list$samp_one_ph_cpnB_vec, s.data.list$samp_two_ph_cpnB_vec)$p.val
      seg.het.data.list$seg_comb_cpnAB_unpaired_pval <- pchisq((sum(log(c(seg.het.data.list$cpnA_unpaired_samp_one_two_pval, seg.het.data.list$cpnB_unpaired_samp_one_two_pval))) * - 2),
                                                                df = 2 * 2,
                                                                lower.tail = FALSE)

      seg.het.data.list$cpnA_diff_amplitude <- abs(mean(as.numeric(s.data.list$samp_one_ph_cpnA_vec)) - mean(as.numeric(s.data.list$samp_two_ph_cpnA_vec)))
      seg.het.data.list$cpnB_diff_amplitude <- abs(mean(as.numeric(s.data.list$samp_one_ph_cpnB_vec)) - mean(as.numeric(s.data.list$samp_two_ph_cpnB_vec)))

      seg.het.data.list$raw_cpn_MSAI <- testPairOfSamplesSegForMSAIWithCPN(seg.name, samp.one, samp.two, asas.ph.raw.cpn.lists) 

    }

    seg.het.data.list[[paste0(analysis.type, "_total_cpn_paired_pval")]] <- tTestWithErrorHandling(s.data.list$samp_one_tot_raw_cpn_vec, s.data.list$samp_two_tot_raw_cpn_vec, paired = TRUE)$p.val

    seg.het.data.list[[paste0(analysis.type, "_samp_one_cpn_event_vec")]] <- getSampLogRSegCPNEventVec(samp.one, s.data.list$samp_one_logr_vec, purity.est.vec, ploidy.est.vec,
                                                       gain.threshold, amp.threshold, loss.threshold, deep.loss.threshold)

    seg.het.data.list[[paste0(analysis.type, "_samp_two_cpn_event_vec")]] <- getSampLogRSegCPNEventVec(samp.two, s.data.list$samp_two_logr_vec, purity.est.vec, ploidy.est.vec,
                                                     gain.threshold, amp.threshold, loss.threshold, deep.loss.threshold)

    seg.het.data.list[[paste0(analysis.type, "_logr_samp_one_vs_two_pval")]] <- tTestWithErrorHandling(s.data.list$samp_one_logr_vec,  s.data.list$samp_two_logr_vec)$p.val;

  }

  return(seg.het.data.list)

}




calculateSegPairwiseCNHeterogeneityMeasures <- function(seg.samp.comp.het.analysis.list, pval.for.diff.amplitude.thresh = 0.001) {

  data.present.vec <- seg.samp.comp.het.analysis.list$data_present_vec

  seg.het.measures.list <- list()
  seg.het.measures.list$data_present_vec <- data.present.vec

  if (all(seg.samp.comp.het.analysis.list$data_present_vec == FALSE)) {

    return(seg.het.measures.list)

  }


  seg.het.measures.list$cpnA_het_size_by_diff <- 0
  seg.het.measures.list$cpnB_het_size_by_diff <- 0 

  seg.het.measures.list$seg_size <- seg.samp.comp.het.analysis.list$seg_size

  if( data.present.vec["hz_phased"] == TRUE ) {

    #Todo make the 0.001 a parameter..
    if (seg.samp.comp.het.analysis.list$cpnA_paired_samp_one_two_pval < pval.for.diff.amplitude.thresh) {

         seg.het.measures.list$cpnA_het_size_by_diff <- seg.samp.comp.het.analysis.list$seg_size * seg.samp.comp.het.analysis.list$cpnA_diff_amplitude

    }#end if (seg.paired.p.vals.vec[1] < 0.001) {

    if (seg.samp.comp.het.analysis.list$cpnB_paired_samp_one_two_pval < pval.for.diff.amplitude.thresh) {

         seg.het.measures.list$cpnB_het_size_by_diff <- seg.samp.comp.het.analysis.list$seg_size * seg.samp.comp.het.analysis.list$cpnB_diff_amplitude

    }#end if (seg.paired.p.vals.vec[2] < 0.001) { 


  }# 

  seg.het.measures.list$cpnA_cpnB_tot_het_size_by_diff <- seg.het.measures.list$cpnA_het_size_by_diff + seg.het.measures.list$cpnB_het_size_by_diff

  return(seg.het.measures.list)
  
}




                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                