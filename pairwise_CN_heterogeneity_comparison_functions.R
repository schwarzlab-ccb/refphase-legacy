

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



checkSegGainOrLossSignificance <- function(ss.baf.lr.baf.df, samp.to.check, ploidy.est.vec, purity.est.vec, threshold, cpn.change.type="gain") {

  min.thresh.raw.exp.samp       <-  ploidy.est.vec[samp.to.check] * (2 ^ (threshold))
  min.thresh.logr.exp.for.samp      <-  calcExpectedLogR(min.thresh.raw.exp.samp, 0, ploidy.est.vec[samp.to.check], purity.est.vec[samp.to.check])
  median.seg.logr                  <- median(ss.baf.lr.baf.df$logR)
  
  if (nrow(ss.baf.lr.baf.df) >= 3) {

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


#Todo harmonise > or >= for minimum SNP number throughout code.
#Todo add in 0.001 threshold as a parameter that can be varied
calculatePairwiseCNHeterogeneity <- function(samp.one,
                                             samp.two,
                                             asas.both.lr.baf.lists,
                                             asas.lr.and.or.baf.lists,
                                             min.cons.segs.df,
                                             ploidy.est.vec,
                                             purity.est.vec,
                                             all.seg.phasing.info.list,
                                             gamma.val,
                                             gain.threshold = log2(2.5 / 2),
                                             loss.threshold = log2(1.5 / 2),
                                             amp.threshold = log2(4 / 2),
                                             deep.loss.threshold = log2(0.5 / 2),
                                             min.num.hz.snp.phasing=5
                                             ){

     cpn.event.names.vec <- c('median_seg_logr','low_ci_logr_vs_exp_logr','high_ci_logr_vs_exp_logr'
                                        ,'min_exp_logr_gain','gain_pval'
                                        ,'min_exp_logr_amp','amp_pval'
                                        ,'min_exp_logr_loss','loss_pval'
                                        ,'min_exp_logr_deep_loss','deep_loss_pval')


	#assuming that there will always be a phasing vector
	s.one.two.comparison <- list()
	chr.start.end.string.vec <- paste(min.cons.segs.df[ , 1], min.cons.segs.df[ , 2], min.cons.segs.df[ , 3], sep=":")
	pw.CN.het.score.vec <- c()
	seg.size.vec <- c()




     #Todo, why are we using seg.name or the index of the seg, choose one throughout code and stick to it 
	for ( i in 1:length(chr.start.end.string.vec) ) {

		seg.name <- chr.start.end.string.vec[i]
    print(seg.name)
          samp.one.two.string <- paste(samp.one, samp.two, sep=":")
          s.one.s.b.lr.baf.df <- asas.both.lr.baf.lists[[samp.one]][[seg.name]]
          s.two.s.b.lr.baf.df <- asas.both.lr.baf.lists[[samp.two]][[seg.name]]

          s.one.s.ao.lr.baf.df <- asas.lr.and.or.baf.lists[[samp.one]][[seg.name]]
          s.two.s.ao.lr.baf.df <- asas.lr.and.or.baf.lists[[samp.two]][[seg.name]]

          s.one.two.comparison[[seg.name]] <- initialiseSampOneTwoSegPairwiseComparisonList(samp.one, samp.two, samp.one.two.string)
          print(nrow(s.one.s.b.lr.baf.df))
          if (nrow(s.one.s.b.lr.baf.df) < min.num.hz.snp.phasing) {

               s.one.s.cpnA.cpnB.df <- calculateSNPsCopyNumberAB(logr.vec = s.two.s.ao.lr.baf.df$logR,
                                                                 baf.vec = s.two.s.ao.lr.baf.df$BAF,
                                                                 tumor.ploidy = ploidy.est.vec[samp.one],
                                                                 tumor.purity = purity.est.vec[samp.one],
                                                                 gamma.val = gamma.val)
               s.two.s.cpnA.cpnB.df <- calculateSNPsCopyNumberAB(logr.vec = s.two.s.ao.lr.baf.df$logR,
                                                                 baf.vec = s.two.s.ao.lr.baf.df$BAF,
                                                                 tumor.ploidy = ploidy.est.vec[samp.two],
                                                                 tumor.purity = purity.est.vec[samp.two],
                                                                 gamma.val = gamma.val)

               seg.paired.p.vals.vec<- c(tTestWithErrorHandling(s.one.s.cpnA.cpnB.df[,1],s.two.s.cpnA.cpnB.df[,1],paired=TRUE)$p.val,
                                         tTestWithErrorHandling(s.one.s.cpnA.cpnB.df[,2],s.two.s.cpnA.cpnB.df[,2],paired=TRUE)$p.val)
               combined.cpnA.cpnB.paired.pval   <- seg.paired.p.vals.vec[1]
               seg.unpaired.p.vals.vec  <- c(tTestWithErrorHandling(s.one.s.cpnA.cpnB.df[,1],s.two.s.cpnA.cpnB.df[,1])$p.val,
                                             tTestWithErrorHandling(s.one.s.cpnA.cpnB.df[,2],s.two.s.cpnA.cpnB.df[,2])$p.val)
               combined.cpnA.cpnB.unpaired.pval <- seg.unpaired.p.vals.vec[1]#pchisq((sum(log(PValsUnPaired))*-2), df=length(PValsUnPaired)*2, lower.tail=F)
               
               # also get a p-value that relates to both samples combined
               total.cpn.paired.pval <- tTestWithErrorHandling(c(phased.s.one.s.cpnA.cpnB.df[, 1] + phased.s.one.s.cpnA.cpnB.df[, 2]),
                                                               c(phased.s.two.s.cpnA.cpnB.df[, 1] + phased.s.two.s.cpnA.cpnB.df[, 2]), paired = TRUE)$p.val


               seg.name.split.vec <- strsplit(chr.start.end.string.vec[i], split=":")
               seg.size <- as.numeric(unlist(seg.name.split.vec)[3]) - as.numeric(unlist(seg.name.split.vec)[2])

               s.one.two.comparison[[seg.name]][[samp.one]]$chr_vec <- s.one.s.b.lr.baf.df$chrom
               s.one.two.comparison[[seg.name]][[samp.one]]$pos_vec <- s.one.s.b.lr.baf.df$pos
               s.one.two.comparison[[seg.name]][[samp.two]]$chr_vec <- s.one.s.b.lr.baf.df$chrom
               s.one.two.comparison[[seg.name]][[samp.two]]$pos_vec <- s.one.s.b.lr.baf.df$pos

               s.one.two.comparison[[seg.name]][[samp.one]]$ph_cpnA_vec <- NA
               s.one.two.comparison[[seg.name]][[samp.one]]$ph_cpnB_vec <- NA
               s.one.two.comparison[[seg.name]][[samp.one]]$mean_ph_cpnA <- NA
               s.one.two.comparison[[seg.name]][[samp.one]]$mean_ph_cpnB <- NA

               s.one.two.comparison[[seg.name]][[samp.one]]$un_ph_cpn_total <- s.one.s.cpnA.cpnB.df[, 1] + s.one.s.cpnA.cpnB.df[, 2]
               s.one.two.comparison[[seg.name]][[samp.one]]$un_ph_cpn_total_mean <- mean(s.one.s.cpnA.cpnB.df[, 1] + s.one.s.cpnA.cpnB.df[, 2])

               s.one.two.comparison[[seg.name]][[samp.two]]$ph_cpnA_vec <- NA
               s.one.two.comparison[[seg.name]][[samp.two]]$ph_cpnB_vec <- NA
               s.one.two.comparison[[seg.name]][[samp.two]]$mean_ph_cpnA <- NA
               s.one.two.comparison[[seg.name]][[samp.two]]$mean_ph_cpnB <- NA

               s.one.two.comparison[[seg.name]][[samp.two]]$un_ph_cpn_total <- s.two.s.cpnA.cpnB.df[, 1] + s.two.s.cpnA.cpnB.df[, 2]
               s.one.two.comparison[[seg.name]][[samp.two]]$un_ph_cpn_total_mean <- mean(s.two.s.cpnA.cpnB.df[, 1] + s.two.s.cpnA.cpnB.df[, 2])

               s.one.two.comparison[[seg.name]][[samp.one.two.string]]$ph_cpnA_het_paired_pval <- NA
               s.one.two.comparison[[seg.name]][[samp.one.two.string]]$ph_cpnB_het_paired_pval <- NA
               s.one.two.comparison[[seg.name]][[samp.one.two.string]]$ph_cpnA_het_unpaired_pval <- NA
               s.one.two.comparison[[seg.name]][[samp.one.two.string]]$ph_cpnB_het_unpaired_pval <- NA

               s.one.two.comparison[[seg.name]][[samp.one.two.string]]$un_ph_cpnA_cpnB_het_comb_paired_pval    <- combined.cpnA.cpnB.paired.pval
               s.one.two.comparison[[seg.name]][[samp.one.two.string]]$un_ph_cpnA_cpnB_het_comb_unpaired_pval  <- combined.cpnA.cpnB.unpaired.pval
               s.one.two.comparison[[seg.name]][[samp.one.two.string]]$un_ph_combinedUnPairedPval      <- total.cpn.paired.pval

               s.one.two.comparison[[seg.name]][[samp.one.two.string]]$ph_cpnA_cpnB_het_comb_paired_pval <- NA
               s.one.two.comparison[[seg.name]][[samp.one.two.string]]$ph_cpnA_cpnB_het_comb_unpaired_pval <- NA
               s.one.two.comparison[[seg.name]][[samp.one.two.string]]$ph_tot_cpn_paired_pval <- NA

               s.one.two.comparison[[seg.name]][[samp.one.two.string]]$ph_cpnA_diff_by_size_het <- NA
               s.one.two.comparison[[seg.name]][[samp.one.two.string]]$ph_cpnB_diff_by_size_het <- NA
               s.one.two.comparison[[seg.name]][[samp.one.two.string]]$comb_cpnA_cpnB_diff_by_size_het <- NA
               s.one.two.comparison[[seg.name]][[samp.one.two.string]]$MSAI <- NA


               gain.samp.one.df <- checkSegGainOrLossSignificance(s.one.s.ao.lr.baf.df, samp.one, ploidy.est.vec, purity.est.vec, gain.threshold, "gain")
               amp.samp.one.df <- checkSegGainOrLossSignificance(s.one.s.ao.lr.baf.df, samp.one, ploidy.est.vec, purity.est.vec, amp.threshold, "gain")
               loss.samp.one.df <- checkSegGainOrLossSignificance(s.one.s.ao.lr.baf.df, samp.one, ploidy.est.vec, purity.est.vec, loss.threshold, "loss")
               deep.loss.samp.one.df <- checkSegGainOrLossSignificance(s.one.s.ao.lr.baf.df, samp.one, ploidy.est.vec, purity.est.vec, deep.loss.threshold, "loss")
               samp.one.cpn.event.vec <- c(gain.samp.one.df[, 3:5, drop=FALSE]
                                        ,gain.samp.one.df[, 2, drop=FALSE], gain.samp.one.df[, 7, drop=FALSE]
                                        ,amp.samp.one.df[, 2,drop=FALSE], amp.samp.one.df[, 7, drop=FALSE]
                                        ,loss.samp.one.df[, 2, drop=FALSE], loss.samp.one.df[, 7, drop=FALSE]
                                        ,deep.loss.samp.one.df[, 2, drop=FALSE], deep.loss.samp.one.df[, 7, drop=FALSE])
               names(samp.one.cpn.event.vec) <- cpn.event.names.vec

               s.one.two.comparison[[seg.name]][[samp.one]]$seg_cpn_event_status_vec <- samp.one.cpn.event.vec

               gain.samp.two.df <- checkSegGainOrLossSignificance(s.one.s.ao.lr.baf.df, samp.two, ploidy.est.vec, purity.est.vec, gain.threshold, "gain")
               amp.samp.two.df <- checkSegGainOrLossSignificance(s.one.s.ao.lr.baf.df, samp.two, ploidy.est.vec, purity.est.vec, amp.threshold, "gain")
               loss.samp.two.df <- checkSegGainOrLossSignificance(s.one.s.ao.lr.baf.df, samp.two, ploidy.est.vec, purity.est.vec, loss.threshold, "loss")
               deep.loss.samp.two.df <- checkSegGainOrLossSignificance(s.one.s.ao.lr.baf.df, samp.two, ploidy.est.vec, purity.est.vec, deep.loss.threshold, "loss")
               samp.two.cpn.event.vec <- c(gain.samp.two.df[, 3:5, drop=FALSE]
                                        ,gain.samp.two.df[, 2, drop=FALSE], gain.samp.two.df[, 7, drop=FALSE]
                                        ,amp.samp.two.df[, 2,drop=FALSE], amp.samp.two.df[, 7, drop=FALSE]
                                        ,loss.samp.two.df[, 2, drop=FALSE], loss.samp.two.df[, 7, drop=FALSE]
                                        ,deep.loss.samp.two.df[, 2, drop=FALSE], deep.loss.samp.two.df[, 7, drop=FALSE])
               names(samp.two.cpn.event.vec) <- cpn.event.names.vec

               s.one.two.comparison[[seg.name]][[samp.two]]$seg_cpn_event_status_vec <- samp.two.cpn.event.vec

               next;

          }#end f (nrow(s.one.s.b.lr.baf.df) < min.num.hz.snp.phasing) {

          if (nrow(s.one.s.b.lr.baf.df) > min.num.hz.snp.phasing) {

               #Todo getting the major and A and B columns can be made into a function...
               s.one.s.cpnA.cpnB.df <- calculateSNPsCopyNumberAB(logr.vec = s.one.s.b.lr.baf.df$logR,
                                                                 baf.vec = s.one.s.b.lr.baf.df$BAF,
                                                                 tumor.ploidy = ploidy.est.vec[samp.one],
                                                                 tumor.purity = purity.est.vec[samp.one],
                                                                 gamma.val = gamma.val)
               s.two.s.cpnA.cpnB.df <- calculateSNPsCopyNumberAB(logr.vec = s.two.s.b.lr.baf.df$logR,
                                                                 baf.vec = s.two.s.b.lr.baf.df$BAF,
                                                                 tumor.ploidy = ploidy.est.vec[samp.two],
                                                                 tumor.purity = purity.est.vec[samp.two],
                                                                 gamma.val = gamma.val)
               s.one.s.max.cpn.vec <- apply(s.one.s.cpnA.cpnB.df, 1, max)
               s.two.s.max.cpn.vec <- apply(s.two.s.cpnA.cpnB.df, 1, max)
               s.one.s.median.max.cpn <- median(s.one.s.max.cpn.vec)
               s.two.s.median.max.cpn <- median(s.two.s.max.cpn.vec)
               s.one.s.is.max.cpn.cpnA.cpnB.logical.df <- ifelse(s.one.s.cpnA.cpnB.df==s.one.s.max.cpn.vec, 1, 0)
               s.two.s.is.max.cpn.cpnA.cpnB.logical.df <- ifelse(s.two.s.cpnA.cpnB.df==s.two.s.max.cpn.vec, 1, 0)

               if (identical(s.one.s.is.max.cpn.cpnA.cpnB.logical.df, s.two.s.is.max.cpn.cpnA.cpnB.logical.df)) {
                    s.one.s.classifier <- 0.05#Why these magic numbers...
                    s.two.s.classifier <- 1#Why these magic numbers...
               }#end if (identical(s.one.s.cpnA.cpnB.is.max.df, s.two.s.cpnA.cpnB.is.max.df)) {

               if ( ! identical(s.one.s.is.max.cpn.cpnA.cpnB.logical.df, s.two.s.is.max.cpn.cpnA.cpnB.logical.df)) {
                    #Testing Major vs minor cpn in samp two versus samp.one
                    s.one.s.classifier <- tTestWithErrorHandling(c(s.two.s.cpnA.cpnB.df[which(s.one.s.is.max.cpn.cpnA.cpnB.logical.df[,1]==1), 1], 
                                                                   s.two.s.cpnA.cpnB.df[which(s.one.s.is.max.cpn.cpnA.cpnB.logical.df[,2]==1), 1]),
                                                                 c(s.two.s.cpnA.cpnB.df[which(s.one.s.is.max.cpn.cpnA.cpnB.logical.df[,1]!=1), 2],
                                                                   s.two.s.cpnA.cpnB.df[which(s.one.s.is.max.cpn.cpnA.cpnB.logical.df[,2]!=1), 2]))$p.value
                    #Tesing Major vs minor copn in samp one versus samp two
                    s.two.s.classifier <- tTestWithErrorHandling(c(s.one.s.cpnA.cpnB.df[which(s.two.s.is.max.cpn.cpnA.cpnB.logical.df[,1]==1), 1], 
                                                                   s.one.s.cpnA.cpnB.df[which(s.two.s.is.max.cpn.cpnA.cpnB.logical.df[,2]==1), 1]),
                                                                 c(s.one.s.cpnA.cpnB.df[which(s.two.s.is.max.cpn.cpnA.cpnB.logical.df[,1]!=1), 2],
                                                                   s.one.s.cpnA.cpnB.df[which(s.two.s.is.max.cpn.cpnA.cpnB.logical.df[,2]!=1), 2]))$p.value
                    if (is.na(s.one.s.classifier)) {

                        s.one.s.classifier <- 1
                        
                    }
                    if (is.na(s.two.s.classifier)) {

                        s.two.s.classifier <- 1

                    }

               }#end if ( ! identical(s.one.s.cpnA.cpnB.is.max.df, s.two.s.cpnA.cpnB.is.max.df)) {

               if (s.one.s.classifier <= s.two.s.classifier) {

                    maj.is.cpnA.idx.vec <- which(s.one.s.is.max.cpn.cpnA.cpnB.logical.df[ , 1] == 1)
                    maj.is.cpnB.idx.vec <- which(s.one.s.is.max.cpn.cpnA.cpnB.logical.df[ , 2] == 1)

               }#end if (s.one.s.classifier <= s.two.s.classifier)  {

               if (s.two.s.classifier > s.two.s.classifier) {

                    maj.is.cpnA.idx.vec <- which(s.one.s.is.max.cpn.cpnA.cpnB.logical.df[ , 1] == 1)
                    maj.is.cpnB.idx.vec <- which(s.one.s.is.max.cpn.cpnA.cpnB.logical.df[ , 2] == 1)

               }#end if (s.two.s.classifier > s.two.s.classifier) {

               if (length(all.seg.phasing.info.list[[seg.name]]$maj_is_cpnA_idx_vec) > 0) {

                    maj.is.cpnA.idx.vec <- all.seg.phasing.info.list[[seg.name]]$maj_is_cpnA_idx_vec
                    
                } else { #end if (length(all.seg.phasing.info.list[[seg.name]]$maj_is_cpnA_idx_vec) > 0)
                
                    maj.is.cpnA.idx.vec <- c()
                }


                if (length(all.seg.phasing.info.list[[seg.name]]$maj_is_cpnB_idx_vec) > 0) {

                     maj.is.cpnA.idx.vec <- all.seg.phasing.info.list[[seg.name]]$maj_is_cpnB_idx_vec
                     
                 } else { #end if (length(all.seg.phasing.info.list[[seg.name]]$maj_is_cpnA_idx_vec) > 0)
                 
                     maj.is.cpnB.idx.vec <- c()
                 }




               maj.is.cpnB.idx.vec <- all.seg.phasing.info.list[[seg.name]]$maj_is_cpnB_idx_vec
               #Now we need to apply the phasing/classifier to both samples

              phased.s.one.s.cpnA.cpnB.df <- s.one.s.cpnA.cpnB.df
              phased.s.two.s.cpnA.cpnB.df <- s.two.s.cpnA.cpnB.df

               if (length(maj.is.cpnA.idx.vec) > 0) {

                  phased.s.one.s.cpnA.cpnB.df[maj.is.cpnA.idx.vec, 1] <- s.one.s.cpnA.cpnB.df[maj.is.cpnA.idx.vec, 1]
                  phased.s.one.s.cpnA.cpnB.df[maj.is.cpnA.idx.vec, 2] <- s.one.s.cpnA.cpnB.df[maj.is.cpnA.idx.vec, 2]

                  phased.s.two.s.cpnA.cpnB.df[maj.is.cpnA.idx.vec, 1] <- s.two.s.cpnA.cpnB.df[maj.is.cpnA.idx.vec, 1]
                  phased.s.two.s.cpnA.cpnB.df[maj.is.cpnA.idx.vec, 2] <- s.two.s.cpnA.cpnB.df[maj.is.cpnA.idx.vec, 2]

               }

               if (length(maj.is.cpnB.idx.vec) > 0) {

                  phased.s.one.s.cpnA.cpnB.df[maj.is.cpnB.idx.vec, 1] <- s.one.s.cpnA.cpnB.df[maj.is.cpnB.idx.vec, 2]
                  phased.s.one.s.cpnA.cpnB.df[maj.is.cpnB.idx.vec, 2] <- s.one.s.cpnA.cpnB.df[maj.is.cpnB.idx.vec, 1]

                  phased.s.two.s.cpnA.cpnB.df[maj.is.cpnB.idx.vec, 1] <- s.two.s.cpnA.cpnB.df[maj.is.cpnB.idx.vec, 2]
                  phased.s.two.s.cpnA.cpnB.df[maj.is.cpnB.idx.vec, 2] <- s.two.s.cpnA.cpnB.df[maj.is.cpnB.idx.vec, 1]

               }



               seg.paired.p.vals.vec <- c(tTestWithErrorHandling(phased.s.one.s.cpnA.cpnB.df[, 1],
                                                                 phased.s.two.s.cpnA.cpnB.df[, 1], paired = TRUE)$p.val,
                                          tTestWithErrorHandling(phased.s.one.s.cpnA.cpnB.df[, 2],
                                                                 phased.s.two.s.cpnA.cpnB.df[, 2], paired = TRUE)$p.val)

               seg.comb.paired.p.val <- pchisq((sum(log(seg.paired.p.vals.vec)) * - 2), df = length(seg.paired.p.vals.vec) * 2, lower.tail = FALSE)

               seg.unpaired.p.vals.vec <- c(tTestWithErrorHandling(phased.s.one.s.cpnA.cpnB.df[, 1],
                                                                   phased.s.two.s.cpnA.cpnB.df[, 1])$p.val,
                                            tTestWithErrorHandling(phased.s.one.s.cpnA.cpnB.df[, 2],
                                                                   phased.s.two.s.cpnA.cpnB.df[, 2])$p.val)

               seg.comb.unpaired.p.val <- pchisq((sum(log(seg.unpaired.p.vals.vec)) * - 2), df = length(seg.unpaired.p.vals.vec) * 2, lower.tail = FALSE)

               # also get a p-value that relates to both samples combined
               total.cpn.paired.pval <- tTestWithErrorHandling(c(phased.s.one.s.cpnA.cpnB.df[, 1] + phased.s.one.s.cpnA.cpnB.df[, 2]),
                                                               c(phased.s.two.s.cpnA.cpnB.df[, 1] + phased.s.two.s.cpnA.cpnB.df[, 2]), paired = TRUE)$p.val

               seg.cpnA.het.size.by.diff <- 0            
               seg.cpnB.het.size.by.diff <- 0

               seg.name.split.vec <- strsplit(chr.start.end.string.vec[i], split=":")
               seg.size <- as.numeric(unlist(seg.name.split.vec)[3]) - as.numeric(unlist(seg.name.split.vec)[2])

               #Todo make the 0.001 a parameter..
               if (seg.paired.p.vals.vec[1] < 0.001) {

                    cpnA.diff.amplitude <- abs(mean(as.numeric(phased.s.one.s.cpnA.cpnB.df[,1])) - mean(as.numeric(phased.s.two.s.cpnA.cpnB.df[,1])))
                    seg.cpnA.het.size.by.diff <- seg.size * cpnA.diff.amplitude

               }#end if (seg.paired.p.vals.vec[1] < 0.001) {

               if (seg.paired.p.vals.vec[2] < 0.001) {

                    cpnB.diff.amplitude <- abs(mean(as.numeric(phased.s.one.s.cpnA.cpnB.df[,2])) - mean(as.numeric(phased.s.two.s.cpnA.cpnB.df[,2])))
                    seg.cpnB.het.size.by.diff <- seg.size * cpnB.diff.amplitude

               }#end if (seg.paired.p.vals.vec[2] < 0.001) { 

               seg.tot.cpnA.cpnB.het <- seg.cpnA.het.size.by.diff + seg.cpnB.het.size.by.diff
               pw.CN.het.score.vec <- c(pw.CN.het.score.vec, seg.tot.cpnA.cpnB.het)
               seg.size.vec <- c(seg.size.vec, seg.size)

               s.one.two.comparison[[seg.name]][[samp.one]]$chr_vec <- s.one.s.b.lr.baf.df$chrom
               s.one.two.comparison[[seg.name]][[samp.one]]$pos_vec <- s.one.s.b.lr.baf.df$chrom

               s.one.two.comparison[[seg.name]][[samp.one]]$ph_cpnA_vec <- as.numeric(phased.s.one.s.cpnA.cpnB.df[,1])
               s.one.two.comparison[[seg.name]][[samp.one]]$ph_cpnB_vec <- as.numeric(phased.s.one.s.cpnA.cpnB.df[,2])
               s.one.two.comparison[[seg.name]][[samp.one]]$mean_ph_cpnA <- mean(as.numeric(phased.s.one.s.cpnA.cpnB.df[,1]))
               s.one.two.comparison[[seg.name]][[samp.one]]$mean_ph_cpnB <- mean(as.numeric(phased.s.one.s.cpnA.cpnB.df[,2]))
               s.one.two.comparison[[seg.name]][[samp.one]]$un_ph_cpn_total <- NA
               s.one.two.comparison[[seg.name]][[samp.one]]$un_ph_cpn_total_mean <- NA

               s.one.two.comparison[[seg.name]][[samp.two]]$ph_cpnA_vec <- as.numeric(phased.s.two.s.cpnA.cpnB.df[,1])
               s.one.two.comparison[[seg.name]][[samp.two]]$ph_cpnB_vec <- as.numeric(phased.s.two.s.cpnA.cpnB.df[,2])
               s.one.two.comparison[[seg.name]][[samp.two]]$mean_ph_cpnA <- mean(as.numeric(phased.s.two.s.cpnA.cpnB.df[,1]))
               s.one.two.comparison[[seg.name]][[samp.two]]$mean_ph_cpnB <- mean(as.numeric(phased.s.two.s.cpnA.cpnB.df[,2]))
               s.one.two.comparison[[seg.name]][[samp.two]]$un_ph_cpn_total <- NA
               s.one.two.comparison[[seg.name]][[samp.two]]$un_ph_cpn_total_mean <- NA

               s.one.two.comparison[[seg.name]][[samp.one.two.string]]$un_ph_cpnA_cpnB_het_comb_paired_pval    <- NA
               s.one.two.comparison[[seg.name]][[samp.one.two.string]]$un_ph_cpnA_cpnB_het_comb_unpaired_pval  <- NA

               s.one.two.comparison[[seg.name]][[samp.one.two.string]]$ph_cpnA_het_paired_pval <- seg.paired.p.vals.vec[1]
               s.one.two.comparison[[seg.name]][[samp.one.two.string]]$ph_cpnB_het_paired_pval <- seg.paired.p.vals.vec[2]
               s.one.two.comparison[[seg.name]][[samp.one.two.string]]$ph_cpnA_het_unpaired_pval <- seg.unpaired.p.vals.vec[1]
               s.one.two.comparison[[seg.name]][[samp.one.two.string]]$ph_cpnB_het_unpaired_pval <- seg.unpaired.p.vals.vec[2]
               s.one.two.comparison[[seg.name]][[samp.one.two.string]]$ph_cpnA_cpnB_het_comb_paired_pval <- seg.comb.paired.p.val
               s.one.two.comparison[[seg.name]][[samp.one.two.string]]$ph_cpnA_cpnB_het_comb_unpaired_pval <- seg.comb.unpaired.p.val
               s.one.two.comparison[[seg.name]][[samp.one.two.string]]$ph_tot_cpn_paired_pval <- total.cpn.paired.pval
               s.one.two.comparison[[seg.name]][[samp.one.two.string]]$ph_cpnA_diff_by_size_het <- seg.cpnA.het.size.by.diff
               s.one.two.comparison[[seg.name]][[samp.one.two.string]]$ph_cpnB_diff_by_size_het <- seg.cpnB.het.size.by.diff
               s.one.two.comparison[[seg.name]][[samp.one.two.string]]$comb_cpnA_cpnB_diff_by_size_het <- seg.cpnA.het.size.by.diff + seg.cpnB.het.size.by.diff
               s.one.two.comparison[[seg.name]][[samp.one.two.string]]$MSAI <- 0


               if (tTestWithErrorHandling(phased.s.one.s.cpnA.cpnB.df[, 1], phased.s.one.s.cpnA.cpnB.df[, 2], paired=TRUE)$p.val < 0.001 &
                    tTestWithErrorHandling(phased.s.two.s.cpnA.cpnB.df[, 1], phased.s.two.s.cpnA.cpnB.df[,2], paired=TRUE)$p.val < 0.001) {

                    if( (mean(phased.s.one.s.cpnA.cpnB.df[, 1] - phased.s.one.s.cpnA.cpnB.df[, 2]) > 0 & 
                         mean(phased.s.two.s.cpnA.cpnB.df[, 1] - phased.s.one.s.cpnA.cpnB.df[, 2]) < 0) | 
                        (mean(phased.s.one.s.cpnA.cpnB.df[, 2] - phased.s.one.s.cpnA.cpnB.df[, 1]) < 0 & 
                         mean(phased.s.two.s.cpnA.cpnB.df[, 2] - phased.s.one.s.cpnA.cpnB.df[, 1]) > 0) ){

                         s.one.two.comparison[[seg.name]][[samp.one.two.string]]$MSAI <- 1

                    }

               }#end if (tTestWithErrorHandling(phased.s.one.s.cpnA.cpnB.df[, 1], phased.s.one.s.cpnA.cpnB.df[, 2], paired=TRUE)$p.val < 0.001 &
                    #tTestWithErrorHandling(phased.s.two.s.cpnA.cpnB.df[, 1], phased.s.two.s.cpnA.cpnB.df[,2], paired=TRUE)$p.val < 0.001) {


               gain.samp.one.df <- checkSegGainOrLossSignificance(s.one.s.ao.lr.baf.df, samp.one, ploidy.est.vec, purity.est.vec, gain.threshold, "gain")
               amp.samp.one.df <- checkSegGainOrLossSignificance(s.one.s.ao.lr.baf.df, samp.one, ploidy.est.vec, purity.est.vec, amp.threshold, "gain")
               loss.samp.one.df <- checkSegGainOrLossSignificance(s.one.s.ao.lr.baf.df, samp.one, ploidy.est.vec, purity.est.vec, loss.threshold, "loss")
               deep.loss.samp.one.df <- checkSegGainOrLossSignificance(s.one.s.ao.lr.baf.df, samp.one, ploidy.est.vec, purity.est.vec, deep.loss.threshold, "loss")
               samp.one.cpn.event.vec <- c(gain.samp.one.df[, 3:5, drop=FALSE]
                                        ,gain.samp.one.df[, 2, drop=FALSE], gain.samp.one.df[, 7, drop=FALSE]
                                        ,amp.samp.one.df[, 2,drop=FALSE], amp.samp.one.df[, 7, drop=FALSE]
                                        ,loss.samp.one.df[, 2, drop=FALSE], loss.samp.one.df[, 7, drop=FALSE]
                                        ,deep.loss.samp.one.df[, 2, drop=FALSE], deep.loss.samp.one.df[, 7, drop=FALSE])
               names(samp.one.cpn.event.vec) <- cpn.event.names.vec

               s.one.two.comparison[[seg.name]][[samp.one]]$seg_cpn_event_status_vec <- samp.one.cpn.event.vec

               gain.samp.two.df <- checkSegGainOrLossSignificance(s.one.s.ao.lr.baf.df, samp.two, ploidy.est.vec, purity.est.vec, gain.threshold, "gain")
               amp.samp.two.df <- checkSegGainOrLossSignificance(s.one.s.ao.lr.baf.df, samp.two, ploidy.est.vec, purity.est.vec, amp.threshold, "gain")
               loss.samp.two.df <- checkSegGainOrLossSignificance(s.one.s.ao.lr.baf.df, samp.two, ploidy.est.vec, purity.est.vec, loss.threshold, "loss")
               deep.loss.samp.two.df <- checkSegGainOrLossSignificance(s.one.s.ao.lr.baf.df, samp.two, ploidy.est.vec, purity.est.vec, deep.loss.threshold, "loss")
               samp.two.cpn.event.df <- c(gain.samp.two.df[, 3:5, drop=FALSE]
                                        ,gain.samp.two.df[, 2, drop=FALSE], gain.samp.two.df[, 7, drop=FALSE]
                                        ,amp.samp.two.df[, 2,drop=FALSE], amp.samp.two.df[, 7, drop=FALSE]
                                        ,loss.samp.two.df[, 2, drop=FALSE], loss.samp.two.df[, 7, drop=FALSE]
                                        ,deep.loss.samp.two.df[, 2, drop=FALSE], deep.loss.samp.two.df[, 7, drop=FALSE])
               samp.two.cpn.event.vec <- cpn.event.names.vec

               s.one.two.comparison[[seg.name]][[samp.two]]$seg_cpn_event_status_vec <- samp.two.cpn.event.vec


          }#end if if (nrow(logRBAFSample1) > min.num.hz.snp.phasing)

	}#end for ( i in 1:length(chr.start.end.string.vec) )

     average.pairwise.het.over.all.segs <- sum(pw.CN.het.score.vec) / sum(seg.size.vec)

     s.one.two.comparison[[samp.one.two.string]]$pairwise_avg_cn_het_score <- average.pairwise.het.over.all.segs

     return(s.one.two.comparison)

}#End calculatePairwiseCNHeterogeneity <- function(samp.one                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            