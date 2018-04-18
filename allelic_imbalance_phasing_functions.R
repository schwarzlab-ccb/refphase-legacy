# options(scipen=999);
# library(compiler)
# library(parallel)
# library(data.table)
# library(mclust)
# library(GenomicRanges)
# library(Rsamtools)
# library(gtools)
# library(diptest)
# library(gplots)





##
##Params:
##
##Returns:
##
getCombinedBAFDipTestPhasing <- function(asas.lr.baf.lists,
                                        min.cons.segs.df,
                                        ploidy.est.vec,
                                        purity.est.vec,
                                        gamma.val=1,
                                        min.num.hz.snp.phasing=5,
                                        min.prop.snp.sep.threshold=0.75,
                                        dip.test.p.val.threshold=0.01,
                                        comb.wilcox.snp.sep.p.val.threshold=0.01) {

  sample.vec <- names(asas.lr.baf.lists)
  sample.to.phase.seg.list <- list()
  new.baf.df.list <- list()
  chr.start.end.string.vec <- paste(min.cons.segs.df[ ,1], min.cons.segs.df[ ,2], min.cons.segs.df[ ,3], sep=":")

  for ( seg.idx in 1:length(chr.start.end.string.vec) ) {

    seg.name <- chr.start.end.string.vec[seg.idx]
    #Test whether phasing is possible for this segment:
    first.ss.lr.baf.df <- asas.lr.baf.lists[[sample.vec[1]]][[seg.idx]]

    if (nrow(first.ss.lr.baf.df) < min.num.hz.snp.phasing) {

      next;#Might actally want this to be a break as they should all have too few snps as it is joint segmentation with the same SNPs...

    }#end if( nrow(ref.ss.lr.baf.df) < min.hz.snp.phasing)

    all.samp.seg.phasing.ability.list <- getAssessmentOfEachSamplesAbilityToPhaseAllOthersForSeg(asas.lr.baf.lists=asas.lr.baf.lists,
                                                                                                sample.vec=sample.vec,
                                                                                                seg.idx=seg.idx,
                                                                                                ploidy.est.vec=ploidy.est.vec,
                                                                                                purity.est.vec=purity.est.vec,
                                                                                                gamma.val=gamma.val,
                                                                                                min.num.hz.snp.phasing=min.num.hz.snp.phasing)
    #save to output list put for this seg:
    sample.to.phase.seg.list[[seg.name]]$all_samp_phasing_comparison_list <- all.samp.seg.phasing.ability.list

    #all samples summary phasing df  
    a.s.sum.phas.df <- data.frame(mean_prop_SNPs_sep_across_samples = sapply(all.samp.seg.phasing.ability.list, function(x) x$mean_prop_SNPs_sep_across_samples),
                                            diptest_for_AI_pval = sapply(all.samp.seg.phasing.ability.list, function(x) x$diptest_for_AI_pval),
                                            ASCAT_called_as_AI = sapply(all.samp.seg.phasing.ability.list, function(x) x$ASCAT_called_as_AI),
                                            comb_pval_SNP_sep_across_samps = sapply(all.samp.seg.phasing.ability.list, function(x) x$comb_wilcox_sep_pval_across_samps),
                                            stringsAsFactors=FALSE)
    rownames(a.s.sum.phas.df) <- sample.vec
    a.s.sum.phas.df <- a.s.sum.phas.df[with(a.s.sum.phas.df, order(diptest_for_AI_pval)), ]
    #possible phasing samples passing the phasing thresholds:
    poss.phas.s.df <- a.s.sum.phas.df[a.s.sum.phas.df$ASCAT_called_as_AI==TRUE | a.s.sum.phas.df$diptest_for_AI_pval < dip.test.p.val.threshold, ,drop=FALSE]

    # first scenario - there is evidence of allelic imbalance (from either ASCAT or dip.test)
    if (nrow(poss.phas.s.df) > 0) {

      #If ASCAT only chose one region then choose this region as the reference
      if (nrow(poss.phas.s.df[poss.phas.s.df$ASCAT_called_as_AI==TRUE, ,drop=FALSE]) == 1) {

        sample.to.phase.seg <- rownames(poss.phas.s.df[poss.phas.s.df$ASCAT_called_as_AI==TRUE, ,drop=FALSE])[1] 
        sample.to.phase.seg.list[[seg.name]]$samp_chosen_to_phase_seg <- sample.to.phase.seg
        sample.to.phase.seg.list[[seg.name]]$maj_is_cpnA_idx_vec <- all.samp.seg.phasing.ability.list[[sample.to.phase.seg]]$maj_is_cpnA_idx_vec
        sample.to.phase.seg.list[[seg.name]]$maj_is_cpnB_idx_vec <- all.samp.seg.phasing.ability.list[[sample.to.phase.seg]]$maj_is_cpnB_idx_vec

      }#end if (length(all.samp.s.ASCAT.called.AI.logical.vec) == 0)

      #If ASCAT calls more than one region to use as reference:
      if (nrow(poss.phas.s.df[poss.phas.s.df$ASCAT_called_as_AI==TRUE, ,drop=FALSE]) > 1) {

        #Choose the maximum phasing region from already sorted by maximum proportion of snps separated dataframe:
        sample.to.phase.seg <- rownames(poss.phas.s.df)[1]


        sample.to.phase.seg.list[[seg.name]]$samp_chosen_to_phase_seg <- sample.to.phase.seg
        sample.to.phase.seg.list[[seg.name]]$maj_is_cpnA_idx_vec <- all.samp.seg.phasing.ability.list[[sample.to.phase.seg]]$maj_is_cpnA_idx_vec
        sample.to.phase.seg.list[[seg.name]]$maj_is_cpnB_idx_vec <- all.samp.seg.phasing.ability.list[[sample.to.phase.seg]]$maj_is_cpnB_idx_vec

        #We assume we can improve the phasing by using a combination from all regions if the phasing from this sample it falls below a threshold:
        if (poss.phas.s.df[sample.to.phase.seg, ,drop=FALSE]$comb_pval_SNP_sep_across_samps < comb.wilcox.snp.sep.p.val.threshold &
            poss.phas.s.df[sample.to.phase.seg, ,drop=FALSE]$mean_prop_SNPs_sep_across_samples < min.prop.snp.sep.threshold) {

          #Call a function to do mixture phasing if we want:
          sample.to.phase.seg.list <- getMixturePhasingResultsForSegment(sample.to.phase.seg.list,
                                                                        sample.to.phase.seg,
                                                                        sample.vec,
                                                                        poss.phas.s.df,
                                                                        asas.lr.baf.lists,
                                                                        seg.idx, ploidy.est.vec, purity.est.vec, gamma.val)

        }#end if(poss.phas.s.df[sample.to.phase.seg,]$comb_wilcox_sep_pval_across_samps < comb.wilcox.snp.sep.p.val.threshold &
         #poss.phas.s.df[sample.to.phase.seg,]$mean_prop_SNPs_sep_across_samples > prop.snp.sep.threshold)

      }#end if (nrow(poss.phas.s.df[poss.phas.s.df$ASCAT_called_as_AI==TRUE,]) > 1)

    }#end if (nrow(poss.phas.s.df) > 0)
    # second scenario - there is no evidence of allelic imbalance (from ASCAT or dip.test) # in this case -randomly assign SNPs as major or minor
    if (nrow(poss.phas.s.df) == 0 ) {

      possible.maj.cpnA.maj.cpnB.idxs.df <- c(all.samp.seg.phasing.ability.list[[1]]$maj_is_cpnA_idx_vec, all.samp.seg.phasing.ability.list[[1]]$maj_is_cpnB_idx_vec) 
      major.random.idxs <- sample(possible.maj.cpnA.maj.cpnB.idxs.df, size=length(possible.maj.cpnA.maj.cpnB.idxs.df) / 2)
      minor.random.idxs <- possible.maj.cpnA.maj.cpnB.idxs.df[ ! possible.maj.cpnA.maj.cpnB.idxs.df %in% major.random.idxs]

      sample.to.phase.seg.list[[seg.name]]$samp_chosen_to_phase_seg <- 'none_snps_randomized'
      sample.to.phase.seg.list[[seg.name]]$maj_is_cpnA_idx_vec <- sort(major.random.idxs)
      sample.to.phase.seg.list[[seg.name]]$maj_is_cpnB_idx_vec <- sort(minor.random.idxs)

    }#end if (nrow(poss.phas.s.df) == 0) {

  }#end for ( seg.idx in 1:length(chr.start.end.vec) )

  return(sample.to.phase.seg.list)
}#end getCombinedBAFDipTestPhasing




##
##Params:
##
##Returns:
##
calculateSNPsRawCopyNumberAB <- function(logr.vec, baf.vec, tumor.ploidy, tumor.purity, gamma.val = 1) {

  #print("funCalcN_withBAF")
  cpnA <- (tumor.purity - 1 + baf.vec * 2 ^ (logr.vec / gamma.val) * ((1 - tumor.purity) * 2 + tumor.purity * tumor.ploidy)) / tumor.purity
  cpnB <- (tumor.purity - 1 - (baf.vec - 1) * 2 ^ (logr.vec / gamma.val) * ((1 - tumor.purity) * 2 + tumor.purity * tumor.ploidy)) / tumor.purity
  #print("Returning from funCalcN_withBAF")
  return(cbind(cpnA,cpnB))

}#end calculateSNPCopyNumberAB <- function(




##
##Params:
##
##Returns:
##
assessPhasingAcrossSelectedSamples <- function(sample.comparison.list,
                                              samp,
                                              asas.lr.baf.lists,
                                              samples.to.test.vec,
                                              ref.ss.maj.is.cpnA.idx.vec,
                                              ref.ss.maj.is.cpnB.idx.vec,
                                              ploidy.est.vec,
                                              purity.est.vec) {

  other.ss.phasing.df <- data.frame(separation_wilcox_pval=rep(1, length(samples.to.test.vec)),
                                    prop_separated=numeric(length(samples.to.test.vec)),
                                    concordant_w_ref=numeric(length(samples.to.test.vec)),
                                    stringsAsFactors=FALSE)

  for ( other.samp in samples.to.test.vec ) {

    oth.ss.lr.baf.df <- asas.lr.baf.lists[[other.samp]][[seg.idx]]
     <- calculateSNPsCopyNumberAB(logr.vec=oth.ss.lr.baf.df$logR, baf.vec=oth.ss.lr.baf.df$BAF, tumor.ploidy=ploidy.est.vec[other.samp], tumor.purity=purity.est.vec[other.samp], gamma.val=gamma.val)

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
                                                                    seg.idx,
                                                                    ploidy.est.vec,
                                                                    purity.est.vec,
                                                                    gamma.val=1,
                                                                    min.num.hz.snp.phasing=5){

  #other sample segment phasing df, used when working out the optimum sample in all regions, for each segment the same dataframe is reused to save time allocating memory...
  print(sample.vec)
  num.other.samps <- length(sample.vec) - 1
  other.ss.phasing.df <- data.frame(separation_wilcox_pval=rep(1, num.other.samps),
                                    prop_separated=numeric(num.other.samps),
                                    concordant_w_ref=numeric(num.other.samps),
                                    stringsAsFactors=FALSE,
                                    drop=FALSE)
  sample.comparison.list <- list()

  for ( samp in sample.vec ) {

      #ref sample seg logr baf df
      ref.ss.lr.baf.df <- asas.lr.baf.lists[[samp]][[seg.idx]]
      if( nrow(ref.ss.lr.baf.df) < min.num.hz.snp.phasing) {

        next;#Might actally want this to be a break as they should all have too few snps as it is joint segmentation with the same SNPs...

      }#end if( nrow(ref.ss.lr.baf.df) < min.hz.snp.phasing)

      #renaming the dataframe rows for the other samples, now we are trying a new reference sample
      rownames(other.ss.phasing.df) <- sample.vec[samp != sample.vec]

      #Get phasing:
      #reference sample segment cpn A cpn B
      print(ploidy.est.vec)
      print(purity.est.vec)
      #Todo getting the major and A and B columns can be made into a function...
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

      samples.to.test.vec <- sample.vec[ ! samp%in%sample.vec]

      sample.comparison.list <- assessPhasingAcrossSelectedSamples(sample.comparison.list = sample.comparison.list,
                                                                  samp = samp,
                                                                  asas.lr.baf.lists = asas.lr.baf.lists,
                                                                  samples.to.test.vec = samples.to.test.vec,
                                                                  ref.ss.maj.is.cpnA.idx.vec = ref.ss.maj.is.cpnA.idx.vec,
                                                                  ref.ss.maj.is.cpnB.idx.vec = ref.ss.maj.is.cpnB.idx.vec,
                                                                  ploidy.est.vec = ploidy.est.vec,
                                                                  purity.est.vec = purity.est.vec)

  }#end for( samp in sample.vec ) 

  return(sample.comparison.list)
}#end of getAssessmentOfEachSamplesAbilityToPhaseAllOthersForSeg <- function(asas.lr.baf.lists, sample.vec, seg.idx)

##Current issues, allowing samples through questionably
##


##
##Params:
##
##Returns:
##
getMixturePhasingResultsForSegment <- function(sample.to.phase.seg.list,
                                              sample.to.phase.seg,
                                              sample.vec,
                                              poss.phas.s.df,
                                              asas.lr.baf.lists,
                                              seg.idx,
                                              ploidy.est.vec,
                                              purity.est.vec,
                                              gamma.val) {

  #Wipe the crrent sample.to.seg.phase list and update to reflect the choice to use a mixture of all the samples:
  sample.to.phase.seg.list[[seg.name]]$samp_chosen_to_phase_seg <- "mixture"
  #We will update the rest at the end of working out the optimum major and minor cpn indices...

  phas.ss.lr.baf.df <- asas.lr.baf.lists[[sample.to.phase.seg]][[seg.idx]]
  phas.ss.cpnA.cpnB.df <- calculateSNPsCopyNumberAB(logr.vec=phas.ss.lr.baf.df$logR, baf.vec=phas.ss.lr.baf.df$BAF, tumor.ploidy=ploidy.est.vec[sample.to.phase.seg], tumor.purity=purity.est.vec[sample.to.phase.seg], gamma.val=gamma.val)
  
  other.samps.vec <- sample.vec[ ! sample.to.phase.seg %in% sample.vec]

  oth.s.cnpA.cpnB.match.ref.vec <- all.samp.seg.phasing.ability.list[sample.to.phase.seg]$oth_samp_cpnA_cpnB_match_ref_vec
  names(oth.s.cnpA.cpnB.match.ref.vec) <- other.samps.vec

  other.passing.samples.vec <- other.samps.vec[other.samps.vec%in%rownames(poss.phas.s.df)]
  all.passing.samples.vec <- c(sample.to.phase.seg, other.passing.samples.vec)

  #check whether should be considered positive or negative
  new.maj.is.cpnA.idx.vec  <- c()
  new.maj.is.cpnB.idx.vec  <- c()
  change.from.ref <- c()

  #Now, for each SNP work out whether it should be major or minor across all samples...?
  for ( SNP.idx in 1:nrow(phas.ss.cpnA.cpnB.df) ) {

    if (SNP.idx%in%all.samp.seg.phasing.ability.list[[sample.to.phase.seg]]$maj_is_cpnA_idx_vec) {
      current.maj.cpn.col.idx <- 1
      current.min.cpn.col.idx <- 2
    }#end if (SNP.idx%in%all.samp.seg.phasing.ability.list[[sample.to.phase.seg]]$maj_is_cpnA_idx_vec)

    if (SNP.idx%in%all.samp.seg.phasing.ability.list[[sample.to.phase.seg]]$maj_is_cpnB_idx_vec) {
      current.maj.cpn.col.idx <- 2
      current.min.cpn.col.idx <- 1
    }# end if (SNP.idx%in%all.samp.seg.phasing.ability.list[[sample.to.phase.seg]]$maj_is_cpnB_idx_vec)
    
    all.s.SNP.sep.scores <- c()

    for ( samp in all.passing.samples.vec ) {

      ss.lr.baf.df <- asas.lr.baf.lists[[samp]][[seg.idx]]
      ss.cpnA.cpnB.df <- calculateSNPsCopyNumberAB(logr.vec=ss.lr.baf.df$logR, baf.vec=ss.lr.baf.df$BAF, tumor.ploidy=ploidy.est.vec[samp], tumor.purity=purity.est.vec[samp], gamma.val=gamma.val)

      ss.SNP.cpnA <- ss.cpnA.cpnB.df[SNP.idx, current.maj.cpn]
      ss.SNP.cpnB <- ss.cpnA.cpnB.df[SNP.idx, current.min.cpn]

      s.cpnAB.match.ref.factor <- oth.s.cnpA.cpnB.match.ref.vec[samp]
      #sample snp sep score
      s.SNP.sep.score <- (ss.SNP.cpnA - ss.SNP.cpnB) * s.cpnAB.match.ref.factor#Couldn't we just absolute this? why do we need multiply it by a factor.

      all.s.SNP.sep.scores <- c(all.s.SNP.sep.scores, s.SNP.sep.score)

    }#end for ( samp in all.passing.samples.vec ) {

    #
    if (mean(all.SNP.sep.score) >= 0) {

      change.from.ref <- c(change.from.ref, 0)

      if (current.maj.cpn.col.idx == 1) {

        new.maj.is.cpnA.idx.vec <- c(new.maj.is.cpnA.idx.vec, SNP.idx)

      } else if (current.maj.cpn.col.idx == 2) {

        new.maj.is.cpnB.idx.vec <- c(new.maj.is.cpnB.idx.vec, SNP.idx)

      }

    }#end if (mean(all.SNP.sep.score) > 0) 


    if (mean(all.SNP.sep.score) < 0) {

      change.from.ref <- c(change.from.ref, 1)

      if (current.maj.cpn.col.idx == 1) {

        new.maj.is.cpnB.idx.vec <- c(new.maj.is.cpnB.idx.vec, SNP.idx)

      } else if (current.maj.cpn.col.idx == 2) {

        new.maj.is.cpnA.idx.vec <- c(new.maj.is.cpnA.idx.vec, SNP.idx)

      }

    }#end if (mean(all.SNP.sep.score) < 0) 

  }#end for ( SNP.idx in 1:nrow(phas.ss.cpnA.cpnB.df) )

  #So have now come up with whether whether cpnA or cpnB are the Major cpn through all samples allowing for MSAI

  sample.to.phase.seg.list[[seg.name]]$maj_is_cpnA_idx_vec <- new.maj.is.cpnA.idx.vec
  sample.to.phase.seg.list[[seg.name]]$maj_is_cpnB_idx_vec <- new.maj.is.cpnB.idx.vec

  #Now check the ability of the combined mixture snp major cpn A or B analysis accross all samples
  sample.to.phase.seg.list[[seg.name]]$all_samp_mixture_phasing_comparison_list <- assessPhasingAcrossSelectedSamples(sample.comparison.list=list(), 
                                                                                                              sample,
                                                                                                              asas.lr.baf.lists,
                                                                                                              sample.vec,
                                                                                                              ref.ss.maj.is.cpnA.idx.vec,
                                                                                                              ref.ss.maj.is.cpnB.idx.vec,
                                                                                                              ploidy.est.vec,
                                                                                                              purity.est.vec)


  return(sample.to.phase.seg.list)
}




#This is used to generate CCF values so needs to be added to the sample.to.phase.seg.list
createAllSegsPhasedBAFDataFrame <- function(asas.lr.baf.lists, sample.to.phase.seg.list) {

  sample.vec <- names(asas.lr.baf.lists)
  chr.start.end.string.vec <- paste(min.cons.segs.df[ ,1], min.cons.segs.df[ ,2], min.cons.segs.df[ ,3], sep=":")
  all.samp.ph.BAF.df.list <- list()

  for ( i in 1:length(chr.start.end.string.vec) ) {

    seg.name <- chr.start.end.string.vec[seg.idx]

    for ( samp in sample.vec ) {

      samp.BAF.df <- asas.lr.baf.lists[[samp]][[seg.name]]
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

      if (nrow(samp.BAF.df[is.na(samp.BAF.df$maj_is_cpnA_idx_vec), ,drop=FALSE])>0)){

        stop("Check error in sample BAF df generation")
        
      }#end if (nrow(sample.BAF.df[is.na(sample.BAF.df$maj_is_cpnA_idx_vec), ,drop=FALSE])>0)){

      samp.BAF.df$ph_median_segmented_BAF <- min(median(samp.BAF.df[samp.BAF.df$maj_is_cpnB_idx_vec==1,]$BAF),
                                           median(samp.BAF.df[samp.BAF.df$maj_is_cpnB_idx_vec==0,]$BAF))

      all.samp.BAF.df.list[[samp]]$ph_samp_BAF_df <- rbind(all.samp.BAF.df.list[[samp]]$ph_samp_BAF_df, samp.BAF.df)

    }#end for ( samp in sample.vec ) {

  }#end for ( i in 1:length(chr.start.end.string.vec) ) {

  return(all.samp.ph.BAF.df.list) 

}#end createSegPhasedBAFTable <- function(....)





                         