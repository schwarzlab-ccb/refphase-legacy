

addOriginalCPNValuesToOutputSegmentation <- function(sample.vec,
                                        seg.name.vec, 
                                        processed.orig.cpn.segs.df,
                                        all.segs.df) {

    min.cons.segs.raw.cpnA.df <- generateMinimumConsistentSegmentation(processed.orig.cpn.segs.df, chroms=1:22, col.to.use="nAraw")
    min.cons.segs.raw.cpnB.df <- generateMinimumConsistentSegmentation(processed.orig.cpn.segs.df, chroms=1:22, col.to.use="nBraw")
    min.cons.segs.nMinor.df <- generateMinimumConsistentSegmentation(processed.orig.cpn.segs.df, chroms=1:22, col.to.use="nMinor")
    min.cons.segs.nMajor.df <- generateMinimumConsistentSegmentation(processed.orig.cpn.segs.df, chroms=1:22, col.to.use="nMajor")

    rownames(min.cons.segs.raw.cpnA.df) <- paste(min.cons.segs.raw.cpnA.df[,1],min.cons.segs.raw.cpnA.df[,2],min.cons.segs.raw.cpnA.df[,3],sep=":") 
    rownames(min.cons.segs.raw.cpnB.df) <- paste(min.cons.segs.raw.cpnA.df[,1],min.cons.segs.raw.cpnA.df[,2],min.cons.segs.raw.cpnA.df[,3],sep=":")
    rownames(min.cons.segs.nMinor.df) <- paste(min.cons.segs.raw.cpnA.df[,1],min.cons.segs.raw.cpnA.df[,2],min.cons.segs.raw.cpnA.df[,3],sep=":") 
    rownames(min.cons.segs.nMajor.df) <- paste(min.cons.segs.raw.cpnA.df[,1],min.cons.segs.raw.cpnA.df[,2],min.cons.segs.raw.cpnA.df[,3],sep=":")

    for ( samp in sample.vec ) {

        for ( seg.name in seg.name.vec ) {

            if (length(sample.vec) > 1) {
                samp.seg.name <- paste0(samp, ":", seg.name)
            } else {
                samp.seg.name <- seg.name
            }
            all.segs.df[samp.seg.name,]$orig_cpnA_raw <- min.cons.segs.raw.cpnB.df[seg.name, samp]
            all.segs.df[samp.seg.name,]$orig_cpnB_raw <- min.cons.segs.raw.cpnA.df[seg.name, samp]
            all.segs.df[samp.seg.name,]$orig_nMinor <- min.cons.segs.nMinor.df[seg.name, samp]
            all.segs.df[samp.seg.name,]$orig_nMajor <- min.cons.segs.nMajor.df[seg.name, samp]
            all.segs.df[samp.seg.name,]$orig_cnTotal <- all.segs.df[samp.seg.name,]$orig_nMinor + all.segs.df[samp.seg.name,]$orig_nMajor


        }#end for ( seg.name in seg.name.vec ) {
    }#end for ( samp in sample.vec ) {

    return(all.segs.df)

}#end addOriginalRawValuesToOutputSegmentation



#Todo refactor and remove this function.... can probably just be done in a for loop in main after the addition of the CCF values..
addPhasingInformationToOutputSegmentation <- function(sample.vec,
                                        seg.name.vec, 
                                        sample.to.phase.seg.list,
                                        all.pw.cn.het.comp.lists,
                                        all.segs.df) {


        #Todo, refactor to make this function work on the basis of a single segment:
        for ( samp in sample.vec ) {

            for ( seg.name in seg.name.vec ) {

                if (length(sample.vec) > 1) {
                    samp.seg.name <- paste0(samp, ":", seg.name)
                } else {
                    samp.seg.name <- seg.name
                }


                if (length(sample.to.phase.seg.list[[seg.name]][["samp_chosen_to_phase_seg"]]) == 0) {

                    all.segs.df[samp.seg.name, ]$samp_used_to_phase <- "too_few_SNPs"

                } else {

                    samp.used.to.phase <-  sample.to.phase.seg.list[[seg.name]][["samp_chosen_to_phase_seg"]]

                    all.segs.df[samp.seg.name, ]$samp_used_to_phase <- samp.used.to.phase

                    #Check here if the comparison has been done in the reverse direction... not sure it has
                    switch.comp.to.ref.vec <- c(all.pw.cn.het.comp.lists[[paste0(samp, ":", samp.used.to.phase)]][[seg.name]]$MSAI,
                                             all.pw.cn.het.comp.lists[[paste0(samp.used.to.phase, ":", samp)]][[seg.name]]$MSAI)
                    switch.comp.to.ref <- switch.comp.to.ref.vec[!sapply(switch.comp.to.ref.vec,is.null)]
                    #If a sample was actually chose as reference it is 
                    if (is.null(switch.comp.to.ref)){

                        all.segs.df[samp.seg.name,]$MSAI <- FALSE 

                    } else if (is.na(switch.comp.to.ref)){

                        all.segs.df[samp.seg.name,]$MSAI <- "missing" 
                        
                    } else if (switch.comp.to.ref == 1){

                        all.segs.df[samp.seg.name,]$MSAI <- TRUE 
                        
                    }  else if(switch.comp.to.ref == 0) {
                        
                        all.segs.df[samp.seg.name,]$MSAI <- FALSE  
                    }
                }#end if (length(sample.to.phase.seg.list[[seg.name]][["samp_chosen_to_phase_seg"]]) == 0) {

                #Need to add the checkf for when there couldn't be any phasng?, for now rely on the the fact it's already there

            }#end for ( seg.name in seg.name.vec ) {
        }#end for ( samp in sample.vec ) {

        return(all.segs.df)
}#end addPhasingInformatonToOutputSegmentation <- function(...



generateFinalSegmentationOutput <- function(seg.name.vec
                                 ,all.samp.ph.BAF.df.list      
                                 ,all.seg.CCF.vals.list
                                 ,all.seg.CCF.diff.from.exp.list
                                 ,ploidy.est.vec
                                 ,purity.est.vec
                                 ,clonal.pval = 0.001
                                 ,diff.from.oth.p.val = 0.01
                                 ,min.CCF.to.consider.val = 0.15
                                 ,min.cpn.to.consider.val = 0.15
                                 ,min.num.other.chrs = 3
                                 ,gain.threshold    = log2(2.5 / 2)
                                 ,loss.threshold    = log2(1.5 / 2)
                                 ,amp.threshold     = log2(4 / 2)
                                 ,deep.loss.threshold= log2(0.5 / 2)
                                 ,tot.cpn.vs.ploidy.pval.thresh = 0.01
                                 ,chrs.to.use.vec = 1:22
                                 ,cpnX.vec = c("cpnA", "cpnB")
                                 ,reduced.seg.flag=FALSE
                                 ,max.cpn.val = 10
                                ) {


    sample.vec <- names(ploidy.est.vec)
    if(length(sample.vec) == 1){
        all.seg.df <- initialiseAllSegFullDataframe(seg.name.vec, sample.vec, make.row.names.seg.name.alone=TRUE)
    } else {
        all.seg.df <- initialiseAllSegFullDataframe(seg.name.vec, sample.vec)
    }

    for ( samp in sample.vec ) {

            for ( seg.name in seg.name.vec ) {

                    if (length(sample.vec) > 1) {
                        samp.seg.name <- paste0(samp, ":", seg.name)
                    } else {
                        samp.seg.name <- seg.name
                    }

                    all.seg.df[samp.seg.name,]$Ploidy <- ploidy.est.vec[samp]
                    all.seg.df[samp.seg.name,]$ACF <- purity.est.vec[samp]

                    chr.start.end.vec <- unlist(strsplit(seg.name, split=":"))
                    ss.ph.baf.df <- all.samp.ph.BAF.df.list[[samp]]
                    ss.baf.lr.baf.df <- ss.ph.baf.df[ss.ph.baf.df$chrom == as.numeric(chr.start.end.vec[1]) &
                                                    ss.ph.baf.df$pos >= as.numeric(chr.start.end.vec[2]) &
                                                    ss.ph.baf.df$pos <= as.numeric(chr.start.end.vec[3]),,drop=FALSE]
                    all.seg.df[samp.seg.name,]$median_BAF  <- ss.baf.lr.baf.df$ph_median_segmented_BAF[1]

                    gain.samp.df <- checkSegGainOrLossSignificance(ss.baf.lr.baf.df, samp, ploidy.est.vec, purity.est.vec, gain.threshold, "gain")
                    amp.samp.df <- checkSegGainOrLossSignificance(ss.baf.lr.baf.df, samp, ploidy.est.vec, purity.est.vec, amp.threshold, "gain")
                    loss.samp.df <- checkSegGainOrLossSignificance(ss.baf.lr.baf.df, samp, ploidy.est.vec, purity.est.vec, loss.threshold, "loss")
                    deep.loss.samp.df <- checkSegGainOrLossSignificance(ss.baf.lr.baf.df, samp, ploidy.est.vec, purity.est.vec, deep.loss.threshold, "loss")
                    all.seg.df[samp.seg.name,]$cnTotal <-  gain.samp.df[, 3]#at the moment these are both set to be median logr   
                    all.seg.df[samp.seg.name,]$median_logr <- gain.samp.df[, 3]#at the moment these are both set to be median logr
                    all.seg.df[samp.seg.name,]$gain_pval <- gain.samp.df[, 7]
                    all.seg.df[samp.seg.name,]$amp_pval <-  amp.samp.df[, 7]
                    all.seg.df[samp.seg.name,]$loss_pval <- loss.samp.df[, 7]
                    all.seg.df[samp.seg.name,]$deep_loss_pval <- deep.loss.samp.df[, 7]
                    
                    cpn.event.vs.ploidy.vec <- as.numeric(c(gain.samp.df[, 7],  amp.samp.df[, 7],  loss.samp.df[, 7], deep.loss.samp.df[, 7]))
                    names(cpn.event.vs.ploidy.vec) <- c("gain", "amp", "loss", "deep_loss")
                    non.na.seg.cpn.event <- cpn.event.vs.ploidy.vec[! is.na(cpn.event.vs.ploidy.vec)]
                    significant.seg.cpn.event <- non.na.seg.cpn.event[non.na.seg.cpn.event <= tot.cpn.vs.ploidy.pval.thresh]
                    #Error check here:
                    if (length(significant.seg.cpn.event) > 1) {
                        #Note need to check that this makes sense
                        print("More than one total cpn event vs ploidy is signficant by thresholds")
                        print(cpn.event.vs.ploidy.vec)
                        all.seg.df[samp.seg.name,]$cpn_event_vs_ploidy <- names(significant.seg.cpn.event)[1]

                    } else if (length(significant.seg.cpn.event) == 1) {

                        all.seg.df[samp.seg.name,]$cpn_event_vs_ploidy <- names(significant.seg.cpn.event)

                    } else {

                        all.seg.df[samp.seg.name,]$cpn_event_vs_ploidy <- "none"

                    }


                    #for now add cpnX initialisation
                    cpnX <- "cpnA"
                    if ( ! seg.name %in% names(all.seg.CCF.vals.list[[cpnX]]) ) {

                        all.seg.df[samp.seg.name,]$samp_used_to_phase <- "too_few_SNPs"
                        next

                    }
                    if ( length(all.seg.CCF.vals.list[[cpnX]][[seg.name]]) == 1 ) {

                        all.seg.df[samp.seg.name,]$samp_used_to_phase <- "too_few_SNPs"
                        next

                    }

                    for ( cpnX in cpnX.vec ) {

                        if(reduced.seg.flag == FALSE){

                            seg.CCF.val.list <- all.seg.CCF.vals.list[[cpnX]][[seg.name]][[samp]]
                            seg.diff.from.others.list <- all.seg.CCF.diff.from.exp.list[[cpnX]][[seg.name]][[samp]]

                        } else {

                            seg.CCF.val.list <- all.seg.CCF.vals.list[[seg.name]]
                            seg.diff.from.others.list <- all.seg.CCF.diff.from.exp.list#Todo, Fill in next... 
                            
                        }
                        
                        #These are currently assigned twice... wiping the original.. definitely needs refactoring:
                        ph.cpnA.vec <- seg.CCF.val.list[[paste0(cpnX, "_orig_ph_obs_cpnA_vec")]]
                        ph.cpnB.vec <- seg.CCF.val.list[[paste0(cpnX, "_orig_ph_obs_cpnB_vec")]]
                        ph.cpnX.vec <- seg.CCF.val.list[[paste0(cpnX, "_orig_ph_obs_", cpnX, "_vec")]]
                        #print(ph.cpnX.vec)
                        all.seg.df[samp.seg.name,]$n.het <- length(ph.cpnA.vec)
                        ss.mean.cpnX.obs <- mean(ph.cpnX.vec)
                        ss.mean.cpnX.obs <- ifelse(ss.mean.cpnX.obs < 0, 0, ss.mean.cpnX.obs)

                        all.seg.df[samp.seg.name,]$nBraw <- ifelse(median(as.numeric(ph.cpnB.vec)) < 0 , 0, median(as.numeric(ph.cpnB.vec)))
                        all.seg.df[samp.seg.name,]$nAraw <- ifelse(median(as.numeric(ph.cpnA.vec)) < 0 , 0, median(as.numeric(ph.cpnA.vec)))
                        all.seg.df[samp.seg.name,]$nMinor <-  round(as.numeric(all.seg.df[samp.seg.name,]$nBraw))
                        all.seg.df[samp.seg.name,]$nMajor <-  round(as.numeric(all.seg.df[samp.seg.name,]$nAraw))
                        all.seg.df[samp.seg.name,]$cnTotal <- as.numeric(all.seg.df[samp.seg.name,]$nMinor) + as.numeric(all.seg.df[samp.seg.name,]$nMajor)
                        #end of double assignment
                        all.seg.df[samp.seg.name, paste0("psc_",cpnX, "_low")] <- median(as.numeric(seg.CCF.val.list[[paste0("poss_sc_",cpnX, "_low")]]))
                        all.seg.df[samp.seg.name, paste0("psc_",cpnX, "_high")] <- median(as.numeric(seg.CCF.val.list[[paste0("poss_sc_",cpnX, "_high")]]))
           
                        ss.frac.poss.sc.cpnX.low.vec <- as.numeric(seg.CCF.val.list[[paste0("frac_", cpnX, "_low_vec")]])
                        ss.frac.poss.sc.cpnX.high.vec <- as.numeric(seg.CCF.val.list[[paste0("frac_", cpnX, "_high_vec")]])

                        all.seg.df[samp.seg.name, paste0("frac_psc_",cpnX,"_low")] <- median(as.numeric(seg.CCF.val.list[[paste0("mean_frac_",cpnX,"_low")]]))
                        all.seg.df[samp.seg.name, paste0("frac_psc_",cpnX,"_high")] <- median(as.numeric(seg.CCF.val.list[[paste0("mean_frac_",cpnX,"_high")]]))

                        #Test for evidence of subclonality within the segment for cpnA
                        mu.ccf.val <- 1 - min.CCF.to.consider.val
                        if (all.seg.df[samp.seg.name, paste0("psc_", cpnX, "_low")] == all.seg.df[samp.seg.name, paste0("psc_", cpnX, "_high")]) { 
                            #print("balanced cpn")
                            all.seg.df[samp.seg.name, paste0("sc_", cpnX, "_high_or_low_dominate")] <- NA
                            all.seg.df[samp.seg.name,  paste0("sc_", cpnX, "_vs_min_CCF_pval")] <- 1
                            all.seg.df[samp.seg.name, paste0("sc_", cpnX, "_vs_min_cpn_pval")] <- 1


                        }  else if ( median(ss.frac.poss.sc.cpnX.low.vec) > median(ss.frac.poss.sc.cpnX.high.vec) ) {
                            cpnA.test.subclonal <- wilcox.test(ss.frac.poss.sc.cpnX.low.vec, mu = mu.ccf.val, alternative = "less", exact = FALSE)
                            all.seg.df[samp.seg.name, paste0("sc_", cpnX, "_high_or_low_dominate")] <- "LOW"
                            all.seg.df[samp.seg.name, paste0("sc_", cpnX, "_vs_min_CCF_pval")] <- cpnA.test.subclonal$p.value

                            print("LOW")
                            print("median(ss.frac.poss.sc.cpnX.low.vec)")
                            print(median(ss.frac.poss.sc.cpnX.low.vec))
                            print("mean(ss.frac.poss.sc.cpnX.low.vec)")
                            print(mean(ss.frac.poss.sc.cpnX.low.vec))
                            print(cpnA.test.subclonal)

                            cpnA.test.subclonal <- wilcox.test(ss.frac.poss.sc.cpnX.low.vec, mu = mu.ccf.val, alternative = "less", exact = FALSE)
                            all.seg.df[samp.seg.name, paste0("sc_", cpnX, "_high_or_low_dominate")] <- "LOW"
                            all.seg.df[samp.seg.name, paste0("sc_", cpnX, "_vs_min_CCF_pval")] <- cpnA.test.subclonal$p.value


                            #Now let's try with min.cpn.to.consider.val
                            mu.cpn.val <- all.seg.df[samp.seg.name, paste0("psc_", cpnX, "_low")] + min.cpn.to.consider.val
                            cpnA.cpn.mu.test.subclonal <- wilcox.test(ph.cpnX.vec, mu = mu.cpn.val, alternative = "greater", exact = FALSE)
                            all.seg.df[samp.seg.name, paste0("sc_", cpnX, "_vs_min_cpn_pval")] <- cpnA.cpn.mu.test.subclonal$p.value




                            
                        } else {#end   if ( median(seg.CCF.val.list[[seg.name]][[samp]]$frac_cpnA_low_vec) > median(seg.CCF.val.list[[seg.name]][[samp]]$frac_cpnA_high_vec) ) {
                            cpnA.test.subclonal <- wilcox.test(ss.frac.poss.sc.cpnX.high.vec, mu = mu.ccf.val, alternative = "less", exact = FALSE)
                            all.seg.df[samp.seg.name, paste0("sc_", cpnX, "_high_or_low_dominate")] <- "HIGH"
                            all.seg.df[samp.seg.name, paste0("sc_", cpnX, "_vs_min_CCF_pval")] <- cpnA.test.subclonal$p.value
                            print("HIGH")
                            print("median(ss.frac.poss.sc.cpnX.high.vec)")
                            print(median(ss.frac.poss.sc.cpnX.high.vec))
                            print("mean(ss.frac.poss.sc.cpnX.high.vec)")
                            print(mean(ss.frac.poss.sc.cpnX.high.vec))
                            print(cpnA.test.subclonal)

                            mu.cpn.val <- all.seg.df[samp.seg.name, paste0("psc_", cpnX, "_high")] - min.cpn.to.consider.val
                            cpnA.cpn.mu.test.subclonal <- wilcox.test(ph.cpnX.vec, mu = mu.cpn.val, alternative = "less", exact = FALSE)
                            all.seg.df[samp.seg.name, paste0("sc_", cpnX, "_vs_min_cpn_pval")] <- cpnA.cpn.mu.test.subclonal$p.value

                        }
                        #readline(prompt="Press [enter] to continue")
                        #We need to refactor this magic number out of this
                        cpnX.pval.vs.other.chroms <- 1e-120
                        all.seg.df[samp.seg.name, paste0("num_other_chrs_same_", cpnX, "_cpn")] <- length(seg.diff.from.others.list[[paste0(cpnX, "_num_chrs_w_match_int_cpn_vec")]])
                        if (all.seg.df[samp.seg.name, paste0("num_other_chrs_same_", cpnX, "_cpn")] >= min.num.other.chrs) {

                            cpnX.pval.vs.other.chroms <- seg.diff.from.others.list[[paste0(cpnX,"_pval_vs_oth_cpn_match_chrs_raw_cpn")]]

                        }
                        all.seg.df[samp.seg.name, paste0("vs_other_chrs_same_", cpnX, "_pval")] <- cpnX.pval.vs.other.chroms 

                        if (all.seg.df[samp.seg.name, paste0("frac_psc_",cpnX,"_low")] > 1 |
                           all.seg.df[samp.seg.name, paste0("frac_psc_",cpnX,"_high")] > 1 |
                           all.seg.df[samp.seg.name, paste0("sc_",cpnX,"_vs_min_CCF_pval")] >= clonal.pval |
                           all.seg.df[samp.seg.name, paste0("vs_other_chrs_same_",cpnX,"_pval")] >= diff.from.oth.p.val) {

                            all.seg.df[samp.seg.name, paste0("within_s_sc_", cpnX, "_present")] <-  FALSE

                        } else {

                            all.seg.df[samp.seg.name, paste0("within_s_sc_", cpnX, "_present")] <- TRUE

                        }


                        if (all.seg.df[samp.seg.name, paste0("frac_psc_",cpnX,"_low")] > 1 |
                           all.seg.df[samp.seg.name, paste0("frac_psc_",cpnX,"_high")] > 1 |
                           all.seg.df[samp.seg.name, paste0("sc_",cpnX,"_vs_min_cpn_pval")] >= clonal.pval |
                           all.seg.df[samp.seg.name, paste0("vs_other_chrs_same_",cpnX,"_pval")] >= diff.from.oth.p.val) {

                            all.seg.df[samp.seg.name, paste0("within_s_sc_", cpnX, "_present_using_min_cpn_diff")] <-  FALSE

                        } else {

                            all.seg.df[samp.seg.name, paste0("within_s_sc_", cpnX, "_present_using_min_cpn_diff")] <- TRUE

                        }

                    }#end for ( cpnX in c("cpnA", "cpnB") ) {

            }#end for ( seg.name in seg.name.vec ) {

    }#end for ( samp in sample.vec ) {

    #Disease level stats filled out here, will separate this out into another function at a later date




    return(all.seg.df)

}#end generateFinalSegmentationOutput <- function(...)




# ##
# ##Params:
# ##
# ##Returns:
# ##
# calculateClonalSubclonalTotalCPNvsPloidy <- function(all.seg.df, sample.vec, seg.name.vec){


#     #So we need to sort out definitions here:
#     #Is there subconal cpn within a single region
#     for ( seg.name in seg.name.vec ) {

#         #Make the all the seg.names
#         samp.seg.name.vec <- paste0(sample.vec, ":", seg.name)

#         #So cpn events vs ploidy can be calculated independently of the phasing so should always be present..

#         #If there is a clonal copy number event relative to ploidy present for this seg, what is it

#         #At the moment we will just follow if all segments passed their own purity and ploidy adjusted thresholds,
#         all.seg.df[samp.seg.name.vec, "across_s_cl_vs_ploidy_present"] <- 
#         #If there are subclonal (or heterogeneous) copy number events relatve to ploidy present for this seg what are they
#         all.seg.df[samp.seg.name.vec, "across_s_sc_vs_ploidy_present"] <- 

#         #For MSAI and any allele specfic events we need phasing to be present, so check it
#         if (all.seg.df[samp.seg.name.vec[1],]$samp_used_to_phase == "too_few_SNPs" | grep(""))

#         #Calcuate cpnA and cpnB across sample comparisons
#         for ( cpnX in c("cpnA", "cpnB") ) {

#             #In any sample is there WITHIN sample subclonal copynumber
#             all.seg.df[samp.seg.name.vec, paste0("sc_", cpnX, "_within_any_samp_present")] <- any(all.seg.df[samp.seg.name.vec, paste0("sc_", cpnX, "_within_any_samp_present")]) 

#             #Is the allele specific INTEGER/clonal within sample copy number different BETWEEN/ACROSS samples
#             #make vector..

#             all.seg.df[samp.seg.name.vec, paste0("het_", cpnX, "withi_samps_present")] <- 

#         }#end or ( cpnX in c("cpnA", "cpnB") ) {

#     }#end for ( seg.name in seg.name.vec ) {

#     return(all.seg.df)
# }




#addSegmentAcrossAllSamplesClonalSubclonalCPNvsPloidy






testFinalOutputDataframe <- function(final.joint.segs.df) {

    #Tests that should hold for the final output:
    potential.msai.vec <- ! is.na(final.joint.segs.df$nAraw) & (final.joint.segs.df$nAraw < final.joint.segs.df$nBraw)
    if ( ! all(as.logical(final.joint.segs.df[potential.msai.vec == TRUE,]$MSAI) == TRUE) ) {

        stop("MSAI column not flled in correctly")

    }

    #Test original adding

    if ( ! all(final.joint.segs.df$orig_nMajor < final.joint.segs.df$orig_nMinor) ) {

        stop("Original nMajor not always greter than orginal nMinor")

    }#end if ( ! all(final.joint.segs.df$orig_nMajor < final.joint.segs.df$orig_nMinor) ) {


}#end testFinalOutputDataframe

initialiseAllSegFullDataframe <- function(seg.name.vec, sample.vec, make.row.names.seg.name.alone=FALSE) { 

        seg.chr.vec <- sapply(seg.name.vec, function(x) as.numeric(unlist(strsplit(x,split=":"))[1]))
        seg.startpos.vec <-  sapply(seg.name.vec, function(x) as.numeric(unlist(strsplit(x,split=":"))[2]))
        seg.endpos.vec <-  sapply(seg.name.vec, function(x) as.numeric(unlist(strsplit(x,split=":"))[3]))
        nr <- length(seg.name.vec) * length(sample.vec)
        sample.rows  <- c()
        for ( samp in sample.vec ) {
            sample.rows <- c(sample.rows, rep(samp, length(seg.name.vec)))
        }
        all.seg.df <- data.frame(
                #New processed calls we might want
                sample_name = sample.rows,
                seg_name = rep(seg.name.vec, length(sample.vec)),
                chr = rep(seg.chr.vec, length(sample.vec)),
                startpos = rep(seg.startpos.vec, length(sample.vec)),
                endpos = rep(seg.endpos.vec, length(sample.vec)),
                n.het = rep(NA, nr),
                cnTotal = rep(NA, nr),
                nMinor = rep(NA, nr),         
                nMajor = rep(NA, nr),
                Ploidy = rep(NA, nr),
                ACF = rep(NA, nr),
                nAraw = rep(NA, nr),
                nBraw = rep(NA, nr),
                within_s_sc_cpnA_present = rep(NA, nr),
                within_s_sc_cpnA_present_using_min_cpn_diff = rep(NA, nr),
                sc_cpnA_high_or_low_dominate = rep(NA, nr),
                psc_cpnA_low = rep(NA, nr),
                psc_cpnA_high = rep(NA, nr),
                within_s_sc_cpnB_present = rep(NA, nr),
                within_s_sc_cpnA_present_using_min_cpn_diff = rep(NA, nr),
                sc_cpnB_high_or_low_dominate = rep(NA, nr),
                psc_cpnB_low = rep(NA, nr),
                psc_cpnB_high = rep(NA, nr),
                #Copy number event calling:
                cpn_event_vs_ploidy = rep(NA, nr),
                median_BAF = rep(NA, nr),
                median_logr = rep(NA, nr),
                mean_BAF = rep(NA, nr),
                mean_logr = rep(NA, nr),
                stdev_BAF = rep(NA, nr),
                stdev_logr = rep(NA, nr),
                num_logr_only = rep(NA, nr),
                num_both_logr_BAF = rep(NA, nr),
                gain_pval = rep(NA, nr),
                amp_pval = rep(NA, nr),
                loss_pval = rep(NA, nr),
                deep_loss_pval = rep(NA, nr),
                #information for subclonal calling:
                ancestral_cpn_state = rep(NA, nr),#if lost or gained ccf different? what i it relative to, ancestral stated medicc
                loss_vs_gain_ccf = rep(NA, nr),
                sc_cpnA_vs_min_cpn_pval = rep(NA, nr),
                sc_cpnA_vs_min_CCF_pval = rep(NA, nr),
                frac_psc_cpnA_low = rep(NA, nr),
                frac_psc_cpnA_high = rep(NA, nr),
                vs_other_chrs_same_cpnA_pval = rep(NA, nr),
                num_other_chrs_same_cpnA_cpn = rep(NA, nr),
                sc_cpnB_vs_min_cpn_pval = rep(NA, nr),
                sc_cpnB_vs_min_CCF_pval = rep(NA, nr),
                frac_psc_cpnB_low = rep(NA, nr),
                frac_psc_cpnB_high = rep(NA, nr),
                vs_other_chrs_same_cpnB_pval = rep(NA, nr),
                num_other_chrs_same_cpnB_cpn = rep(NA, nr),
                #Link back to original segmentation
                orig_seg_name = rep(NA, nr),
                orig_nMinor = rep(NA, nr),
                orig_nMajor = rep(NA, nr),
                orig_cpnA_raw  = rep(NA, nr),
                orig_cpnB_raw = rep(NA, nr),
                orig_cnTotal= rep(NA, nr),
                nMinor_change = rep(NA, nr),
                nMajor_change = rep(NA, nr),
                cnTotal_change = rep(NA, nr),
                #joined
                #joined_arm_or_focal = rep(NA, nr),
                #Allelic imabalance information
                within_s_cl_AI_present = rep(NA, nr),
                within_s_sc_AI_present = rep(NA, nr), 
                within_s_cl_LOH_present =  rep(NA, nr),
                within_s_sc_LOH_present = rep(NA, nr),
                new_within_s_cl_AI_detected = rep(NA, nr),
                new_within_s_sc_AI_detected = rep(NA, nr),
                orig_part_of_arm_or_focal = rep(NA, nr),
                samp_used_to_phase  = rep(NA, nr),
                MSAI = rep(NA, nr),#Todo, decide whether this is used or diff from ref samp
                dominant_allele = rep(NA, nr),
                is_parallel_gain_loss = rep(NA, nr),
                #disease level stats
                any_across_s_MSAI = rep(NA, nr),
                any_within_s_sc_cpnA_present = rep(NA, nr),
                any_within_s_sc_cpnB_present = rep(NA, nr),
                across_s_cl_cpnA = rep(NA, nr),#
                across_s_sc_cpnA = rep(NA, nr),#
                across_s_cl_cpnB = rep(NA, nr),#
                across_s_sc_cpnB = rep(NA, nr),#
                across_s_cl_AI = rep(NA, nr),
                across_s_sc_AI = rep(NA, nr),
                across_s_cl_LOH = rep(NA, nr),
                across_s_sc_LOH = rep(NA, nr),
                across_s_cl_cpn_vs_ploidy_present = rep(NA, nr),
                across_s_sc_cpn_vs_ploidy_present = rep(NA, nr),
                #sample summary data that could go into a different dataframe
                orig_wgii = rep(NA, nr),
                wgii = rep(NA, nr),
                within_s_cl_wgii = rep(NA, nr),
                within_s_sc_wgii = rep(NA, nr),
                logr_data_quality = rep(NA, nr),
                TRACERx_het_measure = rep(NA, nr),
                mphase_het_measure = rep(NA, nr),
                #GD = rep(NA, nr),#don't use until we are actually using GD input for something.
                cpn_fit_quality = rep(NA, nr),
                #Purity adjusted thresholds for calling gains/losses
                purity_adj_logr_thresh_gain = rep(NA, nr),
                purity_adj_logr_thresh_loss = rep(NA, nr),
                purity_adj_logr_thresh_deep_loss = rep(NA, nr),
                purity_adj_logr_thresh_amp = rep(NA, nr),
                stringsAsFactors = FALSE
                )
        if (make.row.names.seg.name.alone == TRUE ) {
            #Todo, change logic here to something more sensible... this should only be allowed when there's only one sample
            row.names.vec <- seg.name.vec
        } else {
            row.names.vec <- c()
            for ( samp in sample.vec ) {
                    row.names.vec <- c(row.names.vec, paste(samp, seg.name.vec, sep=":"))
            }

        }

        rownames(all.seg.df) <- row.names.vec
        return(all.seg.df)

}


# #Original segmentation dataframe:
# all.seg.df <- data.frame(
#         sample_name = sample.rows,
#         seg_name = rep(seg.name.vec, length(sample.vec)),
#         chr = rep(seg.chr.vec, length(sample.vec)),
#         startpos = rep(seg.startpos.vec, length(sample.vec)),
#         endpos = rep(seg.endpos.vec, length(sample.vec)),
#         n.het = rep(NA, nr),
#         cnTotal = rep(NA, nr),
#         nMinor = rep(NA, nr),         
#         nMajor = rep(NA, nr),
#         Ploidy = rep(NA, nr),
#         ACF = rep(NA, nr),
#         nAraw = rep(NA, nr),
#         nBraw = rep(NA, nr),
#         sc_cpnA_present = rep(NA, nr),
#         sc_cpnB_present = rep(NA, nr),
#         poss_sc_cpnA_low = rep(NA, nr),
#         poss_sc_cpnA_high = rep(NA, nr),
#         poss_sc_cpnB_low = rep(NA, nr),
#         poss_sc_cpnB_high = rep(NA, nr),
#         frac_psc_cpnA_low = rep(NA, nr),
#         frac_psc_cpnA_high = rep(NA, nr),
#         frac_psc_cpnB_low = rep(NA, nr),
#         frac_psc_cpnB_high = rep(NA, nr),
#         orig_seg_name = rep(NA, nr),
#         orig_nMinor = rep(NA, nr),
#         orig_nMajor = rep(NA, nr),
#         orig_cpnA_raw  = rep(NA, nr),
#         orig_cpnB_raw = rep(NA, nr),
#         orig_cnTotal= rep(NA, nr),
#         median_BAF = rep(NA, nr),
#         median_logr = rep(NA, nr),
#         gain_pval = rep(NA, nr),
#         amp_pval = rep(NA, nr),
#         loss_pval = rep(NA, nr),
#         deep_loss_pval = rep(NA, nr),
#         sc_cpnA_high_or_low_dominate = rep(NA, nr),
#         sc_cpnA_vs_min_CCF_pval = rep(NA, nr),
#         sc_cpnB_high_or_low_dominate = rep(NA, nr),
#         sc_cpnB_vs_min_CCF_pval = rep(NA, nr),
#         num_other_chrs_same_cpnA_cpn = rep(NA, nr),
#         vs_other_chrs_same_cpnA_pval = rep(NA, nr),
#         num_other_chrs_same_cpnB_cpn = rep(NA, nr),
#         vs_other_chrs_same_cpnB_pval = rep(NA, nr),
#         samp_used_to_phase  = rep(NA, nr),
#         MSAI = rep(NA, nr),#Todo, decide whether this is used or diff from ref samp
#         dominant_allele = rep(NA, nr),
#         cpn_event_vs_ploidy = rep(NA, nr),
#         is_parallel_cpn_event = rep(NA, nr),
#         AI_present = rep(NA, nr),
#         LOH_present = rep(NA, nr),
#         new_AI_present = rep(NA, nr),
#         nMinor_change = rep(NA, nr),
#         nMajor_change = rep(NA, nr),
#         het_cpnA_acc_samps_present = rep(NA, nr),
#         het_cpnB_acc_samps_present = rep(NA, nr),
#         stringsAsFactors = FALSE
#         )








                                                                                                                                                                                                                                                                                                                                          