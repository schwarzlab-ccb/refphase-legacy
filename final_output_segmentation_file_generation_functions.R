

#Todo, don't change this as should be exactly the same as in pairwise_CN_heterogeneity_comparison_functions.R
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




generateFinalSegmentationOutput <- function(seg.name.vec
                                 ,all.samp.ph.BAF.df.list      
                                 ,all.seg.CCF.vals.list
                                 ,ploidy.est.vec
                                 ,purity.est.vec
                                 ,clonal.pval = 0.001
                                 ,diff.from.oth.p.val = 0.01
                                 ,min.CCF.to.consider.val = 0.05
                                 ,pval.to.separate = 0.001
                                 ,min.num.other.chrs = 3
                                 ,gain.threshold    = log2(2.5 / 2)
                                 ,loss.threshold    = log2(1.5 / 2)
                                 ,amp.threshold     = log2(4 / 2)
                                 ,deep.loss.threshold= log2(0.5 / 2)
                                 ,chrs.to.use.vec = 1:22
                                 ,loss.gain.pval = 0.01
                                 ,cpnX.vec = c("cpnA", "cpnB")
                                 ,reduced.seg.flag=FALSE
                                ) {


    mu.ccf.val <- 1 - min.CCF.to.consider.val
    #Note that the min.cons.segs.df passed to this function was generated on the ploidy_adj_tot_cpn column
    #Is ths really necessary? If we're using ploidy adjusted total cpn this might produce slightly smaller joint segs df but  why?
    #Todo, pull this out of this function and put in main()


    # #These really need to passed as parameters to the function... Can't really be refactored for the joined segments....? Could use the median over all
    # #segments that have been included...? I would ignore for the rejoined segments. 
    # min.cons.segs.raw.cpnA.df <- generateMinimumConsistentSegmentation(processed.orig.cpn.segs.df, chroms=1:22, col.to.use="nAraw")
    # min.cons.segs.raw.cpnB.df <- generateMinimumConsistentSegmentation(processed.orig.cpn.segs.df, chroms=1:22, col.to.use="nBraw")
    # #If the rownames are exactly the same, what's the point? original code uses Araw for both...
    # rownames(min.cons.segs.raw.cpnA.df) <- paste(min.cons.segs.raw.cpnA.df[,1],min.cons.segs.raw.cpnA.df[,2],min.cons.segs.raw.cpnA.df[,3],sep=":") 
    # rownames(min.cons.segs.raw.cpnB.df) <- paste(min.cons.segs.raw.cpnA.df[,1],min.cons.segs.raw.cpnA.df[,2],min.cons.segs.raw.cpnA.df[,3],sep=":")

    sample.vec <- names(ploidy.est.vec)
    print(sample.vec)
    #first.sample <- sample.vec[1]
    # #Note we may be missing a segment that fell through, what about those without eno
    # orig.seg.name.vec <- paste(min.cons.segs.ploidy.adj.cpn.df[, 1], min.cons.segs.ploidy.adj.cpn.df[, 2], min.cons.segs.ploidy.adj.cpn.df[, 3], sep=":")
    # seg.name.vec <- orig.seg.name.vec[orig.seg.name.vec %in% names(all.pw.cn.het.comp.lists[[1]])]

    # if (length(orig.seg.name.vec) != length(seg.name.vec)) {
    #         print(orig.seg.name.vec[ ! orig.seg.name.vec %in% seg.name.vec])
    #         print("not present in the all.pw.cn.het.comp.lists[[1]] and excluded.")
    # } else {#end  if (length(orig.seg.name.vec) != length(seg.name.vec)) {

    #         print("All segments present in min.cons.segs.ploidy.adj.cpn.df present in all.pw.cn.het.comp.lists[[1]]")
    # } 
    #seg.name.vec <- orig.seg.name.vec[orig.seg.name.vec %in% names(all.pw.cn.het.comp.lists[[1]])]
    #Now on a per sample basis:

    all.seg.df <- initialiseAllSegFullDataframe(seg.name.vec, sample.vec)

    for ( samp in sample.vec ) {

            for ( seg.name in seg.name.vec ) {

                    samp.seg.name <- paste0(samp, ":", seg.name)

                    chr.start.end.vec <- unlist(strsplit(seg.name,split=":"))
                    all.seg.df[samp.seg.name,]$chr <- chr.start.end.vec[1]
                    all.seg.df[samp.seg.name,]$startpos <- chr.start.end.vec[2]  
                    all.seg.df[samp.seg.name,]$endpos <- chr.start.end.vec[3]
                    all.seg.df[samp.seg.name,]$Ploidy <- ploidy.est.vec[samp]
                    all.seg.df[samp.seg.name,]$ACF <- purity.est.vec[samp]

                    ss.ph.baf.df <- all.samp.ph.BAF.df.list[[samp]]
                    ss.baf.lr.baf.df <- ss.ph.baf.df[ss.ph.baf.df$chrom == as.numeric(chr.start.end.vec[1]) &
                                                    ss.ph.baf.df$pos >= as.numeric(chr.start.end.vec[2]) &
                                                    ss.ph.baf.df$pos <= as.numeric(chr.start.end.vec[3]),,drop=FALSE]
                    all.seg.df[samp.seg.name,]$median_BAF  <- ss.baf.lr.baf.df$ph_median_segmented_BAF[1]
                    #Todo, this is the part that needs to be refactored...
                    #Might as well just calculate the gains and losses here rather than at any point...
                    # if (recalculate.cn.events == TRUE) {

                    gain.samp.df <- checkSegGainOrLossSignificance(ss.baf.lr.baf.df, samp.one, ploidy.est.vec, purity.est.vec, gain.threshold, "gain")
                    amp.samp.df <- checkSegGainOrLossSignificance(ss.baf.lr.baf.df, samp.one, ploidy.est.vec, purity.est.vec, amp.threshold, "gain")
                    loss.samp.df <- checkSegGainOrLossSignificance(ss.baf.lr.baf.df, samp.one, ploidy.est.vec, purity.est.vec, loss.threshold, "loss")
                    deep.loss.samp.df <- checkSegGainOrLossSignificance(ss.baf.lr.baf.df, samp.one, ploidy.est.vec, purity.est.vec, deep.loss.threshold, "loss")

                    all.seg.df[samp.seg.name,]$cnTotal <-  gain.samp.df[, 3]#at the moment these are both set to be median logr   
                    all.seg.df[samp.seg.name,]$median_logr <- gain.samp.df[, 3]#at the moment these are both set to be median logr
                    all.seg.df[samp.seg.name,]$gain_pval <- gain.samp.df[, 2]
                    all.seg.df[samp.seg.name,]$amp_pval <-  amp.samp.df[, 2]
                    all.seg.df[samp.seg.name,]$loss_pval <- loss.samp.df[, 2]
                    all.seg.df[samp.seg.name,]$deep_loss_pval <- deep.loss.samp.df[, 2]

                    # } else {

                    #     css.cpn.event.status.vec <- all.pw.cn.het.comp.lists[[samp]][[seg.name]]$seg_cpn_event_status_vec
                    #     all.seg.df[samp.seg.name,]$cnTotal <- css.cpn.event.status.vec["median_seg_logr"]          
                    #     all.seg.df[samp.seg.name,]$median_logr <- css.cpn.event.status.vec["median_seg_logr"]
                    #     all.seg.df[samp.seg.name,]$gain_pval <- css.cpn.event.status.vec["gain_pval"]
                    #     all.seg.df[samp.seg.name,]$amp_pval <-  css.cpn.event.status.vec["amp_pval"]
                    #     all.seg.df[samp.seg.name,]$loss_pval <- css.cpn.event.status.vec["loss_pval"]
                    #     all.seg.df[samp.seg.name,]$deep_loss_pval <- css.cpn.event.status.vec["deep_loss_pval"]

                    # }

                    # #These can't really be 
                    # all.seg.df[samp.seg.name,]$orig_cpnA_raw <- min.cons.segs.raw.cpnA.df[seg.name, samp]
                    # all.seg.df[samp.seg.name,]$orig_cpnB_raw <- min.cons.segs.raw.cpnB.df[seg.name, samp]
                    #print(seg.name)
                    # seg.info.fss.list <- all.pw.cn.het.comp.lists[[first.sample]][[seg.name]]
                    # seg.info.css.list <- all.pw.cn.het.comp.lists[[samp]][[seg.name]]

                    #So now we need to rewrite this section to also skip when passing the ccfs
                    #print(str(seg.info.fss.list))
                    #So this needs to have phasing information...
                    #Might be good to add phasing information separately to this data frame...

                    #So now the other way this could be checked is by looking to see if it's present in the ccflist
                    if ( ! seg.name %in% names(all.seg.CCF.vals.list) ) {

                        all.seg.df[samp.seg.name,]$samp_used_to_phase <- "too_few_SNPs"
                        next

                    }

                    #print(str(all.seg.CCF.vals.list[[seg.name]]))

                    if ( length(all.seg.CCF.vals.list[[seg.name]]) == 1 ) {

                        all.seg.df[samp.seg.name,]$samp_used_to_phase <- "too_few_SNPs"
                        next

                    }

                    if(reduced.seg.flag == FALSE){
                        seg.CCF.val.list <- all.seg.CCF.vals.list[[seg.name]][[samp]]
                    } else {
                        seg.CCF.val.list <- all.seg.CCF.vals.list[[seg.name]]
                    }

                    print("str(seg.CCF.val.list, max=1)")
                    print(str(seg.CCF.val.list, max=1))


                    for ( cpnX in cpnX.vec ) {

                        #These are currently assigned twice... wiping the original.. definitely needs refactoring:
                        ph.cpnA.vec <- seg.CCF.val.list[[seg.name]][[samp]][[paste0(cpnX,"_orig_ph_obs_cpnA_vec")]]
                        ph.cpnB.vec <- seg.CCF.val.list[[seg.name]][[samp]][[paste0(cpnX,"_orig_ph_obs_cpnB_vec")]]
                        print("head(ph.cpnA.vec)")
                        print(head(ph.cpnA.vec))
                        print(head(ph.cpnB.vec))
                        all.seg.df[samp.seg.name,]$n.het <- length(ph.cpnA.vec)
                        all.seg.df[samp.seg.name,]$nBraw <- ifelse(median(as.numeric(ph.cpnA.vec)) < 0 , 0, median(as.numeric(ph.cpnA.vec)))
                        all.seg.df[samp.seg.name,]$nAraw <- ifelse(median(as.numeric(ph.cpnB.vec)) < 0 , 0, median(as.numeric(ph.cpnB.vec)))
                        all.seg.df[samp.seg.name,]$nMinor <-  round(as.numeric(all.seg.df[samp.seg.name,]$nBraw))
                        all.seg.df[samp.seg.name,]$nMajor <-  round(as.numeric(all.seg.df[samp.seg.name,]$nAraw))
                        all.seg.df[samp.seg.name,]$cnTotal <- as.numeric(all.seg.df[samp.seg.name,]$nMinor) + as.numeric(all.seg.df[samp.seg.name,]$nMajor)
                        #end of double assignment

                        #Todo, we need to refactor
                        print(str(seg.CCF.val.list[[seg.name]][[samp]]))
                        all.seg.df[samp.seg.name, paste0("poss_sc_",cpnX, "_low")] <- median(seg.CCF.val.list[[seg.name]][[samp]][[paste0("poss_sc_",cpnX, "_low")]])
                        all.seg.df[samp.seg.name, paste0("poss_sc_",cpnX, "_high")] <- median(seg.CCF.val.list[[seg.name]][[samp]][[paste0("poss_sc_",cpnX, "_high")]])
                        all.seg.df[samp.seg.name, paste0("frac_psc_",cpnX,"_low")] <- seg.CCF.val.list[[seg.name]][[samp]][[paste0("mean_frac_",cpnX,"_low")]]
                        all.seg.df[samp.seg.name, paste0("frac_psc_",cpnX,"_high")] <- seg.CCF.val.list[[seg.name]][[samp]][[paste0("mean_frac_",cpnX,"_high")]]

                        #Test for evidence of subclonality within the segment for cpnA
                        if (all.seg.df[samp.seg.name, paste0("poss_sc_", cpnX, "_low")] == all.seg.df[samp.seg.name, paste0("poss_sc_", cpnX, "_low")]) {

                            all.seg.df[samp.seg.name, paste0("sc_", cpnX, "_high_or_low_dominate")] <- NA
                            all.seg.df[samp.seg.name,  paste0("sc_", cpnX, "_vs_min_CCF_pval")] <- 1

                        }  else if ( median(seg.CCF.val.list[[seg.name]][[samp]][[paste0("frac_",cpnX,"_low_vec")]]) > median(seg.CCF.val.list[[seg.name]][[samp]][[paste0("frac_",cpnX,"_high_vec")]]) ) {

                            cpnA.test.subclonal <- wilcox.test(seg.CCF.val.list[[seg.name]][[samp]][[paste0("frac_",cpnX,"_low_vec")]], mu = mu.ccf.val, alternative = "less", exact=FALSE)
                            all.seg.df[samp.seg.name, paste0("sc_", cpnX, "_high_or_low_dominate")] <- "HIGH"
                            all.seg.df[samp.seg.name, paste0("sc_", cpnX, "_vs_min_CCF_pval")] <- cpnA.test.subclonal$p.value
                            
                        } else {#end   if ( median(seg.CCF.val.list[[seg.name]][[samp]]$frac_cpnA_low_vec) > median(seg.CCF.val.list[[seg.name]][[samp]]$frac_cpnA_high_vec) ) {

                            cpnA.test.subclonal <- wilcox.test(seg.CCF.val.list[[seg.name]][[samp]][[paste0("frac_",cpnX,"_high_vec")]], mu = mu.ccf.val, alternative = "less", exact=FALSE)
                            all.seg.df[samp.seg.name, paste0("sc_", cpnX, "_high_or_low_dominate")] <- "LOW"
                            all.seg.df[samp.seg.name, paste0("sc_", cpnX, "_vs_min_CCF_pval")] <- cpnA.test.subclonal$p.value
                        }
                        
                        #Check if different to other chromosomes now
                        seg.diff.from.others.list <- all.seg.CCF.vals.list[[seg.name]]$seg_diff_from_expected_list[[samp]]
                        #Todo, make this a default parameter so can be set

                        cpnX.pval.vs.other.chroms <- 1e-120
                        all.seg.df[samp.seg.name, paste0("num_other_chrs_same_", cpnX, "_cpn")] <- length(seg.diff.from.others.list[[paste0(cpnX, "_num_chrs_w_match_int_cpn_vec")]])
                        if (length(seg.diff.from.others.list[[paste0(cpnX,"_num_chrs_w_match_int_cpn")]]) >= min.num.other.chrs) {

                            cpnX.pval.vs.other.chroms <- seg.diff.from.others.list[[paste0(cpnX,"_pval_vs_oth_cpn_match_chrs_raw_cpn")]]

                        }
                        all.seg.df[samp.seg.name, paste0("vs_other_chrs_same_", cpnX, "_pval")] <- cpnX.pval.vs.other.chroms 

                        #Todo, is it true that the hetA or hetB means that there is always some heterogeneity across the samples..., go back through code for edge cases.
                        all.seg.df[samp.seg.name, paste0("het_", cpnX, "_acc_samps_present")] <- all.seg.df[samp.seg.name,paste0("poss_sc_", cpnX, "_low")] != all.seg.df[samp.seg.name,paste0("poss_sc_",cpnX, "_high")]

                        if(all.seg.df[samp.seg.name, paste0("frac_psc_",cpnX,"_low")] > 1 |
                           all.seg.df[samp.seg.name, paste0("frac_psc_",cpnX,"_high")] > 1 |
                           all.seg.df[samp.seg.name, paste0("sc_",cpnX,"_vs_min_CCF_pval")] >= clonal.pval |
                           all.seg.df[samp.seg.name, paste0("vs_other_chrs_same_",cpnX,"_pval")] >= diff.from.oth.p.val) {

                            all.seg.df[samp.seg.name, paste0("sc_",cpnX,"_present")] <-  FALSE

                        } else {

                            all.seg.df[samp.seg.name, paste0("sc_",cpnX,"_present")] <- TRUE

                        }

                    }#end for ( cpnX in c("cpnA", "cpnB") ) {

            }#end for ( seg.name in seg.name.vec ) {

    }#end for ( samp in sample.vec ) {
        
    return(all.seg.df)

}#end generateFinalSegmentationOutput <- function(...)


initialiseAllSegFullDataframe <- function(seg.name.vec, sample.vec) { 

        nr <- length(seg.name.vec) * length(sample.vec)
        sample.rows  <- c()
        for ( samp in sample.vec ) {
            sample.rows <- c(sample.rows, rep(samp, length(seg.name.vec)))
        }
        all.seg.df <- data.frame(
                sample_name = sample.rows,
                seg_name= rep(seg.name.vec, length(sample.vec)),
                chr = rep(NA, nr),
                startpos = rep(NA, nr),
                endpos = rep(NA, nr),
                n.het = rep(NA, nr),
                cnTotal = rep(NA, nr),         
                nMajor = rep(NA, nr),
                nMinor = rep(NA, nr),
                Ploidy = rep(NA, nr),
                ACF = rep(NA, nr),
                nAraw = rep(NA, nr),
                nBraw = rep(NA, nr),
                psi = rep(NA, nr),
                poss_sc_cpnA_low = rep(NA, nr),
                poss_sc_cpnA_high = rep(NA, nr),
                poss_sc_cpnB_low = rep(NA, nr),
                poss_sc_cpnB_high = rep(NA, nr),
                frac_psc_cpnA_low = rep(NA, nr),
                frac_psc_cpnA_high = rep(NA, nr),
                frac_psc_cpnB_low = rep(NA, nr),
                frac_psc_cpnB_high = rep(NA, nr),
                orig_cpnA_raw  = rep(NA, nr),
                orig_cpnB_raw = rep(NA, nr),
                sc_cpnA_present = rep(NA, nr),
                sc_cpnB_present = rep(NA, nr),
                het_cpnA_acc_samps_present = rep(NA, nr),
                het_cpnB_acc_samps_present = rep(NA, nr),
                median_BAF = rep(NA, nr),
                median_logr = rep(NA, nr),
                gain_pval = rep(NA, nr),
                amp_pval = rep(NA, nr),
                loss_pval = rep(NA, nr),
                deep_loss_pval = rep(NA, nr),
                sc_cpnA_high_or_low_dominate = rep(NA, nr),
                sc_cpnA_vs_min_CCF_pval = rep(NA, nr),
                sc_cpnB_high_or_low_dominate = rep(NA, nr),
                sc_cpnB_vs_min_CCF_pval = rep(NA, nr),
                num_other_chrs_same_cpnA_cpn = rep(NA, nr),
                vs_other_chrs_same_cpnA_pval = rep(NA, nr),
                num_other_chrs_same_cpnB_cpn = rep(NA, nr),
                vs_other_chrs_same_cpnB_pval = rep(NA, nr),
                samp_used_to_phase  = rep(NA, nr),
                MSAI = rep(NA, nr),#Todo, decide whether this is used or diff from ref samp
                stringsAsFactors = FALSE
                )
        row.names.vec <- c()
        for ( samp in sample.vec ) {
                row.names.vec <- c(row.names.vec, paste(samp, seg.name.vec, sep=":"))
        }
        rownames(all.seg.df) <- row.names.vec
        return(all.seg.df)

}




                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           