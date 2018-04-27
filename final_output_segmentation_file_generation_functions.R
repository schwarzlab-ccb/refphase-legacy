



generateFinalSegmentationOutput <- function(all.pw.cn.het.comp.lists
                                 ,sample.to.phase.seg.list
                                 ,all.samp.ph.BAF.df.list      
                                 ,all.seg.CCF.vals.list
                                 ,processed.orig.cpn.segs.df
                                 ,min.cons.segs.ploidy.adj.cpn.df
                                 ,ploidy.est.vec
                                 ,purity.est.vec
                                 ,clonal.pval = 0.001
                                 ,diff.from.oth.p.val = 0.01
                                 ,min.CCF.to.consider.val =0.05
                                 ,pval.to.separate = 0.001
                                 ,min.num.other.chrs = 3
                                 ,gain.threshold    = log2(2.5 / 2)
                                 ,loss.threshold    = log2(1.5 / 2)
                                 ,amp.threshold     = log2(4 / 2)
                                 ,deep.loss.threshold= log2(0.5 / 2)
                                 ,chrs.to.use.vec = 1:22
                                 ,ploidy = 'mean'
                                 ,loss.gain.pval =0.01  
                                ) {
    mu.ccf.val <- 1 - min.CCF.to.consider.val
    #Note that the min.cons.segs.df passed to this function was generated on the ploidy_adj_tot_cpn column
    #Is ths really necessary? If we're using ploidy adjusted total cpn this might produce slightly smaller joint segs df but  why?
    #Todo, pull this out of this function and put in main()
    min.cons.segs.raw.cpnA.df <- generateMinimumConsistentSegmentation(processed.orig.cpn.segs.df, chroms=1:22, col.to.use="nAraw")
    min.cons.segs.raw.cpnB.df <- generateMinimumConsistentSegmentation(processed.orig.cpn.segs.df, chroms=1:22, col.to.use="nBraw")
    #If the rownames are exactly the same, what's the point? original code uses Araw for both...
    rownames(min.cons.segs.raw.cpnA.df) <- paste(min.cons.segs.raw.cpnA.df[,1],min.cons.segs.raw.cpnA.df[,2],min.cons.segs.raw.cpnA.df[,3],sep=":") 
    rownames(min.cons.segs.raw.cpnB.df) <- paste(min.cons.segs.raw.cpnA.df[,1],min.cons.segs.raw.cpnA.df[,2],min.cons.segs.raw.cpnA.df[,3],sep=":")

    sample.vec <- names(ploidy.est.vec)
    first.sample <- sample.vec[1]
    #Note we may be missing a segment that fell through, what about those without eno
    orig.seg.name.vec <- paste(min.cons.segs.ploidy.adj.cpn.df[, 1], min.cons.segs.ploidy.adj.cpn.df[, 2], min.cons.segs.ploidy.adj.cpn.df[, 3], sep=":")
    seg.name.vec <- orig.seg.name.vec[orig.seg.name.vec %in% names(all.pw.cn.het.comp.lists[[1]])]

    if (length(orig.seg.name.vec) != length(seg.name.vec)) {
            print(orig.seg.name.vec[ ! orig.seg.name.vec %in% seg.name.vec])
            print("not present in the all.pw.cn.het.comp.lists[[1]] and excluded.")
    } else {#end  if (length(orig.seg.name.vec) != length(seg.name.vec)) {

            print("All segments present in min.cons.segs.ploidy.adj.cpn.df present in all.pw.cn.het.comp.lists[[1]]")
    } 


    #Now on a per sample basis:

    all.samp.new.cpn.vec <- c()
    all.seg.df <- initialiseAllSegFullDataframe(seg.name.vec, sample.vec)
    print(all.seg.df)

    for ( samp in sample.vec ) {

            new.cpn.vec <- c()

            for ( seg.name in seg.name.vec ) {

                    samp.seg.name <- paste0(samp, ":", seg.name)

                    seg.info.fss.list <- all.pw.cn.het.comp.lists[[first.sample]][[seg.name]]
                    seg.info.css.list <- all.pw.cn.het.comp.lists[[samp]][[seg.name]]
                    css.cpn.event.status.vec <- all.pw.cn.het.comp.lists[[samp]][[seg.name]]$seg_cpn_event_status_vec
                    chr.start.end.vec <- unlist(strsplit(seg.name,split=":"))
                    all.seg.df[samp.seg.name,]$chr <- seg.info.fss.list$chr[1]
                    all.seg.df[samp.seg.name,]$startpos <- chr.start.end.vec[2]  
                    all.seg.df[samp.seg.name,]$endpos <- chr.start.end.vec[3]
                    all.seg.df[samp.seg.name,]$cnTotal <- css.cpn.event.status.vec$median_seg_logr          
                    all.seg.df[samp.seg.name,]$Ploidy <- ploidy.est.vec[samp]
                    all.seg.df[samp.seg.name,]$ACF <- purity.est.vec[samp]
                    all.seg.df[samp.seg.name,]$orig_cpnA_raw <- min.cons.segs.raw.cpnA.df[seg.name, samp]
                    all.seg.df[samp.seg.name,]$orig_cpnB_raw <- min.cons.segs.raw.cpnB.df[seg.name, samp]
                    all.seg.df[samp.seg.name,]$median_logr <- css.cpn.event.status.vec["median_seg_logr"]
                    all.seg.df[samp.seg.name,]$gain_pval <- css.cpn.event.status.vec["gain_pval"]
                    all.seg.df[samp.seg.name,]$amp_pval <-  css.cpn.event.status.vec["amp_pval"]
                    all.seg.df[samp.seg.name,]$loss_pval <- css.cpn.event.status.vec["loss_pval"]
                    all.seg.df[samp.seg.name,]$deep_loss_pval <- css.cpn.event.status.vec["deep_loss_pval"]
                    switch.comp.to.ref <- NA

                    if ( is.na(seg.info.fss.list$ph_cpnA_vec[1])  | is.na(seg.info.fss.list$ph_cpnB_vec[1])) {

                            all.seg.df[samp.seg.name,]$samp_used_to_phase <- "too_few_SNPs"

                            next

                    } else {#end if ( is.na(seg.info.fss.list$ph_cpnA_vec[1]) ) {

                            #Look for MSAI
                            samp.used.to.phase <- sample.to.phase.seg.list[[seg.name]]$samp_chosen_to_phase_seg
                            all.seg.df[samp.seg.name,]$samp_used_to_phase <- samp.used.to.phase
                            print(all.pw.cn.het.comp.lists[[paste0(samp, ":", samp.used.to.phase)]][[seg.name]]$MSAI)
                            switch.comp.to.ref <- all.pw.cn.het.comp.lists[[paste0(samp, ":", samp.used.to.phase)]][[seg.name]]$MSAI
                            print("switch.comp.to.ref")
                            print(switch.comp.to.ref)
                            if (length(switch.comp.to.ref) != 0) {
                                if ( ! is.null(switch.comp.to.ref) ){
                                    if (switch.comp.to.ref == 1){
                                        all.seg.df[samp.seg.name,]$MSAI <- TRUE 
                                    }
                                    if(switch.comp.to.ref == 0){
                                        all.seg.df[samp.seg.name,]$MSAI <- FALSE  
                                    }
                             }
                         }
                    }

                    seg.CCF.val.list <- all.seg.CCF.vals.list
                    #Todo refacor to allow 
                    all.seg.df[samp.seg.name,]$n.het <- length(seg.info.css.list$ph_cpnA_vec)
                    all.seg.df[samp.seg.name,]$nBraw <- ifelse(median(as.numeric(seg.info.css.list$ph_cpnB_vec)) < 0 , 0, median(as.numeric(seg.info.css.list$ph_cpnB_vec)))
                    all.seg.df[samp.seg.name,]$nAraw <- ifelse(median(as.numeric(seg.info.css.list$ph_cpnA_vec)) < 0 , 0, median(as.numeric(seg.info.css.list$ph_cpnA_vec)))
                    all.seg.df[samp.seg.name,]$nMinor <-  round(as.numeric(all.seg.df[samp.seg.name,]$nBraw))
                    all.seg.df[samp.seg.name,]$nMajor <-  round(as.numeric(all.seg.df[samp.seg.name,]$nAraw))
                    all.seg.df[samp.seg.name,]$cnTotal <- as.numeric(all.seg.df[samp.seg.name,]$nMinor) + as.numeric(all.seg.df[samp.seg.name,]$nMajor)
                    ss.ph.baf.df <- all.samp.ph.BAF.df.list[[samp]]
                    all.seg.df[samp.seg.name,]$median_BAF  <- ss.ph.baf.df[ss.ph.baf.df$chrom == as.numeric(chr.start.end.vec[1]) &
                                                                            ss.ph.baf.df$pos >= as.numeric(chr.start.end.vec[2]) &
                                                                            ss.ph.baf.df$pos <= as.numeric(chr.start.end.vec[3]),,drop=FALSE]$ph_median_segmented_BAF[1]
 
            
                    for ( cpnX in c("cpnA", "cpnB") ) {

                        all.seg.df[samp.seg.name, paste0("poss_sc_",cpnX, "_low")] <- seg.CCF.val.list[[seg.name]][[paste0("poss_sc_",cpnX, "_low")]]
                        all.seg.df[samp.seg.name, paste0("poss_sc_",cpnX, "_high")] <- seg.CCF.val.list[[seg.name]][[paste0("poss_sc_",cpnX, "_high")]]
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
}                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              