




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
                                )

{
        #Todo add error checking
        #processed.orig.cpn.segs.df looks like this:
        #           sample chr  startpos    endpos n.het cnTotal nMajor nMinor   Ploidy  ACF    nAraw    nBraw patient sample_name tot_raw_cpn max_raw_cpn min_raw_cpn cpn_to_use ploidy_adj_tot_cpn
        # 1 B_MET002_MT-M1   1    762496  23885674 16393       2      1      1 1.885366 0.31 1.082364 1.082364  MET002   MET002:M1    2.164729    1.082364    1.082364   2.164729          0.1993419
        # 2 B_MET002_MT-M1   1  23885742  34666653  6453       1      1      0 1.885366 0.31 1.035798 0.000000  MET002   MET002:M1    1.035798    1.035798    0.000000   1.035798         -0.8641015

        #define 
        
        mu <- 1 - min.CCF.to.consider.val
        #Note that the min.cons.segs.df passed to this function was generated on the ploidy_adj_tot_cpn column
        #Is ths really necessary? If we're using ploidy adjusted total cpn this might produce slightly smaller joint segs df but  why? 
        min.cons.segs.raw.cpnA.df <- generateMinimumConsistentSegmentation(processed.orig.cpn.segs.df, chroms=1:22, col.to.use="nAraw")
        min.cons.segs.raw.cpnB.df <- generateMinimumConsistentSegmentation(processed.orig.cpn.segs.df, chroms=1:22, col.to.use="nBraw")
        #If the rownames are exactly the same, what's the point? original code uses Araw for both...
        rownames(min.cons.segs.raw.cpnA.df) <- paste(min.cons.segs.raw.cpnA.df[,1],min.cons.segs.raw.cpnA.df[,2],min.cons.segs.raw.cpnA.df[,3],sep=":") 
        rownames(min.cons.segs.raw.cpnB.df) <- paste(min.cons.segs.raw.cpnA.df[,1],min.cons.segs.raw.cpnA.df[,2],min.cons.segs.raw.cpnA.df[,3],sep=":")

        sample.vec <- names(ploidy.est.vec)
        first.sample <- sample.vec[1]
        #Note we may be missing a segment that fell through, what about those without eno
        orig.seg.name.vec <- paste(min.cons.segs.ploidy.adj.cpn.df[,1], min.cons.segs.ploidy.adj.cpn.df[,2], min.cons.segs.ploidy.adj.cpn.df[,3], sep=":")
        print("orig.seg.name.vec")
        print(orig.seg.name.vec)
        seg.name.vec <- orig.seg.name.vec[orig.seg.name.vec %in% names(all.pw.cn.het.comp.lists[[1]])]
        print("seg.name.vec")
        print(seg.name.vec)

        if (length(orig.seg.name.vec) != length(seg.name.vec)) {
                print(orig.seg.name.vec[ ! orig.seg.name.vec %in% seg.name.vec])
                print("not present in the all.pw.cn.het.comp.lists[[1]] and excluded.")
        } else {#end  if (length(orig.seg.name.vec) != length(seg.name.vec)) {

                print("All segments present in min.cons.segs.ploidy.adj.cpn.df present in all.pw.cn.het.comp.lists[[1]]")
        } 


        #Now on a per sample basis:

        all.samp.new.cpn.vec <- c()
        all.seg.df <- initialiseAllSegFullDataframe(seg.name.vec)

        for ( samp in sample.vec ) {

                new.cpn.vec <- c()

                for ( seg.name in seg.name.vec ) {

                        samp.seg.name <- paste0(samp, ":", seg.name)

                        seg.info.fss.list <- all.pw.cn.het.comp.lists[[first.sample]][[seg.name]]
                        seg.info.css.list <- all.pw.cn.het.comp.lists[[samp]][[seg.name]]
                        css.cpn.event.status.vec <- all.pw.cn.het.comp.lists[[samp]][[seg.name]]$seg_cpn_event_status_vec
                        chr.start.end.vec <- unlist(strsplit(seg.name,split=":"))

                        all.seg.df[samp.seg.name,]$sample <- samp
                        all.seg.df[samp.seg.name,]$chr <- seg.info.fss.list$chr[1]
                        all.seg.df[samp.seg.name,]$startpos <- chr.start.end.vec[2]
                        all.seg.df[samp.seg.name,]$endpos <- chr.start.end.vec[3]
                        all.seg.df[samp.seg.name,]$cnTotal <- median(seg.info.css.list$cnTotal)             
                        all.seg.df[samp.seg.name,]$Ploidy <- ploidy.est.vec[samp]
                        all.seg.df[samp.seg.name,]$ACF <- purity.est.vec[samp]
                        all.seg.df[samp.seg.name,]$orig_cpnA_raw <- min.cons.segs.raw.cpnA.df[seg.name, samp]
                        all.seg.df[samp.seg.name,]$orig_cpnB_raw <- min.cons.segs.raw.cpnB.df[seg.name, samp]
                        all.seg.df[samp.seg.name,]$median_seg_logr <- css.cpn.event.status.vec["median_seg_logr"]
                        all.seg.df[samp.seg.name,]$gain_pval <- css.cpn.event.status.vec["gain_pval"]
                        all.seg.df[samp.seg.name,]$amp_pval <-  css.cpn.event.status.vec["amp_pval"]
                        all.seg.df[samp.seg.name,]$loss_pval <- css.cpn.event.status.vec["loss_pval"]
                        all.seg.df[samp.seg.name,]$deep_loss_pval <- css.cpn.event.status.vec["deep_loss_pval"]

                        if ( is.na(seg.info.fss.list$ph_cpnA_vec[1]) ) {

                                all.seg.df[samp.seg.name,]$samp_used_to_phase <- "too_few_SNPs"

                                next

                        } else {#end if ( is.na(seg.info.fss.list$ph_cpnA_vec[1]) ) {

                                #Look for MSAI 
                                samp.used.to.phase <- sample.to.phase.seg.list[[seg.name]]$samp_chosen_to_phase_seg
                                switch.comp.to.ref <- all.pw.cn.het.comp.lists[[paste0(samp, ":", samp.used.to.phase)]][[seg.name]]$MSAI

                                if ( length(switch.comp.to.ref) == 0 ) {

                                        switch.comp.to.ref <- all.pw.cn.het.comp.lists[[paste0(samp.used.to.phase, ":", samp)]][[seg.name]]$MSAI

                                        if (  length(switch.comp.to.ref) == 0 ) {

                                                switch.comp.to.ref <- NA

                                        }
                                }

                                all.seg.df[samp.seg.name,]$MSAI <- switch.comp.to.ref 

                        }

                        # css.ccf.list <- seg.CCF.val.list[[seg.name]][[samp]]
                        # all.samps.ccf.list <- seg.CCF.val.list[[seg.name]][[seg.name]]
                        # all.oth.samps.ccf.list <- seg.CCF.val.list[[seg.name]]

                        all.seg.df[samp.seg.name,]$n.het <- length(seg.info.css.list$ph_cpnA_vec)
                        all.seg.df[samp.seg.name,]$nBraw <- ifelse(median(as.numeric(seg.info.css.list$ph_cpnB_vec)) < 0 , 0, median(as.numeric(seg.info.css.list$ph_cpnB_vec)))
                        all.seg.df[samp.seg.name,]$nAraw <- ifelse(median(as.numeric(seg.info.css.list$ph_cpnA_vec)) < 0 , 0, median(as.numeric(seg.info.css.list$ph_cpnA_vec)))
                        all.seg.df[samp.seg.name,]$nMinor <-  round(all.seg.df[samp.seg.name,]$nBraw)
                        all.seg.df[samp.seg.name,]$nMajor <-  round(all.seg.df[samp.seg.name,]$nAraw )
                        all.seg.df[samp.seg.name,]$cnTotal <- all.seg.df[samp.seg.name,]$nMinor + all.seg.df[samp.seg.name,]$nMajor
                        all.seg.df[samp.seg.name,]$poss_sc_cpnA_low <- seg.CCF.val.list[[seg.name]][[samp]]$poss_sc_cpnA_low
                        all.seg.df[samp.seg.name,]$poss_sc_cpnA_high <- seg.CCF.val.list[[seg.name]][[samp]]$poss_sc_cpnA_high
                        all.seg.df[samp.seg.name,]$poss_sc_cpnB_low <- seg.CCF.val.list[[seg.name]][[samp]]$poss_sc_cpnB_low
                        all.seg.df[samp.seg.name,]$poss_sc_cpnB_high <- seg.CCF.val.list[[seg.name]][[samp]]$poss_sc_cpnB_high
                        all.seg.df[samp.seg.name,]$frac_psc_cpnA_low <- seg.CCF.val.list[[seg.name]][[samp]]$mean_frac_cpnA_low 
                        all.seg.df[samp.seg.name,]$frac_psc_cpnA_high <- seg.CCF.val.list[[seg.name]][[samp]]$mean_frac_cpnA_high
                        all.seg.df[samp.seg.name,]$frac_psc_cpnB_low <- seg.CCF.val.list[[seg.name]][[samp]]$mean_frac_cpnB_low
                        all.seg.df[samp.seg.name,]$frac_psc_cpnB_high <- seg.CCF.val.list[[seg.name]][[samp]]$mean_frac_cpnB_high

                        ss.ph.baf.df <- all.samp.ph.BAF.df.list[[samp]]
                        all.seg.df[samp.seg.name,]$median_BAF  <- ss.ph.baf.df[ss.ph.baf.df$chrom == chr & ss.ph.baf.df$pos >= startpos & ss.ph.baf.df$pos <= endpos,,drop=FALSE]$ph_median_segmented_BAF[1]


                        if ( median() ) {

                        }


                }#end for ( seg.name in seg.name.vec ) {

        }#end for ( samp in sample.vec ) {

}#end generateFinalSegmentationOutput <- function(...)




initialiseAllSegFullDataframe <- function(seg.name.vec, sample.vec) { 

        nr <- length(seg.name.vec) * length(sample.vec)
        all.seg.df <- data.frame(
                sample = rep(NA, nr),
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
                frac_psc_cpnB_low = rep(NA, nr),
                orig_cpnA_raw  = rep(NA, nr),
                orig_cpnB_raw = rep(NA, nr),
                sc_cpnA_present = rep(NA, nr),
                sc_cpnB_present = rep(NA, nr),
                het_cpnA_present = rep(NA, nr),
                het_cpnB_present = rep(NA, nr),
                median_BAF = rep(NA, nr),
                median_seg_logr = rep(NA, nr),
                gain_pval = rep(NA, nr),
                amp_pval = rep(NA, nr),
                loss_pval = rep(NA, nr),
                deep_loss_pval = rep(NA, nr),
                small_pval_cpnA = rep(NA, nr),
                small_pval_cpnB = rep(NA, nr),
                other_chrs_same_cpnA_cpn = rep(NA, nr),
                other_chrs_cpnA_pval = rep(NA, nr),
                other_chrs_same_cpnB_cpn = rep(NA, nr),
                other_chrs_cpnB_pval = rep(NA, nr),
                samp_used_to_phase  = rep(NA, nr),
                diff_from_ref_samp  = rep(NA, nr),
                stringsAsFactors=FALSE
                )
        row.names.vec <- c()
        for ( samp in sample.vec ) {
                row.names.vec <- c(row.names.vec, paste(samp, seg.name.vec, sep=":"))
        }
        rownames(all.seg.df) <- row.names.vec
        return(all.seg.df)
}                                                                                               