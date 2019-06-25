  print("Called plotSampleLogRBAFCPNABsAndMSAIAcrossGenome....")
    print("params: samp")
    print(samp)
    seg.name.vec <- names(all.seg.MSAI.by.BAF.ltnpd.present.list)
    chr.vec <- as.numeric(sapply(seg.name.vec, function(x) strsplit(x, split = ":")[[1]][1]))
    start.vec <- as.numeric(sapply(seg.name.vec, function(x) strsplit(x, split = ":")[[1]][2]))
    end.vec <- as.numeric(sapply(seg.name.vec, function(x) strsplit(x, split = ":")[[1]][3]))

    chr.start.end.df <- data.frame(chrom = as.numeric(chr.vec), startpos = as.numeric(start.vec), endpos = as.numeric(end.vec), stringsAsFactors=FALSE)

    # samp.BAF.df <- patient.phased.BAF.df.list[[1]]
    # chr.start.end.df <- samp.BAF.df[, c("chrom", "pos", "pos")]
    # chr.start.end.df <- apply(chr.start.end.df, 2, as.numeric)

    chr.vec <- unique(chr.start.end.df[, "chrom"])

    chroms.length.vec <- getChromLengthsVec(chr.start.end.df)
    names(chroms.length.vec) <- chr.vec
    chrom.segs.to.plot.df <- getCumulativeChromLengthsForPlottingDF(chr.start.end.df, chroms.length.vec, chr.vec)
    rownames(chrom.segs.to.plot.df) <- seg.name.vec

    min.x               <- 0
    min.y               <- -0.5
    max.y               <- 8.5
    max.x               <- as.numeric(max(chrom.segs.to.plot.df[,3]))

    num.samples  <- length(sample.vec)
    first.samp <- sample.vec[1] 
    last.samp <- sample.vec[length(sample.vec)] 
    #sampJoinedSeg <- funJoinAdjacentSegments(copyITHRDS$sampComparisonList) is the equivalent
    cpnA.allele.color = '#de2d2699'
    cpnB.allele.color= '#3182bd99'
    #png(file = paste(PLOT_OUTPUT_DIR, patient, "_BAF_phased", ".png", sep=""), width = 1200, height = 1400)


    track.name.vec <- c("logr only sites", "logr hz sites", "BAF", "cpnA", "cpnB", "mcs allele specific cpn","final allele specific cpn", "original allele specific cpn")
    num.tracks <- length(track.name.vec)
    first.track <- track.name.vec[1]
    last.track <-  track.name.vec[2]

    par(mfrow = c(num.tracks, 1))

  for ( track.name in track.name.vec ) {

      print(track.name)
      #Set out border size for sub plot
      if (track.name == first.track) {

        par(mar = c(2, 3, 4, 3))

      } else if (track.name == last.track) {

        par(mar = c(4, 3, 2, 3))

      } else {

        par(mar = c(3, 3, 3, 3))

      }

      if (track.name == "logr only sites") {
        print("plotting logr only sites")

        min.x               <- 0
        min.y               <- -2
        max.y               <- 3
        max.x               <- as.numeric(max(chrom.segs.to.plot.df[, 3]))
        plot(1
             ,xlim = c(min.x, max.x)
             ,ylim = c(min.y, max.y)
             ,type = 'n'
             ,xaxt = 'n'
             ,yaxs = 'i'
             ,yaxt = 'n'
             ,xaxs = 'i'
             ,ylab = ""
             ,xlab = ""
             ,lwd = 2
             ,main = paste0(samp, " logR from non-hz sites")
             ,bty = "n")
        box()

        seg.name.vec <- names(asas.SNPs.w.both.lr.baf.snp.lists[[samp]]) #names(all.seg.CCF.vals.list)

        #print(seg.name.vec)
        for ( seg.name in seg.name.vec ) {
          #print(seg.name)
          #Will be NA if too few SNPs and can't be plotted, length of NA == 1, do this check to avoid warnings about checking only first element of list.


          chr.start.end.vec <- unlist(strsplit(seg.name, split = ":"))

          #print("print(str(seg.ccf.cpnX.frac.vec))")
          #print(str(seg.ccf.cpnX.frac.vec))
          seg.to.plot.orig.df <- asas.SNPS.lr.and.or.baf.lists.path[[samp]][[seg.name]]

          seg.to.plot.orig.df<-  seg.to.plot.orig.df[ ! seg.to.plot.orig.df$BAF %in% c(NA, 2),]
          seg.to.plot.df <- seg.to.plot.orig.df[, c("chrom", "pos", "pos", "logR", "logR", "logR","logRsegmented")]

          #print(head(seg.to.plot.df))
          if (nrow(seg.to.plot.df) == 0) {

            next

          }#end if(is.na(seg.to.plot.df[1,1])) {

          #toPlot <- cbind(logRBAFSample1$chrom,logRBAFSample1$pos,logRBAFSample1$pos,ABSample1,ABSample2)
          seg.to.plot.df[,1:3] <- getCumulativeChromLengthsForPlottingDF(apply(seg.to.plot.df[,1:3, drop=FALSE], 2, as.numeric), chroms.length.vec, chr.vec)

          seg.to.plot.df <- data.frame(seg.to.plot.df, stringsAsFactors = FALSE)
          colnames(seg.to.plot.df) <- c("chrom", "pos", "pos", "logR", "logR", "logR","logRsegmented")


          points(x = as.numeric(seg.to.plot.df$pos), y = as.numeric(seg.to.plot.df$logR)
                 ,col = "grey"
                 ,cex = .5
                 ,pch = 16)#ifelse(MSAIstatus==1,'#238b4550',ifelse(PValsPaired<0.001&totalPval<0.001,'#3182bd50','#fc927250')),cex=0.75,pch=16)
          par(lend = 1)

          segments(x0 = min(as.numeric(seg.to.plot.df$pos))
                   ,y0 = seg.to.plot.df[, "logRsegmented"]
                   ,x1 = max(as.numeric(seg.to.plot.df$pos))
                   ,y1 = seg.to.plot.df[, "logRsegmented"]
                   ,col = "#de2d26"
                   ,lwd = 3)
          
        }#end   for ( seg.name in seg.name.vec ) {

          axis(side = 2, las = 1)
          abline(h = seq((min.y + 0.5), (max.y - 0.5), by = 1), lty = 'dashed', lwd = 0.75, col = 'gray')
          
          mtext(side = 3
                ,at = getChromMidpointsVec(chrom.segs.to.plot.df[,1:3], unique(chr.start.end.df[,1]))
                ,text = sort(unique(chrom.segs.to.plot.df[,1]))
                ,cex = seq(0.6, 0.4, length.out = length(unique(chrom.segs.to.plot.df[,1])))
                ,line = -1
                #,lwd=0.5
          )
          
          for (i in sort(unique(as.numeric(chrom.segs.to.plot.df[,1]))))
          {
            abline(v = max(as.numeric(chrom.segs.to.plot.df[chrom.segs.to.plot.df[,1] == i, 3]))
                   ,lwd = 0.5)
            #,lty='dashed')
          }#en

          #box()
          

      } else if ( track.name == "logr hz sites") { 
        print("plotting hz logr sites")

          min.x               <- 0
          min.y               <- -2
          max.y               <- 3
          max.x               <- as.numeric(max(chrom.segs.to.plot.df[, 3]))
          plot(1
               ,xlim = c(min.x, max.x)
               ,ylim = c(min.y, max.y)
               ,type = 'n'
               ,xaxt = 'n'
               ,yaxs = 'i'
               ,yaxt = 'n'
               ,xaxs = 'i'
               ,ylab = ""
               ,xlab = ""
               ,lwd = 2
               ,main = paste0(samp, " logr from hz sites")
               ,bty = "n")
          box()

          seg.name.vec <- names(asas.SNPs.w.both.lr.baf.snp.lists[[samp]]) #names(all.seg.CCF.vals.list)

          #print(seg.name.vec)
          for ( seg.name in seg.name.vec ) {
            #print(seg.name)
            #Will be NA if too few SNPs and can't be plotted, length of NA == 1, do this check to avoid warnings about checking only first element of list.


            chr.start.end.vec <- unlist(strsplit(seg.name, split = ":"))

            #print("print(str(seg.ccf.cpnX.frac.vec))")
            #print(str(seg.ccf.cpnX.frac.vec))
            seg.to.plot.orig.df <- asas.SNPs.w.both.lr.baf.snp.lists[[samp]][[seg.name]]

            seg.to.plot.df <- seg.to.plot.orig.df[, c("chrom", "pos", "pos", "logR", "logR", "logR","logRsegmented")]

            #print(head(seg.to.plot.df))
            if (nrow(seg.to.plot.df) == 0) {

              next

            }#end if(is.na(seg.to.plot.df[1,1])) {

            #toPlot <- cbind(logRBAFSample1$chrom,logRBAFSample1$pos,logRBAFSample1$pos,ABSample1,ABSample2)
            seg.to.plot.df[,1:3] <- getCumulativeChromLengthsForPlottingDF(apply(seg.to.plot.df[,1:3, drop=FALSE], 2, as.numeric), chroms.length.vec, chr.vec)

            seg.to.plot.df <- data.frame(seg.to.plot.df, stringsAsFactors = FALSE)
            colnames(seg.to.plot.df) <- c("chrom", "pos", "pos", "logR", "logR", "logR","logRsegmented")


            points(x = as.numeric(seg.to.plot.df$pos), y = as.numeric(seg.to.plot.df$logR)
                   ,col = "grey"
                   ,cex = .5
                   ,pch = 16)#ifelse(MSAIstatus==1,'#238b4550',ifelse(PValsPaired<0.001&totalPval<0.001,'#3182bd50','#fc927250')),cex=0.75,pch=16)
            par(lend = 1)

            segments(x0 = min(as.numeric(seg.to.plot.df$pos))
                     ,y0 = seg.to.plot.df[, "logRsegmented"]
                     ,x1 = max(as.numeric(seg.to.plot.df$pos))
                     ,y1 = seg.to.plot.df[, "logRsegmented"]
                     ,col = "#de2d26"
                     ,lwd = 3)
            
          }#end   for ( seg.name in seg.name.vec ) {

            axis(side = 2, las = 1)
            abline(h = seq((min.y + 0.5), (max.y - 0.5), by = 1), lty = 'dashed', lwd = 0.75, col = 'gray')
            
            mtext(side = 3
                  ,at = getChromMidpointsVec(chrom.segs.to.plot.df[,1:3], unique(chr.start.end.df[,1]))
                  ,text = sort(unique(chrom.segs.to.plot.df[,1]))
                  ,cex = seq(0.6, 0.4, length.out = length(unique(chrom.segs.to.plot.df[,1])))
                  ,line = -1
                  #,lwd=0.5
            )
            
            for (i in sort(unique(as.numeric(chrom.segs.to.plot.df[,1]))))
            {
              abline(v = max(as.numeric(chrom.segs.to.plot.df[chrom.segs.to.plot.df[,1] == i, 3]))
                     ,lwd = 0.5)
              #,lty='dashed')
            }#en

            #box()
            

      } else if (track.name == "BAF") {

        print("plotting BAFs...")

          samp.BAF.df <- patient.phased.BAF.df.list[[samp]]
          samp.BAF.df <- samp.BAF.df[, c("chrom", "pos", "pos", "BAF","BAFsegmented","logR","logRsegmented", "samp","sample_phased_by_this_sample","maj_is_cpnA_logic_vec","maj_is_cpnB_logic_vec","ph_median_segmented_BAF")]
          samp.BAF.df.to.plot <- getCumulativeChromLengthsForPlottingDF(samp.BAF.df, chroms.length.vec, chr.vec)
          samp.BAF.df <- samp.BAF.df.to.plot

          samp.BAF.df$color_to_plot <- ifelse(samp.BAF.df$maj_is_cpnA_logic_vec == 1, cpnA.allele.color, cpnB.allele.color)
          samp.BAF.df$mphase_segmented_BAF <- 0.5

          if(nrow(samp.BAF.df[grepl("none", samp.BAF.df$sample_phased_by_this_sample),]) > 0) {

            samp.BAF.df[grepl("none", samp.BAF.df$sample_phased_by_this_sample),]$ph_median_segmented_BAF <- 0.5


          }

          rownames(samp.BAF.df) <- 1:nrow(samp.BAF.df)


          plot(samp.BAF.df[, 2], samp.BAF.df$BAF,
              col = samp.BAF.df$color_to_plot,
              xlim = c(min.x, max.x),
              pch = 16,
              ylim = c(0, 1),
              cex = 0.5,
              ,xaxt = 'n'
              ,yaxs = 'i'
              ,yaxt = 'n'
              ,xaxs = 'i'
              ,main = samp,
              ylab = "BAF",
              xlab = "",
              las = 1, 
              xaxt = 'n',
              bty = "n")
          
          mtext(side = 3,
                at = getChromMidpointsVec(chrom.segs.to.plot.df, unique(chr.start.end.df[,1])),
                text = sort(unique(chrom.segs.to.plot.df[, 1])),
                cex = seq(0.6, 0.4, length.out = length(unique(chrom.segs.to.plot.df[, 1]))),
                line = -1)

          msai.present.seg.color <- "red"
          non.msai.seg.color <- "yellow"
          problem.seg.color <- "black"

          reference.samp.color <- "green"
          non.reference.samp.color <- "blue"
          no.detection.possible.color <- "grey"
          mixture.phasing.color <- "pink"


          for ( chr in sort(unique(as.numeric(chrom.segs.to.plot.df[, 1])))) {

            #plot the segments of MSAI detection:
              chrom.seg.name.vec <- grep(paste0("^",chr,":"), names(all.seg.MSAI.by.BAF.ltnpd.present.list), value=TRUE)
              #print(chrom.seg.name.vec)

              for (chrom.seg.name in chrom.seg.name.vec) {

                #MSAI seg colors
                if ( length(all.seg.MSAI.by.BAF.ltnpd.present.list[[chrom.seg.name]]$MSAI_present_any_sample_by_thresh) > 0 ) {
                  #Color for the MSAI
                  #print(all.seg.MSAI.by.BAF.ltnpd.present.list[[chrom.seg.name]]$MSAI_present_any_sample_by_thresh %in% TRUE)
                  msai.seg.color <- ifelse(all.seg.MSAI.by.BAF.ltnpd.present.list[[chrom.seg.name]]$MSAI_present_any_sample_by_thresh %in% TRUE, msai.present.seg.color, non.msai.seg.color)

                } else {#end if ( length(all.seg.MSAI.by.BAF.ltnpd.present.list[[chrom.seg.name]]$MSAI_present_any_sample_by_thresh) > 0 ) {

                  msai.seg.color <- problem.seg.color

                }#end else

                #Plot whether MSAI detected
                segments(x0 = chrom.segs.to.plot.df[chrom.seg.name, 2]
                         ,y0 = 0
                         ,x1 = chrom.segs.to.plot.df[chrom.seg.name, 3]
                         ,y1 = 0
                         ,col = msai.seg.color
                         ,lwd = 6)

                #Plot if reference
                if( length(all.seg.MSAI.by.BAF.ltnpd.present.list[[chrom.seg.name]]$MSAI_detection_feasible) > 0 ) {

                  if (all.seg.MSAI.by.BAF.ltnpd.present.list[[chrom.seg.name]]$MSAI_detection_feasible %in% TRUE) {

                    seg.msai.df <- all.seg.MSAI.by.BAF.ltnpd.present.list[[chrom.seg.name]]$MSAI_by_sample_df

                    if (grepl("mixture",all.seg.MSAI.by.BAF.ltnpd.present.list[[chrom.seg.name]]$MSAI_detection_info)){

                      ref.seg.color <- mixture.phasing.color

                    } else if (seg.msai.df[seg.msai.df$sample_name == samp,]$sample_is_ref == TRUE) {

                      ref.seg.color <- reference.samp.color

                    } else {

                      ref.seg.color <- non.reference.samp.color

                    }

                  } else {

                    ref.seg.color <- no.detection.possible.color

                  }

                } else  {

                  ref.seg.color <- problem.seg.color

                }
                segments(x0 = chrom.segs.to.plot.df[chrom.seg.name, 2]
                         ,y0 = 0.05
                         ,x1 = chrom.segs.to.plot.df[chrom.seg.name, 3]
                         ,y1 = 0.05
                         ,col = ref.seg.color
                         ,lwd = 6)
                
              } 

            abline(v = max(as.numeric(chrom.segs.to.plot.df[chrom.segs.to.plot.df[, 1] == chr, 3])), lwd = 0.5)

                points(samp.BAF.df[, 2], samp.BAF.df$BAFsegmented, col = 'black', ylim = c(0 , 1), pch = 15, cex = 0.25)
                points(samp.BAF.df[, 2], samp.BAF.df$ph_median_segmented_BAF, col = '#78c679', ylim = c(0 , 1), pch = 15, cex = 0.5)

                points(samp.BAF.df[, 2], 1 - samp.BAF.df$BAFsegmented, col = 'black', ylim = c(0 , 1), pch = 15, cex = 0.25)
                points(samp.BAF.df[, 2], 1 - samp.BAF.df$ph_median_segmented_BAF, col = '#78c679', ylim = c(0 , 1), pch = 15, cex = 0.5)


          }#end for ( chr in sort(unique(as.numeric(chrom.segs.to.plot.df[, 1])))) {
          box()

      } else if (track.name == "cpnA" | track.name == "cpnB" ) {


        cpnX <- track.name
        print("Plotting ph cpn")
        #"B_MET004_MT-M1"
        Colours <- colorpanel(100, 'steelblue4', '#deebf7', 'steelblue4')

        #filter out the segments that are unphased
        phased.joint.segs.df  <- full.joint.segs.df[full.joint.segs.df$samp_used_to_phase != "too_few_SNPs",]
        #phased.joint.segs.df <- full.joint.segs.df
        # chr.start.end.df <- phased.joint.segs.df[,c("chr", "startpos", "endpos")]
        # chr.start.end.df <- apply(chr.start.end.df, 2, as.numeric)
        # #print("head(chr.start.end.df)")

        # chr.vec <- unique(full.joint.segs.df[, "chr"])
        # #print(head(chr.start.end.df))
        # chroms.length.vec <- getChromLengthsVec(chr.start.end.df)

        # names(chroms.length.vec) <- unique(chr.start.end.df[, 1])

        # #print("chroms.length.vec")
        # #print(chroms.length.vec)
        # chrom.segs.to.plot.df <- getCumulativeChromLengthsForPlottingDF(chr.start.end.df, chroms.length.vec, chr.vec)

        #print("head(chrom.segs.to.plot.df)")
        #print(head(chrom.segs.to.plot.df))
        min.x               <- 0
        min.y               <- -0.5
        max.y               <- 8.5
        max.x               <- as.numeric(max(chrom.segs.to.plot.df[, 3]))

        plot(1
             ,xlim = c(min.x, max.x)
             ,ylim = c(min.y, max.y)
             ,type = 'n'
             ,xaxt = 'n'
             ,yaxs = 'i'
             ,yaxt = 'n'
             ,xaxs = 'i'
             ,ylab = ""
             ,xlab = ""
             ,lwd = 2
             ,main = paste0(samp, " ", cpnX, " allele")
             ,bty = "n")
        box()

        samp.vals.list <- asas.ph.raw.cpn.lists[[samp]]
        samp.one.vals.list <- asas.ph.raw.cpn.lists[[first.samp]]
        seg.name.vec <- phased.joint.segs.df$seg_name #names(all.seg.CCF.vals.list)
        #print(seg.name.vec)
        for ( seg.name in seg.name.vec ) {
          #print(seg.name)
          #Will be NA if too few SNPs and can't be plotted, length of NA == 1, do this check to avoid warnings about checking only first element of list.
          if( length(all.seg.CCF.vals.list[[cpnX]][[seg.name]]) == 1 ) {

            next

          }

          chr.start.end.vec <- unlist(strsplit(seg.name, split = ":"))
          seg.samp.cpnX.vals.list <- samp.vals.list[[seg.name]]
          #print("str(seg.samp.cpnX.vals.list, max=1)")
          #print(str(seg.samp.cpnX.vals.list, max=1))
          seg.samp.one.cpnX.vals.list <- samp.one.vals.list[[seg.name]]
          #print("str(seg.samp.one.cpnX.vals.list, max=1)")
          #print(str(seg.samp.one.cpnX.vals.list, max=1))
          seg.ccf.cpnX.frac.vec   <- all.seg.CCF.vals.list[[cpnX]][[seg.name]][[samp]][[paste0("frac_", cpnX, "_low_vec")]]

          #print("print(str(seg.ccf.cpnX.frac.vec))")
          #print(str(seg.ccf.cpnX.frac.vec))

          seg.to.plot.df <- cbind(seg.samp.one.cpnX.vals.list$chr_vec,
                                  seg.samp.one.cpnX.vals.list$pos_vec,
                                  seg.samp.one.cpnX.vals.list$pos_vec,
                                  seg.samp.cpnX.vals.list[[paste0("ph_", cpnX, "_vec")]],
                                  seg.samp.cpnX.vals.list[[paste0("ph_", cpnX, "_vec")]],
                                  seg.samp.cpnX.vals.list[[paste0("ph_", cpnX, "_vec")]],
                                  seg.ccf.cpnX.frac.vec
                                  )

          if (is.na(seg.to.plot.df[1,1])) {

            next

          }#end if(is.na(seg.to.plot.df[1,1])) {
          
          #toPlot <- cbind(logRBAFSample1$chrom,logRBAFSample1$pos,logRBAFSample1$pos,ABSample1,ABSample2)
          seg.to.plot.df[,1:3] <- getCumulativeChromLengthsForPlottingDF(apply(seg.to.plot.df[,1:3], 2, as.numeric), chroms.length.vec, chr.vec)
          seg.to.plot.df <- data.frame(seg.to.plot.df, stringsAsFactors = FALSE)
          colnames(seg.to.plot.df) <- c('chr','pos','pos','Sample1nA','Sample1nA','Sample1nA','CCF')

          seg.to.plot.df$Sample1nA <- ifelse(as.numeric(seg.to.plot.df$Sample1nA) > max.y, (max.y - 0.25), as.numeric(seg.to.plot.df$Sample1nA))

          colourCode <- ifelse(round(as.numeric(seg.to.plot.df$CCF) * 100) > 100, 100, round(as.numeric(seg.to.plot.df$CCF) * 100))
          colourCode <- ifelse(colourCode == 0, 1, colourCode)


          points(x = as.numeric(seg.to.plot.df$pos), y = as.numeric(seg.to.plot.df$Sample1nA)
                 ,col = Colours[colourCode]
                 ,cex = .5
                 ,pch = 16)#ifelse(MSAIstatus==1,'#238b4550',ifelse(PValsPaired<0.001&totalPval<0.001,'#3182bd50','#fc927250')),cex=0.75,pch=16)
          par(lend = 1)
          segColour <- as.numeric(all.seg.CCF.vals.list[[cpnX]][[seg.name]][[samp]][[paste0("mean_frac_", cpnX, "_low")]])
          segColour <- Colours[(ifelse(segColour > 1, 1, segColour) * 100)]
          
          #So now very quickly we can check wehter that segment has been called significantly subclonal or not:
          #Here we want to implement the two options.. either use the CCF threshold or the CPN threshold



          sc.present <- full.joint.segs.df[paste0(samp, ":", seg.name), paste0("within_samp_sc_", cpnX, "_present_using_min_cpn_diff")]


          if (sc.present == TRUE) {

            ccf.bar.color <- "yellow"

          } else {

            ccf.bar.color <-"#525252" 

          }

          if (seg.name == "6:61888694:170806902") {

            STOP
          }
          segments(x0 = min(as.numeric(seg.to.plot.df$pos))
                   ,y0 = all.seg.CCF.vals.list[[cpnX]][[seg.name]][[samp]][[paste0("poss_sc_", cpnX, "_low")]]
                   ,x1 = max(as.numeric(seg.to.plot.df$pos))
                   ,y1 = all.seg.CCF.vals.list[[cpnX]][[seg.name]][[samp]][[paste0("poss_sc_", cpnX, "_low")]]
                   ,col = "#de2d26"
                   ,lwd = 3)
          
          segments(x0 = min(as.numeric(seg.to.plot.df$pos))
                   ,y0 = all.seg.CCF.vals.list[[cpnX]][[seg.name]][[samp]][[paste0("poss_sc_", cpnX, "_high")]]
                   ,x1 = max(as.numeric(seg.to.plot.df$pos))
                   ,y1 = all.seg.CCF.vals.list[[cpnX]][[seg.name]][[samp]][[paste0("poss_sc_", cpnX, "_high")]]
                   ,col = "#31a354"
                   ,lwd = 3)

          segments(x0 = min(as.numeric(seg.to.plot.df$pos))
                   ,y0 = mean(c(as.numeric(seg.to.plot.df$Sample1nA)))
                   ,x1 = max(as.numeric(seg.to.plot.df$pos)), y1 = mean(c(as.numeric(seg.to.plot.df$Sample1nA)))
                   ,col = ccf.bar.color
                   ,lwd = 3)
          
        }#end   for ( seg.name in seg.name.vec ) {

          axis(side = 2, las = 1)
          abline(h = seq((min.y + 0.5), (max.y - 0.5), by = 1), lty = 'dashed', lwd = 0.75, col = 'gray')
          
          mtext(side = 3
                ,at = getChromMidpointsVec(chrom.segs.to.plot.df[,1:3], unique(chr.start.end.df[,1]))
                ,text = sort(unique(chrom.segs.to.plot.df[,1]))
                ,cex = seq(0.6, 0.4, length.out = length(unique(chrom.segs.to.plot.df[,1])))
                ,line = -1
                #,lwd=0.5
          )
          
          for (i in sort(unique(as.numeric(chrom.segs.to.plot.df[,1]))))
          {
            abline(v = max(as.numeric(chrom.segs.to.plot.df[chrom.segs.to.plot.df[,1] == i, 3]))
                   ,lwd = 0.5)
            #,lty='dashed')
          }#en

          #box()

      } else if ( track.name == "mcs allele specific cpn" ) {

        print("Plotting MCS segs...")

        segs.df <- mcs.segs.df
        cpn.cols.to.plot <- c("int_ph_cpnA","int_ph_cpnB")
        color.vec = c("green", "black")
        print("samp")
        print(samp)
        first.segs.df <- segs.df[segs.df$sample == samp,]
        first.chr.start.end.df <- first.segs.df[,c("chr", "startpos", "endpos")]
        first.chr.start.end.df <- apply(first.chr.start.end.df, 2, as.numeric)
        #print("head(chr.start.end.df)")
        chr.vec <- unique(first.segs.df[, "chr"])
        #print(head(chr.start.end.df))
        first.chroms.length.vec <- getChromLengthsVec(first.chr.start.end.df)
        names(first.chroms.length.vec) <- unique(first.chr.start.end.df[, 1])

        #print("chroms.length.vec")
        #print(chroms.length.vec)
        first.chrom.segs.to.plot.df <- getCumulativeChromLengthsForPlottingDF(first.chr.start.end.df, first.chroms.length.vec, chr.vec)

        plot(1
             ,xlim = c(min.x, max.x)
             ,ylim = c(min.y, max.y)
             ,type = 'n'
             ,xaxt = 'n'
             ,yaxs = 'i'
             ,yaxt = 'n'
             ,xaxs = 'i'
             ,ylab = ""
             ,xlab = ""
             ,lwd = 2
             ,main = paste0(samp, " phased mcs segments allele specific copy number")
             ,bty = "n")
        box()

        samp.seg.df <- segs.df[segs.df$sample == samp,]
        samp.chr.start.end.df <- samp.seg.df[,c("chr", "startpos", "endpos")]
        samp.chr.start.end.df <- apply(samp.chr.start.end.df, 2, as.numeric)
        samp.chroms.length.vec <- getChromLengthsVec(samp.chr.start.end.df)
        names(samp.chroms.length.vec) <- unique(samp.chr.start.end.df[, 1])
        samp.chrom.segs.to.plot.df <- getCumulativeChromLengthsForPlottingDF(samp.chr.start.end.df, samp.chroms.length.vec, chr.vec)

        #print(head(samp.chrom.segs.to.plot.df))

        for ( k in 1:nrow(samp.seg.df) ) {

          for(j in 1:length(cpn.cols.to.plot)) {

            segments(x0 = samp.chrom.segs.to.plot.df[k, "startpos"]
                     ,y0 = samp.seg.df[k, cpn.cols.to.plot[j]] - ((j-1) * 0.05)
                     ,x1 = samp.chrom.segs.to.plot.df[k, "endpos"]
                     ,y1 = samp.seg.df[k, cpn.cols.to.plot[j]] - ((j-1) * 0.05)
                     ,col = color.vec[j]
                     ,lwd = 3
                     ,lend=2)

          }#end for(cpn.col in cpn.cols.to.plot) {

        }#end   for ( seg.name in seg.name.vec ) {

        axis(side = 2, las = 1)
        abline(h = seq((min.y + 0.5), (max.y - 0.5), by = 1), lty = 'dashed', lwd = 0.75, col = 'gray')
        
        mtext(side = 3
              ,at = getChromMidpointsVec(first.chrom.segs.to.plot.df[,1:3], unique(first.chr.start.end.df[,1]))
              ,text = sort(unique(first.chrom.segs.to.plot.df[,1]))
              ,cex = seq(0.6, 0.4, length.out = length(unique(first.chrom.segs.to.plot.df[,1])))
              ,line = -1
              #,lwd=0.5
              )
            
        for (i in sort(unique(as.numeric(first.chrom.segs.to.plot.df[,1]))))
        {
          abline(v = max(as.numeric(first.chrom.segs.to.plot.df[first.chrom.segs.to.plot.df[,1] == i, 3]))
                 ,lwd = 0.5)
          #,lty='dashed')
        }#en

          
      } else if (track.name == "final allele specific cpn") {

        print("Plotting final allele specific segs...")
        cpn.cols.to.plot  <- c("int_ph_cpnA","int_ph_cpnB")
        color.vec = c("green", "black")
        segs.df <- new.joined.segs.df

        first.segs.df <- segs.df[segs.df$sample == samp,]
        first.chr.start.end.df <- first.segs.df[,c("chr", "startpos", "endpos")]
        first.chr.start.end.df <- apply(first.chr.start.end.df, 2, as.numeric)
        #print("head(chr.start.end.df)")
        chr.vec <- unique(first.segs.df[, "chr"])
        #print(head(chr.start.end.df))
        first.chroms.length.vec <- getChromLengthsVec(first.chr.start.end.df)
        names(first.chroms.length.vec) <- unique(first.chr.start.end.df[, 1])

        #print("chroms.length.vec")
        #print(chroms.length.vec)
        first.chrom.segs.to.plot.df <- getCumulativeChromLengthsForPlottingDF(first.chr.start.end.df, first.chroms.length.vec, chr.vec)

        plot(1
             ,xlim = c(min.x, max.x)
             ,ylim = c(min.y, max.y)
             ,type = 'n'
             ,xaxt = 'n'
             ,yaxs = 'i'
             ,yaxt = 'n'
             ,xaxs = 'i'
             ,ylab = ""
             ,xlab = ""
             ,lwd = 2
             ,main = paste0(samp, " final MPhase joined allele specific cpn")
             ,bty = "n")
        box()

        samp.seg.df <- segs.df[segs.df$sample == samp,]
        samp.chr.start.end.df <- samp.seg.df[,c("chr", "startpos", "endpos")]
        samp.chr.start.end.df <- apply(samp.chr.start.end.df, 2, as.numeric)
        samp.chroms.length.vec <- getChromLengthsVec(samp.chr.start.end.df)
        names(samp.chroms.length.vec) <- unique(samp.chr.start.end.df[, 1])
        samp.chrom.segs.to.plot.df <- getCumulativeChromLengthsForPlottingDF(samp.chr.start.end.df, samp.chroms.length.vec, chr.vec)

        #print(head(samp.chrom.segs.to.plot.df))

        for ( k in 1:nrow(samp.seg.df) ) {

          for(j in 1:length(cpn.cols.to.plot)) {

            segments(x0 = samp.chrom.segs.to.plot.df[k, "startpos"]
                     ,y0 = samp.seg.df[k, cpn.cols.to.plot[j]] - ((j-1) * 0.05)
                     ,x1 = samp.chrom.segs.to.plot.df[k, "endpos"]
                     ,y1 = samp.seg.df[k, cpn.cols.to.plot[j]] - ((j-1) * 0.05)
                     ,col = color.vec[j]
                     ,lwd = 3
                     ,lend=2)

          }#end for(cpn.col in cpn.cols.to.plot) {

        }#end   for ( seg.name in seg.name.vec ) {

        axis(side = 2, las = 1)
        abline(h = seq((min.y + 0.5), (max.y - 0.5), by = 1), lty = 'dashed', lwd = 0.75, col = 'gray')
        
        mtext(side = 3
              ,at = getChromMidpointsVec(first.chrom.segs.to.plot.df[,1:3], unique(first.chr.start.end.df[,1]))
              ,text = sort(unique(first.chrom.segs.to.plot.df[,1]))
              ,cex = seq(0.6, 0.4, length.out = length(unique(first.chrom.segs.to.plot.df[,1])))
              ,line = -1
              #,lwd=0.5
              )
            
        for (i in sort(unique(as.numeric(first.chrom.segs.to.plot.df[,1]))))
        {
          abline(v = max(as.numeric(first.chrom.segs.to.plot.df[first.chrom.segs.to.plot.df[,1] == i, 3]))
                 ,lwd = 0.5)
          #,lty='dashed')
        }#en


      } else if (track.name == "original allele specific cpn") {

        segs.df <- orig.segs.df
        cpn.cols.to.plot <- c("int_ph_cpnA","int_ph_cpnB")
        color.vec = c("green", "black")

        first.segs.df <- segs.df[segs.df$sample == samp,]
        first.chr.start.end.df <- first.segs.df[,c("chr", "startpos", "endpos")]
        first.chr.start.end.df <- apply(first.chr.start.end.df, 2, as.numeric)
        #print("head(chr.start.end.df)")
        chr.vec <- unique(first.segs.df[, "chr"])
        #print(head(chr.start.end.df))
        first.chroms.length.vec <- getChromLengthsVec(first.chr.start.end.df)
        names(first.chroms.length.vec) <- unique(first.chr.start.end.df[, 1])

        #print("chroms.length.vec")
        #print(chroms.length.vec)
        first.chrom.segs.to.plot.df <- getCumulativeChromLengthsForPlottingDF(first.chr.start.end.df, first.chroms.length.vec, chr.vec)
        plot(1
             ,xlim = c(min.x, max.x)
             ,ylim = c(min.y, max.y)
             ,type = 'n'
             ,xaxt = 'n'
             ,yaxs = 'i'
             ,yaxt = 'n'
             ,xaxs = 'i'
             ,ylab = ""
             ,xlab = ""
             ,lwd = 2
             ,main = paste0(samp, " original unphased allele specific copy number")
             ,bty = "n")
        box()

        samp.seg.df <- segs.df[segs.df$sample == samp,]
        samp.chr.start.end.df <- samp.seg.df[,c("chr", "startpos", "endpos")]
        samp.chr.start.end.df <- apply(samp.chr.start.end.df, 2, as.numeric)
        samp.chroms.length.vec <- getChromLengthsVec(samp.chr.start.end.df)
        names(samp.chroms.length.vec) <- unique(samp.chr.start.end.df[, 1])
        samp.chrom.segs.to.plot.df <- getCumulativeChromLengthsForPlottingDF(samp.chr.start.end.df, samp.chroms.length.vec, chr.vec)

        #print(head(samp.chrom.segs.to.plot.df))

        for ( k in 1:nrow(samp.seg.df) ) {

          for(j in 1:length(cpn.cols.to.plot)) {

            segments(x0 = samp.chrom.segs.to.plot.df[k, "startpos"]
                     ,y0 = samp.seg.df[k, cpn.cols.to.plot[j]] - ((j-1) * 0.05)
                     ,x1 = samp.chrom.segs.to.plot.df[k, "endpos"]
                     ,y1 = samp.seg.df[k, cpn.cols.to.plot[j]] - ((j-1) * 0.05)
                     ,col = color.vec[j]
                     ,lwd = 3
                     ,lend=2)

          }#end for(cpn.col in cpn.cols.to.plot) {

        }#end   for ( seg.name in seg.name.vec ) {

        axis(side = 2, las = 1)
        abline(h = seq((min.y + 0.5), (max.y - 0.5), by = 1), lty = 'dashed', lwd = 0.75, col = 'gray')
        
        mtext(side = 3
              ,at = getChromMidpointsVec(first.chrom.segs.to.plot.df[,1:3], unique(first.chr.start.end.df[,1]))
              ,text = sort(unique(first.chrom.segs.to.plot.df[,1]))
              ,cex = seq(0.6, 0.4, length.out = length(unique(first.chrom.segs.to.plot.df[,1])))
              ,line = -1
              #,lwd=0.5
              )
            
        for (i in sort(unique(as.numeric(first.chrom.segs.to.plot.df[,1]))))
        {
          abline(v = max(as.numeric(first.chrom.segs.to.plot.df[first.chrom.segs.to.plot.df[,1] == i, 3]))
                 ,lwd = 0.5)
          #,lty='dashed')
        }#en



      }


  }#end tracname vec