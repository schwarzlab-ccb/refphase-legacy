##
##Params:
##
##Returns:
##
funPlotAllRegions_withJoining <- function(all.pw.cn.het.comp.lists, full.joint.segs.df, sample.vec,final.comb.seg.df, all.seg.CCF.vals.list, cpnX){

  Colours <- colorpanel(100,'steelblue4','#deebf7','steelblue4')

  #First get the consensus segmentation, probably want to filter out the smallest segments:
  #filter out the segments that are unphased
  phased.joint.segs.df  <- full.joint.segs.df[full.joint.segs.df$samp_used_to_phase != "too_few_SNPs",]
  #phased.joint.segs.df <- full.joint.segs.df
  chr.start.end.df <- phased.joint.segs.df[,c("chr","startpos","endpos")]
  chr.start.end.df <- apply(chr.start.end.df,2,as.numeric)
  chroms.length.vec <- getChromLengthsVec(chr.start.end.df)
  names(chroms.length.vec) <- 1:22
  chrom.segs.to.plot.df <- getCumulativeChromLengthsForPlottingDF(chr.start.end.df, chroms.length.vec)

  min.x               <- 0
  min.y               <- -0.5
  max.y               <- 8.5
  max.x               <- as.numeric(max(chrom.segs.to.plot.df[,3]))

  num.samples  <- length(sample.vec)
  first.samp <- sample.vec[1] 
  last.samp <- sample.vec[length(sample.vec)] 
  #regionJoinedSeg <- funJoinAdjacentSegments(copyITHRDS$regionComparisonList) is the equivalent
  final.comb.seg.df

  par(mfrow = c(num.samples, 1))

  for ( samp in sample.vec ) {

    #Set out border size for sub plot
    if (samp == first.samp) {

      par(mar = c(2,3,4,3))

    } else if (samp == last.samp) {

      par(mar = c(4,3,2,3))

    } else {

      par(mar = c(3,3,3,3))

    }
    plot(1
         ,xlim=c(min.x,max.x)
         ,ylim=c(min.y,max.y)
         ,type='n'
         ,xaxt='n'
         ,yaxs='i'
         ,yaxt='n'
         ,xaxs='i'
         ,ylab=""
         ,xlab=""
         ,lwd=2
         ,main=paste(samp," ",cpnX," allele",sep="")
         ,bty="n")
    box()

    #"B_MET004_MT-M1"
    samp.vals.list <- all.pw.cn.het.comp.lists[[samp]]
    samp.one.vals.list <- all.pw.cn.het.comp.lists[[first.samp]]

    seg.name.vec <- phased.joint.segs.df$seg_name #names(all.seg.CCF.vals.list)
    for ( seg.name in seg.name.vec ) {

      if( is.na(all.seg.CCF.vals.list[[seg.name]]) ) {

        next

      }

      chr.start.end.vec <- unlist(strsplit(seg.name, split=":"))
      seg.samp.cpnX.vals.list <- samp.vals.list[[seg.name]]
      seg.samp.one.cpnX.vals.list <- samp.one.vals.list[[seg.name]]
      seg.ccf.cpnX.frac.vec   <- all.seg.CCF.vals.list[[seg.name]][[samp]][[paste0("frac_",cpnX,"_low_vec")]]

      
      seg.to.plot.df <- cbind(seg.samp.one.cpnX.vals.list$chr_vec,
                              seg.samp.one.cpnX.vals.list$pos_vec,
                              seg.samp.one.cpnX.vals.list$pos_vec,
                              seg.samp.cpnX.vals.list[[paste0("ph_",cpnX,"_vec")]],
                              seg.samp.cpnX.vals.list[[paste0("ph_",cpnX,"_vec")]],
                              seg.samp.cpnX.vals.list[[paste0("ph_",cpnX,"_vec")]],
                              seg.ccf.cpnX.frac.vec
                              )

      if(is.na(seg.to.plot.df[1,1])) {

        next

      }#end if(is.na(seg.to.plot.df[1,1])) {
      
      #toPlot <- cbind(logRBAFSample1$chrom,logRBAFSample1$pos,logRBAFSample1$pos,ABSample1,ABSample2)
      seg.to.plot.df[,1:3] <- getCumulativeChromLengthsForPlottingDF(apply(seg.to.plot.df[,1:3],2,as.numeric), chroms.length.vec)
      seg.to.plot.df <- data.frame(seg.to.plot.df,stringsAsFactors = FALSE)
      colnames(seg.to.plot.df) <- c('chr','pos','pos','Sample1nA','Sample1nA','Sample1nA','CCF')

      seg.to.plot.df$Sample1nA <- ifelse(as.numeric(seg.to.plot.df$Sample1nA)>max.y,(max.y-0.25),as.numeric(seg.to.plot.df$Sample1nA))

      colourCode <- ifelse(round(as.numeric(seg.to.plot.df$CCF)*100)>100,100,round(as.numeric(seg.to.plot.df$CCF)*100))
      colourCode <- ifelse(colourCode==0,1,colourCode)


      points(x = as.numeric(seg.to.plot.df$pos),y=as.numeric(seg.to.plot.df$Sample1nA)
             ,col=Colours[colourCode]
             ,cex=.5
             ,pch=16)#ifelse(MSAIstatus==1,'#238b4550',ifelse(PValsPaired<0.001&totalPval<0.001,'#3182bd50','#fc927250')),cex=0.75,pch=16)
      par(lend=1)
      segColour <- as.numeric(all.seg.CCF.vals.list[[seg.name]][[samp]][[paste0("mean_frac_",cpnX,"_low")]])
      segColour <- Colours[(ifelse(segColour > 1, 1, segColour)*100)]
      
      #So now very quickly we can check wehter that segment has been called significantly subclonal or not:

      sc.present <- full.joint.segs.df[paste0(samp,":",seg.name), paste0("sc_",cpnX,"_present")]

      if (sc.present == TRUE) {

        ccf.bar.color <- "yellow"

      } else {

        ccf.bar.color <-"#525252" 

      }

      print(all.seg.CCF.vals.list[[seg.name]][[paste0("poss_sc_",cpnX,"_low")]])
      print(all.seg.CCF.vals.list[[seg.name]][[paste0("poss_sc_",cpnX,"_high")]])

      segments(x0=min(as.numeric(seg.to.plot.df$pos))
               ,y0=all.seg.CCF.vals.list[[seg.name]][[samp]][[paste0("poss_sc_",cpnX,"_low")]]
               ,x1=max(as.numeric(seg.to.plot.df$pos))
               ,y1=all.seg.CCF.vals.list[[seg.name]][[samp]][[paste0("poss_sc_",cpnX,"_low")]]
               ,col="#de2d26"
               ,lwd=3)
      
      segments(x0=min(as.numeric(seg.to.plot.df$pos))
               ,y0=all.seg.CCF.vals.list[[seg.name]][[samp]][[paste0("poss_sc_",cpnX,"_high")]]
               ,x1=max(as.numeric(seg.to.plot.df$pos))
               ,y1=all.seg.CCF.vals.list[[seg.name]][[samp]][[paste0("poss_sc_",cpnX,"_high")]]
               ,col="#31a354"
               ,lwd=3)

      segments(x0=min(as.numeric(seg.to.plot.df$pos))
               ,y0=mean(c(as.numeric(seg.to.plot.df$Sample1nA)))
               ,x1=max(as.numeric(seg.to.plot.df$pos)),y1=mean(c(as.numeric(seg.to.plot.df$Sample1nA)))
               ,col=ccf.bar.color
               ,lwd=3)
      


    }#end   for ( seg.name in seg.name.vec ) {

      axis(side = 2, las = 1)
      abline(h=seq((min.y+0.5),(max.y-0.5),by=1),lty='dashed',lwd=0.75,col='gray')
      
      mtext(side=3
            ,at=fun.chrom.mean(chrom.segs.to.plot.df[,1:3])
            ,text=sort(unique(chrom.segs.to.plot.df[,1]))
            ,cex=seq(0.6,0.4,length.out=length(unique(chrom.segs.to.plot.df[,1])))
            ,line=-1
            #,lwd=0.5
      )
      
      for (i in sort(unique(as.numeric(chrom.segs.to.plot.df[,1]))))
      {
        abline(v=max(as.numeric(chrom.segs.to.plot.df[chrom.segs.to.plot.df[,1]==i,3]))
               ,lwd=0.5)
        #,lty='dashed')
      }#en

  }#end for ( samp in sample.vec ) {

}



#These small functions should bve quick to rewrite:
#was fun.chrom.length
getChromLengthsVec <- function(seg){
  # This function returns the length of chromosome 1,...,22
  # seg: first three columns of the minimum conistent region matrix
  chrom.length <- c()
  for(i in 1:22){
    sub.chrom <- subset(seg, seg[, 1] == i)
    chrom.length <- c(chrom.length, sub.chrom[nrow(sub.chrom), 3])
    
  }
  return(chrom.length)
  
}

#Was fun.add.chrom
getCumulativeChromLengthsForPlottingDF <-function(seg, chrom.length){
  # This function adds the length of chromosome 1,...,(i - 1) to chromsomes 2,...,22
  # seg: first three columns of the minimum conistent region matrix
  # chrom.length: output of fun.chrom.length(...)
  if(1 %in% seg[, 1])
    seg.tmp <-  subset(seg, seg[, 1] == 1)
  else
    seg.tmp <- matrix(NA, nr = 0, nc = ncol(seg))  
  for(i in 2:22){
    if(!(i %in% seg[, 1]))
      next()
    sub.chrom <- subset(seg, seg[, 1] == i)
    sub.chrom[, 2] <- sub.chrom[, 2] + sum(chrom.length[1:(i - 1)])
    sub.chrom[, 3] <- sub.chrom[, 3] + sum(chrom.length[1:(i - 1)])
    seg.tmp <- rbind(seg.tmp, sub.chrom)
  }
  return(seg.tmp)
  
}

fun.chrom.mean<-function(seg){
  # This function finds the middle of the positions for each chromosome.
  # seg: Matrix with three columns, containing the positions.
  chrom.mean <- c()
  for(i in 1:22){
    sub.chrom <- subset(seg, seg[, 1] == i)
    chrom.mean <- c(chrom.mean, (min(sub.chrom[, 2]) + max(sub.chrom[, 3]))/2)
    
  }
  return(chrom.mean)
}

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    