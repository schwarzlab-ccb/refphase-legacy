

getSNPDataForTRACERxWGS <- function(snp.file.path) {

  pt.file.paths.df <- read.table(snp.file.path, stringsAsFactors = FALSE, sep = "\t", header = TRUE)
  sample.vec <- pt.file.paths.df$sample


  all.regions.ASCAT.list <- list()
  for ( samp in sample.vec ) {

    print(samp)
    samp.logr.baf.path <- pt.file.paths.df[pt.file.paths.df$sample == samp ,]$samp_logr_baf_file_path
    load(samp.logr.baf.path)#gives tumor.baf_logr
    samp.seg.path <- pt.file.paths.df[pt.file.paths.df$sample == samp ,]$samp_seg_file_path
    load(samp.seg.path)#gives seg

    ASCAT.summary               <- tumor.baf_logr[,c("chrom", "pos")]
    ASCAT.summary$BAF           <- tumor.baf_logr$BAF
    ASCAT.summary$logR          <- tumor.baf_logr[, samp]
    ASCAT.summary$BAFsegmented  <- NA
    ASCAT.summary$logRsegmented <- NA
    ASCAT.summary$sample <- samp

    all.regions.ASCAT.list[[samp]]  <- 


  }



  itx <- iter(sample.vec)
  #for ( samp in sample.vec ) {
  all.samp.SNP.df <- foreach(samp = itx, .combine = rbind) %dopar% {



    
    for ( j in 1:nrow(seg) ) {

        if(nrow(ASCAT.summary[ASCAT.summary$chrom == seg[j ,]$chr 
                            & ASCAT.summary$pos >= seg[j ,]$startpos
                            & ASCAT.summary$pos <= seg[j ,]$endpos,]) > 0) {

          if(seg[j,]$nMajor==seg[j,]$nMinor){
            
            ASCAT.summary[ASCAT.summary$chrom == seg[j ,]$chr
                        & ASCAT.summary$pos >= seg[j ,]$startpos 
                        & ASCAT.summary$pos <= seg[j ,]$endpos ,]$BAFsegmented <- 0.5
            
          } else {
            
            avg.seg.baf <- 0.5 - mean(abs(tumor.baf_logr[tumor.baf_logr$chrom == seg[j ,]$chr
                                                       & tumor.baf_logr$pos >= seg[j ,]$startpos
                                                       & tumor.baf_logr$pos <= seg[j ,]$endpos,]$BAF - 0.5))
            
            ASCAT.summary[ASCAT.summary$chrom == seg[j ,]$chr 
                        & ASCAT.summary$pos >= seg[j ,]$startpos 
                        & ASCAT.summary$pos <= seg[j ,]$endpos ,]$BAFsegmented <- avg.seg.baf
            
          }

          avg.seg.logr <- mean(tumor.baf_logr[tumor.baf_logr$chrom == seg[j ,]$chr 
                             & tumor.baf_logr$pos >= seg[j ,]$startpos 
                             & tumor.baf_logr$pos <= seg[j ,]$endpos, samp])
          
          ASCAT.summary[ASCAT.summary$chrom == seg[j ,]$chr 
                      & ASCAT.summary$pos >= seg[j ,]$startpos 
                      & ASCAT.summary$pos <= seg[j ,]$endpos,]$logRsegmented <- avg.seg.logr

        }#end if(nrow(ASCAT.summary[ASCAT.summary$chrom==seg[j,]$chr & ASCAT.summary$pos...

      }#end for ( j in 1:nrow(seg) ) {
      
      ASCAT.summary

  }# end  for ( samp in sample.vec ) {


    all.samp.SNP.list <- list()
    for ( samp in sample.vec ) {

      samp.w.segmented.info.SNP.df <- all.samp.SNP.df[all.samp.SNP.df$sample==samp,]
      samp.w.segmented.info.SNP.df <- samp.w.segmented.info.SNP.df[, c("chr", "pos", "BAF", "BAFsegmented", "logR", "logRsegmented")]
      all.samp.SNP.list[[samp]] <- samp.w.segmented.info.SNP.df

    }

    return(all.regions.ASCAT.list)
}#end getSNPDataForTRACERxWGS

