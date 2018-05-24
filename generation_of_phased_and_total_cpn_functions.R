





getAllSegAllSamplePhasedRawCPNList <- function(sample.vec, ploidy.est.vec, purity.est.vec, min.cons.segs.df, all.seg.phasing.info.list, asas.both.lr.baf.lists, gamma.val = 1, min.num.hz.snp.phasing = 5) {

    seg.name.vec <- paste(min.cons.segs.df[ , 1], min.cons.segs.df[ , 2], min.cons.segs.df[ , 3], sep=":")
    #Initialise all lists 
    all.samples.phased.cpn.list <- list()
    for ( samp in sample.vec ) {
      all.samples.phased.cpn.list[[samp]] <- list()
      for ( seg.name in seg.name.vec ) {
        all.samples.phased.cpn.list[[samp]][[seg.name]] <- list()
        all.samples.phased.cpn.list[[samp]][[seg.name]]$chr_vec <- NA 
        all.samples.phased.cpn.list[[samp]][[seg.name]]$pos_vec <- NA
        all.samples.phased.cpn.list[[samp]][[seg.name]]$logr_vec <- NA
        all.samples.phased.cpn.list[[samp]][[seg.name]]$ph_cpnA_vec <- NA
        all.samples.phased.cpn.list[[samp]][[seg.name]]$ph_cpnB_vec <- NA
        all.samples.phased.cpn.list[[samp]][[seg.name]]$mean_ph_cpnA <- NA
        all.samples.phased.cpn.list[[samp]][[seg.name]]$mean_ph_cpnB <- NA
        all.samples.phased.cpn.list[[samp]][[seg.name]]$un_ph_cpn_total
        all.samples.phased.cpn.list[[samp]][[seg.name]]$un_ph_cpn_total_mean <- NA
    }#end for ( seg.name in seg.name.vec ) {
  }#end for ( samp in sample.vec ) {
  #Now add all phased copy number usin gthe all.se.phasing.info.list
    for ( samp in sample.vec) {

      for ( seg.name in seg.name.vec ) {
        #We may have the problem here that some segments don't have any snps in them at all, check for this perhaps? see if fails
          samp.b.lr.baf.df <- asas.both.lr.baf.lists[[samp]][[seg.name]]
          if (nrow(samp.b.lr.baf.df) == 0) {
            next
          }#end   if (nrow(samp.b.lr.baf.df) == 0) {

          samp.seg.cpnA.cpnB.df <- calculateSNPsRawCopyNumberAB(logr.vec = samp.b.lr.baf.df$logR,
                                                            baf.vec = samp.b.lr.baf.df$BAF,
                                                            tumor.ploidy = ploidy.est.vec[samp],
                                                            tumor.purity = purity.est.vec[samp],
                                                            gamma.val = gamma.val)
   
        if (nrow(samp.b.lr.baf.df) < min.num.hz.snp.phasing) {

          all.samples.phased.cpn.list[[samp]][[seg.name]]$chr_vec <- samp.b.lr.baf.df$chrom
          all.samples.phased.cpn.list[[samp]][[seg.name]]$pos_vec <- samp.b.lr.baf.df$pos
          all.samples.phased.cpn.list[[samp]][[seg.name]]$logr_vec <- samp.b.lr.baf.df$logR
          all.samples.phased.cpn.list[[samp]][[seg.name]]$un_ph_cpn_total <- samp.seg.cpnA.cpnB.df[, 1] + samp.seg.cpnA.cpnB.df[, 2]
          all.samples.phased.cpn.list[[samp]][[seg.name]]$un_ph_cpn_total_mean <- mean(samp.seg.cpnA.cpnB.df[, 1] + samp.seg.cpnA.cpnB.df[, 2])


        } else {#end if (nrow(samp.b.lr.baf.df) > min.num.hz.snp.phasing) { now applying phasing to those that were phased.

          samp.seg.upper.band.cpn.vals.vec <- apply(samp.seg.cpnA.cpnB.df, 1, max)
          samp.seg.median.upper.band.cpn <- median(samp.seg.upper.band.cpn.vals.vec)
          samp.seg.cpnA.cpnB.is.upper.vec <- ifelse(samp.seg.cpnA.cpnB.df == samp.seg.upper.band.cpn.vals.vec, 1, 0)
          phased.samp.seg.cpnA.cpnB.df <- samp.seg.cpnA.cpnB.df
          #Get indices to phase cpnA
          if (length(all.seg.phasing.info.list[[seg.name]]$maj_is_cpnA_idx_vec) > 0) {

               is.upper.and.cpnA.in.ref.idx.vec <- all.seg.phasing.info.list[[seg.name]]$maj_is_cpnA_idx_vec
               phased.samp.seg.cpnA.cpnB.df[is.upper.and.cpnA.in.ref.idx.vec, 1] <- samp.seg.cpnA.cpnB.df[is.upper.and.cpnA.in.ref.idx.vec, 1]
               phased.samp.seg.cpnA.cpnB.df[is.upper.and.cpnA.in.ref.idx.vec, 2] <- samp.seg.cpnA.cpnB.df[is.upper.and.cpnA.in.ref.idx.vec, 2]

           
          }
          #Get indices to phase cpnB
          if (length(all.seg.phasing.info.list[[seg.name]]$maj_is_cpnB_idx_vec) > 0) {

               is.lower.and.cpnB.in.ref.idx.vec <- all.seg.phasing.info.list[[seg.name]]$maj_is_cpnB_idx_vec
               
               phased.samp.seg.cpnA.cpnB.df[is.lower.and.cpnB.in.ref.idx.vec, 1] <- samp.seg.cpnA.cpnB.df[is.lower.and.cpnB.in.ref.idx.vec, 2]
               phased.samp.seg.cpnA.cpnB.df[is.lower.and.cpnB.in.ref.idx.vec, 2] <- samp.seg.cpnA.cpnB.df[is.lower.and.cpnB.in.ref.idx.vec, 1]
               
          } 
          #Store this information in the list:
          all.samples.phased.cpn.list[[samp]][[seg.name]]$chr_vec <- samp.b.lr.baf.df$chrom
          all.samples.phased.cpn.list[[samp]][[seg.name]]$pos_vec <- samp.b.lr.baf.df$pos
          all.samples.phased.cpn.list[[samp]][[seg.name]]$logr_vec <- samp.b.lr.baf.df$logR
          all.samples.phased.cpn.list[[samp]][[seg.name]]$ph_cpnA_vec <- as.numeric(phased.samp.seg.cpnA.cpnB.df[,1])
          all.samples.phased.cpn.list[[samp]][[seg.name]]$ph_cpnB_vec <- as.numeric(phased.samp.seg.cpnA.cpnB.df[,2])
          all.samples.phased.cpn.list[[samp]][[seg.name]]$mean_ph_cpnA <- mean(as.numeric(phased.samp.seg.cpnA.cpnB.df[,1]))
          all.samples.phased.cpn.list[[samp]][[seg.name]]$mean_ph_cpnB <- mean(as.numeric(phased.samp.seg.cpnA.cpnB.df[,2]))

        }#end if(nrow(samp.b.lr.baf.df) > min.num.hz.snp.phasing)
      }#end for ( seg.name in chr.start.end.string.vec ) {
    }#end for   for ( samp in sample.vec) {
    return(all.samples.phased.cpn.list)
}#end getAllSegAllSamplePhasedRawCPNList







getAllSampLogrAlonePositionsLogrAndTotalCPNDFsList <- function(sample.vec, asas.SNPS.lr.and.or.baf.lists, ploidy.est.vec, purity.est.vec, gamma.val) {

  all.samp.logr.df.list <- list()

  for ( samp in sample.vec ) {

    samp.logr.baf.df <- do.call(rbind, asas.SNPS.lr.and.or.baf.lists[[samp]])

    if ( nrow(samp.logr.baf.df[is.na(samp.logr.baf.df$BAF) ,]) > 0 ) {

      samp.logr.df <- samp.logr.baf.df[ ! is.na(samp.logr.baf.df$BAF) ,]

    } else if ( nrow(samp.logr.baf.df[samp.logr.baf.df$BAF == 2 ,]) > 0) {

      samp.logr.df <- samp.logr.baf.df[ samp.logr.baf.df$BAF == 2 ,]

    } else {

      print("No logr only SNPs, 0 row dataframe will be passed to functions.")
      samp.logr.df <- samp.logr.baf.df[ samp.logr.baf.df$BAF == 2 ,]

    }


    if (nrow(samp.logr.baf.df) > 0 ) {

      mock.cpnA.cpnB.logr.df <- calculateSNPsRawCopyNumberAB(logr.vec = samp.logr.df$logR,
                                                            baf.vec = rep(0.5, nrow(samp.logr.df)),
                                                            tumor.ploidy = ploidy.est.vec[samp],
                                                            tumor.purity = purity.est.vec[samp],
                                                            gamma.val = gamma.val)
      samp.logr.df$tot_raw_un_ph_cpn <- unlist(rowSums(mock.cpnA.cpnB.logr.df))

      #samp.logr.df$tot_raw_un_ph_cpn <- rowSums(mock.cpnA.cpnB.logr.df)

    }

    all.samp.logr.df.list[[samp]] <- samp.logr.df
    
  }


  return(all.samp.logr.df.list)

}#end getAllSampLogrAlonePositionsLogrAndTotalCPNDFsList





#This is used to generate CCF values so needs to be added to the sample.to.phase.seg.list
createAllSegsPhasedBAFDataFrame <- function(asas.lr.baf.lists, sample.to.phase.seg.list, min.num.hz.snp.phasing = 5) {

  sample.vec <- names(asas.lr.baf.lists)
 
  #Now we need to make sure that the asas.lr.baf.lists and sample.to.phase.seg.list 
  if (length(names(sample.to.phase.seg.list)) != length(asas.lr.baf.lists[[1]])) {

    print("Missing or removed segments:")
    print(names(sample.to.phase.seg.list)[ ! names(asas.lr.baf.lists[[1]]) %in% names(sample.to.phase.seg.list)])
    print("Continuing creation of phased BAF dataframe")

  }

  seg.name.vec <- names(sample.to.phase.seg.list)
  all.samp.ph.BAF.df.list <- list()

  for ( seg.name in seg.name.vec ) {

    for ( samp in sample.vec ) {
      #Todo, make sure all segments are referred to by name, define this earlier on. Makes no sense to sometimes be using
      #indices and other times be using the name format, causes confusion and bugs.
      samp.BAF.df <- asas.lr.baf.lists[[samp]][[seg.name]]

      if (nrow(samp.BAF.df) < min.num.hz.snp.phasing) {

        print(seg.name)
        print("Has fewer than 5 heterozygous SNPs, not updating corresponding BAFs in table")
        next;

      }#end if(nrow(samp.BAF.df)<5) {
      
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

      if (nrow(samp.BAF.df[is.na(samp.BAF.df$maj_is_cpnA_logic_vec), ,drop=FALSE]) > 0) {

        print("samp.BAF.df")
        print(samp.BAF.df)
        print("nrow(samp.BAF.df)")
        print(nrow(samp.BAF.df))
        stop("Check error in sample BAF df generation")
        
      }#end if (nrow(samp.BAF.df[is.na(samp.BAF.df$maj_is_cpnA_idx_vec), ,drop=FALSE])>0)){

      samp.BAF.df$ph_median_segmented_BAF <- min(median(samp.BAF.df[samp.BAF.df$maj_is_cpnB_logic_vec == 1 ,]$BAF),
                                           median(samp.BAF.df[samp.BAF.df$maj_is_cpnB_logic_vec == 0 ,]$BAF))

      all.samp.ph.BAF.df.list[[samp]] <- rbind(all.samp.ph.BAF.df.list[[samp]], samp.BAF.df)

    }#end for ( samp in sample.vec ) {

  }#end for ( i in 1:length(chr.start.end.string.vec) ) {

  return(all.samp.ph.BAF.df.list) 

}#end createSegPhasedBAFTable <- function(....)
