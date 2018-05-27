
makeInputPathFileForTRACERxExomeSNPArray <- function(pt.id, sample.vec, release.dir, array.type, save.output.path) {

  pt.dir <- paste(release.dir, pt.id, "/", array.type, "/ASCAT/", sep="")
  #sample.vec <- grep("FAILED",grep(".txt$|.complete$|*EXCLUDE*|*ManualQC*",list.files(pt.dir), value=TRUE, invert=TRUE), invert=TRUE, value=TRUE)

  pt.file.paths.df <- data.frame(patient_id = rep(pt.id, length(sample.vec)),
                sample = sample.vec,
                samp_ascat_bc_path= rep("None", length(sample.vec)),
                stringsAsFactors = FALSE)



  for ( samp in sample.vec ) {
    #We need to check for the existence of directories that look like this: B_MET005_MT-M1_ManualQC
    samp.poss.qc.file.path <- paste(pt.dir, samp, "_ManualQC", "/", samp, ".ascat.bc.RData", sep="")
    samp.regular.file.path <- paste(pt.dir, samp, "/", samp, ".ascat.bc.RData", sep="")

    if (file.exists(samp.poss.qc.file.path) ) {

      pt.file.paths.df[pt.file.paths.df$sample == samp, ]$samp_ascat_bc_path <- samp.poss.qc.file.path


    } else if (file.exists(samp.regular.file.path)) {

      pt.file.paths.df[pt.file.paths.df$sample == samp, ]$samp_ascat_bc_path <- samp.regular.file.path
      
    } else {

      stop("No sample ascat bc object present in either a QC or regular directory")

    }

  }

  write.table(pt.file.paths.df, file = save.output.path , row.names = FALSE, sep = "\t", quote = FALSE)

}

getSNPDataForTRACERxExomeSNPArray <- function(pt.file.paths.df) {

  sample.vec <- pt.file.paths.df$sample
  all.samples.ASCAT.list <- list()
  
  for ( samp in sample.vec ) {

    samp.ASCAT.obj <- pt.file.paths.df[pt.file.paths.df$sample == samp, ]$samp_ascat_bc_path
    samp.ASCAT.obj <- load(samp.ASCAT.obj)
    samp.ASCAT.obj <- get(samp.ASCAT.obj)

    samp.ASCAT.summary               <- samp.ASCAT.obj$SNPpos
    samp.ASCAT.summary$BAF           <- NA
    samp.ASCAT.summary$BAFsegmented  <- NA
    samp.ASCAT.summary$logR          <- NA
    samp.ASCAT.summary$logRsegmented <- NA

    samp.ASCAT.summary[rownames(samp.ASCAT.obj$Tumor_LogR_segmented),]$logRsegmented           <- samp.ASCAT.obj$Tumor_LogR_segmented[,1]
    samp.ASCAT.summary[rownames(samp.ASCAT.obj$Tumor_BAF),]$BAF                                <- samp.ASCAT.obj$Tumor_BAF[[1]]
    samp.ASCAT.summary[rownames(samp.ASCAT.obj$Tumor_BAF_segmented[[1]]),]$BAFsegmented        <- samp.ASCAT.obj$Tumor_BAF_segmented[[1]][,1]
    samp.ASCAT.summary[rownames(samp.ASCAT.obj$Tumor_LogR),]$logR                              <- samp.ASCAT.obj$Tumor_LogR[[1]]


    all.samples.ASCAT.list[[samp]] <- samp.ASCAT.summary

  }#end getSNPDataForTRACERxExomeSNPArray

  return(all.samples.ASCAT.list)

}#end getSNPDataForTRACERxExomeSNPArray <- function(pt.file.paths.df) {



makeInputPathFileForTRACERxWGS <- function(pt.id, release.dir, array.type, save.output.path) {

  #Do we not care about the ASCAT_QC directories, we can come back to these a bit later
  pt.dir <- paste(release.dir, pt.id, "/", array.type, "/ASCAT/", sep = "")
  sample.vec <- grep("FAILED",grep(".txt$|.complete$|*EXCLUDE*|*ManualQC*", list.files(pt.dir), value=TRUE, invert=TRUE), invert=TRUE, value=TRUE)
  pt.file.paths.df <- data.frame(patient_id = rep(pt.id, length(sample.vec)),
                sample = sample.vec,
                samp_logr_baf_file_path = rep("None", length(sample.vec)),
                samp_seg_file_path = rep("None", length(sample.vec)),
                stringsAsFactors = FALSE)


  for ( samp in sample.vec ) {

    samp.logr.baf.path <- paste(pt.dir, samp, "/", samp, ".tumor.baf_logr.RData", sep="")

    pt.file.paths.df[pt.file.paths.df$sample==samp, ]$samp_logr_baf_file_path <- samp.logr.baf.path

    samp.seg.path <- paste(pt.dir, samp, "/", samp, ".ascat.seg.RData", sep="")

    pt.file.paths.df[pt.file.paths.df$sample==samp, ]$samp_seg_file_path <- samp.seg.path

  }#end for ( samp in sample.vec ) {

  write.table(pt.file.paths.df, file = save.output.path ,row.names=FALSE, sep="\t", quote=FALSE)

}#end makeInputPathFileForTRACERxWGS <- function(pt.id, release.dir, array.type, save.output.path) {




getSNPDataForTRACERxWGS <- function(snp.file.path) {

  pt.file.paths.df <- read.table(snp.file.path, stringsAsFactors = FALSE, sep = "\t", header = TRUE)
  sample.vec <- pt.file.paths.df$sample
  all.regions.ASCAT.list <- list()

  for ( samp in sample.vec ) {

    samp.logr.baf.path <- pt.file.paths.df[pt.file.paths.df$sample == samp ,]$samp_logr_baf_file_path
    load(samp.logr.baf.path)#gives tumor.baf_logr
    samp.seg.path <- pt.file.paths.df[pt.file.paths.df$sample == samp ,]$samp_seg_file_path
    load(samp.seg.path)#gives seg

    ASCAT.summary               <- tumor.baf_logr[,c("chrom", "pos")]
    ASCAT.summary$BAF           <- tumor.baf_logr$BAF
    ASCAT.summary$logR          <- tumor.baf_logr[, samp]
    ASCAT.summary$BAFsegmented  <- NA
    ASCAT.summary$logRsegmented <- NA

    

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
      
      all.regions.ASCAT.list[[samp]] <- ASCAT.summary

  }# end  for ( samp in sample.vec ) {

    return(all.regions.ASCAT.list)
}#end getSNPDataForTRACERxWGS




##
##Params:
##
##Returns:
##
errorCheckAndAddPloidyAdjRawCPNToSegDF <- function(pt.id, orig.cpn.segs.df){

  #TODO add checks that that:
  #all samples have the same number of chromosomes
  #There is more than one sample
  #minor is never larger than major cn
  #end value seg end postion is always greater than start
  #chromosomes and segments are sorted correctly
  #No NA values
  #Purity no greater than 1, flag if it is.
  #Check if there are snps in all segements

  orig.cpn.segs.df <- orig.cpn.segs.df[orig.cpn.segs.df$patient==pt.id,,drop=FALSE]

  samples <- unique(orig.cpn.segs.df$sample)

  orig.cpn.segs.df$tot_raw_cpn <- as.numeric(orig.cpn.segs.df$nAraw) + as.numeric(orig.cpn.segs.df$nBraw)
  orig.cpn.segs.df$min_raw_cpn <- apply(cbind(orig.cpn.segs.df$nAraw, orig.cpn.segs.df$nBraw), 1, min)
  orig.cpn.segs.df$max_raw_cpn <- apply(cbind(orig.cpn.segs.df$nAraw, orig.cpn.segs.df$nBraw), 1, max)
  orig.cpn.segs.df$cpn_to_use <- orig.cpn.segs.df$tot_raw_cpn

  orig.cpn.segs.df$ploidy_adj_tot_cpn <- log2(as.numeric(orig.cpn.segs.df$tot_raw_cpn) / as.numeric(orig.cpn.segs.df$Ploidy))

  return(orig.cpn.segs.df)
}



##
##Params:
##
##Returns:
##
getAllSNPsWithLogRatioDataTRACERxFormatASCATBC <- function(pt.id, array.type="exome", release.dir){

  pt.dir <- paste(release.dir, pt.id, "/", array.type, "/ASCAT/", sep="")
  regions <- grep("FAILED",grep(".txt$|.complete$|*EXCLUDE*|*ManualQC*",list.files(pt.dir), value=TRUE, invert=TRUE), invert=TRUE, value=TRUE)
  all.regions.ASCAT.list <- list()

  for ( region in regions ) {

      reg.ASCAT.obj <- load(paste(pt.dir,region,"/",region,".ascat.bc.RData",sep=""))
      reg.ASCAT.obj <- get(reg.ASCAT.obj)

      #let's put together a table with what we need
      reg.ASCAT.summary               <- reg.ASCAT.obj$SNPpos
      reg.ASCAT.summary$BAF           <- NA
      reg.ASCAT.summary$BAFsegmented  <- NA
      reg.ASCAT.summary$logR          <- NA
      reg.ASCAT.summary$logRsegmented <- NA
      
      reg.ASCAT.summary[rownames(reg.ASCAT.obj$Tumor_LogR_segmented),]$logRsegmented           <- reg.ASCAT.obj$Tumor_LogR_segmented[,1]
      reg.ASCAT.summary[rownames(reg.ASCAT.obj$Tumor_BAF) ,]$BAF                                <- reg.ASCAT.obj$Tumor_BAF[[1]]
      reg.ASCAT.summary[rownames(reg.ASCAT.obj$Tumor_BAF_segmented[[1]]),]$BAFsegmented        <- reg.ASCAT.obj$Tumor_BAF_segmented[[1]][,1]
      reg.ASCAT.summary[rownames(reg.ASCAT.obj$Tumor_LogR) ,]$logR                              <- reg.ASCAT.obj$Tumor_LogR[[1]]

      all.regions.ASCAT.list[[region]] <- reg.ASCAT.summary

  }#end for region in regions

  return(all.regions.ASCAT.list)

}


##
##Params:
##
##Returns:
##
getSNPsWithBAFandLogRatioDataTRACERxFormatASCATBC <- function(pt.id, array.type="exome", release.dir){

  pt.dir <- paste(release.dir, pt.id, "/", array.type, "/ASCAT/", sep="")
  regions <- grep("FAILED",grep(".txt$|.complete$|*EXCLUDE*|*ManualQC*",list.files(pt.dir), value=TRUE, invert=TRUE), invert=TRUE, value=TRUE)
  all.regions.ASCAT.list <- list()

  for ( region in regions ) {

      reg.ASCAT.obj <- load(paste(pt.dir,region,"/",region,".ascat.bc.RData",sep=""))
      reg.ASCAT.obj <- get(reg.ASCAT.obj)

      #let's put together a table with what we need
      reg.ASCAT.summary               <- reg.ASCAT.obj$SNPpos
      reg.ASCAT.summary$BAF           <- NA
      reg.ASCAT.summary$BAFsegmented  <- NA
      reg.ASCAT.summary$logR          <- NA
      reg.ASCAT.summary$logRsegmented <- NA
      
      reg.ASCAT.summary[rownames(reg.ASCAT.obj$Tumor_LogR_segmented),]$logRsegmented           <- reg.ASCAT.obj$Tumor_LogR_segmented[,1]
      reg.ASCAT.summary[rownames(reg.ASCAT.obj$Tumor_BAF),]$BAF                                <- reg.ASCAT.obj$Tumor_BAF[[1]]
      reg.ASCAT.summary[rownames(reg.ASCAT.obj$Tumor_BAF_segmented[[1]]),]$BAFsegmented        <- reg.ASCAT.obj$Tumor_BAF_segmented[[1]][,1]
      reg.ASCAT.summary[rownames(reg.ASCAT.obj$Tumor_LogR),]$logR                              <- reg.ASCAT.obj$Tumor_LogR[[1]]

      all.regions.ASCAT.list[[region]] <- reg.ASCAT.summary

  }#end for region in regions

  return(all.regions.ASCAT.list)

}

##
##Params:
##
##Returns:
##
getSNPDataFromTRACERxFormatBAFLogRData <- function(pt.id, release.dir, array.type = "WGS"){

  pt.dir <- paste(release.dir, pt.id, "/", array.type, "/ASCAT/", sep="")
  regions <- grep("FAILED",grep(".txt$|.complete$|*EXCLUDE*|*ManualQC*", list.files(pt.dir), value = TRUE, invert = TRUE), invert = TRUE, value = TRUE)

  for ( i in 1:length(regions) ) {

    reg.logr.baf.path <- paste(DirToConsider, regions[i], "/", regions[i], ".tumor.baf_logr.RData", sep="")
    load(reg.logr.baf.path)#gives tumor.baf_logr obj

    reg.seg.path <- paste(DirToConsider, regions[i], "/", regions[i], ".ascat.seg.RData", sep="")
    load(reg.seg.path)#gives reg.seg.path

    ASCAT.summary               <- tumor.baf_logr[, c("chrom", "pos")]
    ASCAT.summary$BAF           <- tumor.baf_logr$BAF
    ASCAT.summary$logR          <- tumor.baf_logr[, regions[i]]
    ASCAT.summary$BAFsegmented  <- NA
    ASCAT.summary$logRsegmented <- NA

    all.regions.ASCAT.list <- list()

    #Update the fake ASCAT.summaary BAFsegmetned and logRsegmented
    for ( j in 1:nrow(seg) ) {

      if(nrow(ASCAT.summary[ASCAT.summary$chrom == seg[j ,]$chr
                            & ASCAT.summary$pos >= seg[j ,]$startpos
                            & ASCAT.summary$pos <= seg[j ,]$endpos ,]) > 0) {

        if(seg[j,]$nMajor == seg[j,]$nMinor) {
          
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
                                          & tumor.baf_logr$pos <= seg[j ,]$endpos, regions[i]])

        ASCAT.summary[ASCAT.summary$chrom == seg[j ,]$chr
                    & ASCAT.summary$pos >= seg[j ,]$startpos
                    & ASCAT.summary$pos <= seg[j ,]$endpos,]$logRsegmented <- avg.seg.logr

      }#end if(nrow(ASCAT.summary[ASCAT.summary$chrom==seg[j,]$chr & ASCAT.summary$pos...

    }#end for ( j in 1:nrow(seg) ) {
    
    all.regions.ASCAT.list[[regions[i]]] <- ASCAT.summary

  }#end for ( i in 1:length(regions) )

  return(all.regions.ASCAT.list)
}

##
##Params:
##
##Returns:
##
cleanAndErrorCheckMinConsSegs <- function(min.cons.seg.df) {

  #TODO report any rows removed due to being NA
  #Any values updated to from Inf or -Inf
  #Adjust if very skewed values, perhaps add to function?

  min.cons.segs.df[min.cons.segs.df == '-Inf'] <- -3
  min.cons.segs.df[min.cons.segs.df == 'Inf'] <- 3
  #remove any NA rows 
  min.cons.seg.df <- min.cons.seg.df[!is.na(rowSums(min.cons.seg.df[,4:ncol(min.cons.seg.df)])),,drop=FALSE]

  return(min.cons.seg.df)
} #end cleanAndErrorCheckMinConsSegs <- function(....


##
##Params:
##
##Returns:
##
cleanAndErrorCheckSNPsMPhase <- function(all.samp.SNP.list, min.cons.segs.df, allow.logr.only.flag = FALSE) {
  #TODO add check that every single all.samp.SNP.list has the same SNPs..
  #Note error checking of removing NAs from the min.cons.seg.df should already have been performed
  min.cpn.vals.df <- min.cons.segs.df[, 4:ncol(min.cons.segs.df)]
  seg.name.vec <- paste(min.cons.segs.df[, 1], min.cons.segs.df[, 2], min.cons.segs.df[, 3], sep = ":")
  sample.vec <- colnames(min.cpn.vals.df)
  min.cpn.vals.df <- min.cpn.vals.df[!is.na(rowSums(min.cpn.vals.df)), ,drop = FALSE] 
  working.all.samp.SNP.list <- list()
  all.samps.prob.SNPs.vec <- c()

  for ( sample in sample.vec ) {

    samp.prob.SNPs.vec <- c()

    for ( seg.idx in 1:nrow(min.cpn.vals.df) ) {

      samp.SNP.df <- all.samp.SNP.list[[sample]]

      samp.seg.logR.BAF.df <- samp.SNP.df[ samp.SNP.df$chr == min.cons.segs.df[seg.idx, 1]
                                          & samp.SNP.df$pos >= min.cons.segs.df[seg.idx, 2]
                                          & samp.SNP.df$pos <= min.cons.segs.df[seg.idx, 3], ,drop = FALSE]

      #Catch any snp values that have a BAF but no logr in any region, these are likely incorrect
      # samp.prob.SNPs.vec <- c(samp.prob.SNPs.vec, rownames(samp.seg.logR.BAF.df[ ! is.na(samp.seg.logR.BAF.df$BAFsegmented)
      #                                                                            & (is.na(samp.seg.logR.BAF.df$logRsegmented)
      #                                                                               | is.na(samp.seg.logR.BAF.df$logR))
      #                                                                            ,]))
      #Remove any snps that have no logr but no BAF is the flag is false

      if (allow.logr.only.flag == FALSE) {

        samp.prob.SNPs.vec <- c(samp.prob.SNPs.vec, rownames(samp.seg.logR.BAF.df[samp.seg.logR.BAF.df$BAF == 2,]))
        samp.seg.logR.BAF.df <- samp.seg.logR.BAF.df[as.numeric(samp.seg.logR.BAF.df$BAF) != 2, ,drop = FALSE]
        #samp.seg.logR.BAF.df <- samp.seg.logR.BAF.df[! is.na(samp.seg.logR.BAF.df$BAFsegmented), ,drop=FALSE]

      }
      #Todo,
      #There's one more case we need to consider that a only on or some samples might have a logr SNP or BAF SNP as NA but the others don't, if this were the case
      #We'd want to remove these SNPs
      #Todo,
      #Another case is that we might have any of the SNPs removed might be the start or end position of the joint segs, if this is the case it might break downstream
      #analyses, we should stop the analysis if this is the case and say clean up the input data.
     
      colnames(samp.seg.logR.BAF.df)[colnames(samp.seg.logR.BAF.df) == "chr"] <- "chrom"
      working.all.samp.SNP.list[[sample]][[seg.name.vec[seg.idx]]] <- samp.seg.logR.BAF.df


    }#end for ( seg.idx in 1:nrow(min.cpn.vals.df) )

    all.samps.prob.SNPs.vec <- c(all.samps.prob.SNPs.vec, samp.prob.SNPs.vec)

  }#end for ( sample in samples.vec ) {

  all.samps.prob.SNPs.vec <- unique(all.samps.prob.SNPs.vec)

  print("Number of problem SNPs being removed:")
  print(length(all.samps.prob.SNPs.vec))
  for ( sample in sample.vec ) {

    for (seg.idx in 1:nrow(min.cpn.vals.df)) {

      working.all.samp.SNP.list[[sample]][[seg.name.vec[seg.idx]]] <- working.all.samp.SNP.list[[sample]][[seg.name.vec[seg.idx]]][!rownames(working.all.samp.SNP.list[[sample]][[seg.name.vec[seg.idx]]])%in%all.samps.prob.SNPs.vec, ,drop=FALSE] 

    }#end for (seg.idx in 1:nrow(min.cpn.vals.df))

  }#end for ( sample in sample.vec ) {

  return(working.all.samp.SNP.list)

}#end <- function(all.samp.SNP.list, min.cons.seg.df)

##
##Params:
##
##Returns:
##
cleanAndErrorCheckSNPs <- function(all.samp.SNP.list, min.cons.seg.df, allow.logr.only.flag = FALSE) {
  #TODO add check that every single all.samp.SNP.list has the same SNPs..
  #Note error checking of removing NAs from the min.cons.seg.df should already have been performed
  min.cpn.vals.df <- min.cons.seg.df[,4:ncol(min.cons.seg.df)]
  seg.name.vec <- paste(min.cons.seg.df[, 1], min.cons.segs.df[, 2], min.cons.segs.df[, 3], sep=":")
  min.cpn.vals.df <- min.cpn.vals.df[ ! is.na(rowSums(min.cpn.vals.df)), ,drop = FALSE] 
  sample.vec <- colnames(min.cpn.vals.df)
  working.all.samp.SNP.list <- list()
  all.samps.prob.SNPs.vec <- c()

  for ( sample in sample.vec ) {

    samp.prob.SNPs.vec <- c()

    for ( seg.idx in 1:nrow(min.cpn.vals.df) ) {

      samp.SNP.df <- all.samp.SNP.list[[sample]]

      samp.seg.logR.BAF.df <- samp.SNP.df[ samp.SNP.df$chrom == min.cons.seg.df[seg.idx, 1] 
                                         & samp.SNP.df$pos >= min.cons.seg.df[seg.idx, 2] 
                                         & samp.SNP.df$pos <= min.cons.seg.df[seg.idx, 3], ,drop = FALSE]

      #Catch any snp values that have a BAF but no logr in any region, these are likely incorrect
      # samp.prob.SNPs.vec <- c(samp.prob.SNPs.vec, rownames(samp.seg.logR.BAF.df[ ! is.na(samp.seg.logR.BAF.df$BAFsegmented)
      #                                                                            # & (is.na(samp.seg.logR.BAF.df$logRsegmented)
      #                                                                            #    | is.na(samp.seg.logR.BAF.df$logR))
      #                                                                            ,]))
      #Remove any snps that have no logr but no BAF is the flag is false
      if (allow.logr.only.flag == FALSE) {

        samp.prob.SNPs.vec <- c(samp.prob.SNPs.vec, rownames(samp.seg.logR.BAF.df[is.na(samp.seg.logR.BAF.df$BAFsegmented),]))
        samp.seg.logR.BAF.df <- samp.seg.logR.BAF.df[samp.seg.logR.BAF.df$BAF != 2, ,drop = FALSE]
        samp.seg.logR.BAF.df <- samp.seg.logR.BAF.df[! is.na(samp.seg.logR.BAF.df$BAFsegmented), ,drop = FALSE]

      }#end if (allow.logr.only.flag == FALSE) {
      #Todo,
      #There's one more case we need to consider that a only on or some samples might have a logr SNP or BAF SNP as NA but the others don't, if this were the case
      #We'd want to remove these SNPs
      #Todo,
      #Another case is that we might have any of the SNPs removed might be the start or end position of the joint segs, if this is the case it might break downstream
      #analyses, we should stop the analysis if this is the case and say clean up the input data.
      working.all.samp.SNP.list[[sample]][[seg.name.vec[seg.idx]]] <- samp.seg.logR.BAF.df


    }#end for ( seg.idx in 1:nrow(min.cpn.vals.df) )

    all.samps.prob.SNPs.vec <- c(all.samps.prob.SNPs.vec, samp.prob.SNPs.vec)

  }#end for ( sample in samples.vec ) {

  all.samps.prob.SNPs.vec <- unique(all.samps.prob.SNPs.vec)
  print("Number of problem SNPs being removed:")
  print(length(all.samps.prob.SNPs.vec))
  for ( sample in sample.vec ) {

    for (seg.idx in 1:nrow(min.cpn.vals.df)) {

      working.all.samp.SNP.list[[sample]][[seg.name.vec[seg.idx]]] <- working.all.samp.SNP.list[[sample]][[seg.name.vec[seg.idx]]][ ! rownames(working.all.samp.SNP.list[[sample]][[seg.name.vec[seg.idx]]])%in%all.samps.prob.SNPs.vec, ,drop=FALSE] 

    }#end for (seg.idx in 1:nrow(min.cpn.vals.df))

  }#end for ( sample in sample.vec ) {

  return(working.all.samp.SNP.list)

}#end <- function(all.samp.SNP.list, min.cons.seg.df)
                                                        