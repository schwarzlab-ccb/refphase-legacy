
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


  }#end for ( samp in sample.vec ) {

  write.table(pt.file.paths.df, file = save.output.path ,row.names=FALSE, sep="\t", quote=FALSE)

}#end makeInputPathFileForTRACERxWGS <- function(pt.id, release.dir, array.type, save.output.path) {



# getSNPDataForTRACERxWGS <- function(snp.file.path) {

#   pt.file.paths.df <- read.table(snp.file.path, stringsAsFactors = FALSE, sep = "\t", header = TRUE)
#   sample.vec <- pt.file.paths.df$sample
#   all.regions.ASCAT.list <- list()

#   itx <- iter(sample.vec)
#   #for ( samp in sample.vec ) {
#   all.samp.SNP.df <- foreach(samp = itx, .combine = rbind) %dopar% {

#     print(samp)
#     samp.logr.baf.path <- pt.file.paths.df[pt.file.paths.df$sample == samp ,]$samp_logr_baf_file_path
#     load(samp.logr.baf.path)#gives tumor.baf_logr
#     samp.seg.path <- pt.file.paths.df[pt.file.paths.df$sample == samp ,]$samp_seg_file_path
#     load(samp.seg.path)#gives seg

#     ASCAT.summary               <- tumor.baf_logr[,c("chrom", "pos")]
#     ASCAT.summary$BAF           <- tumor.baf_logr$BAF
#     ASCAT.summary$logR          <- tumor.baf_logr[, samp]
#     ASCAT.summary$BAFsegmented  <- NA
#     ASCAT.summary$logRsegmented <- NA
#     ASCAT.summary$sample <- samp
#     colnames(ASCAT.summary)[colnames(ASCAT.summary) == "chrom"] <- chr

#     for ( j in 1:nrow(seg) ) {

#         if(nrow(ASCAT.summary[ASCAT.summary$chr == seg[j ,]$chr 
#                             & ASCAT.summary$pos >= seg[j ,]$startpos
#                             & ASCAT.summary$pos <= seg[j ,]$endpos,]) > 0) {

#           if(seg[j,]$nMajor==seg[j,]$nMinor){
            
#             ASCAT.summary[ASCAT.summary$chr == seg[j ,]$chr
#                         & ASCAT.summary$pos >= seg[j ,]$startpos 
#                         & ASCAT.summary$pos <= seg[j ,]$endpos ,]$BAFsegmented <- 0.5
            
#           } else {
            
#             avg.seg.baf <- 0.5 - mean(abs(tumor.baf_logr[tumor.baf_logr$chr == seg[j ,]$chr
#                                                        & tumor.baf_logr$pos >= seg[j ,]$startpos
#                                                        & tumor.baf_logr$pos <= seg[j ,]$endpos,]$BAF - 0.5))
            
#             ASCAT.summary[ASCAT.summary$chrom == seg[j ,]$chr 
#                         & ASCAT.summary$pos >= seg[j ,]$startpos 
#                         & ASCAT.summary$pos <= seg[j ,]$endpos ,]$BAFsegmented <- avg.seg.baf
            
#           }

#           avg.seg.logr <- mean(tumor.baf_logr[tumor.baf_logr$chr == seg[j ,]$chr 
#                              & tumor.baf_logr$pos >= seg[j ,]$startpos 
#                              & tumor.baf_logr$pos <= seg[j ,]$endpos, samp])
          
#           ASCAT.summary[ASCAT.summary$chr == seg[j ,]$chr 
#                       & ASCAT.summary$pos >= seg[j ,]$startpos 
#                       & ASCAT.summary$pos <= seg[j ,]$endpos,]$logRsegmented <- avg.seg.logr

#         }#end if(nrow(ASCAT.summary[ASCAT.summary$chrom==seg[j,]$chr & ASCAT.summary$pos...

#       }#end for ( j in 1:nrow(seg) ) {
      
#       ASCAT.summary

#   }# end  for ( samp in sample.vec ) {


#     all.samp.SNP.list <- list()
#     for ( samp in sample.vec ) {

#       samp.w.segmented.info.SNP.df <- all.samp.SNP.df[all.samp.SNP.df$sample==samp,]
#       samp.w.segmented.info.SNP.df <- samp.w.segmented.info.SNP.df[, c("chr", "pos", "BAF", "BAFsegmented", "logR", "logRsegmented")]
#       all.samp.SNP.list[[samp]] <- samp.w.segmented.info.SNP.df

#     }

#     return(all.regions.ASCAT.list)
# }#end getSNPDataForTRACERxWGS





getSNPDataForTRACERxWGS <- function(snp.file.path) {

  pt.file.paths.df <- read.table(snp.file.path, stringsAsFactors = FALSE, sep = "\t", header = TRUE)
  sample.vec <- pt.file.paths.df$sample

  all.regions.ASCAT.SNP.list <- list()
  all.regions.ASCAT.seg.df <- NULL
  print("new WGS reader")
  for (samp in sample.vec) {

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
    colnames(ASCAT.summary)[colnames(ASCAT.summary) == "chrom"] <- "chr"

    if (is.null(all.regions.ASCAT.seg.df)) {

      all.regions.ASCAT.seg.df <- seg

    } else {

      all.regions.ASCAT.seg.df <- rbind(all.regions.ASCAT.seg.df, seg)
    }

    all.regions.ASCAT.SNP.list[[samp]] <- ASCAT.summary

  }


  all.regions.ASCAT.seg.df$start_idx <- NA
  all.regions.ASCAT.seg.df$end_idx <- NA
  all.regions.ASCAT.seg.df$avg_seg_BAF <- NA
  all.regions.ASCAT.seg.df$avg_seg_LogR <-NA
  print(head(all.regions.ASCAT.seg.df))
  itx <- iter(all.regions.ASCAT.seg.df, by="row")

  #for ( samp in sample.vec ) {
  all.samp.SNP.df <- foreach(seg = itx, .combine = rbind) %dopar% {

        #print(seg)
        seg.to.return <- seg

        # print(seg$chr)
        # print(seg$startpos)
        # print(seg$endpos)
        ASCAT.summary <- all.regions.ASCAT.SNP.list[[seg$sample]]

        ASCAT.summary$idx <- 1:nrow(ASCAT.summary)

        if(nrow(ASCAT.summary[ASCAT.summary$chr == seg$chr 
                            & ASCAT.summary$pos >= seg$startpos
                            & ASCAT.summary$pos <= seg$endpos,]) > 0) {

          idx.vec <- ASCAT.summary[ASCAT.summary$chr == seg$chr 
                            & ASCAT.summary$pos >= seg$startpos
                            & ASCAT.summary$pos <= seg$endpos,]$idx
          #print(idx.vec)
          seg.to.return$start_idx <- min(idx.vec)
          seg.to.return$end_idx <- max(idx.vec)

          if(seg$nMajor == seg$nMinor){
            
            # ASCAT.summary[ASCAT.summary$chrom == seg$chr
            #             & ASCAT.summary$pos >= seg$startpos 
            #             & ASCAT.summary$pos <= seg$endpos ,]$BAFsegmented <- 0.5
                        seg.to.return$avg_seg_BAF <- 0.5
            
          } else {
            
            avg.seg.baf <- 0.5 - mean(abs(tumor.baf_logr[tumor.baf_logr$chr == seg$chr
                                                       & tumor.baf_logr$pos >= seg$startpos
                                                       & tumor.baf_logr$pos <= seg$endpos,]$BAF - 0.5))
            
            # ASCAT.summary[ASCAT.summary$chrom == seg$chr 
            #             & ASCAT.summary$pos >= seg$startpos 
            #             & ASCAT.summary$pos <= seg$endpos ,]$BAFsegmented <- avg.seg.baf
                          seg.to.return$avg_seg_BAF <- avg.seg.baf

            
          }

          avg.seg.logr <- mean(tumor.baf_logr[tumor.baf_logr$chr == seg$chr 
                             & tumor.baf_logr$pos >= seg$startpos 
                             & tumor.baf_logr$pos <= seg$endpos, samp])
          
          # ASCAT.summary[ASCAT.summary$chrom == seg$chr 
          #             & ASCAT.summary$pos >= seg$startpos 
          #             & ASCAT.summary$pos <= seg$endpos,]$logRsegmented <- avg.seg.logr
                      seg.to.return$avg_seg_LogR <- avg.seg.logr

        }#end if(nrow(ASCAT.summary[ASCAT.summary$chrom==seg[j,]$chr & ASCAT.summary$pos...
      
      seg.to.return

  }# end  for ( samp in sample.vec ) {
  print("all.samp.SNP.df")
  print(head(all.samp.SNP.df))
  all.samp.SNP.df <- all.samp.SNP.df[with(all.samp.SNP.df, order(sample, chr, startpos)), ]
  all.samp.SNP.list <- list()
  for ( samp in sample.vec ) {

      samp.seg.idx.df <- all.samp.SNP.df[all.samp.SNP.df$sample == samp,]    
      samp.SNP.df <- all.regions.ASCAT.SNP.list[[samp]]

    for ( j in 1:nrow(samp.seg.idx.df) ) {

      if(is.na(samp.seg.idx.df[j,]$start_idx)) {

        next

      }
      samp.SNP.df[as.numeric(samp.seg.idx.df[j,]$start_idx):as.numeric(samp.seg.idx.df[j,]$end_idx),]$BAFsegmented <- samp.seg.idx.df[j,]$avg_seg_BAF
      samp.SNP.df[as.numeric(samp.seg.idx.df[j,]$start_idx):as.numeric(samp.seg.idx.df[j,]$end_idx),]$logRsegmented <- samp.seg.idx.df[j,]$avg_seg_LogR

    }

    all.samp.SNP.list[[samp]] <- samp.SNP.df[, c("chr", "pos", "BAF", "BAFsegmented", "logR", "logRsegmented")]
  }

    return(all.samp.SNP.list)
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
        samp.seg.logR.BAF.df <- samp.seg.logR.BAF.df[! is.na(samp.seg.logR.BAF.df$BAFsegmented), ,drop = FALSE]

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
  #min.cpn.vals.df <- min.cpn.vals.df[ ! is.na(rowSums(min.cpn.vals.df)), ,drop = FALSE] 
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
      # if (allow.logr.only.flag == FALSE) {

      #   samp.prob.SNPs.vec <- c(samp.prob.SNPs.vec, rownames(samp.seg.logR.BAF.df[is.na(samp.seg.logR.BAF.df$BAFsegmented),]))
      #   samp.seg.logR.BAF.df <- samp.seg.logR.BAF.df[samp.seg.logR.BAF.df$BAF != 2, ,drop = FALSE]
      #   samp.seg.logR.BAF.df <- samp.seg.logR.BAF.df[! is.na(samp.seg.logR.BAF.df$BAFsegmented), ,drop = FALSE]

      # }#end if (allow.logr.only.flag == FALSE) {
      #Todo,
      #There's one more case we need to consider that a only on or some samples might have a logr SNP or BAF SNP as NA but the others don't, if this were the case
      #We'd want to remove these SNPs
      #Todo,
      #Another case is that we might have any of the SNPs removed might be the start or end position of the joint segs, if this is the case it might break downstream
      #analyses, we should stop the analysis if this is the case and say clean up the input data.
      working.all.samp.SNP.list[[sample]][[seg.name.vec[seg.idx]]] <- samp.seg.logR.BAF.df


    }#end for ( seg.idx in 1:nrow(min.cpn.vals.df) )

    #all.samps.prob.SNPs.vec <- c(all.samps.prob.SNPs.vec, samp.prob.SNPs.vec)

  }#end for ( sample in samples.vec ) {

  # all.samps.prob.SNPs.vec <- unique(all.samps.prob.SNPs.vec)
  # print("Number of problem SNPs being removed:")
  # print(length(all.samps.prob.SNPs.vec))
  # for ( sample in sample.vec ) {

  #   for (seg.idx in 1:nrow(min.cpn.vals.df)) {

  #     working.all.samp.SNP.list[[sample]][[seg.name.vec[seg.idx]]] <- working.all.samp.SNP.list[[sample]][[seg.name.vec[seg.idx]]][ ! rownames(working.all.samp.SNP.list[[sample]][[seg.name.vec[seg.idx]]])%in%all.samps.prob.SNPs.vec, ,drop=FALSE] 

  #   }#end for (seg.idx in 1:nrow(min.cpn.vals.df))

  # }#end for ( sample in sample.vec ) {

  return(working.all.samp.SNP.list)

}#end <- function(all.samp.SNP.list, min.cons.seg.df)



getChromosomeArmGrange <- function(chr, chr.arm, chrom.info.df) {

  if (chr.arm == "p") {

    arm.df <- data.frame(
    chrom  = chr,
    start = chrom.info.df[chrom.info.df$chrom == chr,]$p_start,
    end = chrom.info.df[chrom.info.df$chrom == chr,]$p_end,
    stringsAsFactors = FALSE)

  } else if (chr.arm == "q") {

    arm.df <- data.frame(
    chrom  = chr,
    start = chrom.info.df[chrom.info.df$chrom ==chr,]$q_start,
    end = chrom.info.df[chrom.info.df$chrom == chr,]$q_end,
    stringsAsFactors = FALSE)

  }

  arm.gr <- makeGRangesFromDataFrame(arm.df,
                            keep.extra.columns=FALSE,
                            ignore.strand=FALSE,
                            seqinfo=NULL,
                            seqnames.field=c("seqnames", "seqname",
                                             "chromosome", "chrom",
                                             "chr", "chromosome_name",
                                             "seqid"),
                            start.field=c("Start_Position","start"),
                            end.field=c("End_Position", "end", "stop"),
                            strand.field="strand",
                            starts.in.df.are.0based=FALSE)

  return(arm.gr)

}


##
##Params:
##
##Returns:
##
getAllPatientSamplesAndIntersectionGRangesList <- function(sample.vec, f.joint.segs.df) {

  #Now what about intersection of all samples what proportion of the genome to they cover:
  samp.gr.list <- list()
  overall.intersection.gr <- NULL
  for ( samp in sample.vec ) {

    samp.gr <- makeGRangesFromDataFrame(f.joint.segs.df[f.joint.segs.df$sample==samp 
                      ,],
                            keep.extra.columns=FALSE,
                            ignore.strand=FALSE,
                            seqinfo=NULL,
                            seqnames.field=c("seqnames", "seqname",
                                             "chromosome", "chrom",
                                             "chr", "chromosome_name",
                                             "seqid"),
                            start.field=c("Start_Position","start","startpos"),
                            end.field=c("End_Position", "end", "stop", "endpos"),
                            strand.field="strand",
                            starts.in.df.are.0based=FALSE)
    samp.gr.list[[samp]] <- samp.gr

    if ( is.null(overall.intersection.gr) ) {

      overall.intersection.gr <- samp.gr

    } else {

      overall.intersection.gr <- intersect(overall.intersection.gr, samp.gr)

    }#end ! is.null (overall.intersection.gr)

  }#en for ( samp in sample.vec ) {
  samp.gr.list[["overall_intersection"]] <- overall.intersection.gr


  return(samp.gr.list)
}

##
##Params:
##
##Returns:
##
compareOverlapsBetweenSamples <- function(sample.vec, f.joint.segs.df, hg19.info.df, chrs.to.use.vec=CHRS_TO_USE_VEC) {

  #For each sample how much of the chromosome arm
  #chrs.to.use.vec <- 1:22
  #Subset the chromos
  chr.info.df <- hg19.info.df[hg19.info.df$chrom %in% chrs.to.use.vec,]
  p.name.vec <- chr.info.df$p_name
  q.name.vec <- chr.info.df$q_name
  arm.vec <- c(rbind(p.name.vec, q.name.vec)) 
  p.q.vec <- c("p", "q")

  samp.gr.list <- getAllPatientSamplesAndIntersectionGRangesList(sample.vec, f.joint.segs.df)
  #So now let's make an arm coverage data.frame
  arm.coverage.df <- data.frame(matrix(vector(), length(arm.vec), (length(sample.vec) + 2),
                          dimnames = list(c(), c("arm", sample.vec, "overall_intersection"))),
                          stringsAsFactors = FALSE,
                          check.names = FALSE)
  arm.coverage.df$arm <- arm.vec
  #how 

  #So now we want to take the intersection and compare the sizes of the

  for ( samp in names(samp.gr.list) ) {

    #Get sample granges:
    samp.gr <- samp.gr.list[[samp]]
    #print(samp)
    #print(samp.gr)
    for ( chr in chrs.to.use.vec ) {
      for ( arm in p.q.vec) {

        chr.arm.gr <- getChromosomeArmGrange(chr, arm, chr.info.df)
        chr.arm.overlap.gr <- intersect(samp.gr, chr.arm.gr)
        arm.coverage.df[arm.coverage.df$arm==paste0(chr, arm), samp] <- sum(as.numeric(width(chr.arm.overlap.gr))) / as.numeric(width(chr.arm.gr))

      }#end for ( arm in p.q.vec) {
    }# for ( chr in chrs.to.use.vec ) {

  }#end for ( samp in sample.vec ) {

  arm.coverage.df$diff_from_max <-  apply(arm.coverage.df[,2:ncol(arm.coverage.df)], 1, max) - arm.coverage.df$overall_intersection

  #Report max and union

  return(arm.coverage.df)
}


joinTwoGrangesIntoOne <- function(gr.one, gr.two) {

  #Removee any metadata columns
  mcols(gr.one) <- NULL
  mcols(gr.two) <- NULL

  chr <- unique(as.character(seqnames(gr.one)), as.character(seqnames(gr.two)))
  comb.grs <- c(gr.two, gr.one)
  s.comb.grs <- sort(comb.grs)
  if (length(chr) > 1) {

    print("More than one chromsome given as input, returning input unjoined")
    return(s.comb.grs)

  }

  new.start <- start(s.comb.grs)[1]
  new.end <- end(s.comb.grs)[length(end(s.comb.grs))]

  joined.gr <- GRanges(seqnames = chr, 
             ranges = IRanges(new.start, new.end))

  # gr.list <- GRangesList(gr.one,gr.two)
  # comb.grs <- do.call(getMethod(c, "GenomicRanges"), gr.list)
  return(joined.gr)
  #the preceding segment
}


removeSingleBasePairSegmentsFromIntersectionGranges <- function(all.grs, min.width = 1) {

  #Calculate width:
  all.grs$width <- width(all.grs)
  #Get grs that are one
  width.one.grs <- all.grs[ all.grs$width <= min.width ]

  working.no.ones.grs <- all.grs
  mcols(working.no.ones.grs) <- NULL
  #print("width.one.grs")
  #print(width.one.grs)

  if(length(width.one.grs) == 0) {

    return(working.no.ones.grs)

    
  }
  for ( idx in 1:length(width.one.grs) ) {

    gr.ante.idx <- follow(width.one.grs[idx], working.no.ones.grs)
    gr.post.idx <- precede(width.one.grs[idx], working.no.ones.grs)
    merge.poss.idxs <- c(gr.ante.idx, gr.post.idx)
    non.na.merge.poss.idxs <- merge.poss.idxs[! is.na(merge.poss.idxs)]

    #Check thats
    if ( length(non.na.merge.poss.idxs) == 0 ) {

      print("No ante or post segments")
      next()

    }
    ante.post.grs <- working.no.ones.grs[non.na.merge.poss.idxs]
    nearest.gr <- ante.post.grs[nearest(width.one.grs[idx], ante.post.grs)]
    joined.gr <- joinTwoGrangesIntoOne(nearest.gr, width.one.grs[idx])
    all.idx.vec <- 1:length(working.no.ones.grs)
    hits.obj <- findOverlaps(joined.gr, working.no.ones.grs)
    hit.idx.vec <- subjectHits(hits.obj)
    non.overlapping.idx.vec <- all.idx.vec[ ! all.idx.vec %in% hit.idx.vec]
    non.overlapping.grs <- working.no.ones.grs[non.overlapping.idx.vec]

    working.no.ones.grs <- sort(c(non.overlapping.grs, joined.gr))

  }
  #so now let's subset the granges
  return(working.no.ones.grs)
}


##
##Params:
##
##Returns:
##
getChromosomeArmsGRangeList <- function(current.chr, chrom.info.df) {

  arm.df.list <- list() 
  p.arm.df <- data.frame(
  chrom  = current.chr,
  start = chrom.info.df[chrom.info.df$chrom == current.chr,]$p_start,
  end = chrom.info.df[chrom.info.df$chrom == current.chr,]$p_end,
  stringsAsFactors = FALSE)

  count <- 1

  if(chrom.info.df[chrom.info.df$chrom == current.chr,]$acrocentric == FALSE) {

    arm.df.list[[count]] <- makeGRangesFromDataFrame(p.arm.df,
                             keep.extra.columns=FALSE,
                             ignore.strand=FALSE,
                             seqinfo=NULL,
                             seqnames.field=c("seqnames", "seqname",
                                              "chromosome", "chrom",
                                              "chr", "chromosome_name",
                                              "seqid"),
                             start.field=c("Start_Position","start"),
                             end.field=c("End_Position", "end", "stop"),
                             strand.field="strand",
                             starts.in.df.are.0based=FALSE)

    count <- count + 1
  }

  q.arm.df <- data.frame(
  chrom  = current.chr,
  start = chrom.info.df[chrom.info.df$chrom == current.chr,]$q_start,
  end = chrom.info.df[chrom.info.df$chrom == current.chr,]$q_end,
  stringsAsFactors = FALSE)

  arm.df.list[[count]] <- makeGRangesFromDataFrame(q.arm.df,
                           keep.extra.columns=FALSE,
                           ignore.strand=FALSE,
                           seqinfo=NULL,
                           seqnames.field=c("seqnames", "seqname",
                                            "chromosome", "chrom",
                                            "chr", "chromosome_name",
                                            "seqid"),
                           start.field=c("Start_Position","start"),
                           end.field=c("End_Position", "end", "stop"),
                           strand.field="strand",
                           starts.in.df.are.0based=FALSE)

  return(arm.df.list)
}



## Overlap types
## upstream_overlap: startup & endin
## Q --------------
## S       -------------
##
## downstream_overlap: startin & enddown
## Q       -------------
## S --------------
##
## query_contained_in_subject: startin & endin 
## Q     -----
## S --------------
##
## subject_inside_query: startup & enddown
## Q --------------
## S     -----
#Params:
# query: one or more granges/iranges that you wish to compare to the subject
# subject: one or more grangaes/ranges that you wish to be compared with the query
# output: specifies the output type either "gr" for genomic range or "df" for dataframe 
#Returns a dataframe/grange with following columns/fields:
# seqnames: this is just Grange name for chrom of the input GRange
# start: start of the input GRange
# end: end of the input GRange
# width: width of the input GRange
# strand: strand (mostly useless)
# chrom: seqname of input renamed as chrom
# query_hit_idx: idx of the GRange from input query GRange for this row
# subject_hit_idx: matching idx of the overlapping output Subject segment found for this row
# query_start: genomic start of query
# query_end: genomic end of query
# subject_start: genomic start of matching segment
# subject_end: genomic end of matching segment
# overlap_start: genomic start of overlap
# overlap_end: genomic end of overlap
# overlap_length: length of overlap
# overlap_prop_query: proportion of query Grange covered by overlap
# overlap_prop_subject: proportion of subject/result grange covered by overlap
# overlap_type: what kind of overlap it is (edited)
# 
#Also returns all the metadata in the subject, ie original columns holding non-genomic position information
#Todo add widths for everything
getOverlapGRanges <- function(query, subject, output="gr", ...) {
        require(GenomicRanges); require(IRanges)
        
        ## Input check
        if(!((class(query) == "GRanges" & class(subject) == "GRanges") | (class(query) == "IRanges" & class(subject) == "IRanges"))) {
                stop("Query and subject need to be of same class, either GRanges or IRanges!")
        }
        
        ## Find overlapping ranges
        if(class(query) == "GRanges") {
          seqlengths(query) <- rep(NA, length(seqlengths(query)))
          seqlengths(subject) <- rep(NA, length(seqlengths(subject)))
        }
  overlap.index.mat <- as.matrix(findOverlaps(query, subject, ...))

        query <- query[overlap.index.mat[, 1]]
        subject <- subject[overlap.index.mat[, 2]]
        overlap.mat <- cbind(query_start = start(query), query_end = end(query), subject_start = start(subject), subject_end = end(subject))
  
        ## Pre-queries for overlaps
  startup <- overlap.mat[,"subject_start"] < overlap.mat[,"query_start"]
  enddown <- overlap.mat[,"subject_end"] > overlap.mat[,"query_end"]
  startin <- overlap.mat[,"subject_start"] >= overlap.mat[,"query_start"] & overlap.mat[,"subject_start"] <= overlap.mat[, "query_end"]
  endin <- overlap.mat[,"subject_end"] >= overlap.mat[,"query_start"] & overlap.mat[,"subject_end"] <=  overlap.mat[, "query_end"]

  ## Overlap types
  olup <- startup & endin
  oldown <- startin & enddown
  inside <- startin & endin 
  contained <- startup & enddown
  
        ## Overlap types in one vector
  overlap_type <- rep("", length(overlap.mat[,"query_start"]))
  overlap_type[olup] <- "upstream_overlap"
  overlap_type[oldown] <- "downstream_overlap"
  overlap_type[inside] <- "subject_inside_query" 
  overlap_type[contained] <- "query_contained_in_subject"
  
        ## Overlap positions
  overlap_start <- rep(0, length(overlap.mat[, "query_start"]))
  overlap_end <- rep(0, length(overlap.mat[, "query_start"]))
  overlap_start[olup] <- overlap.mat[, "query_start"][olup]
  overlap_end[olup] <- overlap.mat[, "subject_end"][olup]
  overlap_start[oldown] <- overlap.mat[, "subject_start"][oldown]
  overlap_end[oldown] <- overlap.mat[, "query_end"][oldown]
  overlap_start[inside] <- overlap.mat[, "subject_start"][inside]
  overlap_end[inside] <- overlap.mat[, "subject_end"][inside]
  overlap_start[contained] <- overlap.mat[, "query_start"][contained]
  overlap_end[contained] <- overlap.mat[, "query_end"][contained]

  ## Absolute and relative length of overlaps
  overlap_length <- (overlap_end - overlap_start) + 1
        overlap_prop_query <- overlap_length / width(query) * 100
        overlap_prop_subject <- overlap_length / width(subject) * 100
  
  ## Output type
        oldf <- data.frame(query_hit_idx=overlap.index.mat[,1], subject_hit_idx=overlap.index.mat[,2], overlap.mat, overlap_start, overlap_end, overlap_length, overlap_prop_query, overlap_prop_subject, overlap_type)
        if(class(query) == "GRanges") {
    oldf <- cbind(chrom=as.character(seqnames(query)), oldf)
  }
        if(output=="df") {
                return(cbind(oldf, as.data.frame(mcols(subject))))
        }
        if(output=="gr") {
                if(class(query)=="GRanges") {
                        elementMetadata(query) <- cbind(as.data.frame(elementMetadata(query)), oldf, mcols(subject))
                }
                if(class(query)=="IRanges") {#Todo add the metadata in th eiranges 
                        query <- GRanges(seqnames = Rle(rep("dummy", length(query))), ranges = IRanges(start=oldf[,"query_start"], end=oldf[,"query_end"]), strand = Rle(strand(rep("+", length(query)))), oldf)  
                }
                return(query)
        }
}








# ##
# ##Params:
# ##
# ##Returns:
# ##
# cleanAndErrorCheckSNPsParallel <- function(all.samp.SNP.list, min.cons.seg.df, allow.logr.only.flag = FALSE, num.cores = 8) {
#   #TODO add check that every single all.samp.SNP.list has the same SNPs..
#   #Note error checking of removing NAs from the min.cons.seg.df should already have been performed

#   min.cpn.vals.df <- min.cons.seg.df[,4:ncol(min.cons.seg.df)]
#   seg.name.vec <- paste(min.cons.seg.df[, 1], min.cons.segs.df[, 2], min.cons.segs.df[, 3], sep=":")
#   #min.cpn.vals.df <- min.cpn.vals.df[ ! is.na(rowSums(min.cpn.vals.df)), ,drop = FALSE] 
#   working.all.samp.SNP.list <- list()
#   all.samps.prob.SNPs.vec <- c()

#   #Make vectors for sample vec:
#   samp.seg.vec <- rep(seg.name.vec, length(sample.vec))
#   sample.rep.vec <- c()
#   for ( samp in sample.vec ) {

#     sample.rep.vec <- c(sample.rep.vec, rep(samp, length(seg.name.vec)))

#   }
#   iterable.min.cons.seg.df <- cbind(do.call("rbind", replicate(length(sample.vec), min.cons.seg.df[,1:3], simplify = FALSE)), sample = sample.rep.vec)

#   # cl <- makeCluster(num.cores)
#   # registerDoParallel(cl)
#   itx <- iter(iterable.min.cons.seg.df, by="row")
#   print("Staring parallel assignation of snps to segs:")
#   results.list <- foreach(min.cons.seg.df = itx, .combine = c) %dopar% {

#     samp.SNP.df <- all.samp.SNP.list[[ min.cons.seg.df[, 4]]]

#     samp.seg.logR.BAF.df <- samp.SNP.df[ samp.SNP.df$chrom == min.cons.seg.df[, 1] 
#                                        & samp.SNP.df$pos >= min.cons.seg.df[, 2] 
#                                        & samp.SNP.df$pos <= min.cons.seg.df[, 3], ,drop = FALSE]

#     if (allow.logr.only.flag == FALSE) {

#       samp.prob.SNPs.vec <- rownames(samp.seg.logR.BAF.df[is.na(samp.seg.logR.BAF.df$BAFsegmented),])
#       samp.seg.logR.BAF.df <- samp.seg.logR.BAF.df[samp.seg.logR.BAF.df$BAF != 2, ,drop = FALSE]
#       samp.seg.logR.BAF.df <- samp.seg.logR.BAF.df[! is.na(samp.seg.logR.BAF.df$BAFsegmented), ,drop = FALSE]

#     }
#     samp.seg.logR.BAF.df.list <- list(samp.seg.logR.BAF.df, samp.prob.SNPs.vec)
#     names(samp.seg.logR.BAF.df.list) <- paste0(min.cons.seg.df[, 4], "@", min.cons.seg.df[, 1], ":", min.cons.seg.df[, 2], ":", min.cons.seg.df[, 3])

#     return(samp.seg.logR.BAF.df.list)                               
#   }

#   #Create the list dataobject we're expecting:

#   for ( result.name in names(results.list) ) {

#     samp <- strsplit(result.name, "@")[[1]][1]
#     seg.name <- strsplit(result.name, "@")[[1]][2]

#     working.all.samp.SNP.list[[samp]][[seg.name]] <- results.list[[result.name]][[1]]

#   }
#   #Now get all the problem snps
#   all.samps.prob.SNPs.vec <- sapply(results.list, function(x) x[[1]][2])
#   all.samps.prob.SNPs.vec <- unique(all.samps.prob.SNPs.vec)
#   print("Number of problem SNPs being removed:")
#   print(length(all.samps.prob.SNPs.vec))


#   for ( sample in sample.vec ) {

#     for (seg.idx in 1:nrow(min.cpn.vals.df)) {

#       working.all.samp.SNP.list[[sample]][[seg.name.vec[seg.idx]]] <- working.all.samp.SNP.list[[sample]][[seg.name.vec[seg.idx]]][ ! rownames(working.all.samp.SNP.list[[sample]][[seg.name.vec[seg.idx]]])%in%all.samps.prob.SNPs.vec, ,drop=FALSE] 

#     }#end for (seg.idx in 1:nrow(min.cpn.vals.df))

#   }#end for ( sample in sample.vec ) {

#   return(working.all.samp.SNP.list)

# }#end <- function(all.samp.SNP.list, min.cons.seg.df)
 

# # ##
# ##Params:
# ##
# ##Returns:
# ##
# cleanAndErrorCheckSNPsParallel <- function(all.samp.SNP.list, min.cons.seg.df, allow.logr.only.flag = FALSE, num.cores = 8) {
#   # TODO add check that every single all.samp.SNP.list has the same SNPs..
#   # Note error checking of removing NAs from the min.cons.seg.df should already have been performed

#   min.cpn.vals.df <- min.cons.seg.df[,4:ncol(min.cons.seg.df)]
#   seg.name.vec <- paste(min.cons.seg.df[, 1], min.cons.segs.df[, 2], min.cons.segs.df[, 3], sep=":")
#   #min.cpn.vals.df <- min.cpn.vals.df[ ! is.na(rowSums(min.cpn.vals.df)), ,drop = FALSE] 
#   working.all.samp.SNP.list <- list()
#   all.samps.prob.SNPs.vec <- c()

#   #Make vectors for sample vec:
#   samp.seg.vec <- rep(seg.name.vec, length(sample.vec))
#   sample.rep.vec <- c()
#   for ( samp in sample.vec ) {

#     working.all.samp.SNP.list[[samp]] <- list()
#     sample.rep.vec <- c(sample.rep.vec, rep(samp, length(seg.name.vec)))

#   }

#   iterable.min.cons.seg.df <- cbind(do.call("rbind", replicate(length(sample.vec), min.cons.seg.df[,1:3], simplify = FALSE)), sample = sample.rep.vec)
#   cl <- makeCluster(NUM_CORES)
#   registerDoParallel(cl)
#   itx <- iter(iterable.min.cons.seg.df, by="row")
#   print("Starting parallel assignation of snps to segs:")
#   # Create class which holds multiple results for each loop iteration.
#   # Each loop iteration populates two properties: $result1 and $result2.
#   # For a great tutorial on S3 classes, see: 
#   # http://www.cyclismo.org/tutorial/R/s3Classes.html#creating-an-s3-class
#   MultiResultClass <- function(samp_seg_name=NULL, samp_seg_logR_BAF_df=NULL, samp_prob_SNPs_vec=NULL)
#   {
#     me <- list(
#       samp_seg_name = samp_seg_name,
#       samp_seg_logR_BAF_df = samp_seg_logR_BAF_df,
#       samp_prob_SNPs_vec  = samp_prob_SNPs_vec
#     )
#     # Set the name for the class
#     class(me) <- append(class(me),"MultiResultClass")
#     return(me)
#   }
#   sites.in.each.samp.seg.list <- foreach(min.cons.seg.df = itx, .combine = c) %dopar% {

#     samp.SNP.df <- all.samp.SNP.list[[ min.cons.seg.df[, 4]]]

#     samp.seg.logR.BAF.df <- samp.SNP.df[ samp.SNP.df$chrom == min.cons.seg.df[, 1] 
#                                        & samp.SNP.df$pos >= min.cons.seg.df[, 2] 
#                                        & samp.SNP.df$pos <= min.cons.seg.df[, 3], ,drop = FALSE]

#     if (allow.logr.only.flag == FALSE) {

#       samp.prob.SNPs.vec <- rownames(samp.seg.logR.BAF.df[is.na(samp.seg.logR.BAF.df$BAFsegmented),])
#       samp.seg.logR.BAF.df <- samp.seg.logR.BAF.df[samp.seg.logR.BAF.df$BAF != 2, ,drop = FALSE]
#       samp.seg.logR.BAF.df <- samp.seg.logR.BAF.df[! is.na(samp.seg.logR.BAF.df$BAFsegmented), ,drop = FALSE]

#     } else {

#       samp.prob.SNPs.vec <- c()

#     }

#     result <- MultiResultClass(samp_seg_name = paste0(min.cons.seg.df[, 4], "@", min.cons.seg.df[, 1], ":", min.cons.seg.df[, 2], ":", min.cons.seg.df[, 3]),
#                   samp_seg_logR_BAF_df = samp.seg.logR.BAF.df,
#                   samp_prob_SNPs_vec = samp.prob.SNPs.vec)
#     result.list <- list(result)
#     names(result.list) <- paste0(min.cons.seg.df[, 4], "@", min.cons.seg.df[, 1], ":", min.cons.seg.df[, 2], ":", min.cons.seg.df[, 3])
#     result.list

#   }

#   all.samps.prob.SNPs.vec <- unique(unlist(sapply(sites.in.each.samp.seg.list, function(x) x[["samp_prob_SNPs_vec"]])))
#   print("Number of problem SNPs being removed:")
#   print(length(all.samps.prob.SNPs.vec))

#   if (length(all.samps.prob.SNPs.vec) > 0) {

#     samp.seg.itx <- iter(sites.in.each.samp.seg.list)
#     clean.sites.in.each.samp.seg.list <- foreach(samp.seg.mrc = samp.seg.itx, .combine = c) %dopar% {

#       seg.name <- strsplit(samp.seg.mrc$samp_seg_name, "@")[[1]][2]
#       samp.seg.mrc.list <- sites.in.each.samp.seg.list[grep(seg.name, names(sites.in.each.samp.seg.list), value = TRUE)]
#       samp.seg.prob.snp.vec <- unique(unlist(sapply(samp.seg.mrc.list, function(x) x[["samp_prob_SNPs_vec"]])))
#       samp.seg.logR.BAF.df <-  samp.seg.mrc$samp_seg_logR_BAF_df[ ! rownames(samp.seg.mrc$samp_seg_logR_BAF_df) %in% samp.seg.prob.snp.vec, ,drop = FALSE ]
#       result <- MultiResultClass(samp_seg_name = samp.seg.mrc$samp_seg_name,
#                     samp_seg_logR_BAF_df = samp.seg.logR.BAF.df,
#                     samp_prob_SNPs_vec = c())
#       result.list <- list(result)
#       names(result.list) <- samp.seg.mrc$samp_seg_name
#       result.list
#     }

#   } else {

#     clean.sites.in.each.samp.seg.list <- sites.in.each.samp.seg.list

#   }

#   stopCluster(cl)

#   for ( samp.seg.name in names(clean.sites.in.each.samp.seg.list) ) {

#     samp <- strsplit(samp.seg.name, "@")[[1]][1]
#     seg.name <- strsplit(samp.seg.name, "@")[[1]][2]
#     working.all.samp.SNP.list[[samp]][[seg.name]] <- clean.sites.in.each.samp.seg.list[[samp.seg.name]]$samp_seg_logR_BAF_df

#   }

#   return(working.all.samp.SNP.list)

# }#end cleanAndErrorCheckSNPsParallel <- function(all.samp.SNP.list, min.cons.seg.df,

  ##
    ##Params:
    ##
    ##Returns:
    ##
    readInMEDICCInputFormatParallel <- function(snp.file.df, processed.orig.cpn.segs.df, sample.vec, num.cores = 8) {
      # TODO add check that every single all.samp.SNP.list has the same SNPs..
      # Note error checking of removing NAs from the min.cons.seg.df should already have been performed
      # registerDoParallel(cores = num.cores)
      # print(foreach::getDoParWorkers())
      # print(getDoParName())
      # print(getDoParVersion())

     itx <- iter(sample.vec)
     all.samp.SNP.df <- foreach(samp = itx, .combine = rbind) %dopar% {

      print(samp)
      #print(head(snp.file.df))

      ASCAT.summary <- data.frame(chr = as.numeric(snp.file.df$Chr),
                            pos = as.numeric(snp.file.df$Pos),
                      BAF = as.numeric(snp.file.df[, paste0(gsub("-" , "\\.", samp), "_BAF")]),
                      BAFsegmented = rep(NA, length(snp.file.df$Chr)),
                      logR = as.numeric(snp.file.df[, paste0(gsub("-", "\\.", samp), "_LogR")]),
                      logRsegmented = rep(NA, length(snp.file.df$Chr)),
                      stringsAsFactors = FALSE)
      ASCAT.summary$sample <- samp
      #print(head(ASCAT.summary))

      if(nrow(processed.orig.cpn.segs.df) == 0) {

        stop("processed.orig.cpn.segs.df has no rows...")
      }
      seg <- processed.orig.cpn.segs.df[processed.orig.cpn.segs.df$sample == samp,]
      if(nrow(processed.orig.cpn.segs.df) == 0) {

        stop("processed.orig.cpn.segs.df subsetted by sample has no rows...")

      }
      for ( j in 1:nrow(seg) ) {
        #print(seg[j,])
        #Are there any BAF values that aren't set at 2 (and thus shouldn't be used for allele specific cpn) present in this segment?
        if (nrow(ASCAT.summary[ASCAT.summary$chr == seg[j ,]$chr
                    & ASCAT.summary$pos >= seg[j ,]$startpos 
                    & ASCAT.summary$pos <= seg[j ,]$endpos 
                    & ASCAT.summary$BAF != 2 ,]) > 0) {
          #Is this segment balanced, if so set the BAFsegmented to 0.5
          if (seg[j , ,]$nMajor == seg[j ,]$nMinor) {

            ASCAT.summary[ASCAT.summary$chr == seg[j ,]$chr
                    & ASCAT.summary$pos >= seg[j ,]$startpos
                    & ASCAT.summary$pos <= seg[j ,]$endpos,]$BAFsegmented <- 0.5
            
          } else {#Segment is balanced so we need to calculate the value of the BAF segmented and apply it to all these values 
          

            avg.seg.baf <- 0.5 - mean(abs(as.numeric(ASCAT.summary[ASCAT.summary$chr == seg[j ,]$chr
                                        & ASCAT.summary$pos >= seg[j ,]$startpos
                                        & ASCAT.summary$pos <= seg[j ,]$endpos
                                        & ASCAT.summary$BAF != 2 , ]$BAF)          - 0.5))

            ASCAT.summary[ASCAT.summary$chr == seg[j ,]$chr 
                    & ASCAT.summary$pos >= seg[j ,]$startpos 
                    & ASCAT.summary$pos <= seg[j ,]$endpos, ]$BAFsegmented <- avg.seg.baf

          }
     
        }#end if(nrow(ASCAT.summary[ASCAT.summary$chrom==seg[j,]$chr & ASCAT.summary$pos...
        #We now need to calculate the average logr of any snps that are present
        if (nrow(ASCAT.summary[ASCAT.summary$chr == seg[j ,]$chr
                             & ASCAT.summary$pos >= seg[j ,]$startpos
                             & ASCAT.summary$pos <= seg[j ,]$endpos , ]) > 0) {

          avg.seg.logr <- mean(ASCAT.summary[ASCAT.summary$chr == seg[j ,]$chr
                     & ASCAT.summary$pos >= seg[j ,]$startpos
                     & ASCAT.summary$pos <= seg[j ,]$endpos, ]$logR)

          ASCAT.summary[ASCAT.summary$chr == seg[j ,]$chr
                  & ASCAT.summary$pos >= seg[j ,]$startpos
                  & ASCAT.summary$pos <= seg[j ,]$endpos , ]$logRsegmented <- avg.seg.logr

        }

      }#end for ( j in 1:nrow(seg) ) {

      ASCAT.summary

    }#end for ( samp in sample.vec ) {


    all.samp.SNP.list <- list()
    for ( samp in sample.vec ) {

      samp.w.segmented.info.SNP.df <- all.samp.SNP.df[all.samp.SNP.df$sample==samp,]
      samp.w.segmented.info.SNP.df <- samp.w.segmented.info.SNP.df[, c("chr", "pos", "BAF", "BAFsegmented", "logR", "logRsegmented")]
      all.samp.SNP.list[[samp]] <- samp.w.segmented.info.SNP.df

    }


      return(all.samp.SNP.list)
  }#end cleanAndErrorCheckSNPsParallel <- function(all.samp.SNP.list, min.cons.seg.df



##
##Params:
##
##Returns:
##
cleanAndErrorCheckSNPsParallel <- function(all.samp.SNP.list, min.cons.seg.df, allow.logr.only.flag = FALSE, data.type = "MEDICC") {
  # TODO add check that every single all.samp.SNP.list has the same SNPs..
  # Note error checking of removing NAs from the min.cons.seg.df should already have been performed

  min.cpn.vals.df <- min.cons.seg.df[,4:ncol(min.cons.seg.df)]
  seg.name.vec <- paste(min.cons.seg.df[, 1], min.cons.segs.df[, 2], min.cons.segs.df[, 3], sep=":")
  #min.cpn.vals.df <- min.cpn.vals.df[ ! is.na(rowSums(min.cpn.vals.df)), ,drop = FALSE] 
  working.all.samp.SNP.list <- list()
  all.samps.prob.SNPs.vec <- c()

  #Make vectors for sample vec:
  samp.seg.vec <- rep(seg.name.vec, length(sample.vec))
  sample.rep.vec <- c()
  for ( samp in sample.vec ) {

    working.all.samp.SNP.list[[samp]] <- list()
    sample.rep.vec <- c(sample.rep.vec, rep(samp, length(seg.name.vec)))

  }

  iterable.min.cons.seg.df <- cbind(do.call("rbind", replicate(length(sample.vec), min.cons.seg.df[,1:3], simplify = FALSE)), sample = sample.rep.vec)
  #cl <- makeCluster(NUM_CORES)
  #registerDoParallel(cl)
  itx <- iter(iterable.min.cons.seg.df, by="row")
  print("Starting parallel assignation of snps to segs:")
  # Create class which holds multiple results for each loop iteration.
  # Each loop iteration populates two properties: $result1 and $result2.
  # For a great tutorial on S3 classes, see: 
  # http://www.cyclismo.org/tutorial/R/s3Classes.html#creating-an-s3-class
  MultiResultClass <- function(samp_seg_name=NULL, samp_seg_logR_BAF_df=NULL, samp_prob_SNPs_vec=NULL)
  {
    me <- list(
      samp_seg_name = samp_seg_name,
      samp_seg_logR_BAF_df = samp_seg_logR_BAF_df,
      samp_prob_SNPs_vec  = samp_prob_SNPs_vec
    )
    # Set the name for the class
    class(me) <- append(class(me),"MultiResultClass")
    return(me)
  }
  sites.in.each.samp.seg.list <- foreach(min.cons.seg.df = itx, .combine = c) %dopar% {
  #for( i in 1:nrow(iterable.min.cons.seg.df)) {

    #print(iterable.min.cons.seg.df[i,])
    #min.cons.seg.df <- iterable.min.cons.seg.df[i, ,drop=FALSE]
   # print(min.cons.seg.df)
    samp.SNP.df <- all.samp.SNP.list[[ min.cons.seg.df[, 4]]]

    if( "chr" %in% colnames(samp.SNP.df)) {

          colnames(samp.SNP.df)[colnames(samp.SNP.df) == "chr"] <- "chrom"

    }

    samp.seg.logR.BAF.df <- samp.SNP.df[ samp.SNP.df$chrom == min.cons.seg.df[, 1] 
                                       & samp.SNP.df$pos >= min.cons.seg.df[, 2] 
                                       & samp.SNP.df$pos <= min.cons.seg.df[, 3], ,drop = FALSE]
    samp.prob.SNPs.vec <- c()
    if (data.type == "TRACERx_WGS") {

      samp.prob.SNPs.vec <- c()

    } else if (allow.logr.only.flag == FALSE & data.type == "TRACERx_exome") {

      samp.prob.SNPs.vec <- rownames(samp.seg.logR.BAF.df[is.na(samp.seg.logR.BAF.df$BAFsegmented),])
      samp.seg.logR.BAF.df <- samp.seg.logR.BAF.df[samp.seg.logR.BAF.df$BAF != 2, ,drop = FALSE]
      samp.seg.logR.BAF.df <- samp.seg.logR.BAF.df[! is.na(samp.seg.logR.BAF.df$BAFsegmented), ,drop = FALSE]

    } else if (allow.logr.only.flag == FALSE & data.type == "MEDICC") {

      samp.prob.SNPs.vec <- c(samp.prob.SNPs.vec, rownames(samp.seg.logR.BAF.df[samp.seg.logR.BAF.df$BAF == 2,]))
      samp.seg.logR.BAF.df <- samp.seg.logR.BAF.df[as.numeric(samp.seg.logR.BAF.df$BAF) != 2, ,drop = FALSE]
      samp.seg.logR.BAF.df <- samp.seg.logR.BAF.df[! is.na(samp.seg.logR.BAF.df$BAFsegmented), ,drop = FALSE]
      #colnames(samp.seg.logR.BAF.df)[colnames(samp.seg.logR.BAF.df) == "chr"] <- "chrom"

    } else {

      samp.prob.SNPs.vec <- c()
    
    }

    result <- MultiResultClass(samp_seg_name = paste0(min.cons.seg.df[, 4], "@", min.cons.seg.df[, 1], ":", min.cons.seg.df[, 2], ":", min.cons.seg.df[, 3]),
                  samp_seg_logR_BAF_df = samp.seg.logR.BAF.df,
                  samp_prob_SNPs_vec = samp.prob.SNPs.vec)
    result.list <- list(result)
    names(result.list) <- paste0(min.cons.seg.df[, 4], "@", min.cons.seg.df[, 1], ":", min.cons.seg.df[, 2], ":", min.cons.seg.df[, 3])
    result.list

  }

  all.samps.prob.SNPs.vec <- unique(unlist(sapply(sites.in.each.samp.seg.list, function(x) x[["samp_prob_SNPs_vec"]])))
  print("Number of problem SNPs being removed:")
  print(length(all.samps.prob.SNPs.vec))
  print(head(all.samps.prob.SNPs.vec))

  if (length(all.samps.prob.SNPs.vec) > 0) {

    samp.seg.itx <- iter(sites.in.each.samp.seg.list)
    clean.sites.in.each.samp.seg.list <- foreach(samp.seg.mrc = samp.seg.itx, .combine = c) %dopar% {

      seg.name <- strsplit(samp.seg.mrc$samp_seg_name, "@")[[1]][2]
      samp.seg.mrc.list <- sites.in.each.samp.seg.list[grep(seg.name, names(sites.in.each.samp.seg.list), value = TRUE)]
      samp.seg.prob.snp.vec <- unique(unlist(sapply(samp.seg.mrc.list, function(x) x[["samp_prob_SNPs_vec"]])))
      samp.seg.logR.BAF.df <-  samp.seg.mrc$samp_seg_logR_BAF_df[ ! rownames(samp.seg.mrc$samp_seg_logR_BAF_df) %in% samp.seg.prob.snp.vec, ,drop = FALSE ]
      result <- MultiResultClass(samp_seg_name = samp.seg.mrc$samp_seg_name,
                    samp_seg_logR_BAF_df = samp.seg.logR.BAF.df,
                    samp_prob_SNPs_vec = c())
      result.list <- list(result)
      names(result.list) <- samp.seg.mrc$samp_seg_name
      result.list
    }

  } else {

    clean.sites.in.each.samp.seg.list <- sites.in.each.samp.seg.list

  }

  #stopCluster(cl)

  for ( samp.seg.name in names(clean.sites.in.each.samp.seg.list) ) {

    samp <- strsplit(samp.seg.name, "@")[[1]][1]
    seg.name <- strsplit(samp.seg.name, "@")[[1]][2]
    working.all.samp.SNP.list[[samp]][[seg.name]] <- clean.sites.in.each.samp.seg.list[[samp.seg.name]]$samp_seg_logR_BAF_df

  }

  return(working.all.samp.SNP.list)

}#end cleanAndErrorCheckSNPsParallel <- function(all.samp.SNP.list, min.cons.seg.df,












 # print("HERE")
 # for ( sample in sample.vec ) {

 #   samp.prob.SNPs.vec <- c()

 #   for ( seg.idx in 1:nrow(min.cpn.vals.df) ) {

 #     print(sample)
 #     print(seg.idx)
 #     samp.SNP.df <- all.samp.SNP.list[[sample]]

 #     samp.seg.logR.BAF.df <- samp.SNP.df[ samp.SNP.df$chrom == min.cons.seg.df[seg.idx, 1] 
 #                                        & samp.SNP.df$pos >= min.cons.seg.df[seg.idx, 2] 
 #                                        & samp.SNP.df$pos <= min.cons.seg.df[seg.idx, 3], ,drop = FALSE]

 #     #Catch any snp values that have a BAF but no logr in any region, these are likely incorrect
 #     samp.prob.SNPs.vec <- c(samp.prob.SNPs.vec, rownames(samp.seg.logR.BAF.df[ ! is.na(samp.seg.logR.BAF.df$BAFsegmented)
 #                                                                                & (is.na(samp.seg.logR.BAF.df$logRsegmented)
 #                                                                                   | is.na(samp.seg.logR.BAF.df$logR))
 #                                                                                ,]))
 #     #Remove any snps that have no logr but no BAF is the flag is false
 #     if (allow.logr.only.flag == FALSE) {

 #       samp.prob.SNPs.vec <- c(samp.prob.SNPs.vec, rownames(samp.seg.logR.BAF.df[is.na(samp.seg.logR.BAF.df$BAFsegmented),]))
 #       samp.seg.logR.BAF.df <- samp.seg.logR.BAF.df[samp.seg.logR.BAF.df$BAF != 2, ,drop = FALSE]
 #       samp.seg.logR.BAF.df <- samp.seg.logR.BAF.df[! is.na(samp.seg.logR.BAF.df$BAFsegmented), ,drop = FALSE]

 #     }#end if (allow.logr.only.flag == FALSE) {
 #     #Todo,
 #     #There's one more case we need to consider that a only on or some samples might have a logr SNP or BAF SNP as NA but the others don't, if this were the case
 #     #We'd want to remove these SNPs
 #     #Todo,
 #     #Another case is that we might have any of the SNPs removed might be the start or end position of the joint segs, if this is the case it might break downstream
 #     #analyses, we should stop the analysis if this is the case and say clean up the input data.
 #     working.all.samp.SNP.list[[sample]][[seg.name.vec[seg.idx]]] <- samp.seg.logR.BAF.df


 #   }#end for ( seg.idx in 1:nrow(min.cpn.vals.df) )

 #   all.samps.prob.SNPs.vec <- c(all.samps.prob.SNPs.vec, samp.prob.SNPs.vec)

 # }#end for ( sample in samples.vec ) {


      # all.samp.SNP.list <- list()
    # for ( samp in sample.vec ) {

    #   print(samp)
    #   print(head(snp.file.df))

    #   ASCAT.summary <- data.frame(chr = as.numeric(snp.file.df$Chr),
    #                         pos = as.numeric(snp.file.df$Pos),
    #                   BAF = as.numeric(snp.file.df[, paste0(gsub("-" , "\\.", samp), "_BAF")]),
    #                   BAFsegmented = rep(NA, length(snp.file.df$Chr)),
    #                   logR = as.numeric(snp.file.df[, paste0(gsub("-", "\\.", samp), "_LogR")]),
    #                   logRsegmented = rep(NA, length(snp.file.df$Chr)),
    #                   stringsAsFactors = FALSE)

    #   print(head(ASCAT.summary))

    #   seg <- processed.orig.cpn.segs.df[processed.orig.cpn.segs.df$sample == samp,]

    #   for ( j in 1:nrow(seg) ) {
    #     #print(seg[j,])
    #     #Are there any BAF values that aren't set at 2 (and thus shouldn't be used for allele specific cpn) present in this segment?
    #     if (nrow(ASCAT.summary[ASCAT.summary$chr == seg[j ,]$chr
    #                 & ASCAT.summary$pos >= seg[j ,]$startpos 
    #                 & ASCAT.summary$pos <= seg[j ,]$endpos 
    #                 & ASCAT.summary$BAF != 2 ,]) > 0) {
    #       #Is this segment balanced, if so set the BAFsegmented to 0.5
    #       if (seg[j ,]$nMajor == seg[j ,]$nMinor) {

    #         ASCAT.summary[ASCAT.summary$chr == seg[j ,]$chr
    #                 & ASCAT.summary$pos >= seg[j ,]$startpos
    #                 & ASCAT.summary$pos <= seg[j ,]$endpos,]$BAFsegmented <- 0.5
            
    #       } else {#Segment is balanced so we need to calculate the value of the BAF segmented and apply it to all these values 
          

    #         avg.seg.baf <- 0.5 - mean(abs(as.numeric(ASCAT.summary[ASCAT.summary$chr == seg[j ,]$chr
    #                                     & ASCAT.summary$pos >= seg[j ,]$startpos
    #                                     & ASCAT.summary$pos <= seg[j ,]$endpos
    #                                     & ASCAT.summary$BAF != 2 ,]$BAF)          - 0.5))

    #         ASCAT.summary[ASCAT.summary$chr == seg[j ,]$chr 
    #                 & ASCAT.summary$pos >= seg[j ,]$startpos 
    #                 & ASCAT.summary$pos <= seg[j ,]$endpos ,]$BAFsegmented <- avg.seg.baf

    #       }
     
    #     }#end if(nrow(ASCAT.summary[ASCAT.summary$chrom==seg[j,]$chr & ASCAT.summary$pos...
    #     #We now need to calculate the average logr of any snps that are present
    #     if (nrow(ASCAT.summary[ASCAT.summary$chr == seg[j ,]$chr
    #                          & ASCAT.summary$pos >= seg[j ,]$startpos
    #                          & ASCAT.summary$pos <= seg[j ,]$endpos ,]) > 0) {

    #       avg.seg.logr <- mean(ASCAT.summary[ASCAT.summary$chr == seg[j ,]$chr
    #                  & ASCAT.summary$pos >= seg[j ,]$startpos
    #                  & ASCAT.summary$pos <= seg[j ,]$endpos, ]$logR)

    #       ASCAT.summary[ASCAT.summary$chr == seg[j ,]$chr
    #               & ASCAT.summary$pos >= seg[j ,]$startpos
    #               & ASCAT.summary$pos <= seg[j ,]$endpos ,]$logRsegmented <- avg.seg.logr

    #     }

    #   }#end for ( j in 1:nrow(seg) ) {

    #   all.samp.SNP.list[[samp]] <- ASCAT.summary

    # }#end for ( samp in sample.vec ) {

