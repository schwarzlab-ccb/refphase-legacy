
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

  }#end for ( samp in sample.vec ) {


  write.table(pt.file.paths.df, file = save.output.path , row.names = FALSE, sep = "\t", quote = FALSE)

}#end makeInputPathFileForTRACERxExomeSNPArray

getSNPDataForTRACERxExomeSNPArray <- function(pt.file.paths.df, joint.segs.df) {

  # sample.vec <- pt.file.paths.df$sample
  # all.samples.ASCAT.list <- list()

  # for ( samp in sample.vec ) {

  #   samp.ASCAT.obj <- pt.file.paths.df[pt.file.paths.df$sample == samp, ]$samp_ascat_bc_path
  #   #/camp/project/proj-tracerx-lung/tracerx/release_tx500_180615//B_LTX094//exome/ASCAT/B_LTX094_SU_T1-R3/B_LTX094_SU_T1-R3.ascat.bc.RData
  #   samp.ASCAT.obj <- load(samp.ASCAT.obj)
  #   samp.ASCAT.obj <- get(samp.ASCAT.obj)

  #   samp.ASCAT.summary               <- samp.ASCAT.obj$SNPpos
  #   samp.ASCAT.summary$BAF           <- NA
  #   samp.ASCAT.summary$BAFsegmented  <- NA
  #   samp.ASCAT.summary$logR          <- NA
  #   samp.ASCAT.summary$logRsegmented <- NA

  #   samp.ASCAT.summary[rownames(samp.ASCAT.obj$Tumor_LogR_segmented),]$logRsegmented           <- samp.ASCAT.obj$Tumor_LogR_segmented[,1]
  #   samp.ASCAT.summary[rownames(samp.ASCAT.obj$Tumor_BAF),]$BAF                                <- samp.ASCAT.obj$Tumor_BAF[[1]]
  #   samp.ASCAT.summary[rownames(samp.ASCAT.obj$Tumor_BAF_segmented[[1]]),]$BAFsegmented        <- samp.ASCAT.obj$Tumor_BAF_segmented[[1]][,1]
  #   samp.ASCAT.summary[rownames(samp.ASCAT.obj$Tumor_LogR),]$logR                              <- samp.ASCAT.obj$Tumor_LogR[[1]]


  #   all.samples.ASCAT.list[[samp]] <- samp.ASCAT.summary

  # }#end getSNPDataForTRACERxExomeSNPArray

  # return(all.samples.ASCAT.list)


  # pt.dir <- paste(release.dir, pt.id, "/", array.type, "/ASCAT/", sep="")
  # regions <- grep("FAILED",grep(".txt$|.complete$|*EXCLUDE*|*ManualQC*", list.files(pt.dir), value = TRUE, invert = TRUE), invert = TRUE, value = TRUE)
  regions <- pt.file.paths.df$sample
  all.regions.ASCAT.list <- list()

  samp.ASCAT.obj.path <- pt.file.paths.df[pt.file.paths.df$sample == regions[1], ]$samp_ascat_bc_path
  normal.logr.baf.path <- gsub(".ascat.bc.", ".normal.baf_logr.", samp.ASCAT.obj.path)

  if(file.exists(normal.logr.baf.path)){

      load(normal.logr.baf.path)

  } else {

    normal.logr.baf.path <- gsub("ASCAT", "ASCAT_QC", normal.logr.baf.path)
    load(normal.logr.baf.path)
  }
  #Gives normal.baf_logr
  is.hom.BAF.vec <- as.numeric(normal.baf_logr[, "BAF"]) < 0.3 | as.numeric(normal.baf_logr[, "BAF"]) > 0.7
  is.not.het.BAF.vec <- is.na(normal.baf_logr[, "BAF"]) | is.hom.BAF.vec

  for ( i in 1:length(regions) ) {

    samp.ASCAT.obj.path <- pt.file.paths.df[pt.file.paths.df$sample == regions[i], ]$samp_ascat_bc_path
    reg.logr.baf.path <- gsub(".ascat.bc.", ".tumor.baf_logr.", samp.ASCAT.obj.path)
    #reg.seg.path <- gsub(".bc.", ".seg.", samp.ASCAT.obj.path)

    #reg.logr.baf.path <- paste(DirToConsider, regions[i], "/", regions[i], ".tumor.baf_logr.RData", sep="")
    cat("Attempting to load:", reg.logr.baf.path, "\n")
    if (file.exists(reg.logr.baf.path)) {

      load(reg.logr.baf.path)

    } else {

      reg.logr.baf.path <- gsub("ASCAT", "ASCAT_QC", reg.logr.baf.path)
      load(reg.logr.baf.path)

    }
    #gives tumor.baf_logr obj

    #reg.seg.path <- paste(DirToConsider, regions[i], "/", regions[i], ".ascat.seg.RData", sep="")
    #load(reg.seg.path)#gives reg.seg.path

    seg <- joint.segs.df[joint.segs.df$sample %in% regions[i],]

    ASCAT.summary               <- tumor.baf_logr[, c("chrom", "pos")]
    ASCAT.summary$BAF           <- tumor.baf_logr$BAF
    ASCAT.summary[is.not.het.BAF.vec,]$BAF <- NA#Making homozygous SNPs by ASCAT's definition. read as NAs for now.
    ASCAT.summary$logR          <- tumor.baf_logr[, regions[i]]
    ASCAT.summary$BAFsegmented  <- NA
    ASCAT.summary$logRsegmented <- NA


    #Update the fake ASCAT.summaary BAFsegmetned and logRsegmented
    for ( j in 1:nrow(seg) ) {

      print(seg[j,])

      print(nrow(ASCAT.summary[ASCAT.summary$chrom == seg[j ,]$chr
                            & ASCAT.summary$pos >= seg[j ,]$startpos
                            & ASCAT.summary$pos <= seg[j ,]$endpos ,]))

      if(nrow(ASCAT.summary[ASCAT.summary$chrom == seg[j ,]$chr
                            & ASCAT.summary$pos >= seg[j ,]$startpos
                            & ASCAT.summary$pos <= seg[j ,]$endpos ,]) > 0) {


        if(seg[j,]$nMajor == seg[j,]$nMinor) {

          ASCAT.summary[ASCAT.summary$chrom == seg[j ,]$chr
                        & ASCAT.summary$pos >= seg[j ,]$startpos
                        & ASCAT.summary$pos <= seg[j ,]$endpos ,]$BAFsegmented <- 0.5

                        print("0.5")

        } else {#end if(seg[j,]$nMajor == seg[j,]$nMinor) {

          if(nrow(ASCAT.summary[ASCAT.summary$chrom == seg[j ,]$chr
                            & ASCAT.summary$pos >= seg[j ,]$startpos
                            & ASCAT.summary$pos <= seg[j ,]$endpos
                            & ! is.na(ASCAT.summary$BAF) ,]) > 0 ) {

            avg.seg.baf <- 0.5 - mean(abs(tumor.baf_logr[tumor.baf_logr$chrom == seg[j ,]$chr
                                                       & tumor.baf_logr$pos >= seg[j ,]$startpos
                                                       & tumor.baf_logr$pos <= seg[j ,]$endpos
                                                       & ! is.na(ASCAT.summary$BAF),]$BAF - 0.5))

            ASCAT.summary[ASCAT.summary$chrom == seg[j ,]$chr
                      & ASCAT.summary$pos >= seg[j ,]$startpos
                      & ASCAT.summary$pos <= seg[j ,]$endpos ,]$BAFsegmented <- avg.seg.baf

                      print(avg.seg.baf)
          } else {#end if nrow(ASCAT.summary[ASCAT.summary$chrom == seg[j ,]$chr....


             ASCAT.summary[ASCAT.summary$chrom == seg[j ,]$chr
                      & ASCAT.summary$pos >= seg[j ,]$startpos
                      & ASCAT.summary$pos <= seg[j ,]$endpos ,]$BAFsegmented <- 0.5

          }

        }#end



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


# Read in a set of ascat.bc files and create the requisite all.samp.SNP.list
#
# > x <- readRDS("LTX121_all_samp_SNP_list.RDS")
# > str(x)
# List of 2
#  $ A_LTX121_SU_T1-R1:'data.frame':      938946 obs. of  6 variables:
#   ..$ chrom        : num [1:938946] 1 1 1 1 1 1 1 1 1 1 ...
#   ..$ pos          : int [1:938946] 762273 762496 762497 777237 777318 777335 777361 777456 783071 783147 ...
#   ..$ BAF          : num [1:938946] NA NA NA NA NA NA NA NA NA NA ...
#   ..$ logR         : num [1:938946] -0.049 -0.161 -0.161 -0.542 -0.46 -0.46 -0.46 -0.166 -0.047 -0.013 ...
#   ..$ BAFsegmented : num [1:938946] 0.231 0.231 0.231 0.231 0.231 ...
#   ..$ logRsegmented: num [1:938946] -0.327 -0.327 -0.327 -0.327 -0.327 ...
#
readInASCAT <- function(sample_table, segs, samples) {
  all_snps <- list()
  print(sample_table)

  for (sample in samples) {
    sample_info <- sample_table[sample,]
    cat("Starting to read sample", sample_info$samp_ascat_bc_path, "\n")

    if(grepl("RData", sample_info$samp_ascat_bc_path)) {
      print("ascat.bc object in RData file, loading...")
      load(sample_info$samp_ascat_bc_path)
      ascat_bc <- ascat.bc
    } else {
      ascat_bc <- readRDS(sample_info$samp_ascat_bc_path)
    }

    validate_ascat_bc_input(ascat_bc)

    # First check to see which SNPs are heterozygous
    normal_baf <- ascat_bc$Germline_BAF[[1]]
    is_homozygous <- normal_baf < 0.3 | normal_baf > 0.7
    is_heterozygous <- !(is.na(normal_baf) | is_homozygous)
    cat("# of heterozygous SNPs:", sum(is_heterozygous), "\n")

    # Mask non-heterozygous SNPs with NA
    baf <- ascat_bc$Tumor_BAF[[1]]
    baf[!is_heterozygous] <- NA

    chrom <- chrom2integer(ascat_bc$SNPpos$chr)

    samp_data <- data.frame(chrom = chrom,
                            pos = ascat_bc$SNPpos$pos,
                            BAF = baf,
                            logR = ascat_bc$Tumor_LogR[[1]],
                            logRsegmented = as.numeric(ascat_bc$Tumor_LogR_segmented),
                            snp_pos = paste0(chrom, "_", ascat_bc$SNPpos$pos),
                            stringsAsFactors = FALSE)

    all_snps[[sample]] <- samp_data

  }

  #So, in TRACERx we have the issue that ASCAT, when run with exactly the same input ascat occasionally does not seem to have
  #all the same positions. So, in order to get xphase to run I have to make sure that the same snps/positions are kept for all
  #samples. Given the way M/XPhase works (finding the nearest SNP to a breakpoint) this should be fine even if we get rid of
  #some positions that are the start or end indices of samples' segments.

  all_snps <- ascat_set_bafsegmented(segs, all_snps)
  all_snps <- keep_common_snps(all_snps)

  return(all_snps)

}

# Subset snp data frames to keep only those snps that are present in all of the
# others
keep_common_snps <- function(all_snps) {
  message("Subsetting to shared SNPs only")

  samples <- names(all_snps)

  # Identify all SNPs that are present in every sample
  snp.pos.vec <- all_snps[[1]]$snp_pos
  for (sample in samples[-1]) {
    snp.pos.vec <- snp.pos.vec[snp.pos.vec %in% all_snps[[sample]]$snp_pos]
  }

  # Remove SNP positions that are not in all samples
  for (sample in samples) {
    all_snps[[sample]] <- all_snps[[sample]][ all_snps[[sample]]$snp_pos %in% snp.pos.vec, ]
  }

  return(all_snps)
}

# Read in a single ascat.bc file and create the requisite gl.snp.df that is used to compare the allelic 
# imbalance present in tumour regions to the exact same set of hz snps in the germline. 
#
# > x <- readRDS("LTX121_all_samp_SNP_list.RDS")
# > str(x)
# List of 2
#  $ A_LTX121_SU_T1-R1:'data.frame':      938946 obs. of  6 variables:
#   ..$ chrom        : num [1:938946] 1 1 1 1 1 1 1 1 1 1 ...
#   ..$ pos          : int [1:938946] 762273 762496 762497 777237 777318 777335 777361 777456 783071 783147 ...
#   ..$ BAF          : num [1:938946] NA NA NA NA NA NA NA NA NA NA ...
#   ..$ logR         : num [1:938946] -0.049 -0.161 -0.161 -0.542 -0.46 -0.46 -0.46 -0.166 -0.047 -0.013 ...
#   ..$ BAFsegmented : num [1:938946] 0.231 0.231 0.231 0.231 0.231 ...
#   ..$ logRsegmented: num [1:938946] -0.327 -0.327 -0.327 -0.327 -0.327 ...
#
readInGermlineFromASCAT <- function(ascat_bc_rds_path, snp_pos_vec) {

  cat("Starting to read sample", ascat_bc_rds_path, "\n")

  if(grepl("RData", ascat_bc_rds_path)) {

    print("ascat.bc object in RData file, loading...")
    load(ascat_bc_rds_path)
    ascat_bc <- ascat.bc

  } else {

    ascat_bc <- readRDS(ascat_bc_rds_path)


  }

  # First check to see which SNPs are heterozygous
  normal_baf <- ascat_bc$Germline_BAF[[1]]
  is_homozygous <- normal_baf < 0.3 | normal_baf > 0.7
  is_heterozygous <- !(is.na(normal_baf) | is_homozygous)
  cat("# of heterozygous SNPs:", sum(is_heterozygous), "\n")

  # Mask non-heterozygous SNPs with NA
  baf <- ascat_bc$Germline_BAF[[1]]
  baf[!is_heterozygous] <- NA

  chrom <- chrom2integer(as.character(ascat_bc$SNPpos$chr))
  gl_data_df <- data.frame(chrom = chrom,
                          pos = ascat_bc$SNPpos$pos,
                          BAF = baf,
                          logR = ascat_bc$Germline_LogR[[1]],
                          # FIXME: huge hack. BAFsegmented does not contain all SNPs, so we
                          # substitute the original BAFs here. This is only used for plotting I
                          # think, so it shouldn't alter the CN estimation.
                          BAFsegmented = rep(as.numeric(NA), length(baf)),
                          #samp_data$BAFsegmented = ascat_bc$Tumor_BAF_segmented[is_heterozygous]
                          logRsegmented = as.numeric(ascat_bc$Tumor_LogR_segmented),
                          snp_pos = paste0(chrom, "_", ascat_bc$SNPpos$pos),
                          stringsAsFactors = FALSE)

  final_gl_data_df <- gl_data_df[gl_data_df$snp_pos %in% snp_pos_vec,]

  return(gl_data_df)
}


# Since the ASCAT data files I have don't include a BAFsegmented value, I
# calculate this from scratch. There is one special value:
#
# if ASCAT calls the region as balanced (e.g. 1-1 or 2-2 allele-specific copy numbers)
# -> force BAFsegmented to be 0.5
#
# this is used in the decision of which sample to use as a phasing reference,
# and whether a segment should be phased at all (it won't be if all samples
# have ASCAT calling the region balanced)
ascat_set_bafsegmented <- function(segs, all_snps) {
  message("Calculating segmented BAF values.")

  sample_names <- names(all_snps)
  for (samp_idx in seq_along(all_snps)) {
    samp_segs <- segs[segs$sample == sample_names[[samp_idx]],]
    samp_data <- all_snps[[samp_idx]]
    for (seg_idx in seq_len(nrow(samp_segs))) {

      seg <- samp_segs[seg_idx,]

      # Segments include both the start and end positions
      seg_data_idx <- (samp_data$chrom == seg$chr &
                       samp_data$pos >= seg$startpos &
                       samp_data$pos <= seg$endpos)

      # If ASCAT doesn't call allelic imbalance, then we force BAFsegmented to
      # 0.5 to indicate this. Otherwise, we use the mean BAF as the BAF for the
      # segment
      #
      # Also deal with segments with NA copy number (comes from sequenza)
      if (is.na(seg$nMajor) | is.na(seg$nMinor) | seg$nMajor == seg$nMinor) {
        all_snps[[samp_idx]][seg_data_idx, "BAFsegmented"] <- 0.5
      } else {
        mean_baf <- 0.5 - mean(abs(samp_data[seg_data_idx, "BAF"] - 0.5), na.rm = TRUE)
        all_snps[[samp_idx]][seg_data_idx, "BAFsegmented"] <- mean_baf
      }
    }
  }

  return(all_snps)
}


ascat_set_logrsegmented <- function(segs, all_snps) {
  message("Calculating segmented logR values.")

  sample_names <- names(all_snps)
  for (samp_idx in seq_along(all_snps)) {
    samp_segs <- segs[segs$sample == sample_names[[samp_idx]],]
    samp_data <- all_snps[[samp_idx]]
    for (seg_idx in seq_len(nrow(samp_segs))) {

      seg <- samp_segs[seg_idx,]

      # Segments include both the start and end positions
      seg_data_idx <- (samp_data$chrom == seg$chr &
                       samp_data$pos >= seg$startpos &
                       samp_data$pos <= seg$endpos)

      mean_logr <- mean(samp_data[seg_data_idx, "logR"], na.rm = TRUE)
      all_snps[[samp_idx]][seg_data_idx, "logRsegmented"] <- mean_logr
    }
  }

  return(all_snps)
}


##
##Params:
##
##Returns:
##
readInMEDICCInputFormatParallel <- function(snp.file.df, processed.orig.cpn.segs.df, sample.vec) {
      # TODO add check that every single all.samp.SNP.list has the same SNPs..
      # Note error checking of removing NAs from the min.cons.seg.df should already have been performed
      # registerDoParallel(cores = num.cores)
    print(foreach::getDoParWorkers())
      # print(getDoParName())
      # print(getDoParVersion())
     print("Calling readInMEDICCInputFormatParallel")
     print("Reading in sites data for:")
     print(sample.vec)
     print("Starting parallel work:")

     num.pre.parallel.samps <- length(unique(sample.vec))
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
      print(paste0("Returning: ", samp))
      print(head(ASCAT.summary))
      ASCAT.summary
      # returnable.list <- list(ASCAT.summary)
      # names(returnable.list) <- samp
      # returnable.list


    }#end for ( samp in sample.vec ) {
    print("After parallel read in:")
    print(head(all.samp.SNP.df))
    print(unique(all.samp.SNP.df$sample))
    print(length(unique(all.samp.SNP.df$sample)))
    print(unique(all.samp.SNP.df$sample))

    num.post.parallel.samps <- length(unique(all.samp.SNP.df$sample))

    if (num.post.parallel.samps != num.pre.parallel.samps) {

      print("parallelisation is not working properly, some parallel workers are not returning jobs correctly....")
      print("Either investigate why some cores are not returning results, change to a different node or run in series.")
      print("Running in series now...")
    }

    print("Assigning to list in expected format:")
    all.samp.SNP.list <- list()
    for ( samp in sample.vec ) {
      print(samp)
      samp.w.segmented.info.SNP.df <- all.samp.SNP.df[all.samp.SNP.df$sample==samp ,]
      print(head(samp.w.segmented.info.SNP.df))
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
cleanAndErrorCheckSNPsParallel <- function(all.samp.SNP.list, mcs.seg.to.use.df, allow.logr.only.flag = FALSE, data.type = "MEDICC") {
  # Note error checking of removing NAs from the min.cons.seg.df should already have been performed

  min.cpn.vals.df <- mcs.seg.to.use.df[,4:ncol(mcs.seg.to.use.df)]
  seg.name.vec <- paste(mcs.seg.to.use.df[, 1], mcs.seg.to.use.df[, 2], mcs.seg.to.use.df[, 3], sep=":")
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

  iterable.min.cons.seg.df <- cbind(do.call("rbind", replicate(length(sample.vec), mcs.seg.to.use.df[,1:3], simplify = FALSE)), sample = sample.rep.vec)
  library("gdata")
  iterable.min.cons.seg.df <- drop.levels(iterable.min.cons.seg.df)
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
  sites.in.each.samp.seg.list <- foreach(it.min.cons.seg.df = itx, .combine = c) %dopar% {

    samp.SNP.df <- all.samp.SNP.list[[ it.min.cons.seg.df[, 4]]]

    if( "chr" %in% colnames(samp.SNP.df)) {

          colnames(samp.SNP.df)[colnames(samp.SNP.df) == "chr"] <- "chrom"

    }
    samp.seg.logR.BAF.df <- samp.SNP.df[ samp.SNP.df$chrom == it.min.cons.seg.df[, 1]
                                       & samp.SNP.df$pos >= it.min.cons.seg.df[, 2]
                                       & samp.SNP.df$pos <= it.min.cons.seg.df[, 3], ,drop = FALSE]
    samp.prob.SNPs.vec <- c()
    if (data.type == "TRACERx_WGS") {

      samp.prob.SNPs.vec <- c()

    } else if (allow.logr.only.flag == FALSE & data.type == "TRACERx_exome") {

      samp.prob.SNPs.vec <- rownames(samp.seg.logR.BAF.df[is.na(samp.seg.logR.BAF.df$BAFsegmented), ,drop = FALSE ])
      samp.seg.logR.BAF.df <- samp.seg.logR.BAF.df[samp.seg.logR.BAF.df$BAF != 2, ,drop = FALSE]
      samp.seg.logR.BAF.df <- samp.seg.logR.BAF.df[! is.na(samp.seg.logR.BAF.df$BAFsegmented), ,drop = FALSE]

    } else if (allow.logr.only.flag == FALSE & data.type == "MEDICC") {

      samp.prob.SNPs.vec <- c(samp.prob.SNPs.vec, rownames(samp.seg.logR.BAF.df[samp.seg.logR.BAF.df$BAF == 2, , drop = FALSE ]))
      samp.seg.logR.BAF.df <- samp.seg.logR.BAF.df[as.numeric(samp.seg.logR.BAF.df$BAF) != 2, ,drop = FALSE]
      samp.seg.logR.BAF.df <- samp.seg.logR.BAF.df[! is.na(samp.seg.logR.BAF.df$BAFsegmented), ,drop = FALSE]
      #colnames(samp.seg.logR.BAF.df)[colnames(samp.seg.logR.BAF.df) == "chr"] <- "chrom"

    } else if (allow.logr.only.flag == FALSE & data.type == "ASCAT") {

      samp.prob.SNPs.vec <- rownames(samp.seg.logR.BAF.df)[is.na(samp.seg.logR.BAF.df$BAF)]
      samp.seg.logR.BAF.df <- samp.seg.logR.BAF.df[!is.na(samp.seg.logR.BAF.df$BAF), ,drop = FALSE]

    } else if (allow.logr.only.flag == FALSE & data.type == "sequenza") {

      samp.prob.SNPs.vec <- rownames(samp.seg.logR.BAF.df)[is.na(samp.seg.logR.BAF.df$BAF)]
      samp.seg.logR.BAF.df <- samp.seg.logR.BAF.df[!is.na(samp.seg.logR.BAF.df$BAF), ,drop = FALSE]

    } else {

      samp.prob.SNPs.vec <- c()

    }

    result <- MultiResultClass(samp_seg_name = paste0(it.min.cons.seg.df[, 4], "@", it.min.cons.seg.df[, 1], ":", it.min.cons.seg.df[, 2], ":", it.min.cons.seg.df[, 3]),
                  samp_seg_logR_BAF_df = samp.seg.logR.BAF.df,
                  samp_prob_SNPs_vec = samp.prob.SNPs.vec)
    result.list <- list(result)
    names(result.list) <- paste0(it.min.cons.seg.df[, 4], "@", it.min.cons.seg.df[, 1], ":", it.min.cons.seg.df[, 2], ":", it.min.cons.seg.df[, 3])
    result.list

  }

  #saveRDS(sites.in.each.samp.seg.list, file = )
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


  for ( samp.seg.name in names(clean.sites.in.each.samp.seg.list) ) {

    samp <- strsplit(samp.seg.name, "@")[[1]][1]
    seg.name <- strsplit(samp.seg.name, "@")[[1]][2]
    working.all.samp.SNP.list[[samp]][[seg.name]] <- clean.sites.in.each.samp.seg.list[[samp.seg.name]]$samp_seg_logR_BAF_df

  }

  return(working.all.samp.SNP.list)

}#end cleanAndErrorCheckSNPsParallel <- function(all.samp.SNP.list, min.cons.seg.df,


# -------------------------------------------------------------------------
# Sequenza input format
#
# Read BAF and logR values from seqz files that are generated for sequenza.
# These files are seqz format, a tab-delimited and gzipped file format with the
# following columns:
#
# - chromosome
# - position
# - base.ref
# - depth.normal
# - depth.tumor
# - depth.ratio
# - Af
# - Bf
# - zygosity.normal
# - GC.percent
# - good.reads
# - AB.normal
# - AB.tumor
# - tumor.strand
#
# example:
# $ tabix CB1001-11-D01-X1.gz "chr1:10,000,000-20,000,000;" | head
# chr1    10007515        C       15      5       0.333   1.0     0       hom     32      5       C       .       0
# chr1    10007516        T       15      5       0.333   1.0     0       hom     32      5       T       .       0
# chr1    10007517        T       16      5       0.312   1.0     0       hom     32      5       T       .       0
# chr1    10007518        C       16      5       0.312   1.0     0       hom     32      5       C       .       0
# chr1    10007519        A       15      5       0.333   1.0     0       hom     32      5       A       .       0
# chr1    10007520        C       17      5       0.294   1.0     0       hom     32      5       C       .       0
#
# The object that is returned should have the format:
# > x <- readRDS("LTX121_all_samp_SNP_list.RDS")
# > str(x)
# List of 2
#  $ A_LTX121_SU_T1-R1:'data.frame':      938946 obs. of  6 variables:
#   ..$ chrom        : num [1:938946] 1 1 1 1 1 1 1 1 1 1 ...
#   ..$ pos          : int [1:938946] 762273 762496 762497 777237 777318 777335 777361 777456 783071 783147 ...
#   ..$ BAF          : num [1:938946] NA NA NA NA NA NA NA NA NA NA ...
#   ..$ logR         : num [1:938946] -0.049 -0.161 -0.161 -0.542 -0.46 -0.46 -0.46 -0.166 -0.047 -0.013 ...
#   ..$ BAFsegmented : num [1:938946] 0.231 0.231 0.231 0.231 0.231 ...
#   ..$ logRsegmented: num [1:938946] -0.327 -0.327 -0.327 -0.327 -0.327 ...
read_snps_sequenza <- function(input_files, segs, samples) {

  library(readr)

  all_snps <- list()
  gl_snps <- list()

  for (sample in samples) {

    input_file <- input_files[[sample]]

    # All data that is available in a .seqz file. We only need a subset of it
    # but this is a useful reference in case we need more info at a later date
    #
    # col_types = "ciciidddcddccc"
    # col_names = c("chromosome", "position", "base.ref", "depth.normal",
    #               "depth.tumor", "depth.ratio", "Af", "Bf", "zygosity.normal",
    #               "GC.percent", "good.reads", "AB.normal", "AB.tumor",
    #               "tumor.strand")
    # Using readr
    message("Reading SNP data for sample: ", sample)
    col_types <- "ci---d-dc-----"  # Skip unnecessary columns
    col_names <- c("chromosome", "position", "depth.ratio", "Bf", "zygosity.normal")
    seqz_data <- read_tsv(input_file, col_types = col_types,
                          progress = TRUE, col_names = col_names,
                          skip = 1)

    message("Creating SNP data frame.")
    chrom <- chrom2integer(seqz_data$chromosome)
    pos <- seqz_data$position
    all_snps[[sample]] <- data.frame(chrom = chrom,
                                     pos = pos,
                                     BAF = seqz_data$Bf,
                                     logR = log2(seqz_data$depth.ratio),
                                     zygosity.normal = seqz_data$zygosity.normal,
                                     snp_pos = paste0(chrom, "_", pos),
                                     stringsAsFactors = FALSE)

    # Mask non-hz SNPs with NA
    all_snps[[sample]]$BAF[all_snps[[sample]]$zygosity.normal != "het"] <- NA

  }

  all_snps <- ascat_set_bafsegmented(segs, all_snps)
  all_snps <- ascat_set_logrsegmented(segs, all_snps)
  all_snps <- keep_common_snps(all_snps)

  str(all_snps)

  return(list(all_samples = all_snps, germline = gl_snps))
}


# Convert the internal refphase snps object into an ascat.bc object for use in
# internal rerun of ASCAT
#
# example object format:
# List of 14
#  $ Tumor_LogR          :'data.frame':   675960 obs. of  1 variable:
#   ..$ CB1001-11-D01-X1: num [1:675960] 0.573 0.602 0.622 0.77 0.505 ...
#  $ Tumor_BAF           :'data.frame':   675960 obs. of  1 variable:
#   ..$ BAF: num [1:675960] 0.0435 0 0 0 0 ...
#  $ Tumor_LogR_segmented: num [1:675960, 1] -0.838 -0.838 -0.838 -0.838 -0.838 ...
#   ..- attr(*, "dimnames")=List of 2
#   .. ..$ : chr [1:675960] "1" "2" "3" "4" ...
#   .. ..$ : chr "CB1001-11-D01-X1"
#  $ Tumor_BAF_segmented :List of 1
#   ..$ : num [1:29679, 1] 0.104 0.104 0.104 0.104 0.104 ...
#   .. ..- attr(*, "dimnames")=List of 2
#   .. .. ..$ : chr [1:29679] "250" "274" "278" "314" ...
#   .. .. ..$ : NULL
#  $ Germline_LogR       :'data.frame':   675960 obs. of  1 variable:
#   ..$ CB1001-11-D01-X1: num [1:675960] 0 0 0 0 0 0 0 0 0 0 ...
#  $ Germline_BAF        :'data.frame':   675960 obs. of  1 variable:
#   ..$ BAF: num [1:675960] 0 0 0 0 0 0 0 0 0 0 ...
#  $ SNPpos              :'data.frame':   675960 obs. of  2 variables:
#   ..$ chrs: chr [1:675960] "1" "1" "1" "1" ...
#   ..$ pos : int [1:675960] 13116 13118 13273 13284 13380 138781 138817 138829 139060 139294 ...
#  $ ch                  :List of 22
#   ..$ : int [1:66137] 1 2 3 4 5 6 7 8 9 10 ...
#   ..$ : int [1:51827] 66138 66139 66140 66141 66142 66143 66144 66145 66146 66147 ...
#   ..$ : int [1:8266] 652529 652530 652531 652532 652533 652534 652535 652536 652537 652538 ...
#   ..$ : int [1:15166] 660795 660796 660797 660798 660799 660800 660801 660802 660803 660804 ...
#  $ chr                 :List of 28
#   ..$ : int [1:34703] 1 2 3 4 5 6 7 8 9 10 ...
#   ..$ : int [1:31434] 34704 34705 34706 34707 34708 34709 34710 34711 34712 34713 ...
#   ..$ : int [1:8266] 652529 652530 652531 652532 652533 652534 652535 652536 652537 652538 ...
#   ..$ : int [1:15166] 660795 660796 660797 660798 660799 660800 660801 660802 660803 660804 ...
#  $ chrs                : chr [1:22] "1" "2" "3" "4" ...
#  $ samples             : chr "CB1001-11-D01-X1"
#  $ gender              : chr "XX"
#  $ sexchromosomes      : chr [1:2] "X" "Y"
#  $ failedarrays        : NULL
snps_to_ascat_bc <- function(snps, segs, sample) {

  ascat_bc <- list(Tumor_LogR = data.frame(sample = snps$LogR),
                   Tumor_BAF = data.frame(sample = snps$BAF),
                   Tumor_LogR_segmented = snps$LogR_segmented,
                   Germline_LogR = NA,
                   Germline_BAF = NA,
                   SNPpos = data.frame(chrs = snps$chrom,
                                       pos = snps$pos),
                   ch = NA,
                   chr = NA,
                   chrs = unique(sort(snps$chrom)),
                   samples = sample,
                   gender = "XX",
                   sexchromosomes = c("X", "Y"),
                   failedarrays = NULL)

  ascat_bc$ch <- lapply(ascat_bc$chrs, function(x) { which(ascat_bc$SNPpos$chrs == x) })
  ascat_bc$chr <- split_genome(bc$SNPpos)  # From ascat.R

  #colnames(ascat_bc$Germline_BAF) <- "BAF"
  colnames(ascat_bc$Tumor_BAF) <- "BAF"

  # Validate our created object
  validate_ascat_bc_input(ascat_bc)

  return(ascat_bc)
}
