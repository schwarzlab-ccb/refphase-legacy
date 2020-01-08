
# addOriginalCPNValuesAndCallsVsPloidyToOutputSegmentation <- function(sample.vec,
#                                        seg.name.vec,
#                                        processed.orig.cpn.segs.df,
#                                        all.segs.df) {
#
#    min.cons.segs.raw.cpnA.df <- generateMinimumConsistentSegmentation(processed.orig.cpn.segs.df, chroms = 1:22, col.to.use = "nAraw")
#    min.cons.segs.raw.cpnB.df <- generateMinimumConsistentSegmentation(processed.orig.cpn.segs.df, chroms = 1:22, col.to.use = "nBraw")
#    min.cons.segs.nMinor.df <- generateMinimumConsistentSegmentation(processed.orig.cpn.segs.df, chroms = 1:22, col.to.use = "nMinor")
#    min.cons.segs.nMajor.df <- generateMinimumConsistentSegmentation(processed.orig.cpn.segs.df, chroms = 1:22, col.to.use = "nMajor")
#    #print("nrow(min.cons.segs.raw.cpnA.df)")
#    #print(nrow(min.cons.segs.raw.cpnA.df))
#    rownames(min.cons.segs.raw.cpnA.df) <- paste(min.cons.segs.raw.cpnA.df[,1], min.cons.segs.raw.cpnA.df[,2], min.cons.segs.raw.cpnA.df[,3], sep = ":")
#    rownames(min.cons.segs.raw.cpnB.df) <- paste(min.cons.segs.raw.cpnA.df[,1], min.cons.segs.raw.cpnA.df[,2], min.cons.segs.raw.cpnA.df[,3], sep = ":")
#    rownames(min.cons.segs.nMinor.df) <- paste(min.cons.segs.raw.cpnA.df[,1], min.cons.segs.raw.cpnA.df[,2], min.cons.segs.raw.cpnA.df[,3], sep = ":")
#    rownames(min.cons.segs.nMajor.df) <- paste(min.cons.segs.raw.cpnA.df[,1], min.cons.segs.raw.cpnA.df[,2], min.cons.segs.raw.cpnA.df[,3], sep = ":")
#
#    for ( samp in sample.vec ) {
#
#        for ( seg.name in seg.name.vec ) {
#            #print(seg.name)
#            if (length(sample.vec) > 1) {
#
#                samp.seg.name <- paste0(samp, ":", seg.name)
#
#            } else {
#
#                samp.seg.name <- seg.name
#            }
#            #print(samp.seg.name)
#            all.segs.df[samp.seg.name,]$orig_cpnA_raw <- min.cons.segs.raw.cpnB.df[seg.name, samp]
#            all.segs.df[samp.seg.name,]$orig_cpnB_raw <- min.cons.segs.raw.cpnA.df[seg.name, samp]
#            all.segs.df[samp.seg.name,]$orig_nMinor <- min.cons.segs.nMinor.df[seg.name, samp]
#            all.segs.df[samp.seg.name,]$orig_nMajor <- min.cons.segs.nMajor.df[seg.name, samp]
#            all.segs.df[samp.seg.name,]$orig_cnTotal <- all.segs.df[samp.seg.name,]$orig_nMinor + all.segs.df[samp.seg.name,]$orig_nMajor
#
#        }#end for ( seg.name in seg.name.vec ) {
#
#    }#end for ( samp in sample.vec ) {
#
#    return(all.segs.df)
#
# }#end addOriginalRawValuesToOutputSegmentation


# Todo refactor and remove this function.... can probably just be done in a for loop in main after the addition of the CCF values..
addPhasingInformationToOutputSegmentation <- function(sample.vec,
                                                      seg.name.vec,
                                                      sample.to.phase.seg.list,
                                                      all.pw.cn.het.comp.lists = NULL,
                                                      all.segs.df) {

  #("Segment names: addPhasingInformationToOutputSegmentation")
  #print(seg.name.vec)

  # Todo, refactor to make this function work on the basis of a single segment:
  for (samp in sample.vec) {
    for (seg.name in seg.name.vec) {
      if (length(sample.vec) > 1) {
        samp.seg.name <- paste0(samp, ":", seg.name)
      } else {
        samp.seg.name <- seg.name
      }

      if (length(sample.to.phase.seg.list[[seg.name]][["samp_chosen_to_phase_seg"]]) == 0) {
        all.segs.df[samp.seg.name, ]$samp_used_to_phase <- "too_few_SNPs"
        #message("1) ", samp.seg.name, ": too_few_SNPs")
      } else {
        samp.used.to.phase <- sample.to.phase.seg.list[[seg.name]][["samp_chosen_to_phase_seg"]]

        all.segs.df[samp.seg.name, ]$samp_used_to_phase <- samp.used.to.phase

        # Check here if the comparison has been done in the reverse direction... not sure it has
        switch.comp.to.ref.vec <- c(
          all.pw.cn.het.comp.lists[[paste0(samp, ":", samp.used.to.phase)]][[seg.name]][["comparison_output_list"]]$raw_cpn_MSAI,
          all.pw.cn.het.comp.lists[[paste0(samp.used.to.phase, ":", samp)]][[seg.name]][["comparison_output_list"]]$raw_cpn_MSAI
        )
        switch.comp.to.ref <- switch.comp.to.ref.vec[!sapply(switch.comp.to.ref.vec, is.null)]
        # If a sample was actually chose as reference it is
        if (is.null(switch.comp.to.ref)) {
          all.segs.df[samp.seg.name, ]$raw_cpn_MSAI <- FALSE
        } else if (is.na(switch.comp.to.ref)) {
          all.segs.df[samp.seg.name, ]$raw_cpn_MSAI <- "missing"
        } else if (switch.comp.to.ref == 1) {
          all.segs.df[samp.seg.name, ]$raw_cpn_MSAI <- TRUE
        } else if (switch.comp.to.ref == 0) {
          all.segs.df[samp.seg.name, ]$raw_cpn_MSAI <- FALSE
        }
      }

      # Need to add the checkf for when there couldn't be any phasng?, for now rely on the the fact it's already there
    }
  }

  all.segs.df
}


addNewRawAndIntCPNCallsToFinalSegmentationParallel <- function(seg.name.vec,
                                                               sample.vec,
                                                               all.samp.logr.df.list,
                                                               all.samp.ph.BAF.df.list,
                                                               asas.lr.baf.lists,
                                                               all.seg.CCF.vals.list,
                                                               all.seg.df,
                                                               cpnX.vec = c("cpnA", "cpnB")) {

  # Make vectors for sample vec:
  all.samp.seg.name.vec <- c()
  sample.rep.vec <- c()
  for (samp in sample.vec) {
    all.samp.seg.name.vec <- c(all.samp.seg.name.vec, seg.name.vec)
    sample.rep.vec <- c(sample.rep.vec, rep(samp, (length(seg.name.vec))))
  }
  iterable.samp.seg.df <- cbind(sample = sample.rep.vec, seg_name = all.samp.seg.name.vec, stringsAsFactors = FALSE)
  itx <- iter(iterable.samp.seg.df, by = "row")
  logdebug("Starting addNewRawAndIntCPNCallsToFinalSegmentationParallel:")
  # Create class which holds multiple results for each loop iteration.
  # Each loop iteration populates two properties: $result1 and $result2.
  # For a great tutorial on S3 classes, see:
  # http://www.cyclismo.org/tutorial/R/s3Classes.html#creating-an-s3-class
  MultiResultClass <- function() {
    me <- list(
      nBraw = NA,
      nAraw = NA,
      nMinor = NA,
      nMajor = NA,
      cnTotal = NA,
      n.het = NA,
      num_ph_hz_sites = NA,
      ph_raw_cpnA = NA,
      ph_raw_cpnB = NA,
      ph_int_cpnA = NA,
      ph_int_cpnB = NA,
      ph_hz_total_raw_cpn = NA,
      min_int_cpnX = NA,
      max_int_cpnX = NA,
      ph_hz_total_cpn = NA,
      non_hz_total_cpn = NA,
      median_BAF = NA,
      mean_BAF = NA,
      stdev_BAF = NA,
      median_hz_logr = NA,
      mean_hz_logr = NA,
      stdev_hz_logr = NA,
      num_un_ph_hz_sites = NA,
      num_non_hz_sites = NA,
      median_non_hz_logr = NA,
      mean_non_hz_logr = NA,
      stdev_non_hz_logr = NA,
      median_comb_logr = NA,
      mean_comb_logr = NA,
      stdev_comb_logr = NA,
      non_hz_total_raw_cpn = NA,
      total_raw_cpn = NA,
      total_int_cpn = NA,
      non_hz_total_raw_cpn = NA,
      cnTotal = NA,
      total_raw_cpn = NA
    )
    class(me) <- append(class(me), "MultiResultClass")
    return(me)
  }

  new.raw.int.cpn.calls.list <- foreach(samp.seg.df = itx, .combine = c) %dopar% {
    ret.mrc <- MultiResultClass()
    samp <- samp.seg.df[, 1]
    seg.name <- samp.seg.df[, 2]
    cpnX <- "cpnA"
    # Fill in flag that shows whether segment has been phasing and CCF valculuation process, ie has greater than min num of hz snp positions
    if (length(all.seg.CCF.vals.list[[cpnX]][[seg.name]]) == 1) {
      # Todo remove this NA signalling that a seg has fewer than the minimum number of hz snps, replace with a field that we can check...
      if (is.na(all.seg.CCF.vals.list[[cpnX]][[seg.name]])) {
        seg.has.fewer.than.min.hz.flag <- TRUE
      } else {
        seg.has.fewer.than.min.hz.flag <- FALSE
      }
    } else {
      seg.has.fewer.than.min.hz.flag <- FALSE
    }

    chr.start.end.vec <- unlist(strsplit(seg.name, split = ":"))
    # if the segment passes the minimum number of hz snps it will have more than one sample at  level
    if (seg.has.fewer.than.min.hz.flag == FALSE) {
      seg.CCF.val.list <- all.seg.CCF.vals.list[[cpnX]][[seg.name]][[samp]]
      # These are currently assigned twice... wiping the original.. definitely needs refactoring:
      ph.cpnA.vec <- seg.CCF.val.list[[paste0(cpnX, "_orig_ph_obs_cpnA_vec")]]
      ph.cpnB.vec <- seg.CCF.val.list[[paste0(cpnX, "_orig_ph_obs_cpnB_vec")]]

      ret.mrc$nBraw <- min(ifelse(median(as.numeric(ph.cpnB.vec)) < 0, 0, median(as.numeric(ph.cpnB.vec))), ifelse(median(as.numeric(ph.cpnA.vec)) < 0, 0, median(as.numeric(ph.cpnA.vec))))
      ret.mrc$nAraw <- max(ifelse(median(as.numeric(ph.cpnB.vec)) < 0, 0, median(as.numeric(ph.cpnB.vec))), ifelse(median(as.numeric(ph.cpnA.vec)) < 0, 0, median(as.numeric(ph.cpnA.vec))))
      ret.mrc$nMinor <- min(round(as.numeric(ret.mrc$nBraw)), round(as.numeric(ret.mrc$nAraw)))
      ret.mrc$nMajor <- max(round(as.numeric(ret.mrc$nBraw)), round(as.numeric(ret.mrc$nAraw)))
      ret.mrc$cnTotal <- as.numeric(ret.mrc$nMinor) + as.numeric(ret.mrc$nMajor)
      ret.mrc$n.het <- length(ph.cpnA.vec)
      ret.mrc$num_ph_hz_sites <- length(ph.cpnA.vec)
      ret.mrc$ph_raw_cpnA <- ifelse(median(as.numeric(ph.cpnA.vec)) < 0, 0, median(as.numeric(ph.cpnA.vec)))
      ret.mrc$ph_raw_cpnB <- ifelse(median(as.numeric(ph.cpnB.vec)) < 0, 0, median(as.numeric(ph.cpnB.vec)))
      ret.mrc$ph_int_cpnA <- round(as.numeric(ret.mrc$ph_raw_cpnA))
      ret.mrc$ph_int_cpnB <- round(as.numeric(ret.mrc$ph_raw_cpnB))
      ret.mrc$ph_hz_total_raw_cpn <- as.numeric(ret.mrc$ph_raw_cpnA) + as.numeric(ret.mrc$ph_raw_cpnB)
      ret.mrc$min_int_cpnX <- min(ret.mrc$ph_int_cpnA, ret.mrc$ph_int_cpnB)
      ret.mrc$max_int_cpnX <- max(ret.mrc$ph_int_cpnA, ret.mrc$ph_int_cpnB)
      ret.mrc$ph_hz_total_cpn <- round(ret.mrc$ph_hz_total_raw_cpn)
      ret.mrc$non_hz_total_cpn <- round(ret.mrc$non_hz_total_raw_cpn)
    }
    # Load in all available BAF and LogR data to fill in the rest of the fields
    # Todo, not entirely sure how seg.asas.lr.baf.df is contributing here...
    seg.asas.lr.baf.df <- asas.lr.baf.lists[[samp]][[seg.name]]
    datapoints.w.BAF.df <- all.samp.ph.BAF.df.list[[samp]]
    ss.ph.baf.df <- all.samp.ph.BAF.df.list[[samp]]
    ss.baf.lr.baf.df <- datapoints.w.BAF.df[datapoints.w.BAF.df$chrom == as.numeric(chr.start.end.vec[1]) &
      datapoints.w.BAF.df$pos >= as.numeric(chr.start.end.vec[2]) &
      datapoints.w.BAF.df$pos <= as.numeric(chr.start.end.vec[3]), , drop = FALSE]
    datapoints.w.logr.only.df <- all.samp.logr.df.list[[samp]]

    # Fill in the fields based on those sites that re hz:
    ret.mrc$median_BAF <- ss.baf.lr.baf.df$ph_median_segmented_BAF[1]
    ret.mrc$mean_BAF <- mean(abs(0.5 - ss.baf.lr.baf.df$BAF)) # Todo, check this was the correct change to mirrored
    ret.mrc$stdev_BAF <- sd(abs(0.5 - ss.baf.lr.baf.df$BAF)) # Todo, check this was the correct change to mirrored
    ret.mrc$median_hz_logr <- median(ss.baf.lr.baf.df$logR)
    ret.mrc$mean_hz_logr <- mean(ss.baf.lr.baf.df$logR)
    ret.mrc$stdev_hz_logr <- sd(ss.baf.lr.baf.df$logR)

    if (seg.has.fewer.than.min.hz.flag == TRUE) {
      ret.mrc$num_un_ph_hz_sites <- nrow(ss.baf.lr.baf.df)
    }

    # Fields that only depend on the non_hz_logr data
    if (nrow(datapoints.w.logr.only.df) > 0) {
      ss.dp.w.logr.only.df <- datapoints.w.logr.only.df[datapoints.w.logr.only.df$chrom == as.numeric(chr.start.end.vec[1]) &
        datapoints.w.logr.only.df$pos >= as.numeric(chr.start.end.vec[2]) &
        datapoints.w.logr.only.df$pos <= as.numeric(chr.start.end.vec[3]), , drop = FALSE]

      ret.mrc$num_non_hz_sites <- nrow(ss.dp.w.logr.only.df)

      if (nrow(ss.dp.w.logr.only.df) > 0) {
        ret.mrc$median_non_hz_logr <- median(ss.dp.w.logr.only.df$logR)
        ret.mrc$mean_non_hz_logr <- mean(ss.dp.w.logr.only.df$logR)
        ret.mrc$stdev_non_hz_logr <- sd(ss.dp.w.logr.only.df$logR)
        # Cmombined logr measures
        ret.mrc$median_comb_logr <- median(c(ss.baf.lr.baf.df$logR, ss.dp.w.logr.only.df$logR))
        ret.mrc$mean_comb_logr <- mean(c(ss.baf.lr.baf.df$logR, ss.dp.w.logr.only.df$logR))
        ret.mrc$stdev_comb_logr <- sd(c(ss.baf.lr.baf.df$logR, ss.dp.w.logr.only.df$logR))
        ret.mrc$non_hz_total_raw_cpn <- median(ss.dp.w.logr.only.df$tot_raw_un_ph_cpn)
      }

      if (seg.has.fewer.than.min.hz.flag == FALSE) {
        ret.mrc$total_raw_cpn <- median(c(ss.dp.w.logr.only.df$tot_raw_un_ph_cpn, (as.numeric(ph.cpnA.vec) + as.numeric(ph.cpnB.vec))))
        ret.mrc$total_int_cpn <- round(median(c(ss.dp.w.logr.only.df$tot_raw_un_ph_cpn, (as.numeric(ph.cpnA.vec) + as.numeric(ph.cpnB.vec)))))
      } else {
        # Todo, get the unphased SNP's total cpn here... we don't calculate this currently as it's a heterozygous SNP that falls below the threshold
        # for ccf calculation.
        ret.mrc$total_raw_cpn <- median(ss.dp.w.logr.only.df$tot_raw_un_ph_cpn)
        ret.mrc$total_int_cpn <- round(median(ss.dp.w.logr.only.df$tot_raw_un_ph_cpn))
      }
    } else if (nrow(datapoints.w.logr.only.df) == 0 & seg.has.fewer.than.min.hz.flag == TRUE) {


      # Nothing we can add here for now: as there's no way of accessing unphased hz snps
      # if( nrow(ss.baf.lr.baf.df) > 0 ) {

      # We could check this obkect for total copy number but... we haven't calculated unph total copy number in the object
      # and it would be irritating to have to call that function here. Sort out at a later date, likely  a v minor subset...

      # }
      # ret.mrc$total_raw_cpn <- median(as.numeric(ph.cpnA.vec) + median(as.numeric(ph.cpnB.vec)))
      # ret.mrc$total_int_cpn <- round(median(as.numeric(ph.cpnA.vec) + as.numeric(ph.cpnB.vec)))
    } else if (nrow(datapoints.w.logr.only.df) == 0 & seg.has.fewer.than.min.hz.flag == FALSE) {
      ret.mrc$total_raw_cpn <- median(as.numeric(ph.cpnA.vec) + median(as.numeric(ph.cpnB.vec)))
      ret.mrc$total_int_cpn <- round(median(as.numeric(ph.cpnA.vec) + as.numeric(ph.cpnB.vec)))
    }

    # Now let's calculate the int toal copy number that we're going to use:
    # Let's set the total copy number

    if (nrow(datapoints.w.logr.only.df) > 0 & seg.has.fewer.than.min.hz.flag == TRUE) {
      ret.mrc$non_hz_total_raw_cpn <- median(ss.dp.w.logr.only.df$tot_raw_un_ph_cpn)
      ret.mrc$cnTotal <- round(median(ss.dp.w.logr.only.df$tot_raw_un_ph_cpn))
    } else if (nrow(datapoints.w.logr.only.df) == 0 & seg.has.fewer.than.min.hz.flag == TRUE) {

      # Note this doesn't guarantee that there will be phased cpn available. Rework at a later date
      # ret.mrc$total_raw_cpn <- median(c(as.numeric(ph.cpnA.vec) + as.numeric(ph.cpnB.vec)))
    }


    returnable.list <- list(ret.mrc)
    names(returnable.list) <- paste0(samp.seg.df[, 1], ":", samp.seg.df[, 2])
    returnable.list
  }

  for (samp in sample.vec) {
    for (seg.name in seg.name.vec) {
      if (length(sample.vec) > 1) {
        samp.seg.name <- paste0(samp, ":", seg.name)
      } else {
        samp.seg.name <- seg.name
      }

      if (!samp.seg.name %in% rownames(all.seg.df)) {
        logwarn("Seg name not in all.seg.df:")
        logwarn(samp.seg.name)
        next
      }

      ret.mrc <- new.raw.int.cpn.calls.list[[samp.seg.name]]

      for (j in 1:length(ret.mrc)) {
        all.seg.df[samp.seg.name, names(ret.mrc)[j]] <- ret.mrc[names(ret.mrc)[j]]
      }
    }
  }

  all.seg.df
}


addNewRawAndIntCPNCallsToFinalSegmentation <- function(seg.name.vec,
                                                       sample.vec,
                                                       all.samp.logr.df.list,
                                                       all.samp.ph.BAF.df.list,
                                                       asas.lr.baf.lists,
                                                       all.seg.CCF.vals.list,
                                                       all.seg.df,
                                                       cpnX.vec = c("cpnA", "cpnB")) {

  for (samp in sample.vec) {
    for (seg.name in seg.name.vec) {
      for (cpnX in cpnX.vec) {
        if (length(sample.vec) > 1) {
          samp.seg.name <- paste0(samp, ":", seg.name)
        } else {
          samp.seg.name <- seg.name
        }

        if (!samp.seg.name %in% rownames(all.seg.df)) {
          logwarn("Seg name not in all.seg.df:")
          logwarn(samp.seg.name)
          next
        }
        # Fill in flag that shows whether segment has been phasing and CCF valculuation process, ie has greater than min num of hz snp positions
        if (length(all.seg.CCF.list[[cpnX]][[seg.name]]) == 1) {
          # Todo remove this NA signalling that a seg has fewer than the minimum number of hz snps, replace with a field that we can check...
          if (is.na(all.seg.CCF.list[[cpnX]][[seg.name]])) {
            seg.has.fewer.than.min.hz.flag <- TRUE
          } else {
            seg.has.fewer.than.min.hz.flag <- FALSE
          }
        } else {
          seg.has.fewer.than.min.hz.flag <- FALSE
        }

        chr.start.end.vec <- unlist(strsplit(seg.name, split = ":"))
        # if the segment passes the minimum number of hz snps it will have more than one sample at  level
        if (seg.has.fewer.than.min.hz.flag == FALSE) {
          seg.CCF.val.list <- all.seg.CCF.vals.list[[cpnX]][[seg.name]][[samp]]
          # These are currently assigned twice... wiping the original.. definitely needs refactoring:
          ph.cpnA.vec <- seg.CCF.val.list[[paste0(cpnX, "_orig_ph_obs_cpnA_vec")]]
          ph.cpnB.vec <- seg.CCF.val.list[[paste0(cpnX, "_orig_ph_obs_cpnB_vec")]]

          all.seg.df[samp.seg.name, ]$nBraw <- min(ifelse(median(as.numeric(ph.cpnB.vec)) < 0, 0, median(as.numeric(ph.cpnB.vec))), ifelse(median(as.numeric(ph.cpnA.vec)) < 0, 0, median(as.numeric(ph.cpnA.vec))))
          all.seg.df[samp.seg.name, ]$nAraw <- max(ifelse(median(as.numeric(ph.cpnB.vec)) < 0, 0, median(as.numeric(ph.cpnB.vec))), ifelse(median(as.numeric(ph.cpnA.vec)) < 0, 0, median(as.numeric(ph.cpnA.vec))))
          all.seg.df[samp.seg.name, ]$nMinor <- min(round(as.numeric(all.seg.df[samp.seg.name, ]$nBraw)), round(as.numeric(all.seg.df[samp.seg.name, ]$nAraw)))
          all.seg.df[samp.seg.name, ]$nMajor <- max(round(as.numeric(all.seg.df[samp.seg.name, ]$nBraw)), round(as.numeric(all.seg.df[samp.seg.name, ]$nAraw)))
          all.seg.df[samp.seg.name, ]$cnTotal <- as.numeric(all.seg.df[samp.seg.name, ]$nMinor) + as.numeric(all.seg.df[samp.seg.name, ]$nMajor)
          all.seg.df[samp.seg.name, ]$n.het <- length(ph.cpnA.vec)
          all.seg.df[samp.seg.name, ]$num_ph_hz_sites <- length(ph.cpnA.vec)
          all.seg.df[samp.seg.name, ]$ph_raw_cpnA <- ifelse(median(as.numeric(ph.cpnA.vec)) < 0, 0, median(as.numeric(ph.cpnA.vec)))
          all.seg.df[samp.seg.name, ]$ph_raw_cpnB <- ifelse(median(as.numeric(ph.cpnB.vec)) < 0, 0, median(as.numeric(ph.cpnB.vec)))
          all.seg.df[samp.seg.name, ]$ph_int_cpnA <- round(as.numeric(all.seg.df[samp.seg.name, ]$ph_raw_cpnA))
          all.seg.df[samp.seg.name, ]$ph_int_cpnB <- round(as.numeric(all.seg.df[samp.seg.name, ]$ph_raw_cpnB))
          all.seg.df[samp.seg.name, ]$ph_hz_total_raw_cpn <- as.numeric(all.seg.df[samp.seg.name, ]$ph_raw_cpnA) + as.numeric(all.seg.df[samp.seg.name, ]$ph_raw_cpnB)
          all.seg.df[samp.seg.name, ]$min_int_cpnX <- min(all.seg.df[samp.seg.name, ]$ph_int_cpnA, all.seg.df[samp.seg.name, ]$ph_int_cpnB)
          all.seg.df[samp.seg.name, ]$max_int_cpnX <- max(all.seg.df[samp.seg.name, ]$ph_int_cpnA, all.seg.df[samp.seg.name, ]$ph_int_cpnB)

          all.seg.df[samp.seg.name, ]$ph_hz_total_cpn <- round(all.seg.df[samp.seg.name, ]$ph_hz_total_raw_cpn)
          all.seg.df[samp.seg.name, ]$non_hz_total_cpn <- round(all.seg.df[samp.seg.name, ]$non_hz_total_raw_cpn)
        }
        # Load in all available BAF and LogR data to fill in the rest of the fields
        # Todo, not entirely sure how seg.asas.lr.baf.df is contributing here...
        seg.asas.lr.baf.df <- asas.lr.baf.lists[[samp]][[seg.name]]
        datapoints.w.BAF.df <- all.samp.ph.BAF.df.list[[samp]]
        ss.ph.baf.df <- all.samp.ph.BAF.df.list[[samp]]
        ss.baf.lr.baf.df <- datapoints.w.BAF.df[datapoints.w.BAF.df$chrom == as.numeric(chr.start.end.vec[1]) &
          datapoints.w.BAF.df$pos >= as.numeric(chr.start.end.vec[2]) &
          datapoints.w.BAF.df$pos <= as.numeric(chr.start.end.vec[3]), , drop = FALSE]
        datapoints.w.logr.only.df <- all.samp.logr.df.list[[samp]]

        all.seg.df[samp.seg.name, ]$median_BAF <- ss.baf.lr.baf.df$ph_median_segmented_BAF[1]
        all.seg.df[samp.seg.name, ]$mean_BAF <- mean(abs(0.5 - ss.baf.lr.baf.df$BAF)) # Todo, check this was the correct change to mirrored
        all.seg.df[samp.seg.name, ]$stdev_BAF <- sd(abs(0.5 - ss.baf.lr.baf.df$BAF)) # Todo, check this was the correct change to mirrored
        all.seg.df[samp.seg.name, ]$median_hz_logr <- median(ss.baf.lr.baf.df$logR)
        all.seg.df[samp.seg.name, ]$mean_hz_logr <- mean(ss.baf.lr.baf.df$logR)
        all.seg.df[samp.seg.name, ]$stdev_hz_logr <- sd(ss.baf.lr.baf.df$logR)

        if (isTRUE(seg.has.fewer.than.min.hz.flag)) {
          all.seg.df[samp.seg.name, ]$num_un_ph_hz_sites <- nrow(ss.baf.lr.baf.df)
        }

        # Fields that only depend on the non_hz_logr data
        if (nrow(datapoints.w.logr.only.df) > 0) {
          ss.dp.w.logr.only.df <- datapoints.w.logr.only.df[datapoints.w.logr.only.df$chrom == as.numeric(chr.start.end.vec[1]) &
            datapoints.w.logr.only.df$pos >= as.numeric(chr.start.end.vec[2]) &
            datapoints.w.logr.only.df$pos <= as.numeric(chr.start.end.vec[3]), , drop = FALSE]

          all.seg.df[samp.seg.name, ]$num_non_hz_sites <- nrow(ss.dp.w.logr.only.df)

          if (nrow(ss.dp.w.logr.only.df) > 0) {
            all.seg.df[samp.seg.name, ]$median_non_hz_logr <- median(ss.dp.w.logr.only.df$logR)
            all.seg.df[samp.seg.name, ]$mean_non_hz_logr <- mean(ss.dp.w.logr.only.df$logR)
            all.seg.df[samp.seg.name, ]$stdev_non_hz_logr <- sd(ss.dp.w.logr.only.df$logR)
            # Cmombined logr measures
            all.seg.df[samp.seg.name, ]$median_comb_logr <- median(c(ss.baf.lr.baf.df$logR, ss.dp.w.logr.only.df$logR))
            all.seg.df[samp.seg.name, ]$mean_comb_logr <- mean(c(ss.baf.lr.baf.df$logR, ss.dp.w.logr.only.df$logR))
            all.seg.df[samp.seg.name, ]$stdev_comb_logr <- sd(c(ss.baf.lr.baf.df$logR, ss.dp.w.logr.only.df$logR))
            all.seg.df[samp.seg.name, ]$non_hz_total_raw_cpn <- median(ss.dp.w.logr.only.df$tot_raw_un_ph_cpn)
          }

          if (isFALSE(seg.has.fewer.than.min.hz.flag)) {
            all.seg.df[samp.seg.name, ]$total_raw_cpn <- median(c(ss.dp.w.logr.only.df$tot_raw_un_ph_cpn, (as.numeric(ph.cpnA.vec) + as.numeric(ph.cpnB.vec))))
            all.seg.df[samp.seg.name, ]$total_int_cpn <- round(median(c(ss.dp.w.logr.only.df$tot_raw_un_ph_cpn, (as.numeric(ph.cpnA.vec) + as.numeric(ph.cpnB.vec)))))
          } else {
            # Todo, get the unphased SNP's total cpn here... we don't calculate this currently as it's a heterozygous SNP that falls below the threshold
            # for ccf calculation.
            all.seg.df[samp.seg.name, ]$total_raw_cpn <- median(ss.dp.w.logr.only.df$tot_raw_un_ph_cpn)
            all.seg.df[samp.seg.name, ]$total_int_cpn <- round(median(ss.dp.w.logr.only.df$tot_raw_un_ph_cpn))
          }
        } else if (nrow(datapoints.w.logr.only.df) == 0 & isTRUE(seg.has.fewer.than.min.hz.flag)) { # end  if(nrow(datapoints.w.logr.only.df) > 0) {

          # all.seg.df[samp.seg.name,]$total_raw_cpn <- median(as.numeric(ph.cpnA.vec) + median(as.numeric(ph.cpnB.vec)))
          # all.seg.df[samp.seg.name,]$total_int_cpn <- round(median(as.numeric(ph.cpnA.vec) + as.numeric(ph.cpnB.vec)))
        } else if (nrow(datapoints.w.logr.only.df) == 0 & isFALSE(seg.has.fewer.than.min.hz.flag)) {
          all.seg.df[samp.seg.name, ]$total_raw_cpn <- median(as.numeric(ph.cpnA.vec) + median(as.numeric(ph.cpnB.vec)))
          all.seg.df[samp.seg.name, ]$total_int_cpn <- round(median(as.numeric(ph.cpnA.vec) + as.numeric(ph.cpnB.vec)))
        }

        # Now let's calculate the int toal copy number that we're going to use:
        # Let's set the total copy number
        if (nrow(datapoints.w.logr.only.df) > 0 & seg.has.fewer.than.min.hz.flag == TRUE) {
          all.seg.df[samp.seg.name, ]$non_hz_total_raw_cpn <- median(ss.dp.w.logr.only.df$tot_raw_un_ph_cpn)
          all.seg.df[samp.seg.name, ]$cnTotal <- round(median(ss.dp.w.logr.only.df$tot_raw_un_ph_cpn))
        } else if (nrow(datapoints.w.logr.only.df) == 0 & seg.has.fewer.than.min.hz.flag == FALSE) {
          all.seg.df[samp.seg.name, ]$total_raw_cpn <- median(c(as.numeric(ph.cpnA.vec) + as.numeric(ph.cpnB.vec)))
        }
      }
    }
  }

  all.seg.df
}



retestSubclonalAmps <- function(all.seg.df, all.samp.ph.BAF.df.list, all.samp.logr.df.list, gain.threshold, amp.threshold, p.val.for.converting.gain.to.amp = 0.01) {
  non.na.cpn.event.df <- all.seg.df[!is.na(all.seg.df$cpn_event_vs_ploidy), ]
  amp.all.seg.df <- non.na.cpn.event.df[non.na.cpn.event.df$cpn_event_vs_ploidy == "amp", ]

  if (nrow(amp.all.seg.df) > 0) {
    # So now we want to go check if there are any sample that has a gain called in the same seg
    for (amp.idx in 1:nrow(amp.all.seg.df)) {
      current.amp.samp <- amp.all.seg.df[amp.idx, ]$sample_name
      current.amp.seg.name <- amp.all.seg.df[amp.idx, ]$seg_name
      gain.segs.matching.amp.df <- all.seg.df[all.seg.df$seg_name == current.amp.seg.name & all.seg.df$cpn_event_vs_ploidy == "gain", ]

      if (nrow(gain.segs.matching.amp.df) > 0) {

        # Now we now there are gains in other samples let's load in the amp logr data to enable a comparison:
        c.a.chr.start.end.vec <- unlist(strsplit(current.amp.seg.name, split = ":"))

        c.a.samp.datapoints.w.BAF.df <- all.samp.ph.BAF.df.list[[current.amp.samp]]
        c.a.baf.lr.baf.df <- c.a.samp.datapoints.w.BAF.df[c.a.samp.datapoints.w.BAF.df$chrom == as.numeric(c.a.chr.start.end.vec[1]) &
          c.a.samp.datapoints.w.BAF.df$pos >= as.numeric(c.a.chr.start.end.vec[2]) &
          c.a.samp.datapoints.w.BAF.df$pos <= as.numeric(c.a.chr.start.end.vec[3]), , drop = FALSE]
        c.a.samp.datapoints.w.logr.only.df <- all.samp.logr.df.list[[current.amp.samp]]

        if (nrow(c.a.samp.datapoints.w.logr.only.df) > 0) {
          c.a.dp.w.logr.only.df <- c.a.samp.datapoints.w.logr.only.df[c.a.samp.datapoints.w.logr.only.df$chrom == as.numeric(chr.start.end.vec[1]) &
            c.a.samp.datapoints.w.logr.only.df$pos >= as.numeric(chr.start.end.vec[2]) &
            c.a.samp.datapoints.w.logr.only.df$pos <= as.numeric(chr.start.end.vec[3]), , drop = FALSE]
          if (nrow(c.a.dp.w.logr.only.df) > 0) {
            c.a.comb.baf.lr.baf.df <- rbind(
              c.a.baf.lr.baf.df[, c("chrom", "pos", "BAF", "BAFsegmented", "logR", "logRsegmented")],
              c.a.dp.w.logr.only.df[, c("chrom", "pos", "BAF", "BAFsegmented", "logR", "logRsegmented")]
            )
          } else {
            c.a.comb.baf.lr.baf.df <- c.a.baf.lr.baf.df
          }
        }

        # now let's check each of these matching gains vs the amp used to select them
        for (gain.idx in 1:nrow(gain.segs.matching.amp.df)) {
          current.gain.samp <- gain.segs.matching.amp.df[gain.idx, ]$sample_name
          current.gain.seg.name <- gain.segs.matching.amp.df[gain.idx, ]$seg_name
          c.g.chr.start.end.vec <- unlist(strsplit(current.gain.seg.name, split = ":"))
          c.g.samp.datapoints.w.BAF.df <- all.samp.ph.BAF.df.list[[current.gain.samp]]
          c.g.baf.lr.baf.df <- c.g.samp.datapoints.w.BAF.df[c.g.samp.datapoints.w.BAF.df$chrom == as.numeric(c.g.chr.start.end.vec[1]) &
            c.g.samp.datapoints.w.BAF.df$pos >= as.numeric(c.g.chr.start.end.vec[2]) &
            c.g.samp.datapoints.w.BAF.df$pos <= as.numeric(c.g.chr.start.end.vec[3]), , drop = FALSE]
          c.g.samp.datapoints.w.logr.only.df <- all.samp.logr.df.list[[current.gain.samp]]

          if (nrow(c.g.samp.datapoints.w.logr.only.df) > 0) {
            c.g.dp.w.logr.only.df <- c.g.samp.datapoints.w.logr.only.df[c.g.samp.datapoints.w.logr.only.df$chrom == as.numeric(chr.start.end.vec[1]) &
              c.g.samp.datapoints.w.logr.only.df$pos >= as.numeric(chr.start.end.vec[2]) &
              c.g.samp.datapoints.w.logr.only.df$pos <= as.numeric(chr.start.end.vec[3]), , drop = FALSE]
            if (nrow(c.g.dp.w.logr.only.df) > 0) {
              c.g.comb.baf.lr.baf.df <- rbind(
                c.g.baf.lr.baf.df[, c("chrom", "pos", "BAF", "BAFsegmented", "logR", "logRsegmented")],
                c.g.dp.w.logr.only.df[, c("chrom", "pos", "BAF", "BAFsegmented", "logR", "logRsegmented")]
              )
            } else {
              ss.comb.baf.lr.baf.df <- ss.baf.lr.baf.df
            }
          }

          # Now test vs the amp logr:
          if (tTestWithErrorHandling(c.a.comb.baf.lr.baf.df$logR, c.g.comb.baf.lr.baf.df$logR)$p.val > p.val.for.converting.gain.to.amp) {
            all.seg.df[paste0(current.gain.samp, ":", current.gain.seg.name), ]$gain_conv_to_amp <- "amp"
          }
        }
      }
    }
  } else {
    logdebug("No amps")
  }

  return(all.seg.df)
}


addOriginalSegmentInfoToFinalSegmentation <- function(seg.name.vec,
                                                      sample.vec,
                                                      processed.orig.cpn.segs.df,
                                                      gain.threshold,
                                                      loss.threshold,
                                                      amp.threshold,
                                                      deep.loss.threshold,
                                                      all.seg.df) {

  # We need to work out which segment this one originally came from in the input.
  # Data we need for this?

  # The original segmentation...

  for (samp in sample.vec) {
    for (seg.name in seg.name.vec) {
      if (length(sample.vec) > 1) {
        samp.seg.name <- paste0(samp, ":", seg.name)
      } else {
        samp.seg.name <- seg.name
      }

      chr.start.end.vec <- unlist(strsplit(seg.name, split = ":"))

      chr.samp.seg.df <- processed.orig.cpn.segs.df[processed.orig.cpn.segs.df$sample == samp
      & as.numeric(processed.orig.cpn.segs.df$chr) == as.numeric(chr.start.end.vec[1]), ,
      drop = FALSE
      ]

      cpn.call.state.vec <- rep("neutral", nrow(chr.samp.seg.df))
      tot.raw.cpn.vec <- chr.samp.seg.df$nAraw + chr.samp.seg.df$nBraw
      log.tot.raw.over.ploidy.cpn.vec <- log2(tot.raw.cpn.vec / chr.samp.seg.df$Ploidy)

      cpn.call.state.vec[log.tot.raw.over.ploidy.cpn.vec >= log2(2.5 / 2)] <- "gain"
      cpn.call.state.vec[log.tot.raw.over.ploidy.cpn.vec <= log2(1.5 / 2)] <- "loss"
      cpn.call.state.vec[tot.raw.cpn.vec >= ((2 * chr.samp.seg.df$Ploidy[1]) + 1)] <- "amp"
      cpn.call.state.vec[tot.raw.cpn.vec <= 0.5] <- "deep_loss"

      chr.samp.seg.df$orig_cpn_call_vs_ploidy <- cpn.call.state.vec


      matching.chr.samp.seg.df <- chr.samp.seg.df[chr.samp.seg.df$startpos == as.numeric(chr.start.end.vec[2])
      | chr.samp.seg.df$endpos == as.numeric(chr.start.end.vec[3])
      | (chr.samp.seg.df$startpos >= as.numeric(chr.start.end.vec[2]) & chr.samp.seg.df$endpos <= as.numeric(chr.start.end.vec[3]))
      | (as.numeric(chr.start.end.vec[2]) >= chr.samp.seg.df$startpos & as.numeric(chr.start.end.vec[3]) <= chr.samp.seg.df$endpos), , drop = FALSE ]

      # add a column defining gain and loss etc...
      if (nrow(matching.chr.samp.seg.df) == 0) {
        logwarn(chr.samp.seg.df)
        logwarn(samp)
        logwarn(seg.name)
        logwarn(matching.chr.samp.seg.df)
        logwarn("Error, minimum consistent region segmentation doesn't map back to any original segments")
        # Which segment does it overlap most with...
        chr.samp.seg.range.df <- chr.samp.seg.df[, c("chr", "startpos", "endpos")]
        colnames(chr.samp.seg.range.df) <- c("Chromosome", "start", "end")
        logdebug("chr.samp.seg.range.df")
        print(chr.samp.seg.range.df)
        ref.gen.range <- makeGRangesFromDataFrame(chr.samp.seg.range.df)
        seg.test.df <- data.frame(
          Chromosome = chr.start.end.vec[1],
          start = as.numeric(chr.start.end.vec[2]),
          end = as.numeric(chr.start.end.vec[3]),
          stringsAsFactors = FALSE
        )
        test.gen.range <- makeGRangesFromDataFrame(seg.test.df)
        #print("test.gen.range")
        #print(test.gen.range)
        # find overlaps:
        hit.obj <- findOverlaps(ref.gen.range, test.gen.range)
        overlaps <- pintersect(ref.gen.range[queryHits(hit.obj)], test.gen.range[subjectHits(hit.obj)])
        overlap.df <- as.data.frame(overlaps)
        percent.overlap.vec <- width(overlaps) / width(test.gen.range[subjectHits(hit.obj)])
        overlap.df$percent_overlap <- as.numeric(percent.overlap.vec)
        # Should probably do percent overlaps



        if (nrow(overlap.df) > 0) {
          overlap.df <- overlap.df[with(overlap.df, order(-width)), ]

          #print("chr.samp.seg.df")
          #print(chr.samp.seg.df)
          #print("overlap.df")
          #print(overlap.df)

          most.overlapping.seg.df <- chr.samp.seg.df[chr.samp.seg.df$chr %in% overlap.df$seqnames &
            chr.samp.seg.df$startpos %in% overlap.df$start &
            chr.samp.seg.df$endpos %in% overlap.df$end, , drop = FALSE]

          all.seg.df[samp.seg.name, ]$orig_seg_name <- paste(most.overlapping.seg.df[1, ]$chr, most.overlapping.seg.df[1, ]$startpos, most.overlapping.seg.df[1, ]$endpos, sep = ":")
          all.seg.df[samp.seg.name, ]$orig_cpn_call_vs_ploidy <- most.overlapping.seg.df[1, ]$orig_cpn_call_vs_ploidy

          if (all.seg.df[samp.seg.name, ]$orig_seg_name == "NA:NA:NA") {
            nearest.obj <- nearest(test.gen.range, ref.gen.range)
            dist <- distanceToNearest(test.gen.range, ref.gen.range, select = "arbitrary")
            nearest.grange <- ref.gen.range[nearest.obj[1], ]

            nearest.gr.df <- as.data.frame(nearest.grange)

            nearest.seg.df <- chr.samp.seg.df[chr.samp.seg.df$chr %in% nearest.gr.df$seqnames &
              chr.samp.seg.df$startpos %in% nearest.gr.df$start &
              chr.samp.seg.df$endpos %in% nearest.gr.df$end, , drop = FALSE]


            all.seg.df[samp.seg.name, ]$orig_seg_name <- paste(nearest.seg.df[1, ]$chr, nearest.seg.df[1, ]$startpos, nearest.seg.df[1, ]$endpos, sep = ":")
            all.seg.df[samp.seg.name, ]$orig_cpn_call_vs_ploidy <- nearest.seg.df[1, ]$orig_cpn_call_vs_ploidy

            # Let's do the nearest thing...
          }
        } else {
          stop("No seg exactly matches and no segment overlaps, error in the minimum consistent segmentation that needs to be investigated.")
        }
      } else if (nrow(matching.chr.samp.seg.df) > 2) {
        print(chr.samp.seg.df)
        print(samp)
        print(seg.name)
        print(matching.chr.samp.seg.df)
        stop("Error, minimum consistent region segmentation maps back to more than two original segments")
      } else if (nrow(matching.chr.samp.seg.df) == 2) {

        # Defualt to larger segment... that contains
        reordered.matching.chr.samp.seg.df <- matching.chr.samp.seg.df[with(matching.chr.samp.seg.df, order(-n.het)), ]
        all.seg.df[samp.seg.name, ]$orig_seg_name <- paste(reordered.matching.chr.samp.seg.df[1, ]$chr, reordered.matching.chr.samp.seg.df[1, ]$startpos, reordered.matching.chr.samp.seg.df[1, ]$endpos, sep = ":")
        all.seg.df[samp.seg.name, ]$orig_cpn_call_vs_ploidy <- reordered.matching.chr.samp.seg.df[1, ]$orig_cpn_call_vs_ploidy
      } else {
        all.seg.df[samp.seg.name, ]$orig_cpn_call_vs_ploidy <- matching.chr.samp.seg.df[1, ]$orig_cpn_call_vs_ploidy
        all.seg.df[samp.seg.name, ]$orig_seg_name <- paste(matching.chr.samp.seg.df$chr, matching.chr.samp.seg.df$startpos, matching.chr.samp.seg.df$endpos, sep = ":")
      }
      # Ok so now let's get the original name from something looking like this:
      #           sample chr  startpos    endpos n.het cnTotal nMajor nMinor   Ploidy  ACF     nAraw       nBraw patient sample_name tot_raw_cpn min_raw_cpn max_raw_cpn cpn_to_use ploidy_adj_tot_cpn
      # 1 B_MET048_MT-M1   1    762496   4834651  4329       2      1      1 2.293931 0.83 0.8747692 0.874769163  MET048   MET048:M1   1.7495383 0.874769163   0.8747692  1.7495383         -0.3908477
    }
  }

  all.seg.df
}



addWithinSampleChangeToCpnAndAllelicImbalanceCalls <- function(all.seg.df, sample.vec, seg.name.vec) {
  for (samp in sample.vec) {
    for (seg.name in seg.name.vec) {
      if (length(sample.vec) > 1) {
        samp.seg.name <- paste0(samp, ":", seg.name)
      } else {
        samp.seg.name <- seg.name
      }
      # In some cases we may not be able to do this as the they are
      # Add where the new nMinor differs
      if (is.na(all.seg.df[samp.seg.name, ]$nMinor) == TRUE) {
        # If couldn't phase then there won't be any change to the nMinmor, should probably change this to
        # print("skipping to next seg.")
        next
      }

      # Changes to copy number states...
      all.seg.df[samp.seg.name, ]$nMinor_change <- all.seg.df[samp.seg.name, ]$nMinor != all.seg.df[samp.seg.name, ]$orig_nMinor
      # Add where the new nMajor differs
      all.seg.df[samp.seg.name, ]$nMajor_change <- all.seg.df[samp.seg.name, ]$nMajor != all.seg.df[samp.seg.name, ]$orig_nMajor
      # Add where the total copy number differs
      all.seg.df[samp.seg.name, ]$cnTotal_change <- all.seg.df[samp.seg.name, ]$cnTotal != all.seg.df[samp.seg.name, ]$orig_cnTotal
      all.seg.df[samp.seg.name, ]$within_samp_cl_AI_present <- (as.numeric(all.seg.df[samp.seg.name, ]$nMinor) != as.numeric(all.seg.df[samp.seg.name, ]$nMajor))
      all.seg.df[samp.seg.name, ]$within_samp_cl_LOH_present <- min(as.numeric(all.seg.df[samp.seg.name, ]$nMinor), as.numeric(all.seg.df[samp.seg.name, ]$nMajor)) == 0

      # Todo, Should this be exor?
      all.seg.df[samp.seg.name, ]$within_samp_sc_AI_present <- (all.seg.df[samp.seg.name, ]$nMinor == all.seg.df[samp.seg.name, ]$nMajor) &
        (all.seg.df[samp.seg.name, ]$within_samp_sc_cpnA_present == TRUE | all.seg.df[samp.seg.name, ]$within_samp_sc_cpnB_present == TRUE)


      all.seg.df[samp.seg.name, ]$within_samp_sc_LOH_present <- (all.seg.df[samp.seg.name, ]$nMinor == all.seg.df[samp.seg.name, ]$nMajor) &
        ((all.seg.df[samp.seg.name, ]$within_samp_sc_cpnA_present == TRUE & all.seg.df[samp.seg.name, ]$psc_cpnA_low == 0 & all.seg.df[samp.seg.name, ]$sc_cpnA_high_or_low_dominate == "LOW") |
          (all.seg.df[samp.seg.name, ]$within_samp_sc_cpnB_present == TRUE & all.seg.df[samp.seg.name, ]$psc_cpnB_low == 0 & all.seg.df[samp.seg.name, ]$sc_cpnB_high_or_low_dominate == "LOW"))

      # Add where new within sample clonal AI has changed
      all.seg.df[samp.seg.name, ]$new_within_samp_cl_AI_detected <- (all.seg.df[samp.seg.name, ]$orig_nMajor == all.seg.df[samp.seg.name, ]$orig_nMinor) &
        (all.seg.df[samp.seg.name, ]$nMinor != all.seg.df[samp.seg.name, ]$nMajor)
      # Add where new within sample clonal LOH has changed
      all.seg.df[samp.seg.name, ]$new_within_samp_cl_LOH_detected <- (min(all.seg.df[samp.seg.name, ]$orig_nMinor, all.seg.df[samp.seg.name, ]$orig_nMajor) != 0) &
        (min(all.seg.df[samp.seg.name, ]$nMinor, all.seg.df[samp.seg.name, ]$nMinor) == 0)
    }
  }

  all.seg.df
}


# Todo, factor in MSAI nto determining whether these events are clonal...? Or do this outside the loop
determineAcrossSampleClonalityOfLogrEventsVsPloidy <- function(sample.vec,
                                                               seg.name.vec,
                                                               all.seg.df,
                                                               all.samp.logr.df.list,
                                                               all.samp.ph.BAF.df.list,
                                                               gain.threshold,
                                                               loss.threshold,
                                                               amp.threshold,
                                                               deep.loss.threshold,
                                                               tot.cpn.vs.ploidy.pval.thresh) {

  # This vector is used to record a summary of the calls vs ploidy if there subclonal instances
  # We can think about adding MSAI to this at a later date, for the moment if MSAI s detected as being recorded
  # We will make sure everything is recorded as subclonal.
  cpn.event.name.vec <- c("deep_loss", "loss", "neutral", "gain", "amp")

  for (seg.name in seg.name.vec) {
    sample.seg.vec <- paste0(sample.vec, ":", seg.name)
    cpn.events.vs.ploidy.vec <- all.seg.df[sample.seg.vec, ]$cpn_event_vs_ploidy
    uniq.cpn.events.vec <- unique(all.seg.df[sample.seg.vec, ]$cpn_event_vs_ploidy)

    if (length(uniq.cpn.events.vec) == 1) {
      if (is.na(uniq.cpn.events.vec[1])) {
        all.seg.df[sample.seg.vec, ]$across_samp_cl_cpn_vs_ploidy_present <- NA
        all.seg.df[sample.seg.vec, ]$across_samp_sc_cpn_vs_ploidy_present <- NA
      } else {
        all.seg.df[sample.seg.vec, ]$across_samp_cl_cpn_vs_ploidy_present <- uniq.cpn.events.vec[1]
        all.seg.df[sample.seg.vec, ]$across_samp_sc_cpn_vs_ploidy_present <- "no_subclonal_events"
      }
      # So what are the colkumn names that we need for this.
    } else { # So now we know that some subclonal events are present, if there is an NA present in the vector now it must be an error

      if (length(uniq.cpn.events.vec[is.na(uniq.cpn.events.vec)]) > 0) {
        logwarn("WARNING: Shouldn't be a mixture of NA and cpn event vs ploidy results for the same segment")
        logwarn("This may occur if one of the segments had logr values that were all the same while others segments logr values varied,this would lead to tTestWithErrorHandling returning NA for the case with constant data and a p.value for the others. If the segment is large, there may be a more serious bug.")
        logwarn(seg.name)
        logwarn(uniq.cpn.events.vec)
        logwarn(all.seg.df[sample.seg.vec, c(
          "num_ph_hz_sites",
          "num_un_ph_hz_sites",
          "num_non_hz_sites",
          "cpn_event_vs_ploidy"
        )])
        logwarn("Setting all instances of this segment's cpn_event_vs_ploidy and associated columns to NA")
        all.seg.df[sample.seg.vec, c(
          "gain_hz_logr_pval",
          "amp_hz_logr_pval",
          "loss_hz_logr_pval",
          "deep_loss_hz_logr_pval",
          "gain_non_hz_logr_pval",
          "amp_non_hz_logr_pval",
          "loss_non_hz_pval",
          "deep_loss_non_hz_pval",
          "gain_comb_logr_pval",
          "amp_comb_logr_pval",
          "loss_comb_logr_pval",
          "deep_loss_comb_logr_pval",
          "cpn_event_vs_ploidy"
        )] <- NA
      } else if (all(is.na(all.seg.df[sample.seg.vec, ]$MSAI)) | all(as.logical(all.seg.df[sample.seg.vec, ]$MSAI)) == FALSE) {

        # So now to make this easy, would be good to have a sterotyped way of reporting the cpn events. e.g. deep_loss:0:loss:0:neutral:0:gain:0:amp:0
        # Do we want to report the samples themselves? no, too much detail can easily look them up...
        cpn.event.count.vec <- rep(0, length(cpn.event.name.vec))
        names(cpn.event.count.vec) <- cpn.event.name.vec
        for (cpn.event in cpn.event.name.vec) {
          cpn.event.count.vec[cpn.event] <- length(cpn.events.vs.ploidy.vec[cpn.events.vs.ploidy.vec == cpn.event])
        } # end for ( cpn.event in cpn.event.name.vec ) {
        # This interleaves the names and counts into a single vector that's collapsed into a single string.
        cpn.event.count.string.df <- cbind(names(cpn.event.count.vec), rep("=", length(cpn.event.count.vec)), cpn.event.count.vec)
        paste(apply(cpn.event.count.string.df, 1, paste, collapse = ""), collapse = ";")

        all.seg.df[sample.seg.vec, ]$across_samp_sc_cpn_vs_ploidy_present <- paste(apply(cpn.event.count.string.df, 1, paste, collapse = ""), collapse = ";")
        all.seg.df[sample.seg.vec, ]$across_samp_cl_cpn_vs_ploidy_present <- "subclonal_events_present"
      }
    }
    # If the first value is NA then all others will be and we should move on
  }

  all.seg.df
}


addWithinSampleClonalityAndCCFToFinalSegmentation <- function(seg.name.vec,
                                                              all.seg.CCF.vals.list,
                                                              all.seg.CCF.diff.from.exp.list,
                                                              all.seg.df,
                                                              sample.vec,
                                                              clonal.pval = 0.001,
                                                              diff.from.oth.p.val = 0.01,
                                                              min.CCF.to.consider.val = 0.15,
                                                              min.cpn.to.consider.val = 0.15,
                                                              min.num.other.chrs = 3,
                                                              chrs.to.use.vec = 1:22,
                                                              cpnX.vec = c("cpnA", "cpnB"),
                                                              reduced.seg.flag = FALSE,
                                                              max.cpn.val = 10) {
  logdebug("Segment names: addWithinSampleClonalityAndCCFToFinalSegmentation")
  #print(seg.name.vec)

  for (samp in sample.vec) {
    for (seg.name in seg.name.vec) {
      if (length(sample.vec) > 1) {
        samp.seg.name <- paste0(samp, ":", seg.name)
      } else {
        samp.seg.name <- seg.name
      }

      if (!seg.name %in% names(all.seg.CCF.vals.list[[cpnX.vec[1]]])) {
        all.seg.df[samp.seg.name, ]$samp_used_to_phase <- "too_few_SNPs"
        #message("2) ", samp.seg.name, ": too_few_SNPs")
        next
      }
      #str(all.seg.CCF.vals.list[[cpnX.vec[1]]])
      if (length(all.seg.CCF.vals.list[[cpnX.vec[1]]][[seg.name]]) == 1) {
        all.seg.df[samp.seg.name, ]$samp_used_to_phase <- "too_few_SNPs"
        #message("3) ", samp.seg.name, ": too_few_SNPs")
        #str(all.seg.CCF.vals.list[[cpnX.vec[1]]][[seg.name]])
        next
      }

      for (cpnX in cpnX.vec) {
        seg.CCF.val.list <- all.seg.CCF.vals.list[[cpnX]][[seg.name]][[samp]]
        seg.diff.from.others.list <- all.seg.CCF.diff.from.exp.list[[cpnX]][[seg.name]][[samp]]
        ph.cpnX.vec <- seg.CCF.val.list[[paste0(cpnX, "_orig_ph_obs_", cpnX, "_vec")]]


        all.seg.df[samp.seg.name, paste0("psc_", cpnX, "_low")] <- median(as.numeric(seg.CCF.val.list[[paste0("poss_sc_", cpnX, "_low")]]))
        all.seg.df[samp.seg.name, paste0("psc_", cpnX, "_high")] <- median(as.numeric(seg.CCF.val.list[[paste0("poss_sc_", cpnX, "_high")]]))

        ss.frac.poss.sc.cpnX.low.vec <- as.numeric(seg.CCF.val.list[[paste0("frac_", cpnX, "_low_vec")]])
        ss.frac.poss.sc.cpnX.high.vec <- as.numeric(seg.CCF.val.list[[paste0("frac_", cpnX, "_high_vec")]])

        all.seg.df[samp.seg.name, paste0("frac_psc_", cpnX, "_low")] <- median(as.numeric(seg.CCF.val.list[[paste0("mean_frac_", cpnX, "_low")]]))
        all.seg.df[samp.seg.name, paste0("frac_psc_", cpnX, "_high")] <- median(as.numeric(seg.CCF.val.list[[paste0("mean_frac_", cpnX, "_high")]]))

        # Test for evidence of subclonality within the segment for cpnA
        mu.ccf.val <- 1 - min.CCF.to.consider.val

        cpnX.mean.and.median.flag <- TRUE
        if (all.seg.df[samp.seg.name, paste0("psc_", cpnX, "_low")] == all.seg.df[samp.seg.name, paste0("psc_", cpnX, "_high")]) {
          # print("balanced cpn")
          all.seg.df[samp.seg.name, paste0("sc_", cpnX, "_high_or_low_dominate")] <- NA
          all.seg.df[samp.seg.name, paste0("sc_", cpnX, "_vs_min_CCF_pval")] <- 1
          all.seg.df[samp.seg.name, paste0("sc_", cpnX, "_vs_min_cpn_pval")] <- 1
        } else if (median(ss.frac.poss.sc.cpnX.low.vec) > median(ss.frac.poss.sc.cpnX.high.vec)) {
          cpnA.test.subclonal <- wilcoxTestWithErrorHandling(ss.frac.poss.sc.cpnX.low.vec, mu = mu.ccf.val, alternative = "less", exact = FALSE)
          all.seg.df[samp.seg.name, paste0("sc_", cpnX, "_high_or_low_dominate")] <- "LOW"
          all.seg.df[samp.seg.name, paste0("sc_", cpnX, "_vs_min_CCF_pval")] <- cpnA.test.subclonal$p.value

          # Now let's try with min.cpn.to.consider.val
          mu.cpn.val <- all.seg.df[samp.seg.name, paste0("psc_", cpnX, "_low")] + min.cpn.to.consider.val
          cpnA.cpn.mu.test.subclonal <- wilcoxTestWithErrorHandling(ph.cpnX.vec, mu = mu.cpn.val, alternative = "greater", exact = FALSE)
          all.seg.df[samp.seg.name, paste0("sc_", cpnX, "_vs_min_cpn_pval")] <- cpnA.cpn.mu.test.subclonal$p.value

          cpnX.mean.and.median.flag <- mean(ph.cpnX.vec) > mu.cpn.val & median(ph.cpnX.vec) > mu.cpn.val
        } else { # end   if ( median(seg.CCF.val.list[[seg.name]][[samp]]$frac_cpnA_low_vec) > median(seg.CCF.val.list[[seg.name]][[samp]]$frac_cpnA_high_vec) ) {
          cpnA.test.subclonal <- wilcoxTestWithErrorHandling(ss.frac.poss.sc.cpnX.high.vec, mu = mu.ccf.val, alternative = "less", exact = FALSE)
          all.seg.df[samp.seg.name, paste0("sc_", cpnX, "_high_or_low_dominate")] <- "HIGH"
          all.seg.df[samp.seg.name, paste0("sc_", cpnX, "_vs_min_CCF_pval")] <- cpnA.test.subclonal$p.value

          mu.cpn.val <- all.seg.df[samp.seg.name, paste0("psc_", cpnX, "_high")] - min.cpn.to.consider.val
          cpnA.cpn.mu.test.subclonal <- wilcoxTestWithErrorHandling(ph.cpnX.vec, mu = mu.cpn.val, alternative = "less", exact = FALSE)
          all.seg.df[samp.seg.name, paste0("sc_", cpnX, "_vs_min_cpn_pval")] <- cpnA.cpn.mu.test.subclonal$p.value

          cpnX.mean.and.median.flag <- mean(ph.cpnX.vec) < mu.cpn.val & median(ph.cpnX.vec) < mu.cpn.val
        }
        # readline(prompt="Press [enter] to continue")
        # We need to refactor this magic number out of this
        cpnX.pval.vs.other.chroms <- 1e-120
        abs.diff.from.oth.match.cpn.chrs.flag <- TRUE
        directionality.vs.oth.cpn.match.chrs.raw.cpn <- TRUE
        all.seg.df[samp.seg.name, paste0("num_other_chrs_same_", cpnX, "_cpn")] <- length(seg.diff.from.others.list[[paste0(cpnX, "_num_chrs_w_match_int_cpn_vec")]])


        if (all.seg.df[samp.seg.name, paste0("num_other_chrs_same_", cpnX, "_cpn")] >= min.num.other.chrs) {
          cpnX.pval.vs.other.chroms <- seg.diff.from.others.list[[paste0(cpnX, "_pval_vs_oth_cpn_match_chrs_raw_cpn")]]
          mean.oth.cpn.match.chrs.raw.cpn <- seg.diff.from.others.list[[paste0(cpnX, "_mean_oth_cpn_match_chrs_raw_cpn")]]
          median.oth.cpn.match.chrs.raw.cpn <- seg.diff.from.others.list[[paste0(cpnX, "_median_oth_cpn_match_chrs_raw_cpn")]]
          abs.diff.from.oth.match.cpn.chrs.flag <- abs(mean(ph.cpnX.vec) - mean.oth.cpn.match.chrs.raw.cpn) >= 0.04 & abs(median(ph.cpnX.vec) - median.oth.cpn.match.chrs.raw.cpn) >= 0.04


          if (all.seg.df[samp.seg.name, paste0("psc_", cpnX, "_high")] <= all.seg.df[samp.seg.name, paste0("ph_int_", cpnX)]) {
            if (all.seg.df[samp.seg.name, paste0("sc_", cpnX, "_high_or_low_dominate")] %in% "HIGH") {
              directionality.vs.oth.cpn.match.chrs.raw.cpn <- mean(ph.cpnX.vec) < mean.oth.cpn.match.chrs.raw.cpn
            } else if (all.seg.df[samp.seg.name, paste0("sc_", cpnX, "_high_or_low_dominate")] %in% "LOW") {
              directionality.vs.oth.cpn.match.chrs.raw.cpn <- mean(ph.cpnX.vec) > mean.oth.cpn.match.chrs.raw.cpn
            }
          } else if (all.seg.df[samp.seg.name, paste0("psc_", cpnX, "_high")] - all.seg.df[samp.seg.name, paste0("ph_int_", cpnX)] == 1) {
            if (all.seg.df[samp.seg.name, paste0("sc_", cpnX, "_high_or_low_dominate")] %in% "HIGH") {
              directionality.vs.oth.cpn.match.chrs.raw.cpn <- mean(ph.cpnX.vec) > mean.oth.cpn.match.chrs.raw.cpn
            } else if (all.seg.df[samp.seg.name, paste0("sc_", cpnX, "_high_or_low_dominate")] %in% "LOW") {
              directionality.vs.oth.cpn.match.chrs.raw.cpn <- mean(ph.cpnX.vec) < mean.oth.cpn.match.chrs.raw.cpn
            }
          }
        }

        all.seg.df[samp.seg.name, paste0("vs_other_chrs_same_", cpnX, "_pval")] <- cpnX.pval.vs.other.chroms

        if (all.seg.df[samp.seg.name, paste0("frac_psc_", cpnX, "_low")] > 1 |
          all.seg.df[samp.seg.name, paste0("frac_psc_", cpnX, "_high")] > 1 |
          all.seg.df[samp.seg.name, paste0("sc_", cpnX, "_vs_min_CCF_pval")] >= clonal.pval |
          all.seg.df[samp.seg.name, paste0("vs_other_chrs_same_", cpnX, "_pval")] >= diff.from.oth.p.val) {
          all.seg.df[samp.seg.name, paste0("within_samp_sc_", cpnX, "_present")] <- FALSE
        } else {
          all.seg.df[samp.seg.name, paste0("within_samp_sc_", cpnX, "_present")] <- TRUE
        }

        if (all.seg.df[samp.seg.name, paste0("frac_psc_", cpnX, "_low")] > 1 |
          all.seg.df[samp.seg.name, paste0("frac_psc_", cpnX, "_high")] > 1 |
          all.seg.df[samp.seg.name, paste0("sc_", cpnX, "_vs_min_cpn_pval")] >= clonal.pval |
          all.seg.df[samp.seg.name, paste0("sc_", cpnX, "_vs_min_CCF_pval")] >= clonal.pval |
          all.seg.df[samp.seg.name, paste0("vs_other_chrs_same_", cpnX, "_pval")] >= diff.from.oth.p.val) {
          all.seg.df[samp.seg.name, paste0("within_samp_sc_", cpnX, "_present_using_min_cpn_diff")] <- FALSE
        } else {
          if (cpnX.mean.and.median.flag == TRUE & abs.diff.from.oth.match.cpn.chrs.flag == TRUE & directionality.vs.oth.cpn.match.chrs.raw.cpn == TRUE) {
            all.seg.df[samp.seg.name, paste0("within_samp_sc_", cpnX, "_present_using_min_cpn_diff")] <- TRUE
          } else {
            all.seg.df[samp.seg.name, paste0("within_samp_sc_", cpnX, "_present_using_min_cpn_diff")] <- FALSE
          }
        }
      }
    }
  }

  all.seg.df
}


addLogrDefinedCPNEventsToOutputSegmentationParallel <- function(sample.vec,
                                                                seg.name.vec,
                                                                all.seg.df,
                                                                ploidy.est.vec,
                                                                purity.est.vec,
                                                                all.samp.logr.df.list,
                                                                all.samp.ph.BAF.df.list,
                                                                gain.threshold,
                                                                loss.threshold,
                                                                amp.threshold,
                                                                deep.loss.threshold,
                                                                tot.cpn.vs.ploidy.pval.thresh) {

  checkSegGainOrLossSignificance <- function(logr.vec, samp.to.check, ploidy.est.vec, purity.est.vec, threshold, cpn.change.type = "gain", min.num.for.logr.test = 3) {
    min.thresh.raw.exp.samp <- ploidy.est.vec[samp.to.check] * (2^(threshold))
    min.thresh.logr.exp.for.samp <- calcExpectedLogR(min.thresh.raw.exp.samp, 0, ploidy.est.vec[samp.to.check], purity.est.vec[samp.to.check])
    median.seg.logr <- median(logr.vec)

    if (length(logr.vec) >= min.num.for.logr.test) {
      gen.t.test.results <- tTestWithErrorHandling(logr.vec, mu = min.thresh.logr.exp.for.samp)

      if (cpn.change.type == "gain") {
        specific.t.test.alt <- "greater"
      } else if (cpn.change.type == "loss") {
        specific.t.test.alt <- "less"
      }

      specific.t.test.result <- tTestWithErrorHandling(logr.vec, mu = min.thresh.logr.exp.for.samp, alternative = specific.t.test.alt)

      low.obs.logr.conf.int <- specific.t.test.result$conf.int[1]
      high.obs.logr.conf.int <- specific.t.test.result$conf.int[2]
      gen.t.test.pval <- gen.t.test.results$p.value
      alt.t.test.pval <- specific.t.test.result$p.value
    }

    # Todo, refactor this magic number...
    if (length(logr.vec) < 3) {
      low.obs.logr.conf.int <- NA
      high.obs.logr.conf.int <- NA
      gen.t.test.pval <- NA
      alt.t.test.pval <- NA
    }

    out <- cbind(min.thresh.raw.exp.samp,
      min.thresh.logr.exp.for.samp,
      median.seg.logr,
      low.obs.logr.conf.int,
      high.obs.logr.conf.int,
      gen.t.test.pval,
      alt.t.test.pval,
      stringsAsFactors = FALSE
    )

    out
  }


  calcExpectedLogR <- function(exp.intA, exp.intB, tumor.ploidy, tumor.purity) {
    exp.logr <- log2((2 * (1 - tumor.purity) + tumor.purity * (exp.intA + exp.intB)) / (2 * (1 - tumor.purity) + (tumor.purity * tumor.ploidy)))

    return(exp.logr)
  }


  tTestWithErrorHandling <- function(x, y = NULL, alternative = c("two.sided", "less", "greater"),
                                     mu = 0, paired = FALSE, var.equal = FALSE, conf.level = 0.95, errors.not.permitted = "") {
    test.results <- tryCatch(
      {
        test.results <- t.test(x = x, y = y, alternative = alternative, mu = mu, paired = paired, var.equal = var.equal, conf.level = conf.level)

        if (is.na(test.results[["p.value"]]) | is.nan(test.results[["p.value"]])) {
          test.results[["stat"]] <- 1
          test.results[["conf.int"]] <- c(NA, NA)
          test.results[["p.value"]] <- 1
          test.results[["p.val"]] <- 1
        }
        return(test.results)
      },
      error = function(err) {
        message(err)
        write(paste(err$message, "\n"), stderr())
        write("\n", stderr())
        if (grepl("NULL", err$message) & "NULL" %in% errors.not.permitted) {
          stop("t.test ERROR with a NULL")
        } else if (grepl("enough", err$message) & "enough" %in% errors.not.permitted) {
          stop("t.test ERROR with not enough")
        } else if (grepl("constant", err$message) & "constant" %in% errors.not.permitted) {
          stop("t.test ERROR with constant data")
        } else if (grepl("enough|constant|NULL", err$message) & "other" %in% errors.not.permitted) {
          stop(err$message)
        }
        test.results <- list()
        test.results[["stat"]] <- 1
        test.results[["conf.int"]] <- c(NA, NA)
        test.results[["p.value"]] <- 1
        test.results[["p.val"]] <- 1
        return(test.results)
      }
    )

    test.results
  }

  # Make vectors for sample vec:
  all.samp.seg.name.vec <- c()
  sample.rep.vec <- c()
  for (samp in sample.vec) {
    all.samp.seg.name.vec <- c(all.samp.seg.name.vec, seg.name.vec)
    sample.rep.vec <- c(sample.rep.vec, rep(samp, (length(seg.name.vec))))
  }
  iterable.samp.seg.df <- cbind(sample = sample.rep.vec, seg_name = all.samp.seg.name.vec)
  itx <- iter(iterable.samp.seg.df, by = "row")
  logdebug("Starting parallel addLogrDefinedCPNEventsToOutputSegmentation:")
  MultiResultClass <- function() {
    me <- list(
      gain_hz_logr_pval = NA,
      amp_hz_logr_pval = NA,
      loss_hz_logr_pval = NA,
      deep_loss_hz_logr_pval = NA,
      gain_non_hz_logr_pval = NA,
      amp_non_hz_logr_pval = NA,
      loss_non_hz_pval = NA,
      deep_loss_non_hz_pval = NA,
      gain_comb_logr_pval = NA,
      amp_comb_logr_pval = NA,
      loss_comb_logr_pval = NA,
      deep_loss_comb_logr_pval = NA,
      cpn_event_vs_ploidy = NA
    )
    class(me) <- append(class(me), "MultiResultClass")
    return(me)
  }

  cpn.events.vs.ploidy.by.samp.seg.list <- foreach(samp.seg.df = itx, .combine = c, .errorhandling = "pass") %dopar% {
    ret.mrc <- MultiResultClass()
    samp <- samp.seg.df[, 1]
    seg.name <- samp.seg.df[, 2]
    chr.start.end.vec <- unlist(strsplit(seg.name, split = ":"))


    datapoints.w.BAF.df <- all.samp.ph.BAF.df.list[[samp]]
    ss.baf.lr.baf.df <- datapoints.w.BAF.df[datapoints.w.BAF.df$chrom == as.numeric(chr.start.end.vec[1]) &
      datapoints.w.BAF.df$pos >= as.numeric(chr.start.end.vec[2]) &
      datapoints.w.BAF.df$pos <= as.numeric(chr.start.end.vec[3]), , drop = FALSE]

    gain.hz.only.samp.df <- checkSegGainOrLossSignificance(ss.baf.lr.baf.df$logR, samp, ploidy.est.vec, purity.est.vec, gain.threshold, "gain", min.num.for.logr.test = 3)
    amp.hz.only.samp.df <- checkSegGainOrLossSignificance(ss.baf.lr.baf.df$logR, samp, ploidy.est.vec, purity.est.vec, amp.threshold, "gain", min.num.for.logr.test = 3)
    loss.hz.only.samp.df <- checkSegGainOrLossSignificance(ss.baf.lr.baf.df$logR, samp, ploidy.est.vec, purity.est.vec, loss.threshold, "loss", min.num.for.logr.test = 3)
    deep.loss.hz.only.samp.df <- checkSegGainOrLossSignificance(ss.baf.lr.baf.df$logR, samp, ploidy.est.vec, purity.est.vec, deep.loss.threshold, "loss", min.num.for.logr.test = 3)
    ret.mrc$gain_hz_logr_pval <- amp.hz.only.samp.df[, 7]
    ret.mrc$amp_hz_logr_pval <- amp.hz.only.samp.df[, 7]
    ret.mrc$loss_hz_logr_pval <- loss.hz.only.samp.df[, 7]
    ret.mrc$deep_loss_hz_logr_pval <- deep.loss.hz.only.samp.df[, 7]

    datapoints.w.logr.only.df <- all.samp.logr.df.list[[samp]]

    if (nrow(datapoints.w.logr.only.df) > 0) {
      ss.dp.w.logr.only.df <- datapoints.w.logr.only.df[datapoints.w.logr.only.df$chrom == as.numeric(chr.start.end.vec[1]) &
        datapoints.w.logr.only.df$pos >= as.numeric(chr.start.end.vec[2]) &
        datapoints.w.logr.only.df$pos <= as.numeric(chr.start.end.vec[3]), , drop = FALSE]

      gain.non.hz.samp.df <- checkSegGainOrLossSignificance(ss.dp.w.logr.only.df$logR, samp, ploidy.est.vec, purity.est.vec, gain.threshold, "gain", min.num.for.logr.test = 3)
      amp.non.hz.samp.df <- checkSegGainOrLossSignificance(ss.dp.w.logr.only.df$logR, samp, ploidy.est.vec, purity.est.vec, amp.threshold, "gain", min.num.for.logr.test = 3)
      loss.non.hz.samp.df <- checkSegGainOrLossSignificance(ss.dp.w.logr.only.df$logR, samp, ploidy.est.vec, purity.est.vec, loss.threshold, "loss", min.num.for.logr.test = 3)
      deep.loss.non.hz.samp.df <- checkSegGainOrLossSignificance(ss.dp.w.logr.only.df$logR, samp, ploidy.est.vec, purity.est.vec, deep.loss.threshold, "loss", min.num.for.logr.test = 3)
      ret.mrc$gain_non_hz_logr_pval <- gain.non.hz.samp.df[, 7]
      ret.mrc$amp_non_hz_logr_pval <- amp.non.hz.samp.df[, 7]
      ret.mrc$loss_non_hz_pval <- loss.non.hz.samp.df[, 7]
      ret.mrc$deep_loss_non_hz_pval <- deep.loss.non.hz.samp.df[, 7]

      # now we can simply add these data points to the dataframe that contained the datapoints that had both logr and BAF, do we need to re-sort? leave for now.
      if (nrow(ss.dp.w.logr.only.df) > 0) {
        ss.comb.baf.lr.baf.df <- rbind(ss.baf.lr.baf.df[, c("chrom", "pos", "BAF", "BAFsegmented", "logR", "logRsegmented")],
          ss.dp.w.logr.only.df[, c("chrom", "pos", "BAF", "BAFsegmented", "logR", "logRsegmented")],
          stringsAsFactors = FALSE
        )
      } else {
        ss.comb.baf.lr.baf.df <- ss.baf.lr.baf.df
      }
    } else {
      ss.comb.baf.lr.baf.df <- ss.baf.lr.baf.df
    }

    gain.comb.logr.samp.df <- checkSegGainOrLossSignificance(ss.comb.baf.lr.baf.df$logR, samp, ploidy.est.vec, purity.est.vec, gain.threshold, "gain", min.num.for.logr.test = 3)
    amp.comb.logr.samp.df <- checkSegGainOrLossSignificance(ss.comb.baf.lr.baf.df$logR, samp, ploidy.est.vec, purity.est.vec, amp.threshold, "gain", min.num.for.logr.test = 3)
    loss.comb.logr.samp.df <- checkSegGainOrLossSignificance(ss.comb.baf.lr.baf.df$logR, samp, ploidy.est.vec, purity.est.vec, loss.threshold, "loss", min.num.for.logr.test = 3)
    deep.loss.comb.logr.samp.df <- checkSegGainOrLossSignificance(ss.comb.baf.lr.baf.df$logR, samp, ploidy.est.vec, purity.est.vec, deep.loss.threshold, "loss", min.num.for.logr.test = 3)
    ret.mrc$gain_comb_logr_pval <- gain.comb.logr.samp.df[, 7]
    ret.mrc$amp_comb_logr_pval <- amp.comb.logr.samp.df[, 7]
    ret.mrc$loss_comb_logr_pval <- loss.comb.logr.samp.df[, 7]
    ret.mrc$deep_loss_comb_logr_pval <- deep.loss.comb.logr.samp.df[, 7]

    cpn.event.vs.ploidy.vec <- as.numeric(c(gain.comb.logr.samp.df[, 7], amp.comb.logr.samp.df[, 7], loss.comb.logr.samp.df[, 7], deep.loss.comb.logr.samp.df[, 7]))

    names(cpn.event.vs.ploidy.vec) <- c("gain", "amp", "loss", "deep_loss")
    non.na.seg.cpn.event <- cpn.event.vs.ploidy.vec[!is.na(cpn.event.vs.ploidy.vec)]
    if (length(non.na.seg.cpn.event) == 0) {
      ret.mrc$cpn_event_vs_ploidy <- NA # had less than min.num of datapoints to examine for a logr event.
    } else {
      significant.seg.cpn.event <- non.na.seg.cpn.event[non.na.seg.cpn.event <= tot.cpn.vs.ploidy.pval.thresh]
      # Error check here:
      if (length(significant.seg.cpn.event) > 1) {
        if ("amp" %in% names(significant.seg.cpn.event)) {
          ret.mrc$cpn_event_vs_ploidy <- "amp"
        } else if ("deep_loss" %in% names(significant.seg.cpn.event)) {
          ret.mrc$cpn_event_vs_ploidy <- "deep_loss"
        } else {
          logwarn(significant.seg.cpn.event)
          stop("Error, more than one cpn event vs ploidy is signficant and neither are amp or deep_loss, likely bug present")
        }
      } else if (length(significant.seg.cpn.event) == 1) {
        ret.mrc$cpn_event_vs_ploidy <- names(significant.seg.cpn.event)
      } else {
        ret.mrc$cpn_event_vs_ploidy <- "neutral"
      }
    }

    returnable.list <- list(ret.mrc)
    names(returnable.list) <- paste0(samp.seg.df[, 1], ":", samp.seg.df[, 2])
    returnable.list
  }

  for (samp in sample.vec) {
    for (seg.name in seg.name.vec) {
      if (length(sample.vec) > 1) {
        samp.seg.name <- paste0(samp, ":", seg.name)
      } else {
        samp.seg.name <- seg.name
      }


      ret.mrc <- cpn.events.vs.ploidy.by.samp.seg.list[[samp.seg.name]]

      for (j in 1:length(ret.mrc)) {
        all.seg.df[samp.seg.name, names(ret.mrc)[j]] <- ret.mrc[names(ret.mrc)[j]]
      }
    }
  }

  all.seg.df
}


addLogrDefinedCPNEventsToOutputSegmentation <- function(sample.vec,
                                                        seg.name.vec,
                                                        all.seg.df,
                                                        ploidy.est.vec,
                                                        purity.est.vec,
                                                        all.samp.logr.df.list,
                                                        all.samp.ph.BAF.df.list,
                                                        gain.threshold,
                                                        loss.threshold,
                                                        amp.threshold,
                                                        deep.loss.threshold,
                                                        tot.cpn.vs.ploidy.pval.thresh) {
  for (samp in sample.vec) {
    for (seg.name in seg.name.vec) {
      if (length(sample.vec) > 1) {
        samp.seg.name <- paste0(samp, ":", seg.name)
      } else {
        samp.seg.name <- seg.name
      }

      chr.start.end.vec <- unlist(strsplit(seg.name, split = ":"))
      datapoints.w.BAF.df <- all.samp.ph.BAF.df.list[[samp]]
      ss.baf.lr.baf.df <- datapoints.w.BAF.df[datapoints.w.BAF.df$chrom == as.numeric(chr.start.end.vec[1]) &
        datapoints.w.BAF.df$pos >= as.numeric(chr.start.end.vec[2]) &
        datapoints.w.BAF.df$pos <= as.numeric(chr.start.end.vec[3]), , drop = FALSE]

      gain.hz.only.samp.df <- checkSegGainOrLossSignificance(ss.baf.lr.baf.df$logR, samp, ploidy.est.vec, purity.est.vec, gain.threshold, "gain", min.num.for.logr.test = 3)
      amp.hz.only.samp.df <- checkSegGainOrLossSignificance(ss.baf.lr.baf.df$logR, samp, ploidy.est.vec, purity.est.vec, amp.threshold, "gain", min.num.for.logr.test = 3)
      loss.hz.only.samp.df <- checkSegGainOrLossSignificance(ss.baf.lr.baf.df$logR, samp, ploidy.est.vec, purity.est.vec, loss.threshold, "loss", min.num.for.logr.test = 3)
      deep.loss.hz.only.samp.df <- checkSegGainOrLossSignificance(ss.baf.lr.baf.df$logR, samp, ploidy.est.vec, purity.est.vec, deep.loss.threshold, "loss", min.num.for.logr.test = 3)
      all.seg.df[samp.seg.name, ]$gain_hz_logr_pval <- amp.hz.only.samp.df[, 7]
      all.seg.df[samp.seg.name, ]$amp_hz_logr_pval <- amp.hz.only.samp.df[, 7]
      all.seg.df[samp.seg.name, ]$loss_hz_logr_pval <- loss.hz.only.samp.df[, 7]
      all.seg.df[samp.seg.name, ]$deep_loss_hz_logr_pval <- deep.loss.hz.only.samp.df[, 7]


      datapoints.w.logr.only.df <- all.samp.logr.df.list[[samp]]
      if (nrow(datapoints.w.logr.only.df) > 0) {
        ss.dp.w.logr.only.df <- datapoints.w.logr.only.df[datapoints.w.logr.only.df$chrom == as.numeric(chr.start.end.vec[1]) &
          datapoints.w.logr.only.df$pos >= as.numeric(chr.start.end.vec[2]) &
          datapoints.w.logr.only.df$pos <= as.numeric(chr.start.end.vec[3]), , drop = FALSE]

        gain.non.hz.samp.df <- checkSegGainOrLossSignificance(ss.dp.w.logr.only.df$logR, samp, ploidy.est.vec, purity.est.vec, gain.threshold, "gain", min.num.for.logr.test = 3)
        amp.non.hz.samp.df <- checkSegGainOrLossSignificance(ss.dp.w.logr.only.df$logR, samp, ploidy.est.vec, purity.est.vec, amp.threshold, "gain", min.num.for.logr.test = 3)
        loss.non.hz.samp.df <- checkSegGainOrLossSignificance(ss.dp.w.logr.only.df$logR, samp, ploidy.est.vec, purity.est.vec, loss.threshold, "loss", min.num.for.logr.test = 3)
        deep.loss.non.hz.samp.df <- checkSegGainOrLossSignificance(ss.dp.w.logr.only.df$logR, samp, ploidy.est.vec, purity.est.vec, deep.loss.threshold, "loss", min.num.for.logr.test = 3)
        all.seg.df[samp.seg.name, ]$gain_non_hz_logr_pval <- gain.non.hz.samp.df[, 7]
        all.seg.df[samp.seg.name, ]$amp_non_hz_logr_pval <- amp.non.hz.samp.df[, 7]
        all.seg.df[samp.seg.name, ]$loss_non_hz_pval <- loss.non.hz.samp.df[, 7]
        all.seg.df[samp.seg.name, ]$deep_loss_non_hz_pval <- deep.loss.non.hz.samp.df[, 7]

        # now we can simply add these data points to the dataframe that contained the datapoints that had both logr and BAF, do we need to re-sort? leave for now.
        if (nrow(ss.dp.w.logr.only.df) > 0) {
          ss.comb.baf.lr.baf.df <- rbind(ss.baf.lr.baf.df[, c("chrom", "pos", "BAF", "BAFsegmented", "logR", "logRsegmented")],
            ss.dp.w.logr.only.df[, c("chrom", "pos", "BAF", "BAFsegmented", "logR", "logRsegmented")],
            stringsAsFactors = FALSE
          )
        } else {
          ss.comb.baf.lr.baf.df <- ss.baf.lr.baf.df
        }
      } else {
        ss.comb.baf.lr.baf.df <- ss.baf.lr.baf.df
      }

      gain.comb.logr.samp.df <- checkSegGainOrLossSignificance(ss.comb.baf.lr.baf.df$logR, samp, ploidy.est.vec, purity.est.vec, gain.threshold, "gain", min.num.for.logr.test = 3)
      amp.comb.logr.samp.df <- checkSegGainOrLossSignificance(ss.comb.baf.lr.baf.df$logR, samp, ploidy.est.vec, purity.est.vec, amp.threshold, "gain", min.num.for.logr.test = 3)
      loss.comb.logr.samp.df <- checkSegGainOrLossSignificance(ss.comb.baf.lr.baf.df$logR, samp, ploidy.est.vec, purity.est.vec, loss.threshold, "loss", min.num.for.logr.test = 3)
      deep.loss.comb.logr.samp.df <- checkSegGainOrLossSignificance(ss.comb.baf.lr.baf.df$logR, samp, ploidy.est.vec, purity.est.vec, deep.loss.threshold, "loss", min.num.for.logr.test = 3)
      all.seg.df[samp.seg.name, ]$gain_comb_logr_pval <- gain.comb.logr.samp.df[, 7]
      all.seg.df[samp.seg.name, ]$amp_comb_logr_pval <- amp.comb.logr.samp.df[, 7]
      all.seg.df[samp.seg.name, ]$loss_comb_logr_pval <- loss.comb.logr.samp.df[, 7]
      all.seg.df[samp.seg.name, ]$deep_loss_comb_logr_pval <- deep.loss.comb.logr.samp.df[, 7]

      cpn.event.vs.ploidy.vec <- as.numeric(c(gain.comb.logr.samp.df[, 7], amp.comb.logr.samp.df[, 7], loss.comb.logr.samp.df[, 7], deep.loss.comb.logr.samp.df[, 7]))

      names(cpn.event.vs.ploidy.vec) <- c("gain", "amp", "loss", "deep_loss")
      non.na.seg.cpn.event <- cpn.event.vs.ploidy.vec[!is.na(cpn.event.vs.ploidy.vec)]
      if (length(non.na.seg.cpn.event) == 0) {
        all.seg.df[samp.seg.name, ]$cpn_event_vs_ploidy <- NA # had less than min.num of datapoints to examine for a logr event.
      } else {
        significant.seg.cpn.event <- non.na.seg.cpn.event[non.na.seg.cpn.event <= tot.cpn.vs.ploidy.pval.thresh]
        if (length(significant.seg.cpn.event) > 1) {
          if ("amp" %in% names(significant.seg.cpn.event)) {
            all.seg.df[samp.seg.name, ]$cpn_event_vs_ploidy <- "amp"
          } else if ("deep_loss" %in% names(significant.seg.cpn.event)) {
            all.seg.df[samp.seg.name, ]$cpn_event_vs_ploidy <- "deep_loss"
          } else {
            logwarn(significant.seg.cpn.event)
            stop("Error, more than one cpn event vs ploidy is signficant and neither are amp or deep_loss, likely bug present")
          }
        } else if (length(significant.seg.cpn.event) == 1) {
          all.seg.df[samp.seg.name, ]$cpn_event_vs_ploidy <- names(significant.seg.cpn.event)
        } else {
          all.seg.df[samp.seg.name, ]$cpn_event_vs_ploidy <- "neutral"
        }
      }
    }
  }

  all.seg.df
}


addMSAIDetectionByLTNPDToOutputSegmentation <- function(sample.vec,
                                                        seg.name.vec,
                                                        all.seg.MSAI.by.BAF.ltnpd.present.list,
                                                        all.segs.df) {
  for (seg.name in seg.name.vec) {
    seg.MSAI.info.list <- all.seg.MSAI.by.BAF.ltnpd.present.list[[seg.name]]
    samp.seg.name.vec <- paste0(sample.vec, ":", seg.name)

    if (seg.MSAI.info.list$MSAI_detection_feasible == FALSE) {

      # so now get all the relevant segments and mark them as not being possible to mark
      all.segs.df[samp.seg.name.vec, ]$MSAI <- NA
      all.segs.df[samp.seg.name.vec, ]$any_across_samp_MSAI <- NA
    } else if (seg.MSAI.info.list$MSAI_detection_feasible == TRUE) {
      if (seg.MSAI.info.list$MSAI_present_any_sample_by_thresh == FALSE) {
        all.segs.df[samp.seg.name.vec, ]$any_across_samp_MSAI <- FALSE
        all.segs.df[samp.seg.name.vec, ]$MSAI <- FALSE
        all.segs.df[samp.seg.name.vec, ]$dominant_allele <- "cpnA"
      } else if (seg.MSAI.info.list$MSAI_present_any_sample_by_thresh == TRUE) {
        all.segs.df[samp.seg.name.vec, ]$any_across_samp_MSAI <- TRUE
        # So now regernate the samp.seg.name.vec to be sure it's in the same order as the MSAI df
        seg.MSAI.df <- seg.MSAI.info.list$MSAI_by_sample_df
        samp.seg.name.vec <- paste0(seg.MSAI.df$sample_name, ":", seg.name)
        all.segs.df[samp.seg.name.vec, ]$MSAI <- seg.MSAI.df$MSAI_pval_passed_thresh
        all.segs.df[samp.seg.name.vec, ]$dominant_allele <- seg.MSAI.df$poss_dominant_allele
      }
    }
  }

  all.segs.df
}


addParallelEventsByMSAIToOutputSegmentation <- function(all.segs.df) {

  # So first just subset by anything with MSAI calls:
  matched.MSAI.loss.segs.df <- all.segs.df[(!is.na(all.segs.df$any_across_samp_MSAI))
  & all.segs.df$any_across_samp_MSAI == TRUE
  & all.segs.df$cpn_event_vs_ploidy %in% c("loss", "deep_loss"), ]

  matched.MSAI.gain.segs.df <- all.segs.df[(!is.na(all.segs.df$any_across_samp_MSAI))
  & all.segs.df$any_across_samp_MSAI == TRUE
  & all.segs.df$cpn_event_vs_ploidy %in% c("gain", "amp"), ]

  matched.MSAI.LOH.segs.df <- all.segs.df[(!is.na(all.segs.df$any_across_samp_MSAI))
  & all.segs.df$any_across_samp_MSAI == TRUE
  & all.segs.df$min_int_cpnX == 0, ]

  matched.MSAI.LOH.or.loss.segs.df <- all.segs.df[(!is.na(all.segs.df$any_across_samp_MSAI))
  & all.segs.df$any_across_samp_MSAI == TRUE
  & (all.segs.df$min_int_cpnX == 0 | all.segs.df$cpn_event_vs_ploidy %in% c("loss", "deep_loss")), ]

  MSAI.gain.seg.names.vec <- unique(matched.MSAI.gain.segs.df$seg_name)
  MSAI.loss.seg.names.vec <- unique(matched.MSAI.loss.segs.df$seg_name)
  MSAI.LOH.seg.names.vec <- unique(matched.MSAI.LOH.segs.df$seg_name)
  MSAI.LOH.or.loss.seg.names.vec <- unique(matched.MSAI.LOH.or.loss.segs.df$seg_name)

  # Todo add strict AMP definition here....
  for (seg.name in MSAI.gain.seg.names.vec) {
    seg.MSAI.gain.df <- matched.MSAI.gain.segs.df[matched.MSAI.gain.segs.df$seg_name == seg.name, ]

    if ("cpnA" %in% seg.MSAI.gain.df$dominant_allele & "cpnB" %in% seg.MSAI.gain.df$dominant_allele) {
      all.segs.df[rownames(seg.MSAI.gain.df), ]$is_parallel_gain_and_or_amp <- TRUE
    }
  }

  for (seg.name in MSAI.loss.seg.names.vec) {
    seg.MSAI.loss.df <- matched.MSAI.loss.segs.df[matched.MSAI.loss.segs.df$seg_name == seg.name, ]

    if ("cpnA" %in% seg.MSAI.loss.df$dominant_allele & "cpnB" %in% seg.MSAI.loss.df$dominant_allele) {
      all.segs.df[rownames(seg.MSAI.loss.df), ]$is_parallel_loss <- TRUE
    }
  }

  for (seg.name in MSAI.LOH.seg.names.vec) {
    seg.MSAI.LOH.df <- matched.MSAI.LOH.segs.df[matched.MSAI.LOH.segs.df$seg_name == seg.name, ]

    if ("cpnA" %in% seg.MSAI.LOH.df$dominant_allele & "cpnB" %in% seg.MSAI.LOH.df$dominant_allele) {
      all.segs.df[rownames(seg.MSAI.LOH.df), ]$is_parallel_LOH <- TRUE
    }
  }

  for (seg.name in MSAI.LOH.or.loss.seg.names.vec) {
    seg.MSAI.LOH.or.loss.df <- matched.MSAI.LOH.or.loss.segs.df[matched.MSAI.LOH.or.loss.segs.df$seg_name == seg.name, ]

    if ("cpnA" %in% seg.MSAI.LOH.or.loss.df$dominant_allele & "cpnB" %in% seg.MSAI.LOH.or.loss.df$dominant_allele) {
      all.segs.df[rownames(seg.MSAI.LOH.or.loss.df), ]$is_parallel_loss_and_or_LOH <- TRUE
    }
  }

  # Todo, Subclonal parallel loss/gains..
  # Todo, parallel amps /deep losses?
  # Todo, use the conversion to parallel amps?
  all.segs.df
}


addAcrossSampleAlleleSpecificClonalityCalls <- function(all.seg.df, sample.vec, seg.name.vec) {

  for (seg.name in seg.name.vec) {
    samp.seg.name.vec <- paste0(sample.vec, ":", seg.name)
    # all.samps.seg.df <- all.seg.df[samp.seg.name.vec,]
    if (is.na(all.seg.df[samp.seg.name.vec, ]$within_samp_sc_cpnA_present[1]) | is.na(all.seg.df[samp.seg.name.vec, ]$within_samp_sc_cpnB_present[1])) {
      next
    } else {
      all.seg.df[samp.seg.name.vec, ]$any_within_samp_sc_cpnA <- any(all.seg.df[samp.seg.name.vec, ]$within_samp_sc_cpnA_present)
      all.seg.df[samp.seg.name.vec, ]$any_within_samp_sc_cpnB <- any(all.seg.df[samp.seg.name.vec, ]$within_samp_sc_cpnB_present)

      # Now do we want to bring in MSAI here? It's not really clonal if there is any instance of MSAI affecting any of the segments... Todo, figure this out.
      all.seg.df[samp.seg.name.vec, ]$across_samp_cl_AI <- all(all.seg.df[samp.seg.name.vec, ]$ph_int_cpnA != all.seg.df[samp.seg.name.vec, ]$ph_int_cpnB)
      # across_samp_sc_AI = rep(NA, nr),
      all.seg.df[samp.seg.name.vec, ]$across_samp_sc_AI <- any(all.seg.df[samp.seg.name.vec, ]$ph_int_cpnA != all.seg.df[samp.seg.name.vec, ]$ph_int_cpnB) &
        (!all.seg.df[samp.seg.name.vec, ]$across_samp_cl_AI)
      # across_samp_cl_LOH = rep(NA, nr),
      all.seg.df[samp.seg.name.vec, ]$across_samp_cl_LOH <- all(all.seg.df[samp.seg.name.vec, ]$min_int_cpnX == 0)
      # across_samp_sc_LOH = rep(NA, nr),
      all.seg.df[samp.seg.name.vec, ]$across_samp_sc_LOH <- any(all.seg.df[samp.seg.name.vec, ]$min_int_cpnX == 0) &
        (!all.seg.df[samp.seg.name.vec, ]$across_samp_cl_LOH)
    }
  }

  all.seg.df
}


addPanCanColumnDefinitionToSegDF <- function(all.seg.df) {
  for (seg.name in seg.name.vec) {
    samp.seg.name.vec <- paste0(sample.vec, ":", seg.name)

    ss.df <- all.seg.df[samp.seg.name.vec, ]

    # TODO FIX THIS, currently reduces the amount of MSAI that's not correct. For now use other columns...
    if (is.na(all.seg.df[samp.seg.name.vec, ]$within_samp_sc_cpnA_present[1]) | is.na(all.seg.df[samp.seg.name.vec, ]$within_samp_sc_cpnB_present[1])) {
      next
    } else {

      # Todo refactor the strict definitions so they work correctly...
      all.seg.df[samp.seg.name.vec, ]$pt_seg_parallel_MSAI_strict_gain <- any(ss.df$is_parallel_gain_and_or_amp, na.rm = TRUE)
      all.seg.df[samp.seg.name.vec, ]$pt_seg_parallel_MSAI_gain_or_amp <- any(ss.df$is_parallel_gain_and_or_amp, na.rm = TRUE)
      all.seg.df[samp.seg.name.vec, ]$pt_seg_parallel_MSAI_strict_amp <- any(ss.df$is_parallel_gain_and_or_amp, na.rm = TRUE)

      all.seg.df[samp.seg.name.vec, ]$pt_seg_parallel_MSAI_loss <- any(ss.df$is_parallel_loss, na.rm = TRUE)
      all.seg.df[samp.seg.name.vec, ]$pt_seg_parallel_MSAI_LOH <- any(ss.df$is_parallel_LOH, na.rm = TRUE)
      all.seg.df[samp.seg.name.vec, ]$pt_seg_parallel_MSAI_loss_or_LOH <- any(ss.df$is_parallel_loss_and_or_LOH, na.rm = TRUE)
      # TODO revisit this! and change definition so that it's EITHER clonal LOH OR clonal loss vs ploidy not some combination...
      all.seg.df[samp.seg.name.vec, ]$pt_seg_cl_loss_LOH_no_MSAI <- (!any(ss.df$any_across_samp_MSAI, na.rm = TRUE)) & all(ss.df$cpn_event_vs_ploidy %in% c("loss", "deep_loss") | ss.df$min_int_cpnX == 0)
      all.seg.df[samp.seg.name.vec, ]$pt_seg_cl_gain_amp_no_MSAI <- (!any(ss.df$any_across_samp_MSAI, na.rm = TRUE)) & all(ss.df$cpn_event_vs_ploidy %in% c("gain", "amp"))
      all.seg.df[samp.seg.name.vec, ]$pt_seg_sc_loss_LOH_incl_MSAI <- (any(ss.df$any_across_samp_MSAI, na.rm = TRUE) & all(ss.df$cpn_event_vs_ploidy %in% c("loss", "amp"))) | (!(all(ss.df$cpn_event_vs_ploidy %in% c("loss", "deep_loss") | ss.df$min_int_cpnX == 0))) & (any(ss.df$cpn_event_vs_ploidy %in% c("loss", "deep_loss") | ss.df$min_int_cpnX == 0))
      all.seg.df[samp.seg.name.vec, ]$pt_seg_sc_gain_amp_incl_MSAI <- (any(ss.df$any_across_samp_MSAI, na.rm = TRUE) & all(ss.df$cpn_event_vs_ploidy %in% c("gain", "amp"))) | ((!all(ss.df$cpn_event_vs_ploidy %in% c("gain", "amp"))) & any(ss.df$cpn_event_vs_ploidy %in% c("gain", "amp")))
    }
  }

  all.seg.df
}


testFinalOutputDataframe <- function(final.joint.segs.df) {

  # Tests that should hold for the final output:
  potential.msai.vec <- !is.na(final.joint.segs.df$nAraw) & (final.joint.segs.df$nAraw < final.joint.segs.df$nBraw)
  if (!all(as.logical(final.joint.segs.df[potential.msai.vec == TRUE, ]$MSAI) == TRUE)) {
    stop("MSAI column not flled in correctly")
  }

  # Test original adding
  if (!all(final.joint.segs.df$orig_nMajor >= final.joint.segs.df$orig_nMinor)) {
    stop("Original nMajor not always greater than orginal nMinor")
  }
}


initialiseFullOutputDataframe <- function(seg.name.vec, ploidy.est.vec, purity.est.vec, make.row.names.seg.name.alone = FALSE) {
  sample.vec <- names(purity.est.vec)
  seg.chr.vec <- sapply(seg.name.vec, function(x) as.numeric(unlist(strsplit(x, split = ":"))[1]))
  seg.startpos.vec <- sapply(seg.name.vec, function(x) as.numeric(unlist(strsplit(x, split = ":"))[2]))
  seg.endpos.vec <- sapply(seg.name.vec, function(x) as.numeric(unlist(strsplit(x, split = ":"))[3]))
  nr <- length(seg.name.vec) * length(sample.vec)
  sample.rows.vec <- c()
  sample.ploidy.vec <- c()
  sample.purity.vec <- c()
  for (samp in sample.vec) {
    sample.rows.vec <- c(sample.rows.vec, rep(samp, length(seg.name.vec)))
    sample.ploidy.vec <- c(sample.ploidy.vec, rep(ploidy.est.vec[samp], length(seg.name.vec)))
    sample.purity.vec <- c(sample.purity.vec, rep(purity.est.vec[samp], length(seg.name.vec)))
  }
  all.seg.df <- data.frame(
    # New processed calls we might want
    # Strict ASCAT naming convention:
    sample = sample.rows.vec,
    chr = rep(seg.chr.vec, length(sample.vec)),
    startpos = rep(seg.startpos.vec, length(sample.vec)),
    endpos = rep(seg.endpos.vec, length(sample.vec)),
    n.het = rep(0, nr),
    cnTotal = rep(NA, nr),
    nMinor = rep(NA, nr),
    nMajor = rep(NA, nr),
    Ploidy = sample.ploidy.vec,
    ACF = sample.purity.vec,
    nAraw = rep(NA, nr),
    nBraw = rep(NA, nr),
    # MPhase phased copy number calls:
    patient_id = rep(NA, nr),
    data_type = rep(NA, nr),
    analysis_group = rep(NA, nr),
    sample_name = sample.rows.vec,
    seg_name = rep(seg.name.vec, length(sample.vec)),
    chrom = rep(seg.chr.vec, length(sample.vec)),
    start_pos = rep(seg.startpos.vec, length(sample.vec)),
    end_pos = rep(seg.endpos.vec, length(sample.vec)),
    num_ph_hz_sites = rep(0, nr),
    num_un_ph_hz_sites = rep(0, nr),
    num_non_hz_sites = rep(0, nr),
    total_int_cpn = rep(NA, nr),
    ph_int_cpnA = rep(NA, nr),
    ph_int_cpnB = rep(NA, nr),
    ph_hz_total_raw_cpn = rep(NA, nr),
    non_hz_total_raw_cpn = rep(NA, nr),
    total_raw_cpn = rep(NA, nr),
    ph_raw_cpnA = rep(NA, nr),
    ph_raw_cpnB = rep(NA, nr),
    min_int_cpnX = rep(NA, nr),
    max_int_cpnX = rep(NA, nr),
    sample_ploidy = sample.ploidy.vec,
    sample_purity = sample.purity.vec,
    # Subclonal copy number calls
    within_samp_sc_cpnA_present = rep(NA, nr),
    within_samp_sc_cpnA_present_using_min_cpn_diff = rep(NA, nr),
    sc_cpnA_high_or_low_dominate = rep(NA, nr),
    psc_cpnA_low = rep(NA, nr),
    psc_cpnA_high = rep(NA, nr),
    within_samp_sc_cpnB_present = rep(NA, nr),
    within_samp_sc_cpnB_present_using_min_cpn_diff = rep(NA, nr),
    sc_cpnB_high_or_low_dominate = rep(NA, nr),
    psc_cpnB_low = rep(NA, nr),
    psc_cpnB_high = rep(NA, nr),
    # Copy number event calling:
    ph_hz_total_cpn = rep(NA, nr),
    non_hz_total_cpn = rep(NA, nr),
    # Chosen copynumber event vs ploidy based on combined non_hz and hz logr pval
    cpn_event_vs_ploidy = rep(NA, nr),
    gain_conv_to_amp = rep(NA, nr),
    # BAF info
    median_BAF = rep(NA, nr),
    mean_BAF = rep(NA, nr),
    stdev_BAF = rep(NA, nr),
    # hz logr values
    median_hz_logr = rep(NA, nr),
    mean_hz_logr = rep(NA, nr),
    stdev_hz_logr = rep(NA, nr),
    # non hz values
    median_non_hz_logr = rep(NA, nr),
    mean_non_hz_logr = rep(NA, nr),
    stdev_non_hz_logr = rep(NA, nr),
    # combined hz values
    median_comb_logr = rep(NA, nr),
    mean_comb_logr = rep(NA, nr),
    stdev_comb_logr = rep(NA, nr),
    # pvalues for tests hz logr
    gain_hz_logr_pval = rep(NA, nr),
    amp_hz_logr_pval = rep(NA, nr),
    loss_hz_logr_pval = rep(NA, nr),
    deep_loss_hz_logr_pval = rep(NA, nr),
    # pvales for tests with non-hz logr
    gain_non_hz_logr_pval = rep(NA, nr),
    amp_non_hz_logr_pval = rep(NA, nr),
    loss_non_hz_pval = rep(NA, nr),
    deep_loss_non_hz_pval = rep(NA, nr),
    # pvales for tests with non-hz logr
    gain_comb_logr_pval = rep(NA, nr),
    amp_comb_logr_pval = rep(NA, nr),
    loss_comb_logr_pval = rep(NA, nr),
    deep_loss_comb_logr_pval = rep(NA, nr),
    # information for subclonal calling:
    ancestral_cpn_state = rep(NA, nr), # if lost or gained ccf different? what i it relative to, ancestral stated medicc
    loss_vs_gain_ccf = rep(NA, nr),
    sc_cpnA_vs_min_cpn_pval = rep(NA, nr),
    sc_cpnA_vs_min_CCF_pval = rep(NA, nr),
    frac_psc_cpnA_low = rep(NA, nr),
    frac_psc_cpnA_high = rep(NA, nr),
    vs_other_chrs_same_cpnA_pval = rep(NA, nr),
    num_other_chrs_same_cpnA_cpn = rep(NA, nr),
    sc_cpnB_vs_min_cpn_pval = rep(NA, nr),
    sc_cpnB_vs_min_CCF_pval = rep(NA, nr),
    frac_psc_cpnB_low = rep(NA, nr),
    frac_psc_cpnB_high = rep(NA, nr),
    vs_other_chrs_same_cpnB_pval = rep(NA, nr),
    num_other_chrs_same_cpnB_cpn = rep(NA, nr),
    # Link back to original segmentation
    orig_seg_name = rep(NA, nr),
    orig_cpn_call_vs_ploidy = rep(NA, nr),
    orig_nMinor = rep(NA, nr),
    orig_nMajor = rep(NA, nr),
    orig_cpnA_raw = rep(NA, nr),
    orig_cpnB_raw = rep(NA, nr),
    orig_cnTotal = rep(NA, nr),
    nMinor_change = rep(NA, nr),
    nMajor_change = rep(NA, nr),
    cnTotal_change = rep(NA, nr),
    # joined
    # joined_arm_or_focal = rep(NA, nr),
    # Allelic imabalance information
    within_samp_cl_AI_present = rep(NA, nr),
    within_samp_sc_AI_present = rep(NA, nr),
    within_samp_cl_LOH_present = rep(NA, nr),
    within_samp_sc_LOH_present = rep(NA, nr),
    new_within_samp_cl_AI_detected = rep(NA, nr),
    new_within_samp_cl_LOH_detected = rep(NA, nr),
    orig_part_of_arm_or_focal = rep(NA, nr),
    samp_used_to_phase = rep(NA, nr),
    MSAI = rep(NA, nr), # Means tha
    raw_cpn_MSAI = rep(NA, nr), # Todo, decide whether this is used or diff from ref samp
    dominant_allele = rep(NA, nr),
    is_parallel_gain_and_or_amp = rep(NA, nr),
    is_parallel_loss = rep(NA, nr),
    is_parallel_LOH = rep(NA, nr),
    is_parallel_loss_and_or_LOH = rep(NA, nr),
    # disease level stats
    any_across_samp_MSAI = rep(NA, nr),
    any_within_samp_sc_cpnA = rep(NA, nr),
    any_within_samp_sc_cpnB = rep(NA, nr),
    # across_samp_cl_cpnA = rep(NA, nr),#
    # across_samp_sc_cpnA = rep(NA, nr),#
    # across_samp_cl_cpnB = rep(NA, nr),#
    # across_samp_sc_cpnB = rep(NA, nr),#
    across_samp_cl_AI = rep(NA, nr),
    across_samp_sc_AI = rep(NA, nr),
    across_samp_cl_LOH = rep(NA, nr),
    across_samp_sc_LOH = rep(NA, nr),
    across_samp_cl_cpn_vs_ploidy_present = rep(NA, nr),
    across_samp_sc_cpn_vs_ploidy_present = rep(NA, nr),
    # Patient plotting
    pt_seg_parallel_MSAI_strict_gain = rep(NA, nr),
    pt_seg_parallel_MSAI_gain_or_amp = rep(NA, nr),
    pt_seg_parallel_MSAI_strict_amp = rep(NA, nr),
    pt_seg_parallel_MSAI_loss = rep(NA, nr),
    pt_seg_parallel_MSAI_LOH = rep(NA, nr),
    pt_seg_parallel_MSAI_loss_or_LOH = rep(NA, nr),
    pt_seg_cl_loss_LOH_no_MSAI = rep(NA, nr),
    pt_seg_cl_gain_amp_no_MSAI = rep(NA, nr),
    pt_seg_sc_loss_LOH_incl_MSAI = rep(NA, nr),
    pt_seg_sc_gain_amp_incl_MSAI = rep(NA, nr),
    # sample summary data that could go into a different dataframe
    orig_WGII = rep(NA, nr),
    mphase_WGII = rep(NA, nr),
    within_samp_cl_WGII = rep(NA, nr),
    within_samp_sc_WGII = rep(NA, nr),
    logr_data_quality = rep(NA, nr),
    TRACERx_het_measure = rep(NA, nr),
    mphase_het_measure = rep(NA, nr),
    # GD = rep(NA, nr),#don't use until we are actually using GD input for something.
    cpn_fit_quality = rep(NA, nr),
    # Purity adjusted thresholds for calling gains/losses
    purity_adj_logr_thresh_gain = rep(NA, nr),
    purity_adj_logr_thresh_loss = rep(NA, nr),
    purity_adj_logr_thresh_deep_loss = rep(NA, nr),
    purity_adj_logr_thresh_amp = rep(NA, nr),
    stringsAsFactors = FALSE
  )

  if (make.row.names.seg.name.alone == TRUE) {
    # Todo, change logic here to something more sensible... this should only be allowed when there's only one sample
    row.names.vec <- seg.name.vec
  } else {
    row.names.vec <- c()
    for (samp in sample.vec) {
      row.names.vec <- c(row.names.vec, paste(samp, seg.name.vec, sep = ":"))
    }
  }

  rownames(all.seg.df) <- row.names.vec
  all.seg.df
}


addOriginalSegmentInfoToFinalSegmentationNew <- function(all.seg.df, all.orig.samp.segs.anno.df) {

  # Really assuming order is correct here!
  all.seg.df$overlap_w_n_orig_samps <- all.orig.samp.segs.anno.df$overlap_w_n_samps
  all.seg.df$num_segs_maps_to_in_orig_samp <- all.orig.samp.segs.anno.df$num_segs_maps_to_in_samp
  all.seg.df$in_intersection_w_orig_segs <- all.orig.samp.segs.anno.df$in_intersection
  all.seg.df$in_disjoint_set_w_orig_segs <- all.orig.samp.segs.anno.df$in_disjoin_only
  all.seg.df$not_present_in_orig_segs <- all.orig.samp.segs.anno.df$in_arm_setdiff
  all.seg.df$orig_samp_seg_name <- all.orig.samp.segs.anno.df$uniq_id
  all.seg.df$orig_cpnA_raw <- all.orig.samp.segs.anno.df$nAraw
  all.seg.df$orig_cpnB_raw <- all.orig.samp.segs.anno.df$nBraw
  all.seg.df$orig_nMinor <- all.orig.samp.segs.anno.df$nMinor
  all.seg.df$orig_nMajor <- all.orig.samp.segs.anno.df$nMajor
  all.seg.df$orig_cnTotal <- all.orig.samp.segs.anno.df$cnTotal
  all.seg.df$orig_cpn_call_vs_ploidy <- all.orig.samp.segs.anno.df$cpn_call
  all.seg.df$orig_seg_in_whole_arm_event <- all.orig.samp.segs.anno.df$arm_event_present
  all.seg.df$orig_seg_is_arm_or_focal_event <- all.orig.samp.segs.anno.df$arm_or_focal_event
  all.seg.df$orig_arm_seg_present_in <- all.orig.samp.segs.anno.df$within_chr_arm
  all.seg.df$chr_arm <- all.orig.samp.segs.anno.df$chr_arm

  all.seg.df
}
