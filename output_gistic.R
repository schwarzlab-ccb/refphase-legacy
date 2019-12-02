# Now for each segment in the patient for the gains we want the maximum copynumber and for losses the minimum
collapseSegsByPatientAndMinMaxMean <- function(all.segs.df,
                                               markers.to.use.as.default = "all_sites",
                                               cpn.to.use = "total_raw_cpn",
                                               min.max.mean.to.use = "mean") {

  # There should be no NA segmens left after the filterling
  if (!min.max.mean.to.use %in% c("mean", "max", "min")) {
    print("min.max.mean.to.use")
    print(min.max.mean.to.use)
    stop("Error in function to use to calculate copy number.")
  }
  if (!cpn.to.use %in% c("total_raw_cpn", "ph_hz_total_raw_cpn", "non_hz_total_raw_cpn", "min_raw_cpn_for_GISTIC_loss", "seg_cn_for_gistic")) {
    print("cpn.to.use")
    print(cpn.to.use)
    stop("Error copy number column chosen.")
  }

  function.to.apply <- match.fun(min.max.mean.to.use)
  all.segs.df$samp_seg_name <- paste0(all.segs.df$sample, "@", all.segs.df$seg_name)
  all.segs.df$patient_seg_name <- paste0(all.segs.df$patient, "@", all.segs.df$seg_name)
  all.segs.raw.cpn.df <- all.segs.df[, c("total_raw_cpn", "ph_hz_total_raw_cpn", "non_hz_total_raw_cpn", "min_raw_cpn_for_GISTIC_loss", "seg_cn_for_gistic")]
  raw.cpn.segs.by.pt.df <- aggregate(all.segs.raw.cpn.df, by = list(all.segs.df$patient_seg_name), FUN = function.to.apply)
  cpn.to.use.vec <- raw.cpn.segs.by.pt.df[, cpn.to.use]


  # all.segs.df$patient_seg_name <- paste0(all.segs.df$patient, "@", all.segs.df$segName)
  # Use these idxs to get the markers present.
  all.segs.sites.df <- all.segs.df[match(raw.cpn.segs.by.pt.df$Group.1, all.segs.df$patient_seg_name), c("num_ph_hz_sites", "num_un_ph_hz_sites", "num_non_hz_sites")]
  all.sites.vec <- rowSums(all.segs.sites.df, na.rm = TRUE)
  hz.sites.vec <- rowSums(all.segs.sites.df[, c("num_ph_hz_sites", "num_un_ph_hz_sites")], na.rm = TRUE)
  non.hz.sites.vec <- rowSums(all.segs.sites.df[, c("num_non_hz_sites"), drop = FALSE], na.rm = TRUE)

  # Now ploidy to calculate the CN...
  all.samps.ploidy.vec <- all.segs.df[match(all.segs.df$patient_seg_name, raw.cpn.segs.by.pt.df$Group.1), "sample_ploidy"]
  processed.cpn.to.use.vec <- as.numeric(cpn.to.use.vec)
  processed.cpn.to.use.vec[processed.cpn.to.use.vec <= 0] <- 0.00001

  depth.ratio.vec <- log2(processed.cpn.to.use.vec) - 1
  # depth.ratio.vec <- log2(as.numeric(all.samps.ploidy.vec) * as.numeric(cpn.to.use.vec)) - 1
  if (markers.to.use.as.default == "all_sites") {
    default.sites.vec <- all.sites.vec
  } else if (markers.to.use.as.default == "hz_sites") {
    default.sites.vec <- hz.sites.vec
  } else if (markers.to.use.as.default == "non_hz_sites") {
    default.sites.vec <- non.hz.sites.vec
  } else {
    print("markers.to.use.as.default")
    print(markers.to.use.as.default)
    stop("Error in markers to use as default selection.")
  }
  # So now we want to make a GISTIC dataframe:
  # (1)  Sample           (sample name)
  # (2)  Chromosome  (chromosome number)
  # (3)  Start Position  (segment start position, in bases)
  # (4)  End Position   (segment end position, in bases)
  # (5)  Num markers      (number of markers in segment)
  # (6)  Seg.CN       (log2() -1 of copy number)
  sample.vec <- sapply(raw.cpn.segs.by.pt.df$Group.1, function(x) strsplit(split = "@", x)[[1]][[1]])
  seg.name.vec <- sapply(raw.cpn.segs.by.pt.df$Group.1, function(x) strsplit(split = "@", x)[[1]][[2]])
  chr.vec <- as.numeric(sapply(seg.name.vec, function(x) strsplit(split = ":", x)[[1]][[1]]))
  start.pos.vec <- as.numeric(sapply(seg.name.vec, function(x) strsplit(split = ":", x)[[1]][[2]]))
  end.pos.vec <- as.numeric(sapply(seg.name.vec, function(x) strsplit(split = ":", x)[[1]][[3]]))
  print(length(sample.vec))
  print(length(chr.vec))
  print(length(start.pos.vec))
  print(length(default.sites.vec))
  print(length(depth.ratio.vec))

  GISTIC.extended.input.df <- data.frame(
    Sample = sample.vec,
    Chromosome = as.numeric(chr.vec),
    Start_Position = as.numeric(start.pos.vec),
    End_Position = as.numeric(end.pos.vec),
    Num_Markers = as.numeric(default.sites.vec),
    Seg.CN = as.numeric(depth.ratio.vec),
    stringsAsFactors = FALSE
  )

  GISTIC.extended.input.df <- GISTIC.extended.input.df[with(GISTIC.extended.input.df, order(Sample, Chromosome, Start_Position)), ]

  return(GISTIC.extended.input.df)
}
