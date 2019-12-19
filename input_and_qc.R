

cleanAndErrorCheckSNPs <- function(all.samp.SNP.list, min.cons.seg.df, allow.logr.only.flag = FALSE) {
  # TODO add check that every single all.samp.SNP.list has the same SNPs..
  # Note error checking of removing NAs from the min.cons.seg.df should already have been performed
  min.cpn.vals.df <- min.cons.seg.df[, 4:ncol(min.cons.seg.df)]
  seg.name.vec <- paste(min.cons.seg.df[, 1], min.cons.segs.df[, 2], min.cons.segs.df[, 3], sep = ":")
  sample.vec <- colnames(min.cpn.vals.df)
  working.all.samp.SNP.list <- list()
  all.samps.prob.SNPs.vec <- c()

  for (sample in sample.vec) {
    samp.prob.SNPs.vec <- c()

    for (seg.idx in 1:nrow(min.cpn.vals.df)) {
      samp.SNP.df <- all.samp.SNP.list[[sample]]

      samp.seg.logR.BAF.df <- samp.SNP.df[ samp.SNP.df$chrom == min.cons.seg.df[seg.idx, 1]
      & samp.SNP.df$pos >= min.cons.seg.df[seg.idx, 2]
      & samp.SNP.df$pos <= min.cons.seg.df[seg.idx, 3], , drop = FALSE]

      # Todo,
      # There's one more case we need to consider that a only on or some samples might have a logr SNP or BAF SNP as NA but the others don't, if this were the case
      # We'd want to remove these SNPs
      # Todo,
      # Another case is that we might have any of the SNPs removed might be the start or end position of the joint segs, if this is the case it might break downstream
      # analyses, we should stop the analysis if this is the case and say clean up the input data.
      working.all.samp.SNP.list[[sample]][[seg.name.vec[seg.idx]]] <- samp.seg.logR.BAF.df
    }

    # all.samps.prob.SNPs.vec <- c(all.samps.prob.SNPs.vec, samp.prob.SNPs.vec)
  }

  # all.samps.prob.SNPs.vec <- unique(all.samps.prob.SNPs.vec)
  # print("Number of problem SNPs being removed:")
  # print(length(all.samps.prob.SNPs.vec))
  # for ( sample in sample.vec ) {

  #   for (seg.idx in 1:nrow(min.cpn.vals.df)) {

  #     working.all.samp.SNP.list[[sample]][[seg.name.vec[seg.idx]]] <- working.all.samp.SNP.list[[sample]][[seg.name.vec[seg.idx]]][ ! rownames(working.all.samp.SNP.list[[sample]][[seg.name.vec[seg.idx]]])%in%all.samps.prob.SNPs.vec, ,drop=FALSE]

  #   }#end for (seg.idx in 1:nrow(min.cpn.vals.df))

  # }#end for ( sample in sample.vec ) {

  return(working.all.samp.SNP.list)
}


getChromosomeArmGrange <- function(chr, chr.arm, chrom.info.df) {
  if (chr.arm == "p") {
    arm.df <- data.frame(
      chrom = chr,
      start = chrom.info.df[chrom.info.df$chrom == chr, ]$p_start,
      end = chrom.info.df[chrom.info.df$chrom == chr, ]$p_end,
      stringsAsFactors = FALSE
    )
  } else if (chr.arm == "q") {
    arm.df <- data.frame(
      chrom = chr,
      start = chrom.info.df[chrom.info.df$chrom == chr, ]$q_start,
      end = chrom.info.df[chrom.info.df$chrom == chr, ]$q_end,
      stringsAsFactors = FALSE
    )
  }

  arm.gr <- makeGRangesFromDataFrame(arm.df,
    keep.extra.columns = FALSE,
    ignore.strand = FALSE,
    seqinfo = NULL,
    seqnames.field = c(
      "seqnames", "seqname",
      "chromosome", "chrom",
      "chr", "chromosome_name",
      "seqid"
    ),
    start.field = c("Start_Position", "start"),
    end.field = c("End_Position", "end", "stop"),
    strand.field = "strand",
    starts.in.df.are.0based = FALSE
  )

  return(arm.gr)
}


getAllPatientSamplesAndIntersectionGRangesList <- function(sample.vec, f.joint.segs.df) {

  # Now what about intersection of all samples what proportion of the genome to they cover:
  samp.gr.list <- list()
  overall.intersection.gr <- NULL
  for (samp in sample.vec) {
    samp.gr <- makeGRangesFromDataFrame(f.joint.segs.df[f.joint.segs.df$sample == samp, ],
      keep.extra.columns = FALSE,
      ignore.strand = FALSE,
      seqinfo = NULL,
      seqnames.field = c(
        "seqnames", "seqname",
        "chromosome", "chrom",
        "chr", "chromosome_name",
        "seqid"
      ),
      start.field = c("Start_Position", "start", "startpos"),
      end.field = c("End_Position", "end", "stop", "endpos"),
      strand.field = "strand",
      starts.in.df.are.0based = FALSE
    )
    samp.gr.list[[samp]] <- samp.gr

    if (is.null(overall.intersection.gr)) {
      overall.intersection.gr <- samp.gr
    } else {
      overall.intersection.gr <- intersect(overall.intersection.gr, samp.gr)
    }
  }

  samp.gr.list[["overall_intersection"]] <- overall.intersection.gr

  return(samp.gr.list)
}


compareOverlapsBetweenSamples <- function(sample.vec, f.joint.segs.df, hg19.info.df, chrs.to.use.vec = CHRS_TO_USE_VEC) {

  # For each sample how much of the chromosome arm
  # chrs.to.use.vec <- 1:22
  # Subset the chromos
  chr.info.df <- hg19.info.df[hg19.info.df$chrom %in% chrs.to.use.vec, ]
  p.name.vec <- chr.info.df$p_name
  q.name.vec <- chr.info.df$q_name
  arm.vec <- c(rbind(p.name.vec, q.name.vec))
  p.q.vec <- c("p", "q")

  samp.gr.list <- getAllPatientSamplesAndIntersectionGRangesList(sample.vec, f.joint.segs.df)
  # So now let's make an arm coverage data.frame
  arm.coverage.df <- data.frame(matrix(vector(), length(arm.vec), (length(sample.vec) + 2),
    dimnames = list(c(), c("arm", sample.vec, "overall_intersection"))
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
  )
  arm.coverage.df$arm <- arm.vec

  # So now we want to take the intersection and compare the sizes of the

  for (samp in names(samp.gr.list)) {

    # Get sample granges:
    samp.gr <- samp.gr.list[[samp]]
    for (chr in chrs.to.use.vec) {
      for (arm in p.q.vec) {
        chr.arm.gr <- getChromosomeArmGrange(chr, arm, chr.info.df)
        chr.arm.overlap.gr <- intersect(samp.gr, chr.arm.gr)
        arm.coverage.df[arm.coverage.df$arm == paste0(chr, arm), samp] <- sum(as.numeric(width(chr.arm.overlap.gr))) / as.numeric(width(chr.arm.gr))
      }
    }
  }

  arm.coverage.df$diff_from_max <- apply(arm.coverage.df[, 2:ncol(arm.coverage.df)], 1, max) - arm.coverage.df$overall_intersection

  # Report max and union

  return(arm.coverage.df)
}


joinTwoGrangesIntoOne <- function(gr.one, gr.two) {

  # Removee any metadata columns
  mcols(gr.one) <- NULL
  mcols(gr.two) <- NULL

  chr <- unique(as.character(seqnames(gr.one)), as.character(seqnames(gr.two)))
  comb.grs <- c(gr.two, gr.one)
  s.comb.grs <- sort(comb.grs)
  if (length(chr) > 1) {
    message("More than one chromsome given as input, returning input unjoined")
    return(s.comb.grs)
  }

  new.start <- start(s.comb.grs)[1]
  new.end <- end(s.comb.grs)[length(end(s.comb.grs))]

  joined.gr <- GRanges(
    seqnames = chr,
    ranges = IRanges(new.start, new.end)
  )

  return(joined.gr)
}


removeSingleBasePairSegmentsFromIntersectionGranges <- function(all.grs, min.width = 1) {

  # Calculate width:
  all.grs$width <- width(all.grs)
  # Get grs that are one
  width.one.grs <- all.grs[ all.grs$width <= min.width ]

  working.no.ones.grs <- all.grs
  mcols(working.no.ones.grs) <- NULL

  if (length(width.one.grs) == 0) {
    return(working.no.ones.grs)
  }

  for (idx in 1:length(width.one.grs)) {
    gr.ante.idx <- follow(width.one.grs[idx], working.no.ones.grs)
    gr.post.idx <- precede(width.one.grs[idx], working.no.ones.grs)
    merge.poss.idxs <- c(gr.ante.idx, gr.post.idx)
    non.na.merge.poss.idxs <- merge.poss.idxs[!is.na(merge.poss.idxs)]

    # Check thats
    if (length(non.na.merge.poss.idxs) == 0) {
      message("No ante or post segments")
      next()
    }
    ante.post.grs <- working.no.ones.grs[non.na.merge.poss.idxs]
    nearest.gr <- ante.post.grs[nearest(width.one.grs[idx], ante.post.grs)]
    joined.gr <- joinTwoGrangesIntoOne(nearest.gr, width.one.grs[idx])
    all.idx.vec <- 1:length(working.no.ones.grs)
    hits.obj <- findOverlaps(joined.gr, working.no.ones.grs)
    hit.idx.vec <- subjectHits(hits.obj)
    non.overlapping.idx.vec <- all.idx.vec[ !all.idx.vec %in% hit.idx.vec]
    non.overlapping.grs <- working.no.ones.grs[non.overlapping.idx.vec]

    working.no.ones.grs <- sort(c(non.overlapping.grs, joined.gr))
  }

  return(working.no.ones.grs)
}


getChromosomeArmsGRangeList <- function(current.chr, chrom.info.df) {
  arm.df.list <- list()
  p.arm.df <- data.frame(
    chrom = current.chr,
    start = chrom.info.df[chrom.info.df$chrom == current.chr, ]$p_start,
    end = chrom.info.df[chrom.info.df$chrom == current.chr, ]$p_end,
    stringsAsFactors = FALSE
  )

  count <- 1

  if (chrom.info.df[chrom.info.df$chrom == current.chr, ]$acrocentric == FALSE) {
    arm.df.list[[count]] <- makeGRangesFromDataFrame(p.arm.df,
      keep.extra.columns = FALSE,
      ignore.strand = FALSE,
      seqinfo = NULL,
      seqnames.field = c(
        "seqnames", "seqname",
        "chromosome", "chrom",
        "chr", "chromosome_name",
        "seqid"
      ),
      start.field = c("Start_Position", "start"),
      end.field = c("End_Position", "end", "stop"),
      strand.field = "strand",
      starts.in.df.are.0based = FALSE
    )

    count <- count + 1
  }

  q.arm.df <- data.frame(
    chrom = current.chr,
    start = chrom.info.df[chrom.info.df$chrom == current.chr, ]$q_start,
    end = chrom.info.df[chrom.info.df$chrom == current.chr, ]$q_end,
    stringsAsFactors = FALSE
  )

  arm.df.list[[count]] <- makeGRangesFromDataFrame(q.arm.df,
    keep.extra.columns = FALSE,
    ignore.strand = FALSE,
    seqinfo = NULL,
    seqnames.field = c(
      "seqnames", "seqname",
      "chromosome", "chrom",
      "chr", "chromosome_name",
      "seqid"
    ),
    start.field = c("Start_Position", "start"),
    end.field = c("End_Position", "end", "stop"),
    strand.field = "strand",
    starts.in.df.are.0based = FALSE
  )

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
# Params:
# query: one or more granges/iranges that you wish to compare to the subject
# subject: one or more grangaes/ranges that you wish to be compared with the query
# output: specifies the output type either "gr" for genomic range or "df" for dataframe
# Returns a dataframe/grange with following columns/fields:
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
# Also returns all the metadata in the subject, ie original columns holding non-genomic position information
# Todo add widths for everything
getOverlapGRanges <- function(query, subject, output = "gr", ...) {
  require(GenomicRanges)
  require(IRanges)

  ## Input check
  if (!((class(query) == "GRanges" & class(subject) == "GRanges") | (class(query) == "IRanges" & class(subject) == "IRanges"))) {
    stop("Query and subject need to be of same class, either GRanges or IRanges!")
  }

  ## Find overlapping ranges
  if (class(query) == "GRanges") {
    seqlengths(query) <- rep(NA, length(seqlengths(query)))
    seqlengths(subject) <- rep(NA, length(seqlengths(subject)))
  }
  overlap.index.mat <- as.matrix(findOverlaps(query, subject, ...))

  query <- query[overlap.index.mat[, 1]]
  subject <- subject[overlap.index.mat[, 2]]
  overlap.mat <- cbind(query_start = start(query), query_end = end(query), subject_start = start(subject), subject_end = end(subject))

  ## Pre-queries for overlaps
  startup <- overlap.mat[, "subject_start"] < overlap.mat[, "query_start"]
  enddown <- overlap.mat[, "subject_end"] > overlap.mat[, "query_end"]
  startin <- overlap.mat[, "subject_start"] >= overlap.mat[, "query_start"] & overlap.mat[, "subject_start"] <= overlap.mat[, "query_end"]
  endin <- overlap.mat[, "subject_end"] >= overlap.mat[, "query_start"] & overlap.mat[, "subject_end"] <= overlap.mat[, "query_end"]

  ## Overlap types
  olup <- startup & endin
  oldown <- startin & enddown
  inside <- startin & endin
  contained <- startup & enddown

  ## Overlap types in one vector
  overlap_type <- rep("", length(overlap.mat[, "query_start"]))
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
  oldf <- data.frame(query_hit_idx = overlap.index.mat[, 1], subject_hit_idx = overlap.index.mat[, 2], overlap.mat, overlap_start, overlap_end, overlap_length, overlap_prop_query, overlap_prop_subject, overlap_type)
  if (class(query) == "GRanges") {
    oldf <- cbind(chrom = as.character(seqnames(query)), oldf)
  }
  if (output == "df") {
    return(cbind(oldf, as.data.frame(mcols(subject))))
  }
  if (output == "gr") {
    if (class(query) == "GRanges") {
      elementMetadata(query) <- cbind(as.data.frame(elementMetadata(query)), oldf, mcols(subject))
    }
    if (class(query) == "IRanges") { # Todo add the metadata in th eiranges
      query <- GRanges(seqnames = Rle(rep("dummy", length(query))), ranges = IRanges(start = oldf[, "query_start"], end = oldf[, "query_end"]), strand = Rle(strand(rep("+", length(query)))), oldf)
    }
    return(query)
  }
}


# Subset snp data frames to keep only those snps that are present in all of the
# others
keep_common_snps <- function(all_snps) {
  message("Subsetting to shared SNPs only")

  samples <- names(all_snps)

  # Identify all SNPs that are present in every sample
  snp.pos.vec <- row.names(all_snps[[1]])
  for (sample in samples[-1]) {
    snp.pos.vec <- snp.pos.vec[snp.pos.vec %in% row.names(all_snps[[sample]])]
  }

  # Remove SNP positions that are not in all samples
  for (sample in samples) {
    all_snps[[sample]] <- all_snps[[sample]][snp.pos.vec, ]
  }

  return(all_snps)
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
snps_set_bafsegmented <- function(segs, all_snps) {
  message("Calculating segmented BAF values.")
  warning("Centering BAF values at 0.475 rather than 0.5 to account for reference bias. This value was chosen based on mean BAF across all samples in this cohort."

  sample_names <- names(all_snps)
  for (samp_idx in seq_along(all_snps)) {
    samp_segs <- segs[segs$sample == sample_names[[samp_idx]], ]
    samp_data <- all_snps[[samp_idx]]
    for (seg_idx in seq_len(nrow(samp_segs))) {
      seg <- samp_segs[seg_idx, ]

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
        center_baf <- 0.475  # not 0.5, to account for reference bias
        mean_baf <- center_baf - mean(abs(samp_data[seg_data_idx, "BAF"] - center_baf), na.rm = TRUE)
        all_snps[[samp_idx]][seg_data_idx, "BAFsegmented"] <- mean_baf
      }
    }
  }

  return(all_snps)
}


snps_set_logrsegmented <- function(segs, all_snps) {
  message("Calculating segmented logR values.")

  sample_names <- names(all_snps)
  for (samp_idx in seq_along(all_snps)) {
    samp_segs <- segs[segs$sample == sample_names[[samp_idx]], ]
    samp_data <- all_snps[[samp_idx]]
    for (seg_idx in seq_len(nrow(samp_segs))) {
      seg <- samp_segs[seg_idx, ]

      # Segments include both the start and end positions
      seg_data_idx <- (samp_data$chrom == seg$chr &
        samp_data$pos >= seg$startpos &
        samp_data$pos <= seg$endpos)

      mean_logr <- mean(samp_data[seg_data_idx, "logR"], na.rm = TRUE)
      all_snps[[samp_idx]][seg_data_idx, "logRsegmented"] <- mean_logr
    }

    # Deal with NAs
    logr_na <- is.na(all_snps[[samp_idx]]$logRsegmented)
    if (any(logr_na)) {
      message("Found NA values in logR segmented, setting to 0.0")
      all_snps[[samp_idx]][logr_na, "logRsegmented"] <- 0.0
    }
  }

  return(all_snps)
}


cleanAndErrorCheckSNPsParallel <- function(all.samp.SNP.list, mcs.seg.to.use.df, samples, allow.logr.only.flag = FALSE, data.type = "MEDICC") {
  message("Clean and error check SNPs (parallel version)")
  # Note error checking of removing NAs from the min.cons.seg.df should already
  # have been performed

  min.cpn.vals.df <- mcs.seg.to.use.df[, 4:ncol(mcs.seg.to.use.df)]
  seg.name.vec <- paste(mcs.seg.to.use.df[, "Chromosome"], mcs.seg.to.use.df[, "Start_Position"], mcs.seg.to.use.df[, "End_Position"], sep = ":")
  working.all.samp.SNP.list <- list()
  all.samps.prob.SNPs.vec <- c()

  # Make vectors for sample vec:
  samp.seg.vec <- rep(seg.name.vec, length(samples))
  sample.rep.vec <- c()
  for (samp in samples) {
    working.all.samp.SNP.list[[samp]] <- list()
    sample.rep.vec <- c(sample.rep.vec, rep(samp, length(seg.name.vec)))
  }

  iterable.min.cons.seg.df <- cbind(do.call("rbind", replicate(length(samples), mcs.seg.to.use.df[, c("Chromosome", "Start_Position", "End_Position")], simplify = FALSE)), sample = sample.rep.vec)
  library("gdata")
  iterable.min.cons.seg.df <- drop.levels(iterable.min.cons.seg.df)
  itx <- iter(iterable.min.cons.seg.df, by = "row")
  message("Starting parallel assignation of snps to segs:")
  # Create class which holds multiple results for each loop iteration.
  # Each loop iteration populates two properties: $result1 and $result2.
  # For a great tutorial on S3 classes, see:
  # http://www.cyclismo.org/tutorial/R/s3Classes.html#creating-an-s3-class
  MultiResultClass <- function(samp_seg_name = NULL, samp_seg_logR_BAF_df = NULL, samp_prob_SNPs_vec = NULL) {
    me <- list(
      samp_seg_name = samp_seg_name,
      samp_seg_logR_BAF_df = samp_seg_logR_BAF_df,
      samp_prob_SNPs_vec = samp_prob_SNPs_vec
    )
    # Set the name for the class
    class(me) <- append(class(me), "MultiResultClass")
    return(me)
  }
  sites.in.each.samp.seg.list <- foreach(it.min.cons.seg.df = itx, .combine = c) %dopar% {
    samp.SNP.df <- all.samp.SNP.list[[ it.min.cons.seg.df[, "sample"]]]

    if ("chr" %in% colnames(samp.SNP.df)) {
      colnames(samp.SNP.df)[colnames(samp.SNP.df) == "chr"] <- "chrom"
    }
    samp.seg.logR.BAF.df <- samp.SNP.df[ samp.SNP.df$chrom == it.min.cons.seg.df[, "Chromosome"] & samp.SNP.df$pos >= it.min.cons.seg.df[, "Start_Position"] & samp.SNP.df$pos <= it.min.cons.seg.df[, "End_Position"], , drop = FALSE]
    samp.prob.SNPs.vec <- c()

    if (data.type == "TRACERx_WGS") {
      samp.prob.SNPs.vec <- c()
    } else if (allow.logr.only.flag == FALSE & data.type == "TRACERx_exome") {
      samp.prob.SNPs.vec <- rownames(samp.seg.logR.BAF.df[is.na(samp.seg.logR.BAF.df$BAFsegmented), , drop = FALSE ])
      samp.seg.logR.BAF.df <- samp.seg.logR.BAF.df[samp.seg.logR.BAF.df$BAF != 2, , drop = FALSE]
      samp.seg.logR.BAF.df <- samp.seg.logR.BAF.df[!is.na(samp.seg.logR.BAF.df$BAFsegmented), , drop = FALSE]
    } else if (allow.logr.only.flag == FALSE & data.type == "MEDICC") {
      samp.prob.SNPs.vec <- c(samp.prob.SNPs.vec, rownames(samp.seg.logR.BAF.df[samp.seg.logR.BAF.df$BAF == 2, , drop = FALSE ]))
      samp.seg.logR.BAF.df <- samp.seg.logR.BAF.df[as.numeric(samp.seg.logR.BAF.df$BAF) != 2, , drop = FALSE]
      samp.seg.logR.BAF.df <- samp.seg.logR.BAF.df[!is.na(samp.seg.logR.BAF.df$BAFsegmented), , drop = FALSE]
      # colnames(samp.seg.logR.BAF.df)[colnames(samp.seg.logR.BAF.df) == "chr"] <- "chrom"
    } else if (allow.logr.only.flag == FALSE & data.type == "ASCAT") {
      samp.prob.SNPs.vec <- rownames(samp.seg.logR.BAF.df)[is.na(samp.seg.logR.BAF.df$BAF)]
      samp.seg.logR.BAF.df <- samp.seg.logR.BAF.df[!is.na(samp.seg.logR.BAF.df$BAF), , drop = FALSE]
    } else if (allow.logr.only.flag == FALSE & data.type == "sequenza") {
      samp.prob.SNPs.vec <- rownames(samp.seg.logR.BAF.df)[is.na(samp.seg.logR.BAF.df$BAF)]
      samp.seg.logR.BAF.df <- samp.seg.logR.BAF.df[!is.na(samp.seg.logR.BAF.df$BAF), , drop = FALSE]
    } else {
      samp.prob.SNPs.vec <- c()
    }

    result <- MultiResultClass(
      samp_seg_name = paste0(it.min.cons.seg.df[, "sample"], "@", it.min.cons.seg.df[, "Chromosome"], ":", it.min.cons.seg.df[, "Start_Position"], ":", it.min.cons.seg.df[, "End_Position"]),
      samp_seg_logR_BAF_df = samp.seg.logR.BAF.df,
      samp_prob_SNPs_vec = samp.prob.SNPs.vec
    )
    result.list <- list(result)
    names(result.list) <- paste0(it.min.cons.seg.df[, "sample"], "@", it.min.cons.seg.df[, "Chromosome"], ":", it.min.cons.seg.df[, "Start_Position"], ":", it.min.cons.seg.df[, "End_Position"])
    result.list
  }

  all.samps.prob.SNPs.vec <- unique(unlist(sapply(sites.in.each.samp.seg.list, function(x) x[["samp_prob_SNPs_vec"]])))
  message("Number of SNPs with NA BAF = ", length(all.samps.prob.SNPs.vec), ". These will be removed.")

  if (length(all.samps.prob.SNPs.vec) > 0) {
    samp.seg.itx <- iter(sites.in.each.samp.seg.list)

    clean.sites.in.each.samp.seg.list <- foreach(samp.seg.mrc = samp.seg.itx, .combine = c) %dopar% {
      seg.name <- strsplit(samp.seg.mrc$samp_seg_name, "@")[[1]][2]
      samp.seg.mrc.list <- sites.in.each.samp.seg.list[grep(seg.name, names(sites.in.each.samp.seg.list), value = TRUE)]
      samp.seg.prob.snp.vec <- unique(unlist(sapply(samp.seg.mrc.list, function(x) x[["samp_prob_SNPs_vec"]])))
      samp.seg.logR.BAF.df <- samp.seg.mrc$samp_seg_logR_BAF_df[ !rownames(samp.seg.mrc$samp_seg_logR_BAF_df) %in% samp.seg.prob.snp.vec, , drop = FALSE ]
      result <- MultiResultClass(
        samp_seg_name = samp.seg.mrc$samp_seg_name,
        samp_seg_logR_BAF_df = samp.seg.logR.BAF.df,
        samp_prob_SNPs_vec = c()
      )
      result.list <- list(result)
      names(result.list) <- samp.seg.mrc$samp_seg_name
      result.list
    }
  } else {
    clean.sites.in.each.samp.seg.list <- sites.in.each.samp.seg.list
  }

  for (samp.seg.name in names(clean.sites.in.each.samp.seg.list)) {
    samp <- strsplit(samp.seg.name, "@")[[1]][1]
    seg.name <- strsplit(samp.seg.name, "@")[[1]][2]
    working.all.samp.SNP.list[[samp]][[seg.name]] <- clean.sites.in.each.samp.seg.list[[samp.seg.name]]$samp_seg_logR_BAF_df
  }

  return(working.all.samp.SNP.list)
}
