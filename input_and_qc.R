

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
  # message("Reading in ascat-format snp data from the following files")
  # print(sample_table)

  for (sample in samples) {
    sample_info <- sample_table[sample, ]
    message("Starting to read sample", sample_info$samp_ascat_bc_path)

    if (grepl("RData", sample_info$samp_ascat_bc_path)) {
      message("ascat.bc object in RData file, loading...")
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
    message("# of heterozygous SNPs: ", sum(is_heterozygous))

    # Mask non-heterozygous SNPs with NA
    baf <- ascat_bc$Tumor_BAF[[1]]
    baf[!is_heterozygous] <- NA

    chrom <- chrom2integer(ascat_bc$SNPpos$chr)

    samp_data <- data.frame(
      chrom = chrom,
      pos = ascat_bc$SNPpos$pos,
      BAF = baf,
      logR = ascat_bc$Tumor_LogR[[1]],
      germline_zygosity = ifelse(is_heterozygous, "het", "hom"),
      logRsegmented = as.numeric(ascat_bc$Tumor_LogR_segmented),
      snp_pos = paste0(chrom, "_", ascat_bc$SNPpos$pos),
      row.names = paste0(chrom, "_", ascat_bc$SNPpos$pos),
      stringsAsFactors = FALSE
    )

    all_snps[[sample]] <- samp_data
  }

  # So, in TRACERx we have the issue that ASCAT, when run with exactly the same
  # input ascat occasionally does not seem to have all the same positions. So,
  # in order to get xphase to run I have to make sure that the same
  # snps/positions are kept for all samples. Given the way M/XPhase works
  # (finding the nearest SNP to a breakpoint) this should be fine even if we
  # get rid of some positions that are the start or end indices of samples'
  # segments.
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
ascat_set_bafsegmented <- function(segs, all_snps) {
  message("Calculating segmented BAF values.")

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


readInMEDICCInputFormatParallel <- function(snp.file.df, processed.orig.cpn.segs.df, sample.vec) {
  # TODO add check that every single all.samp.SNP.list has the same SNPs..
  # Note error checking of removing NAs from the min.cons.seg.df should already have been performed
  # registerDoParallel(cores = num.cores)
  # print(foreach::getDoParWorkers())
  # print(getDoParName())
  # print(getDoParVersion())
  message("Calling readInMEDICCInputFormatParallel")
  message("Reading in sites data for:")
  print(sample.vec)
  message("Starting parallel work:")

  num.pre.parallel.samps <- length(unique(sample.vec))
  itx <- iter(sample.vec)
  all.samp.SNP.df <- foreach(samp = itx, .combine = rbind) %dopar% {
    ASCAT.summary <- data.frame(
      chr = as.numeric(snp.file.df$Chr),
      pos = as.numeric(snp.file.df$Pos),
      BAF = as.numeric(snp.file.df[, paste0(gsub("-", "\\.", samp), "_BAF")]),
      BAFsegmented = rep(NA, length(snp.file.df$Chr)),
      logR = as.numeric(snp.file.df[, paste0(gsub("-", "\\.", samp), "_LogR")]),
      logRsegmented = rep(NA, length(snp.file.df$Chr)),
      stringsAsFactors = FALSE
    )
    ASCAT.summary$sample <- samp

    if (nrow(processed.orig.cpn.segs.df) == 0) {
      stop("processed.orig.cpn.segs.df has no rows...")
    }
    seg <- processed.orig.cpn.segs.df[processed.orig.cpn.segs.df$sample == samp, ]
    if (nrow(processed.orig.cpn.segs.df) == 0) {
      stop("processed.orig.cpn.segs.df subsetted by sample has no rows...")
    }

    for (j in 1:nrow(seg)) {
      # Are there any BAF values that aren't set at 2 (and thus shouldn't be
      # used for allele specific cpn) present in this segment?
      if (nrow(ASCAT.summary[ASCAT.summary$chr == seg[j, ]$chr
      & ASCAT.summary$pos >= seg[j, ]$startpos
      & ASCAT.summary$pos <= seg[j, ]$endpos
      & ASCAT.summary$BAF != 2, ]) > 0) {

        # Is this segment balanced, if so set the BAFsegmented to 0.5
        if (seg[j, , ]$nMajor == seg[j, ]$nMinor) {
          ASCAT.summary[ASCAT.summary$chr == seg[j, ]$chr
          & ASCAT.summary$pos >= seg[j, ]$startpos
          & ASCAT.summary$pos <= seg[j, ]$endpos, ]$BAFsegmented <- 0.5
        } else {

          avg.seg.baf <- 0.5 - mean(abs(as.numeric(ASCAT.summary[ASCAT.summary$chr == seg[j, ]$chr
          & ASCAT.summary$pos >= seg[j, ]$startpos
          & ASCAT.summary$pos <= seg[j, ]$endpos
          & ASCAT.summary$BAF != 2, ]$BAF) - 0.5))

          ASCAT.summary[ASCAT.summary$chr == seg[j, ]$chr
          & ASCAT.summary$pos >= seg[j, ]$startpos
          & ASCAT.summary$pos <= seg[j, ]$endpos, ]$BAFsegmented <- avg.seg.baf
        }
      }

      # We now need to calculate the average logr of any snps that are present
      if (nrow(ASCAT.summary[ASCAT.summary$chr == seg[j, ]$chr
      & ASCAT.summary$pos >= seg[j, ]$startpos
      & ASCAT.summary$pos <= seg[j, ]$endpos, ]) > 0) {
        avg.seg.logr <- mean(ASCAT.summary[ASCAT.summary$chr == seg[j, ]$chr
        & ASCAT.summary$pos >= seg[j, ]$startpos
        & ASCAT.summary$pos <= seg[j, ]$endpos, ]$logR)

        ASCAT.summary[ASCAT.summary$chr == seg[j, ]$chr
        & ASCAT.summary$pos >= seg[j, ]$startpos
        & ASCAT.summary$pos <= seg[j, ]$endpos, ]$logRsegmented <- avg.seg.logr
      }
    }
    ASCAT.summary
  }

  num.post.parallel.samps <- length(unique(all.samp.SNP.df$sample))

  if (num.post.parallel.samps != num.pre.parallel.samps) {
    message("parallelisation is not working properly, some parallel workers are not returning jobs correctly....")
    message("Either investigate why some cores are not returning results, change to a different node or run in series.")
    message("Running in series now...")
  }

  message("Assigning to list in expected format:")
  all.samp.SNP.list <- list()
  for (samp in sample.vec) {
    message(samp)
    samp.w.segmented.info.SNP.df <- all.samp.SNP.df[all.samp.SNP.df$sample == samp, ]
    print(head(samp.w.segmented.info.SNP.df))
    samp.w.segmented.info.SNP.df <- samp.w.segmented.info.SNP.df[, c("chr", "pos", "BAF", "BAFsegmented", "logR", "logRsegmented")]
    all.samp.SNP.list[[samp]] <- samp.w.segmented.info.SNP.df
  }

  return(all.samp.SNP.list)
}


cleanAndErrorCheckSNPsParallel <- function(all.samp.SNP.list, mcs.seg.to.use.df, samples, allow.logr.only.flag = FALSE, data.type = "MEDICC") {
  # Note error checking of removing NAs from the min.cons.seg.df should already
  # have been performed

  min.cpn.vals.df <- mcs.seg.to.use.df[, 4:ncol(mcs.seg.to.use.df)]
  seg.name.vec <- paste(mcs.seg.to.use.df[, 1], mcs.seg.to.use.df[, 2], mcs.seg.to.use.df[, 3], sep = ":")
  # min.cpn.vals.df <- min.cpn.vals.df[ ! is.na(rowSums(min.cpn.vals.df)), ,drop = FALSE]
  working.all.samp.SNP.list <- list()
  all.samps.prob.SNPs.vec <- c()

  # Make vectors for sample vec:
  samp.seg.vec <- rep(seg.name.vec, length(samples))
  sample.rep.vec <- c()
  for (samp in samples) {
    working.all.samp.SNP.list[[samp]] <- list()
    sample.rep.vec <- c(sample.rep.vec, rep(samp, length(seg.name.vec)))
  }

  iterable.min.cons.seg.df <- cbind(do.call("rbind", replicate(length(samples), mcs.seg.to.use.df[, 1:3], simplify = FALSE)), sample = sample.rep.vec)
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
    samp.SNP.df <- all.samp.SNP.list[[ it.min.cons.seg.df[, 4]]]

    if ("chr" %in% colnames(samp.SNP.df)) {
      colnames(samp.SNP.df)[colnames(samp.SNP.df) == "chr"] <- "chrom"
    }
    samp.seg.logR.BAF.df <- samp.SNP.df[ samp.SNP.df$chrom == it.min.cons.seg.df[, 1]
    & samp.SNP.df$pos >= it.min.cons.seg.df[, 2]
    & samp.SNP.df$pos <= it.min.cons.seg.df[, 3], , drop = FALSE]
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
      samp_seg_name = paste0(it.min.cons.seg.df[, 4], "@", it.min.cons.seg.df[, 1], ":", it.min.cons.seg.df[, 2], ":", it.min.cons.seg.df[, 3]),
      samp_seg_logR_BAF_df = samp.seg.logR.BAF.df,
      samp_prob_SNPs_vec = samp.prob.SNPs.vec
    )
    result.list <- list(result)
    names(result.list) <- paste0(it.min.cons.seg.df[, 4], "@", it.min.cons.seg.df[, 1], ":", it.min.cons.seg.df[, 2], ":", it.min.cons.seg.df[, 3])
    result.list
  }

  all.samps.prob.SNPs.vec <- unique(unlist(sapply(sites.in.each.samp.seg.list, function(x) x[["samp_prob_SNPs_vec"]])))
  message("Number of problem SNPs being removed:", length(all.samps.prob.SNPs.vec))
  print(head(all.samps.prob.SNPs.vec))

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
read_snps_sequenza <- function(input_files, segs, samples, het_only = TRUE) {
  library(readr)

  all_snps <- list()

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
    # Skip unnecessary columns
    col_types <- "cic--dddc--c--"
    col_names <- c("chromosome", "position", "base.ref", "depth.ratio", "Af", "Bf", "zygosity.normal", "AB.normal")
    seqz_data <- read_tsv(input_file,
      col_types = col_types,
      progress = TRUE, col_names = col_names,
      skip = 1
    )

    # Only keep heterozygous positions
    if (isTRUE(het_only)) {
      germline_het <- seqz_data$zygosity.normal == "het"
      seqz_data <- seqz_data[germline_het, ]
    }

    message("Creating SNP data frame.")
    chrom <- chrom2integer(seqz_data$chromosome)
    pos <- seqz_data$position

    # Bf is always < Af, so we need to check the AB.normal to see if we need to flip it.
    # This is probably not the perfect thing to do. Look into this topic more.
    flip_baf <- nchar(seqz_data$AB.normal) > 0 & seqz_data$base.ref != substr(seqz_data$AB.normal, 1, 1)
    baf <- seqz_data$Bf
    baf[flip_baf] <- seqz_data$Af[flip_baf]

    all_snps[[sample]] <- data.frame(
      chrom = chrom,
      pos = pos,
      BAF = baf,
      logR = log2(seqz_data$depth.ratio),
      germline_zygosity = seqz_data$zygosity.normal,
      snp_pos = paste0(chrom, "_", pos),
      row.names = paste0(chrom, "_", pos),
      stringsAsFactors = FALSE
    )

    # Mask non-hz SNPs with NA
    all_snps[[sample]]$BAF[all_snps[[sample]]$germline_zygosity != "het"] <- NA
  }

  all_snps <- ascat_set_bafsegmented(segs, all_snps)
  all_snps <- ascat_set_logrsegmented(segs, all_snps)
  all_snps <- keep_common_snps(all_snps)

  str(all_snps)

  return(all_snps)
}


# Convert the internal refphase snps object into an ascat.bc object for use in
# internal rerun of ASCAT
#
# > str(ascat.bc)
# List of 14
# $ Tumor_LogR          :'data.frame':   938946 obs. of  1 variable:
# ..$ A_LTX121_SU_T1-R1: num [1:938946] -0.0542 -0.2103 -0.2103 -0.4974 -0.3978 ...
# $ Tumor_BAF           :'data.frame':   938946 obs. of  1 variable:
# ..$ BAF: num [1:938946] 1 NA NA NA NA ...
# $ Tumor_LogR_segmented: num [1:938946, 1] -0.335 -0.335 -0.335 -0.335 -0.335 ...
# ..- attr(*, "dimnames")=List of 2
# .. ..$ : chr [1:938946] "1" "2" "3" "4" ...
# .. ..$ : chr "A_LTX121_SU_T1-R1"
# $ Tumor_BAF_segmented :List of 1
# ..$ : num [1:33666, 1] 0.23 0.23 0.23 0.23 0.23 ...
# .. ..- attr(*, "dimnames")=List of 2
# .. .. ..$ : chr [1:33666] "41" "105" "123" "135" ...
# .. .. ..$ : NULL
# $ Germline_LogR       :'data.frame':   938946 obs. of  1 variable:
# ..$ A_LTX121_SU_T1-R1: num [1:938946] 0 NA NA NA NA NA NA NA NA NA ...
# $ Germline_BAF        :'data.frame':   938946 obs. of  1 variable:
# ..$ BAF: num [1:938946] 1.01 NA NA NA NA ...
# $ SNPpos              :'data.frame':   938946 obs. of  2 variables:
# ..$ chrom: chr [1:938946] "1" "1" "1" "1" ...
# ..$ pos  : int [1:938946] 762273 762496 762497 777237 777318 777335 777361 777456 783071 783147 ...
# $ ch                  :List of 23
# ..$ : int [1:92902] 1 2 3 4 5 6 7 8 9 10 ...
# ..$ : int [1:68365] 92903 92904 92905 92906 92907 92908 92909 92910 92911 92912 ...
# ..$ : int [1:53503] 161268 161269 161270 161271 161272 161273 161274 161275 161276 161277 ...
# ..$ : int [1:22578] 899048 899049 899050 899051 899052 899053 899054 899055 899056 899057 ...
# ..$ : int [1:17321] 921626 921627 921628 921629 921630 921631 921632 921633 921634 921635 ...
# $ chr                 :List of 28
# ..$ : int [1:50594] 1 2 3 4 5 6 7 8 9 10 ...
# ..$ : int [1:42308] 50595 50596 50597 50598 50599 50600 50601 50602 50603 50604 ...
# ..$ : int [1:68365] 92903 92904 92905 92906 92907 92908 92909 92910 92911 92912 ...
# ..$ : int [1:53503] 161268 161269 161270 161271 161272 161273 161274 161275 161276 161277 ...
# ..$ : int [1:17321] 921626 921627 921628 921629 921630 921631 921632 921633 921634 921635 ...
# $ chrs                : chr [1:23] "1" "2" "3" "4" ...
# $ samples             : chr "A_LTX121_SU_T1-R1"
# $ gender              : chr "XX"
# $ sexchromosomes      : chr [1:2] "X" "Y"
# $ failedarrays        : NULL
snps_to_ascat_bc <- function(snps, sample_name) {
  bafseg <- matrix(snps$BAFsegmented, ncol = 1, dimnames = list(as.character(seq_along(snps$BAFsegmented))))
  # Remove positions with NA in BAF
  bafseg <- bafseg[!is.na(snps$BAF), , drop = FALSE]

  ascat_bc <- list(
    Tumor_LogR = data.frame(sample = snps$logR),
    Tumor_BAF = data.frame(sample = snps$BAF),
    Tumor_LogR_segmented = matrix(snps$logRsegmented,
      ncol = 1,
      dimnames = list(as.character(seq_along(snps$logRsegmented)))
    ),
    Tumor_BAF_segmented = list(bafseg),
    Germline_LogR = data.frame(),
    Germline_BAF = data.frame(),
    SNPpos = data.frame(
      chrs = as.character(snps$chrom),
      pos = snps$pos
    ),
    ch = NA,
    chr = NA,
    chrs = as.character(unique(sort(snps$chrom))),
    samples = c(sample_name),
    gender = c("XX"),
    sexchromosomes = c("23", "24"),
    failedarrays = NULL
  )

  ascat_bc$ch <- lapply(ascat_bc$chrs, function(x) {
    which(ascat_bc$SNPpos$chrs == x)
  })
  ascat_bc$chr <- split_genome(ascat_bc$SNPpos) # From ascat.R

  # colnames(ascat_bc$Germline_BAF) <- "BAF"
  colnames(ascat_bc$Tumor_BAF) <- "BAF"
  colnames(ascat_bc$Tumor_LogR) <- sample_name

  return(ascat_bc)
}
