# A collection of functions to ensure that our input data matches what the
# script expects.
#
# It is broken down into functions to validate the segmentation, the input
# baf/logR data, and both of those datasets together.

# -------------------------------------------------------------------------
# SNP validation (logR, BAF and more)
# -------------------------------------------------------------------------

validate_snp_all_samples <- function(tumor_snps) {
  logdebug("Validating SNP data for all samples.")

  samples <- names(tumor_snps)
  # Make sure the tumor_snps list has names
  stopifnot(!is.null(samples))

  # We need at least two samples
  stopifnot(length(tumor_snps) >= 2)

  # Check individual SNP input files
  for (sample in samples) {
    logdebug(paste0("Checking SNP input file: ", sample))
    validate_snp_input(tumor_snps[[sample]])
  }

  # Some checks across all SNP input files

  # -------------------------------------------------------------------------
  # 1) Positions are identical in all input files (both for Tumor and Germline)
  stopifnot(all(sapply(tumor_snps[-1], FUN = function(x, y) {
    identical(rownames(x), rownames(y))
  }, tumor_snps[[1]])))

  # Checks 2 and 3 are not necessary if the above check passes.

  # -------------------------------------------------------------------------
  # 2) All the same chromosomes are represented in every SNP file (covered by
  # the above but do anyway)

  # -------------------------------------------------------------------------
  # 3) The total number of lines/rows in each dataframe is the same for each
  # SNP file - check on a chromosome by chromosome basis as well (overkill but
  # won’t hurt)

  # -------------------------------------------------------------------------
  # 4) All SNPs that are “BAF” SNPs ie between 0.3 and 0.7 in the germline have
  # a BAF value in every SNP file.

  # From the original input code. Can probably be skipped.
  # Checks to see if all samples have the same number of snps (NOT that those
  # snps are identical), and that there are at least 2 samples
  if (length(unique(unlist(lapply(tumor_snps, nrow)))) != 1 | length(tumor_snps) < 2) {
    stop("Not all sites are present in all samples' sites information read in or too few samples are present at this stage, this needs to be fixed before running MPhase.")
  }

  logdebug("Validation passed.")
}

validate_snp_input <- function(snp_data) {

  # -------------------------------------------------------------------------
  # 1) Are the columns in the expected order (check for types of data and
  # column names)
  expected_cols <- list(
    chrom = integer(),
    pos = integer(),
    BAF = numeric(),
    logR = numeric(),
    BAFsegmented = numeric(),
    logRsegmented = numeric()
  )
  # snp_pos = character())

  # Make sure we have (at least) the columns that we need
  stopifnot(all(names(expected_cols) %in% names(snp_data)))

  # Check types
  for (expected_col in names(expected_cols)) {
    stopifnot(class(snp_data[, expected_col]) == class(expected_cols[[expected_col]]))
  }

  # -------------------------------------------------------------------------
  # 2) Are any chromosome names not in the format 1-24 (check for NAs)
  stopifnot(all(snp_data$chrom %in% 1:24))

  # -------------------------------------------------------------------------
  # 3) Are all positions of SNPs positive integers
  # stopifnot(is.integer(snp_data$pos))  # Already checked above
  stopifnot(all(snp_data$pos >= 0))

  # -------------------------------------------------------------------------
  # 4) Are all chromosomes and positions in the correct (1-24.. 1-200000)
  # sorted order

  # -------------------------------------------------------------------------
  # 5) Are there any NA values in the LogR columns for any regions
  stopifnot(all(!is.na(snp_data$logR)))

  # -------------------------------------------------------------------------
  # 6) Are there any Inf or -Inf values in the LogR columns
  stopifnot(all(is.finite(snp_data$logR)))

  # -------------------------------------------------------------------------
  # 7) Given the chosen reference genome are the positions for each chromosome
  # possible

  # -------------------------------------------------------------------------
  # 8) Are there any duplicate positions
  stopifnot(anyDuplicated(snp_data[, c("chrom", "pos")]) == 0)

  # -------------------------------------------------------------------------
  # 9) Are there positions present for all expected chromosomes if none present
  # for any autosome error
  stopifnot(all(1:22 %in% snp_data$chrom))

  # else note that 23/24 are missing
  if (!all(c(23, 24) %in% snp_data$chrom)) {
    logwarn("No segment in either sex chromosome (23 or 24)")
  }
}

# Validate raw ASCAT bc file
validate_ascat_bc_input <- function(bc) {
  # -------------------------------------------------------------------------
  # 1) Are the columns in the expected order (check for types of data and
  # column names)
  #   > x <- readRDS("ascat_bc.rds")
  # > lapply(x, class)
  # $Tumor_LogR [1] "data.frame"
  # $Tumor_BAF [1] "data.frame"
  # $Tumor_LogR_segmented [1] "matrix"
  # $Tumor_BAF_segmented [1] "list"
  # $Germline_LogR [1] "data.frame"
  # $Germline_BAF [1] "data.frame"
  # $SNPpos [1] "data.frame"
  # $ch [1] "list"
  # $chr [1] "list"
  # $chrs [1] "character"
  # $samples [1] "character"
  # $gender [1] "character"
  # $sexchromosomes [1] "character"
  # $failedarrays [1] "NULL"
  expected_cols <- list(
    Tumor_LogR = data.frame(),
    Tumor_BAF = data.frame(),
    Tumor_LogR_segmented = matrix(),
    Tumor_BAF_segmented = list(),
    Germline_LogR = data.frame(),
    Germline_BAF = data.frame(),
    SNPpos = data.frame(),
    ch = list(),
    chr = list(),
    chrs = character(),
    samples = character(),
    gender = character(),
    sexchromosomes = character()
  )

  # We don't care about the type of this one too much, it just needs to be present
  # failedarrays = character())

  # Do we need all of the columns?
  stopifnot(all(names(expected_cols) %in% names(bc)))
  stopifnot("failedarrays" %in% names(bc))

  # # Check types
  for (expected_col in names(expected_cols)) {
    if (class(bc[[expected_col]]) != class(expected_cols[[expected_col]])) {
      logwarn(paste0(expected_col, " doesn't have the correct class. Expected: ", class(expected_cols[[expected_col]]), " Found: ", class(bc[[expected_col]])))
      stopifnot(FALSE)
    }
  }

  # -------------------------------------------------------------------------
  # 2) Are any chromosome names not in the format 1-24 (check for NAs)
  stopifnot(all(bc$chrs %in% as.character(1:24)))
  stopifnot(all(bc$SNPpos$chrs %in% as.character(1:24)))


  # -------------------------------------------------------------------------
  # 3) Are all positions of SNPs positive integers
  stopifnot(is.integer(bc$SNPpos$pos))
  stopifnot(bc$SNPpos$pos >= 0)

  # -------------------------------------------------------------------------
  # 4) Are all chromosomes and positions in the correct (1-24.. 1-200000)
  # sorted order

  # -------------------------------------------------------------------------
  # 5) Are there any NA values in the LogR columns for any regions
  logwarn(all(!is.na(bc[["Tumor_LogR"]][, 1])))
  logwarn(all(!is.na(bc[["Tumor_LogR_segmented"]][, 1])))
  logwarn(all(!is.na(bc[["Tumor_BAF_segmented"]][[1]][, 1])))

  # -------------------------------------------------------------------------
  # 6) Are there any Inf or -Inf values in the LogR columns
  logwarn(all(is.finite(bc[["Tumor_LogR"]][, 1])))
  logwarn(all(is.finite(bc[["Tumor_LogR_segmented"]][, 1])))
  logwarn(all(is.finite(bc[["Tumor_BAF_segmented"]][[1]][, 1])))

  # -------------------------------------------------------------------------
  # 7) Given the chosen reference genome are the positions for each chromosome
  # possible

  # -------------------------------------------------------------------------
  # 8) Are there any duplicate positions
  stopifnot(anyDuplicated(bc$SNPpos) == 0)

  # -------------------------------------------------------------------------
  # 9) Are there positions present for all expected chromosomes if none present
  # for any autosome error, else note that 23/24 are missing
  stopifnot(as.character(1:22) %in% unique(bc$SNPpos$chrs))
}

# -------------------------------------------------------------------------
# Segmentation validation
# -------------------------------------------------------------------------

# Call this on the table with segmentations for all samples
validate_seg_all_samples <- function(segs) {
  # Run validate_seg_input() on each sample separately
  samples <- unique(segs[, "sample"])
  for (sample in samples) {
    validate_seg_input(segs[segs[, "sample"] == sample, ])
  }
}

#
validate_seg_input <- function(seg) {
  # -------------------------------------------------------------------------
  # 1) Are the columns in the expected order (check for types of data and
  # column names)
  expected_cols <- list(
    sample = character(),
    chr = integer(),
    startpos = integer(),
    endpos = integer(),
    n.het = integer(),
    cnTotal = integer(),
    nMajor = integer(),
    nMinor = integer(),
    Ploidy = numeric(),
    ACF = numeric(),
    nAraw = numeric(),
    nBraw = numeric(),
    callState = character()
  )

  # Do we need all of the columns?
  stopifnot(identical(names(seg), names(expected_cols)))

  # Check types
  for (expected_col in names(expected_cols)) {
    if (expected_col %in% c("nAraw", "nBraw", "callState")) {
      next
    }
    # message("column: ", expected_col, "\tsegmentation: ", class(seg[,expected_col]), "\texpected: ", class(expected_cols[[expected_col]]))
    stopifnot(class(seg[, expected_col]) == class(expected_cols[[expected_col]]))
  }

  # -------------------------------------------------------------------------
  # 2) Are any chromosome names not in the format 1-24 (check for NAs)

  # No NA's
  stopifnot(all(!is.na(seg[, "chr"])))

  # Must be in 1-24
  stopifnot(seg[, "chr"] %in% 1:24)

  # Must not be missing a chromosome in 1-22
  stopifnot(all(1:22 %in% seg[, "chr"]))

  # -------------------------------------------------------------------------
  # 3) Within a chromosome are all segment start positions <= end positions
  stopifnot(with(seg, all(startpos <= endpos)))

  # -------------------------------------------------------------------------
  # 4) Are there any single SNP segments ie start position == end position
  stopifnot(with(seg, all(startpos != endpos)))

  # -------------------------------------------------------------------------
  # 5) Are all segments (by chromosomes and positions) in the correct (1-24..
  # 1-200000) sorted order
  # FIXME
  # stopifnot(!is.unsorted(order(seg[,c("chr", "startpos", "endpos")])))

  # -------------------------------------------------------------------------
  # 6) Do any segments within the same chromosome overlap
  # FIXME: not implemented

  # 7) Do all chromosomes have at least one segment present, if none for an
  # autosome error, if none for a sex chromosome flag
  # NOTE: already checked autosomes in step 2)
  if (!any(c(23, 24) %in% seg[, "chr"])) {
    logwarn("No segment in either sex chromosome (23 or 24)")
  }

  # 8) Given the chosen reference genome are the positions for each chromosome
  # possible
  # FIXME: not implemented

  logdebug("Segment input validation passed.")
}
