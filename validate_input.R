# A collection of functions to ensure that our input data matches what the
# script expects.
# 
# It is broken down into functions to validate the segmentation, the input
# baf/logR data, and both of those datasets together.

# -------------------------------------------------------------------------
# SNP validation (logR, BAF and more)
# -------------------------------------------------------------------------

validate_snp_all_samples <- function(tumor_snps, germline_snps) {
  samples <- names(tumor_snps)
  # Make sure the tumor_snps list has names
  stopifnot(!is_null(samples))

  # Check individual SNP input files
  for (sample in samples) {
    message("Checking SNP input file: ", sample)
    validate_snp_input(tumor_snps[[sample]])
  }

  # Some checks across all SNP input files

  # -------------------------------------------------------------------------
  # 1) Positions are identical in all input files (both for Tumor and Germline)

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

  
}

validate_snp_input <- function(tumor_snp) {

  # -------------------------------------------------------------------------
  # 1) Are the columns in the expected order (check for types of data and
  # column names)
  # expected_cols <- list(sample = character(),
  #                       chr = integer(),
  #                       startpos = integer(),
  #                       endpos = integer(),
  #                       n.het = integer(),
  #                       cnTotal = integer(),
  #                       nMajor = integer(),
  #                       nMinor = integer(),
  #                       Ploidy = numeric(),
  #                       ACF = numeric(),
  #                       nAraw = numeric(),
  #                       nBraw = numeric(),
  #                       callState = character())

  # # Do we need all of the columns?
  # stopifnot(identical(names(seg), names(expected_cols)))

  # # Check types
  # for (expected_col in names(expected_cols))
  #   stopifnot(class(seg[,expected_col]) == class(expected_cols[[expected_col]]))

  # -------------------------------------------------------------------------
  # 2) Are any chromosome names not in the format 1-24 (check for NAs)

  # -------------------------------------------------------------------------
  # 3) Are all positions of SNPs positive integers

  # -------------------------------------------------------------------------
  # 4) Are all chromosomes and positions in the correct (1-24.. 1-200000)
  # sorted order

  # -------------------------------------------------------------------------
  # 5) Are there any NA values in the LogR columns for any regions

  # -------------------------------------------------------------------------
  # 6) Are there any Inf or -Inf values in the LogR columns

  # -------------------------------------------------------------------------
  # 7) Given the chosen reference genome are the positions for each chromosome
  # possible

  # -------------------------------------------------------------------------
  # 8) Are there any duplicate positions
  
  # -------------------------------------------------------------------------
  # 9) Are there positions present for all expected chromosomes if none present
  # for any autosome error, else note that 23/24 are missing 
  

}

# Validate raw ASCAT bc file
validate_ascat_bc_input <- function(bc) {
  # -------------------------------------------------------------------------
  # 1) Are the columns in the expected order (check for types of data and
  # column names)
  # expected_cols <- list(sample = character(),
  #                       chr = integer(),
  #                       startpos = integer(),
  #                       endpos = integer(),
  #                       n.het = integer(),
  #                       cnTotal = integer(),
  #                       nMajor = integer(),
  #                       nMinor = integer(),
  #                       Ploidy = numeric(),
  #                       ACF = numeric(),
  #                       nAraw = numeric(),
  #                       nBraw = numeric(),
  #                       callState = character())

  # # Do we need all of the columns?
  # stopifnot(identical(names(seg), names(expected_cols)))

  # # Check types
  # for (expected_col in names(expected_cols))
  #   stopifnot(class(seg[,expected_col]) == class(expected_cols[[expected_col]]))

  # -------------------------------------------------------------------------
  # 2) Are any chromosome names not in the format 1-24 (check for NAs)
  

  # -------------------------------------------------------------------------
  # 3) Are all positions of SNPs positive integers
  stopifnot(is.integer(bc$SNPpos$pos))
  stopifnot(bc$SNPpos$pos >= 0)

  # -------------------------------------------------------------------------
  # 4) Are all chromosomes and positions in the correct (1-24.. 1-200000)
  # sorted order

  # -------------------------------------------------------------------------
  # 5) Are there any NA values in the LogR columns for any regions
  stopifnot(!any(is.na(bc[["Tumor_LogR"]][,1])))
  stopifnot(!any(is.na(bc[["Tumor_LogR_segmented"]][,1])))
  stopifnot(!any(is.na(bc[["Tumor_BAF_segmented"]][[1]][,1])))

  # -------------------------------------------------------------------------
  # 6) Are there any Inf or -Inf values in the LogR columns
  #stopifnot(any(is.finite(bc[["Tumor_LogR"]][,1])))
  #stopifnot(any(is.finite(bc[["Tumor_LogR_segmented"]][,1])))
  #stopifnot(any(is.finite(bc[["Tumor_BAF_segmented"]][[1]][,1])))

  # -------------------------------------------------------------------------
  # 7) Given the chosen reference genome are the positions for each chromosome
  # possible

  # -------------------------------------------------------------------------
  # 8) Are there any duplicate positions
  stopifnot(anyDuplicated(bc$SNPpos) == 0)
  
  # -------------------------------------------------------------------------
  # 9) Are there positions present for all expected chromosomes if none present
  # for any autosome error, else note that 23/24 are missing 
}

# -------------------------------------------------------------------------
# Segmentation validation
# -------------------------------------------------------------------------

# Call this on the table with segmentations for all samples
validate_seg_all_samples <- function(segs) {
  # Run validate_seg_input() on each sample separately
  samples <- unique(segs[,"sample"])
  for (sample in samples) {
    validate_seg_input(segs[segs[,"sample"] == sample,])
  }
}

# 
validate_seg_input <- function(seg) {
  # -------------------------------------------------------------------------
  # 1) Are the columns in the expected order (check for types of data and
  # column names)
  expected_cols <- list(sample = character(),
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
                        callState = character())

  # Do we need all of the columns?
  stopifnot(identical(names(seg), names(expected_cols)))

  # Check types
  for (expected_col in names(expected_cols)) {
    if (expected_col %in% c("nAraw", "nBraw", "callState"))
      next
    message("column: ", expected_col, "\tsegmentation: ", class(seg[,expected_col]), "\texpected: ", class(expected_cols[[expected_col]]))
    stopifnot(class(seg[,expected_col]) == class(expected_cols[[expected_col]]))
  }

  # -------------------------------------------------------------------------
  # 2) Are any chromosome names not in the format 1-24 (check for NAs)

  # No NA's
  stopifnot(all(!is.na(seg[,"chr"])))

  # Must be in 1-24
  stopifnot(seg[,"chr"] %in% 1:24)

  # Must not be missing a chromosome in 1-22
  stopifnot(all(1:22 %in% seg[,"chr"]))

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
  #stopifnot(!is.unsorted(order(seg[,c("chr", "startpos", "endpos")])))

  # -------------------------------------------------------------------------
  # 6) Do any segments within the same chromosome overlap
  # FIXME: not implemented

  # 7) Do all chromosomes have at least one segment present, if none for an
  # autosome error, if none for a sex chromosome flag
  # NOTE: already checked autosomes in step 2)
  if (!any(c(23, 24) %in% seg[,"chr"]))
    warning("No segment in either sex chromosome (23 or 24)")

  # 8) Given the chosen reference genome are the positions for each chromosome
  # possible
  # FIXME: not implemented

  message("Segment input validation passed.")
  return(TRUE)

}
