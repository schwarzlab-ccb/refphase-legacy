# A collection of functions to ensure that our input data matches what the
# script expects.
# 
# It is broken down into functions to validate the segmentation, the input
# baf/logR data, and both of those datasets together.

# Wrapper for all validation functions 
validate_snp_input <- function(snp_data) {
  
}

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
  for (expected_col in names(expected_cols))
    stopifnot(class(seg[,expected_col]) == class(expected_cols[[expected_col]]))

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
  stopifnot(with(seg, any(startpos == endpos)))

  # -------------------------------------------------------------------------
  # 5) Are all segments (by chromosomes and positions) in the correct (1-24..
  # 1-200000) sorted order
  stopifnot(!is.unsorted(order(seg[,c("chr", "startpos", "endpos")])))

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
