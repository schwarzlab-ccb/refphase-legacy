joinSegAndReturnNewArmGranges <- function(arm.idx.vec, curr.gr.idx, adj.idx, adj.segs.mat, sorted.arm.grs) {
  kept.arm.idxs <- arm.idx.vec[ !arm.idx.vec %in% c(curr.gr.idx, adj.segs.mat[adj.idx, "gr_idx"]) ]

  if (length(kept.arm.idxs) > 0) {
    kept.arm.chr.vec <- as.numeric(as.character(seqnames(sorted.arm.grs[kept.arm.idxs])))
    kept.arm.start.vec <- start(sorted.arm.grs[kept.arm.idxs])
    kept.arm.end.vec <- end(sorted.arm.grs[kept.arm.idxs])
  } else {
    kept.arm.chr.vec <- c()
    kept.arm.start.vec <- c()
    kept.arm.end.vec <- c()
  }

  new.seg.chr <- adj.segs.mat[adj.idx, "chrom"]
  new.seg.start <- min(adj.segs.mat[adj.idx, "startpos"], start(sorted.arm.grs[curr.gr.idx]))
  new.seg.end <- max(adj.segs.mat[adj.idx, "endpos"], end(sorted.arm.grs[curr.gr.idx]))

  new.arm.chr.vec <- c(new.seg.chr, kept.arm.chr.vec)
  new.arm.start.vec <- c(new.seg.start, kept.arm.start.vec)
  new.arm.end.vec <- c(new.seg.end, kept.arm.end.vec)

  new.arm.grs <- GRanges(
    seqnames = new.arm.chr.vec,
    ranges = IRanges(new.arm.start.vec, new.arm.end.vec),
    strand = rep("*", length(new.arm.chr.vec))
  )

  return(new.arm.grs)
}


# Try to join using phasing
attemptToJoinSegbyPhasedCPN <- function(curr.seg.mat, adj.segs.mat, samp.arm.sites.df, arm.idx.vec, curr.gr.idx, sorted.arm.grs, cpnX.string, MIN_PH_CPN) {
  new.unsorted.arm.grs <- sorted.arm.grs

  for (k in 1:nrow(adj.segs.mat)) {
    if (curr.seg.mat[, "num_ph_cpn"] >= MIN_PH_CPN & adj.segs.mat[k, "num_ph_cpn"] >= MIN_PH_CPN) {
      cur.ph.cpn.vec <- samp.arm.sites.df[samp.arm.sites.df$pos >= curr.seg.mat[, "startpos"]
      & samp.arm.sites.df$pos <= curr.seg.mat[, "endpos"], cpnX.string]

      adj.ph.cpn.vec <- samp.arm.sites.df[samp.arm.sites.df$pos >= adj.segs.mat[k, "startpos"]
      & samp.arm.sites.df$pos <= adj.segs.mat[k, "endpos"], cpnX.string]

      ph.p.val <- wilcoxTestWithErrorHandling(cur.ph.cpn.vec, adj.ph.cpn.vec)$p.val

      if (ph.p.val > PVAL_THRESH_TO_JOIN_PH_SEGS) {
        # Join segment and make a new set of granges for the arm
        new.unsorted.arm.grs <- joinSegAndReturnNewArmGranges(
          arm.idx.vec = arm.idx.vec,
          curr.gr.idx = curr.gr.idx,
          adj.idx = k,
          adj.segs.mat = adj.segs.mat,
          sorted.arm.grs = sorted.arm.grs
        )
        return(new.unsorted.arm.grs)
      }
    }
  }

  return(new.unsorted.arm.grs)
}


# If join return
attemptToJoinSegbyTotalCPN <- function(curr.seg.mat, adj.segs.mat, samp.arm.sites.df, arm.idx.vec, curr.gr.idx, sorted.arm.grs, NUM_TOT_CPN) {
  new.unsorted.arm.grs <- sorted.arm.grs

  for (k in 1:nrow(adj.segs.mat)) {

    # If joining by phased copy number is possible don't join
    if (curr.seg.mat[, "num_tot_cpn"] >= NUM_TOT_CPN & adj.segs.mat[k, "num_tot_cpn"] >= NUM_TOT_CPN) {
      join.by.phased.cpn.not.possible <- curr.seg.mat[, "num_ph_cpn"] < MIN_PH_CPN | adj.segs.mat[k, "num_ph_cpn"] < MIN_PH_CPN

      if (isTRUE(join.by.phased.cpn.not.possible)) {
        cur.tot.cpn.vec <- rle(samp.arm.sites.df[samp.arm.sites.df$pos >= curr.seg.mat[, "startpos"]
        & samp.arm.sites.df$pos <= curr.seg.mat[, "endpos"], "tot_raw_un_ph_cpn"])$values

        adj.tot.cpn.vec <- rle(samp.arm.sites.df[samp.arm.sites.df$pos >= adj.segs.mat[k, "startpos"]
        & samp.arm.sites.df$pos <= adj.segs.mat[k, "endpos"], "tot_raw_un_ph_cpn"])$values

        tot.p.val <- wilcoxTestWithErrorHandling(cur.tot.cpn.vec, adj.tot.cpn.vec)$p.val

        if (tot.p.val > PVAL_THRESH_TO_JOIN_TOT_SEGS) {
          new.unsorted.arm.grs <- joinSegAndReturnNewArmGranges(
            arm.idx.vec = arm.idx.vec,
            curr.gr.idx = curr.gr.idx,
            adj.idx = k,
            adj.segs.mat = adj.segs.mat,
            sorted.arm.grs = sorted.arm.grs
          )

          return(new.unsorted.arm.grs)
        }
      }
    }
  }

  return(new.unsorted.arm.grs)
}


attemptToJoinSegbyMeanTotalCPNAlone <- function(curr.seg.mat, adj.segs.mat, samp.arm.sites.df, arm.idx.vec, curr.gr.idx, sorted.arm.grs, NUM_TOT_CPN, CPN_AVG_THRESH) {
  new.unsorted.arm.grs <- sorted.arm.grs

  # Allow joining by means...
  for (k in 1:nrow(adj.segs.mat)) {
    if (curr.seg.mat[, "num_tot_cpn"] > 0 & curr.seg.mat[, "num_tot_cpn"] < NUM_TOT_CPN & adj.segs.mat[k, "num_tot_cpn"] > 0) {
      curr.seg.mean <- curr.seg.mat[, "mean_tot_cpn"]
      curr.seg.median <- curr.seg.mat[, "median_tot_cpn"]

      adj.seg.mean <- adj.segs.mat[k, "mean_tot_cpn"]
      adj.seg.median <- adj.segs.mat[k, "median_tot_cpn"]

      if (abs(curr.seg.mean - adj.seg.mean) <= CPN_AVG_THRESH | abs(curr.seg.median - adj.seg.median) <= CPN_AVG_THRESH) {
        new.unsorted.arm.grs <- joinSegAndReturnNewArmGranges(
          arm.idx.vec = arm.idx.vec,
          curr.gr.idx = curr.gr.idx,
          adj.idx = k,
          adj.segs.mat = adj.segs.mat,
          sorted.arm.grs = sorted.arm.grs
        )

        return(new.unsorted.arm.grs)
      }
    }
  }

  return(new.unsorted.arm.grs)
}


joinBlankSegsToMostConfidentCall <- function(curr.seg.mat, adj.segs.mat, samp.arm.sites.df, arm.idx.vec, curr.gr.idx, sorted.arm.grs) {
  new.unsorted.arm.grs <- sorted.arm.grs
  adj.segs.mat <- adj.segs.mat[order(adj.segs.mat[, "num_ph_cpn"], adj.segs.mat[, "num_tot_cpn"]), , drop = FALSE]
  adj.segs.mat <- adj.segs.mat[1, , drop = FALSE]

  # Allow joining by means...
  for (k in 1:nrow(adj.segs.mat)) {
    if (curr.seg.mat[, "num_tot_cpn"] == 0) {
      new.unsorted.arm.grs <- joinSegAndReturnNewArmGranges(
        arm.idx.vec = arm.idx.vec,
        curr.gr.idx = curr.gr.idx,
        adj.idx = k,
        adj.segs.mat = adj.segs.mat,
        sorted.arm.grs = sorted.arm.grs
      )

      return(new.unsorted.arm.grs)
    }
  }

  return(new.unsorted.arm.grs)
}


# Get all segment names  sorted by whether they're on a p por q arm...

wilcoxTestWithErrorHandling <- function(x, y = NULL,
                                        alternative = c("two.sided", "less", "greater"),
                                        mu = 0, paired = FALSE, exact = FALSE, correct = TRUE,
                                        conf.int = FALSE, conf.level = 0.95, errors.not.permitted = "") {
  test.results <- tryCatch(
    {
      test.results <- wilcox.test(
        x = x,
        y = y,
        alternative = alternative,
        mu = mu,
        paired = paired,
        exact = exact,
        correct = correct,
        conf.int = conf.int,
        conf.level = conf.level
      )

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
      test.results[["p.val"]] <- 1
      return(test.results)
    }
  )

  if (is.na(test.results$p.val)) {
    test.results$p.val <- 1
  }

  return(test.results)
}


creatingArmJoiningInfoMatrixAndGranges <- function(arm.grs, samp.arm.sites.df, cpnX.string) {
  proc.arm.grs <- arm.grs
  proc.arm.grs$seg_width <- GenomicRanges::width(proc.arm.grs)
  joining.vals.mat.col.names.vec <- c(
    "chrom", "startpos", "endpos", "gr_idx",
    "num_ph_cpn", "mean_ph_cpn", "median_ph_cpn",
    "num_tot_cpn", "mean_tot_cpn", "median_tot_cpn"
  )
  joining.vals.mat <- matrix(ncol = 10, nrow = length(proc.arm.grs), 0)
  colnames(joining.vals.mat) <- joining.vals.mat.col.names.vec

  for (idx in 1:length(proc.arm.grs)) {
    # Get the idx of the grange in the sorted by size arm grs list
    joining.vals.mat[idx, "chrom"] <- as.numeric(as.character(GenomicRanges::seqnames(proc.arm.grs[idx])))
    joining.vals.mat[idx, "startpos"] <- GenomicRanges::start(proc.arm.grs[idx])
    joining.vals.mat[idx, "endpos"] <- GenomicRanges::end(proc.arm.grs[idx])
    seg.sites.info.df <- samp.arm.sites.df[samp.arm.sites.df$chrom == as.numeric(as.character(GenomicRanges::seqnames(proc.arm.grs[idx])))
    & samp.arm.sites.df$pos >= GenomicRanges::start(proc.arm.grs[idx])
    & samp.arm.sites.df$pos <= GenomicRanges::end(proc.arm.grs[idx]), ]

    seg.ph.cpn.sites.df <- seg.sites.info.df[!is.na(seg.sites.info.df[, cpnX.string]), ]
    num.ph.cpn <- length(seg.ph.cpn.sites.df[, cpnX.string])

    if (num.ph.cpn > 0) {
      joining.vals.mat[idx, "num_ph_cpn"] <- num.ph.cpn
      joining.vals.mat[idx, "mean_ph_cpn"] <- mean(seg.ph.cpn.sites.df[, cpnX.string])
      joining.vals.mat[idx, "median_ph_cpn"] <- median(seg.ph.cpn.sites.df[, cpnX.string])
    }

    # using run length encodings here as we have the issue that many logr
    # values, due to being 1000G positiosn are often very close together and
    # thus are, in reality, a single datapoint...
    # Could be a problem with very low coverage samples...?todo, check this
    num.tot.cpn <- length(rle(seg.sites.info.df$tot_raw_un_ph_cpn)$values)

    if (num.tot.cpn > 0) {
      joining.vals.mat[idx, "num_tot_cpn"] <- num.tot.cpn
      joining.vals.mat[idx, "mean_tot_cpn"] <- mean(base::rle(seg.sites.info.df$tot_raw_un_ph_cpn)$values)
      joining.vals.mat[idx, "median_tot_cpn"] <- median(base::rle(seg.sites.info.df$tot_raw_un_ph_cpn)$values)
    }
  }

  proc.arm.grs$num_ph_cpn <- joining.vals.mat[, "num_ph_cpn"]
  proc.arm.grs$num_tot_cpn <- joining.vals.mat[, "num_tot_cpn"]
  sorted.arm.grs <- proc.arm.grs[ order(proc.arm.grs$num_ph_cpn, proc.arm.grs$num_tot_cpn, proc.arm.grs$seg_width) ]
  joining.vals.mat <- joining.vals.mat[match(GenomicRanges::start(sorted.arm.grs), joining.vals.mat[, "startpos"]), , drop = FALSE]
  joining.vals.mat[, "gr_idx"] <- 1:nrow(joining.vals.mat)

  return(list(joining.vals.mat, sorted.arm.grs))
}


attemptToJoinOneSegmentWithinChromArm <- function(current.arm.grs, joining.vals.mat, samp.arm.sites.df, sorted.arm.grs, cpnX.string, MIN_PH_CPN, NUM_TOT_CPN, CPN_AVG_THRESH) {

  # Edge case where this only one segment:
  if (length(current.arm.grs) == 1) {
    return(current.arm.grs)
  }

  arm.idx.vec <- 1:length(current.arm.grs)
  old.length.current.arm.grs <- length(current.arm.grs)
  curr.ord.idx <- 1
  while (curr.ord.idx <= length(current.arm.grs) & length(current.arm.grs) == old.length.current.arm.grs) {
    curr.gr <- current.arm.grs[curr.ord.idx]
    precede.and.follow.gr.idx.vec <- c(follow(curr.gr, current.arm.grs), precede(curr.gr, current.arm.grs))
    sort.non.na.gr.pre.and.fol.idx.vec <- sort(precede.and.follow.gr.idx.vec[ !is.na(precede.and.follow.gr.idx.vec)])
    curr.seg.mat <- joining.vals.mat[joining.vals.mat[, "gr_idx"] == curr.ord.idx, , drop = FALSE]

    adj.segs.mat <- joining.vals.mat[sort.non.na.gr.pre.and.fol.idx.vec, , drop = FALSE]

    new.unsorted.arm.grs <- attemptToJoinSegbyPhasedCPN(
      curr.seg.mat = curr.seg.mat,
      adj.segs.mat = adj.segs.mat,
      samp.arm.sites.df = samp.arm.sites.df,
      arm.idx.vec = arm.idx.vec,
      curr.gr.idx = curr.ord.idx,
      sorted.arm.grs = current.arm.grs,
      cpnX.string = cpnX.string,
      MIN_PH_CPN = MIN_PH_CPN
    )

    if (length(new.unsorted.arm.grs) != length(current.arm.grs)) {
      current.arm.grs <- new.unsorted.arm.grs
      next
    }
    new.unsorted.arm.grs <- attemptToJoinSegbyTotalCPN(
      curr.seg.mat = curr.seg.mat,
      adj.segs.mat = adj.segs.mat,
      samp.arm.sites.df = samp.arm.sites.df,
      arm.idx.vec = arm.idx.vec,
      curr.gr.idx = curr.ord.idx,
      sorted.arm.grs = current.arm.grs,
      NUM_TOT_CPN = NUM_TOT_CPN
    )
    if (length(new.unsorted.arm.grs) != length(current.arm.grs)) {
      current.arm.grs <- new.unsorted.arm.grs
      next
    }

    new.unsorted.arm.grs <- attemptToJoinSegbyMeanTotalCPNAlone(
      curr.seg.mat = curr.seg.mat,
      adj.segs.mat = adj.segs.mat,
      samp.arm.sites.df = samp.arm.sites.df,
      arm.idx.vec = arm.idx.vec,
      curr.gr.idx = curr.ord.idx,
      sorted.arm.grs = current.arm.grs,
      NUM_TOT_CPN = NUM_TOT_CPN,
      CPN_AVG_THRESH = CPN_AVG_THRESH
    )

    if (length(new.unsorted.arm.grs) != length(current.arm.grs)) {
      current.arm.grs <- new.unsorted.arm.grs
      next
    }

    new.unsorted.arm.grs <- joinBlankSegsToMostConfidentCall(
      curr.seg.mat = curr.seg.mat,
      adj.segs.mat = adj.segs.mat,
      samp.arm.sites.df = samp.arm.sites.df,
      arm.idx.vec = arm.idx.vec,
      curr.gr.idx = curr.ord.idx,
      sorted.arm.grs = current.arm.grs
    )

    if (length(new.unsorted.arm.grs) != length(current.arm.grs)) {
      current.arm.grs <- new.unsorted.arm.grs
      next
    }

    curr.ord.idx <- curr.ord.idx + 1
  }

  return(current.arm.grs)
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
  # Input check
  if (!((class(query) == "GRanges" & class(subject) == "GRanges") | (class(query) == "IRanges" & class(subject) == "IRanges"))) {
    stop("Query and subject need to be of same class, either GRanges or IRanges!")
  }

  # Find overlapping ranges
  if (class(query) == "GRanges") {
    seqlengths(query) <- rep(NA, length(seqlengths(query)))
    seqlengths(subject) <- rep(NA, length(seqlengths(subject)))
  }
  overlap.index.mat <- as.matrix(findOverlaps(query, subject, ...))

  query <- query[overlap.index.mat[, 1]]
  subject <- subject[overlap.index.mat[, 2]]
  overlap.mat <- cbind(query_start = start(query), query_end = end(query), subject_start = start(subject), subject_end = end(subject))

  # Pre-queries for overlaps
  startup <- overlap.mat[, "subject_start"] < overlap.mat[, "query_start"]
  enddown <- overlap.mat[, "subject_end"] > overlap.mat[, "query_end"]
  startin <- overlap.mat[, "subject_start"] >= overlap.mat[, "query_start"] & overlap.mat[, "subject_start"] <= overlap.mat[, "query_end"]
  endin <- overlap.mat[, "subject_end"] >= overlap.mat[, "query_start"] & overlap.mat[, "subject_end"] <= overlap.mat[, "query_end"]

  # Overlap types
  olup <- startup & endin
  oldown <- startin & enddown
  inside <- startin & endin
  contained <- startup & enddown

  # Overlap types in one vector
  overlap_type <- rep("", length(overlap.mat[, "query_start"]))
  overlap_type[olup] <- "upstream_overlap"
  overlap_type[oldown] <- "downstream_overlap"
  overlap_type[inside] <- "subject_inside_query"
  overlap_type[contained] <- "query_contained_in_subject"

  # Overlap positions
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

  # Absolute and relative length of overlaps
  overlap_length <- (overlap_end - overlap_start) + 1
  overlap_prop_query <- overlap_length / width(query) * 100
  overlap_prop_subject <- overlap_length / width(subject) * 100

  # Output type
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
    if (class(query) == "IRanges") {
      # Todo add the metadata in th eiranges
      query <- GRanges(seqnames = Rle(rep("dummy", length(query))), ranges = IRanges(start = oldf[, "query_start"], end = oldf[, "query_end"]), strand = Rle(strand(rep("+", length(query)))), oldf)
    }
    return(query)
  }
}
