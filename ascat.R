# ASCAT version 2.3, 15/07/2014
# author: Peter Van Loo
# PCF and ASPCF: Gro Nilsen
# GC correction: Jiqiu Cheng
# Updated  2015-02-11: Lowered threshold for setting ACF to 1, from 0.2 to 0.15, added functionality to read files from R objects both to ascat.loadData and GCcorrection


# plots SNP array data
# input: an ASCAT object (e.g. from ASCAT.ASPCF) and plots the SNP array data before and after segmentation
ascat.plotSegmentedData <- function(ASCATobj) {
  for (arraynr in seq_len(ncol(ASCATobj$Tumor_LogR))) {
    Select_nonNAs <- rownames(ASCATobj$Tumor_BAF_segmented[[arraynr]])
    AllIDs <- 1:dim(ASCATobj$Tumor_LogR)[1]
    names(AllIDs) <- rownames(ASCATobj$Tumor_LogR)
    HetIDs <- AllIDs[Select_nonNAs]
    png(filename = paste(ASCATobj$samples[arraynr], ".ASPCF.png", sep = ""), width = 2000, height = 1000, res = 200)
    par(mar = c(0.5, 5, 5, 0.5), mfrow = c(2, 1), cex = 0.4, cex.main = 3, cex.axis = 2)
    r <- ASCATobj$Tumor_LogR_segmented[rownames(ASCATobj$Tumor_BAF_segmented[[arraynr]]), arraynr]
    beta <- ASCATobj$Tumor_BAF_segmented[[arraynr]][, ]
    plot(c(1, length(r)), c(-1, 1), type = "n", xaxt = "n", main = paste(colnames(ASCATobj$Tumor_BAF)[arraynr], ", LogR", sep = ""), xlab = "", ylab = "")
    points(ASCATobj$Tumor_LogR[rownames(ASCATobj$Tumor_BAF_segmented[[arraynr]]), arraynr], col = "red", pch = ifelse(dim(ASCATobj$Tumor_LogR)[1] > 100000, ".", 20))
    points(r, col = "green")
    abline(v = 0.5, lty = 1, col = "lightgrey")
    chrk_tot_len <- 0
    for (j in seq_len(length(ASCATobj$ch))) {
      chrk <- intersect(ASCATobj$ch[[j]], HetIDs)
      chrk_tot_len_prev <- chrk_tot_len
      chrk_tot_len <- chrk_tot_len + length(chrk)
      vpos <- chrk_tot_len
      tpos <- (chrk_tot_len + chrk_tot_len_prev) / 2
      text(tpos, 1, ASCATobj$chrs[j], pos = 1, cex = 2)
      abline(v = vpos + 0.5, lty = 1, col = "lightgrey")
    }
    plot(c(1, length(beta)), c(0, 1), type = "n", xaxt = "n", main = paste(colnames(ASCATobj$Tumor_BAF)[arraynr], ", BAF", sep = ""), xlab = "", ylab = "")
    points(ASCATobj$Tumor_BAF[rownames(ASCATobj$Tumor_BAF_segmented[[arraynr]]), arraynr], col = "red", pch = ifelse(dim(ASCATobj$Tumor_LogR)[1] > 100000, ".", 20))
    points(beta, col = "green")
    points(1 - beta, col = "green")
    abline(v = 0.5, lty = 1, col = "lightgrey")
    chrk_tot_len <- 0
    for (j in seq_len(length(ASCATobj$ch))) {
      chrk <- intersect(ASCATobj$ch[[j]], HetIDs)
      chrk_tot_len_prev <- chrk_tot_len
      chrk_tot_len <- chrk_tot_len + length(chrk)
      vpos <- chrk_tot_len
      tpos <- (chrk_tot_len + chrk_tot_len_prev) / 2
      text(tpos, 1, ASCATobj$chrs[j], pos = 1, cex = 2)
      abline(v = vpos + 0.5, lty = 1, col = "lightgrey")
    }
    dev.off()
  }
}


# the ASCAT main function, calculating the allele-specific copy numbers
# input: (i) an ASCAT object from ascat.aspcf.
# (ii) gamma: technology parameter, compaction of Log R profiles (expected decrease in case of deletion in diploid sample, 100 % aberrant cells; 1 in ideal case, 0.55 of Illumina 109K arrays)
# (iii) rho_manual: optional argument to override ASCAT optimization and supply rho and psi parameters (not recommended) (purity, use same format (e.g. 0.45), multiple samples, input vector)
# (iv) psi_manual: optional argument to override ASCAT optimization and supply rho and psi parameters (not recommended) (ploidy, use same format (e.g. 2.45), multiple samples, input vector)
# output: an ASCAT output object, containing:
# 1. nA: copy number of the A allele
# 2. nB: copy number of the B allele
# 3. aberrantcellfraction: the aberrant cell fraction of all arrays
# 4. ploidy: the ploidy of all arrays
# 5. psi: ploidy parameter, an estimate of the tumour's ploidy for all arrays
# 6. goodnessOfFit: the goodness of fit for each copy number solution
# 7. failedarrays: arrays on which ASCAT analysis failed
# 8. nonaberrantarrays = arrays on which ASCAT analysis indicates that they so virtually no aberrations
# 9. segments: an array containing the copy number segments of each sample (not including failed arrays)
# 10. segments_raw: an array containing the copy number segments of each sample without any rounding applied
# note: for copy number only probes, nA contains the copy number value and nB = 0.
ascat.runAscat <- function(ASCATobj, gamma = 0.55, rho_manual = NA, psi_manual = NA) {
  goodarrays <- NULL
  res <- vector("list", ncol(ASCATobj$Tumor_LogR))
  for (arraynr in seq_len(ncol(ASCATobj$Tumor_LogR))) {
    logdebug(paste0("Sample ", ASCATobj$samples[arraynr], " (", arraynr, "/", length(ASCATobj$samples), ")"))
    lrr <- ASCATobj$Tumor_LogR[, arraynr]
    names(lrr) <- rownames(ASCATobj$Tumor_LogR)
    baf <- ASCATobj$Tumor_BAF[, arraynr]
    names(baf) <- rownames(ASCATobj$Tumor_BAF)
    lrrsegm <- ASCATobj$Tumor_LogR_segmented[, arraynr]
    names(lrrsegm) <- rownames(ASCATobj$Tumor_LogR_segmented)
    bafsegm <- ASCATobj$Tumor_BAF_segmented[[arraynr]][, ]
    names(bafsegm) <- rownames(ASCATobj$Tumor_BAF_segmented[[arraynr]])
    failedqualitycheck <- F
    if (ASCATobj$samples[arraynr] %in% ASCATobj$failedarrays) {
      failedqualitycheck <- T
    }
    if (all(is.na(rho_manual))) {
      res[[arraynr]] <- runASCAT(
        lrr, baf, lrrsegm, bafsegm, ASCATobj$gender[arraynr], ASCATobj$SNPpos, ASCATobj$ch, ASCATobj$chrs, ASCATobj$sexchromosomes, failedqualitycheck,
        paste(ASCATobj$samples[arraynr], ".sunrise.png", sep = ""), paste(ASCATobj$samples[arraynr], ".ASCATprofile.png", sep = ""),
        paste(ASCATobj$samples[arraynr], ".rawprofile.png", sep = ""), paste(ASCATobj$samples[arraynr], ".aberrationreliability.png", sep = ""),
        gamma
      )
    } else {
      res[[arraynr]] <- runASCAT(
        lrr, baf, lrrsegm, bafsegm, ASCATobj$gender[arraynr], ASCATobj$SNPpos, ASCATobj$ch, ASCATobj$chrs, ASCATobj$sexchromosomes, failedqualitycheck,
        paste(ASCATobj$samples[arraynr], ".sunrise.png", sep = ""), paste(ASCATobj$samples[arraynr], ".ASCATprofile.png", sep = ""),
        paste(ASCATobj$samples[arraynr], ".rawprofile.png", sep = ""), paste(ASCATobj$samples[arraynr], ".aberrationreliability.png", sep = ""),
        gamma, rho_manual, psi_manual
      )
    }
    if (!is.na(res[[arraynr]]$rho)) {
      goodarrays[length(goodarrays) + 1] <- arraynr
    }
  }

  if (length(goodarrays) > 0) {
    n1 <- matrix(nrow = dim(ASCATobj$Tumor_LogR)[1], ncol = length(goodarrays))
    n2 <- matrix(nrow = dim(ASCATobj$Tumor_LogR)[1], ncol = length(goodarrays))
    rownames(n1) <- rownames(ASCATobj$Tumor_LogR)
    rownames(n2) <- rownames(ASCATobj$Tumor_LogR)
    colnames(n1) <- colnames(ASCATobj$Tumor_LogR)[goodarrays]
    colnames(n2) <- colnames(ASCATobj$Tumor_LogR)[goodarrays]
    for (i in 1:length(goodarrays)) {
      n1[, i] <- res[[goodarrays[i]]]$nA
      n2[, i] <- res[[goodarrays[i]]]$nB
    }

    tp <- vector(length = length(goodarrays))
    psi <- vector(length = length(goodarrays))
    ploidy <- vector(length = length(goodarrays))
    goodnessOfFit <- vector(length = length(goodarrays))
    naarrays <- NULL
    for (i in 1:length(goodarrays)) {
      tp[i] <- res[[goodarrays[i]]]$rho
      psi[i] <- res[[goodarrays[i]]]$psi
      ploidy[i] <- mean(res[[goodarrays[i]]]$nA + res[[goodarrays[i]]]$nB, na.rm = T)
      goodnessOfFit[i] <- res[[goodarrays[i]]]$goodnessOfFit
      if (res[[goodarrays[i]]]$nonaberrant) {
        naarrays <- c(naarrays, ASCATobj$samples[goodarrays[i]])
      }
    }
    fa <- colnames(ASCATobj$Tumor_LogR)[-goodarrays]
    names(tp) <- colnames(n1)
    names(ploidy) <- colnames(n1)
    names(psi) <- colnames(n1)
    names(goodnessOfFit) <- colnames(n1)

    seg <- NULL
    for (i in 1:length(goodarrays)) {
      segje <- res[[goodarrays[i]]]$seg
      seg <- rbind(seg, cbind(
        ASCATobj$samples[goodarrays[i]], as.vector(ASCATobj$SNPpos[segje[, 1], 1]),
        ASCATobj$SNPpos[segje[, 1], 2],
        ASCATobj$SNPpos[segje[, 2], 2], segje[, 3], segje[, 4]
      ))
    }
    colnames(seg) <- c("sample", "chr", "startpos", "endpos", "nMajor", "nMinor")
    seg <- data.frame(seg, stringsAsFactors = F)
    seg[, 3] <- as.numeric(seg[, 3])
    seg[, 4] <- as.numeric(seg[, 4])
    seg[, 5] <- as.numeric(seg[, 5])
    seg[, 6] <- as.numeric(seg[, 6])

    seg_raw <- NULL
    for (i in 1:length(goodarrays)) {
      segje <- res[[goodarrays[i]]]$seg_raw
      seg_raw <- rbind(seg_raw, cbind(
        ASCATobj$samples[goodarrays[i]], as.vector(ASCATobj$SNPpos[segje[, 1], 1]),
        ASCATobj$SNPpos[segje[, 1], 2],
        ASCATobj$SNPpos[segje[, 2], 2], segje[, 3], segje[, 4:ncol(segje)]
      ))
    }
    colnames(seg_raw) <- c("sample", "chr", "startpos", "endpos", "nMajor", "nMinor", "nAraw", "nBraw")

    seg_raw <- data.frame(seg_raw, stringsAsFactors = F)
    seg_raw[, 3] <- as.numeric(seg_raw[, 3])
    seg_raw[, 4] <- as.numeric(seg_raw[, 4])
    seg_raw[, 5] <- as.numeric(seg_raw[, 5])
    seg_raw[, 6] <- as.numeric(seg_raw[, 6])
    seg_raw[, 7] <- as.numeric(seg_raw[, 7])
    seg_raw[, 8] <- as.numeric(seg_raw[, 8])
    # TBKW edit
    nArawfull <- res[[arraynr]]$nArawfull
    nBrawfull <- res[[arraynr]]$nBrawfull
  }
  else {
    n1 <- NULL
    n2 <- NULL
    tp <- NULL
    ploidy <- NULL
    psi <- NULL
    goodnessOfFit <- NULL
    fa <- colnames(ASCATobj$Tumor_LogR)
    naarrays <- NULL
    seg <- NULL
    seg_raw <- NULL
    # TBKWedit
    nArawfull <- NULL
    nBrawfull <- NULL
  }
  return(list(
    nA = n1, nB = n2, aberrantcellfraction = tp, ploidy = ploidy, psi = psi, goodnessOfFit = goodnessOfFit,
    failedarrays = fa, nonaberrantarrays = naarrays, segments = seg, segments_raw = seg_raw, nArawfull = nArawfull, nBrawfull = nBrawfull
  ))
}


# helper function to split the genome into parts
split_genome <- function(SNPpos) {

  # look for gaps of more than 5Mb (arbitrary treshold to account for big centremeres or other gaps) and chromosome borders
  bigHoles <- which(diff(SNPpos[, 2]) >= 5000000) + 1
  chrBorders <- which(SNPpos[1:(dim(SNPpos)[1] - 1), 1] != SNPpos[2:(dim(SNPpos)[1]), 1]) + 1

  holes <- unique(sort(c(bigHoles, chrBorders)))

  startseg <- c(1, holes)
  endseg <- c(holes - 1, dim(SNPpos)[1])

  chr <- list()
  for (i in 1:length(startseg)) {
    chr[[i]] <- startseg[i]:endseg[i]
  }

  return(chr)
}


# function to make segments of constant LRR and BAF
# this function is more general and does not depend on specifically ASPCF output
# it can also handle segmention performed on LRR and BAF separately
make_segments <- function(r, b) {
  m <- matrix(ncol = 2, nrow = length(b))
  m[, 1] <- r
  m[, 2] <- b
  m <- as.matrix(na.omit(m))
  pcf_segments <- matrix(ncol = 3, nrow = dim(m)[1])
  colnames(pcf_segments) <- c("r", "b", "length")
  index <- 0
  previousb <- -1
  previousr <- 1E10
  for (i in 1:dim(m)[1]) {
    if (m[i, 2] != previousb || m[i, 1] != previousr) {
      index <- index + 1
      count <- 1
      pcf_segments[index, "r"] <- m[i, 1]
      pcf_segments[index, "b"] <- m[i, 2]
    }
    else {
      count <- count + 1
    }
    pcf_segments[index, "length"] <- count
    previousb <- m[i, 2]
    previousr <- m[i, 1]
  }
  pcf_segments <- as.matrix(na.omit(pcf_segments))[, ]
  return(pcf_segments)
}

# function to create the distance matrix (distance for a range of ploidy and tumor percentage values)
# input: segmented LRR and BAF and the value for gamma
create_distance_matrix <- function(segments, gamma) {
  s <- segments
  psi_pos <- seq(1, 6, 0.05)
  #rho_pos <- seq(0.1, 1.05, 0.01)
  rho_pos <- seq(0.1, 1.00, 0.01)
  d <- matrix(nrow = length(psi_pos), ncol = length(rho_pos))
  rownames(d) <- psi_pos
  colnames(d) <- rho_pos
  dmin <- 1E20
  for (i in 1:length(psi_pos)) {
    psi <- psi_pos[i]
    for (j in 1:length(rho_pos)) {
      rho <- rho_pos[j]
      nA <- (rho - 1 - (s[, "b"] - 1) * 2^(s[, "r"] / gamma) * ((1 - rho) * 2 + rho * psi)) / rho
      nB <- (rho - 1 + s[, "b"] * 2^(s[, "r"] / gamma) * ((1 - rho) * 2 + rho * psi)) / rho
      # choose the minor allele
      nMinor <- NULL
      if (sum(nA, na.rm = T) < sum(nB, na.rm = T)) {
        nMinor <- nA
      }
      else {
        nMinor <- nB
      }
      d[i, j] <- sum(abs(nMinor - pmax(round(nMinor), 0))^2 * s[, "length"] * ifelse(s[, "b"] == 0.5, 0.05, 1), na.rm = T)
    }
  }
  return(d)
}


# the ASCAT main function
# lrr: (unsegmented) log R, in genomic sequence (all probes), with probe IDs
# baf: (unsegmented) B Allele Frequency, in genomic sequence (all probes), with probe IDs
# lrrsegmented: log R, segmented, in genomic sequence (all probes), with probe IDs
# bafsegmented: B Allele Frequency, segmented, in genomic sequence (only probes heterozygous in germline), with probe IDs
# chromosomes: a list containing c vectors, where c is the number of chromosomes and every vector contains all probe numbers per chromosome
# chrnames: a vector containing the names for the chromosomes (e.g. c(1:22,"X"))
# sexchromosomes: a vector containing the names for the sex chromosomes
# failedqualitycheck: did the sample fail any previous quality check or not?
# distancepng: if NA: distance is plotted, if filename is given, the plot is written to a .png file
# copynumberprofilespng: if NA: possible copy number profiles are plotted, if filename is given, the plot is written to a .png file
# nonroundedprofilepng: if NA: copy number profile before rounding is plotted (total copy number as well as the copy number of the minor allele), if filename is given, the plot is written to a .png file
# aberrationreliabilitypng: if NA: aberration reliability score is plotted, if filename is given, the plot is written to a .png file
# gamma: technology parameter, compaction of Log R profiles (expected decrease in case of deletion in diploid sample, 100 % aberrant cells; 1 in ideal case, 0.55 of Illumina 109K arrays)
runASCAT <- function(lrr, baf, lrrsegmented, bafsegmented, gender, SNPpos, chromosomes, chrnames, sexchromosomes, failedqualitycheck = F,
                     distancepng = NA, copynumberprofilespng = NA, nonroundedprofilepng = NA, aberrationreliabilitypng = NA, gamma = 0.55,
                     rho_manual = NA, psi_manual = NA) {
  ch <- chromosomes
  chrs <- chrnames
  b <- bafsegmented
  raw_b <- baf[as.numeric(attributes(bafsegmented)$names)]
  raw_b <- ifelse(raw_b > 0.5, raw_b, 1 - raw_b)
  raw_b <- 1 - raw_b
  #message("length(raw_b): ", length(raw_b))
  r <- lrrsegmented[names(bafsegmented)]

  library(RColorBrewer)

  SNPposhet <- SNPpos[names(bafsegmented), ]
  autoprobes <- !(SNPposhet[, 1] %in% sexchromosomes)

  b2 <- b[autoprobes]
  r2 <- r[autoprobes]

  s <- make_segments(r2, b2)
  d <- create_distance_matrix(s, gamma) # Write out to get manual psi/rho ### NJB
  sampleid <- sub(".aberrationreliability.png$", "", aberrationreliabilitypng) ### Get sample ID
  save(d, file = paste(sampleid, "distance_matrix.RData", sep = ".")) ## Saving distance matrix for future optimization of rhp & psi

  if (!is.na(distancepng)) {
    png(filename = distancepng, width = 1000, height = 1000, res = 1000 / 7)
  }

  par(mar = c(5, 5, 0.5, 0.5), cex = 0.75, cex.lab = 2, cex.axis = 2)

  hmcol <- rev(colorRampPalette(brewer.pal(10, "RdBu"))(256))
  ploidy_points <- as.numeric(rownames(d))
  purity_points <- as.numeric(colnames(d))
  image(x = ploidy_points, y = purity_points, z = log(d), col = hmcol, xlab = "Ploidy", ylab = "Aberrant cell fraction")

  #axis(1, at = ploidy_points, label = ploidy_points)
  #axis(2, at = purity_points, label = purity_points)

  #axis(1, at = seq(0, 1, by = 1 / 5), label = seq(1, 6, by = 1))
  ##axis(2, at = seq(0, 1 / 1.05, by = 1 / 3 / 1.05), label = seq(0.1, 1, by = 3 / 10))
  #axis(2, at = seq(0, 1, by = 1 / 3 ), label = seq(0.1, 1, by = 3 / 10))

  TheoretMaxdist <- sum(rep(0.25, dim(s)[1]) * s[, "length"] * ifelse(s[, "b"] == 0.5, 0.05, 1), na.rm = T)

  # flag the sample as non-aberrant if necessary
  nonaberrant <- F
  MINABB <- 0.03
  MINABBREGION <- 0.005

  percentAbb <- sum(ifelse(s[, "b"] == 0.5, 0, 1) * s[, "length"]) / sum(s[, "length"])
  maxsegAbb <- max(ifelse(s[, "b"] == 0.5, 0, s[, "length"])) / sum(s[, "length"])
  if (percentAbb <= MINABB & maxsegAbb <= MINABBREGION) {
    nonaberrant <- T
  }

  MINPLOIDY <- 1.5
  MAXPLOIDY <- 5.5
  MINRHO <- 0.05  # TBKW: set in XPhase this is set in  Changed this to 0.15 from 0.2, 2015-02-11. Got several ok samples that could go lower from TRACERx ### Changed again to 0.10, 2015-05-05
  #MAXRHO <- 0.95
  MAXRHO <- 0.85

  # Tweak the limits if we have passed in a range of manual rho values
  if (length(rho_manual) > 1) {
    MINRHO <- max(MINRHO, rho_manual[[1]])
    MAXRHO <- min(MAXRHO, rho_manual[[2]])
    # Reset rho_manual
    rho_manual <- NA
  }

  # Tweak the limits if we have passed in a range of manual psi/plidy values
  if (length(psi_manual) > 1) {
    MINPLOIDY <- max(MINPLOIDY, psi_manual[[1]])
    MAXPLOIDY <- min(MAXPLOIDY, psi_manual[[2]])
    # Reset psi_manual
    psi_manual <- NA
  }

  # Highlight the region that we consider to contain reasonable purity and
  # ploidy estimates
  rect(xleft = MINPLOIDY, ybottom = MINRHO, xright = MAXPLOIDY, ytop = MAXRHO,
       col = adjustcolor("yellow", 0.1), border = NA)

  MINGOODNESSOFFIT <- 80
  MINPERCZERO <- 0.02
  MINPERCZEROABB <- 0.1
  MINPERCODDEVEN <- 0.05
  MINPLOIDYSTRICT <- 1.7
  MAXPLOIDYSTRICT <- 2.3

  nropt <- 0
  localmin <- NULL
  optima <- list()

  if (!failedqualitycheck && is.na(rho_manual)) {

    ### NJB edit 2015-05-05, Saving alternative solutions
    inv.d <- 1 / d
    alt.ploidy <- cbind(as.numeric(rownames(inv.d)), apply(inv.d, 1, max))
    alt.cellularity <- cbind(as.numeric(colnames(inv.d)), apply(inv.d, 2, max))
    alt.p.indx <- which(diff(sign(diff(alt.ploidy[, 2]))) == -2) + 1
    alt.c.indx <- apply(inv.d[rownames(alt.ploidy[alt.p.indx, , drop = F]), , drop = F], 1, function(x) {
      which(x == max(x))
    })
    if (length(ascat.runAscat) > 0) {
      alt.sol <- cbind(alt.ploidy[alt.p.indx, 1, drop = F], alt.cellularity[alt.c.indx, , drop = F])
      colnames(alt.sol) <- c("Ploidy", "Cellularity", "Goodness_of_fit")
      alt.sol <- alt.sol[order(alt.sol[, 3], decreasing = T), , drop = F]
      rownames(alt.sol) <- 1:nrow(alt.sol)
      alt.sol[, 3] <- (1 - ((1 / alt.sol[, 3]) / TheoretMaxdist)) * 100
      write.table(alt.sol, file = paste(sampleid, "alternative_solutions.txt", sep = "."), quote = F, row.names = F, sep = "\t")
    }
    #### end edit

    # zeroth step, take global min inside the range of purity and ploidy values
    # that we want
    d_allowed <- d[as.numeric(rownames(d)) >= MINPLOIDY & as.numeric(rownames(d)) <= MAXPLOIDY,
                   as.numeric(colnames(d)) >= MINRHO & as.numeric(colnames(d)) <= MAXRHO]
    m <- min(d_allowed)
    min_index <- which(d_allowed == m, arr.ind = TRUE)
    i <- min_index[1, "row"]
    j <- min_index[1, "col"]
    psi <- as.numeric(rownames(d_allowed)[i])
    rho <- as.numeric(colnames(d_allowed)[j])
    nA <- (rho - 1 - (s[, "b"] - 1) * 2^(s[, "r"] / gamma) * ((1 - rho) * 2 + rho * psi)) / rho
    nB <- (rho - 1 + s[, "b"] * 2^(s[, "r"] / gamma) * ((1 - rho) * 2 + rho * psi)) / rho

    logdebug(paste("Absolute best m:", m))
    logdebug(paste("Absolute best psi:", psi, MINPLOIDY, MAXPLOIDY))
    logdebug(paste("Absolute best rho:", rho, MINRHO, MAXRHO))
    # ploidy is recalculated based on results, to avoid bias (due to differences in normalization of LogR)
    ploidy <- sum((nA + nB) * s[, "length"]) / sum(s[, "length"])
    percentzero <- (sum((round(nA) == 0) * s[, "length"]) + sum((round(nB) == 0) * s[, "length"])) / sum(s[, "length"])
    goodnessOfFit <- (1 - m / TheoretMaxdist) * 100
    logdebug(paste("ploidy:", ploidy, MINPLOIDY, MAXPLOIDY))
    logdebug(paste("percentzero:", percentzero, MINPERCZERO))
    logdebug(paste("goodnessOfFit:", goodnessOfFit, MINGOODNESSOFFIT))
    logdebug(paste("nonaberrant:", nonaberrant))
    #if (!nonaberrant & ploidy > MINPLOIDY & ploidy < MAXPLOIDY & rho >= MINRHO & goodnessOfFit > MINGOODNESSOFFIT & percentzero > MINPERCZERO) {
    #if (!nonaberrant & ploidy > MINPLOIDY & ploidy < MAXPLOIDY & rho >= MINRHO & goodnessOfFit > MINGOODNESSOFFIT) {
    if (!nonaberrant) {
      min_index <- which(d == m, arr.ind = TRUE)
      i <- min_index[1, "row"]
      j <- min_index[1, "col"]
      nropt <- nropt + 1
      optima[[nropt]] <- c(m, i, j, ploidy, goodnessOfFit)
      localmin[nropt] <- m
    }

    logdebug(paste0("0) Number of optimal solutions selected: ", nropt))

    # first, try with all filters
    # if no solution, drop the percentzero > MINPERCZERO filter (allow non-aberrant solutions - but limit the ploidy options)
    if (nropt == 0) {
      for (i in 4:(nrow(d) - 3)) {
        for (j in 4:(ncol(d) - 3)) {
          m <- d[i, j]
          seld <- d[(i - 3):(i + 3), (j - 3):(j + 3)]
          seld[4, 4] <- max(seld)
          if (min(seld) > m) {
            psi <- as.numeric(rownames(d)[i])
            rho <- as.numeric(colnames(d)[j])
            nA <- (rho - 1 - (s[, "b"] - 1) * 2^(s[, "r"] / gamma) * ((1 - rho) * 2 + rho * psi)) / rho
            nB <- (rho - 1 + s[, "b"] * 2^(s[, "r"] / gamma) * ((1 - rho) * 2 + rho * psi)) / rho

            # ploidy is recalculated based on results, to avoid bias (due to differences in normalization of LogR)
            ploidy <- sum((nA + nB) * s[, "length"]) / sum(s[, "length"])

            percentzero <- (sum((round(nA) == 0) * s[, "length"]) + sum((round(nB) == 0) * s[, "length"])) / sum(s[, "length"])

            goodnessOfFit <- (1 - m / TheoretMaxdist) * 100

            if (!nonaberrant & ploidy > MINPLOIDY & ploidy < MAXPLOIDY & rho >= MINRHO & goodnessOfFit > MINGOODNESSOFFIT & percentzero > MINPERCZERO) {
              nropt <- nropt + 1
              optima[[nropt]] <- c(m, i, j, ploidy, goodnessOfFit)
              localmin[nropt] <- m
            }
          }
        }
      }
    }

    logdebug(paste0("1) Number of optimal solutions selected: ", nropt))

    # if no solution, drop the percentzero > MINPERCZERO filter (allow non-aberrant solutions - but limit the ploidy options)
    if (nropt == 0) {
      for (i in 4:(dim(d)[1] - 3)) {
        for (j in 4:(dim(d)[2] - 3)) {
          m <- d[i, j]
          seld <- d[(i - 3):(i + 3), (j - 3):(j + 3)]
          seld[4, 4] <- max(seld)
          if (min(seld) > m) {
            psi <- as.numeric(rownames(d)[i])
            rho <- as.numeric(colnames(d)[j])
            nA <- (rho - 1 - (s[, "b"] - 1) * 2^(s[, "r"] / gamma) * ((1 - rho) * 2 + rho * psi)) / rho
            nB <- (rho - 1 + s[, "b"] * 2^(s[, "r"] / gamma) * ((1 - rho) * 2 + rho * psi)) / rho

            # ploidy is recalculated based on results, to avoid bias (due to differences in normalization of LogR)
            ploidy <- sum((nA + nB) * s[, "length"]) / sum(s[, "length"])

            perczeroAbb <- (sum((round(nA) == 0) * s[, "length"] * ifelse(s[, "b"] == 0.5, 0, 1)) + sum((round(nB) == 0) * s[, "length"] * ifelse(s[, "b"] == 0.5, 0, 1))) / sum(s[, "length"] * ifelse(s[, "b"] == 0.5, 0, 1))
            # the next can happen if BAF is a flat line at 0.5
            if (is.na(perczeroAbb)) {
              perczeroAbb <- 0
            }

            goodnessOfFit <- (1 - m / TheoretMaxdist) * 100

            if (ploidy > MINPLOIDYSTRICT & ploidy < MAXPLOIDYSTRICT & rho >= MINRHO & goodnessOfFit > MINGOODNESSOFFIT & perczeroAbb > MINPERCZEROABB) {
              nropt <- nropt + 1
              optima[[nropt]] <- c(m, i, j, ploidy, goodnessOfFit)
              localmin[nropt] <- m
            }
          }
        }
      }
    }
    logdebug(paste0("2) Number of optimal solutions selected: ", nropt))

    # if still no solution, allow solutions with 100% aberrant cells (include the borders with rho = 1), but in first instance, keep the percentzero > 0.01 filter
    if (nropt == 0) {
      # first, include borders
      cold <- which(as.numeric(colnames(d)) > 1)
      d[, cold] <- 1E20
      for (i in 4:(dim(d)[1] - 3)) {
        for (j in 4:(dim(d)[2] - 3)) {
          m <- d[i, j]
          seld <- d[(i - 3):(i + 3), (j - 3):(j + 3)]
          seld[4, 4] <- max(seld)
          if (min(seld) > m) {
            psi <- as.numeric(rownames(d)[i])
            rho <- as.numeric(colnames(d)[j])
            nA <- (rho - 1 - (s[, "b"] - 1) * 2^(s[, "r"] / gamma) * ((1 - rho) * 2 + rho * psi)) / rho
            nB <- (rho - 1 + s[, "b"] * 2^(s[, "r"] / gamma) * ((1 - rho) * 2 + rho * psi)) / rho

            # ploidy is recalculated based on results, to avoid bias (due to differences in normalization of LogR)
            ploidy <- sum((nA + nB) * s[, "length"]) / sum(s[, "length"])

            percentzero <- (sum((round(nA) == 0) * s[, "length"]) + sum((round(nB) == 0) * s[, "length"])) / sum(s[, "length"])
            percOddEven <- sum((round(nA) %% 2 == 0 & round(nB) %% 2 == 1 | round(nA) %% 2 == 1 & round(nB) %% 2 == 0) * s[, "length"]) / sum(s[, "length"])
            perczeroAbb <- (sum((round(nA) == 0) * s[, "length"] * ifelse(s[, "b"] == 0.5, 0, 1)) + sum((round(nB) == 0) * s[, "length"] * ifelse(s[, "b"] == 0.5, 0, 1))) / sum(s[, "length"] * ifelse(s[, "b"] == 0.5, 0, 1))
            if (is.na(perczeroAbb)) {
              perczeroAbb <- 0
            }

            goodnessOfFit <- (1 - m / TheoretMaxdist) * 100

            if (!nonaberrant & ploidy > MINPLOIDY & ploidy < MAXPLOIDY & rho >= MINRHO & goodnessOfFit > MINGOODNESSOFFIT &
              (perczeroAbb > MINPERCZEROABB | percentzero > MINPERCZERO | percOddEven > MINPERCODDEVEN)) {
              nropt <- nropt + 1
              optima[[nropt]] <- c(m, i, j, ploidy, goodnessOfFit)
              localmin[nropt] <- m
            }
          }
        }
      }
    }
    logdebug(paste0("3) Number of optimal solutions selected: ", nropt))

    # if still no solution, drop the percentzero > MINPERCENTZERO filter, but strict ploidy borders
    if (nropt == 0) {
      for (i in 4:(dim(d)[1] - 3)) {
        for (j in 4:(dim(d)[2] - 3)) {
          m <- d[i, j]
          seld <- d[(i - 3):(i + 3), (j - 3):(j + 3)]
          seld[4, 4] <- max(seld)
          if (min(seld) > m) {
            psi <- as.numeric(rownames(d)[i])
            rho <- as.numeric(colnames(d)[j])
            nA <- (rho - 1 - (s[, "b"] - 1) * 2^(s[, "r"] / gamma) * ((1 - rho) * 2 + rho * psi)) / rho
            nB <- (rho - 1 + s[, "b"] * 2^(s[, "r"] / gamma) * ((1 - rho) * 2 + rho * psi)) / rho

            # ploidy is recalculated based on results, to avoid bias (due to differences in normalization of LogR)
            ploidy <- sum((nA + nB) * s[, "length"]) / sum(s[, "length"])

            perczeroAbb <- (sum((round(nA) == 0) * s[, "length"] * ifelse(s[, "b"] == 0.5, 0, 1)) + sum((round(nB) == 0) * s[, "length"] * ifelse(s[, "b"] == 0.5, 0, 1))) / sum(s[, "length"] * ifelse(s[, "b"] == 0.5, 0, 1))
            # the next can happen if BAF is a flat line at 0.5
            if (is.na(perczeroAbb)) {
              perczeroAbb <- 0
            }

            goodnessOfFit <- (1 - m / TheoretMaxdist) * 100

            if (ploidy > MINPLOIDYSTRICT & ploidy < MAXPLOIDYSTRICT & rho >= MINRHO & goodnessOfFit > MINGOODNESSOFFIT) {
              nropt <- nropt + 1
              optima[[nropt]] <- c(m, i, j, ploidy, goodnessOfFit)
              localmin[nropt] <- m
            }
          }
        }
      }
    }
    logdebug(paste0("4) Number of optimal solutions selected: ", nropt))
  }

  if (!is.na(rho_manual)) {
    rho <- rho_manual
    psi <- psi_manual

    nA <- (rho - 1 - (s[, "b"] - 1) * 2^(s[, "r"] / gamma) * ((1 - rho) * 2 + rho * psi)) / rho
    nB <- (rho - 1 + s[, "b"] * 2^(s[, "r"] / gamma) * ((1 - rho) * 2 + rho * psi)) / rho

    # ploidy is recalculated based on results, to avoid bias (due to differences in normalization of LogR)
    ploidy <- sum((nA + nB) * s[, "length"]) / sum(s[, "length"])

    nMinor <- NULL
    if (sum(nA, na.rm = T) < sum(nB, na.rm = T)) {
      nMinor <- nA
    }
    else {
      nMinor <- nB
    }
    m <- sum(abs(nMinor - pmax(round(nMinor), 0))^2 * s[, "length"] * ifelse(s[, "b"] == 0.5, 0.05, 1), na.rm = T)

    goodnessOfFit <- (1 - m / TheoretMaxdist) * 100

    nropt <- 1
    optima[[1]] <- c(m, rho, psi, ploidy, goodnessOfFit)
    localmin[1] <- m
  }


  if (nropt > 0) {
    if (is.na(rho_manual)) {
      optlim <- sort(localmin)[1]
      for (i in 1:length(optima)) {
        if (optima[[i]][1] == optlim) {
          psi_opt1 <- as.numeric(rownames(d)[optima[[i]][2]])
          rho_opt1 <- as.numeric(colnames(d)[optima[[i]][3]])
          if (rho_opt1 > 1) {
            rho_opt1 <- 1
          }
          ploidy_opt1 <- optima[[i]][4]
          goodnessOfFit_opt1 <- optima[[i]][5]
          #points((psi_opt1 - 1) / 5, (rho_opt1 - 0.1) / 0.95, col = "green", pch = "X", cex = 2)
          points(psi_opt1, rho_opt1, col = "green", pch = 4, cex = 2)
        }
      }
    } else {
      rho_opt1 <- optima[[1]][2]
      psi_opt1 <- optima[[1]][3]
      ploidy_opt1 <- optima[[1]][4]
      goodnessOfFit_opt1 <- optima[[1]][5]
      #points((psi_opt1 - 1) / 5, (rho_opt1 - 0.1) / 0.95, col = "green", pch = "X", cex = 2)
      points(psi_opt1, rho_opt1, col = "green", pch = 4, cex = 2)
    }
  }

  if (!is.na(distancepng)) {
    dev.off()
  }


  if (nropt > 0) {
    if (!is.na(nonroundedprofilepng)) {
      png(filename = nonroundedprofilepng, width = 2000, height = 500, res = 200)
    }
    else {
      windows(10, 5)
    }

    par(mar = c(0.5, 5, 5, 0.5), cex = 0.4, cex.main = 3, cex.axis = 2.5)

    rho <- rho_opt1
    psi <- psi_opt1

    SNPposhet <- SNPpos[names(bafsegmented), ]
    haploidchrs <- unique(c(substring(gender, 1, 1), substring(gender, 2, 2)))
    if (substring(gender, 1, 1) == substring(gender, 2, 2)) {
      haploidchrs <- setdiff(haploidchrs, substring(gender, 1, 1))
    }
    diploidprobes <- !(SNPposhet[, 1] %in% haploidchrs)
    nullchrs <- setdiff(sexchromosomes, unique(c(substring(gender, 1, 1), substring(gender, 2, 2))))
    nullprobes <- SNPposhet[, 1] %in% nullchrs

    nAfull <- ifelse(diploidprobes,
      (rho - 1 - (b - 1) * 2^(r / gamma) * ((1 - rho) * 2 + rho * psi)) / rho,
      ifelse(nullprobes, 0,
        ifelse(b < 0.5, (rho - 1 + ((1 - rho) * 2 + rho * psi) * 2^(r / gamma)) / rho, 0)
      )
    )
    nBfull <- ifelse(diploidprobes,
      (rho - 1 + b * 2^(r / gamma) * ((1 - rho) * 2 + rho * psi)) / rho,
      ifelse(nullprobes, 0,
        ifelse(b < 0.5, 0, (rho - 1 + ((1 - rho) * 2 + rho * psi) * 2^(r / gamma)) / rho)
      )
    )
    # TBKW edit...
    #message("length(raw_b): ", length(raw_b))
    #message("length(b): ", length(b))
    nArawfull <- ifelse(diploidprobes,
      (rho - 1 - (raw_b - 1) * 2^(r / gamma) * ((1 - rho) * 2 + rho * psi)) / rho,
      ifelse(nullprobes, 0,
        ifelse(raw_b < 0.5, (rho - 1 + ((1 - rho) * 2 + rho * psi) * 2^(r / gamma)) / rho, 0)
      )
    )
    nBrawfull <- ifelse(diploidprobes,
      (rho - 1 + raw_b * 2^(r / gamma) * ((1 - rho) * 2 + rho * psi)) / rho,
      ifelse(nullprobes, 0,
        ifelse(raw_b < 0.5, 0, (rho - 1 + ((1 - rho) * 2 + rho * psi) * 2^(r / gamma)) / rho)
      )
    )

    # TBKW edit...
    save(nA, file = paste0(sampleid, ".nA.Rdata"))
    save(nB, file = paste0(sampleid, ".nB.Rdata"))
    save(nArawfull, file = paste0(sampleid, ".nArawfull.Rdata"))
    save(nBrawfull, file = paste0(sampleid, ".nBrawfull.Rdata"))
    save(nBfull, file = paste0(sampleid, ".nBfull.Rdata"))
    save(nAfull, file = paste0(sampleid, ".nAfull.Rdata"))

    nA <- pmax(round(nAfull), 0)
    nB <- pmax(round(nBfull), 0)

    maintitle <- paste("Ploidy: ", sprintf("%1.2f", ploidy_opt1), ", aberrant cell fraction: ", sprintf("%2.0f", rho_opt1 * 100), "%, goodness of fit: ", sprintf("%2.1f", goodnessOfFit_opt1), "%", ifelse(nonaberrant, ", non-aberrant", ""), sep = "")
    plot(c(1, length(nAfull)), c(0, 5), type = "n", xaxt = "n", main = maintitle, xlab = "", ylab = "")
    points(nBfull, col = "blue", pch = "|")
    points(nAfull + nBfull, col = "purple", pch = "|")
    # don't ask me why, but the "|" sign is not centered, so the lines may need to be shifted..
    abline(v = 0, lty = 1, col = "lightgrey")
    chrk_tot_len <- 0
    for (i in 1:length(ch)) {
      chrk <- ch[[i]]
      chrk_hetero <- intersect(names(lrr)[chrk], names(bafsegmented))
      chrk_tot_len_prev <- chrk_tot_len
      chrk_tot_len <- chrk_tot_len + length(chrk_hetero)
      vpos <- chrk_tot_len
      tpos <- (chrk_tot_len + chrk_tot_len_prev) / 2
      text(tpos, 5, chrs[i], pos = 1, cex = 2)
      abline(v = vpos, lty = 1, col = "lightgrey")
    }

    if (!is.na(nonroundedprofilepng)) {
      dev.off()
    }

    rho <- rho_opt1
    psi <- psi_opt1

    diploidprobes <- !(SNPpos[, 1] %in% haploidchrs)
    nullprobes <- SNPpos[, 1] %in% nullchrs

    # this replaces an occurrence of unique that caused problems
    tlr2 <- rle(lrrsegmented)
    tlr <- tlr2$values
    tlrstart <- c(1, cumsum(tlr2$lengths) + 1)
    tlrstart <- tlrstart[1:(length(tlrstart) - 1)]
    tlrend <- cumsum(tlr2$lengths)
    seg <- NULL
    for (i in 1:length(tlr)) {
      logR <- tlr[i]
      # pr = which(lrrsegmented==logR) # this was a problem
      pr <- tlrstart[i]:tlrend[i]
      start <- min(pr)
      end <- max(pr)
      bafke <- bafsegmented[intersect(names(lrrsegmented)[pr], names(bafsegmented))][1]
      # if bafke is NA, this means that we are dealing with a germline homozygous stretch with a copy number change within it.
      # in this case, nA and nB are irrelevant, just their sum matters
      if (is.na(bafke)) {
        bafke <- 0
      }

      nAraw <- ifelse(diploidprobes[start],
        (rho - 1 - (bafke - 1) * 2^(logR / gamma) * ((1 - rho) * 2 + rho * psi)) / rho,
        ifelse(nullprobes[start], 0,
          (rho - 1 + ((1 - rho) * 2 + rho * psi) * 2^(logR / gamma)) / rho
        )
      )
      nBraw <- ifelse(diploidprobes[start], (rho - 1 + bafke * 2^(logR / gamma) * ((1 - rho) * 2 + rho * psi)) / rho, 0)
      # correct for negative values:
      if (nAraw + nBraw < 0) {
        nAraw <- 0
        nBraw <- 0
      }
      else if (nAraw < 0) {
        nBraw <- nAraw + nBraw
        nAraw <- 0
      }
      else if (nBraw < 0) {
        nAraw <- nAraw + nBraw
        nBraw <- 0
      }
      # when evidence for odd copy number in segments of BAF = 0.5, assume a deviation..
      limitround <- 0.5
      nA <- ifelse(bafke == 0.5,
        ifelse(nAraw + nBraw > round(nAraw) + round(nBraw) + limitround,
          round(nAraw) + 1,
          ifelse(nAraw + nBraw < round(nAraw) + round(nBraw) - limitround,
            round(nAraw),
            round(nAraw)
          )
        ),
        round(nAraw)
      )
      nB <- ifelse(bafke == 0.5,
        ifelse(nAraw + nBraw > round(nAraw) + round(nBraw) + limitround,
          round(nBraw),
          ifelse(nAraw + nBraw < round(nAraw) + round(nBraw) - limitround,
            round(nBraw) - 1,
            round(nBraw)
          )
        ),
        round(nBraw)
      )
      if (is.null(seg)) {
        seg <- t(as.matrix(c(start, end, nA, nB)))
        seg_raw <- t(as.matrix(c(start, end, nA, nB, nAraw, nBraw)))
      }
      else {
        seg <- rbind(seg, c(start, end, nA, nB))
        seg_raw <- rbind(seg_raw, c(start, end, nA, nB, nAraw, nBraw))
      }
    }
    colnames(seg) <- c("start", "end", "nA", "nB")
    colnames(seg_raw) <- c("start", "end", "nA", "nB", "nAraw", "nBraw")

    # build helping vector
    chrhelp <- vector(length = length(lrrsegmented))
    for (chrnr in 1:length(ch)) {
      chrke <- ch[[chrnr]]
      chrhelp[chrke] <- chrnr
    }

    # every repeat joins 2 ends. 20 repeats will join about 1 million ends..
    for (rep in 1:20) {
      seg2 <- seg
      seg <- NULL
      skipnext <- F
      for (i in 1:dim(seg2)[1]) {
        if (!skipnext) {
          if (i != dim(seg2)[1] && seg2[i, "nA"] == seg2[i + 1, "nA"] && seg2[i, "nB"] == seg2[i + 1, "nB"] &&
            chrhelp[seg2[i, "end"]] == chrhelp[seg2[i + 1, "start"]]) {
            segline <- c(seg2[i, "start"], seg2[i + 1, "end"], seg2[i, 3:4])
            skipnext <- T
          }
          else {
            segline <- seg2[i, ]
          }

          if (is.null(seg)) {
            seg <- t(as.matrix(segline))
          }
          else {
            seg <- rbind(seg, segline)
          }
        }
        else {
          skipnext <- F
        }
      }
      colnames(seg) <- colnames(seg2)
    }
    rownames(seg) <- NULL

    nMajor <- vector(length = length(lrrsegmented))
    names(nMajor) <- names(lrrsegmented)
    nMinor <- vector(length = length(lrrsegmented))
    names(nMinor) <- names(lrrsegmented)

    for (i in 1:dim(seg)[1]) {
      nMajor[seg[i, "start"]:seg[i, "end"]] <- seg[i, "nA"]
      nMinor[seg[i, "start"]:seg[i, "end"]] <- seg[i, "nB"]
    }

    n1all <- vector(length = length(lrrsegmented))
    names(n1all) <- names(lrrsegmented)
    n2all <- vector(length = length(lrrsegmented))
    names(n2all) <- names(lrrsegmented)

    # note: any of these can have length 0
    NAprobes <- which(is.na(lrr))
    heteroprobes <- setdiff(which(names(lrrsegmented) %in% names(bafsegmented)), NAprobes)
    homoprobes <- setdiff(setdiff(which(!is.na(baf)), heteroprobes), NAprobes)
    CNprobes <- setdiff(which(is.na(baf)), NAprobes)

    n1all[NAprobes] <- NA
    n2all[NAprobes] <- NA
    n1all[CNprobes] <- nMajor[CNprobes] + nMinor[CNprobes]
    n2all[CNprobes] <- 0
    heteroprobes2 <- names(lrrsegmented)[heteroprobes]
    n1all[heteroprobes] <- ifelse(baf[heteroprobes2] <= 0.5, nMajor[heteroprobes], nMinor[heteroprobes])
    n2all[heteroprobes] <- ifelse(baf[heteroprobes2] > 0.5, nMajor[heteroprobes], nMinor[heteroprobes])
    n1all[homoprobes] <- ifelse(baf[homoprobes] <= 0.5, nMajor[homoprobes] + nMinor[homoprobes], 0)
    n2all[homoprobes] <- ifelse(baf[homoprobes] > 0.5, nMajor[homoprobes] + nMinor[homoprobes], 0)

    # plot ASCAT profile
    if (!is.na(copynumberprofilespng)) {
      png(filename = copynumberprofilespng, width = 2000, height = 500, res = 200)
    } else {
      windows(10, 2.5)
    }

    par(mar = c(0.5, 5, 5, 0.5), cex = 0.4, cex.main = 3, cex.axis = 2.5)

    nA2 <- n1all[heteroprobes]
    nB2 <- n2all[heteroprobes]
    nA <- ifelse(nA2 > nB2, nA2, nB2)
    nB <- ifelse(nA2 > nB2, nB2, nA2)
    maintitle <- paste("Ploidy: ", sprintf("%1.2f", ploidy_opt1), ", aberrant cell fraction: ", sprintf("%2.0f", rho_opt1 * 100), "%, goodness of fit: ", sprintf("%2.1f", goodnessOfFit_opt1), "%", ifelse(nonaberrant, ", non-aberrant", ""), sep = "")
    plot(c(1, length(nAfull)), c(0, 5), type = "n", xaxt = "n", main = maintitle, xlab = "", ylab = "")
    points(nA - 0.1, col = "red", pch = "|")
    points(nB + 0.1, col = "green", pch = "|")
    # don't ask me why, but the "|" sign is not centered, so the lines may need to be shifted..
    abline(v = 0, lty = 1, col = "lightgrey")
    chrk_tot_len <- 0
    for (i in 1:length(ch)) {
      chrk <- ch[[i]]
      chrk_hetero <- intersect(names(lrr)[chrk], names(bafsegmented))
      chrk_tot_len_prev <- chrk_tot_len
      chrk_tot_len <- chrk_tot_len + length(chrk_hetero)
      vpos <- chrk_tot_len
      tpos <- (chrk_tot_len + chrk_tot_len_prev) / 2
      text(tpos, 5, chrs[i], pos = 1, cex = 2)
      abline(v = vpos, lty = 1, col = "lightgrey")
    }

    if (!is.na(copynumberprofilespng)) {
      dev.off()
    }

    if (!is.na(aberrationreliabilitypng)) {
      png(filename = aberrationreliabilitypng, width = 2000, height = 500, res = 200)
    } else {
      windows(10, 2.5)
    }

    par(mar = c(0.5, 5, 5, 0.5), cex = 0.4, cex.main = 3, cex.axis = 2.5)

    diploidprobes <- !(SNPposhet[, 1] %in% haploidchrs)
    nullprobes <- SNPposhet[, 1] %in% nullchrs

    rBacktransform <- ifelse(diploidprobes,
      gamma * log((rho * (nA + nB) + (1 - rho) * 2) / ((1 - rho) * 2 + rho * psi), 2),
      # the value for nullprobes is arbitrary (but doesn't matter, as these are not plotted anyway because BAF=0.5)
      ifelse(nullprobes, -10, gamma * log((rho * (nA + nB) + (1 - rho)) / ((1 - rho) * 2 + rho * psi), 2))
    )

    bBacktransform <- ifelse(diploidprobes,
      (1 - rho + rho * nB) / (2 - 2 * rho + rho * (nA + nB)),
      ifelse(nullprobes, 0.5, 0)
    )

    rConf <- ifelse(abs(rBacktransform) > 0.15, pmin(100, pmax(0, 100 * (1 - abs(rBacktransform - r) / abs(r)))), NA)
    bConf <- ifelse(diploidprobes & bBacktransform != 0.5, pmin(100, pmax(0, ifelse(b == 0.5, 100, 100 * (1 - abs(bBacktransform - b) / abs(b - 0.5))))), NA)
    confidence <- ifelse(is.na(rConf), bConf, ifelse(is.na(bConf), rConf, (rConf + bConf) / 2))
    maintitle <- paste("Aberration reliability score (%), average: ", sprintf("%2.0f", mean(confidence, na.rm = T)), "%", sep = "")
    plot(c(1, length(nAfull)), c(0, 100), type = "n", xaxt = "n", main = maintitle, xlab = "", ylab = "")
    points(confidence, col = "blue", pch = "|")
    abline(v = 0, lty = 1, col = "lightgrey")
    chrk_tot_len <- 0
    for (i in 1:length(ch)) {
      chrk <- ch[[i]]
      chrk_hetero <- intersect(names(lrr)[chrk], names(bafsegmented))
      chrk_tot_len_prev <- chrk_tot_len
      chrk_tot_len <- chrk_tot_len + length(chrk_hetero)
      vpos <- chrk_tot_len
      tpos <- (chrk_tot_len + chrk_tot_len_prev) / 2
      text(tpos, 5, chrs[i], pos = 1, cex = 2)
      abline(v = vpos, lty = 1, col = "lightgrey")
    }

    if (!is.na(aberrationreliabilitypng)) {
      dev.off()
    }

    return(list(
      rho = rho_opt1, psi = psi_opt1, goodnessOfFit = goodnessOfFit_opt1,
      nonaberrant = nonaberrant, nA = n1all, nB = n2all, seg = seg,
      seg_raw = seg_raw, nArawfull = data.frame(nArawfull = nArawfull, stringsAsFactors = FALSE),
      nBrawfull = data.frame(nBrawfull = nBrawfull, stringsAsFactors = FALSE)
    ))
  } else {
    return(list(rho = NA, psi = NA, goodnessOfFit = NA, nonaberrant = F, nA = NA, nB = NA, seg = NA, seg_raw = NA))
  }
}


aspcfpart <- function(logRpart, allBflip, a, b, d, sd1, sd2, N, kmin, gamma) {
  from <- max(c(1, a))
  usefrom <- max(c(1, a + d))
  useto <- min(c(N, b - d))

  N <- length(logRpart)
  y1 <- logRpart
  y2 <- allBflip

  # Check that vectors are long enough to run algorithm:
  if (N < 2 * kmin) {
    breakpts <- 0
    return(list(breakpts = breakpts))
  }

  # Find initSum, initKvad, initAve for segment y[1..kmin]
  initSum1 <- sum(y1[1:kmin])
  initKvad1 <- sum(y1[1:kmin]^2)
  initAve1 <- initSum1 / kmin
  initSum2 <- sum(y2[1:kmin])
  initKvad2 <- sum(y2[1:kmin]^2)
  initAve2 <- initSum2 / kmin

  # Define vector of best costs
  bestCost <- rep(0, N)
  cost1 <- (initKvad1 - initSum1 * initAve1) / sd1^2
  cost2 <- (initKvad2 - initSum2 * initAve2) / sd2^2
  bestCost[kmin] <- cost1 + cost2

  # Define vector of best splits
  bestSplit <- rep(0, N)

  # Define vector of best averages
  bestAver1 <- rep(0, N)
  bestAver2 <- rep(0, N)
  bestAver1[kmin] <- initAve1
  bestAver2[kmin] <- initAve2

  # Initialize
  Sum1 <- rep(0, N)
  Sum2 <- rep(0, N)
  Kvad1 <- rep(0, N)
  Kvad2 <- rep(0, N)
  Aver1 <- rep(0, N)
  Aver2 <- rep(0, N)
  Cost <- rep(0, N)

  # We have to treat the region y(1..2*kmin-1) separately, as it
  # cannot be split into two full segments
  kminP1 <- kmin + 1
  for (k in (kminP1):(2 * kmin - 1)) {
    Sum1[kminP1:k] <- Sum1[kminP1:k] + y1[k]
    Aver1[kminP1:k] <- Sum1[kminP1:k] / ((k - kmin):1)
    Kvad1[kminP1:k] <- Kvad1[kminP1:k] + y1[k]^2
    Sum2[kminP1:k] <- Sum2[kminP1:k] + y2[k]
    Aver2[kminP1:k] <- Sum2[kminP1:k] / ((k - kmin):1)
    Kvad2[kminP1:k] <- Kvad2[kminP1:k] + y2[k]^2


    bestAver1[k] <- (initSum1 + Sum1[kminP1]) / k
    bestAver2[k] <- (initSum2 + Sum2[kminP1]) / k
    cost1 <- ((initKvad1 + Kvad1[kminP1]) - k * bestAver1[k]^2) / sd1^2
    cost2 <- ((initKvad2 + Kvad2[kminP1]) - k * bestAver2[k]^2) / sd2^2

    bestCost[k] <- cost1 + cost2
  }

  for (n in (2 * kmin):N) {
    nMkminP1 <- n - kmin + 1

    Sum1[kminP1:n] <- Sum1[kminP1:n] + y1[n]
    Aver1[kminP1:n] <- Sum1[kminP1:n] / ((n - kmin):1)
    Kvad1[kminP1:n] <- Kvad1[kminP1:n] + (y1[n])^2

    cost1 <- (Kvad1[kminP1:nMkminP1] - Sum1[kminP1:nMkminP1] * Aver1[kminP1:nMkminP1]) / sd1^2

    Sum2[kminP1:n] <- Sum2[kminP1:n] + y2[n]
    Aver2[kminP1:n] <- Sum2[kminP1:n] / ((n - kmin):1)
    Kvad2[kminP1:n] <- Kvad2[kminP1:n] + (y2[n])^2
    cost2 <- (Kvad2[kminP1:nMkminP1] - Sum2[kminP1:nMkminP1] * Aver2[kminP1:nMkminP1]) / sd2^2

    Cost[kminP1:nMkminP1] <- bestCost[kmin:(n - kmin)] + cost1 + cost2

    Pos <- which.min(Cost[kminP1:nMkminP1]) + kmin
    cost <- Cost[Pos] + gamma

    aver1 <- Aver1[Pos]
    aver2 <- Aver2[Pos]
    totAver1 <- (Sum1[kminP1] + initSum1) / n
    totCost1 <- ((Kvad1[kminP1] + initKvad1) - n * totAver1 * totAver1) / sd1^2
    totAver2 <- (Sum2[kminP1] + initSum2) / n
    totCost2 <- ((Kvad2[kminP1] + initKvad2) - n * totAver2 * totAver2) / sd2^2
    totCost <- totCost1 + totCost2

    if (totCost < cost) {
      Pos <- 1
      cost <- totCost
      aver1 <- totAver1
      aver2 <- totAver2
    }
    bestCost[n] <- cost
    bestAver1[n] <- aver1
    bestAver2[n] <- aver2
    bestSplit[n] <- Pos - 1
  }

  # Trace back
  n <- N
  breakpts <- n
  while (n > 0) {
    breakpts <- c(bestSplit[n], breakpts)
    n <- bestSplit[n]
  }

  breakpts <- breakpts + from - 1
  breakpts <- breakpts[breakpts >= usefrom & breakpts <= useto]

  return(list(breakpts = breakpts))
}


div <- function(a, b, c) {
  if (nargs() < 3) {
    c <- 0
  }

  if (b > 0) {
    v <- a / b
  } else {
    v <- c
  }

  return(v)
}


psi <- function(x, z) {
  xwin <- x
  xwin[x < -z] <- -z
  xwin[x > z] <- z
  return(xwin)
}
