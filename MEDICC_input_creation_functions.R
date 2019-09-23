#Needs MCS segmentation for the set of samples passed to MEDICC.
makeMEDICCInputListsFromSegmentationDF <- function(segs.to.make.MEDICC.input.df, sample.vec, chrs.to.use.vec, chr.col, start.col, end.col, cpnA.col, cpnB.col) {

	#To do this we need to make it match the
	MEDICC.FASTA.file.list <- list()
	for (chr in CHRS_TO_USE_VEC) {

		MEDICC.FASTA.file.list[[chr]] <- list()
		MEDICC.FASTA.file.list[[chr]]$chrom <- chr
		MEDICC.FASTA.file.list[[chr]][["diploid"]]$samp <- ">diploid"

		#select by
		chr.samp.segs.df <- segs.to.make.MEDICC.input.df[segs.to.make.MEDICC.input.df$sample == sample.vec[1]
													& segs.to.make.MEDICC.input.df[,chr.col] == chr, , drop = FALSE]
		#print(head(chr.samp.segs.df))
		MEDICC.FASTA.file.list[[chr]][["diploid"]]$cpnA <- paste(rep("1", nrow(chr.samp.segs.df)), collapse="")
		MEDICC.FASTA.file.list[[chr]][["diploid"]]$cpnB <- paste(rep("1", nrow(chr.samp.segs.df)), collapse="")


		for (samp in sample.vec) {

			chr.samp.segs.df <- segs.to.make.MEDICC.input.df[segs.to.make.MEDICC.input.df$sample == samp
														& segs.to.make.MEDICC.input.df[, chr.col] == chr, , drop = FALSE]

			#So now let's format the copy number so it will fit MEDICC properly
			MEDICC.FASTA.file.list[[chr]][[samp]]$samp <- paste0(">",samp)
			chr.samp.cpnA.vec <- chr.samp.segs.df[, cpnA.col]
			string.chr.samp.cpnA.vec <- as.character(chr.samp.cpnA.vec)
			string.chr.samp.cpnA.vec[which(chr.samp.cpnA.vec >= 10)] <- "a"
			chr.samp.cpnB.vec <- chr.samp.segs.df[, cpnB.col]
			string.chr.samp.cpnB.vec <- as.character(chr.samp.cpnB.vec)
			string.chr.samp.cpnB.vec[which(chr.samp.cpnB.vec >= 10)] <- "a"
			#So now let's select the data frame
			MEDICC.FASTA.file.list[[chr]][[samp]]$chrom_vec <- rep(chr, nrow(chr.samp.segs.df))
			MEDICC.FASTA.file.list[[chr]][[samp]]$startpos_vec <- chr.samp.segs.df[,  start.col]
			MEDICC.FASTA.file.list[[chr]][[samp]]$endpos_vec <- chr.samp.segs.df[, end.col]
			MEDICC.FASTA.file.list[[chr]][[samp]]$cpnA_vec <- chr.samp.segs.df[, cpnA.col]
			MEDICC.FASTA.file.list[[chr]][[samp]]$cpnB_vec <- chr.samp.segs.df[, cpnB.col]
			MEDICC.FASTA.file.list[[chr]][[samp]]$cpnA_string <- paste(string.chr.samp.cpnA.vec, collapse = "")
			MEDICC.FASTA.file.list[[chr]][[samp]]$cpnB_string <-  paste(string.chr.samp.cpnB.vec, collapse = "")

		}#end for (samp in sample.vec) {

	}#for (chr in CHRS_TO_USE_VEC) {

	return(MEDICC.FASTA.file.list)

}#End makeMEDICCInputListsFromSegmentationDF


# Create input files in the exact format that MEDICC expects them
# Includes:
# _desc.txt file:
# chr_1 CB1001_chr1_A.fasta CB1001_chr1_B.fasta
# chr_2 CB1001_chr2_A.fasta CB1001_chr2_B.fasta
# chr_3 CB1001_chr3_A.fasta CB1001_chr3_B.fasta
# chr_4 CB1001_chr4_A.fasta CB1001_chr4_B.fasta
# chr_5 CB1001_chr5_A.fasta CB1001_chr5_B.fasta
# chr_6 CB1001_chr6_A.fasta CB1001_chr6_B.fasta
#
# and a .fasta file for each haplotype: (e.g. CB1001_chr1_A.fasta)
# >diploid
# 1111
# >CB1001-11-D01-X1
# 2111
# >CB1001-16-D01-X1
# 2211
# >CB1001-15-D01-X1
# 2211
# >CB1001-17-D01-X1
# 2211
# >CB1001-18-D01-X1
# 2211
createMEDICCInputFiles <- function(medicc_fasta, samples, output_dir) {
  message("Saving MEDICC input files.")

  if (!dir.exists(output_dir))
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # Create the desc.txt file
  desc_file <- file.path(output_dir, "desc.txt")
  for (chromosome in medicc_fasta) {
    chrom <- paste0("chr", chromosome$chrom)
    fasta_file_a <- file.path(output_dir, paste0(chrom, "_A.fasta"))
    fasta_file_b <- file.path(output_dir, paste0(chrom, "_B.fasta"))
    # We don't need full paths to the input .fasta files, they are assumed to
    # be in the same dir as the desc.txt file itself
    #cat(chrom, fasta_file_a, fasta_file_b, "\n", file=desc_file, append=TRUE)
    cat(chrom, paste0(chrom, "_A.fasta"), paste0(chrom, "_B.fasta"), "\n", file=desc_file, append=TRUE)

    cat(chromosome[["diploid"]]$samp, "\n", file=fasta_file_a, append = TRUE)
    cat(chromosome[["diploid"]]$cpnA, "\n", file=fasta_file_a, append = TRUE)
    cat(chromosome[["diploid"]]$samp, "\n", file=fasta_file_b, append = TRUE)
    cat(chromosome[["diploid"]]$cpnB, "\n", file=fasta_file_b, append = TRUE)

    for (sample in samples) {
      cat(chromosome[[sample]]$samp, "\n", file=fasta_file_a, append = TRUE)
      cat(chromosome[[sample]]$cpnA_string, "\n", file=fasta_file_a, append = TRUE)
      cat(chromosome[[sample]]$samp, "\n", file=fasta_file_b, append = TRUE)
      cat(chromosome[[sample]]$cpnB_string, "\n", file=fasta_file_b, append = TRUE)
    }
  }
}



#Calculates the WGII, GII, FLOH and FLOH scores based on the integer copy number calls fed into it, must take a dataframe of with column names:
# c("sample","chrom", "startpos", "endpos", "n.het", "cnTotal", "nMajor","nMinor", "Ploidy", "ACF", "nAraw", "nBraw")
#But Ploidy and ACF may be blank as not acutally used by the the funciton.
calculateWGII <- function(seg.df, threshold = 0.6, include.sex.chrom = FALSE) {

  colnames(seg.df) <- c("sample","chrom", "startpos", "endpos", "n.het", "cnTotal", "nMajor","nMinor", "Ploidy", "ACF", "nAraw", "nBraw")
  cat("calculateWGII: segmentation data frame:\n")
  #str(seg.df)
  #First column is meant to be the sample name
  seg.df[, 1] <- as.character(seg.df[, 1])
  sample.vec <- unique(seg.df[, 1])
  num.samples <- length(sample.vec)
  output.mat <- matrix(NA, nrow = num.samples, ncol = 4)
  colnames(output.mat) <- c("GII", "wGII", "FLOH", "wFLOH")
  rownames(output.mat) <- sample.vec
  chrom.length.vec <- setNames(rep(NA, length(unique(seg.df[,2]))), unique(seg.df[,2]))

  if( ! all(seg.df[,"nMinor"] <= seg.df[,"nMajor"]) ) {

    cat("Warning!! nB  not always <= nA!!  -- Correcting for internal use (only!)\n") # In case ASCAT people change the algorithm
    seg.df <- transform(seg.df, nMinor = ifelse( nMinor > nMajor , nMajor, nMinor),
                                nMajor = ifelse( nMinor > nMajor , nMinor, nMajor))

  }#end if( ! all(seg.df[,"nMinor"] <= seg.df[,"nMajor"]) ) {
  cat("Setting sample ploidy according to major proportion copy number. Ploidy < 2 re-set at 2\n")
  ploidy.mat <- calculatePloidy(seg.df)
  #print("ploidy.mat")
  #print(ploidy.mat)
  ploidy.vec <- setNames(ploidy.mat[,2], rownames(ploidy.mat))
  ploidy.vec[ploidy.vec < 2] <- 2

  chrom.vec <- unique(seg.df[, 2])
  for(chrom in chrom.vec) {

    chrom.length.vec[chrom] <- max(as.numeric(seg.df[seg.df[, 2] %in% chrom, 4])) - min(as.numeric(seg.df[seg.df[,2] %in% chrom, 3]))


  }#end for(i in names(chrom.length)){ # restrict to the part that can be measured
  chrom.length.vec <- setNames(as.numeric(chrom.length.vec), names(chrom.length.vec))

  for (samp in sample.vec) {

    samp.seg.df <- seg.df[seg.df[, 1] %in% samp,]
    #cat("Segmentation for:", samp, "\n")
    #print(str(samp.seg.df))

    samp.seg.SCNA.df <- samp.seg.df[samp.seg.df[, 6] < (ploidy.vec[samp] - threshold) | samp.seg.df[, 6] > (ploidy.vec[samp] + threshold),] # Restrict to aberrant segments
    samp.seg.LOH.df <- samp.seg.df[samp.seg.df[, 8] == 0 & samp.seg.df[, 7] != 0,] # Restrict to LOH segments

    output.mat[samp, 1] <- sum(as.numeric(samp.seg.SCNA.df[, 4]) - as.numeric(samp.seg.SCNA.df[, 3])) / sum(as.numeric(chrom.length.vec)) # GII
    output.mat[samp, 3] <- sum(as.numeric(samp.seg.LOH.df[, 4] - as.numeric(samp.seg.LOH.df[, 3])) / sum(as.numeric(chrom.length.vec))) # FLOH

    wgii.vec <- vector()
    wfloh.vec <- vector()

    if ( ! include.sex.chrom == TRUE ) {

       chrom.length.vec <- chrom.length.vec[ ! names(chrom.length.vec) %in% c("X", "Y", "x", "y", 23, 24)]

    }#end if ( ! include.sex.chrom == TRUE ) {
    #print("chrom.vec")
    #print(chrom.vec)
    for (chrom in chrom.vec) {

      chrom.samp.seg.SCNA.df <- samp.seg.SCNA.df[samp.seg.SCNA.df[, 2] %in% chrom,]
      chrom.samp.seg.LOH.df <- samp.seg.LOH.df[samp.seg.LOH.df[, 2] %in% chrom,]

      #cat("Chrom df:", samp, chrom, "\n")
      #cat("SCNA\n")
      #print(str(chrom.samp.seg.SCNA.df))
      #cat("LOH\n")
      #print(str(chrom.samp.seg.LOH.df))

      wgii.vec <- c(wgii.vec, sum(chrom.samp.seg.SCNA.df[, 4] - chrom.samp.seg.SCNA.df[, 3]) / chrom.length.vec[chrom])
      wfloh.vec <- c(wfloh.vec, sum(chrom.samp.seg.LOH.df[,4] - chrom.samp.seg.LOH.df[,3]) / chrom.length.vec[chrom])

      #print("==== wgii.vec")
      #print(wgii.vec)
      #print("==== wfloh.vec")
      #print(wfloh.vec)

    }#end for (chrom in chrom.vec) {
    output.mat[samp, 2] <- sum(wgii.vec) / length(chrom.length.vec)
    output.mat[samp, 4] <- sum(wfloh.vec) / length(chrom.length.vec)

  }#end for (samp in sample.vec)

  #cat("WGII output matrix\n")
  #str(output.mat)
  #print(output.mat)

  return(output.mat)

}#calculateWGII


calculatePloidy  <- function(seg.df, tot.cpn.col.idx = 6, calculate.raw.ploidy = TRUE) {

  seg.df[, 1]  <- as.character(seg.df[, 1])
  colnames(seg.df) <- c("sample","chrom", "startpos", "endpos", "n.het", "cnTotal", "nMajor","nMinor", "Ploidy", "ACF", "nAraw", "nBraw")
  sample.vec <- unique(seg.df$sample)
  num.samples <- length(sample.vec)
  out.ploidy.mat <- matrix(nrow = num.samples, ncol = 14)
  row.names(out.ploidy.mat) <- sample.vec
  colnames(out.ploidy.mat) <- c("Major", "wMajor", "MeanRaw", 0:10)  #wMajor is weighted by chromosome, MeanRaw is mean of the raw CN values from ASCAT 2.3

  for (samp in sample.vec) {

    #Integer ploidy
    samp.seg.df <- seg.df[seg.df$sample %in% samp,]
    bp.size.by.tot.cpn.vec  <- c()
    uniq.tot.cpn.val.vec <- unique(samp.seg.df[, tot.cpn.col.idx])

    for (tot.cpn.val in uniq.tot.cpn.val.vec) {

      samp.seg.w.cpn.val.df <- samp.seg.df[samp.seg.df[, tot.cpn.col.idx] %in% tot.cpn.val,]
      bp.size.by.tot.cpn.vec  <- c(bp.size.by.tot.cpn.vec , setNames(sum(as.numeric(samp.seg.w.cpn.val.df$endpos - samp.seg.w.cpn.val.df$startpos)), tot.cpn.val))

    }#end for (tot.cpn.state in uniq.tot.cpn.val.vec) {
    summ.tot.cpn.by.bp.size <- (bp.size.by.tot.cpn.vec / sum(as.numeric(bp.size.by.tot.cpn.vec)))
    sort.tot.cpn.by.bp.size <- sort(summ.tot.cpn.by.bp.size, decreasing = TRUE)
    out.ploidy.mat[samp, "Major"] <- as.numeric(names(sort.tot.cpn.by.bp.size)[1])
    out.ploidy.mat[samp, 4:13] <- sort.tot.cpn.by.bp.size[colnames(out.ploidy.mat)[4:13]]

    #Weighted ploidy:
    w.major.vec <- c()
    chrom.vec <- unique(seg.df[, 2])
    for (chrom in chrom.vec) {

      chrom.samp.seg.df <- samp.seg.df[samp.seg.df[, 2] %in% chrom,]
      #Get the sizes of the base pairs
      chr.tot.cpn.val.bp.per.seg.vec <- setNames(as.numeric(chrom.samp.seg.df[, 4]) - as.numeric(chrom.samp.seg.df[, 3]), chrom.samp.seg.df[, tot.cpn.col.idx])
      #Sum up all the seg sizes by the tot cpn value that they had
      chr.tot.cpn.bp.vec <- sort(tapply(chr.tot.cpn.val.bp.per.seg.vec, names(chr.tot.cpn.val.bp.per.seg.vec),sum), decreasing = TRUE)
      #Take the top value and add it to
      w.major.vec <- c(w.major.vec, setNames(names(chr.tot.cpn.bp.vec[1]), chrom))

    }#end for (chrom in chrom.vec) {
    out.ploidy.mat[samp , "wMajor"] <- as.numeric(names(sort(table(w.major.vec))))[1]

    #Raw ploidy:

    if (calculate.raw.ploidy == TRUE) {

      samp.tot.raw.cpn.vec <- as.numeric(rowSums(samp.seg.df[, c("nAraw", "nBraw")]))
      samp.seg.size.vec <- as.numeric(samp.seg.df[, 4]) - as.numeric(samp.seg.df[, 3])

      out.ploidy.mat[samp, "MeanRaw"] <- sum(as.numeric(samp.tot.raw.cpn.vec) * (samp.seg.size.vec/sum(samp.seg.size.vec)))

    }#end if (calculate.raw.ploidy == TRUE) {

  }#end for (samp in sample.vec) {

  if (calculate.raw.ploidy == FALSE) {

    out.ploidy.mat <- out.ploidy.mat[, -c(3) ]

  }

  return(out.ploidy.mat)

}#end calculatePloidy  <- function(seg.df, cpn.col.idx = 6) {
