annotate_joint_segs <- function(segs, snps, patient.phased.BAF.df.list, conf) {
  logwarn("TODO: check this function (annotate_joint_segs()) to be sure it works properly without germline information.")
  sample.vec <- unique(segs$sample)

  anno.seg.df <- segs

  # So determine which SNPs are hz:
  snp.file.df <- snps

  snp.file.df$snp_pos <- paste0(snp.file.df$chrom, "_", snp.file.df$pos)
  snp.file.df$snp_idx <- seq_len(nrow(snp.file.df))

  # #CB
  # > head(snp.file.df)
  #   chrom    pos BAF logR BAFsegmented logRsegmented  snp_pos snp_idx
  # 1     1  13116  NA    0           NA    -0.8384372  1_13116       1
  # 2     1  13118  NA    0           NA    -0.8384372  1_13118       2
  # 3     1  13273  NA    0           NA    -0.8384372  1_13273       3
  # 4     1  13284  NA    0           NA    -0.8384372  1_13284       4
  # 5     1  13380  NA    0           NA    -0.8384372  1_13380       5
  # 6     1 138781  NA    0           NA    -0.8384372 1_138781       6

  # #LTX121
  # > head(snp.file.df)
  #   chrom    pos BAF logR BAFsegmented logRsegmented  snp_pos
  # 1     1 762273  NA    0           NA    -0.3346204 1_762273
  # 2     1 762496  NA   NA           NA    -0.3346204 1_762496
  # 3     1 762497  NA   NA           NA    -0.3346204 1_762497
  # 4     1 777237  NA   NA           NA    -0.3346204 1_777237
  # 5     1 777318  NA   NA           NA    -0.3346204 1_777318
  # 6     1 777335  NA   NA           NA    -0.3346204 1_777335
  # Differences can be seen here and are noted when the gl.snp.df has been read in


  # 	#get germline column
  gl.col.name <- "BAF"
  # 	logr.only.snp.df <- snp.file.df[is.na(snp.file.df[,gl.col.name]),]
  actual.snp.file.df <- snp.file.df[!is.na(snp.file.df[, gl.col.name]), ]
  #
  #  gl.homo.snp.df <- data.frame()
  #  gl.hz.snp.df <- data.frame()
  #  if (conf$input_format == "ASCAT") {
  #    gl.homo.snp.df <- actual.snp.file.df[! is.na(actual.snp.file.df[,gl.col.name]) & (as.numeric(actual.snp.file.df[,gl.col.name]) > 0.7 | as.numeric(actual.snp.file.df[,gl.col.name]) < 0.3) ,]
  #    gl.hz.snp.df <- actual.snp.file.df[! is.na(actual.snp.file.df[,gl.col.name]) & (as.numeric(actual.snp.file.df[,gl.col.name]) <= 0.7 &  as.numeric(actual.snp.file.df[,gl.col.name]) >= 0.3),]
  #  } else {
  gl.homo.snp.df <- actual.snp.file.df[actual.snp.file.df$germline_zygosity != "het", ]
  gl.hz.snp.df <- actual.snp.file.df[actual.snp.file.df$germline_zygosity == "het", ]
  logr.only.snp.df <- snp.file.df[is.na(snp.file.df[, gl.col.name]), ]

  #  }
  #
  # 	gl.hz.snp.df$hz_snp_idx <- seq_len(nrow(gl.hz.snp.df))

  # Note anno seg is sorted.all.bp.all.samps.segs.anno.df
  anno.seg.df$start_pos <- paste0(anno.seg.df$Chromosome, "_", anno.seg.df$Start_Position)
  anno.seg.df$start_is_hz <- anno.seg.df$start_pos %in% gl.hz.snp.df$snp_pos
  anno.seg.df$start_is_homo <- anno.seg.df$start_pos %in% gl.homo.snp.df$snp_pos
  anno.seg.df$start_is_logr <- anno.seg.df$start_pos %in% logr.only.snp.df$snp_pos
  anno.seg.df$end_pos <- paste0(anno.seg.df$Chromosome, "_", anno.seg.df$End_Position)
  anno.seg.df$end_is_hz <- anno.seg.df$end_pos %in% gl.hz.snp.df$snp_pos
  anno.seg.df$end_is_homo <- anno.seg.df$end_pos %in% gl.homo.snp.df$snp_pos
  anno.seg.df$end_is_logr <- anno.seg.df$end_pos %in% logr.only.snp.df$snp_pos
  # make a quick summary column
  anno.seg.df$start_is_pos <- anno.seg.df$start_is_hz | anno.seg.df$start_is_homo | anno.seg.df$start_is_logr
  anno.seg.df$end_is_pos <- anno.seg.df$end_is_hz | anno.seg.df$end_is_homo | anno.seg.df$end_is_logr
  anno.seg.df$start_snp_idx <- NA
  anno.seg.df$end_snp_idx <- NA
  anno.seg.df[anno.seg.df$start_is_pos, ]$start_snp_idx <- snp.file.df[match(anno.seg.df[anno.seg.df$start_is_pos, ]$start_pos, snp.file.df$snp_pos), ]$snp_idx
  anno.seg.df[anno.seg.df$end_is_pos, ]$end_snp_idx <- snp.file.df[match(anno.seg.df[anno.seg.df$end_is_pos, ]$end_pos, snp.file.df$snp_pos), ]$snp_idx
  # potential column names of use: in_intersection in_disjoin_only in_arm_setdiff
  # #For no now let's only take those those segments that are in the intersection
  # int.anno.seg.df <- anno.seg.df[anno.seg.df$in_intersection%in%TRUE & anno.seg.df$full_sample %in% sample.vec[1],]

  # Let's try adding in the disjoin...
  int.anno.seg.df <- anno.seg.df[anno.seg.df$in_intersection | anno.seg.df$in_disjoin_only & anno.seg.df$full_sample %in% sample.vec[1], ]

  # > head(int.anno.seg.df[,ncol(int.anno.seg.df):(ncol(int.anno.seg.df)-12)])
  #                                      end_snp_idx start_snp_idx end_is_pos start_is_pos end_is_logr end_is_homo end_is_hz    end_pos start_is_logr start_is_homo start_is_hz  start_pos chr_arm
  # CB1001-11-D01-X1:1:138781:16421937          7928             6       TRUE         TRUE       FALSE       FALSE      TRUE 1_16421937          TRUE         FALSE       FALSE   1_138781      1p
  # CB1001-11-D01-X1:1:16421938:16490693          NA            NA      FALSE        FALSE       FALSE       FALSE     FALSE 1_16490693         FALSE         FALSE       FALSE 1_16421938      1p
  # CB1001-11-D01-X1:1:16490694:16535960          NA          7973      FALSE         TRUE       FALSE       FALSE     FALSE 1_16535960         FALSE         FALSE        TRUE 1_16490694      1p
  # CB1001-11-D01-X1:1:16535961:16948870        8287          7979       TRUE         TRUE        TRUE       FALSE     FALSE 1_16948870          TRUE         FALSE       FALSE 1_16535961      1p
  # CB1001-11-D01-X1:1:16948871:17232700          NA            NA      FALSE        FALSE       FALSE       FALSE     FALSE 1_17232700         FALSE         FALSE       FALSE 1_16948871      1p
  # CB1001-11-D01-X1:1:17232701:17243971        8487          8477       TRUE         TRUE       FALSE       FALSE      TRUE 1_17243971         FALSE         FALSE        TRUE 1_17232701      1p

  # > head(int.anno.seg.df[,ncol(int.anno.seg.df):(ncol(int.anno.seg.df)-12)])
  #                                         end_snp_idx start_snp_idx end_is_pos start_is_pos end_is_logr end_is_homo end_is_hz     end_pos start_is_logr start_is_homo start_is_hz   start_pos chr_arm
  # A_LTX121_SU_T1-R1:1:762273:26136282           17836             1       TRUE         TRUE        TRUE       FALSE     FALSE  1_26136282          TRUE         FALSE       FALSE    1_762273      1p
  # A_LTX121_SU_T1-R1:1:26136283:26136297            NA            NA      FALSE        FALSE       FALSE       FALSE     FALSE  1_26136297         FALSE         FALSE       FALSE  1_26136283      1p
  # A_LTX121_SU_T1-R1:1:26136298:121310061        50594         17837       TRUE         TRUE        TRUE       FALSE     FALSE 1_121310061          TRUE         FALSE       FALSE  1_26136298      1p
  # A_LTX121_SU_T1-R1:1:144052752:145096660       50645         50595       TRUE         TRUE        TRUE       FALSE     FALSE 1_145096660          TRUE         FALSE       FALSE 1_144052752      1q
  # A_LTX121_SU_T1-R1:1:145096661:145109582          NA            NA      FALSE        FALSE       FALSE       FALSE     FALSE 1_145109582         FALSE         FALSE       FALSE 1_145096661      1q
  # A_LTX121_SU_T1-R1:1:145109583:145477446       50734         50646       TRUE         TRUE        TRUE       FALSE     FALSE 1_145477446         FALSE         FALSE        TRUE 1_145109583      1q

  if (conf$genome == "hg19") {
    info.path <- system.file("extdata", "ucsc_hg19_chrominfo.txt", package = "refphase")
  } else if (conf$genome == "hg38") {
    info.path <- system.file("extdata", "ucsc_hg36_chrominfo.txt", package = "refphase")
  } else {
    stop("Genome version not recognized: ", conf$genome)
  }
  chrom.info.df <- read.table(info.path, stringsAsFactors = FALSE, sep = "\t", header = TRUE)

  # Let's just add the arm classificaitons to all the snp positions
  snp.file.df$arm <- NA
  for (ch in unique(int.anno.seg.df$chrom)) {
    # get start and end of p
    p.start <- chrom.info.df[chrom.info.df$chrom %in% ch, ]$p_start
    p.end <- chrom.info.df[chrom.info.df$chrom %in% ch, ]$p_end
    q.start <- chrom.info.df[chrom.info.df$chrom %in% ch, ]$q_start
    q.end <- chrom.info.df[chrom.info.df$chrom %in% ch, ]$q_end

    if (nrow(snp.file.df[snp.file.df$chrom %in% ch & snp.file.df$pos >= p.start & snp.file.df$pos <= p.end, ]) > 0) {
      snp.file.df[snp.file.df$chrom %in% ch & snp.file.df$pos >= p.start & snp.file.df$pos <= p.end, ]$arm <- paste0(ch, "p")
    }

    if (nrow(snp.file.df[snp.file.df$chrom %in% ch & snp.file.df$pos >= q.start & snp.file.df$pos <= q.end, ]) > 0) {
      snp.file.df[snp.file.df$chrom %in% ch & snp.file.df$pos >= q.start & snp.file.df$pos <= q.end, ]$arm <- paste0(ch, "q")
    }
  } # end for (chr in unique(int.anno.seg.df$chrom)) {

  # #IN CB there are some SNPs that do not have an arm assiged:
  # >  length(unique(snp.file.df$arm))
  # [1] 41
  # > unique(snp.file.df$arm
  # + )
  #  [1] "1p"  "1q"  "2p"  "2q"  "3p"  "3q"  "4p"  "4q"  "5p"  "5q"
  #  "6p"  "6q"  "7p"  "7q"  "8p"  "8q"  "9p"  "9q"  "10p" "10q" "11p"
  #  "11q" "12p" "12q" "13q" "14q" "15q" "16p" "16q" "17p" "17q" "18p"
  #  "18q" "19p" "19q" "20p" "20q" "21p" "21q" "22q" NA

  # vs LTX121

  # > length(unique(snp.file.df$arm))
  # [1] 40
  # > unique(snp.file.df$arm)
  #  [1] "1p"  "1q"  "2p"  "2q"  "3p"  "3q"  "4p"  "4q"  "5p"  "5q"
  #    "6p"  "6q"  "7p"  "7q"  "8p"  "8q"  "9p"  "9q"  "10p" "10q" "11p"
  #     "11q" "12p" "12q" "13q" "14q" "15q" "16p" "16q" "17p" "17q" "18p"
  #      "18q" "19p" "19q" "20p" "20q" "21p" "21q" "22q"

  # So we need to write a function that finds the nearest SNP that is
  # 1) in the same chromosome arm
  # 2) not already the the start or end of another segment
  # 3) Not within another segment...
  # Remove any of the SNPs that are already starting or ending segments...
  all.snp.grs <- makeGRangesFromDataFrame(snp.file.df[ !snp.file.df$snp_pos %in% int.anno.seg.df$start_pos & !snp.file.df$snp_pos %in% int.anno.seg.df$end_pos, ],
    keep.extra.columns = FALSE,
    ignore.strand = FALSE,
    seqinfo = NULL,
    seqnames.field = c("chrom"),
    start.field = c("pos"),
    end.field = c("pos"),
    strand.field = "strand",
    starts.in.df.are.0based = FALSE
  )
  # CB
  # > head(all.snp.grs)
  # GRanges object with 6 ranges and 0 metadata columns:
  #     seqnames    ranges strand
  #        <Rle> <IRanges>  <Rle>
  #   1        1     13116      *
  #   2        1     13118      *
  #   3        1     13273      *
  #   4        1     13284      *
  #   5        1     13380      *
  #   7        1    138817      *
  #   -------
  # LTX121
  # > head(all.snp.grs)
  # GRanges object with 6 ranges and 0 metadata columns:
  #     seqnames    ranges strand
  #        <Rle> <IRanges>  <Rle>
  #   2        1    762496      *
  #   3        1    762497      *
  #   4        1    777237      *
  #   5        1    777318      *
  #   6        1    777335      *
  #   7        1    777361      *

  # Get indices of SNPs with missing start or end idxs
  int.anno.seg.df$strand <- "*"
  all.int.grs <- makeGRangesFromDataFrame(int.anno.seg.df,
    keep.extra.columns = FALSE,
    ignore.strand = FALSE,
    seqinfo = NULL,
    seqnames.field = c("Chromosome"),
    start.field = c("Start_Position"),
    end.field = c("End_Position"),
    strand.field = "strand",
    starts.in.df.are.0based = FALSE
  )

  # > head(all.int.grs)
  # GRanges object with 6 ranges and 0 metadata columns:
  #                                        seqnames            ranges strand
  #                                           <Rle>         <IRanges>  <Rle>
  #     CB1001-11-D01-X1:1:138781:16421937        1   138781-16421937      *
  #   CB1001-11-D01-X1:1:16421938:16490693        1 16421938-16490693      *
  #   CB1001-11-D01-X1:1:16490694:16535960        1 16490694-16535960      *
  #   CB1001-11-D01-X1:1:16535961:16948870        1 16535961-16948870      *
  #   CB1001-11-D01-X1:1:16948871:17232700        1 16948871-17232700      *
  #   CB1001-11-D01-X1:1:17232701:17243971        1 17232701-17243971      *

  # > head(all.int.grs)
  # GRanges object with 6 ranges and 0 metadata columns:
  #                                           seqnames              ranges strand
  #                                              <Rle>           <IRanges>  <Rle>
  #       A_LTX121_SU_T1-R1:1:762273:26136282        1     762273-26136282      *
  #     A_LTX121_SU_T1-R1:1:26136283:26136297        1   26136283-26136297      *
  #    A_LTX121_SU_T1-R1:1:26136298:121310061        1  26136298-121310061      *
  #   A_LTX121_SU_T1-R1:1:144052752:145096660        1 144052752-145096660      *
  #   A_LTX121_SU_T1-R1:1:145096661:145109582        1 145096661-145109582      *
  #   A_LTX121_SU_T1-R1:1:145109583:145477446        1 145109583-145477446      *
  #

  # Calculate all of these in one step, much faster than doing it one at a time
  before <- follow(all.int.grs, all.snp.grs)
  after <- precede(all.int.grs, all.snp.grs)
  nearby <- nearest(all.int.grs, all.snp.grs)

  logdebug("Setting nearest preceding/following SNP")
  lacking.snp.idx.vec <- which(!(int.anno.seg.df$start_is_pos & int.anno.seg.df$end_is_pos))
  int.anno.seg.df$nearest_following_snp <- NA
  int.anno.seg.df$nearest_preceding_snp <- NA
  for (i in seq_along(all.int.grs)) {
    # Check the seqnames part as SNPs really need to be in the same need to be in the same chromosome arm...?
    if (!i %in% lacking.snp.idx.vec) {
      next
    }
    #logdebug(paste0(" ", i))

    curr.arm <- int.anno.seg.df[i, ]$chr_arm
    i_grs <- all.int.grs[i]
    # Precede is wrong here follow is correct, check the granges documentation.
    before_i <- before[[i]]
    if (!is.na(before_i)) {
      if (snp.file.df[before_i, ]$arm == curr.arm) {
        int.anno.seg.df[i, ]$nearest_preceding_snp <- start(all.snp.grs[before_i])
      }
    }

    # follow is wrong here precede is correct, check the granges documentation.
    after_i <- after[[i]]
    if (!is.na(after_i)) {
      if (snp.file.df[after_i, ]$arm == curr.arm) {
        int.anno.seg.df[i, ]$nearest_following_snp <- start(all.snp.grs[after_i])
      }
    }

    # Any remaining examples are going to be the first or only
    if (is.na(int.anno.seg.df[i, ]$nearest_following_snp) & is.na(int.anno.seg.df[i, ]$nearest_preceding_snp)) {
      near_snp_start <- start(all.snp.grs[nearby[[i]]])
      int.anno.seg.df[i, ]$nearest_following_snp <- near_snp_start
      int.anno.seg.df[i, ]$nearest_preceding_snp <- near_snp_start
    }
  }
  logdebug("done")

  int.anno.seg.df[is.na(int.anno.seg.df$start_snp_idx) | is.na(int.anno.seg.df$end_snp_idx), c("Start_Position", "End_Position", "start_snp_idx", "end_snp_idx", "nearest_preceding_snp", "nearest_following_snp") ]
  int.anno.seg.df[int.anno.seg.df$start_is_pos %in% FALSE, ]$Start_Position <- int.anno.seg.df[int.anno.seg.df$start_is_pos %in% FALSE, ]$nearest_following_snp
  int.anno.seg.df[is.na(int.anno.seg.df$Start_Position), ]$Start_Position <- int.anno.seg.df[is.na(int.anno.seg.df$Start_Position), ]$nearest_preceding_snp
  int.anno.seg.df[int.anno.seg.df$end_is_pos %in% FALSE, ]$End_Position <- int.anno.seg.df[int.anno.seg.df$end_is_pos %in% FALSE, ]$nearest_preceding_snp
  int.anno.seg.df[is.na(int.anno.seg.df$End_Position), ]$End_Position <- int.anno.seg.df[is.na(int.anno.seg.df$End_Position), ]$nearest_following_snp

  f.int.anno.seg.df <- int.anno.seg.df
  for (i in 2:nrow(f.int.anno.seg.df)) {
    if (f.int.anno.seg.df[ (i - 1), ]$chr_arm == f.int.anno.seg.df[i, ]$chr_arm) {
      if (f.int.anno.seg.df[ (i - 1), ]$End_Position > f.int.anno.seg.df[i, ]$Start_Position) {
        preceding.end.pos <- f.int.anno.seg.df[ (i - 1), ]$End_Position
        f.int.anno.seg.df[ (i - 1), ]$End_Position <- f.int.anno.seg.df[i, ]$Start_Position
        f.int.anno.seg.df[i, ]$Start_Position <- preceding.end.pos
      } # end if(f.int.anno.seg.df[ (i - 1),]$End_Position > f.int.anno.seg.df[i, ]$Start_Position) {
    } # end if(f.int.anno.seg.df[ (i - 1),]$chr_arm == f.int.anno.seg.df[i,]$chr_arm) {
  } # end for (i in 2:nrow(f.int.anno.seg.df)) {

  f.int.anno.seg.df$start_pos <- paste0(f.int.anno.seg.df$Chromosome, "_", f.int.anno.seg.df$Start_Position)
  f.int.anno.seg.df$start_is_hz <- f.int.anno.seg.df$start_pos %in% gl.hz.snp.df$snp_pos
  f.int.anno.seg.df$start_is_homo <- f.int.anno.seg.df$start_pos %in% gl.homo.snp.df$snp_pos
  f.int.anno.seg.df$start_is_logr <- f.int.anno.seg.df$start_pos %in% logr.only.snp.df$snp_pos
  f.int.anno.seg.df$end_pos <- paste0(f.int.anno.seg.df$Chromosome, "_", f.int.anno.seg.df$End_Position)
  f.int.anno.seg.df$end_is_hz <- f.int.anno.seg.df$end_pos %in% gl.hz.snp.df$snp_pos
  f.int.anno.seg.df$end_is_homo <- f.int.anno.seg.df$end_pos %in% gl.homo.snp.df$snp_pos
  f.int.anno.seg.df$end_is_logr <- f.int.anno.seg.df$end_pos %in% logr.only.snp.df$snp_pos

  # make a quick summary column
  f.int.anno.seg.df$start_is_pos <- f.int.anno.seg.df$start_is_hz | f.int.anno.seg.df$start_is_homo | f.int.anno.seg.df$start_is_logr
  f.int.anno.seg.df$end_is_pos <- f.int.anno.seg.df$end_is_hz | f.int.anno.seg.df$end_is_homo | f.int.anno.seg.df$end_is_logr
  # So now let's assign these new snps to be the starts and ends...
  f.int.anno.seg.df[f.int.anno.seg.df$start_is_pos %in% TRUE, ]$start_snp_idx <- snp.file.df[match(f.int.anno.seg.df[f.int.anno.seg.df$start_is_pos %in% TRUE, ]$start_pos, snp.file.df$snp_pos), ]$snp_idx
  f.int.anno.seg.df[f.int.anno.seg.df$end_is_pos %in% TRUE, ]$end_snp_idx <- snp.file.df[match(f.int.anno.seg.df[f.int.anno.seg.df$end_is_pos %in% TRUE, ]$end_pos, snp.file.df$snp_pos), ]$snp_idx

  # So we're going to need the indices of the BAF snps, only needs to be done once as will be the same for all samples
  phased.baf.info.df <- patient.phased.BAF.df.list[[1]]
  phased.baf.info.df$snp_pos <- paste0(phased.baf.info.df$chrom, "_", phased.baf.info.df$pos)
  # Now we need to know their idx in the original data frame
  phased.baf.info.df$snp_idx <- snp.file.df[match(phased.baf.info.df$snp_pos, snp.file.df$snp_pos), ]$snp_idx
  # The variable was renamed here to allow debugging but is just continues.
  all.samp.j.seg.df <- f.int.anno.seg.df

  return(list(
    all.samp.j.seg.df = all.samp.j.seg.df,
    phased.baf.info.df = phased.baf.info.df,
    gl.hz.snp.df = gl.hz.snp.df
  ))
}
