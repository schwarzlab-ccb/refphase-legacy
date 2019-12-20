# -------------------------------------------------------------------------
# Sequenza input format
#
# Read BAF and logR values from seqz files that are generated for sequenza.
# These files are seqz format, a tab-delimited and gzipped file format with the
# following columns:
#
# * chromosome
# * position
# * base.ref
# * depth.normal
# * depth.tumor
# * depth.ratio
# * Af
# * Bf
# * zygosity.normal
# * GC.percent
# * good.reads
# * AB.normal
# * AB.tumor
# * tumor.strand
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

  for (samp in samples) {
    input_file <- input_files[[samp]]

    # All data that is available in a .seqz file. We only need a subset of it
    # but this is a useful reference in case we need more info at a later date
    #
    # col_types = "ciciidddcddccc"
    # col_names = c("chromosome", "position", "base.ref", "depth.normal",
    #               "depth.tumor", "depth.ratio", "Af", "Bf", "zygosity.normal",
    #               "GC.percent", "good.reads", "AB.normal", "AB.tumor",
    #               "tumor.strand")
    # Using readr
    message("Reading SNP data for sample: ", samp)
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

    chrom <- chrom2integer(seqz_data$chromosome)
    pos <- seqz_data$position

    # Bf is always < Af, so we need to check the AB.normal to see if we need to flip it.
    # This is probably not the perfect thing to do. Look into this topic more.
    flip_baf <- nchar(seqz_data$AB.normal) > 0 & seqz_data$base.ref != substr(seqz_data$AB.normal, 1, 1)
    baf <- seqz_data$Bf
    baf[flip_baf] <- seqz_data$Af[flip_baf]

    warning("Hack to adjust for BAF reference bias (mean BAF is ~ 0.475)")
    baf <- baf + (0.5 - 0.475)
    baf[baf > 1.0] <- 1.0

    all_snps[[samp]] <- data.frame(
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
    all_snps[[samp]]$BAF[all_snps[[samp]]$germline_zygosity != "het"] <- NA

    # Look into reference bias
    #library(data.table)
    #baf_tmp <- all_snps[[samp]]$BAF[!is.na(all_snps[[samp]]$BAF)]
    #message("Mean of BAF: ", mean(baf_tmp))
    #message("Median of BAF: ", median(baf_tmp))
    #message("Rolling mean of BAF (window 100 hz positions):")
    #print(frollmean(baf_tmp, 100))
  }

  all_snps <- snps_set_bafsegmented(segs, all_snps)
  all_snps <- snps_set_logrsegmented(segs, all_snps)
  all_snps <- keep_common_snps(all_snps)

  all_snps
}
