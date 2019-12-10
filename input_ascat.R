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
read_snps_ascat <- function(sample_table, segs, samples) {
  all_snps <- list()

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
  all_snps <- snps_set_bafsegmented(segs, all_snps)
  all_snps <- keep_common_snps(all_snps)

  all_snps
}
