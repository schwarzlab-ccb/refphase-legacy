# Read input data: segmentation and snps with associated BAF/logR and more
read_input <- function(seg_file, snp_file, input_format) {
  init.f.joint.segs.df <- read.table(seg.file.path, stringsAsFactors = FALSE, sep = "\t", header = TRUE)

  # This function will stop the script if the input is not valid
  if (conf$validate_input) {
    validate_seg_all_samples(init.f.joint.segs.df)
  }

  f.joint.segs.df <- reformatInitialSegmentationInput(f.joint.segs.df = init.f.joint.segs.df)

  sample.vec <- unique(f.joint.segs.df$sample)
  ploidy.est.vec <- f.joint.segs.df[ !duplicated(f.joint.segs.df$sample), "Ploidy"]
  names(ploidy.est.vec) <- unique(f.joint.segs.df$sample)
  purity.est.vec <- f.joint.segs.df[ !duplicated(f.joint.segs.df$sample), "ACF"]
  names(purity.est.vec) <- unique(f.joint.segs.df$sample)

  if (input_format == "ASCAT") {
    # Parse a bunch of rds files containing ascat.bc objects
    pt.snp.file.paths.df <- read.delim(snp.file.path,
      stringsAsFactors = FALSE,
      header = TRUE,
      row.names = "sample"
    )

    if (nchar(conf$snp_data_prefix) > 0) {
      pt.snp.file.paths.df$samp_ascat_bc_path <- file.path(conf$snp_data_prefix, pt.snp.file.paths.df$samp_ascat_bc_path)
    }

    all.samp.SNP.list <- read_snps_ascat(
      sample_table = pt.snp.file.paths.df,
      segs = f.joint.segs.df,
      samples = sample.vec
    )

    # TBKW edit, I'm assuming that these ascat objects will also contain the germline values and am pulling
    # them out for this "xphase" version that requires the germline BAFs to compare the allelic imbalance in each tumour
    # region segment vs the equivalent segment in the germline. This comparison is done to try and counter the effect
    # of noise and "reference bias" that might result in false allelic imbalance. The correction that is performed to avert this
    # is primitive as of 16/08/2019 and will require rewriting as well as testing. It should be noted the presence of
    # bias observed by allelic imbalance is from separate work I've performed and eye-balling the data, it needs to be
    # shown in the data...

    # Making the assumption here that all regions that were analysed had
    # gl.snp.df <- readInGermlineFromASCAT(ascat_bc_rds_path = pt.snp.file.paths.df[1,]$samp_ascat_bc_path,
    # 									snp_pos_vec = paste0(all.samp.SNP.list[[1]]$chrom
    # 														,"_"
    # 														,all.samp.SNP.list[[1]]$pos))
  } else if (input_format == "sequenza") {
    pt.snp.file.paths.df <- read.delim(snp.file.path,
      stringsAsFactors = FALSE,
      header = TRUE,
      row.names = "sample"
    )

    message("Start parsing sequenza input.\n")

    input_files <- pt.snp.file.paths.df$seqz_file
    if (nchar(conf$snp_data_prefix) > 0) {
      input_files <- file.path(conf$snp_data_prefix, pt.snp.file.paths.df$seqz_file)
    }
    names(input_files) <- rownames(pt.snp.file.paths.df)

    # Read in SNPS, keeping homozygous positions
    all.samp.SNP.list <- read_snps_sequenza(
      input_files = input_files,
      segs = f.joint.segs.df,
      samples = sample.vec,
      het_only = FALSE
    )

    message("Finished parsing sequenza input.\n")
  } else {
    stop(paste0("Incorrect snp data input format (", input_format, ")"))
  }

  # Immediately validate input SNP information (BAF, logR and more)
  if (conf$validate_input) {
    validate_snp_all_samples(all.samp.SNP.list)
  }

  # So we need the centromere postions for each arm:
  # Add in p_end and q_start

  # TBKW  20190818 edit so here is where the joint segmentation stuff has been dumped into a function to allow multiple
  # calls. Absolutely needs refactoring.
  sorted.all.bp.all.samps.segs.anno.df <- createFullyAnnotatedJointSegmentation(f.joint.segs.df = f.joint.segs.df, samples = sample.vec)
  # TBKW edit for xphase version
  anno.seg.df <- sorted.all.bp.all.samps.segs.anno.df
  # sorted.all.bp.all.samps.segs.anno.df.path <- paste0(PLOT_OUTPUT_DIR, patient, "_sorted_all_bp_all_samps_segs_anno_df.RDS")
  # saveRDS(sorted.all.bp.all.samps.segs.anno.df, file = sorted.all.bp.all.samps.segs.anno.df.path)
  # Now we want to save this annotated version:
  cpn.chr.start.end.df <- sorted.all.bp.all.samps.segs.anno.df[sorted.all.bp.all.samps.segs.anno.df$full_sample == sorted.all.bp.all.samps.segs.anno.df[1, ]$full_sample, 1:3]
  mock.cpn.matrix <- matrix(data = 0, nrow = nrow(cpn.chr.start.end.df), ncol = length(sample.vec))
  colnames(mock.cpn.matrix) <- sample.vec
  min.cons.segs.df <- cbind(cpn.chr.start.end.df, mock.cpn.matrix)

  # Minimum consistent segmentation is finished! Let's get the SNPs and their LogR and BAF values amd assign to segments.
  # All SNPs in these lists must have both a logr and baf value: all sample all seg logr baf lists is asas.lr.baf
  message("cleanAndErrorCheckSNPsParallel first")
  asas.SNPs.w.both.lr.baf.snp.lists <- cleanAndErrorCheckSNPsParallel(all.samp.SNP.list,
    mcs.seg.to.use.df = min.cons.segs.df,
    samples = sample.vec,
    allow.logr.only.flag = FALSE,
    data.type = input_format
  )

  if (length(unique(sapply(asas.SNPs.w.both.lr.baf.snp.lists, length))) != 1) {
    stop("Not all segments are present in all samples' snp assignations, likely something has gone wrong, is the parallelisation working correctly?")
  }

  # SNPs that only have LogR are permitted in these lists, but it also contains all the SNPs with both BAF and LogR
  message("cleanAndErrorCheckSNPsParallel second")
  asas.SNPS.lr.and.or.baf.lists <- cleanAndErrorCheckSNPsParallel(all.samp.SNP.list,
    mcs.seg.to.use.df = min.cons.segs.df,
    samples = sample.vec,
    allow.logr.only.flag = TRUE,
    data.type = input_format
  )

  if (length(unique(sapply(asas.SNPS.lr.and.or.baf.lists, length))) != 1) {
    stop("Not all segments are present in all samples' snp assignations, likely something has gone wrong, is the parallelisation working correctly?")
  }

  if (length(asas.SNPs.w.both.lr.baf.snp.lists) <= 1) {
    print(str(asas.SNPs.w.both.lr.baf.snp.lists, max = 1))
    stop("Only one or fewer samples present after reading in data.")
  }

  return(list(
    all.samp.SNP.list = all.samp.SNP.list,
    asas.SNPS.lr.and.or.baf.lists = asas.SNPS.lr.and.or.baf.lists,
    asas.SNPs.w.both.lr.baf.snp.lists = asas.SNPs.w.both.lr.baf.snp.lists,
    f.joint.segs.df = f.joint.segs.df,
    min.cons.segs.df = min.cons.segs.df,
    ploidy.est.vec = ploidy.est.vec,
    pt.snp.file.paths.df = pt.snp.file.paths.df,
    purity.est.vec = purity.est.vec,
    sorted.all.bp.all.samps.segs.anno.df = sorted.all.bp.all.samps.segs.anno.df,
    anno.seg.df = anno.seg.df
  ))
}


# Read a ascat.bc object that has been saved by an earlier ASCAT run
read_ascat_bc <- function(filename) {
  if (grepl("\\.rdata$", tolower(filename))) {
    load(filename)
  } else if (grepl("\\.rds$", tolower(filename))) {
    ascat.bc <- readRDS(filename)
  } else {
    stop("Input file format not recognized: ", filename)
  }

  ascat.bc
}
