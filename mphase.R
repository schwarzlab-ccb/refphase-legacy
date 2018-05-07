#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
patient <- args[1]
snp.file.path <- args[2]
seg.file.path <- args[3]
output.directory.path <- args[4]

#####Get relevant libraries:
#Todo check this is the only required library.
library(diptest)


###Source relevant scripts
source("./data_input_and_QC_functions.R")
source("./minimum_consistent_segmentation_functions.R")
source("./allelic_imbalance_phasing_functions.R")
source("./pairwise_CN_heterogeneity_comparison_functions.R")
source("./cancer_cell_fraction_calculation_functions.R")
source("./final_output_segmentation_file_generation_functions.R")
source("./join_segments_functions.R")
source("./CN_heterogeneity_plotting_functions.R")


###############This is past the input stage and into the analysis:
#Set constants for all analyses here... eventually these can be command line set
MIN_NUM_HZ_SNP_PHASING <- 5
PLOT_OUTPUT_DIR <- NULL
GAMMA_VAL <- 1
MAX_CPN_FOR_SUBCLONALITY <- 10
MIN_CCF_TO_CONSIDER <- 0.15
MIN_NUM_CHRS <- 3
PVAL_THRESH_VS_MIN_CCF <- 0.001
PVAL_THRESH_VS_OTH_CHRS_SAME_CPN <- 0.01
GAIN_LOGR_THRESH <- log2(2.5 / 2)
LOSS_LOGR_THRESH <- log2(1.5 / 2)
AMP_LOGR_THRESH <- log2(4 / 2)
DEEP_LOSS_LOGR_THRESH <- log2(0.5 / 2)
CHRS_TO_USE_VEC <- 1:22
CONF_INT_FOR_CPNX_BAF_VS_EXP <- 0.9999
PVAL_THRESH_TO_JOIN_PH_SEGS <- 0.001
MIN_PROP_SNP_SEP_THRESH <- 0.75
PVAL_COMB_WILCOX_PH_THRESH <- 0.01
PVAL_DIPTEST_AI_THRESH <- 0.01
ALLOW_DIPTEST_ALONE_PH_FLAG <- FALSE

if (is.null(output.directory.path)) {

	stop("No output path provided, aborting...")

}


#so now we need to read in the MSAI data:
if ( ! file.exists(snp.file.path)) {

	stop("No file at the snp.file.path provided in the second argument")

}

if ( ! file.exists(seg.file.path)) {

	stop("No file at the seg.file.path provided in the third argument")

}

#For purposes of debugging
snp.file.path <- "/home/camp/watkint/working/watkint/data/brastianos/release_v2/B_MET048/exome/MSAI/MET048_joint_snps.txt"

seg.file.path <- "/home/camp/watkint/working/watkint/data/brastianos/release_v2/B_MET048/exome/MSAI/MET048_joint_segs_all_ascat_info.txt"

snp.file.df <- read.table(snp.file.path, stringsAsFactors=FALSE, sep="\t", header=TRUE)
#Note this should be the full oen
seg.file.df <- read.table(seg.file.path, stringsAsFactors=FALSE, sep="\t", header=TRUE)
#Designed to handle the MSAI nput files
seg.file.df$sample <- gsub("\\.", "-", gsub("_LogR", "", seg.file.df$sample))
seg.file.df$callState <- NULL
seg.file.df$patient <- patient
seg.file.df$sample_name <- paste(patient,sapply(f.joint.segs.df$sample, function(x) strsplit(x, split="-")[[1]][2]), sep=":")

processed.orig.cpn.segs.df <- processOriginalCopyNumberSegments(patient, seg.file.df)

samples <- unique(orig.cpn.segs.df$sample)

purity.estimates <- processed.orig.cpn.segs.df[ ! duplicated(processed.orig.cpn.segs.df$sample), "ACF"]
names(purity.estimates) <- samples

ploidy.estimates <- processed.orig.cpn.segs.df[ ! duplicated(processed.orig.cpn.segs.df$sample), "Ploidy"]
names(ploidy.estimates) <- samples

#Todo accept precomputed mininmum consistent segmentation
min.cons.segs.df <- generateMinimumConsistentSegmentation(processed.orig.cpn.segs.df, chroms=1:22, col.to.use="ploidy_adj_tot_cpn")
print(head(min.cons.segs.df))

min.cons.segs.df <- cleanAndErrorCheckMinConsSegs(min.cons.segs.df)
#colnames(min.cons.segs.df)[1:3] <- c("chrom", "startpos", "endpos")
print(head(min.cons.segs.df))


#Ok so now we want to load in the snps from this kind of file:
#              Name Chr Position ref alt Scaled_CADD_score B_MET048_MT.M1_LogR B_MET048_MT.M1_tot B_MET048_MT.M1_var B_MET048_MT.M1_BAF B_MET048_T1.R1_LogR B_MET048_T1.R1_tot B_MET048_T1.R1_var B_MET048_T1.R1_BAF
# 1 chr1:762496:G:G   1   762496   G   G                 0               0.509                  0                  0                  2               0.108                  0                  0                  2

# samples in sample.vec look like this:
# [1] "B_MET048_MT-M1" "B_MET048_T1-R1"

#So we need to tmake a fata frame for each of the 
all.samp.SNP.list <- list()
for ( samp in sample.vec ) {

	ASCAT.summary <- data.frame(chr = snp.file.df$Chr,
		                					pos = snp.file.df$Pos,
		   									BAF = snp.file.df[, paste0(gsub(samp, "-" ,"\\.","_BAF")])],
		   									BAFsegmented = rep(NA, length(snp.file.df$Chr)),
		   									logr = snp.file.df[, paste0(gsub(samp, "-", "\\.", "_LogR"))],
		   									logRsegmented = rep(NA, length(snp.file.df$Chr)),
		   									stringsAsFactors = FALSE)

	#So I guess we want the mnimum consistent segmentation? or the original?
	#           sample chr  startpos    endpos n.het cnTotal nMajor nMinor   Ploidy  ACF     nAraw       nBraw patient sample_name tot_raw_cpn max_raw_cpn min_raw_cpn cpn_to_use ploidy_adj_tot_cpn
	# 1 B_MET048_MT-M1   1    762496   4834651  4329       2      1      1 2.293931 0.83 0.8747692 0.874769163  MET048   MET048:M1   1.7495383   0.8747692 0.874769163  1.7495383         -0.3908477
	# 2 B_MET048_MT-M1   1   4847589  34684306 18442       1      1      0 2.293931 0.83 0.8888665 0.000000000  MET048   MET048:M1   0.8888665   0.8888665 0.000000000  0.8888665         -1.3677832

	seg <- processed.orig.cpn.segs.df[processed.orig.cpn.segs.df$sample==sample,]

	#Update the fake ASCAT.summaary BAFsegmetned and logRsegmented
	for ( j in 1:nrow(seg) ) {

	  if(nrow(ASCAT.summary[ASCAT.summary$chrom==seg[j,]$chr & ASCAT.summary$pos >= seg[j,]$startpos & ASCAT.summary$pos <= seg[j,]$endpos,]) > 0) {

	    if(seg[j,]$nMajor==seg[j,]$nMinor){
	      
	      ASCAT.summary[ASCAT.summary$chrom==seg[j,]$chr & ASCAT.summary$pos >= seg[j,]$startpos & ASCAT.summary$pos <= seg[j,]$endpos,]$BAFsegmented <- 0.5
	      
	    } else {
	      
	      avg.seg.baf <- 0.5 - mean(abs(tumor.baf_logr[tumor.baf_logr$chrom==seg[j,]$chr & tumor.baf_logr$pos >= seg[j,]$startpos & tumor.baf_logr$pos <= seg[j,]$endpos,]$BAF - 0.5))
	      ASCAT.summary[ASCAT.summary$chrom==seg[j,]$chr & ASCAT.summary$pos >= seg[j,]$startpos & ASCAT.summary$pos <= seg[j,]$endpos,]$BAFsegmented <- avg.seg.baf
	      
	    }

	    avg.seg.logr <- mean(tumor.baf_logr[tumor.baf_logr$chrom==seg[j,]$chr & tumor.baf_logr$pos >= seg[j,]$startpos & tumor.baf_logr$pos <= seg[j,]$endpos, regions[i]])
	    ASCAT.summary[ASCAT.summary$chrom==seg[j,]$chr & ASCAT.summary$pos >= seg[j,]$startpos & ASCAT.summary$pos <= seg[j,]$endpos,]$logRsegmented <- avg.seg.logr

	  }#end if(nrow(ASCAT.summary[ASCAT.summary$chrom==seg[j,]$chr & ASCAT.summary$pos...

	}#end for ( j in 1:nrow(seg) ) {

	all.samp.SNP.list[[samp]] <- ASCAT.summary	                					
}#end for ( samp in sample.vec ) {


asas.SNPs.w.both.lr.baf.snp.lists <- cleanAndErrorCheckSNPs(all.samp.SNP.list, min.cons.segs.df, allow.logr.only.flag = FALSE)

asas.SNPS.lr.and.or.baf.lists <- cleanAndErrorCheckSNPs(all.samp.SNP.list, min.cons.segs.df, allow.logr.only.flag = TRUE)


#Todo, further refactor code so that everything is done at the level of a segment
#Todo, add options to phase using baf dip test not just ascat regions
#Todo, reduce the number of parameters per function
#Todo, Add function definitions documentaiton
#Todo, clarify the loops and next statements that are controlling what happens when there are less than 5 hz snps in a sample.
#Todo, check if the mixture phasing is ever being selected... in the example with B_MET002 no mixture phasing called.

all.seg.phasing.info.list <- getCombinedBAFDipTestPhasing(asas.lr.baf.lists = asas.SNPs.w.both.lr.baf.snp.lists,
                                            min.cons.segs.df = min.cons.segs.df,
                                            ploidy.est.vec = ploidy.est.vec,
                                            purity.est.vec = purity.est.vec,
                                            gamma.val = GAMMA_VAL,
                                            min.num.hz.snp.phasing = MIN_NUM_HZ_SNP_PHASING,
                                            min.prop.snp.sep.threshold = MIN_PROP_SNP_SEP_THRESH,
                                            dip.test.p.val.threshold = PVAL_DIPTEST_AI_THRESH,
                                            comb.wilcox.snp.sep.p.val.threshold = PVAL_COMB_WILCOX_PH_THRESH,
                                            allow.dip.test.without.ASCAT.to.phase = ALLOW_DIPTEST_ALONE_PH_FLAG
                                            )


#Todo unsure of what to do with the baf logr lists which have positions that only have logr here, at the moment exclude them but how
#does that affect later run
patient.phased.BAF.df.list <- createAllSegsPhasedBAFDataFrame(asas.lr.baf.lists = asas.SNPs.w.both.lr.baf.snp.lists,
                                                         sample.to.phase.seg.list = all.seg.phasing.info.list,
                                                         min.num.hz.snp.phasing = MIN_NUM_HZ_SNP_PHASING
                                                         )


#Todo check wha happens to the missing rows, are these segments that contained less than 5 snps and thus didn't have something created?

sample.vec <- names(asas.SNPs.w.both.lr.baf.snp.lists)

all.pw.cn.het.comp.lists <- list()
for (i in 1:(length(sample.vec) - 1)) {

  for(j in (i + 1):length(sample.vec)) {

    samp.one <- sample.vec[i]
    samp.two <- sample.vec[j]
    samp.one.two <- paste(samp.one, ":", samp.two, sep="")

    samp.one.two.comparison.list <- calculatePairwiseCNHeterogeneity(samp.one = samp.one,
                                                                  samp.two = samp.two,
                                                                  asas.both.lr.baf.lists = asas.SNPs.w.both.lr.baf.snp.lists,
                                                                  asas.lr.and.or.baf.lists = asas.SNPS.lr.and.or.baf.lists,
                                                                  min.cons.segs.df = min.cons.segs.df,
                                                                  ploidy.est.vec = ploidy.est.vec,
                                                                  purity.est.vec = purity.est.vec,
                                                                  all.seg.phasing.info.list = all.seg.phasing.info.list,
                                                                  gamma.val = GAMMA_VAL,
                                                                  min.num.hz.snp.phasing = MIN_NUM_HZ_SNP_PHASING,
                                                                  gain.threshold = GAIN_LOGR_THRESH,
                                                                  loss.threshold = LOSS_LOGR_THRESH,
                                                                  amp.threshold = AMP_LOGR_THRESH,
                                                                  deep.loss.threshold = DEEP_LOSS_LOGR_THRESH
                                                                  )
    #Todo, no need to keep reassigning here...
    for ( seg.name in names(samp.one.two.comparison.list) ) {

      all.pw.cn.het.comp.lists[[samp.one]][[seg.name]] <- samp.one.two.comparison.list[[seg.name]][[samp.one]]
      all.pw.cn.het.comp.lists[[samp.two]][[seg.name]] <- samp.one.two.comparison.list[[seg.name]][[samp.two]]
      all.pw.cn.het.comp.lists[[samp.one.two]][[seg.name]] <- samp.one.two.comparison.list[[seg.name]][[samp.one.two]]

    }

    all.pw.cn.het.comp.lists[[samp.one.two]]$pairwise_avg_cn_het_score <- samp.one.two.comparison.list[[seg.name]][[samp.one.two]]$pairwise_avg_cn_het_score

  }#end for (i in 1:(length(sample.vec) - 1)) {

}#end for(j in (i + 1):length(sample.vec)) {



#Let's get the CCF
seg.name.vec <- names(all.pw.cn.het.comp.lists[[1]])
#Check to make sure that all segments are present in all samples... might not be necessary todo, check.
for ( samp in sample.vec ) {

      seg.name.vec <- intersect(seg.name.vec, names(all.pw.cn.het.comp.lists[[samp]]))

}#end for ( samp in sample.vec ) {   


source("/home/camp/watkint/working/watkint/scripts/multiregioncopynum/cancer_cell_fraction_calculation_functions.R")
all.seg.CCF.list <- list()
all.seg.int.ob.xp.cpn.df.list <- list()
all.seg.CCF.diff.from.exp.list <- list()
all.seg.clonal.in.all.samps.list <- list()
all.seg.greater.than.two.diff.list <- list()
#next assess for each chromosome and region, what are the the raw values underlying CN calls, these will be used to help test for subclonality.
samps.int.cpn.raw.cpn.by.chr.list <- generateAllSamplesByChromIntCPNUnderlyingRawValuesLists(all.pw.cn.het.comp.lists, sample.vec, seg.name.vec)

for (cpnX in c("cpnA", "cpnB")) {

  #Here is where we're going to have iniatalise on all the lists
  all.seg.CCF.list[[cpnX]] <- list()
  all.seg.int.ob.xp.cpn.df.list[[cpnX]] <- list()
  all.seg.CCF.diff.from.exp.list[[cpnX]] <- list()
  all.seg.clonal.in.all.samps.list[[cpnX]] <- list()
  all.seg.greater.than.two.diff.list[[cpnX]] <- list()

  for ( seg.name in seg.name.vec ) {

    all.seg.CCF.list[[cpnX]][[seg.name]] <- list()
    all.seg.int.ob.xp.cpn.df.list[[cpnX]][[seg.name]] <- list()
    all.seg.CCF.diff.from.exp.list[[cpnX]][[seg.name]] <- list()
    all.seg.clonal.in.all.samps.list[[cpnX]][[seg.name]] <- list()
    all.seg.greater.than.two.diff.list[[cpnX]][[seg.name]] <- list()

  }

  for (seg.name in seg.name.vec) {

    if (is.na(all.pw.cn.het.comp.lists[[1]][[seg.name]]$ph_cpnA_vec[1])) {

      all.seg.CCF.list[[cpnX]][[seg.name]] <- NA
      all.seg.int.ob.xp.cpn.df.list[[cpnX]][[seg.name]] <- NA
      all.seg.CCF.diff.from.exp.list[[cpnX]][[seg.name]] <- NA
      all.seg.clonal.in.all.samps.list[[cpnX]][[seg.name]] <- NA
      all.seg.greater.than.two.diff.list[[cpnX]][[seg.name]] <- NA
      next 

    }

    #does both A and B, todo, refactor so only does one
    s.samps.int.ob.xp.cpn.df <- getSegObsExpCPNBySampleDF(all.pw.cn.het.comp.lists = all.pw.cn.het.comp.lists,
                                                          sample.vec = sample.vec,
                                                          seg.name = seg.name,
                                                          cpnX.vec = cpnX,
                                                          reduced.seg.flag = FALSE) 

    #does only cpnX
    seg.cpnX.diff.from.exp.list <- calculateWhetherSegCPNDiffFromIntegerExpectedAndAllMatchingSampCPN(all.pw.cn.het.comp.lists = all.pw.cn.het.comp.lists,
                                                                                                      s.samps.int.ob.xp.cpn.df = s.samps.int.ob.xp.cpn.df,
                                                                                                      seg.name = seg.name,
                                                                                                      samps.int.cpn.raw.cpn.by.chr.list = samps.int.cpn.raw.cpn.by.chr.list,
                                                                                                      cpnX = cpnX,
                                                                                                      chrs.can.use.for.ccf.test = CHRS_TO_USE_VEC,
                                                                                                      reduced.seg.flag = FALSE)
    #does only cpnX
    is.segment.clonal.in.all.samps.flag <- determineCCFClonalityOfSegment(s.samps.int.ob.xp.cpn.df = s.samps.int.ob.xp.cpn.df,
                                                                          seg.name = seg.name,
                                                                          all.samp.ph.BAF.df.list = patient.phased.BAF.df.list,
                                                                          purity.est.vec = purity.est.vec,
                                                                          cpnX = cpnX,
                                                                          min.CCF.to.consider.val = MIN_CCF_TO_CONSIDER,
                                                                          confidence.interval.threshold = CONF_INT_FOR_CPNX_BAF_VS_EXP)

    #does only cpnX
    is.there.cpnX.difference.greater.than.two.flag <- calculateTwoValueDiffFlag(s.samps.int.ob.xp.cpn.df = s.samps.int.ob.xp.cpn.df,
                                                                seg.cpnX.diff.from.exp.list = seg.cpnX.diff.from.exp.list,
                                                                cpnX = cpnX,
                                                                min.CCF.to.consider.val = MIN_CCF_TO_CONSIDER,
                                                                confidence.interval.threshold = CONF_INT_FOR_CPNX_BAF_VS_EXP,
                                                                min.num.other.chrs = MIN_NUM_CHRS)

    #does only cpnX
    seg.psc.low.high.cpn.and.use.preset.ccf.flag.vec <- getPossibleLowAndHighSubclonalCPNStates(seg.name = seg.name,
                                                                                                cpnX = cpnX,
                                                                                                s.samps.int.ob.xp.cpn.df = s.samps.int.ob.xp.cpn.df,
                                                                                                all.clonal.flag = is.segment.clonal.in.all.samps.flag,
                                                                                                two.value.diff.flag = is.there.cpnX.difference.greater.than.two.flag)
    seg.psc.cpnX.low  <- seg.psc.low.high.cpn.and.use.preset.ccf.flag.vec[1]
    seg.psc.cpnX.high  <- seg.psc.low.high.cpn.and.use.preset.ccf.flag.vec[2]
    seg.use.default.ccf <- as.logical(seg.psc.low.high.cpn.and.use.preset.ccf.flag.vec[3])

    #update segment for each sample 
    for ( samp  in sample.vec ) {

      all.seg.CCF.list[[cpnX]][[seg.name]][[samp]] <- list()
      ss.ph.cpnX.vals.vec <- all.pw.cn.het.comp.lists[[samp]][[seg.name]][[paste0("ph_", cpnX, "_vec")]]

      seg.samp.CCF.psc.cpn.df <- getSegCPNXMeanAndAllRawCCFsDF(ss.mean.cpnX.obs = s.samps.int.ob.xp.cpn.df[samp, paste0("mean_obs_", cpnX)],
                                                              ss.ph.cpnX.vals.vec = ss.ph.cpnX.vals.vec,
                                                              seg.psc.cpnX.low = seg.psc.cpnX.low,
                                                              seg.psc.cpnX.high = seg.psc.cpnX.high,
                                                              seg.use.default.ccf.flag = seg.use.default.ccf,
                                                              all.clonal.flag = is.segment.clonal.in.all.samps.flag,
                                                              max.cpn.val = MAX_CPN_FOR_SUBCLONALITY)

      all.seg.CCF.list[[cpnX]][[seg.name]][[samp]][[paste0("poss_sc_", cpnX, "_low")]] <- seg.samp.CCF.psc.cpn.df$poss_sc_cpnX_low
      all.seg.CCF.list[[cpnX]][[seg.name]][[samp]][[paste0("poss_sc_", cpnX, "_high")]] <- seg.samp.CCF.psc.cpn.df$poss_sc_cpnX_high
      all.seg.CCF.list[[cpnX]][[seg.name]][[samp]][[paste0("mean_frac_", cpnX, "_low")]] <- seg.samp.CCF.psc.cpn.df$mean_frac_cpnX_low
      all.seg.CCF.list[[cpnX]][[seg.name]][[samp]][[paste0("mean_frac_", cpnX, "_high")]] <- seg.samp.CCF.psc.cpn.df$mean_frac_cpnX_high
      all.seg.CCF.list[[cpnX]][[seg.name]][[samp]][[paste0(cpnX, "_orig_ph_obs_cpnA_vec")]]   <- all.pw.cn.het.comp.lists[[samp]][[seg.name]][[paste0("ph_cpnA_vec")]]
      all.seg.CCF.list[[cpnX]][[seg.name]][[samp]][[paste0(cpnX, "_orig_ph_obs_cpnB_vec")]]   <- all.pw.cn.het.comp.lists[[samp]][[seg.name]][[paste0("ph_cpnB_vec")]]
      all.seg.CCF.list[[cpnX]][[seg.name]][[samp]][[paste0("clean_ph_obs_", cpnX, "_vec")]]   <- seg.samp.CCF.psc.cpn.df$clean_ph_obs_cpnX
      all.seg.CCF.list[[cpnX]][[seg.name]][[samp]][[paste0("frac_", cpnX, "_low_vec")]]  <- seg.samp.CCF.psc.cpn.df$frac_cpnX_low 
      all.seg.CCF.list[[cpnX]][[seg.name]][[samp]][[paste0("frac_", cpnX, "_high_vec")]]  <- seg.samp.CCF.psc.cpn.df$frac_cpnX_high

    }#end  for ( samp  in sample.vec ) {

    all.seg.int.ob.xp.cpn.df.list[[cpnX]][[seg.name]] <- s.samps.int.ob.xp.cpn.df
    all.seg.clonal.in.all.samps.list[[cpnX]][[seg.name]] <- is.segment.clonal.in.all.samps.flag
    all.seg.greater.than.two.diff.list[[cpnX]][[seg.name]] <- is.there.cpnX.difference.greater.than.two.flag
    all.seg.CCF.diff.from.exp.list[[cpnX]][[seg.name]] <- seg.cpnX.diff.from.exp.list

  }#for (seg.name %in% seg.name.vec 
    #for the moment these stays here so we don't wipe it, the list should be updated twice so that A and B are stored within it

}#end for (cpnX in c("cpnA", "cpnB")) {

#Todo remove the dataframe initialisation from the function and modify it to work on a segment by segment basis.
#Todo rename clonal.pval to be more meaningful...
#Todo Check [1] "B_MET019_MT-M2:1:762496:53608176" to see if LOW/HIGH fraction working properly.

#Ok so our issue is we're always getting the low because the CCFs are not being recalculated correctly.
#We need to record whether a segment is clonal within a sample or not

#Let's initialise the first joint output:

#empty.unjoined.mcr.segs.df <- initialiseAllSegFullDataframe(seg.name.vec, sample.vec, make.row.names.seg.name.alone = FALSE)

#So no
full.joint.segs.df <- generateFinalSegmentationOutput(seg.name.vec = names(all.pw.cn.het.comp.lists[[1]])
                                ,all.samp.ph.BAF.df.list = patient.phased.BAF.df.list     
                                ,all.seg.CCF.vals.list = all.seg.CCF.list
                                ,all.seg.CCF.diff.from.exp.list = all.seg.CCF.diff.from.exp.list
                                ,ploidy.est.vec =  ploidy.est.vec
                                ,purity.est.vec = purity.est.vec
                                ,clonal.pval = PVAL_THRESH_VS_MIN_CCF
                                ,diff.from.oth.p.val = PVAL_THRESH_VS_OTH_CHRS_SAME_CPN
                                ,min.CCF.to.consider.val = MIN_CCF_TO_CONSIDER
                                ,min.num.other.chrs = MIN_NUM_CHRS
                                ,gain.threshold    = GAIN_LOGR_THRESH
                                ,loss.threshold    = LOSS_LOGR_THRESH
                                ,amp.threshold     = AMP_LOGR_THRESH
                                ,deep.loss.threshold = DEEP_LOSS_LOGR_THRESH
                                ,chrs.to.use.vec = CHRS_TO_USE_VEC
                                ,cpnX.vec = c("cpnA", "cpnB")
                                ,reduced.seg.flag = FALSE
                                ,max.cpn.val = MAX_CPN_FOR_SUBCLONALITY
                                )

#These are kept in separate functions as the original segmentaiton was refactored so it could run when segments are rejoined at the end of the script.
full.joint.segs.w.phasing.df <- addPhasingInformationToOutputSegmentation(sample.vec = names(purity.est.vec),
                                        seg.name.vec = names(all.pw.cn.het.comp.lists[[1]]), 
                                        sample.to.phase.seg.list = all.seg.phasing.info.list,
                                        all.pw.cn.het.comp.lists = all.pw.cn.het.comp.lists,
                                        all.segs.df = full.joint.segs.df)

full.joint.segs.w.phasing.and.raw.df <- addOriginalCPNValuesToOutputSegmentation(sample.vec,
                                        seg.name.vec = names(all.pw.cn.het.comp.lists[[1]]), 
                                        processed.orig.cpn.segs.df,
                                        full.joint.segs.w.phasing.df)


#Todo plot how much of the original segmentation actually comes all the way through the code.
full.joint.segs.df <- full.joint.segs.w.phasing.and.raw.df




source("/home/camp/watkint/working/watkint/scripts/multiregioncopynum/CN_heterogeneity_plotting_functions.R")


png(file = paste(PLOT_OUTPUT_DIR, "CCF_cpnA", ".png", sep = ""), width = 1200, height = 1400)

plotAllSampleCCFsAcrossGenome(all.pw.cn.het.comp.lists = all.pw.cn.het.comp.lists,
                              full.joint.segs.df = full.joint.segs.w.phasing.df,
                              sample.vec = sample.vec,
                              all.seg.CCF.vals.list = all.seg.CCF.list,
                              cpnX = "cpnA")

dev.off()

png(file = paste(PLOT_OUTPUT_DIR, "CCF_cpnB", ".png", sep = ""), width = 1200, height = 1400)

plotAllSampleCCFsAcrossGenome(all.pw.cn.het.comp.lists = all.pw.cn.het.comp.lists,
                              full.joint.segs.df = full.joint.segs.w.phasing.df,
                              sample.vec = sample.vec,
                              all.seg.CCF.vals.list = all.seg.CCF.list,
                              cpnX = "cpnB")
dev.off()

#Now let's plot baf
png(file = paste(PLOT_OUTPUT_DIR, "BAF_phased", ".png", sep=""), width = 1200, height = 1400)

plotAllSampleBAFsAcrossGenome(sample.vec, patient.phased.BAF.df.list, cpnA.allele.color = '#de2d2699', cpnB.allele.color= '#3182bd99')

dev.off()

#Now for these sample by sample comparisons the names of PNGs should contain the two samples being being compared really so we'll have to the plotting within the functions

plotOldNewComparisonOfAllCPN(sample.vec, full.joint.segs.df, plot.output.dir = PLOT_OUTPUT_DIR)


plotAllSamplePairwiseComparisonSegmentCPN(sample.vec, full.joint.segs.df, plot.output.dir = PLOT_OUTPUT_DIR)
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 