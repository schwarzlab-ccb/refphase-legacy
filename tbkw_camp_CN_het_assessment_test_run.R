
options(scipen=999);
library(compiler)
library(parallel)
library(data.table)
library(mclust)
library(GenomicRanges)
library(Rsamtools)
library(gtools)
library(diptest)
library(gplots)


#Debugging
patient <- "MET048"
paramsFile <- paste0("/camp/lab/swantonc/working/watkint/data/brastianos/release_v2/B_",patient,"/exome/commandHistory/B_",patient,".exome.params.txt")
designFile <- paste0("/camp/lab/swantonc/working/watkint/data/brastianos/release_v2/B_",patient,"/exome/commandHistory/B_",patient,".exome.design.txt")
params.patient <- read.delim(paramsFile,header=FALSE,sep="\t",stringsAsFactors=F);
for (r in 1:nrow(params.patient)) {
    print(params.patient[r,1])
    print(params.patient[r,2])
    assign(x=params.patient[r,1],value=params.patient[r,2]);
}
msubCount <- as.numeric(msubCount);
nProcs    <- as.numeric(nProcs);
mapq      <- as.numeric(mapq);
baseq     <- as.numeric(baseq);
overwrite <- as.logical(overwrite);
varscanDir <- as.character(varscanDir)
patientWorkDir <- as.character(patientWorkDir)
gatk.genomeFasta <- as.character(gatk.genomeFasta)
gatk.genomeIndex <- as.character(gatk.genomeIndex)
bam.location <- as.character(bamProcDir)
scriptDir <- as.character(scriptDir)
dataDir <- as.character(dataDir)
platypusDir <- as.character(platypusDir)
msaiDir <- as.character(msaiDir)

temp.scriptDir <- "/camp/lab/swantonc/working/watkint/scripts/tracerx-pipeline/"
CN.heterogeneity.dir <- paste(patientWorkDir,'CN_heterogeneity', sep="/")

#NOTE, we may not want this test of msai seg finishing to be the limit as we could rescue new samples. #TODO change this later.
msai.complete.file.path <- paste0(msaiDir, "/msai_results/", patient, "_msai_segs.txt")
print(msai.complete.file.path)
ascat.complete.file.path <- paste(ascatDir,patient,".ASCAT.complete",sep='')
print(ascat.complete.file.path)
#what directory are we looking at for msai regions:

joint.segs.file.path <- paste(msaiDir, patient, "_joint_segs.txt", sep='')
print(joint.segs.file.path)

joint.snps.file.path <- paste(msaiDir, patient, "_joint_snps.txt", sep='')
print(joint.snps.file.path)

full.joint.segs.file.path <- paste(msaiDir, patient, "_joint_segs_all_ascat_info.txt", sep='')

#How do we know that the reference phasing has run? Need the path to the ref phasing complete file.

CN.heterogeneity.out.complete <- paste(CN.heterogeneity.dir,"/",patient,".CN.heterogeneity.complete",sep='');

f.joint.segs.df <- read.table(full.joint.segs.file.path, stringsAsFactors=FALSE, sep="\t", header=TRUE)
    f.joint.segs.df$sample <- gsub("\\.", "-", gsub("_LogR", "", f.joint.segs.df$sample))
    f.joint.segs.df$callState <- NULL
    f.joint.segs.df$patient <- patient
    f.joint.segs.df$sample_name <- paste(patient,sapply(f.joint.segs.df$sample, function(x) strsplit(x, split="-")[[1]][2]), sep=":")


#The patient is supposed to match what's in the segs:

#full.patient <- paste("B_", patient, sep="")
array.type <- "exome"
pt.id <- patient
full.pt <- paste("B_",patient, sep="")
orig.cpn.segs.df <- f.joint.segs.df

source("/home/camp/watkint/working/watkint/scripts/multiregioncopynum/data_input_and_QC_functions.R")
source("/home/camp/watkint/working/watkint/scripts/multiregioncopynum/minimum_consistent_segmentation_functions.R")


print(head(orig.cpn.segs.df))
print(pt.id)
processed.orig.cpn.segs.df <- processOriginalCopyNumberSegments(pt.id, orig.cpn.segs.df)
print(head(processed.orig.cpn.segs.df))

ploidy.est.vec <- processed.orig.cpn.segs.df[!duplicated(processed.orig.cpn.segs.df$sample),'Ploidy']
names(ploidy.est.vec) <- unique(processed.orig.cpn.segs.df$sample)
purity.est.vec <- processed.orig.cpn.segs.df[!duplicated(processed.orig.cpn.segs.df$sample),'ACF']
names(purity.est.vec) <- unique(processed.orig.cpn.segs.df$sample)



#Todo accept precomputed mininmum consistent segmentation
min.cons.segs.df <- generateMinimumConsistentSegmentation(processed.orig.cpn.segs.df, chroms=1:22, col.to.use="ploidy_adj_tot_cpn")
print(head(min.cons.segs.df))


min.cons.segs.df <- cleanAndErrorCheckMinConsSegs(min.cons.segs.df)
print(head(min.cons.segs.df))
#Todo might need to add column names so looks like below:?

#Clean and check the SNP data, removing NAs from the segments and SNPs files.
#      Chr     Start       End B_MET002_MT-M1 B_MET002_T1-R1
# [1,]   1    762496  23885674      0.1993419     0.05387792


#This needs to be refactored further to take a dataframe of region names and ascat.bc object paths
#We should also rework the WGS to be put into ascat.bc objects or work out another data format but move on from here...
if(array.type=="WGS"){

  all.samp.SNP.list <- getSNPDataFromTRACERxFormatBAFLogRData(pt.id = full.pt, release.dir = workDir, array.type = array.type)

} else if(array.type%in%c("snp_array", "exome")){
  
  all.samp.SNP.list <- getSNPsWithBAFandLogRatioDataTRACERxFormatASCATBC(pt.id = full.pt, release.dir = workDir, array.type = array.type)

}
str(all.samp.SNP.list)

#all sample all seg logr baf lists is asas.lr.baf
#All SNPs in these lists must have both a logr and baf value
asas.SNPs.w.both.lr.baf.snp.lists <- cleanAndErrorCheckSNPs(all.samp.SNP.list, min.cons.segs.df, allow.logr.only.flag = FALSE)
#Todo, is it necessary to keep both lists of SNPs like this?
#SNPs that only have LogR are permitted in these lists, but it also contains all the SNPs with both BAF and LogR
asas.SNPS.lr.and.or.baf.lists <- cleanAndErrorCheckSNPs(all.samp.SNP.list, min.cons.segs.df, allow.logr.only.flag = TRUE)

#####Get relevant libraries:
#Todo check this is the only required library.
library(diptest)



###Source relevant scripts
source("/home/camp/watkint/working/watkint/scripts/multiregioncopynum/data_input_and_QC_functions.R")
source("/home/camp/watkint/working/watkint/scripts/multiregioncopynum/minimum_consistent_segmentation_functions.R")
source("/home/camp/watkint/working/watkint/scripts/multiregioncopynum/allelic_imbalance_phasing_functions.R")
source("/home/camp/watkint/working/watkint/scripts/multiregioncopynum/pairwise_CN_heterogeneity_comparison_functions.R")
source("/home/camp/watkint/working/watkint/scripts/multiregioncopynum/cancer_cell_fraction_calculation_functions.R")
source("/home/camp/watkint/working/watkint/scripts/multiregioncopynum/final_output_segmentation_file_generation_functions.R")
source("/home/camp/watkint/working/watkint/scripts/multiregioncopynum/join_segments_functions.R")
source("/home/camp/watkint/working/watkint/scripts/multiregioncopynum/CN_heterogeneity_plotting_functions.R")


###############This is past the input stage and into the analysis:
#Set constants for all analyses here... eventually these can be command line set
MIN_NUM_HZ_SNP_PHASING <- 5
PLOT_OUTPUT_DIR <- "/home/camp/watkint/working/watkint/scripts/multiregioncopynum_plots/"
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
  #Here is where we're going to have iniatalise on all the lists for either cpnA or cpnB
  all.seg.CCF.list[[cpnX]] <- list()
  all.seg.int.ob.xp.cpn.df.list[[cpnX]] <- list()
  all.seg.CCF.diff.from.exp.list[[cpnX]] <- list()
  all.seg.clonal.in.all.samps.list[[cpnX]] <- list()
  all.seg.greater.than.two.diff.list[[cpnX]] <- list()

  for ( seg.name in seg.name.vec ) {
    #Now we need to initialise lists for each segment for whichever cpnA or cpnB we are working on.
    all.seg.CCF.list[[cpnX]][[seg.name]] <- list()
    all.seg.int.ob.xp.cpn.df.list[[cpnX]][[seg.name]] <- list()
    all.seg.CCF.diff.from.exp.list[[cpnX]][[seg.name]] <- list()
    all.seg.clonal.in.all.samps.list[[cpnX]][[seg.name]] <- list()
    all.seg.greater.than.two.diff.list[[cpnX]][[seg.name]] <- list()

  }

  for (seg.name in seg.name.vec) {
    #Check if the phasing wasn't performed for this segment, if so set all values to NA and move on.
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
source("/home/camp/watkint/working/watkint/scripts/multiregioncopynum/final_output_segmentation_file_generation_functions.R")

#So no
full.joint.segs.df <- generateFinalSegmentationOutput(seg.name.vec = names(all.pw.cn.het.comp.lists[[1]])
                                ,all.samp.ph.BAF.df.list = patient.phased.BAF.df.list     
                                ,all.seg.CCF.vals.list = all.seg.CCF.list
                                ,all.seg.CCF.diff.from.exp.list = all.seg.CCF.diff.from.exp.list
                                ,ploidy.est.vec =  ploidy.est.vec
                                ,purity.est.vec = purity.est.vec
                                ,clonal.pval = PVAL_THRESH_VS_MIN_CCF#Todo median data is sometimes > the threshold but it still comes out highly significant... check this, looking at patient ME048
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






#Joining segment stuff here will have to be rewritten to match the functions I've rewritten



# #Join segments together for each sample... for now let's just ignore the samples that have fewer than 5 SNPs so they fall out entirely
# filt.segs.vec <- c()
# for ( seg.name in seg.name.vec ) {

#   if (length(all.pw.cn.het.comp.lists[[sample.vec[1]]][[seg.name]]$chr_vec) >= MIN_NUM_HZ_SNP_PHASING &
#       length(all.pw.cn.het.comp.lists[[sample.vec[1]]][[seg.name]]$ph_cpnA_vec) > 0 &
#       length(all.pw.cn.het.comp.lists[[sample.vec[1]]][[seg.name]]$ph_cpnB_vec) > 0) {
    
#     filt.segs.vec <- c(filt.segs.vec, seg.name)
    
#   }

# }#end for ( seg.name in seg.name.vec ) {

# #Now rejoin the segs for cpnA and cpnB separately for each sample so the
# new.joined.cpnA.cpnB.segs.list <- list()

# for ( samp in sample.vec ) {

#     all.segs.chr.pos.df <- do.call(rbind, strsplit(filt.segs.vec, split=":"))
#     all.segs.chr.pos.df <- apply(all.segs.chr.pos.df, 2, as.numeric)
#     rownames(all.segs.chr.pos.df) <- filt.segs.vec
#     new.joined.cpnA.cpnB.segs.list[[samp]] <- list()

#     for ( cpnX in c("cpnA", "cpnB")  ) {

#       samp.cpnX.segs.list <- list()

#       for ( chr in unique(all.segs.chr.pos.df[, 1]) ) {
#           #We've got an issue here that we sometimes have segments that are too small to be joined together, i.e. they don't have enough SNPs for the Wilcox tests
#           #This is because the number of het snps is under the minimum number for SNP phasing, so the only segments that can be joined by this method are those that
#           #are next to one another and have more than this minimum number of snps... We should only pass "runs" like these to the joining function...
#           #And flag the segments that couldn't be joined.
#           chr.seg.df <- all.segs.chr.pos.df[all.segs.chr.pos.df[, 1] %in% chr, ,drop = FALSE]
#           chr.joined.seg.list <- joinAdjacentSegsWithinSampleChromosome(samp, chr.seg.df, cpnX.string = cpnX, all.pw.cn.het.comp.lists, all.seg.CCF.vals.list, pval.to.separate = PVAL_THRESH_TO_JOIN_PH_SEGS)
#           samp.cpnX.segs.list  <- c(samp.cpnX.segs.list ,chr.joined.seg.list)

#       }#end for ( chrs in unique(all.segs.chr.pos.df[, 1]) ) {

#       new.joined.cpnA.cpnB.segs.list[[samp]][[paste0(cpnX, "_joined_segs_list")]] <- samp.cpnX.segs.list

#     }#end  for ( cpnX in c("cpnA", "cpnB")  ) {
# }#end for ( samp in sample.vec ) {


# #Now for each sample we need to introduce the breakpoints recaclulate gains/losses + CCFs based on the new segment sizes...
# #Todo, make all CCF calculations performed in the final segmentation and the CCF scripts the same by function calls, there's currently duplicated code there.
# new.cpn.segs.list <- list()
# for ( samp in sample.vec ) { 

#   new.cpn.segs.list[[samp]] <- list()
#   seg.cpnA.cnpB.diff.exp.list <- list()
#   seg.cpnA.cnpB.diff.exp.list[[samp]] <- list()

#   for ( cpnX in c("cpnA", "cpnB") ) {

#     seg.name.vec <- unique(names(new.joined.cpnA.cpnB.segs.list[[samp]][[paste0(cpnX, "_joined_segs_list")]]))


#     samps.int.cpn.raw.cpn.by.chr.list <- generateAllSamplesByChromIntCPNUnderlyingRawValuesLists(all.pw.cn.het.comp.lists = new.joined.cpnA.cpnB.segs.list[[samp]][[paste0(cpnX, "_joined_segs_list")]],
#                                                                                                   sample.vec = samp,
#                                                                                                   seg.name.vec = seg.name.vec,
#                                                                                                   cpnX.vec = c(cpnX),
#                                                                                                   reduced.seg.flag = TRUE)

#     for ( seg.name in seg.name.vec ) {

#       s.samps.int.ob.xp.cpn.df <- getSegObsExpCPNBySampleDF(all.pw.cn.het.comp.lists = new.joined.cpnA.cpnB.segs.list[[samp]][[paste0(cpnX, "_joined_segs_list")]],
#                                                             sample.vec = samp,
#                                                             seg.name = seg.name,
#                                                             cpnX.vec = c(cpnX),
#                                                             reduced.seg.flag = TRUE)

#       seg.cpnA.cnpB.diff.exp.list <- calculateWhetherSegCPNDiffFromIntegerExpectedAndAllMatchingSampCPN(samp.pvals.diff.from.exp.list = seg.cpnA.cnpB.diff.exp.list,
#                                                                                                         all.pw.cn.het.comp.lists = new.joined.cpnA.cpnB.segs.list[[samp]][[paste0(cpnX, "_joined_segs_list")]],
#                                                                                                         s.samps.int.ob.xp.cpn.df,
#                                                                                                         seg.name,
#                                                                                                         chrs.for.ccf.comp.vec,
#                                                                                                         samps.int.cpn.raw.cpn.by.chr.list,
#                                                                                                         cpnA.or.cpnB = cpnX,
#                                                                                                         reduced.seg.flag = TRUE)

#     }#end for ( seg.name in seg.name.vec ) {


#     #Todo we need pull out segment clonality calculation so that it can be passed to generate final segments


#     #functions to update everything? 
#     #Spout generate generate final segmenation ouput so it works per segment, not all wrapped up.
#     #So now we want to recalcuate everything that went into the final segmentation but for this new set of segments, except reduced seg = TRUE
#     cpnX.joint.segs.df <- generateFinalSegmentationOutput(seg.name.vec = seg.name.vec
#                                 ,all.samp.ph.BAF.df.list = patient.phased.BAF.df.list     
#                                 ,all.seg.CCF.vals.list = new.joined.cpnA.cpnB.segs.list[[samp]][[paste0(cpnX, "_joined_segs_list")]]
#                                 ,all.seg.CCF.diff.from.exp.list = seg.cpnA.cnpB.diff.exp.list
#                                 ,ploidy.est.vec =  ploidy.est.vec[samp]
#                                 ,purity.est.vec = purity.est.vec[samp]
#                                 ,clonal.pval = PVAL_THRESH_VS_MIN_CCF
#                                 ,diff.from.oth.p.val = PVAL_THRESH_VS_OTH_CHRS_SAME_CPN
#                                 ,min.CCF.to.consider.val = MIN_CCF_TO_CONSIDER
#                                 ,min.num.other.chrs = MIN_NUM_CHRS
#                                 ,gain.threshold    = GAIN_LOGR_THRESH
#                                 ,loss.threshold    = LOSS_LOGR_THRESH
#                                 ,amp.threshold     = AMP_LOGR_THRESH
#                                 ,deep.loss.threshold = DEEP_LOSS_LOGR_THRESH 
#                                 ,chrs.to.use.vec = CHRS_TO_USE_VEC
#                                 ,cpnX.vec = c(cpnX)
#                                 ,reduced.seg.flag = TRUE
#                                 ,max.cpn.val = MAX_CPN_FOR_SUBCLONALITY
#                                 )
#     new.cpn.segs.list[[samp]][[paste0("new_", cpnX, "_joined_segs_df")]] <- cpnX.joint.segs.df

#   }#end for ( cpnX in c("cpnA", "cpnB") ) {

#   ###So this can be a function call to get the combined dtartas and ends: that allow us to initialise the dataframe

#   #NOW take the A and B joined segments, and reintroduce breakpoints so that A and B have the same breakpoints. However we want to maintain the same CCFs that we got from the 
#   #larger segments before the reintroduction of breakpoints, we want to make sure that no recalculation needs to occur within the finalsegmentation output? or this populate?

#   new.samp.seg.name.vec <- getCPNSegNamesForSampleJoinedAB(new.cpn.segs.list, samp) 
#   #intialise the dataframe to fill in:
#   final.samp.comb.joined.cpnA.cpnB.seg.df <- initialiseAllSegFullDataframe(seg.name.vec = new.samp.seg.name.vec, samp, make.row.names.seg.name.alone = TRUE)

#   #so now we want to populate the dataframe
#   #Todo, may need to reclaculate CCFs etc based on these new boundaries? Or do we just keep the old ones that were likely more accurate as the came from larger uninterruped sections of cpnA
#   #and cpnB

#   for ( seg.name in new.samp.seg.name.vec) {


#     for ( cpnX in c("cpnA", "cpnB") ) {
      
#       #Todo, bad idea to keep replacing this variable.. do it by side effect of function...?
#       final.samp.comb.joined.cpnA.cpnB.seg.df <- populateOneSegOfCPNSampFinalJoinedDF(seg.name = seg.name,
#                                                                                       samp = samp,
#                                                                                       cpnX = cpnX,
#                                                                                       new.cpn.segs.list = new.cpn.segs.list,
#                                                                                       final.comb.seg.df = final.samp.comb.joined.cpnA.cpnB.seg.df)

#     }#end for ( cpnX in c("cpnA", "cpnB") ) {

#   }#for ( seg.name in new.seg.name.vec) {

# }#end for ( samp in sample.vec ) { 





                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 