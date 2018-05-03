
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
patient <- "MET018"
paramsFile <- "/camp/lab/swantonc/working/watkint/data/brastianos/release_v2/B_MET018/exome/commandHistory/B_MET018.exome.params.txt"
designFile <- "/camp/lab/swantonc/working/watkint/data/brastianos/release_v2/B_MET018/exome/commandHistory/B_MET018.exome.design.txt"
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
print(CN.heterogeneity.out.complete)

#source("/home/camp/watkint/working/watkint/scripts/modified_ascat/CN_heterogeneity_package_scripts/minimum_consistent_segmentation_functions.R")
#source("/home/camp/watkint/working/watkint/scripts/modified_ascat/CN_heterogeneity_package_scripts/data_input_and_QC_functions.R")
#source("/home/camp/watkint/working/watkint/scripts/modified_ascat/CN_heterogeneity_package_scripts/allelic_imbalance_phasing_functions.R")

f.joint.segs.df <- read.table(full.joint.segs.file.path, stringsAsFactors=FALSE, sep="\t", header=TRUE)
    f.joint.segs.df$sample <- gsub("\\.", "-", gsub("_LogR", "", f.joint.segs.df$sample))
    f.joint.segs.df$callState <- NULL
    f.joint.segs.df$patient <- patient
    f.joint.segs.df$sample_name <- paste(patient,sapply(f.joint.segs.df$sample, function(x) strsplit(x, split="-")[[1]][2]), sep=":")


#The patient is supposed to match what's in the segs:

source("/home/camp/watkint/working/watkint/scripts/multiregioncopynum/data_input_and_QC_functions.R")
#full.patient <- paste("B_", patient, sep="")
array.type <- "exome"
pt.id <- patient
full.pt <- paste("B_",patient, sep="")
orig.cpn.segs.df <- f.joint.segs.df

print(head(orig.cpn.segs.df))
print(pt.id)
processed.orig.cpn.segs.df <- processOriginalCopyNumberSegments(pt.id, orig.cpn.segs.df)
print(head(processed.orig.cpn.segs.df))

ploidy.est.vec <- processed.orig.cpn.segs.df[!duplicated(processed.orig.cpn.segs.df$sample),'Ploidy']
names(ploidy.est.vec) <- unique(processed.orig.cpn.segs.df$sample)
purity.est.vec <- processed.orig.cpn.segs.df[!duplicated(processed.orig.cpn.segs.df$sample),'ACF']
names(purity.est.vec) <- unique(processed.orig.cpn.segs.df$sample)




#Todo accept precomputed mininmum consistent segmentation
source("/home/camp/watkint/working/watkint/scripts/multiregioncopynum/minimum_consistent_segmentation_functions.R")
min.cons.segs.df <- generateMinimumConsistentSegmentation(processed.orig.cpn.segs.df, chroms=1:22, col.to.use="ploidy_adj_tot_cpn")
print(head(min.cons.segs.df))


source("/home/camp/watkint/working/watkint/scripts/multiregioncopynum/data_input_and_QC_functions.R")
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
asas.SNPs.w.both.lr.baf.snp.lists <- cleanAndErrorCheckSNPs(all.samp.SNP.list, min.cons.segs.df, allow.logr.only.flag=FALSE)
#Todo, is it necessary to keep both lists of SNPs like this?
#SNPs that only have LogR are permitted in these lists, but it also contains all the SNPs with both BAF and LogR
source("/home/camp/watkint/working/watkint/scripts/multiregioncopynum/data_input_and_QC_functions.R")
asas.SNPS.lr.and.or.baf.lists <- cleanAndErrorCheckSNPs(all.samp.SNP.list, min.cons.segs.df, allow.logr.only.flag=TRUE)


source("/home/camp/watkint/working/watkint/scripts/multiregioncopynum/allelic_imbalance_phasing_functions.R")



###############This is past the input stage and into the analysis:
#Set constants for all analyses here... eventually these can be command line set
MIN_NUM_HZ_SNP_PHASING <- 5


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
                                            gamma.val = 1)

source("/home/camp/watkint/working/watkint/scripts/multiregioncopynum/allelic_imbalance_phasing_functions.R")

#Todo unsure of what to do with the baf logr lists which have positions that only have logr here, at the moment exclude them but how
#does that affect later run
patient.phased.BAF.df.list <- createAllSegsPhasedBAFDataFrame(asas.lr.baf.lists = asas.SNPs.w.both.lr.baf.snp.lists,
                                                         sample.to.phase.seg.list = all.seg.phasing.info.list
                                                         )

#Add in plotting functions


#Todo check wha happens to the missing rows, are these segments that contained less than 5 snps and thus didn't have something created?
#what's going on.

source("/home/camp/watkint/working/watkint/scripts/multiregioncopynum/pairwise_CN_heterogeneity_comparison_functions.R")
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
                                                                  gamma.val = 1,
                                                                  min.num.hz.snp.phasing = MIN_NUM_HZ_SNP_PHASING
                                                                  )
    #Todo, no need o keep reassigning here...
    for ( seg.name in names(samp.one.two.comparison.list) ) {

      all.pw.cn.het.comp.lists[[samp.one]][[seg.name]] <- samp.one.two.comparison.list[[seg.name]][[samp.one]]
      all.pw.cn.het.comp.lists[[samp.two]][[seg.name]] <- samp.one.two.comparison.list[[seg.name]][[samp.two]]
      all.pw.cn.het.comp.lists[[samp.one.two]][[seg.name]] <- samp.one.two.comparison.list[[seg.name]][[samp.one.two]]

    }

    all.pw.cn.het.comp.lists[[samp.one.two]]$pairwise_avg_cn_het_score <- samp.one.two.comparison.list[[seg.name]][[samp.one.two]]$pairwise_avg_cn_het_score

  }#end for (i in 1:(length(sample.vec) - 1)) {

}#end for(j in (i + 1):length(sample.vec)) {

#Todo, add in plotting functions

#Let's get the CCF

source("/home/camp/watkint/working/watkint/scripts/multiregioncopynum/cancer_cell_fraction_calculation_functions.R")

seg.name.vec <- names(all.pw.cn.het.comp.lists[[1]])
#Check to make sure that all segments are present in all samples... why is this necessary
#Todo, check why this bug is ocurring
for ( samp in sample.vec ) {

      seg.name.vec <- intersect(seg.name.vec, names(all.pw.cn.het.comp.lists[[samp]]))

}#end for ( samp in sample.vec ) {   


seg.r.cpn.same.diff.across.samps.list <- generateAllSegsRoundedCPNSameOrDifferentAcrossSamplesLists(all.pw.cn.het.comp.lists, sample.vec, seg.name.vec)
#next assess for each chromosome and region, what are the different raw-values
samps.int.cpn.raw.cpn.by.chr.list <- generateAllSamplesByChromIntCPNUnderlyingRawValuesLists(all.pw.cn.het.comp.lists, sample.vec, seg.name.vec)

all.seg.CCF.vals.list <- list()
all.seg.CCF.diff.from.exp.list <- list()
for (seg.name in seg.name.vec) {

  #all.samp.ph.BAF.df.list used to be in the all.pw.cn.het.comp.lists but was refactored out, need to add a call or pass it to this function.
  seg.CCF.vals <- getIndividualCCFestimates_withPvals(seg.name = seg.name,
                                          purity.est.vec = purity.est.vec,
                                          all.pw.cn.het.comp.lists = all.pw.cn.het.comp.lists,
                                          all.samp.ph.BAF.df.list = patient.phased.BAF.df.list,
                                          samps.int.cpn.raw.cpn.by.chr.list = samps.int.cpn.raw.cpn.by.chr.list,
                                          max.cpn.val = 10,
                                          min.CCF.to.consider.val = 0.15,
                                          chrs.can.use.for.ccf.comp = c(1:22),
                                          confidence.interval.threshold = 0.9999,
                                          min.num.other.chrs = 3)
  print(str(seg.CCF.vals, max=1))
  if(  ! is.na(seg.CCF.vals) ){
    all.seg.CCF.diff.from.exp.list[[seg.name]] <- seg.CCF.vals[[seg.name]]$seg_diff_from_expected_list
    seg.CCF.vals$seg_diff_from_expected_list  <- NULL  
  } else {
    all.seg.CCF.diff.from.exp.list[[seg.name]] <- NA
  } 
    
  #deleting eg_diff_from_expected_list from seg.CCF.vals, I guess eventually we'd like to refactor these out into different functions
  all.seg.CCF.vals.list[[seg.name]] <- seg.CCF.vals
}
#Todo, fix the need to remove upper layer of list and rename
all.seg.CCF.vals.list <- unlist(all.seg.CCF.vals.list, recursive=FALSE)
names(all.seg.CCF.vals.list) <- seg.name.vec



#Todo, rewrite the plotting functions for the CCFs..
source("/home/camp/watkint/working/watkint/scripts/multiregioncopynum/final_output_segmentation_file_generation_functions.R")

#Todo, rewrite functions to get the ASCAT style copy number segmentation
full.joint.segs.df <- generateFinalSegmentationOutput(seg.name.vec = names(all.pw.cn.het.comp.lists[[1]])
                                ,all.samp.ph.BAF.df.list = patient.phased.BAF.df.list     
                                ,all.seg.CCF.vals.list = all.seg.CCF.vals.list
                                ,all.seg.CCF.diff.from.exp.list = all.seg.CCF.diff.from.exp.list
                                ,ploidy.est.vec =  ploidy.est.vec
                                ,purity.est.vec = purity.est.vec
                                ,clonal.pval = 0.001
                                ,diff.from.oth.p.val = 0.01
                                ,min.CCF.to.consider.val = 0.05
                                ,pval.to.separate = 0.001
                                ,min.num.other.chrs = 3
                                ,gain.threshold    = log2(2.5 / 2)
                                ,loss.threshold    = log2(1.5 / 2)
                                ,amp.threshold     = log2(4 / 2)
                                ,deep.loss.threshold= log2(0.5 / 2)
                                ,chrs.to.use.vec = 1:22
                                ,loss.gain.pval = 0.01  
                                )



full.joint.segs.w.phasing.df <- addPhasingInformationToOutputSegmentation(sample.vec = names(purity.est.vec),
                                        seg.name.vec = names(all.pw.cn.het.comp.lists[[1]]), 
                                        sample.to.phase.seg.list = all.seg.phasing.info.list,
                                        all.pw.cn.het.comp.lists = all.pw.cn.het.comp.lists,
                                        all.segs.df = full.joint.segs.df)


full.joint.segs.w.phasing.and.raw.df <- addOriginalCPNValuesToOutputSegmentation(sample.vec,
                                        seg.name.vec=names(all.pw.cn.het.comp.lists[[1]]), 
                                        processed.orig.cpn.segs.df,
                                        full.joint.segs.w.phasing.df)


#Todo plot how much of the original segmentation actually comes all the way through the code.



#Join segments together for each sample... for now let's just ignore the samples that have fewer than 5 SNPs...
source("/home/camp/watkint/working/watkint/scripts/multiregioncopynum/join_segments_functions.R")
#Pull out all pairwise comparisons
new.segs.list <- list()
filt.segs.vec <- c()
print(length(seg.name.vec))
for ( seg.name in seg.name.vec ) {

  if (length(all.pw.cn.het.comp.lists[[sample.vec[1]]][[seg.name]]$chr_vec) >= MIN_NUM_HZ_SNP_PHASING &
      length(all.pw.cn.het.comp.lists[[sample.vec[1]]][[seg.name]]$ph_cpnA_vec) > 0 &
      length(all.pw.cn.het.comp.lists[[sample.vec[1]]][[seg.name]]$ph_cpnB_vec) > 0) {
    
    filt.segs.vec <- c(filt.segs.vec, seg.name)
    
  }

}#end for ( seg.name in seg.name.vec ) {
print(length(seg.name.vec))

#Now rejoin the segs for cpnA and cpnB separately for each sample.
new.joined.cpnA.cpnB.segs.list <- list()

for ( samp in sample.vec ) {
    print(samp)
    all.segs.chr.pos.df <- do.call(rbind, strsplit(filt.segs.vec, split=":"))
    all.segs.chr.pos.df <- apply(all.segs.chr.pos.df, 2, as.numeric)
    rownames(all.segs.chr.pos.df) <- filt.segs.vec
    new.joined.cpnA.cpnB.segs.list[[samp]] <- list()

    for ( cpnX in c("cpnA", "cpnB")  ) {
      print(cpnX)
      samp.cpnX.segs.list <- list()

      for ( chr in unique(all.segs.chr.pos.df[,1]) ) {
          print(chr)
          #We've got an issue here that we sometimes have segments that are too small to be joined together, i.e. they don't have enough SNPs for the Wilcox tests
          #This is because the number of het snps is under the minimum number for SNP phasing, so the only segments that can be joined by this method are those that
          #are next to one another and have more than this minimum number of snps... We should only pass "runs" like these to the joining function...
          #And flag the segments that couldn't be joined.
          chr.seg.df <- all.segs.chr.pos.df[all.segs.chr.pos.df[,1] %in% chr, ,drop = FALSE]
          chr.joined.seg.list <- joinAdjacentSegsWithinSampleChromosome(samp, chr.seg.df, cpnX.string=cpnX, all.pw.cn.het.comp.lists, all.seg.CCF.vals.list, pval.to.separate = 0.001)
          samp.cpnX.segs.list  <- c(samp.cpnX.segs.list ,chr.joined.seg.list)

          # if(chr==2){
          #   print(names(chr.joined.seg.list))
          #   STOP
          # }
      }#end or ( chrs in all.segs.chr.pos.df[,1] ) {
        #This stupidity is to try and match list nesting...
        new.joined.cpnA.cpnB.segs.list[[samp]][[paste0(cpnX,"_joined_segs_list")]] <- samp.cpnX.segs.list
        print(samp.cpnX.segs.list)
    }#end  for ( cpnX in c("cpnA", "cpnB")  ) {s

}#end for ( samp in sample.vec ) {

print(str(new.joined.cpnA.cpnB.segs.list, max=1))
print(str(new.joined.cpnA.cpnB.segs.list[[1]], max=1))
#So next we need to reassess cpn changes in these larger joined segments.



source("/home/camp/watkint/working/watkint/scripts/multiregioncopynum/final_output_segmentation_file_generation_functions.R")
source("/home/camp/watkint/working/watkint/scripts/multiregioncopynum/cancer_cell_fraction_calculation_functions.R")

new.comb.cpn.seg.mat       <- c()

chrs.can.use.for.ccf.comp  <- c(1:22)
new.cpn.segs.list <- list()
for ( samp in sample.vec ) { 

  new.cpn.segs.list[[samp]] <- list()


  #get sample segs, need to do unique as we have repeated segments in each sample massively...., need to go back and find out why and remove
  seg.cpnA.cnpB.diff.exp.list <- list()
  seg.cpnA.cnpB.diff.exp.list[[samp]] <- list()
  for ( cpnX in c("cpnA", "cpnB") ) {

    seg.name.vec <- unique(names(new.joined.cpnA.cpnB.segs.list[[samp]][[paste0(cpnX, "_joined_segs_list")]]))

    samps.int.cpn.raw.cpn.by.chr.list <- generateAllSamplesByChromIntCPNUnderlyingRawValuesLists(all.pw.cn.het.comp.lists=new.joined.cpnA.cpnB.segs.list[[samp]][[paste0(cpnX, "_joined_segs_list")]], samp, seg.name.vec, cpnX.vec=c(cpnX), reduced.seg.flag=TRUE)

    print(str(samps.int.cpn.raw.cpn.by.chr.list, max=1))

    for ( seg.name in seg.name.vec ) {
      print("##################")
      print(seg.name)
      s.samps.int.ob.xp.cpn.df <- getSegObsExpCPNBySampleDF(all.pw.cn.het.comp.lists=new.joined.cpnA.cpnB.segs.list[[samp]][[paste0(cpnX, "_joined_segs_list")]], sample.vec=samp, seg.name=seg.name, cpnX.vec=c(cpnX), reduced.seg.flag=TRUE)

      chr <- strsplit(seg.name, split=":")[1]
      chrs.for.ccf.comp.vec <- chrs.can.use.for.ccf.comp[ ! chrs.can.use.for.ccf.comp%in% chr] #Going to be used to assess CCFs across the each sample

      seg.cpnA.cnpB.diff.exp.list <- calculateWhetherSegCPNDiffFromIntegerExpectedAndAllMatchingSampCPN(samp.pvals.diff.from.exp.list=seg.cpnA.cnpB.diff.exp.list,
                                                                                                        all.pw.cn.het.comp.lists=new.joined.cpnA.cpnB.segs.list[[samp]][[paste0(cpnX, "_joined_segs_list")]],
                                                                                                        s.samps.int.ob.xp.cpn.df,
                                                                                                        seg.name,
                                                                                                        chrs.for.ccf.comp.vec,
                                                                                                        samps.int.cpn.raw.cpn.by.chr.list,
                                                                                                        cpnA.or.cpnB=cpnX,
                                                                                                        reduced.seg.flag=TRUE)

    }

    #print(str(seg.cpnA.cnpB.diff.exp.list, max=1))
    #So 
    print(seg.name.vec)
    #So now we want to recalcuate everything that went into the final segmentation but for this new set of segments
    cpnX.joint.segs.df <- generateFinalSegmentationOutput(seg.name.vec = seg.name.vec
                                ,all.samp.ph.BAF.df.list = patient.phased.BAF.df.list     
                                ,all.seg.CCF.vals.list = new.joined.cpnA.cpnB.segs.list[[samp]][[paste0(cpnX, "_joined_segs_list")]]
                                ,all.seg.CCF.diff.from.exp.list = seg.cpnA.cnpB.diff.exp.list
                                ,ploidy.est.vec =  ploidy.est.vec[samp]
                                ,purity.est.vec = purity.est.vec[samp]
                                ,clonal.pval = 0.001
                                ,diff.from.oth.p.val = 0.01
                                ,min.CCF.to.consider.val =0.05
                                ,pval.to.separate = 0.001
                                ,min.num.other.chrs = 3
                                ,gain.threshold    = log2(2.5 / 2)
                                ,loss.threshold    = log2(1.5 / 2)
                                ,amp.threshold     = log2(4 / 2)
                                ,deep.loss.threshold= log2(0.5 / 2)
                                ,chrs.to.use.vec = 1:22
                                ,loss.gain.pval = 0.01
                                ,cpnX.vec = c(cpnX)
                                ,reduced.seg.flag=TRUE 
                                )
    new.cpn.segs.list[[samp]][[paste0("new_",cpnX,"_joined_segs_df")]] <- cpnX.joint.segs.df

  }#end for ( cpnX in c("cpnA", "cpnB") ) {

  #Now here we want to make a merged version
  samp.cpnA.df <- new.cpn.segs.list[[samp]][[paste0("new_cpnA_joined_segs_df")]]
  samp.cpnB.df <- new.cpn.segs.list[[samp]][[paste0("new_cpnB_joined_segs_df")]]
  unq.seg.name.vec <- unique(c(samp.cpnA.df$seg_name, samp.cpnB.df$seg_name))
  comb.segs.df <- do.call(rbind, strsplit(unq.seg.name.vec, split=":")) 
  comb.segs.df      <- comb.segs.df[order(as.numeric(comb.segs.df[,1]),as.numeric(comb.segs.df[,2]),as.numeric(comb.segs.df[,3])),,drop=FALSE]
  comb.segs.df  <- data.frame(comb.segs.df, stringsAsFactors=FALSE)
  comb.segs.df$seg_name <- paste(comb.segs.df[,1], comb.segs.df[,2], comb.segs.df[,3], sep=":")
  colnames(comb.segs.df)[1:3] <- c("chr", "start", "end")
  rownames(comb.segs.df) <- comb.segs.df$seg_name

  #Now check which names should be removed:

  for ( i in 2:nrow(comb.segs.df) ) {
    #Note not sure if we actually need to check that is less than nrow..., try without
    if (as.numeric(comb.segs.df[i,]$start) < as.numeric(comb.segs.df[ (i - 1) ,]$end) &
        comb.segs.df[i,]$chr == comb.segs.df[ (i - 1) , ]$chr &
        i <= nrow(comb.segs.df)) {

        #remove the offending row
        comb.segs.df <- comb.segs.df[ -i ,]

    }#end if (as.numeric(comb.segs.df[i,]$start...))

  }#end for ( i in 2:nrow(comb.segs.df) ) {

  #So now we can initialise the actual final joint segmentation matrix

  #function(seg.name.vec, sample.vec) {
  new.seg.name.vec <- rownames(comb.segs.df)
  final.comb.seg.df <- initialiseAllSegFullDataframe(seg.name.vec = new.seg.name.vec, samp, make.row.names.seg.name.alone=TRUE)

  #So as we don't know how many rows there are going to be in this
  #Todo really hese aren't accurate CCFs etc.. the values would actually change with the size of the segments..., this whole section 
  #needs a rethink..
  #Todo the start and end positions of the segments haven't really been sorted out, need to be redone.
  for ( seg.name in new.seg.name.vec) {

    for ( cpnX in c("cpnA", "cpnB") ) {

      raw.col.name <- ifelse(cpnX == "cpnA", "nAraw", "nBraw")

      samp.cpnX.df <- new.cpn.segs.list[[samp]][[paste0("new_",cpnX,"_joined_segs_df")]]

      if ( seg.name %in% samp.cpnX.df$seg_name ){

        #Now we need to update the dataframe if this worked.. can put into anot her function if we want
        #Todo add he additional columns that have been added since the update.
        final.comb.seg.df[seg.name, raw.col.name] <- samp.cpnX.df[seg.name, raw.col.name]
        final.comb.seg.df[seg.name, paste0("poss_sc_",cpnX,"_low")]  <- samp.cpnX.df[seg.name, paste0("poss_sc_",cpnX,"_low")]
        final.comb.seg.df[seg.name, paste0("poss_sc_",cpnX,"_high")]  <- samp.cpnX.df[seg.name, paste0("poss_sc_",cpnX,"_high")]
        final.comb.seg.df[seg.name, paste0("frac_psc_",cpnX,"_low")]  <- samp.cpnX.df[seg.name, paste0("frac_psc_",cpnX,"_low")]
        final.comb.seg.df[seg.name, paste0("frac_psc_",cpnX,"_high")]  <- samp.cpnX.df[seg.name, paste0("frac_psc_",cpnX,"_high")]
        final.comb.seg.df[seg.name, paste0("frac_psc_",cpnX,"_high")]  <- samp.cpnX.df[seg.name, paste0("frac_psc_",cpnX,"_high")]
        #Todo check that this is the right column:
        final.comb.seg.df[seg.name ,paste0("sc_",cpnX,"_present")] <- samp.cpnX.df[seg.name ,paste0("sc_",cpnX,"_present")]
        final.comb.seg.df[seg.name ,paste0("orig_",cpnX,"_raw")] <- samp.cpnX.df[seg.name ,paste0("orig_",cpnX,"_raw")]
        final.comb.seg.df[seg.name,"n.het"] <- samp.cpnX.df[seg.name, "n.het"]
        final.comb.seg.df[seg.name, paste0("sc_",cpnX,"_high_or_low_dominate")] <- samp.cpnX.df[seg.name, paste0("sc_",cpnX,"_high_or_low_dominate")]
        final.comb.seg.df[seg.name, paste0("sc_",cpnX,"_vs_min_CCF_pval")] <- samp.cpnX.df[seg.name, paste0("sc_",cpnX,"_vs_min_CCF_pval")]

        final.comb.seg.df[seg.name, paste0("num_other_chrs_same_",cpnX,"_cpn")] <- samp.cpnX.df[seg.name, paste0("num_other_chrs_same_",cpnX,"_cpn")]
        final.comb.seg.df[seg.name, paste0("vs_other_chrs_same_",cpnX,"_pval")] <- samp.cpnX.df[seg.name, paste0("vs_other_chrs_same_",cpnX,"_pval")]
        #Todo potentially add in the gains and losses

      } else {

        seg.chr <- as.numeric(unlist(strsplit(seg.name,split=":"))[1])
        seg.start <- as.numeric(unlist(strsplit(seg.name,split=":"))[2])
        seg.end <- as.numeric(unlist(strsplit(seg.name,split=":"))[3])
        subset.samp.cpnX.df <- samp.cpnX.df[samp.cpnX.df$chr == seg.chr, ,drop = FALSE]
        subset.samp.cpnX.df <- subset.samp.cpnX.df[as.numeric(subset.samp.cpnX.df$startpos) <= seg.start &
                                                    as.numeric(subset.samp.cpnX.df$endpos) >= seg.start, ,drop=FALSE]
        #Todo refacor so that these this updating is actually part of a function tow which you pass two dataframes?
        final.comb.seg.df[seg.name, raw.col.name] <- subset.samp.cpnX.df[, raw.col.name]
        final.comb.seg.df[seg.name, paste0("poss_sc_",cpnX,"_low")]  <- subset.samp.cpnX.df[, paste0("poss_sc_",cpnX,"_low")]
        final.comb.seg.df[seg.name, paste0("poss_sc_",cpnX,"_high")]  <- subset.samp.cpnX.df[, paste0("poss_sc_",cpnX,"_high")]
        final.comb.seg.df[seg.name, paste0("frac_psc_",cpnX,"_low")]  <- subset.samp.cpnX.df[, paste0("frac_psc_",cpnX,"_low")]
        final.comb.seg.df[seg.name, paste0("frac_psc_",cpnX,"_high")]  <- subset.samp.cpnX.df[, paste0("frac_psc_",cpnX,"_high")]
        final.comb.seg.df[seg.name, paste0("frac_psc_",cpnX,"_high")]  <- subset.samp.cpnX.df[, paste0("frac_psc_",cpnX,"_high")]
        #Todo check that this is the right column:
        final.comb.seg.df[seg.name ,paste0("sc_",cpnX,"_present")] <- subset.samp.cpnX.df[ ,paste0("sc_",cpnX,"_present")]
        final.comb.seg.df[seg.name ,paste0("orig_",cpnX,"_raw")] <- subset.samp.cpnX.df[ ,paste0("orig_",cpnX,"_raw")]
        final.comb.seg.df[seg.name,"n.het"] <- subset.samp.cpnX.df[, "n.het"]
        final.comb.seg.df[seg.name, paste0("sc_",cpnX,"_high_or_low_dominate")] <- subset.samp.cpnX.df[, paste0("sc_",cpnX,"_high_or_low_dominate")]
        final.comb.seg.df[seg.name, paste0("sc_",cpnX,"_vs_min_CCF_pval")] <- subset.samp.cpnX.df[, paste0("sc_",cpnX,"_vs_min_CCF_pval")]

        final.comb.seg.df[seg.name, paste0("num_other_chrs_same_",cpnX,"_cpn")] <- subset.samp.cpnX.df[, paste0("num_other_chrs_same_",cpnX,"_cpn")]
        final.comb.seg.df[seg.name, paste0("vs_other_chrs_same_",cpnX,"_pval")] <- subset.samp.cpnX.df[, paste0("vs_other_chrs_same_",cpnX,"_pval")]

      }

    }#end for ( cpnX in c("cpnA", "cpnB") ) {

  }#for ( seg.name in new.seg.name.vec) {
  
}#end for ( samp in sample.vec ) { 


source("/home/camp/watkint/working/watkint/scripts/multiregioncopynum/CN_heterogeneity_plotting_functions.R")
pdf("test.pdf")
funPlotAllRegions_withJoining(all.pw.cn.het.comp.lists, full.joint.segs.df=full.joint.segs.w.phasing.df, sample.vec, final.comb.seg.df, all.seg.CCF.vals.list, cpnX="cpnA")
par()
funPlotAllRegions_withJoining(all.pw.cn.het.comp.lists, full.joint.segs.df=full.joint.segs.w.phasing.df, sample.vec, final.comb.seg.df, all.seg.CCF.vals.list, cpnX="cpnB")
dev.off()


final.comb.seg.df

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         