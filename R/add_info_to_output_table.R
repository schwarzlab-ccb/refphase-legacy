# Another monster function to move code out of the main refphase script.
# Appends additional columns to the output table that contain useful
# information.

add_info_to_output_table <- function(asas.ph.raw.cpn.lists, ploidy.est.vec,
                                     purity.est.vec,
                                     all.samp.non.hz.logr.and.tot.cpn.df.list,
                                     patient.phased.BAF.df.list,
                                     asas.SNPs.w.both.lr.baf.snp.lists,
                                     all.seg.CCF.list,
                                     all.seg.CCF.diff.from.exp.list,
                                     all.seg.phasing.info.list,
                                     all.seg.MSAI.by.BAF.ltnpd.present.list,
                                     constants,
                                     sorted.all.bp.all.samps.segs.anno.df,
                                     sample.vec,
                                     seg.name.vec) {

  # Let's initialise the first joint output:
  logdebug("Initialising minimum consistent segmentation (MCS) output dataframe")
  empty.full.joint.segs.df <- initialiseFullOutputDataframe(
    seg.name.vec = names(asas.ph.raw.cpn.lists[[1]]),
    ploidy.est.vec = ploidy.est.vec,
    purity.est.vec = purity.est.vec,
    make.row.names.seg.name.alone = FALSE
  )

  # This adds the raw cn for both segments that could be phased and those that couldn't
  logdebug("Adding new copy number calls to MCS output based on the phased raw copy number")
  new.cpn.full.joint.segs.df <- addNewRawAndIntCPNCallsToFinalSegmentationParallel(
    seg.name.vec = names(asas.ph.raw.cpn.lists[[1]]),
    sample.vec = names(purity.est.vec),
    all.samp.logr.df.list = all.samp.non.hz.logr.and.tot.cpn.df.list,
    all.samp.ph.BAF.df.list = patient.phased.BAF.df.list,
    all.seg.CCF.vals.list = all.seg.CCF.list,
    asas.lr.baf.lists = asas.SNPs.w.both.lr.baf.snp.lists,
    all.seg.df = empty.full.joint.segs.df,
    cpnX.vec = c("cpnA", "cpnB")
  )

  # print("Mapping MCS segments back to original input and copy number states...")
  full.joint.segs.w.old.cpn.vs.ploidy.df <- addOriginalSegmentInfoToFinalSegmentationNew(new.cpn.full.joint.segs.df, sorted.all.bp.all.samps.segs.anno.df)

  logdebug("Adding CCF information and subclonal cpn calls to MCS output")
  clonality.and.ccf.joint.segs.df <- addWithinSampleClonalityAndCCFToFinalSegmentation(
    seg.name.vec = names(asas.ph.raw.cpn.lists[[1]]),
    all.seg.CCF.vals.list = all.seg.CCF.list,
    all.seg.CCF.diff.from.exp.list = all.seg.CCF.diff.from.exp.list,
    all.seg.df = full.joint.segs.w.old.cpn.vs.ploidy.df # new.cpn.full.joint.segs.df
    , sample.vec = names(purity.est.vec),
    clonal.pval = constants$PVAL_THRESH_VS_MIN_CCF,
    diff.from.oth.p.val = constants$PVAL_THRESH_VS_OTH_CHRS_SAME_CPN,
    min.CCF.to.consider.val = constants$MIN_CCF_TO_CONSIDER,
    min.cpn.to.consider.val = constants$MIN_CPN_TO_CONSIDER,
    min.num.other.chrs = constants$MIN_NUM_CHRS,
    chrs.to.use.vec = constants$CHRS_TO_USE_VEC,
    cpnX.vec = c("cpnA", "cpnB"),
    max.cpn.val = constants$MAX_CPN_FOR_SUBCLONALITY
  )

  logdebug("Adding total copy number estimates to MCS output for segments with only logr information")
  full.joint.segs.df.w.cpn.vs.ploidy.df <- addLogrDefinedCPNEventsToOutputSegmentationParallel(
    sample.vec = names(purity.est.vec),
    seg.name.vec = names(asas.ph.raw.cpn.lists[[1]]),
    all.seg.df = clonality.and.ccf.joint.segs.df,
    ploidy.est.vec = ploidy.est.vec,
    purity.est.vec = purity.est.vec,
    all.samp.logr.df.list = all.samp.non.hz.logr.and.tot.cpn.df.list,
    all.samp.ph.BAF.df.list = patient.phased.BAF.df.list,
    gain.threshold = constants$GAIN_LOGR_THRESH,
    loss.threshold = constants$LOSS_LOGR_THRESH,
    amp.threshold = constants$AMP_LOGR_THRESH,
    deep.loss.threshold = constants$DEEP_LOSS_LOGR_THRESH,
    tot.cpn.vs.ploidy.pval.thresh = constants$PVAL_THRESH_VS_LOGR_FOR_CPN_EVENT
  )

  logdebug("Adding cpn event vs ploidy by logr measures to MCS output")
  full.joint.segs.logr.vs.ploidy.clonality.df <- determineAcrossSampleClonalityOfLogrEventsVsPloidy(
    sample.vec = names(purity.est.vec),
    seg.name.vec = names(asas.ph.raw.cpn.lists[[1]]),
    all.seg.df = full.joint.segs.df.w.cpn.vs.ploidy.df,
    all.samp.logr.df.list = all.samp.non.hz.logr.and.tot.cpn.df.list,
    all.samp.ph.BAF.df.list = patient.phased.BAF.df.list,
    gain.threshold = constants$GAIN_LOGR_THRESH,
    loss.threshold = constants$LOSS_LOGR_THRESH,
    amp.threshold = constants$AMP_LOGR_THRESH,
    deep.loss.threshold = constants$DEEP_LOSS_LOGR_THRESH,
    tot.cpn.vs.ploidy.pval.thresh = constants$PVAL_THRESH_VS_LOGR_FOR_CPN_EVENT
  )

  logdebug("Testing whether any gains should be reclassified as amplifications")
  # full.joint.check.amps.vs.ploidy.clonality.df <- retestSubclonalAmps(all.seg.df = full.joint.segs.logr.vs.ploidy.clonality.df,
  # 																	all.samp.logr.df.list = all.samp.non.hz.logr.and.tot.cpn.df.list,
  # 																	all.samp.ph.BAF.df.list = patient.phased.BAF.df.list,
  # 																	gain.threshold = GAIN_LOGR_THRESH,
  # 																	amp.threshold = AMP_LOGR_THRESH,
  # 																	p.val.for.converting.gain.to.amp = PVAL_THRESH_CONVERT_GAIN_TO_AMP)

  # These are kept in separate functions as the original segmentaiton was refactored so it could run when segments are rejoined at the end of the script.
  logdebug("Adding phasing information and MSAI detection by raw cpn calls to MCS output")
  full.joint.segs.w.phasing.df <- addPhasingInformationToOutputSegmentation(
    sample.vec = names(purity.est.vec),
    seg.name.vec = names(asas.ph.raw.cpn.lists[[1]]),
    sample.to.phase.seg.list = all.seg.phasing.info.list,
    all.pw.cn.het.comp.lists = NULL,
    all.segs.df = full.joint.segs.logr.vs.ploidy.clonality.df
  )

  # #Now we want to add in the allelic imbalance clonal and subclonal inrormation to the final data frame.
  # print("Adding across sample clonal/subclonal summary to MCS output...")
  # full.joint.segs.w.phasing.and.raw.df <- addOriginalCPNValuesAndCallsVsPloidyToOutputSegmentation(sample.vec = sample.vec,
  #                                         seg.name.vec = names(asas.ph.raw.cpn.lists[[1]]),
  #                                         processed.orig.cpn.segs.df = processed.orig.cpn.segs.df,
  #                                         all.segs.df = full.joint.segs.w.phasing.df)

  logdebug("Adding MSAI detected by LTNPD method to the MCS output")
  full.joint.segs.w.ltnpd.MSAI.df <- addMSAIDetectionByLTNPDToOutputSegmentation(
    sample.vec = sample.vec,
    seg.name.vec = names(asas.ph.raw.cpn.lists[[1]]),
    all.seg.MSAI.by.BAF.ltnpd.present.list = all.seg.MSAI.by.BAF.ltnpd.present.list,
    all.segs.df = full.joint.segs.w.phasing.df
  )

  logdebug("Defining any parallel gain or loss events detected via MSAI and addingto MCS output")
  full.joint.segs.w.parallel.cpn.events.df <- addParallelEventsByMSAIToOutputSegmentation(all.segs.df = full.joint.segs.w.ltnpd.MSAI.df)

  logdebug("Adding comparison between original cpn calls and MPhase copy number calling to MCS output")
  full.joint.segs.within.and.changed.cpn.info.df <- addWithinSampleChangeToCpnAndAllelicImbalanceCalls(
    all.seg.df = full.joint.segs.w.parallel.cpn.events.df,
    sample.vec = sample.vec,
    seg.name.vec = unique(seg.name.vec)
  )

  logdebug("Adding across sample estimates of allelic imbalance and LOH clonality to MCS output")
  full.joint.segs.w.across.samp.AI.stats.df <- addAcrossSampleAlleleSpecificClonalityCalls(
    all.seg.df = full.joint.segs.within.and.changed.cpn.info.df,
    sample.vec = sample.vec,
    seg.name.vec = seg.name.vec
  )

  # Legacy mostly, trying to make definitions that Emilia was doing and put them ion here to ensure they were correct
  full.joint.segs.w.across.plotting.cols.df <- addPanCanColumnDefinitionToSegDF(full.joint.segs.w.across.samp.AI.stats.df, seg.name.vec = seg.name.vec, sample.vec = sample.vec)
  full.joint.segs.df <- full.joint.segs.w.across.plotting.cols.df

  full.joint.segs.df
}
