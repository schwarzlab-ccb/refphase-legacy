ccf_from_raw_cpn <- function(asas.ph.raw.cpn.lists, sample.vec, seg.name.vec) {
  message("*******************************************************************")
  message("Calculating cancer cell fractions (CCFs) using phased raw cpn")
  message("*******************************************************************")

  samps.int.cpn.raw.cpn.by.chr.list <- generateAllSamplesByChromIntCPNUnderlyingRawValuesLists(
    asas.ph.raw.cpn.lists = asas.ph.raw.cpn.lists,
    sample.vec = sample.vec,
    seg.name.vec = seg.name.vec,
    cpnX.vec = c("cpnA", "cpnB"),
    reduced.seg.flag = FALSE
  )
  if (!is.na(conf$debug_dir)) {
    saveRDS(samps.int.cpn.raw.cpn.by.chr.list, file = samps.int.cpn.raw.cpn.by.chr.list.path)
  }

  all.seg.CCF.list <- list()
  all.seg.CCF.diff.from.exp.list <- list()

  # These four variables are never used beyond this next code block
  all.seg.int.ob.xp.cpn.df.list <- list()
  all.seg.clonal.in.all.samps.list <- list()
  all.seg.greater.than.two.diff.list <- list()

  for (cpnX in c("cpnA", "cpnB")) {
    all.seg.CCF.list[[cpnX]] <- list()
    all.seg.int.ob.xp.cpn.df.list[[cpnX]] <- list()
    all.seg.CCF.diff.from.exp.list[[cpnX]] <- list()
    all.seg.clonal.in.all.samps.list[[cpnX]] <- list()
    all.seg.greater.than.two.diff.list[[cpnX]] <- list()

    for (seg.name in seg.name.vec) {
      all.seg.CCF.list[[cpnX]][[seg.name]] <- list()
      all.seg.int.ob.xp.cpn.df.list[[cpnX]][[seg.name]] <- list()
      all.seg.CCF.diff.from.exp.list[[cpnX]][[seg.name]] <- list()
      all.seg.clonal.in.all.samps.list[[cpnX]][[seg.name]] <- list()
      all.seg.greater.than.two.diff.list[[cpnX]][[seg.name]] <- list()
    }

    for (seg.name in seg.name.vec) {
      if (is.na(asas.ph.raw.cpn.lists[[1]][[seg.name]]$ph_cpnA_vec[1])) {
        all.seg.CCF.list[[cpnX]][[seg.name]] <- NA
        all.seg.int.ob.xp.cpn.df.list[[cpnX]][[seg.name]] <- NA
        all.seg.CCF.diff.from.exp.list[[cpnX]][[seg.name]] <- NA
        all.seg.clonal.in.all.samps.list[[cpnX]][[seg.name]] <- NA
        all.seg.greater.than.two.diff.list[[cpnX]][[seg.name]] <- NA
        next
      }

      # does both A and B, todo, refactor so only does one
      s.samps.int.ob.xp.cpn.df <- getSegObsExpCPNBySampleDF(
        asas.ph.raw.cpn.lists = asas.ph.raw.cpn.lists,
        sample.vec = sample.vec,
        seg.name = seg.name,
        cpnX.vec = cpnX,
        reduced.seg.flag = FALSE
      )

      # does only cpnX
      seg.cpnX.diff.from.exp.list <- calculateWhetherSegCPNDiffFromIntegerExpectedAndAllMatchingSampCPN(
        asas.ph.raw.cpn.lists = asas.ph.raw.cpn.lists,
        s.samps.int.ob.xp.cpn.df = s.samps.int.ob.xp.cpn.df,
        seg.name = seg.name,
        samps.int.cpn.raw.cpn.by.chr.list = samps.int.cpn.raw.cpn.by.chr.list,
        cpnX = cpnX,
        chrs.can.use.for.ccf.test = CHRS_TO_USE_VEC,
        reduced.seg.flag = FALSE
      )
      # does only cpnX
      is.segment.clonal.in.all.samps.flag <- determineCCFClonalityOfSegment(
        s.samps.int.ob.xp.cpn.df = s.samps.int.ob.xp.cpn.df,
        seg.name = seg.name,
        asas.ph.raw.cpn.lists = asas.ph.raw.cpn.lists,
        all.samp.ph.BAF.df.list = patient.phased.BAF.df.list,
        purity.est.vec = purity.est.vec,
        cpnX = cpnX,
        min.CCF.to.consider.val = MIN_CCF_TO_CONSIDER,
        confidence.interval.threshold = CONF_INT_FOR_CPNX_BAF_VS_EXP
      )

      # does only cpnX
      is.there.cpnX.difference.greater.than.two.flag <- calculateTwoValueDiffFlag(
        s.samps.int.ob.xp.cpn.df = s.samps.int.ob.xp.cpn.df,
        seg.cpnX.diff.from.exp.list = seg.cpnX.diff.from.exp.list,
        cpnX = cpnX,
        min.CCF.to.consider.val = MIN_CCF_TO_CONSIDER,
        confidence.interval.threshold = CONF_INT_FOR_CPNX_BAF_VS_EXP,
        min.num.other.chrs = MIN_NUM_CHRS
      )

      # does only cpnX
      seg.psc.low.high.cpn.and.use.preset.ccf.flag.vec <- getPossibleLowAndHighSubclonalCPNStates(
        seg.name = seg.name,
        cpnX = cpnX,
        s.samps.int.ob.xp.cpn.df = s.samps.int.ob.xp.cpn.df,
        all.clonal.flag = is.segment.clonal.in.all.samps.flag,
        two.value.diff.flag = is.there.cpnX.difference.greater.than.two.flag,
        min.CCF.to.consider.val = MIN_CCF_TO_CONSIDER
      )


      seg.psc.cpnX.low <- seg.psc.low.high.cpn.and.use.preset.ccf.flag.vec[1]
      seg.psc.cpnX.high <- seg.psc.low.high.cpn.and.use.preset.ccf.flag.vec[2]
      seg.use.default.ccf <- as.logical(seg.psc.low.high.cpn.and.use.preset.ccf.flag.vec[3])

      # update segment for each sample
      for (samp in sample.vec) {
        all.seg.CCF.list[[cpnX]][[seg.name]][[samp]] <- list()
        # ss.ph.cpnX.vals.vec <- all.pw.cn.het.comp.lists[[samp]][[seg.name]][[paste0("ph_", cpnX, "_vec")]]
        ss.ph.cpnX.vals.vec <- asas.ph.raw.cpn.lists[[samp]][[seg.name]][[paste0("ph_", cpnX, "_vec")]]

        seg.samp.CCF.psc.cpn.df <- getSegCPNXMeanAndAllRawCCFsDF(
          ss.mean.cpnX.obs = s.samps.int.ob.xp.cpn.df[samp, paste0("mean_obs_", cpnX)],
          ss.ph.cpnX.vals.vec = ss.ph.cpnX.vals.vec,
          seg.psc.cpnX.low = seg.psc.cpnX.low,
          seg.psc.cpnX.high = seg.psc.cpnX.high,
          seg.use.default.ccf.flag = seg.use.default.ccf,
          all.clonal.flag = is.segment.clonal.in.all.samps.flag,
          max.cpn.val = MAX_CPN_FOR_SUBCLONALITY
        )

        all.seg.CCF.list[[cpnX]][[seg.name]][[samp]][[paste0("poss_sc_", cpnX, "_low")]] <- seg.samp.CCF.psc.cpn.df$poss_sc_cpnX_low
        all.seg.CCF.list[[cpnX]][[seg.name]][[samp]][[paste0("poss_sc_", cpnX, "_high")]] <- seg.samp.CCF.psc.cpn.df$poss_sc_cpnX_high
        all.seg.CCF.list[[cpnX]][[seg.name]][[samp]][[paste0("mean_frac_", cpnX, "_low")]] <- seg.samp.CCF.psc.cpn.df$mean_frac_cpnX_low
        all.seg.CCF.list[[cpnX]][[seg.name]][[samp]][[paste0("mean_frac_", cpnX, "_high")]] <- seg.samp.CCF.psc.cpn.df$mean_frac_cpnX_high
        all.seg.CCF.list[[cpnX]][[seg.name]][[samp]][[paste0(cpnX, "_orig_ph_obs_cpnA_vec")]] <- asas.ph.raw.cpn.lists[[samp]][[seg.name]][[paste0("ph_cpnA_vec")]]
        all.seg.CCF.list[[cpnX]][[seg.name]][[samp]][[paste0(cpnX, "_orig_ph_obs_cpnB_vec")]] <- asas.ph.raw.cpn.lists[[samp]][[seg.name]][[paste0("ph_cpnB_vec")]]
        all.seg.CCF.list[[cpnX]][[seg.name]][[samp]][[paste0("clean_ph_obs_", cpnX, "_vec")]] <- seg.samp.CCF.psc.cpn.df$clean_ph_obs_cpnX
        all.seg.CCF.list[[cpnX]][[seg.name]][[samp]][[paste0("frac_", cpnX, "_low_vec")]] <- seg.samp.CCF.psc.cpn.df$frac_cpnX_low
        all.seg.CCF.list[[cpnX]][[seg.name]][[samp]][[paste0("frac_", cpnX, "_high_vec")]] <- seg.samp.CCF.psc.cpn.df$frac_cpnX_high
      } # end  for ( samp  in sample.vec ) {

      all.seg.int.ob.xp.cpn.df.list[[cpnX]][[seg.name]] <- s.samps.int.ob.xp.cpn.df
      all.seg.clonal.in.all.samps.list[[cpnX]][[seg.name]] <- is.segment.clonal.in.all.samps.flag
      all.seg.greater.than.two.diff.list[[cpnX]][[seg.name]] <- is.there.cpnX.difference.greater.than.two.flag
      all.seg.CCF.diff.from.exp.list[[cpnX]][[seg.name]] <- seg.cpnX.diff.from.exp.list
    }
  }

  # Save CCF data:
  if (!is.na(conf$debug_dir)) {
    saveRDS(all.seg.CCF.list, file = all.seg.CCF.list.path)
    saveRDS(all.seg.int.ob.xp.cpn.df.list, file = all.seg.int.ob.xp.cpn.df.list.path)
    saveRDS(all.seg.CCF.diff.from.exp.list, file = all.seg.CCF.diff.from.exp.list.path)
    saveRDS(all.seg.clonal.in.all.samps.list, file = all.seg.clonal.in.all.samps.list.path)
    saveRDS(all.seg.greater.than.two.diff.list, file = all.seg.greater.than.two.diff.list.path)
  }

  list(all.seg.CCF.list = all.seg.CCF.list, all.seg.CCF.diff.from.exp.list = all.seg.CCF.diff.from.exp.list)
}
