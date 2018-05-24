##
##Params:
##
##Returns:
##
populateOneSegOfCPNSampFinalJoinedDF <- function(seg.name, samp, cpnX, new.cpn.segs.list, final.comb.seg.df) {

  raw.col.name <- ifelse(cpnX == "cpnA", "nAraw", "nBraw")

  samp.cpnX.df <- new.cpn.segs.list[[samp]][[paste0("new_", cpnX, "_joined_segs_df")]]

  if ( seg.name %in% samp.cpnX.df$seg_name ){

    #Now we need to update the dataframe if this worked.. can put into anot her function if we want
    #Todo add he additional columns that have been added since the update.
    final.comb.seg.df[seg.name, raw.col.name] <- samp.cpnX.df[seg.name, raw.col.name]
    final.comb.seg.df[seg.name, paste0("poss_sc_", cpnX, "_low")]  <- samp.cpnX.df[seg.name, paste0("poss_sc_", cpnX, "_low")]
    final.comb.seg.df[seg.name, paste0("poss_sc_", cpnX, "_high")]  <- samp.cpnX.df[seg.name, paste0("poss_sc_", cpnX, "_high")]
    final.comb.seg.df[seg.name, paste0("frac_psc_", cpnX, "_low")]  <- samp.cpnX.df[seg.name, paste0("frac_psc_", cpnX, "_low")]
    final.comb.seg.df[seg.name, paste0("frac_psc_", cpnX, "_high")]  <- samp.cpnX.df[seg.name, paste0("frac_psc_", cpnX, "_high")]
    final.comb.seg.df[seg.name, paste0("frac_psc_", cpnX, "_high")]  <- samp.cpnX.df[seg.name, paste0("frac_psc_", cpnX, "_high")]
    #Todo check that this is the right column:
    final.comb.seg.df[seg.name ,paste0("sc_", cpnX, "_present")] <- samp.cpnX.df[seg.name, paste0("sc_", cpnX, "_present")]
    final.comb.seg.df[seg.name ,paste0("orig_", cpnX, "_raw")] <- samp.cpnX.df[seg.name, paste0("orig_", cpnX, "_raw")]
    final.comb.seg.df[seg.name, "n.het"] <- samp.cpnX.df[seg.name, "n.het"]
    final.comb.seg.df[seg.name, paste0("sc_", cpnX, "_high_or_low_dominate")] <- samp.cpnX.df[seg.name, paste0("sc_" ,cpnX, "_high_or_low_dominate")]
    final.comb.seg.df[seg.name, paste0("sc_", cpnX, "_vs_min_CCF_pval")] <- samp.cpnX.df[seg.name, paste0("sc_", cpnX, "_vs_min_CCF_pval")]

    final.comb.seg.df[seg.name, paste0("num_other_chrs_same_", cpnX, "_cpn")] <- samp.cpnX.df[seg.name, paste0("num_other_chrs_same_", cpnX, "_cpn")]
    final.comb.seg.df[seg.name, paste0("vs_other_chrs_same_", cpnX, "_pval")] <- samp.cpnX.df[seg.name, paste0("vs_other_chrs_same_", cpnX, "_pval")]
    #Todo potentially add in the gains and losses

  } else {#end if ( seg.name %in% samp.cpnX.df$seg_name ){

    seg.chr <- as.numeric(unlist(strsplit(seg.name, split = ":"))[1])
    seg.start <- as.numeric(unlist(strsplit(seg.name, split = ":"))[2])
    seg.end <- as.numeric(unlist(strsplit(seg.name, split = ":"))[3])
    subset.samp.cpnX.df <- samp.cpnX.df[samp.cpnX.df$chr == seg.chr, ,drop = FALSE]
    subset.samp.cpnX.df <- subset.samp.cpnX.df[as.numeric(subset.samp.cpnX.df$startpos) <= seg.start &
                                                as.numeric(subset.samp.cpnX.df$endpos) >= seg.start, ,drop = FALSE]
    #Todo refacor so that these this updating is actually part of a function tow which you pass two dataframes?
    final.comb.seg.df[seg.name, raw.col.name] <- subset.samp.cpnX.df[, raw.col.name]
    final.comb.seg.df[seg.name, paste0("poss_sc_", cpnX, "_low")]  <- subset.samp.cpnX.df[, paste0("poss_sc_", cpnX, "_low")]
    final.comb.seg.df[seg.name, paste0("poss_sc_", cpnX, "_high")]  <- subset.samp.cpnX.df[, paste0("poss_sc_", cpnX, "_high")]
    final.comb.seg.df[seg.name, paste0("frac_psc_", cpnX, "_low")]  <- subset.samp.cpnX.df[, paste0("frac_psc_", cpnX, "_low")]
    final.comb.seg.df[seg.name, paste0("frac_psc_", cpnX, "_high")]  <- subset.samp.cpnX.df[, paste0("frac_psc_", cpnX, "_high")]
    final.comb.seg.df[seg.name, paste0("frac_psc_", cpnX, "_high")]  <- subset.samp.cpnX.df[, paste0("frac_psc_", cpnX, "_high")]
    #Todo check that this is the right column:
    final.comb.seg.df[seg.name ,paste0("sc_", cpnX, "_present")] <- subset.samp.cpnX.df[ ,paste0("sc_", cpnX, "_present")]
    final.comb.seg.df[seg.name ,paste0("orig_", cpnX, "_raw")] <- subset.samp.cpnX.df[ ,paste0("orig_", cpnX, "_raw")]
    final.comb.seg.df[seg.name,"n.het"] <- subset.samp.cpnX.df[, "n.het"]
    final.comb.seg.df[seg.name, paste0("sc_", cpnX, "_high_or_low_dominate")] <- subset.samp.cpnX.df[, paste0("sc_", cpnX, "_high_or_low_dominate")]
    final.comb.seg.df[seg.name, paste0("sc_", cpnX, "_vs_min_CCF_pval")] <- subset.samp.cpnX.df[, paste0("sc_", cpnX, "_vs_min_CCF_pval")]

    final.comb.seg.df[seg.name, paste0("num_other_chrs_same_", cpnX, "_cpn")] <- subset.samp.cpnX.df[, paste0("num_other_chrs_same_", cpnX, "_cpn")]
    final.comb.seg.df[seg.name, paste0("vs_other_chrs_same_", cpnX, "_pval")] <- subset.samp.cpnX.df[, paste0("vs_other_chrs_same_", cpnX, "_pval")]

  }#end if ( ! seg.name %in% samp.cpnX.df$seg_name ){

  return(final.comb.seg.df)

}#end populateCPNSampFinalJoinedDF <- function(...)



getCPNSegNamesForSampleJoinedAB <- function(new.cpn.segs.list, samp) {

  samp.cpnA.df <- new.cpn.segs.list[[samp]][[paste0("new_cpnA_joined_segs_df")]]
  samp.cpnB.df <- new.cpn.segs.list[[samp]][[paste0("new_cpnB_joined_segs_df")]]
  unq.seg.name.vec <- unique(c(samp.cpnA.df$seg_name, samp.cpnB.df$seg_name))
  comb.segs.df <- do.call(rbind, strsplit(unq.seg.name.vec, split=":")) 
  comb.segs.df      <- comb.segs.df[order(as.numeric(comb.segs.df[, 1]), as.numeric(comb.segs.df[, 2]), as.numeric(comb.segs.df[, 3])), ,drop = FALSE]
  comb.segs.df  <- data.frame(comb.segs.df, stringsAsFactors = FALSE)
  comb.segs.df$seg_name <- paste(comb.segs.df[, 1], comb.segs.df[, 2], comb.segs.df[, 3], sep = ":")
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

  new.seg.name.vec <- rownames(comb.segs.df)

  return(new.seg.name.vec)
}#end getCPNSegNamesForSampleJoinedAB(...)


##Params:
##
##Returns:
##
getSegForJoiningDF <- function(cpnX.string, samp, seg.name, all.pw.cn.het.comp.lists, all.seg.CCF.vals.list) {


    num.rows <- length(all.pw.cn.het.comp.lists[[samp]][[seg.name]][[paste0("ph_cpnA_vec")]])

    seg.to.join.df <- data.frame(orig_ph_cpnA_vec = all.pw.cn.het.comp.lists[[samp]][[seg.name]][[paste0("ph_cpnA_vec")]],
                orig_ph_cpnB_vec = all.pw.cn.het.comp.lists[[samp]][[seg.name]][[paste0("ph_cpnB_vec")]],
               chr = all.pw.cn.het.comp.lists[[samp]][[seg.name]]$chr_vec,
               pos = all.pw.cn.het.comp.lists[[samp]][[seg.name]]$pos_vec,
               low_cpnX_frac = all.seg.CCF.vals.list[[seg.name]][[samp]][[paste0("frac_", cpnX.string, "_low_vec")]],
               high_cpnX_frac = all.seg.CCF.vals.list[[seg.name]][[samp]][[paste0("frac_", cpnX.string, "_high_vec")]], 
                cpnX_frac = all.seg.CCF.vals.list[[seg.name]][[samp]][[paste0("frac_", cpnX.string, "_high_vec")]],
               low_sc_cpnX = rep(all.seg.CCF.vals.list[[seg.name]][[samp]][[paste0("poss_sc_", cpnX.string, "_low")]], num.rows),
               high_sc_cpnX = rep(all.seg.CCF.vals.list[[seg.name]][[samp]][[paste0("poss_sc_", cpnX.string, "_high")]], num.rows),
               stringsAsFactors = FALSE)

    return(seg.to.join.df)
}

updateJoinedSegList <- function(samp, cpnX.string, samp.cpnX.joined.segs.list, seg.name, seg.df) {

    start <- as.numeric(unlist(strsplit(seg.name,split=":"))[2])
    end <- as.numeric(unlist(strsplit(seg.name,split=":"))[3])

    if (end < start){
      
      print(seg.name)
      
    }

    samp.cpnX.joined.segs.list[[seg.name]][[paste0("ph_", cpnX.string, "_vec")]] <- seg.df$cpnX
    samp.cpnX.joined.segs.list[[seg.name]][["chr_vec"]] <- seg.df$chr
    samp.cpnX.joined.segs.list[[seg.name]][["pos_vec"]] <- seg.df$pos
    #Added the original ph cpn so that this can be kept track of as well
    samp.cpnX.joined.segs.list[[seg.name]][[paste0(cpnX.string,"_orig_ph_obs_cpnA_vec")]] <- seg.df$orig_ph_cpnA_vec
    samp.cpnX.joined.segs.list[[seg.name]][[paste0(cpnX.string,"_orig_ph_obs_cpnB_vec")]] <- seg.df$orig_ph_cpnB_vec
    ######
    samp.cpnX.joined.segs.list[[seg.name]][[paste0("frac_", cpnX.string, "_low_vec")]] <- seg.df$low_cpnX_frac
    samp.cpnX.joined.segs.list[[seg.name]][[paste0("frac_", cpnX.string, "_high_vec")]] <- seg.df$high_cpnX_frac
    samp.cpnX.joined.segs.list[[seg.name]][[paste0(cpnX.string, "_frac")]] <- seg.df$cpnX_frac
    samp.cpnX.joined.segs.list[[seg.name]][[paste0("poss_sc_", cpnX.string, "_low")]] <- seg.df$low_sc_cpnX
    samp.cpnX.joined.segs.list[[seg.name]][[paste0("poss_sc_", cpnX.string, "_high")]] <- seg.df$high_sc_cpnX

    return(samp.cpnX.joined.segs.list)

}#end updateJoinedSegList

wilcoxTestWithErrorHandling <- function(x, y = NULL,
            alternative = c("two.sided", "less", "greater"),
            mu = 0, paired = FALSE, exact = FALSE, correct = TRUE,
            conf.int = FALSE, conf.level = 0.95, errors.not.permitted="")
             {
    test.results <- tryCatch({wilcox.test(x=x, y=y, alternative=alternative, mu=mu, paired=paired, exact=exact, correct=correct, conf.int = conf.int, conf.level = conf.level)},
       error = function(err) {
          message(err)
          write(paste(err$message, "\n"), stderr())
          write("\n", stderr())
          if(grepl("NULL", err$message) & "NULL" %in% errors.not.permitted){
            stop("t.test ERROR with a NULL")
          } else if (grepl("enough", err$message) & "enough" %in% errors.not.permitted){
            stop("t.test ERROR with not enough")
          } else if (grepl("constant", err$message) & "constant" %in% errors.not.permitted){
            stop("t.test ERROR with constant data")
          } else if(grepl("enough|constant|NULL", err$message) & "other" %in% errors.not.permitted){
            stop(err$message)
          }
          # cat("\n")
          # paste("Input for wilcox test")
          # print("x")
          # print(x)
          # print("y")
          # print(y)
          # print("mu")
          # print(mu)
          # print("paired")
          # print(paired)
          # print("conf.level")
          # print(conf.level)
          test.results <- list()
          test.results[["stat"]] <- 1
          test.results[["conf.int"]] <- c(NA, NA)
          test.results[["p.val"]] <- 1
          return(test.results)
        })

        if (is.na(test.results$p.val)) {

            test.results$p.val <- 1
        }

        return(test.results)
}#end wilcoxTestWithErrorHandling
  

shouldSegsBeJoined <- function(s.curr.df, s.next.df, pval.to.separate) {

    cpnX.vec.s.curr.vs.s.next.test <- wilcoxTestWithErrorHandling(s.curr.df$cpnX, s.next.df$cpnX)
    frac.vec.s.curr.vs.s.next.test <- wilcoxTestWithErrorHandling(s.curr.df$cpnX_frac, s.next.df$cpnX_frac)

    if (cpnX.vec.s.curr.vs.s.next.test$p.val > pval.to.separate & frac.vec.s.curr.vs.s.next.test$p.val > pval.to.separate) {

        return(TRUE)

    } else {

        return(FALSE)
    }

}#end shouldSegsBeJoined


##
##Params:
##
##Returns:
##
joinAdjacentSegsWithinSampleChromosome <- function(samp, chr.seg.df, cpnX.string, all.pw.cn.het.comp.lists, all.seg.CCF.vals.list, pval.to.separate) {
    chr.seg.names.vec <- rownames(chr.seg.df)
    s.curr.name <- chr.seg.names.vec[1] 
    s.curr.df <- getSegForJoiningDF(cpnX.string, samp, seg.name = s.curr.name, all.pw.cn.het.comp.lists, all.seg.CCF.vals.list)
    
    #samp.cpnX.segs.list 
    samp.cpnX.segs.list <- list()
    #Doing this so it's compatible with generate final segmentation function... may want to refactor.
    if (length(chr.seg.names.vec) == 1) {

        samp.cpnX.segs.list <- updateJoinedSegList(samp, cpnX.string, samp.cpnX.joined.segs.list=samp.cpnX.segs.list, s.curr.name, s.curr.df)

    } else {

        s.next.name <- chr.seg.names.vec[2]
        s.next.df <- getSegForJoiningDF(cpnX.string, samp, seg.name=s.next.name, all.pw.cn.het.comp.lists, all.seg.CCF.vals.list)

    }# if (length(chr.seg.names.vec) > 1)
    #Check if less than 3 segments
    if (length(chr.seg.names.vec) == 2) {

            k <- 1

            if (shouldSegsBeJoined(s.curr.df, s.next.df, pval.to.separate) == TRUE) {

                new.seg.name <- paste(unlist(strsplit(s.curr.name,split=":"))[1],unlist(strsplit(s.curr.name,split=":"))[2],unlist(strsplit(s.next.name,split=":"))[3],sep=":")
                s.curr.df <- rbind(s.curr.df, s.next.df)
                s.curr.name <- new.seg.name
                samp.cpnX.segs.list <- updateJoinedSegList(samp, cpnX.string, samp.cpnX.joined.segs.list=samp.cpnX.segs.list, new.seg.name, s.curr.df)

            } else {

                samp.cpnX.segs.list <- updateJoinedSegList(samp, cpnX.string, samp.cpnX.joined.segs.list=samp.cpnX.segs.list, s.curr.name, s.curr.df)
                samp.cpnX.segs.list <- updateJoinedSegList(samp, cpnX.string, samp.cpnX.joined.segs.list=samp.cpnX.segs.list, s.next.name, s.next.df)

            }#end if (cpnX.vec.s.one.vs.s.two.test$p.val <= pval.to.separate &

    }# end if (length(chr.seg.names.vec) == 2) {
    k <- 3
    while (k <= length(chr.seg.names.vec)) {
        if (shouldSegsBeJoined(s.curr.df, s.next.df, pval.to.separate) == TRUE) {

            new.seg.name <- paste(unlist(strsplit(s.curr.name,split=":"))[1],unlist(strsplit(s.curr.name,split=":"))[2],unlist(strsplit(s.next.name,split=":"))[3],sep=":")
            s.curr.df <- rbind(s.curr.df, s.next.df)
            #Get next seg to consider
            s.curr.name <- new.seg.name
            s.next.name <- chr.seg.names.vec[k]
            s.next.df <- getSegForJoiningDF(cpnX.string, samp, seg.name=s.next.name, all.pw.cn.het.comp.lists, all.seg.CCF.vals.list)
            k <- k + 1

            if (k > length(chr.seg.names.vec)) {

                if (shouldSegsBeJoined(s.curr.df, s.next.df, pval.to.separate) == TRUE) {
                    new.seg.name <- paste(unlist(strsplit(s.curr.name,split=":"))[1],unlist(strsplit(s.curr.name,split=":"))[2],unlist(strsplit(s.next.name,split=":"))[3],sep=":")
                    s.curr.df <- rbind(s.curr.df, s.next.df)
                    samp.cpnX.segs.list <- updateJoinedSegList(samp, cpnX.string, samp.cpnX.joined.segs.list=samp.cpnX.segs.list, new.seg.name, s.curr.df)
                    #Get next seg to consider
                    s.curr.name <- s.next.name
                    s.next.name <- chr.seg.names.vec[k]
                    s.next.df <- getSegForJoiningDF(cpnX.string, samp, seg.name=s.next.name, all.pw.cn.het.comp.lists, all.seg.CCF.vals.list)
                } else {#end if (shouldSegsBeJoined(s.curr.df, s.next.df) == TRUE)
                    #Keep the two segs under consideration and prepare to exit the while
                    samp.cpnX.segs.list <- updateJoinedSegList(samp, cpnX.string, samp.cpnX.joined.segs.list=samp.cpnX.segs.list, s.curr.name, s.curr.df)
                    samp.cpnX.segs.list <- updateJoinedSegList(samp, cpnX.string, samp.cpnX.joined.segs.list=samp.cpnX.segs.list, s.next.name, s.next.df)
                }#end if (shouldSegsBeJoined(s.curr.df, s.next.df) == FALSE)
            }#end if (k > length(chr.seg.names.vec)) {
        }#end if (shouldSegsBeJoined(s.curr.df, s.next.df) == TRUE) 

        if (shouldSegsBeJoined(s.curr.df, s.next.df, pval.to.separate) == FALSE & k <= length(chr.seg.names.vec)) {
            #keep seg as is and move to the next one for comparisons
            samp.cpnX.segs.list <- updateJoinedSegList(samp, cpnX.string, samp.cpnX.joined.segs.list=samp.cpnX.segs.list, s.curr.name, s.curr.df)
            #Update current seg and get next one
            s.curr.name <- s.next.name
            s.curr.df <- s.next.df
            s.next.name <- chr.seg.names.vec[k]
            s.next.df <- getSegForJoiningDF(cpnX.string, samp, seg.name=s.next.name, all.pw.cn.het.comp.lists, all.seg.CCF.vals.list)
            k <- k + 1

            if (k > length(chr.seg.names.vec)) {

                if (shouldSegsBeJoined(s.curr.df, s.next.df, pval.to.separate) == TRUE) {
                    #segs can be joined
                    new.seg.name <- paste(unlist(strsplit(s.curr.name,split=":"))[1],unlist(strsplit(s.curr.name,split=":"))[2],unlist(strsplit(s.next.name,split=":"))[3],sep=":")
                    s.curr.df <- rbind(s.curr.df, s.next.df)
                    samp.cpnX.segs.list <- updateJoinedSegList(samp, cpnX.string, samp.cpnX.joined.segs.list=samp.cpnX.segs.list, new.seg.name, s.curr.df)
                    s.curr.name <- s.next.name
                    s.next.name <- chr.seg.names.vec[k]
                    s.next.df <- getSegForJoiningDF(cpnX.string, samp, seg.name=s.next.name, all.pw.cn.het.comp.lists, all.seg.CCF.vals.list)
                } else {
                  
                    samp.cpnX.segs.list <- updateJoinedSegList(samp, cpnX.string, samp.cpnX.joined.segs.list=samp.cpnX.segs.list, s.curr.name, s.curr.df)
                    samp.cpnX.segs.list <- updateJoinedSegList(samp, cpnX.string, samp.cpnX.joined.segs.list=samp.cpnX.segs.list, s.next.name, s.next.df)
                }#end  if (shouldSegsBeJoined(s.curr.df, s.next.df, pval.to.separate) == TRUE) {
            }#end  if (k > length(chr.seg.names.vec)) {
        }#end if (shouldSegsBeJoined(s.curr.df, s.next.df) == FALSE & k <= length(chr.seg.names.vec)) {
    }# end while (k <= length(chr.seg.names.vec)) {
    return(samp.cpnX.segs.list)
}#end joinAdjacentSegsWithinSampleChromosome