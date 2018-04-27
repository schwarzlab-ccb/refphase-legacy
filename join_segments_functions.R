
##Params:
##
##Returns:
##
getSegForJoiningDF <- function(cpnX.string, samp, seg.name, all.pw.cn.het.comp.lists, all.seg.CCF.vals.list) {

    seg.to.join.df <- data.frame(cpnX = all.pw.cn.het.comp.lists[[samp]][[seg.name]][[paste0("ph_", cpnX.string, "_vec")]],
               chr = all.pw.cn.het.comp.lists[[samp]][[seg.name]]$chr_vec,
               pos = all.pw.cn.het.comp.lists[[samp]][[seg.name]]$pos_vec,
               low_cpnX_frac = all.seg.CCF.vals.list[[seg.name]][[samp]][[paste0("frac_", cpnX.string, "_low_vec")]],
               high_cpnX_frac = all.seg.CCF.vals.list[[seg.name]][[samp]][[paste0("frac_", cpnX.string, "_high_vec")]], 
               cpnX_frac = all.pw.cn.het.comp.lists[[seg.name]][[samp]][[paste0("ph_", cpnX.string, "_vec")]],
               low_sc_cpnX = rep(all.seg.CCF.vals.list[[seg.name]][[paste0("poss_sc_", cpnX.string, "_low")]]),
               high_sc_cpnX = rep(all.seg.CCF.vals.list[[seg.name]][[paste0("poss_sc_", cpnX.string, "_high")]]),
               stringsAsFactors = FALSE)

    return(seg.join.df)
}

updateJoinedSegList <- function(cpnX.string, samp.cpnX.joined.segs.list, seg.name, seg.df)

    samp.cpnX.joined.segs.list[[seg.name]][[paste0("ph_", cpnX.string, "_vec")]] <- seg.df$cpnX
    samp.cpnX.joined.segs.list[[seg.name]][["chr_vec"]] <- seg.df$chr
    samp.cpnX.joined.segs.list[[seg.name]][["pos_vec"]] <- seg.df$pos
    samp.cpnX.joined.segs.list[[seg.name]][[paste0("frac_", cpnX.string, "_low_vec")]] <- seg.df$low_cpnX_frac
    samp.cpnX.joined.segs.list[[seg.name]][[paste0("frac_", cpnX.string, "_high_vec")]] <- seg.df$high_cpnX_frac
    samp.cpnX.joined.segs.list[[seg.name]][[paste0(cpnX.string, "_frac")]] <- seg.df$cpnX_frac
    samp.cpnX.joined.segs.list[[seg.name]][[paste0("poss_sc_", cpnX.string, "_low")]] <- seg.df$low_sc_cpnX
    samp.cpnX.joined.segs.list[[seg.name]][[paste0("poss_sc_", cpnX.string, "_high")]] <- seg.df$high_sc_cpnX

    return(samp.cpnX.joined.segs.list)

}#end updateJoinedSegList

wilcoxTestWithErrorHandling <- function(x,y) {
    out <- wilcox.test(x,y)
    if(is.na(out$p.value))
    {
      out$p.value <- 1
    }
    return(out)
}#end wilcoxTestWithErrorHandling
  

shouldSegsBeJoined <- function(s.curr.df, s.next.df){

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
joinAdjacentSegsWithinSampleChromosome <- function(chr.seg.df, pval.to.separate) {

    #Initialise variables:
    segs.cpnX.used.vec <- c()
    samp.cpnX.segs.list <- list()

    chr.seg.names.vec <- rownames(chr.seg.df)
    s.curr.name <-chr.seg.names.vec[1] 
    s.curr.df <- getSegForJoiningDF(cpnX.string, samp, seg.name = s.curr.name, all.pw.cn.het.comp.lists, all.seg.CCF.vals.list)
    #Call 
    if (length(chr.seg.names.vec) == 1) {

        samp.cpnX.segs.list <- updateJoinedSegList(cpnX.string, samp.cpnX.joined.segs.list=samp.cpnX.segs.list, s.curr.name, s.curr.df)
        segs.cpnX.used.vec <- c(segs.cpnX.used.vec, seg.curr.name)

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
                samp.cpnX.segs.list <- updateJoinedSegList(cpnX.string, samp.cpnX.joined.segs.list=samp.cpnX.segs.list, new.seg.name, s.curr.df)
                segs.cpnX.used.vec <- c(segs.cpnX.used.vec, new.seg.name)

            } else {

                samp.cpnX.segs.list <- updateJoinedSegList(cpnX.string, samp.cpnX.joined.segs.list=samp.cpnX.segs.list, s.curr.name, s.curr.df)
                samp.cpnX.segs.list <- updateJoinedSegList(cpnX.string, samp.cpnX.joined.segs.list=samp.cpnX.segs.list, s.next.name, s.next.df)
                segs.cpnX.used.vec <- c(segs.cpnX.used.vec, s.curr.name, s.next.name)

            }#end if (cpnX.vec.s.one.vs.s.two.test$p.val <= pval.to.separate &
    }# end if (length(chr.seg.names.vec) == 2) {

    k <- 3

    while (k <= length(chr.seg.names.vec)) {

        if (shouldSegsBeJoined(s.curr.df, s.next.df, pval.to.separate) == TRUE) {

            new.seg.name <- paste(unlist(strsplit(s.curr.name,split=":"))[1],unlist(strsplit(s.curr.name,split=":"))[2],unlist(strsplit(s.next.name,split=":"))[3],sep=":")
            s.curr.df <- rbind(s.curr.df, s.next.df)
            samp.cpnX.segs.list <- updateJoinedSegList(cpnX.string, samp.cpnX.joined.segs.list=samp.cpnX.segs.list, new.seg.name, s.curr.df)
            segs.cpnX.used.vec <- c(segs.cpnX.used.vec, new.seg.name)
            #Get next seg to consider
            s.next.name <- chr.seg.names.vec[k]
            s.next.df <- getSegForJoiningDF(cpnX.string, samp, seg.name=s.next.name, all.pw.cn.het.comp.lists, all.seg.CCF.vals.list)
            k <- k + 1

            if (k > length(chr.seg.names.vec)) {

                if (shouldSegsBeJoined(s.curr.df, s.next.df, pval.to.separate) == TRUE) {

                    new.seg.name <- paste(unlist(strsplit(s.curr.name,split=":"))[1],unlist(strsplit(s.curr.name,split=":"))[2],unlist(strsplit(s.next.name,split=":"))[3],sep=":")
                    s.curr.df <- rbind(s.curr.df, s.next.df)
                    samp.cpnX.segs.list <- updateJoinedSegList(cpnX.string, samp.cpnX.joined.segs.list=samp.cpnX.segs.list, new.seg.name, s.curr.df)
                    segs.cpnX.used.vec <- c(segs.cpnX.used.vec, new.seg.name)
                    #Get next seg to consider
                    s.next.name <- chr.seg.names.vec[k]
                    s.next.df <- getSegForJoiningDF(cpnX.string, samp, seg.name=s.next.name, all.pw.cn.het.comp.lists, all.seg.CCF.vals.list)

                } else {#end if (shouldSegsBeJoined(s.curr.df, s.next.df) == TRUE)
                    #Keep the two segs under consideration and prepare to exit the while
                    samp.cpnX.segs.list <- updateJoinedSegList(cpnX.string, samp.cpnX.joined.segs.list=samp.cpnX.segs.list, s.curr.name, s.curr.df)
                    samp.cpnX.segs.list <- updateJoinedSegList(cpnX.string, samp.cpnX.joined.segs.list=samp.cpnX.segs.list, s.next.name, s.next.df)
                    segs.cpnX.used.vec <- c(segs.cpnX.used.vec, s.curr.name, s.next.name)

                }#end if (shouldSegsBeJoined(s.curr.df, s.next.df) == FALSE)
            }#end if (k > length(chr.seg.names.vec)) {
        }#end if (shouldSegsBeJoined(s.curr.df, s.next.df) == TRUE) 

        if (shouldSegsBeJoined(s.curr.df, s.next.df, pval.to.separate) == FALSE & k <= length(chr.seg.names.vec)) {

            #keep seg as is and move to the next one for comparisons
            samp.cpnX.segs.list <- updateJoinedSegList(cpnX.string, samp.cpnX.joined.segs.list=samp.cpnX.segs.list, s.curr.name, s.curr.df)
            segs.cpnX.used.vec <- c(segs.cpnX.used.vec, s.curr.name)
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
                    samp.cpnX.segs.list <- updateJoinedSegList(cpnX.string, samp.cpnX.joined.segs.list=samp.cpnX.segs.list, new.seg.name, s.curr.df)
                    segs.cpnX.used.vec <- c(segs.cpnX.used.vec, new.seg.name)
                    s.next.name <- chr.seg.names.vec[k]
                    s.next.df <- getSegForJoiningDF(cpnX.string, samp, seg.name=s.next.name, all.pw.cn.het.comp.lists, all.seg.CCF.vals.list)

                } else {

                    samp.cpnX.segs.list <- updateJoinedSegList(cpnX.string, samp.cpnX.joined.segs.list=samp.cpnX.segs.list, s.curr.name, s.curr.df)
                    samp.cpnX.segs.list <- updateJoinedSegList(cpnX.string, samp.cpnX.joined.segs.list=samp.cpnX.segs.list, s.next.name, s.next.df)
                    segs.cpnX.used.vec <- c(segs.cpnX.used.vec, s.curr.name, s.next.name)

                }#end  if (shouldSegsBeJoined(s.curr.df, s.next.df, pval.to.separate) == TRUE) {
            }#end  if (k > length(chr.seg.names.vec)) {
        }#end if (shouldSegsBeJoined(s.curr.df, s.next.df) == FALSE & k <= length(chr.seg.names.vec)) {
    }# end while (k <= length(chr.seg.names.vec)) {

    return(object)
}#end joinAdjacentSegsWithinSampleChromosome


##
##Params:
##
##Returns:
##
joinAdjacentSegsWithCCF <- function(all.pw.cn.het.comp.lists, all.seg.CCF.vals.list, pval.to.separate, cpnX.string){

    #Pull out all pairwise comparisons
    sample.vec <- grep(":", names(all.pw.cn.het.comp.lists), invert=TRUE, value=TRUE)
    seg.name.vec <- names(all.pw.cn.het.comp.lists[[sample.vec[1]]])
    new.segs.list <- list()
    filt.segs.vec <- c()

    for ( seg.name in seg.name.vec ) {

        if (is.na(all.pw.cn.het.comp.lists[[sample.vec[1]]][[seg.name]]$ph_cpnB_vec)) {

            filt.segs.vec <- c(filt.segs.vec, seg.name)

        }#end if (is.na(all.pw.cn.het.comp.lists[[sample.vec[1]]][[seg.name]]$ph_cpnB_vec)) {s

    }#end for ( seg.name in seg.name.vec ) {

    #So there should be a function call 
    samp.one <- sample.vec[1]

    for ( samp in sample.vec ) {

        all.segs.chr.pos.df <- do.call(rbind, strsplit(filt.segs.vec, split=":"))
        all.segs.chr.pos.df <- apply(all.segs.chr.pos.df, 2, as.numeric)
        rownames(all.segs.chr.pos.df) <- filt.segs.vec

        samp.cpnX.segs.list <- list()

        for ( chr in all.segs.chr.pos.df[,1] ) {

            chr.seg.df <- all.segs.chr.pos.df[all.segs.chr.pos.df[,1]%in%chr, ,drop = FALSE]

            joinAdjacentSegsWithinSampleChromosome <- 
            #We can insert logic here to do chromosome arms if we want

        }#end or ( chrs in all.segs.chr.pos.df[,1] ) {

    }#end for ( samp in sample.vec ) {

    return(object)
}
