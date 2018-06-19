




generateMinimumConsistentSegmentation <- function(orig.cpn.segs.df, chroms.vec = 1:22, mc.cores = 2, col.to.use = "cpn_to_use") {
  require(parallel)
  orig.cpn.segs.mat <- as.matrix(orig.cpn.segs.df)

  orig.cpn.segs.mat <- apply(orig.cpn.segs.mat,2,function(x)gsub('\\s+', '',x))#remove any lingering whitespace from matrix conversion
  all.new.seg.chr.vec <- c()
  all.new.seg.start.vec <- c()
  all.new.seg.end.vec <- c()


  for ( i in 1:length(chroms.vec) ) {

    chr.cpn.segs.mat <- orig.cpn.segs.mat[as.numeric(orig.cpn.segs.mat[,2]) == as.numeric(chroms.vec[i]),,drop=FALSE]
    chr.seg.start.vec <- unique(as.numeric(chr.cpn.segs.mat[,3]))
    chr.seg.end.vec <- unique(as.numeric(chr.cpn.segs.mat[,4]))
    chr.seg.start.vec <- chr.seg.start.vec[ ! chr.seg.start.vec %in% chr.seg.end.vec]#Added this to remove Inf ends sometimes resulted
    joint.chr.start.vec <- c()
    joint.chr.end.vec <- c()
    current.end <- 0
    end.next <- 0
    while (current.end < max(chr.seg.start.vec) & 
            length(chr.seg.start.vec[chr.seg.start.vec > current.end]) > 0) {
      ##Added this condition to remove warning messages from a numeric(0) passed to min in the loop comparison, replaced this:
                #& max(chr.seg.end.vec) > min(chr.seg.start.vec[chr.seg.start.vec > current.end])#potential bug with v low end positions? test later.
      current.start <- as.numeric(min(chr.seg.start.vec[chr.seg.start.vec > current.end]))
      current.end <- as.numeric(min(chr.seg.end.vec[chr.seg.end.vec > current.start]))
      joint.chr.start.vec <- c(joint.chr.start.vec, current.start)
      joint.chr.end.vec <- c(joint.chr.end.vec, current.end)

    } #end while(current.end < max....)

    all.new.seg.chr.vec <- c(all.new.seg.chr.vec, rep(chroms.vec[i], length(joint.chr.start.vec)))
    all.new.seg.start.vec <- c(all.new.seg.start.vec, joint.chr.start.vec)
    all.new.seg.end.vec <- c(all.new.seg.end.vec, joint.chr.end.vec)

  } #end for ( i in 1:length(chroms.vec) )

  chr.pos.mat <- matrix(NA, length(all.new.seg.start.vec), ncol = 3)
  colnames(chr.pos.mat) <- c("chr", "start", "end")
  chr.pos.mat <- cbind(as.numeric(all.new.seg.chr.vec), as.numeric(all.new.seg.start.vec), as.numeric(all.new.seg.end.vec))


  joint.seg.cpn.chr.mat.list <- mclapply(unique(chr.pos.mat[,1]), function(x){

      t(apply(chr.pos.mat[chr.pos.mat[,1] == x,, drop=FALSE], 1, mapOriginalCopyNumberToJointSegmentation, orig.cpn.segs.mat = orig.cpn.segs.mat, col.to.use=col.to.use))

    }#end function within mclapply..
  )#end mclapply

  joint.seg.cpn.mat <- do.call(rbind, joint.seg.cpn.chr.mat.list)
  final.joint.seg.mat <- cbind(chr.pos.mat, joint.seg.cpn.mat)
  return(final.joint.seg.mat)

} #end generateMinimumConsistentSegmentation <- function(...


#Todo could fail if there are segments present in some but not others...
#Check for this.
mapOriginalCopyNumberToJointSegmentation <- function(chr.start.end.vec, orig.cpn.segs.mat, col.to.use) {

  joint.seg.cpn.by.samp.vec <- rep(NA, length(unique(orig.cpn.segs.mat[, 1])))
  names(joint.seg.cpn.by.samp.vec) <- unique(orig.cpn.segs.mat[, 1])
  orig.cpn.chr.segs.mat <- orig.cpn.segs.mat[as.numeric(orig.cpn.segs.mat[,2]) == chr.start.end.vec[1],]
  red.cpn.chr.segs.mat <- cbind(orig.cpn.chr.segs.mat[,c(2:5)], orig.cpn.chr.segs.mat[, col.to.use])
                                                             #orig start                          #joint end                     #orig end                          #joint end
  overlapping.orig.cpn.chr.seg.idxs <- which(( (as.numeric(red.cpn.chr.segs.mat[, 2]) <= chr.start.end.vec[3] &  as.numeric(red.cpn.chr.segs.mat[, 3]) >= chr.start.end.vec[3]) #orig.start <= joint.end & orig.end >= joint.end#PICK UP SEGS OVERLAPPING TO RIGHT BOUNDARY
                            #orig start                          #joint start                    #orig end                         #joint start
               |   (as.numeric(red.cpn.chr.segs.mat[, 2]) <= chr.start.end.vec[2] &  as.numeric(red.cpn.chr.segs.mat[, 3]) >= chr.start.end.vec[2]) #orig.start <= joint.start & orig.end >= joint.start#PICK UP SEGS OVERLAPPING LEFT BOUNDARY
                            #orig start                           #joint start                   #orig end                         #joint end
               |   (as.numeric(red.cpn.chr.segs.mat[, 2]) >= chr.start.end.vec[2] &  as.numeric(red.cpn.chr.segs.mat[, 3]) <= chr.start.end.vec[3]  ) #orig.start >= joint.start & orig.end <= joint.start#SEGS W/IN (So equal to) JOINT segs
  ))
  
  overlap.orig.chr.segs.w.joint.mat <- orig.cpn.chr.segs.mat[overlapping.orig.cpn.chr.seg.idxs,, drop=FALSE]
  #Gets the all indices of the run of rows of the same sample that overlapped in the original segmentation with the joint segment 
  sample.run.idxs.logical.vec <- duplicated(overlap.orig.chr.segs.w.joint.mat[, 1], fromLast=TRUE) == TRUE | duplicated(overlap.orig.chr.segs.w.joint.mat[, 1]) == TRUE
  #Get rows of samples that don't have multiple overlapping segments with this joint segment
  singleton.samp.overlap.chr.segs.w.joint.mat <- overlap.orig.chr.segs.w.joint.mat[! sample.run.idxs.logical.vec, ,drop=FALSE]
  #Record the copy number values of the samples wthich mtacht single segment and that don't require extra work in the vector
  joint.seg.cpn.by.samp.vec[ singleton.samp.overlap.chr.segs.w.joint.mat[, 1] ] <- as.numeric(singleton.samp.overlap.chr.segs.w.joint.mat[, col.to.use])
  #For joint segs where there are multiple matching cpn segments from the same sample in the original dataframe
  if (length(which(sample.run.idxs.logical.vec)) >= 1) {

    #Pull out the only the samples with multiple matching original segments to this joint segment
    samples.w.multi.overlap.w.joint.mat <- overlap.orig.chr.segs.w.joint.mat[ sample.run.idxs.logical.vec , ]
    samples.w.m.overlap.cpn.val.vec <- rep(NA, length(unique(samples.w.multi.overlap.w.joint.mat[,1])))
    names(samples.w.m.overlap.cpn.val.vec) <- unique(samples.w.multi.overlap.w.joint.mat[,1])

    for (sample in unique(samples.w.multi.overlap.w.joint.mat[,1])) {

      matching.sample.idxs <- which(samples.w.multi.overlap.w.joint.mat[,1] == sample)
      seg.overlap.sizes <- c()

      for (idx in matching.sample.idxs){

          seg.overlap.sizes  <- c(seg.overlap.sizes, pmax(c(chr.start.end.vec[2], chr.start.end.vec[3])) - 
                                  pmin(c(as.numeric(samples.w.multi.overlap.w.joint.mat[idx, 3]), as.numeric(samples.w.multi.overlap.w.joint.mat[idx,4]))))

      }#end for(idx in matching.sample.idxs)

      samples.w.m.overlap.cpn.val.vec[sample] <- samples.w.multi.overlap.w.joint.mat[matching.sample.idxs, col.to.use][which(seg.overlap.sizes==max(seg.overlap.sizes))[1]]

    }#end for(sample in unique(samples.w.multi.overlap.w.joint.mat[,1])) {

    joint.seg.cpn.by.samp.vec[names(samples.w.m.overlap.cpn.val.vec)]  <- as.numeric(samples.w.m.overlap.cpn.val.vec)

  }#end if (length(which(sample.run.idxs.logical.vec)) > 1)...

  return(joint.seg.cpn.by.samp.vec)               

} #end mapOriginalCopyNumberToJointSegmentation <- function(...



#We will now make minimum consistent segments that will cover both arms of the chromosomes:
##
##Params:
##
##Returns:
##Todo make version of the genome that's used a input option, for now stick with the hg19
getChromosomeArmsGRangeList <- function(current.chr, chrom.info.df) {

  arm.df.list <- list() 
  p.arm.df <- data.frame(
  chrom  = current.chr,
  start = chrom.info.df[chrom.info.df$chrom == current.chr,]$p_start,
  end = chrom.info.df[chrom.info.df$chrom == current.chr,]$p_end,
  stringsAsFactors = FALSE)

  count <- 1

  if(chrom.info.df[chrom.info.df$chrom == current.chr,]$acrocentric == FALSE) {

    arm.df.list[[count]] <- makeGRangesFromDataFrame(p.arm.df,
                             keep.extra.columns=FALSE,
                             ignore.strand=FALSE,
                             seqinfo=NULL,
                             seqnames.field=c("seqnames", "seqname",
                                              "chromosome", "chrom",
                                              "chr", "chromosome_name",
                                              "seqid"),
                             start.field=c("Start_Position","start"),
                             end.field=c("End_Position", "end", "stop"),
                             strand.field="strand",
                             starts.in.df.are.0based=FALSE)

    count <- count + 1
  }



  q.arm.df <- data.frame(
  chrom  = current.chr,
  start = chrom.info.df[chrom.info.df$chrom == current.chr,]$q_start,
  end = chrom.info.df[chrom.info.df$chrom == current.chr,]$q_end,
  stringsAsFactors = FALSE)

  arm.df.list[[count]] <- makeGRangesFromDataFrame(q.arm.df,
                           keep.extra.columns=FALSE,
                           ignore.strand=FALSE,
                           seqinfo=NULL,
                           seqnames.field=c("seqnames", "seqname",
                                            "chromosome", "chrom",
                                            "chr", "chromosome_name",
                                            "seqid"),
                           start.field=c("Start_Position","start"),
                           end.field=c("End_Position", "end", "stop"),
                           strand.field="strand",
                           starts.in.df.are.0based=FALSE)



  return(arm.df.list)
}
