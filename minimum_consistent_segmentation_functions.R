

reformatInitialSegmentationInput <- function(f.joint.segs.df) {

  f.joint.segs.df$sample <- gsub("\\.", "-", gsub("_LogR", "", f.joint.segs.df$sample))
  f.joint.segs.df$callState <- NULL
  f.joint.segs.df$patient <- patient
  f.joint.segs.df$sample_name <- paste(patient,sapply(f.joint.segs.df$sample, function(x) strsplit(x, split = "-")[[1]][2]), sep = ":")
  cpn.call.state.vec <- rep("neutral", nrow(f.joint.segs.df))
  tot.raw.cpn.vec <- f.joint.segs.df$nAraw + f.joint.segs.df$nBraw
  log.tot.raw.over.ploidy.cpn.vec <- log2(tot.raw.cpn.vec / f.joint.segs.df$Ploidy)
  cpn.call.state.vec[log.tot.raw.over.ploidy.cpn.vec >= log2(2.5 / 2)] <- "gain"
  cpn.call.state.vec[log.tot.raw.over.ploidy.cpn.vec <= log2(1.5 / 2)] <- "loss"
  cpn.call.state.vec[tot.raw.cpn.vec >= ((2 * f.joint.segs.df$Ploidy) + 1)] <- "amp"
  cpn.call.state.vec[tot.raw.cpn.vec <= 0.5] <- "deep_loss"
  f.joint.segs.df$cpn_call <- cpn.call.state.vec
  f.joint.segs.df$arm_or_focal_event <- "none"
  f.joint.segs.df$arm_event_present <- "none"
  f.joint.segs.df$within_chr_arm <- "not_in_arm"
  f.joint.segs.df$uniq_id <- paste0(f.joint.segs.df$sample, "@", f.joint.segs.df$chr, ":", f.joint.segs.df$startpos, ":", f.joint.segs.df$endpos)
  rownames(f.joint.segs.df) <- 1:nrow(f.joint.segs.df)
  return(f.joint.segs.df)

}#end reformatInitialSegmentationInput <- function(f.joint.segs.df) {



#TBKW As of 20190818 this function is used 
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







#Created as of 20190818 to allow multiple  calls of the segmentation
#It absolutely needs breaking up into multiple smaller functions
#And various design decisions regarding whether we should maintain the capability of covering
#the entire genome in copy number segments should be mulled over, currently, as they are 
#are clearly marked and easily dropped I would be in favour of keeping them.
createFullyAnnotatedJointSegmentation <- function(f.joint.segs.df) {



  all.samp.grs <- makeGRangesFromDataFrame(f.joint.segs.df,
                               keep.extra.columns=FALSE,
                               ignore.strand=FALSE,
                               seqinfo=NULL,
                               seqnames.field=c("seqnames", "seqname",
                                                "chromosome", "chrom",
                                                "chr", "chromosome_name",
                                                "seqid"),
                               start.field=c("Start_Position","start", "startpos"),
                               end.field=c("End_Position", "end", "stop", "endpos"),
                               strand.field="strand",
                               starts.in.df.are.0based=FALSE)

  all.samp.w.metadata.grs <- makeGRangesFromDataFrame(f.joint.segs.df,
                          keep.extra.columns=TRUE,
                          ignore.strand=FALSE,
                          seqinfo=NULL,
                          seqnames.field=c("seqnames", "seqname",
                                           "chromosome", "chrom",
                                           "chr", "chromosome_name",
                                           "seqid"),
                          start.field=c("Start_Position","start","startpos"),
                          end.field=c("End_Position", "end", "stop", "endpos"),
                          strand.field="strand",
                          starts.in.df.are.0based=FALSE)

  print("Introducing breakpoints at starts and ends of centromeres in original segmentation to allow calculation of arm level events in the original...")
  new.all.samp.w.metadata.grs <- sort(all.samp.w.metadata.grs)
  for (chr in CHRS_TO_USE_VEC) {

    p.end.gr <- makeGRangesFromDataFrame(chrom.info.df[chrom.info.df$chrom == chr,],
                          keep.extra.columns=FALSE,
                          ignore.strand=FALSE,
                          seqinfo=NULL,
                          seqnames.field=c("chrom"),
                          start.field=c("p_end"),
                          end.field=c("p_end"),
                          strand.field="strand",
                          starts.in.df.are.0based=FALSE)
    #print(p.end.gr)
    q.start.gr <- makeGRangesFromDataFrame(chrom.info.df[chrom.info.df$chrom == chr,],
                            keep.extra.columns=FALSE,
                            ignore.strand=FALSE,
                            seqinfo=NULL,
                            seqnames.field=c("chrom"),
                            start.field=c("q_start"),
                            end.field=c("q_start"),
                            strand.field="strand",
                            starts.in.df.are.0based=FALSE)
    #print(q.start.gr)

    bps.to.introduce.grs.list <- list(p.end.gr, q.start.gr)
    for (new.bp.gr in bps.to.introduce.grs.list) {

      #DEBUG
      #new.bp.gr <- bps.to.introduce.grs.list[[i]]

      overlapping.seg.grs <- getOverlapGRanges(new.bp.gr, new.all.samp.w.metadata.grs)
      #So go through each one, split into two,add to the dataset, remove the original, sort
      if(length(overlapping.seg.grs) > 0) {

        for ( i in 1:length(overlapping.seg.grs) ) {

          new.seg.one.template <- new.all.samp.w.metadata.grs[overlapping.seg.grs[i]$subject_hit_idx]
          new.seg.one.start <- overlapping.seg.grs[i]$subject_start
          new.seg.one.end <- start(new.bp.gr)

          start(new.seg.one.template) <- new.seg.one.start
          end(new.seg.one.template) <- new.seg.one.end
          new.seg.one.template$uniq_id <- paste0(new.seg.one.template$sample,
                              "@",
                              as.character(seqnames(new.seg.one.template)),
                              ":",
                              start(new.seg.one.template),
                              ":",
                              end(new.seg.one.template))

          new.seg.two.template <- new.all.samp.w.metadata.grs[overlapping.seg.grs[i]$subject_hit_idx]
          new.seg.two.start  <- start(new.bp.gr) + 1
          new.seg.two.end <- overlapping.seg.grs[i]$subject_end

          start(new.seg.two.template) <- new.seg.two.start
          end(new.seg.two.template) <- new.seg.two.end
          new.seg.two.template$uniq_id <- paste0(new.seg.two.template$sample,
                              "@",
                              as.character(seqnames(new.seg.two.template)),
                              ":",
                              start(new.seg.two.template),
                              ":",
                              end(new.seg.two.template))
          #Remove the original, so get the subject idx
          all.orig.idx.vec <- 1:length(new.all.samp.w.metadata.grs)
          idx.of.orig.segment <- overlapping.seg.grs[i]$subject_hit_idx
          new.all.samp.w.metadata.grs <- new.all.samp.w.metadata.grs[ ! all.orig.idx.vec %in% idx.of.orig.segment ]

          #So now add the new segs to the dataset and sort:
          new.all.samp.w.metadata.grs <- sort(c(new.all.samp.w.metadata.grs, new.seg.one.template, new.seg.two.template))
        }#end for ( i in 1:length(overlapping.seg.grs) ) {
      }#end if(length(overlapping.seg.grs) > 0) {
    }#end     for (new.bp.gr in bps.to.introduce.grs.list) {
  }#end for ( chr in CHRS_TO_USE_VEC ) {


  #So actually defining the segments that are lost or gained makes sense here..
  print("Calculating arm level gain, loss, neutral, and LOH events in the original segmentation...")
  all.samp.w.metadata.grs <- new.all.samp.w.metadata.grs
  gain.loss.name.list <- list(c("neutral"), c("gain", "amp"), c("loss", "deep_loss"),  c("LOH"))
  names(gain.loss.name.list) <- c("neutral", "gain", "loss", "LOH")
  p.or.q.vec <- c("p", "q")
  for ( chr in CHRS_TO_USE_VEC ) {
    chr.gr.list <- getChromosomeArmsGRangeList(chr, chrom.info.df)
    for ( i in 1:length(chr.gr.list) ) {

      curr.chr.arm.gr <- chr.gr.list[[i]]
      all.samp.overlap.grs <- getOverlapGRanges(curr.chr.arm.gr, all.samp.w.metadata.grs)

      for (samp in sample.vec) {

        arm.level.events.present <- list()
        for (j in 1:length(gain.loss.name.list)) {
          cpn.event.vec <- gain.loss.name.list[[j]]
          #cpn.event.vec <- c("loss", "deep_loss")
          #First does the total area of the event add up to 80% of an arm:
          if ( ! names(gain.loss.name.list)[j] %in% "LOH") {
            samp.gain.amp.grs <- all.samp.overlap.grs[ all.samp.overlap.grs$cpn_call %in% cpn.event.vec &
                                 all.samp.overlap.grs$sample == samp
                                ]

            samp.all.cpn.grs <- all.samp.overlap.grs[ all.samp.overlap.grs$sample == samp ]

          } else {

            samp.gain.amp.grs <- all.samp.overlap.grs[ all.samp.overlap.grs$nMinor == 0 &
                                   all.samp.overlap.grs$sample == samp
                                  ]

          }

          samp.arm.all.uniq.ids.vec <- samp.all.cpn.grs$uniq_id
          if( length(all.samp.w.metadata.grs[all.samp.w.metadata.grs$uniq_id %in% samp.arm.all.uniq.ids.vec]) == 0) {

            next

          }

          all.samp.w.metadata.grs[all.samp.w.metadata.grs$uniq_id %in% samp.arm.all.uniq.ids.vec]$within_chr_arm <- paste0(chr, p.or.q.vec[i])
          if (sum(as.numeric(samp.gain.amp.grs$overlap_prop_query))/100 >= 0.74 & sum(samp.gain.amp.grs[ samp.gain.amp.grs$overlap_prop_query >= 0.2 ]$overlap_prop_query)/100 >= 0.5){

              #So all segments in the arm should be annotated as being in an arm level event:
              #All the larger segments
              arm.level.events.present[[names(gain.loss.name.list)[j]]] <- list()
              #So now

              arm.level.events.present[[names(gain.loss.name.list)[j]]]$arm_event_uniq_ids <- samp.gain.amp.grs[ samp.gain.amp.grs$overlap_prop_query >= 0.2 & samp.gain.amp.grs$cpn_call %in% cpn.event.vec ]$uniq_id
              arm.level.events.present[[names(gain.loss.name.list)[j]]]$focal_event_uniq_ids <- samp.all.cpn.grs[ (samp.all.cpn.grs$overlap_prop_query < 0.2 & samp.all.cpn.grs$cpn_call %in% cpn.event.vec)
                                    | (! (samp.all.cpn.grs$cpn_call %in% cpn.event.vec)) ]$uniq_id
              #so anything that isn't an arm and was gained will be focal
            }

          }
          if (length(arm.level.events.present) >= 1) {

            if (names(arm.level.events.present)[1] != "LOH") {

              all.samp.w.metadata.grs[all.samp.w.metadata.grs$uniq_id %in% arm.level.events.present[[1]]$arm_event_uniq_ids ]$arm_or_focal_event  <- "arm"
              if (length(arm.level.events.present[[1]]$focal_event_uniq_ids) > 0) {

                all.samp.w.metadata.grs[all.samp.w.metadata.grs$uniq_id %in% arm.level.events.present[[1]]$focal_event_uniq_ids  ]$arm_or_focal_event  <- "focal"

              }

            }
            all.samp.w.metadata.grs[all.samp.w.metadata.grs$uniq_id %in% samp.arm.all.uniq.ids.vec ]$arm_event_present <- paste(names(arm.level.events.present), collapse="_")

          } else {

            if (length(all.samp.w.metadata.grs[all.samp.w.metadata.grs$uniq_id %in% samp.arm.all.uniq.ids.vec ]) > 0) {

              all.samp.w.metadata.grs[all.samp.w.metadata.grs$uniq_id %in% samp.arm.all.uniq.ids.vec ]$arm_or_focal_event <- "focal"

            }
          }
          #Next is there one segment with at least 50% of the arm in it:
        }#end for (cpn.event.vec in gain.loss.name.list) {
      }#end for (samp in sample.vec) {
  }#end for ( chr in CHRS_TO_USE_VEC ) {

  all.samp.w.metadata.df <- as.data.frame(all.samp.w.metadata.grs)
  all.samp.w.metadata.df <- all.samp.w.metadata.df[with(all.samp.w.metadata.df, order(sample, seqnames, start)), ]
  #print(all.samp.w.metadata.df)
  #Let's generate the minimum consistent segmentation now:

  all.samp.w.metadata.conv.df <- all.samp.w.metadata.df[, grep("strand|width", colnames(all.samp.w.metadata.df), value = TRUE, invert=TRUE)]
  colnames(all.samp.w.metadata.conv.df)[colnames(all.samp.w.metadata.conv.df)=="seqnames"] <- "chr"
  colnames(all.samp.w.metadata.conv.df)[colnames(all.samp.w.metadata.conv.df)=="start"] <- "startpos"
  colnames(all.samp.w.metadata.conv.df)[colnames(all.samp.w.metadata.conv.df)=="end"] <- "endpos"
  #Reorder
  correct.order.cols <- c( "sample","chr","startpos","endpos","n.het","cnTotal","nMajor","nMinor","Ploidy","ACF","nAraw","nBraw","patient","sample_name","cpn_call","arm_or_focal_event","arm_event_present","within_chr_arm","uniq_id")
  f.joint.segs.df.new <- all.samp.w.metadata.conv.df[, correct.order.cols]

  print("Getting proportion of chromosome arms in selected genome covered by segments from each (and all) samples:")
  overlaps.between.samples.df <- compareOverlapsBetweenSamples(sample.vec, f.joint.segs.df, chrom.info.df, chrs.to.use.vec= CHRS_TO_USE_VEC)
  print(overlaps.between.samples.df)
  #So the question is here are we going to allow single base bair segments?
  print("Calculating intersections by GRanges...")
  intersection.grs.list <- getAllPatientSamplesAndIntersectionGRangesList(sample.vec, f.joint.segs.df)
  intersection.grs <- intersection.grs.list[["overall_intersection"]]
  non.ones.all.samp.intersection.grs <- removeSingleBasePairSegmentsFromIntersectionGranges(intersection.grs)

  #disjoin gets rid of all metadata
  print("Calculating disjoint set by GRanges...")
  dj.all.samp.grs <- disjoin(all.samp.grs)
  print("Removing segments that are smaller by minimum width by merging with the closest segment...")#Todo revisit this potentially
  non.ones.dj.all.samp.grs <- removeSingleBasePairSegmentsFromIntersectionGranges(dj.all.samp.grs)


  #Now we need to make sure that every bp of the relevant genome's non-acrocentric chromosome arms has a
  #corresponding segment, dummy or otherwise
  print("Creating set of segments that cover all bp in the non-acrocentric selected chromosomes' arms of the selected genome...")
  all.bp.chr.grs <- NULL
  for ( chr in CHRS_TO_USE_VEC ) {
    print(chr)
    curr.chr.grs <- non.ones.dj.all.samp.grs[seqnames(non.ones.dj.all.samp.grs) == chr]
    #This is done to stop intersect merging segments that start and end right after each other:
    alternating.idx.vec <- (seq(1, length(curr.chr.grs)) %% 2) + 1#plus one to match 1-based indexing
    plus.minus.vec <- c("+", "-")
    strand(curr.chr.grs) <- plus.minus.vec[alternating.idx.vec]
    chr.gr.list <- getChromosomeArmsGRangeList(chr, chrom.info.df)
    #So do this for the chromosome p and q arms (if p arm is present)
    for ( i in 1:length(chr.gr.list) ) {

      curr.chr.arm.gr <- chr.gr.list[[i]]
      curr.chr.arm.gr.plus <- curr.chr.arm.gr
      curr.chr.arm.gr.minus <- curr.chr.arm.gr
      strand(curr.chr.arm.gr.plus) <- "+"
      strand(curr.chr.arm.gr.minus) <- "-"
      curr.chr.arm.gr.strands <- c(curr.chr.arm.gr.plus, curr.chr.arm.gr.minus)
      #intersect taking into account strand to prevent merging of adjacent segments:
      overlapping.grs <- intersect(curr.chr.grs, curr.chr.arm.gr.strands, ignore.strand = FALSE)
      #Reset starand as it isn't used for antyhing else
      strand(overlapping.grs) <- "*"
      #No need to take into account strand as there is no way two segments would be merged by the setdiff operation.
      non.overlapping.grs <- setdiff(curr.chr.arm.gr, curr.chr.grs, ignore.strand = TRUE)
      sorted.all.bp.arm.grs <- sort(c(overlapping.grs, non.overlapping.grs))

      if ( is.null(all.bp.chr.grs) ) {

        all.bp.chr.grs <- sorted.all.bp.arm.grs

      } else {

        all.bp.chr.grs <- c(all.bp.chr.grs, sorted.all.bp.arm.grs)
      }
    }
  }


  print("Mapping the segments covering all bp of all arms back to input segments where appropriate....")
  ptm <- proc.time()
  itx <- iter(as.data.frame(all.bp.chr.grs), by="row")
  results <- foreach(curr.all.bp.chr.df = itx, .combine = c) %dopar% {

   # for (i in 1:length(all.bp.chr.grs)) {

   #  curr.all.bp.chr.df <- all.bp.chr.grs[i]
    require(GenomicRanges)
    curr.all.bp.chr.gr <- makeGRangesFromDataFrame(curr.all.bp.chr.df,
                            keep.extra.columns=TRUE,
                            ignore.strand=FALSE,
                            seqinfo=NULL,
                            seqnames.field=c("seqnames", "seqname",
                                             "chromosome", "chrom",
                                             "chr", "chromosome_name",
                                             "seqid"),
                            start.field=c("Start_Position","start","startpos"),
                            end.field=c("End_Position", "end", "stop", "endpos"),
                            strand.field="strand",
                            starts.in.df.are.0based=FALSE)

    this.seg.name <- paste0(as.character(seqnames(curr.all.bp.chr.gr)), ":", start(curr.all.bp.chr.gr), ":", end(curr.all.bp.chr.gr))
    this.seg.samp.seg.names.vec <- paste(sample.vec, ":", this.seg.name, sep="")
    info.cols.vec <- c("overlap_w_n_samps","num_segs_maps_to_in_samp" ,"in_intersection", "in_disjoin_only", "in_arm_setdiff")

    all.samp.this.seg.bp.df <- data.frame(matrix(FALSE, length(sample.vec), length(info.cols.vec),
                            dimnames = list(this.seg.samp.seg.names.vec, info.cols.vec)),
                            stringsAsFactors = FALSE,
                            check.names = FALSE)

    names(this.seg.samp.seg.names.vec) <- sample.vec
    all.samp.overlap.grs <- getOverlapGRanges(curr.all.bp.chr.gr, all.samp.w.metadata.grs)

    if (length(all.samp.overlap.grs) > 0) {

      all.samp.overlap.df <- as.data.frame(all.samp.overlap.grs, row.names = NULL)[0,][1:length(sample.vec),]
      rownames(all.samp.overlap.df) <- this.seg.samp.seg.names.vec

      matching.grs.sample.vec <- unique(all.samp.overlap.grs$sample)
      non.matching.sample.vec <- sample.vec[ ! sample.vec %in% matching.grs.sample.vec]
      all.samp.this.seg.bp.df[, "overlap_w_n_samps"] <- length(matching.grs.sample.vec)

    } else {

      all.samp.overlap.df <- as.data.frame(all.samp.overlap.grs)[1,]
      matching.grs.sample.vec <- c()
      non.matching.sample.vec <- sample.vec
    }

    if (all(sample.vec %in% matching.grs.sample.vec)) {

      all.samp.this.seg.bp.df[, "in_intersection"] <- TRUE

    } else if ( any(sample.vec %in% matching.grs.sample.vec) ) {#So now this means there's at least one segment missing

      all.samp.this.seg.bp.df[, "in_disjoin_only"] <- TRUE

    } else {

      all.samp.this.seg.bp.df[, "in_arm_setdiff"] <- TRUE

    }

    for ( match.samp in matching.grs.sample.vec ) {

      match.samp.overlap.grs <- all.samp.overlap.grs[all.samp.overlap.grs$sample == match.samp]
      all.samp.this.seg.bp.df[this.seg.samp.seg.names.vec[match.samp], "num_segs_maps_to_in_samp"] <- length(match.samp.overlap.grs)
      max.samp.overlap.gr <- match.samp.overlap.grs[ match.samp.overlap.grs$overlap_prop_query == max(match.samp.overlap.grs$overlap_prop_query) ]
      all.samp.overlap.df[this.seg.samp.seg.names.vec[match.samp],] <- as.data.frame(max.samp.overlap.gr, row.names=NULL)[1,]

    }

    if(length(non.matching.sample.vec) > 0) {

      all.samp.this.seg.bp.df[this.seg.samp.seg.names.vec[non.matching.sample.vec], "num_segs_maps_to_in_samp"] <- 0

    }


    #Creare a test dataframe
    return(list(cbind(all.samp.this.seg.bp.df, all.samp.overlap.df)))

  }
  #stopCluster(cl)
  time.passed <- proc.time() - ptm
  print("time.passed")
  print(time.passed)

  all.bp.all.samps.segs.df <- do.call(rbind, results)

  full_sample <- sapply(rownames(all.bp.all.samps.segs.df), function(x) strsplit(x, split = ":")[[1]][1])
  Chromosome <- as.numeric(sapply(rownames(all.bp.all.samps.segs.df), function(x) strsplit(x, split = ":")[[1]][2]))
  Start_Position <- as.numeric(sapply(rownames(all.bp.all.samps.segs.df), function(x) strsplit(x, split = ":")[[1]][3]))
  End_Position <- as.numeric(sapply(rownames(all.bp.all.samps.segs.df), function(x) strsplit(x, split = ":")[[1]][4]))

  #Let's make a pretend set of copy number columns that we can apply to the

  all.bp.all.samps.segs.anno.df <- cbind(Chromosome, Start_Position, End_Position, full_sample, all.bp.all.samps.segs.df)
  sorted.all.bp.all.samps.segs.anno.df <- all.bp.all.samps.segs.anno.df[with(all.bp.all.samps.segs.anno.df, order(full_sample, Chromosome, Start_Position)), ]
  sorted.all.bp.all.samps.segs.anno.df[1:50, 1:9]



  ###Add the chromosoem arms to the segments to make later joininbg etc much easier...
  print("Adding chromosome arm annotation to all segments original and new...")
  all.samp.all.bp.anno.grs <- makeGRangesFromDataFrame(sorted.all.bp.all.samps.segs.anno.df[,1:4],
                          keep.extra.columns=TRUE,
                          ignore.strand=FALSE,
                          seqinfo=NULL,
                          seqnames.field=c("Chromosome"),
                          start.field=c("Start_Position"),
                          end.field=c("End_Position"),
                          strand.field="*",
                          starts.in.df.are.0based=FALSE)
  all.samp.all.bp.anno.grs$chr_arm <- "none"

  for ( chr in CHRS_TO_USE_VEC ) {

    chr.gr.list <- getChromosomeArmsGRangeList(chr, chrom.info.df)
    #print(chr)
    if (length(chr.gr.list) == 2 ) {

      p.or.q.vec <- c("p", "q")

    } else if (length(chr.gr.list) == 1 ) {

      p.or.q.vec <- c("q")

    } else {

      stop("Error, no matching arms")

    }

    for ( i in 1:length(chr.gr.list) ) {

      curr.chr.arm.gr <- chr.gr.list[[i]]
      all.samp.overlap.grs <- getOverlapGRanges(all.samp.all.bp.anno.grs, curr.chr.arm.gr)
      all.samp.all.bp.anno.grs[ all.samp.overlap.grs$query_hit_idx ]$chr_arm <- paste0(chr, p.or.q.vec[i])

    }#end for ( i in 1:length(chr.gr.list) ) {

  }#end   for ( chr in CHRS_TO_USE_VEC ) {
  sorted.all.bp.all.samps.segs.anno.df$chr_arm <- all.samp.all.bp.anno.grs$chr_arm
  #So now we have chromosome arms on everything...

  return(sorted.all.bp.all.samps.segs.anno.df)

}#end createFullyAnnotatedJointSegmentation

#TBKW As of 20190818 this function is NOT used 
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

#TBKW As of 20190818 this function is NOT used 
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


