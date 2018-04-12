#################################### End of function definitions ##########################################################################
#####################################################################################################################################################

funEvaluateCopyITH <- function(TRACERxID
                               ,segForEvaluation
                               ,tracerxReleaseDir="/farm/tracerx/lung/release/"
                               ,arrayType
                               ,ploidy='mean'
                               ,plotting=FALSE
                               ,use.raw=TRUE
                               ,plottingDir="~/Google_Drive/Manuscripts/MSAI_paper/tempFigures/"
                               ,gammaToUse
                               ,methodToPhase=NULL
                               ,pValThresh = 0.01
)
{
  # select the sample of interest
  PatientSeg <- segForEvaluation[segForEvaluation$LTXID%in%TRACERxID,,drop=FALSE]
  patientCPNtable <- PatientSeg
  if(length(unique(patientCPNtable$sample))==1)
  {
    return(NA)
    #next;
  }
  regions <- unique(patientCPNtable$regionName)
  
  
  patientCPNtable$rawCopyNum  <- as.numeric(patientCPNtable$nBraw)+as.numeric(patientCPNtable$nAraw)
  patientCPNtable$minRaw      <- apply(cbind(patientCPNtable$nBraw,patientCPNtable$nAraw),1,min)
  patientCPNtable$maxRaw      <- apply(cbind(patientCPNtable$nBraw,patientCPNtable$nAraw),1,max)
  patientCPNtable$cpnToUse    <- patientCPNtable$cnTotal
  patientCPNtable$cpnPloidyCorrected <- NA
  
  purityEstimates <- patientCPNtable[!duplicated(patientCPNtable$sample),'ACF']
  names(purityEstimates) <- unique(patientCPNtable$sample )
  
  ploidyEstimates <- patientCPNtable[!duplicated(patientCPNtable$sample),'Ploidy']
  names(ploidyEstimates) <- unique(patientCPNtable$sample )
  
  
  if(use.raw){
    patientCPNtable$cpnToUse <-patientCPNtable$rawCopyNum 
  }
  
  if(ploidy=='main'){
    samp.ploidy<- calc.ploidy(patientCPNtable)
    patientCPNtable$cpnPloidyCorrected <- log2(patientCPNtable$cpnToUse/samp.ploidy[match(patientCPNtable$sample,rownames(samp.ploidy)),2])
    
  }
  if(ploidy=='mean'){
    patientCPNtable$cpnPloidyCorrected <- log2(patientCPNtable$cpnToUse/patientCPNtable$Ploidy)
    
  } 
  
  patientmin.cons.CPN     <- min.cons.fn(patientCPNtable, chromosomes=1:22,colToUse = "cpnPloidyCorrected") 
  patientmin.cons.CPN[patientmin.cons.CPN == '-Inf'] <- -3
  patientmin.cons.CPN[patientmin.cons.CPN == 'Inf'] <- 3
  
  # first, do this the standard way
  PatientClonalDissectionAmpsAndDels <- freq_clonalcn.fn(min.cons.seg = patientmin.cons.CPN
                                                         ,sample.names = colnames(patientmin.cons.CPN[,4:ncol(patientmin.cons.CPN)])
                                                         ,type = 'all'
                                                         ,gain.threshold = log2(2.5/2)
                                                         ,loss.threshold = log2(1.5/2)
                                                         ,min.gain.threshold = log2(2.5/2)
                                                         ,min.loss.threshold = log2(1.5/2)
  )
  
  # next, do it, but with the actual logR and BAF values. 
  # perform a t-test to explore whether regions are significantly different from one another. 
  # first, get the logR and BAF for the patient
  fullName <- grep(TRACERxID,list.files(tracerxReleaseDir),value=TRUE)
  print(tracerxReleaseDir)
  if(length(fullName)==0)
  {
    stop('cannot find TRACERxID in tracerxReleaseDir: You probably need to mount drive!')
  }
  DirToConsider <- paste(tracerxReleaseDir,"/",fullName,sep="")
  DirToConsider <- paste(DirToConsider,"/",arrayType,"/ASCAT/",sep="")
  
  #ASCATobj  <- getRegionASCATobj(DirToConsider = DirToConsider)
  ASCATdata <- getASCAT_BAF_logR(SampleID = fullName
                                 ,tracerxReleaseDir = tracerxReleaseDir
                                 ,arrayType = arrayType
                                 ,returnNumSNPs = FALSE
                                 ,everything = TRUE
                                 ,performPlotting = FALSE
  )
  min.cons.samples <- patientmin.cons.CPN[,4:ncol(patientmin.cons.CPN)]
  patientmin.cons.CPN <- patientmin.cons.CPN[!is.na(rowSums(min.cons.samples)),,drop=FALSE]
  min.cons.samples <- min.cons.samples[!is.na(rowSums(min.cons.samples)),,drop=FALSE]
  
  regionList <- colnames(min.cons.samples)
  
  bigList <-list()
  AllSnpsToRemove <- c()
  for (region in regionList)
  {
    snpsToRemove <- c()
    for (x in 1:nrow(min.cons.samples))
    {
      logRBAF <-  ASCATdata[[region]][ASCATdata[[region]]$chrom==patientmin.cons.CPN[x,1]&ASCATdata[[region]]$pos>=(patientmin.cons.CPN[x,2])&ASCATdata[[region]]$pos<=patientmin.cons.CPN[x,3],,drop=FALSE]
      snpsToRemove <- c(snpsToRemove,rownames(logRBAF[is.na(logRBAF$BAFsegmented),]))
      logRBAF <-  logRBAF[!is.na(logRBAF$BAFsegmented),,drop=FALSE]
      
      bigList[[region]][[x]] <-logRBAF
      
    }
    
    AllSnpsToRemove <- c(AllSnpsToRemove,snpsToRemove)
  }
  AllSnpsToRemove <- unique(AllSnpsToRemove)
  
  
  bigListLogR <-list()
  AllSnpsToRemove <- c()
  for (region in regionList)
  {
    snpsToRemove <- c()
    for (x in 1:nrow(min.cons.samples))
    {
      logRBAF <-  ASCATdata[[region]][ASCATdata[[region]]$chrom==patientmin.cons.CPN[x,1]&ASCATdata[[region]]$pos>=patientmin.cons.CPN[x,2]&ASCATdata[[region]]$pos<=patientmin.cons.CPN[x,3],,drop=FALSE]
      #snpsToRemove <- c(snpsToRemove,rownames(logRBAF[is.na(logRBAF$BAFsegmented),]))
      #logRBAF <-  logRBAF[!is.na(logRBAF$BAFsegmented),,drop=FALSE]
      
      bigListLogR[[region]][[x]] <-logRBAF
      
    }
    
    #AllSnpsToRemove <- c(AllSnpsToRemove,snpsToRemove)
  }
  #AllSnpsToRemove <- unique(AllSnpsToRemove)
  
  # remove any offending SNPs (which are NA in one region)
  for (region in regionList)
  {
    for (x in 1:nrow(min.cons.samples))
    {
      bigList[[region]][[x]] <- bigList[[region]][[x]][!rownames(bigList[[region]][[x]])%in%AllSnpsToRemove,,drop=FALSE]
    }
  }
  
  regionList <- names(bigList)
  
  if(length(methodToPhase)>0)
  {
    
    outPutOfPhasing  <- performAIphasing(regionList = regionList
                                         ,bigList = bigList
                                         ,min.cons.seg = patientmin.cons.CPN
                                         ,plotting = FALSE
                                         ,method = methodToPhase
                                         ,purityEstimates = purityEstimates
                                         ,ploidyEstimates = ploidyEstimates
                                         ,gammaToUse = gammaToUse
    )
    whichRegionToPhaseFrom <- outPutOfPhasing$whichRegionToPhaseFrom
    newBAFtable            <- outPutOfPhasing$newBAFtable
  }
  
  if(length(methodToPhase)==0)
  {
    whichRegionToPhaseFrom <- NULL
    
  }
  
  
  regionComparisonList <- list()
  HeterogeneityScores  <- c()
  if(plotting==TRUE)
  {
    pdf(paste(plottingDir,TRACERxID,"_",methodToPhase,"_CPNHeterogeneity.pdf",sep=""),width = 12,height=8)
  }
  for(w in 1:(length(regionList)-1)){
    for(j in (w+1):length(regionList)){
      
      regionOne <- regionList[w]
      regionTwo <- regionList[j]
      
      regionComOutPut <- calcDiffBetweenTwoRegions_withLogR( regionOne = regionOne
                                                    ,regionTwo = regionTwo
                                                    ,bigList = bigList
                                                    ,bigListLogR = bigListLogR
                                                    ,gammaToUse = gammaToUse 
                                                    ,min.cons.seg = patientmin.cons.CPN
                                                    ,plotting = plotting
                                                    ,minYaxis = (-4)
                                                    ,maxYaxis = 4
                                                    ,purityEstimates=purityEstimates
                                                    ,ploidyEstimates=ploidyEstimates
                                                    ,whichRegionToPhaseFrom = whichRegionToPhaseFrom
      )
      
      
      
      # regionComOutPut <- calcDiffBetweenTwoRegions( regionOne = regionOne
      #                                               ,regionTwo = regionTwo
      #                                               ,bigList = bigList
      #                                               ,min.cons.seg = patientmin.cons.CPN
      #                                               ,plotting = plotting
      #                                               ,minYaxis = (-4)
      #                                               ,maxYaxis = 4
      #                                               ,purityEstimates=purityEstimates
      #                                               ,ploidyEstimates=ploidyEstimates
      #                                               , whichRegionToPhaseFrom = NULL
      # )
      #regionComOutPut$`L_LTX064_SU_T1-R1:L_LTX064_SU_T1-R2`
      for (chrSeg in names(regionComOutPut)[1:(length(names(regionComOutPut))-1)])
      {
        regionComparisonList[[regionOne]][[chrSeg]] <- regionComOutPut[[chrSeg]][[regionOne]]
        regionComparisonList[[regionTwo]][[chrSeg]] <- regionComOutPut[[chrSeg]][[regionTwo]]
        regionComparisonList[[paste(regionOne,regionTwo,sep=":")]][[chrSeg]] <- regionComOutPut[[chrSeg]][[paste(regionOne,regionTwo,sep=":")]]
        
      }
      regionComparisonList[[paste(regionOne,regionTwo,sep=":")]]$HeterogeneityScore <- regionComOutPut[[paste(regionOne,regionTwo,sep=":")]]$HeterogeneityScore
      
      
      HeterogeneityScores <- c(HeterogeneityScores,regionComparisonList[[paste(regionOne,regionTwo,sep=":")]]$HeterogeneityScore)
      names(HeterogeneityScores)[length(HeterogeneityScores)] <- paste(regionOne,regionTwo,sep=":")
      
    }}
  if(plotting==TRUE)
  {
    dev.off()
  }
  
  outPut <- list()
  outPut$HeterogeneityScores <- HeterogeneityScores
  outPut$regionComparisonList <- regionComparisonList
  
  # calculate trunk and branch gain and loss
  
  ITHcalcs <- getProportionTrunkBranchAberrant(regionComparisonList,0.001, TRACERxID=TRACERxID)
  
  outPut$simpleSummary    <- ITHcalcs$simpleSummary
  outPut$simpleSummaryNew <- ITHcalcs$simpleSummaryNew
  outPut$minCopyLists     <- ITHcalcs$minCopy
  outPut$newBAFtable      <- newBAFtable
  outPut$whichRegionToPhaseFrom      <- whichRegionToPhaseFrom
  
  
  return(outPut)
}

#########################################################################################Functions called by funEvaluateCopyITH
############################################################################################################################################################
min.cons.fn<- function(segments, chromosomes = 1:22, mc.cores=2,colToUse = "cpnPloidyCorrected"){
  print("Called min.cons.fn")
  #assumes the following columns (and order)
  # sample, chr, startpos, endpos
  require(parallel)
  segments <- as.matrix(segments)
  chrom.all<- c()
  start.all<- c()
  end.all<- c()
  for(i in chromosomes){
    #   cat(i)
    s.chrom<- subset(segments, as.numeric(segments[,2]) == i)
    start<- unique(as.numeric(s.chrom[, 3]))
    end<-  unique(as.numeric(s.chrom[, 4]))
    end1<- 0
    start.final<- c()
    end.final<- c()
    #   print(i)
    count<- 1
    while(end1<  max(start)&  max(end) >  min(start[start >  end1])){
      #This will give a warning... Poor coding, should fix
      start1<- min(start[start>  end1])
      end1<- min(end[end>  start1])
      start.final<- c(start.final, start1)
      end.final<- c(end.final, end1)
      count<- count + 1
    }
    chrom.all<- c(chrom.all, rep(i, length(start.final)))
    start.all<- c(start.all, start.final)
    end.all<- c(end.all, end.final)
  }
  mat.pos<- matrix(NA, length(start.all), nc = 3)
  colnames(mat.pos)<- c("Chrom", "Start", "End")
  mat.pos<- cbind(as.numeric(chrom.all), as.numeric(start.all), as.numeric(end.all))
  #print(head(mat.pos))
  # cat('\nAssembling matrix....\n')
  # mat.tmp.cn<- t(apply(mat.pos, 1,fun.mapcopy.pos, seg = segments))
  tmp.out <- mclapply(unique(mat.pos[,1]), function(x){
    t(apply(mat.pos[mat.pos[,1] %in% x,,drop=F], 1, fun.mapcopy.pos, seg = segments,colToUse=colToUse))}, mc.cores=mc.cores)
  mat.tmp.cn <- do.call(rbind, tmp.out)
  #print(head(mat.tmp.cpn))
  #print(head(mat.pos))
  # print(str(mat.pos))
  #print(str(mat.tmp.cn))
  mat.final.cn<- cbind(mat.pos, mat.tmp.cn)
  colnames(mat.final.cn)<- c("Chr", "Start", "End", as.character(unique(segments[,1])))
  print("Returning from min.cons.fn")
  return(mat.final.cn)
}


fun.mapcopy.pos<-function(chrom.pos, seg,colToUse = "cpnPloidyCorrected") {
  print("Called fun.mapcopy.pos")
  vect.final <- rep(NA, length(unique(seg[, 1])))
  names(vect.final) <- unique(seg[, 1])
  sub  <- subset(seg, as.numeric(seg[, 2]) == chrom.pos[1])
  sub2 <- cbind(sub[,c(2:5)],sub[,colToUse])
  ww <- which(((as.numeric(sub2[, 2]) <= chrom.pos[3] &  as.numeric(sub2[, 3]) >= chrom.pos[ 3]  )
               |   (as.numeric(sub2[, 2]) <= chrom.pos[2] &  as.numeric(sub2[, 3]) >= chrom.pos[2])  
               |   (as.numeric(sub2[, 2]) >= chrom.pos[2]  &  as.numeric(sub2[, 3]) <= chrom.pos[3]  )
  ))
  
  sub.tmp <- sub[ww, ,drop=F]
  
  dups <- duplicated(sub.tmp[, 1], fromLast = T) == T | duplicated(sub.tmp[, 1]) == T
  
  n.dup <- sub.tmp[!dups, ,drop=F]
  vect.final[n.dup[, 1]] <- as.numeric(n.dup[, colToUse]) 
  
  if(length(which(dups)) >= 1){
    dup <- sub.tmp[dups, ]
    val.dup <- rep(NA, length(unique(dup[, 1])))
    names(val.dup) <- unique(dup[, 1])
    for(kk in unique(dup[, 1])){
      w.tmp <- which(dup[, 1] == kk)
      tmp <- c()
      for(ii in w.tmp)
        tmp <- c(tmp, pmax(c(chrom.pos[ 2], chrom.pos[3])) - pmin(c(as.numeric(dup[ii, 3]), as.numeric(dup[ii, 4]))))
      val.dup[kk] <- dup[w.tmp,  colToUse][which(tmp == max(tmp))[1]]      
    } 
    
    vect.final[names(val.dup)] <- as.numeric(val.dup) 
  }   
  print("Reurning from fun.mapcopy.pos")
  return(vect.final)               
}


calc.ploidy <-function(seg, cnCol=6,check.names=FALSE){
  print("Called calc.ploidy")
  #20140311 : added wMajor, weighted ploidy by chromosome
  #20151209 : Added the option to call mean raw copy if it exists
  #seg = ASCAT segmented output (or any segmented output), with total CN in column 6
  seg[,1] <- as.character(seg[,1])
  if(check.names){
    seg <- check.names.fn(seg)
  }
  samples <- unique(seg[,1])
  out.ploidy <- matrix(nrow=length(samples), ncol=14)
  rownames(out.ploidy) <- samples
  colnames(out.ploidy) <- c('Major','wMajor','MeanRaw',0:10) #wMajor is weighted by chromosome, MeanRaw is mean of the raw CN values from ASCAT 2.3
  for(j in samples){
    sample.seg <- seg[seg[,1] %in% j,]
    ploidy <- vector()
    for(k in unique(sample.seg[,cnCol])){
      tmp <- sample.seg[sample.seg[,cnCol] %in% k,]
      ploidy <- c(ploidy, setNames(sum(tmp[,4]-tmp[,3]), k))
    }
    ploidy <- (ploidy/sum(ploidy))[order(ploidy/sum(ploidy),decreasing=T)]
    out.ploidy[j,1] <- as.numeric(names(ploidy)[1])
    out.ploidy[j,3:13] <- ploidy[colnames(out.ploidy)[3:13]]
    wMajor <- c()
    for(k in unique(sample.seg[,2])){
      tmp <- sample.seg[sample.seg[,2] %in% k,]
      tmp <- setNames(tmp[,4]-tmp[,3], tmp[,cnCol])
      wMajor <- c(wMajor,setNames(as.numeric(names(sort(tapply(tmp, names(tmp),sum), decreasing=T)[1])), k))
    }
    out.ploidy[j,2] <- as.numeric(names(sort(table(wMajor), decreasing=T)[1]))
    if(any(grepl('nAraw', colnames(seg)))){
      a <- rowSums(sample.seg[,c('nAraw','nBraw')])
      b <- sample.seg[,4]-sample.seg[,3]
      out.ploidy[j,'MeanRaw'] <- sum(a * b/sum(b))
    }
  }
  if(!any(grepl('nAraw', colnames(seg)))){
    out.ploidy <- out.ploidy[,-c(3)]
  }
  print("Returning from calc.ploidy")
  return(out.ploidy)
}


freq_clonalcn.fn<-function(min.cons.seg, sample.names, type, gain.threshold=log2(2.5/2), loss.threshold=log2(1.5/2),min.gain.threshold=gain.threshold, min.loss.threshold=loss.threshold, min.samples=2){ 
  print("Called freq_clonalcn.fn")
  samples <- sort(unique(substr(sample.names,3,8)))
  gain.min.cons <- loss.min.cons <- gain.trunk.min.cons <- gain.branch.min.cons <-loss.trunk.min.cons <-loss.branch.min.cons <- cbind(min.cons.seg[,1:3],matrix(0, nrow=nrow(min.cons.seg), ncol=length(samples)))
  colnames(loss.min.cons) <-colnames(gain.min.cons) <-colnames(gain.trunk.min.cons)<-colnames(gain.branch.min.cons) <- colnames(loss.trunk.min.cons) <-colnames(loss.branch.min.cons) <- c(colnames(gain.trunk.min.cons)[1:3], sort(samples))
  
  for(i in 4:ncol(gain.trunk.min.cons)){
    min.cons.samples <- min.cons.seg[,substr(colnames(min.cons.seg),3,8) %in% colnames(gain.trunk.min.cons)[i],drop=F]
    
    clonalgain <- t(apply(min.cons.samples,1,function(x){
      tmp <- x
      tmp[tmp>=gain.threshold] <- 1
      tmp <- tmp[!is.na(tmp)]
      if(any(tmp >= 1)){
        tmp[tmp>= min.gain.threshold] <- 1
      }
      tmp[tmp < 1] <- 0
      c(sum(tmp),length(tmp))
    }))
    clonalloss <- t(apply(min.cons.samples,1,function(x){
      tmp <- x
      tmp[tmp<= loss.threshold] <- -1
      tmp <- tmp[!is.na(tmp)]
      if(any(tmp <= -1)){
        tmp[tmp<= min.loss.threshold] <- -1
      }
      tmp[tmp > -1] <- 0
      c(abs(sum(tmp)), length(tmp))
    }))
    
    gain.branch.min.cons[clonalgain[,1] >0 & clonalgain[,1] < clonalgain[,2],i] <- 1
    loss.branch.min.cons[clonalloss[,1] >0 & clonalloss[,1] < clonalloss[,2],i] <- -1
    
    gain.trunk.min.cons[clonalgain[,1] == clonalgain[,2],i] <- 1
    loss.trunk.min.cons[clonalloss[,1] == clonalloss[,2],i] <- -1
    
    
    gain.min.cons[clonalgain[,1] >0,i] <- 1
    loss.min.cons[clonalloss[,1] >0,i] <- -1
    gain.min.cons[clonalgain[,2] < min.samples,i] <- gain.trunk.min.cons[clonalgain[,2] < min.samples,i] <- gain.branch.min.cons[clonalgain[,2] < min.samples,i] <- NA
    loss.min.cons[clonalloss[,2] < min.samples,i] <- loss.trunk.min.cons[clonalloss[,2] < min.samples,i] <-loss.branch.min.cons[clonalloss[,2] < min.samples,i] <- NA
  }
  
  if(type=='all' ){
    outlist <- list(gain.min.cons, loss.min.cons, gain.trunk.min.cons, loss.trunk.min.cons,gain.branch.min.cons, loss.branch.min.cons)
    names(outlist) <- c('all_gain','all_loss','trunk_gain','trunk_loss','branch_gain','branch_loss')
  }
  if(type=='clonal' ){
    outlist <- list(gain.trunk.min.cons, loss.trunk.min.cons)
  }
  if(type=='subclonal' ){
    outlist <- list(gain.branch.min.cons, loss.branch.min.cons)
  }
  print("Returning from freq_clonalcn.fn")
  return(outlist)
}



getRegionASCATobj <- function(DirToConsider,arrayType="exome")
{
  print("Called getRegionASCATobj")
  outObject <- list()
  regions <- grep("FAILED",grep(".txt$|.complete$|*EXCLUDE*|*ManualQC*",list.files(DirToConsider),value=TRUE,invert=TRUE),invert=TRUE,value=TRUE)
  #print(regions)
  if(arrayType%in%c('exome','snp_array'))
  {
    for (region in regions)
    {
      ASCATobjTemp <- load(paste(DirToConsider,region,"/",region,".ascat.bc.RData",sep=""))
      ASCATobjTemp <- get(ASCATobjTemp)
      
      outObject[[region]] <- ASCATobjTemp
      
    }
  }
  if(arrayType%in%c('WGS'))
  {
    for (region in regions)
    {
      ASCATobjTemp <- load(paste(DirToConsider,region,"/",region,".tumor.baf_logr.RData",sep=""))
      ASCATobjTemp <- get(ASCATobjTemp)
      
      outObject[[region]] <- ASCATobjTemp
      
    }
  }
  
  #print(names(outObject))
  print("Return getRegionASCATobj")
  return(outObject)
}


getASCAT_BAF_logR <- function(SampleID 
                              ,tracerxReleaseDir="/farm/tracerx/lung/release/" 
                              ,arrayType="snp_array"
                              ,chr=6 ,start=29909037 ,stop =31324964 ,padding=15000000 ,performPlotting=TRUE ,returnNumSNPs=FALSE ,everything =FALSE)
{
  print("Called getASCAT_BAF_logR")
  DirToConsider <- paste(tracerxReleaseDir,SampleID,sep="")
  DirToConsider <- paste(DirToConsider,"/",arrayType,"/ASCAT/",sep="")
  
  if(arrayType%in%c('WGS'))
  {
   
      
      #load("/home/camp/watkint/working/watkint/data/yates/release/Y_BTX001/WGS/CN_heterogeneity/regionASCATsummary.RData")
      
      #logr.df <- read.table("/home/camp/watkint/working/watkint/data/yates/battenberg_logr_custom/Y_BTX001/Y_BTX001_new_snp_df.txt", stringsAsFactors=FALSE, sep="\t", header=TRUE)
      #logr.df$chr_pos <- paste(logr.df$Chr, "_", logr.df$Position, sep="")
      
      #return(regionASCATsummary)
      
      
      regions <- grep("FAILED",grep(".txt$|.complete$|*EXCLUDE*|*ManualQC*",list.files(DirToConsider),value=TRUE,invert=TRUE),invert=TRUE,value=TRUE)
      
      #So now we have the regons we want to load in the LogR and BAF data
      #We also want to load in the segmentation data so we can spoof it
      regionASCATsummary <- list()
      regionSNPsummary <- list()
      
      for ( i in 1:length(regions) ) {
        print("regions[i]")
        print(regions[i])
        reg.logr.baf.path <- paste(DirToConsider, regions[i], "/", regions[i], ".tumor.baf_logr.RData", sep="")
        load(reg.logr.baf.path)
        #tumor.baf_logr
        print("tumor.baf_logr")
        print(head(tumor.baf_logr))
        
        reg.seg.path <- paste(DirToConsider, regions[i], "/", regions[i], ".ascat.seg.RData", sep="")
        load(reg.seg.path)
        print("seg")
        print(head(seg))
        #seg
        
        ASCATsummary               <- tumor.baf_logr[,c("chrom", "pos")]
        ASCATsummary$BAF           <- tumor.baf_logr$BAF
        ASCATsummary$logR          <- tumor.baf_logr[,regions[i]]
        ASCATsummary$BAFsegmented  <- NA
        ASCATsummary$logRsegmented <- NA
        print("nrow(seg)")
        print(nrow(seg))
        for ( j in 1:nrow(seg) ) {
          print("j")
          print(j)
          if(nrow(ASCATsummary[ASCATsummary$chrom==seg[j,]$chr & ASCATsummary$pos >= seg[j,]$startpos & ASCATsummary$pos <= seg[j,]$endpos,]) > 0) {
            if(seg[j,]$nMajor==seg[j,]$nMinor){
              
              ASCATsummary[ASCATsummary$chrom==seg[j,]$chr & ASCATsummary$pos >= seg[j,]$startpos & ASCATsummary$pos <= seg[j,]$endpos,]$BAFsegmented <- 0.5
              
            } else {
              
              avg.seg.baf <- 0.5 - mean(abs(tumor.baf_logr[tumor.baf_logr$chrom==seg[j,]$chr & tumor.baf_logr$pos >= seg[j,]$startpos & tumor.baf_logr$pos <= seg[j,]$endpos,]$BAF - 0.5))
              ASCATsummary[ASCATsummary$chrom==seg[j,]$chr & ASCATsummary$pos >= seg[j,]$startpos & ASCATsummary$pos <= seg[j,]$endpos,]$BAFsegmented <- avg.seg.baf
              
            }
            avg.seg.logr <- mean(tumor.baf_logr[tumor.baf_logr$chrom==seg[j,]$chr & tumor.baf_logr$pos >= seg[j,]$startpos & tumor.baf_logr$pos <= seg[j,]$endpos, regions[i]])
            ASCATsummary[ASCATsummary$chrom==seg[j,]$chr & ASCATsummary$pos >= seg[j,]$startpos & ASCATsummary$pos <= seg[j,]$endpos,]$logRsegmented <- avg.seg.logr
          }
        }#nrow(seg)
        
        regionASCATsummary[[regions[i]]] <- ASCATsummary
        print("str(regionASCATsummary[[regions[i]]])")
        print(str(regionASCATsummary[[regions[i]]]))
      }#length(regions)
      #save(regionASCATsummary, file="/home/camp/watkint/working/watkint/data/yates/release/Y_BTX001/WGS/CN_heterogeneity/regionASCATsummary.RData")
      
      
      # i <- 1 
      # region.names <- names(regionASCATsummary)
      # for ( i in 1:length(regionASCATsummary) ) {
      
      #   region.df <- regionASCATsummary[[i]]
      #   region <- region.names[i]
      #   print(region)
      #   mod.region <- paste(gsub("-", ".", region), "_LogR", sep="")
      #   print(mod.region)
      #   #region.df$chr_pos <- paste(region.df$chrom, "_", region.df$pos, sep="")
      #   region.df$logr <- logr.df[, mod.region]
      #   #region.df$chr_pos <- NULL
      #   regionASCATsummary[[region]] <- region.df
      # }
      # print(str(regionASCATsummary))
      
      
      return(regionASCATsummary)
    }#if(arrayType=="WGS")
  
  if(arrayType%in%c('exome','snp_array'))
  {
  getRegionASCATobj <- function(DirToConsider)
  {
    outObject <- list()
    regions <- grep("FAILED",grep(".txt$|.complete$|*EXCLUDE*|*ManualQC*",list.files(DirToConsider),value=TRUE,invert=TRUE),invert=TRUE,value=TRUE)
    #print(regions)
    for (region in regions)
    {
      ASCATobjTemp <- load(paste(DirToConsider,region,"/",region,".ascat.bc.RData",sep=""))
      ASCATobjTemp <- get(ASCATobjTemp)
      
      outObject[[region]] <- ASCATobjTemp
      
    }
    #print(names(outObject))
    return(outObject)
  }
  
  regionASCATobj <- getRegionASCATobj(DirToConsider)
  regionASCATsummary <- list()
  regionSNPsummary   <- list()
  for (region in names(regionASCATobj))
  {
    
    ascatRegionBC <- regionASCATobj[[region]]
    
    #let's put together a table with what we need
    ASCATsummary               <- ascatRegionBC$SNPpos
    ASCATsummary$BAF           <- NA
    ASCATsummary$BAFsegmented  <- NA
    ASCATsummary$logR          <- NA
    ASCATsummary$logRsegmented <- NA
    
    ASCATsummary[rownames(ascatRegionBC$Tumor_LogR_segmented),]$logRsegmented           <- ascatRegionBC$Tumor_LogR_segmented[,1]
    ASCATsummary[rownames(ascatRegionBC$Tumor_BAF),]$BAF                                <- ascatRegionBC$Tumor_BAF[[1]]
    ASCATsummary[rownames(ascatRegionBC$Tumor_BAF_segmented[[1]]),]$BAFsegmented        <- ascatRegionBC$Tumor_BAF_segmented[[1]][,1]
    ASCATsummary[rownames(ascatRegionBC$Tumor_LogR),]$logR                              <- ascatRegionBC$Tumor_LogR[[1]]
    
    if(!everything)
    {
      ASCATsummary <- ASCATsummary[ASCATsummary$chrom%in%6&ASCATsummary$pos>=(start-padding)&ASCATsummary$pos<=(stop+padding),,drop=FALSE]
    }
    
    if(performPlotting)
    {
      par(mfrow=c(2,1))
      par(mar=c(2,2,2,2))
      plot(ASCATsummary$pos,ASCATsummary$logRsegmented,ylim=c(-1,1))
      points(ASCATsummary$pos,ASCATsummary$logR,ylim=c(-1,1))
      
      abline(v=start,lty='dashed',col='red')
      abline(v=stop,lty='dashed',col='red')
      
      plot(ASCATsummary$pos,ASCATsummary$BAFsegmented,ylim=c(0,1))
      points(ASCATsummary$pos[!is.na(ASCATsummary$BAFsegmented)],ASCATsummary$BAF[!is.na(ASCATsummary$BAFsegmented)],ylim=c(0,1))
      
      abline(v=start,lty='dashed',col='red')
      abline(v=stop,lty='dashed',col='red')
    }
    regionASCATsummary[[region]] <- ASCATsummary
    
    
  }
  
  
  HLA_genePos <- list()
  HLA_genePos[["hla_a"]] <- c(6,29909037,29913661)
  HLA_genePos[["hla_b"]] <- c(6,31321649,31324964) 
  HLA_genePos[["hla_c"]] <- c(6,31236526,31239869)
  
  HLA_A_logRBAF <- regionASCATsummary[[1]][regionASCATsummary[[1]]$pos>=HLA_genePos[['hla_a']][2]&  regionASCATsummary[[1]]$pos<=HLA_genePos[['hla_a']][3],]
  HLA_A_logRBAF <- HLA_A_logRBAF[!is.na(HLA_A_logRBAF$BAFsegmented),,drop=FALSE]
  nrAsnps       <- nrow(HLA_A_logRBAF)
  
  HLA_B_logRBAF <- regionASCATsummary[[1]][regionASCATsummary[[1]]$pos>=HLA_genePos[['hla_b']][2]&  regionASCATsummary[[1]]$pos<=HLA_genePos[['hla_b']][3],]
  HLA_B_logRBAF <- HLA_B_logRBAF[!is.na(HLA_B_logRBAF$BAFsegmented),,drop=FALSE]
  nrBsnps       <- nrow(HLA_B_logRBAF)
  
  HLA_C_logRBAF <- regionASCATsummary[[1]][regionASCATsummary[[1]]$pos>=HLA_genePos[['hla_c']][2]&  regionASCATsummary[[1]]$pos<=HLA_genePos[['hla_c']][3],]
  HLA_C_logRBAF <- HLA_B_logRBAF[!is.na(HLA_B_logRBAF$BAFsegmented),,drop=FALSE]
  nrCsnps       <- nrow(HLA_C_logRBAF)
  
  
  outPut  <- cbind(nrAsnps,nrBsnps,nrCsnps)
  rownames(outPut) <- SampleID
  
  print("Returning from getASCAT_BAF_logR")
  if(returnNumSNPs)
  {
    return(outPut)
  }
  
  if(!returnNumSNPs)
  {
    return(regionASCATsummary)
  }
  
  }
  
}



performAIphasing  <- function(regionList
                              ,bigList
                              ,min.cons.seg
                              ,plotting=TRUE
                              ,gammaToUse = 1
                              ,method=c('sig.wilcox.test')
                              ,purityEstimates
                              ,ploidyEstimates
                              ,dipTestPvalThresh = 0.01
)
{
  
  
  if(method=='sig.wilcox.test')
  {
    whichRegionToPhaseFrom <- list()
    names_ComparisonOutPut <- apply(min.cons.seg[,1:3],1,PasteVector,sep=":")
    
    
    for (x in 1:length(names_ComparisonOutPut))
    {
      cat(x)
      regionCompairson  <- list()
      maxPhasingAbility <- c()
      showsAI           <- c()
      for (region in regionList)
      {
        # region of Interest
        logRBAFSample1 <- bigList[[region]][[x]] #ASCATdata[[regionOne]][ASCATdata[[regionOne]]$chrom==gain.min.cons[x,1]&ASCATdata[[regionOne]]$pos>=gain.min.cons[x,2]&ASCATdata[[regionOne]]$pos<=gain.min.cons[x,3],,drop=FALSE]
        if(nrow(logRBAFSample1)<5)
        {
          next;
        }
        
        
        AndBSample1 <- funCalcN_withBAF(logRSites = logRBAFSample1$logR,bafSites = logRBAFSample1$BAF,tumorPloidy = ploidyEstimates[region],tumorPurity = purityEstimates[region],gamma = gammaToUse)
        Sample1Max <- apply(AndBSample1,1,max)
        tmpSample1 <- ifelse(AndBSample1>=Sample1Max,1,0)
        MajorColumn1 <- which(tmpSample1[,1]==1)
        
        MajorColumn2    <- which(tmpSample1[,2]==1)
        #check whether significant in this sample (it will be)
        RefRegionTtest  <- tTestWithErrorHandling(c(AndBSample1[MajorColumn1,1],AndBSample1[MajorColumn2,2]),c(c(AndBSample1[MajorColumn1,2],AndBSample1[MajorColumn2,1])),paired=TRUE)
        showAI          <- as.numeric(tTestWithErrorHandling(c(AndBSample1[MajorColumn1,1],AndBSample1[MajorColumn2,2]),c(c(AndBSample1[MajorColumn1,2],AndBSample1[MajorColumn2,1])))$p.val)
        showsAI         <- c(showsAI,showAI)
        #check other regions
        OtherRegionPvals  <- c()
        PropSNPsSeparate  <- c()
        ConcordantWithRef <- c()
        for (otherRegion in regionList[!regionList%in%region])
        {
          logRBAFOtherRegion <- bigList[[otherRegion]][[x]] #ASCATdata[[regionOne]][ASCATdata[[regionOne]]$chrom==gain.min.cons[x,1]&ASCATdata[[regionOne]]$pos>=gain.min.cons[x,2]&ASCATdata[[regionOne]]$pos<=gain.min.cons[x,3],,drop=FALSE]
          AndBOtherRegion    <- funCalcN_withBAF(logRSites = logRBAFOtherRegion$logR,bafSites = logRBAFOtherRegion$BAF,tumorPloidy = ploidyEstimates[otherRegion],tumorPurity = purityEstimates[otherRegion],gamma = gammaToUse)
          OtherRegionTtest   <- wilcox.test(c(AndBOtherRegion[MajorColumn1,1],AndBOtherRegion[MajorColumn2,2]),c(c(AndBOtherRegion[MajorColumn1,2],AndBOtherRegion[MajorColumn2,1])),paired=TRUE)
          OtherRegionPvals   <- c(OtherRegionPvals,OtherRegionTtest$p.value)
          
          regionSNPseparate_sameAsReference    <- length(which(c(AndBOtherRegion[MajorColumn1,1],AndBOtherRegion[MajorColumn2,2]) -c(c(AndBOtherRegion[MajorColumn1,2],AndBOtherRegion[MajorColumn2,1]))>0))/nrow(AndBOtherRegion)
          regionSNPseparate_diffFromReference  <- length(which(c(AndBOtherRegion[MajorColumn1,1],AndBOtherRegion[MajorColumn2,2]) -c(c(AndBOtherRegion[MajorColumn1,2],AndBOtherRegion[MajorColumn2,1]))<0))/nrow(AndBOtherRegion)
          
          if(regionSNPseparate_sameAsReference>=regionSNPseparate_diffFromReference)
          {
            PropSNPsSeparate <- c(PropSNPsSeparate,regionSNPseparate_sameAsReference)
            ConcordantWithRef <- c(ConcordantWithRef,1)
          }
          
          if(regionSNPseparate_sameAsReference<regionSNPseparate_diffFromReference)
          {
            PropSNPsSeparate <- c(PropSNPsSeparate,regionSNPseparate_diffFromReference)
            ConcordantWithRef <- c(ConcordantWithRef,-1)
          }
          #PropSNPsSeparate   <- c(PropSNPsSeparate,length(which(c(AndBOtherRegion[MajorColumn1,1],AndBOtherRegion[MajorColumn2,2]) -c(c(AndBOtherRegion[MajorColumn1,2],AndBOtherRegion[MajorColumn2,1]))>0))/nrow(AndBOtherRegion))
          
          
          
        }
        combinedPvals <- pchisq((sum(log(OtherRegionPvals))*-2), df=length(OtherRegionPvals)*2, lower.tail=F)
        regionCompairson[[region]]$RefRegionPval        <- RefRegionTtest$p.value
        regionCompairson[[region]]$RefUnpaired          <- showAI
        regionCompairson[[region]]$combinedPval         <- combinedPvals
        regionCompairson[[region]]$otherPvals           <- OtherRegionPvals
        regionCompairson[[region]]$otherSNPsSeparate    <- PropSNPsSeparate
        regionCompairson[[region]]$ConcordantWithRef    <- ConcordantWithRef
        regionCompairson[[region]]$MajorColumn1         <- MajorColumn1
        regionCompairson[[region]]$MajorColumn2         <- MajorColumn2
        
        maxPhasingAbility <- c(maxPhasingAbility,mean(PropSNPsSeparate))
        
      }
      if(nrow(logRBAFSample1)<5)
      {
        next;
      }
      names(maxPhasingAbility) <- regionList
      names(showsAI) <- regionList
      
      maxPhasingAbility <- maxPhasingAbility[showsAI<0.05]
      
      if(length(maxPhasingAbility)==0)
      {
        next;
      }
      
      # choose a reference
      regionToChose <- names(which(maxPhasingAbility==max(maxPhasingAbility)))[1]
      
      segName <- names_ComparisonOutPut[[x]]
      whichRegionToPhaseFrom[[segName]]$regionToUse <- regionToChose
      whichRegionToPhaseFrom[[segName]]$MajorColumn1 <- regionCompairson[[regionToChose]]$MajorColumn1
      whichRegionToPhaseFrom[[segName]]$MajorColumn2 <- regionCompairson[[regionToChose]]$MajorColumn2
      whichRegionToPhaseFrom[[segName]]$regionComparisons <- regionCompairson
      
      
      
    }
    
    return(whichRegionToPhaseFrom)  
  }
  
  if(method=='combinedBAF')
  {
    whichRegionToPhaseFrom <- list()
    names_ComparisonOutPut <- apply(min.cons.seg[,1:3],1,PasteVector,sep=":")
    
    
    for (x in 1:length(names_ComparisonOutPut))
    {
      cat(x)
      regionCompairson  <- list()
      maxPhasingAbility <- c()
      showsAI           <- c()
      for (region in regionList)
      {
        # region of Interest
        logRBAFSample1 <- bigList[[region]][[x]] #ASCATdata[[regionOne]][ASCATdata[[regionOne]]$chrom==gain.min.cons[x,1]&ASCATdata[[regionOne]]$pos>=gain.min.cons[x,2]&ASCATdata[[regionOne]]$pos<=gain.min.cons[x,3],,drop=FALSE]
        if(nrow(logRBAFSample1)<5)
        {
          next;
        }
        
        AndBSample1 <- funCalcN_withBAF(logRSites = logRBAFSample1$logR,bafSites = logRBAFSample1$BAF,tumorPloidy = ploidyEstimates[region],tumorPurity = purityEstimates[region],gamma = gammaToUse)
        Sample1Max <- apply(AndBSample1,1,max)
        tmpSample1 <- ifelse(AndBSample1>=Sample1Max,1,0)
        MajorColumn1    <- which(tmpSample1[,1]==1)
        MajorColumn2    <- which(tmpSample1[,2]==1)
        
        #check whether significant in this sample (it will be)
        RefRegionTtest  <- tTestWithErrorHandling(c(AndBSample1[MajorColumn1,1],AndBSample1[MajorColumn2,2]),c(c(AndBSample1[MajorColumn1,2],AndBSample1[MajorColumn2,1])),paired=TRUE)
        showAI          <- as.numeric(tTestWithErrorHandling(c(AndBSample1[MajorColumn1,1],AndBSample1[MajorColumn2,2]),c(c(AndBSample1[MajorColumn1,2],AndBSample1[MajorColumn2,1])))$p.val)
        showsAI         <- c(showsAI,showAI)
        #check other regions
        OtherRegionPvals  <- c()
        PropSNPsSeparate  <- c()
        ConcordantWithRef <- c()
        for (otherRegion in regionList[!regionList%in%region])
        {
          logRBAFOtherRegion <- bigList[[otherRegion]][[x]] #ASCATdata[[regionOne]][ASCATdata[[regionOne]]$chrom==gain.min.cons[x,1]&ASCATdata[[regionOne]]$pos>=gain.min.cons[x,2]&ASCATdata[[regionOne]]$pos<=gain.min.cons[x,3],,drop=FALSE]
          AndBOtherRegion    <- funCalcN_withBAF(logRSites = logRBAFOtherRegion$logR,bafSites = logRBAFOtherRegion$BAF,tumorPloidy = ploidyEstimates[otherRegion],tumorPurity = purityEstimates[otherRegion],gamma = gammaToUse)
          OtherRegionTtest   <- wilcox.test(c(AndBOtherRegion[MajorColumn1,1],AndBOtherRegion[MajorColumn2,2]),c(c(AndBOtherRegion[MajorColumn1,2],AndBOtherRegion[MajorColumn2,1])),paired=TRUE)
          OtherRegionPvals   <- c(OtherRegionPvals,OtherRegionTtest$p.value)
          
          regionSNPseparate_sameAsReference    <- length(which(c(AndBOtherRegion[MajorColumn1,1],AndBOtherRegion[MajorColumn2,2]) -c(c(AndBOtherRegion[MajorColumn1,2],AndBOtherRegion[MajorColumn2,1]))>0))/nrow(AndBOtherRegion)
          regionSNPseparate_diffFromReference  <- length(which(c(AndBOtherRegion[MajorColumn1,1],AndBOtherRegion[MajorColumn2,2]) -c(c(AndBOtherRegion[MajorColumn1,2],AndBOtherRegion[MajorColumn2,1]))<0))/nrow(AndBOtherRegion)
          
          if(regionSNPseparate_sameAsReference>=regionSNPseparate_diffFromReference)
          {
            PropSNPsSeparate <- c(PropSNPsSeparate,regionSNPseparate_sameAsReference)
            ConcordantWithRef <- c(ConcordantWithRef,1)
          }
          
          if(regionSNPseparate_sameAsReference<regionSNPseparate_diffFromReference)
          {
            PropSNPsSeparate <- c(PropSNPsSeparate,regionSNPseparate_diffFromReference)
            ConcordantWithRef <- c(ConcordantWithRef,-1)
          }
          #PropSNPsSeparate   <- c(PropSNPsSeparate,length(which(c(AndBOtherRegion[MajorColumn1,1],AndBOtherRegion[MajorColumn2,2]) -c(c(AndBOtherRegion[MajorColumn1,2],AndBOtherRegion[MajorColumn2,1]))>0))/nrow(AndBOtherRegion))
          
          
          
        }
        combinedPvals <- pchisq((sum(log(OtherRegionPvals))*-2), df=length(OtherRegionPvals)*2, lower.tail=F)
        regionCompairson[[region]]$RefRegionPval        <- RefRegionTtest$p.value
        regionCompairson[[region]]$RefUnpaired          <- showAI
        regionCompairson[[region]]$combinedPval         <- combinedPvals
        regionCompairson[[region]]$otherPvals           <- OtherRegionPvals
        regionCompairson[[region]]$otherSNPsSeparate    <- PropSNPsSeparate
        regionCompairson[[region]]$ConcordantWithRef    <- ConcordantWithRef
        regionCompairson[[region]]$MajorColumn1         <- MajorColumn1
        regionCompairson[[region]]$MajorColumn2         <- MajorColumn2
        
        maxPhasingAbility <- c(maxPhasingAbility,mean(PropSNPsSeparate))
        
      }
      if(nrow(logRBAFSample1)<5)
      {
        next;
      }
      names(maxPhasingAbility) <- regionList
      names(showsAI) <- regionList
      
      maxPhasingAbility <- maxPhasingAbility[showsAI<0.05]
      
      if(length(maxPhasingAbility)==0)
      {
        next;
      }
      
      # choose a reference
      regionToChose <- names(which(maxPhasingAbility==max(maxPhasingAbility)))[1]
      
      # choose a reference
      if(mean(regionCompairson[[regionToChose]]$otherSNPsSeparate)>0.9)
        #make sure the region chosen 
      {
        
        segName <- names_ComparisonOutPut[[x]]
        whichRegionToPhaseFrom[[segName]]$regionToUse <- regionToChose
        whichRegionToPhaseFrom[[segName]]$MajorColumn1 <- regionCompairson[[regionToChose]]$MajorColumn1
        whichRegionToPhaseFrom[[segName]]$MajorColumn2 <- regionCompairson[[regionToChose]]$MajorColumn2
        whichRegionToPhaseFrom[[segName]]$regionComparisons <- regionCompairson
        
      }
      
      # make sure the region chosen 
      # now if the region isn't >0.9 across all regions
      if(regionCompairson[[regionToChose]]$combinedPval<0.01&mean(regionCompairson[[regionToChose]]$otherSNPsSeparate)<0.9)
      {
        # use a combination of all regions
        logRBAFSample1 <- bigList[[regionToChose]][[x]] #ASCATdata[[regionOne]][ASCATdata[[regionOne]]$chrom==gain.min.cons[x,1]&ASCATdata[[regionOne]]$pos>=gain.min.cons[x,2]&ASCATdata[[regionOne]]$pos<=gain.min.cons[x,3],,drop=FALSE]
        AndBSample1 <- funCalcN_withBAF(logRSites = logRBAFSample1$logR,bafSites = logRBAFSample1$BAF,tumorPloidy = ploidyEstimates[regionToChose],tumorPurity = purityEstimates[regionToChose],gamma = gammaToUse)
        
        nA <- c(AndBSample1[regionCompairson[[regionToChose]]$MajorColumn1,1],AndBSample1[regionCompairson[[regionToChose]]$MajorColumn2,2])
        nB <- c(AndBSample1[regionCompairson[[regionToChose]]$MajorColumn1,2],AndBSample1[regionCompairson[[regionToChose]]$MajorColumn2,1])
        
        OrderToConsider <- order(nA-nB)
        #barplot(c(nA-nB)[OrderToConsider],las=1)
        otherRegions  <- regionList[!regionList%in%regionToChose]
        
        OtherRegionConcordance        <- regionCompairson[[regionToChose]]$ConcordantWithRef
        names(OtherRegionConcordance) <- otherRegions
        RefRegion <- c(1)
        names(RefRegion) <- regionToChose
        OtherRegionConcordance <- c(RefRegion,OtherRegionConcordance)
        
        otherRegions  <- otherRegions[otherRegions%in%c(names(showsAI[showsAI<0.01]))]
        allSigRegions <- c(regionToChose,otherRegions)
        
        #check whether should be considered positive or negative
        MajorColumn1  <- c()
        MajorColumn2  <- c()
        meanSNPscores <- c()
        listScores    <- list()
        changeFromRef <- c()
        
        for (SNP in 1:nrow(AndBSample1))
        {
          columnCurrently <- NA
          if(SNP%in%regionCompairson[[regionToChose]]$MajorColumn1)
          {
            MajorCurrently <- 1
            MinorCurrently <- 2
          }
          
          if(SNP%in%regionCompairson[[regionToChose]]$MajorColumn2)
          {
            MajorCurrently <- 2
            MinorCurrently <-1
            
          }
          
          SNPscores <- c()
          for (region in allSigRegions)
          {
            logRBAFSample1 <- bigList[[region]][[x]] #ASCATdata[[regionOne]][ASCATdata[[regionOne]]$chrom==gain.min.cons[x,1]&ASCATdata[[regionOne]]$pos>=gain.min.cons[x,2]&ASCATdata[[regionOne]]$pos<=gain.min.cons[x,3],,drop=FALSE]
            AndBSample1 <- funCalcN_withBAF(logRSites = logRBAFSample1$logR,bafSites = logRBAFSample1$BAF,tumorPloidy = ploidyEstimates[region],tumorPurity = purityEstimates[region],gamma = gammaToUse)
            
            nA <- AndBSample1[SNP,MajorCurrently]
            nB <- AndBSample1[SNP,MinorCurrently]
            multFactor <- OtherRegionConcordance[region]
            
            DiffScore <- (nA-nB)*multFactor
            
            SNPscores <- c(SNPscores,DiffScore)
            listScores[[region]]$diff <- c(listScores[[region]]$diff,DiffScore)
            
          }
          
          if(mean(SNPscores)>0)
          {
            changeFromRef <- c(changeFromRef,0)
            #MajorToUse <- columnCurrently
            if(MajorCurrently==1)
            {
              MajorColumn1 <- c(MajorColumn1,SNP)
            }
            
            if(MajorCurrently==2)
            {
              MajorColumn2 <- c(MajorColumn2,SNP)
            }
            
            
          }
          
          if(mean(SNPscores)<0)
          {
            changeFromRef <- c(changeFromRef,1)
            
            if(MajorCurrently==1)
            {
              MajorColumn2 <- c(MajorColumn2,SNP)
            }
            
            if(MajorCurrently==2)
            {
              MajorColumn1 <- c(MajorColumn1,SNP)
            }
          }
          
          meanSNPscores <- c(meanSNPscores,mean(SNPscores))
          
        }
        
        # barplot(listScores$`L_LTX062_SU_T1-R1`$diff)
        # barplot(listScores$`L_LTX062_SU_T1-R2`$diff)
        # barplot(listScores$`L_LTX062_SU_T1-R3`$diff)
        # barplot(listScores$`L_LTX062_SU_T1-R4`$diff)
        # 
        
        segName <- names_ComparisonOutPut[[x]]
        whichRegionToPhaseFrom[[segName]]$regionToUse <- "mixture"
        whichRegionToPhaseFrom[[segName]]$MajorColumn1 <- MajorColumn1
        whichRegionToPhaseFrom[[segName]]$MajorColumn2 <- MajorColumn2
        
        
        #finally, check the ability of the combined to separate A and B in each region
        OtherRegionPvals  <- c()
        PropSNPsSeparate  <- c()
        ConcordantWithRef <- c()
        for (otherRegion in regionList)
        {
          logRBAFOtherRegion <- bigList[[otherRegion]][[x]] #ASCATdata[[regionOne]][ASCATdata[[regionOne]]$chrom==gain.min.cons[x,1]&ASCATdata[[regionOne]]$pos>=gain.min.cons[x,2]&ASCATdata[[regionOne]]$pos<=gain.min.cons[x,3],,drop=FALSE]
          AndBOtherRegion    <- funCalcN_withBAF(logRSites = logRBAFOtherRegion$logR,bafSites = logRBAFOtherRegion$BAF,tumorPloidy = ploidyEstimates[otherRegion],tumorPurity = purityEstimates[otherRegion],gamma = gammaToUse)
          OtherRegionTtest   <- wilcox.test(c(AndBOtherRegion[MajorColumn1,1],AndBOtherRegion[MajorColumn2,2]),c(c(AndBOtherRegion[MajorColumn1,2],AndBOtherRegion[MajorColumn2,1])),paired=TRUE)
          OtherRegionPvals   <- c(OtherRegionPvals,OtherRegionTtest$p.value)
          
          regionSNPseparate_sameAsReference    <- length(which(c(AndBOtherRegion[MajorColumn1,1],AndBOtherRegion[MajorColumn2,2]) -c(c(AndBOtherRegion[MajorColumn1,2],AndBOtherRegion[MajorColumn2,1]))>0))/nrow(AndBOtherRegion)
          regionSNPseparate_diffFromReference  <- length(which(c(AndBOtherRegion[MajorColumn1,1],AndBOtherRegion[MajorColumn2,2]) -c(c(AndBOtherRegion[MajorColumn1,2],AndBOtherRegion[MajorColumn2,1]))<0))/nrow(AndBOtherRegion)
          
          if(regionSNPseparate_sameAsReference>=regionSNPseparate_diffFromReference)
          {
            PropSNPsSeparate <- c(PropSNPsSeparate,regionSNPseparate_sameAsReference)
            ConcordantWithRef <- c(ConcordantWithRef,1)
          }
          
          if(regionSNPseparate_sameAsReference<regionSNPseparate_diffFromReference)
          {
            PropSNPsSeparate <- c(PropSNPsSeparate,regionSNPseparate_diffFromReference)
            ConcordantWithRef <- c(ConcordantWithRef,-1)
          }
          #PropSNPsSeparate   <- c(PropSNPsSeparate,length(which(c(AndBOtherRegion[MajorColumn1,1],AndBOtherRegion[MajorColumn2,2]) -c(c(AndBOtherRegion[MajorColumn1,2],AndBOtherRegion[MajorColumn2,1]))>0))/nrow(AndBOtherRegion))
          
          
          
        }
        
        
        names(PropSNPsSeparate) <- regionList
        names(ConcordantWithRef) <- regionList
        names(OtherRegionPvals)  <- regionList
        whichRegionToPhaseFrom[[segName]]$PropSNPsSeparate <- PropSNPsSeparate
        whichRegionToPhaseFrom[[segName]]$ConcordantWithRef <- ConcordantWithRef
        whichRegionToPhaseFrom[[segName]]$RegionPvals <- OtherRegionPvals
        whichRegionToPhaseFrom[[segName]]$RegionOrder <- regionList
        whichRegionToPhaseFrom[[segName]]$regionComparisons <- regionCompairson
        
        
        
        
      }
      
    }
    
    return(whichRegionToPhaseFrom)  
  }
  
  if(method=='combinedBAFDipTest')
  {
    whichRegionToPhaseFrom <- list()
    newBAFtable            <- list()
    
    
    names_ComparisonOutPut <- apply(min.cons.seg[,1:3],1,PasteVector,sep=":")
    
    
    for (x in 1:length(names_ComparisonOutPut))
    {
      cat(x)
      segName <- names_ComparisonOutPut[[x]]
      regionCompairson  <- list()
      maxPhasingAbility <- c()
      showsAI           <- c()
      ASCATcallsAI      <- c()
      for (region in regionList)
      {
        # region of Interest
        logRBAFSample1 <- bigList[[region]][[x]] #ASCATdata[[regionOne]][ASCATdata[[regionOne]]$chrom==gain.min.cons[x,1]&ASCATdata[[regionOne]]$pos>=gain.min.cons[x,2]&ASCATdata[[regionOne]]$pos<=gain.min.cons[x,3],,drop=FALSE]
        if(nrow(logRBAFSample1)<5)
        {
          next;
        }
        

        AndBSample1 <- funCalcN_withBAF(logRSites = logRBAFSample1$logR,bafSites = logRBAFSample1$BAF,tumorPloidy = ploidyEstimates[region],tumorPurity = purityEstimates[region],gamma = gammaToUse)
        Sample1Max <- apply(AndBSample1,1,max)
        tmpSample1 <- ifelse(AndBSample1>=Sample1Max,1,0)
        MajorColumn1    <- which(tmpSample1[,1]==1)
        MajorColumn2    <- which(tmpSample1[,2]==1)

        
        #check whether significant in this sample (it will be)
        
        # check whether significantly different from 0.5
        # tTestWithErrorHandling(logRBAFSample1$BAF-0.5)
        RefRegionTtest  <- dip.test(logRBAFSample1$BAF)
        showAI          <- as.numeric(dip.test(logRBAFSample1$BAF)$p.val)
        showsAI         <- c(showsAI,showAI)
        ASCATAI         <- ifelse(logRBAFSample1$BAFsegmented[1]!=0.5,TRUE,FALSE)
        ASCATcallsAI    <- c(ASCATcallsAI,ASCATAI)
        #check other regions
        OtherRegionPvals  <- c()
        PropSNPsSeparate  <- c()
        ConcordantWithRef <- c()
        for (otherRegion in regionList[!regionList%in%region])
        {
         
      
          logRBAFOtherRegion <- bigList[[otherRegion]][[x]] #ASCATdata[[regionOne]][ASCATdata[[regionOne]]$chrom==gain.min.cons[x,1]&ASCATdata[[regionOne]]$pos>=gain.min.cons[x,2]&ASCATdata[[regionOne]]$pos<=gain.min.cons[x,3],,drop=FALSE]
          AndBOtherRegion    <- funCalcN_withBAF(logRSites = logRBAFOtherRegion$logR,bafSites = logRBAFOtherRegion$BAF,tumorPloidy = ploidyEstimates[otherRegion],tumorPurity = purityEstimates[otherRegion],gamma = gammaToUse)

          #print(MajorColumn1)
          MajorColumn1 <- MajorColumn1[MajorColumn1 <= nrow(AndBOtherRegion)]
          newMajorColumn1 <- MajorColumn1[MajorColumn1 <= nrow(AndBOtherRegion)]
          MajorColumn2 <- MajorColumn2[MajorColumn2 <= nrow(AndBOtherRegion)]
          newMajorColumn2 <- MajorColumn2
          OtherRegionTtest   <- wilcox.test(c(AndBOtherRegion[MajorColumn1,1],AndBOtherRegion[MajorColumn2,2]),c(c(AndBOtherRegion[MajorColumn1,2],AndBOtherRegion[MajorColumn2,1])),paired=TRUE)

          OtherRegionPvals   <- c(OtherRegionPvals,OtherRegionTtest$p.value)
          
          regionSNPseparate_sameAsReference    <- length(which(c(AndBOtherRegion[MajorColumn1,1],AndBOtherRegion[MajorColumn2,2]) -c(c(AndBOtherRegion[MajorColumn1,2],AndBOtherRegion[MajorColumn2,1]))>0))/nrow(AndBOtherRegion)
          regionSNPseparate_diffFromReference  <- length(which(c(AndBOtherRegion[MajorColumn1,1],AndBOtherRegion[MajorColumn2,2]) -c(c(AndBOtherRegion[MajorColumn1,2],AndBOtherRegion[MajorColumn2,1]))<0))/nrow(AndBOtherRegion)
          
          if(regionSNPseparate_sameAsReference>=regionSNPseparate_diffFromReference)
          {
            PropSNPsSeparate <- c(PropSNPsSeparate,regionSNPseparate_sameAsReference)
            ConcordantWithRef <- c(ConcordantWithRef,1)
          }
          
          if(regionSNPseparate_sameAsReference<regionSNPseparate_diffFromReference)
          {
            PropSNPsSeparate <- c(PropSNPsSeparate,regionSNPseparate_diffFromReference)
            ConcordantWithRef <- c(ConcordantWithRef,-1)
          }
          #PropSNPsSeparate   <- c(PropSNPsSeparate,length(which(c(AndBOtherRegion[MajorColumn1,1],AndBOtherRegion[MajorColumn2,2]) -c(c(AndBOtherRegion[MajorColumn1,2],AndBOtherRegion[MajorColumn2,1]))>0))/nrow(AndBOtherRegion))
          
          
          
        }
        combinedPvals <- pchisq((sum(log(OtherRegionPvals))*-2), df=length(OtherRegionPvals)*2, lower.tail=F)
        regionCompairson[[region]]$RefRegionPval        <- RefRegionTtest$p.value
        regionCompairson[[region]]$RefUnpaired          <- showAI
        regionCompairson[[region]]$combinedPval         <- combinedPvals
        regionCompairson[[region]]$otherPvals           <- OtherRegionPvals
        regionCompairson[[region]]$otherSNPsSeparate    <- PropSNPsSeparate
        regionCompairson[[region]]$ConcordantWithRef    <- ConcordantWithRef
        regionCompairson[[region]]$MajorColumn1         <- MajorColumn1
        regionCompairson[[region]]$MajorColumn2         <- MajorColumn2
        regionCompairson[[region]]$ASCAT_AI             <- ASCATAI
        
        maxPhasingAbility <- c(maxPhasingAbility,mean(PropSNPsSeparate))
        
      }
      if(nrow(logRBAFSample1)<5)
      {
        next;
      }
      names(maxPhasingAbility) <- regionList
      names(showsAI) <- regionList
      names(ASCATcallsAI) <- regionList
      
      maxPhasingAbility <- maxPhasingAbility[showsAI<dipTestPvalThresh|ASCATcallsAI==TRUE]
      showsAI           <- showsAI[showsAI<dipTestPvalThresh|ASCATcallsAI==TRUE]
      
      showsAI_sig       <- showsAI[showsAI<dipTestPvalThresh]
      ASCATcallsAI      <- ASCATcallsAI[ASCATcallsAI==TRUE]
      
      # consider the different cases and use appropriate thresholds. 
      
      # first scenario - there is no evidence of allelic imbalance (from ASCAT or dip.test)
      # in this case -randomly assign SNPs as major or minor
      if(length(maxPhasingAbility)==0)
      {
        
        # now, none of them showed any evidence of BAF separtion. Therefore, make sure it's just a random selection
        #segName <- names_ComparisonOutPut[[x]]
        
        #TOM EDIT######
        #possibleColumns <- c(regionCompairson[[1]]$MajorColumn1,regionCompairson[[1]]$MajorColumn2)
        possibleColumns <- c(newMajorColumn1,newMajorColumn2)
        #TOM EDIT#########End
        
        
        MajorRandom <- sample(possibleColumns,size=length(possibleColumns)/2)
        MinorRandom <- possibleColumns[!possibleColumns%in%MajorRandom]
        
        whichRegionToPhaseFrom[[segName]]$regionToUse <- 'randomSelection'
        whichRegionToPhaseFrom[[segName]]$MajorColumn1 <- sort(MajorRandom)
        whichRegionToPhaseFrom[[segName]]$MajorColumn2 <- sort(MinorRandom)
        whichRegionToPhaseFrom[[segName]]$regionComparisons <- regionCompairson
        # next;
      }
      
      # second scenario - there is evidence of allelic imbalance (from either ASCAT or dip.test [or])
      if(length(maxPhasingAbility)>0)
      {
        maxPhasingAbility <- maxPhasingAbility[order(showsAI)]
        
        # choose a reference
        if(length(ASCATcallsAI)==1)
        {
          regionToChose <- names(maxPhasingAbility[names(ASCATcallsAI)])
          # just use this region as the one to choose.  
          #segName <- names_ComparisonOutPut[[x]]
          whichRegionToPhaseFrom[[segName]]$regionToUse <- regionToChose
          whichRegionToPhaseFrom[[segName]]$MajorColumn1 <- regionCompairson[[regionToChose]]$MajorColumn1
          whichRegionToPhaseFrom[[segName]]$MajorColumn2 <- regionCompairson[[regionToChose]]$MajorColumn2
          whichRegionToPhaseFrom[[segName]]$regionComparisons <- regionCompairson
          
          
        }
        
        # choose a reference
        if(length(ASCATcallsAI)>1)
        {
          regionToChose <- names(which(maxPhasingAbility[ASCATcallsAI]==max(maxPhasingAbility[ASCATcallsAI])))[1]
          # choose a reference
          
          # now if the region isn't >0.75 across all regions, should we combine all
          if(regionCompairson[[regionToChose]]$combinedPval<0.01&mean(regionCompairson[[regionToChose]]$otherSNPsSeparate)<0.75)
          {
            # use a combination of all regions
            logRBAFSample1 <- bigList[[regionToChose]][[x]] #ASCATdata[[regionOne]][ASCATdata[[regionOne]]$chrom==gain.min.cons[x,1]&ASCATdata[[regionOne]]$pos>=gain.min.cons[x,2]&ASCATdata[[regionOne]]$pos<=gain.min.cons[x,3],,drop=FALSE]
            AndBSample1 <- funCalcN_withBAF(logRSites = logRBAFSample1$logR,bafSites = logRBAFSample1$BAF,tumorPloidy = ploidyEstimates[regionToChose],tumorPurity = purityEstimates[regionToChose],gamma = gammaToUse)
            
            nA <- c(AndBSample1[regionCompairson[[regionToChose]]$MajorColumn1,1],AndBSample1[regionCompairson[[regionToChose]]$MajorColumn2,2])
            nB <- c(AndBSample1[regionCompairson[[regionToChose]]$MajorColumn1,2],AndBSample1[regionCompairson[[regionToChose]]$MajorColumn2,1])
            
            OrderToConsider <- order(nA-nB)
            #barplot(c(nA-nB)[OrderToConsider],las=1)
            otherRegions  <- regionList[!regionList%in%regionToChose]
            
            OtherRegionConcordance        <- regionCompairson[[regionToChose]]$ConcordantWithRef
            names(OtherRegionConcordance) <- otherRegions
            RefRegion <- c(1)
            names(RefRegion) <- regionToChose
            OtherRegionConcordance <- c(RefRegion,OtherRegionConcordance)
            
            otherRegions  <- otherRegions[otherRegions%in%c(names(showsAI[showsAI<0.01]))]
            allSigRegions <- c(regionToChose,otherRegions)
            
            #check whether should be considered positive or negative
            MajorColumn1  <- c()
            MajorColumn2  <- c()
            meanSNPscores <- c()
            listScores    <- list()
            changeFromRef <- c()
            
            for (SNP in 1:nrow(AndBSample1))
            {
              columnCurrently <- NA
              if(SNP%in%regionCompairson[[regionToChose]]$MajorColumn1)
              {
                MajorCurrently <- 1
                MinorCurrently <- 2
              }
              
              if(SNP%in%regionCompairson[[regionToChose]]$MajorColumn2)
              {
                MajorCurrently <- 2
                MinorCurrently <-1
                
              }
              
              SNPscores <- c()
              for (region in allSigRegions)
              {
                logRBAFSample1 <- bigList[[region]][[x]] #ASCATdata[[regionOne]][ASCATdata[[regionOne]]$chrom==gain.min.cons[x,1]&ASCATdata[[regionOne]]$pos>=gain.min.cons[x,2]&ASCATdata[[regionOne]]$pos<=gain.min.cons[x,3],,drop=FALSE]
                AndBSample1 <- funCalcN_withBAF(logRSites = logRBAFSample1$logR,bafSites = logRBAFSample1$BAF,tumorPloidy = ploidyEstimates[region],tumorPurity = purityEstimates[region],gamma = gammaToUse)
                
                nA <- AndBSample1[SNP,MajorCurrently]
                nB <- AndBSample1[SNP,MinorCurrently]
                multFactor <- OtherRegionConcordance[region]
                
                DiffScore <- (nA-nB)*multFactor
                
                SNPscores <- c(SNPscores,DiffScore)
                listScores[[region]]$diff <- c(listScores[[region]]$diff,DiffScore)
                
              }
              
              if(mean(SNPscores)>0)
              {
                changeFromRef <- c(changeFromRef,0)
                #MajorToUse <- columnCurrently
                if(MajorCurrently==1)
                {
                  MajorColumn1 <- c(MajorColumn1,SNP)
                }
                
                if(MajorCurrently==2)
                {
                  MajorColumn2 <- c(MajorColumn2,SNP)
                }
                
                
              }
              
              if(mean(SNPscores)<0)
              {
                changeFromRef <- c(changeFromRef,1)
                
                if(MajorCurrently==1)
                {
                  MajorColumn2 <- c(MajorColumn2,SNP)
                }
                
                if(MajorCurrently==2)
                {
                  MajorColumn1 <- c(MajorColumn1,SNP)
                }
              }
              
              meanSNPscores <- c(meanSNPscores,mean(SNPscores))
              
            }
            
            # barplot(listScores$`L_LTX062_SU_T1-R1`$diff)
            # barplot(listScores$`L_LTX062_SU_T1-R2`$diff)
            # barplot(listScores$`L_LTX062_SU_T1-R3`$diff)
            # barplot(listScores$`L_LTX062_SU_T1-R4`$diff)
            
            #segName <- names_ComparisonOutPut[[x]]
            whichRegionToPhaseFrom[[segName]]$regionToUse <- "mixture"
            whichRegionToPhaseFrom[[segName]]$MajorColumn1 <- MajorColumn1
            whichRegionToPhaseFrom[[segName]]$MajorColumn2 <- MajorColumn2
            
            
            #finally, check the ability of the combined to separate A and B in each region
            OtherRegionPvals  <- c()
            PropSNPsSeparate  <- c()
            ConcordantWithRef <- c()
            for (otherRegion in regionList)
            {
              logRBAFOtherRegion <- bigList[[otherRegion]][[x]] #ASCATdata[[regionOne]][ASCATdata[[regionOne]]$chrom==gain.min.cons[x,1]&ASCATdata[[regionOne]]$pos>=gain.min.cons[x,2]&ASCATdata[[regionOne]]$pos<=gain.min.cons[x,3],,drop=FALSE]
              AndBOtherRegion    <- funCalcN_withBAF(logRSites = logRBAFOtherRegion$logR,bafSites = logRBAFOtherRegion$BAF,tumorPloidy = ploidyEstimates[otherRegion],tumorPurity = purityEstimates[otherRegion],gamma = gammaToUse)
              OtherRegionTtest   <- wilcox.test(c(AndBOtherRegion[MajorColumn1,1],AndBOtherRegion[MajorColumn2,2]),c(c(AndBOtherRegion[MajorColumn1,2],AndBOtherRegion[MajorColumn2,1])),paired=TRUE)
              OtherRegionPvals   <- c(OtherRegionPvals,OtherRegionTtest$p.value)
              
              regionSNPseparate_sameAsReference    <- length(which(c(AndBOtherRegion[MajorColumn1,1],AndBOtherRegion[MajorColumn2,2]) -c(c(AndBOtherRegion[MajorColumn1,2],AndBOtherRegion[MajorColumn2,1]))>0))/nrow(AndBOtherRegion)
              regionSNPseparate_diffFromReference  <- length(which(c(AndBOtherRegion[MajorColumn1,1],AndBOtherRegion[MajorColumn2,2]) -c(c(AndBOtherRegion[MajorColumn1,2],AndBOtherRegion[MajorColumn2,1]))<0))/nrow(AndBOtherRegion)
              
              if(regionSNPseparate_sameAsReference>=regionSNPseparate_diffFromReference)
              {
                PropSNPsSeparate <- c(PropSNPsSeparate,regionSNPseparate_sameAsReference)
                ConcordantWithRef <- c(ConcordantWithRef,1)
              }
              
              if(regionSNPseparate_sameAsReference<regionSNPseparate_diffFromReference)
              {
                PropSNPsSeparate <- c(PropSNPsSeparate,regionSNPseparate_diffFromReference)
                ConcordantWithRef <- c(ConcordantWithRef,-1)
              }
              #PropSNPsSeparate   <- c(PropSNPsSeparate,length(which(c(AndBOtherRegion[MajorColumn1,1],AndBOtherRegion[MajorColumn2,2]) -c(c(AndBOtherRegion[MajorColumn1,2],AndBOtherRegion[MajorColumn2,1]))>0))/nrow(AndBOtherRegion))
              
              
              
            }
            
            
            names(PropSNPsSeparate) <- regionList
            names(ConcordantWithRef) <- regionList
            names(OtherRegionPvals)  <- regionList
            whichRegionToPhaseFrom[[segName]]$PropSNPsSeparate <- PropSNPsSeparate
            whichRegionToPhaseFrom[[segName]]$ConcordantWithRef <- ConcordantWithRef
            whichRegionToPhaseFrom[[segName]]$RegionPvals <- OtherRegionPvals
            whichRegionToPhaseFrom[[segName]]$RegionOrder <- regionList
            whichRegionToPhaseFrom[[segName]]$regionComparisons <- regionCompairson
            
            
            
            
          }
          
          if(length(whichRegionToPhaseFrom[[segName]]$regionToUse)==0)
          {
            
            segName <- names_ComparisonOutPut[[x]]
            whichRegionToPhaseFrom[[segName]]$regionToUse <- regionToChose
            whichRegionToPhaseFrom[[segName]]$MajorColumn1 <- regionCompairson[[regionToChose]]$MajorColumn1
            whichRegionToPhaseFrom[[segName]]$MajorColumn2 <- regionCompairson[[regionToChose]]$MajorColumn2
            whichRegionToPhaseFrom[[segName]]$regionComparisons <- regionCompairson
            
          }
          
          
          #if(whichRegionToPhaseFrom)
          
          
        }
        
        
        if(length(whichRegionToPhaseFrom[[segName]]$regionToUse)==0)
        {
          
          # now, none of them showed any evidence of BAF separtion. Therefore, make sure it's just a random selection
          segName <- names_ComparisonOutPut[[x]]
          possibleColumns <- c(regionCompairson[[1]]$MajorColumn1,regionCompairson[[1]]$MajorColumn2)
          MajorRandom <- sample(possibleColumns,size=length(possibleColumns)/2)
          MinorRandom <- possibleColumns[!possibleColumns%in%MajorRandom]
          
          whichRegionToPhaseFrom[[segName]]$regionToUse <- 'randomSelection'
          whichRegionToPhaseFrom[[segName]]$MajorColumn1 <- sort(MajorRandom)
          whichRegionToPhaseFrom[[segName]]$MajorColumn2 <- sort(MinorRandom)
          whichRegionToPhaseFrom[[segName]]$regionComparisons <- regionCompairson
          #next;
        }
        
      }
      
      # next, let's put together a new BAF table (might as well have logR etc in there too)
      for (region in regionList)
      {
        
        
        segName <- names_ComparisonOutPut[[x]]
        regionBAFtable            <- bigList[[region]][[x]] #ASCATdata[[regionOne]][ASCATdata[[regionOne]]$chrom==gain.min.cons[x,1]&ASCATdata[[regionOne]]$pos>=gain.min.cons[x,2]&ASCATdata[[regionOne]]$pos<=gain.min.cons[x,3],,drop=FALSE]
        regionBAFtable$region     <- region
        regionBAFtable$phasedFrom <- whichRegionToPhaseFrom[[segName]]$regionToUse
        regionBAFtable$nA         <- 0
        regionBAFtable$nB         <- 0
        
        if(length(whichRegionToPhaseFrom[[segName]]$MajorColumn1)>1)
        {
          regionBAFtable[whichRegionToPhaseFrom[[segName]]$MajorColumn1,]$nA <- 1
        }
        if(length(whichRegionToPhaseFrom[[segName]]$MajorColumn2)>1)
        {
          regionBAFtable[whichRegionToPhaseFrom[[segName]]$MajorColumn2,]$nB <- 1
        }
        
        if(nrow(regionBAFtable[is.na(regionBAFtable$nA),,drop=FALSE])>0)
        {
          stop('check error')
        }
        
        # add a new BAFsegmented
        #which(regionBAFtable$nB==1)
        
        regionBAFtable$newBAFsegmented <- min(median(regionBAFtable[regionBAFtable$nB==1,]$BAF),median(regionBAFtable[regionBAFtable$nB==0,]$BAF))
        
        newBAFtable[[region]]$regionBAFtable <- rbind(newBAFtable[[region]]$regionBAFtable,regionBAFtable)
      }
      
    }
    
    
    fullOutPut <- list()
    fullOutPut$whichRegionToPhaseFrom <- whichRegionToPhaseFrom
    fullOutPut$newBAFtable <- newBAFtable
    return(fullOutPut)
    #return(whichRegionToPhaseFrom)  
  }
  
}




PasteVector<-function(v,sep=""){
  print("Called PasteVector")
  
  vt <- v[1];
  if(length(v) > 1){
    for(g in 2:length(v)){
      vt <- paste(vt,v[g],sep=sep)
      
    }
  }
  vt <- paste(vt," EnD",sep="");
  out.v <- sub(" EnD","",vt);
  out.v <- sub("NA , ","",out.v);
  out.v <- sub(" , NA","",out.v);
  out.v <- sub(" , NA , "," , ",out.v);
  print("Returning from PasteVector")
  return(out.v);
  
}


funCalcN_withBAF <- function(logRSites,bafSites,tumorPloidy,tumorPurity,gamma)
{
  #print("funCalcN_withBAF")
  nA <- (tumorPurity-1+bafSites*2^(logRSites/gamma)*((1-tumorPurity)*2+tumorPurity*tumorPloidy))/tumorPurity
  nB <- (tumorPurity-1-(bafSites-1)*2^(logRSites/gamma)*((1-tumorPurity)*2+tumorPurity*tumorPloidy))/tumorPurity
  #print("Returning from funCalcN_withBAF")
  return(cbind(nA,nB))
  
}



calcDiffBetweenTwoRegions <- function(regionOne
                                      , regionTwo
                                      , bigList
                                      , min.cons.seg
                                      , gammaToUse=gammaToUse
                                      , plotting=TRUE
                                      , minYaxis=(-2)
                                      , maxYaxis=2
                                      , purityEstimates
                                      , ploidyEstimates
                                      #, phasing='singleRegion'
                                      , whichRegionToPhaseFrom=NULL
                                      , gainThreshold    = log2(2.5/2)
                                      , lossThreshold    = log2(1.5/2)
                                      , ampThreshold     = log2(4/2)
                                      , deepLossThreshold= log2(0.5/2)
)
{
  
  print("Called calcDiffBetweenTwoRegions")
  # Make sure any copy numbers that are too high are still included
  #layout(rbind(1,2))
  
  if(length(whichRegionToPhaseFrom)==0)
  {
    ComparisonOutPut <- list()
    names_ComparisonOutPut <- apply(min.cons.seg[,1:3],1,PasteVector,sep=":")
    HeterogeneityScore <- c()
    SizesOfSegments    <- c()
    for (x in 1:length(names_ComparisonOutPut))
    {
      
      
      
      
      #cat(x)
      # regionToConsider
      logRBAFSample1 <- bigList[[regionOne]][[x]] #ASCATdata[[regionOne]][ASCATdata[[regionOne]]$chrom==gain.min.cons[x,1]&ASCATdata[[regionOne]]$pos>=gain.min.cons[x,2]&ASCATdata[[regionOne]]$pos<=gain.min.cons[x,3],,drop=FALSE]
      logRBAFSample2 <- bigList[[regionTwo]][[x]] #ASCATdata[[regionTwo]][ASCATdata[[regionTwo]]$chrom==gain.min.cons[x,1]&ASCATdata[[regionTwo]]$pos>=gain.min.cons[x,2]&ASCATdata[[regionTwo]]$pos<=gain.min.cons[x,3],,drop=FALSE]
      #logRBAFSample1 <- logRBAFSample1[!is.na(logRBAFSample1$BAFsegmented),,drop=FALSE]
      #logRBAFSample2 <- logRBAFSample2[!is.na(logRBAFSample2$BAFsegmented),,drop=FALSE]
      
      
      
      
      
      if(nrow(logRBAFSample1)<5)
      {
        
        AndBSample1 <- funCalcN_withBAF(logRSites = logRBAFSample1$logR,bafSites = logRBAFSample1$BAF,tumorPloidy = ploidyEstimates[regionOne],tumorPurity = purityEstimates[regionOne],gamma = gammaToUse)
        AndBSample2 <- funCalcN_withBAF(logRSites = logRBAFSample2$logR,bafSites = logRBAFSample2$BAF,tumorPloidy = ploidyEstimates[regionTwo],tumorPurity = purityEstimates[regionTwo],gamma = gammaToUse)
        
        ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionOne]]$chr <- NA
        ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionOne]]$pos <- NA
        ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionOne]]$nA <- NA
        ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionOne]]$nB <- NA
        ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionOne]]$nA_mean <- NA
        ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionOne]]$nB_mean <- NA
        ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$nA_PairedPval        <- NA
        ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$nB_PairedPval        <- NA
        ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$nB_unPairedPval      <- NA
        ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$nA_unPairedPval      <- NA
        ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$combinedPairedPval   <- NA
        ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$combinedUnPairedPval <- NA
        ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$nA_heterogeneity          <- NA
        ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$nB_heterogeneity          <- NA
        ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$combinedHeterogeneity     <- NA
        
        
        #outPut <- cbind(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)
        #return(outPut)
        next;
      }
      
      # make a decision regarding which samples decides A and B
      # this should really be done on all samples together to ensure consistency
      
      AndBSample1 <- funCalcN_withBAF(logRSites = logRBAFSample1$logR,bafSites = logRBAFSample1$BAF,tumorPloidy = ploidyEstimates[regionOne],tumorPurity = purityEstimates[regionOne],gamma = gammaToUse)
      AndBSample2 <- funCalcN_withBAF(logRSites = logRBAFSample2$logR,bafSites = logRBAFSample2$BAF,tumorPloidy = ploidyEstimates[regionTwo],tumorPurity = purityEstimates[regionTwo],gamma = gammaToUse)
      
      Sample1MedianMax <- median(apply(AndBSample1,1,max))
      Sample2MedianMax <- median(apply(AndBSample2,1,max))
      
      Sample1Max <- apply(AndBSample1,1,max)
      Sample2Max <- apply(AndBSample2,1,max)
      tmpSample1 <- ifelse(AndBSample1==Sample1Max,1,0)
      tmpSample2 <- ifelse(AndBSample2==Sample2Max,1,0)
      
      if(identical(tmpSample1,tmpSample2))
      {
        Sample1Classifier <- 0.05
        Sample2Classifier <- 1
      }
      # check if one is significant using the other classification. If so, use that one
      
      if(!identical(tmpSample1,tmpSample2))
      {
        Sample1Classifier <- tTestWithErrorHandling(c(AndBSample2[which(tmpSample1[,1]==1),1],AndBSample2[which(tmpSample1[,2]==1),2]),c(AndBSample2[which(tmpSample1[,1]!=1),2],AndBSample2[which(tmpSample1[,2]!=1),2]))$p.value
        Sample2Classifier <- tTestWithErrorHandling(c(AndBSample1[which(tmpSample2[,1]==1),1],AndBSample1[which(tmpSample2[,2]==1),2]),c(AndBSample1[which(tmpSample2[,1]!=1),2],AndBSample1[which(tmpSample2[,2]!=1),2]))$p.value
        
        
        
        if(is.na(Sample1Classifier))
        {
          Sample1Classifier <- 1
        }
        if(is.na(Sample2Classifier))
        {
          Sample2Classifier <- 1
        }
      }
      
      
      # Major <- apply(cbind(AndBSample1,AndBSample2),1,max)
      # tmpSample1 <- ifelse(AndBSample1==Major,1,0)
      # tmpSample2 <- ifelse(AndBSample2==Major,1,0)
      
      ABSample1 <- AndBSample1
      ABSample2 <- AndBSample2
      
      if(Sample1Classifier<=Sample2Classifier)
      {
        MajorColumn1 <- which(tmpSample1[,1]==1)
        MajorColumn2 <- which(tmpSample1[,2]==1)
      }
      if(Sample1Classifier>Sample2Classifier)
      {
        MajorColumn1 <- which(tmpSample2[,1]==1)
        MajorColumn2 <- which(tmpSample2[,2]==1)
      }
      #Sample1Decides <-  which(cbind(rowSums(tmpSample1),rowSums(tmpSample2))[,1]==1)
      #Sample2Decides <-  which(cbind(rowSums(tmpSample1),rowSums(tmpSample2))[,2]==1)
      
      #MajorColumn1   <- c(Sample1Decides[which(tmpSample1[Sample1Decides,1]==1)],Sample2Decides[which(tmpSample2[Sample2Decides,1]==1)])
      #MajorColumn2   <- c(Sample1Decides[which(tmpSample1[Sample1Decides,2]==1)],Sample2Decides[which(tmpSample2[Sample2Decides,2]==1)])
      
      ABSample1[MajorColumn1,1] <- AndBSample1[MajorColumn1,1]
      ABSample1[MajorColumn2,1] <- AndBSample1[MajorColumn2,2]
      ABSample1[MajorColumn1,2] <- AndBSample1[MajorColumn1,2]
      ABSample1[MajorColumn2,2] <- AndBSample1[MajorColumn2,1]
      
      ABSample2[MajorColumn1,1] <- AndBSample2[MajorColumn1,1]
      ABSample2[MajorColumn2,1] <- AndBSample2[MajorColumn2,2]
      ABSample2[MajorColumn1,2] <- AndBSample2[MajorColumn1,2]
      ABSample2[MajorColumn2,2] <- AndBSample2[MajorColumn2,1]
      
      
      
      PValsPaired          <- c(tTestWithErrorHandling(ABSample1[,1],ABSample2[,1],paired=TRUE)$p.val,tTestWithErrorHandling(ABSample1[,2],ABSample2[,2],paired=TRUE)$p.val)
      combinedPvalPaired   <- pchisq((sum(log(PValsPaired))*-2), df=length(PValsPaired)*2, lower.tail=F)
      PValsUnPaired        <- c(tTestWithErrorHandling(ABSample1[,1],ABSample2[,1])$p.val,tTestWithErrorHandling(ABSample1[,2],ABSample2[,2])$p.val)
      combinedPvalUnPaired <- pchisq((sum(log(PValsUnPaired))*-2), df=length(PValsUnPaired)*2, lower.tail=F)
      
      # also get a p-value that relates to both samples combined
      totalCPNpval         <- tTestWithErrorHandling(c(ABSample1[,1]+ABSample1[,2]),c(ABSample2[,1]+ABSample2[,2]),paired=TRUE)$p.val
      
      nA_heterogeneity       <- 0
      nB_heterogeneity       <- 0
      SizeOfSegment           <- strsplit(names_ComparisonOutPut[x],split=":")
      SizeOfSegment           <- as.numeric(unlist(SizeOfSegment)[3])-as.numeric(unlist(SizeOfSegment)[2])
      
      if(PValsPaired[1]<0.001)
      {
        nAamplitudeOfDifference <- abs(mean(as.numeric(ABSample1[,1]))-mean(as.numeric(ABSample2[,1])))
        SizeOfSegment           <- strsplit(names_ComparisonOutPut[x],split=":")
        SizeOfSegment           <- as.numeric(unlist(SizeOfSegment)[3])-as.numeric(unlist(SizeOfSegment)[2])
        nA_heterogeneity        <- nAamplitudeOfDifference*SizeOfSegment
      }
      
      if(PValsPaired[2]<0.001)
      {
        nBamplitudeOfDifference <- abs(mean(as.numeric(ABSample1[,2]))-mean(as.numeric(ABSample2[,2])))
        SizeOfSegment           <- strsplit(names_ComparisonOutPut[x],split=":")
        SizeOfSegment           <- as.numeric(unlist(SizeOfSegment)[3])-as.numeric(unlist(SizeOfSegment)[2])
        nB_heterogeneity        <- nBamplitudeOfDifference*SizeOfSegment
      }
      
      total_heterogeneity <- nA_heterogeneity+nB_heterogeneity
      HeterogeneityScore  <- c(HeterogeneityScore,total_heterogeneity)
      SizesOfSegments     <- c(SizesOfSegments,SizeOfSegment)
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionOne]]$chr     <- logRBAFSample1$chrom
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionOne]]$pos     <- logRBAFSample1$pos
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionOne]]$nA      <- as.numeric(ABSample1[,1])
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionOne]]$nB      <- as.numeric(ABSample1[,2])
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionOne]]$nA_mean <- mean(as.numeric(ABSample1[,1]))
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionOne]]$nB_mean <- mean(as.numeric(ABSample1[,2]))
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionTwo]]$nA      <- as.numeric(ABSample2[,1])
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionTwo]]$nB      <- as.numeric(ABSample2[,2])
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionTwo]]$nA_mean <- mean(as.numeric(ABSample2[,1]))
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionTwo]]$nB_mean <- mean(as.numeric(ABSample2[,2]))
      
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$nA_PairedPval             <- PValsPaired[1]
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$nB_PairedPval             <- PValsPaired[2]
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$nB_unPairedPval           <- PValsUnPaired[2]
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$nA_unPairedPval           <- PValsUnPaired[1]
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$combinedPairedPval        <- combinedPvalPaired
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$combinedUnPairedPval      <- combinedPvalUnPaired
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$combinedUnPairedPval      <- combinedPvalUnPaired
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$totalCPNPairedPval        <- totalCPNpval
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$nA_heterogeneity          <- nA_heterogeneity
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$nB_heterogeneity          <- nB_heterogeneity
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$combinedHeterogeneity     <- nB_heterogeneity+nA_heterogeneity
      # classify the segment as showing MSAI
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$msAI                      <- 0
      # if(PValsPaired[1]<0.001&PValsPaired[2]<0.001)
      # {
      
      if(tTestWithErrorHandling(ABSample1[,1],ABSample1[,2],paired=TRUE)$p.val<0.001&tTestWithErrorHandling(ABSample2[,1],ABSample2[,2],paired=TRUE)$p.val<0.001)
      {
        if((mean(ABSample1[,1]-ABSample1[,2])>0&mean(ABSample2[,1]-ABSample2[,2])<0)|(mean(ABSample1[,1]-ABSample1[,2])<0&mean(ABSample2[,1]-ABSample2[,2])>0))
        {
          ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$msAI                      <- 1
        }
        
      }
      
      
      checkGainSignificance <- function(logRBAFRegion,regionToCheck,ploidyEstimates,gainThreshold)
      {
        minGainRawExpRegion       <-  ploidyEstimates[regionToCheck]*(2^(gainThreshold))
        minGainLogRExpRegion      <-  funCalcExpectedLogR(minGainRawExpRegion,0,ploidyEstimates[regionToCheck],purityEstimates[regionToCheck])
        logRObs                   <-  mean(logRBAFRegion$logR)
        
        if(nrow(logRBAFRegion)>=3)
        {
          Ttest        <- tTestWithErrorHandling(logRBAFRegion$logR,mu = minGainLogRExpRegion)
          TtestGreater <- tTestWithErrorHandling(logRBAFRegion$logR,mu = minGainLogRExpRegion,alternative = 'greater')
          #TtestLess    <- tTestWithErrorHandling(logRBAF$logR,mu = minLossLogRExp,alternative = 'less')
          logRObs_low  <- Ttest$conf.int[1]
          logRObs_high <- Ttest$conf.int[2]
          PVal         <- Ttest$p.value
          PvalGain     <- TtestGreater$p.value
          
        }
        if(nrow(logRBAFRegion)<3)
        {
          logRObs_low  <- NA
          logRObs_high <- NA
          PVal         <- NA
          PvalGain     <- NA
          
          
        }
        
        out <- cbind(minGainRawExpRegion,minGainLogRExpRegion,logRObs,logRObs_low,logRObs_high,PVal,PvalGain)  
        return(out)
        
        
      }
      checkLossSignificance <- function(logRBAFRegion,regionToCheck,ploidyEstimates,lossThreshold)
      {
        minLossRawExpRegion       <-  ploidyEstimates[regionToCheck]*(2^(lossThreshold))
        minLossLogRExpRegion      <-  funCalcExpectedLogR(minLossRawExpRegion,0,ploidyEstimates[regionToCheck],purityEstimates[regionToCheck])
        logRObs                   <- mean(logRBAFRegion$logR)
        
        if(nrow(logRBAFRegion)>=3)
        {
          Ttest        <- tTestWithErrorHandling(logRBAFRegion$logR,mu = minLossLogRExpRegion)
          TtestLess <- tTestWithErrorHandling(logRBAFRegion$logR,mu = minLossLogRExpRegion,alternative = 'less')
          #TtestLess    <- tTestWithErrorHandling(logRBAF$logR,mu = minLossLogRExp,alternative = 'less')
          logRObs_low  <- Ttest$conf.int[1]
          logRObs_high <- Ttest$conf.int[2]
          PVal         <- Ttest$p.value
          PvalLoss     <- TtestLess$p.value
          
        }
        # if(nrow(logRBAF)<3)
        # {
        #   logRObs_low  <- NA
        #   logRObs_high <- NA
        #   PVal         <- NA
        #   PvalLoss     <- NA
        
        
        # }
        
        out <- cbind(minLossRawExpRegion,minLossLogRExpRegion,logRObs,logRObs_low,logRObs_high,PVal,PvalLoss)  
        return(out)
        
        
      }
      
      
      GainRegionOne     <- checkGainSignificance(logRBAFSample1,regionToCheck = regionOne,ploidyEstimates,gainThreshold = gainThreshold)
      AmpRegionOne      <- checkGainSignificance(logRBAFSample1,regionToCheck = regionOne,ploidyEstimates,gainThreshold = ampThreshold)
      LossRegionOne     <- checkLossSignificance(logRBAFSample1,regionToCheck = regionOne,ploidyEstimates,lossThreshold = lossThreshold)
      DeepLossRegionOne <- checkLossSignificance(logRBAFSample1,regionToCheck = regionOne,ploidyEstimates,lossThreshold = deepLossThreshold)
      
      
      # put together into One table
      regionOneGainOrLoss <- c(GainRegionOne[,3:5,drop=FALSE]
                               ,GainRegionOne[,2,drop=FALSE],GainRegionOne[,7,drop=FALSE]
                               ,AmpRegionOne[,2,drop=FALSE],AmpRegionOne[,7,drop=FALSE]
                               ,LossRegionOne[,2,drop=FALSE],LossRegionOne[,7,drop=FALSE]
                               ,DeepLossRegionOne[,2,drop=FALSE],DeepLossRegionOne[,7,drop=FALSE]
      )
      names(regionOneGainOrLoss) <- c('logRObs','logRObs_low','logRObs_high'
                                      ,'minGainExpLogR','gainPval'
                                      ,'minAmpExpLogR','ampPval'
                                      ,'minLossExpLogR','lossPval'
                                      ,'minDeepLossExpLogR','DeepLossPval')
      
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionOne]]$segStatus     <- regionOneGainOrLoss
      
      
      GainRegionTwo     <- checkGainSignificance(logRBAFSample2,regionToCheck = regionTwo,ploidyEstimates,gainThreshold = gainThreshold)
      AmpRegionTwo      <- checkGainSignificance(logRBAFSample2,regionToCheck = regionTwo,ploidyEstimates,gainThreshold = ampThreshold)
      LossRegionTwo     <- checkLossSignificance(logRBAFSample2,regionToCheck = regionTwo,ploidyEstimates,lossThreshold = lossThreshold)
      DeepLossRegionTwo <- checkLossSignificance(logRBAFSample2,regionToCheck = regionTwo,ploidyEstimates,lossThreshold = deepLossThreshold)
      
      
      # put together into Two table
      regionTwoGainOrLoss <- c(GainRegionTwo[,3:5,drop=FALSE]
                               ,GainRegionTwo[,2,drop=FALSE],GainRegionTwo[,7,drop=FALSE]
                               ,AmpRegionTwo[,2,drop=FALSE],AmpRegionTwo[,7,drop=FALSE]
                               ,LossRegionTwo[,2,drop=FALSE],LossRegionTwo[,7,drop=FALSE]
                               ,DeepLossRegionTwo[,2,drop=FALSE],DeepLossRegionTwo[,7,drop=FALSE]
      )
      names(regionTwoGainOrLoss) <- c('logRObs','logRObs_low','logRObs_high'
                                      ,'minGainExpLogR','gainPval'
                                      ,'minAmpExpLogR','ampPval'
                                      ,'minLossExpLogR','lossPval'
                                      ,'minDeepLossExpLogR','DeepLossPval')
      
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionTwo]]$segStatus     <- regionTwoGainOrLoss
      
      
      
      # # for it to be classified as MSAI the difference between A1-A2 must be positive and B1-B2 negative (or the other way round)
      # if((mean(as.numeric(ABSample1[,1]))-mean(as.numeric(ABSample2[,1]))>0&mean(as.numeric(ABSample1[,2]))-mean(as.numeric(ABSample2[,2]))<0)|(mean(as.numeric(ABSample1[,1]))-mean(as.numeric(ABSample2[,1]))<0&mean(as.numeric(ABSample1[,2]))-mean(as.numeric(ABSample2[,2]))>0))
      # {
      #   # we also must observe significant difference in A1 and B1 as well as A2 and B2. And this significant difference must be opposite. 
      #   
      # }
      # }
      
      
      
      
    }
    
    averageHeterogeneity <- sum(HeterogeneityScore)/sum(SizesOfSegments)
    
    #ComparisonOutPut$`1:762273:11831543`$`L_LTX064_SU_T1-R6:L_LTX064_SU_T1-R7`$nA_heterogeneity
    #ComparisonOutPut$`1:762273:11831543`$`L_LTX064_SU_T1-R6:L_LTX064_SU_T1-R7`$nB_heterogeneity
    ComparisonOutPut[[paste(regionOne,regionTwo,sep=":")]]$HeterogeneityScore <- averageHeterogeneity 
    
  }
  
  if(length(whichRegionToPhaseFrom)!=0)
  {
    ComparisonOutPut <- list()
    names_ComparisonOutPut <- apply(min.cons.seg[,1:3],1,PasteVector,sep=":")
    HeterogeneityScore <- c()
    SizesOfSegments    <- c()
    for (x in 1:length(names_ComparisonOutPut))
    {
      
      #cat(x)
      # regionToConsider
      NameToUse <- names_ComparisonOutPut[x]
      
      
      logRBAFSample1 <- bigList[[regionOne]][[x]] #ASCATdata[[regionOne]][ASCATdata[[regionOne]]$chrom==gain.min.cons[x,1]&ASCATdata[[regionOne]]$pos>=gain.min.cons[x,2]&ASCATdata[[regionOne]]$pos<=gain.min.cons[x,3],,drop=FALSE]
      logRBAFSample2 <- bigList[[regionTwo]][[x]] #ASCATdata[[regionTwo]][ASCATdata[[regionTwo]]$chrom==gain.min.cons[x,1]&ASCATdata[[regionTwo]]$pos>=gain.min.cons[x,2]&ASCATdata[[regionTwo]]$pos<=gain.min.cons[x,3],,drop=FALSE]
      #logRBAFSample1 <- logRBAFSample1[!is.na(logRBAFSample1$BAFsegmented),,drop=FALSE]
      #logRBAFSample2 <- logRBAFSample2[!is.na(logRBAFSample2$BAFsegmented),,drop=FALSE]
      
      if(nrow(logRBAFSample1)<5)
      {
        
        AndBSample1 <- funCalcN_withBAF(logRSites = logRBAFSample1$logR,bafSites = logRBAFSample1$BAF,tumorPloidy = ploidyEstimates[regionOne],tumorPurity = purityEstimates[regionOne],gamma = gammaToUse)
        AndBSample2 <- funCalcN_withBAF(logRSites = logRBAFSample2$logR,bafSites = logRBAFSample2$BAF,tumorPloidy = ploidyEstimates[regionTwo],tumorPurity = purityEstimates[regionTwo],gamma = gammaToUse)
        
        ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionOne]]$chr <- NA
        ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionOne]]$pos <- NA
        ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionOne]]$nA <- NA
        ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionOne]]$nB <- NA
        ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionOne]]$nA_mean <- NA
        ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionOne]]$nB_mean <- NA
        ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$nA_PairedPval        <- NA
        ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$nB_PairedPval        <- NA
        ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$nB_unPairedPval      <- NA
        ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$nA_unPairedPval      <- NA
        ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$combinedPairedPval   <- NA
        ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$combinedUnPairedPval <- NA
        ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$nA_heterogeneity          <- NA
        ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$nB_heterogeneity          <- NA
        ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$combinedHeterogeneity     <- NA
        
        
        #outPut <- cbind(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)
        #return(outPut)
        next;
      }
      
      # make a decision regarding which samples decides A and B
      # this should really be done on all samples together to ensure consistency
      
      AndBSample1 <- funCalcN_withBAF(logRSites = logRBAFSample1$logR,bafSites = logRBAFSample1$BAF,tumorPloidy = ploidyEstimates[regionOne],tumorPurity = purityEstimates[regionOne],gamma = gammaToUse)
      AndBSample2 <- funCalcN_withBAF(logRSites = logRBAFSample2$logR,bafSites = logRBAFSample2$BAF,tumorPloidy = ploidyEstimates[regionTwo],tumorPurity = purityEstimates[regionTwo],gamma = gammaToUse)
      
      Sample1MedianMax <- median(apply(AndBSample1,1,max))
      Sample2MedianMax <- median(apply(AndBSample2,1,max))
      
      Sample1Max <- apply(AndBSample1,1,max)
      Sample2Max <- apply(AndBSample2,1,max)
      tmpSample1 <- ifelse(AndBSample1==Sample1Max,1,0)
      tmpSample2 <- ifelse(AndBSample2==Sample2Max,1,0)
      
      if(identical(tmpSample1,tmpSample2))
      {
        Sample1Classifier <- 0.05
        Sample2Classifier <- 1
      }
      # check if one is significant using the other classification. If so, use that one
      
      if(!identical(tmpSample1,tmpSample2))
      {
        Sample1Classifier <- tTestWithErrorHandling(c(AndBSample2[which(tmpSample1[,1]==1),1],AndBSample2[which(tmpSample1[,2]==1),2]),c(AndBSample2[which(tmpSample1[,1]!=1),2],AndBSample2[which(tmpSample1[,2]!=1),2]))$p.value
        Sample2Classifier <- tTestWithErrorHandling(c(AndBSample1[which(tmpSample2[,1]==1),1],AndBSample1[which(tmpSample2[,2]==1),2]),c(AndBSample1[which(tmpSample2[,1]!=1),2],AndBSample1[which(tmpSample2[,2]!=1),2]))$p.value
        
        
        
        if(is.na(Sample1Classifier))
        {
          Sample1Classifier <- 1
        }
        if(is.na(Sample2Classifier)) 
        {
          Sample2Classifier <- 1
        }
      }
      
      
      # Major <- apply(cbind(AndBSample1,AndBSample2),1,max)
      # tmpSample1 <- ifelse(AndBSample1==Major,1,0)
      # tmpSample2 <- ifelse(AndBSample2==Major,1,0)
      
      ABSample1 <- AndBSample1
      ABSample2 <- AndBSample2
      
      if(Sample1Classifier<=Sample2Classifier)
      {
        MajorColumn1 <- which(tmpSample1[,1]==1)
        MajorColumn2 <- which(tmpSample1[,2]==1)
      }
      if(Sample1Classifier>Sample2Classifier)
      {
        MajorColumn1 <- which(tmpSample2[,1]==1)
        MajorColumn2 <- which(tmpSample2[,2]==1)
      }
      
      if(length(whichRegionToPhaseFrom[[NameToUse]]$MajorColumn1)>0)
      {
        MajorColumn1 <- whichRegionToPhaseFrom[[NameToUse]]$MajorColumn1
        MajorColumn2 <- whichRegionToPhaseFrom[[NameToUse]]$MajorColumn2
        
      }
      
      #Sample1Decides <-  which(cbind(rowSums(tmpSample1),rowSums(tmpSample2))[,1]==1)
      #Sample2Decides <-  which(cbind(rowSums(tmpSample1),rowSums(tmpSample2))[,2]==1)
      
      #MajorColumn1   <- c(Sample1Decides[which(tmpSample1[Sample1Decides,1]==1)],Sample2Decides[which(tmpSample2[Sample2Decides,1]==1)])
      #MajorColumn2   <- c(Sample1Decides[which(tmpSample1[Sample1Decides,2]==1)],Sample2Decides[which(tmpSample2[Sample2Decides,2]==1)])
      
      ABSample1[MajorColumn1,1] <- AndBSample1[MajorColumn1,1]
      ABSample1[MajorColumn2,1] <- AndBSample1[MajorColumn2,2]
      ABSample1[MajorColumn1,2] <- AndBSample1[MajorColumn1,2]
      ABSample1[MajorColumn2,2] <- AndBSample1[MajorColumn2,1]
      
      ABSample2[MajorColumn1,1] <- AndBSample2[MajorColumn1,1]
      ABSample2[MajorColumn2,1] <- AndBSample2[MajorColumn2,2]
      ABSample2[MajorColumn1,2] <- AndBSample2[MajorColumn1,2]
      ABSample2[MajorColumn2,2] <- AndBSample2[MajorColumn2,1]
      
      
      
      PValsPaired          <- c(tTestWithErrorHandling(ABSample1[,1],ABSample2[,1],paired=TRUE)$p.val,tTestWithErrorHandling(ABSample1[,2],ABSample2[,2],paired=TRUE)$p.val)
      combinedPvalPaired   <- pchisq((sum(log(PValsPaired))*-2), df=length(PValsPaired)*2, lower.tail=F)
      PValsUnPaired        <- c(tTestWithErrorHandling(ABSample1[,1],ABSample2[,1])$p.val,tTestWithErrorHandling(ABSample1[,2],ABSample2[,2])$p.val)
      combinedPvalUnPaired <- pchisq((sum(log(PValsUnPaired))*-2), df=length(PValsUnPaired)*2, lower.tail=F)
      
      # also get a p-value that relates to both samples combined
      totalCPNpval         <- tTestWithErrorHandling(c(ABSample1[,1]+ABSample1[,2]),c(ABSample2[,1]+ABSample2[,2]),paired=TRUE)$p.val
      
      nA_heterogeneity       <- 0
      nB_heterogeneity       <- 0
      SizeOfSegment           <- strsplit(names_ComparisonOutPut[x],split=":")
      SizeOfSegment           <- as.numeric(unlist(SizeOfSegment)[3])-as.numeric(unlist(SizeOfSegment)[2])
      
      if(PValsPaired[1]<0.001)
      {
        nAamplitudeOfDifference <- abs(mean(as.numeric(ABSample1[,1]))-mean(as.numeric(ABSample2[,1])))
        SizeOfSegment           <- strsplit(names_ComparisonOutPut[x],split=":")
        SizeOfSegment           <- as.numeric(unlist(SizeOfSegment)[3])-as.numeric(unlist(SizeOfSegment)[2])
        nA_heterogeneity        <- nAamplitudeOfDifference*SizeOfSegment
      }
      
      if(PValsPaired[2]<0.001)
      {
        nBamplitudeOfDifference <- abs(mean(as.numeric(ABSample1[,2]))-mean(as.numeric(ABSample2[,2])))
        SizeOfSegment           <- strsplit(names_ComparisonOutPut[x],split=":")
        SizeOfSegment           <- as.numeric(unlist(SizeOfSegment)[3])-as.numeric(unlist(SizeOfSegment)[2])
        nB_heterogeneity        <- nBamplitudeOfDifference*SizeOfSegment
      }
      
      total_heterogeneity <- nA_heterogeneity+nB_heterogeneity
      HeterogeneityScore  <- c(HeterogeneityScore,total_heterogeneity)
      SizesOfSegments     <- c(SizesOfSegments,SizeOfSegment)
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionOne]]$chr     <- logRBAFSample1$chrom
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionOne]]$pos     <- logRBAFSample1$pos
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionOne]]$nA      <- as.numeric(ABSample1[,1])
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionOne]]$nB      <- as.numeric(ABSample1[,2])
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionOne]]$nA_mean <- mean(as.numeric(ABSample1[,1]))
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionOne]]$nB_mean <- mean(as.numeric(ABSample1[,2]))
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionTwo]]$nA      <- as.numeric(ABSample2[,1])
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionTwo]]$nB      <- as.numeric(ABSample2[,2])
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionTwo]]$nA_mean <- mean(as.numeric(ABSample2[,1]))
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionTwo]]$nB_mean <- mean(as.numeric(ABSample2[,2]))
      
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$nA_PairedPval             <- PValsPaired[1]
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$nB_PairedPval             <- PValsPaired[2]
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$nB_unPairedPval           <- PValsUnPaired[2]
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$nA_unPairedPval           <- PValsUnPaired[1]
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$combinedPairedPval        <- combinedPvalPaired
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$combinedUnPairedPval      <- combinedPvalUnPaired
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$combinedUnPairedPval      <- combinedPvalUnPaired
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$totalCPNPairedPval        <- totalCPNpval
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$nA_heterogeneity          <- nA_heterogeneity
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$nB_heterogeneity          <- nB_heterogeneity
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$combinedHeterogeneity     <- nB_heterogeneity+nA_heterogeneity
      # classify the segment as showing MSAI
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$msAI                      <- 0
      # if(PValsPaired[1]<0.001&PValsPaired[2]<0.001)
      # {
      
      if(tTestWithErrorHandling(ABSample1[,1],ABSample1[,2],paired=TRUE)$p.val<0.001&tTestWithErrorHandling(ABSample2[,1],ABSample2[,2],paired=TRUE)$p.val<0.001)
      {
        if((mean(ABSample1[,1]-ABSample1[,2])>0&mean(ABSample2[,1]-ABSample2[,2])<0)|(mean(ABSample1[,1]-ABSample1[,2])<0&mean(ABSample2[,1]-ABSample2[,2])>0))
        {
          ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$msAI                      <- 1
        }
        
      }
      
      
      
      checkGainSignificance <- function(logRBAFRegion,regionToCheck,ploidyEstimates,gainThreshold)
      {
        minGainRawExpRegion       <-  ploidyEstimates[regionToCheck]*(2^(gainThreshold))
        minGainLogRExpRegion      <-  funCalcExpectedLogR(minGainRawExpRegion,0,ploidyEstimates[regionToCheck],purityEstimates[regionToCheck])
        logRObs                   <- mean(logRBAFRegion$logR)
        
        if(nrow(logRBAFRegion)>=3)
        {
          Ttest        <- tTestWithErrorHandling(logRBAFRegion$logR,mu = minGainLogRExpRegion)
          TtestGreater <- tTestWithErrorHandling(logRBAFRegion$logR,mu = minGainLogRExpRegion,alternative = 'greater')
          #TtestLess    <- tTestWithErrorHandling(logRBAF$logR,mu = minLossLogRExp,alternative = 'less')
          logRObs_low  <- Ttest$conf.int[1]
          logRObs_high <- Ttest$conf.int[2]
          PVal         <- Ttest$p.value
          PvalGain     <- TtestGreater$p.value
          
        }
        if(nrow(logRBAFRegion)<3)
        {
          logRObs_low  <- NA
          logRObs_high <- NA
          PVal         <- NA
          PvalGain     <- NA
          
          
        }
        
        out <- cbind(minGainRawExpRegion,minGainLogRExpRegion,logRObs,logRObs_low,logRObs_high,PVal,PvalGain)  
        return(out)
        
        
      }
      checkLossSignificance <- function(logRBAFRegion,regionToCheck,ploidyEstimates,lossThreshold)
      {
        minLossRawExpRegion       <-  ploidyEstimates[regionToCheck]*(2^(lossThreshold))
        minLossLogRExpRegion      <-  funCalcExpectedLogR(minLossRawExpRegion,0,ploidyEstimates[regionToCheck],purityEstimates[regionToCheck])
        logRObs                   <- mean(logRBAFRegion$logR)
        
        if(nrow(logRBAFRegion)>=3)
        {
          Ttest        <- tTestWithErrorHandling(logRBAFRegion$logR,mu = minLossLogRExpRegion)
          TtestLess <- tTestWithErrorHandling(logRBAFRegion$logR,mu = minLossLogRExpRegion,alternative = 'less')
          #TtestLess    <- tTestWithErrorHandling(logRBAF$logR,mu = minLossLogRExp,alternative = 'less')
          logRObs_low  <- Ttest$conf.int[1]
          logRObs_high <- Ttest$conf.int[2]
          PVal         <- Ttest$p.value
          PvalLoss     <- TtestLess$p.value
          
        }
        if(nrow(logRBAFRegion)<3)
        {
          logRObs_low  <- NA
          logRObs_high <- NA
          PVal         <- NA
          PvalLoss     <- NA
          
          
        }
        
        out <- cbind(minLossRawExpRegion,minLossLogRExpRegion,logRObs,logRObs_low,logRObs_high,PVal,PvalLoss)  
        return(out)
        
        
      }
      
      
      GainRegionOne     <- checkGainSignificance(logRBAFSample1,regionToCheck = regionOne,ploidyEstimates,gainThreshold = gainThreshold)
      AmpRegionOne      <- checkGainSignificance(logRBAFSample1,regionToCheck = regionOne,ploidyEstimates,gainThreshold = ampThreshold)
      LossRegionOne     <- checkLossSignificance(logRBAFSample1,regionToCheck = regionOne,ploidyEstimates,lossThreshold = lossThreshold)
      DeepLossRegionOne <- checkLossSignificance(logRBAFSample1,regionToCheck = regionOne,ploidyEstimates,lossThreshold = deepLossThreshold)
      
      
      # put together into One table
      regionOneGainOrLoss <- c(GainRegionOne[,3:5,drop=FALSE]
                               ,GainRegionOne[,2,drop=FALSE],GainRegionOne[,7,drop=FALSE]
                               ,AmpRegionOne[,2,drop=FALSE],AmpRegionOne[,7,drop=FALSE]
                               ,LossRegionOne[,2,drop=FALSE],LossRegionOne[,7,drop=FALSE]
                               ,DeepLossRegionOne[,2,drop=FALSE],DeepLossRegionOne[,7,drop=FALSE]
      )
      names(regionOneGainOrLoss) <- c('logRObs','logRObs_low','logRObs_high'
                                      ,'minGainExpLogR','gainPval'
                                      ,'minAmpExpLogR','ampPval'
                                      ,'minLossExpLogR','lossPval'
                                      ,'minDeepLossExpLogR','DeepLossPval')
      
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionOne]]$segStatus     <- regionOneGainOrLoss
      
      
      GainRegionTwo     <- checkGainSignificance(logRBAFSample2,regionToCheck = regionTwo,ploidyEstimates,gainThreshold = gainThreshold)
      AmpRegionTwo      <- checkGainSignificance(logRBAFSample2,regionToCheck = regionTwo,ploidyEstimates,gainThreshold = ampThreshold)
      LossRegionTwo     <- checkLossSignificance(logRBAFSample2,regionToCheck = regionTwo,ploidyEstimates,lossThreshold = lossThreshold)
      DeepLossRegionTwo <- checkLossSignificance(logRBAFSample2,regionToCheck = regionTwo,ploidyEstimates,lossThreshold = deepLossThreshold)
      
      
      # put together into Two table
      regionTwoGainOrLoss <- c(GainRegionTwo[,3:5,drop=FALSE]
                               ,GainRegionTwo[,2,drop=FALSE],GainRegionTwo[,7,drop=FALSE]
                               ,AmpRegionTwo[,2,drop=FALSE],AmpRegionTwo[,7,drop=FALSE]
                               ,LossRegionTwo[,2,drop=FALSE],LossRegionTwo[,7,drop=FALSE]
                               ,DeepLossRegionTwo[,2,drop=FALSE],DeepLossRegionTwo[,7,drop=FALSE]
      )
      names(regionTwoGainOrLoss) <- c('logRObs','logRObs_low','logRObs_high'
                                      ,'minGainExpLogR','gainPval'
                                      ,'minAmpExpLogR','ampPval'
                                      ,'minLossExpLogR','lossPval'
                                      ,'minDeepLossExpLogR','DeepLossPval')
      
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionTwo]]$segStatus     <- regionTwoGainOrLoss
      
      # # for it to be classified as MSAI the difference between A1-A2 must be positive and B1-B2 negative (or the other way round)
      # if((mean(as.numeric(ABSample1[,1]))-mean(as.numeric(ABSample2[,1]))>0&mean(as.numeric(ABSample1[,2]))-mean(as.numeric(ABSample2[,2]))<0)|(mean(as.numeric(ABSample1[,1]))-mean(as.numeric(ABSample2[,1]))<0&mean(as.numeric(ABSample1[,2]))-mean(as.numeric(ABSample2[,2]))>0))
      # {
      #   # we also must observe significant difference in A1 and B1 as well as A2 and B2. And this significant difference must be opposite. 
      #   
      # }
      # }
      
      
      
      
    }
    
    averageHeterogeneity <- sum(HeterogeneityScore)/sum(SizesOfSegments)
    
    #ComparisonOutPut$`1:762273:11831543`$`L_LTX064_SU_T1-R6:L_LTX064_SU_T1-R7`$nA_heterogeneity
    #ComparisonOutPut$`1:762273:11831543`$`L_LTX064_SU_T1-R6:L_LTX064_SU_T1-R7`$nB_heterogeneity
    ComparisonOutPut[[paste(regionOne,regionTwo,sep=":")]]$HeterogeneityScore <- averageHeterogeneity 
    
  }
  
  
  if(plotting==TRUE)
  {
    
    par(mfrow=c(4,1))
    
    #major Allele####
    chrom.length.copy   <- fun.chrom.length(min.cons.seg[,1:3])
    names(chrom.length.copy)   <- 1:22
    chrom.segs          <- fun.add.chrom(min.cons.seg[,1:3],chrom.length.copy)
    min.x               <- 0
    min.y               <- minYaxis
    max.y               <- maxYaxis
    max.x               <- as.numeric(max(chrom.segs[,'End']))
    par(mar=c(2,5,5,5))
    plot(1
         ,xlim=c(min.x,max.x)
         ,ylim=c(min.y,max.y)
         ,type='n'
         ,xaxt='n'
         ,yaxs='i'
         ,yaxt='n'
         ,xaxs='i'
         ,ylab=""
         ,xlab=""
         ,lwd=2
         ,main=paste("Major allele ", regionOne," vs ",regionTwo,sep="")
         ,bty="n")
    box()
    
    for (x in 1:length(names_ComparisonOutPut))
    {
      ChrCoords <- unlist(strsplit(names_ComparisonOutPut[x],split=":"))
      toPlot <- cbind(ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionOne]]$chr
                      ,ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionOne]]$pos
                      ,ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionOne]]$pos
                      ,ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionOne]]$nA
                      ,ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionOne]]$nB
                      ,ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionTwo]]$nA
                      ,ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionTwo]]$nB)
      if(is.na(toPlot[1,1]))
      {
        next;
      }
      
      #toPlot <- cbind(logRBAFSample1$chrom,logRBAFSample1$pos,logRBAFSample1$pos,ABSample1,ABSample2)
      toPlot[,1:3] <- fun.add.chrom(apply(toPlot[,1:3],2,as.numeric),chrom.length.copy)
      toPlot <- data.frame(toPlot,stringsAsFactors = FALSE)
      colnames(toPlot) <- c('chr','pos','pos','Sample1nA','Sample1nB','Sample2nA','Sample2nB')
      PValsPaired <- as.numeric(ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$nA_PairedPval)
      totalPval   <-  ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$totalCPNPairedPval#as.numeric(ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$nA_PairedPval)
      
      
      MSAIstatus  <- as.numeric(ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$msAI)
      points(x = as.numeric(toPlot$pos),y=as.numeric(toPlot$Sample1nA)-as.numeric(toPlot$Sample2nA),col=ifelse(MSAIstatus==1,'#238b4550',ifelse(PValsPaired<0.001&totalPval<0.001,'#3182bd50','#fc927250')),cex=0.75,pch=16)
      par(lend=1)
      segments(x0=min(as.numeric(toPlot$pos))
               ,y0=mean(c(as.numeric(toPlot$Sample1nA)-as.numeric(toPlot$Sample2nA)))
               ,x1=max(as.numeric(toPlot$pos)),y1=mean(c(as.numeric(toPlot$Sample1nA)-as.numeric(toPlot$Sample2nA)))
               ,col=ifelse(MSAIstatus==1,'#00441b',ifelse(PValsPaired<0.001&totalPval<0.001,'#08306b','#a50f15')),lwd=3)
    }
    
    axis(side =2,las=1)
    abline(h=seq((minYaxis+1),(maxYaxis-1),by=1),lty='dashed',lwd=2,col='gray')
    
    mtext(side=3
          ,at=fun.chrom.mean(chrom.segs[,1:3])
          ,text=sort(unique(chrom.segs[,1]))
          ,cex=seq(0.6,0.4,length.out=length(unique(chrom.segs[,1])))
          ,line=-1
          #,lwd=0.5
    )
    
    for (i in sort(unique(as.numeric(chrom.segs[,1]))))
    {
      abline(v=max(as.numeric(chrom.segs[chrom.segs[,1]==i,3]))
             ,lwd=0.5)
      #,lty='dashed')
    }
    
    #minor Allele ####
    par(mar=c(5,5,2,5))
    chrom.length.copy   <- fun.chrom.length(min.cons.seg[,1:3])
    names(chrom.length.copy)   <- 1:22
    chrom.segs          <- fun.add.chrom(min.cons.seg[,1:3],chrom.length.copy)
    min.x               <- 0
    min.y               <- minYaxis
    max.y               <- maxYaxis
    max.x               <- as.numeric(max(chrom.segs[,'End']))
    
    plot(1
         ,xlim=c(min.x,max.x)
         ,ylim=c(min.y,max.y)
         ,type='n'
         ,xaxt='n'
         ,yaxs='i'
         ,yaxt='n'
         ,xaxs='i'
         ,ylab=""
         ,xlab=""
         ,lwd=2
         ,main=paste("Minor allele ", regionOne," vs ",regionTwo,sep="")
         ,bty="n")
    box()
    
    
    for (x in 1:length(names_ComparisonOutPut))
    {
      ChrCoords <- unlist(strsplit(names_ComparisonOutPut[x],split=":"))
      toPlot <- cbind(ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionOne]]$chr
                      ,ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionOne]]$pos
                      ,ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionOne]]$pos
                      ,ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionOne]]$nA
                      ,ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionOne]]$nB
                      ,ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionTwo]]$nA
                      ,ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionTwo]]$nB)
      if(is.na(toPlot[1,1]))
      {
        next;
      }
      
      #toPlot <- cbind(logRBAFSample1$chrom,logRBAFSample1$pos,logRBAFSample1$pos,ABSample1,ABSample2)
      toPlot[,1:3] <- fun.add.chrom(apply(toPlot[,1:3],2,as.numeric),chrom.length.copy)
      toPlot <- data.frame(toPlot,stringsAsFactors = FALSE)
      colnames(toPlot) <- c('chr','pos','pos','Sample1nA','Sample1nB','Sample2nA','Sample2nB')
      PValsPaired <- as.numeric(ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$nB_PairedPval)
      totalPval   <-  ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$totalCPNPairedPval#as.numeric(ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$nA_PairedPval)
      MSAIstatus  <- as.numeric(ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$msAI)
      points(x = as.numeric(toPlot$pos),y=as.numeric(toPlot$Sample1nB)-as.numeric(toPlot$Sample2nB),col=ifelse(MSAIstatus==1,'#238b4550',ifelse(PValsPaired<0.001&totalPval<0.001,'#3182bd50','#fc927250')),cex=0.75,pch=16)
      par(lend=1)
      segments(x0=min(as.numeric(toPlot$pos))
               ,y0=mean(c(as.numeric(toPlot$Sample1nB)-as.numeric(toPlot$Sample2nB)))
               ,x1=max(as.numeric(toPlot$pos)),y1=mean(c(as.numeric(toPlot$Sample1nB)-as.numeric(toPlot$Sample2nB)))
               ,col=ifelse(MSAIstatus==1,'#00441b',ifelse(PValsPaired<0.001&totalPval<0.001,'#08306b','#a50f15')),lwd=3)
      
      
      
      #  points(x = as.numeric(toPlot$pos),y=as.numeric(toPlot$Sample1nB)-as.numeric(toPlot$Sample2nB),col=ifelse(PValsPaired<0.001,'#3182bd50','#fb6a4a50'),cex=0.75,pch=16)
      # segments(x0=min(as.numeric(toPlot$pos))
      #          ,y0=mean(c(as.numeric(toPlot$Sample1nB)-as.numeric(toPlot$Sample2nB)))
      #          ,x1=max(as.numeric(toPlot$pos)),y1=mean(c(as.numeric(toPlot$Sample1nB)-as.numeric(toPlot$Sample2nB)))
      #          ,col=ifelse(PValsPaired<0.001,'#08519c','#a50f15'),lwd=3)
    }
    
    
    axis(side =2,las=1)
    
    mtext(side=3
          ,at=fun.chrom.mean(chrom.segs[,1:3])
          ,text=sort(unique(chrom.segs[,1]))
          ,cex=seq(0.6,0.4,length.out=length(unique(chrom.segs[,1])))
          ,line=-1
          #,lwd=0.5
    )
    
    for (i in sort(unique(as.numeric(chrom.segs[,1]))))
    {
      abline(v=max(as.numeric(chrom.segs[chrom.segs[,1]==i,3]))
             ,lwd=0.5)
      #,lty='dashed')
    }
    abline(h=seq((minYaxis+1),(maxYaxis-1),by=1),lty='dashed',lwd=2,col='gray')
    
    
    
    
    
  }
  print("Returning from calcDiffBetweenTwoRegions")
  invisible(ComparisonOutPut)
  
}


calcDiffBetweenTwoRegions_withLogR <- function(regionOne
                                      , regionTwo
                                      , bigList
                                      , bigListLogR
                                      , min.cons.seg
                                      , gammaToUse=gammaToUse
                                      , plotting=TRUE
                                      , minYaxis=(-2)
                                      , maxYaxis=2
                                      , purityEstimates
                                      , ploidyEstimates
                                      #, phasing='singleRegion'
                                      , whichRegionToPhaseFrom=NULL
                                      , gainThreshold    = log2(2.5/2)
                                      , lossThreshold    = log2(1.5/2)
                                      , ampThreshold     = log2(4/2)
                                      , deepLossThreshold= log2(0.5/2)
)
{
  
  print("Called calcDiffBetweenTwoRegions_withLogR")
  # Make sure any copy numbers that are too high are still included
  #layout(rbind(1,2))
  
  if(length(whichRegionToPhaseFrom)==0)
  {
    ComparisonOutPut <- list()
    names_ComparisonOutPut <- apply(min.cons.seg[,1:3],1,PasteVector,sep=":")
    HeterogeneityScore <- c()
    SizesOfSegments    <- c()
    for (x in 1:length(names_ComparisonOutPut))
    {
      
      
      
      
      #cat(x)
      # regionToConsider
      logRBAFSample1 <- bigList[[regionOne]][[x]] #ASCATdata[[regionOne]][ASCATdata[[regionOne]]$chrom==gain.min.cons[x,1]&ASCATdata[[regionOne]]$pos>=gain.min.cons[x,2]&ASCATdata[[regionOne]]$pos<=gain.min.cons[x,3],,drop=FALSE]
      logRBAFSample2 <- bigList[[regionTwo]][[x]] #ASCATdata[[regionTwo]][ASCATdata[[regionTwo]]$chrom==gain.min.cons[x,1]&ASCATdata[[regionTwo]]$pos>=gain.min.cons[x,2]&ASCATdata[[regionTwo]]$pos<=gain.min.cons[x,3],,drop=FALSE]
      #logRBAFSample1 <- logRBAFSample1[!is.na(logRBAFSample1$BAFsegmented),,drop=FALSE]
      #logRBAFSample2 <- logRBAFSample2[!is.na(logRBAFSample2$BAFsegmented),,drop=FALSE]
      
      
      
      
      
      if(nrow(logRBAFSample1)<5)
      {
        
        AndBSample1 <- funCalcN_withBAF(logRSites = logRBAFSample1$logR,bafSites = logRBAFSample1$BAF,tumorPloidy = ploidyEstimates[regionOne],tumorPurity = purityEstimates[regionOne],gamma = gammaToUse)
        AndBSample2 <- funCalcN_withBAF(logRSites = logRBAFSample2$logR,bafSites = logRBAFSample2$BAF,tumorPloidy = ploidyEstimates[regionTwo],tumorPurity = purityEstimates[regionTwo],gamma = gammaToUse)
        
        ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionOne]]$chr <- NA
        ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionOne]]$pos <- NA
        ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionOne]]$nA <- NA
        ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionOne]]$nB <- NA
        ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionOne]]$nA_mean <- NA
        ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionOne]]$nB_mean <- NA
        ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$nA_PairedPval        <- NA
        ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$nB_PairedPval        <- NA
        ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$nB_unPairedPval      <- NA
        ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$nA_unPairedPval      <- NA
        ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$combinedPairedPval   <- NA
        ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$combinedUnPairedPval <- NA
        ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$nA_heterogeneity          <- NA
        ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$nB_heterogeneity          <- NA
        ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$combinedHeterogeneity     <- NA
        
        
        #outPut <- cbind(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)
        #return(outPut)
        next;
      }
      
      # make a decision regarding which samples decides A and B
      # this should really be done on all samples together to ensure consistency
      
      AndBSample1 <- funCalcN_withBAF(logRSites = logRBAFSample1$logR,bafSites = logRBAFSample1$BAF,tumorPloidy = ploidyEstimates[regionOne],tumorPurity = purityEstimates[regionOne],gamma = gammaToUse)
      AndBSample2 <- funCalcN_withBAF(logRSites = logRBAFSample2$logR,bafSites = logRBAFSample2$BAF,tumorPloidy = ploidyEstimates[regionTwo],tumorPurity = purityEstimates[regionTwo],gamma = gammaToUse)
      
      Sample1MedianMax <- median(apply(AndBSample1,1,max))
      Sample2MedianMax <- median(apply(AndBSample2,1,max))
      
      Sample1Max <- apply(AndBSample1,1,max)
      Sample2Max <- apply(AndBSample2,1,max)
      tmpSample1 <- ifelse(AndBSample1==Sample1Max,1,0)
      tmpSample2 <- ifelse(AndBSample2==Sample2Max,1,0)
      
      if(identical(tmpSample1,tmpSample2))
      {
        Sample1Classifier <- 0.05
        Sample2Classifier <- 1
      }
      # check if one is significant using the other classification. If so, use that one
      
      if(!identical(tmpSample1,tmpSample2))
      {
        Sample1Classifier <- tTestWithErrorHandling(c(AndBSample2[which(tmpSample1[,1]==1),1],AndBSample2[which(tmpSample1[,2]==1),2]),c(AndBSample2[which(tmpSample1[,1]!=1),2],AndBSample2[which(tmpSample1[,2]!=1),2]))$p.value
        Sample2Classifier <- tTestWithErrorHandling(c(AndBSample1[which(tmpSample2[,1]==1),1],AndBSample1[which(tmpSample2[,2]==1),2]),c(AndBSample1[which(tmpSample2[,1]!=1),2],AndBSample1[which(tmpSample2[,2]!=1),2]))$p.value
        
        
        
        if(is.na(Sample1Classifier))
        {
          Sample1Classifier <- 1
        }
        if(is.na(Sample2Classifier))
        {
          Sample2Classifier <- 1
        }
      }
      
      
      # Major <- apply(cbind(AndBSample1,AndBSample2),1,max)
      # tmpSample1 <- ifelse(AndBSample1==Major,1,0)
      # tmpSample2 <- ifelse(AndBSample2==Major,1,0)
      
      ABSample1 <- AndBSample1
      ABSample2 <- AndBSample2
      
      if(Sample1Classifier<=Sample2Classifier)
      {
        MajorColumn1 <- which(tmpSample1[,1]==1)
        MajorColumn2 <- which(tmpSample1[,2]==1)
      }
      if(Sample1Classifier>Sample2Classifier)
      {
        MajorColumn1 <- which(tmpSample2[,1]==1)
        MajorColumn2 <- which(tmpSample2[,2]==1)
      }
      #Sample1Decides <-  which(cbind(rowSums(tmpSample1),rowSums(tmpSample2))[,1]==1)
      #Sample2Decides <-  which(cbind(rowSums(tmpSample1),rowSums(tmpSample2))[,2]==1)
      
      #MajorColumn1   <- c(Sample1Decides[which(tmpSample1[Sample1Decides,1]==1)],Sample2Decides[which(tmpSample2[Sample2Decides,1]==1)])
      #MajorColumn2   <- c(Sample1Decides[which(tmpSample1[Sample1Decides,2]==1)],Sample2Decides[which(tmpSample2[Sample2Decides,2]==1)])
      
      ABSample1[MajorColumn1,1] <- AndBSample1[MajorColumn1,1]
      ABSample1[MajorColumn2,1] <- AndBSample1[MajorColumn2,2]
      ABSample1[MajorColumn1,2] <- AndBSample1[MajorColumn1,2]
      ABSample1[MajorColumn2,2] <- AndBSample1[MajorColumn2,1]
      
      ABSample2[MajorColumn1,1] <- AndBSample2[MajorColumn1,1]
      ABSample2[MajorColumn2,1] <- AndBSample2[MajorColumn2,2]
      ABSample2[MajorColumn1,2] <- AndBSample2[MajorColumn1,2]
      ABSample2[MajorColumn2,2] <- AndBSample2[MajorColumn2,1]
      
      
      
      PValsPaired          <- c(tTestWithErrorHandling(ABSample1[,1],ABSample2[,1],paired=TRUE)$p.val,tTestWithErrorHandling(ABSample1[,2],ABSample2[,2],paired=TRUE)$p.val)
      combinedPvalPaired   <- pchisq((sum(log(PValsPaired))*-2), df=length(PValsPaired)*2, lower.tail=F)
      PValsUnPaired        <- c(tTestWithErrorHandling(ABSample1[,1],ABSample2[,1])$p.val,tTestWithErrorHandling(ABSample1[,2],ABSample2[,2])$p.val)
      combinedPvalUnPaired <- pchisq((sum(log(PValsUnPaired))*-2), df=length(PValsUnPaired)*2, lower.tail=F)
      
      # also get a p-value that relates to both samples combined
      totalCPNpval         <- tTestWithErrorHandling(c(ABSample1[,1]+ABSample1[,2]),c(ABSample2[,1]+ABSample2[,2]),paired=TRUE)$p.val
      
      nA_heterogeneity       <- 0
      nB_heterogeneity       <- 0
      SizeOfSegment           <- strsplit(names_ComparisonOutPut[x],split=":")
      SizeOfSegment           <- as.numeric(unlist(SizeOfSegment)[3])-as.numeric(unlist(SizeOfSegment)[2])
      
      if(PValsPaired[1]<0.001)
      {
        nAamplitudeOfDifference <- abs(mean(as.numeric(ABSample1[,1]))-mean(as.numeric(ABSample2[,1])))
        SizeOfSegment           <- strsplit(names_ComparisonOutPut[x],split=":")
        SizeOfSegment           <- as.numeric(unlist(SizeOfSegment)[3])-as.numeric(unlist(SizeOfSegment)[2])
        nA_heterogeneity        <- nAamplitudeOfDifference*SizeOfSegment
      }
      
      if(PValsPaired[2]<0.001)
      {
        nBamplitudeOfDifference <- abs(mean(as.numeric(ABSample1[,2]))-mean(as.numeric(ABSample2[,2])))
        SizeOfSegment           <- strsplit(names_ComparisonOutPut[x],split=":")
        SizeOfSegment           <- as.numeric(unlist(SizeOfSegment)[3])-as.numeric(unlist(SizeOfSegment)[2])
        nB_heterogeneity        <- nBamplitudeOfDifference*SizeOfSegment
      }
      
      total_heterogeneity <- nA_heterogeneity+nB_heterogeneity
      HeterogeneityScore  <- c(HeterogeneityScore,total_heterogeneity)
      SizesOfSegments     <- c(SizesOfSegments,SizeOfSegment)
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionOne]]$chr     <- logRBAFSample1$chrom
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionOne]]$pos     <- logRBAFSample1$pos
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionOne]]$nA      <- as.numeric(ABSample1[,1])
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionOne]]$nB      <- as.numeric(ABSample1[,2])
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionOne]]$nA_mean <- mean(as.numeric(ABSample1[,1]))
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionOne]]$nB_mean <- mean(as.numeric(ABSample1[,2]))
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionTwo]]$nA      <- as.numeric(ABSample2[,1])
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionTwo]]$nB      <- as.numeric(ABSample2[,2])
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionTwo]]$nA_mean <- mean(as.numeric(ABSample2[,1]))
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionTwo]]$nB_mean <- mean(as.numeric(ABSample2[,2]))
      
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$nA_PairedPval             <- PValsPaired[1]
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$nB_PairedPval             <- PValsPaired[2]
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$nB_unPairedPval           <- PValsUnPaired[2]
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$nA_unPairedPval           <- PValsUnPaired[1]
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$combinedPairedPval        <- combinedPvalPaired
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$combinedUnPairedPval      <- combinedPvalUnPaired
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$combinedUnPairedPval      <- combinedPvalUnPaired
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$totalCPNPairedPval        <- totalCPNpval
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$nA_heterogeneity          <- nA_heterogeneity
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$nB_heterogeneity          <- nB_heterogeneity
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$combinedHeterogeneity     <- nB_heterogeneity+nA_heterogeneity
      # classify the segment as showing MSAI
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$msAI                      <- 0
      # if(PValsPaired[1]<0.001&PValsPaired[2]<0.001)
      # {
      
      if(tTestWithErrorHandling(ABSample1[,1],ABSample1[,2],paired=TRUE)$p.val<0.001&tTestWithErrorHandling(ABSample2[,1],ABSample2[,2],paired=TRUE)$p.val<0.001)
      {
        if((mean(ABSample1[,1]-ABSample1[,2])>0&mean(ABSample2[,1]-ABSample2[,2])<0)|(mean(ABSample1[,1]-ABSample1[,2])<0&mean(ABSample2[,1]-ABSample2[,2])>0))
        {
          ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$msAI                      <- 1
        }
        
      }
      
      
      checkGainSignificance <- function(logRBAFRegion,regionToCheck,ploidyEstimates,gainThreshold)
      {
        minGainRawExpRegion       <-  ploidyEstimates[regionToCheck]*(2^(gainThreshold))
        minGainLogRExpRegion      <-  funCalcExpectedLogR(minGainRawExpRegion,0,ploidyEstimates[regionToCheck],purityEstimates[regionToCheck])
        logRObs                   <-  mean(logRBAFRegion$logR)
        
        if(nrow(logRBAFRegion)>=3)
        {
          Ttest        <- tTestWithErrorHandling(logRBAFRegion$logR,mu = minGainLogRExpRegion)
          TtestGreater <- tTestWithErrorHandling(logRBAFRegion$logR,mu = minGainLogRExpRegion,alternative = 'greater')
          #TtestLess    <- tTestWithErrorHandling(logRBAF$logR,mu = minLossLogRExp,alternative = 'less')
          logRObs_low  <- Ttest$conf.int[1]
          logRObs_high <- Ttest$conf.int[2]
          PVal         <- Ttest$p.value
          PvalGain     <- TtestGreater$p.value
          
        }
        if(nrow(logRBAFRegion)<3)
        {
          logRObs_low  <- NA
          logRObs_high <- NA
          PVal         <- NA
          PvalGain     <- NA
          
          
        }
        
        out <- cbind(minGainRawExpRegion,minGainLogRExpRegion,logRObs,logRObs_low,logRObs_high,PVal,PvalGain)  
        return(out)
        
        
      }
      checkLossSignificance <- function(logRBAFRegion,regionToCheck,ploidyEstimates,lossThreshold)
      {
        minLossRawExpRegion       <-  ploidyEstimates[regionToCheck]*(2^(lossThreshold))
        minLossLogRExpRegion      <-  funCalcExpectedLogR(minLossRawExpRegion,0,ploidyEstimates[regionToCheck],purityEstimates[regionToCheck])
        logRObs                   <- mean(logRBAFRegion$logR)
        
        if(nrow(logRBAFRegion)>=3)
        {
          Ttest        <- tTestWithErrorHandling(logRBAFRegion$logR,mu = minLossLogRExpRegion)
          TtestLess <- tTestWithErrorHandling(logRBAFRegion$logR,mu = minLossLogRExpRegion,alternative = 'less')
          #TtestLess    <- tTestWithErrorHandling(logRBAF$logR,mu = minLossLogRExp,alternative = 'less')
          logRObs_low  <- Ttest$conf.int[1]
          logRObs_high <- Ttest$conf.int[2]
          PVal         <- Ttest$p.value
          PvalLoss     <- TtestLess$p.value
          
        }
        # if(nrow(logRBAF)<3)
        # {
        #   logRObs_low  <- NA
        #   logRObs_high <- NA
        #   PVal         <- NA
        #   PvalLoss     <- NA
        
        
        # }
        
        out <- cbind(minLossRawExpRegion,minLossLogRExpRegion,logRObs,logRObs_low,logRObs_high,PVal,PvalLoss)  
        return(out)
        
        
      }
      
      
      GainRegionOne     <- checkGainSignificance(logRBAFSample1,regionToCheck = regionOne,ploidyEstimates,gainThreshold = gainThreshold)
      AmpRegionOne      <- checkGainSignificance(logRBAFSample1,regionToCheck = regionOne,ploidyEstimates,gainThreshold = ampThreshold)
      LossRegionOne     <- checkLossSignificance(logRBAFSample1,regionToCheck = regionOne,ploidyEstimates,lossThreshold = lossThreshold)
      DeepLossRegionOne <- checkLossSignificance(logRBAFSample1,regionToCheck = regionOne,ploidyEstimates,lossThreshold = deepLossThreshold)
      
      
      # put together into One table
      regionOneGainOrLoss <- c(GainRegionOne[,3:5,drop=FALSE]
                               ,GainRegionOne[,2,drop=FALSE],GainRegionOne[,7,drop=FALSE]
                               ,AmpRegionOne[,2,drop=FALSE],AmpRegionOne[,7,drop=FALSE]
                               ,LossRegionOne[,2,drop=FALSE],LossRegionOne[,7,drop=FALSE]
                               ,DeepLossRegionOne[,2,drop=FALSE],DeepLossRegionOne[,7,drop=FALSE]
      )
      names(regionOneGainOrLoss) <- c('logRObs','logRObs_low','logRObs_high'
                                      ,'minGainExpLogR','gainPval'
                                      ,'minAmpExpLogR','ampPval'
                                      ,'minLossExpLogR','lossPval'
                                      ,'minDeepLossExpLogR','DeepLossPval')
      
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionOne]]$segStatus     <- regionOneGainOrLoss
      
      
      GainRegionTwo     <- checkGainSignificance(logRBAFSample2,regionToCheck = regionTwo,ploidyEstimates,gainThreshold = gainThreshold)
      AmpRegionTwo      <- checkGainSignificance(logRBAFSample2,regionToCheck = regionTwo,ploidyEstimates,gainThreshold = ampThreshold)
      LossRegionTwo     <- checkLossSignificance(logRBAFSample2,regionToCheck = regionTwo,ploidyEstimates,lossThreshold = lossThreshold)
      DeepLossRegionTwo <- checkLossSignificance(logRBAFSample2,regionToCheck = regionTwo,ploidyEstimates,lossThreshold = deepLossThreshold)
      
      
      # put together into Two table
      regionTwoGainOrLoss <- c(GainRegionTwo[,3:5,drop=FALSE]
                               ,GainRegionTwo[,2,drop=FALSE],GainRegionTwo[,7,drop=FALSE]
                               ,AmpRegionTwo[,2,drop=FALSE],AmpRegionTwo[,7,drop=FALSE]
                               ,LossRegionTwo[,2,drop=FALSE],LossRegionTwo[,7,drop=FALSE]
                               ,DeepLossRegionTwo[,2,drop=FALSE],DeepLossRegionTwo[,7,drop=FALSE]
      )
      names(regionTwoGainOrLoss) <- c('logRObs','logRObs_low','logRObs_high'
                                      ,'minGainExpLogR','gainPval'
                                      ,'minAmpExpLogR','ampPval'
                                      ,'minLossExpLogR','lossPval'
                                      ,'minDeepLossExpLogR','DeepLossPval')
      
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionTwo]]$segStatus     <- regionTwoGainOrLoss
      
      
      
      # # for it to be classified as MSAI the difference between A1-A2 must be positive and B1-B2 negative (or the other way round)
      # if((mean(as.numeric(ABSample1[,1]))-mean(as.numeric(ABSample2[,1]))>0&mean(as.numeric(ABSample1[,2]))-mean(as.numeric(ABSample2[,2]))<0)|(mean(as.numeric(ABSample1[,1]))-mean(as.numeric(ABSample2[,1]))<0&mean(as.numeric(ABSample1[,2]))-mean(as.numeric(ABSample2[,2]))>0))
      # {
      #   # we also must observe significant difference in A1 and B1 as well as A2 and B2. And this significant difference must be opposite. 
      #   
      # }
      # }
      
      
      
      
    }
    
    averageHeterogeneity <- sum(HeterogeneityScore)/sum(SizesOfSegments)
    
    #ComparisonOutPut$`1:762273:11831543`$`L_LTX064_SU_T1-R6:L_LTX064_SU_T1-R7`$nA_heterogeneity
    #ComparisonOutPut$`1:762273:11831543`$`L_LTX064_SU_T1-R6:L_LTX064_SU_T1-R7`$nB_heterogeneity
    ComparisonOutPut[[paste(regionOne,regionTwo,sep=":")]]$HeterogeneityScore <- averageHeterogeneity 
    
  }
  
  if(length(whichRegionToPhaseFrom)!=0)
  {
    ComparisonOutPut <- list()
    names_ComparisonOutPut <- apply(min.cons.seg[,1:3],1,PasteVector,sep=":")
    HeterogeneityScore <- c()
    SizesOfSegments    <- c()
    for (x in 1:length(names_ComparisonOutPut))
    {
      
      #cat(x)
      # regionToConsider
      NameToUse <- names_ComparisonOutPut[x]
      
      
      logRBAFSample1 <- bigList[[regionOne]][[x]] #ASCATdata[[regionOne]][ASCATdata[[regionOne]]$chrom==gain.min.cons[x,1]&ASCATdata[[regionOne]]$pos>=gain.min.cons[x,2]&ASCATdata[[regionOne]]$pos<=gain.min.cons[x,3],,drop=FALSE]
      logRBAFSample2 <- bigList[[regionTwo]][[x]] #ASCATdata[[regionTwo]][ASCATdata[[regionTwo]]$chrom==gain.min.cons[x,1]&ASCATdata[[regionTwo]]$pos>=gain.min.cons[x,2]&ASCATdata[[regionTwo]]$pos<=gain.min.cons[x,3],,drop=FALSE]
      logRSample1 <- bigListLogR[[regionOne]][[x]] #ASCATdata[[regionOne]][ASCATdata[[regionOne]]$chrom==gain.min.cons[x,1]&ASCATdata[[regionOne]]$pos>=gain.min.cons[x,2]&ASCATdata[[regionOne]]$pos<=gain.min.cons[x,3],,drop=FALSE]
      logRSample2 <- bigListLogR[[regionTwo]][[x]] #ASCATdata[[regionTwo]][ASCATdata[[regionTwo]]$chrom==gain.min.cons[x,1]&ASCATdata[[regionTwo]]$pos>=gain.min.cons[x,2]&ASCATdata[[regionTwo]]$pos<=gain.min.cons[x,3],,drop=FALSE]
      
      #logRSample3 <- bigListLogR[[regionThree]][[x]] #ASCATdata[[regionOne]][ASCATdata[[regionOne]]$chrom==gain.min.cons[x,1]&ASCATdata[[regionOne]]$pos>=gain.min.cons[x,2]&ASCATdata[[regionOne]]$pos<=gain.min.cons[x,3],,drop=FALSE]
      #logRSample4 <- bigListLogR[[regionFourth]][[x]] #ASCATdata[[regionTwo]][ASCATdata[[regionTwo]]$chrom==gain.min.cons[x,1]&ASCATdata[[regionTwo]]$pos>=gain.min.cons[x,2]&ASCATdata[[regionTwo]]$pos<=gain.min.cons[x,3],,drop=FALSE]
      
      # plot(logRSample1$logR,col=ifelse(is.na(logRSample1$BAFsegmented),'#bdbdbd99','#de2d26'))
      # abline(h=unique(logRSample1$logRsegmented))
      # plot(logRSample2$logR,col=ifelse(is.na(logRSample2$BAFsegmented),'#bdbdbd99','#de2d26'))
      # abline(h=unique(logRSample2$logRsegmented))
      # 
      # plot(logRSample3$logR,col=ifelse(is.na(logRSample3$BAFsegmented),'#bdbdbd99','#de2d26'))
      # abline(h=unique(logRSample3$logRsegmented))
      # plot(logRSample4$logR,col=ifelse(is.na(logRSample3$BAFsegmented),'#bdbdbd99','#de2d26'))
      # abline(h=unique(logRSample4$logRsegmented))
      # 
      # test <-logRSample1[!is.na(logRSample1$BAF),,drop=FALSE]
      
      #logRBAFSample1 <- logRBAFSample1[!is.na(logRBAFSample1$BAFsegmented),,drop=FALSE]
      #logRBAFSample2 <- logRBAFSample2[!is.na(logRBAFSample2$BAFsegmented),,drop=FALSE]
      
      if(nrow(logRBAFSample1)<5)
      {
        NameSplit <- unlist(strsplit(NameToUse,split=":"))
        
        ABSample1 <- funCalcN_withBAF(logRSites = logRSample1$logR,bafSites = rep(0.5,nrow(logRSample1)),tumorPloidy = ploidyEstimates[regionOne],tumorPurity = purityEstimates[regionOne],gamma = gammaToUse)
        ABSample2 <- funCalcN_withBAF(logRSites = logRSample2$logR,bafSites = rep(0.5,nrow(logRSample2)),tumorPloidy = ploidyEstimates[regionTwo],tumorPurity = purityEstimates[regionTwo],gamma = gammaToUse)
        
        PValsPaired          <- c(tTestWithErrorHandling(ABSample1[,1],ABSample2[,1],paired=TRUE)$p.val,tTestWithErrorHandling(ABSample1[,2],ABSample2[,2],paired=TRUE)$p.val)
        combinedPvalPaired   <- PValsPaired[1]
        PValsUnPaired        <- c(tTestWithErrorHandling(ABSample1[,1],ABSample2[,1])$p.val,tTestWithErrorHandling(ABSample1[,2],ABSample2[,2])$p.val)
        combinedPvalUnPaired <- PValsUnPaired[1]#pchisq((sum(log(PValsUnPaired))*-2), df=length(PValsUnPaired)*2, lower.tail=F)
        
        # also get a p-value that relates to both samples combined
        totalCPNpval         <- tTestWithErrorHandling(c(ABSample1[,1]+ABSample1[,2]),c(ABSample2[,1]+ABSample2[,2]),paired=TRUE)$p.val
        
        nA_heterogeneity       <- 0
        nB_heterogeneity       <- 0
        SizeOfSegment           <- strsplit(names_ComparisonOutPut[x],split=":")
        SizeOfSegment           <- as.numeric(unlist(SizeOfSegment)[3])-as.numeric(unlist(SizeOfSegment)[2])
        
        if(PValsPaired[1]<0.001)
        {
          nAamplitudeOfDifference <- abs(mean(as.numeric(ABSample1[,1]))-mean(as.numeric(ABSample2[,1])))
          SizeOfSegment           <- strsplit(names_ComparisonOutPut[x],split=":")
          SizeOfSegment           <- as.numeric(unlist(SizeOfSegment)[3])-as.numeric(unlist(SizeOfSegment)[2])
          nA_heterogeneity        <- nAamplitudeOfDifference*SizeOfSegment
        }
        
        if(PValsPaired[2]<0.001)
        {
          nBamplitudeOfDifference <- abs(mean(as.numeric(ABSample1[,2]))-mean(as.numeric(ABSample2[,2])))
          SizeOfSegment           <- strsplit(names_ComparisonOutPut[x],split=":")
          SizeOfSegment           <- as.numeric(unlist(SizeOfSegment)[3])-as.numeric(unlist(SizeOfSegment)[2])
          nB_heterogeneity        <- nBamplitudeOfDifference*SizeOfSegment
        }
        
        
        total_heterogeneity <- nA_heterogeneity+nB_heterogeneity
        HeterogeneityScore  <- c(HeterogeneityScore,total_heterogeneity)
        SizesOfSegments     <- c(SizesOfSegments,SizeOfSegment)
        ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionOne]]$chr     <- logRSample1$chrom
        ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionOne]]$pos     <- logRSample1$pos
        ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionOne]]$nA      <- NA #as.numeric(ABSample1[,1])
        ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionOne]]$nB      <- NA #as.numeric(ABSample1[,2])
        ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionOne]]$nA_mean <- NA #mean(as.numeric(ABSample1[,1]))
        ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionOne]]$nB_mean <- NA 
        ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionOne]]$cnTotal <- ABSample1[,1]+ABSample1[,2]
        ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionOne]]$cnTotal_mean <- mean(ABSample1[,1]+ABSample1[,2])##mean(as.numeric(ABSample1[,2]))
        ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionTwo]]$nA      <- NA #as.numeric(ABSample2[,1])
        ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionTwo]]$nB      <- NA #as.numeric(ABSample2[,2])
        ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionTwo]]$nA_mean <- NA #mean(as.numeric(ABSample2[,1]))
        ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionTwo]]$nB_mean <- NA #mean(as.numeric(ABSample2[,2]))
        ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionTwo]]$cnTotal <- ABSample2[,1]+ABSample2[,2]
        ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionTwo]]$cnTotal_mean <- mean(ABSample2[,1]+ABSample2[,2])#mean(as.numeric(ABSample2[,2]))
        
        ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$nA_PairedPval             <- NA #PValsPaired[1]
        ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$nB_PairedPval             <- NA #PValsPaired[2]
        ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$nB_unPairedPval           <- NA #PValsUnPaired[2]
        ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$nA_unPairedPval           <- NA #PValsUnPaired[1]
        ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$combinedPairedPval        <- combinedPvalPaired
        ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$combinedUnPairedPval      <- combinedPvalUnPaired
        ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$combinedUnPairedPval      <- combinedPvalUnPaired
        ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$totalCPNPairedPval        <- totalCPNpval
        ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$nA_heterogeneity          <- NA #nA_heterogeneity
        ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$nB_heterogeneity          <- NA #nB_heterogeneity
        ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$combinedHeterogeneity     <- NA #nB_heterogeneity+nA_heterogeneity
        # classify the segment as showing MSAI
        ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$msAI                      <- NA
        # if(PValsPaired[1]<0.001&PValsPaired[2]<0.001)
        # {
        
        checkGainSignificance <- function(logRBAFRegion,regionToCheck,ploidyEstimates,gainThreshold)
        {
          minGainRawExpRegion       <-  ploidyEstimates[regionToCheck]*(2^(gainThreshold))
          minGainLogRExpRegion      <-  funCalcExpectedLogR(minGainRawExpRegion,0,ploidyEstimates[regionToCheck],purityEstimates[regionToCheck])
          logRObs                   <- median(logRBAFRegion$logR)
          
          if(nrow(logRBAFRegion)>=2)
          {
            Ttest        <- tTestWithErrorHandling(logRBAFRegion$logR,mu = minGainLogRExpRegion)
            TtestGreater <- tTestWithErrorHandling(logRBAFRegion$logR,mu = minGainLogRExpRegion,alternative = 'greater')
            #TtestLess    <- tTestWithErrorHandling(logRBAF$logR,mu = minLossLogRExp,alternative = 'less')
            logRObs_low  <- Ttest$conf.int[1]
            logRObs_high <- Ttest$conf.int[2]
            PVal         <- Ttest$p.value
            PvalGain     <- TtestGreater$p.value
            
          }
          if(nrow(logRBAFRegion)<2)
          {
            logRObs_low  <- NA
            logRObs_high <- NA
            PVal         <- NA
            PvalGain     <- NA
            
            
          }
          
          out <- cbind(minGainRawExpRegion,minGainLogRExpRegion,logRObs,logRObs_low,logRObs_high,PVal,PvalGain)  
          return(out)
          
          
        }
        checkLossSignificance <- function(logRBAFRegion,regionToCheck,ploidyEstimates,lossThreshold)
        {
          minLossRawExpRegion       <-  ploidyEstimates[regionToCheck]*(2^(lossThreshold))
          minLossLogRExpRegion      <-  funCalcExpectedLogR(minLossRawExpRegion,0,ploidyEstimates[regionToCheck],purityEstimates[regionToCheck])
          logRObs                   <- median(logRBAFRegion$logR)
          
          if(nrow(logRBAFRegion)>=2)
          {
            Ttest        <- tTestWithErrorHandling(logRBAFRegion$logR,mu = minLossLogRExpRegion)
            TtestLess <- tTestWithErrorHandling(logRBAFRegion$logR,mu = minLossLogRExpRegion,alternative = 'less')
            #TtestLess    <- tTestWithErrorHandling(logRBAF$logR,mu = minLossLogRExp,alternative = 'less')
            logRObs_low  <- Ttest$conf.int[1]
            logRObs_high <- Ttest$conf.int[2]
            PVal         <- Ttest$p.value
            PvalLoss     <- TtestLess$p.value
            
          }
          if(nrow(logRBAFRegion)<2)
          {
            logRObs_low  <- NA
            logRObs_high <- NA
            PVal         <- NA
            PvalLoss     <- NA
            
            
          }
          
          out <- cbind(minLossRawExpRegion,minLossLogRExpRegion,logRObs,logRObs_low,logRObs_high,PVal,PvalLoss)  
          return(out)
          
          
        }
        
        
        GainRegionOne     <- checkGainSignificance(logRSample1,regionToCheck = regionOne,ploidyEstimates,gainThreshold = gainThreshold)
        AmpRegionOne      <- checkGainSignificance(logRSample1,regionToCheck = regionOne,ploidyEstimates,gainThreshold = ampThreshold)
        LossRegionOne     <- checkLossSignificance(logRSample1,regionToCheck = regionOne,ploidyEstimates,lossThreshold = lossThreshold)
        DeepLossRegionOne <- checkLossSignificance(logRSample1,regionToCheck = regionOne,ploidyEstimates,lossThreshold = deepLossThreshold)
        
        # put together into One table
        regionOneGainOrLoss <- c(GainRegionOne[,3:5,drop=FALSE]
                                 ,GainRegionOne[,2,drop=FALSE],GainRegionOne[,7,drop=FALSE]
                                 ,AmpRegionOne[,2,drop=FALSE],AmpRegionOne[,7,drop=FALSE]
                                 ,LossRegionOne[,2,drop=FALSE],LossRegionOne[,7,drop=FALSE]
                                 ,DeepLossRegionOne[,2,drop=FALSE],DeepLossRegionOne[,7,drop=FALSE]
        )
        names(regionOneGainOrLoss) <- c('logRObs','logRObs_low','logRObs_high'
                                        ,'minGainExpLogR','gainPval'
                                        ,'minAmpExpLogR','ampPval'
                                        ,'minLossExpLogR','lossPval'
                                        ,'minDeepLossExpLogR','DeepLossPval')
        
        ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionOne]]$segStatus     <- regionOneGainOrLoss
        
        
        GainRegionTwo     <- checkGainSignificance(logRSample2,regionToCheck = regionTwo,ploidyEstimates,gainThreshold = gainThreshold)
        AmpRegionTwo      <- checkGainSignificance(logRSample2,regionToCheck = regionTwo,ploidyEstimates,gainThreshold = ampThreshold)
        LossRegionTwo     <- checkLossSignificance(logRSample2,regionToCheck = regionTwo,ploidyEstimates,lossThreshold = lossThreshold)
        DeepLossRegionTwo <- checkLossSignificance(logRSample2,regionToCheck = regionTwo,ploidyEstimates,lossThreshold = deepLossThreshold)
        
        
        # put together into Two table
        regionTwoGainOrLoss <- c(GainRegionTwo[,3:5,drop=FALSE]
                                 ,GainRegionTwo[,2,drop=FALSE],GainRegionTwo[,7,drop=FALSE]
                                 ,AmpRegionTwo[,2,drop=FALSE],AmpRegionTwo[,7,drop=FALSE]
                                 ,LossRegionTwo[,2,drop=FALSE],LossRegionTwo[,7,drop=FALSE]
                                 ,DeepLossRegionTwo[,2,drop=FALSE],DeepLossRegionTwo[,7,drop=FALSE]
        )
        names(regionTwoGainOrLoss) <- c('logRObs','logRObs_low','logRObs_high'
                                        ,'minGainExpLogR','gainPval'
                                        ,'minAmpExpLogR','ampPval'
                                        ,'minLossExpLogR','lossPval'
                                        ,'minDeepLossExpLogR','DeepLossPval')
        
        ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionTwo]]$segStatus     <- regionTwoGainOrLoss
        
        
        
        #outPut <- cbind(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)
        #return(outPut)
        next;
      }
      
      # make a decision regarding which samples decides A and B
      # this should really be done on all samples together to ensure consistency
      
      AndBSample1 <- funCalcN_withBAF(logRSites = logRBAFSample1$logR,bafSites = logRBAFSample1$BAF,tumorPloidy = ploidyEstimates[regionOne],tumorPurity = purityEstimates[regionOne],gamma = gammaToUse)
      AndBSample2 <- funCalcN_withBAF(logRSites = logRBAFSample2$logR,bafSites = logRBAFSample2$BAF,tumorPloidy = ploidyEstimates[regionTwo],tumorPurity = purityEstimates[regionTwo],gamma = gammaToUse)
      
      Sample1MedianMax <- median(apply(AndBSample1,1,max))
      Sample2MedianMax <- median(apply(AndBSample2,1,max))
      
      Sample1Max <- apply(AndBSample1,1,max)
      Sample2Max <- apply(AndBSample2,1,max)
      tmpSample1 <- ifelse(AndBSample1==Sample1Max,1,0)
      tmpSample2 <- ifelse(AndBSample2==Sample2Max,1,0)
      
      if(identical(tmpSample1,tmpSample2))
      {
        Sample1Classifier <- 0.05
        Sample2Classifier <- 1
      }
      # check if one is significant using the other classification. If so, use that one
      
      if(!identical(tmpSample1,tmpSample2))
      {
        Sample1Classifier <- tTestWithErrorHandling(c(AndBSample2[which(tmpSample1[,1]==1),1],AndBSample2[which(tmpSample1[,2]==1),2]),c(AndBSample2[which(tmpSample1[,1]!=1),2],AndBSample2[which(tmpSample1[,2]!=1),2]))$p.value
        Sample2Classifier <- tTestWithErrorHandling(c(AndBSample1[which(tmpSample2[,1]==1),1],AndBSample1[which(tmpSample2[,2]==1),2]),c(AndBSample1[which(tmpSample2[,1]!=1),2],AndBSample1[which(tmpSample2[,2]!=1),2]))$p.value
        
        
        
        if(is.na(Sample1Classifier))
        {
          Sample1Classifier <- 1
        }
        if(is.na(Sample2Classifier)) 
        {
          Sample2Classifier <- 1
        }
      }
      
      
      # Major <- apply(cbind(AndBSample1,AndBSample2),1,max)
      # tmpSample1 <- ifelse(AndBSample1==Major,1,0)
      # tmpSample2 <- ifelse(AndBSample2==Major,1,0)
      
      ABSample1 <- AndBSample1
      ABSample2 <- AndBSample2
      
      if(Sample1Classifier<=Sample2Classifier)
      {
        MajorColumn1 <- which(tmpSample1[,1]==1)
        MajorColumn2 <- which(tmpSample1[,2]==1)
      }
      if(Sample1Classifier>Sample2Classifier)
      {
        MajorColumn1 <- which(tmpSample2[,1]==1)
        MajorColumn2 <- which(tmpSample2[,2]==1)
      }
      
      if(length(whichRegionToPhaseFrom[[NameToUse]]$MajorColumn1)>0)
      {
        MajorColumn1 <- whichRegionToPhaseFrom[[NameToUse]]$MajorColumn1
        MajorColumn2 <- whichRegionToPhaseFrom[[NameToUse]]$MajorColumn2
        
      }
      
      #Sample1Decides <-  which(cbind(rowSums(tmpSample1),rowSums(tmpSample2))[,1]==1)
      #Sample2Decides <-  which(cbind(rowSums(tmpSample1),rowSums(tmpSample2))[,2]==1)
      
      #MajorColumn1   <- c(Sample1Decides[which(tmpSample1[Sample1Decides,1]==1)],Sample2Decides[which(tmpSample2[Sample2Decides,1]==1)])
      #MajorColumn2   <- c(Sample1Decides[which(tmpSample1[Sample1Decides,2]==1)],Sample2Decides[which(tmpSample2[Sample2Decides,2]==1)])
      
      ABSample1[MajorColumn1,1] <- AndBSample1[MajorColumn1,1]
      ABSample1[MajorColumn2,1] <- AndBSample1[MajorColumn2,2]
      ABSample1[MajorColumn1,2] <- AndBSample1[MajorColumn1,2]
      ABSample1[MajorColumn2,2] <- AndBSample1[MajorColumn2,1]
      
      ABSample2[MajorColumn1,1] <- AndBSample2[MajorColumn1,1]
      ABSample2[MajorColumn2,1] <- AndBSample2[MajorColumn2,2]
      ABSample2[MajorColumn1,2] <- AndBSample2[MajorColumn1,2]
      ABSample2[MajorColumn2,2] <- AndBSample2[MajorColumn2,1]
      
      
      
      PValsPaired          <- c(tTestWithErrorHandling(ABSample1[,1],ABSample2[,1],paired=TRUE)$p.val,tTestWithErrorHandling(ABSample1[,2],ABSample2[,2],paired=TRUE)$p.val)
      combinedPvalPaired   <- pchisq((sum(log(PValsPaired))*-2), df=length(PValsPaired)*2, lower.tail=F)
      PValsUnPaired        <- c(tTestWithErrorHandling(ABSample1[,1],ABSample2[,1])$p.val,tTestWithErrorHandling(ABSample1[,2],ABSample2[,2])$p.val)
      combinedPvalUnPaired <- pchisq((sum(log(PValsUnPaired))*-2), df=length(PValsUnPaired)*2, lower.tail=F)
      
      # also get a p-value that relates to both samples combined
      totalCPNpval         <- tTestWithErrorHandling(c(ABSample1[,1]+ABSample1[,2]),c(ABSample2[,1]+ABSample2[,2]),paired=TRUE)$p.val
      
      nA_heterogeneity       <- 0
      nB_heterogeneity       <- 0
      SizeOfSegment           <- strsplit(names_ComparisonOutPut[x],split=":")
      SizeOfSegment           <- as.numeric(unlist(SizeOfSegment)[3])-as.numeric(unlist(SizeOfSegment)[2])
      
      if(PValsPaired[1]<0.001)
      {
        nAamplitudeOfDifference <- abs(mean(as.numeric(ABSample1[,1]))-mean(as.numeric(ABSample2[,1])))
        SizeOfSegment           <- strsplit(names_ComparisonOutPut[x],split=":")
        SizeOfSegment           <- as.numeric(unlist(SizeOfSegment)[3])-as.numeric(unlist(SizeOfSegment)[2])
        nA_heterogeneity        <- nAamplitudeOfDifference*SizeOfSegment
      }
      
      if(PValsPaired[2]<0.001)
      {
        nBamplitudeOfDifference <- abs(mean(as.numeric(ABSample1[,2]))-mean(as.numeric(ABSample2[,2])))
        SizeOfSegment           <- strsplit(names_ComparisonOutPut[x],split=":")
        SizeOfSegment           <- as.numeric(unlist(SizeOfSegment)[3])-as.numeric(unlist(SizeOfSegment)[2])
        nB_heterogeneity        <- nBamplitudeOfDifference*SizeOfSegment
      }
      
      total_heterogeneity <- nA_heterogeneity+nB_heterogeneity
      HeterogeneityScore  <- c(HeterogeneityScore,total_heterogeneity)
      SizesOfSegments     <- c(SizesOfSegments,SizeOfSegment)
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionOne]]$chr     <- logRBAFSample1$chrom
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionOne]]$pos     <- logRBAFSample1$pos
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionOne]]$nA      <- as.numeric(ABSample1[,1])
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionOne]]$nB      <- as.numeric(ABSample1[,2])
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionOne]]$nA_mean <- mean(as.numeric(ABSample1[,1]))
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionOne]]$nB_mean <- mean(as.numeric(ABSample1[,2]))
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionTwo]]$nA      <- as.numeric(ABSample2[,1])
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionTwo]]$nB      <- as.numeric(ABSample2[,2])
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionTwo]]$nA_mean <- mean(as.numeric(ABSample2[,1]))
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionTwo]]$nB_mean <- mean(as.numeric(ABSample2[,2]))
      
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$nA_PairedPval             <- PValsPaired[1]
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$nB_PairedPval             <- PValsPaired[2]
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$nB_unPairedPval           <- PValsUnPaired[2]
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$nA_unPairedPval           <- PValsUnPaired[1]
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$combinedPairedPval        <- combinedPvalPaired
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$combinedUnPairedPval      <- combinedPvalUnPaired
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$combinedUnPairedPval      <- combinedPvalUnPaired
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$totalCPNPairedPval        <- totalCPNpval
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$nA_heterogeneity          <- nA_heterogeneity
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$nB_heterogeneity          <- nB_heterogeneity
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$combinedHeterogeneity     <- nB_heterogeneity+nA_heterogeneity
      # classify the segment as showing MSAI
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$msAI                      <- 0
      # if(PValsPaired[1]<0.001&PValsPaired[2]<0.001)
      # {
      
      if(tTestWithErrorHandling(ABSample1[,1],ABSample1[,2],paired=TRUE)$p.val<0.001&tTestWithErrorHandling(ABSample2[,1],ABSample2[,2],paired=TRUE)$p.val<0.001)
      {
        if((mean(ABSample1[,1]-ABSample1[,2])>0&mean(ABSample2[,1]-ABSample2[,2])<0)|(mean(ABSample1[,1]-ABSample1[,2])<0&mean(ABSample2[,1]-ABSample2[,2])>0))
        {
          ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$msAI                      <- 1
        }
        
      }
      
      
      
      checkGainSignificance <- function(logRBAFRegion,regionToCheck,ploidyEstimates,gainThreshold)
      {
        minGainRawExpRegion       <-  ploidyEstimates[regionToCheck]*(2^(gainThreshold))
        minGainLogRExpRegion      <-  funCalcExpectedLogR(minGainRawExpRegion,0,ploidyEstimates[regionToCheck],purityEstimates[regionToCheck])
        logRObs                   <- median(logRBAFRegion$logR)
        
        if(nrow(logRBAFRegion)>=3)
        {
          Ttest        <- tTestWithErrorHandling(logRBAFRegion$logR,mu = minGainLogRExpRegion)
          TtestGreater <- tTestWithErrorHandling(logRBAFRegion$logR,mu = minGainLogRExpRegion,alternative = 'greater')
          #TtestLess    <- tTestWithErrorHandling(logRBAF$logR,mu = minLossLogRExp,alternative = 'less')
          logRObs_low  <- Ttest$conf.int[1]
          logRObs_high <- Ttest$conf.int[2]
          PVal         <- Ttest$p.value
          PvalGain     <- TtestGreater$p.value
          
        }
        if(nrow(logRBAFRegion)<3)
        {
          logRObs_low  <- NA
          logRObs_high <- NA
          PVal         <- NA
          PvalGain     <- NA
          
          
        }
        
        out <- cbind(minGainRawExpRegion,minGainLogRExpRegion,logRObs,logRObs_low,logRObs_high,PVal,PvalGain)  
        return(out)
        
        
      }
      checkLossSignificance <- function(logRBAFRegion,regionToCheck,ploidyEstimates,lossThreshold)
      {
        minLossRawExpRegion       <-  ploidyEstimates[regionToCheck]*(2^(lossThreshold))
        minLossLogRExpRegion      <-  funCalcExpectedLogR(minLossRawExpRegion,0,ploidyEstimates[regionToCheck],purityEstimates[regionToCheck])
        logRObs                   <- median(logRBAFRegion$logR)
        
        if(nrow(logRBAFRegion)>=3)
        {
          Ttest        <- tTestWithErrorHandling(logRBAFRegion$logR,mu = minLossLogRExpRegion)
          TtestLess <- tTestWithErrorHandling(logRBAFRegion$logR,mu = minLossLogRExpRegion,alternative = 'less')
          #TtestLess    <- tTestWithErrorHandling(logRBAF$logR,mu = minLossLogRExp,alternative = 'less')
          logRObs_low  <- Ttest$conf.int[1]
          logRObs_high <- Ttest$conf.int[2]
          PVal         <- Ttest$p.value
          PvalLoss     <- TtestLess$p.value
          
        }
        if(nrow(logRBAFRegion)<3)
        {
          logRObs_low  <- NA
          logRObs_high <- NA
          PVal         <- NA
          PvalLoss     <- NA
          
          
        }
        
        out <- cbind(minLossRawExpRegion,minLossLogRExpRegion,logRObs,logRObs_low,logRObs_high,PVal,PvalLoss)  
        return(out)
        
        
      }
      
      
      GainRegionOne     <- checkGainSignificance(logRSample1,regionToCheck = regionOne,ploidyEstimates,gainThreshold = gainThreshold)
      AmpRegionOne      <- checkGainSignificance(logRSample1,regionToCheck = regionOne,ploidyEstimates,gainThreshold = ampThreshold)
      LossRegionOne     <- checkLossSignificance(logRSample1,regionToCheck = regionOne,ploidyEstimates,lossThreshold = lossThreshold)
      DeepLossRegionOne <- checkLossSignificance(logRSample1,regionToCheck = regionOne,ploidyEstimates,lossThreshold = deepLossThreshold)
      
      # put together into One table
      regionOneGainOrLoss <- c(GainRegionOne[,3:5,drop=FALSE]
                               ,GainRegionOne[,2,drop=FALSE],GainRegionOne[,7,drop=FALSE]
                               ,AmpRegionOne[,2,drop=FALSE],AmpRegionOne[,7,drop=FALSE]
                               ,LossRegionOne[,2,drop=FALSE],LossRegionOne[,7,drop=FALSE]
                               ,DeepLossRegionOne[,2,drop=FALSE],DeepLossRegionOne[,7,drop=FALSE]
      )
      names(regionOneGainOrLoss) <- c('logRObs','logRObs_low','logRObs_high'
                                      ,'minGainExpLogR','gainPval'
                                      ,'minAmpExpLogR','ampPval'
                                      ,'minLossExpLogR','lossPval'
                                      ,'minDeepLossExpLogR','DeepLossPval')
      
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionOne]]$segStatus     <- regionOneGainOrLoss
      
      
      GainRegionTwo     <- checkGainSignificance(logRSample2,regionToCheck = regionTwo,ploidyEstimates,gainThreshold = gainThreshold)
      AmpRegionTwo      <- checkGainSignificance(logRSample2,regionToCheck = regionTwo,ploidyEstimates,gainThreshold = ampThreshold)
      LossRegionTwo     <- checkLossSignificance(logRSample2,regionToCheck = regionTwo,ploidyEstimates,lossThreshold = lossThreshold)
      DeepLossRegionTwo <- checkLossSignificance(logRSample2,regionToCheck = regionTwo,ploidyEstimates,lossThreshold = deepLossThreshold)
      
      
      # put together into Two table
      regionTwoGainOrLoss <- c(GainRegionTwo[,3:5,drop=FALSE]
                               ,GainRegionTwo[,2,drop=FALSE],GainRegionTwo[,7,drop=FALSE]
                               ,AmpRegionTwo[,2,drop=FALSE],AmpRegionTwo[,7,drop=FALSE]
                               ,LossRegionTwo[,2,drop=FALSE],LossRegionTwo[,7,drop=FALSE]
                               ,DeepLossRegionTwo[,2,drop=FALSE],DeepLossRegionTwo[,7,drop=FALSE]
      )
      names(regionTwoGainOrLoss) <- c('logRObs','logRObs_low','logRObs_high'
                                      ,'minGainExpLogR','gainPval'
                                      ,'minAmpExpLogR','ampPval'
                                      ,'minLossExpLogR','lossPval'
                                      ,'minDeepLossExpLogR','DeepLossPval')
      
      ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionTwo]]$segStatus     <- regionTwoGainOrLoss
      
      # # for it to be classified as MSAI the difference between A1-A2 must be positive and B1-B2 negative (or the other way round)
      # if((mean(as.numeric(ABSample1[,1]))-mean(as.numeric(ABSample2[,1]))>0&mean(as.numeric(ABSample1[,2]))-mean(as.numeric(ABSample2[,2]))<0)|(mean(as.numeric(ABSample1[,1]))-mean(as.numeric(ABSample2[,1]))<0&mean(as.numeric(ABSample1[,2]))-mean(as.numeric(ABSample2[,2]))>0))
      # {
      #   # we also must observe significant difference in A1 and B1 as well as A2 and B2. And this significant difference must be opposite. 
      #   
      # }
      # }
      
      
      
      
    }
    
    averageHeterogeneity <- sum(HeterogeneityScore)/sum(SizesOfSegments)
    
    #ComparisonOutPut$`1:762273:11831543`$`L_LTX064_SU_T1-R6:L_LTX064_SU_T1-R7`$nA_heterogeneity
    #ComparisonOutPut$`1:762273:11831543`$`L_LTX064_SU_T1-R6:L_LTX064_SU_T1-R7`$nB_heterogeneity
    ComparisonOutPut[[paste(regionOne,regionTwo,sep=":")]]$HeterogeneityScore <- averageHeterogeneity 
    
  }
  
  
  if(plotting==TRUE)
  {
    
    par(mfrow=c(4,1))
    
    #major Allele####
    chrom.length.copy   <- fun.chrom.length(min.cons.seg[,1:3])
    names(chrom.length.copy)   <- 1:22
    chrom.segs          <- fun.add.chrom(min.cons.seg[,1:3],chrom.length.copy)
    min.x               <- 0
    min.y               <- minYaxis
    max.y               <- maxYaxis
    max.x               <- as.numeric(max(chrom.segs[,'End']))
    par(mar=c(2,5,5,5))
    plot(1
         ,xlim=c(min.x,max.x)
         ,ylim=c(min.y,max.y)
         ,type='n'
         ,xaxt='n'
         ,yaxs='i'
         ,yaxt='n'
         ,xaxs='i'
         ,ylab=""
         ,xlab=""
         ,lwd=2
         ,main=paste("Major allele ", regionOne," vs ",regionTwo,sep="")
         ,bty="n")
    box()
    
    for (x in 1:length(names_ComparisonOutPut))
    {
      ChrCoords <- unlist(strsplit(names_ComparisonOutPut[x],split=":"))
      toPlot <- cbind(ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionOne]]$chr
                      ,ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionOne]]$pos
                      ,ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionOne]]$pos
                      ,ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionOne]]$nA
                      ,ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionOne]]$nB
                      ,ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionTwo]]$nA
                      ,ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionTwo]]$nB)
      if(is.na(toPlot[1,1]))
      {
        next;
      }
      
      #toPlot <- cbind(logRBAFSample1$chrom,logRBAFSample1$pos,logRBAFSample1$pos,ABSample1,ABSample2)
      toPlot[,1:3] <- fun.add.chrom(apply(toPlot[,1:3],2,as.numeric),chrom.length.copy)
      toPlot <- data.frame(toPlot,stringsAsFactors = FALSE)
      colnames(toPlot) <- c('chr','pos','pos','Sample1nA','Sample1nB','Sample2nA','Sample2nB')
      PValsPaired <- as.numeric(ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$nA_PairedPval)
      totalPval   <-  ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$totalCPNPairedPval#as.numeric(ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$nA_PairedPval)
      
      
      MSAIstatus  <- as.numeric(ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$msAI)
      points(x = as.numeric(toPlot$pos),y=as.numeric(toPlot$Sample1nA)-as.numeric(toPlot$Sample2nA),col=ifelse(MSAIstatus==1,'#238b4550',ifelse(PValsPaired<0.001&totalPval<0.001,'#3182bd50','#fc927250')),cex=0.75,pch=16)
      par(lend=1)
      segments(x0=min(as.numeric(toPlot$pos))
               ,y0=mean(c(as.numeric(toPlot$Sample1nA)-as.numeric(toPlot$Sample2nA)))
               ,x1=max(as.numeric(toPlot$pos)),y1=mean(c(as.numeric(toPlot$Sample1nA)-as.numeric(toPlot$Sample2nA)))
               ,col=ifelse(MSAIstatus==1,'#00441b',ifelse(PValsPaired<0.001&totalPval<0.001,'#08306b','#a50f15')),lwd=3)
    }
    
    axis(side =2,las=1)
    abline(h=seq((minYaxis+1),(maxYaxis-1),by=1),lty='dashed',lwd=2,col='gray')
    
    mtext(side=3
          ,at=fun.chrom.mean(chrom.segs[,1:3])
          ,text=sort(unique(chrom.segs[,1]))
          ,cex=seq(0.6,0.4,length.out=length(unique(chrom.segs[,1])))
          ,line=-1
          #,lwd=0.5
    )
    
    for (i in sort(unique(as.numeric(chrom.segs[,1]))))
    {
      abline(v=max(as.numeric(chrom.segs[chrom.segs[,1]==i,3]))
             ,lwd=0.5)
      #,lty='dashed')
    }
    
    #minor Allele ####
    par(mar=c(5,5,2,5))
    chrom.length.copy   <- fun.chrom.length(min.cons.seg[,1:3])
    names(chrom.length.copy)   <- 1:22
    chrom.segs          <- fun.add.chrom(min.cons.seg[,1:3],chrom.length.copy)
    min.x               <- 0
    min.y               <- minYaxis
    max.y               <- maxYaxis
    max.x               <- as.numeric(max(chrom.segs[,'End']))
    
    plot(1
         ,xlim=c(min.x,max.x)
         ,ylim=c(min.y,max.y)
         ,type='n'
         ,xaxt='n'
         ,yaxs='i'
         ,yaxt='n'
         ,xaxs='i'
         ,ylab=""
         ,xlab=""
         ,lwd=2
         ,main=paste("Minor allele ", regionOne," vs ",regionTwo,sep="")
         ,bty="n")
    box()
    
    
    for (x in 1:length(names_ComparisonOutPut))
    {
      ChrCoords <- unlist(strsplit(names_ComparisonOutPut[x],split=":"))
      toPlot <- cbind(ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionOne]]$chr
                      ,ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionOne]]$pos
                      ,ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionOne]]$pos
                      ,ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionOne]]$nA
                      ,ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionOne]]$nB
                      ,ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionTwo]]$nA
                      ,ComparisonOutPut[[names_ComparisonOutPut[x]]][[regionTwo]]$nB)
      if(is.na(toPlot[1,1]))
      {
        next;
      }
      
      #toPlot <- cbind(logRBAFSample1$chrom,logRBAFSample1$pos,logRBAFSample1$pos,ABSample1,ABSample2)
      toPlot[,1:3] <- fun.add.chrom(apply(toPlot[,1:3],2,as.numeric),chrom.length.copy)
      toPlot <- data.frame(toPlot,stringsAsFactors = FALSE)
      colnames(toPlot) <- c('chr','pos','pos','Sample1nA','Sample1nB','Sample2nA','Sample2nB')
      PValsPaired <- as.numeric(ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$nB_PairedPval)
      totalPval   <-  ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$totalCPNPairedPval#as.numeric(ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$nA_PairedPval)
      MSAIstatus  <- as.numeric(ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$msAI)
      points(x = as.numeric(toPlot$pos),y=as.numeric(toPlot$Sample1nB)-as.numeric(toPlot$Sample2nB),col=ifelse(MSAIstatus==1,'#238b4550',ifelse(PValsPaired<0.001&totalPval<0.001,'#3182bd50','#fc927250')),cex=0.75,pch=16)
      par(lend=1)
      segments(x0=min(as.numeric(toPlot$pos))
               ,y0=mean(c(as.numeric(toPlot$Sample1nB)-as.numeric(toPlot$Sample2nB)))
               ,x1=max(as.numeric(toPlot$pos)),y1=mean(c(as.numeric(toPlot$Sample1nB)-as.numeric(toPlot$Sample2nB)))
               ,col=ifelse(MSAIstatus==1,'#00441b',ifelse(PValsPaired<0.001&totalPval<0.001,'#08306b','#a50f15')),lwd=3)
      
      
      
      #  points(x = as.numeric(toPlot$pos),y=as.numeric(toPlot$Sample1nB)-as.numeric(toPlot$Sample2nB),col=ifelse(PValsPaired<0.001,'#3182bd50','#fb6a4a50'),cex=0.75,pch=16)
      # segments(x0=min(as.numeric(toPlot$pos))
      #          ,y0=mean(c(as.numeric(toPlot$Sample1nB)-as.numeric(toPlot$Sample2nB)))
      #          ,x1=max(as.numeric(toPlot$pos)),y1=mean(c(as.numeric(toPlot$Sample1nB)-as.numeric(toPlot$Sample2nB)))
      #          ,col=ifelse(PValsPaired<0.001,'#08519c','#a50f15'),lwd=3)
    }
    
    
    axis(side =2,las=1)
    
    mtext(side=3
          ,at=fun.chrom.mean(chrom.segs[,1:3])
          ,text=sort(unique(chrom.segs[,1]))
          ,cex=seq(0.6,0.4,length.out=length(unique(chrom.segs[,1])))
          ,line=-1
          #,lwd=0.5
    )
    
    for (i in sort(unique(as.numeric(chrom.segs[,1]))))
    {
      abline(v=max(as.numeric(chrom.segs[chrom.segs[,1]==i,3]))
             ,lwd=0.5)
      #,lty='dashed')
    }
    abline(h=seq((minYaxis+1),(maxYaxis-1),by=1),lty='dashed',lwd=2,col='gray')
    
    
    
    
    
  }
  print("Returning from calcDiffBetweenTwoRegions")
  invisible(ComparisonOutPut)
  
}



funCalcExpectedLogR <- function(exp_intA,exp_intB,tumorPloidy,tumorPurity)
{
  print("Called funCalcExpectedLogR")
  expectedLogR <- log2((2*(1-tumorPurity)+tumorPurity*(exp_intA+exp_intB))/(2*(1-tumorPurity)+(tumorPurity*tumorPloidy)))
  print("Returning from funCalcExpectedLogR")
  return(expectedLogR)
}

funCalcExpectedBAF <- function(exp_intA,exp_intB,tumorPloidy,tumorPurity)
{
  print("Called funCalcExpectedBAF")
  expectedBAF <- (1-tumorPurity+(tumorPurity*exp_intB))/(2-2*tumorPurity+tumorPurity*(exp_intA+exp_intB))
  print("Returning from funCalcExpectedBAF")
  return(expectedBAF)
}


fun.chrom.length<- function(seg){
  # This function returns the length of chromosome 1,...,22
  # seg: first three columns of the minimum conistent region matrix
  print("Called fun.chrom.length")
  chrom.length <- c()
  for(i in 1:22){
    sub.chrom <- subset(seg, seg[, 1] == i)
    chrom.length <- c(chrom.length, sub.chrom[nrow(sub.chrom), 3])
    
  }
  print("Returning from fun.chrom.length")
  return(chrom.length)
  
}


fun.add.chrom<-function(seg, chrom.length){
  print("Called fun.add.chrom")
  # This function adds the length of chromosome 1,...,(i - 1) to chromsomes 2,...,22
  # seg: first three columns of the minimum conistent region matrix
  # chrom.length: output of fun.chrom.length(...)
  if(1 %in% seg[, 1])
    seg.tmp <-  subset(seg, seg[, 1] == 1)
  else
    seg.tmp <- matrix(NA, nr = 0, nc = ncol(seg))  
  for(i in 2:22){
    if(!(i %in% seg[, 1]))
      next()
    sub.chrom <- subset(seg, seg[, 1] == i)
    sub.chrom[, 2] <- sub.chrom[, 2] + sum(chrom.length[1:(i - 1)])
    sub.chrom[, 3] <- sub.chrom[, 3] + sum(chrom.length[1:(i - 1)])
    seg.tmp <- rbind(seg.tmp, sub.chrom)
  }
  print("Returning from fun.add.chrom")
  return(seg.tmp)
  
}

fun.chrom.mean<-function(seg){
  print("Called fun.chrom.mean")
  # This function finds the middle of the positions for each chromosome.
  # seg: Matrix with three columns, containing the positions.
  chrom.mean <- c()
  for(i in 1:22){
    sub.chrom <- subset(seg, seg[, 1] == i)
    chrom.mean <- c(chrom.mean, (min(sub.chrom[, 2]) + max(sub.chrom[, 3]))/2)
    
  }
  print("Returning from fun.chrom.mean")
  return(chrom.mean)
}


###########################################################################Only  funPlotAllRegions_withJoining plotting scripts:
funPlotAllRegions_withJoining <- function(copyITHRDS)
{
  #str(copyITHRDS)
  #str(copyITHRDS$newBAFtable)
  
  # first, get the min.cons.regions
  names_ComparisonOutPut <- names(copyITHRDS[[2]][[1]])
  min.cons.seg <- names(copyITHRDS[[2]][[1]])
  min.cons.seg <- do.call(rbind,strsplit(min.cons.seg,split=":"))
  min.cons.seg <- apply(min.cons.seg,2,as.numeric)
  # next plot each sample individually
  chrom.length.copy   <- fun.chrom.length(min.cons.seg[,1:3])
  names(chrom.length.copy)   <- 1:22
  chrom.segs          <- fun.add.chrom(min.cons.seg[,1:3],chrom.length.copy)
  min.x               <- 0
  min.y               <- -0.5
  max.y               <- 8.5
  max.x               <- as.numeric(max(chrom.segs[,3]))
  
  # how many samples are we dealing with
  Regions     <- grep(":",names(copyITHRDS$regionComparisonList),invert=TRUE,value=TRUE)
  numRegions  <- length(Regions)
  firstRegion <- Regions[1] 
  lastRegion <- Regions[length(Regions)] 
  regionJoinedSeg <- funJoinAdjacentSegments(copyITHRDS$regionComparisonList)
  
  
  par(mfrow=c(numRegions,1))
  
  for (region in Regions)
  {
    if(region ==firstRegion)
    {
      par(mar=c(2,3,4,3))
    }
    
    if(region ==lastRegion)
    {
      par(mar=c(4,3,2,3))
    }
    
    if(!region %in%c(firstRegion,lastRegion))
    {
      par(mar=c(3,3,3,3))
    }
    plot(1
         ,xlim=c(min.x,max.x)
         ,ylim=c(min.y,max.y)
         ,type='n'
         ,xaxt='n'
         ,yaxs='i'
         ,yaxt='n'
         ,xaxs='i'
         ,ylab=""
         ,xlab=""
         ,lwd=2
         ,main=paste(region," A and B alleles",sep="")
         ,bty="n")
    box()
    
    # get the regionSeg
    regionSeg <- regionJoinedSeg[[region]]
    regionSegsToUse_nA <- names(regionSeg$segs_nA)
    regionSegsToUse_nB <- names(regionSeg$segs_nB)
    
    for (x in 1:length(regionSegsToUse_nA))
    {
      ChrCoords <- unlist(strsplit(regionSegsToUse_nA[x],split=":"))
      
      toPlot <- cbind(regionSeg$segs_nA[[regionSegsToUse_nA[x]]]$chr
                      ,regionSeg$segs_nA[[regionSegsToUse_nA[x]]]$pos
                      ,regionSeg$segs_nA[[regionSegsToUse_nA[x]]]$pos
                      ,regionSeg$segs_nA[[regionSegsToUse_nA[x]]]$nA
                      ,regionSeg$segs_nA[[regionSegsToUse_nA[x]]]$nA
                      ,regionSeg$segs_nA[[regionSegsToUse_nA[x]]]$nA
                      ,regionSeg$segs_nA[[regionSegsToUse_nA[x]]]$nA)
      # toPlot <- cbind(copyITHRDS$regionComparisonList[[firstRegion]][[regionSegsToUse_nA[x]]]$chr
      #                 ,copyITHRDS$regionComparisonList[[firstRegion]][[regionSegsToUse_nA[x]]]$pos
      #                 ,copyITHRDS$regionComparisonList[[firstRegion]][[regionSegsToUse_nA[x]]]$pos
      #                 ,copyITHRDS$regionComparisonList[[region]][[regionSegsToUse_nA[x]]]$nA
      #                 ,copyITHRDS$regionComparisonList[[region]][[regionSegsToUse_nA[x]]]$nA
      #                 ,copyITHRDS$regionComparisonList[[region]][[regionSegsToUse_nA[x]]]$nA
      #                 ,copyITHRDS$regionComparisonList[[region]][[regionSegsToUse_nA[x]]]$nA)
      
      
      if(is.na(toPlot[1,1]))
      {
        next;
      }
      
      #toPlot <- cbind(logRBAFSample1$chrom,logRBAFSample1$pos,logRBAFSample1$pos,ABSample1,ABSample2)
      toPlot[,1:3] <- fun.add.chrom(apply(toPlot[,1:3],2,as.numeric),chrom.length.copy)
      toPlot <- data.frame(toPlot,stringsAsFactors = FALSE)
      colnames(toPlot) <- c('chr','pos','pos','Sample1nA','Sample1nB','Sample2nA','Sample2nB')
      
      toPlot$Sample1nA <- ifelse(as.numeric(toPlot$Sample1nA)>max.y,(max.y-0.25),as.numeric(toPlot$Sample1nA))
      toPlot$Sample1nA <- ifelse(as.numeric(toPlot$Sample1nA)<min.y,(min.y+0.25),as.numeric(toPlot$Sample1nA))
      
      toPlot$Sample1nB <- ifelse(as.numeric(toPlot$Sample1nB)>max.y,(max.y-0.25),as.numeric(toPlot$Sample1nB))
      toPlot$Sample1nB <- ifelse(as.numeric(toPlot$Sample1nB)<min.y,(min.y+0.25),as.numeric(toPlot$Sample1nB))
      
      #PValsPaired <- as.numeric(ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$nA_PairedPval)
      #totalPval   <-  ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$totalCPNPairedPval#as.numeric(ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$nA_PairedPval)
      
      
      #MSAIstatus  <- as.numeric(ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$msAI)
      points(x = as.numeric(toPlot$pos),y=as.numeric(toPlot$Sample1nA),col='#2ca25f50'
             ,cex=.5
             ,pch=16)#ifelse(MSAIstatus==1,'#238b4550',ifelse(PValsPaired<0.001&totalPval<0.001,'#3182bd50','#fc927250')),cex=0.75,pch=16)
      par(lend=1)
      segments(x0=min(as.numeric(toPlot$pos))
               ,y0=mean(c(as.numeric(toPlot$Sample1nA)))
               ,x1=max(as.numeric(toPlot$pos)),y1=mean(c(as.numeric(toPlot$Sample1nA)))
               ,col='#00441b'
               ,lwd=3)#ifelse(MSAIstatus==1,'#00441b',ifelse(PValsPaired<0.001&totalPval<0.001,'#08306b','#a50f15')),lwd=3)
      
      # points(x = as.numeric(toPlot$pos),y=as.numeric(toPlot$Sample1nB),col='#4292c650'
      #        ,cex=0.75
      #        ,pch=16)#ifelse(MSAIstatus==1,'#238b4550',ifelse(PValsPaired<0.001&totalPval<0.001,'#3182bd50','#fc927250')),cex=0.75,pch=16)
      # par(lend=1)
      # segments(x0=min(as.numeric(toPlot$pos))
      #          ,y0=mean(c(as.numeric(toPlot$Sample1nB)))
      #          ,x1=max(as.numeric(toPlot$pos)),y1=mean(c(as.numeric(toPlot$Sample1nB)))
      #          ,col='#08306b'
      #          ,lwd=3)#ifelse(MSAIstatus==1,'#00441b',ifelse(PValsPaired<0.001&totalPval<0.001,'#08306b','#a50f15')),lwd=3)
      
    }
    
    for (x in 1:length(regionSegsToUse_nB))
    {
      ChrCoords <- unlist(strsplit(regionSegsToUse_nB[x],split=":"))
      
      toPlot <- cbind(regionSeg$segs_nB[[regionSegsToUse_nB[x]]]$chr
                      ,regionSeg$segs_nB[[regionSegsToUse_nB[x]]]$pos
                      ,regionSeg$segs_nB[[regionSegsToUse_nB[x]]]$pos
                      ,regionSeg$segs_nB[[regionSegsToUse_nB[x]]]$nB
                      ,regionSeg$segs_nB[[regionSegsToUse_nB[x]]]$nB
                      ,regionSeg$segs_nB[[regionSegsToUse_nB[x]]]$nB
                      ,regionSeg$segs_nB[[regionSegsToUse_nB[x]]]$nB)
      # toPlot <- cbind(copyITHRDS$regionComparisonList[[firstRegion]][[regionSegsToUse_nB[x]]]$chr
      #                 ,copyITHRDS$regionComparisonList[[firstRegion]][[regionSegsToUse_nB[x]]]$pos
      #                 ,copyITHRDS$regionComparisonList[[firstRegion]][[regionSegsToUse_nB[x]]]$pos
      #                 ,copyITHRDS$regionComparisonList[[region]][[regionSegsToUse_nB[x]]]$nA
      #                 ,copyITHRDS$regionComparisonList[[region]][[regionSegsToUse_nB[x]]]$nA
      #                 ,copyITHRDS$regionComparisonList[[region]][[regionSegsToUse_nB[x]]]$nA
      #                 ,copyITHRDS$regionComparisonList[[region]][[regionSegsToUse_nB[x]]]$nA)
      
      
      if(is.na(toPlot[1,1]))
      {
        next;
      }
      
      #toPlot <- cbind(logRBAFSample1$chrom,logRBAFSample1$pos,logRBAFSample1$pos,ABSample1,ABSample2)
      toPlot[,1:3] <- fun.add.chrom(apply(toPlot[,1:3],2,as.numeric),chrom.length.copy)
      toPlot <- data.frame(toPlot,stringsAsFactors = FALSE)
      colnames(toPlot) <- c('chr','pos','pos','Sample1nA','Sample1nB','Sample2nA','Sample2nB')
      
      toPlot$Sample1nA <- ifelse(as.numeric(toPlot$Sample1nA)>max.y,(max.y-0.25),as.numeric(toPlot$Sample1nA))
      toPlot$Sample1nA <- ifelse(as.numeric(toPlot$Sample1nA)<min.y,(min.y+0.25),as.numeric(toPlot$Sample1nA))
      
      toPlot$Sample1nB <- ifelse(as.numeric(toPlot$Sample1nB)>max.y,(max.y-0.25),as.numeric(toPlot$Sample1nB))
      toPlot$Sample1nB <- ifelse(as.numeric(toPlot$Sample1nB)<min.y,(min.y+0.25),as.numeric(toPlot$Sample1nB))
      
      #PValsPaired <- as.numeric(ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$nA_PairedPval)
      #totalPval   <-  ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$totalCPNPairedPval#as.numeric(ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$nA_PairedPval)
      
      
      #MSAIstatus  <- as.numeric(ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$msAI)
      # points(x = as.numeric(toPlot$pos),y=as.numeric(toPlot$Sample1nA),col='#2ca25f50'
      #        ,cex=0.75
      #        ,pch=16)#ifelse(MSAIstatus==1,'#238b4550',ifelse(PValsPaired<0.001&totalPval<0.001,'#3182bd50','#fc927250')),cex=0.75,pch=16)
      # par(lend=1)
      # segments(x0=min(as.numeric(toPlot$pos))
      #          ,y0=mean(c(as.numeric(toPlot$Sample1nA)))
      #          ,x1=max(as.numeric(toPlot$pos)),y1=mean(c(as.numeric(toPlot$Sample1nA)))
      #          ,col='#00441b'
      #          ,lwd=3)#ifelse(MSAIstatus==1,'#00441b',ifelse(PValsPaired<0.001&totalPval<0.001,'#08306b','#a50f15')),lwd=3)
      
      points(x = as.numeric(toPlot$pos),y=as.numeric(toPlot$Sample1nB),col='#4292c650'
             ,cex=0.5
             ,pch=16)#ifelse(MSAIstatus==1,'#238b4550',ifelse(PValsPaired<0.001&totalPval<0.001,'#3182bd50','#fc927250')),cex=0.75,pch=16)
      par(lend=1)
      segments(x0=min(as.numeric(toPlot$pos))
               ,y0=mean(c(as.numeric(toPlot$Sample1nB)))
               ,x1=max(as.numeric(toPlot$pos)),y1=mean(c(as.numeric(toPlot$Sample1nB)))
               ,col='#08306b'
               ,lwd=3)#ifelse(MSAIstatus==1,'#00441b',ifelse(PValsPaired<0.001&totalPval<0.001,'#08306b','#a50f15')),lwd=3)
      
    }
    
    axis(side =2,las=1)
    abline(h=seq((min.y+0.5),(max.y-0.5),by=1),lty='dashed',lwd=0.75,col='gray')
    
    mtext(side=3
          ,at=fun.chrom.mean(chrom.segs[,1:3])
          ,text=sort(unique(chrom.segs[,1]))
          ,cex=seq(0.6,0.4,length.out=length(unique(chrom.segs[,1])))
          ,line=-1
          #,lwd=0.5
    )
    
    for (i in sort(unique(as.numeric(chrom.segs[,1]))))
    {
      abline(v=max(as.numeric(chrom.segs[chrom.segs[,1]==i,3]))
             ,lwd=0.5)
      #,lty='dashed')
    }
  }
  
  newBAFtable <- copyITHRDS$newBAFtable
  # next, plot the BAF profiles and see whether these are improved with rescue phasing
  for (region in Regions)
  {
    print(region)
    test <- newBAFtable[[region]]$regionBAFtable
    #str(newBAFtable)
    #str(test)
    test$col <- ifelse(test$nA==1,'#de2d2699','#3182bd99')
    firstThreeColumns <- apply(cbind(test$chrom,test$pos,test$pos),2,as.numeric)
    chrom.length.copy   <- fun.chrom.length(firstThreeColumns[,1:3])
    chrom.segs          <- fun.add.chrom(firstThreeColumns[,1:3],chrom.length.copy)
    
    
    test[test$phasedFrom%in%'randomSelection',]$newBAFsegmented <- 0.5
    rownames(test) <- 1:nrow(test)
    plot(chrom.segs[,2],test$BAF,col=test$col,pch=16,ylim=c(0,1),cex=0.5,main=region,ylab='BAF',xlab='',las=1,xaxt='n')
    mtext(side=3
          ,at=fun.chrom.mean(chrom.segs[,1:3])
          ,text=sort(unique(chrom.segs[,1]))
          ,cex=seq(0.6,0.4,length.out=length(unique(chrom.segs[,1])))
          ,line=-1
          #,lwd=0.5
    )
    for (i in sort(unique(as.numeric(chrom.segs[,1]))))
    {
      abline(v=max(as.numeric(chrom.segs[chrom.segs[,1]==i,3]))
             ,lwd=0.5)
      #,lty='dashed')
    }
    
    # par(lend=2)
    points(chrom.segs[,2],test$BAFsegmented,col='black',ylim=c(0,1),pch=15,cex=0.25)
    points(chrom.segs[,2],test$newBAFsegmented,col='#78c679',ylim=c(0,1),pch=15,cex=0.5)
    #par(lend=2)
    points(chrom.segs[,2],1-test$BAFsegmented,col='black',ylim=c(0,1),pch=15,cex=0.25)
    points(chrom.segs[,2],1-test$newBAFsegmented,col='#78c679',ylim=c(0,1),pch=15,cex=0.5)
    
    #test2 <- test[as.numeric(rownames(test))>22000&as.numeric(rownames(test))<25000&test$BAFsegmented==0.5,,drop=FALSE]
    #test <- test[!test$phasedFrom%in%'randomSelection',,drop=FALSE]
    #test <- test[test$BAFsegmented==0.5,,drop=FALSE]
    # plot(rownames(test[test$nA==1,,drop=FALSE]),test$BAF[test$nA==1],col=ifelse(test[test$nA==1,]$phasedFrom=='randomSelection','#bdbdbd99','#de2d2699'),pch=16,ylim=c(0,1))
    # points(rownames(test[test$nA==1,,drop=FALSE]),1-test$BAFsegmented[test$nA==1],col=ifelse(test[test$nA==1,]$phasedFrom=='randomSelection','#bdbdbd99','#31a354'),pch=16,ylim=c(0,1))
    # points(rownames(test[test$nA==0,,drop=FALSE]),test$BAF[test$nA==0],col=ifelse(test[test$nA==0,]$phasedFrom=='randomSelection','#bdbdbd99','#3182bd99'),pch=16)
    #points(rownames(test[test$nA==0,,drop=FALSE]),test$BAFsegmented[test$nA==0],col=ifelse(test[test$nA==0,]$phasedFrom=='randomSelection','#bdbdbd99','#31a354'),pch=16,ylim=c(0,1))
  }
  
  
  
}

######################################## Only helpers to funPlotAllRegions_withJoining

funJoinAdjacentSegments <- function(regionComparisonList
                                    ,PvalToSeparate=0.005)
{
  print("Called funJoinAdjacentSegments")
  regionsToConsider <- grep(":",names(regionComparisonList),invert=TRUE,value=TRUE)
  segsToConsider    <- names(regionComparisonList[[1]])
  
  newList <- list()
  segsToUse <- c()
  for (seg in segsToConsider)
  {
    if(is.na(regionComparisonList[[1]][[seg]]$nB))
    {
      next;
    }
    segsToUse <- c(segsToUse,seg)
  }
  segsToConsider <- segsToUse
  
  
  regionOne <- regionsToConsider[1]
  for (region in regionsToConsider)
  {
    
    #first check segs for each chromosome
    segTable <- do.call(rbind,strsplit(segsToConsider,split=":"))
    segTable <- apply(segTable,2,as.numeric)
    rownames(segTable) <- segsToConsider
    
    regionSegs_nA <- list()
    regionSegs_nB <- list()
    for (chr in unique(segTable[,1]))
    {
      segTableSub <- segTable[segTable[,1]%in%chr,,drop=FALSE]
      segsToConsiderSubs <- rownames(segTableSub)
      
      # start with nB
      segOneB <- regionComparisonList[[region]][[segsToConsiderSubs[1]]]$nB
      segOneBpos <- regionComparisonList[[regionOne]][[segsToConsiderSubs[1]]]$pos
      segOneBchr <- regionComparisonList[[regionOne]][[segsToConsiderSubs[1]]]$chr
      
      segTwoB <- regionComparisonList[[region]][[segsToConsiderSubs[2]]]$nB
      segTwoBpos <- regionComparisonList[[regionOne]][[segsToConsiderSubs[2]]]$pos
      segTwoBchr <- regionComparisonList[[regionOne]][[segsToConsiderSubs[2]]]$chr
      
      segOneName <- segsToConsiderSubs[1]
      segTwoName <- segsToConsiderSubs[2]
      
      segsThatAreUsed <- c()
      segsThatAreUsedB  <- c()
      #check if less than 3 segments
      if (length(segsToConsiderSubs)==2)
      {
        k <- 1
        if(tTestWithErrorHandling(segOneB,segTwoB)$p.val>PvalToSeparate)
        {
          segName <- paste(unlist(strsplit(segOneName,split=":"))[1],unlist(strsplit(segOneName,split=":"))[2],unlist(strsplit(segTwoName,split=":"))[3],sep=":")
          segB    <- c(segOneB,segTwoB)
          
          
          #regionSegs[[segName]]$nB <- segB  
          segOneBpos <- c(segOneBpos,segTwoBpos)
          segOneBchr <- c(segOneBchr,segTwoBchr)
          segOneB    <- segB
          segOneName <- segName
          regionSegs_nB[[segOneName]]$nB <- segOneB
          regionSegs_nB[[segOneName]]$chr <- segOneBchr 
          regionSegs_nB[[segOneName]]$pos <- segOneBpos 
          segsThatAreUsed <- c(segsThatAreUsed,segOneName)
          
        }
        
        if(tTestWithErrorHandling(segOneB,segTwoB)$p.val<=PvalToSeparate)
        {
          regionSegs_nB[[segOneName]]$nB <- segOneB 
          regionSegs_nB[[segOneName]]$chr <- segOneBchr 
          regionSegs_nB[[segOneName]]$pos <- segOneBpos 
          
          regionSegs_nB[[segTwoName]]$nB <- segTwoB 
          regionSegs_nB[[segTwoName]]$chr <- segTwoBchr 
          regionSegs_nB[[segTwoName]]$pos <- segTwoBpos 
          segsThatAreUsed <- c(segsThatAreUsed,segOneName,segTwoName)
          
        }
        
      }
      
      if (length(segsToConsiderSubs)==1)
      {
        #regionSegs[[segName]]$nB <- segB  
        # segOneBpos <- c(segOneBpos,segTwoBpos)
        # segOneBchr <- c(segOneBchr,segTwoBchr)
        # segOneB    <- segB
        # segOneName <- segName
        regionSegs_nB[[segOneName]]$nB <- segOneB
        regionSegs_nB[[segOneName]]$chr <- segOneBchr 
        regionSegs_nB[[segOneName]]$pos <- segOneBpos 
        segsThatAreUsed <- c(segsThatAreUsed,segOneName)
        
        
      }
      
      k <- 3
      while(k <= length(segsToConsiderSubs))
      {
        #print(segsToConsiderSubs[k])
        if(tTestWithErrorHandling(segOneB,segTwoB)$p.val>PvalToSeparate)
        {
          #seg one and two can be joined 
          
          #combine the segments
          segName <- paste(unlist(strsplit(segOneName,split=":"))[1],unlist(strsplit(segOneName,split=":"))[2],unlist(strsplit(segTwoName,split=":"))[3],sep=":")
          segB    <- c(segOneB,segTwoB)
          
          
          #regionSegs[[segName]]$nB <- segB  
          segOneBpos <- c(segOneBpos,segTwoBpos)
          segOneBchr <- c(segOneBchr,segTwoBchr)
          segOneB    <- segB
          segOneName <- segName
          
          
          segTwoB   <-  regionComparisonList[[region]][[segsToConsiderSubs[k]]]$nB
          segTwoBchr   <-  regionComparisonList[[regionOne]][[segsToConsiderSubs[k]]]$chr
          segTwoBpos   <-  regionComparisonList[[regionOne]][[segsToConsiderSubs[k]]]$pos
          segTwoName <- segsToConsiderSubs[k]
          
          k <- k+1
          
          # check whether k exceeds number of segments
          if (k>length(segsToConsiderSubs))
          {
            if(tTestWithErrorHandling(segOneB,segTwoB)$p.val>PvalToSeparate)
            {
              segName <- paste(unlist(strsplit(segOneName,split=":"))[1],unlist(strsplit(segOneName,split=":"))[2],unlist(strsplit(segTwoName,split=":"))[3],sep=":")
              segB    <- c(segOneB,segTwoB)
              
              
              #regionSegs[[segName]]$nB <- segB  
              segOneBpos <- c(segOneBpos,segTwoBpos)
              segOneBchr <- c(segOneBchr,segTwoBchr)
              segOneB    <- segB
              segOneName <- segName
              regionSegs_nB[[segOneName]]$nB <- segOneB
              regionSegs_nB[[segOneName]]$chr <- segOneBchr 
              regionSegs_nB[[segOneName]]$pos <- segOneBpos 
              segsThatAreUsed <- c(segsThatAreUsed,segOneName)
              
            }
            
            if(tTestWithErrorHandling(segOneB,segTwoB)$p.val<=PvalToSeparate)
            {
              regionSegs_nB[[segOneName]]$nB <- segOneB 
              regionSegs_nB[[segOneName]]$chr <- segOneBchr 
              regionSegs_nB[[segOneName]]$pos <- segOneBpos 
              
              regionSegs_nB[[segTwoName]]$nB <- segTwoB 
              regionSegs_nB[[segTwoName]]$chr <- segTwoBchr 
              regionSegs_nB[[segTwoName]]$pos <- segTwoBpos 
              segsThatAreUsed <- c(segsThatAreUsed,segOneName,segTwoName)
              
            }
            
          }
          
        }
        
        if(tTestWithErrorHandling(segOneB,segTwoB)$p.val<=PvalToSeparate&k<=length(segsToConsiderSubs))
        {
          segsThatAreUsed <- c(segsThatAreUsed,segOneName)
          #print(k)
          # now add this to the regionSeg
          regionSegs_nB[[segOneName]]$nB <- segOneB 
          regionSegs_nB[[segOneName]]$chr <- segOneBchr 
          regionSegs_nB[[segOneName]]$pos <- segOneBpos 
          
          
          segOneB  <- segTwoB
          segOneName <- segTwoName
          segOneBpos <- segTwoBpos
          segOneBchr <- segTwoBchr
          
          segTwoB  <- regionComparisonList[[region]][[segsToConsiderSubs[k]]]$nB
          segTwoBchr   <-  regionComparisonList[[regionOne]][[segsToConsiderSubs[k]]]$chr
          segTwoBpos   <-  regionComparisonList[[regionOne]][[segsToConsiderSubs[k]]]$pos
          segTwoName  <- segsToConsiderSubs[k]
          
          k <- k+1
          if (k>length(segsToConsiderSubs))
          {
            if(tTestWithErrorHandling(segOneB,segTwoB)$p.val>PvalToSeparate)
            {
              segName <- paste(unlist(strsplit(segOneName,split=":"))[1],unlist(strsplit(segOneName,split=":"))[2],unlist(strsplit(segTwoName,split=":"))[3],sep=":")
              segB    <- c(segOneB,segTwoB)
              
              
              #regionSegs[[segName]]$nB <- segB  
              segOneBpos <- c(segOneBpos,segTwoBpos)
              segOneBchr <- c(segOneBchr,segTwoBchr)
              segOneB    <- segB
              segOneName <- segName
              regionSegs_nB[[segOneName]]$nB <- segOneB
              regionSegs_nB[[segOneName]]$chr <- segOneBchr 
              regionSegs_nB[[segOneName]]$pos <- segOneBpos 
              segsThatAreUsed <- c(segsThatAreUsed,segOneName)
              
            }
            
            if(tTestWithErrorHandling(segOneB,segTwoB)$p.val<=PvalToSeparate)
            {
              regionSegs_nB[[segOneName]]$nB <- segOneB 
              regionSegs_nB[[segOneName]]$chr <- segOneBchr 
              regionSegs_nB[[segOneName]]$pos <- segOneBpos 
              
              regionSegs_nB[[segTwoName]]$nB <- segTwoB 
              regionSegs_nB[[segTwoName]]$chr <- segTwoBchr 
              regionSegs_nB[[segTwoName]]$pos <- segTwoBpos 
              segsThatAreUsed <- c(segsThatAreUsed,segOneName,segTwoName)
              
            }
            
          }
          
        }
      }     
      
      # move to nA
      segOneA <- regionComparisonList[[region]][[segsToConsiderSubs[1]]]$nA
      segOneApos <- regionComparisonList[[regionOne]][[segsToConsiderSubs[1]]]$pos
      segOneAchr <- regionComparisonList[[regionOne]][[segsToConsiderSubs[1]]]$chr
      
      segTwoA <- regionComparisonList[[region]][[segsToConsiderSubs[2]]]$nA
      segTwoApos <- regionComparisonList[[regionOne]][[segsToConsiderSubs[2]]]$pos
      segTwoAchr <- regionComparisonList[[regionOne]][[segsToConsiderSubs[2]]]$chr
      
      segOneName <- segsToConsiderSubs[1]
      segTwoName <- segsToConsiderSubs[2]
      
      #segsThatAreUsedB <- c()
      if (length(segsToConsiderSubs)==2)
      {
        k <- 1
        if(tTestWithErrorHandling(segOneA,segTwoA)$p.val>PvalToSeparate)
        {
          segName <- paste(unlist(strsplit(segOneName,split=":"))[1],unlist(strsplit(segOneName,split=":"))[2],unlist(strsplit(segTwoName,split=":"))[3],sep=":")
          segB    <- c(segOneA,segTwoA)
          
          
          #regionSegs[[segName]]$nA <- segB  
          segOneApos <- c(segOneApos,segTwoApos)
          segOneAchr <- c(segOneAchr,segTwoAchr)
          segOneA    <- segB
          segOneName <- segName
          regionSegs_nA[[segOneName]]$nA <- segOneA
          regionSegs_nA[[segOneName]]$chr <- segOneAchr 
          regionSegs_nA[[segOneName]]$pos <- segOneApos 
          segsThatAreUsed <- c(segsThatAreUsed,segOneName)
          
        }
        
        if(tTestWithErrorHandling(segOneA,segTwoA)$p.val<=PvalToSeparate)
        {
          regionSegs_nA[[segOneName]]$nA <- segOneA 
          regionSegs_nA[[segOneName]]$chr <- segOneAchr 
          regionSegs_nA[[segOneName]]$pos <- segOneApos 
          
          regionSegs_nA[[segTwoName]]$nA <- segTwoA 
          regionSegs_nA[[segTwoName]]$chr <- segTwoAchr 
          regionSegs_nA[[segTwoName]]$pos <- segTwoApos 
          segsThatAreUsed <- c(segsThatAreUsed,segOneName,segTwoName)
          
        }
        
      }
      
      if (length(segsToConsiderSubs)==1)
      {
        #regionSegs[[segName]]$nB <- segB  
        #segOneBpos <- c(segOneBpos,segTwoApos)
        #segOneBchr <- c(segOneBchr,segTwoAchr)
        #segOneB    <- segB
        #segOneName <- segName
        regionSegs_nA[[segOneName]]$nA <- segOneA
        regionSegs_nA[[segOneName]]$chr <- segOneAchr 
        regionSegs_nA[[segOneName]]$pos <- segOneApos 
        segsThatAreUsed <- c(segsThatAreUsed,segOneName)
        
        
      }
      
      
      k <- 3
      while(k <= length(segsToConsiderSubs))
        #while(k < 14)
      {
        if(tTestWithErrorHandling(segOneA,segTwoA)$p.val>PvalToSeparate)
        {
          #combine the segments
          segName <- paste(unlist(strsplit(segOneName,split=":"))[1],unlist(strsplit(segOneName,split=":"))[2],unlist(strsplit(segTwoName,split=":"))[3],sep=":")
          segA    <- c(segOneA,segTwoA)
          
          
          #regionSegs[[segName]]$nA <- segA  
          segOneApos <- c(segOneApos,segTwoApos)
          segOneAchr <- c(segOneAchr,segTwoAchr)
          segOneA    <- segA
          segOneName <- segName
          
          
          segTwoA   <-  regionComparisonList[[region]][[segsToConsiderSubs[k]]]$nA
          segTwoAchr   <-  regionComparisonList[[regionOne]][[segsToConsiderSubs[k]]]$chr
          segTwoApos   <-  regionComparisonList[[regionOne]][[segsToConsiderSubs[k]]]$pos
          segTwoName <- segsToConsiderSubs[k]
          
          k <- k+1
          
          # check whether k exceeds number of segments
          if (k>length(segsToConsiderSubs))
          {
            if(tTestWithErrorHandling(segOneA,segTwoA)$p.val>PvalToSeparate)
            {
              segName <- paste(unlist(strsplit(segOneName,split=":"))[1],unlist(strsplit(segOneName,split=":"))[2],unlist(strsplit(segTwoName,split=":"))[3],sep=":")
              segA    <- c(segOneA,segTwoA)
              
              
              #regionSegs[[segName]]$nA <- segA  
              segOneApos <- c(segOneApos,segTwoApos)
              segOneAchr <- c(segOneAchr,segTwoAchr)
              segOneA    <- segA
              segOneName <- segName
              regionSegs_nA[[segOneName]]$nA <- segOneA
              regionSegs_nA[[segOneName]]$chr <- segOneAchr 
              regionSegs_nA[[segOneName]]$pos <- segOneApos 
              segsThatAreUsedB <- c(segsThatAreUsedB,segOneName)
              
            }
            
            if(tTestWithErrorHandling(segOneA,segTwoA)$p.val<=PvalToSeparate)
            {
              regionSegs_nA[[segOneName]]$nA <- segOneA 
              regionSegs_nA[[segOneName]]$chr <- segOneAchr 
              regionSegs_nA[[segOneName]]$pos <- segOneApos 
              
              regionSegs_nA[[segTwoName]]$nA <- segTwoA 
              regionSegs_nA[[segTwoName]]$chr <- segTwoAchr 
              regionSegs_nA[[segTwoName]]$pos <- segTwoApos 
              segsThatAreUsedB <- c(segsThatAreUsedB,segOneName,segTwoName)
              
            }
            
          }
        }
        
        if(tTestWithErrorHandling(segOneA,segTwoA)$p.val<=PvalToSeparate&k<=length(segsToConsiderSubs))
        {
          #print(k)
          # now add this to the regionSeg
          regionSegs_nA[[segOneName]]$nA <- segOneA 
          regionSegs_nA[[segOneName]]$chr <- segOneAchr 
          regionSegs_nA[[segOneName]]$pos <- segOneApos 
          segsThatAreUsedB <- c(segsThatAreUsedB,segOneName)
          
          
          segOneA  <- segTwoA
          segOneName <- segTwoName
          segOneApos <- segTwoApos
          segOneAchr <- segTwoAchr
          
          segTwoA  <- regionComparisonList[[region]][[segsToConsiderSubs[k]]]$nA
          segTwoAchr   <-  regionComparisonList[[regionOne]][[segsToConsiderSubs[k]]]$chr
          segTwoApos   <-  regionComparisonList[[regionOne]][[segsToConsiderSubs[k]]]$pos
          segTwoName  <- segsToConsiderSubs[k]
          
          k <- k+1
          if (k>length(segsToConsiderSubs))
          {
            if(tTestWithErrorHandling(segOneA,segTwoA)$p.val>PvalToSeparate)
            {
              segName <- paste(unlist(strsplit(segOneName,split=":"))[1],unlist(strsplit(segOneName,split=":"))[2],unlist(strsplit(segTwoName,split=":"))[3],sep=":")
              segA    <- c(segOneA,segTwoA)
              
              
              #regionSegs[[segName]]$nA <- segA  
              segOneApos <- c(segOneApos,segTwoApos)
              segOneAchr <- c(segOneAchr,segTwoAchr)
              segOneA    <- segA
              segOneName <- segName
              regionSegs_nA[[segOneName]]$nA <- segOneA
              regionSegs_nA[[segOneName]]$chr <- segOneAchr 
              regionSegs_nA[[segOneName]]$pos <- segOneApos 
              segsThatAreUsedB <- c(segsThatAreUsedB,segOneName)
              
            }
            
            if(tTestWithErrorHandling(segOneA,segTwoA)$p.val<=PvalToSeparate)
            {
              regionSegs_nA[[segOneName]]$nA <- segOneA 
              regionSegs_nA[[segOneName]]$chr <- segOneAchr 
              regionSegs_nA[[segOneName]]$pos <- segOneApos 
              
              regionSegs_nA[[segTwoName]]$nA <- segTwoA 
              regionSegs_nA[[segTwoName]]$chr <- segTwoAchr 
              regionSegs_nA[[segTwoName]]$pos <- segTwoApos 
              segsThatAreUsedB <- c(segsThatAreUsedB,segOneName,segTwoName)
              
            }
            
          }
        }
        
      }
      
      
    }
    
    newList[[region]]$segs_nA <- regionSegs_nA
    newList[[region]]$segs_nB <- regionSegs_nB
    
  }
  print("Returing from funJoinAdjacentSegments")
  return(newList)
  
}


funJoinAdjacentSegments_withCCF <- function(regionComparisonList
                                            ,CCFestimates
                                            ,PvalToSeparate=0.005)
{
  
  regionsToConsider <- grep(":",names(regionComparisonList),invert=TRUE,value=TRUE)
  segsToConsider    <- names(regionComparisonList[[1]])
  for (region in regionsToConsider)
  {
    segsToConsider <- intersect(segsToConsider,names(regionComparisonList[[region]]))
  }
  newWilcox.test <- function(x,y)
  {
    out <- wilcox.test(x,y)
    if(is.na(out$p.value))
    {
      out$p.value <- 1
    }
    return(out)
  }
  
  newList <- list()
  segsToUse <- c()
  for (seg in segsToConsider)
  {
    if(is.na(regionComparisonList[[1]][[seg]]$nB))
    {
      next;
    }
    segsToUse <- c(segsToUse,seg)
  }
  segsToConsider <- segsToUse
  
  
  regionOne <- regionsToConsider[1]
  for (region in regionsToConsider)
  {
    
    #first check segs for each chromosome
    segTable <- do.call(rbind,strsplit(segsToConsider,split=":"))
    segTable <- apply(segTable,2,as.numeric)
    rownames(segTable) <- segsToConsider
    #yourrighthere
    regionSegs_nA <- list()
    regionSegs_nB <- list()
    for (chr in unique(segTable[,1]))
    {
      segTableSub <- segTable[segTable[,1]%in%chr,,drop=FALSE]
      segsToConsiderSubs <- rownames(segTableSub)
      
      # start with nB
      {
        segOneB   <- regionComparisonList[[region]][[segsToConsiderSubs[1]]]$nB
        segOneBpos <- regionComparisonList[[regionOne]][[segsToConsiderSubs[1]]]$pos
        segOneBchr <- regionComparisonList[[regionOne]][[segsToConsiderSubs[1]]]$chr
        segOneB1frac <- CCFestimates[[segsToConsiderSubs[1]]][[segsToConsiderSubs[[1]]]][[region]]$fracB1_all
        segOneB2frac <- CCFestimates[[segsToConsiderSubs[1]]][[segsToConsiderSubs[[1]]]][[region]]$fracB2_all
        segOneBfrac <- CCFestimates[[segsToConsiderSubs[1]]][[segsToConsiderSubs[[1]]]][[region]]$fracB1_all
        segOneB1 <- rep(CCFestimates[[segsToConsiderSubs[1]]][[segsToConsiderSubs[[1]]]]$nB1,length(segOneB))
        segOneB2 <- rep(CCFestimates[[segsToConsiderSubs[1]]][[segsToConsiderSubs[[1]]]]$nB2,length(segOneB))
        
        segsThatAreUsed <- c()
        segsThatAreUsedB  <- c()
        if (length(segsToConsiderSubs)==1)
        {
          #regionSegs[[segName]]$nB <- segB  
          # segOneBpos <- c(segOneBpos,segTwoBpos)
          # segOneBchr <- c(segOneBchr,segTwoBchr)
          # segOneB    <- segB
          segOneName <- segsToConsiderSubs
          regionSegs_nB[[segOneName]]$nB <- segOneB
          regionSegs_nB[[segOneName]]$nBfrac <- segOneBfrac
          regionSegs_nB[[segOneName]]$nB1frac <- segOneB1frac
          regionSegs_nB[[segOneName]]$nB2frac <- segOneB2frac
          regionSegs_nB[[segOneName]]$nB1 <- segOneB1
          regionSegs_nB[[segOneName]]$nB2 <- segOneB2
          regionSegs_nB[[segOneName]]$chr <- segOneBchr 
          regionSegs_nB[[segOneName]]$pos <- segOneBpos 
          segsThatAreUsed <- c(segsThatAreUsed,segOneName)
          #next;
          
        }
        
        #note - add B1 and B2
        
        if(length(segsToConsiderSubs)>1)
        {
        
        segTwoB <- regionComparisonList[[region]][[segsToConsiderSubs[2]]]$nB
        segTwoBpos <- regionComparisonList[[regionOne]][[segsToConsiderSubs[2]]]$pos
        segTwoBchr <- regionComparisonList[[regionOne]][[segsToConsiderSubs[2]]]$chr
        segTwoBfrac <- CCFestimates[[segsToConsiderSubs[2]]][[segsToConsiderSubs[[2]]]][[region]]$fracB1_all
        segTwoB1frac <- CCFestimates[[segsToConsiderSubs[2]]][[segsToConsiderSubs[[2]]]][[region]]$fracB1_all
        segTwoB2frac <- CCFestimates[[segsToConsiderSubs[2]]][[segsToConsiderSubs[[2]]]][[region]]$fracB2_all
        segTwoB1 <- rep(CCFestimates[[segsToConsiderSubs[2]]][[segsToConsiderSubs[[2]]]]$nB1,length(segTwoB))
        segTwoB2 <- rep(CCFestimates[[segsToConsiderSubs[2]]][[segsToConsiderSubs[[2]]]]$nB2,length(segTwoB))
        
        
        segOneName <- segsToConsiderSubs[1]
        segTwoName <- segsToConsiderSubs[2]
        }
        
        #check if less than 3 segments
        if (length(segsToConsiderSubs)==2)
        {
          k <- 1
          if(newWilcox.test(segOneB,segTwoB)$p.val>PvalToSeparate&newWilcox.test(segOneBfrac,segTwoBfrac)$p.val>PvalToSeparate)
          {
            segName <- paste(unlist(strsplit(segOneName,split=":"))[1],unlist(strsplit(segOneName,split=":"))[2],unlist(strsplit(segTwoName,split=":"))[3],sep=":")
            segB    <- c(segOneB,segTwoB)
            segBfrac <- c(segOneBfrac,segTwoBfrac)
            segB1frac <- c(segOneB1frac,segTwoB1frac)
            segB2frac <- c(segOneB1frac,segTwoB2frac)
            segB1    <- c(segOneB1,segTwoB1)
            segB2    <- c(segOneB2,segTwoB2)
            
            #regionSegs[[segName]]$nB <- segB  
            segOneBpos <- c(segOneBpos,segTwoBpos)
            segOneBchr <- c(segOneBchr,segTwoBchr)
            segOneB    <- segB
            segOneBfrac    <- segBfrac
            segOneB1frac    <- segB1frac
            segOneB2frac    <- segB2frac
            segOneB1    <- segB1
            segOneB2    <- segB2
            segOneName <- segName
            regionSegs_nB[[segOneName]]$nB <- segOneB
            regionSegs_nB[[segOneName]]$nBfrac <- segOneBfrac
            regionSegs_nB[[segOneName]]$nB1frac <- segOneB1frac
            regionSegs_nB[[segOneName]]$nB2frac <- segOneB2frac
            regionSegs_nB[[segOneName]]$nB1 <- segOneB1
            regionSegs_nB[[segOneName]]$nB2 <- segOneB2
            regionSegs_nB[[segOneName]]$chr <- segOneBchr 
            regionSegs_nB[[segOneName]]$pos <- segOneBpos 
            segsThatAreUsed <- c(segsThatAreUsed,segOneName)
            
          }
          
          if(newWilcox.test(segOneB,segTwoB)$p.val<=PvalToSeparate|newWilcox.test(segOneBfrac,segTwoBfrac)$p.val<=PvalToSeparate)
          {
            regionSegs_nB[[segOneName]]$nB <- segOneB
            regionSegs_nB[[segOneName]]$nBfrac <- segOneBfrac
            regionSegs_nB[[segOneName]]$nB1frac <- segOneB1frac
            regionSegs_nB[[segOneName]]$nB2frac <- segOneB2frac
            regionSegs_nB[[segOneName]]$nB1 <- segOneB1
            regionSegs_nB[[segOneName]]$nB2 <- segOneB2
            
            regionSegs_nB[[segOneName]]$chr <- segOneBchr 
            regionSegs_nB[[segOneName]]$pos <- segOneBpos 
            
            regionSegs_nB[[segTwoName]]$nB <- segTwoB 
            regionSegs_nB[[segTwoName]]$nBfrac <- segTwoBfrac
            regionSegs_nB[[segTwoName]]$nB1frac <- segTwoB1frac
            regionSegs_nB[[segTwoName]]$nB2frac <- segTwoB2frac
            regionSegs_nB[[segTwoName]]$nB1 <- segTwoB1 
            regionSegs_nB[[segTwoName]]$nB2 <- segTwoB2 
            
            regionSegs_nB[[segTwoName]]$chr <- segTwoBchr 
            regionSegs_nB[[segTwoName]]$pos <- segTwoBpos 
            segsThatAreUsed <- c(segsThatAreUsed,segOneName,segTwoName)
            
          }
          
        }
        
        
        
        k <- 3
        while(k <= length(segsToConsiderSubs))
        {
          #print(segsToConsiderSubs[k])
          if(newWilcox.test(segOneB,segTwoB)$p.val>PvalToSeparate&newWilcox.test(segOneBfrac,segTwoBfrac)$p.val>PvalToSeparate)
          {
            #seg one and two can be joined 
            
            #combine the segments
            segName <- paste(unlist(strsplit(segOneName,split=":"))[1],unlist(strsplit(segOneName,split=":"))[2],unlist(strsplit(segTwoName,split=":"))[3],sep=":")
            segB    <- c(segOneB,segTwoB)
            segBfrac <- c(segOneBfrac,segTwoBfrac)
            segB1frac <- c(segOneBfrac,segTwoB1frac)
            segB2frac <- c(segOneBfrac,segTwoB2frac)
            segB1 <- c(segOneB1,segTwoB1)
            segB2 <- c(segOneB2,segTwoB2)
            
            #regionSegs[[segName]]$nB <- segB  
            segOneBpos <- c(segOneBpos,segTwoBpos)
            segOneBchr <- c(segOneBchr,segTwoBchr)
            segOneB    <- segB
            segOneBfrac    <- segBfrac
            segOne1Bfrac    <- segB1frac
            segOne2Bfrac    <- segB2frac
            segOneB1    <- segB1
            segOneB2    <- segB2
            segOneName <- segName
            
            
            segTwoB      <-  regionComparisonList[[region]][[segsToConsiderSubs[k]]]$nB
            segTwoBchr   <-  regionComparisonList[[regionOne]][[segsToConsiderSubs[k]]]$chr
            segTwoBpos   <-  regionComparisonList[[regionOne]][[segsToConsiderSubs[k]]]$pos
            segTwoBfrac  <- CCFestimates[[segsToConsiderSubs[k]]][[segsToConsiderSubs[[k]]]][[region]]$fracB1_all
            segTwoB1frac  <- CCFestimates[[segsToConsiderSubs[k]]][[segsToConsiderSubs[[k]]]][[region]]$fracB1_all
            segTwoB2frac  <- CCFestimates[[segsToConsiderSubs[k]]][[segsToConsiderSubs[[k]]]][[region]]$fracB2_all
            segTwoB1  <- rep(CCFestimates[[segsToConsiderSubs[k]]][[segsToConsiderSubs[[k]]]]$nB1,length(segTwoB))
            segTwoB2  <- rep(CCFestimates[[segsToConsiderSubs[k]]][[segsToConsiderSubs[[k]]]]$nB2,length(segTwoB))
            
            segTwoName   <- segsToConsiderSubs[k]
            
            k <- k+1
            
            # check whether k exceeds number of segments
            if (k>length(segsToConsiderSubs))
            {
              if(newWilcox.test(segOneB,segTwoB)$p.val>PvalToSeparate&newWilcox.test(segOneBfrac,segTwoBfrac)$p.val>PvalToSeparate)
              {
                segName <- paste(unlist(strsplit(segOneName,split=":"))[1],unlist(strsplit(segOneName,split=":"))[2],unlist(strsplit(segTwoName,split=":"))[3],sep=":")
                segB    <- c(segOneB,segTwoB)
                segBfrac <- c(segOneBfrac,segTwoBfrac)
                segB1frac <- c(segOneBfrac,segTwoB1frac)
                segB2frac <- c(segOneBfrac,segTwoB2frac)
                segB1 <- c(segOneB1,segTwoB1)
                segB2 <- c(segOneB2,segTwoB2)
                
                #regionSegs[[segName]]$nB <- segB  
                segOneBpos <- c(segOneBpos,segTwoBpos)
                segOneBchr <- c(segOneBchr,segTwoBchr)
                segOneB    <- segB
                segOneBfrac <- segBfrac
                segOneB1frac <- segB1frac
                segOneB2frac <- segB2frac
                segOneB1 <- segB1
                segOneB2 <- segB2
                
                
                
                segOneName <- segName
                regionSegs_nB[[segOneName]]$nB <- segOneB
                regionSegs_nB[[segOneName]]$nBfrac <- segOneBfrac
                regionSegs_nB[[segOneName]]$nB1frac <- segOneB1frac
                regionSegs_nB[[segOneName]]$nB2frac <- segOneB2frac
                regionSegs_nB[[segOneName]]$nB1 <- segOneB1
                regionSegs_nB[[segOneName]]$nB2 <- segOneB2
                regionSegs_nB[[segOneName]]$chr <- segOneBchr 
                regionSegs_nB[[segOneName]]$pos <- segOneBpos 
                segsThatAreUsed <- c(segsThatAreUsed,segOneName)
                
              }
              
              if(newWilcox.test(segOneB,segTwoB)$p.val<=PvalToSeparate|newWilcox.test(segOneBfrac,segTwoBfrac)$p.val<=PvalToSeparate)
              {
                regionSegs_nB[[segOneName]]$nB <- segOneB 
                regionSegs_nB[[segOneName]]$nBfrac <- segOneBfrac
                regionSegs_nB[[segOneName]]$nB1frac <- segOneB1frac
                regionSegs_nB[[segOneName]]$nB2frac <- segOneB2frac
                regionSegs_nB[[segOneName]]$nB1 <- segOneB1
                regionSegs_nB[[segOneName]]$nB2 <- segOneB2
                regionSegs_nB[[segOneName]]$chr <- segOneBchr 
                regionSegs_nB[[segOneName]]$pos <- segOneBpos 
                
                regionSegs_nB[[segTwoName]]$nB <- segTwoB
                regionSegs_nB[[segTwoName]]$nBfrac <- segTwoBfrac
                regionSegs_nB[[segTwoName]]$nB1frac <- segTwoB1frac
                regionSegs_nB[[segTwoName]]$nB2frac <- segTwoB2frac
                regionSegs_nB[[segTwoName]]$nB1 <- segTwoB1
                regionSegs_nB[[segTwoName]]$nB2 <- segTwoB2
                regionSegs_nB[[segTwoName]]$chr <- segTwoBchr 
                regionSegs_nB[[segTwoName]]$pos <- segTwoBpos 
                segsThatAreUsed <- c(segsThatAreUsed,segOneName,segTwoName)
                
              }
              
            }
            
          }
          
          if((newWilcox.test(segOneB,segTwoB)$p.val<=PvalToSeparate|newWilcox.test(segOneBfrac,segTwoBfrac)$p.val<=PvalToSeparate)&k<=length(segsToConsiderSubs))
          {
            segsThatAreUsed <- c(segsThatAreUsed,segOneName)
            #print(k)
            # now add this to the regionSeg
            regionSegs_nB[[segOneName]]$nB <- segOneB 
            regionSegs_nB[[segOneName]]$nBfrac <- segOneBfrac 
            regionSegs_nB[[segOneName]]$nB1frac <- segOneB1frac 
            regionSegs_nB[[segOneName]]$nB2frac <- segOneB2frac 
            regionSegs_nB[[segOneName]]$nB1 <- segOneB1
            regionSegs_nB[[segOneName]]$nB2 <- segOneB2
            regionSegs_nB[[segOneName]]$chr <- segOneBchr 
            regionSegs_nB[[segOneName]]$pos <- segOneBpos 
            
            
            segOneB  <- segTwoB
            segOneBfrac  <- segTwoBfrac
            segOneB1frac  <- segTwoB1frac
            segOneB2frac  <- segTwoB2frac
            segOneB1  <- segTwoB1
            segOneB2  <- segTwoB2
            segOneName <- segTwoName
            segOneBpos <- segTwoBpos
            segOneBchr <- segTwoBchr
            
            segTwoB  <- regionComparisonList[[region]][[segsToConsiderSubs[k]]]$nB
            segTwoBfrac  <-  CCFestimates[[segsToConsiderSubs[k]]][[segsToConsiderSubs[[k]]]][[region]]$fracB1_all#regionComparisonList[[region]][[segsToConsiderSubs[k]]]$nBfrac
            segTwoB1frac  <-  CCFestimates[[segsToConsiderSubs[k]]][[segsToConsiderSubs[[k]]]][[region]]$fracB1_all#regionComparisonList[[region]][[segsToConsiderSubs[k]]]$nBfrac
            segTwoB2frac  <-  CCFestimates[[segsToConsiderSubs[k]]][[segsToConsiderSubs[[k]]]][[region]]$fracB2_all#regionComparisonList[[region]][[segsToConsiderSubs[k]]]$nBfrac
            segTwoB1  <-  rep(CCFestimates[[segsToConsiderSubs[k]]][[segsToConsiderSubs[[k]]]]$nB1,length(segTwoB))
            segTwoB2  <-  rep(CCFestimates[[segsToConsiderSubs[k]]][[segsToConsiderSubs[[k]]]]$nB2,length(segTwoB))#regionComparisonList[[region]][[segsToConsiderSubs[k]]]$nBfrac
            
            segTwoBchr   <-  regionComparisonList[[regionOne]][[segsToConsiderSubs[k]]]$chr
            segTwoBpos   <-  regionComparisonList[[regionOne]][[segsToConsiderSubs[k]]]$pos
            segTwoName  <- segsToConsiderSubs[k]
            
            k <- k+1
            if (k>length(segsToConsiderSubs))
            {
              if(newWilcox.test(segOneB,segTwoB)$p.val>PvalToSeparate&newWilcox.test(segOneBfrac,segTwoBfrac)$p.val>PvalToSeparate)
              {
                segName <- paste(unlist(strsplit(segOneName,split=":"))[1],unlist(strsplit(segOneName,split=":"))[2],unlist(strsplit(segTwoName,split=":"))[3],sep=":")
                segB    <- c(segOneB,segTwoB)
                segBfrac    <- c(segOneBfrac,segTwoBfrac)
                segB1frac    <- c(segOneBfrac,segTwoB1frac)
                segB2frac    <- c(segOneBfrac,segTwoB2frac)
                segB1    <- c(segOneB1,segTwoB1)
                segB2    <- c(segOneB2,segTwoB2)
                
                
                #regionSegs[[segName]]$nB <- segB  
                segOneBpos <- c(segOneBpos,segTwoBpos)
                segOneBchr <- c(segOneBchr,segTwoBchr)
                segOneB    <- segB
                segOneBfrac    <- segBfrac
                segOneB1frac    <- segB1frac
                segOneB2frac    <- segB2frac
                segOneB1    <- segB1
                segOneB2    <- segB2
                segOneName <- segName
                regionSegs_nB[[segOneName]]$nB <- segOneB
                regionSegs_nB[[segOneName]]$nBfrac <- segOneBfrac
                regionSegs_nB[[segOneName]]$nB1frac <- segOneB1frac
                regionSegs_nB[[segOneName]]$nB2frac <- segOneB2frac
                regionSegs_nB[[segOneName]]$nB1 <- segOneB1
                regionSegs_nB[[segOneName]]$nB2 <- segOneB2
                regionSegs_nB[[segOneName]]$chr <- segOneBchr 
                regionSegs_nB[[segOneName]]$pos <- segOneBpos 
                segsThatAreUsed <- c(segsThatAreUsed,segOneName)
                
              }
              
              if(newWilcox.test(segOneB,segTwoB)$p.val<=PvalToSeparate|newWilcox.test(segOneBfrac,segTwoBfrac)$p.val<=PvalToSeparate)
              {
                regionSegs_nB[[segOneName]]$nB <- segOneB 
                regionSegs_nB[[segOneName]]$nBfrac <- segOneBfrac
                regionSegs_nB[[segOneName]]$nB1frac <- segOneB1frac
                regionSegs_nB[[segOneName]]$nB2frac <- segOneB2frac
                regionSegs_nB[[segOneName]]$nB1 <- segOneB1
                regionSegs_nB[[segOneName]]$nB2 <- segOneB2
                regionSegs_nB[[segOneName]]$chr <- segOneBchr 
                regionSegs_nB[[segOneName]]$pos <- segOneBpos 
                
                regionSegs_nB[[segTwoName]]$nB <- segTwoB 
                regionSegs_nB[[segTwoName]]$nBfrac <- segTwoBfrac 
                regionSegs_nB[[segTwoName]]$nB1frac <- segTwoB1frac 
                regionSegs_nB[[segTwoName]]$nB2frac <- segTwoB2frac 
                regionSegs_nB[[segTwoName]]$nB1 <- segTwoB1
                regionSegs_nB[[segTwoName]]$nB2 <- segTwoB2
                regionSegs_nB[[segTwoName]]$chr <- segTwoBchr 
                regionSegs_nB[[segTwoName]]$pos <- segTwoBpos 
                segsThatAreUsed <- c(segsThatAreUsed,segOneName,segTwoName)
                
              }
              
            }
            
          }
        } 
      }
      
      # move to nA
      {
        segOneA   <- regionComparisonList[[region]][[segsToConsiderSubs[1]]]$nA
        segOneApos <- regionComparisonList[[regionOne]][[segsToConsiderSubs[1]]]$pos
        segOneAchr <- regionComparisonList[[regionOne]][[segsToConsiderSubs[1]]]$chr
        segOneAfrac <- CCFestimates[[segsToConsiderSubs[1]]][[segsToConsiderSubs[[1]]]][[region]]$fracA1_all
        segOneA1frac <- CCFestimates[[segsToConsiderSubs[1]]][[segsToConsiderSubs[[1]]]][[region]]$fracA1_all
        segOneA2frac <- CCFestimates[[segsToConsiderSubs[1]]][[segsToConsiderSubs[[1]]]][[region]]$fracA2_all
        segOneA1 <- rep(CCFestimates[[segsToConsiderSubs[1]]][[segsToConsiderSubs[[1]]]]$nA1,length(segOneA))
        segOneA2 <- rep(CCFestimates[[segsToConsiderSubs[1]]][[segsToConsiderSubs[[1]]]]$nA2,length(segOneA))
        
        segsThatAreUsed <- c()
        segsThatAreUsedA  <- c()
        
        if (length(segsToConsiderSubs)==1)
        {
          #regionSegs[[segName]]$nA <- segA  
          # segOneApos <- c(segOneApos,segTwoApos)
          # segOneAchr <- c(segOneAchr,segTwoAchr)
          # segOneA    <- segA
          segOneName <- segsToConsiderSubs
          regionSegs_nA[[segOneName]]$nA <- segOneA
          regionSegs_nA[[segOneName]]$nAfrac <- segOneAfrac
          regionSegs_nA[[segOneName]]$nA1frac <- segOneA1frac
          regionSegs_nA[[segOneName]]$nA2frac <- segOneA2frac
          regionSegs_nA[[segOneName]]$nA1 <- segOneA1
          regionSegs_nA[[segOneName]]$nA2 <- segOneA2
          regionSegs_nA[[segOneName]]$chr <- segOneAchr 
          regionSegs_nA[[segOneName]]$pos <- segOneApos 
          segsThatAreUsed <- c(segsThatAreUsed,segOneName)
          #next;
          
        }
        
        if(length(segsToConsiderSubs)>1)
        {
        segTwoA <- regionComparisonList[[region]][[segsToConsiderSubs[2]]]$nA
        segTwoApos <- regionComparisonList[[regionOne]][[segsToConsiderSubs[2]]]$pos
        segTwoAchr <- regionComparisonList[[regionOne]][[segsToConsiderSubs[2]]]$chr
        segTwoAfrac <- CCFestimates[[segsToConsiderSubs[2]]][[segsToConsiderSubs[[2]]]][[region]]$fracA1_all
        segTwoA1frac <- CCFestimates[[segsToConsiderSubs[2]]][[segsToConsiderSubs[[2]]]][[region]]$fracA1_all
        segTwoA2frac <- CCFestimates[[segsToConsiderSubs[2]]][[segsToConsiderSubs[[2]]]][[region]]$fracA2_all
        segTwoA1 <- rep(CCFestimates[[segsToConsiderSubs[2]]][[segsToConsiderSubs[[2]]]]$nA1,length(segTwoA))
        segTwoA2 <- rep(CCFestimates[[segsToConsiderSubs[2]]][[segsToConsiderSubs[[2]]]]$nA2,length(segTwoA))
        
        
        segOneName <- segsToConsiderSubs[1]
        segTwoName <- segsToConsiderSubs[2]
        }
        #check if less than 3 segments
        if (length(segsToConsiderSubs)==2)
        {
          k <- 1
          if(newWilcox.test(segOneA,segTwoA)$p.val>PvalToSeparate&newWilcox.test(segOneAfrac,segTwoAfrac)$p.val>PvalToSeparate)
          {
            segName <- paste(unlist(strsplit(segOneName,split=":"))[1],unlist(strsplit(segOneName,split=":"))[2],unlist(strsplit(segTwoName,split=":"))[3],sep=":")
            segA    <- c(segOneA,segTwoA)
            segAfrac <- c(segOneAfrac,segTwoAfrac)
            segA1frac <- c(segOneA1frac,segTwoA1frac)
            segA2frac <- c(segOneA2frac,segTwoA2frac)
            segA1    <- c(segOneA1,segTwoA1)
            segA2    <- c(segOneA2,segTwoA2)
            
            #regionSegs[[segName]]$nA <- segA  
            segOneApos <- c(segOneApos,segTwoApos)
            segOneAchr <- c(segOneAchr,segTwoAchr)
            segOneA    <- segA
            segOneAfrac    <- segAfrac
            segOneA1frac    <- segA1frac
            segOneA2frac    <- segA2frac
            segOneA1    <- segA1
            segOneA2    <- segA2
            segOneName <- segName
            regionSegs_nA[[segOneName]]$nA <- segOneA
            regionSegs_nA[[segOneName]]$nAfrac <- segOneAfrac
            regionSegs_nA[[segOneName]]$nA1frac <- segOneA1frac
            regionSegs_nA[[segOneName]]$nA2frac <- segOneA2frac
            regionSegs_nA[[segOneName]]$nA1 <- segOneA1
            regionSegs_nA[[segOneName]]$nA2 <- segOneA2
            regionSegs_nA[[segOneName]]$chr <- segOneAchr 
            regionSegs_nA[[segOneName]]$pos <- segOneApos 
            segsThatAreUsed <- c(segsThatAreUsed,segOneName)
            
          }
          
          if(newWilcox.test(segOneA,segTwoA)$p.val<=PvalToSeparate|newWilcox.test(segOneAfrac,segTwoAfrac)$p.val<=PvalToSeparate)
          {
            regionSegs_nA[[segOneName]]$nA <- segOneA
            regionSegs_nA[[segOneName]]$nAfrac <- segOneAfrac
            regionSegs_nA[[segOneName]]$nA1frac <- segOneA1frac
            regionSegs_nA[[segOneName]]$nA2frac <- segOneA2frac
            regionSegs_nA[[segOneName]]$nA1 <- segOneA1
            regionSegs_nA[[segOneName]]$nA2 <- segOneA2
            
            regionSegs_nA[[segOneName]]$chr <- segOneAchr 
            regionSegs_nA[[segOneName]]$pos <- segOneApos 
            
            regionSegs_nA[[segTwoName]]$nA <- segTwoA 
            regionSegs_nA[[segTwoName]]$nAfrac <- segTwoAfrac
            regionSegs_nA[[segTwoName]]$nA1frac <- segTwoA1frac
            regionSegs_nA[[segTwoName]]$nA2frac <- segTwoA2frac
            regionSegs_nA[[segTwoName]]$nA1 <- segTwoA1 
            regionSegs_nA[[segTwoName]]$nA2 <- segTwoA2 
            
            regionSegs_nA[[segTwoName]]$chr <- segTwoAchr 
            regionSegs_nA[[segTwoName]]$pos <- segTwoApos 
            segsThatAreUsed <- c(segsThatAreUsed,segOneName,segTwoName)
            
          }
          
        }
        
        
        
        k <- 3
        while(k <= length(segsToConsiderSubs))
        {
          #print(segsToConsiderSubs[k])
          if(newWilcox.test(segOneA,segTwoA)$p.val>PvalToSeparate&newWilcox.test(segOneAfrac,segTwoAfrac)$p.val>PvalToSeparate)
          {
            #seg one and two can be joined 
            
            #combine the segments
            segName <- paste(unlist(strsplit(segOneName,split=":"))[1],unlist(strsplit(segOneName,split=":"))[2],unlist(strsplit(segTwoName,split=":"))[3],sep=":")
            segA    <- c(segOneA,segTwoA)
            segAfrac <- c(segOneAfrac,segTwoAfrac)
            segA1frac <- c(segOneA1frac,segTwoA1frac)
            segA2frac <- c(segOneA2frac,segTwoA2frac)
            segA1 <- c(segOneA1,segTwoA1)
            segA2 <- c(segOneA2,segTwoA2)
            
            #regionSegs[[segName]]$nA <- segA  
            segOneApos <- c(segOneApos,segTwoApos)
            segOneAchr <- c(segOneAchr,segTwoAchr)
            segOneA    <- segA
            segOneAfrac    <- segAfrac
            segOneA1frac    <- segA1frac
            segOneA2frac    <- segA2frac
            segOneA1    <- segA1
            segOneA2    <- segA2
            segOneName <- segName
            
            
            segTwoA      <-  regionComparisonList[[region]][[segsToConsiderSubs[k]]]$nA
            segTwoAchr   <-  regionComparisonList[[regionOne]][[segsToConsiderSubs[k]]]$chr
            segTwoApos   <-  regionComparisonList[[regionOne]][[segsToConsiderSubs[k]]]$pos
            segTwoAfrac  <- CCFestimates[[segsToConsiderSubs[k]]][[segsToConsiderSubs[[k]]]][[region]]$fracA1_all
            segTwoA1frac  <- CCFestimates[[segsToConsiderSubs[k]]][[segsToConsiderSubs[[k]]]][[region]]$fracA1_all
            segTwoA2frac  <- CCFestimates[[segsToConsiderSubs[k]]][[segsToConsiderSubs[[k]]]][[region]]$fracA2_all
            segTwoA1  <- rep(CCFestimates[[segsToConsiderSubs[k]]][[segsToConsiderSubs[[k]]]]$nA1,length(segTwoA))
            segTwoA2  <- rep(CCFestimates[[segsToConsiderSubs[k]]][[segsToConsiderSubs[[k]]]]$nA2,length(segTwoA))
            
            segTwoName   <- segsToConsiderSubs[k]
            
            k <- k+1
            
            # check whether k exceeds number of segments
            if (k>length(segsToConsiderSubs))
            {
              if(newWilcox.test(segOneA,segTwoA)$p.val>PvalToSeparate&newWilcox.test(segOneAfrac,segTwoAfrac)$p.val>PvalToSeparate)
              {
                segName <- paste(unlist(strsplit(segOneName,split=":"))[1],unlist(strsplit(segOneName,split=":"))[2],unlist(strsplit(segTwoName,split=":"))[3],sep=":")
                segA    <- c(segOneA,segTwoA)
                segAfrac <- c(segOneAfrac,segTwoAfrac)
                segA1frac <- c(segOneA1frac,segTwoA1frac)
                segA2frac <- c(segOneA2frac,segTwoA2frac)
                segA1 <- c(segOneA1,segTwoA1)
                segA2 <- c(segOneA2,segTwoA2)
                
                #regionSegs[[segName]]$nA <- segA  
                segOneApos <- c(segOneApos,segTwoApos)
                segOneAchr <- c(segOneAchr,segTwoAchr)
                segOneA    <- segA
                segOneAfrac <- segAfrac
                segOneA1frac <- segA1frac
                segOneA2frac <- segA2frac
                segOneA1 <- segA1
                segOneA2 <- segA2
                
                
                
                segOneName <- segName
                regionSegs_nA[[segOneName]]$nA <- segOneA
                regionSegs_nA[[segOneName]]$nAfrac <- segOneAfrac
                regionSegs_nA[[segOneName]]$nA1frac <- segOneA1frac
                regionSegs_nA[[segOneName]]$nA2frac <- segOneA2frac
                regionSegs_nA[[segOneName]]$nA1 <- segOneA1
                regionSegs_nA[[segOneName]]$nA2 <- segOneA2
                regionSegs_nA[[segOneName]]$chr <- segOneAchr 
                regionSegs_nA[[segOneName]]$pos <- segOneApos 
                segsThatAreUsed <- c(segsThatAreUsed,segOneName)
                
              }
              
              if(newWilcox.test(segOneA,segTwoA)$p.val<=PvalToSeparate|newWilcox.test(segOneAfrac,segTwoAfrac)$p.val<=PvalToSeparate)
              {
                regionSegs_nA[[segOneName]]$nA <- segOneA 
                regionSegs_nA[[segOneName]]$nAfrac <- segOneAfrac
                regionSegs_nA[[segOneName]]$nA1frac <- segOneA1frac
                regionSegs_nA[[segOneName]]$nA2frac <- segOneA2frac
                regionSegs_nA[[segOneName]]$nA1 <- segOneA1
                regionSegs_nA[[segOneName]]$nA2 <- segOneA2
                regionSegs_nA[[segOneName]]$chr <- segOneAchr 
                regionSegs_nA[[segOneName]]$pos <- segOneApos 
                
                regionSegs_nA[[segTwoName]]$nA <- segTwoA
                regionSegs_nA[[segTwoName]]$nAfrac <- segTwoAfrac
                regionSegs_nA[[segTwoName]]$nA1frac <- segTwoA1frac
                regionSegs_nA[[segTwoName]]$nA2frac <- segTwoA2frac
                regionSegs_nA[[segTwoName]]$nA1 <- segTwoA1
                regionSegs_nA[[segTwoName]]$nA2 <- segTwoA2
                regionSegs_nA[[segTwoName]]$chr <- segTwoAchr 
                regionSegs_nA[[segTwoName]]$pos <- segTwoApos 
                segsThatAreUsed <- c(segsThatAreUsed,segOneName,segTwoName)
                
              }
              
            }
            
          }
          
          if((newWilcox.test(segOneA,segTwoA)$p.val<=PvalToSeparate|newWilcox.test(segOneAfrac,segTwoAfrac)$p.val<=PvalToSeparate)&k<=length(segsToConsiderSubs))
          {
            segsThatAreUsed <- c(segsThatAreUsed,segOneName)
            #print(k)
            # now add this to the regionSeg
            regionSegs_nA[[segOneName]]$nA <- segOneA 
            regionSegs_nA[[segOneName]]$nAfrac <- segOneAfrac 
            regionSegs_nA[[segOneName]]$nA1frac <- segOneA1frac 
            regionSegs_nA[[segOneName]]$nA2frac <- segOneA2frac 
            regionSegs_nA[[segOneName]]$nA1 <- segOneA1
            regionSegs_nA[[segOneName]]$nA2 <- segOneA2
            regionSegs_nA[[segOneName]]$chr <- segOneAchr 
            regionSegs_nA[[segOneName]]$pos <- segOneApos 
            
            
            segOneA  <- segTwoA
            segOneAfrac  <- segTwoAfrac
            segOneA1frac  <- segTwoA1frac
            segOneA2frac  <- segTwoA2frac
            segOneA1  <- segTwoA1
            segOneA2  <- segTwoA2
            segOneName <- segTwoName
            segOneApos <- segTwoApos
            segOneAchr <- segTwoAchr
            
            segTwoA  <- regionComparisonList[[region]][[segsToConsiderSubs[k]]]$nA
            segTwoAfrac  <-  CCFestimates[[segsToConsiderSubs[k]]][[segsToConsiderSubs[[k]]]][[region]]$fracA1_all#regionComparisonList[[region]][[segsToConsiderSubs[k]]]$nAfrac
            segTwoA1frac  <-  CCFestimates[[segsToConsiderSubs[k]]][[segsToConsiderSubs[[k]]]][[region]]$fracA1_all#regionComparisonList[[region]][[segsToConsiderSubs[k]]]$nAfrac
            segTwoA2frac  <-  CCFestimates[[segsToConsiderSubs[k]]][[segsToConsiderSubs[[k]]]][[region]]$fracA2_all#regionComparisonList[[region]][[segsToConsiderSubs[k]]]$nAfrac
            segTwoA1  <-  rep(CCFestimates[[segsToConsiderSubs[k]]][[segsToConsiderSubs[[k]]]]$nA1,length(segTwoA))
            segTwoA2  <-  rep(CCFestimates[[segsToConsiderSubs[k]]][[segsToConsiderSubs[[k]]]]$nA2,length(segTwoA))#regionComparisonList[[region]][[segsToConsiderSubs[k]]]$nAfrac
            
            segTwoAchr   <-  regionComparisonList[[regionOne]][[segsToConsiderSubs[k]]]$chr
            segTwoApos   <-  regionComparisonList[[regionOne]][[segsToConsiderSubs[k]]]$pos
            segTwoName  <- segsToConsiderSubs[k]
            
            k <- k+1
            if (k>length(segsToConsiderSubs))
            {
              if(newWilcox.test(segOneA,segTwoA)$p.val>PvalToSeparate&newWilcox.test(segOneAfrac,segTwoAfrac)$p.val>PvalToSeparate)
              {
                segName <- paste(unlist(strsplit(segOneName,split=":"))[1],unlist(strsplit(segOneName,split=":"))[2],unlist(strsplit(segTwoName,split=":"))[3],sep=":")
                segA    <- c(segOneA,segTwoA)
                segAfrac    <- c(segOneAfrac,segTwoAfrac)
                segA1frac    <- c(segOneA1frac,segTwoA1frac)
                segA2frac    <- c(segOneA2frac,segTwoA2frac)
                segA1    <- c(segOneA1,segTwoA1)
                segA2    <- c(segOneA2,segTwoA2)
                
                
                #regionSegs[[segName]]$nA <- segA  
                segOneApos <- c(segOneApos,segTwoApos)
                segOneAchr <- c(segOneAchr,segTwoAchr)
                segOneA    <- segA
                segOneAfrac    <- segAfrac
                segOneA1frac    <- segA1frac
                segOneA2frac    <- segA2frac
                segOneA1    <- segA1
                segOneA2    <- segA2
                segOneName <- segName
                regionSegs_nA[[segOneName]]$nA <- segOneA
                regionSegs_nA[[segOneName]]$nAfrac <- segOneAfrac
                regionSegs_nA[[segOneName]]$nA1frac <- segOneA1frac
                regionSegs_nA[[segOneName]]$nA2frac <- segOneA2frac
                regionSegs_nA[[segOneName]]$nA1 <- segOneA1
                regionSegs_nA[[segOneName]]$nA2 <- segOneA2
                regionSegs_nA[[segOneName]]$chr <- segOneAchr 
                regionSegs_nA[[segOneName]]$pos <- segOneApos 
                segsThatAreUsed <- c(segsThatAreUsed,segOneName)
                
              }
              
              if(newWilcox.test(segOneA,segTwoA)$p.val<=PvalToSeparate|newWilcox.test(segOneAfrac,segTwoAfrac)$p.val<=PvalToSeparate)
              {
                regionSegs_nA[[segOneName]]$nA <- segOneA 
                regionSegs_nA[[segOneName]]$nAfrac <- segOneAfrac
                regionSegs_nA[[segOneName]]$nA1frac <- segOneA1frac
                regionSegs_nA[[segOneName]]$nA2frac <- segOneA2frac
                regionSegs_nA[[segOneName]]$nA1 <- segOneA1
                regionSegs_nA[[segOneName]]$nA2 <- segOneA2
                regionSegs_nA[[segOneName]]$chr <- segOneAchr 
                regionSegs_nA[[segOneName]]$pos <- segOneApos 
                
                regionSegs_nA[[segTwoName]]$nA <- segTwoA 
                regionSegs_nA[[segTwoName]]$nAfrac <- segTwoAfrac 
                regionSegs_nA[[segTwoName]]$nA1frac <- segTwoA1frac 
                regionSegs_nA[[segTwoName]]$nA2frac <- segTwoA2frac 
                regionSegs_nA[[segTwoName]]$nA1 <- segTwoA1
                regionSegs_nA[[segTwoName]]$nA2 <- segTwoA2
                regionSegs_nA[[segTwoName]]$chr <- segTwoAchr 
                regionSegs_nA[[segTwoName]]$pos <- segTwoApos 
                segsThatAreUsed <- c(segsThatAreUsed,segOneName,segTwoName)
                
              }
              
            }
            
          }
        } 
      }
    }
    
    newList[[region]]$segs_nA <- regionSegs_nA
    newList[[region]]$segs_nB <- regionSegs_nB
    
  }
  
  return(newList)
  
}


getASCATstyleCOPYseg <- function(copyITHRDS
                                 ,CCFestimates
                                 ,patientCPNtable
                                 ,purityEstimates
                                 ,ploidyEstimates
                                 ,clonalPvalue=0.001
                                 ,diffFromOtherPvalToUse =0.01
                                 #,mu=0.9
                                 ,plotting=TRUE
                                 ,minCCFtoConsider=0.05
                                 ,PvalToSeparate=0.001
                                 ,minOtherChrs=3
                                , gainThreshold    = log2(2.5/2)
                                , lossThreshold    = log2(1.5/2)
                                , ampThreshold     = log2(4/2)
                                , deepLossThreshold= log2(0.5/2)
                                , chrsToConsider=1:22
                                ,ploidy='mean'
                                ,LossGainPval =0.01
                                
                                )

{
  #we'll produce a few different outputs
  # 1 minimum consistent region matrix for each sample, with A and B, and CCF, and gain and loss, p-values etc recorded
  # ideally we also include the original ASCAT estimate here too. 
  # 2 joint segmented matrix for each sample for A and B separately (without CCF)
  
  today <- Sys.Date()
  today <- as.Date(today,format="%Y-%m-%d")
  dateToUse <- PasteVector(unlist(strsplit(as.character(today),split="-")),sep="")
  
  # select the sample of interest
  if(length(unique(patientCPNtable$sample))==1)
  {
    return(NA)
    #next;
  }
  regions <- unique(patientCPNtable$regionName)
  
  
  patientCPNtable$rawCopyNum  <- as.numeric(patientCPNtable$nBraw)+as.numeric(patientCPNtable$nAraw)
  patientCPNtable$minRaw      <- apply(cbind(patientCPNtable$nBraw,patientCPNtable$nAraw),1,min)
  patientCPNtable$maxRaw      <- apply(cbind(patientCPNtable$nBraw,patientCPNtable$nAraw),1,max)
  patientCPNtable$cpnToUse    <- patientCPNtable$cnTotal
  patientCPNtable$cpnPloidyCorrected <- NA
  
  purityEstimates <- patientCPNtable[!duplicated(patientCPNtable$sample),'ACF']
  names(purityEstimates) <- unique(patientCPNtable$sample )
  
  ploidyEstimates <- patientCPNtable[!duplicated(patientCPNtable$sample),'Ploidy']
  names(ploidyEstimates) <- unique(patientCPNtable$sample) 
  patientCPNtable$cpnToUse <-patientCPNtable$rawCopyNum 
  
  
  if(ploidy=='main'){
    samp.ploidy<- calc.ploidy(patientCPNtable)
    patientCPNtable$cpnPloidyCorrected <- log2(patientCPNtable$cpnToUse/samp.ploidy[match(patientCPNtable$sample,rownames(samp.ploidy)),2])
    
  }
  if(ploidy=='mean'){
    patientCPNtable$cpnPloidyCorrected <- log2(patientCPNtable$cpnToUse/patientCPNtable$Ploidy)
    
  } 
  
  patientmin.cons.CPN     <- min.cons.fn(patientCPNtable, chromosomes=1:22,colToUse = "cpnPloidyCorrected") 
  patientmin.cons.CPN[patientmin.cons.CPN == '-Inf'] <- -3
  patientmin.cons.CPN[patientmin.cons.CPN == 'Inf'] <- 3
  
  # first, do this the standard way
  PatientClonalDissectionAmpsAndDels <- freq_clonalcn.fn(min.cons.seg = patientmin.cons.CPN
                                                         ,sample.names = colnames(patientmin.cons.CPN[,4:ncol(patientmin.cons.CPN)])
                                                         ,type = 'all'
                                                         ,gain.threshold = log2(2.5/2)
                                                         ,loss.threshold = log2(1.5/2)
                                                         ,min.gain.threshold = log2(2.5/2)
                                                         ,min.loss.threshold = log2(1.5/2)
  )
  
  
  #define mu based on minCCFtoConsider
  mu <- 1-minCCFtoConsider
  
  # let's start with the easy one
  regions  <- grep(":",names(copyITHRDS$regionComparisonList),invert=TRUE,value=TRUE)
  segNamesToConsider <- names(copyITHRDS$regionComparisonList[[1]])
  #segNamesToConsider <- c()
  for (region in regions)
  {
    segNamesToConsider <- intersect(segNamesToConsider,names(copyITHRDS$regionComparisonList[[region]]))
  }
  segNamesToConsider <- unique(segNamesToConsider)
  #intersect(segNamesToConsider,names(CCFestimates))
  regionOne <- regions[1]
  
  # to permit a direct comparison of whether regions change in copy number, also perform a minimum consisten region analysis
  patientmin.cons.CPN_aRaw     <- min.cons.fn(patientCPNtable, chromosomes=1:22,colToUse = c('nAraw')) 
  patientmin.cons.CPN_bRaw     <- min.cons.fn(patientCPNtable, chromosomes=1:22,colToUse = c('nBraw')) 
  
  # let's also produce the nicolai style output here. 
  rownames(patientmin.cons.CPN_bRaw)<- rownames(patientmin.cons.CPN_aRaw) <- paste(patientmin.cons.CPN_aRaw[,1],patientmin.cons.CPN_aRaw[,2],patientmin.cons.CPN_aRaw[,3],sep=":")  
  
  
 
  
  
  
  
  # do this region at a time
  fullTable <- c()
  for (region in regions)
  {
    newCOPYtable <- c()
    for (segName in segNamesToConsider)
    {
      segDetailsRegionOne  <- copyITHRDS$regionComparisonList[[regionOne]][[segName]]
      segDetailsSpecRegion <- copyITHRDS$regionComparisonList[[region]][[segName]]
      segStatus            <- copyITHRDS$regionComparisonList[[region]][[segName]]$segStatus
      
      if(is.na(segDetailsRegionOne$nA[1]))
      {
        #stop('check this')
        sample               <- region
        segNameSplit         <- unlist(strsplit(segName,split=":"))
        chr                  <- as.numeric(segDetailsRegionOne$chr[1])
        startpos             <- as.numeric(segNameSplit[2])
        endpos               <- as.numeric(segNameSplit[3])
        n.het                <- 0#length(segDetailsSpecRegion$pos)
        cnTotal              <- median(segDetailsSpecRegion$cnTotal)#NA#ifelse(median(as.numeric(segDetailsSpecRegion$nA))<0,0,median(as.numeric(segDetailsSpecRegion$nA))+ifelse(median(as.numeric(segDetailsSpecRegion$nB))<0,0,median(as.numeric(segDetailsSpecRegion$nB))))
        nMajor               <- NA#round(ifelse(median(as.numeric(segDetailsSpecRegion$nA))<0,0,median(as.numeric(segDetailsSpecRegion$nA))))
        nMinor               <- NA#round(ifelse(median(as.numeric(segDetailsSpecRegion$nB))<0,0,median(as.numeric(segDetailsSpecRegion$nB))))
        Ploidy               <- ploidyEstimates[region]
        ACF                  <- purityEstimates[region]
        nAraw                <- NA#ifelse(median(as.numeric(segDetailsSpecRegion$nA))<0,0,median(as.numeric(segDetailsSpecRegion$nA)))
        nBraw                <- NA#ifelse(median(as.numeric(segDetailsSpecRegion$nB))<0,0,median(as.numeric(segDetailsSpecRegion$nB)))
        psi                  <- NA
        date                 <- dateToUse
        
        #additional features1
        nA1                 <-  NA#as.numeric(ccfDetailsAll$nA1)
        nA2                 <-  NA#as.numeric(ccfDetailsAll$nA2)
        nB1                 <-  NA#as.numeric(ccfDetailsAll$nB1)
        nB2                 <-  NA#as.numeric(ccfDetailsAll$nB2)
        fracA1              <-  NA#as.numeric(ccfDetailsSpecRegion$fracA1_mean)
        fracA2              <-  NA#as.numeric(ccfDetailsSpecRegion$fracA2_mean)
        fracB1              <-  NA#as.numeric(ccfDetailsSpecRegion$fracB1_mean)
        fracB2              <-  NA#as.numeric(ccfDetailsSpecRegion$fracB2_mean)
        origAraw <- patientmin.cons.CPN_aRaw[segName,region]
        origBraw <- patientmin.cons.CPN_bRaw[segName,region]
        
        subclonalA <- subclonalB <- heterogeneousA<-heterogeneousB <- obsBAF <- smallPvalA <- smallPvalB <- A_OtherChrsWithState<- A_OtherChrPval <-B_OtherChrsWithState <- B_OtherChrPval <- NA
        regionUsedToPhase <- 'tooFewSNPs'
        diffFromRef <- NA
        
        
        
        #additional features2
        logRObs            <- segStatus["logRObs"]
        gainPval           <- segStatus["gainPval"]
        ampPval            <- segStatus["ampPval"]
        lossPval           <- segStatus["lossPval"]
        DeepLossPval       <- segStatus["DeepLossPval"]
        
        rowOut <- cbind(sample,chr,startpos,endpos,n.het,cnTotal,nMajor,nMinor,Ploidy,ACF,nAraw,nBraw,psi,date,nA1,nA2,nB1,nB2,fracA1,fracA2,fracB1,fracB2,subclonalA,subclonalB,heterogeneousA,heterogeneousB,origAraw,origBraw
                        ,logRObs
                        ,obsBAF
                        ,gainPval
                        ,ampPval
                        ,lossPval
                        ,DeepLossPval
                        ,smallPvalA
                        ,smallPvalB
                        ,A_OtherChrsWithState
                        ,A_OtherChrPval
                        ,B_OtherChrsWithState
                        ,B_OtherChrPval
                        ,segName
                        ,regionUsedToPhase
                        ,diffFromRef)
        
        newCOPYtable <- rbind(newCOPYtable,rowOut)
        
      }
      
      
      if(!is.na(segDetailsRegionOne$nA[1]))
      {
      # check whether the segment shows evidence of MSAI
      regionUsedToPhase   <- copyITHRDS$whichRegionToPhaseFrom[[segName]]$regionToUse
      diffFromRef         <- copyITHRDS$regionComparisonList[[paste(regionUsedToPhase,region,sep=":")]][[segName]]$msAI
      if(length(diffFromRef)==0)
      {
        diffFromRef         <- copyITHRDS$regionComparisonList[[paste(region,regionUsedToPhase,sep=":")]][[segName]]$msAI
        
        if(length(diffFromRef)==0)
        {
        diffFromRef <- NA
        }
      }
      
      ccfDetailsSpecRegion <- CCFestimates[[segName]][[segName]][[region]]
      ccfDetailsAll        <- CCFestimates[[segName]][[segName]]
      ccfDetailsSpecRegion_wothers <- CCFestimates[[segName]]$diffFromExpected[[region]]
      
      # essentials
      # changing from mean to median 04/04/2018
      sample               <- region
      chr                  <- as.numeric(segDetailsRegionOne$chr[1])
      startpos             <- min(as.numeric(segDetailsRegionOne$pos))
      endpos               <- max(as.numeric(segDetailsRegionOne$pos))
      n.het                <- length(segDetailsSpecRegion$nA)
      cnTotal              <- ifelse(median(as.numeric(segDetailsSpecRegion$nA))<0,0,median(as.numeric(segDetailsSpecRegion$nA))+ifelse(median(as.numeric(segDetailsSpecRegion$nB))<0,0,median(as.numeric(segDetailsSpecRegion$nB))))
      nMajor               <- round(ifelse(median(as.numeric(segDetailsSpecRegion$nA))<0,0,median(as.numeric(segDetailsSpecRegion$nA))))
      nMinor               <- round(ifelse(median(as.numeric(segDetailsSpecRegion$nB))<0,0,median(as.numeric(segDetailsSpecRegion$nB))))
      Ploidy               <- ploidyEstimates[region]
      ACF                  <- purityEstimates[region]
      nAraw                <- ifelse(median(as.numeric(segDetailsSpecRegion$nA))<0,0,median(as.numeric(segDetailsSpecRegion$nA)))
      nBraw                <- ifelse(median(as.numeric(segDetailsSpecRegion$nB))<0,0,median(as.numeric(segDetailsSpecRegion$nB)))
      psi                  <- NA
      date                 <- dateToUse
      
      #additional features1
      nA1                 <-  as.numeric(ccfDetailsAll$nA1)
      nA2                 <-  as.numeric(ccfDetailsAll$nA2)
      nB1                 <-  as.numeric(ccfDetailsAll$nB1)
      nB2                 <-  as.numeric(ccfDetailsAll$nB2)
      fracA1              <-  as.numeric(ccfDetailsSpecRegion$fracA1_mean)
      fracA2              <-  as.numeric(ccfDetailsSpecRegion$fracA2_mean)
      fracB1              <-  as.numeric(ccfDetailsSpecRegion$fracB1_mean)
      fracB2              <-  as.numeric(ccfDetailsSpecRegion$fracB2_mean)
      
      #additional features2
      logRObs            <- segStatus["logRObs"]
      gainPval           <- segStatus["gainPval"]
      ampPval            <- segStatus["ampPval"]
      lossPval           <- segStatus["lossPval"]
      DeepLossPval       <- segStatus["DeepLossPval"]
      
      #also get values for the expected BAF
      newBAFtable<- copyITHRDS$newBAFtable[[region]]$regionBAFtable 
      BAFsegment  <- newBAFtable[newBAFtable$chrom==chr&newBAFtable$pos>=startpos&newBAFtable$pos<=endpos,,drop=FALSE]
      obsBAF          <-  median(BAFsegment$newBAFsegmented)
      explogR_nA1nB1  <-  funCalcExpectedLogR(exp_intA = nA1,exp_intB = nB1,tumorPloidy =Ploidy,tumorPurity = ACF)
      explogR_nA1nB2  <-  funCalcExpectedLogR(exp_intA = nA1,exp_intB = nB2,tumorPloidy =Ploidy,tumorPurity = ACF)
      explogR_nA2nB1  <-  funCalcExpectedLogR(exp_intA = nA2,exp_intB = nB1,tumorPloidy =Ploidy,tumorPurity = ACF)
      explogR_nA2nB2  <-  funCalcExpectedLogR(exp_intA = nA2,exp_intB = nB2,tumorPloidy =Ploidy,tumorPurity = ACF)
      
      expBAF_nA1nB1  <-  funCalcExpectedBAF(exp_intA = nA1,exp_intB = nB1,tumorPloidy =Ploidy,tumorPurity = ACF)
      expBAF_nA1nB2  <-  funCalcExpectedBAF(exp_intA = nA1,exp_intB = nB2,tumorPloidy =Ploidy,tumorPurity = ACF)
      expBAF_nA2nB1  <-  funCalcExpectedBAF(exp_intA = nA2,exp_intB = nB1,tumorPloidy =Ploidy,tumorPurity = ACF)
      expBAF_nA2nB2  <-  funCalcExpectedBAF(exp_intA = nA2,exp_intB = nB2,tumorPloidy =Ploidy,tumorPurity = ACF)
      
      
      
      #testfor Subclonality
      # check if the larger fraction is significantly less than mu
      
      
      if(median(ccfDetailsSpecRegion$fracA1_all)>median(ccfDetailsSpecRegion$fracA2_all))
      {
        TtestAsubclonal <- wilcox.test(ccfDetailsSpecRegion$fracA1_all,mu=mu,alternative = 'less')
      }
      
      if(median(ccfDetailsSpecRegion$fracA1_all)<=median(ccfDetailsSpecRegion$fracA2_all))
      {
        TtestAsubclonal <- wilcox.test(ccfDetailsSpecRegion$fracA2_all,mu=mu,alternative = 'less')
      }
      if(median(ccfDetailsSpecRegion$fracB1_all)>median(ccfDetailsSpecRegion$fracB2_all))
      {
        TtestBsubclonal <- wilcox.test(ccfDetailsSpecRegion$fracB1_all,mu=mu,alternative = 'less')
      }
      
      if(median(ccfDetailsSpecRegion$fracB1_all)<=median(ccfDetailsSpecRegion$fracB2_all))
      {
        TtestBsubclonal <- wilcox.test(ccfDetailsSpecRegion$fracB2_all,mu=mu,alternative = 'less')
      }
      
      # also check whether evidence the segment is different from others
      nA_PvalFromOthers <- 1e-120
      if(length(ccfDetailsSpecRegion_wothers$nA_nrMatchedAChrs)>=minOtherChrs)
      {
        nA_PvalFromOthers <- ccfDetailsSpecRegion_wothers$nApval_fromOtherChrs
      }
      nB_PvalFromOthers <- 1e-120
      if( length(ccfDetailsSpecRegion_wothers$nB_nrMatchedAChrs)>=minOtherChrs)
      {
        nB_PvalFromOthers <- ccfDetailsSpecRegion_wothers$nBpval_fromOtherChrs
      }
      
      A_OtherChrsWithState <- length(ccfDetailsSpecRegion_wothers$nA_nrMatchedAChrs)
      A_OtherChrPval       <- nA_PvalFromOthers
      
      B_OtherChrsWithState <- length(ccfDetailsSpecRegion_wothers$nB_nrMatchedAChrs)
      B_OtherChrPval       <- nB_PvalFromOthers
      
      # determine whether any evidence for subclonal A or B
      subclonalA         <- TRUE
      subclonalB         <- TRUE
      heterogeneousA     <- TRUE
      heterogeneousB     <- TRUE
      smallPvalA         <- TtestAsubclonal$p.value
      smallPvalB         <- TtestBsubclonal$p.value
      
      
      if(nA1==nA2)
      {
        subclonalA         <- FALSE
        heterogeneousA     <- FALSE
      }
      if(nB1==nB2)
      {
        subclonalB         <- FALSE
        heterogeneousB     <- FALSE
      }
      if(fracA1>=1|fracA2>=1|smallPvalA>=clonalPvalue|nA_PvalFromOthers>=diffFromOtherPvalToUse)
      {
        subclonalA <- FALSE
      }
      if(fracB1>=1|fracB2>=1|smallPvalB>=clonalPvalue|nB_PvalFromOthers>=diffFromOtherPvalToUse)
      {
        subclonalB <- FALSE
      }
      
      origAraw <- patientmin.cons.CPN_aRaw[segName,region]
      origBraw <- patientmin.cons.CPN_bRaw[segName,region]
      
      
      # now put together into a table
      rowOut <- cbind(sample,chr,startpos,endpos,n.het,cnTotal,nMajor,nMinor,Ploidy,ACF,nAraw,nBraw,psi,date,nA1,nA2,nB1,nB2,fracA1,fracA2,fracB1,fracB2,subclonalA,subclonalB,heterogeneousA,heterogeneousB,origAraw,origBraw
                      ,logRObs
                      ,obsBAF
                      ,gainPval
                      ,ampPval
                      ,lossPval
                      ,DeepLossPval
                      ,smallPvalA
                      ,smallPvalB
                      ,A_OtherChrsWithState
                      ,A_OtherChrPval
                      ,B_OtherChrsWithState
                      ,B_OtherChrPval
                      ,segName
                      ,regionUsedToPhase
                      ,diffFromRef)
      
      
    
      newCOPYtable <- rbind(newCOPYtable,rowOut)
      
      }
      
    }
    fullTable <- rbind(fullTable,newCOPYtable)
  }
  
  rownames(fullTable) <- 1:nrow(fullTable)
  fullTable <- data.frame(fullTable,stringsAsFactors = FALSE)
  
  #let's add our new calls to this one
  allZeros <- PatientClonalDissectionAmpsAndDels$all_gain
  allZeros[,4] <- 0
  rownames(allZeros) <- paste(allZeros[,1],allZeros[,2],allZeros[,3],sep=":")
  
  PatientClonalDissectionAmpsAndDels$all_gain_new <- PatientClonalDissectionAmpsAndDels$all_loss_new <- PatientClonalDissectionAmpsAndDels$trunk_gain_new <- PatientClonalDissectionAmpsAndDels$trunk_loss_new <- PatientClonalDissectionAmpsAndDels$branch_gain_new <- PatientClonalDissectionAmpsAndDels$branch_loss_new <- allZeros 

  # check each segment whether it is gained or lost
  for (segName in unique(fullTable$segName))
  {
    
    
    segRegions <- fullTable[fullTable$segName%in%segName,,drop=FALSE]
    if(TRUE%in%is.na(as.numeric(segRegions$gainPval)))
    {
      next
    }
    
    if(min(as.numeric(segRegions$gainPval),as.numeric(segRegions$lossPval))>LossGainPval)
    {
      next;
    }

    if(min(as.numeric(segRegions$gainPval))<=LossGainPval)
    {
      PatientClonalDissectionAmpsAndDels$all_gain_new[segName,4] <- 1
    }
    if(min(as.numeric(segRegions$lossPval))<=LossGainPval)
    {
      PatientClonalDissectionAmpsAndDels$all_loss_new[segName,4] <- -1
    }
    
    #is it ubiquitously gained?
    if(max(as.numeric(segRegions$gainPval))<=LossGainPval)
    {
      PatientClonalDissectionAmpsAndDels$trunk_gain_new[segName,4] <- 1
      next;
    }
    
    #is it ubiquitously loss?
    if(max(as.numeric(segRegions$lossPval))<=LossGainPval)
    {
      PatientClonalDissectionAmpsAndDels$trunk_loss_new[segName,4] <- -1
      next;
    }
    
    #is it heterogeneous gain?
    if(min(as.numeric(segRegions$gainPval))<=LossGainPval)
    {
      PatientClonalDissectionAmpsAndDels$branch_gain_new[segName,4] <- 1
    }
    
    #is it heterogeneous loss?
    if(min(as.numeric(segRegions$lossPval))<=LossGainPval)
    {
      PatientClonalDissectionAmpsAndDels$branch_loss_new[segName,4] <- -1
    }
    
    
  }
  
  # oldTrunk <- fullTable[fullTable$segName%in%rownames(PatientClonalDissectionAmpsAndDels$trunk_gain_new[PatientClonalDissectionAmpsAndDels$trunk_gain[,4]==1&PatientClonalDissectionAmpsAndDels$trunk_gain_new[,4]!=1,,drop=FALSE]),,drop=FALSE]
  # newTrunk <- fullTable[fullTable$segName%in%rownames(PatientClonalDissectionAmpsAndDels$trunk_gain_new[PatientClonalDissectionAmpsAndDels$trunk_gain[,4]!=1&PatientClonalDissectionAmpsAndDels$trunk_gain_new[,4]==1,,drop=FALSE]),,drop=FALSE]
  # 
  # oldTrunkLoss <- fullTable[fullTable$segName%in%rownames(PatientClonalDissectionAmpsAndDels$trunk_loss_new[PatientClonalDissectionAmpsAndDels$trunk_loss[,4]==(-1)&PatientClonalDissectionAmpsAndDels$trunk_loss_new[,4]!=(-1),,drop=FALSE]),,drop=FALSE]
  # newTrunkLoss <- fullTable[fullTable$segName%in%rownames(PatientClonalDissectionAmpsAndDels$trunk_loss_new[PatientClonalDissectionAmpsAndDels$trunk_loss[,4]!=(-1)&PatientClonalDissectionAmpsAndDels$trunk_loss_new[,4]==(-1),,drop=FALSE]),,drop=FALSE]
  # 
  # 
  # newTrunk <- newTrunk[order(newTrunk$segName),,drop=FALSE]
  # oldTrunk <- oldTrunk[order(oldTrunk$segName),,drop=FALSE]
  # 
  # 
  # newTrunkLoss <- newTrunkLoss[order(newTrunkLoss$segName),,drop=FALSE]
  # oldTrunkLoss <- oldTrunkLoss[order(oldTrunkLoss$segName),,drop=FALSE]
  
  
  if(plotting==TRUE)
  {
    
    test <- fullTable[fullTable$subclonalA==TRUE,,drop=FALSE]
    test <- fullTable[fullTable$subclonalA==FALSE,,drop=FALSE]
    
    par(mfrow=c(2,2))
    for (region in regions)
    {
      regionTable <- fullTable[fullTable$sample%in%region,,drop=FALSE]
      #defineYmax and Xmax
      AmaxVal <- min(c(5,ceiling(max(c(as.numeric(regionTable$nAraw),as.numeric(regionTable$origAraw)),na.rm = TRUE))))
      BmaxVal <- min(c(5,ceiling(max(c(as.numeric(regionTable$nBraw),as.numeric(regionTable$origBraw)),na.rm=TRUE))))
      plot( ifelse(as.numeric(regionTable$nAraw)>=AmaxVal,AmaxVal,as.numeric(regionTable$nAraw))
            ,ifelse(as.numeric(regionTable$origAraw)>=AmaxVal,AmaxVal,as.numeric(regionTable$origAraw))
            ,pch=16,cex=log(as.numeric(regionTable$n.het),5),col=ifelse(regionTable$subclonalA=="TRUE",'#de2d2699','#3182bd99'),las=1
            ,xlab="new Araw "
            ,ylab="ASCAT Araw"
            ,xlim=c(0,AmaxVal)
            ,ylim=c(0,AmaxVal)
            ,main=region)
      
      abline(0,1)
      abline(h=seq(0,20,by=1),lty='dashed')
      abline(v=seq(0,20,by=1),lty='dashed')
      plot(ifelse(as.numeric(regionTable$nBraw)>=BmaxVal,BmaxVal,as.numeric(regionTable$nBraw))
           ,ifelse(as.numeric(regionTable$origBraw)>=BmaxVal,BmaxVal,as.numeric(regionTable$origBraw))
           ,pch=16,cex=log(as.numeric(regionTable$n.het),10),col=ifelse(regionTable$subclonalB=="TRUE",'#de2d2699','#3182bd99'),las=1
           ,xlab="new Braw "
           ,ylab="ASCAT Braw"
           ,xlim=c(0,BmaxVal)
           ,ylim=c(0,BmaxVal)
           ,main=region)
      abline(h=seq(0,20,by=1),lty='dashed')
      abline(v=seq(0,20,by=1),lty='dashed')
      abline(0,1)
      
    }
    
    # next, compare each and every region to each other
    for (x in 1:(length(regions)-1))
    {
      for (y in (x+1):length(regions))
      {
        regionX <- regions[x]
        regionY <- regions[y]
        
        
        
        regionXTable <- fullTable[fullTable$sample%in%regionX,,drop=FALSE]
        regionYTable <- fullTable[fullTable$sample%in%regionY,,drop=FALSE]
        
        
        AmaxVal <- min(c(5,ceiling(max(c(as.numeric(regionXTable$nAraw),as.numeric(regionYTable$nAraw)),na.rm = TRUE))))
        BmaxVal <- min(c(5,ceiling(max(c(as.numeric(regionXTable$nBraw),as.numeric(regionYTable$nBraw)),na.rm = TRUE))))
        regionXTable$nAraw <- ifelse(as.numeric(regionXTable$nAraw)>AmaxVal,AmaxVal,as.numeric(regionXTable$nAraw))
        regionYTable$nAraw <- ifelse(as.numeric(regionYTable$nAraw)>AmaxVal,AmaxVal,as.numeric(regionYTable$nAraw))
        regionXTable$nBraw <- ifelse(as.numeric(regionXTable$nBraw)>BmaxVal,BmaxVal,as.numeric(regionXTable$nBraw))
        regionYTable$nBraw <- ifelse(as.numeric(regionYTable$nBraw)>BmaxVal,BmaxVal,as.numeric(regionYTable$nBraw))
        
        plot(regionXTable$nAraw
             ,regionYTable$nAraw
             ,pch=16,cex=log(as.numeric(regionTable$n.het),10),col=ifelse(regionXTable$subclonalA=="TRUE"|regionYTable$subclonalA=="TRUE",'#de2d2699','#3182bd99'),las=1
             ,xlab=paste(regionX,"Araw")
             ,ylab=paste(regionY,"Araw")
             ,main=paste(regionX,regionY)
             ,xlim=c(0,AmaxVal)
             ,ylim=c(0,AmaxVal))
        abline(0,1)
        abline(h=seq(0,20,by=1),lty='dashed')
        abline(v=seq(0,20,by=1),lty='dashed')
        
        
        # regionXTable <- fullTable[fullTable$sample%in%regionX,,drop=FALSE]
        # regionYTable <- fullTable[fullTable$sample%in%regionY,,drop=FALSE]
        plot(regionXTable$nBraw,regionYTable$nBraw,pch=16,cex=log(as.numeric(regionTable$n.het),10),col=ifelse(regionXTable$subclonalB=="TRUE"|regionYTable$subclonalB=="TRUE",'#de2d2699','#3182bd99'),las=1
             ,xlab=paste(regionX,"Braw")
             ,ylab=paste(regionY,"Braw")
             ,main=paste(regionX,regionY)
             ,xlim=c(0,BmaxVal)
             ,ylim=c(0,BmaxVal))
        abline(0,1)
        abline(h=seq(0,20,by=1),lty='dashed')
        abline(v=seq(0,20,by=1),lty='dashed')
        
        
        plot(regionXTable$fracA2,regionYTable$fracA2,pch=16,cex=log(as.numeric(regionTable$n.het),10),col=ifelse(regionXTable$subclonalA=="TRUE"|regionYTable$subclonalA=="TRUE",'#de2d2699','#3182bd99'),las=1
             ,xlab=paste(regionX,"fracA2")
             ,ylab=paste(regionY,"fracA2")
             ,main=paste(regionX,regionY)
             ,xlim=c(0,1.5)
             ,ylim=c(0,1.5))
        abline(0,1)
        abline(h=seq(0,20,by=1),lty='dashed')
        abline(v=seq(0,20,by=1),lty='dashed')
        
        plot(regionXTable$fracB2,regionYTable$fracB2,pch=16,cex=log(as.numeric(regionTable$n.het),10),col=ifelse(regionXTable$subclonalB=="TRUE"|regionYTable$subclonalB=="TRUE",'#de2d2699','#3182bd99'),las=1
             ,xlab=paste(regionX,"fracB2")
             ,ylab=paste(regionY,"fracB2")
             ,main=paste(regionX,regionY)
             ,xlim=c(0,1.5)
             ,ylim=c(0,1.5))
        abline(0,1)
        abline(h=seq(0,20,by=1),lty='dashed')
        abline(v=seq(0,20,by=1),lty='dashed')
        
        
        
        
        
      }
    }
  }
  
  # next, extract the segmented matrix (with joining)
  Regions     <- grep(":",names(copyITHRDS$regionComparisonList),invert=TRUE,value=TRUE)
  numRegions  <- length(Regions)
  firstRegion <- Regions[1] 
  lastRegion      <- Regions[length(Regions)] 
  #regionJoinedSeg <- funJoinAdjacentSegments(copyITHRDS$regionComparisonList,PvalToSeparate = PvalToSeparate)
  regionJoinedSegCCF <- funJoinAdjacentSegments_withCCF(regionComparisonList = copyITHRDS$regionComparisonList
                                                        ,CCFestimates = CCFestimates
                                                        ,PvalToSeparate = PvalToSeparate
  )
  
  
  #create a new comparison matrix with these segs
  regionChrList <- list()
  # next assess for each chromosome and region, what are the different raw-values
  for (region in regions)
  {
    
    #regionChrList[[region]] <- list()
    
    for (segName in segNamesToConsider)
    {
      intNumber <- round(median(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nA)))
      chrNumber <- unlist(strsplit(segName,split=":"))[1]
      regionChrList[[region]][[as.character(intNumber)]][[chrNumber]]$nAvals <- c(regionChrList[[region]][[as.character(intNumber)]][[chrNumber]]$nAvals,copyITHRDS$regionComparisonList[[region]][[segName]]$nA)
    }
    
    for (segName in segNamesToConsider)
    {
      intNumber <- round(median(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nB)))
      chrNumber <- unlist(strsplit(segName,split=":"))[1]
      regionChrList[[region]][[as.character(intNumber)]][[chrNumber]]$nBvals <- c(regionChrList[[region]][[as.character(intNumber)]][[chrNumber]]$nBvals,copyITHRDS$regionComparisonList[[region]][[segName]]$nB)
    }
    
  }
  
  segMatCopyCombinedCCF          <- c()
  segMatCopy_nA_allCCF           <- c()
  segMatCopy_nB_allCCF           <- c()
  for (region in names(regionJoinedSegCCF))
  {
    segMatCopy_nA_CCF <- c()
    segMatCopy_nB_CCF <- c()
    for (segName in unique(names(regionJoinedSegCCF[[region]][['segs_nA']])))
    {
      
      nAraw    <- median(as.numeric(regionJoinedSegCCF[[region]][['segs_nA']][[segName]]$nA))
      chr      <- mean(as.numeric(regionJoinedSegCCF[[region]][['segs_nA']][[segName]]$chr))
      startpos <- min(as.numeric(regionJoinedSegCCF[[region]][['segs_nA']][[segName]]$pos))
      endpos   <- max(as.numeric(regionJoinedSegCCF[[region]][['segs_nA']][[segName]]$pos))
      sample   <- region
      
      segDetailsSpecRegion <- regionJoinedSegCCF[[region]]$segs_nA[[segName]]
      
      # essentials
      sample               <- region
      chrName              <- as.numeric(segDetailsSpecRegion$chr[1])
      startpos             <- min(as.numeric(segDetailsSpecRegion$pos))
      endpos               <- max(as.numeric(segDetailsSpecRegion$pos))
      n.het                <- length(segDetailsSpecRegion$nA)
      cnTotal              <- NA #ifelse(as.numeric(segDetailsSpecRegion$nA_mean)<0,0,as.numeric(segDetailsSpecRegion$nA_mean))+ifelse(as.numeric(segDetailsSpecRegion$nB_mean)<0,0,as.numeric(segDetailsSpecRegion$nB_mean))
      nMajor               <- NA #round(ifelse(as.numeric(segDetailsSpecRegion$nA_mean)<0,0,as.numeric(segDetailsSpecRegion$nA_mean)))
      nMinor               <- NA #round(ifelse(as.numeric(segDetailsSpecRegion$nB_mean)<0,0,as.numeric(segDetailsSpecRegion$nB_mean)))
      Ploidy               <- ploidyEstimates[region]
      ACF                  <- purityEstimates[region]
      nAraw                <- ifelse(median(segDetailsSpecRegion$nA)<0,0,median(segDetailsSpecRegion$nA))
      nBraw                <- NA #ifelse(segDetailsSpecRegion$nB_mean<0,0,segDetailsSpecRegion$nB_mean)
      psi                  <- NA
      date                 <- dateToUse
      
      #additional features
      nA1                 <-  median(as.numeric(segDetailsSpecRegion$nA1))
      nA2                 <-  median(as.numeric(segDetailsSpecRegion$nA2))# as.numeric(ccfDetailsAll$nA2)
      nB1                 <-  NA#as.numeric(ccfDetailsAll$nB1)
      nB2                 <-  NA#as.numeric(ccfDetailsAll$nB2)
      fracA1              <-  median(as.numeric(segDetailsSpecRegion$nAfrac)) #as.numeric(ccfDetailsSpecRegion$fracA1_mean)
      fracA2              <-  median(as.numeric(segDetailsSpecRegion$nA2frac))
      #fracA2              <-  ifelse(1-fracA1<0,0,1-fracA1)
      #as.numeric(ccfDetailsSpecRegion$fracA2_mean)
      fracB1              <-  NA #as.numeric(ccfDetailsSpecRegion$fracB1_mean)
      fracB2              <-  NA #¢as.numeric(ccfDetailsSpecRegion$fracB2_mean)
      
      
      #additional features2
      #also get values for the expected BAF
      newBAFtable<- copyITHRDS$newBAFtable[[region]]$regionBAFtable 
      BAFsegment  <- newBAFtable[newBAFtable$chrom==chr&newBAFtable$pos>=startpos&newBAFtable$pos<=endpos,,drop=FALSE]
      
      checkGainSignificance <- function(logRBAFRegion,regionToCheck,ploidyEstimates,purityEstimates,gainThreshold)
      {
        minGainRawExpRegion       <-  ploidyEstimates[regionToCheck]*(2^(gainThreshold))
        minGainLogRExpRegion      <-  funCalcExpectedLogR(minGainRawExpRegion,0,ploidyEstimates[regionToCheck],purityEstimates[regionToCheck])
        logRObs                   <-  mean(logRBAFRegion$logR)
        
        if(nrow(logRBAFRegion)>=3)
        {
          Ttest        <- tTestWithErrorHandling(logRBAFRegion$logR,mu = minGainLogRExpRegion)
          TtestGreater <- tTestWithErrorHandling(logRBAFRegion$logR,mu = minGainLogRExpRegion,alternative = 'greater')
          #TtestLess    <- tTestWithErrorHandling(logRBAF$logR,mu = minLossLogRExp,alternative = 'less')
          logRObs_low  <- Ttest$conf.int[1]
          logRObs_high <- Ttest$conf.int[2]
          PVal         <- Ttest$p.value
          PvalGain     <- TtestGreater$p.value
          
        }
        if(nrow(logRBAFRegion)<3)
        {
          logRObs_low  <- NA
          logRObs_high <- NA
          PVal         <- NA
          PvalGain     <- NA
          
          
        }
        
        out <- cbind(minGainRawExpRegion,minGainLogRExpRegion,logRObs,logRObs_low,logRObs_high,PVal,PvalGain)  
        return(out)
        
        
      }
      checkLossSignificance <- function(logRBAFRegion,regionToCheck,ploidyEstimates,purityEstimates,lossThreshold)
      {
        minLossRawExpRegion       <-  ploidyEstimates[regionToCheck]*(2^(lossThreshold))
        minLossLogRExpRegion      <-  funCalcExpectedLogR(minLossRawExpRegion,0,ploidyEstimates[regionToCheck],purityEstimates[regionToCheck])
        logRObs                   <- mean(logRBAFRegion$logR)
        
        if(nrow(logRBAFRegion)>=3)
        {
          Ttest        <- tTestWithErrorHandling(logRBAFRegion$logR,mu = minLossLogRExpRegion)
          TtestLess <- tTestWithErrorHandling(logRBAFRegion$logR,mu = minLossLogRExpRegion,alternative = 'less')
          #TtestLess    <- tTestWithErrorHandling(logRBAF$logR,mu = minLossLogRExp,alternative = 'less')
          logRObs_low  <- Ttest$conf.int[1]
          logRObs_high <- Ttest$conf.int[2]
          PVal         <- Ttest$p.value
          PvalLoss     <- TtestLess$p.value
          
        }
        # if(nrow(logRBAF)<3)
        # {
        #   logRObs_low  <- NA
        #   logRObs_high <- NA
        #   PVal         <- NA
        #   PvalLoss     <- NA
        
        
        # }
        
        out <- cbind(minLossRawExpRegion,minLossLogRExpRegion,logRObs,logRObs_low,logRObs_high,PVal,PvalLoss)  
        return(out)
        
        
      }
      

      GainRegionOne     <- checkGainSignificance(logRBAFRegion = BAFsegment
                                                 ,regionToCheck = region
                                                 ,ploidyEstimates
                                                 ,purityEstimates
                                                 ,gainThreshold = gainThreshold)
      AmpRegionOne      <- checkGainSignificance(BAFsegment,regionToCheck = region,ploidyEstimates,purityEstimates,gainThreshold = ampThreshold)
      LossRegionOne     <- checkLossSignificance(BAFsegment,regionToCheck = region,ploidyEstimates,purityEstimates,lossThreshold = lossThreshold)
      DeepLossRegionOne <- checkLossSignificance(BAFsegment,regionToCheck = region,ploidyEstimates,purityEstimates,lossThreshold = deepLossThreshold)
      
      
      # put together into One table
      segStatus <- c(GainRegionOne[,3:5,drop=FALSE]
                               ,GainRegionOne[,2,drop=FALSE],GainRegionOne[,7,drop=FALSE]
                               ,AmpRegionOne[,2,drop=FALSE],AmpRegionOne[,7,drop=FALSE]
                               ,LossRegionOne[,2,drop=FALSE],LossRegionOne[,7,drop=FALSE]
                               ,DeepLossRegionOne[,2,drop=FALSE],DeepLossRegionOne[,7,drop=FALSE]
      )
      names(segStatus) <- c('logRObs','logRObs_low','logRObs_high'
                                      ,'minGainExpLogR','gainPval'
                                      ,'minAmpExpLogR','ampPval'
                                      ,'minLossExpLogR','lossPval'
                                      ,'minDeepLossExpLogR','DeepLossPval')
      
      logRObs            <- median(as.numeric(BAFsegment$BAF))
      gainPval           <- segStatus["gainPval"]
      ampPval            <- segStatus["ampPval"]
      lossPval           <- segStatus["lossPval"]
      DeepLossPval       <- segStatus["DeepLossPval"]
      
      obsBAF          <-  median(BAFsegment$newBAFsegmented)

      # A1 ttest
      if(length(unique(segDetailsSpecRegion$nAfrac))==1)
      {
        if(unique(segDetailsSpecRegion$nAfrac)==1)
        {
          TtestA_subclonal <- list()
          TtestA_subclonal$p.value <- 1
        }
        
      }
      
      if(length(unique(segDetailsSpecRegion$nAfrac))>1)
      {
        # check whether to use A1 or A2
        if(fracA1>fracA2)
        {
          TtestA_subclonal   <- wilcox.test(segDetailsSpecRegion$nA1frac,mu=mu,alternative = "less")
          
          
        }
        if(fracA1<fracA2)
        {
          TtestA_subclonal   <- wilcox.test(segDetailsSpecRegion$nA2frac,mu=mu,alternative = "less")
          
        }
        
      }
      
      # also check whether it is significantly different from other chromosomes
      if(fracA1>fracA2)
      {
         #choose nA1 as integer reference
         rawAs   <- segDetailsSpecRegion$nA
         intNumA <- nA1
         
         fullAs <- c()
         nrMatchedA_Chrs <- c()
         for (chr in chrsToConsider)
         {
           if(length(regionChrList[[region]][[as.character(intNumA)]][[as.character(chr)]]$nAvals)==0)
           {
             next;
           }
           if(length(regionChrList[[region]][[as.character(intNumA)]][[as.character(chr)]]$nAvals)!=0)
           {
             nrMatchedA_Chrs <- c(nrMatchedA_Chrs,as.character(chr))
             fullAs <- c(fullAs,regionChrList[[region]][[as.character(intNumA)]][[as.character(chr)]]$nAvals)
           }
           
         }
         TtestAs             <- tTestWithErrorHandling(rawAs,fullAs)
         TtestDiffFromIntA   <- tTestWithErrorHandling(rawAs,mu=intNumA)
         nrMatchedA_Chrs     <- length(nrMatchedA_Chrs)
         
     }
      
      
      if(fracA1<=fracA2)
      {
        #choose nA1 as integer reference
        rawAs   <- segDetailsSpecRegion$nA
        intNumA <- nA2
        
        fullAs <- c()
        nrMatchedA_Chrs <- c()
        for (chr in chrsToConsider)
        {
          if(length(regionChrList[[region]][[as.character(intNumA)]][[as.character(chr)]]$nAvals)==0)
          {
            next;
          }
          if(length(regionChrList[[region]][[as.character(intNumA)]][[as.character(chr)]]$nAvals)!=0)
          {
            nrMatchedA_Chrs <- c(nrMatchedA_Chrs,as.character(chr))
            fullAs <- c(fullAs,regionChrList[[region]][[as.character(intNumA)]][[as.character(chr)]]$nAvals)
          }
          
        }
        TtestAs             <- tTestWithErrorHandling(rawAs,fullAs)
        TtestDiffFromIntA   <- tTestWithErrorHandling(rawAs,mu=intNumA)
        nrMatchedA_Chrs     <- length(nrMatchedA_Chrs)
        
      }  
      
      
      # determine whether any evidence for subclonal A or B
      subclonalA         <- TRUE
      subclonalB         <- NA
      heterogeneousA     <- TRUE
      heterogeneousB     <- NA
      diffFromOtherPval  <- 1e-100
     
      if(nrMatchedA_Chrs>=minOtherChrs)
      {
        diffFromOtherPval <- TtestAs$p.value
      }
      # if(nB1==nB2)
      # {
      #   subclonalB         <- FALSE
      #   heterogeneousB     <- FALSE
      # }
      smallerPval <- TtestA_subclonal$p.value
      
      if(fracA1>=mu|fracA2>=mu|smallerPval>=clonalPvalue|diffFromOtherPval>diffFromOtherPvalToUse)
      {
        subclonalA <- FALSE
      }
      if(nA1==nA2)
      {
        subclonalA         <- FALSE
        heterogeneousA     <- FALSE
      }
      
      origAraw <- NA
      origBraw <- NA
      smallPvalA <- smallerPval
      smallPvalB <- NA
      A_OtherChrsWithState <- nrMatchedA_Chrs
      A_OtherChrPval       <- diffFromOtherPval
      B_OtherChrsWithState <- NA
      B_OtherChrPval       <- NA
      
      chr <- chrName
      rowOut <- cbind(sample,chr,startpos,endpos,n.het,cnTotal,nMajor,nMinor,Ploidy,ACF,nAraw,nBraw,psi,date,nA1,nA2,nB1,nB2,fracA1,fracA2,fracB1,fracB2,subclonalA,subclonalB,heterogeneousA,heterogeneousB,origAraw,origBraw,segName
                      ,logRObs
                      ,obsBAF
                      ,gainPval
                      ,ampPval
                      ,lossPval
                      ,DeepLossPval
                      ,smallPvalA
                      ,smallPvalB
                      ,A_OtherChrsWithState
                      ,A_OtherChrPval
                      ,B_OtherChrsWithState
                      ,B_OtherChrPval)
      
      # segOut   <- cbind(sample,chr,startpos,endpos,segName,nAraw)
      segMatCopy_nA_CCF <- rbind(segMatCopy_nA_CCF,rowOut)
    }
    
    for (segName in unique(names(regionJoinedSegCCF[[region]][['segs_nB']])))
    {
      
      nBraw    <- median(as.numeric(regionJoinedSegCCF[[region]][['segs_nB']][[segName]]$nB))
      chrName      <- mean(as.numeric(regionJoinedSegCCF[[region]][['segs_nB']][[segName]]$chr))
      startpos <- min(as.numeric(regionJoinedSegCCF[[region]][['segs_nB']][[segName]]$pos))
      endpos   <- max(as.numeric(regionJoinedSegCCF[[region]][['segs_nB']][[segName]]$pos))
      sample   <- region
      
      
      
      segDetailsSpecRegion <- regionJoinedSegCCF[[region]]$segs_nB[[segName]]
      
      # essentials
      sample               <- region
      chr                  <- as.numeric(segDetailsSpecRegion$chr[1])
      startpos             <- min(as.numeric(segDetailsSpecRegion$pos))
      endpos               <- max(as.numeric(segDetailsSpecRegion$pos))
      n.het                <- length(segDetailsSpecRegion$nB)
      cnTotal              <- NA #ifelse(as.numeric(segDetailsSpecRegion$nA_mean)<0,0,as.numeric(segDetailsSpecRegion$nA_mean))+ifelse(as.numeric(segDetailsSpecRegion$nB_mean)<0,0,as.numeric(segDetailsSpecRegion$nB_mean))
      nMajor               <- NA #round(ifelse(as.numeric(segDetailsSpecRegion$nA_mean)<0,0,as.numeric(segDetailsSpecRegion$nA_mean)))
      nMinor               <- NA #round(ifelse(as.numeric(segDetailsSpecRegion$nB_mean)<0,0,as.numeric(segDetailsSpecRegion$nB_mean)))
      Ploidy               <- ploidyEstimates[region]
      ACF                  <- purityEstimates[region]
      nAraw                <- NA #ifelse(mean(segDetailsSpecRegion$nA)<0,0,mean(segDetailsSpecRegion$nA))
      nBraw                <- ifelse(median(segDetailsSpecRegion$nB)<0,0,median(segDetailsSpecRegion$nB))#NA #ifelse(segDetailsSpecRegion$nB_mean<0,0,segDetailsSpecRegion$nB_mean)
      psi                  <- NA
      date                 <- dateToUse
      
      #additional features
      nA1                 <-  NA#mean(as.numeric(segDetailsSpecRegion$nA1))
      nA2                 <-  NA#mean(as.numeric(segDetailsSpecRegion$nA2))# as.numeric(ccfDetailsAll$nA2)
      nB1                 <-  median(as.numeric(segDetailsSpecRegion$nB1))#NA#as.numeric(ccfDetailsAll$nB1)
      nB2                 <-  median(as.numeric(segDetailsSpecRegion$nB2))#NA#as.numeric(ccfDetailsAll$nB2)
      fracA1              <-  NA#median(as.numeric(segDetailsSpecRegion$nAfrac)) #as.numeric(ccfDetailsSpecRegion$fracA1_mean)
      fracA2              <-  NA#ifelse(1-fracA1<0,0,1-fracA1)
      #as.numeric(ccfDetailsSpecRegion$fracA2_mean)
      fracB1              <-  median(as.numeric(segDetailsSpecRegion$nB1frac)) #as.numeric(ccfDetailsSpecRegion$fracA1_mean)#NA #as.numeric(ccfDetailsSpecRegion$fracB1_mean)
      fracB2              <-  median(as.numeric(segDetailsSpecRegion$nB2frac))#ifelse(1-fracB1<0,0,1-fracB1)#NA #¢as.numeric(ccfDetailsSpecRegion$fracB2_mean)
      
      
      #additional features2
      #also get values for the expected BAF
      newBAFtable<- copyITHRDS$newBAFtable[[region]]$regionBAFtable 
      BAFsegment  <- newBAFtable[newBAFtable$chrom==chr&newBAFtable$pos>=startpos&newBAFtable$pos<=endpos,,drop=FALSE]
      
      checkGainSignificance <- function(logRBAFRegion,regionToCheck,ploidyEstimates,purityEstimates,gainThreshold)
      {
        minGainRawExpRegion       <-  ploidyEstimates[regionToCheck]*(2^(gainThreshold))
        minGainLogRExpRegion      <-  funCalcExpectedLogR(minGainRawExpRegion,0,ploidyEstimates[regionToCheck],purityEstimates[regionToCheck])
        logRObs                   <-  mean(logRBAFRegion$logR)
        
        if(nrow(logRBAFRegion)>=3)
        {
          Ttest        <- tTestWithErrorHandling(logRBAFRegion$logR,mu = minGainLogRExpRegion)
          TtestGreater <- tTestWithErrorHandling(logRBAFRegion$logR,mu = minGainLogRExpRegion,alternative = 'greater')
          #TtestLess    <- tTestWithErrorHandling(logRBAF$logR,mu = minLossLogRExp,alternative = 'less')
          logRObs_low  <- Ttest$conf.int[1]
          logRObs_high <- Ttest$conf.int[2]
          PVal         <- Ttest$p.value
          PvalGain     <- TtestGreater$p.value
          
        }
        if(nrow(logRBAFRegion)<3)
        {
          logRObs_low  <- NA
          logRObs_high <- NA
          PVal         <- NA
          PvalGain     <- NA
          
          
        }
        
        out <- cbind(minGainRawExpRegion,minGainLogRExpRegion,logRObs,logRObs_low,logRObs_high,PVal,PvalGain)  
        return(out)
        
        
      }
      checkLossSignificance <- function(logRBAFRegion,regionToCheck,ploidyEstimates,purityEstimates,lossThreshold)
      {
        minLossRawExpRegion       <-  ploidyEstimates[regionToCheck]*(2^(lossThreshold))
        minLossLogRExpRegion      <-  funCalcExpectedLogR(minLossRawExpRegion,0,ploidyEstimates[regionToCheck],purityEstimates[regionToCheck])
        logRObs                   <- mean(logRBAFRegion$logR)
        
        if(nrow(logRBAFRegion)>=3)
        {
          Ttest        <- tTestWithErrorHandling(logRBAFRegion$logR,mu = minLossLogRExpRegion)
          TtestLess <- tTestWithErrorHandling(logRBAFRegion$logR,mu = minLossLogRExpRegion,alternative = 'less')
          #TtestLess    <- tTestWithErrorHandling(logRBAF$logR,mu = minLossLogRExp,alternative = 'less')
          logRObs_low  <- Ttest$conf.int[1]
          logRObs_high <- Ttest$conf.int[2]
          PVal         <- Ttest$p.value
          PvalLoss     <- TtestLess$p.value
          
        }
        # if(nrow(logRBAF)<3)
        # {
        #   logRObs_low  <- NA
        #   logRObs_high <- NA
        #   PVal         <- NA
        #   PvalLoss     <- NA
        
        
        # }
        
        out <- cbind(minLossRawExpRegion,minLossLogRExpRegion,logRObs,logRObs_low,logRObs_high,PVal,PvalLoss)  
        return(out)
        
        
      }
      
      
      GainRegionOne     <- checkGainSignificance(BAFsegment,regionToCheck = region,ploidyEstimates,purityEstimates,gainThreshold = gainThreshold)
      AmpRegionOne      <- checkGainSignificance(BAFsegment,regionToCheck = region,ploidyEstimates,purityEstimates,gainThreshold = ampThreshold)
      LossRegionOne     <- checkLossSignificance(BAFsegment,regionToCheck = region,ploidyEstimates,purityEstimates,lossThreshold = lossThreshold)
      DeepLossRegionOne <- checkLossSignificance(BAFsegment,regionToCheck = region,ploidyEstimates,purityEstimates,lossThreshold = deepLossThreshold)
      
      
      # put together into One table
      segStatus <- c(GainRegionOne[,3:5,drop=FALSE]
                     ,GainRegionOne[,2,drop=FALSE],GainRegionOne[,7,drop=FALSE]
                     ,AmpRegionOne[,2,drop=FALSE],AmpRegionOne[,7,drop=FALSE]
                     ,LossRegionOne[,2,drop=FALSE],LossRegionOne[,7,drop=FALSE]
                     ,DeepLossRegionOne[,2,drop=FALSE],DeepLossRegionOne[,7,drop=FALSE]
      )
      names(segStatus) <- c('logRObs','logRObs_low','logRObs_high'
                            ,'minGainExpLogR','gainPval'
                            ,'minAmpExpLogR','ampPval'
                            ,'minLossExpLogR','lossPval'
                            ,'minDeepLossExpLogR','DeepLossPval')
      
      logRObs            <- median(as.numeric(BAFsegment$BAF))
      gainPval           <- segStatus["gainPval"]
      ampPval            <- segStatus["ampPval"]
      lossPval           <- segStatus["lossPval"]
      DeepLossPval       <- segStatus["DeepLossPval"]
      
      obsBAF          <-  median(BAFsegment$newBAFsegmented)
      
      
      # # B1 ttest
      if(length(unique(segDetailsSpecRegion$nBfrac))==1)
      {
        if(unique(segDetailsSpecRegion$nBfrac)==1)
        {
          TtestB_subclonal <- list()
          TtestB_subclonal$p.value <- 1
        }
        
      }
      # 
      if(length(unique(segDetailsSpecRegion$nBfrac))>1)
      {
        if(fracB1>fracB2)
        {
          TtestB_subclonal   <- wilcox.test(segDetailsSpecRegion$nB1frac,mu=mu,alternative = "less")
          
        }
        
        if(fracB2>fracB1)
        {
          TtestB_subclonal   <- wilcox.test(segDetailsSpecRegion$nB1frac,mu=mu,alternative = "less")
          
        }
        
      }
      
      if(fracB1>fracB2)
      {
        #choose nA1 as integer reference
        rawBs   <- segDetailsSpecRegion$nB
        intNumB <- nB1
        
        fullBs <- c()
        nrMatchedB_Chrs <- c()
        for (chr in chrsToConsider)
        {
          if(length(regionChrList[[region]][[as.character(intNumB)]][[as.character(chr)]]$nBvals)==0)
          {
            next;
          }
          if(length(regionChrList[[region]][[as.character(intNumB)]][[as.character(chr)]]$nBvals)!=0)
          {
            nrMatchedB_Chrs <- c(nrMatchedB_Chrs,as.character(chr))
            fullBs <- c(fullBs,regionChrList[[region]][[as.character(intNumB)]][[as.character(chr)]]$nBvals)
          }
          
        }
        TtestBs             <- tTestWithErrorHandling(rawBs,fullBs)
        TtestDiffFromIntB   <- tTestWithErrorHandling(rawBs,mu=intNumB)
        nrMatchedB_Chrs     <- length(nrMatchedB_Chrs)
        
      }
      
      
      if(fracB1<=fracB2)
      {
        #choose nA1 as integer reference
        rawBs   <- segDetailsSpecRegion$nB
        intNumB <- nB2
        
        fullBs <- c()
        nrMatchedB_Chrs <- c()
        for (chr in chrsToConsider)
        {
          if(length(regionChrList[[region]][[as.character(intNumB)]][[as.character(chr)]]$nBvals)==0)
          {
            next;
          }
          if(length(regionChrList[[region]][[as.character(intNumB)]][[as.character(chr)]]$nBvals)!=0)
          {
            nrMatchedB_Chrs <- c(nrMatchedB_Chrs,as.character(chr))
            fullBs <- c(fullBs,regionChrList[[region]][[as.character(intNumB)]][[as.character(chr)]]$nBvals)
          }
          
        }
        TtestBs             <- tTestWithErrorHandling(rawBs,fullBs)
        TtestDiffFromIntB   <- tTestWithErrorHandling(rawBs,mu=intNumB)
        nrMatchedB_Chrs     <- length(nrMatchedB_Chrs)
        
      }
      
      # # A1 ttest
      # if(length(unique(segDetailsSpecRegion$nBfrac))==1)
      # {
      #   if(unique(segDetailsSpecRegion$nAfrac)==1)
      #   {
      #     TtestA1_subclonal <- list()
      #     TtestA1_subclonal$p.value <- 1
      #   }
      #   
      # }
      # 
      # if(length(unique(segDetailsSpecRegion$nAfrac))>1)
      # {
      #   TtestA1_subclonal   <- wilcox.test(segDetailsSpecRegion$nAfrac,mu=mu,alternative = "less")
      # }
      # 
      # 
      # #A2 ttest
      # if(length(unique(ifelse(c(1-segDetailsSpecRegion$nAfrac)<0,0,c(1-segDetailsSpecRegion$nAfrac))))==1)
      # {
      #   if(unique(ifelse(c(1-segDetailsSpecRegion$nAfrac)<0,0,c(1-segDetailsSpecRegion$nAfrac)))%in%c(1,0))
      #   {
      #     TtestA2_subclonal <- list()
      #     TtestA2_subclonal$p.value <- 1
      #   }
      #   
      # }
      
      # determine whether any evidence for subclonal A or B
      subclonalA         <- NA
      subclonalB         <- TRUE
      heterogeneousA     <- NA
      heterogeneousB     <- TRUE
      
      # if(nA1==nA2)
      # {
      #   subclonalA         <- FALSE
      #   heterogeneousA     <- FALSE
      # }
     
      # if(fracA1>=1|fracA2>=1|TtestA1_subclonal$p.value>=subclonalPvalue|TtestA2_subclonal$p.value>=subclonalPvalue|TtestA1_subclonal$p.value>=subclonalPvalue)
      # {
      #   subclonalA <- FALSE
      # }
      smallerPval <- min(TtestB_subclonal$p.value)
      TtestDiffOtherChrs <- 1e-100
      
      if(nrMatchedB_Chrs>=minOtherChrs)
      {
        TtestDiffOtherChrs <- TtestBs$p.value
      }
      
      
      if(fracB1>=mu|fracB2>=mu|smallerPval>=clonalPvalue|TtestDiffOtherChrs>=diffFromOtherPvalToUse)
      {
        subclonalB <- FALSE
      }
      
      if(nB1==nB2)
      {
        subclonalB         <- FALSE
        heterogeneousB     <- FALSE
      }
      
         
      origAraw <- NA
      origBraw <- NA
      smallPvalA <- NA#smallerPval
      smallPvalB <- smallerPval
      A_OtherChrsWithState <- #nrMatchedA_Chrs
      A_OtherChrPval       <- #diffFromOtherPval
      B_OtherChrsWithState <- nrMatchedB_Chrs
      B_OtherChrPval       <- TtestDiffOtherChrs
      
      chr <- chrName
      rowOut <- cbind(sample,chr,startpos,endpos,n.het,cnTotal,nMajor,nMinor,Ploidy,ACF,nAraw,nBraw,psi,date,nA1,nA2,nB1,nB2,fracA1,fracA2,fracB1,fracB2,subclonalA,subclonalB,heterogeneousA,heterogeneousB,origAraw,origBraw,segName
                      ,logRObs
                      ,obsBAF
                      ,gainPval
                      ,ampPval
                      ,lossPval
                      ,DeepLossPval
                      ,smallPvalA
                      ,smallPvalB
                      ,A_OtherChrsWithState
                      ,A_OtherChrPval
                      ,B_OtherChrsWithState
                      ,B_OtherChrPval)
      
      
      # segOut   <- cbind(sample,chr,startpos,endpos,segName,nAraw)
      segMatCopy_nB_CCF <- rbind(segMatCopy_nB_CCF,rowOut)
    }
    
    # next, combine the two
    rownames(segMatCopy_nA_CCF) <- NULL -> rownames(segMatCopy_nB_CCF)
    segMatCopy_nA_CCF<- data.frame(segMatCopy_nA_CCF,stringsAsFactors=FALSE)
    segMatCopy_nB_CCF<- data.frame(segMatCopy_nB_CCF,stringsAsFactors=FALSE)
    rownames(segMatCopy_nA_CCF) <- unique(names(regionJoinedSegCCF[[region]][['segs_nA']]))
    rownames(segMatCopy_nB_CCF) <- unique(names(regionJoinedSegCCF[[region]][['segs_nB']]))
    
    segNamesToConsider <- unique(c(segMatCopy_nA_CCF$segName,segMatCopy_nB_CCF$segName))
    segNamesTable      <- do.call(rbind,strsplit(segNamesToConsider,split=":"))
    segNamesTable      <- segNamesTable[order(as.numeric(segNamesTable[,1]),as.numeric(segNamesTable[,2]),as.numeric(segNamesTable[,3])),,drop=FALSE]
    
    segNamesTable      <- data.frame(segNamesTable,stringsAsFactors = FALSE)
    
    segNamesTable$segName <- paste(segNamesTable[,1],segNamesTable[,2],segNamesTable[,3],sep=":")
    colnames(segNamesTable)[1:3] <- c('chr','startpos','endpos')
    rownames(segNamesTable) <- segNamesTable$segName
    
    segNamesTableNew <- segNamesTable
    # check which ones should be removed
    for (i in 2:nrow(segNamesTable))
    {
      if((as.numeric(segNamesTable[i,]$startpos)<as.numeric(segNamesTable[(i-1),]$endpos))&(segNamesTable[i,]$chr==segNamesTable[(i-1),]$chr)&i<=nrow(segNamesTable))
      {
        segNamesTable <- segNamesTable[-i,]
      }
      
    }
    
    
    # next, seg the appropriate values for nA and nB
    segNamesTable$sample  <- region
    segNamesTable$n.het   <- NA
    segNamesTable$cnTotal <- NA
    segNamesTable$cnMajor <- NA
    segNamesTable$cnMinor <- NA
    segNamesTable$Ploidy  <- ploidyEstimates[region]
    segNamesTable$ACF     <- purityEstimates[region]
    segNamesTable$nAraw <- NA
    segNamesTable$nBraw <- NA
    segNamesTable$psi <- NA
    segNamesTable$date <- dateToUse
    
    segNamesTable$nA1<- NA
    segNamesTable$nA2 <- NA
    segNamesTable$nB1 <- NA
    segNamesTable$nB2 <- NA
    
    
    segNamesTable$fracA1 <- NA
    segNamesTable$fracA2 <- NA
    segNamesTable$fracB1 <- NA
    segNamesTable$fracB2 <- NA
    
    segNamesTable$subclonalA <- NA
    segNamesTable$subclonalB <- NA
    segNamesTable$heterogeneousA <- NA
    segNamesTable$heterogeneousB <- NA
    
    segNamesTable$origAraw <- NA
    segNamesTable$origBraw <- NA
    
    #segNamesTable$segName <- NA
    segNamesTable$logRObs <- NA
    segNamesTable$obsBAF <- NA
    segNamesTable$gainPval <- NA
    segNamesTable$ampPval <- NA
    segNamesTable$lossPval <- NA
    segNamesTable$DeepLossPval <- NA
    segNamesTable$smallPvalA <- NA
    segNamesTable$smallPvalB <- NA
    segNamesTable$A_OtherChrsWithState <- NA
    segNamesTable$A_OtherChrPval <- NA
    segNamesTable$B_OtherChrsWithState <- NA
    segNamesTable$B_OtherChrPval <- NA
    
    
    
    
    for (segName in segNamesTable$segName)
    {
      #do we have this for segMatCopy_nA
      if(segName%in%segMatCopy_nA_CCF$segName)
      {
        segNamesTable[segName,]$nAraw <- segMatCopy_nA_CCF[segName,]$nAraw
        segNamesTable[segName,]$nA1 <- segMatCopy_nA_CCF[segName,]$nA1
        segNamesTable[segName,]$nA2 <- segMatCopy_nA_CCF[segName,]$nA2
        segNamesTable[segName,]$fracA1 <- segMatCopy_nA_CCF[segName,]$fracA1
        segNamesTable[segName,]$fracA2 <- segMatCopy_nA_CCF[segName,]$fracA2
        segNamesTable[segName,]$subclonalA <- segMatCopy_nA_CCF[segName,]$subclonalA
        segNamesTable[segName,]$heterogeneousA <- segMatCopy_nA_CCF[segName,]$heterogeneousA
        segNamesTable[segName,]$origAraw <- segMatCopy_nA_CCF[segName,]$origAraw
        segNamesTable[segName,]$n.het <- segMatCopy_nA_CCF[segName,]$n.het
        segNamesTable[segName,]$smallPvalA <- segMatCopy_nA_CCF[segName,]$smallPvalA
        segNamesTable[segName,]$A_OtherChrsWithState <- segMatCopy_nA_CCF[segName,]$A_OtherChrsWithState
        segNamesTable[segName,]$A_OtherChrPval <- segMatCopy_nA_CCF[segName,]$A_OtherChrPval
        
        segNamesTable[segName,]$logRObs <- segMatCopy_nA_CCF[segName,]$logRObs
        segNamesTable[segName,]$obsBAF <- segMatCopy_nA_CCF[segName,]$obsBAF
        segNamesTable[segName,]$gainPval <- segMatCopy_nA_CCF[segName,]$gainPval
        segNamesTable[segName,]$ampPval <- segMatCopy_nA_CCF[segName,]$ampPval
        segNamesTable[segName,]$lossPval <- segMatCopy_nA_CCF[segName,]$lossPval
        segNamesTable[segName,]$DeepLossPval <- segMatCopy_nA_CCF[segName,]$DeepLossPval
        
        
        
        
      }
      if(segName%in%segMatCopy_nB_CCF$segName)
      {
        segNamesTable[segName,]$nBraw <- segMatCopy_nB_CCF[segName,]$nBraw
        segNamesTable[segName,]$nB1 <- segMatCopy_nB_CCF[segName,]$nB1
        segNamesTable[segName,]$nB2 <- segMatCopy_nB_CCF[segName,]$nB2
        segNamesTable[segName,]$fracB1 <- segMatCopy_nB_CCF[segName,]$fracB1
        segNamesTable[segName,]$fracB2 <- segMatCopy_nB_CCF[segName,]$fracB2
        segNamesTable[segName,]$subclonalB <- segMatCopy_nB_CCF[segName,]$subclonalB
        segNamesTable[segName,]$heterogeneousB <- segMatCopy_nB_CCF[segName,]$heterogeneousB
        segNamesTable[segName,]$origBraw <- segMatCopy_nB_CCF[segName,]$origBraw
        segNamesTable[segName,]$n.het <- segMatCopy_nB_CCF[segName,]$n.het
        segNamesTable[segName,]$smallPvalB <- segMatCopy_nB_CCF[segName,]$smallPvalB
        
        segNamesTable[segName,]$B_OtherChrsWithState <- segMatCopy_nB_CCF[segName,]$B_OtherChrsWithState
        segNamesTable[segName,]$B_OtherChrPval <- segMatCopy_nB_CCF[segName,]$B_OtherChrPval
        
        segNamesTable[segName,]$logRObs <- segMatCopy_nB_CCF[segName,]$logRObs
        segNamesTable[segName,]$obsBAF <- segMatCopy_nB_CCF[segName,]$obsBAF
        segNamesTable[segName,]$gainPval <- segMatCopy_nB_CCF[segName,]$gainPval
        segNamesTable[segName,]$ampPval <- segMatCopy_nB_CCF[segName,]$ampPval
        segNamesTable[segName,]$lossPval <- segMatCopy_nB_CCF[segName,]$lossPval
        segNamesTable[segName,]$DeepLossPval <- segMatCopy_nB_CCF[segName,]$DeepLossPval
        
        
      }
      if(!segName%in%segMatCopy_nB_CCF$segName)
      {
        
        #segName
        segChr   <- as.numeric(unlist(strsplit(segName,split=":"))[1])
        segStart <- as.numeric(unlist(strsplit(segName,split=":"))[2])
        segStop  <- as.numeric(unlist(strsplit(segName,split=":"))[3])
        
        segMatCopy_nBsubset <- segMatCopy_nB_CCF[segMatCopy_nB_CCF$chr==segChr,,drop=FALSE]
        segMatCopy_nBsubset$startpos <- as.numeric(do.call(rbind,(strsplit(segMatCopy_nBsubset$segName,split=":")))[,2])
        segMatCopy_nBsubset$endpos <- as.numeric(do.call(rbind,(strsplit(segMatCopy_nBsubset$segName,split=":")))[,3])
        
        segNamesTable[segName,]$nBraw <- as.numeric(segMatCopy_nBsubset[as.numeric(segMatCopy_nBsubset$startpos)<=segStart&as.numeric(segMatCopy_nBsubset$endpos)>=segStart,]$nBraw)
        segNamesTable[segName,]$nB1 <- as.numeric(segMatCopy_nBsubset[as.numeric(segMatCopy_nBsubset$startpos)<=segStart&as.numeric(segMatCopy_nBsubset$endpos)>=segStart,]$nB1)
        segNamesTable[segName,]$nB2 <- as.numeric(segMatCopy_nBsubset[as.numeric(segMatCopy_nBsubset$startpos)<=segStart&as.numeric(segMatCopy_nBsubset$endpos)>=segStart,]$nB2)
        segNamesTable[segName,]$fracB1 <- as.numeric(segMatCopy_nBsubset[as.numeric(segMatCopy_nBsubset$startpos)<=segStart&as.numeric(segMatCopy_nBsubset$endpos)>=segStart,]$fracB1)
        segNamesTable[segName,]$fracB2 <- as.numeric(segMatCopy_nBsubset[as.numeric(segMatCopy_nBsubset$startpos)<=segStart&as.numeric(segMatCopy_nBsubset$endpos)>=segStart,]$fracB2)
        segNamesTable[segName,]$subclonalB <- segMatCopy_nBsubset[as.numeric(segMatCopy_nBsubset$startpos)<=segStart&as.numeric(segMatCopy_nBsubset$endpos)>=segStart,]$subclonalB
        segNamesTable[segName,]$heterogeneousB <- segMatCopy_nBsubset[as.numeric(segMatCopy_nBsubset$startpos)<=segStart&as.numeric(segMatCopy_nBsubset$endpos)>=segStart,]$heterogeneousB
        segNamesTable[segName,]$origBraw <- as.numeric(segMatCopy_nBsubset[as.numeric(segMatCopy_nBsubset$startpos)<=segStart&as.numeric(segMatCopy_nBsubset$endpos)>=segStart,]$origBraw)
        segNamesTable[segName,]$smallPvalB <- as.numeric(segMatCopy_nBsubset[as.numeric(segMatCopy_nBsubset$startpos)<=segStart&as.numeric(segMatCopy_nBsubset$endpos)>=segStart,]$smallPvalB)
        segNamesTable[segName,]$B_OtherChrsWithState <- as.numeric(segMatCopy_nBsubset[as.numeric(segMatCopy_nBsubset$startpos)<=segStart&as.numeric(segMatCopy_nBsubset$endpos)>=segStart,]$B_OtherChrsWithState)
        segNamesTable[segName,]$B_OtherChrPval  <- as.numeric(segMatCopy_nBsubset[as.numeric(segMatCopy_nBsubset$startpos)<=segStart&as.numeric(segMatCopy_nBsubset$endpos)>=segStart,]$B_OtherChrPval)
        
        
        
        
      }
      if(!segName%in%segMatCopy_nA_CCF$segName)
      {
        #segName
        segChr   <- as.numeric(unlist(strsplit(segName,split=":"))[1])
        segStart <- as.numeric(unlist(strsplit(segName,split=":"))[2])
        segStop  <- as.numeric(unlist(strsplit(segName,split=":"))[3])
        
        segMatCopy_nAsubset <- segMatCopy_nA_CCF[segMatCopy_nA_CCF$chr==segChr,,drop=FALSE]
        segMatCopy_nAsubset$startpos <- as.numeric(do.call(rbind,(strsplit(segMatCopy_nAsubset$segName,split=":")))[,2])
        segMatCopy_nAsubset$endpos <- as.numeric(do.call(rbind,(strsplit(segMatCopy_nAsubset$segName,split=":")))[,3])
        
        segNamesTable[segName,]$nAraw <- as.numeric(segMatCopy_nAsubset[as.numeric(segMatCopy_nAsubset$startpos)<=segStart&as.numeric(segMatCopy_nAsubset$endpos)>=segStart,]$nAraw)
        segNamesTable[segName,]$nA1 <- as.numeric(segMatCopy_nAsubset[as.numeric(segMatCopy_nAsubset$startpos)<=segStart&as.numeric(segMatCopy_nAsubset$endpos)>=segStart,]$nA1)
        segNamesTable[segName,]$nA2 <- as.numeric(segMatCopy_nAsubset[as.numeric(segMatCopy_nAsubset$startpos)<=segStart&as.numeric(segMatCopy_nAsubset$endpos)>=segStart,]$nA2)
        segNamesTable[segName,]$fracA1 <- as.numeric(segMatCopy_nAsubset[as.numeric(segMatCopy_nAsubset$startpos)<=segStart&as.numeric(segMatCopy_nAsubset$endpos)>=segStart,]$fracA1)
        segNamesTable[segName,]$fracA2 <- as.numeric(segMatCopy_nAsubset[as.numeric(segMatCopy_nAsubset$startpos)<=segStart&as.numeric(segMatCopy_nAsubset$endpos)>=segStart,]$fracA2)
        segNamesTable[segName,]$subclonalA <- segMatCopy_nAsubset[as.numeric(segMatCopy_nAsubset$startpos)<=segStart&as.numeric(segMatCopy_nAsubset$endpos)>=segStart,]$subclonalA
        segNamesTable[segName,]$heterogeneousA <- segMatCopy_nAsubset[as.numeric(segMatCopy_nAsubset$startpos)<=segStart&as.numeric(segMatCopy_nAsubset$endpos)>=segStart,]$heterogeneousA
        segNamesTable[segName,]$origAraw <- as.numeric(segMatCopy_nAsubset[as.numeric(segMatCopy_nAsubset$startpos)<=segStart&as.numeric(segMatCopy_nAsubset$endpos)>=segStart,]$origAraw)
        segNamesTable[segName,]$smallPvalA <- segMatCopy_nAsubset[as.numeric(segMatCopy_nAsubset$startpos)<=segStart&as.numeric(segMatCopy_nAsubset$endpos)>=segStart,]$smallPvalA
        segNamesTable[segName,]$A_OtherChrsWithState <- segMatCopy_nAsubset[as.numeric(segMatCopy_nAsubset$startpos)<=segStart&as.numeric(segMatCopy_nAsubset$endpos)>=segStart,]$A_OtherChrsWithState
        segNamesTable[segName,]$A_OtherChrPval <- as.numeric(segMatCopy_nAsubset[as.numeric(segMatCopy_nAsubset$startpos)<=segStart&as.numeric(segMatCopy_nAsubset$endpos)>=segStart,]$A_OtherChrPval)
        
        #segNamesTable[segName,]$nAraw <- segMatCopy_nA[do.call(rbind,strsplit(segMatCopy_nAB$segName,split=":"))[,1]==unlist(strsplit(segName,split=":"))[1]&as.numeric(do.call(rbind,strsplit(segMatCopy_nA$segName,split=":"))[,2])<=as.numeric(unlist(strsplit(segName,split=":"))[2])&as.numeric(do.call(rbind,strsplit(segMatCopy_nA$segName,split=":"))[,3])>=as.numeric(unlist(strsplit(segName,split=":"))[3]),]$nAraw
        
        
      }
      
      
    }
    segNamesTable$sample <- region
    segMatCopyCombinedCCF   <- rbind(segMatCopyCombinedCCF,segNamesTable)
    segMatCopy_nA_allCCF        <- rbind(segMatCopy_nA_allCCF,segMatCopy_nA_CCF)
    segMatCopy_nB_allCCF        <- rbind(segMatCopy_nB_allCCF,segMatCopy_nB_CCF)
  }
  
  
  toReturn <- list()
  toReturn$fullCopyTable         <-   fullTable
  toReturn$segMatCopy_nA_allCCF  <-   segMatCopy_nA_allCCF
  toReturn$segMatCopy_nB_allCCF  <-   segMatCopy_nB_allCCF
  toReturn$segMatCopyCombinedCCF <-   segMatCopyCombinedCCF
  toReturn$regionJoinedSegCCF    <-   regionJoinedSegCCF
  toReturn$PatientClonalDissectionAmpsAndDels <- PatientClonalDissectionAmpsAndDels
  return(toReturn)
  
}


#put together copy numberTable
getTableForCancerStudy <- function(releaseDir
                                   ,cancerCode="MET")
{
  allTumors <- grep(cancerCode,list.files(releaseDir),value=TRUE)
  fullSegTable  <- c()
  copyAsegTable <- c()
  copyBsegTable <- c()
  copyCombinedSegTable <- c()
  missingData <- c()
  
  for (tumor in allTumors)
  {
    allOutPut <- list.files(paste(releaseDir,tumor,"/exome/CN_heterogeneity/",sep=""))
    
    if(length(grep("_segmentation_output.RDS",allOutPut,value=TRUE))==0)
    {
      missingData <- c(missingData,tumor)
      next;
    }
    segListFile <- paste(releaseDir,tumor,"/exome/CN_heterogeneity/",grep("_segmentation_output.RDS",allOutPut,value=TRUE),sep="")
    
    if(!file.exists(segListFile))
    {
      missingData <- c(missingData,tumor)
      next;
    }
    # first load the segmentation_ouput
    segList   <- readRDS(segListFile)
    fullSegTable  <- rbind(fullSegTable,segList$fullCopyTable)
    copyAsegTable <- rbind(copyAsegTable,segList$segMatCopy_nA_allCCF)
    copyBsegTable <- rbind(copyBsegTable,segList$segMatCopy_nB_allCCF)
    copyCombinedSegTable <- rbind(copyCombinedSegTable,segList$segMatCopyCombinedCCF)
  }
  
  plot(apply(cbind(as.numeric(fullSegTable$nAraw),as.numeric(fullSegTable$nBraw)),1,max),as.numeric(fullSegTable$origAraw),pch=16,col='#bdbdbd99',xlim=c(0,200),ylim=c(0,200),xlab='new Araw',ylab='original Araw')
  abline(0,1)
  
  testTable <- fullSegTable[apply(cbind(as.numeric(fullSegTable$nAraw),as.numeric(fullSegTable$nBraw)),1,max)>20&as.numeric(fullSegTable$origAraw)<10,,drop=FALSE]
  
}

##############################################################################End of plotting script ##########################
##############################################################################################################################


getCCFractionNeeded <- function (obsValue,lowExpValue, highExpValue)
{
  if(obsValue>highExpValue)
  {
    if(obsValue/highExpValue!=Inf)
    {
      return(c(0,obsValue/highExpValue))
    }
    
    if(obsValue/highExpValue==Inf)
    {
      return(c(1,0))
    }
  }
  
  if(obsValue<lowExpValue)
  {
    return(c(1,0))
  }
  
  # because my brain isn't working, i'm going to do a brute force 2digit CCF calc
  CCFtoUse <- 0.01
  estValue <- ((CCFtoUse * lowExpValue))+((1-CCFtoUse) * highExpValue) 
  estValue <- abs(obsValue-estValue)
  
  for (CCFvalue in seq(0.0,1.5,by=.01))
  {
    CestValue <- ((CCFvalue * lowExpValue))+((1-CCFvalue) * highExpValue) 
    CestValue <- abs(obsValue-CestValue)
    if(CestValue<estValue)
    {
      estValue <- CestValue
      CCFtoUse <- CCFvalue
    }
  }
  
  
  return(c(CCFtoUse,(1-CCFtoUse)))
}

getSegCCFEstimate           <- function(copyITHRDS
                                        ,purityEstimates
                                        ,ploidyEstimates
                                        ,highestTrustedAmp=5
                                        ,pValThresh=0.01
                                        ,PvalToSeparate=0.001
                                        ,gainPValToUse = 'gainPval'
                                        ,lossPValToUse = 'lossPval'
                                        ,minCCFtoConsider=0.15
                                        ,plotting=FALSE
                                        ,confInt=0.99999
                                        #,plotHeterogeneity
)
{
  #need to obtain an ancestral state for each copy number 
  # (initially, simply do this by voting (will fail for subclonal GD events)).
  
  regions  <- names(purityEstimates)
  regions  <- grep(":",names(copyITHRDS$regionComparisonList),invert=TRUE,value=TRUE)
  segNames <- names(copyITHRDS$regionComparisonList[[names(ploidyEstimates)[1]]])
  for (region in regions)
  {
    segNames <- intersect(segNames,names(copyITHRDS$regionComparisonList[[region]]))
    
  }
  #print(regions)
  #print(segNames)
  # first, check for each segName if the rounded numbers are the same 
  # - if they are the same, then we can easily calculate the CCF for each region
  sameInteger <- c()
  diffInteger <- c()
  sameIntegerA <- c()
  diffIntegerA <- c()
  sameIntegerB <- c()
  diffIntegerB <- c()
  for (segName in segNames)
  {
    nA <- c()
    nB <- c()
    
    for (region in regions)
    {
      nA <- c(nA,round(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean))
      nB <- c(nB,round(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean))
    }
    
    if(length(unique(nA))==1&length(unique(nB))==1)
    {
      sameInteger <- c(sameInteger,segName)
    }
    
    if(length(unique(nA))!=1|length(unique(nB))!=1)
    {
      diffInteger <- c(diffInteger,segName)
      
    }
    if(length(unique(nA))==1)
    {
      sameIntegerA <- c(sameIntegerA,segName)
    }
    if(length(unique(nA))!=1)
    {
      diffIntegerA <- c(diffIntegerA,segName)
    }
    
    if(length(unique(nB))==1)
    {
      sameIntegerB <- c(sameIntegerB,segName)
    }
    if(length(unique(nB))!=1)
    {
      diffIntegerB <- c(diffIntegerB,segName)
    }
    
    
  }
  
  #print(sameInteger)
  
  #
  segCCFvalues  <- list()
  maxCPN        <- 5
  
  getIndividualCCFestimates <- function(segName
                                        ,regions
                                        ,copyITHRDS
                                        ,maxCPN=5
                                        ,minCCFtoConsider=0.15
                                        ,confInt=0.9999
                                        
  )
  {
    
    segCCFvalues <- list()
    
    strsplitSegName <- unlist(strsplit(segName,split=":"))
    chrName       <- strsplitSegName[1]
    postStartName <- strsplitSegName[2]
    postEndName   <- strsplitSegName[3]
    
    obsAvals <- c()
    obsBvals <- c()
    expAvals <- c()
    expBvals <- c()
    
    for (region in regions)
    {
      nAobs   <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean) #copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean
      nBobs   <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean)
      nAexp   <- max(0,round(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean))
      nBexp   <- max(0,round(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean))
      
      obsAvals <- c(obsAvals,nAobs)
      obsBvals <- c(obsBvals,nBobs)
      expAvals <- c(expAvals,nAexp)
      expBvals <- c(expBvals,nBexp)
    }
    
    names(obsAvals) <- names(obsBvals) <- names(expAvals) <- names(expBvals) <- regions
    
    
    # first, do this for Avals
    if(length(unique(obsAvals))==1)
    {
      # we can determine the low and high values easily 
      nA1  <- 1 # (this should really be ploidy)
      nA2  <- unique(obsAvals)
      
      segCCFvalues[[segName]]$nA1   <- nA1
      segCCFvalues[[segName]]$nA2   <- nA2
      
      # now go region by region, and check the fraction of each (which fraction to count as the CCF depends on whether it is a gain or a loss - determine this below)
      for (region in regions)
      {
        nAobs           <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean) #copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean
        nA1             <- nA1#floor(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean))
        nA2             <- nA2#ceiling(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean))
        
        fracsToUse      <- getCCFractionNeeded(nAobs,nA1,nA2)
        fracA1  <-       0
        fracA2  <-       1
        
        segCCFvalues[[segName]][[region]]$fracA1_mean  <- fracA1
        segCCFvalues[[segName]][[region]]$fracA2_mean  <- fracA2
        segCCFvalues[[segName]][[region]]$nAobs_mean   <- nAobs
        #segCCFvalues[[segName]][[region]]$CCF_mean     <- fracA1
        
        # lets also add the raw values (just in case)
        nAobs_all       <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nA<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nA) #copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean
        nAobs_all       <- ifelse(nAobs_all>=maxCPN,maxCPN,nAobs_all)
        nA1             <- nA1#floor(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean))
        nA2             <- nA2#ceiling(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean))
        fracsToUse     <- sapply(nBobs_all,getCCFractionNeeded,nA1,nA2)
        fracsToUse[fracsToUse=="Inf"] <- 1
        
        fracA1_all      <-apply(fracsToUse,2,max)
        # fracifelse(as.numeric(nB2)-as.numeric(nBobs_all)<0
        #  ,abs(as.numeric(nB2)-as.numeric(nBobs_all))+1
        #  ,as.numeric(nB2)-as.numeric(nBobs_all))
        fracA2_all       <-apply(fracsToUse,2,min)
        
        
        #tTestWithErrorHandling(fracA1_all,mu = 1)
        
        segCCFvalues[[segName]][[region]]$fracA1_all  <- fracA1_all
        segCCFvalues[[segName]][[region]]$fracA2_all  <- fracA2_all
        segCCFvalues[[segName]][[region]]$nAobs_all   <- nAobs_all
        #segCCFvalues[[segName]][[region]]$fracA2_all  <- fracA2
        
        
      }
      
      #next;
    }
    
    # let's check whether the ceiling and floor lead to only 2 options
    if(length(unique(c(ceiling(obsAvals),floor(obsAvals))))==2)
    {
      # first check whether all regions are clonal
      # check whether only one state when rounding
      if(length(unique(round(obsAvals)))==1)
      {
        
        # check whether all regions can be classified as clonal
        roundedIntger <- round(obsAvals)[1]
        allClonal <- TRUE
        nonClonalRegion <- c()
        for (region in regions)
        {
          Test <- tTestWithErrorHandling(copyITHRDS$regionComparisonList[[region]][[segName]]$nA,mu = roundedIntger,conf.level = confInt)
          if(Test$p.value<(1-confInt))
          {
            allClonal <- FALSE
            nonClonalRegion <- c(nonClonalRegion,region)
          }
        }
        
        
        ttestClonal <-TRUE
        nonClonalRegionBAF <- c()
        for (region in regions)
        {
          
          segSplit <- unlist(strsplit(segName,split=":"))
          test <- copyITHRDS$newBAFtable[[region]]$regionBAFtable[copyITHRDS$newBAFtable[[region]]$regionBAFtable$pos>=as.numeric(segSplit[2])&copyITHRDS$newBAFtable[[region]]$regionBAFtable$pos<=as.numeric(segSplit[3])&copyITHRDS$newBAFtable[[region]]$regionBAFtable$chrom==as.numeric(segSplit[1]),,drop=FALSE]
          nAraw   <- round(mean(copyITHRDS$regionComparisonList[[region]][[segName]]$nA))
          nBraw   <- round(mean(copyITHRDS$regionComparisonList[[region]][[segName]]$nB))
          
          #get the SD of the sample
          # SD <- sd(c(test[test$nA==1,'BAF'],1-test[test$nB==1,'BAF']))
          # Number <- length(c(test[test$nA==1,'BAF'],1-test[test$nB==1,'BAF']))
          # mu <- funCalcExpectedBAF(expBvals[region],expAvals[region],ploidyEstimates[region],purityEstimates[region])
          # #tTestWithErrorHandling(c(test[test$nA==1,'BAF'],1-test[test$nB==1,'BAF']),rnorm(n = Number,sd = 0.5,mean=mu))
          #tTestWithErrorHandling(c(test[test$nA==1,'BAF'],1-test[test$nB==1,'BAF']),mu = funCalcExpectedBAF(expBvals[region],expAvals[region],ploidyEstimates[region],purityEstimates[region]))
          tTestSpec <- tTestWithErrorHandling(c(test[test$nA==1,'BAF'],1-test[test$nB==1,'BAF']),mu = funCalcExpectedBAF(expBvals[region],expAvals[region],ploidyEstimates[region],purityEstimates[region]))
          #tTest <- c(tTest,tTestSpec$p.value)
          if(tTestSpec$p.value<(1-confInt))
          {
            ttestClonal <- FALSE
            nonClonalRegionBAF <- c(nonClonalRegionBAF,region)
          }
          
          
        }
        
        
        # finally, check what the CCF is of each region if you assume there is subclonality
        maxCCF <- 0
        nA1  <- min(unique(floor(obsAvals)))
        nA2  <- max(unique(ceiling(obsAvals)))
        
        # now go region by region, and check the fraction of each (which fraction to count as the CCF depends on whether it is a gain or a loss - determine this below)
        for (region in regions)
        {
          nAobs           <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean) #copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean
          nA1             <- nA1#floor(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean))
          nA2             <- nA2#ceiling(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean))
          
          fracsToUse      <- getCCFractionNeeded(nAobs,nA1,nA2)
          minCCF <- min(fracsToUse)
          if(minCCF>maxCCF)
          {
            maxCCF <- minCCF
          }
        }
        
        # next, determine whether to explore subclonality
        if(allClonal==TRUE|ttestClonal==TRUE|maxCCF<minCCFtoConsider)
        {
          allClonal <- TRUE
        }
        
        
        if(allClonal)
        {
          nA1  <- round(obsAvals)[1]
          nA2  <- round(obsAvals)[1]
          
          segCCFvalues[[segName]]$nA1   <- nA1
          segCCFvalues[[segName]]$nA2   <- nA2
          
          # now go region by region, and check the fraction of each (which fraction to count as the CCF depends on whether it is a gain or a loss - determine this below)
          for (region in regions)
          {
            nAobs           <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean) #copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean
            nA1             <- nA1#floor(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean))
            nA2             <- nA2#ceiling(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean))
            fracsToUse      <- getCCFractionNeeded(nAobs,nA1,nA2)
            
            fracA1  <- max(fracsToUse)
            # ifelse(as.numeric(nA2)-as.numeric(nAobs)<0
            #   ,abs(as.numeric(nA2)-as.numeric(nAobs))+1
            #   ,as.numeric(nA2)-as.numeric(nAobs))
            fracA2  <- min(fracsToUse)
            
            segCCFvalues[[segName]][[region]]$fracA1_mean  <- fracA1
            segCCFvalues[[segName]][[region]]$fracA2_mean  <- fracA2
            segCCFvalues[[segName]][[region]]$nAobs_mean   <- nAobs
            #segCCFvalues[[segName]][[region]]$CCF_mean     <- fracA1
            
            # lets also add the raw values (just in case)
            nAobs_all      <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nA<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nA) #copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean
            nAobs_all       <- ifelse(nAobs_all>=maxCPN,maxCPN,nAobs_all)
            nA1             <- nA1#floor(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean))
            nA2             <- nA2#ceiling(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean))
            fracsToUse     <- sapply(nAobs_all,getCCFractionNeeded,nA1,nA2)
            fracsToUse[fracsToUse=="Inf"] <- 1
            
            
            fracA1_all      <-apply(fracsToUse,2,max)
            # fracifelse(as.numeric(nA2)-as.numeric(nAobs_all)<0
            #  ,abs(as.numeric(nA2)-as.numeric(nAobs_all))+1
            #  ,as.numeric(nA2)-as.numeric(nAobs_all))
            fracA2_all       <-apply(fracsToUse,2,min)
            segCCFvalues[[segName]][[region]]$fracA1_all  <- fracA1_all
            segCCFvalues[[segName]][[region]]$fracA2_all  <- fracA2_all
            segCCFvalues[[segName]][[region]]$nAobs_all   <- nAobs_all
            #segCCFvalues[[segName]][[region]]$fracA2_all  <- fracA2
            
            
          }
        }
        
        if(!allClonal)
        {
          # we can determine the low and high values easily 
          nA1  <- min(unique(floor(obsAvals)))
          nA2  <- max(unique(ceiling(obsAvals)))
          
          segCCFvalues[[segName]]$nA1   <- nA1
          segCCFvalues[[segName]]$nA2   <- nA2
          
          # now go region by region, and check the fraction of each (which fraction to count as the CCF depends on whether it is a gain or a loss - determine this below)
          for (region in regions)
          {
            nAobs           <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean) #copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean
            nA1             <- nA1#floor(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean))
            nA2             <- nA2#ceiling(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean))
            
            fracsToUse      <- getCCFractionNeeded(nAobs,nA1,nA2)
            fracA1  <-       fracsToUse[1]#ifelse(as.numeric(nA2)-as.numeric(nAobs)<0
            # ,abs(as.numeric(nA2)-as.numeric(nAobs))+1
            # ,as.numeric(nA2)-as.numeric(nAobs))
            fracA2  <- fracsToUse[2]#ifelse(1-fracA1<0,0,1-fracA1)
            
            segCCFvalues[[segName]][[region]]$fracA1_mean  <- fracA1
            segCCFvalues[[segName]][[region]]$fracA2_mean  <- fracA2
            segCCFvalues[[segName]][[region]]$nAobs_mean   <- nAobs
            #segCCFvalues[[segName]][[region]]$CCF_mean     <- fracA1
            
            # lets also add the raw values (just in case)
            nAobs_all      <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nA<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nA) #copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean
            nAobs_all       <- ifelse(nAobs_all>=maxCPN,maxCPN,nAobs_all)
            nA1             <- nA1#floor(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean))
            nA2             <- nA2#ceiling(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean))
            fracsToUse     <- sapply(nAobs_all,getCCFractionNeeded,nA1,nA2)
            fracsToUse[fracsToUse=="Inf"] <- 1
            
            
            fracA1_all      <-fracsToUse[1,] 
            # fracifelse(as.numeric(nA2)-as.numeric(nAobs_all)<0
            #  ,abs(as.numeric(nA2)-as.numeric(nAobs_all))+1
            #  ,as.numeric(nA2)-as.numeric(nAobs_all))
            fracA2_all       <-fracsToUse[2,] 
            
            #tTestWithErrorHandling(fracA1_all,mu = 1)
            
            segCCFvalues[[segName]][[region]]$fracA1_all  <- fracA1_all
            segCCFvalues[[segName]][[region]]$fracA2_all  <- fracA2_all
            segCCFvalues[[segName]][[region]]$nAobs_all   <- nAobs_all
            #segCCFvalues[[segName]][[region]]$fracA2_all  <- fracA2
            
            
          }
        }
        
        
        
      }
      
      if(length(unique(round(obsAvals)))!=1)
      {
        # we can determine the low and high values easily 
        nA1  <- min(unique(floor(obsAvals)))
        nA2  <- max(unique(ceiling(obsAvals)))
        
        segCCFvalues[[segName]]$nA1   <- nA1
        segCCFvalues[[segName]]$nA2   <- nA2
        
        # now go region by region, and check the fraction of each (which fraction to count as the CCF depends on whether it is a gain or a loss - determine this below)
        for (region in regions)
        {
          nAobs           <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean) #copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean
          nA1             <- nA1#floor(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean))
          nA2             <- nA2#ceiling(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean))
          
          fracsToUse      <- getCCFractionNeeded(nAobs,nA1,nA2)
          fracA1  <-       fracsToUse[1]#ifelse(as.numeric(nA2)-as.numeric(nAobs)<0
          # ,abs(as.numeric(nA2)-as.numeric(nAobs))+1
          # ,as.numeric(nA2)-as.numeric(nAobs))
          fracA2  <- fracsToUse[2]#ifelse(1-fracA1<0,0,1-fracA1)
          
          segCCFvalues[[segName]][[region]]$fracA1_mean  <- fracA1
          segCCFvalues[[segName]][[region]]$fracA2_mean  <- fracA2
          segCCFvalues[[segName]][[region]]$nAobs_mean   <- nAobs
          #segCCFvalues[[segName]][[region]]$CCF_mean     <- fracA1
          
          # lets also add the raw values (just in case)
          nAobs_all      <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nA<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nA) #copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean
          nAobs_all       <- ifelse(nAobs_all>=maxCPN,maxCPN,nAobs_all)
          nA1             <- nA1#floor(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean))
          nA2             <- nA2#ceiling(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean))
          fracsToUse     <- sapply(nAobs_all,getCCFractionNeeded,nA1,nA2)
          fracsToUse[fracsToUse=="Inf"] <- 1
          
          
          fracA1_all      <-fracsToUse[1,] 
          # fracifelse(as.numeric(nA2)-as.numeric(nAobs_all)<0
          #  ,abs(as.numeric(nA2)-as.numeric(nAobs_all))+1
          #  ,as.numeric(nA2)-as.numeric(nAobs_all))
          fracA2_all       <-fracsToUse[2,] 
          
          #tTestWithErrorHandling(fracA1_all,mu = 1)
          
          segCCFvalues[[segName]][[region]]$fracA1_all  <- fracA1_all
          segCCFvalues[[segName]][[region]]$fracA2_all  <- fracA2_all
          segCCFvalues[[segName]][[region]]$nAobs_all   <- nAobs_all
          #segCCFvalues[[segName]][[region]]$fracA2_all  <- fracA2
          
          
        }
      }
      #next;
    }
    
    # let's check whether the ceiling and floor lead to only 2 options
    if(length(unique(c(ceiling(obsAvals),floor(obsAvals))))>2)
    {
      # check whether only one state when rounding
      if(length(unique(round(obsAvals)))==1)
      {
        
        # check whether all regions can be classified as clonal
        roundedIntger <- round(obsAvals)[1]
        allClonal <- TRUE
        nonClonalRegion <- c()
        for (region in regions)
        {
          Test <- tTestWithErrorHandling(copyITHRDS$regionComparisonList[[region]][[segName]]$nA,mu = roundedIntger,conf.level = confInt)
          if(Test$p.value<(1-confInt))
          {
            allClonal <- FALSE
            nonClonalRegion <- c(nonClonalRegion,region)
          }
        }
        
        ttestClonal <-TRUE
        nonClonalRegionBAF <- c()
        for (region in regions)
        {
          
          segSplit <- unlist(strsplit(segName,split=":"))
          test <- copyITHRDS$newBAFtable[[region]]$regionBAFtable[copyITHRDS$newBAFtable[[region]]$regionBAFtable$pos>=as.numeric(segSplit[2])&copyITHRDS$newBAFtable[[region]]$regionBAFtable$pos<=as.numeric(segSplit[3])&copyITHRDS$newBAFtable[[region]]$regionBAFtable$chrom==as.numeric(segSplit[1]),,drop=FALSE]
          nAraw   <- round(mean(copyITHRDS$regionComparisonList[[region]][[segName]]$nA))
          nBraw   <- round(mean(copyITHRDS$regionComparisonList[[region]][[segName]]$nB))
          
          tTestSpec <- tTestWithErrorHandling(c(test[test$nA==1,'BAF'],1-test[test$nB==1,'BAF']),mu = funCalcExpectedBAF(expBvals[region],expAvals[region],ploidyEstimates[region],purityEstimates[region]))
          #tTest <- c(tTest,tTestSpec$p.value)
          if(tTestSpec$p.value<(1-confInt))
          {
            ttestClonal <- FALSE
            nonClonalRegionBAF <- c(nonClonalRegionBAF,region)
          }
          
          
        }
        
        
        # finally, check what the CCF is of each region if you assume there is subclonality
        maxCCF <- 0
        nA1  <- min(unique(floor(obsAvals)))
        nA2  <- max(unique(round(obsAvals)))
        if(abs(min(max(ceiling(obsAvals))-obsAvals))<abs(max(min(floor(obsAvals))-obsAvals)))
        {
          nA1  <- round(obsAvals)[1]
          nA2  <- max(unique(ceiling(obsAvals)))
        }
        
        if(abs(min(max(ceiling(obsAvals))-obsAvals))>=abs(max(min(floor(obsAvals))-obsAvals)))
        {
          nA1  <- min(unique(floor(obsAvals)))
          nA2  <- max(unique(round(obsAvals)))
        }
        # now go region by region, and check the fraction of each (which fraction to count as the CCF depends on whether it is a gain or a loss - determine this below)
        for (region in regions)
        {
          nAobs           <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean) #copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean
          nA1             <- nA1#floor(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean))
          nA2             <- nA2#ceiling(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean))
          
          fracsToUse      <- getCCFractionNeeded(nAobs,nA1,nA2)
          minCCF <- min(fracsToUse)
          if(minCCF>maxCCF)
          {
            maxCCF <- minCCF
          }
        }
        
        # next, determine whether to explore subclonality
        if(allClonal==TRUE|ttestClonal==TRUE|maxCCF<minCCFtoConsider)
        {
          allClonal <- TRUE
        }
        
        
        if(allClonal)
        {
          nA1  <- round(obsAvals)[1]
          nA2  <- round(obsAvals)[1]
          
          segCCFvalues[[segName]]$nA1   <- nA1
          segCCFvalues[[segName]]$nA2   <- nA2
          
          # now go region by region, and check the fraction of each (which fraction to count as the CCF depends on whether it is a gain or a loss - determine this below)
          for (region in regions)
          {
            nAobs           <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean) #copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean
            nA1             <- nA1#floor(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean))
            nA2             <- nA2#ceiling(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean))
            fracsToUse      <- getCCFractionNeeded(nAobs,nA1,nA2)
            
            fracA1  <- max(fracsToUse)
            # ifelse(as.numeric(nA2)-as.numeric(nAobs)<0
            #   ,abs(as.numeric(nA2)-as.numeric(nAobs))+1
            #   ,as.numeric(nA2)-as.numeric(nAobs))
            fracA2  <- min(fracsToUse)
            
            segCCFvalues[[segName]][[region]]$fracA1_mean  <- fracA1
            segCCFvalues[[segName]][[region]]$fracA2_mean  <- fracA2
            segCCFvalues[[segName]][[region]]$nAobs_mean   <- nAobs
            #segCCFvalues[[segName]][[region]]$CCF_mean     <- fracA1
            
            # lets also add the raw values (just in case)
            nAobs_all      <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nA<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nA) #copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean
            nAobs_all       <- ifelse(nAobs_all>=maxCPN,maxCPN,nAobs_all)
            nA1             <- nA1#floor(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean))
            nA2             <- nA2#ceiling(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean))
            fracsToUse     <- sapply(nAobs_all,getCCFractionNeeded,nA1,nA2)
            fracsToUse[fracsToUse=="Inf"] <- 1
            
            
            fracA1_all      <-apply(fracsToUse,2,max)
            # fracifelse(as.numeric(nA2)-as.numeric(nAobs_all)<0
            #  ,abs(as.numeric(nA2)-as.numeric(nAobs_all))+1
            #  ,as.numeric(nA2)-as.numeric(nAobs_all))
            fracA2_all       <-apply(fracsToUse,2,min)
            segCCFvalues[[segName]][[region]]$fracA1_all  <- fracA1_all
            segCCFvalues[[segName]][[region]]$fracA2_all  <- fracA2_all
            segCCFvalues[[segName]][[region]]$nAobs_all   <- nAobs_all
            #segCCFvalues[[segName]][[region]]$fracA2_all  <- fracA2
            
            
          }
        }
        
        if(!allClonal)
        {
          # do we twist up or down
          if(abs(min(max(ceiling(obsAvals))-obsAvals))<abs(max(min(floor(obsAvals))-obsAvals)))
          {
            # take the ceiling 
            # we can determine the low and high values easily 
            nA1  <- round(obsAvals)[1]
            nA2  <- max(unique(ceiling(obsAvals)))
            
            segCCFvalues[[segName]]$nA1   <- nA1
            segCCFvalues[[segName]]$nA2   <- nA2
            
            # now go region by region, and check the fraction of each (which fraction to count as the CCF depends on whether it is a gain or a loss - determine this below)
            for (region in regions)
            {
              nAobs           <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean) #copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean
              nA1             <- nA1#floor(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean))
              nA2             <- nA2#ceiling(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean))
              fracsToUse      <- getCCFractionNeeded(nAobs,nA1,nA2)
              
              fracA1  <- fracsToUse[1]
              # ifelse(as.numeric(nA2)-as.numeric(nAobs)<0
              #   ,abs(as.numeric(nA2)-as.numeric(nAobs))+1
              #   ,as.numeric(nA2)-as.numeric(nAobs))
              fracA2  <- fracsToUse[2]
              
              segCCFvalues[[segName]][[region]]$fracA1_mean  <- fracA1
              segCCFvalues[[segName]][[region]]$fracA2_mean  <- fracA2
              segCCFvalues[[segName]][[region]]$nAobs_mean   <- nAobs
              #segCCFvalues[[segName]][[region]]$CCF_mean     <- fracA1
              
              # lets also add the raw values (just in case)
              nAobs_all      <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nA<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nA) #copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean
              nAobs_all       <- ifelse(nAobs_all>=maxCPN,maxCPN,nAobs_all)
              nA1             <- nA1#floor(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean))
              nA2             <- nA2#ceiling(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean))
              fracsToUse     <- sapply(nAobs_all,getCCFractionNeeded,nA1,nA2)
              fracsToUse[fracsToUse=="Inf"] <- 1
              
              
              fracA1_all      <-fracsToUse[1,] 
              # fracifelse(as.numeric(nA2)-as.numeric(nAobs_all)<0
              #  ,abs(as.numeric(nA2)-as.numeric(nAobs_all))+1
              #  ,as.numeric(nA2)-as.numeric(nAobs_all))
              fracA2_all       <-fracsToUse[2,] 
              segCCFvalues[[segName]][[region]]$fracA1_all  <- fracA1_all
              segCCFvalues[[segName]][[region]]$fracA2_all  <- fracA2_all
              segCCFvalues[[segName]][[region]]$nAobs_all   <- nAobs_all
              #segCCFvalues[[segName]][[region]]$fracA2_all  <- fracA2
              
              
            }
            
            #next;
            
          }
          if(abs(min(max(ceiling(obsAvals))-obsAvals))>=abs(max(min(floor(obsAvals))-obsAvals)))
          {
            # take the floor 
            # we can determine the low and high values easily 
            nA1  <- min(unique(floor(obsAvals)))
            nA2  <- max(unique(round(obsAvals)))
            
            segCCFvalues[[segName]]$nA1   <- nA1
            segCCFvalues[[segName]]$nA2   <- nA2
            
            # now go region by region, and check the fraction of each (which fraction to count as the CCF depends on whether it is a gain or a loss - determine this below)
            for (region in regions)
            {
              nAobs           <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean) #copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean
              nA1             <- nA1#floor(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean))
              nA2             <- nA2#ceiling(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean))
              fracsToUse      <- getCCFractionNeeded(nAobs,nA1,nA2)
              
              fracA1  <- fracsToUse[1]
              # ifelse(as.numeric(nA2)-as.numeric(nAobs)<0
              #   ,abs(as.numeric(nA2)-as.numeric(nAobs))+1
              #   ,as.numeric(nA2)-as.numeric(nAobs))
              fracA2  <- fracsToUse[2]
              #ifelse(1-fracA1<0,0,1-fracA1)
              
              segCCFvalues[[segName]][[region]]$fracA1_mean  <- fracA1
              segCCFvalues[[segName]][[region]]$fracA2_mean  <- fracA2
              segCCFvalues[[segName]][[region]]$nAobs_mean   <- nAobs
              #segCCFvalues[[segName]][[region]]$CCF_mean     <- fracA1
              
              # lets also add the raw values (just in case)
              nAobs_all      <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nA<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nA) #copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean
              nAobs_all       <- ifelse(nAobs_all>=maxCPN,maxCPN,nAobs_all)
              nA1             <- nA1#floor(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean))
              nA2             <- nA2#ceiling(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean))
              fracsToUse     <- sapply(nAobs_all,getCCFractionNeeded,nA1,nA2)
              fracsToUse[fracsToUse=="Inf"] <- 1
              
              
              fracA1_all      <-fracsToUse[1,] 
              # fracifelse(as.numeric(nA2)-as.numeric(nAobs_all)<0
              #  ,abs(as.numeric(nA2)-as.numeric(nAobs_all))+1
              #  ,as.numeric(nA2)-as.numeric(nAobs_all))
              fracA2_all       <-fracsToUse[2,] 
              
              segCCFvalues[[segName]][[region]]$fracA1_all  <- fracA1_all
              segCCFvalues[[segName]][[region]]$fracA2_all  <- fracA2_all
              segCCFvalues[[segName]][[region]]$nAobs_all   <- nAobs_all
              #segCCFvalues[[segName]][[region]]$fracA2_all  <- fracA2
              
              
            }
            
            #next;
            
          }
        }
        
        
        
      }
      
      #check whether there are only two states when rounding
      if(length(unique(round(obsAvals)))==2)
      {
        
        # we can determine the low and high values easily
        nA1  <- min(round(obsAvals))
        nA2  <- max(round(obsAvals))
        
        segCCFvalues[[segName]]$nA1   <- nA1
        segCCFvalues[[segName]]$nA2   <- nA2
        
        # now go region by region, and check the fraction of each (which fraction to count as the CCF depends on whether it is a gain or a loss - determine this below)
        for (region in regions)
        {
          nAobs           <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean) #copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean
          nA1             <- nA1#floor(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean))
          nA2             <- nA2#ceiling(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean))
          fracsToUse      <- getCCFractionNeeded(nAobs,nA1,nA2)
          fracA1  <- fracsToUse[1] 
          # ifelse(as.numeric(nA2)-as.numeric(nAobs)<0
          # ,abs(as.numeric(nA2)-as.numeric(nAobs))+1
          # ,as.numeric(nA2)-as.numeric(nAobs))
          fracA2  <- fracsToUse[2]
          # ifelse(1-fracA1<0,0,1-fracA1)
          
          segCCFvalues[[segName]][[region]]$fracA1_mean  <- fracA1
          segCCFvalues[[segName]][[region]]$fracA2_mean  <- fracA2
          segCCFvalues[[segName]][[region]]$nAobs_mean   <- nAobs
          #segCCFvalues[[segName]][[region]]$CCF_mean     <- fracA1
          
          # lets also add the raw values (just in case)
          nAobs_all      <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nA<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nA) #copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean
          nAobs_all       <- ifelse(nAobs_all>=maxCPN,maxCPN,nAobs_all)
          nA1             <- nA1#floor(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean))
          nA2             <- nA2#ceiling(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean))
          fracsToUse     <- sapply(nAobs_all,getCCFractionNeeded,nA1,nA2)
          fracsToUse[fracsToUse=="Inf"] <- 1
          
          
          fracA1_all      <-fracsToUse[1,] 
          # fracifelse(as.numeric(nA2)-as.numeric(nAobs_all)<0
          #  ,abs(as.numeric(nA2)-as.numeric(nAobs_all))+1
          #  ,as.numeric(nA2)-as.numeric(nAobs_all))
          fracA2_all       <-fracsToUse[2,] 
          # ifelse(1-fracA1_all<0,0,1-fracA1_all)
          
          segCCFvalues[[segName]][[region]]$fracA1_all  <- fracA1_all
          segCCFvalues[[segName]][[region]]$fracA2_all  <- fracA2_all
          segCCFvalues[[segName]][[region]]$nAobs_all   <- nAobs_all
          #segCCFvalues[[segName]][[region]]$fracA2_all  <- fracA2
          
          
        }
        
        #next;
      }
      
      
      
      
      # check whether an amplification greater than 8
      
      #if(min(expAvals))
      
      # posValues  <- unique(sort(expAvals))
      # lowValue   <- posValues[1]
      # highValue  <- posValues[2]
      # 
      # # choose the max and min that maximise the CCF across the tumour regions (in other words, minimize heterogeneity where possible)
      # posFracA1 <- c()
      # posFracA2 <- c()
      # valuesChosen <- c()
      # for (x in 1:(length(posValues)-1))
      # {
      #   for (y in (x+1):length(posValues))
      #   {
      #     
      #     posLowValue  <- posValues[x]
      #     posHighValue <- posValues[y]
      #     
      #     CCFvalues <- c()
      #     for (i in 1:length(obsAvals))
      #     {
      #       CCFvalues <-c(CCFvalues,(getCCFractionNeeded(obsAvals[i],posLowValue,posHighValue))[2])
      #     }
      #     valuesChosen <- rbind(valuesChosen,c(posLowValue,posHighValue))
      #     posFracA1 <- rbind(posFracA1,CCFvalues)
      #   }
      # }
      
      if(length(unique(round(obsAvals)))>2)
      {
        # we can determine the low and high values easily 
        nA1  <- min(unique(floor(obsAvals+minCCFtoConsider)))
        nA2  <- max(unique(ceiling(obsAvals-minCCFtoConsider)))
        
        segCCFvalues[[segName]]$nA1   <- nA1
        segCCFvalues[[segName]]$nA2   <- nA2
        
        # now go region by region, and check the fraction of each (which fraction to count as the CCF depends on whether it is a gain or a loss - determine this below)
        for (region in regions)
        {
          nAobs           <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean) #copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean
          nA1             <- nA1#floor(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean))
          nA2             <- nA2#ceiling(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean))
          
          fracsToUse      <- getCCFractionNeeded(nAobs,nA1,nA2)
          fracA1  <- fracsToUse[1] 
          # ifelse(as.numeric(nA2)-as.numeric(nAobs)<0
          # ,abs(as.numeric(nA2)-as.numeric(nAobs))+1
          # ,as.numeric(nA2)-as.numeric(nAobs))
          fracA2  <- fracsToUse[2]
          # ifelse(1-fracA1<0,0,1-fracA1)
          
          segCCFvalues[[segName]][[region]]$fracA1_mean  <- fracA1
          segCCFvalues[[segName]][[region]]$fracA2_mean  <- fracA2
          segCCFvalues[[segName]][[region]]$nAobs_mean   <- nAobs
          #segCCFvalues[[segName]][[region]]$CCF_mean     <- fracA1
          
          # lets also add the raw values (just in case)
          nAobs_all      <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nA<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nA) #copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean
          nAobs_all       <- ifelse(nAobs_all>=maxCPN,maxCPN,nAobs_all)
          nA1             <- nA1#floor(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean))
          nA2             <- nA2#ceiling(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean))
          fracsToUse     <- sapply(nAobs_all,getCCFractionNeeded,nA1,nA2)
          fracsToUse[fracsToUse=="Inf"] <- 1
          
          fracA1_all      <-as.numeric(fracsToUse[1,])
          # fracifelse(as.numeric(nB2)-as.numeric(nBobs_all)<0
          #  ,abs(as.numeric(nB2)-as.numeric(nBobs_all))+1
          #  ,as.numeric(nB2)-as.numeric(nBobs_all))
          fracA2_all      <-as.numeric(fracsToUse[2,])
          
          segCCFvalues[[segName]][[region]]$fracA1_all  <- fracA1_all
          segCCFvalues[[segName]][[region]]$fracA2_all  <- fracA2_all
          segCCFvalues[[segName]][[region]]$nAobs_all   <- nAobs_all
          #segCCFvalues[[segName]][[region]]$fracA2_all  <- fracA2
          
          
        }
      }
      
      
      
      #next;
    }
    
    # next, do this for Bvals
    if(length(unique(obsBvals))==1)
    {
      # we can determine the low and high values easily 
      nB1  <- 1 # (this should really be ploidy)
      nB2  <- unique(obsBvals)
      
      segCCFvalues[[segName]]$nB1   <- nB1
      segCCFvalues[[segName]]$nB2   <- nB2
      
      # now go region by region, and check the fraction of each (which fraction to count as the CCF depends on whether it is a gain or a loss - determine this below)
      for (region in regions)
      {
        nBobs           <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean) #copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean
        nB1             <- nB1#floor(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean))
        nB2             <- nB2#ceiling(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean))
        
        fracsToUse      <- getCCFractionNeeded(nBobs,nB1,nB2)
        fracB1  <-       0
        fracB2  <-       1
        
        segCCFvalues[[segName]][[region]]$fracB1_mean  <- fracB1
        segCCFvalues[[segName]][[region]]$fracB2_mean  <- fracB2
        segCCFvalues[[segName]][[region]]$nBobs_mean   <- nBobs
        #segCCFvalues[[segName]][[region]]$CCF_mean     <- fracB1
        
        # lets also add the raw values (just in case)
        nBobs_all      <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nB<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nB) #copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean
        nBobs_all       <- ifelse(nBobs_all>=maxCPN,maxCPN,nBobs_all)
        
        nB1             <- nB1#floor(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean))
        nB2             <- nB2#ceiling(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean))
        fracsToUse     <- sapply(nBobs_all,getCCFractionNeeded,nB1,nB2)
        fracsToUse[fracsToUse=="Inf"] <- 1
        
        fracB1_all      <-apply(fracsToUse,2,max)
        # fracifelse(as.numeric(nB2)-as.numeric(nBobs_all)<0
        #  ,abs(as.numeric(nB2)-as.numeric(nBobs_all))+1
        #  ,as.numeric(nB2)-as.numeric(nBobs_all))
        fracB2_all       <-apply(fracsToUse,2,min)
        
        # tTestWithErrorHandling(fracB1_all,mu = 1)
        
        segCCFvalues[[segName]][[region]]$fracB1_all  <- fracB1_all
        segCCFvalues[[segName]][[region]]$fracB2_all  <- fracB2_all
        segCCFvalues[[segName]][[region]]$nBobs_all   <- nBobs_all
        #segCCFvalues[[segName]][[region]]$fracB2_all  <- fracB2
        
        
      }
      
      # next;
    }
    
    # let's check whether the ceiling and floor lead to only 2 options
    if(length(unique(c(ceiling(obsBvals),floor(obsBvals))))==2)
    {
      
      if(length(unique(round(obsBvals)))==1)
      {
        
        # check whether all regions can be classified as clonal
        roundedIntger <- round(obsBvals)[1]
        allClonal <- TRUE
        nonClonalRegion <- c()
        for (region in regions)
        {
          Test <- tTestWithErrorHandling(copyITHRDS$regionComparisonList[[region]][[segName]]$nB,mu = roundedIntger,conf.level = confInt)
          if(Test$p.value<(1-confInt))
          {
            allClonal <- FALSE
            nonClonalRegion <- c(nonClonalRegion,region)
          }
        }
        
        
        # also, check whether significant different BAF distribution
        
        ttestClonal <-TRUE
        nonClonalRegionBAF <- c()
        for (region in regions)
        {
          
          segSplit <- unlist(strsplit(segName,split=":"))
          test <- copyITHRDS$newBAFtable[[region]]$regionBAFtable[copyITHRDS$newBAFtable[[region]]$regionBAFtable$pos>=as.numeric(segSplit[2])&copyITHRDS$newBAFtable[[region]]$regionBAFtable$pos<=as.numeric(segSplit[3])&copyITHRDS$newBAFtable[[region]]$regionBAFtable$chrom==as.numeric(segSplit[1]),,drop=FALSE]
          #nAraw   <- round(mean(copyITHRDS$regionComparisonList[[region]][[segName]]$nA))
          #nBraw   <- round(mean(copyITHRDS$regionComparisonList[[region]][[segName]]$nB))
          
          tTestSpec <- tTestWithErrorHandling(c(test[test$nA==1,'BAF'],1-test[test$nB==1,'BAF']),mu = funCalcExpectedBAF(expBvals[region],expAvals[region],ploidyEstimates[region],purityEstimates[region]))
          #tTest <- c(tTest,tTestSpec$p.value)
          if(tTestSpec$p.value<(1-confInt))
          {
            ttestClonal <- FALSE
            nonClonalRegionBAF <- c(nonClonalRegionBAF,region)
          }
          
          
        }
        
        
        # finally, check what the CCF is of each region if you assume there is subclonality
        maxCCF <- 0
        nB1  <- min(unique(floor(obsBvals)))
        nB2  <- max(unique(ceiling(obsBvals)))
        
        # now go region by region, and check the fraction of each (which fraction to count as the CCF depends on whether it is a gain or a loss - determine this below)
        for (region in regions)
        {
          nBobs           <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean) #copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean
          nB1             <- nB1#floor(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean))
          nB2             <- nB2#ceiling(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean))
          
          fracsToUse      <- getCCFractionNeeded(nBobs,nB1,nB2)
          minCCF <- min(fracsToUse)
          if(minCCF>maxCCF)
          {
            maxCCF <- minCCF
          }
        }
        
        # next, determine whether to explore subclonality
        if(allClonal==TRUE|ttestClonal==TRUE|maxCCF<minCCFtoConsider)
        {
          allClonal <- TRUE
        }
        
        
        
        if(allClonal)
        {
          nB1  <- round(obsBvals)[1]
          nB2  <- round(obsBvals)[1]
          
          segCCFvalues[[segName]]$nB1   <- nB1
          segCCFvalues[[segName]]$nB2   <- nB2
          
          # now go region by region, and check the fraction of each (which fraction to count as the CCF depends on whether it is a gain or a loss - determine this below)
          for (region in regions)
          {
            nBobs           <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean) #copyITHRDS$regionComparisonList[[region]][[segnBme]]$nB_mean
            nB1             <- nB1#floor(as.numeric(copyITHRDS$regionComparisonList[[region]][[segnBme]]$nB_mean))
            nB2             <- nB2#ceiling(as.numeric(copyITHRDS$regionComparisonList[[region]][[segnBme]]$nB_mean))
            fracsToUse      <- getCCFractionNeeded(nBobs,nB1,nB2)
            
            fracB1  <- max(fracsToUse)
            fracB2  <- min(fracsToUse)
            
            segCCFvalues[[segName]][[region]]$fracB1_mean  <- fracB1
            segCCFvalues[[segName]][[region]]$fracB2_mean  <- fracB2
            segCCFvalues[[segName]][[region]]$nBobs_mean   <- nBobs
            #segCCFvalues[[segnBme]][[region]]$CCF_mean     <- fracB1
            
            # lets also add the raw values (just in case)
            nBobs_all      <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nB<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nB) #copyITHRDS$regionComparisonList[[region]][[segnBme]]$nB_mean
            nBobs_all       <- ifelse(nBobs_all>=maxCPN,maxCPN,nBobs_all)
            nB1             <- nB1#floor(as.numeric(copyITHRDS$regionComparisonList[[region]][[segnBme]]$nB_mean))
            nB2             <- nB2#ceiling(as.numeric(copyITHRDS$regionComparisonList[[region]][[segnBme]]$nB_mean))
            fracsToUse     <- sapply(nBobs_all,getCCFractionNeeded,nB1,nB2)
            fracsToUse[fracsToUse=="Inf"] <- 1
            
            fracB1_all      <-apply(fracsToUse,2,max)
            # fracifelse(as.numeric(nB2)-as.numeric(nBobs_all)<0
            #  ,abs(as.numeric(nB2)-as.numeric(nBobs_all))+1
            #  ,as.numeric(nB2)-as.numeric(nBobs_all))
            fracB2_all       <-apply(fracsToUse,2,min)
            segCCFvalues[[segName]][[region]]$fracB1_all  <- fracB1_all
            segCCFvalues[[segName]][[region]]$fracB2_all  <- fracB2_all
            segCCFvalues[[segName]][[region]]$nBobs_all   <- nBobs_all
            #segCCFvalues[[segnBme]][[region]]$fracB2_all  <- fracB2
            
            
          }
        }
        
        if(!allClonal)
        {
          # we can determine the low and high values easily 
          # we can determine the low and high values easily 
          nB1  <- min(unique(floor(obsBvals)))
          nB2  <- max(unique(ceiling(obsBvals)))
          
          segCCFvalues[[segName]]$nB1   <- nB1
          segCCFvalues[[segName]]$nB2   <- nB2
          
          # now go region by region, and check the fraction of each (which fraction to count as the CCF depends on whether it is a gain or a loss - determine this below)
          for (region in regions)
          {
            nBobs           <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean) #copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean
            nB1             <- nB1#floor(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean))
            nB2             <- nB2#ceiling(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean))
            
            fracsToUse      <- getCCFractionNeeded(nBobs,nB1,nB2)
            fracB1  <-       fracsToUse[1]#ifelse(as.numeric(nB2)-as.numeric(nBobs)<0
            # ,abs(as.numeric(nB2)-as.numeric(nBobs))+1
            # ,as.numeric(nB2)-as.numeric(nBobs))
            fracB2  <- fracsToUse[2]#ifelse(1-fracB1<0,0,1-fracB1)
            
            segCCFvalues[[segName]][[region]]$fracB1_mean  <- fracB1
            segCCFvalues[[segName]][[region]]$fracB2_mean  <- fracB2
            segCCFvalues[[segName]][[region]]$nBobs_mean   <- nBobs
            #segCCFvalues[[segName]][[region]]$CCF_mean     <- fracB1
            
            # lets also add the raw values (just in case)
            nBobs_all      <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nB<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nB) #copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean
            nBobs_all       <- ifelse(nBobs_all>=maxCPN,maxCPN,nBobs_all)
            
            nB1             <- nB1#floor(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean))
            nB2             <- nB2#ceiling(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean))
            fracsToUse     <- sapply(nBobs_all,getCCFractionNeeded,nB1,nB2)
            fracsToUse[fracsToUse=="Inf"] <- 1
            
            
            fracB1_all      <-fracsToUse[1,] 
            # fracifelse(as.numeric(nB2)-as.numeric(nBobs_all)<0
            #  ,abs(as.numeric(nB2)-as.numeric(nBobs_all))+1
            #  ,as.numeric(nB2)-as.numeric(nBobs_all))
            fracB2_all       <-fracsToUse[2,] 
            
            #tTestWithErrorHandling(fracB1_all,mu = 1)
            
            segCCFvalues[[segName]][[region]]$fracB1_all  <- fracB1_all
            segCCFvalues[[segName]][[region]]$fracB2_all  <- fracB2_all
            segCCFvalues[[segName]][[region]]$nBobs_all   <- nBobs_all
            #segCCFvalues[[segName]][[region]]$fracB2_all  <- fracB2
            
            
          }
        }
        
        
        
        
      }
      
      if(length(unique(round(obsBvals)))!=1)
      {
        # we can determine the low and high values easily 
        nB1  <- min(unique(floor(obsBvals)))
        nB2  <- max(unique(ceiling(obsBvals)))
        
        segCCFvalues[[segName]]$nB1   <- nB1
        segCCFvalues[[segName]]$nB2   <- nB2
        
        # now go region by region, and check the fraction of each (which fraction to count as the CCF depends on whether it is a gain or a loss - determine this below)
        for (region in regions)
        {
          nBobs           <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean) #copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean
          nB1             <- nB1#floor(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean))
          nB2             <- nB2#ceiling(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean))
          
          fracsToUse      <- getCCFractionNeeded(nBobs,nB1,nB2)
          fracB1  <-       fracsToUse[1]#ifelse(as.numeric(nB2)-as.numeric(nBobs)<0
          # ,abs(as.numeric(nB2)-as.numeric(nBobs))+1
          # ,as.numeric(nB2)-as.numeric(nBobs))
          fracB2  <- fracsToUse[2]#ifelse(1-fracB1<0,0,1-fracB1)
          
          segCCFvalues[[segName]][[region]]$fracB1_mean  <- fracB1
          segCCFvalues[[segName]][[region]]$fracB2_mean  <- fracB2
          segCCFvalues[[segName]][[region]]$nBobs_mean   <- nBobs
          #segCCFvalues[[segName]][[region]]$CCF_mean     <- fracB1
          
          # lets also add the raw values (just in case)
          nBobs_all      <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nB<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nB) #copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean
          nBobs_all       <- ifelse(nBobs_all>=maxCPN,maxCPN,nBobs_all)
          
          nB1             <- nB1#floor(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean))
          nB2             <- nB2#ceiling(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean))
          fracsToUse     <- sapply(nBobs_all,getCCFractionNeeded,nB1,nB2)
          fracsToUse[fracsToUse=="Inf"] <- 1
          
          
          
          fracB1_all      <-fracsToUse[1,] 
          # fracifelse(as.numeric(nB2)-as.numeric(nBobs_all)<0
          #  ,abs(as.numeric(nB2)-as.numeric(nBobs_all))+1
          #  ,as.numeric(nB2)-as.numeric(nBobs_all))
          fracB2_all       <-fracsToUse[2,] 
          
          #tTestWithErrorHandling(fracB1_all,mu = 1)
          
          segCCFvalues[[segName]][[region]]$fracB1_all  <- fracB1_all
          segCCFvalues[[segName]][[region]]$fracB2_all  <- fracB2_all
          segCCFvalues[[segName]][[region]]$nBobs_all   <- nBobs_all
          #segCCFvalues[[segName]][[region]]$fracB2_all  <- fracB2
          
          
        }
      }
      
      #next;
    }
    
    # let's check whether the ceiling and floor lead to only 2 options
    if(length(unique(c(ceiling(obsBvals),floor(obsBvals))))>2)
    {
      # check whether only one state when rounding
      if(length(unique(round(obsBvals)))==1)
      {
        roundedIntger <- round(obsBvals)[1]
        allClonal <- TRUE
        nonClonalRegion <- c()
        for (region in regions)
        {
          Test <- tTestWithErrorHandling(copyITHRDS$regionComparisonList[[region]][[segName]]$nB,mu = roundedIntger,conf.level = confInt)
          if(Test$p.value<(1-confInt))
          {
            allClonal <- FALSE
            nonClonalRegion <- c(nonClonalRegion,region)
          }
        }
        
        
        ttestClonal <-TRUE
        nonClonalRegionBAF <- c()
        for (region in regions)
        {
          
          segSplit <- unlist(strsplit(segName,split=":"))
          test <- copyITHRDS$newBAFtable[[region]]$regionBAFtable[copyITHRDS$newBAFtable[[region]]$regionBAFtable$pos>=as.numeric(segSplit[2])&copyITHRDS$newBAFtable[[region]]$regionBAFtable$pos<=as.numeric(segSplit[3])&copyITHRDS$newBAFtable[[region]]$regionBAFtable$chrom==as.numeric(segSplit[1]),,drop=FALSE]
          nAraw   <- round(mean(copyITHRDS$regionComparisonList[[region]][[segName]]$nA))
          nBraw   <- round(mean(copyITHRDS$regionComparisonList[[region]][[segName]]$nB))
          
          tTestSpec <- tTestWithErrorHandling(c(test[test$nA==1,'BAF'],1-test[test$nB==1,'BAF']),mu = funCalcExpectedBAF(expBvals[region],expAvals[region],ploidyEstimates[region],purityEstimates[region]))
          #tTest <- c(tTest,tTestSpec$p.value)
          if(tTestSpec$p.value<(1-confInt))
          {
            ttestClonal <- FALSE
            nonClonalRegionBAF <- c(nonClonalRegionBAF,region)
          }
          
          
        }
        
        
        # finally, check what the CCF is of each region if you assume there is subclonality
        maxCCF <- 0
        nB1  <- min(unique(floor(obsBvals)))
        nB2  <- max(unique(ceiling(obsBvals)))
        
        # now go region by region, and check the fraction of each (which fraction to count as the CCF depends on whether it is a gain or a loss - determine this below)
        for (region in regions)
        {
          nBobs           <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean) #copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean
          nB1             <- nB1#floor(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean))
          nB2             <- nB2#ceiling(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean))
          
          fracsToUse      <- getCCFractionNeeded(nBobs,nB1,nB2)
          minCCF <- min(fracsToUse)
          if(minCCF>maxCCF)
          {
            maxCCF <- minCCF
          }
        }
        
        # next, determine whether to explore subclonality
        if(allClonal==TRUE|ttestClonal==TRUE|maxCCF<minCCFtoConsider)
        {
          allClonal <- TRUE
        }
        
        
        if(allClonal)
        {
          nB1  <- round(obsBvals)[1]
          nB2  <- round(obsBvals)[1]
          
          segCCFvalues[[segName]]$nB1   <- nB1
          segCCFvalues[[segName]]$nB2   <- nB2
          
          # now go region by region, and check the fraction of each (which fraction to count as the CCF depends on whether it is a gain or a loss - determine this below)
          for (region in regions)
          {
            nBobs           <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean) #copyITHRDS$regionComparisonList[[region]][[segnBme]]$nB_mean
            nB1             <- nB1#floor(as.numeric(copyITHRDS$regionComparisonList[[region]][[segnBme]]$nB_mean))
            nB2             <- nB2#ceiling(as.numeric(copyITHRDS$regionComparisonList[[region]][[segnBme]]$nB_mean))
            fracsToUse      <- getCCFractionNeeded(nBobs,nB1,nB2)
            
            fracB1  <- max(fracsToUse)
            # ifelse(as.numeric(nB2)-as.numeric(nBobs)<0
            #   ,abs(as.numeric(nB2)-as.numeric(nBobs))+1
            #   ,as.numeric(nB2)-as.numeric(nBobs))
            fracB2  <- min(fracsToUse)
            
            segCCFvalues[[segName]][[region]]$fracB1_mean  <- fracB1
            segCCFvalues[[segName]][[region]]$fracB2_mean  <- fracB2
            segCCFvalues[[segName]][[region]]$nBobs_mean   <- nBobs
            #segCCFvalues[[segnBme]][[region]]$CCF_mean     <- fracB1
            
            # lets also add the raw values (just in case)
            nBobs_all      <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nB<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nB) #copyITHRDS$regionComparisonList[[region]][[segnBme]]$nB_mean
            nBobs_all       <- ifelse(nBobs_all>=maxCPN,maxCPN,nBobs_all)
            nB1             <- nB1#floor(as.numeric(copyITHRDS$regionComparisonList[[region]][[segnBme]]$nB_mean))
            nB2             <- nB2#ceiling(as.numeric(copyITHRDS$regionComparisonList[[region]][[segnBme]]$nB_mean))
            fracsToUse     <- sapply(nBobs_all,getCCFractionNeeded,nB1,nB2)
            fracsToUse[fracsToUse=="Inf"] <- 1
            
            
            fracB1_all      <-apply(fracsToUse,2,max)
            # fracifelse(as.numeric(nB2)-as.numeric(nBobs_all)<0
            #  ,abs(as.numeric(nB2)-as.numeric(nBobs_all))+1
            #  ,as.numeric(nB2)-as.numeric(nBobs_all))
            fracB2_all       <-apply(fracsToUse,2,min)
            segCCFvalues[[segName]][[region]]$fracB1_all  <- fracB1_all
            segCCFvalues[[segName]][[region]]$fracB2_all  <- fracB2_all
            segCCFvalues[[segName]][[region]]$nBobs_all   <- nBobs_all
            #segCCFvalues[[segnBme]][[region]]$fracB2_all  <- fracB2
            
            
          }
        }
        
        if(!allClonal)
        {
          # do we twist up or down
          if(abs(min(max(ceiling(obsBvals))-obsBvals))<abs(max(min(floor(obsBvals))-obsBvals)))
          {
            # take the ceiling 
            # we can determine the low and high values easily 
            nB1  <- round(obsBvals)[1]
            nB2  <- max(unique(ceiling(obsBvals)))
            
            segCCFvalues[[segName]]$nB1   <- nB1
            segCCFvalues[[segName]]$nB2   <- nB2
            
            # now go region by region, and check the fraction of each (which fraction to count as the CCF depends on whether it is a gain or a loss - determine this below)
            for (region in regions)
            {
              nBobs           <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean) #copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean
              nB1             <- nB1#floor(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean))
              nB2             <- nB2#ceiling(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean))
              fracsToUse      <- getCCFractionNeeded(nBobs,nB1,nB2)
              
              fracB1  <- fracsToUse[1]
              # ifelse(as.numeric(nB2)-as.numeric(nBobs)<0
              #   ,abs(as.numeric(nB2)-as.numeric(nBobs))+1
              #   ,as.numeric(nB2)-as.numeric(nBobs))
              fracB2  <- fracsToUse[2]
              
              segCCFvalues[[segName]][[region]]$fracB1_mean  <- fracB1
              segCCFvalues[[segName]][[region]]$fracB2_mean  <- fracB2
              segCCFvalues[[segName]][[region]]$nBobs_mean   <- nBobs
              #segCCFvalues[[segName]][[region]]$CCF_mean     <- fracB1
              
              # lets also add the raw values (just in case)
              nBobs_all      <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nB<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nB) #copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean
              nBobs_all       <- ifelse(nBobs_all>=maxCPN,maxCPN,nBobs_all)
              nB1             <- nB1#floor(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean))
              nB2             <- nB2#ceiling(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean))
              fracsToUse     <- sapply(nBobs_all,getCCFractionNeeded,nB1,nB2)
              fracsToUse[fracsToUse=="Inf"] <- 1
              
              
              fracB1_all      <-fracsToUse[1,] 
              # fracifelse(as.numeric(nB2)-as.numeric(nBobs_all)<0
              #  ,abs(as.numeric(nB2)-as.numeric(nBobs_all))+1
              #  ,as.numeric(nB2)-as.numeric(nBobs_all))
              fracB2_all       <-fracsToUse[2,] 
              segCCFvalues[[segName]][[region]]$fracB1_all  <- fracB1_all
              segCCFvalues[[segName]][[region]]$fracB2_all  <- fracB2_all
              segCCFvalues[[segName]][[region]]$nBobs_all   <- nBobs_all
              #segCCFvalues[[segName]][[region]]$fracB2_all  <- fracB2
              
              
            }
            
            #next;
            
          }
          if(abs(min(max(ceiling(obsBvals))-obsBvals))>=abs(max(min(floor(obsBvals))-obsBvals)))
          {
            # take the floor 
            # we can determine the low and high values easily 
            nB1  <- min(unique(floor(obsBvals)))
            nB2  <- max(unique(round(obsBvals)))
            
            segCCFvalues[[segName]]$nB1   <- nB1
            segCCFvalues[[segName]]$nB2   <- nB2
            
            # now go region by region, and check the fraction of each (which fraction to count as the CCF depends on whether it is a gain or a loss - determine this below)
            for (region in regions)
            {
              nBobs           <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean) #copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean
              nB1             <- nB1#floor(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean))
              nB2             <- nB2#ceiling(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean))
              fracsToUse      <- getCCFractionNeeded(nBobs,nB1,nB2)
              
              fracB1  <- fracsToUse[1]
              # ifelse(as.numeric(nB2)-as.numeric(nBobs)<0
              #   ,abs(as.numeric(nB2)-as.numeric(nBobs))+1
              #   ,as.numeric(nB2)-as.numeric(nBobs))
              fracB2  <- fracsToUse[2]
              #ifelse(1-fracB1<0,0,1-fracB1)
              
              segCCFvalues[[segName]][[region]]$fracB1_mean  <- fracB1
              segCCFvalues[[segName]][[region]]$fracB2_mean  <- fracB2
              segCCFvalues[[segName]][[region]]$nBobs_mean   <- nBobs
              #segCCFvalues[[segName]][[region]]$CCF_mean     <- fracB1
              
              # lets also add the raw values (just in case)
              nBobs_all      <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nB<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nB) #copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean
              nBobs_all       <- ifelse(nBobs_all>=maxCPN,maxCPN,nBobs_all)
              
              nB1             <- nB1#floor(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean))
              nB2             <- nB2#ceiling(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean))
              fracsToUse     <- sapply(nBobs_all,getCCFractionNeeded,nB1,nB2)
              fracsToUse[fracsToUse=="Inf"] <- 1
              
              
              fracB1_all      <-fracsToUse[1,] 
              # fracifelse(as.numeric(nB2)-as.numeric(nBobs_all)<0
              #  ,abs(as.numeric(nB2)-as.numeric(nBobs_all))+1
              #  ,as.numeric(nB2)-as.numeric(nBobs_all))
              fracB2_all       <-fracsToUse[2,] 
              
              segCCFvalues[[segName]][[region]]$fracB1_all  <- fracB1_all
              segCCFvalues[[segName]][[region]]$fracB2_all  <- fracB2_all
              segCCFvalues[[segName]][[region]]$nBobs_all   <- nBobs_all
              #segCCFvalues[[segName]][[region]]$fracB2_all  <- fracB2
              
              
            }
            
            # next;
            
          }
        }
        
      }
      
      #check whether there are only two states when rounding
      if(length(unique(round(obsBvals)))==2)
      {
        
        # we can determine the low and high values easily
        nB1  <- min(round(obsBvals))
        nB2  <- max(round(obsBvals))
        
        segCCFvalues[[segName]]$nB1   <- nB1
        segCCFvalues[[segName]]$nB2   <- nB2
        
        # now go region by region, and check the fraction of each (which fraction to count as the CCF depends on whether it is a gain or a loss - determine this below)
        for (region in regions)
        {
          nBobs           <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean) #copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean
          nB1             <- nB1#floor(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean))
          nB2             <- nB2#ceiling(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean))
          fracsToUse      <- getCCFractionNeeded(nBobs,nB1,nB2)
          fracB1  <- fracsToUse[1] 
          # ifelse(as.numeric(nB2)-as.numeric(nBobs)<0
          # ,abs(as.numeric(nB2)-as.numeric(nBobs))+1
          # ,as.numeric(nB2)-as.numeric(nBobs))
          fracB2  <- fracsToUse[2]
          # ifelse(1-fracB1<0,0,1-fracB1)
          
          segCCFvalues[[segName]][[region]]$fracB1_mean  <- fracB1
          segCCFvalues[[segName]][[region]]$fracB2_mean  <- fracB2
          segCCFvalues[[segName]][[region]]$nBobs_mean   <- nBobs
          #segCCFvalues[[segName]][[region]]$CCF_mean     <- fracB1
          
          # lets also add the raw values (just in case)
          nBobs_all      <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nB<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nB) #copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean
          nBobs_all       <- ifelse(nBobs_all>=maxCPN,maxCPN,nBobs_all)
          
          nB1             <- nB1#floor(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean))
          nB2             <- nB2#ceiling(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean))
          fracsToUse     <- sapply(nBobs_all,getCCFractionNeeded,nB1,nB2)
          fracsToUse[fracsToUse=="Inf"] <- 1
          
          
          fracB1_all      <-fracsToUse[1,] 
          # fracifelse(as.numeric(nB2)-as.numeric(nBobs_all)<0
          #  ,abs(as.numeric(nB2)-as.numeric(nBobs_all))+1
          #  ,as.numeric(nB2)-as.numeric(nBobs_all))
          fracB2_all       <-fracsToUse[2,] 
          # ifelse(1-fracB1_all<0,0,1-fracB1_all)
          
          segCCFvalues[[segName]][[region]]$fracB1_all  <- fracB1_all
          segCCFvalues[[segName]][[region]]$fracB2_all  <- fracB2_all
          segCCFvalues[[segName]][[region]]$nBobs_all   <- nBobs_all
          #segCCFvalues[[segName]][[region]]$fracB2_all  <- fracB2
          
          
        }
        
        #next;
      }
      
      
      
      
      # check whether an amplification greater than 8
      
      #if(min(expAvals))
      
      # posValues  <- unique(sort(expAvals))
      # lowValue   <- posValues[1]
      # highValue  <- posValues[2]
      # 
      # # choose the max and min that maximise the CCF across the tumour regions (in other words, minimize heterogeneity where possible)
      # posFracA1 <- c()
      # posFracA2 <- c()
      # valuesChosen <- c()
      # for (x in 1:(length(posValues)-1))
      # {
      #   for (y in (x+1):length(posValues))
      #   {
      #     
      #     posLowValue  <- posValues[x]
      #     posHighValue <- posValues[y]
      #     
      #     CCFvalues <- c()
      #     for (i in 1:length(obsBvals))
      #     {
      #       CCFvalues <-c(CCFvalues,(getCCFractionNeeded(obsBvals[i],posLowValue,posHighValue))[2])
      #     }
      #     valuesChosen <- rbind(valuesChosen,c(posLowValue,posHighValue))
      #     posFracA1 <- rbind(posFracA1,CCFvalues)
      #   }
      # }
      
      if(length(unique(round(obsBvals)))>2)
      {
        # we can determine the low and high values easily 
        nB1  <- min(unique(floor(obsBvals+minCCFtoConsider)))
        nB2  <- max(unique(ceiling(obsBvals-minCCFtoConsider)))
        
        segCCFvalues[[segName]]$nB1   <- nB1
        segCCFvalues[[segName]]$nB2   <- nB2
        
        # now go region by region, and check the fraction of each (which fraction to count as the CCF depends on whether it is a gain or a loss - determine this below)
        for (region in regions)
        {
          nBobs           <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean) #copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean
          nB1             <- nB1#floor(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean))
          nB2             <- nB2#ceiling(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean))
          
          fracsToUse      <- getCCFractionNeeded(nBobs,nB1,nB2)
          fracB1  <- fracsToUse[1] 
          # ifelse(as.numeric(nB2)-as.numeric(nBobs)<0
          # ,abs(as.numeric(nB2)-as.numeric(nBobs))+1
          # ,as.numeric(nB2)-as.numeric(nBobs))
          fracB2  <- fracsToUse[2]
          # ifelse(1-fracB1<0,0,1-fracB1)
          
          segCCFvalues[[segName]][[region]]$fracB1_mean  <- fracB1
          segCCFvalues[[segName]][[region]]$fracB2_mean  <- fracB2
          segCCFvalues[[segName]][[region]]$nBobs_mean   <- nBobs
          #segCCFvalues[[segName]][[region]]$CCF_mean     <- fracB1
          
          # lets also add the raw values (just in case)
          nBobs_all      <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nB<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nB) #copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean
          nBobs_all       <- ifelse(nBobs_all>=maxCPN,maxCPN,nBobs_all)
          
          nB1             <- nB1#floor(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean))
          nB2             <- nB2#ceiling(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean))
          
          fracsToUse     <- sapply(nBobs_all,getCCFractionNeeded,nB1,nB2)
          fracsToUse[fracsToUse=="Inf"] <- 1
          
          
          fracB1_all      <-fracsToUse[1,] 
          # fracifelse(as.numeric(nB2)-as.numeric(nBobs_all)<0
          #  ,abs(as.numeric(nB2)-as.numeric(nBobs_all))+1
          #  ,as.numeric(nB2)-as.numeric(nBobs_all))
          fracB2_all       <-fracsToUse[2,] 
          # ifelse(1-fracB1_all<0,0,1-fracB1_all)
          
          
          # fracB1_all  <- ifelse(as.numeric(nB2)-as.numeric(nBobs_all)<0
          #                       ,abs(as.numeric(nB2)-as.numeric(nBobs_all))+1
          #                       ,as.numeric(nB2)-as.numeric(nBobs_all))
          # fracB2_all  <- ifelse(1-fracB1_all<0,0,1-fracB1_all)
          
          segCCFvalues[[segName]][[region]]$fracB1_all  <- fracB1_all
          segCCFvalues[[segName]][[region]]$fracB2_all  <- fracB2_all
          segCCFvalues[[segName]][[region]]$nBobs_all   <- nBobs_all
          #segCCFvalues[[segName]][[region]]$fracB2_all  <- fracB2
          
          
        }
      }
      
      
      
      #next;
    }
    
    
    return(segCCFvalues)
    
  }
  
  # also create a table with the output
  simpleOutPut <- c()
  for (segName in segNames)
  {
    testOut <- getIndividualCCFestimates(segName,regions,copyITHRDS,minCCFtoConsider = minCCFtoConsider,confInt = confInt)
    segCCFvalues[[segName]] <- testOut
    
    rowToAdd <- cbind(segName,as.numeric(testOut[[segName]]$nA1),as.numeric(testOut[[segName]]$nA2))
    
    for (region in regions)
    {
      rowToAdd <- cbind(rowToAdd
                        ,as.numeric(testOut[[segName]][[region]]$fracA1_mean)
                        ,as.numeric(testOut[[segName]][[region]]$fracA2_mean))
    }
    simpleOutPut <- rbind(simpleOutPut,rowToAdd)
    
  }
  #print('got this far')
  #print(segNames)
  
  # first plot the segments (with same integer across regions) and those with different across regions
  CCFtables            <- list()
  CCFtables$fracA1     <- c()
  CCFtables$fracA2     <- c()
  CCFtables$fracB1     <- c()
  CCFtables$fracB2     <- c()
  CCFtables$A_CCF      <- c()
  CCFtables$B_CCF      <- c()
  
  for (segName in segNames)
  {
    
    regionFracA1 <- c()
    regionFracA2 <- c()
    regionFracB1 <- c()
    regionFracB2 <- c()
    
    regionA_CCF <- c()
    regionB_CCF <- c()
    for (region in regions)
    {
      regionFracA1 <- c(regionFracA1,segCCFvalues[[segName]][[segName]][[region]]$fracA1_mean)
      regionFracA2 <- c(regionFracA2,segCCFvalues[[segName]][[segName]][[region]]$fracA2_mean)
      regionFracB1 <- c(regionFracB1,segCCFvalues[[segName]][[segName]][[region]]$fracB1_mean)
      regionFracB2 <- c(regionFracB2,segCCFvalues[[segName]][[segName]][[region]]$fracB2_mean)
      
      regionA_CCF  <- c(regionA_CCF,max(segCCFvalues[[segName]][[segName]][[region]]$fracA2_mean,segCCFvalues[[segName]][[segName]][[region]]$fracA1_mean))
      regionB_CCF  <- c(regionB_CCF,max(segCCFvalues[[segName]][[segName]][[region]]$fracB2_mean,segCCFvalues[[segName]][[segName]][[region]]$fracB1_mean))
    }
    
    CCFtables$fracA1 <- rbind(CCFtables$fracA1,regionFracA1)
    CCFtables$fracA2 <- rbind(CCFtables$fracA2,regionFracA2)
    
    CCFtables$fracB1 <- rbind(CCFtables$fracB1,regionFracB1)
    CCFtables$fracB2 <- rbind(CCFtables$fracB2,regionFracB2)
    
    CCFtables$A_CCF <- rbind(CCFtables$A_CCF,regionA_CCF)
    CCFtables$B_CCF <- rbind(CCFtables$B_CCF,regionB_CCF)
    
    
  }
  
  rownames(CCFtables$fracA1) <- segNames
  rownames(CCFtables$fracA2) <- segNames
  rownames(CCFtables$fracB1) <- segNames
  rownames(CCFtables$fracB2) <- segNames
  rownames(CCFtables$A_CCF)  <- segNames
  rownames(CCFtables$B_CCF)  <- segNames
  
  colnames(CCFtables$fracA1) <- regions
  colnames(CCFtables$fracA2) <- regions
  colnames(CCFtables$fracB1) <- regions
  colnames(CCFtables$fracB2) <- regions
  colnames(CCFtables$A_CCF)  <- regions
  colnames(CCFtables$B_CCF)  <- regions
  
  if(plotting==TRUE)
  {
    #first, determine what is heterogeneous and what is classified as homogeneous.
    # next, simply consider all segNames and any where there is a difference in heterogeneity (according to new scores)
    regionsToConsider <- grep(":",names(copyITHRDS$regionComparisonList),invert=TRUE,value=TRUE)
    segsToConsider    <- names(copyITHRDS$regionComparisonList[[1]])
    
    min.gain.trunk           <- data.frame(matrix(NA,ncol=length(regionsToConsider),nrow=length(segNames)))
    rownames(min.gain.trunk) <- segNames
    colnames(min.gain.trunk) <- regionsToConsider
    min.gain <- min.loss <- min.loss.trunk  <- min.loss.branch <- min.gain.branch <- min.gain.trunk
    
    # go region by region
    for (region in regionsToConsider)
    {
      regionOut <- copyITHRDS$regionComparisonList[[region]]
      for (seg in segNames)
      {
        regionSegStatus <- regionOut[[seg]]$segStatus
        if(length(regionSegStatus)==0)
        {
          next;
        }
        min.gain[seg,region] <- ifelse(as.numeric(regionSegStatus[gainPValToUse])<pValThresh,1,0)
        min.loss[seg,region] <- ifelse(as.numeric(regionSegStatus[lossPValToUse])<pValThresh,1,0)
        
      }
      
    }
    
    
    # remove the NAs
    min.gain <- min.gain[!is.na(rowSums(min.gain)),,drop=FALSE]
    min.loss <- min.loss[!is.na(rowSums(min.loss)),,drop=FALSE]
    
    trunkGain  <- names(which(rowSums(min.gain)==length(regionsToConsider)))
    trunkLoss  <- names(which(rowSums(min.loss)==length(regionsToConsider)))
    branchGain <- names(which(rowSums(min.gain)<length(regionsToConsider)&rowSums(min.gain)>0))
    branchLoss <- names(which(rowSums(min.loss)<length(regionsToConsider)&rowSums(min.loss)>0))
    
    #perform a test to see whether there is a significant difference
    #between the genome of the different regions
    revertToTrunkGain <- c()
    for (seg in branchGain)
    {
      
      GainRegions   <-     names(which(unlist(c(min.gain[seg,])>0)))
      NoGainRegions <-     names(which(unlist(c(min.gain[seg,])==0)))
      
      # make the assumption that the sample with the biggest p-value will be more likely to be similar to the others
      referenceGain <- GainRegions[1]
      refRegionOut <- copyITHRDS$regionComparisonList[[referenceGain]]
      refRegionPval <- as.numeric(refRegionOut[[seg]]$segStatus[gainPValToUse])
      
      
      for(gainRegion in GainRegions)
      {
        
        testRegionOut <- copyITHRDS$regionComparisonList[[gainRegion]]
        if( as.numeric(testRegionOut[[seg]]$segStatus[gainPValToUse])>refRegionPval)
        {
          referenceGain <- gainRegion
          refRegionOut <- copyITHRDS$regionComparisonList[[referenceGain]]
          refRegionPval <- as.numeric(refRegionOut[[seg]]$segStatus[gainPValToUse])
          
        }
        
      }
      
      #next, check whether each NoGainRegion is significantly different from the reference region
      
      reconsider <- c()
      for (noGainRegion in NoGainRegions)
      {
        combinedName <- c(referenceGain,noGainRegion)
        combinedName <- paste(sort(combinedName)[1],sort(combinedName)[2],sep=":")
        if(copyITHRDS$regionComparisonList[[combinedName]][[seg]]$combinedPairedPval>PvalToSeparate)
        {
          reconsider <- c(reconsider,noGainRegion)
        }
        
      }
      
      if(identical(reconsider,NoGainRegions))
      {
        revertToTrunkGain <- c(revertToTrunkGain,seg)
      }
      
      
      
      
    }
    
    revertToTrunkLoss <- c()
    for (seg in branchLoss)
    {
      
      LossRegions   <-     names(which(unlist(c(min.loss[seg,])>0)))
      NoLossRegions <-     names(which(unlist(c(min.loss[seg,])==0)))
      
      # make the assumption that the sample with the smallest p-value will be more likely to be similar to the others
      referenceLoss <- LossRegions[1]
      refRegionOut <- copyITHRDS$regionComparisonList[[referenceLoss]]
      refRegionPval <- as.numeric(refRegionOut[[seg]]$segStatus[lossPValToUse])
      
      
      for(lossRegion in LossRegions)
      {
        
        testRegionOut <- copyITHRDS$regionComparisonList[[lossRegion]]
        if( as.numeric(testRegionOut[[seg]]$segStatus[lossPValToUse])>refRegionPval)
        {
          referenceLoss <- lossRegion
          refRegionOut <- copyITHRDS$regionComparisonList[[referenceLoss]]
          refRegionPval <- as.numeric(refRegionOut[[seg]]$segStatus[lossPValToUse])
          
        }
        
      }
      
      #next, check whether each NoLossRegion is significantly different from the reference region
      
      reconsider <- c()
      for (noLossRegion in NoLossRegions)
      {
        combinedName <- c(referenceLoss,noLossRegion)
        combinedName <- paste(sort(combinedName)[1],sort(combinedName)[2],sep=":")
        if(copyITHRDS$regionComparisonList[[combinedName]][[seg]]$combinedPairedPval>PvalToSeparate)
        {
          reconsider <- c(reconsider,noLossRegion)
        }
        
      }
      
      if(identical(reconsider,NoLossRegions))
      {
        revertToTrunkLoss <- c(revertToTrunkLoss,seg)
      }
      
      
      
      
    }
    
    min.gain[revertToTrunkGain,] <- 1
    min.loss[revertToTrunkLoss,] <- 1
    
    trunkGain  <- names(which(rowSums(min.gain)==length(regionsToConsider)))
    trunkLoss  <- names(which(rowSums(min.loss)==length(regionsToConsider)))
    branchGain <- names(which(rowSums(min.gain)<length(regionsToConsider)&rowSums(min.gain)>0))
    branchLoss <- names(which(rowSums(min.loss)<length(regionsToConsider)&rowSums(min.loss)>0))
    
    
    hist(CCFtables$fracA1[sameIntegerA,],breaks=50,probability = TRUE,ylim=c(0,10),col='#3182bd99',xlab='CCF estimate of ancestral A-allele state',main='CCF ancestral A of ubuqitous and heterogeneous events',xlim=c(0,1.25))
    hist(CCFtables$fracA1[diffIntegerA,],breaks=50,probability = TRUE,ylim=c(0,10),col='#de2d2699',add=TRUE)
    
    hist(CCFtables$fracB1[sameIntegerB,],breaks=50,probability = TRUE,ylim=c(0,10),col='#3182bd99',xlab='CCF estimate of ancestral B-allele state',main='CCF ancestral B of ubuqitous and heterogeneous events',xlim=c(0,1.25))
    hist(CCFtables$fracB1[diffIntegerB,],breaks=50,probability = TRUE,ylim=c(0,10),col='#de2d2699',add=TRUE)
    
    # great a grid for boxplots for each region
    plot(0,xlim=c(0.5,length(regions)+1),ylim=c(-0.15,1.5)
         ,xaxt='n'
         ,yaxt='n'
         ,ylab="CCF"
         ,xlab=""
         ,type='n'
         ,bty='n'
         ,main='A-allele difference')
    
    for (i in 1:length(regions))
    {
      Test <- tTestWithErrorHandling(CCFtables$A_CCF[sameIntegerA,i],CCFtables$A_CCF[diffIntegerA,i])
      beeswarm(CCFtables$A_CCF[sameIntegerA,i],at=i-0.2,cex=1.5,pch=16,col='#de2d2699',method='center',corralWidth = 0.2,corral='wrap',add=TRUE)
      medianLine <- median(CCFtables$A_CCF[sameIntegerA,i])
      segments(x0 = i-0.2-0.15,y0 = medianLine,x1 = i-0.2+0.15,y1 = medianLine,lwd = 5,col = 'black')
      
      beeswarm(CCFtables$A_CCF[diffIntegerA,i],at=i+0.2,cex=1.5,pch=16,col='#3182bd99',method='center',corralWidth = 0.2,corral='wrap',add=TRUE)
      medianLine <- median(CCFtables$A_CCF[diffIntegerA,i])
      segments(x0 = i+0.2-0.15,y0 = medianLine,x1 = i+0.2+0.15,y1 = medianLine,lwd = 5,col = 'black')
      
      symbolToPlot <- ifelse(as.numeric(Test$p.value)<0.001,"***", ifelse(as.numeric(Test$p.value)<0.01,"**",ifelse(as.numeric(Test$p.value)<0.05,"*","ns")))
      text(i,1.25,labels = symbolToPlot,cex=1.5)
      
    }
    
    axis(side = 2,at = seq(0,1.5,by=0.2),las=1)
    axis(side = 1,at = seq(1,length(regions)),labels = regions,las=2)
    
    # great a grid for boxplots for each region
    plot(0,xlim=c(0.5,length(regions)+1),ylim=c(-0.15,1.5)
         ,xaxt='n'
         ,yaxt='n'
         ,ylab="frac A1"
         ,xlab=""
         ,type='n'
         ,bty='n'
         ,main='A-allele difference')
    
    for (i in 1:length(regions))
    {
      Test <- tTestWithErrorHandling(CCFtables$fracA1[sameIntegerA,i],CCFtables$fracA1[diffIntegerA,i])
      beeswarm(CCFtables$fracA1[sameIntegerA,i],at=i-0.2,cex=1.5,pch=16,col='#de2d2699',method='center',corralWidth = 0.2,corral='wrap',add=TRUE)
      medianLine <- median(CCFtables$fracA1[sameIntegerA,i])
      segments(x0 = i-0.2-0.15,y0 = medianLine,x1 = i-0.2+0.15,y1 = medianLine,lwd = 5,col = 'black')
      
      beeswarm(CCFtables$fracA1[diffIntegerA,i],at=i+0.2,cex=1.5,pch=16,col='#3182bd99',method='center',corralWidth = 0.2,corral='wrap',add=TRUE)
      medianLine <- median(CCFtables$fracA1[diffIntegerA,i])
      segments(x0 = i+0.2-0.15,y0 = medianLine,x1 = i+0.2+0.15,y1 = medianLine,lwd = 5,col = 'black')
      
      symbolToPlot <- ifelse(as.numeric(Test$p.value)<0.001,"***", ifelse(as.numeric(Test$p.value)<0.01,"**",ifelse(as.numeric(Test$p.value)<0.05,"*","ns")))
      text(i,1.25,labels = symbolToPlot,cex=1.5)
      
    }
    
    axis(side = 2,at = seq(0,1.5,by=0.2),las=1)
    axis(side = 1,at = seq(1,length(regions)),labels = regions,las=2)
    
    
    # great a grid for boxplots for each region
    plot(0,xlim=c(0.5,length(regions)+1),ylim=c(-0.15,1.5)
         ,xaxt='n'
         ,yaxt='n'
         ,ylab="frac A2"
         ,xlab=""
         ,type='n'
         ,bty='n'
         ,main='A-allele difference')
    
    for (i in 1:length(regions))
    {
      Test <- tTestWithErrorHandling(CCFtables$fracA2[sameIntegerA,i],CCFtables$fracA2[diffIntegerA,i])
      beeswarm(CCFtables$fracA2[sameIntegerA,i],at=i-0.2,cex=1.5,pch=16,col='#de2d2699',method='center',corralWidth = 0.2,corral='wrap',add=TRUE)
      medianLine <- median(CCFtables$fracA2[sameIntegerA,i])
      segments(x0 = i-0.2-0.15,y0 = medianLine,x1 = i-0.2+0.15,y1 = medianLine,lwd = 5,col = 'black')
      
      beeswarm(CCFtables$fracA2[diffIntegerA,i],at=i+0.2,cex=1.5,pch=16,col='#3182bd99',method='center',corralWidth = 0.2,corral='wrap',add=TRUE)
      medianLine <- median(CCFtables$fracA2[diffIntegerA,i])
      segments(x0 = i+0.2-0.15,y0 = medianLine,x1 = i+0.2+0.15,y1 = medianLine,lwd = 5,col = 'black')
      
      symbolToPlot <- ifelse(as.numeric(Test$p.value)<0.001,"***", ifelse(as.numeric(Test$p.value)<0.01,"**",ifelse(as.numeric(Test$p.value)<0.05,"*","ns")))
      text(i,1.25,labels = symbolToPlot,cex=1.5)
      
    }
    
    axis(side = 2,at = seq(0,1.5,by=0.2),las=1)
    axis(side = 1,at = seq(1,length(regions)),labels = regions,las=2)
    
    # great a grid for boxplots for each region
    plot(0,xlim=c(0.5,length(regions)+1),ylim=c(-0.15,1.5)
         ,xaxt='n'
         ,yaxt='n'
         ,ylab="CCF"
         ,xlab=""
         ,type='n'
         ,bty='n'
         ,main='B-allele difference')
    
    for (i in 1:length(regions))
    {
      Test <- tTestWithErrorHandling(CCFtables$B_CCF[sameIntegerB,i],CCFtables$B_CCF[diffIntegerB,i])
      beeswarm(CCFtables$B_CCF[sameIntegerB,i],at=i-0.2,cex=1.5,pch=16,col='#de2d2699',method='center',corralWidth = 0.2,corral='wrap',add=TRUE)
      medianLine <- median(CCFtables$B_CCF[sameIntegerB,i])
      segments(x0 = i-0.2-0.15,y0 = medianLine,x1 = i-0.2+0.15,y1 = medianLine,lwd = 5,col = 'black')
      
      beeswarm(CCFtables$B_CCF[diffIntegerB,i],at=i+0.2,cex=1.5,pch=16,col='#3182bd99',method='center',corralWidth = 0.2,corral='wrap',add=TRUE)
      medianLine <- median(CCFtables$B_CCF[diffIntegerB,i])
      segments(x0 = i+0.2-0.15,y0 = medianLine,x1 = i+0.2+0.15,y1 = medianLine,lwd = 5,col = 'black')
      
      symbolToPlot <- ifelse(as.numeric(Test$p.value)<0.001,"***", ifelse(as.numeric(Test$p.value)<0.01,"**",ifelse(as.numeric(Test$p.value)<0.05,"*","ns")))
      text(i,1.25,labels = symbolToPlot,cex=1.5)
      
    }
    
    axis(side = 2,at = seq(0,1.5,by=0.2),las=1)
    axis(side = 1,at = seq(1,length(regions)),labels = regions,las=2)
    
    
    # great a grid for boxplots for each region
    plot(0,xlim=c(0.5,length(regions)+1),ylim=c(-0.15,1.5)
         ,xaxt='n'
         ,yaxt='n'
         ,ylab="frac B1"
         ,xlab=""
         ,type='n'
         ,bty='n'
         ,main='B-allele difference')
    
    for (i in 1:length(regions))
    {
      Test <- tTestWithErrorHandling(CCFtables$fracB1[sameIntegerB,i],CCFtables$fracB1[diffIntegerB,i])
      beeswarm(CCFtables$fracB1[sameIntegerB,i],at=i-0.2,cex=1.5,pch=16,col='#de2d2699',method='center',corralWidth = 0.2,corral='wrap',add=TRUE)
      medianLine <- median(CCFtables$fracB1[sameIntegerB,i])
      segments(x0 = i-0.2-0.15,y0 = medianLine,x1 = i-0.2+0.15,y1 = medianLine,lwd = 5,col = 'black')
      
      beeswarm(CCFtables$fracB1[diffIntegerB,i],at=i+0.2,cex=1.5,pch=16,col='#3182bd99',method='center',corralWidth = 0.2,corral='wrap',add=TRUE)
      medianLine <- median(CCFtables$fracB1[diffIntegerB,i])
      segments(x0 = i+0.2-0.15,y0 = medianLine,x1 = i+0.2+0.15,y1 = medianLine,lwd = 5,col = 'black')
      
      symbolToPlot <- ifelse(as.numeric(Test$p.value)<0.001,"***", ifelse(as.numeric(Test$p.value)<0.01,"**",ifelse(as.numeric(Test$p.value)<0.05,"*","ns")))
      text(i,1.25,labels = symbolToPlot,cex=1.5)
      
    }
    
    axis(side = 2,at = seq(0,1.5,by=0.2),las=1)
    axis(side = 1,at = seq(1,length(regions)),labels = regions,las=2)
    
    
    plot(0,xlim=c(0.5,length(regions)+1),ylim=c(-0.15,1.5)
         ,xaxt='n'
         ,yaxt='n'
         ,ylab="frac B2"
         ,xlab=""
         ,type='n'
         ,bty='n'
         ,main='B-allele difference')
    
    for (i in 1:length(regions))
    {
      Test <- tTestWithErrorHandling(CCFtables$fracB2[sameIntegerB,i],CCFtables$fracB2[diffIntegerB,i])
      beeswarm(CCFtables$fracB2[sameIntegerB,i],at=i-0.2,cex=1.5,pch=16,col='#de2d2699',method='center',corralWidth = 0.2,corral='wrap',add=TRUE)
      medianLine <- median(CCFtables$fracB1[sameIntegerB,i])
      segments(x0 = i-0.2-0.15,y0 = medianLine,x1 = i-0.2+0.15,y1 = medianLine,lwd = 5,col = 'black')
      
      beeswarm(CCFtables$fracB2[diffIntegerB,i],at=i+0.2,cex=1.5,pch=16,col='#3182bd99',method='center',corralWidth = 0.2,corral='wrap',add=TRUE)
      medianLine <- median(CCFtables$fracB2[diffIntegerB,i])
      segments(x0 = i+0.2-0.15,y0 = medianLine,x1 = i+0.2+0.15,y1 = medianLine,lwd = 5,col = 'black')
      
      symbolToPlot <- ifelse(as.numeric(Test$p.value)<0.001,"***", ifelse(as.numeric(Test$p.value)<0.01,"**",ifelse(as.numeric(Test$p.value)<0.05,"*","ns")))
      text(i,1.25,labels = symbolToPlot,cex=1.5)
      
    }
    
    axis(side = 2,at = seq(0,1.5,by=0.2),las=1)
    axis(side = 1,at = seq(1,length(regions)),labels = regions,las=2)
    
    
    
    # next, consider the same, but with heterogeneous and homogeneous gains and losses
    
    # great a grid for boxplots for each region
    plot(0,xlim=c(0.5,length(regions)+1),ylim=c(-0.15,1.5)
         ,xaxt='n'
         ,yaxt='n'
         ,ylab="CCF"
         ,xlab=""
         ,type='n'
         ,bty='n'
         ,main='A-allele difference trunk vs. branch gain')
    
    for (i in 1:length(regions))
    {
      Test <- tTestWithErrorHandling(CCFtables$A_CCF[trunkGain,i],CCFtables$A_CCF[branchGain,i])
      beeswarm(CCFtables$A_CCF[trunkGain,i],at=i-0.2,cex=1.5,pch=16,col='#de2d2699',method='center',corralWidth = 0.2,corral='wrap',add=TRUE)
      medianLine <- median(CCFtables$A_CCF[trunkGain,i])
      segments(x0 = i-0.2-0.15,y0 = medianLine,x1 = i-0.2+0.15,y1 = medianLine,lwd = 5,col = 'black')
      
      beeswarm(CCFtables$A_CCF[branchGain,i],at=i+0.2,cex=1.5,pch=16,col='#3182bd99',method='center',corralWidth = 0.2,corral='wrap',add=TRUE)
      medianLine <- median(CCFtables$A_CCF[branchGain,i])
      segments(x0 = i+0.2-0.15,y0 = medianLine,x1 = i+0.2+0.15,y1 = medianLine,lwd = 5,col = 'black')
      
      symbolToPlot <- ifelse(as.numeric(Test$p.value)<0.001,"***", ifelse(as.numeric(Test$p.value)<0.01,"**",ifelse(as.numeric(Test$p.value)<0.05,"*","ns")))
      text(i,1.25,labels = symbolToPlot,cex=1.5)
      
    }
    
    axis(side = 2,at = seq(0,1.5,by=0.2),las=1)
    axis(side = 1,at = seq(1,length(regions)),labels = regions,las=2)
    
    
    # great a grid for boxplots for each region
    plot(0,xlim=c(0.5,length(regions)+1),ylim=c(-0.15,1.5)
         ,xaxt='n'
         ,yaxt='n'
         ,ylab="CCF"
         ,xlab=""
         ,type='n'
         ,bty='n'
         ,main='B-allele difference trunk vs. branch gain')
    CCFtables$B_CCF[CCFtables$B_CCF=='Inf'] <- 1
    
    for (i in 1:length(regions))
    {
      Test <- tTestWithErrorHandling(CCFtables$B_CCF[trunkGain,i],CCFtables$B_CCF[branchGain,i])
      beeswarm(CCFtables$B_CCF[trunkGain,i],at=i-0.2,cex=1.5,pch=16,col='#de2d2699',method='center',corralWidth = 0.2,corral='wrap',add=TRUE)
      medianLine <- median(CCFtables$B_CCF[trunkGain,i])
      segments(x0 = i-0.2-0.15,y0 = medianLine,x1 = i-0.2+0.15,y1 = medianLine,lwd = 5,col = 'black')
      
      beeswarm(CCFtables$B_CCF[branchGain,i],at=i+0.2,cex=1.5,pch=16,col='#3182bd99',method='center',corralWidth = 0.2,corral='wrap',add=TRUE)
      medianLine <- median(CCFtables$B_CCF[branchGain,i])
      segments(x0 = i+0.2-0.15,y0 = medianLine,x1 = i+0.2+0.15,y1 = medianLine,lwd = 5,col = 'black')
      
      symbolToPlot <- ifelse(as.numeric(Test$p.value)<0.001,"***", ifelse(as.numeric(Test$p.value)<0.01,"**",ifelse(as.numeric(Test$p.value)<0.05,"*","ns")))
      text(i,1.25,labels = symbolToPlot,cex=1.5)
      
    }
    
    axis(side = 2,at = seq(0,1.5,by=0.2),las=1)
    axis(side = 1,at = seq(1,length(regions)),labels = regions,las=2)
    
    
    
    # great a grid for boxplots for each region
    plot(0,xlim=c(0.5,length(regions)+1),ylim=c(-0.15,1.5)
         ,xaxt='n'
         ,yaxt='n'
         ,ylab="CCF"
         ,xlab=""
         ,type='n'
         ,bty='n'
         ,main='A-allele difference trunk vs. branch loss')
    
    for (i in 1:length(regions))
    {
      Test <- tTestWithErrorHandling(CCFtables$A_CCF[trunkLoss,i],CCFtables$A_CCF[branchLoss,i])
      beeswarm(CCFtables$A_CCF[trunkLoss,i],at=i-0.2,cex=1.5,pch=16,col='#de2d2699',method='center',corralWidth = 0.2,corral='wrap',add=TRUE)
      medianLine <- median(CCFtables$A_CCF[trunkLoss,i])
      segments(x0 = i-0.2-0.15,y0 = medianLine,x1 = i-0.2+0.15,y1 = medianLine,lwd = 5,col = 'black')
      
      beeswarm(CCFtables$A_CCF[branchLoss,i],at=i+0.2,cex=1.5,pch=16,col='#3182bd99',method='center',corralWidth = 0.2,corral='wrap',add=TRUE)
      medianLine <- median(CCFtables$A_CCF[branchLoss,i])
      segments(x0 = i+0.2-0.15,y0 = medianLine,x1 = i+0.2+0.15,y1 = medianLine,lwd = 5,col = 'black')
      
      symbolToPlot <- ifelse(as.numeric(Test$p.value)<0.001,"***", ifelse(as.numeric(Test$p.value)<0.01,"**",ifelse(as.numeric(Test$p.value)<0.05,"*","ns")))
      text(i,1.25,labels = symbolToPlot,cex=1.5)
      
    }
    
    axis(side = 2,at = seq(0,1.5,by=0.2),las=1)
    axis(side = 1,at = seq(1,length(regions)),labels = regions,las=2)
    
    
    # great a grid for boxplots for each region
    plot(0,xlim=c(0.5,length(regions)+1),ylim=c(-0.15,1.5)
         ,xaxt='n'
         ,yaxt='n'
         ,ylab="CCF"
         ,xlab=""
         ,type='n'
         ,bty='n'
         ,main='B-allele difference trunk vs. branch loss')
    
    for (i in 1:length(regions))
    {
      Test <- tTestWithErrorHandling(CCFtables$B_CCF[trunkLoss,i],CCFtables$B_CCF[branchLoss,i])
      beeswarm(CCFtables$B_CCF[trunkLoss,i],at=i-0.2,cex=1.5,pch=16,col='#de2d2699',method='center',corralWidth = 0.2,corral='wrap',add=TRUE)
      medianLine <- median(CCFtables$B_CCF[trunkLoss,i])
      segments(x0 = i-0.2-0.15,y0 = medianLine,x1 = i-0.2+0.15,y1 = medianLine,lwd = 5,col = 'black')
      
      beeswarm(CCFtables$B_CCF[branchLoss,i],at=i+0.2,cex=1.5,pch=16,col='#3182bd99',method='center',corralWidth = 0.2,corral='wrap',add=TRUE)
      medianLine <- median(CCFtables$B_CCF[branchLoss,i])
      segments(x0 = i+0.2-0.15,y0 = medianLine,x1 = i+0.2+0.15,y1 = medianLine,lwd = 5,col = 'black')
      
      symbolToPlot <- ifelse(as.numeric(Test$p.value)<0.001,"***", ifelse(as.numeric(Test$p.value)<0.01,"**",ifelse(as.numeric(Test$p.value)<0.05,"*","ns")))
      text(i,1.25,labels = symbolToPlot,cex=1.5)
      
    }
    
    
    axis(side = 2,at = seq(0,1.5,by=0.2),las=1)
    axis(side = 1,at = seq(1,length(regions)),labels = regions,las=2)
    
    
    heatmap.2(CCFtables$fracA1
              ,dendrogram = 'column'
              #,Rowv = FALSE
              #,RowSideColors=cols[match(most.likely.cluster[rownames(heatmap.ccf)], c(1:no.optima))]
              ,col= colorpanel(100,'white','steelblue4')
              #,margins=c(20,10)
              ,trace='none'
              ,key=FALSE
    )
    
    
    heatmap.2(CCFtables$fracB1
              ,dendrogram = 'column'
              #,Rowv = FALSE
              #,RowSideColors=cols[match(most.likely.cluster[rownames(heatmap.ccf)], c(1:no.optima))]
              ,col= colorpanel(100,'white','steelblue4')
              #,margins=c(20,10)
              ,trace='none'
              ,key=FALSE
    )
    
    
    #finally plot across the genome, and for each segment shade according to how clonal the segment is
    funPlotCCFAcrossGenome(copyITHRDS,segCCFvalues)
    
  }
  # plot simply boxplots showing the 
  
  
  
  
  return(segCCFvalues)
}

getSegNewCCFEstimate           <- function(copyITHRDS
                                           ,purityEstimates
                                           ,ploidyEstimates
                                           ,highestTrustedAmp=10
                                           ,pValThresh=0.01
                                           ,PvalToSeparate=0.001
                                           ,gainPValToUse = 'gainPval'
                                           ,lossPValToUse = 'lossPval'
                                           ,minCCFtoConsider=0.15
                                           ,plotting=FALSE
                                           ,confInt=0.99999
                                           #,plotHeterogeneity
)
{
  #need to obtain an ancestral state for each copy number 
  # (initially, simply do this by voting (will fail for subclonal GD events)).
  
  regions  <- names(purityEstimates)
  regions  <- grep(":",names(copyITHRDS$regionComparisonList),invert=TRUE,value=TRUE)
  segNames <- names(copyITHRDS$regionComparisonList[[names(ploidyEstimates)[1]]])
  for (region in regions)
  {
    segNames <- intersect(segNames,names(copyITHRDS$regionComparisonList[[region]]))
  }
  #print(regions)
  #print(segNames)
  # first, check for each segName if the rounded numbers are the same 
  # - if they are the same, then we can easily calculate the CCF for each region
  sameInteger <- c()
  diffInteger <- c()
  sameIntegerA <- c()
  diffIntegerA <- c()
  sameIntegerB <- c()
  diffIntegerB <- c()
  for (segName in segNames)
  {
    nA <- c()
    nB <- c()
    
    for (region in regions)
    {
      nA <- c(nA,round(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean))
      nB <- c(nB,round(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean))
    }
    
    if(length(unique(nA))==1&length(unique(nB))==1)
    {
      sameInteger <- c(sameInteger,segName)
    }
    
    if(length(unique(nA))!=1|length(unique(nB))!=1)
    {
      diffInteger <- c(diffInteger,segName)
      
    }
    if(length(unique(nA))==1)
    {
      sameIntegerA <- c(sameIntegerA,segName)
    }
    if(length(unique(nA))!=1)
    {
      diffIntegerA <- c(diffIntegerA,segName)
    }
    
    if(length(unique(nB))==1)
    {
      sameIntegerB <- c(sameIntegerB,segName)
    }
    if(length(unique(nB))!=1)
    {
      diffIntegerB <- c(diffIntegerB,segName)
    }
    
    
  }
  
  regionChrList <- list()
  # next assess for each chromosome and region, what are the different raw-values
  for (region in regions)
  {
    
    #regionChrList[[region]] <- list()
    
    for (segName in segNames)
    {
      intNumber <- round(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean)
      chrNumber <- unlist(strsplit(segName,split=":"))[1]
      regionChrList[[region]][[as.character(intNumber)]][[chrNumber]]$nAvals <- c(regionChrList[[region]][[as.character(intNumber)]][[chrNumber]]$nAvals,copyITHRDS$regionComparisonList[[region]][[segName]]$nA)
    }
    
    for (segName in segNames)
    {
      intNumber <- round(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean)
      chrNumber <- unlist(strsplit(segName,split=":"))[1]
      regionChrList[[region]][[as.character(intNumber)]][[chrNumber]]$nBvals <- c(regionChrList[[region]][[as.character(intNumber)]][[chrNumber]]$nBvals,copyITHRDS$regionComparisonList[[region]][[segName]]$nB)
    }
    
  }
  
  
  #print(sameInteger)
  
  #
  segCCFvalues  <- list()
  maxCPN        <- highestTrustedAmp
  
  getIndividualCCFestimates_withPvals <- function(segName
                                                  ,regions
                                                  ,copyITHRDS
                                                  ,regionChrList
                                                  ,maxCPN=5
                                                  ,minCCFtoConsider=0.15
                                                  ,chrsToConsider =c(1:22)
                                                  ,confInt=0.9999
                                                  ,minOtherChrs = 3
                                                  
  )
  {
    
    segCCFvalues <- list()
    
    strsplitSegName <- unlist(strsplit(segName,split=":"))
    chrName       <- strsplitSegName[1]
    postStartName <- strsplitSegName[2]
    postEndName   <- strsplitSegName[3]
    
    # check if we have A and B for this segment
    if(is.na(copyITHRDS$regionComparisonList[[region]][[segName]]$nA[1]))
    {
      return(NA)
    }
    
    
    obsAvals <- c()
    obsBvals <- c()
    expAvals <- c()
    expBvals <- c()
    
    for (region in regions)
    {
      nAobs   <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean) #copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean
      nBobs   <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean)
      nAexp   <- max(0,round(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean))
      nBexp   <- max(0,round(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean))
      
      obsAvals <- c(obsAvals,nAobs)
      obsBvals <- c(obsBvals,nBobs)
      expAvals <- c(expAvals,nAexp)
      expBvals <- c(expBvals,nBexp)
    }
    
    names(obsAvals) <- names(obsBvals) <- names(expAvals) <- names(expBvals) <- regions
    
    #for each one determine whether it is significantly greater or less than average or than the integer
    
    chrsToConsider <- chrsToConsider[!chrsToConsider%in%chrName]
    regionPvalsDiffFromExp <- list()
    for (region in regions)
    {
      intNumA <- round(mean(copyITHRDS$regionComparisonList[[region]][[segName]]$nA))
      intNumB <- round(mean(copyITHRDS$regionComparisonList[[region]][[segName]]$nB))
      rawAs   <- copyITHRDS$regionComparisonList[[region]][[segName]]$nA
      rawBs   <- copyITHRDS$regionComparisonList[[region]][[segName]]$nB
      
      fullAs <- c()
      nrMatchedA_Chrs <- c()
      for (chr in chrsToConsider)
      {
        if(length(regionChrList[[region]][[as.character(intNumA)]][[as.character(chr)]]$nAvals)==0)
        {
          next;
        }
        if(length(regionChrList[[region]][[as.character(intNumA)]][[as.character(chr)]]$nAvals)!=0)
        {
          nrMatchedA_Chrs <- c(nrMatchedA_Chrs,as.character(chr))
          fullAs <- c(fullAs,regionChrList[[region]][[as.character(intNumA)]][[as.character(chr)]]$nAvals)
        }
        
      }
      TtestAs             <- tTestWithErrorHandling(rawAs,fullAs)
      
      # also check whether it differs significantly from the integer number
      TtestDiffFromIntA   <- tTestWithErrorHandling(rawAs,mu=intNumA)
      
      #boxplot(rawAs,fullAs)
      
      
      fullBs <- c()
      nrMatchedB_Chrs <- c()
      
      for (chr in chrsToConsider)
      {
        if(length(regionChrList[[region]][[as.character(intNumB)]][[as.character(chr)]]$nBvals)==0)
        {
          next;
        }
        if(length(regionChrList[[region]][[as.character(intNumB)]][[as.character(chr)]]$nBvals)!=0)
        {
          nrMatchedB_Chrs <- c(nrMatchedB_Chrs,as.character(chr))
          fullBs <- c(fullBs,regionChrList[[region]][[as.character(intNumB)]][[as.character(chr)]]$nBvals)
        }
        
        
      }
      TtestBs <- tTestWithErrorHandling(rawBs,fullBs)
      TtestBsFromInteger <- tTestWithErrorHandling(rawBs,mu=intNumB)
      
      #also just check if it is greater or less than the integer
      greaterThanOtherAs <- median(rawAs) > median(fullAs)
      greaterThanIntAs   <- median(rawAs) > intNumA
      
      greaterThanOtherBs <- median(rawBs) > median(fullBs)
      greaterThanIntBs   <- median(rawBs) > intNumB
      
      regionPvalsDiffFromExp[[region]]$nA_intger             <- intNumA
      regionPvalsDiffFromExp[[region]]$nApval_fromOtherChrs  <- TtestAs$p.value
      regionPvalsDiffFromExp[[region]]$nApval_fromInteger    <- TtestDiffFromIntA$p.value
      regionPvalsDiffFromExp[[region]]$nA_nrMatchedAChrs     <- nrMatchedA_Chrs
      regionPvalsDiffFromExp[[region]]$greaterThanMatchedAs     <- greaterThanOtherAs
      regionPvalsDiffFromExp[[region]]$greaterThanIntAs     <- greaterThanIntAs
      
      
      regionPvalsDiffFromExp[[region]]$nB_intger             <- intNumB
      regionPvalsDiffFromExp[[region]]$nBpval_fromOtherChrs  <- TtestBs$p.value
      regionPvalsDiffFromExp[[region]]$nBpval_fromIntger     <- TtestBsFromInteger$p.value
      regionPvalsDiffFromExp[[region]]$nB_nrMatchedAChrs     <- nrMatchedB_Chrs
      regionPvalsDiffFromExp[[region]]$greaterThanMatchedBs     <- greaterThanOtherBs
      regionPvalsDiffFromExp[[region]]$greaterThanIntBs     <- greaterThanIntBs
      
    }
    
    #add these things to the segCCFvalues
    # first, do this for Avals
    if(length(unique(obsAvals))==1)
    {
      # we can determine the low and high values easily 
      nA1  <- 1 # (this should really be ploidy)
      nA2  <- unique(obsAvals)
      
      segCCFvalues[[segName]]$nA1   <- nA1
      segCCFvalues[[segName]]$nA2   <- nA2
      
      # now go region by region, and check the fraction of each (which fraction to count as the CCF depends on whether it is a gain or a loss - determine this below)
      for (region in regions)
      {
        nAobs           <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean) #copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean
        nA1             <- nA1#floor(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean))
        nA2             <- nA2#ceiling(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean))
        
        fracsToUse      <- getCCFractionNeeded(nAobs,nA1,nA2)
        fracA1  <-       0
        fracA2  <-       1
        
        segCCFvalues[[segName]][[region]]$fracA1_mean  <- fracA1
        segCCFvalues[[segName]][[region]]$fracA2_mean  <- fracA2
        segCCFvalues[[segName]][[region]]$nAobs_mean   <- nAobs
        #segCCFvalues[[segName]][[region]]$CCF_mean     <- fracA1
        
        # lets also add the raw values (just in case)
        nAobs_all       <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nA<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nA) #copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean
        nAobs_all       <- ifelse(nAobs_all>=maxCPN,maxCPN,nAobs_all)
        nA1             <- nA1#floor(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean))
        nA2             <- nA2#ceiling(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean))
        fracsToUse     <- sapply(nBobs_all,getCCFractionNeeded,nA1,nA2)
        fracsToUse[fracsToUse=="Inf"] <- 1
        
        fracA1_all      <-apply(fracsToUse,2,max)
        # fracifelse(as.numeric(nB2)-as.numeric(nBobs_all)<0
        #  ,abs(as.numeric(nB2)-as.numeric(nBobs_all))+1
        #  ,as.numeric(nB2)-as.numeric(nBobs_all))
        fracA2_all       <-apply(fracsToUse,2,min)
        
        
        #tTestWithErrorHandling(fracA1_all,mu = 1)
        
        segCCFvalues[[segName]][[region]]$fracA1_all  <- fracA1_all
        segCCFvalues[[segName]][[region]]$fracA2_all  <- fracA2_all
        segCCFvalues[[segName]][[region]]$nAobs_all   <- nAobs_all
        #segCCFvalues[[segName]][[region]]$fracA2_all  <- fracA2
        
        
      }
      
      #next;
    }
    
    # let's check whether the ceiling and floor lead to only 2 options
    if(length(unique(c(ceiling(obsAvals),floor(obsAvals))))==2)
    {
      # first check whether all regions are clonal
      # check whether only one state when rounding
      if(length(unique(round(obsAvals)))==1)
      {
        
        # check whether all regions can be classified as clonal
        roundedIntger <- round(obsAvals)[1]
        allClonal <- TRUE
        nonClonalRegion <- c()
        for (region in regions)
        {
          Test <- tTestWithErrorHandling(copyITHRDS$regionComparisonList[[region]][[segName]]$nA,mu = roundedIntger,conf.level = confInt)
          if(Test$p.value<(1-confInt))
          {
            allClonal <- FALSE
            nonClonalRegion <- c(nonClonalRegion,region)
          }
        }
        
        
        ttestClonal <-TRUE
        nonClonalRegionBAF <- c()
        for (region in regions)
        {
          
          segSplit <- unlist(strsplit(segName,split=":"))
          test <- copyITHRDS$newBAFtable[[region]]$regionBAFtable[copyITHRDS$newBAFtable[[region]]$regionBAFtable$pos>=as.numeric(segSplit[2])&copyITHRDS$newBAFtable[[region]]$regionBAFtable$pos<=as.numeric(segSplit[3])&copyITHRDS$newBAFtable[[region]]$regionBAFtable$chrom==as.numeric(segSplit[1]),,drop=FALSE]
          nAraw   <- round(mean(copyITHRDS$regionComparisonList[[region]][[segName]]$nA))
          nBraw   <- round(mean(copyITHRDS$regionComparisonList[[region]][[segName]]$nB))
          
          #get the SD of the sample
          # SD <- sd(c(test[test$nA==1,'BAF'],1-test[test$nB==1,'BAF']))
          # Number <- length(c(test[test$nA==1,'BAF'],1-test[test$nB==1,'BAF']))
          # mu <- funCalcExpectedBAF(expBvals[region],expAvals[region],ploidyEstimates[region],purityEstimates[region])
          # #tTestWithErrorHandling(c(test[test$nA==1,'BAF'],1-test[test$nB==1,'BAF']),rnorm(n = Number,sd = 0.5,mean=mu))
          #tTestWithErrorHandling(c(test[test$nA==1,'BAF'],1-test[test$nB==1,'BAF']),mu = funCalcExpectedBAF(expBvals[region],expAvals[region],ploidyEstimates[region],purityEstimates[region]))
          tTestSpec <- tTestWithErrorHandling(c(test[test$nA==1,'BAF'],1-test[test$nB==1,'BAF']),mu = funCalcExpectedBAF(expBvals[region],expAvals[region],ploidyEstimates[region],purityEstimates[region]))
          #tTest <- c(tTest,tTestSpec$p.value)
          if(tTestSpec$p.value<(1-confInt))
          {
            ttestClonal <- FALSE
            nonClonalRegionBAF <- c(nonClonalRegionBAF,region)
          }
          
          
        }
        
        
        # finally, check what the CCF is of each region if you assume there is subclonality
        maxCCF <- 0
        nA1  <- min(unique(floor(obsAvals)))
        nA2  <- max(unique(ceiling(obsAvals)))
        
        # now go region by region, and check the fraction of each (which fraction to count as the CCF depends on whether it is a gain or a loss - determine this below)
        for (region in regions)
        {
          nAobs           <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean) #copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean
          nA1             <- nA1#floor(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean))
          nA2             <- nA2#ceiling(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean))
          
          fracsToUse      <- getCCFractionNeeded(nAobs,nA1,nA2)
          minCCF <- min(fracsToUse)
          if(minCCF>maxCCF)
          {
            maxCCF <- minCCF
          }
        }
        
        # next, determine whether to explore subclonality
        if(allClonal==TRUE|ttestClonal==TRUE|maxCCF<minCCFtoConsider)
        {
          allClonal <- TRUE
        }
        
        
        if(allClonal)
        {
          nA1  <- round(obsAvals)[1]
          nA2  <- round(obsAvals)[1]
          
          segCCFvalues[[segName]]$nA1   <- nA1
          segCCFvalues[[segName]]$nA2   <- nA2
          
          # now go region by region, and check the fraction of each (which fraction to count as the CCF depends on whether it is a gain or a loss - determine this below)
          for (region in regions)
          {
            nAobs           <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean) #copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean
            nA1             <- nA1#floor(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean))
            nA2             <- nA2#ceiling(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean))
            fracsToUse      <- getCCFractionNeeded(nAobs,nA1,nA2)
            
            fracA1  <- max(fracsToUse)
            # ifelse(as.numeric(nA2)-as.numeric(nAobs)<0
            #   ,abs(as.numeric(nA2)-as.numeric(nAobs))+1
            #   ,as.numeric(nA2)-as.numeric(nAobs))
            fracA2  <- min(fracsToUse)
            
            segCCFvalues[[segName]][[region]]$fracA1_mean  <- fracA1
            segCCFvalues[[segName]][[region]]$fracA2_mean  <- fracA2
            segCCFvalues[[segName]][[region]]$nAobs_mean   <- nAobs
            #segCCFvalues[[segName]][[region]]$CCF_mean     <- fracA1
            
            # lets also add the raw values (just in case)
            nAobs_all      <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nA<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nA) #copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean
            nAobs_all       <- ifelse(nAobs_all>=maxCPN,maxCPN,nAobs_all)
            nA1             <- nA1#floor(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean))
            nA2             <- nA2#ceiling(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean))
            fracsToUse     <- sapply(nAobs_all,getCCFractionNeeded,nA1,nA2)
            fracsToUse[fracsToUse=="Inf"] <- 1
            
            
            fracA1_all      <-apply(fracsToUse,2,max)
            # fracifelse(as.numeric(nA2)-as.numeric(nAobs_all)<0
            #  ,abs(as.numeric(nA2)-as.numeric(nAobs_all))+1
            #  ,as.numeric(nA2)-as.numeric(nAobs_all))
            fracA2_all       <-apply(fracsToUse,2,min)
            segCCFvalues[[segName]][[region]]$fracA1_all  <- fracA1_all
            segCCFvalues[[segName]][[region]]$fracA2_all  <- fracA2_all
            segCCFvalues[[segName]][[region]]$nAobs_all   <- nAobs_all
            #segCCFvalues[[segName]][[region]]$fracA2_all  <- fracA2
            
            
          }
        }
        
        if(!allClonal)
        {
          # we can determine the low and high values easily 
          nA1  <- min(unique(floor(obsAvals)))
          nA2  <- max(unique(ceiling(obsAvals)))
          
          segCCFvalues[[segName]]$nA1   <- nA1
          segCCFvalues[[segName]]$nA2   <- nA2
          
          # now go region by region, and check the fraction of each (which fraction to count as the CCF depends on whether it is a gain or a loss - determine this below)
          for (region in regions)
          {
            nAobs           <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean) #copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean
            nA1             <- nA1#floor(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean))
            nA2             <- nA2#ceiling(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean))
            
            fracsToUse      <- getCCFractionNeeded(nAobs,nA1,nA2)
            fracA1  <-       fracsToUse[1]#ifelse(as.numeric(nA2)-as.numeric(nAobs)<0
            # ,abs(as.numeric(nA2)-as.numeric(nAobs))+1
            # ,as.numeric(nA2)-as.numeric(nAobs))
            fracA2  <- fracsToUse[2]#ifelse(1-fracA1<0,0,1-fracA1)
            
            segCCFvalues[[segName]][[region]]$fracA1_mean  <- fracA1
            segCCFvalues[[segName]][[region]]$fracA2_mean  <- fracA2
            segCCFvalues[[segName]][[region]]$nAobs_mean   <- nAobs
            #segCCFvalues[[segName]][[region]]$CCF_mean     <- fracA1
            
            # lets also add the raw values (just in case)
            nAobs_all      <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nA<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nA) #copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean
            nAobs_all       <- ifelse(nAobs_all>=maxCPN,maxCPN,nAobs_all)
            nA1             <- nA1#floor(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean))
            nA2             <- nA2#ceiling(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean))
            fracsToUse     <- sapply(nAobs_all,getCCFractionNeeded,nA1,nA2)
            fracsToUse[fracsToUse=="Inf"] <- 1
            
            
            fracA1_all      <-fracsToUse[1,] 
            # fracifelse(as.numeric(nA2)-as.numeric(nAobs_all)<0
            #  ,abs(as.numeric(nA2)-as.numeric(nAobs_all))+1
            #  ,as.numeric(nA2)-as.numeric(nAobs_all))
            fracA2_all       <-fracsToUse[2,] 
            
            #tTestWithErrorHandling(fracA1_all,mu = 1)
            
            segCCFvalues[[segName]][[region]]$fracA1_all  <- fracA1_all
            segCCFvalues[[segName]][[region]]$fracA2_all  <- fracA2_all
            segCCFvalues[[segName]][[region]]$nAobs_all   <- nAobs_all
            #segCCFvalues[[segName]][[region]]$fracA2_all  <- fracA2
            
            
          }
        }
        
        
        
      }
      
      if(length(unique(round(obsAvals)))!=1)
      {
        # we can determine the low and high values easily 
        nA1  <- min(unique(floor(obsAvals)))
        nA2  <- max(unique(ceiling(obsAvals)))
        
        segCCFvalues[[segName]]$nA1   <- nA1
        segCCFvalues[[segName]]$nA2   <- nA2
        
        # now go region by region, and check the fraction of each (which fraction to count as the CCF depends on whether it is a gain or a loss - determine this below)
        for (region in regions)
        {
          nAobs           <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean) #copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean
          nA1             <- nA1#floor(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean))
          nA2             <- nA2#ceiling(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean))
          
          fracsToUse      <- getCCFractionNeeded(nAobs,nA1,nA2)
          fracA1  <-       fracsToUse[1]#ifelse(as.numeric(nA2)-as.numeric(nAobs)<0
          # ,abs(as.numeric(nA2)-as.numeric(nAobs))+1
          # ,as.numeric(nA2)-as.numeric(nAobs))
          fracA2  <- fracsToUse[2]#ifelse(1-fracA1<0,0,1-fracA1)
          
          segCCFvalues[[segName]][[region]]$fracA1_mean  <- fracA1
          segCCFvalues[[segName]][[region]]$fracA2_mean  <- fracA2
          segCCFvalues[[segName]][[region]]$nAobs_mean   <- nAobs
          #segCCFvalues[[segName]][[region]]$CCF_mean     <- fracA1
          
          # lets also add the raw values (just in case)
          nAobs_all      <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nA<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nA) #copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean
          nAobs_all       <- ifelse(nAobs_all>=maxCPN,maxCPN,nAobs_all)
          nA1             <- nA1#floor(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean))
          nA2             <- nA2#ceiling(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean))
          fracsToUse     <- sapply(nAobs_all,getCCFractionNeeded,nA1,nA2)
          fracsToUse[fracsToUse=="Inf"] <- 1
          
          
          fracA1_all      <-fracsToUse[1,] 
          # fracifelse(as.numeric(nA2)-as.numeric(nAobs_all)<0
          #  ,abs(as.numeric(nA2)-as.numeric(nAobs_all))+1
          #  ,as.numeric(nA2)-as.numeric(nAobs_all))
          fracA2_all       <-fracsToUse[2,] 
          
          #tTestWithErrorHandling(fracA1_all,mu = 1)
          
          segCCFvalues[[segName]][[region]]$fracA1_all  <- fracA1_all
          segCCFvalues[[segName]][[region]]$fracA2_all  <- fracA2_all
          segCCFvalues[[segName]][[region]]$nAobs_all   <- nAobs_all
          #segCCFvalues[[segName]][[region]]$fracA2_all  <- fracA2
          
          
        }
      }
      #next;
    }
    
    # let's check whether the ceiling and floor lead to only 2 options
    if(length(unique(c(ceiling(obsAvals),floor(obsAvals))))>2)
    {
      # check whether only one state when rounding
      if(length(unique(round(obsAvals)))==1)
      {
        
        # check whether all regions can be classified as clonal
        roundedIntger <- round(obsAvals)[1]
        allClonal <- TRUE
        nonClonalRegion <- c()
        for (region in regions)
        {
          Test <- tTestWithErrorHandling(copyITHRDS$regionComparisonList[[region]][[segName]]$nA,mu = roundedIntger,conf.level = confInt)
          if(Test$p.value<(1-confInt))
          {
            allClonal <- FALSE
            nonClonalRegion <- c(nonClonalRegion,region)
          }
        }
        
        ttestClonal <-TRUE
        nonClonalRegionBAF <- c()
        for (region in regions)
        {
          
          segSplit <- unlist(strsplit(segName,split=":"))
          test <- copyITHRDS$newBAFtable[[region]]$regionBAFtable[copyITHRDS$newBAFtable[[region]]$regionBAFtable$pos>=as.numeric(segSplit[2])&copyITHRDS$newBAFtable[[region]]$regionBAFtable$pos<=as.numeric(segSplit[3])&copyITHRDS$newBAFtable[[region]]$regionBAFtable$chrom==as.numeric(segSplit[1]),,drop=FALSE]
          nAraw   <- round(mean(copyITHRDS$regionComparisonList[[region]][[segName]]$nA))
          nBraw   <- round(mean(copyITHRDS$regionComparisonList[[region]][[segName]]$nB))
          
          tTestSpec <- tTestWithErrorHandling(c(test[test$nA==1,'BAF'],1-test[test$nB==1,'BAF']),mu = funCalcExpectedBAF(expBvals[region],expAvals[region],ploidyEstimates[region],purityEstimates[region]))
          #tTest <- c(tTest,tTestSpec$p.value)
          if(tTestSpec$p.value<(1-confInt))
          {
            ttestClonal <- FALSE
            nonClonalRegionBAF <- c(nonClonalRegionBAF,region)
          }
          
          
        }
        
        
        # finally, check what the CCF is of each region if you assume there is subclonality
        maxCCF <- 0
        nA1  <- min(unique(floor(obsAvals)))
        nA2  <- max(unique(round(obsAvals)))
        if(abs(min(max(ceiling(obsAvals))-obsAvals))<abs(max(min(floor(obsAvals))-obsAvals)))
        {
          nA1  <- round(obsAvals)[1]
          nA2  <- max(unique(ceiling(obsAvals)))
        }
        
        if(abs(min(max(ceiling(obsAvals))-obsAvals))>=abs(max(min(floor(obsAvals))-obsAvals)))
        {
          nA1  <- min(unique(floor(obsAvals)))
          nA2  <- max(unique(round(obsAvals)))
        }
        # now go region by region, and check the fraction of each (which fraction to count as the CCF depends on whether it is a gain or a loss - determine this below)
        for (region in regions)
        {
          nAobs           <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean) #copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean
          nA1             <- nA1#floor(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean))
          nA2             <- nA2#ceiling(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean))
          
          fracsToUse      <- getCCFractionNeeded(nAobs,nA1,nA2)
          minCCF <- min(fracsToUse)
          if(minCCF>maxCCF)
          {
            maxCCF <- minCCF
          }
        }
        
        # next, determine whether to explore subclonality
        if(allClonal==TRUE|ttestClonal==TRUE|maxCCF<minCCFtoConsider)
        {
          allClonal <- TRUE
        }
        
        
        if(allClonal)
        {
          nA1  <- round(obsAvals)[1]
          nA2  <- round(obsAvals)[1]
          
          segCCFvalues[[segName]]$nA1   <- nA1
          segCCFvalues[[segName]]$nA2   <- nA2
          
          # now go region by region, and check the fraction of each (which fraction to count as the CCF depends on whether it is a gain or a loss - determine this below)
          for (region in regions)
          {
            nAobs           <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean) #copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean
            nA1             <- nA1#floor(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean))
            nA2             <- nA2#ceiling(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean))
            fracsToUse      <- getCCFractionNeeded(nAobs,nA1,nA2)
            
            fracA1  <- max(fracsToUse)
            # ifelse(as.numeric(nA2)-as.numeric(nAobs)<0
            #   ,abs(as.numeric(nA2)-as.numeric(nAobs))+1
            #   ,as.numeric(nA2)-as.numeric(nAobs))
            fracA2  <- min(fracsToUse)
            
            segCCFvalues[[segName]][[region]]$fracA1_mean  <- fracA1
            segCCFvalues[[segName]][[region]]$fracA2_mean  <- fracA2
            segCCFvalues[[segName]][[region]]$nAobs_mean   <- nAobs
            #segCCFvalues[[segName]][[region]]$CCF_mean     <- fracA1
            
            # lets also add the raw values (just in case)
            nAobs_all      <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nA<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nA) #copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean
            nAobs_all       <- ifelse(nAobs_all>=maxCPN,maxCPN,nAobs_all)
            nA1             <- nA1#floor(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean))
            nA2             <- nA2#ceiling(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean))
            fracsToUse     <- sapply(nAobs_all,getCCFractionNeeded,nA1,nA2)
            fracsToUse[fracsToUse=="Inf"] <- 1
            
            
            fracA1_all      <-apply(fracsToUse,2,max)
            # fracifelse(as.numeric(nA2)-as.numeric(nAobs_all)<0
            #  ,abs(as.numeric(nA2)-as.numeric(nAobs_all))+1
            #  ,as.numeric(nA2)-as.numeric(nAobs_all))
            fracA2_all       <-apply(fracsToUse,2,min)
            segCCFvalues[[segName]][[region]]$fracA1_all  <- fracA1_all
            segCCFvalues[[segName]][[region]]$fracA2_all  <- fracA2_all
            segCCFvalues[[segName]][[region]]$nAobs_all   <- nAobs_all
            #segCCFvalues[[segName]][[region]]$fracA2_all  <- fracA2
            
            
          }
        }
        
        if(!allClonal)
        {
          # should we consider 1 value difference between the two or is it worth considering more
          
          
          twoValueDiff <- TRUE
          if(abs(min(max(ceiling(obsAvals))-obsAvals))>=abs(max(min(floor(obsAvals))-obsAvals)))
          {
            nA1  <- min(unique(floor(obsAvals)))
            nA2  <- max(unique(round(obsAvals)))
            
            # is this choice of nA2, too conservative
            regionsWhichMightBeGreater <- names(which(obsAvals>=nA2))
            for (region in regionsWhichMightBeGreater)
            {
             
              maxCCF      <-  abs(nA2 - obsAvals[region])
              diffFromInt <- (regionPvalsDiffFromExp[[region]]$nApval_fromInteger<1-confInt)
              diffFromOther <- TRUE
              if (length(regionPvalsDiffFromExp[[region]]$nA_nrMatchedAChrs)>=minOtherChrs)
              {
                diffFromOther <- (regionPvalsDiffFromExp[[region]]$nApval_fromOtherChrs<1-confInt)
                
              }
          
              if(diffFromInt==TRUE&diffFromOther==TRUE&maxCCF>=minCCFtoConsider)
              {
                twoValueDiff <- FALSE
              }
              
              
            }
            
          }
          
          if(abs(min(max(ceiling(obsAvals))-obsAvals))<abs(max(min(floor(obsAvals))-obsAvals)))
          {
            nA1  <- round(obsAvals)[1]
            nA2  <- max(unique(ceiling(obsAvals)))
            
            # is this choice of nA1, too conservative
            regionsWhichMightBeLess <- names(which(obsAvals<nA1))
            for (region in regionsWhichMightBeLess)
            {
              
              maxCCF      <-  abs(nA1 - obsAvals[region])
              diffFromInt <- (regionPvalsDiffFromExp[[region]]$nApval_fromInteger<1-confInt)
              diffFromOther <- TRUE
              if (length(regionPvalsDiffFromExp[[region]]$nA_nrMatchedAChrs)>=minOtherChrs)
              {
                diffFromOther <- (regionPvalsDiffFromExp[[region]]$nApval_fromOtherChrs<1-confInt)
                
              }
              
              if(diffFromInt==TRUE&diffFromOther==TRUE&maxCCF>=minCCFtoConsider)
              {
                twoValueDiff <- FALSE
              }
              
              
            }
            
          }

          if(twoValueDiff==TRUE)
          {
            # do we twist up or down
            if(abs(min(max(ceiling(obsAvals))-obsAvals))<abs(max(min(floor(obsAvals))-obsAvals)))
            {
              # take the ceiling 
              # we can determine the low and high values easily 
              nA1  <- round(obsAvals)[1]
              nA2  <- max(unique(ceiling(obsAvals)))
              
              segCCFvalues[[segName]]$nA1   <- nA1
              segCCFvalues[[segName]]$nA2   <- nA2
              
              # now go region by region, and check the fraction of each (which fraction to count as the CCF depends on whether it is a gain or a loss - determine this below)
              for (region in regions)
              {
                nAobs           <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean) #copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean
                nA1             <- nA1#floor(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean))
                nA2             <- nA2#ceiling(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean))
                fracsToUse      <- getCCFractionNeeded(nAobs,nA1,nA2)
                
                fracA1  <- fracsToUse[1]
                # ifelse(as.numeric(nA2)-as.numeric(nAobs)<0
                #   ,abs(as.numeric(nA2)-as.numeric(nAobs))+1
                #   ,as.numeric(nA2)-as.numeric(nAobs))
                fracA2  <- fracsToUse[2]
                
                segCCFvalues[[segName]][[region]]$fracA1_mean  <- fracA1
                segCCFvalues[[segName]][[region]]$fracA2_mean  <- fracA2
                segCCFvalues[[segName]][[region]]$nAobs_mean   <- nAobs
                #segCCFvalues[[segName]][[region]]$CCF_mean     <- fracA1
                
                # lets also add the raw values (just in case)
                nAobs_all      <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nA<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nA) #copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean
                nAobs_all       <- ifelse(nAobs_all>=maxCPN,maxCPN,nAobs_all)
                nA1             <- nA1#floor(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean))
                nA2             <- nA2#ceiling(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean))
                fracsToUse     <- sapply(nAobs_all,getCCFractionNeeded,nA1,nA2)
                fracsToUse[fracsToUse=="Inf"] <- 1
                
                
                fracA1_all      <-fracsToUse[1,] 
                # fracifelse(as.numeric(nA2)-as.numeric(nAobs_all)<0
                #  ,abs(as.numeric(nA2)-as.numeric(nAobs_all))+1
                #  ,as.numeric(nA2)-as.numeric(nAobs_all))
                fracA2_all       <-fracsToUse[2,] 
                segCCFvalues[[segName]][[region]]$fracA1_all  <- fracA1_all
                segCCFvalues[[segName]][[region]]$fracA2_all  <- fracA2_all
                segCCFvalues[[segName]][[region]]$nAobs_all   <- nAobs_all
                #segCCFvalues[[segName]][[region]]$fracA2_all  <- fracA2
                
                
              }
              
              #next;
              
            }
            if(abs(min(max(ceiling(obsAvals))-obsAvals))>=abs(max(min(floor(obsAvals))-obsAvals)))
            {
              # take the floor 
              # we can determine the low and high values easily 
              nA1  <- min(unique(floor(obsAvals)))
              nA2  <- max(unique(round(obsAvals)))
              
              segCCFvalues[[segName]]$nA1   <- nA1
              segCCFvalues[[segName]]$nA2   <- nA2
              
              # now go region by region, and check the fraction of each (which fraction to count as the CCF depends on whether it is a gain or a loss - determine this below)
              for (region in regions)
              {
                nAobs           <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean) #copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean
                nA1             <- nA1#floor(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean))
                nA2             <- nA2#ceiling(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean))
                fracsToUse      <- getCCFractionNeeded(nAobs,nA1,nA2)
                
                fracA1  <- fracsToUse[1]
                # ifelse(as.numeric(nA2)-as.numeric(nAobs)<0
                #   ,abs(as.numeric(nA2)-as.numeric(nAobs))+1
                #   ,as.numeric(nA2)-as.numeric(nAobs))
                fracA2  <- fracsToUse[2]
                #ifelse(1-fracA1<0,0,1-fracA1)
                
                segCCFvalues[[segName]][[region]]$fracA1_mean  <- fracA1
                segCCFvalues[[segName]][[region]]$fracA2_mean  <- fracA2
                segCCFvalues[[segName]][[region]]$nAobs_mean   <- nAobs
                #segCCFvalues[[segName]][[region]]$CCF_mean     <- fracA1
                
                # lets also add the raw values (just in case)
                nAobs_all      <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nA<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nA) #copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean
                nAobs_all       <- ifelse(nAobs_all>=maxCPN,maxCPN,nAobs_all)
                nA1             <- nA1#floor(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean))
                nA2             <- nA2#ceiling(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean))
                fracsToUse     <- sapply(nAobs_all,getCCFractionNeeded,nA1,nA2)
                fracsToUse[fracsToUse=="Inf"] <- 1
                
                
                fracA1_all      <-fracsToUse[1,] 
                # fracifelse(as.numeric(nA2)-as.numeric(nAobs_all)<0
                #  ,abs(as.numeric(nA2)-as.numeric(nAobs_all))+1
                #  ,as.numeric(nA2)-as.numeric(nAobs_all))
                fracA2_all       <-fracsToUse[2,] 
                
                segCCFvalues[[segName]][[region]]$fracA1_all  <- fracA1_all
                segCCFvalues[[segName]][[region]]$fracA2_all  <- fracA2_all
                segCCFvalues[[segName]][[region]]$nAobs_all   <- nAobs_all
                #segCCFvalues[[segName]][[region]]$fracA2_all  <- fracA2
                
                
              }
              
              #next;
              
            }
          }
          
          if(twoValueDiff==FALSE)
          {
            # do we twist up or down
            if(abs(min(max(ceiling(obsAvals))-obsAvals))<abs(max(min(floor(obsAvals))-obsAvals)))
            {
              # take the ceiling 
              # we can determine the low and high values easily 
              nA1  <- min(unique(floor(obsAvals)))#[1]
              nA2  <- max(unique(ceiling(obsAvals)))
              
              segCCFvalues[[segName]]$nA1   <- nA1
              segCCFvalues[[segName]]$nA2   <- nA2
              
              # now go region by region, and check the fraction of each (which fraction to count as the CCF depends on whether it is a gain or a loss - determine this below)
              for (region in regions)
              {
                nAobs           <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean) #copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean
                nA1             <- nA1#floor(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean))
                nA2             <- nA2#ceiling(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean))
                fracsToUse      <- getCCFractionNeeded(nAobs,nA1,nA2)
                
                fracA1  <- fracsToUse[1]
                # ifelse(as.numeric(nA2)-as.numeric(nAobs)<0
                #   ,abs(as.numeric(nA2)-as.numeric(nAobs))+1
                #   ,as.numeric(nA2)-as.numeric(nAobs))
                fracA2  <- fracsToUse[2]
                
                segCCFvalues[[segName]][[region]]$fracA1_mean  <- fracA1
                segCCFvalues[[segName]][[region]]$fracA2_mean  <- fracA2
                segCCFvalues[[segName]][[region]]$nAobs_mean   <- nAobs
                #segCCFvalues[[segName]][[region]]$CCF_mean     <- fracA1
                
                # lets also add the raw values (just in case)
                nAobs_all      <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nA<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nA) #copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean
                nAobs_all       <- ifelse(nAobs_all>=maxCPN,maxCPN,nAobs_all)
                nA1             <- nA1#floor(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean))
                nA2             <- nA2#ceiling(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean))
                fracsToUse     <- sapply(nAobs_all,getCCFractionNeeded,nA1,nA2)
                fracsToUse[fracsToUse=="Inf"] <- 1
                
                
                fracA1_all      <-fracsToUse[1,] 
                # fracifelse(as.numeric(nA2)-as.numeric(nAobs_all)<0
                #  ,abs(as.numeric(nA2)-as.numeric(nAobs_all))+1
                #  ,as.numeric(nA2)-as.numeric(nAobs_all))
                fracA2_all       <-fracsToUse[2,] 
                segCCFvalues[[segName]][[region]]$fracA1_all  <- fracA1_all
                segCCFvalues[[segName]][[region]]$fracA2_all  <- fracA2_all
                segCCFvalues[[segName]][[region]]$nAobs_all   <- nAobs_all
                #segCCFvalues[[segName]][[region]]$fracA2_all  <- fracA2
                
                
              }
              
              #next;
              
            }
            if(abs(min(max(ceiling(obsAvals))-obsAvals))>=abs(max(min(floor(obsAvals))-obsAvals)))
            {
              # take the floor 
              # we can determine the low and high values easily 
              nA1  <- min(unique(floor(obsAvals)))#[1]
              nA2  <- max(unique(ceiling(obsAvals)))
              
              
              segCCFvalues[[segName]]$nA1   <- nA1
              segCCFvalues[[segName]]$nA2   <- nA2
              
              # now go region by region, and check the fraction of each (which fraction to count as the CCF depends on whether it is a gain or a loss - determine this below)
              for (region in regions)
              {
                nAobs           <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean) #copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean
                nA1             <- nA1#floor(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean))
                nA2             <- nA2#ceiling(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean))
                fracsToUse      <- getCCFractionNeeded(nAobs,nA1,nA2)
                
                fracA1  <- fracsToUse[1]
                # ifelse(as.numeric(nA2)-as.numeric(nAobs)<0
                #   ,abs(as.numeric(nA2)-as.numeric(nAobs))+1
                #   ,as.numeric(nA2)-as.numeric(nAobs))
                fracA2  <- fracsToUse[2]
                #ifelse(1-fracA1<0,0,1-fracA1)
                
                segCCFvalues[[segName]][[region]]$fracA1_mean  <- fracA1
                segCCFvalues[[segName]][[region]]$fracA2_mean  <- fracA2
                segCCFvalues[[segName]][[region]]$nAobs_mean   <- nAobs
                #segCCFvalues[[segName]][[region]]$CCF_mean     <- fracA1
                
                # lets also add the raw values (just in case)
                nAobs_all      <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nA<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nA) #copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean
                nAobs_all       <- ifelse(nAobs_all>=maxCPN,maxCPN,nAobs_all)
                nA1             <- nA1#floor(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean))
                nA2             <- nA2#ceiling(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean))
                fracsToUse     <- sapply(nAobs_all,getCCFractionNeeded,nA1,nA2)
                fracsToUse[fracsToUse=="Inf"] <- 1
                
                
                fracA1_all      <-fracsToUse[1,] 
                # fracifelse(as.numeric(nA2)-as.numeric(nAobs_all)<0
                #  ,abs(as.numeric(nA2)-as.numeric(nAobs_all))+1
                #  ,as.numeric(nA2)-as.numeric(nAobs_all))
                fracA2_all       <-fracsToUse[2,] 
                
                segCCFvalues[[segName]][[region]]$fracA1_all  <- fracA1_all
                segCCFvalues[[segName]][[region]]$fracA2_all  <- fracA2_all
                segCCFvalues[[segName]][[region]]$nAobs_all   <- nAobs_all
                #segCCFvalues[[segName]][[region]]$fracA2_all  <- fracA2
                
                
              }
              
              #next;
              
            }
          }

        }
        
        
        
      }
      
      #check whether there are only two states when rounding
      if(length(unique(round(obsAvals)))==2)
      {
        
        # we can determine the low and high values easily
        nA1  <- min(round(obsAvals))
        nA2  <- max(round(obsAvals))
        
        segCCFvalues[[segName]]$nA1   <- nA1
        segCCFvalues[[segName]]$nA2   <- nA2
        
        # now go region by region, and check the fraction of each (which fraction to count as the CCF depends on whether it is a gain or a loss - determine this below)
        for (region in regions)
        {
          nAobs           <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean) #copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean
          nA1             <- nA1#floor(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean))
          nA2             <- nA2#ceiling(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean))
          fracsToUse      <- getCCFractionNeeded(nAobs,nA1,nA2)
          fracA1  <- fracsToUse[1] 
          # ifelse(as.numeric(nA2)-as.numeric(nAobs)<0
          # ,abs(as.numeric(nA2)-as.numeric(nAobs))+1
          # ,as.numeric(nA2)-as.numeric(nAobs))
          fracA2  <- fracsToUse[2]
          # ifelse(1-fracA1<0,0,1-fracA1)
          
          segCCFvalues[[segName]][[region]]$fracA1_mean  <- fracA1
          segCCFvalues[[segName]][[region]]$fracA2_mean  <- fracA2
          segCCFvalues[[segName]][[region]]$nAobs_mean   <- nAobs
          #segCCFvalues[[segName]][[region]]$CCF_mean     <- fracA1
          
          # lets also add the raw values (just in case)
          nAobs_all      <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nA<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nA) #copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean
          nAobs_all       <- ifelse(nAobs_all>=maxCPN,maxCPN,nAobs_all)
          nA1             <- nA1#floor(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean))
          nA2             <- nA2#ceiling(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean))
          fracsToUse     <- sapply(nAobs_all,getCCFractionNeeded,nA1,nA2)
          fracsToUse[fracsToUse=="Inf"] <- 1
          
          
          fracA1_all      <-fracsToUse[1,] 
          # fracifelse(as.numeric(nA2)-as.numeric(nAobs_all)<0
          #  ,abs(as.numeric(nA2)-as.numeric(nAobs_all))+1
          #  ,as.numeric(nA2)-as.numeric(nAobs_all))
          fracA2_all       <-fracsToUse[2,] 
          # ifelse(1-fracA1_all<0,0,1-fracA1_all)
          
          segCCFvalues[[segName]][[region]]$fracA1_all  <- fracA1_all
          segCCFvalues[[segName]][[region]]$fracA2_all  <- fracA2_all
          segCCFvalues[[segName]][[region]]$nAobs_all   <- nAobs_all
          #segCCFvalues[[segName]][[region]]$fracA2_all  <- fracA2
          
          
        }
        
        #next;
      }
      
      
      
      
      # check whether an amplification greater than 8
      
      #if(min(expAvals))
      
      # posValues  <- unique(sort(expAvals))
      # lowValue   <- posValues[1]
      # highValue  <- posValues[2]
      # 
      # # choose the max and min that maximise the CCF across the tumour regions (in other words, minimize heterogeneity where possible)
      # posFracA1 <- c()
      # posFracA2 <- c()
      # valuesChosen <- c()
      # for (x in 1:(length(posValues)-1))
      # {
      #   for (y in (x+1):length(posValues))
      #   {
      #     
      #     posLowValue  <- posValues[x]
      #     posHighValue <- posValues[y]
      #     
      #     CCFvalues <- c()
      #     for (i in 1:length(obsAvals))
      #     {
      #       CCFvalues <-c(CCFvalues,(getCCFractionNeeded(obsAvals[i],posLowValue,posHighValue))[2])
      #     }
      #     valuesChosen <- rbind(valuesChosen,c(posLowValue,posHighValue))
      #     posFracA1 <- rbind(posFracA1,CCFvalues)
      #   }
      # }
      
      if(length(unique(round(obsAvals)))>2)
      {
        # we can determine the low and high values easily 
        nA1  <- min(unique(floor(obsAvals+minCCFtoConsider)))
        nA2  <- max(unique(ceiling(obsAvals-minCCFtoConsider)))
        
        segCCFvalues[[segName]]$nA1   <- nA1
        segCCFvalues[[segName]]$nA2   <- nA2
        
        # now go region by region, and check the fraction of each (which fraction to count as the CCF depends on whether it is a gain or a loss - determine this below)
        for (region in regions)
        {
          nAobs           <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean) #copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean
          nA1             <- nA1#floor(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean))
          nA2             <- nA2#ceiling(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean))
          
          fracsToUse      <- getCCFractionNeeded(nAobs,nA1,nA2)
          fracA1  <- fracsToUse[1] 
          # ifelse(as.numeric(nA2)-as.numeric(nAobs)<0
          # ,abs(as.numeric(nA2)-as.numeric(nAobs))+1
          # ,as.numeric(nA2)-as.numeric(nAobs))
          fracA2  <- fracsToUse[2]
          # ifelse(1-fracA1<0,0,1-fracA1)
          
          segCCFvalues[[segName]][[region]]$fracA1_mean  <- fracA1
          segCCFvalues[[segName]][[region]]$fracA2_mean  <- fracA2
          segCCFvalues[[segName]][[region]]$nAobs_mean   <- nAobs
          #segCCFvalues[[segName]][[region]]$CCF_mean     <- fracA1
          
          # lets also add the raw values (just in case)
          nAobs_all      <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nA<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nA) #copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean
          nAobs_all       <- ifelse(nAobs_all>=maxCPN,maxCPN,nAobs_all)
          nA1             <- nA1#floor(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean))
          nA2             <- nA2#ceiling(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nA_mean))
          fracsToUse     <- sapply(nAobs_all,getCCFractionNeeded,nA1,nA2)
          fracsToUse[fracsToUse=="Inf"] <- 1
          
          fracA1_all      <-as.numeric(fracsToUse[1,])
          # fracifelse(as.numeric(nB2)-as.numeric(nBobs_all)<0
          #  ,abs(as.numeric(nB2)-as.numeric(nBobs_all))+1
          #  ,as.numeric(nB2)-as.numeric(nBobs_all))
          fracA2_all      <-as.numeric(fracsToUse[2,])
          
          segCCFvalues[[segName]][[region]]$fracA1_all  <- fracA1_all
          segCCFvalues[[segName]][[region]]$fracA2_all  <- fracA2_all
          segCCFvalues[[segName]][[region]]$nAobs_all   <- nAobs_all
          #segCCFvalues[[segName]][[region]]$fracA2_all  <- fracA2
          
          
        }
      }
      
      
      
      #next;
    }
    
    # next, do this for Bvals
    if(length(unique(obsBvals))==1)
    {
      # we can determine the low and high values easily 
      nB1  <- 1 # (this should really be ploidy)
      nB2  <- unique(obsBvals)
      
      segCCFvalues[[segName]]$nB1   <- nB1
      segCCFvalues[[segName]]$nB2   <- nB2
      
      # now go region by region, and check the fraction of each (which fraction to count as the CCF depends on whether it is a gain or a loss - determine this below)
      for (region in regions)
      {
        nBobs           <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean) #copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean
        nB1             <- nB1#floor(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean))
        nB2             <- nB2#ceiling(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean))
        
        fracsToUse      <- getCCFractionNeeded(nBobs,nB1,nB2)
        fracB1  <-       0
        fracB2  <-       1
        
        segCCFvalues[[segName]][[region]]$fracB1_mean  <- fracB1
        segCCFvalues[[segName]][[region]]$fracB2_mean  <- fracB2
        segCCFvalues[[segName]][[region]]$nBobs_mean   <- nBobs
        #segCCFvalues[[segName]][[region]]$CCF_mean     <- fracB1
        
        # lets also add the raw values (just in case)
        nBobs_all      <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nB<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nB) #copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean
        nBobs_all       <- ifelse(nBobs_all>=maxCPN,maxCPN,nBobs_all)
        
        nB1             <- nB1#floor(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean))
        nB2             <- nB2#ceiling(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean))
        fracsToUse     <- sapply(nBobs_all,getCCFractionNeeded,nB1,nB2)
        fracsToUse[fracsToUse=="Inf"] <- 1
        
        fracB1_all      <-apply(fracsToUse,2,max)
        # fracifelse(as.numeric(nB2)-as.numeric(nBobs_all)<0
        #  ,abs(as.numeric(nB2)-as.numeric(nBobs_all))+1
        #  ,as.numeric(nB2)-as.numeric(nBobs_all))
        fracB2_all       <-apply(fracsToUse,2,min)
        
        # tTestWithErrorHandling(fracB1_all,mu = 1)
        
        segCCFvalues[[segName]][[region]]$fracB1_all  <- fracB1_all
        segCCFvalues[[segName]][[region]]$fracB2_all  <- fracB2_all
        segCCFvalues[[segName]][[region]]$nBobs_all   <- nBobs_all
        #segCCFvalues[[segName]][[region]]$fracB2_all  <- fracB2
        
        
      }
      
      # next;
    }
    
    # let's check whether the ceiling and floor lead to only 2 options
    if(length(unique(c(ceiling(obsBvals),floor(obsBvals))))==2)
    {
      
      if(length(unique(round(obsBvals)))==1)
      {
        
        # check whether all regions can be classified as clonal
        roundedIntger <- round(obsBvals)[1]
        allClonal <- TRUE
        nonClonalRegion <- c()
        for (region in regions)
        {
          Test <- tTestWithErrorHandling(copyITHRDS$regionComparisonList[[region]][[segName]]$nB,mu = roundedIntger,conf.level = confInt)
          if(Test$p.value<(1-confInt))
          {
            allClonal <- FALSE
            nonClonalRegion <- c(nonClonalRegion,region)
          }
        }
        
        
        # also, check whether significant different BAF distribution
        
        ttestClonal <-TRUE
        nonClonalRegionBAF <- c()
        for (region in regions)
        {
          
          segSplit <- unlist(strsplit(segName,split=":"))
          test <- copyITHRDS$newBAFtable[[region]]$regionBAFtable[copyITHRDS$newBAFtable[[region]]$regionBAFtable$pos>=as.numeric(segSplit[2])&copyITHRDS$newBAFtable[[region]]$regionBAFtable$pos<=as.numeric(segSplit[3])&copyITHRDS$newBAFtable[[region]]$regionBAFtable$chrom==as.numeric(segSplit[1]),,drop=FALSE]
          #nAraw   <- round(mean(copyITHRDS$regionComparisonList[[region]][[segName]]$nA))
          #nBraw   <- round(mean(copyITHRDS$regionComparisonList[[region]][[segName]]$nB))
          
          tTestSpec <- tTestWithErrorHandling(c(test[test$nA==1,'BAF'],1-test[test$nB==1,'BAF']),mu = funCalcExpectedBAF(expBvals[region],expAvals[region],ploidyEstimates[region],purityEstimates[region]))
          #tTest <- c(tTest,tTestSpec$p.value)
          if(tTestSpec$p.value<(1-confInt))
          {
            ttestClonal <- FALSE
            nonClonalRegionBAF <- c(nonClonalRegionBAF,region)
          }
          
          
        }
        
        
        # finally, check what the CCF is of each region if you assume there is subclonality
        maxCCF <- 0
        nB1  <- min(unique(floor(obsBvals)))
        nB2  <- max(unique(ceiling(obsBvals)))
        
        # now go region by region, and check the fraction of each (which fraction to count as the CCF depends on whether it is a gain or a loss - determine this below)
        for (region in regions)
        {
          nBobs           <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean) #copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean
          nB1             <- nB1#floor(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean))
          nB2             <- nB2#ceiling(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean))
          
          fracsToUse      <- getCCFractionNeeded(nBobs,nB1,nB2)
          minCCF <- min(fracsToUse)
          if(minCCF>maxCCF)
          {
            maxCCF <- minCCF
          }
        }
        
        # next, determine whether to explore subclonality
        if(allClonal==TRUE|ttestClonal==TRUE|maxCCF<minCCFtoConsider)
        {
          allClonal <- TRUE
        }
        
        
        
        if(allClonal)
        {
          nB1  <- round(obsBvals)[1]
          nB2  <- round(obsBvals)[1]
          
          segCCFvalues[[segName]]$nB1   <- nB1
          segCCFvalues[[segName]]$nB2   <- nB2
          
          # now go region by region, and check the fraction of each (which fraction to count as the CCF depends on whether it is a gain or a loss - determine this below)
          for (region in regions)
          {
            nBobs           <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean) #copyITHRDS$regionComparisonList[[region]][[segnBme]]$nB_mean
            nB1             <- nB1#floor(as.numeric(copyITHRDS$regionComparisonList[[region]][[segnBme]]$nB_mean))
            nB2             <- nB2#ceiling(as.numeric(copyITHRDS$regionComparisonList[[region]][[segnBme]]$nB_mean))
            fracsToUse      <- getCCFractionNeeded(nBobs,nB1,nB2)
            
            fracB1  <- max(fracsToUse)
            fracB2  <- min(fracsToUse)
            
            segCCFvalues[[segName]][[region]]$fracB1_mean  <- fracB1
            segCCFvalues[[segName]][[region]]$fracB2_mean  <- fracB2
            segCCFvalues[[segName]][[region]]$nBobs_mean   <- nBobs
            #segCCFvalues[[segnBme]][[region]]$CCF_mean     <- fracB1
            
            # lets also add the raw values (just in case)
            nBobs_all      <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nB<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nB) #copyITHRDS$regionComparisonList[[region]][[segnBme]]$nB_mean
            nBobs_all       <- ifelse(nBobs_all>=maxCPN,maxCPN,nBobs_all)
            nB1             <- nB1#floor(as.numeric(copyITHRDS$regionComparisonList[[region]][[segnBme]]$nB_mean))
            nB2             <- nB2#ceiling(as.numeric(copyITHRDS$regionComparisonList[[region]][[segnBme]]$nB_mean))
            fracsToUse     <- sapply(nBobs_all,getCCFractionNeeded,nB1,nB2)
            fracsToUse[fracsToUse=="Inf"] <- 1
            
            fracB1_all      <-apply(fracsToUse,2,max)
            # fracifelse(as.numeric(nB2)-as.numeric(nBobs_all)<0
            #  ,abs(as.numeric(nB2)-as.numeric(nBobs_all))+1
            #  ,as.numeric(nB2)-as.numeric(nBobs_all))
            fracB2_all       <-apply(fracsToUse,2,min)
            segCCFvalues[[segName]][[region]]$fracB1_all  <- fracB1_all
            segCCFvalues[[segName]][[region]]$fracB2_all  <- fracB2_all
            segCCFvalues[[segName]][[region]]$nBobs_all   <- nBobs_all
            #segCCFvalues[[segnBme]][[region]]$fracB2_all  <- fracB2
            
            
          }
        }
        
        if(!allClonal)
        {
          # we can determine the low and high values easily 
          # we can determine the low and high values easily 
          nB1  <- min(unique(floor(obsBvals)))
          nB2  <- max(unique(ceiling(obsBvals)))
          
          segCCFvalues[[segName]]$nB1   <- nB1
          segCCFvalues[[segName]]$nB2   <- nB2
          
          # now go region by region, and check the fraction of each (which fraction to count as the CCF depends on whether it is a gain or a loss - determine this below)
          for (region in regions)
          {
            nBobs           <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean) #copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean
            nB1             <- nB1#floor(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean))
            nB2             <- nB2#ceiling(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean))
            
            fracsToUse      <- getCCFractionNeeded(nBobs,nB1,nB2)
            fracB1  <-       fracsToUse[1]#ifelse(as.numeric(nB2)-as.numeric(nBobs)<0
            # ,abs(as.numeric(nB2)-as.numeric(nBobs))+1
            # ,as.numeric(nB2)-as.numeric(nBobs))
            fracB2  <- fracsToUse[2]#ifelse(1-fracB1<0,0,1-fracB1)
            
            segCCFvalues[[segName]][[region]]$fracB1_mean  <- fracB1
            segCCFvalues[[segName]][[region]]$fracB2_mean  <- fracB2
            segCCFvalues[[segName]][[region]]$nBobs_mean   <- nBobs
            #segCCFvalues[[segName]][[region]]$CCF_mean     <- fracB1
            
            # lets also add the raw values (just in case)
            nBobs_all      <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nB<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nB) #copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean
            nBobs_all       <- ifelse(nBobs_all>=maxCPN,maxCPN,nBobs_all)
            
            nB1             <- nB1#floor(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean))
            nB2             <- nB2#ceiling(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean))
            fracsToUse     <- sapply(nBobs_all,getCCFractionNeeded,nB1,nB2)
            fracsToUse[fracsToUse=="Inf"] <- 1
            
            
            fracB1_all      <-fracsToUse[1,] 
            # fracifelse(as.numeric(nB2)-as.numeric(nBobs_all)<0
            #  ,abs(as.numeric(nB2)-as.numeric(nBobs_all))+1
            #  ,as.numeric(nB2)-as.numeric(nBobs_all))
            fracB2_all       <-fracsToUse[2,] 
            
            #tTestWithErrorHandling(fracB1_all,mu = 1)
            
            segCCFvalues[[segName]][[region]]$fracB1_all  <- fracB1_all
            segCCFvalues[[segName]][[region]]$fracB2_all  <- fracB2_all
            segCCFvalues[[segName]][[region]]$nBobs_all   <- nBobs_all
            #segCCFvalues[[segName]][[region]]$fracB2_all  <- fracB2
            
            
          }
        }
        
        
        
        
      }
      
      if(length(unique(round(obsBvals)))!=1)
      {
        # we can determine the low and high values easily 
        nB1  <- min(unique(floor(obsBvals)))
        nB2  <- max(unique(ceiling(obsBvals)))
        
        segCCFvalues[[segName]]$nB1   <- nB1
        segCCFvalues[[segName]]$nB2   <- nB2
        
        # now go region by region, and check the fraction of each (which fraction to count as the CCF depends on whether it is a gain or a loss - determine this below)
        for (region in regions)
        {
          nBobs           <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean) #copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean
          nB1             <- nB1#floor(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean))
          nB2             <- nB2#ceiling(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean))
          
          fracsToUse      <- getCCFractionNeeded(nBobs,nB1,nB2)
          fracB1  <-       fracsToUse[1]#ifelse(as.numeric(nB2)-as.numeric(nBobs)<0
          # ,abs(as.numeric(nB2)-as.numeric(nBobs))+1
          # ,as.numeric(nB2)-as.numeric(nBobs))
          fracB2  <- fracsToUse[2]#ifelse(1-fracB1<0,0,1-fracB1)
          
          segCCFvalues[[segName]][[region]]$fracB1_mean  <- fracB1
          segCCFvalues[[segName]][[region]]$fracB2_mean  <- fracB2
          segCCFvalues[[segName]][[region]]$nBobs_mean   <- nBobs
          #segCCFvalues[[segName]][[region]]$CCF_mean     <- fracB1
          
          # lets also add the raw values (just in case)
          nBobs_all      <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nB<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nB) #copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean
          nBobs_all       <- ifelse(nBobs_all>=maxCPN,maxCPN,nBobs_all)
          
          nB1             <- nB1#floor(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean))
          nB2             <- nB2#ceiling(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean))
          fracsToUse     <- sapply(nBobs_all,getCCFractionNeeded,nB1,nB2)
          fracsToUse[fracsToUse=="Inf"] <- 1
          
          
          
          fracB1_all      <-fracsToUse[1,] 
          # fracifelse(as.numeric(nB2)-as.numeric(nBobs_all)<0
          #  ,abs(as.numeric(nB2)-as.numeric(nBobs_all))+1
          #  ,as.numeric(nB2)-as.numeric(nBobs_all))
          fracB2_all       <-fracsToUse[2,] 
          
          #tTestWithErrorHandling(fracB1_all,mu = 1)
          
          segCCFvalues[[segName]][[region]]$fracB1_all  <- fracB1_all
          segCCFvalues[[segName]][[region]]$fracB2_all  <- fracB2_all
          segCCFvalues[[segName]][[region]]$nBobs_all   <- nBobs_all
          #segCCFvalues[[segName]][[region]]$fracB2_all  <- fracB2
          
          
        }
      }
      
      #next;
    }
    
    # let's check whether the ceiling and floor lead to only 2 options
    if(length(unique(c(ceiling(obsBvals),floor(obsBvals))))>2)
    {
      # check whether only one state when rounding
      if(length(unique(round(obsBvals)))==1)
      {
        roundedIntger <- round(obsBvals)[1]
        allClonal <- TRUE
        nonClonalRegion <- c()
        for (region in regions)
        {
          Test <- tTestWithErrorHandling(copyITHRDS$regionComparisonList[[region]][[segName]]$nB,mu = roundedIntger,conf.level = confInt)
          if(Test$p.value<(1-confInt))
          {
            allClonal <- FALSE
            nonClonalRegion <- c(nonClonalRegion,region)
          }
        }
        
        
        ttestClonal <-TRUE
        nonClonalRegionBAF <- c()
        for (region in regions)
        {
          
          segSplit <- unlist(strsplit(segName,split=":"))
          test <- copyITHRDS$newBAFtable[[region]]$regionBAFtable[copyITHRDS$newBAFtable[[region]]$regionBAFtable$pos>=as.numeric(segSplit[2])&copyITHRDS$newBAFtable[[region]]$regionBAFtable$pos<=as.numeric(segSplit[3])&copyITHRDS$newBAFtable[[region]]$regionBAFtable$chrom==as.numeric(segSplit[1]),,drop=FALSE]
          nAraw   <- round(mean(copyITHRDS$regionComparisonList[[region]][[segName]]$nA))
          nBraw   <- round(mean(copyITHRDS$regionComparisonList[[region]][[segName]]$nB))
          
          tTestSpec <- tTestWithErrorHandling(c(test[test$nA==1,'BAF'],1-test[test$nB==1,'BAF']),mu = funCalcExpectedBAF(expBvals[region],expAvals[region],ploidyEstimates[region],purityEstimates[region]))
          #tTest <- c(tTest,tTestSpec$p.value)
          if(tTestSpec$p.value<(1-confInt))
          {
            ttestClonal <- FALSE
            nonClonalRegionBAF <- c(nonClonalRegionBAF,region)
          }
          
          
        }
        
        
        # finally, check what the CCF is of each region if you assume there is subclonality
        maxCCF <- 0
        nB1  <- min(unique(floor(obsBvals)))
        nB2  <- max(unique(ceiling(obsBvals)))
        
        # now go region by region, and check the fraction of each (which fraction to count as the CCF depends on whether it is a gain or a loss - determine this below)
        for (region in regions)
        {
          nBobs           <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean) #copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean
          nB1             <- nB1#floor(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean))
          nB2             <- nB2#ceiling(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean))
          
          fracsToUse      <- getCCFractionNeeded(nBobs,nB1,nB2)
          minCCF <- min(fracsToUse)
          if(minCCF>maxCCF)
          {
            maxCCF <- minCCF
          }
        }
        
        # next, determine whether to explore subclonality
        if(allClonal==TRUE|ttestClonal==TRUE|maxCCF<minCCFtoConsider)
        {
          allClonal <- TRUE
        }
        
        
        if(allClonal)
        {
          nB1  <- round(obsBvals)[1]
          nB2  <- round(obsBvals)[1]
          
          segCCFvalues[[segName]]$nB1   <- nB1
          segCCFvalues[[segName]]$nB2   <- nB2
          
          # now go region by region, and check the fraction of each (which fraction to count as the CCF depends on whether it is a gain or a loss - determine this below)
          for (region in regions)
          {
            nBobs           <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean) #copyITHRDS$regionComparisonList[[region]][[segnBme]]$nB_mean
            nB1             <- nB1#floor(as.numeric(copyITHRDS$regionComparisonList[[region]][[segnBme]]$nB_mean))
            nB2             <- nB2#ceiling(as.numeric(copyITHRDS$regionComparisonList[[region]][[segnBme]]$nB_mean))
            fracsToUse      <- getCCFractionNeeded(nBobs,nB1,nB2)
            
            fracB1  <- max(fracsToUse)
            # ifelse(as.numeric(nB2)-as.numeric(nBobs)<0
            #   ,abs(as.numeric(nB2)-as.numeric(nBobs))+1
            #   ,as.numeric(nB2)-as.numeric(nBobs))
            fracB2  <- min(fracsToUse)
            
            segCCFvalues[[segName]][[region]]$fracB1_mean  <- fracB1
            segCCFvalues[[segName]][[region]]$fracB2_mean  <- fracB2
            segCCFvalues[[segName]][[region]]$nBobs_mean   <- nBobs
            #segCCFvalues[[segnBme]][[region]]$CCF_mean     <- fracB1
            
            # lets also add the raw values (just in case)
            nBobs_all      <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nB<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nB) #copyITHRDS$regionComparisonList[[region]][[segnBme]]$nB_mean
            nBobs_all       <- ifelse(nBobs_all>=maxCPN,maxCPN,nBobs_all)
            nB1             <- nB1#floor(as.numeric(copyITHRDS$regionComparisonList[[region]][[segnBme]]$nB_mean))
            nB2             <- nB2#ceiling(as.numeric(copyITHRDS$regionComparisonList[[region]][[segnBme]]$nB_mean))
            fracsToUse     <- sapply(nBobs_all,getCCFractionNeeded,nB1,nB2)
            fracsToUse[fracsToUse=="Inf"] <- 1
            
            
            fracB1_all      <-apply(fracsToUse,2,max)
            # fracifelse(as.numeric(nB2)-as.numeric(nBobs_all)<0
            #  ,abs(as.numeric(nB2)-as.numeric(nBobs_all))+1
            #  ,as.numeric(nB2)-as.numeric(nBobs_all))
            fracB2_all       <-apply(fracsToUse,2,min)
            segCCFvalues[[segName]][[region]]$fracB1_all  <- fracB1_all
            segCCFvalues[[segName]][[region]]$fracB2_all  <- fracB2_all
            segCCFvalues[[segName]][[region]]$nBobs_all   <- nBobs_all
            #segCCFvalues[[segnBme]][[region]]$fracB2_all  <- fracB2
            
            
          }
        }
        
        if(!allClonal)
        {
          # do we twist up or down
          if(abs(min(max(ceiling(obsBvals))-obsBvals))<abs(max(min(floor(obsBvals))-obsBvals)))
          {
            # take the ceiling 
            # we can determine the low and high values easily 
            nB1  <- round(obsBvals)[1]
            nB2  <- max(unique(ceiling(obsBvals)))
            
            segCCFvalues[[segName]]$nB1   <- nB1
            segCCFvalues[[segName]]$nB2   <- nB2
            
            # now go region by region, and check the fraction of each (which fraction to count as the CCF depends on whether it is a gain or a loss - determine this below)
            for (region in regions)
            {
              nBobs           <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean) #copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean
              nB1             <- nB1#floor(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean))
              nB2             <- nB2#ceiling(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean))
              fracsToUse      <- getCCFractionNeeded(nBobs,nB1,nB2)
              
              fracB1  <- fracsToUse[1]
              # ifelse(as.numeric(nB2)-as.numeric(nBobs)<0
              #   ,abs(as.numeric(nB2)-as.numeric(nBobs))+1
              #   ,as.numeric(nB2)-as.numeric(nBobs))
              fracB2  <- fracsToUse[2]
              
              segCCFvalues[[segName]][[region]]$fracB1_mean  <- fracB1
              segCCFvalues[[segName]][[region]]$fracB2_mean  <- fracB2
              segCCFvalues[[segName]][[region]]$nBobs_mean   <- nBobs
              #segCCFvalues[[segName]][[region]]$CCF_mean     <- fracB1
              
              # lets also add the raw values (just in case)
              nBobs_all      <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nB<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nB) #copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean
              nBobs_all       <- ifelse(nBobs_all>=maxCPN,maxCPN,nBobs_all)
              nB1             <- nB1#floor(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean))
              nB2             <- nB2#ceiling(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean))
              fracsToUse     <- sapply(nBobs_all,getCCFractionNeeded,nB1,nB2)
              fracsToUse[fracsToUse=="Inf"] <- 1
              
              
              fracB1_all      <-fracsToUse[1,] 
              # fracifelse(as.numeric(nB2)-as.numeric(nBobs_all)<0
              #  ,abs(as.numeric(nB2)-as.numeric(nBobs_all))+1
              #  ,as.numeric(nB2)-as.numeric(nBobs_all))
              fracB2_all       <-fracsToUse[2,] 
              segCCFvalues[[segName]][[region]]$fracB1_all  <- fracB1_all
              segCCFvalues[[segName]][[region]]$fracB2_all  <- fracB2_all
              segCCFvalues[[segName]][[region]]$nBobs_all   <- nBobs_all
              #segCCFvalues[[segName]][[region]]$fracB2_all  <- fracB2
              
              
            }
            
            #next;
            
          }
          if(abs(min(max(ceiling(obsBvals))-obsBvals))>=abs(max(min(floor(obsBvals))-obsBvals)))
          {
            # take the floor 
            # we can determine the low and high values easily 
            nB1  <- min(unique(floor(obsBvals)))
            nB2  <- max(unique(round(obsBvals)))
            
            segCCFvalues[[segName]]$nB1   <- nB1
            segCCFvalues[[segName]]$nB2   <- nB2
            
            # now go region by region, and check the fraction of each (which fraction to count as the CCF depends on whether it is a gain or a loss - determine this below)
            for (region in regions)
            {
              nBobs           <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean) #copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean
              nB1             <- nB1#floor(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean))
              nB2             <- nB2#ceiling(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean))
              fracsToUse      <- getCCFractionNeeded(nBobs,nB1,nB2)
              
              fracB1  <- fracsToUse[1]
              # ifelse(as.numeric(nB2)-as.numeric(nBobs)<0
              #   ,abs(as.numeric(nB2)-as.numeric(nBobs))+1
              #   ,as.numeric(nB2)-as.numeric(nBobs))
              fracB2  <- fracsToUse[2]
              #ifelse(1-fracB1<0,0,1-fracB1)
              
              segCCFvalues[[segName]][[region]]$fracB1_mean  <- fracB1
              segCCFvalues[[segName]][[region]]$fracB2_mean  <- fracB2
              segCCFvalues[[segName]][[region]]$nBobs_mean   <- nBobs
              #segCCFvalues[[segName]][[region]]$CCF_mean     <- fracB1
              
              # lets also add the raw values (just in case)
              nBobs_all      <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nB<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nB) #copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean
              nBobs_all       <- ifelse(nBobs_all>=maxCPN,maxCPN,nBobs_all)
              
              nB1             <- nB1#floor(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean))
              nB2             <- nB2#ceiling(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean))
              fracsToUse     <- sapply(nBobs_all,getCCFractionNeeded,nB1,nB2)
              fracsToUse[fracsToUse=="Inf"] <- 1
              
              
              fracB1_all      <-fracsToUse[1,] 
              # fracifelse(as.numeric(nB2)-as.numeric(nBobs_all)<0
              #  ,abs(as.numeric(nB2)-as.numeric(nBobs_all))+1
              #  ,as.numeric(nB2)-as.numeric(nBobs_all))
              fracB2_all       <-fracsToUse[2,] 
              
              segCCFvalues[[segName]][[region]]$fracB1_all  <- fracB1_all
              segCCFvalues[[segName]][[region]]$fracB2_all  <- fracB2_all
              segCCFvalues[[segName]][[region]]$nBobs_all   <- nBobs_all
              #segCCFvalues[[segName]][[region]]$fracB2_all  <- fracB2
              
              
            }
            
            # next;
            
          }
        }
        
      }
      
      #check whether there are only two states when rounding
      if(length(unique(round(obsBvals)))==2)
      {
        
        # we can determine the low and high values easily
        nB1  <- min(round(obsBvals))
        nB2  <- max(round(obsBvals))
        
        segCCFvalues[[segName]]$nB1   <- nB1
        segCCFvalues[[segName]]$nB2   <- nB2
        
        # now go region by region, and check the fraction of each (which fraction to count as the CCF depends on whether it is a gain or a loss - determine this below)
        for (region in regions)
        {
          nBobs           <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean) #copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean
          nB1             <- nB1#floor(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean))
          nB2             <- nB2#ceiling(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean))
          fracsToUse      <- getCCFractionNeeded(nBobs,nB1,nB2)
          fracB1  <- fracsToUse[1] 
          # ifelse(as.numeric(nB2)-as.numeric(nBobs)<0
          # ,abs(as.numeric(nB2)-as.numeric(nBobs))+1
          # ,as.numeric(nB2)-as.numeric(nBobs))
          fracB2  <- fracsToUse[2]
          # ifelse(1-fracB1<0,0,1-fracB1)
          
          segCCFvalues[[segName]][[region]]$fracB1_mean  <- fracB1
          segCCFvalues[[segName]][[region]]$fracB2_mean  <- fracB2
          segCCFvalues[[segName]][[region]]$nBobs_mean   <- nBobs
          #segCCFvalues[[segName]][[region]]$CCF_mean     <- fracB1
          
          # lets also add the raw values (just in case)
          nBobs_all      <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nB<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nB) #copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean
          nBobs_all       <- ifelse(nBobs_all>=maxCPN,maxCPN,nBobs_all)
          
          nB1             <- nB1#floor(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean))
          nB2             <- nB2#ceiling(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean))
          fracsToUse     <- sapply(nBobs_all,getCCFractionNeeded,nB1,nB2)
          fracsToUse[fracsToUse=="Inf"] <- 1
          
          
          fracB1_all      <-fracsToUse[1,] 
          # fracifelse(as.numeric(nB2)-as.numeric(nBobs_all)<0
          #  ,abs(as.numeric(nB2)-as.numeric(nBobs_all))+1
          #  ,as.numeric(nB2)-as.numeric(nBobs_all))
          fracB2_all       <-fracsToUse[2,] 
          # ifelse(1-fracB1_all<0,0,1-fracB1_all)
          
          segCCFvalues[[segName]][[region]]$fracB1_all  <- fracB1_all
          segCCFvalues[[segName]][[region]]$fracB2_all  <- fracB2_all
          segCCFvalues[[segName]][[region]]$nBobs_all   <- nBobs_all
          #segCCFvalues[[segName]][[region]]$fracB2_all  <- fracB2
          
          
        }
        
        #next;
      }
      
      
      
      
      # check whether an amplification greater than 8
      
      #if(min(expAvals))
      
      # posValues  <- unique(sort(expAvals))
      # lowValue   <- posValues[1]
      # highValue  <- posValues[2]
      # 
      # # choose the max and min that maximise the CCF across the tumour regions (in other words, minimize heterogeneity where possible)
      # posFracA1 <- c()
      # posFracA2 <- c()
      # valuesChosen <- c()
      # for (x in 1:(length(posValues)-1))
      # {
      #   for (y in (x+1):length(posValues))
      #   {
      #     
      #     posLowValue  <- posValues[x]
      #     posHighValue <- posValues[y]
      #     
      #     CCFvalues <- c()
      #     for (i in 1:length(obsBvals))
      #     {
      #       CCFvalues <-c(CCFvalues,(getCCFractionNeeded(obsBvals[i],posLowValue,posHighValue))[2])
      #     }
      #     valuesChosen <- rbind(valuesChosen,c(posLowValue,posHighValue))
      #     posFracA1 <- rbind(posFracA1,CCFvalues)
      #   }
      # }
      
      if(length(unique(round(obsBvals)))>2)
      {
        # we can determine the low and high values easily 
        nB1  <- min(unique(floor(obsBvals+minCCFtoConsider)))
        nB2  <- max(unique(ceiling(obsBvals-minCCFtoConsider)))
        
        segCCFvalues[[segName]]$nB1   <- nB1
        segCCFvalues[[segName]]$nB2   <- nB2
        
        # now go region by region, and check the fraction of each (which fraction to count as the CCF depends on whether it is a gain or a loss - determine this below)
        for (region in regions)
        {
          nBobs           <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean) #copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean
          nB1             <- nB1#floor(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean))
          nB2             <- nB2#ceiling(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean))
          
          fracsToUse      <- getCCFractionNeeded(nBobs,nB1,nB2)
          fracB1  <- fracsToUse[1] 
          # ifelse(as.numeric(nB2)-as.numeric(nBobs)<0
          # ,abs(as.numeric(nB2)-as.numeric(nBobs))+1
          # ,as.numeric(nB2)-as.numeric(nBobs))
          fracB2  <- fracsToUse[2]
          # ifelse(1-fracB1<0,0,1-fracB1)
          
          segCCFvalues[[segName]][[region]]$fracB1_mean  <- fracB1
          segCCFvalues[[segName]][[region]]$fracB2_mean  <- fracB2
          segCCFvalues[[segName]][[region]]$nBobs_mean   <- nBobs
          #segCCFvalues[[segName]][[region]]$CCF_mean     <- fracB1
          
          # lets also add the raw values (just in case)
          nBobs_all      <- ifelse(copyITHRDS$regionComparisonList[[region]][[segName]]$nB<0,0,copyITHRDS$regionComparisonList[[region]][[segName]]$nB) #copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean
          nBobs_all       <- ifelse(nBobs_all>=maxCPN,maxCPN,nBobs_all)
          
          nB1             <- nB1#floor(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean))
          nB2             <- nB2#ceiling(as.numeric(copyITHRDS$regionComparisonList[[region]][[segName]]$nB_mean))
          
          fracsToUse     <- sapply(nBobs_all,getCCFractionNeeded,nB1,nB2)
          fracsToUse[fracsToUse=="Inf"] <- 1
          
          
          fracB1_all      <-fracsToUse[1,] 
          # fracifelse(as.numeric(nB2)-as.numeric(nBobs_all)<0
          #  ,abs(as.numeric(nB2)-as.numeric(nBobs_all))+1
          #  ,as.numeric(nB2)-as.numeric(nBobs_all))
          fracB2_all       <-fracsToUse[2,] 
          # ifelse(1-fracB1_all<0,0,1-fracB1_all)
          
          
          # fracB1_all  <- ifelse(as.numeric(nB2)-as.numeric(nBobs_all)<0
          #                       ,abs(as.numeric(nB2)-as.numeric(nBobs_all))+1
          #                       ,as.numeric(nB2)-as.numeric(nBobs_all))
          # fracB2_all  <- ifelse(1-fracB1_all<0,0,1-fracB1_all)
          
          segCCFvalues[[segName]][[region]]$fracB1_all  <- fracB1_all
          segCCFvalues[[segName]][[region]]$fracB2_all  <- fracB2_all
          segCCFvalues[[segName]][[region]]$nBobs_all   <- nBobs_all
          #segCCFvalues[[segName]][[region]]$fracB2_all  <- fracB2
          
          
        }
      }
      
      
      
      #next;
    }
    
    segCCFvalues$diffFromExpected <- regionPvalsDiffFromExp
    return(segCCFvalues)
    
  }
  
  # also create a table with the output
  #simpleOutPut <- c()
  for (segName in segNames)
  {
    testOut <- getIndividualCCFestimates_withPvals(segName = segName
                                                   ,regions = regions
                                                   ,copyITHRDS = copyITHRDS
                                                   ,regionChrList = regionChrList
                                                   ,maxCPN = maxCPN
                                                   ,minCCFtoConsider = minCCFtoConsider
                                                   ,confInt = confInt)
    segCCFvalues[[segName]] <- testOut
    
    #rowToAdd <- cbind(segName,as.numeric(testOut[[segName]]$nA1),as.numeric(testOut[[segName]]$nA2))
    
    # for (region in regions)
    # {
    #   rowToAdd <- cbind(rowToAdd
    #                     ,as.numeric(testOut[[segName]][[region]]$fracA1_mean)
    #                     ,as.numeric(testOut[[segName]][[region]]$fracA2_mean))
    # }
    # simpleOutPut <- rbind(simpleOutPut,rowToAdd)
    
  }

  return(segCCFvalues)
}

funPlotCCFAcrossGenome  <- function(copyITHRDS,segCCFvalues)
{
  
  Colours <- colorpanel(100,'steelblue4','#deebf7','steelblue4')
  
  # first, get the min.cons.regions
  names_ComparisonOutPut <- names(copyITHRDS[[2]][[1]])
  min.cons.seg <- names(copyITHRDS[[2]][[1]])
  min.cons.seg <- do.call(rbind,strsplit(min.cons.seg,split=":"))
  min.cons.seg <- apply(min.cons.seg,2,as.numeric)
  # next plot each sample individually
  chrom.length.copy   <- fun.chrom.length(min.cons.seg[,1:3])
  names(chrom.length.copy)   <- 1:22
  chrom.segs          <- fun.add.chrom(min.cons.seg[,1:3],chrom.length.copy)
  min.x               <- 0
  min.y               <- -0.5
  max.y               <- 5.5
  max.x               <- as.numeric(max(chrom.segs[,3]))
  
  # how many samples are we dealing with
  Regions     <- grep(":",names(copyITHRDS$regionComparisonList),invert=TRUE,value=TRUE)
  numRegions  <- length(Regions)
  firstRegion <- Regions[1] 
  lastRegion <- Regions[length(Regions)] 
  #regionJoinedSeg <- funJoinAdjacentSegments(copyITHRDS$regionComparisonList)
  regionJoinedSeg <- copyITHRDS$regionComparisonList
  
  par(mfrow=c(numRegions,1))
  
  for (region in Regions)
  {
    if(region ==firstRegion)
    {
      par(mar=c(2,3,4,3))
    }
    
    if(region ==lastRegion)
    {
      par(mar=c(4,3,2,3))
    }
    
    if(!region %in%c(firstRegion,lastRegion))
    {
      par(mar=c(3,3,3,3))
    }
    plot(1
         ,xlim=c(min.x,max.x)
         ,ylim=c(min.y,max.y)
         ,type='n'
         ,xaxt='n'
         ,yaxs='i'
         ,yaxt='n'
         ,xaxs='i'
         ,ylab=""
         ,xlab=""
         ,lwd=2
         ,main=paste(region," A allele",sep="")
         ,bty="n")
    box()
    
    # get the regionSeg
    #copyITHRDS
    
    #regionSeg          <- regionJoinedSeg[[region]]
    #regionSegsToUse_nA <- names(regionSeg$segs_nA)
    #regionSegsToUse_nB <- names(regionSeg$segs_nB)
    regionSeg  <- copyITHRDS$regionComparisonList[[region]]
    regionOneSeg <- copyITHRDS$regionComparisonList[[firstRegion]]
    
    for (x in 1:length(names(segCCFvalues)))
    {
      ChrCoords <- unlist(strsplit(names(segCCFvalues)[x],split=":"))
      segName   <- names(segCCFvalues)[x]
      regionSegsToUse_nA  <- regionSeg[[segName]]
      regionOneSegToConsider <- regionOneSeg[[segName]]
      if(is.na(regionOneSegToConsider$nA))
      {
        next;
      }
      allCCFvals           <- segCCFvalues[[segName]][[segName]][[region]]$fracA1_all
      
      
      toPlot <- cbind(regionOneSegToConsider$chr
                      ,regionOneSegToConsider$pos
                      ,regionOneSegToConsider$pos
                      ,regionSegsToUse_nA$nA
                      ,regionSegsToUse_nA$nA
                      ,regionSegsToUse_nA$nA
                      ,allCCFvals
      )
      # toPlot <- cbind(copyITHRDS$regionComparisonList[[firstRegion]][[regionSegsToUse_nA[x]]]$chr
      #                 ,copyITHRDS$regionComparisonList[[firstRegion]][[regionSegsToUse_nA[x]]]$pos
      #                 ,copyITHRDS$regionComparisonList[[firstRegion]][[regionSegsToUse_nA[x]]]$pos
      #                 ,copyITHRDS$regionComparisonList[[region]][[regionSegsToUse_nA[x]]]$nA
      #                 ,copyITHRDS$regionComparisonList[[region]][[regionSegsToUse_nA[x]]]$nA
      #                 ,copyITHRDS$regionComparisonList[[region]][[regionSegsToUse_nA[x]]]$nA
      #                 ,copyITHRDS$regionComparisonList[[region]][[regionSegsToUse_nA[x]]]$nA)
      
      
      if(is.na(toPlot[1,1]))
      {
        next;
      }
      
      #toPlot <- cbind(logRBAFSample1$chrom,logRBAFSample1$pos,logRBAFSample1$pos,ABSample1,ABSample2)
      toPlot[,1:3] <- fun.add.chrom(apply(toPlot[,1:3],2,as.numeric),chrom.length.copy)
      toPlot <- data.frame(toPlot,stringsAsFactors = FALSE)
      colnames(toPlot) <- c('chr','pos','pos','Sample1nA','Sample1nA','Sample1nA','CCF')
      
      toPlot$Sample1nA <- ifelse(as.numeric(toPlot$Sample1nA)>max.y,(max.y-0.25),as.numeric(toPlot$Sample1nA))
      # toPlot$Sample1nA <- ifelse(as.numeric(toPlot$Sample1nA)<min.y,(min.y+0.25),as.numeric(toPlot$Sample1nA))
      # 
      # toPlot$Sample1nB <- ifelse(as.numeric(toPlot$Sample1nB)>max.y,(max.y-0.25),as.numeric(toPlot$Sample1nB))
      # toPlot$Sample1nB <- ifelse(as.numeric(toPlot$Sample1nB)<min.y,(min.y+0.25),as.numeric(toPlot$Sample1nB))
      # 
      #PValsPaired <- as.numeric(ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$nA_PairedPval)
      #totalPval   <-  ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$totalCPNPairedPval#as.numeric(ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$nA_PairedPval)
      
      colourCode <- ifelse(round(as.numeric(toPlot$CCF)*100)>100,100,round(as.numeric(toPlot$CCF)*100))
      colourCode <- ifelse(colourCode==0,1,colourCode)
      
      
      #MSAIstatus  <- as.numeric(ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$msAI)
      points(x = as.numeric(toPlot$pos),y=as.numeric(toPlot$Sample1nA)
             ,col=Colours[colourCode]
             ,cex=.5
             ,pch=16)#ifelse(MSAIstatus==1,'#238b4550',ifelse(PValsPaired<0.001&totalPval<0.001,'#3182bd50','#fc927250')),cex=0.75,pch=16)
      par(lend=1)
      segColour <- as.numeric(segCCFvalues[[segName]][[segName]][[region]]$fracA1_mean)
      segColour <- Colours[(ifelse(segColour>1,1,segColour)*100)]
      
      segments(x0=min(as.numeric(toPlot$pos))
               ,y0=mean(c(as.numeric(toPlot$Sample1nA)))
               ,x1=max(as.numeric(toPlot$pos)),y1=mean(c(as.numeric(toPlot$Sample1nA)))
               ,col="#525252"
               ,lwd=3)
      
      segments(x0=min(as.numeric(toPlot$pos))
               ,y0=segCCFvalues[[segName]][[segName]]$nA1
               ,x1=max(as.numeric(toPlot$pos))
               ,y1=segCCFvalues[[segName]][[segName]]$nA1
               ,col="#de2d26"
               ,lwd=3)
      
      segments(x0=min(as.numeric(toPlot$pos))
               ,y0=segCCFvalues[[segName]][[segName]]$nA2
               ,x1=max(as.numeric(toPlot$pos))
               ,y1=segCCFvalues[[segName]][[segName]]$nA2
               ,col="#31a354"
               ,lwd=3)
      
      #ifelse(MSAIstatus==1,'#00441b',ifelse(PValsPaired<0.001&totalPval<0.001,'#08306b','#a50f15')),lwd=3)
      
      # points(x = as.numeric(toPlot$pos),y=as.numeric(toPlot$Sample1nB),col='#4292c650'
      #        ,cex=0.75
      #        ,pch=16)#ifelse(MSAIstatus==1,'#238b4550',ifelse(PValsPaired<0.001&totalPval<0.001,'#3182bd50','#fc927250')),cex=0.75,pch=16)
      # par(lend=1)
      # segments(x0=min(as.numeric(toPlot$pos))
      #          ,y0=mean(c(as.numeric(toPlot$Sample1nB)))
      #          ,x1=max(as.numeric(toPlot$pos)),y1=mean(c(as.numeric(toPlot$Sample1nB)))
      #          ,col='#08306b'
      #          ,lwd=3)#ifelse(MSAIstatus==1,'#00441b',ifelse(PValsPaired<0.001&totalPval<0.001,'#08306b','#a50f15')),lwd=3)
      
    }
    
    axis(side =2,las=1)
    abline(h=seq((min.y+0.5),(max.y-0.5),by=1),lty='dashed',lwd=0.75,col='gray')
    
    mtext(side=3
          ,at=fun.chrom.mean(chrom.segs[,1:3])
          ,text=sort(unique(chrom.segs[,1]))
          ,cex=seq(0.6,0.4,length.out=length(unique(chrom.segs[,1])))
          ,line=-1
          #,lwd=0.5
    )
    
    for (i in sort(unique(as.numeric(chrom.segs[,1]))))
    {
      abline(v=max(as.numeric(chrom.segs[chrom.segs[,1]==i,3]))
             ,lwd=0.5)
      #,lty='dashed')
    }
  }
  
  
  par(mfrow=c(numRegions,1))
  
  for (region in Regions)
  {
    if(region ==firstRegion)
    {
      par(mar=c(2,3,4,3))
    }
    
    if(region ==lastRegion)
    {
      par(mar=c(4,3,2,3))
    }
    
    if(!region %in%c(firstRegion,lastRegion))
    {
      par(mar=c(3,3,3,3))
    }
    plot(1
         ,xlim=c(min.x,max.x)
         ,ylim=c(min.y,max.y)
         ,type='n'
         ,xaxt='n'
         ,yaxs='i'
         ,yaxt='n'
         ,xaxs='i'
         ,ylab=""
         ,xlab=""
         ,lwd=2
         ,main=paste(region," B allele",sep="")
         ,bty="n")
    box()
    
    # get the regionSeg
    #copyITHRDS
    
    regionSeg          <- regionJoinedSeg[[region]]
    #regionSegsToUse_nA <- names(regionSeg$segs_nA)
    #regionSegsToUse_nB <- names(regionSeg$segs_nB)
    regionSeg  <- copyITHRDS$regionComparisonList[[region]]
    regionOneSeg <- copyITHRDS$regionComparisonList[[firstRegion]]
    
    for (x in 1:length(names(segCCFvalues)))
    {
      ChrCoords <- unlist(strsplit(names(segCCFvalues)[x],split=":"))
      segName   <- names(segCCFvalues)[x]
      regionSegsToUse_nA  <- regionSeg[[segName]]
      regionOneSegToConsider <- regionOneSeg[[segName]]
      
      if(is.na(regionSegsToUse_nA$nA))
      {
        next;
      }
      
      allCCFvals           <- segCCFvalues[[segName]][[segName]][[region]]$fracB1_all
      
      
      toPlot <- cbind(regionOneSegToConsider$chr
                      ,regionOneSegToConsider$pos
                      ,regionOneSegToConsider$pos
                      ,regionSegsToUse_nA$nB
                      ,regionSegsToUse_nA$nB
                      ,regionSegsToUse_nA$nB
                      ,allCCFvals
      )
      # toPlot <- cbind(copyITHRDS$regionComparisonList[[firstRegion]][[regionSegsToUse_nA[x]]]$chr
      #                 ,copyITHRDS$regionComparisonList[[firstRegion]][[regionSegsToUse_nA[x]]]$pos
      #                 ,copyITHRDS$regionComparisonList[[firstRegion]][[regionSegsToUse_nA[x]]]$pos
      #                 ,copyITHRDS$regionComparisonList[[region]][[regionSegsToUse_nA[x]]]$nA
      #                 ,copyITHRDS$regionComparisonList[[region]][[regionSegsToUse_nA[x]]]$nA
      #                 ,copyITHRDS$regionComparisonList[[region]][[regionSegsToUse_nA[x]]]$nA
      #                 ,copyITHRDS$regionComparisonList[[region]][[regionSegsToUse_nA[x]]]$nA)
      
      
      if(is.na(toPlot[1,1]))
      {
        next;
      }
      
      #toPlot <- cbind(logRBAFSample1$chrom,logRBAFSample1$pos,logRBAFSample1$pos,ABSample1,ABSample2)
      toPlot[,1:3] <- fun.add.chrom(apply(toPlot[,1:3],2,as.numeric),chrom.length.copy)
      toPlot <- data.frame(toPlot,stringsAsFactors = FALSE)
      colnames(toPlot) <- c('chr','pos','pos','Sample1nA','Sample1nA','Sample1nA','CCF')
      
      toPlot$Sample1nA <- ifelse(as.numeric(toPlot$Sample1nA)>max.y,(max.y-0.25),as.numeric(toPlot$Sample1nA))
      # toPlot$Sample1nA <- ifelse(as.numeric(toPlot$Sample1nA)<min.y,(min.y+0.25),as.numeric(toPlot$Sample1nA))
      # 
      # toPlot$Sample1nB <- ifelse(as.numeric(toPlot$Sample1nB)>max.y,(max.y-0.25),as.numeric(toPlot$Sample1nB))
      # toPlot$Sample1nB <- ifelse(as.numeric(toPlot$Sample1nB)<min.y,(min.y+0.25),as.numeric(toPlot$Sample1nB))
      # 
      #PValsPaired <- as.numeric(ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$nA_PairedPval)
      #totalPval   <-  ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$totalCPNPairedPval#as.numeric(ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$nA_PairedPval)
      #Colours <- colorpanel(100,'white','steelblue4')
      
      colourCode <- ifelse(round(as.numeric(toPlot$CCF)*100)>100,100,round(as.numeric(toPlot$CCF)*100))
      colourCode <- ifelse(colourCode==0,1,colourCode)
      
      
      #MSAIstatus  <- as.numeric(ComparisonOutPut[[names_ComparisonOutPut[x]]][[paste(regionOne,regionTwo,sep=":")]]$msAI)
      points(x = as.numeric(toPlot$pos),y=as.numeric(toPlot$Sample1nA)
             ,col=Colours[colourCode]
             ,cex=.5
             ,pch=16)#ifelse(MSAIstatus==1,'#238b4550',ifelse(PValsPaired<0.001&totalPval<0.001,'#3182bd50','#fc927250')),cex=0.75,pch=16)
      par(lend=1)
      segColour <- as.numeric(segCCFvalues[[segName]][[segName]][[region]]$fracB1_mean)
      segColour <- Colours[(ifelse(segColour>1,1,segColour)*100)]
      
      segments(x0=min(as.numeric(toPlot$pos))
               ,y0=mean(c(as.numeric(toPlot$Sample1nA)))
               ,x1=max(as.numeric(toPlot$pos)),y1=mean(c(as.numeric(toPlot$Sample1nA)))
               ,col="#525252"
               ,lwd=3)#ifelse(MSAIstatus==1,'#00441b',ifelse(PValsPaired<0.001&totalPval<0.001,'#08306b','#a50f15')),lwd=3)
      
      segments(x0=min(as.numeric(toPlot$pos))
               ,y0=segCCFvalues[[segName]][[segName]]$nB1
               ,x1=max(as.numeric(toPlot$pos))
               ,y1=segCCFvalues[[segName]][[segName]]$nB1
               ,col="#de2d26"
               ,lwd=3)
      
      segments(x0=min(as.numeric(toPlot$pos))
               ,y0=segCCFvalues[[segName]][[segName]]$nB2
               ,x1=max(as.numeric(toPlot$pos))
               ,y1=segCCFvalues[[segName]][[segName]]$nB2
               ,col="#31a354"
               ,lwd=3)
      
      # points(x = as.numeric(toPlot$pos),y=as.numeric(toPlot$Sample1nB),col='#4292c650'
      #        ,cex=0.75
      #        ,pch=16)#ifelse(MSAIstatus==1,'#238b4550',ifelse(PValsPaired<0.001&totalPval<0.001,'#3182bd50','#fc927250')),cex=0.75,pch=16)
      # par(lend=1)
      # segments(x0=min(as.numeric(toPlot$pos))
      #          ,y0=mean(c(as.numeric(toPlot$Sample1nB)))
      #          ,x1=max(as.numeric(toPlot$pos)),y1=mean(c(as.numeric(toPlot$Sample1nB)))
      #          ,col='#08306b'
      #          ,lwd=3)#ifelse(MSAIstatus==1,'#00441b',ifelse(PValsPaired<0.001&totalPval<0.001,'#08306b','#a50f15')),lwd=3)
      
    }
    
    axis(side =2,las=1)
    abline(h=seq((min.y+0.5),(max.y-0.5),by=1),lty='dashed',lwd=0.75,col='gray')
    
    mtext(side=3
          ,at=fun.chrom.mean(chrom.segs[,1:3])
          ,text=sort(unique(chrom.segs[,1]))
          ,cex=seq(0.6,0.4,length.out=length(unique(chrom.segs[,1])))
          ,line=-1
          #,lwd=0.5
    )
    
    for (i in sort(unique(as.numeric(chrom.segs[,1]))))
    {
      abline(v=max(as.numeric(chrom.segs[chrom.segs[,1]==i,3]))
             ,lwd=0.5)
      #,lty='dashed')
    }
  }
  
}

getProportionTrunkBranchAberrant <- function(regionComparisonList
                                             ,PvalToSeparate=0.005
                                             ,pValThresh=0.05
                                             ,gainPValToUse = 'gainPval'
                                             ,lossPValToUse = 'lossPval'
                                             ,TRACERxID)
{
  print("Called getProportionTrunkBranchAberrant")
  regionsToConsider <- grep(":",names(regionComparisonList),invert=TRUE,value=TRUE)
  segsToConsider    <- names(regionComparisonList[[1]])
  
  min.gain.trunk           <- data.frame(matrix(NA,ncol=length(regionsToConsider),nrow=length(segsToConsider)))
  rownames(min.gain.trunk) <- segsToConsider
  colnames(min.gain.trunk) <- regionsToConsider
  min.gain <- min.loss <- min.loss.trunk  <- min.loss.branch <- min.gain.branch <- min.gain.trunk
  
  # go region by region
  for (region in regionsToConsider)
  {
    regionOut <- regionComparisonList[[region]]
    for (seg in segsToConsider)
    {
      regionSegStatus <- regionOut[[seg]]$segStatus
      if(length(regionSegStatus)==0)
      {
        next;
      }
      min.gain[seg,region] <- ifelse(as.numeric(regionSegStatus[gainPValToUse])<pValThresh,1,0)
      min.loss[seg,region] <- ifelse(as.numeric(regionSegStatus[lossPValToUse])<pValThresh,1,0)
      
    }
    
  }
  
  # remove the NAs
  min.gain <- min.gain[!is.na(rowSums(min.gain)),,drop=FALSE]
  min.loss <- min.loss[!is.na(rowSums(min.loss)),,drop=FALSE]
  
  trunkGain  <- names(which(rowSums(min.gain)==length(regionsToConsider)))
  trunkLoss  <- names(which(rowSums(min.loss)==length(regionsToConsider)))
  branchGain <- names(which(rowSums(min.gain)<length(regionsToConsider)&rowSums(min.gain)>0))
  branchLoss <- names(which(rowSums(min.loss)<length(regionsToConsider)&rowSums(min.loss)>0))
  
  
  genomeSize    <- do.call(rbind,strsplit(rownames(min.gain),split=":"))
  genomeSize    <- sum(as.numeric(genomeSize[,3])-as.numeric(genomeSize[,2]))
  
  #get a simple summary gain
  GenomeGain <- do.call(rbind,strsplit(names(which(rowSums(min.gain)>0)),split=":"))
  GenomeGain <- (sum(as.numeric(GenomeGain[,3]) - as.numeric(GenomeGain[,2])))
  propGenomeGain  <- GenomeGain/genomeSize
  
  #get a simple summary loss
  GenomeLoss <- do.call(rbind,strsplit(names(which(rowSums(min.loss)>0)),split=":"))
  GenomeLoss <- (sum(as.numeric(GenomeLoss[,3]) - as.numeric(GenomeLoss[,2])))
  propGenomeLoss <- GenomeLoss/genomeSize
  
  # [Note the gain and loss could add up to more than 100% as a segment can be both subject to gain and loss in different samples]
  GenomeAberrant     <- GenomeLoss+GenomeGain
  propGenomeAberrant <- propGenomeGain+propGenomeLoss
  
  GenomeTrunkGain      <- do.call(rbind,strsplit(trunkGain,split=":"))
  GenomeTrunkGain      <- sum(as.numeric(GenomeTrunkGain[,3])-as.numeric(GenomeTrunkGain[,2]))
  propTrunkGain        <- GenomeTrunkGain/GenomeGain
  
  GenomeTrunkLoss      <- do.call(rbind,strsplit(trunkLoss,split=":"))
  GenomeTrunkLoss      <- sum(as.numeric(GenomeTrunkLoss[,3])-as.numeric(GenomeTrunkLoss[,2]))
  propTrunkLoss        <- GenomeTrunkLoss/GenomeLoss
  
  propTrunkAberrant    <- (GenomeTrunkLoss+GenomeTrunkGain)/GenomeAberrant
  
  
  GenomeBranchGain      <- do.call(rbind,strsplit(branchGain,split=":"))
  GenomeBranchGain      <- sum(as.numeric(GenomeBranchGain[,3])-as.numeric(GenomeBranchGain[,2]))
  propBranchGain        <- GenomeBranchGain/GenomeGain
  
  GenomeBranchLoss      <- do.call(rbind,strsplit(branchLoss,split=":"))
  GenomeBranchLoss      <- sum(as.numeric(GenomeBranchLoss[,3])-as.numeric(GenomeBranchLoss[,2]))
  propBranchLoss        <- GenomeBranchLoss/GenomeLoss
  
  propBranchAberrant    <- (GenomeBranchLoss+GenomeBranchGain)/GenomeAberrant
  
  simpleSummary <- cbind(GenomeAberrant,GenomeGain,GenomeLoss,GenomeTrunkGain,GenomeTrunkLoss,GenomeBranchGain,GenomeBranchLoss,propTrunkGain,propTrunkLoss,propTrunkAberrant,propBranchGain,propBranchLoss,propBranchAberrant)
  rownames(simpleSummary) <- TRACERxID
  #outPut$sampleSummary <- simpleSummary
  
  minCopyLists  <- list()
  minCopyLists$min.gain        <- min.gain
  minCopyLists$min.loss        <- min.loss
  
  
  
  #perform a test to see whether there is a significant difference
  #between the genome of the different regions
  revertToTrunkGain <- c()
  for (seg in branchGain)
  {
    
    GainRegions   <-     names(which(unlist(c(min.gain[seg,])>0)))
    NoGainRegions <-     names(which(unlist(c(min.gain[seg,])==0)))
    
    # make the assumption that the sample with the biggest p-value will be more likely to be similar to the others
    referenceGain <- GainRegions[1]
    refRegionOut <- regionComparisonList[[referenceGain]]
    refRegionPval <- as.numeric(refRegionOut[[seg]]$segStatus[gainPValToUse])
    
    
    for(gainRegion in GainRegions)
    {
      
      testRegionOut <- regionComparisonList[[gainRegion]]
      if( as.numeric(testRegionOut[[seg]]$segStatus[gainPValToUse])>refRegionPval)
      {
        referenceGain <- gainRegion
        refRegionOut <- regionComparisonList[[referenceGain]]
        refRegionPval <- as.numeric(refRegionOut[[seg]]$segStatus[gainPValToUse])
        
      }
      
    }
    
    #next, check whether each NoGainRegion is significantly different from the reference region
    
    reconsider <- c()
    for (noGainRegion in NoGainRegions)
    {
      combinedName <- c(referenceGain,noGainRegion)
      combinedName <- paste(sort(combinedName)[1],sort(combinedName)[2],sep=":")
      
      if(length(regionComparisonList[[combinedName]][[seg]]$combinedPairedPval)!=0){
        if(regionComparisonList[[combinedName]][[seg]]$combinedPairedPval>PvalToSeparate)
        {
          reconsider <- c(reconsider,noGainRegion)
        }
      }
      
      
    }
    
    if(identical(reconsider,NoGainRegions))
    {
      revertToTrunkGain <- c(revertToTrunkGain,seg)
    }
    
    
    
    
  }
  
  revertToTrunkLoss <- c()
  for (seg in branchLoss)
  {
    
    LossRegions   <-     names(which(unlist(c(min.loss[seg,])>0)))
    NoLossRegions <-     names(which(unlist(c(min.loss[seg,])==0)))
    
    # make the assumption that the sample with the smallest p-value will be more likely to be similar to the others
    referenceLoss <- LossRegions[1]
    refRegionOut <- regionComparisonList[[referenceLoss]]
    refRegionPval <- as.numeric(refRegionOut[[seg]]$segStatus[lossPValToUse])
    
    
    for(lossRegion in LossRegions)
    {
      
      testRegionOut <- regionComparisonList[[lossRegion]]
      if( as.numeric(testRegionOut[[seg]]$segStatus[lossPValToUse])>refRegionPval)
      {
        referenceLoss <- lossRegion
        refRegionOut <- regionComparisonList[[referenceLoss]]
        refRegionPval <- as.numeric(refRegionOut[[seg]]$segStatus[lossPValToUse])
        
      }
      
    }
    
    #next, check whether each NoLossRegion is significantly different from the reference region
    
    reconsider <- c()
    for (noLossRegion in NoLossRegions)
    {
      combinedName <- c(referenceLoss,noLossRegion)
      combinedName <- paste(sort(combinedName)[1],sort(combinedName)[2],sep=":")
      
      if(length(regionComparisonList[[combinedName]][[seg]]$combinedPairedPval)!=0){
        if(regionComparisonList[[combinedName]][[seg]]$combinedPairedPval>PvalToSeparate)
        {
          reconsider <- c(reconsider,noLossRegion)
        }
      }
      
      
    }
    
    if(identical(reconsider,NoLossRegions))
    {
      revertToTrunkLoss <- c(revertToTrunkLoss,seg)
    }
    
    
    
    
  }
  
  min.gain[revertToTrunkGain,] <- 1
  min.loss[revertToTrunkLoss,] <- 1
  
  trunkGain  <- names(which(rowSums(min.gain)==length(regionsToConsider)))
  trunkLoss  <- names(which(rowSums(min.loss)==length(regionsToConsider)))
  branchGain <- names(which(rowSums(min.gain)<length(regionsToConsider)&rowSums(min.gain)>0))
  branchLoss <- names(which(rowSums(min.loss)<length(regionsToConsider)&rowSums(min.loss)>0))
  
  
  genomeSize    <- do.call(rbind,strsplit(rownames(min.gain),split=":"))
  genomeSize    <- sum(as.numeric(genomeSize[,3])-as.numeric(genomeSize[,2]))
  
  #get a simple summary gain
  GenomeGain <- do.call(rbind,strsplit(names(which(rowSums(min.gain)>0)),split=":"))
  GenomeGain <- (sum(as.numeric(GenomeGain[,3]) - as.numeric(GenomeGain[,2])))
  propGenomeGain  <- GenomeGain/genomeSize
  
  #get a simple summary loss
  GenomeLoss <- do.call(rbind,strsplit(names(which(rowSums(min.loss)>0)),split=":"))
  GenomeLoss <- (sum(as.numeric(GenomeLoss[,3]) - as.numeric(GenomeLoss[,2])))
  propGenomeLoss <- GenomeLoss/genomeSize
  
  # [Note the gain and loss could add up to more than 100% as a segment can be both subject to gain and loss in different samples]
  GenomeAberrant     <- GenomeLoss+GenomeGain
  propGenomeAberrant <- propGenomeGain+propGenomeLoss
  
  GenomeTrunkGain      <- do.call(rbind,strsplit(trunkGain,split=":"))
  GenomeTrunkGain      <- sum(as.numeric(GenomeTrunkGain[,3])-as.numeric(GenomeTrunkGain[,2]))
  propTrunkGain        <- GenomeTrunkGain/GenomeGain
  
  GenomeTrunkLoss      <- do.call(rbind,strsplit(trunkLoss,split=":"))
  GenomeTrunkLoss      <- sum(as.numeric(GenomeTrunkLoss[,3])-as.numeric(GenomeTrunkLoss[,2]))
  propTrunkLoss        <- GenomeTrunkLoss/GenomeLoss
  
  propTrunkAberrant    <- (GenomeTrunkLoss+GenomeTrunkGain)/GenomeAberrant
  
  
  GenomeBranchGain      <- do.call(rbind,strsplit(branchGain,split=":"))
  GenomeBranchGain      <- sum(as.numeric(GenomeBranchGain[,3])-as.numeric(GenomeBranchGain[,2]))
  propBranchGain        <- GenomeBranchGain/GenomeGain
  
  GenomeBranchLoss      <- do.call(rbind,strsplit(branchLoss,split=":"))
  GenomeBranchLoss      <- sum(as.numeric(GenomeBranchLoss[,3])-as.numeric(GenomeBranchLoss[,2]))
  propBranchLoss        <- GenomeBranchLoss/GenomeLoss
  
  propBranchAberrant    <- (GenomeBranchLoss+GenomeBranchGain)/GenomeAberrant
  
  simpleSummaryNew <- cbind(GenomeAberrant,GenomeGain,GenomeLoss,GenomeTrunkGain,GenomeTrunkLoss,GenomeBranchGain,GenomeBranchLoss,propTrunkGain,propTrunkLoss,propTrunkAberrant,propBranchGain,propBranchLoss,propBranchAberrant)
  rownames(simpleSummaryNew) <- TRACERxID
  
  # next, see if there are any gains or losses that should in fact be heterogeneous
  
  revertToBranchGain <- c()
  for (seg in trunkGain)
  {
    
    GainRegions   <-     names(which(unlist(c(min.gain[seg,])>0)))
    NoGainRegions <-     names(which(unlist(c(min.gain[seg,])==0)))
    
    # make the assumption that the sample with the biggest p-value will be more likely to be similar to the others
    biggestGain     <- GainRegions[1]
    biggestGainOut  <- regionComparisonList[[biggestGain]]
    biggestGainPval <- as.numeric(biggestGainOut[[seg]]$segStatus[gainPValToUse])
    smallestGain     <- GainRegions[2]
    smallestGainOut  <- regionComparisonList[[smallestGain]]
    smallestGainPval <- as.numeric(smallestGainOut[[seg]]$segStatus[gainPValToUse])
    
    for(gainRegion in GainRegions)
    {
      
      testRegionOut <- regionComparisonList[[gainRegion]]
      if( as.numeric(testRegionOut[[seg]]$segStatus[gainPValToUse])>smallestGainPval)
      {
        smallestGain <- gainRegion
        smallestGainOut <- regionComparisonList[[smallestGain]]
        smallestGainPval <- as.numeric(smallestGainOut[[seg]]$segStatus[gainPValToUse])
        
      }
      
      if( as.numeric(testRegionOut[[seg]]$segStatus[gainPValToUse])<biggestGainPval)
      {
        biggestGain <- gainRegion
        biggestGainOut <- regionComparisonList[[biggestGain]]
        biggestGainPval <- as.numeric(biggestGainOut[[seg]]$segStatus[gainPValToUse])
        
      }
      
    }
    
    #next, check whether smallestGain is significantly different from biggestGain
    
    #reconsider <- c()
    combinedName <- c(biggestGain,smallestGain)
    combinedName <- paste(sort(combinedName)[1],sort(combinedName)[2],sep=":")
    if(!is.na(as.numeric(regionComparisonList[[combinedName]][[seg]]$nA_PairedPval)))
    {
      if(as.numeric(regionComparisonList[[combinedName]][[seg]]$nA_PairedPval)<PvalToSeparate)
      {
        revertToBranchGain <- c(revertToBranchGain,seg)
      }
    }

    
    
  }
  
  
  revertToBranchLoss <- c()
  for (seg in trunkLoss)
  {
    
    LossRegions   <-     names(which(unlist(c(min.loss[seg,])>0)))
    NoLossRegions <-     names(which(unlist(c(min.loss[seg,])==0)))
    
    # make the assumption that the sample with the biggest p-value will be more likely to be similar to the others
    biggestLoss     <- LossRegions[1]
    biggestLossOut  <- regionComparisonList[[biggestLoss]]
    biggestLossPval <- as.numeric(biggestLossOut[[seg]]$segStatus[lossPValToUse])
    smallestLoss     <- LossRegions[2]
    smallestLossOut  <- regionComparisonList[[smallestLoss]]
    smallestLossPval <- as.numeric(smallestLossOut[[seg]]$segStatus[lossPValToUse])
    
    for(lossRegion in LossRegions)
    {
      
      testRegionOut <- regionComparisonList[[lossRegion]]
      if( as.numeric(testRegionOut[[seg]]$segStatus[lossPValToUse])>smallestLossPval)
      {
        smallestLoss <- lossRegion
        smallestLossOut <- regionComparisonList[[smallestLoss]]
        smallestLossPval <- as.numeric(smallestLossOut[[seg]]$segStatus[lossPValToUse])
        
      }
      
      if( as.numeric(testRegionOut[[seg]]$segStatus[lossPValToUse])<biggestLossPval)
      {
        biggestLoss <- lossRegion
        biggestLossOut <- regionComparisonList[[biggestLoss]]
        biggestLossPval <- as.numeric(biggestLossOut[[seg]]$segStatus[lossPValToUse])
        
      }
      
    }
    
    #next, check whether smallestLoss is significantly different from biggestLoss
    
    #reconsider <- c()
    combinedName <- c(biggestLoss,smallestLoss)
    combinedName <- paste(sort(combinedName)[1],sort(combinedName)[2],sep=":")
    
   
    if(!is.na(regionComparisonList[[combinedName]][[seg]]$nB_PairedPval)&length(regionComparisonList[[combinedName]][[seg]]$nB_PairedPval)!=0){
      if(as.numeric(regionComparisonList[[combinedName]][[seg]]$nB_PairedPval)<PvalToSeparate)
      {
        revertToBranchLoss <- c(revertToBranchLoss,seg)
      }
    }
    
    
    
  }
  
  
  min.gain[revertToBranchGain,] <- 0.5
  min.loss[revertToBranchLoss,] <- 0.5
  
  trunkGain  <- names(which(rowSums(min.gain)==length(regionsToConsider)))
  trunkLoss  <- names(which(rowSums(min.loss)==length(regionsToConsider)))
  branchGain <- names(which(rowSums(min.gain)<length(regionsToConsider)&rowSums(min.gain)>0))
  branchLoss <- names(which(rowSums(min.loss)<length(regionsToConsider)&rowSums(min.loss)>0))
  
  
  genomeSize    <- do.call(rbind,strsplit(rownames(min.gain),split=":"))
  genomeSize    <- sum(as.numeric(genomeSize[,3])-as.numeric(genomeSize[,2]))
  
  #get a simple summary gain
  GenomeGain <- do.call(rbind,strsplit(names(which(rowSums(min.gain)>0)),split=":"))
  GenomeGain <- (sum(as.numeric(GenomeGain[,3]) - as.numeric(GenomeGain[,2])))
  propGenomeGain  <- GenomeGain/genomeSize
  
  #get a simple summary loss
  GenomeLoss <- do.call(rbind,strsplit(names(which(rowSums(min.loss)>0)),split=":"))
  GenomeLoss <- (sum(as.numeric(GenomeLoss[,3]) - as.numeric(GenomeLoss[,2])))
  propGenomeLoss <- GenomeLoss/genomeSize
  
  # [Note the gain and loss could add up to more than 100% as a segment can be both subject to gain and loss in different samples]
  GenomeAberrant     <- GenomeLoss+GenomeGain
  propGenomeAberrant <- propGenomeGain+propGenomeLoss
  
  GenomeTrunkGain      <- do.call(rbind,strsplit(trunkGain,split=":"))
  GenomeTrunkGain      <- sum(as.numeric(GenomeTrunkGain[,3])-as.numeric(GenomeTrunkGain[,2]))
  propTrunkGain        <- GenomeTrunkGain/GenomeGain
  
  GenomeTrunkLoss      <- do.call(rbind,strsplit(trunkLoss,split=":"))
  GenomeTrunkLoss      <- sum(as.numeric(GenomeTrunkLoss[,3])-as.numeric(GenomeTrunkLoss[,2]))
  propTrunkLoss        <- GenomeTrunkLoss/GenomeLoss
  
  propTrunkAberrant    <- (GenomeTrunkLoss+GenomeTrunkGain)/GenomeAberrant
  
  
  GenomeBranchGain      <- do.call(rbind,strsplit(branchGain,split=":"))
  GenomeBranchGain      <- sum(as.numeric(GenomeBranchGain[,3])-as.numeric(GenomeBranchGain[,2]))
  propBranchGain        <- GenomeBranchGain/GenomeGain
  
  GenomeBranchLoss      <- do.call(rbind,strsplit(branchLoss,split=":"))
  GenomeBranchLoss      <- sum(as.numeric(GenomeBranchLoss[,3])-as.numeric(GenomeBranchLoss[,2]))
  propBranchLoss        <- GenomeBranchLoss/GenomeLoss
  
  propBranchAberrant    <- (GenomeBranchLoss+GenomeBranchGain)/GenomeAberrant
  
  simpleSummaryLatest <- cbind(GenomeAberrant,GenomeGain,GenomeLoss,GenomeTrunkGain,GenomeTrunkLoss,GenomeBranchGain,GenomeBranchLoss,propTrunkGain,propTrunkLoss,propTrunkAberrant,propBranchGain,propBranchLoss,propBranchAberrant)
  rownames(simpleSummaryLatest) <- TRACERxID
  
  
  returnOut <- list()
  returnOut$simpleSummary    <- simpleSummary  
  returnOut$simpleSummaryNew <- simpleSummaryNew
  returnOut$simpleSummaryLatest <- simpleSummaryLatest
  returnOut$minCopy          <- minCopyLists
  print("Returning from getProportionTrunkBranchAberrant")
  return(returnOut)
}

#
tTestWithErrorHandling <- function(x, y = NULL, alternative = c("two.sided", "less", "greater"),
                                   mu = 0, paired = FALSE, var.equal = FALSE, conf.level = 0.95,errors.not.permitted=""){
  
  test.results <- tryCatch({
    return(t.test(x=x, y=y, alternative=alternative, mu=mu, paired=paired, var.equal=var.equal, conf.level=conf.level))},
    error = function(err) {
      message(err)
      write(err$message, stderr())
      write("\n", stderr())
      if(grepl("NULL", err$message) & "NULL"%in%errors.not.permitted){
        stop("t.test ERROR with a NULL")
      } else if (grepl("enough", err$message) & "enough"%in%errors.not.permitted){
        stop("t.test ERROR with not enough")
      } else if (grepl("constant", err$message) & "constant"%in%errors.not.permitted){
        stop("t.test ERROR with constant data")
      } else if(grepl("enough|constant|NULL", err$message) & "other"%in%errors.not.permitted){
        stop(err$message)
      }
      cat("\n")
      paste("Input for t test")
      print("x")
      print(x)
      print("y")
      print(y)
      print("mu")
      print(mu)
      print("paired")
      print(paired)
      print("conf.level")
      print(conf.level)
      test.results <- list()
      test.results[["stat"]] <- 1
      test.results[["conf.int"]] <- c(NA, NA)
      test.results[["p.value"]] <- 1
      return(test.results)
    })
  return(test.results)
}
