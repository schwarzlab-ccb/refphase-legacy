source('/camp/lab/swantonc/working/rosentr/scripts/multi/boop.dir.R')

newQCfile       <- "/camp/lab/swantonc/working/mcgrann/projects/allQCs2019-01-21.txt"
#newQCfile       <- ""
fullQCsamples   <- read.table("/camp/lab/swantonc/working/mcgrann/projects/allQCs2018-10-05.txt",sep="\t",stringsAsFactors = FALSE,header=TRUE)
newQCsamples   <- read.table("/camp/lab/swantonc/working/mcgrann/projects/allQCs2018-12-05.txt",sep="\t",stringsAsFactors = FALSE,header=TRUE)

CPNfailSamples    <- c('B_LTX605','M_LTX483')

test <- newQCsamples[!newQCsamples$Sample%in%fullQCsamples$Sample,,drop=FALSE]


tracerx.samples  <- unique(newQCsamples$Patient)
tracerx.script   <- "/camp/lab/swantonc/working/mcgrann/tracerx-exome-pipeline/tracerx.mphase.R"
scratchDir            <- '/camp/project/proj-tracerx-lung/tctProjects/lungTx/scratchFolder/mPHASE_release_20181204/'
releaseDir            <- "/camp/project/proj-tracerx-lung/tracerx/release_tx500_180615/"
releaseCore           <- "/camp/project/proj-tracerx-lung/tracerx/"
RUN.script            <- TRUE
nameToUse             <- "mphase"
samplesToRun          <- c()

# submit <- c('B_LTX288','B_LTX484')
# test <- unique(complete.mPHASEfirst)[!unique(complete.mPHASEfirst)%in%submit]
# submit <- c('U_LTX752', 'U_LTX851', 'U_LTX903', 'U_LTX977')
# submit <- c('M_LTX532')
# submit <- unique(complete.mPHASEfirst)
# submit <- c('B_LTX503')
bigTable <- bigTable[!duplicated(bigTable[,1]),,drop=FALSE]

submitThis <- c(unique(bigTable[,1]))
#submitThis <- c('U_LTX851')
mPHASEmissing <- missingOnes[missingOnes$mPHASE_complete%in%FALSE,,]
submitThis <- submitThis[submitThis%in%mPHASEmissing$TRACERxID]
#submitThis <- c('U_LTX231')
# submitThis <- c('U_LTX979','F_LTX523','F_LTX523','U_LTX852')
# submitThis <- c('U_LTX852')

for(sample in submitThis){
  
  print(sample)
  
  # has this patient been run before
  # tmp1 <- fullQCsamples[fullQCsamples$Patient%in%sample,,drop=FALSE]
  # tmp2 <- newQCsamples[newQCsamples$Patient%in%sample,,drop=FALSE]
  # if(nrow(tmp1)>=1)
  # {
  #   rownames(tmp1) <- 1:nrow(tmp1)
  # }
  # 
  # rownames(tmp2) <- 1:nrow(tmp2)
  # 
  # if(identical(tmp1,tmp2))
  # {
  #   #next;
  # }
  # samplesToRun <- c(samplesToRun,sample)
  specRelease <- paste(releaseCore,bigTable[bigTable[,1]%in%sample,2],"/",sep="")
  
  full.patient         <- sample
  patient              <- unlist(strsplit(full.patient,split="_"))[2]
  designFile           <- paste(specRelease,full.patient,  "/exome/commandHistory/", full.patient, ".exome.design.txt",sep="")
  paramsFile           <- paste(specRelease,full.patient,  "/exome/commandHistory/", full.patient, ".exome.params.txt",sep="")
  save.to.scratch      <- TRUE
  full.QCfile          <- "/camp/project/proj-tracerx-lung/tracerx/release_tx500_180615/sampleStatus/tx450QCs_2018-12-13_gw.txt"
  allow.dip.test.phasing.flag <- FALSE
  gamma.val           <- 1
  num.cores           <- 1
  overwrite.flag      <- TRUE
  data.format.to.use  <- "TRACERx_exome"
  
  # cmdDir
  cmdDir <- paste(scratchDir,full.patient,'/',sep='')
  boop.dir(cmdDir)
  setwd(cmdDir)
  
  # sbatch file
  time.and.date <- gsub(":", "-", gsub(" +", "_", date()))
  
  if(RUN.script){
    
    save.name     <- paste('msub.', time.and.date, '.', full.patient, ".", nameToUse,".sh", sep = "")
    system(paste('touch ', save.name, sep = ''))
    
    msubOut   <- paste(cmdDir,full.patient, '.', time.and.date, ".",nameToUse,".stdout",sep='')
    msubErr   <- paste(cmdDir,full.patient, '.', time.and.date, ".",nameToUse,".stderr",sep='')
    
    # shebang
    bashInterpreter   <- "#!/bin/bash"
    txt <- paste(bashInterpreter, '\n', '\n', sep = '')
    write.table(txt, save.name, quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE) 
    
    # ml required
    txt <- 'module load BEDTools/2.26.0-foss-2016b'
    write.table(txt, save.name, quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
    
    txt <- 'module load SAMtools/1.3.1-foss-2016b'
    write.table(txt, save.name, quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
    
    txt <- 'module load R/3.3.1-foss-2016b-bioc-3.3-libX11-1.6.3'
    write.table(txt, save.name, quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
    
    # txt <- 'module load novoalign/3.07.00'
    # write.table(txt, save.name, quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
    
    txt <- 'module load TracerX-Picard-GATK/0.1-Java-1.7.0_80'
    write.table(txt, save.name, quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
    
    # txt <- 'module load Singularity/2.3.1-GCC-5.4.0-2.26'
    # write.table(txt, save.name, quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
    
    # txt <- 'module load Sambamba/0.6.6'
    # write.table(txt, save.name, quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
    
    # txt <- 'module load BWA/0.7.15-foss-2016b'
    # write.table(txt, save.name, quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
    
    
    #txt <- paste("Rscript /camp/lab/swantonc/working/rosentr/scripts/lohhla/HLAlossFullScript_updated__user_opt_tmp_RR.R ", 
    #""
    txt <- paste("Rscript ", tracerx.script, 
                 " ", patient,
                 " ", designFile,
                 " ", paramsFile,
                 " ", save.to.scratch,
                 " ", newQCfile,
                 " ", allow.dip.test.phasing.flag,
                 " ", gamma.val,
                 " ", num.cores,
                 " ", overwrite.flag,
                 " ", data.format.to.use
                 ,sep="")
   
           
    write.table(txt, save.name, quote = F, row.names = F, col.names = F, append = TRUE)
    
    msubOptions <- ' --partition=TracerX --mem=48GB'
    msubOptions <- ' '
    
    # and submit
    msubStd <- paste("-o ",msubOut," -e ",msubErr,sep="")
    
    cmd <- paste('sbatch -J ', full.patient, msubOptions, ' ', msubStd, ' ', save.name, sep = '')
    write.table(cmd, file = paste('submit.', time.and.date, '.sh', sep = ''), sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    #check whether submitting from camp (or on home computer)
    checkComputer <- "hostname"
    whichComputer <- system(checkComputer,intern=TRUE)
    if(whichComputer=="C02RN1TVG8WL")
    {
      cat('\nsubmitting camp job\n')
      newCMD <- paste("/usr/bin/ssh camp ", paste('sbatch -J ', full.patient, msubOptions, ' ', msubStd, ' ', cmdDir, save.name, sep = ''),sep="")
      cat('\n')
      cat(newCMD)
      system(newCMD)
   
    }
    if(whichComputer!="C02RN1TVG8WL")
    {
      system(cmd)
      
    }
    
    
  }
  
}



