#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
print("args")
print(args)
patient <- args[1]
allow.dip.test.phasing.flag <- as.logical(args[2])
output.directory.base.path <- args[3]
data.format.to.use <- args[4]#Default for this is MEDICC, but can specify TRACERx_exome, TRACERx_snp_array or TRACERx_WGS
snp.file.path <- args[5]#Points towards the snp files, of whatever format they are given that should be used based on data.format to use.
seg.file.path <- args[6]#Pintts towards the seg files, of whatever format they should be determined byh the data.format.to.use flag.
overwrite.flag  <- as.logical(args[7])
gamma.val <- args[8]
num.cores <- args[9]
intersection.or.whole.arm <- args[10]#Can either be intersection or whole_arm


#function taken from https://github.com/molgenis/molgenis-pipelines/wiki/How-to-source-another_file.R-from-within-your-R-script
LocationOfThisScript  <-  function() { # Function LocationOfThisScript returns the location of this .R script (may be needed to source other files in same dir)
	this.file = NULL
	# This file may be 'sourced'
	for (i in - (1:sys.nframe())) {
		if (identical(sys.function(i), base::source)) this.file = (normalizePath(sys.frame(i)$ofile))
	}#end for (i in - (1:sys.nframe())) {

	if (!is.null(this.file)) {

	 return(dirname(this.file))

	}#end if (!is.null(this.file)) {

	# But it may also be called from the command line
	cmd.args <-  commandArgs(trailingOnly = FALSE)
	cmd.args.trailing <- commandArgs(trailingOnly = TRUE)
	cmd.args <- cmd.args[seq.int(from = 1, length.out = length(cmd.args) - length(cmd.args.trailing))]
	res <- gsub("^(?:--file=(.*)|.*)$", "\\1", cmd.args)

	# If multiple --file arguments are given, R uses the last one
	res <- tail(res[res != ""], 1)
	if (0 < length(res)) {

		return(dirname(res))

	}#end if (0 < length(res)) {

	# Both are not the case. Maybe we are in an R GUI?
	return(NULL)
}#end LocationOfThisScript  <-  function() 

current.dir <-  LocationOfThisScript()

# #Debugging option parameters:
# patient <- "MET048"
# allow.dip.test.phasing.flag  <- FALSE
# output.directory.base.path <- "/home/camp/watkint/working/watkint/data/brastianos/release_v2/B_MET048/exome/CN_het"
# data.format.to.use <- "MEDICC"
# snp.file.path <- "/home/camp/watkint/working/watkint/data/brastianos/release_v2/B_MET048/exome/MSAI/MET048_joint_snps.txt"
# seg.file.path <- "/home/camp/watkint/working/watkint/data/brastianos/release_v2/B_MET048/exome/MSAI/MET048_joint_segs_all_ascat_info.txt"
# ##########
# mphase.script.path <- "/home/camp/watkint/working/watkint/scripts/multiregioncopynum/mphase.R"
# #MEDICC format for TRACERx example
# medicc.mphase.command <- paste("Rscript", mphase.script.path, patient, allow.dip.test.phasing.flag, data.format.to.use, snp.file.path, seg.file.path, sep=" ")

# #Debugging option parameters:
# patient <- "MET073"
# allow.dip.test.phasing.flag  <- FALSE
# output.directory.base.path <- paste("/home/camp/watkint/working/watkint/data/brastianos/release_v2/","B_", patient, "/exome/CN_het_TRACERx_no_diptest/", sep = "")
# data.format.to.use <- "TRACERx_exome"
# snp.file.path <- paste0(output.directory.base.path, patient, "_TRACERx_input_paths_df.txt")
# seg.file.path <- paste("/home/camp/watkint/working/watkint/data/brastianos/release_v2/", "B_", patient, "/exome/MSAI/", patient, "_joint_segs_all_ascat_info.txt", sep = "")
# overwrite.flag <- FALSE
# gamma.val <- 1
# ##########
#mphase.script.path <- "/home/camp/watkint/working/watkint/scripts/multiregioncopynum/mphase.R"
#TRACERx format for TRACERx example
#tracerx.mphase.command <- paste("Rscript", mphase.script.path, patient, allow.dip.test.phasing.flag, data.format.to.use, snp.file.path, seg.file.path, sep= " ")

# #Debugging option parameters:
# patient <- "ETX014"
# allow.dip.test.phasing.flag  <- FALSE
# # M_ETX014/exome/CN_het_MEDICC_no_diptest
# output.directory.base.path <- paste("/home/camp/watkint/working/watkint/data/murugaesu/release/","M_", patient, "/exome/CN_het_TRACERx_no_diptest/", sep = "")
# data.format.to.use <- "TRACERx_exome"
# snp.file.path <- paste0(output.directory.base.path, patient, "_TRACERx_input_paths_df.txt")
# seg.file.path <- paste("/home/camp/watkint/working/watkint/data/murugaesu/release/", "M_", patient, "/exome/MSAI/", patient, "_joint_segs_all_ascat_info.txt", sep = "")
# ##########
# mphase.script.path <- "/home/camp/watkint/working/watkint/scripts/multiregioncopynum/mphase.R"
# #TRACERx format for TRACERx example
# tracerx.mphase.command <- paste("Rscript", mphase.script.path, patient, allow.dip.test.phasing.flag, data.format.to.use, snp.file.path, seg.file.path, sep= " ")



# # #Debugging option parameters:
# # patient <- "MET048"
# allow.dip.test.phasing.flag  <- FALSE
# output.directory.base.path <- "/Users/watkint/Desktop/MET048_mPhase_data/CN_het_output/"
# data.format.to.use <- "MEDICC"
# snp.file.path <- "/Users/watkint/Desktop/MET048_mPhase_data/MET048_joint_snps.txt"
# seg.file.path <- "/Users/watkint/Desktop/MET048_mPhase_data/MET048_joint_segs_all_ascat_info.txt"
# ##########
# mphase.script.path <- "/Users/watkint/code/multiregioncopynum/mphase.R"
# #MEDICC format for TRACERx example
# medicc.mphase.command <- paste("Rscript", mphase.script.path, patient, allow.dip.test.phasing.flag, data.format.to.use, snp.file.path, seg.file.path, sep=" ")


# patient <- "PTX001"
# allow.dip.test.phasing.flag <- FALSE
# output.directory.base.path <- "/home/camp/watkint/working/watkint/data/gundem/release/G_PTX001/WGS/CN_het_MEDICC_no_diptest/"
# data.format.to.use <- "MEDICC"
# snp.file.path <- "/home/camp/watkint/working/watkint/data/gundem/release/G_PTX001/WGS/MSAI/PTX001_joint_snps.txt"
# seg.file.path <- "/home/camp/watkint/working/watkint/data/gundem/release/G_PTX001/WGS/MSAI/PTX001_joint_segs_all_ascat_info.txt"


# patient <- "BTX001"
# allow.dip.test.phasing.flag <- FALSE
# output.directory.base.path <- "/home/camp/watkint/working/watkint/data/yates/release/Y_BTX001/WGS//CN_het_TRACERx_WGS_no_diptest/"
# data.format.to.use <- "TRACERx_WGS"
# snp.file.path <- "/home/camp/watkint/working/watkint/data/yates/release/Y_BTX001/WGS//CN_het_TRACERx_WGS_no_diptest/Y_BTX001_input_snp_files.txt"
# seg.file.path <- "/home/camp/watkint/working/watkint/data/yates/release/Y_BTX001/WGS/MSAI/BTX001_joint_segs_all_ascat_info.txt"
#overwrite.flag <- FALSE


# patient <- "KTX167"
# allow.dip.test.phasing.flag <- FALSE
# output.directory.base.path <-  "/home/camp/watkint/working/watkint/data/turajlic_wgs/release/T_KTX167/WGS//CN_het_MEDICC_no_diptest/"
# data.format.to.use <-  "MEDICC"
# snp.file.path <- "/home/camp/watkint/working/watkint/data/turajlic_wgs/release/T_KTX167/WGS/MSAI/KTX167_joint_snps.txt"
# seg.file.path <- "/home/camp/watkint/working/watkint/data/turajlic_wgs/release/T_KTX167/WGS/MSAI/KTX167_joint_segs_all_ascat_info.txt"

#"patient"
# patient <- "MTX005"
# allow.dip.test.phasing.flag <- FALSE
# output.directory.base.path <- paste0("/home/camp/watkint/working/watkint/data/martinez/release_v3/R_",patient,"/snp_array//CN_het_MEDICC_no_diptest/")
# data.format.to.use <- "MEDICC"
# snp.file.path <- paste0("/home/camp/watkint/working/watkint/data/martinez/release_v3/R_",patient,"/snp_array/MSAI/",patient,"_new_all_snps.txt")
# seg.file.path <- paste0("/home/camp/watkint/working/watkint/data/martinez/release_v3/R_",patient,"/snp_array/MSAI/",patient,"_joint_segs_all_ascat_info.txt")
# overwrite.flag <- TRUE
# gamma.val <- 0.5



# patient <- "PTX002"
# allow.dip.test.phasing.flag <- FALSE
# output.directory.base.path <- "/home/camp/watkint/working/watkint/data/gundem/release/G_PTX002/WGS/CN_het_TRACERx_WGS_no_diptest/"
# data.format.to.use <- "TRACERx_WGS"
# snp.file.path <- "/home/camp/watkint/working/watkint/data/gundem/release/G_PTX002/WGS/CN_het_TRACERx_WGS_no_diptest/G_PTX002_input_snp_files.txt"
# seg.file.path <- "/home/camp/watkint/working/watkint/data/gundem/release/G_PTX002/WGS/MSAI/PTX002_joint_segs_all_ascat_info.txt"
# overwrite.flag <- FALSE
# gamma.val <- 1


# patient  <-  "STX009"
# allow.dip.test.phasing.flag <- FALSE
# output.directory.base.path <- "/home/camp/watkint/working/watkint/data/sottoriva_snp/release_v3/C_STX009/snp_array//CN_het_MEDICC_no_diptest/"
# data.format.to.use <- "MEDICC"
# snp.file.path <- "/home/camp/watkint/working/watkint/data/sottoriva_snp/release_v3/C_STX009/snp_array/MSAI/STX009_final_all_snps.txt"
# seg.file.path <- "/home/camp/watkint/working/watkint/data/sottoriva_snp/release_v3/C_STX009/snp_array/MSAI/STX009_joint_segs_all_ascat_info.txt"
# overwrite.flag <- TRUE
# gamma.val <- 0.5
# num.cores <- 12



# patient  <-  "STX009"
# allow.dip.test.phasing.flag <- FALSE
# output.directory.base.path <- "/home/camp/watkint/working/watkint/data/sottoriva_snp/release_v3/C_STX009/snp_array//CN_het_MEDICC_no_diptest/"
# data.format.to.use <- "MEDICC"
# snp.file.path <- "/home/camp/watkint/working/watkint/data/sottoriva_snp/release_v3/C_STX009/snp_array/MSAI/STX009_final_all_snps.txt"
# seg.file.path <- "/home/camp/watkint/working/watkint/data/sottoriva_snp/release_v3/C_STX009/snp_array/MSAI/STX009_joint_segs_all_ascat_info.txt"
# overwrite.flag <- FALSE
# gamma.val <- 0.5

# patient  <-  "LTX021"
# allow.dip.test.phasing.flag <- FALSE
# output.directory.base.path <- "/Users/watkint/Desktop/MET048_mPhase_data/CN_het_TRACERx_no_diptest/"
# data.format.to.use <- "TRACERx_exome"
# snp.file.path <- "/Users/watkint/Desktop/MET048_mPhase_data/A_LTX021/exome/CN_het_TRACERx_no_diptest/LTX021_TRACERx_input_paths_df.txt"
# seg.file.path <- "/Users/watkint/Desktop/MET048_mPhase_data/A_LTX021/exome/MSAI/LTX021_joint_segs_all_ascat_info.txt"
# overwrite.flag <- TRUE
# gamma.val <- 1
# num.cores <- 8
# intersection.or.whole.arm <- "intersection"

# patient  <-  "LTX021"
# allow.dip.test.phasing.flag <- FALSE
# output.directory.base.path <- "/camp/project/proj-tracerx-lung/tctProjects/watkint/TRACERx100_lung_pan_can_analysis/release/A_LTX021/exome/CN_het_TRACERx_no_diptest/"
# data.format.to.use <- "TRACERx_exome"
# snp.file.path <- "/camp/project/proj-tracerx-lung/tctProjects/watkint/TRACERx100_lung_pan_can_analysis/release/A_LTX021/exome/CN_het_TRACERx_no_diptest/LTX021_TRACERx_input_paths_df.txt"
# seg.file.path <- "/camp/project/proj-tracerx-lung/tctProjects/watkint/TRACERx100_lung_pan_can_analysis/release/A_LTX021/exome/MSAI/LTX021_joint_segs_all_ascat_info.txt"
# overwrite.flag <- TRUE
# gamma.val <- 1
# num.cores <- 12
# intersection.or.whole.arm <- "intersection"

# patient  <-  "LTX046"
# allow.dip.test.phasing.flag <- FALSE
# output.directory.base.path <- "/camp/project/proj-tracerx-lung/tctProjects/watkint/TRACERx100_lung_pan_can_analysis/release/B_LTX046/exome/CN_het_TRACERx_no_diptest/"
# data.format.to.use <- "TRACERx_exome"
# snp.file.path <- "/camp/project/proj-tracerx-lung/tctProjects/watkint/TRACERx100_lung_pan_can_analysis/release/B_LTX046/exome/CN_het_TRACERx_no_diptest/LTX046_TRACERx_input_paths_df.txt"
# seg.file.path <- "/camp/project/proj-tracerx-lung/tctProjects/watkint/TRACERx100_lung_pan_can_analysis/release/B_LTX046/exome/MSAI/LTX046_joint_segs_all_ascat_info.txt"
# overwrite.flag <- TRUE
# gamma.val <- 1
# num.cores <- 24


# patient  <-  "LTX041"
# allow.dip.test.phasing.flag <- FALSE
# output.directory.base.path <- "/camp/project/proj-tracerx-lung/tctProjects/watkint/TRACERx100_lung_pan_can_analysis/release/L_LTX041/exome/CN_het_TRACERx_no_diptest/"
# data.format.to.use <- "TRACERx_exome"
# snp.file.path <- "/camp/project/proj-tracerx-lung/tctProjects/watkint/TRACERx100_lung_pan_can_analysis/release/L_LTX041/exome/CN_het_TRACERx_no_diptest/LTX041_TRACERx_input_paths_df.txt"
# seg.file.path <- "/camp/project/proj-tracerx-lung/tctProjects/watkint/TRACERx100_lung_pan_can_analysis/release/L_LTX041/exome/MSAI/LTX041_joint_segs_all_ascat_info.txt"
# overwrite.flag <- TRUE
# gamma.val <- 1
# num.cores <- 24



# patient  <-  "LTX181"
# allow.dip.test.phasing.flag <- FALSE
# output.directory.base.path <- "/camp/project/proj-tracerx-lung/tctProjects/watkint/TRACERx100_lung_pan_can_analysis/release/M_LTX181/exome/CN_het_TRACERx_no_diptest/"
# data.format.to.use <- "TRACERx_exome"
# snp.file.path <- "/camp/project/proj-tracerx-lung/tctProjects/watkint/TRACERx100_lung_pan_can_analysis/release/M_LTX181/exome/CN_het_TRACERx_no_diptest/LTX181_TRACERx_input_paths_df.txt"
# seg.file.path <- "/camp/project/proj-tracerx-lung/tctProjects/watkint/TRACERx100_lung_pan_can_analysis/release/M_LTX181/exome/MSAI/LTX181_joint_segs_all_ascat_info.txt"
# overwrite.flag <- TRUE
# gamma.val <- 1
# num.cores <- 24


# patient <- "MTX002"
# allow.dip.test.phasing.flag <- FALSE
# output.directory.base.path <- "/home/camp/watkint/working/watkint/data/martinez/release_v3/R_MTX002/snp_array//CN_het_MEDICC_no_diptest/"
# data.format.to.use<- "MEDICC"
# snp.file.path<- "/home/camp/watkint/working/watkint/data/martinez/release_v3/R_MTX002/snp_array/MSAI/MTX002_final_all_snps.txt"
# seg.file.path<- "/home/camp/watkint/working/watkint/data/martinez/release_v3/R_MTX002/snp_array/MSAI/MTX002_joint_segs_all_ascat_info.txt"
# overwrite.flag<-TRUE
# gamma.val <-0.5
# num.cores <- 12

# patient <- "BTX003"
# allow.dip.test.phasing.flag <- FALSE
# output.directory.base.path <- "/camp/lab/swantonc/working/watkint/data/loi/release/S_BTX003/exome//CN_het_TRACERx_no_diptest/"
# data.format.to.use <- "TRACERx_exome"
# snp.file.path <- "/camp/lab/swantonc/working/watkint/data/loi/release/S_BTX003/exome//CN_het_TRACERx_no_diptest//BTX003_TRACERx_input_paths_df.txt"
# seg.file.path <- "/camp/lab/swantonc/working/watkint/data/loi/release/S_BTX003/exome//MSAI/BTX003_joint_segs_all_ascat_info.txt"
# overwrite.flag <- TRUE
# gamma.val <- 1
# num.cores <- 12


# patient <- "BTX004"
# allow.dip.test.phasing.flag <- FALSE
# output.directory.base.path <- "/home/camp/watkint/working/watkint/data/yates/release/Y_BTX004/WGS//CN_het_TRACERx_WGS_no_diptest/"
# data.format.to.use <- "TRACERx_WGS"
# snp.file.path <- "/home/camp/watkint/working/watkint/data/yates/release/Y_BTX004/WGS//CN_het_TRACERx_WGS_no_diptest/Y_BTX004_input_snp_files.txt"
# seg.file.path <- "/home/camp/watkint/working/watkint/data/yates/release/Y_BTX004/WGS/MSAI/BTX004_joint_segs_all_ascat_info.txt"
# overwrite.flag <- TRUE
# gamma.val <- 1
# num.cores <- 12


# patient <- "RLX007"
# allow.dip.test.phasing.flag <- FALSE
# output.directory.base.path <- "/camp/lab/swantonc/working/lime/projects/PanCan/renal/gerlinger/release/G_RLX007/exome//CN_het_TRACERx_no_diptest/"
# data.format.to.use <- "TRACERx_exome"
# snp.file.path <- "/camp/lab/swantonc/working/lime/projects/PanCan/renal/gerlinger/release/G_RLX007/exome//CN_het_TRACERx_no_diptest//RLX007_TRACERx_input_paths_df.txt"
# seg.file.path <- "/camp/lab/swantonc/working/lime/projects/PanCan/renal/gerlinger/release/G_RLX007/exome//MSAI/RLX007_joint_segs_all_ascat_info.txt"
# overwrite.flag <- TRUE
# gamma.val <- 1
# num.cores <- 12


# patient <- "KTX128"
# allow.dip.test.phasing.flag <- FALSE
# output.directory.base.path <- "/camp/project/proj-tracerx-lung/tctProjects/lime/projects/PanCan/renal/turajlic_txpipeline/release/E_KTX128/exome//CN_het_TRACERx_no_diptest/"
# data.format.to.use <- "TRACERx_exome"
# snp.file.path <- "/camp/project/proj-tracerx-lung/tctProjects/lime/projects/PanCan/renal/turajlic_txpipeline/release/E_KTX128/exome//CN_het_TRACERx_no_diptest//KTX128_TRACERx_input_paths_df.txt"
# seg.file.path <- "/camp/project/proj-tracerx-lung/tctProjects/lime/projects/PanCan/renal/turajlic_txpipeline/release/E_KTX128/exome//MSAI/KTX128_joint_segs_all_ascat_info.txt"
# overwrite.flag <- TRUE
# gamma.val <- 1
# num.cores <- 12


# patient <- "KTX096"
# allow.dip.test.phasing.flag <- FALSE
# output.directory.base.path <- "/camp/project/proj-tracerx-lung/tctProjects/lime/projects/PanCan/renal/turajlic_txpipeline/release/E_KTX096/exome//CN_het_TRACERx_no_diptest/"
# data.format.to.use <- "TRACERx_exome"
# snp.file.path <- "/camp/project/proj-tracerx-lung/tctProjects/lime/projects/PanCan/renal/turajlic_txpipeline/release/E_KTX096/exome//CN_het_TRACERx_no_diptest//KTX096_TRACERx_input_paths_df.txt"
# seg.file.path <- "/camp/project/proj-tracerx-lung/tctProjects/lime/projects/PanCan/renal/turajlic_txpipeline/release/E_KTX096/exome//MSAI/KTX096_joint_segs_all_ascat_info.txt"
# overwrite.flag <- TRUE
# gamma.val <- 1
# num.cores <- 12


# patient <- "PTX006"
# allow.dip.test.phasing.flag <- FALSE
# output.directory.base.path <- "/home/camp/watkint/working/watkint/data/gundem/release/G_PTX006/WGS//CN_het_TRACERx_WGS_no_diptest/"
# data.format.to.use <- "TRACERx_WGS"
# snp.file.path <- "/home/camp/watkint/working/watkint/data/gundem/release/G_PTX006/WGS//CN_het_TRACERx_WGS_no_diptest/G_PTX006_input_snp_files.txt"
# seg.file.path <- "/home/camp/watkint/working/watkint/data/gundem/release/G_PTX006/WGS/MSAI/PTX006_joint_segs_all_ascat_info.txt"
# overwrite.flag <- TRUE
# gamma.val <- 1
# num.cores <- 12

# patient <- "PTX007"
# allow.dip.test.phasing.flag <- FALSE
# output.directory.base.path <- "/home/camp/watkint/working/watkint/data/gundem/release/G_PTX007/WGS//CN_het_TRACERx_WGS_no_diptest/"
# data.format.to.use <- "TRACERx_WGS"
# snp.file.path <- "/home/camp/watkint/working/watkint/data/gundem/release/G_PTX007/WGS//CN_het_TRACERx_WGS_no_diptest/G_PTX007_input_snp_files.txt"
# seg.file.path <- "/home/camp/watkint/working/watkint/data/gundem/release/G_PTX007/WGS/MSAI/PTX007_joint_segs_all_ascat_info.txt"
# overwrite.flag <- TRUE
# gamma.val <- 1
# num.cores <- 12


# patient <- "KTX139"
# allow.dip.test.phasing.flag <- FALSE
# output.directory.base.path <- "/camp/project/proj-tracerx-lung/tctProjects/lime/projects/PanCan/renal/turajlic_txpipeline/release/E_KTX139/exome//CN_het_TRACERx_no_diptest/"
# data.format.to.use <- "TRACERx_exome"
# snp.file.path <- "/camp/project/proj-tracerx-lung/tctProjects/lime/projects/PanCan/renal/turajlic_txpipeline/release/E_KTX139/exome//CN_het_TRACERx_no_diptest//KTX139_TRACERx_input_paths_df.txt"
# seg.file.path <- "/camp/project/proj-tracerx-lung/tctProjects/lime/projects/PanCan/renal/turajlic_txpipeline/release/E_KTX139/exome//MSAI/KTX139_joint_segs_all_ascat_info.txt"
# overwrite.flag <- TRUE
# gamma.val <- 1
# num.cores <- 12


# patient <- "KTX146"
# allow.dip.test.phasing.flag <- FALSE
# output.directory.base.path <- "/camp/project/proj-tracerx-lung/tctProjects/lime/projects/PanCan/renal/turajlic_txpipeline/release/E_KTX139/exome//CN_het_TRACERx_no_diptest/"
# data.format.to.use <- "TRACERx_exome"
# snp.file.path <- "/camp/project/proj-tracerx-lung/tctProjects/lime/projects/PanCan/renal/turajlic_txpipeline/release/E_KTX146/exome//CN_het_TRACERx_no_diptest//KTX146_TRACERx_input_paths_df.txt"
# seg.file.path <- "/camp/project/proj-tracerx-lung/tctProjects/lime/projects/PanCan/renal/turajlic_txpipeline/release/E_KTX146/exome//MSAI/KTX146_joint_segs_all_ascat_info.txt"
# overwrite.flag <- TRUE
# gamma.val <- 1
# num.cores <- 12


# patient <- "PTX011"
# allow.dip.test.phasing.flag <- FALSE
# output.directory.base.path <- "/home/camp/watkint/working/watkint/data/gundem/release/G_PTX011/WGS//CN_het_TRACERx_WGS_no_diptest/"
# data.format.to.use <- "TRACERx_WGS"
# snp.file.path <- "/home/camp/watkint/working/watkint/data/gundem/release/G_PTX011/WGS//CN_het_TRACERx_WGS_no_diptest/G_PTX011_input_snp_files.txt"
# seg.file.path <- "/home/camp/watkint/working/watkint/data/gundem/release/G_PTX011/WGS/MSAI/PTX011_joint_segs_all_ascat_info.txt"
# overwrite.flag <- TRUE
# gamma.val <- 1
# num.cores <- 12


# patient <- "KTX167"
# allow.dip.test.phasing.flag <- FALSE
# output.directory.base.path <- "/home/camp/watkint/working/watkint/data/turajlic_wgs/release/T_KTX167/WGS//CN_het_TRACERx_no_diptest/"
# data.format.to.use <- "TRACERx_WGS"
# snp.file.path <- "/home/camp/watkint/working/watkint/data/turajlic_wgs/release/T_KTX167/WGS//CN_het_TRACERx_no_diptest/T_KTX167_input_snp_files.txt"
# seg.file.path <- "/home/camp/watkint/working/watkint/data/turajlic_wgs/release/T_KTX167/WGS/MSAI/KTX167_joint_segs_all_ascat_info.txt"
# overwrite.flag <- TRUE
# gamma.val <- 1
# num.cores <- 16


# patient <- "KTX169"
# allow.dip.test.phasing.flag <- FALSE
# output.directory.base.path <- "/home/camp/watkint/working/watkint/data/turajlic_wgs/release/T_KTX169/WGS//CN_het_TRACERx_no_diptest/"
# data.format.to.use <- "TRACERx_WGS"
# snp.file.path <- "/home/camp/watkint/working/watkint/data/turajlic_wgs/release/T_KTX169/WGS//CN_het_TRACERx_no_diptest/T_KTX169_input_snp_files.txt"
# seg.file.path <- "/home/camp/watkint/working/watkint/data/turajlic_wgs/release/T_KTX169/WGS/MSAI/KTX169_joint_segs_all_ascat_info.txt"
# overwrite.flag <- TRUE
# gamma.val <- 1


# patient  <-  "STX007"
# allow.dip.test.phasing.flag <- FALSE
# output.directory.base.path <- "/home/camp/watkint/working/watkint/data/sottoriva_snp/release_v3/C_STX007/snp_array//CN_het_MEDICC_no_diptest/"
# data.format.to.use <- "MEDICC"
# snp.file.path <- "/home/camp/watkint/working/watkint/data/sottoriva_snp/release_v3/C_STX007/snp_array/MSAI/STX007_final_all_snps.txt"
# seg.file.path <- "/home/camp/watkint/working/watkint/data/sottoriva_snp/release_v3/C_STX007/snp_array/MSAI/STX007_joint_segs_all_ascat_info.txt"
# overwrite.flag <- TRUE
# gamma.val <- 0.5
# num.cores <- 12

# patient <- "KTX027"
# allow.dip.test.phasing.flag <- FALSE
# output.directory.base.path <- "/home/camp/watkint/working/watkint/data/turajlic_wgs/release/T_KTX027/WGS//CN_het_TRACERx_no_diptest/"
# data.format.to.use <- "TRACERx_WGS"
# snp.file.path <- "/home/camp/watkint/working/watkint/data/turajlic_wgs/release/T_KTX027/WGS//CN_het_TRACERx_no_diptest/T_KTX027_input_snp_files.txt"
# seg.file.path <- "/home/camp/watkint/working/watkint/data/turajlic_wgs/release/T_KTX027/WGS/MSAI/KTX027_joint_segs_all_ascat_info.txt"
# overwrite.flag <- TRUE
# gamma.val <- 1
# num.cores <- 12


# patient  <-  "LTX121"
# allow.dip.test.phasing.flag <- FALSE
# output.directory.base.path <- "/camp/project/proj-tracerx-lung/tctProjects/watkint/TRACERx100_lung_pan_can_analysis/release/A_LTX121/exome/CN_het_TRACERx_no_diptest/"
# data.format.to.use <- "TRACERx_exome"
# snp.file.path <- "/camp/project/proj-tracerx-lung/tctProjects/watkint/TRACERx100_lung_pan_can_analysis/release/A_LTX121/exome/CN_het_TRACERx_no_diptest/LTX121_TRACERx_input_paths_df.txt"
# seg.file.path <- "/camp/project/proj-tracerx-lung/tctProjects/watkint/TRACERx100_lung_pan_can_analysis/release/A_LTX121/exome/MSAI/LTX121_joint_segs_all_ascat_info.txt"
# overwrite.flag <- TRUE
# gamma.val <- 1
# num.cores <- 12
# intersection.or.whole.arm <- "intersection"


patient <- "MTX006"
allow.dip.test.phasing.flag <- FALSE
output.directory.base.path <- paste0("/home/camp/watkint/working/watkint/data/martinez/release_v3/R_",patient,"/snp_array//CN_het_MEDICC_no_diptest/")
data.format.to.use <- "MEDICC"
snp.file.path <- paste0("/home/camp/watkint/working/watkint/data/martinez/release_v3/R_",patient,"/snp_array/MSAI/",patient,"_new_all_snps.txt")
seg.file.path <- paste0("/home/camp/watkint/working/watkint/data/martinez/release_v3/R_",patient,"/snp_array/MSAI/",patient,"_joint_segs_all_ascat_info.txt")
overwrite.flag <- TRUE
gamma.val <- 0.5
num.cores <- 12

print("patient")
print(patient)
print("allow.dip.test.phasing.flag")
print(allow.dip.test.phasing.flag)
print("output.directory.base.path")
print(output.directory.base.path)
print("data.format.to.use)")
print(data.format.to.use)
print("snp.file.path")
print(snp.file.path)
print("seg.file.path")
print(seg.file.path)
print("overwrite.flag")
print(overwrite.flag)
print("gamma.val")
print(gamma.val)
print("num.cores")
print(num.cores)

#We need to be able to select the intersection of all regions or a full GISTIC compatible union with defined regions that will always be covered
#What proportion of the chosen genome, in this case hg19, arms are covered by the segments that have been inputted.
#We also need to add in the wgii calculations


###Source relevant scripts, best to do this at the top to check for errors...
if( is.null(current.dir) ) {

	current.dir <- ""

} else {

	current.dir <- paste0(current.dir, "/")

}
source(paste0(current.dir, "utility_MPhase_functions.R"))
source(paste0(current.dir, "data_input_and_QC_functions_v2.R"))
source(paste0(current.dir, "minimum_consistent_segmentation_functions.R"))
source(paste0(current.dir, "consensus_phasing_from_multiple_samples_functions_v2.R"))
source(paste0(current.dir, "generation_of_phased_and_total_cpn_functions_v2.R"))
source(paste0(current.dir, "MSAI_detection_functions.R"))
source(paste0(current.dir, "pairwise_cpn_comparison_and_heterogeneity_functions.R"))
source(paste0(current.dir, "cancer_cell_fraction_calculation_functions_v2.R"))
source(paste0(current.dir, "final_output_segmentation_file_generation_functions_v3.R"))
source(paste0(current.dir, "segment_joining_functions.R"))
source(paste0(current.dir, "cpn_heterogeneity_plotting_functions_v2.R"))

#####Get relevant libraries:
#Todo check this is the only required library.
suppressWarnings(suppressMessages(library(diptest)))
suppressWarnings(suppressMessages(library(gplots)))
suppressWarnings(suppressMessages(library(gdata)))
suppressWarnings(suppressMessages(library(GenomicRanges)))
suppressWarnings(suppressMessages(library(doParallel)))

###############This is past the input stage and into the analysis:
#Set constants for all analyses here... eventually these can be command line set
MIN_NUM_HZ_SNP_PHASING <- 5
#PLOT_OUTPUT_DIR <- NULL
OVERWRITE_FLAG <- as.logical(overwrite.flag)
PLOT_OUTPUT_DIR <- paste(output.directory.base.path, "/", sep="")
if ( ! file.exists(PLOT_OUTPUT_DIR) ) {

	print(PLOT_OUTPUT_DIR)
	print("directory doesn't exist")

}
NUM_CORES <- as.numeric(num.cores)


registerDoParallel(cores = NUM_CORES)
print(foreach::getDoParWorkers())
print(getDoParName())
print(getDoParVersion())


INTERSECTION_OR_ARM <- intersection.or.whole.arm
GAMMA_VAL <- as.numeric(gamma.val)
MAX_CPN_FOR_SUBCLONALITY <- 10
MIN_CCF_TO_CONSIDER <- 0.15
MIN_CPN_TO_CONSIDER <- 0.15
MIN_NUM_CHRS <- 3
PVAL_THRESH_VS_MIN_CCF <- 0.001
PVAL_THRESH_VS_OTH_CHRS_SAME_CPN <- 0.01
GAIN_LOGR_THRESH <- log2(2.5 / 2)
LOSS_LOGR_THRESH <- log2(1.5 / 2)
AMP_LOGR_THRESH <- log2(4 / 2)
DEEP_LOSS_LOGR_THRESH <- log2(0.5 / 2)
PVAL_THRESH_VS_LOGR_FOR_CPN_EVENT <- 0.01
CHRS_TO_USE_VEC <- 1:22
CONF_INT_FOR_CPNX_BAF_VS_EXP <- 0.9999
PW_CPNA_CPNB_DIFF_AMPLITUDE_THRESH <- 0.001
PVAL_THRESH_TO_JOIN_PH_SEGS <- 0.001
MIN_PROP_SNP_SEP_THRESH <- 0.75
PVAL_COMB_WILCOX_PH_THRESH <- 0.01
PVAL_DIPTEST_AI_THRESH <- 0.01
ALLOW_DIPTEST_ALONE_PH_FLAG <- allow.dip.test.phasing.flag
MSAI_BAF_LTNPD_PVAL_THRESHOLD <- 0.001
MSAI_CPNX_PVAL_THRESHOLD <- 0.001
PVAL_THRESH_CONVERT_GAIN_TO_AMP  <- 0.01

if (is.null(output.directory.base.path)) {

	stop("No output path provided, aborting...")

}

#so now we need to read in the MSAI data:
if ( ! file.exists(snp.file.path)) {

	stop("No file at the snp.file.path provided in the second argument")

}

if ( ! file.exists(seg.file.path)) {

	stop("No file at the seg.file.path provided in the third argument")

}

###Load in the relevant genome file, #todo make this an option that can be set.
hg19.info.path <- paste0(current.dir, "ucsc_hg19_chrominfo.txt")
hg19.info.df <- read.table(hg19.info.path, stringsAsFactors = FALSE, sep = "\t", header = TRUE)
chrom.info.df <- hg19.info.df 

#Paths for saving/loading objects depending on OVERWRITE_FLAG
#Reading in data object paths
asas.SNPS.lr.and.or.baf.lists.path <- paste0(PLOT_OUTPUT_DIR, patient, "_asas_SNPS_lr_and_or_baf_lists_lists.RDS")
asas.SNPs.w.both.lr.baf.snp.lists.path <- paste0(PLOT_OUTPUT_DIR, patient, "_asas_SNPs_w_both_lr_baf_snp_lists.RDS")
min.cons.segs.df.write.path <- paste0(PLOT_OUTPUT_DIR, patient, "_min_cons_segs_df.txt")
sorted.all.bp.all.samps.segs.anno.df.path <- paste0(PLOT_OUTPUT_DIR, patient, "_sorted_all_bp_all_samps_segs_anno.df.txt")
min.cons.segs.df.path <- paste0(PLOT_OUTPUT_DIR, patient, "_min_cons_segs_df.RDS")
processed.orig.cpn.segs.df.path <- paste0(PLOT_OUTPUT_DIR, patient, "_processed_orig_cpn_segs_df.RDS")

ploidy.est.vec.path <- paste0(PLOT_OUTPUT_DIR, patient, "_ploidy_est_vec.RDS")
purity.est.vec.path <- paste0(PLOT_OUTPUT_DIR, patient, "_purity_est_vec.RDS")
#Analysis object paths
all.seg.phasing.info.list.path <- paste0(PLOT_OUTPUT_DIR, patient, "_all_seg_phasing_info_list.RDS")
patient.phased.BAF.df.list.path <- paste0(PLOT_OUTPUT_DIR, patient, "_phased_BAF_df_list.RDS")
all.seg.MSAI.by.BAF.ltnpd.present.list.path <- paste0(PLOT_OUTPUT_DIR, patient, "_all_seg_MSAI_by_BAF_ltnpd_present_list.RDS")
all.samp.non.hz.logr.and.tot.cpn.df.list.path <- paste0(PLOT_OUTPUT_DIR, patient, "_all_samp_non_hz_logr_and_tot_cpn_df_list.RDS")
asas.ph.raw.cpn.lists.path <- paste0(PLOT_OUTPUT_DIR, patient, "_asas_ph_raw_cpn_lists.RDS")
samps.int.cpn.raw.cpn.by.chr.list.path <- paste0(PLOT_OUTPUT_DIR, patient, "_samps_int_cpn_raw_cpn_by_chr_list.RDS")
all.seg.CCF.list.path <- paste0(PLOT_OUTPUT_DIR, patient, "_all_seg_CCF_list.RDS")
all.seg.int.ob.xp.cpn.df.list.path <- paste0(PLOT_OUTPUT_DIR, patient, "_all_seg_int_ob_xp_cpn_df_list.RDS")
all.seg.CCF.diff.from.exp.list.path <- paste0(PLOT_OUTPUT_DIR, patient, "_all_seg_CCF_diff_from_exp_list.RDS")
all.seg.clonal.in.all.samps.list.path <- paste0(PLOT_OUTPUT_DIR, patient, "_all_seg_clonal_in_all_samps_list.RDS")
all.seg.greater.than.two.diff.list.path <- paste0(PLOT_OUTPUT_DIR, patient, "_all_seg_greater_than_two_diff_list.RDS")
print("Finished generating paths...")
#Load in the data: 
#Note that the afunctions to create the single dataframes that it accepts are in the data_input_and_QC_functions.R
if (OVERWRITE_FLAG == TRUE) {

	f.joint.segs.df <- read.table(seg.file.path, stringsAsFactors = FALSE, sep = "\t", header = TRUE)
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


	sample.vec <- unique(f.joint.segs.df$sample)
	ploidy.est.vec <- f.joint.segs.df[ ! duplicated(f.joint.segs.df$sample), 'Ploidy']
	names(ploidy.est.vec) <- unique(f.joint.segs.df$sample)
	purity.est.vec <- f.joint.segs.df[ ! duplicated(f.joint.segs.df$sample), 'ACF']
	names(purity.est.vec) <- unique(f.joint.segs.df$sample)

	if (data.format.to.use == "MEDICC") {

		print("data format is MEDICC")
		snp.file.df <- read.table(snp.file.path, stringsAsFactors = FALSE, sep = "\t", header = TRUE)
		all.samp.SNP.list <- readInMEDICCInputFormatParallel(snp.file.df = snp.file.df,
															processed.orig.cpn.segs.df = f.joint.segs.df,
															sample.vec = sample.vec,
															num.cores = NUM_CORES)

	} else if (data.format.to.use == "TRACERx_exome" | data.format.to.use == "TRACERx_snp_array") {

		if (data.format.to.use == "TRACERx_exome") {

			array.type <- "exome"

		} else if (data.format.to.use == "TRACERx_snp_array") {

			array.type <- "snp_array"
		}
		#So we want to read in a file that looks like this:	#patient_id sample_id ascat_bc_file_path 
		pt.snp.file.paths.df <- read.table(snp.file.path, stringsAsFactors = FALSE, sep = "\t", header = TRUE)
		print("Reading in SNP data...")
		all.samp.SNP.list <- getSNPDataForTRACERxExomeSNPArray(pt.snp.file.paths.df)



	} else if (data.format.to.use == "TRACERx_WGS") {

		array.type <- "WGS"
		print("Reading in WGS SNP data...")
		all.samp.SNP.list <- getSNPDataForTRACERxWGS(snp.file.path)


	} else {

		print(data.format.to.use)
		stop("Incorrect data.format.to.use input to mphase")

	}

	#So we need the centromere postions for each arm:
	#Add in p_end and q_start


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
	for ( chr in CHRS_TO_USE_VEC ) {

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

		}#end 		for (new.bp.gr in bps.to.introduce.grs.list) {

	}#end for ( chr in CHRS_TO_USE_VEC ) {


	#So actually defining the segments that are lost or gained makes sense here..
	print("Calculating arm level gain, loss, \"neutral\", and LOH events in the original segmentation...")
	all.samp.w.metadata.grs <- new.all.samp.w.metadata.grs
	gain.loss.name.list <- list(c("neutral"), c("gain", "amp"), c("loss", "deep_loss"),  c("LOH"))
	gain.loss.name.list
	names(gain.loss.name.list) <- c("neutral", "gain", "loss", "LOH")
	p.or.q.vec <- c("p", "q")
	for ( chr in CHRS_TO_USE_VEC ) {
		#chr <- 21
		chr.gr.list <- getChromosomeArmsGRangeList(chr, chrom.info.df)
		#print(chr)
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

	}#end	for ( chr in CHRS_TO_USE_VEC ) {

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
	#corresponding segment, dummy or otherise 
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

	 # 	curr.all.bp.chr.df <- all.bp.chr.grs[i]

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

	cpn.chr.start.end.df <- sorted.all.bp.all.samps.segs.anno.df[sorted.all.bp.all.samps.segs.anno.df$full_sample==sorted.all.bp.all.samps.segs.anno.df[1,]$full_sample,1:3]
	mock.cpn.matrix <- matrix(data=0, nrow=nrow(cpn.chr.start.end.df), ncol = length(sample.vec)) 
	colnames(mock.cpn.matrix) <- sample.vec
	min.cons.segs.df <- cbind(cpn.chr.start.end.df, mock.cpn.matrix)

	#All SNPs in these lists must have both a logr and baf value: all sample all seg logr baf lists is asas.lr.baf
	ptm <- proc.time()
	print("cleanAndErrorCheckSNPsParallel first")
	source(paste0(current.dir, "data_input_and_QC_functions_v2.R"))
	asas.SNPs.w.both.lr.baf.snp.lists <- cleanAndErrorCheckSNPsParallel(all.samp.SNP.list,
																		min.cons.segs.df,
																		allow.logr.only.flag = FALSE,
																		data.type = data.format.to.use)
	print(asas.SNPs.w.both.lr.baf.snp.lists.path)
	saveRDS(asas.SNPs.w.both.lr.baf.snp.lists, file = asas.SNPs.w.both.lr.baf.snp.lists.path)
	time.taken <- proc.time() - ptm
	print("time.taken")
	print(time.taken)
	##
	#SNPs that only have LogR are permitted in these lists, but it also contains all the SNPs with both BAF and LogR
	ptm <- proc.time()
	print("cleanAndErrorCheckSNPsParallel second")
	asas.SNPS.lr.and.or.baf.lists <- cleanAndErrorCheckSNPsParallel(all.samp.SNP.list,
																	min.cons.segs.df,
																	allow.logr.only.flag = TRUE, 
																	data.type = data.format.to.use)
	print(asas.SNPS.lr.and.or.baf.lists.path)
	saveRDS(asas.SNPS.lr.and.or.baf.lists, file = asas.SNPS.lr.and.or.baf.lists.path)
	time.taken <- proc.time() - ptm
	print("time.taken")
	print(time.taken)

	if ( length(asas.SNPs.w.both.lr.baf.snp.lists) <= 1 ) {

		print(str(asas.SNPs.w.both.lr.baf.snp.lists, max = 1))
		stop("Only one or fewer samples present after reading in data.")

	}

	#save:
	print("Saving snps files....")
	print(asas.SNPS.lr.and.or.baf.lists.path)
	saveRDS(asas.SNPS.lr.and.or.baf.lists, file = asas.SNPS.lr.and.or.baf.lists.path)
	print(asas.SNPs.w.both.lr.baf.snp.lists.path)
	saveRDS(asas.SNPs.w.both.lr.baf.snp.lists, file = asas.SNPs.w.both.lr.baf.snp.lists.path)
	print(min.cons.segs.df.path)
	saveRDS(min.cons.segs.df, file = min.cons.segs.df.path)
	print(min.cons.segs.df.write.path)
	saveRDS(ploidy.est.vec, file = ploidy.est.vec.path)
	saveRDS(purity.est.vec, file = purity.est.vec.path)
	saveRDS(f.joint.segs.df, file = processed.orig.cpn.segs.df.path)

	saveRDS(sorted.all.bp.all.samps.segs.anno.df, file = sorted.all.bp.all.samps.segs.anno.df.path)
	
} else if (OVERWRITE_FLAG == FALSE) {#end if (OVERWRITE_FLAG == TRUE) {


	#load original data:
	asas.SNPS.lr.and.or.baf.lists <- readRDS(file = asas.SNPS.lr.and.or.baf.lists.path)
	asas.SNPs.w.both.lr.baf.snp.lists <- readRDS(file = asas.SNPs.w.both.lr.baf.snp.lists.path)
	min.cons.segs.df <- readRDS(file = min.cons.segs.df.path)
	f.joint.segs.df <- readRDS(file = processed.orig.cpn.segs.df.path)
	ploidy.est.vec <- readRDS(file = ploidy.est.vec.path)
	purity.est.vec <- readRDS(file = purity.est.vec.path)

}#end else if (OVERWRITE_FLAG == FALSE


#if we're making a 
# Now we need to make sure that the asas.lr.baf.lists and sample.to.phase.seg.list 
# if (length(names(sample.to.phase.seg.list)) != length(asas.lr.baf.lists[[1]])) {

#   print("Missing or removed segments:")
#   print(names(sample.to.phase.seg.list)[ ! names(asas.lr.baf.lists[[1]]) %in% names(sample.to.phase.seg.list)])
#   print("Continuing creation of phased BAF dataframe")

# }

#Todo, Add function definitions documentaiton
#Todo, clarify the loops and next statements that are controlling what happens when there are less than 5 hz snps in a sample.
#Todo, check if the mixture phasing is ever being selected... in the example with B_MET002 no mixture phasing called.
seg.name.vec <- paste(min.cons.segs.df[, 1], min.cons.segs.df[, 2], min.cons.segs.df[, 3], sep = ":")
sample.vec <- names(asas.SNPs.w.both.lr.baf.snp.lists)

if (OVERWRITE_FLAG == TRUE) {

	print("Getting phasing to apply to all samples...")
	source(paste0(current.dir, "consensus_phasing_from_multiple_samples_functions_v2.R"))
	all.seg.phasing.info.list <- getCombinedBAFDipTestPhasing(asas.lr.baf.lists = asas.SNPs.w.both.lr.baf.snp.lists,
	                                            seg.name.vec = seg.name.vec,
	                                            ploidy.est.vec = ploidy.est.vec,
	                                            purity.est.vec = purity.est.vec,
	                                            gamma.val = GAMMA_VAL,
	                                            min.num.hz.snp.phasing = MIN_NUM_HZ_SNP_PHASING,
	                                            min.prop.snp.sep.threshold = MIN_PROP_SNP_SEP_THRESH,
	                                            dip.test.p.val.threshold = PVAL_DIPTEST_AI_THRESH,
	                                            comb.wilcox.snp.sep.p.val.threshold = PVAL_COMB_WILCOX_PH_THRESH,
	                                            allow.dip.test.without.ASCAT.to.phase = ALLOW_DIPTEST_ALONE_PH_FLAG
	                	                            )
	print("Getting phasing to apply to all samples using getCombinedBAFDipTestPhasingParallel...")
	source(paste0(current.dir, "consensus_phasing_from_multiple_samples_functions_v2.R"))
	# all.seg.phasing.info.list <- getCombinedBAFDipTestPhasingParallel(asas.lr.baf.lists = asas.SNPs.w.both.lr.baf.snp.lists,
	#                                             seg.name.vec = seg.name.vec,
	#                                             ploidy.est.vec = ploidy.est.vec,
	#                                             purity.est.vec = purity.est.vec,
	#                                             gamma.val = GAMMA_VAL,
	#                                             min.num.hz.snp.phasing = MIN_NUM_HZ_SNP_PHASING,
	#                                             min.prop.snp.sep.threshold = MIN_PROP_SNP_SEP_THRESH,
	#                                             dip.test.p.val.threshold = PVAL_DIPTEST_AI_THRESH,
	#                                             comb.wilcox.snp.sep.p.val.threshold = PVAL_COMB_WILCOX_PH_THRESH,
	#                                             allow.dip.test.without.ASCAT.to.phase = ALLOW_DIPTEST_ALONE_PH_FLAG
	#                                             )
	
	saveRDS(all.seg.phasing.info.list, file = all.seg.phasing.info.list.path)


	#Todo unsure of what to do with the baf logr lists which have positions that only have logr here, at the moment exclude them but how
	#does that affect later run
	print("Generating phased BAF dataframes...")
	patient.phased.BAF.df.list <- createAllSegsPhasedBAFDataFrame(asas.lr.baf.lists = asas.SNPs.w.both.lr.baf.snp.lists,
                                                         sample.to.phase.seg.list = all.seg.phasing.info.list,
                                                         min.num.hz.snp.phasing = MIN_NUM_HZ_SNP_PHASING
                                                         )
	if(length(patient.phased.BAF.df.list) <= 1) {

		print("Error in generating phased BAF dataframes")
	}
	saveRDS(patient.phased.BAF.df.list, file = patient.phased.BAF.df.list.path)


	#Detect MSAI in the original BAF before it is processed alongside logr to get the Raw allele specifc cpn to be used in the rest of the analysis.
	print("Performing LTNPD MSAI detection on the BAF...")
	all.seg.MSAI.by.BAF.ltnpd.present.list <- getAllSegsMSAIDetectionFromBAFListLTNPD(seg.name.vec = seg.name.vec,
																				  sample.vec  = sample.vec,
																				  asas.SNPS.lr.and.or.baf.lists = asas.SNPS.lr.and.or.baf.lists,
																				  all.seg.phasing.info.list = all.seg.phasing.info.list,
																				  MSAI.BAF.LTPND.pval.thresh = MSAI_BAF_LTNPD_PVAL_THRESHOLD)
	saveRDS(all.seg.MSAI.by.BAF.ltnpd.present.list, file = all.seg.MSAI.by.BAF.ltnpd.present.list.path)

	#Todo add another method of detecting MSAI using only the BAF and a paired T test.
	# all.seg.MSAI.by.BAF.ltnpd.present.list <- getAllSegsMirroredSubclonalAllelicImbalanceDetectionFromBAFList(seg.name.vec = seg.name.vec,
	# 																										   sample.vec  = sample.vec,
	# 																										   asas.SNPS.lr.and.or.baf.lists = asas.SNPS.lr.and.or.baf.lists,
	# 																										   all.seg.phasing.info.list = all.seg.phasing.info.list,
	# 																										   MSAI.BAF.LTPND.pval.thresh = MSAI_BAF_LTNPD_PVAL_THRESHOLD)

	print("Generating datapoints with LogR total cpn dataframes...")
	#This is where the SNP filtering occurs. We can fix this by adding by passing an option of which data type is being analysed...
	all.samp.non.hz.logr.and.tot.cpn.df.list <- getAllSampLogrAlonePositionsLogrAndTotalCPNDFsList(sample.vec = sample.vec,
																				asas.SNPS.lr.and.or.baf.lists = asas.SNPS.lr.and.or.baf.lists,
																				ploidy.est.vec = ploidy.est.vec,
																				purity.est.vec = purity.est.vec,
																				gamma.val = GAMMA_VAL,
																				data.format.to.use = data.format.to.use)
	saveRDS(all.samp.non.hz.logr.and.tot.cpn.df.list, file = all.samp.non.hz.logr.and.tot.cpn.df.list.path)



	#Add this the unphaseable stuff too maybe? Not necessary.
	print("Applying phasing to all samples to get phased raw cpn...")
	asas.ph.raw.cpn.lists <- getAllSegAllSamplePhasedRawCPNList(sample.vec = sample.vec,
																ploidy.est.vec = ploidy.est.vec,
																purity.est.vec = purity.est.vec,
																min.cons.segs.df = min.cons.segs.df,
																all.seg.phasing.info.list = all.seg.phasing.info.list,
																asas.both.lr.baf.lists = asas.SNPs.w.both.lr.baf.snp.lists,
																gamma.val = GAMMA_VAL,
																min.num.hz.snp.phasing = MIN_NUM_HZ_SNP_PHASING)
	saveRDS(asas.ph.raw.cpn.lists, file = asas.ph.raw.cpn.lists.path)


	#Let's get the CCF
	print("Checking all minimum consistent segments are present in each sample...")
	seg.name.vec <- names(asas.ph.raw.cpn.lists[[1]])
	#Check to make sure that all segments are present in all samples... might not be necessary todo, check.
	for ( samp in sample.vec ) {

	      seg.name.vec <- intersect(seg.name.vec, names(asas.ph.raw.cpn.lists[[samp]]))

	}#end for ( samp in sample.vec ) { 

	if (length(seg.name.vec) == 0) {

		stop("No shared segments for some reason, go back and check")

	}


	print("Calculating CCFs using phased raw cpn...")
	all.seg.CCF.list <- list()
	all.seg.int.ob.xp.cpn.df.list <- list()
	all.seg.CCF.diff.from.exp.list <- list()
	all.seg.clonal.in.all.samps.list <- list()
	all.seg.greater.than.two.diff.list <- list()

	samps.int.cpn.raw.cpn.by.chr.list <- generateAllSamplesByChromIntCPNUnderlyingRawValuesLists(asas.ph.raw.cpn.lists = asas.ph.raw.cpn.lists,
																								 sample.vec = sample.vec,
																								 seg.name.vec = seg.name.vec,
																								 cpnX.vec=c("cpnA", "cpnB"),
																								 reduced.seg.flag = FALSE)
	saveRDS(samps.int.cpn.raw.cpn.by.chr.list, file = samps.int.cpn.raw.cpn.by.chr.list.path)


	for (cpnX in c("cpnA", "cpnB")) {
	  print(cpnX)
	  #Here is where we're going to have iniatalise on all the lists
	  all.seg.CCF.list[[cpnX]] <- list()
	  all.seg.int.ob.xp.cpn.df.list[[cpnX]] <- list()
	  all.seg.CCF.diff.from.exp.list[[cpnX]] <- list()
	  all.seg.clonal.in.all.samps.list[[cpnX]] <- list()
	  all.seg.greater.than.two.diff.list[[cpnX]] <- list()

	  for ( seg.name in seg.name.vec ) {

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

	    #does both A and B, todo, refactor so only does one
	    s.samps.int.ob.xp.cpn.df <- getSegObsExpCPNBySampleDF(asas.ph.raw.cpn.lists = asas.ph.raw.cpn.lists,
	                                                          sample.vec = sample.vec,
	                                                          seg.name = seg.name,
	                                                          cpnX.vec = cpnX,
	                                                          reduced.seg.flag = FALSE) 

	    #does only cpnX
	    seg.cpnX.diff.from.exp.list <- calculateWhetherSegCPNDiffFromIntegerExpectedAndAllMatchingSampCPN(asas.ph.raw.cpn.lists = asas.ph.raw.cpn.lists,
	                                                                                                      s.samps.int.ob.xp.cpn.df = s.samps.int.ob.xp.cpn.df,
	                                                                                                      seg.name = seg.name,
	                                                                                                      samps.int.cpn.raw.cpn.by.chr.list = samps.int.cpn.raw.cpn.by.chr.list,
	                                                                                                      cpnX = cpnX,
	                                                                                                      chrs.can.use.for.ccf.test = CHRS_TO_USE_VEC,
	                                                                                                      reduced.seg.flag = FALSE)
	    #does only cpnX
	    is.segment.clonal.in.all.samps.flag <- determineCCFClonalityOfSegment(s.samps.int.ob.xp.cpn.df = s.samps.int.ob.xp.cpn.df,
	                                                                          seg.name = seg.name,
	                                                                          asas.ph.raw.cpn.lists = asas.ph.raw.cpn.lists,
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
	                                                                                                two.value.diff.flag = is.there.cpnX.difference.greater.than.two.flag,
	                                                                                                min.CCF.to.consider.val = MIN_CCF_TO_CONSIDER)


	    seg.psc.cpnX.low  <- seg.psc.low.high.cpn.and.use.preset.ccf.flag.vec[1]
	    seg.psc.cpnX.high  <- seg.psc.low.high.cpn.and.use.preset.ccf.flag.vec[2]
	    seg.use.default.ccf <- as.logical(seg.psc.low.high.cpn.and.use.preset.ccf.flag.vec[3])

	    #update segment for each sample 
	    for ( samp  in sample.vec ) {

	      all.seg.CCF.list[[cpnX]][[seg.name]][[samp]] <- list()
	      #ss.ph.cpnX.vals.vec <- all.pw.cn.het.comp.lists[[samp]][[seg.name]][[paste0("ph_", cpnX, "_vec")]]
	      ss.ph.cpnX.vals.vec <- asas.ph.raw.cpn.lists[[samp]][[seg.name]][[paste0("ph_", cpnX, "_vec")]]

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
	      all.seg.CCF.list[[cpnX]][[seg.name]][[samp]][[paste0(cpnX, "_orig_ph_obs_cpnA_vec")]]   <- asas.ph.raw.cpn.lists[[samp]][[seg.name]][[paste0("ph_cpnA_vec")]]
	      all.seg.CCF.list[[cpnX]][[seg.name]][[samp]][[paste0(cpnX, "_orig_ph_obs_cpnB_vec")]]   <- asas.ph.raw.cpn.lists[[samp]][[seg.name]][[paste0("ph_cpnB_vec")]]
	      all.seg.CCF.list[[cpnX]][[seg.name]][[samp]][[paste0("clean_ph_obs_", cpnX, "_vec")]]   <- seg.samp.CCF.psc.cpn.df$clean_ph_obs_cpnX
	      all.seg.CCF.list[[cpnX]][[seg.name]][[samp]][[paste0("frac_", cpnX, "_low_vec")]]  <- seg.samp.CCF.psc.cpn.df$frac_cpnX_low 
	      all.seg.CCF.list[[cpnX]][[seg.name]][[samp]][[paste0("frac_", cpnX, "_high_vec")]]  <- seg.samp.CCF.psc.cpn.df$frac_cpnX_high

	    }#end  for ( samp  in sample.vec ) {

	    all.seg.int.ob.xp.cpn.df.list[[cpnX]][[seg.name]] <- s.samps.int.ob.xp.cpn.df
	    all.seg.clonal.in.all.samps.list[[cpnX]][[seg.name]] <- is.segment.clonal.in.all.samps.flag
	    all.seg.greater.than.two.diff.list[[cpnX]][[seg.name]] <- is.there.cpnX.difference.greater.than.two.flag
	    all.seg.CCF.diff.from.exp.list[[cpnX]][[seg.name]] <- seg.cpnX.diff.from.exp.list

	  }#for (seg.name %in% seg.name.vec 

	}#end for (cpnX in c("cpnA", "cpnB")) {

	#Save CCF data:
	saveRDS(all.seg.CCF.list, file = all.seg.CCF.list.path)
	saveRDS(all.seg.int.ob.xp.cpn.df.list, file = all.seg.int.ob.xp.cpn.df.list.path)
	saveRDS(all.seg.CCF.diff.from.exp.list, file = all.seg.CCF.diff.from.exp.list.path)
	saveRDS(all.seg.clonal.in.all.samps.list, file = all.seg.clonal.in.all.samps.list.path)
	saveRDS(all.seg.greater.than.two.diff.list, file = all.seg.greater.than.two.diff.list.path)


} else if (OVERWRITE_FLAG == FALSE) {

	all.seg.phasing.info.list <- readRDS(file = all.seg.phasing.info.list.path)
	patient.phased.BAF.df.list <- readRDS(file = patient.phased.BAF.df.list.path)
	all.seg.MSAI.by.BAF.ltnpd.present.list <- readRDS(file = all.seg.MSAI.by.BAF.ltnpd.present.list.path)
	all.samp.non.hz.logr.and.tot.cpn.df.list <- readRDS(file = all.samp.non.hz.logr.and.tot.cpn.df.list.path)
	asas.ph.raw.cpn.lists <- readRDS(file = asas.ph.raw.cpn.lists.path)
	samps.int.cpn.raw.cpn.by.chr.list <- readRDS(file = samps.int.cpn.raw.cpn.by.chr.list.path)
	all.seg.CCF.list <- readRDS(file = all.seg.CCF.list.path)
	all.seg.int.ob.xp.cpn.df.list <- readRDS(file = all.seg.int.ob.xp.cpn.df.list.path)
	all.seg.CCF.diff.from.exp.list <- readRDS(file = all.seg.CCF.diff.from.exp.list.path)
	all.seg.clonal.in.all.samps.list <- readRDS(file = all.seg.clonal.in.all.samps.list.path)
	all.seg.greater.than.two.diff.list <- readRDS(file = all.seg.greater.than.two.diff.list.path)

}
#Todo remove the dataframe initialisation from the function and modify it to work on a segment by segment basis.

#Let's initialise the first joint output:
print("Initialising minimum consistent segmentation (MCS) output dataframe...")
empty.full.joint.segs.df <- initialiseFullOutputDataframe(seg.name.vec = names(asas.ph.raw.cpn.lists[[1]]),
														ploidy.est.vec = ploidy.est.vec,
														purity.est.vec = purity.est.vec,
														make.row.names.seg.name.alone = FALSE) 


#This adds the raw cn for both segments that could be phased and those that couldn't
print("Adding new copy number calls to MCS output based on the phased raw copy number...")
# new.cpn.full.joint.segs.df <- addNewRawAndIntCPNCallsToFinalSegmentation(seg.name.vec = names(asas.ph.raw.cpn.lists[[1]]),
# 																		sample.vec = names(purity.est.vec),
# 																		all.samp.logr.df.list = all.samp.non.hz.logr.and.tot.cpn.df.list,
# 																		all.samp.ph.BAF.df.list = patient.phased.BAF.df.list,
# 																		all.seg.CCF.vals.list = all.seg.CCF.list,
# 																		asas.lr.baf.lists = asas.SNPs.w.both.lr.baf.snp.lists,
# 																		all.seg.df = empty.full.joint.segs.df,
# 																		cpnX.vec = c("cpnA", "cpnB")
# 																		)

source(paste0(current.dir, "/final_output_segmentation_file_generation_functions_v3.R")) 
source(paste0(current.dir, "/utility_MPhase_functions.R"))
new.cpn.full.joint.segs.df <- addNewRawAndIntCPNCallsToFinalSegmentationParallel(seg.name.vec = names(asas.ph.raw.cpn.lists[[1]]),
																		sample.vec = names(purity.est.vec),
																		all.samp.logr.df.list = all.samp.non.hz.logr.and.tot.cpn.df.list,
																		all.samp.ph.BAF.df.list = patient.phased.BAF.df.list,
																		all.seg.CCF.vals.list = all.seg.CCF.list,
																		asas.lr.baf.lists = asas.SNPs.w.both.lr.baf.snp.lists,
																		all.seg.df = empty.full.joint.segs.df,
																		cpnX.vec = c("cpnA", "cpnB")
																		)


addOriginalSegmentInfoToFinalSegmentationNew <- function(all.seg.df, all.orig.samp.segs.anno.df) {

	#Really assuming order is correct here!
	all.seg.df$overlap_w_n_orig_samps <- all.orig.samp.segs.anno.df$overlap_w_n_samps
	all.seg.df$num_segs_maps_to_in_orig_samp <-  all.orig.samp.segs.anno.df$num_segs_maps_to_in_samp
	all.seg.df$in_intersection_w_orig_segs <-  all.orig.samp.segs.anno.df$in_intersection
	all.seg.df$in_disjoint_set_w_orig_segs <-  all.orig.samp.segs.anno.df$in_disjoin_only
	all.seg.df$not_present_in_orig_segs <-  all.orig.samp.segs.anno.df$in_arm_setdiff
	all.seg.df$orig_samp_seg_name <- all.orig.samp.segs.anno.df$uniq_id
	all.seg.df$orig_nMinor <- all.orig.samp.segs.anno.df$nMinor
	all.seg.df$orig_nMajor <- all.orig.samp.segs.anno.df$nMajor
	all.seg.df$orig_cnTotal <- all.orig.samp.segs.anno.df$cnTotal
	all.seg.df$orig_cpn_call_vs_ploidy <- all.orig.samp.segs.anno.df$cpn_call
	all.seg.df$orig_seg_in_whole_arm_event <- all.orig.samp.segs.anno.df$arm_event_present
	all.seg.df$orig_seg_is_arm_or_focal_event <- all.orig.samp.segs.anno.df$arm_or_focal_event
	all.seg.df$orig_arm_seg_present_in <- all.orig.samp.segs.anno.df$within_chr_arm

	return(all.seg.df)

}


#Ok so now say we want to find out how many clonal and subclonal focal events there are, we need to:

#Make focal original id column:

#first identify clonal focal events

#Now that we know which ones are clonal focal in terms of position
#Subset by those that are also the same cpn vs ploidy all the way through
#Now we have our set of clonal ammp+gain / deep_loss+loss / LOH/ LOH combined with loss





#Now we also need  to repeat this exact same process for few and fewer samples
#and split between gain/loss/ and LOH




print("Mapping MCS segments back to original input and copy number states...")
full.joint.segs.w.old.cpn.vs.ploidy.df <- addOriginalSegmentInfoToFinalSegmentationNew(new.cpn.full.joint.segs.df, sorted.all.bp.all.samps.segs.anno.df)





print("Adding CCF information and subclonal cpn calls to MCS output...")
clonality.and.ccf.joint.segs.df <- addWithinSampleClonalityAndCCFToFinalSegmentation(seg.name.vec = names(asas.ph.raw.cpn.lists[[1]])
                                ,all.seg.CCF.vals.list = all.seg.CCF.list
                                ,all.seg.CCF.diff.from.exp.list = all.seg.CCF.diff.from.exp.list
                                ,all.seg.df = full.joint.segs.w.old.cpn.vs.ploidy.df#new.cpn.full.joint.segs.df
                                ,sample.vec = names(purity.est.vec)
                                ,clonal.pval = PVAL_THRESH_VS_MIN_CCF
                                ,diff.from.oth.p.val = PVAL_THRESH_VS_OTH_CHRS_SAME_CPN
                                ,min.CCF.to.consider.val = MIN_CCF_TO_CONSIDER
                                ,min.cpn.to.consider.val = MIN_CPN_TO_CONSIDER
                                ,min.num.other.chrs = MIN_NUM_CHRS
                                ,chrs.to.use.vec = CHRS_TO_USE_VEC
                                ,cpnX.vec = c("cpnA", "cpnB")
                                ,max.cpn.val = MAX_CPN_FOR_SUBCLONALITY
                                )

print("Adding total copy number estimates to MCS output for segments with only logr information...")
source(paste0(current.dir, "/final_output_segmentation_file_generation_functions_v3.R"))
source(paste0(current.dir, "/utility_MPhase_functions.R"))
full.joint.segs.df.w.cpn.vs.ploidy.df <- addLogrDefinedCPNEventsToOutputSegmentationParallel(sample.vec = names(purity.est.vec),
                                                        							seg.name.vec = names(asas.ph.raw.cpn.lists[[1]]),
                                                        							all.seg.df = clonality.and.ccf.joint.segs.df,
                                                        							ploidy.est.vec = ploidy.est.vec,
                                                        							purity.est.vec = purity.est.vec,
                                                        							all.samp.logr.df.list = all.samp.non.hz.logr.and.tot.cpn.df.list,
                                                        							all.samp.ph.BAF.df.list = patient.phased.BAF.df.list,     
                                                       								gain.threshold = GAIN_LOGR_THRESH,
                                                        							loss.threshold = LOSS_LOGR_THRESH,
                                                        							amp.threshold = AMP_LOGR_THRESH,  
                                                        							deep.loss.threshold = DEEP_LOSS_LOGR_THRESH,
                                                        							tot.cpn.vs.ploidy.pval.thresh = PVAL_THRESH_VS_LOGR_FOR_CPN_EVENT)

print(head(full.joint.segs.df.w.cpn.vs.ploidy.df[,"cpn_event_vs_ploidy"]))
print("Adding cpn event vs ploidy by logr measures to MCS output...")
source(paste0(current.dir, "/final_output_segmentation_file_generation_functions_v3.R"))
full.joint.segs.logr.vs.ploidy.clonality.df <- determineAcrossSampleClonalityOfLogrEventsVsPloidy(sample.vec = names(purity.est.vec),
                                                        							seg.name.vec = names(asas.ph.raw.cpn.lists[[1]]),
                                                        							all.seg.df = full.joint.segs.df.w.cpn.vs.ploidy.df,
                                                        							all.samp.logr.df.list = all.samp.non.hz.logr.and.tot.cpn.df.list,
                                                        							all.samp.ph.BAF.df.list = patient.phased.BAF.df.list,     
                                                       								gain.threshold = GAIN_LOGR_THRESH,
                                                        							loss.threshold = LOSS_LOGR_THRESH,
                                                        							amp.threshold = AMP_LOGR_THRESH,  
                                                        							deep.loss.threshold = DEEP_LOSS_LOGR_THRESH,
                                                        							tot.cpn.vs.ploidy.pval.thresh = PVAL_THRESH_VS_LOGR_FOR_CPN_EVENT)


print("Testing whether any gains should be reclassified as amplifications...")
# full.joint.check.amps.vs.ploidy.clonality.df <- retestSubclonalAmps(all.seg.df = full.joint.segs.logr.vs.ploidy.clonality.df,
# 																	all.samp.logr.df.list = all.samp.non.hz.logr.and.tot.cpn.df.list,
# 																	all.samp.ph.BAF.df.list = patient.phased.BAF.df.list,   
# 																	gain.threshold = GAIN_LOGR_THRESH,
# 																	amp.threshold = AMP_LOGR_THRESH,
# 																	p.val.for.converting.gain.to.amp = PVAL_THRESH_CONVERT_GAIN_TO_AMP)

print(head(full.joint.segs.logr.vs.ploidy.clonality.df[,c("across_samp_cl_cpn_vs_ploidy_present","across_samp_sc_cpn_vs_ploidy_present")]))


#These are kept in separate functions as the original segmentaiton was refactored so it could run when segments are rejoined at the end of the script.
print("Adding phasing information and MSAI detection by raw cpn calls to MCS output...")
full.joint.segs.w.phasing.df <- addPhasingInformationToOutputSegmentation(sample.vec = names(purity.est.vec),
                                        seg.name.vec = names(asas.ph.raw.cpn.lists[[1]]), 
                                        sample.to.phase.seg.list = all.seg.phasing.info.list,
                                        all.pw.cn.het.comp.lists = NULL,
                                        all.segs.df = full.joint.segs.logr.vs.ploidy.clonality.df)


# #Now we want to add in the allelic imbalance clonal and subclonal inrormation to the final data frame.
# print("Adding across sample clonal/subclonal summary to MCS output...")
# full.joint.segs.w.phasing.and.raw.df <- addOriginalCPNValuesAndCallsVsPloidyToOutputSegmentation(sample.vec = sample.vec,
#                                         seg.name.vec = names(asas.ph.raw.cpn.lists[[1]]), 
#                                         processed.orig.cpn.segs.df = processed.orig.cpn.segs.df,
#                                         all.segs.df = full.joint.segs.w.phasing.df)

print("Adding MSAI detected by LTNPD method to the MCS output...")
full.joint.segs.w.ltnpd.MSAI.df <- addMSAIDetectionByLTNPDToOutputSegmentation(sample.vec = sample.vec,
											seg.name.vec = names(asas.ph.raw.cpn.lists[[1]]),
											all.seg.MSAI.by.BAF.ltnpd.present.list = all.seg.MSAI.by.BAF.ltnpd.present.list,
											all.segs.df = full.joint.segs.w.phasing.df)#full.joint.segs.w.phasing.and.raw.df)


print(head(full.joint.segs.w.ltnpd.MSAI.df[( ! is.na(full.joint.segs.w.ltnpd.MSAI.df$any_across_samp_MSAI)) &
											full.joint.segs.w.ltnpd.MSAI.df$any_across_samp_MSAI == TRUE,
											c("MSAI", "any_across_samp_MSAI")]))


print("Defining any parallel gain or loss events detected via MSAI and addingto MCS output...")
full.joint.segs.w.parallel.cpn.events.df <- addParallelEventsByMSAIToOutputSegmentation(all.segs.df = full.joint.segs.w.ltnpd.MSAI.df)


print(head(full.joint.segs.w.parallel.cpn.events.df[(! is.na(full.joint.segs.w.parallel.cpn.events.df$MSAI)),
													c("any_across_samp_MSAI", "MSAI",
														"cpn_event_vs_ploidy",
														"is_parallel_gain_and_or_amp",
														"is_parallel_loss",
														"is_parallel_LOH",
														"is_parallel_loss_and_or_LOH")]))

print("Adding comparison between original cpn calls and MPhase copy number calling to MCS output...")
full.joint.segs.within.and.changed.cpn.info.df <- addWithinSampleChangeToCpnAndAllelicImbalanceCalls(full.joint.segs.w.parallel.cpn.events.df)


print("Adding across sample estimates of allelic imbalance and LOH clonality to MCS output...")
full.joint.segs.w.across.samp.AI.stats.df <- addAcrossSampleAlleleSpecificClonalityCalls(all.seg.df = full.joint.segs.within.and.changed.cpn.info.df,
																						sample.vec = sample.vec,
																						seg.name.vec = seg.name.vec)

print(head(full.joint.segs.w.across.samp.AI.stats.df[!is.na(full.joint.segs.w.across.samp.AI.stats.df$within_samp_sc_cpnA_present), c("any_within_samp_sc_cpnA",
"any_within_samp_sc_cpnB",
"across_samp_cl_AI",
"across_samp_sc_AI",
"across_samp_cl_LOH",
"across_samp_sc_LOH")]))


full.joint.segs.w.across.plotting.cols.df <- addPanCanColumnDefinitionToSegDF(full.joint.segs.w.across.samp.AI.stats.df)

print(full.joint.segs.w.across.plotting.cols.df[, c("pt_seg_cl_loss_LOH_no_MSAI", "pt_seg_sc_loss_LOH_incl_MSAI", "pt_seg_cl_gain_amp_no_MSAI", "pt_seg_sc_gain_amp_incl_MSAI")])

full.joint.segs.df <- full.joint.segs.w.across.plotting.cols.df
#Remove segments without any datapoints.
#Todo, merge these int
#print("Removing any segments that do not have underlying datapoints, this is an artefact of the minimum consistent segmentation performed by MPhase...")
#full.joint.segs.df <- full.joint.segs.df[full.joint.segs.df$num_un_ph_hz_sites != 0 | full.joint.segs.df$num_ph_hz_sites != 0  | full.joint.segs.df$num_non_hz_sites != 0,]
#print("Removing segments that do not map back to any original segments")
#full.joint.segs.df <- full.joint.segs.df[full.joint.segs.df$orig_seg_name != "NA:NA:NA",]



full.joint.segs.df.path <- paste0(PLOT_OUTPUT_DIR, patient, "_full_joint_segs_df.RDS")
saveRDS(full.joint.segs.df, file=full.joint.segs.df.path)



shouldSegsWithCPNXValuesBeJoined <- function(s.curr.df, s.next.df, pval.to.separate) {

    cpnX.vec.s.curr.vs.s.next.test <- wilcoxTestWithErrorHandling(s.curr.df$cpnX, s.next.df$cpnX)
    frac.vec.s.curr.vs.s.next.test <- wilcoxTestWithErrorHandling(s.curr.df$cpnX_frac, s.next.df$cpnX_frac)

    if (cpnX.vec.s.curr.vs.s.next.test$p.val > pval.to.separate & frac.vec.s.curr.vs.s.next.test$p.val > pval.to.separate) {

        return(TRUE)

    } else {

        return(FALSE)
    }

}#end shouldSegsBeJoined

# all.orig.samp.segs.anno.df


# #So now we want to make a 
# #Get the chromosome positions:

# #So now for all the segments that are chosen we want to return a dataframe with the following columns:

# #So now we want to add what arm every segment is part of
# #Add whether that arm is a signficant

# #    Arm X..Genes Amp.frequency Amp.z.score Amp.q.value Del.Frequency Del.z.score Del.q.value
# # 1   1p     1300          0.31      3.0200    6.47e-03          0.23      -0.723
# #sig.hist.arm.amp.df <- hist.arm.amp.df[ hist.arm.amp.df[,"Amp.q.value"] <= 0.05 , ]


# sig.hist.arm.del.df <- getGISTICSigArms(hist.type, "DEL")
# getGISTICSigFocals()


# getGISTICSigFocals()

# # sorted.all.bp.all.samps.segs.anno.df <- makeGRangesFromDataFrame(sorted.all.bp.all.samps.segs.anno.df,
# # 	                         keep.extra.columns=FALSE,
# # 	                         ignore.strand=FALSE,
# # 	                         seqinfo=NULL,
# # 	                         seqnames.field=c("seqnames", "seqname",
# # 	                                          "chromosome", "chrom",
# # 	                                          "chr", "chromosome_name",
# # 	                                          "seqid"),
# # 	                         start.field=c("Start_Position","start", "startpos"),
# # 	                         end.field=c("End_Position", "end", "stop", "endpos"),
# # 	                         strand.field="strand",
# 	                         starts.in.df.are.0based=FALSE)

# cpnX.vec <- c("cpnB", "cpnA")
# p.or.q.vec <- c("p", "q")
# for ( chr in CHRS_TO_USE_VEC ) {
# 	#chr <- 21
# 	chr.gr.list <- getChromosomeArmsGRangeList(chr, chrom.info.df)

# 		for (samp in sample.vec) {

# 			for (cpnX in cpnX.vec) {

# 				for (i in 1:length(chr.gr.list)) {

# 				#Make a cnpA set of Granges that match this arm but have all info for both alleles
# 				#Make a cpnB set of Granges that match this arm but have all info for both alleles
# 				getOverlapGRanges

# 				}#end for (i in 1:length(chr.gr.list)) {


# 			}#for (cpnX in cpnX.vec) {


# 		}#end for (samp in sample.vec) {

# }# for ( chr in CHRS_TO_USE_VEC ) {

#Pseudocode the joining:


#Need to add arms to all segments:

# #All possible useful input:
# patient.phased.BAF.df.list
# > str(patient.phased.BAF.df.list[[1]], max=1)
# 'data.frame':	33192 obs. of  11 variables:
#  $ chrom                       : chr  "1" "1" "1" "1" ...
#  $ pos                         : int  871020 879317 880390 881627 892460 894573 896064 897730 902108 906272 ...
#  $ BAF                         : num  0.257 0.248 0.2 0.866 0.143 ...
#  $ BAFsegmented                : num  0.23 0.23 0.23 0.23 0.23 ...
#  $ logR                        : num  -0.271 -0.36 -0.385 -0.208 -0.571 ...
#  $ logRsegmented               : num  -0.335 -0.335 -0.335 -0.335 -0.335 ...
#  $ samp                        : chr  "A_LTX121_SU_T1-R1" "A_LTX121_SU_T1-R1" "A_LTX121_SU_T1-R1" "A_LTX121_SU_T1-R1" ...
#  $ sample_phased_by_this_sample: chr  "A_LTX121_SU_T1-R1" "A_LTX121_SU_T1-R1" "A_LTX121_SU_T1-R1" "A_LTX121_SU_T1-R1" ...
#  $ maj_is_cpnA_logic_vec       : num  0 0 0 1 0 1 0 0 0 1 ...
#  $ maj_is_cpnB_logic_vec       : num  1 1 1 0 1 0 1 1 1 0 ...
#  $ ph_median_segmented_BAF     : num  0.207 0.207 0.207 0.207 0.207 ...
# all.seg.MSAI.by.BAF.ltnpd.present.list
# all.samp.non.hz.logr.and.tot.cpn.df.list
# > str(all.samp.non.hz.logr.and.tot.cpn.df.list[[1]], max=1)
# 'data.frame':	921620 obs. of  7 variables:
#  $ chrom            : chr  "1" "1" "1" "1" ...
#  $ pos              : int  762273 762496 762497 777237 777318 777335 777361 777456 783071 783147 ...
#  $ BAF              : num  1 NA NA NA NA ...
#  $ BAFsegmented     : num  NA NA NA NA NA NA NA NA NA NA ...
#  $ logR             : num  -0.0542 -0.2103 -0.2103 -0.4974 -0.3978 ...
#  $ logRsegmented    : num  -0.335 -0.335 -0.335 -0.335 -0.335 ...
#  $ tot_raw_un_ph_cpn: num  2.98 2.47 2.47 1.67 1.93 ...
# asas.ph.raw.cpn.lists
# > str(asas.ph.raw.cpn.lists[[1]], max=1)
# List of 418
#  $ 1:1:762272            :List of 9
#  $ 1:762273:26136282     :List of 9
# all.seg.CCF.list
# > str(all.seg.CCF.list[[1]] ,max=1)
# List of 418
#  $ 1:1:762272            : logi NA
#  $ 1:762273:26136282     :List of 2
#  $ 1:26136283:26136297   : logi NA
#  > str(all.seg.CCF.list[[1]][2])
#  List of 1
#   $ 1:762273:26136282:List of 2
#    ..$ A_LTX121_SU_T1-R1:List of 9
#    .. ..$ poss_sc_cpnA_low         : num [1:755] 2 2 2 2 2 2 2 2 2 2 ...
#    .. ..$ poss_sc_cpnA_high        : num [1:755] 3 3 3 3 3 3 3 3 3 3 ...
#    .. ..$ mean_frac_cpnA_low       : num [1:755] 0.83 0.83 0.83 0.83 0.83 0.83 0.83 0.83 0.83 0.83 ...
#    .. ..$ mean_frac_cpnA_high      : num [1:755] 0.17 0.17 0.17 0.17 0.17 0.17 0.17 0.17 0.17 0.17 ...
#    .. ..$ cpnA_orig_ph_obs_cpnA_vec: num [1:755] 2.19 2.03 2.17 2.88 1.98 ...
#    .. ..$ cpnA_orig_ph_obs_cpnB_vec: num [1:755] 0.10186 0.000399 -0.205837 -0.400462 -0.501425 ...
#    .. ..$ clean_ph_obs_cpnA_vec    : num [1:755] 2.19 2.03 2.17 2.88 1.98 ...
#    .. ..$ frac_cpnA_low_vec        : num [1:755] 0.81 0.97 0.83 0.12 1 0.69 0.68 0.88 0.86 1 ...
#    .. ..$ frac_cpnA_high_vec       : num [1:755] 0.19 0.03 0.17 0.88 0 ...
#    ..$ A_LTX121_SU_T1-R2:List of 9
#    .. ..$ poss_sc_cpnA_low         : num [1:755] 2 2 2 2 2 2 2 2 2 2 ...
#    .. ..$ poss_sc_cpnA_high        : num [1:755] 3 3 3 3 3 3 3 3 3 3 ...
#    .. ..$ mean_frac_cpnA_low       : num [1:755] 0.74 0.74 0.74 0.74 0.74 0.74 0.74 0.74 0.74 0.74 ...
#    .. ..$ mean_frac_cpnA_high      : num [1:755] 0.26 0.26 0.26 0.26 0.26 0.26 0.26 0.26 0.26 0.26 ...
#    .. ..$ cpnA_orig_ph_obs_cpnA_vec: num [1:755] 1.56 2.26 2.61 2.16 1.74 ...
#    .. ..$ cpnA_orig_ph_obs_cpnB_vec: num [1:755] 0.32 -0.258 -0.22 -0.531 0.212 ...
#    .. ..$ clean_ph_obs_cpnA_vec    : num [1:755] 1.56 2.26 2.61 2.16 1.74 ...
#    .. ..$ frac_cpnA_low_vec        : num [1:755] 1 0.74 0.39 0.84 1 1 0.73 0.38 0.44 1 ...
#    .. ..$ frac_cpnA_high_vec       : num [1:755] 0 0.26 0.61 0.16 0 0 0.27 0.62 0.56 0 ...
# > str(all.seg.CCF.list[[2]][2])
# List of 1
# $ 1:762273:26136282:List of 2
#  ..$ A_LTX121_SU_T1-R1:List of 9
#  .. ..$ poss_sc_cpnB_low         : num [1:755] 1 1 1 1 1 1 1 1 1 1 ...
#  .. ..$ poss_sc_cpnB_high        : num [1:755] 1 1 1 1 1 1 1 1 1 1 ...
#  .. ..$ mean_frac_cpnB_low       : num [1:755] 0 0 0 0 0 0 0 0 0 0 ...
#  .. ..$ mean_frac_cpnB_high      : num [1:755] 1 1 1 1 1 1 1 1 1 1 ...
#  .. ..$ cpnB_orig_ph_obs_cpnA_vec: num [1:755] 2.19 2.03 2.17 2.88 1.98 ...
#  .. ..$ cpnB_orig_ph_obs_cpnB_vec: num [1:755] 0.10186 0.000399 -0.205837 -0.400462 -0.501425 ...
#  .. ..$ clean_ph_obs_cpnB_vec    : num [1:755] 0.10186 0.000399 0 0 0 ...
#  .. ..$ frac_cpnB_low_vec        : num [1:755] 0 0 0 0 0 0 0 0 0 0 ...
#  .. ..$ frac_cpnB_high_vec       : num [1:755] 1 1 1 1 1 1 1 1 1 1 ...
#  ..$ A_LTX121_SU_T1-R2:List of 9
#  .. ..$ poss_sc_cpnB_low         : num [1:755] 1 1 1 1 1 1 1 1 1 1 ...
#  .. ..$ poss_sc_cpnB_high        : num [1:755] 1 1 1 1 1 1 1 1 1 1 ...
#  .. ..$ mean_frac_cpnB_low       : num [1:755] 0 0 0 0 0 0 0 0 0 0 ...
#  .. ..$ mean_frac_cpnB_high      : num [1:755] 1 1 1 1 1 1 1 1 1 1 ...
#  .. ..$ cpnB_orig_ph_obs_cpnA_vec: num [1:755] 1.56 2.26 2.61 2.16 1.74 ...
#  .. ..$ cpnB_orig_ph_obs_cpnB_vec: num [1:755] 0.32 -0.258 -0.22 -0.531 0.212 ...
#  .. ..$ clean_ph_obs_cpnB_vec    : num [1:755] 0.32 0 0 0 0.212 ...
#  .. ..$ frac_cpnB_low_vec        : num [1:755] 0 0 0 0 0 0 0 0 0 0 ...
#  .. ..$ frac_cpnB_high_vec       : num [1:755] 1 1 1 1 1 1 1 1 1 1 ...
# asas.SNPs.w.both.lr.baf.snp.lists
# > str(asas.SNPs.w.both.lr.baf.snp.lists[[1]][2], max=1)
# List of 1
#  $ 1:762273:26136282:'data.frame':	755 obs. of  6 variables:
# > str(asas.SNPs.w.both.lr.baf.snp.lists[[1]][2])
# List of 1
#  $ 1:762273:26136282:'data.frame':	755 obs. of  6 variables:
#   ..$ chrom        : chr [1:755] "1" "1" "1" "1" ...
#   ..$ pos          : int [1:755] 871020 879317 880390 881627 892460 894573 896064 897730 902108 906272 ...
#   ..$ BAF          : num [1:755] 0.257 0.248 0.2 0.866 0.143 ...
#   ..$ BAFsegmented : num [1:755] 0.23 0.23 0.23 0.23 0.23 ...
#   ..$ logR         : num [1:755] -0.271 -0.36 -0.385 -0.208 -0.571 ...
#   ..$ logRsegmented: num [1:755] -0.335 -0.335 -0.335 -0.335 -0.335 ...
# asas.SNPS.lr.and.or.baf.lists
# > str(asas.SNPS.lr.and.or.baf.lists[[1]][2])
# List of 1
#  $ 1:762273:26136282:'data.frame':	17836 obs. of  6 variables:
#   ..$ chrom        : chr [1:17836] "1" "1" "1" "1" ...
#   ..$ pos          : int [1:17836] 762273 762496 762497 777237 777318 777335 777361 777456 783071 783147 ...
#   ..$ BAF          : num [1:17836] 1 NA NA NA NA ...
#   ..$ BAFsegmented : num [1:17836] NA NA NA NA NA NA NA NA NA NA ...
#   ..$ logR         : num [1:17836] -0.0542 -0.2103 -0.2103 -0.4974 -0.3978 ...
#   ..$ logRsegmented: num [1:17836] -0.335 -0.335 -0.335 -0.335 -0.335 ...


#So we want:
# chr
# pos
# hz
# non_hz
# ph_cpn
# tot_pn
# BAF
# logR


joinSegmentsWithinSampleChromosomeArm <- function(arm.segs, arm.samp.ph.tot.cpn.df, cpnX) {


	#Original input segments:

	#Make a cpnA copy of the segments:
	
	#Make a cpnB copy of the segments:

	#Make a cpnA results dataframe and list:

	#Make a cpnA results dataframe and list:

	#For each cpnA or cpnB

	#Sort segments by size

	#Choose the smallest segment

	#Determine what kind of segment this smallest segment is

	#Get the ante segment, if there is one

	#Determine what kind of segment this preceding segment is

	#Get the  pst segment, if there is one.

	#Determine what kind of segment post segment is.

	#So what are the possible cases:

	#So each seg should have a vector that determines its identity and course of action.

	#num_ph_hz,
	#num_unph_hz,
	#num_non_hz,
	#non_cons_logr,
	#has_ph_cpn,

	#If they both have phased copy number -> test for difference using wilcox test
	#									  -> test for MAD between means and see if lower than thesthold
	#Join if nonsignificant/below threshold.
	#else leave separate

	#If one has phased copy number and the other doesn't? -> Test total copy number vs total copy number using wilcox test
	#													  -> Test for mad between means of distribution
	#Join if nonsignificant/below threshold

	#If neither has phased copy number? -> Test total copy number 
	#									-> 

	#If nether has any data behind it?



}



#We are at the level of an arm:

#Start with smallest segment:

#Out of the preceding and subsequent segment attempt the following with the smallest first:

#There are a few cases:

#The segment we are trying to join is:
	#has allele specific cpn -> join if difference is non significant by significance threshold and mean cpn threshold
	#only has total cpn -> test vs total cpn either side, if non-significant join to that segment and extend allele specific calls
	#					-> if different, look at which segment was closest in total copy number to it:
	#						-> extend the minor allele of the segment that is closest to it total cpn, make up the difference in using the major (default cpnA), mark that this was done
	#has no information associated, or below a threahold (3 values?)
	#					-> join to the largest segment as this is likely to be most accurate.




#What do we need to make








#Ok how to join the segments:
#Order of choices:
#First allele specific
#Then total copy number
#Then just joinin
#Do it by arm


#we'te going to need the lists that contain the information about the segments


#Work within a chromosome arm:

#Check whether any changes made for this segment (is this the first time we've seen it, if we see it twice then we know we've finished with the arm)
#Get smallest segment

#Segment with logr and baf? how can we join these?
#Answer: they are completely separate so should be consudered differently:



#What to we do if the segments have logr and no BAF?
#First check whether total copy number differs from neighbors (smallest first) -> join if doesn't and maintain allele specific
#If not then we don't have allele specific, try and keep lower allele count from cloest segment total
#Then increase the other to match the total cn and mark that this was performed.

#Segment with no logr or baf? Joint to largest segment next to them.


#Add new segment and sort the set




png(file = paste(PLOT_OUTPUT_DIR, patient, "_BAF_phased", ".png", sep=""), width = 1200, height = 1400)

print("Plotting phased BAF across the genome...")
plotAllSampleBAFsAcrossGenome(sample.vec, patient.phased.BAF.df.list, cpnA.allele.color = '#de2d2699', cpnB.allele.color= '#3182bd99')

dev.off()

####Plot all genome wide plots
print("Plotting CCFs for cpnA and cpnB across the genome using cpn definition of subclonality...")
png(file = paste(PLOT_OUTPUT_DIR, patient, "_CCF_w_cpn_thresh_cpnA", ".png", sep = ""), width = 1200, height = 1400)

plotAllSampleCCFsAcrossGenome(asas.ph.raw.cpn.lists = asas.ph.raw.cpn.lists,
                              full.joint.segs.df = full.joint.segs.w.phasing.df,
                              sample.vec = sample.vec,
                              all.seg.CCF.vals.list = all.seg.CCF.list,
                              cpnX = "cpnA",
                              use.cpn.or.ccf.thresh = "CPN")

dev.off()

png(file = paste(PLOT_OUTPUT_DIR, patient, "_CCF_w_cpn_thresh_cpnB", ".png", sep = ""), width = 1200, height = 1400)

plotAllSampleCCFsAcrossGenome(asas.ph.raw.cpn.lists = asas.ph.raw.cpn.lists,
                              full.joint.segs.df = full.joint.segs.w.phasing.df,
                              sample.vec = sample.vec,
                              all.seg.CCF.vals.list = all.seg.CCF.list,
                              cpnX = "cpnB",
                              use.cpn.or.ccf.thresh = "CPN")
dev.off()



STOP


####Plot all genome wide plots
print("Plotting CCFs for cpnA and cpnB across the genome using CCF definition of subclonality...")
png(file = paste(PLOT_OUTPUT_DIR, patient, "_CCF_w_CCF_thresh_cpnA", ".png", sep = ""), width = 1200, height = 1400)

plotAllSampleCCFsAcrossGenome(asas.ph.raw.cpn.lists = asas.ph.raw.cpn.lists,
                              full.joint.segs.df = full.joint.segs.w.phasing.df,
                              sample.vec = sample.vec,
                              all.seg.CCF.vals.list = all.seg.CCF.list,
                              cpnX = "cpnA",
                              use.cpn.or.ccf.thresh = "CCF")

dev.off()

png(file = paste(PLOT_OUTPUT_DIR, patient, "_CCF_w_CCF_thresh_cpnB", ".png", sep = ""), width = 1200, height = 1400)

plotAllSampleCCFsAcrossGenome(asas.ph.raw.cpn.lists = asas.ph.raw.cpn.lists,
                              full.joint.segs.df = full.joint.segs.w.phasing.df,
                              sample.vec = sample.vec,
                              all.seg.CCF.vals.list = all.seg.CCF.list,
                              cpnX = "cpnB",
                              use.cpn.or.ccf.thresh = "CCF")
dev.off()



#Plot the region comparison across genome plots, will need their own directories:
old.vs.new.raw.cpn.dir <- paste(PLOT_OUTPUT_DIR, "old_raw_cpn_vs_MPhase_raw_cpn/", sep="")
if( ! file.exists(old.vs.new.raw.cpn.dir)) {

	dir.create(old.vs.new.raw.cpn.dir)

}#end if( ! file.exists(by.chrom.CCF.plot.dir)) {

print("Plotting comparison of old vs new MPhase raw cpn calls..")
plotOldNewComparisonOfAllCPN(sample.vec, full.joint.segs.df, plot.output.dir = old.vs.new.raw.cpn.dir)

pairwise.raw.cpn.dir <- paste(PLOT_OUTPUT_DIR, "pairwise_raw_cpn_comparisons/", sep = "")
if( ! file.exists(pairwise.raw.cpn.dir)) {

	dir.create(pairwise.raw.cpn.dir)

}#end if( ! file.exists(by.chrom.CCF.plot.dir)) {

print("Plotting pairwise sample comparison between MPhase raw cpn calls across the whole genome...")
plotAllSamplePairwiseComparisonSegmentCPN(sample.vec, full.joint.segs.df, plot.output.dir = pairwise.raw.cpn.dir)



#Plot by chromosome plots:

uniq.chrs <- unique(full.joint.segs.w.phasing.df$chr)
by.chrom.CCF.plot.dir <- paste(PLOT_OUTPUT_DIR, "CCF_by_chrom/", sep = "")
if( ! file.exists(by.chrom.CCF.plot.dir)) {

	dir.create(by.chrom.CCF.plot.dir)

}#end if( ! file.exists(by.chrom.CCF.plot.dir)) {
print("Plotting by chromosome CCF estimates for cpnA and cpnB using CCF subclonality thresholds...")
png(file = paste(by.chrom.CCF.plot.dir, "/", patient, "_CCF_w_CCF_thresh_cpnA_", "%02d", ".png", sep = ""), width = 1200, height = 1400)
for ( chrom in uniq.chrs ) {

	plotAllSampleCCFsAcrossGenome(asas.ph.raw.cpn.lists = asas.ph.raw.cpn.lists,
                              full.joint.segs.df = full.joint.segs.w.phasing.df[full.joint.segs.w.phasing.df$chr == chrom,] ,
                              sample.vec = sample.vec,
                              all.seg.CCF.vals.list = all.seg.CCF.list,
                              cpnX = "cpnA",
                              use.cpn.or.ccf.thresh = "CCF")


}
dev.off()

png(file = paste(by.chrom.CCF.plot.dir, "/", patient, "_CCF_w_CCF_thresh_cpnB_", "%02d", ".png", sep = ""), width = 1200, height = 1400)
for ( chrom in uniq.chrs ) {

	plotAllSampleCCFsAcrossGenome(asas.ph.raw.cpn.lists = asas.ph.raw.cpn.lists,
                              full.joint.segs.df = full.joint.segs.w.phasing.df[full.joint.segs.w.phasing.df$chr == chrom,] ,
                              sample.vec = sample.vec,
                              all.seg.CCF.vals.list = all.seg.CCF.list,
                              cpnX = "cpnB",
                              use.cpn.or.ccf.thresh = "CCF")


}
dev.off()

print("Plotting by chromosome CCF estimates for cpnA and cpnB using cpn subclonality thresholds...")
png(file = paste(by.chrom.CCF.plot.dir, "/", patient, "_CCF_w_cpn_thresh_cpnA_","%02d", ".png", sep = ""), width = 1200, height = 1400)
for ( chrom in uniq.chrs ) {

	plotAllSampleCCFsAcrossGenome(asas.ph.raw.cpn.lists = asas.ph.raw.cpn.lists,
                              full.joint.segs.df = full.joint.segs.w.phasing.df[full.joint.segs.w.phasing.df$chr == chrom,] ,
                              sample.vec = sample.vec,
                              all.seg.CCF.vals.list = all.seg.CCF.list,
                              cpnX = "cpnA",
                              use.cpn.or.ccf.thresh = "CPN")


}
dev.off()



png(file = paste(by.chrom.CCF.plot.dir, "/", patient, "_CCF_w_cpn_thresh_cpnB_", "%02d", ".png", sep = ""), width = 1200, height = 1400)
for ( chrom in uniq.chrs ) {

	plotAllSampleCCFsAcrossGenome(asas.ph.raw.cpn.lists = asas.ph.raw.cpn.lists,
                              full.joint.segs.df = full.joint.segs.w.phasing.df[full.joint.segs.w.phasing.df$chr == chrom,] ,
                              sample.vec = sample.vec,
                              all.seg.CCF.vals.list = all.seg.CCF.list,
                              cpnX = "cpnB",
                              use.cpn.or.ccf.thresh = "CPN")


}
dev.off()

by.chrom.BAF.plot.dir <- paste(PLOT_OUTPUT_DIR, "BAF_by_chrom/", sep="")
if( ! file.exists(by.chrom.BAF.plot.dir)) {

	dir.create(by.chrom.BAF.plot.dir)

}#end if( ! file.exists(by.chrom.CCF.plot.dir)) {

#Now let's plot baf
print("Plotting by chromosome phased BAF...")
png(file = paste(by.chrom.BAF.plot.dir, "/", patient, "_BAF_phased_", "%02d", ".png", sep=""), width = 1200, height = 1400)
for ( chrom in uniq.chrs ) {

	chr.pt.ph.BAF.df.list <- patient.phased.BAF.df.list

	for ( samp in sample.vec ) {

		chr.pt.ph.BAF.df.list[[samp]] <- chr.pt.ph.BAF.df.list[[samp]][chr.pt.ph.BAF.df.list[[samp]]$chrom == chrom,]

	}#end for ( samp in sample.vec ) {

	plotAllSampleBAFsAcrossGenome(sample.vec, chr.pt.ph.BAF.df.list, cpnA.allele.color = '#de2d2699', cpnB.allele.color= '#3182bd99')


}#end for ( chrom in uniq.chrs ) {
dev.off()

#testFinalOutputDataframe(full.joint.segs.df)

