
#So let's get the
source("/camp/lab/swantonc/working/watkint/scripts/TRACERx_421/SCNA_utility_functions.R")
full.pt.vec <- getReleasePatientVec(shorten=FALSE, release.path=TX_421_DIR, pattern="LTX")
all.pts.mph.df <- checkMPhaseOutput(full.patient.vec = full.pt.vec, release.path = TX_421_DIR)

#full.pt.vec <- "A_LTX210"

for (full.patient in full.pt.vec) {

	short.pt <- substring(full.patient, 3, 8)

	patient <- short.pt
	dip.test.flag <- FALSE
	output.dir <- paste0("/camp/project/proj-tracerx-lung/tctProjects/watkint/tx421_copy_number_qc/xphase_runs/", full.patient, "/")
	system(paste0("mkdir ",output.dir))

	medicc.output.dir <- paste0("/camp/project/proj-tracerx-lung/tctProjects/watkint/tx421_copy_number_qc/xphase_runs/", full.patient, "/MEDICC_output/")
	input.format <- "ASCAT"
	snp.file.path <- paste0("/camp/project/proj-tracerx-lung/tracerx/release_tx421_20190411/",full.patient,"/exome/MPhase/",short.pt,"_TRACERx_input_paths_df.txt")
	seg.file.path <- paste0("/camp/project/proj-tracerx-lung/tracerx/release_tx421_20190411/",full.patient,"/exome/MPhase/",short.pt,"_joint_segs_all_ascat_info.txt")
	overwrite.flag <- TRUE 
	gamma.val <- 1
	num.cores <- 1
	seg.type <- "intersection"
	chosen.genome <- "hg19"
	seed.val <- 2134234

	script.path <- "/camp/lab/swantonc/working/watkint/scripts/TRACERx_421/multiregioncopynum/xphase.R"
	outputfile <- paste(output.dir, full.patient, "_stdout.txt", sep = "") 
	errorfile <- paste(output.dir, full.patient, "_stderr.txt", sep = "")

	#So now we also want to create an xphase command:
	xphase.command <- paste0("Rscript ",script.path ,
				   " --patient=",patient,
				   " --dip-test=",dip.test.flag,
				   " --output-dir=",output.dir,
				   " --medicc-output-dir=",medicc.output.dir,
				   " --input-format=",input.format,
				   " --snp-data=",snp.file.path,
				   " --seg-data=",seg.file.path,
				   " --overwrite=",overwrite.flag,
				   " --gamma=",gamma.val,
				   " --cores=",num.cores,
				   " --segment-type=",seg.type,
				   " --genome=",chosen.genome,
				   " --seed=",seed.val) 

	pipe.to.sbatch <- paste("sbatch --wrap=","\"", xphase.command,"\"", " -o ", outputfile," -e ",errorfile, " -t 1-0 -c 1 --mem=12G ", sep="")
	system(pipe.to.sbatch)

#STOP
}









# option_list <- list(
#     make_option(c("--patient"), default="CB1002",
#                 help="Patient id"),
#     make_option(c("--dip-test"), default=FALSE,
#                 help="Use dip test (default=TRUE)"),
#     make_option(c("--output-dir"), default="../output/CB1002/mphase",
#                 help="Output directory"),
#     make_option(c("--medicc-output-dir"), default="../output/CB1002/mphase",
#                 help="Output directory"),
#     make_option(c("--input-format"), default="ASCAT",
#                 help="Input data format (default=MEDICC)"),
#     make_option(c("--snp-data"), default="../input/nb-schulte//CB1002_input_paths_df.txt",
#                 help="SNP data"),
#     make_option(c("--seg-data"), default="../input/nb-schulte//CB1002.segmentation",
#                 help="Segmentation data output from ASCAT"),
#     make_option(c("--overwrite"), default=TRUE,
#                 help="Overwrite output (default=TRUE)"),
#     make_option(c("--gamma"), default=1.0,
#                 help="Gamma (for seq=1, for arrays~=0.5 (default=1.0)"),
#     make_option(c("-j", "--cores"), default=1,
#                 help="CPU cores (default=1)"),
#     make_option(c("--segment-type"), default="intersection",
#                 help="Segment type, intersection or whole_arm (default=1)"),
#     make_option(c("--genome"), default="hg38",
#                 help="Human genome version (default=hg19, also supports hg38)"),
#     make_option(c("--seed"), default=2134234,
#                 help="Random seed")
#     )
