# MultiRegionCopyNum

Aneuploidy and chromosomal instability are a hallmark of cancer and result in somatic copy number alteration (SCNAs). Many bioinformatics tools have been written to analyse single tumour regions to estimate the aneuploid state of the cancer genome and the allele-specific copy number state at each locus in the genome. The scripts in this repository constitute a bioinformatics tool described as MPhase that is designed to estimate the allele-specific copy number of multiple samples from the same tumour using phasing data obtained from all samples of a patient's disease.

In order to do this MPhase takes as input for each sample of the same patient's disease investigated:

estimate of purity 
estimate of ploidy
allele specific copy number segmentation 
A set of homozygous and heterozygous SNP positions 
Log ratio values for each position
B-allele frequency values each position

This current version is designed to run from ASCAT output of these values but as long as the data is packaged.

These values are packaged into the following input files:

A tumour sample BAF and LogR file that is an output from ASCAT.
  chrom    pos      BAF            LogR
1     1 762273 1.013648            -0.133
2     1 762496       NA            -0.244

The snp_file_path parameter will be a path to a file containg the paths to all snp files of the above format for each sample:

patient_id	sample	samp_snp_path
patient	sample	path_to_snp_file


A combined segmentation file for all samples from the original ASCAT segmentation containing purity and ploidy estimate:

sample	chr	startpos	endpos	n.het	cnTotal	nMajor	nMinor	Ploidy	ACF	nAraw	nBraw
sample	1	762273	16803458	11708	3	2	1	3.1614395468455	0.14	1.98176736346188	0.616233906001798


The Mphase tool itself is called with the following parameters:

patient = a string representing the patient ID
diptest_phasing_flag = should be set to FALSE as this allows a more lax phasing still in development
output_directory_path =  path to the directory that the output shoud be placed in
data_format_to_use = this should be set to "MEDICC"
snp_file_path = path to the snp file that contains the paths of all snp files to be used in the analysis
seg_file_path =  path to the combined segmentation files
overwrite_flag = "TRUE" or "FALSE" to indicate whether prexisting files in the output 


An example run would be:

Rscript path_to_mphase_dir/Mphase.R patient_001 FALSE output_directory_path MEDICC snp_file_path seg_file_path TRUE











