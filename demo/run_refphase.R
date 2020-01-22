library(refphase)
library(optparse)

option_list <- list(
  make_option(c("--patient"),
    default = "CB1002",
    help = "Patient id"
  ),
  make_option(c("--dip-test"),
    default = FALSE, action = "store_true",
    help = "Use dip test [default = %default]"
  ),
  make_option(c("--output-dir"),
    default = "../output/CB1002/mphase",
    help = "Output directory"
  ),
  make_option(c("--medicc-output-dir"),
    default = NA, type = "character",
    help = "Output directory for MEDICC input files"
  ),
  make_option(c("--input-format"),
    default = "ASCAT",
    help = "Input data format [default = %default]"
  ),
  make_option(c("--snp-data"),
    default = "../input/nb-schulte//CB1002_input_paths_df.txt",
    help = "SNP data"
  ),
  make_option(c("--snp-data-prefix"),
    default = "",
    help = "Directory prefix for the paths to the SNP data in the file above"
  ),
  make_option(c("--seg-data"),
    default = "../input/nb-schulte//CB1002.segmentation",
    help = "Segmentation data output from ASCAT"
  ),
  make_option(c("--gamma"),
    default = 1.0, type = "numeric",
    help = "Gamma (for seq=1, for arrays~=0.5 [default = %default]"
  ),
  make_option(c("-j", "--cores"),
    default = 1, type = "integer",
    help = "CPU cores [default = %default]"
  ),
  make_option(c("--segment-type"),
    default = "intersection",
    help = "Segment type, intersection or whole_arm [default = %default]"
  ),
  make_option(c("--genome"),
    default = "hg38",
    help = "Human genome version [default = %default] (supports hg19 and hg38)"
  ),
  make_option("--no-normalize-logr",
    default = FALSE, action = "store_true",
    help = "Disable normalization of counts before calculating logR (seqz input only) [default = %default]"
  ),
  make_option(c("--overwrite"),
    default = TRUE, action = "store_true",
    help = "Overwrite output [default = %default]"
  ),
  make_option(c("--no-overwrite"),
    action = "store_false", dest = "overwrite",
    help = "Do not overwrite previous output"
  ),
  make_option("--reestimate-purity-ploidy",
    default = FALSE, action = "store_true",
    help = "Re-estimate purity and ploidy of each sample using ASCAT [default = %default]"
  ),
  make_option("--cached-dir",
    default = NA,
    help = "Set this to load intermediate data from the given directory"
  ),
  make_option("--debug-dir",
    default = NA, type="character",
    help = "Set this to save intermediate data to the given directory for later debugging"
  ),
  make_option(c("--validate-input"),
    default = TRUE, action = "store_true", dest = "validate_input",
    help = "Validate input files and abort if input is not valid"
  ),
  make_option(c("--no-validate-input"),
    action = "store_false", dest = "validate_input",
    help = "Do not validate input"
  ),
  make_option("--logging", type="character", default="DEBUG",
    help = "Set the logging level (DEBUG, INFO, WARN, STOP) [default = %default]"
  ),
  make_option(c("--seed"),
    default = 2134234,
    help = "Random seed"
  )
)

conf <- parse_args(OptionParser(option_list = option_list),
		   convert_hyphens_to_underscores = TRUE)

rp <- refphase(conf)

