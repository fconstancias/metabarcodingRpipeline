#!/usr/bin/env Rscript

rm(list = ls())

cat("\n############################################################\n")
cat("Starting \n")
cat("############################################################\n\n")

## ------------------------------------------------------------------------

require("optparse")

## ------------------------------------------------------------------------
suppressPackageStartupMessages(library("optparse"))


option_list = list(
  make_option(c("-i", "--input_directory"), type="character", default=NULL, 
              help="Path of the input directory containing reads in their respective run sub-directories \n
              e.g., -i raw [contains raw/run1 and raw/run2]\n
              N.B.: sample name is extracted from .fastq.gz samples  before the first '_' e.g., XXX-XX \n
              sample names must be unambigious and unique e.g., sample-1 = sample-11, sample-01 != sample-10", metavar="character"),
  make_option(c("-r", "--primer_removed_directory"), type="character", default = NULL, 
              help="Directory containing primer trimmed/demultiplexed samples", metavar="character"),
  make_option(c("-o", "--output_directory"), type="character", default= "dada2", 
              help="Name of the output directory", metavar="character"),
  make_option(c("-c", "--cut_file_pattern"), type="character", default= NULL, 
              help="", metavar="character"),
  make_option(c("--nbases"), type="numeric", default= 20000000, 
              help="nbases for error learning", metavar="character"),
  make_option(c("--trunclen"), type="character", default = "c(260,250)", 
              help="Nucleotide position to truncate the Fwd and Rev reads at [if using -V V4 or V3V4, this parameter is already set]", metavar="numeric"),
  make_option(c("--maxee"), type="character", default = "c(3,4)", 
              help="Maximum expected error for Fwd and Rev reads [if using -V V4 or V3V4, this parameter is already set]", metavar="numeric"),
  make_option(c("--minLen"), type="numeric", default = 225, 
              help="Minimum read length [if using -V V3 or V3V4, this parameter is already set]", metavar="numeric"),
  make_option(c("--maxLen"), type="numeric", default = 475, 
              help="MAximum read length [if using -V V3 or V3V4, this parameter is already set]", metavar="numeric"),
  make_option(c("-T", "--slots"), type="numeric", default = 6, 
              help="Number of threads to perform the analyses", metavar="numeric"),
  make_option(c("-f", "--fun_dir"), type="character", default= NULL, 
              help="Directory containing the R functions", metavar="character")
  
); 
opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser)

## ------------------------------------------------------------------------
# parse_args(opt_parser, args = c("--help"))
# 
# if(is.null(opt$input_directory)) {
#   print_help(opt_parser)
#   stop("You must give an input directory (-i)")
# }
# 
# if(is.null(opt$atropos_binary)) {
#   print_help(opt_parser)
#   stop("You must provide path to atropos program (-a)")
# }

source(opt$fun_dir)

opt$trunclen <- as.vector(strsplit(opt$trunclen, ","))
opt$maxee <- as.vector(strsplit(opt$maxee, ","))

## ------------------------------------------------------------------------
cat(paste0('\n# Input directory: ',opt$input_directory,'.\n'))
cat(paste0('\n# Directory containing read without PCR primers: ',opt$primer_removed_director,'.\n'))
cat(paste0('\n# trunclen: ',opt$trunclen,'.\n'))
cat(paste0('\n# maxee: ',opt$maxee,'.\n'))
cat(paste0('\n# slots: ',opt$slots,'.\n'))
cat(paste0('\n# maxLen: ',opt$maxLen,'.\n'))

## ------------------------------------------------------------------------

run_dada2_filter_denoise_merge_reads(raw_files_path = opt$input_directory,
                                     cut_dir = opt$primer_removed_directory,
                                     cut_file_pattern = opt$cut_file_pattern,
                                     output  = opt$output_directory,
                                     trunclen = opt$trunclen,
                                     maxLen = opt$maxLen,
                                     maxee = opt$maxee,
                                     truncQ = 6,
                                     minLen = opt$minLen, #250 #350 in initial analysis which makes sense according to read/ASV distrib and 16S V3V4 distrib
                                     nthreads = opt$slots,
                                     nbases = opt$nbases,
                                     pool = "pseudo",
                                     minover = 12,
                                     filt_dir = "02_dada2_filtered_denoised_merged",
                                     filt_pattern = c("_filtered.fastq.gz"),
                                     seed_value = 123,
                                     print_plot = FALSE)
## ------------------------------------------------------------------------


cat("\n\nHere some Info for the record keeping:\n\n")
print(sessionInfo())

cat("Parameters Used:\n\n")

cat(paste0('\n# Input directory: ',opt$input_directory,'.\n'))
cat(paste0('\n# Directory containing read without PCR primers: ',opt$primer_removed_director,'.\n'))
cat(paste0('\n# trunclen: ',opt$trunclen,'.\n'))
cat(paste0('\n# trim_length: ',opt$trim_length,'.\n'))
cat(paste0('\n# maxee: ',opt$maxee,'.\n'))
cat(paste0('\n# slots: ',opt$slots,'.\n'))
cat(paste0('\n# maxLen: ',opt$maxLen,'.\n'))

cat("\n############################################################\n")
cat("All done \n")
cat("############################################################\n\n")

print(Sys.time())

## ------------------------------------------------------------------------
