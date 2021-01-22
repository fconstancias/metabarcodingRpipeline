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
  make_option(c("-p", "--prop_sample"), type="numeric", default = 20, 
              help="Proportion of samples to consider for quality plot", metavar="character"),
  make_option(c("-c", "--cut_file_pattern"), type="character", default= c("_primersout_R1_.fastq.gz","_primersout_R2_.fastq.gz"), 
              help="Pattern of primer removed files", metavar="character"),
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

## ------------------------------------------------------------------------
cat(paste0('\n# Input directory: ',opt$input_directory,'.\n'))

## ------------------------------------------------------------------------

run_dada2_qplot(raw_files_path = opt$input_directory,
                cut_dir = opt$primer_removed_directory,
                cut_file_pattern = opt$cut_file_pattern,
                output  = opt$output_directory,
                prop.sample = opt$prop_sample)

## ------------------------------------------------------------------------


cat("\n\nHere some Info for the record keeping:\n\n")
print(sessionInfo())

cat("Parameters Used:\n\n")

cat(paste0('\n# Input directory: ',opt$input_directory,'.\n'))

cat("\n############################################################\n")
cat("All done \n")
cat("############################################################\n\n")

print(Sys.time())

## ------------------------------------------------------------------------
