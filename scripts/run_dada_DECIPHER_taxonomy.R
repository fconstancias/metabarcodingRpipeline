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
  make_option(c("-o", "--output_directory"), type="character", default= "dada2", 
              help="Name of the output directory", metavar="character"),
  make_option(c("-t", "--tax_method"), type="character", default= "dada", 
              help="Name of the output directory", metavar="character"),
  make_option(c("--db"), type="character", default= NULL, 
              help="path of database", metavar="character"),
  make_option(c("--db_species"), type="character", default= NULL, 
              help="path of species level database", metavar="character"),
  make_option(c("--reverse_comp"), type="character", default= TRUE, 
              help="Name of the output directory", metavar="character"),
  make_option(c("--tax_threshold"), type="numeric", default = 60, 
              help="ASV of length outside the range will be discarded [i.e., insilco size exclusion of ASV - if using -V V3 or V3V4, this parameter is already set]", metavar="numeric"),
  make_option(c("--output_bootstraps"), type="character", default = TRUE, 
              help=" ", metavar="numeric"),
  make_option(c("--allow_multiple_assignments"), type="character", default = TRUE, 
              help=" ", metavar="character"),
  make_option(c("-T", "--slots"), type="numeric", default = 6, 
              help="Number of threads to perform the analyses", metavar="character"),
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
cat(paste0('\n# Output directory: ',opt$output_directory,'.\n'))

cat(paste0('\n# tax assignment method: ',opt$tax_method,'.\n'))
cat(paste0('\n# database used : ',opt$db,'.\n'))
cat(paste0('\n# database used for species level assignments: ',opt$db_species,'.\n'))

cat(paste0('\n# Reverse complement? ',opt$reverse_comp,'.\n'))
cat(paste0('\n# Export bootstrap values? ',opt$output_bootstraps,'.\n'))
cat(paste0('\n# Allow multiple species level assignments? ',opt$allow_multiple_assignments,'.\n'))


cat(paste0('\n# Threads: ',opt$slots,'.\n'))

## ------------------------------------------------------------------------

run_dada_DECIPHER_taxonomy(raw_files_path = opt$input_directory,
                           taxa_dir = "04_dada2_taxonomy",
                           method = opt$tax_method, # method = "DECIPHER" or "dada" 
                           threshold = opt$tax_threshold,#60,  # used for DECIPHER and dada2 if outputBootstraps = FALSE
                           tryRC = as.logical(opt$reverse_comp),
                           db = opt$db, # "~/db/DADA2/silva_nr99_v138_train_set.fa.gz",
                           db_species = opt$db_species,#"~/db/DADA2/silva_species_assignment_v138.fa.gz",
                           nthreads = opt$slots,
                           outputBootstraps = as.logical(opt$output_bootstraps), #3 only for dada2 method# outputBootstraps <- TRUE 
                           allowMultiple = as.logical(opt$allow_multiple_assignments),
                           merged_run_dir = "03_dada2_merged_runs_chimera_removed",
                           output = opt$output_directory,
                           seed_value = 123)

## ------------------------------------------------------------------------


cat("\n\nHere some Info for the record keeping:\n\n")
print(sessionInfo())

cat("Parameters Used:\n\n")

## ------------------------------------------------------------------------
cat(paste0('\n# Input directory: ',opt$input_directory,'.\n'))
cat(paste0('\n# Output directory: ',opt$output_directory,'.\n'))

cat(paste0('\n# tax assignment method: ',opt$tax_method,'.\n'))
cat(paste0('\n# database used : ',opt$db,'.\n'))
cat(paste0('\n# database used for species level assignments: ',opt$db_species,'.\n'))

cat(paste0('\n# Reverse complement? ',opt$reverse_comp,'.\n'))
cat(paste0('\n# Export bootstrap values? ',opt$output_bootstraps,'.\n'))
cat(paste0('\n# Allow multiple species level assignments? ',opt$allow_multiple_assignments,'.\n'))


cat(paste0('\n# Threads: ',opt$slots,'.\n'))


cat(paste0('\n# Threads: ',opt$slots,'.\n'))
cat("\n############################################################\n")
cat("All done \n")
cat("############################################################\n\n")

print(Sys.time())

## ------------------------------------------------------------------------



