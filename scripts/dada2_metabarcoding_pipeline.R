#!/usr/bin/env Rscript

rm(list = ls())

cat("\n############################################################\n")
cat("Starting \n")
cat("############################################################\n\n")

## ------------------------------------------------------------------------

require("optparse")
source("https://raw.githubusercontent.com/fconstancias/metabarcodingRpipeline/master/scripts/functions.R")

## ------------------------------------------------------------------------
suppressPackageStartupMessages(library("optparse"))

option_list = list(
  make_option(c("-i", "--input_directory"), type="character", default=NULL, 
              help="Path of the input directory containing raw _R1_ and _R2_ raw reads in their respective run sub-directories \n
              e.g., -i raw [contains raw/run1 and raw/run2]\n
              N.B.: sample name is extracted from .fastq.gz samples  before the first '_' e.g., XXX-XX \n
              sample names must be unambigious and unique e.g., sample-1 = sample-11, sample-01 != sample-10", metavar="character"),
  make_option(c("-a", "--atropos_binary"), type="character", default = "atropos", 
              help="Path of atropos program [used for primer removal]", metavar="character"),
  make_option(c("-o", "--output_directory"), type="character", default= "dada2", 
              help="Name of the output directory", metavar="character"),
  make_option(c("-V", "--pipeline"), type="character", default = "V4", 
              help="V4 or V3V4 will use default primers and parameters as used in the FBT lab [primers, trunc, maxee, overlap, expected length, ...]", metavar="character"),
  make_option(c("-t", "--tax_method"), type="character", default="dada", 
              help="User can specify using dada2 (=dada) or DECIPHER for taxonomic assignments [default dada]", metavar="character"),
  make_option(c("--tax_threshold"), type="numeric", default= 60, 
              help="Thershold for taxonomic assignments [if --tax_metod dada: the minimum bootstrap confidence for assigning a taxonomic level. if --tax_method DECIPHER: Numeric specifying the confidence at which to truncate the output taxonomic classifications. ]", metavar="numeric"),
  make_option(c("--metadata"), type="character", default=NULL, 
              help="Path to excel document containing metadata [Sample identifier column should be sample_name]", metavar="character"),
  make_option(c("--database"), type="character", default = "~/db/DADA2/silva_nr_v138_train_set.fa.gz", 
              help="Path to the taxonomic database", metavar="character"),
  make_option(c("--database_for_species_assignments"), type="character", default = "~/db/DADA2/silva_species_assignment_v138.fa.gz", 
              help="Path to the speies-level taxonomic database [only for --tax_metod  dada]", metavar="character"),
  make_option(c("--phylo"), type="logical", default = "FALSE", 
              help="Compute phylogenetic tree from the ASV sequence ?", metavar="character"),
  make_option(c("--PRIMER_F"), type="character", default = NULL, 
              help="Sequence of the gene specific Fwd primer to be removed with atropos [if using -V V4 or V3V4, this parameter is already set]", metavar="character"),
  make_option(c("--PRIMER_R"), type="character", default = NULL, 
              help="Sequence of the gene specific Rev primer to be removed with atropos [if using -V V4 or V3V4, this parameter is already set]", metavar="character"),
  make_option(c("--minover"), type="numeric", default = 15, 
              help="Minimum overlap for merginf R1 and R2 reads [if using -V V4 or V3V4, this parameter is already set]", metavar="numeric"),
  make_option(c("--trunclen"), type="character", default=c(260,250), 
              help="Nucleotide position to truncate the Fwd and Rev reads at [if using -V V4 or V3V4, this parameter is already set]", metavar="numeric"),
  make_option(c("--trim_length"), type="character", default = c(240,400), 
              help="ASV of length outside the range will be discarded [i.e., insilco size exclusion of ASV - if using -V V3 or V3V4, this parameter is already set]", metavar="numeric"),
  make_option(c("--maxee"), type="character", default = c(3,4), 
              help="Maximum expected error for Fwd and Rev reads [if using -V V4 or V3V4, this parameter is already set]", metavar="numeric"),
  make_option(c("--minLen"), type="numeric", default = 100, 
              help="Minimul read length [if using -V V3 or V3V4, this parameter is already set]", metavar="numeric"),
  make_option(c("-T", "--slots"), type="numeric", default = 4, 
              help="Number of threads to perform the analyses", metavar="numeric")
  
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

run_phylo <- as.logical(opt$phylo)

## ------------------------------------------------------------------------
cat(paste0('\n# Input directory: ',opt$input_directory,'.\n'))
cat(paste0('\n# Tax_method: ',opt$tax_method,'.\n'))
cat(paste0('\n# database: ',opt$database,'.\n'))
cat(paste0('\n# database_for_species_assignments: ',opt$database_for_species_assignments,'.\n'))
## ------------------------------------------------------------------------

run_16S_pipe(raw_files_path = opt$input_directory,
             atropos = opt$atropos_binary,
             out_dir = opt$output_directory,
             V = opt$pipeline,
             tax_method = opt$tax_method,
             metadata = opt$metadata,
             db = opt$database,
             db_species = opt$database_for_species_assignments,
             run_phylo = run_phylo,
             tax_threshold = opt$tax_threshold,
             PRIMER_F = opt$PRIMER_F,
             PRIMER_R = opt$PRIMER_R,
             minover = opt$minover,
             trunclen = opt$trunclen,
             trim_length = opt$trim_length,
             maxee = opt$maxee,
             minLen = opt$minLen,
             SLOTS = opt$slots)

## ------------------------------------------------------------------------


cat("\n\nHere some Info for the record keeping:\n\n")
print(sessionInfo())

cat("Parameters Used:\n\n")

cat(paste0('\n# Input directory: ',opt$input_directory,'.\n'))
cat(paste0('\n# Tax_method: ',opt$tax_method,'.\n'))
cat(paste0('\n# database: ',opt$database,'.\n'))
cat(paste0('\n# database_for_species_assignments: ',opt$database_for_species_assignments,'.\n'))

cat("\n############################################################\n")
cat("All done \n")
cat("############################################################\n\n")

print(Sys.time())

## ------------------------------------------------------------------------
