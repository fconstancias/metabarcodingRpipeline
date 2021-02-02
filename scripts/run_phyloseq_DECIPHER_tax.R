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
  make_option(c("-p", "--phyloseq_path"), type="character", default = NULL, 
              help="Path of the input phyloseq object", metavar="character"),
  make_option(c("-o", "--output"), type="character", default= "~/", 
              help="Name of the output directory - If FALSE do not export - which does not make much sense in that context", metavar="character"),
  make_option(c("--tax_threshold"), type="numeric", default = 60, 
              help="Confidence threshold", metavar="numeric"),
  make_option(c("--db"), type="character", default= NULL, 
              help="path of database", metavar="character"),
  make_option(c("--db_species"), type="character", default= NULL, 
              help="path of species level database", metavar="character"),
  make_option(c("--reverse_comp"), type="character", default= TRUE, 
              help="Reverse complement sequences for taxonomic assignments?", metavar="character"),
  make_option(c("-T", "--slots"), type="numeric", default = 6, 
              help="Number of threads to perform the analyses", metavar="character"),
  make_option(c("-f", "--fun_dir"), type="character", default= NULL, 
              help="Directory containing the R functions", metavar="character"),
  make_option(c("--seed"), type="numeric", default= 123, 
              help="seed number for random generator", metavar="character")
  
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
cat(paste0('\n# Input phyloseq object: ',opt$phyloseq_path,'.\n'))
cat(paste0('\n# Output directory: ',opt$export,'.\n'))

cat(paste0('\n# Confidence threshold: ',opt$tax_threshold,'.\n'))
cat(paste0('\n# database used : ',opt$db,'.\n'))
cat(paste0('\n# database used for species level assignments: ',opt$db_species,'.\n'))

cat(paste0('\n# Reverse complement? ',opt$reverse_comp,'.\n'))

cat(paste0('\n# Threads: ',opt$slots,'.\n'))

cat(paste0('\n# Seed for reproducibility: ',opt$seed,'.\n'))
## ------------------------------------------------------------------------


phyloseq_DECIPHER_tax(physeq = opt$phyloseq_path, # readRDS("data/processed/physeq_update_11_1_21.RDS") %>% subset_taxa(taxa_sums(physeq) > 100000) -> physeq
                      threshold = opt$tax_threshold, # 60 (very high),  50 (high), PB = 10
                      db = opt$db, # db ="~/db/DADA2/SILVA_SSU_r132_March2018.RData" db ="~/db/DADA2/Databases_pbHITdb_pbHITdb_v1.0.0_20180919_IDTAXA.RData" db ="~/db/DADA2/GTDB_r95-mod_August2020.RData" db ="~/db/DADA2/SILVA_SSU_r138_2019.RData" db ="~/db/DADA2/Fungal_LSU_v11_March2018.RData"
                      nthreads = opt$slots,
                      tryRC = s.logical(opt$reverse_comp),
                      export = opt$output,
                      seed_value = opt$seed,
                      bootlabel = "_Confidence",
                      tax_ranks = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
                      return = FALSE)

## ------------------------------------------------------------------------


cat("\n\nHere some Info for the record keeping:\n\n")
print(sessionInfo())

cat("Parameters Used:\n\n")

## ------------------------------------------------------------------------
cat(paste0('\n# Input phyloseq object: ',opt$phyloseq_path,'.\n'))
cat(paste0('\n# Output directory: ',opt$export,'.\n'))

cat(paste0('\n# Bootstrap threshold: ',opt$tax_threshold,'.\n'))
cat(paste0('\n# database used : ',opt$db,'.\n'))
cat(paste0('\n# database used for species level assignments: ',opt$db_species,'.\n'))

cat(paste0('\n# Reverse complement? ',opt$reverse_comp,'.\n'))

cat(paste0('\n# Threads: ',opt$slots,'.\n'))

cat(paste0('\n# Seed for reproducibility: ',opt$seed,'.\n'))

cat("\n############################################################\n")
cat("All done \n")
cat("############################################################\n\n")

print(Sys.time())

## ------------------------------------------------------------------------



