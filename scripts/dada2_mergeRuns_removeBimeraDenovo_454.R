#' @title ...
#' @param .
#' @param ..
#' @author Florentin Constancias
#' @note .
#' @note .
#' @note .
#' @return .
#' @export
#' @examples
#'
#'
#'
#'
#'
#'

run_dada2_mergeRuns_removeBimeraDenovo <- function(seqtab = NULL,
                                                   track = NULL,
                                                   merged_run_dir = "dada2/03_dada2_merged_runs_chimera_removed",
                                                   chimera_method = "consensus",
                                                   trim_length,
                                                   nthreads = 6,
                                                   collapseNoMis = FALSE,
                                                   minOverlap = 20, # https://github.com/benjjneb/dada2/issues/518
                                                   filt_dir = "dada2/02_dada2_filtered_denoised_merged",
                                                   export = TRUE,
                                                   seed_value = 123,
                                                   return = TRUE){
  ## ------------------------------------------------------------------------
  require(tidyverse); require(dada2); require(phyloseq)
  cat(paste0('\n##',"You are using DADA2 version ", packageVersion('dada2'),'\n'))
  cat(paste0('\n##',"You are using tidyverse version ", packageVersion('tidyverse'),'\n\n'))
  cat(paste0('\n##',"You are using phyloseq version ", packageVersion('phyloseq'),'\n\n'))
  
  cat('################################\n\n')
  
  
  ## ------------------------------------------------------------------------
  
  if(is.null(seqtab)  & is.null(track) ){
    
    merged_run_path <- file.path(merged_run_dir)
    
    dir.create(merged_run_path, showWarnings = TRUE, recursive = TRUE)
    
    seqtables <- sort(list.files(filt_dir,
                                 pattern = glob2rx("*_seqtab.rds*"),
                                 full.names = TRUE,
                                 recursive = TRUE))
    
    track <- sort(list.files(filt_dir,
                             pattern = glob2rx("*_track_analysis.tsv*"),
                             full.names = TRUE,
                             recursive = TRUE))
    
    summary <- map(track, 
                   read_tsv) %>% bind_rows()
    
    cat(str_c('\n# removeBimeraDenovo running on ', seqtables, ' file '))
    
    # https://github.com/benjjneb/dada2/issues/345
    
    
    if(seqtables %>% length > 1){
      list.df <- map(seqtables, readRDS)
      
      st.all <- mergeSequenceTables(tables = list.df)
      
      cat('\n# mergeSequenceTables done\n')
      
    }else{
      st.all = seqtables %>%
        readRDS()
      
      cat('\n# only one SequenceTable, no merging to do\n')
      
    }
  }
  ## ------------------------------------------------------------------------
  
  if(!is.null(seqtab) == TRUE & !is.null(track) == TRUE){
    
    summary <- track %>% bind_rows()
    
    if(seqtab %>% length > 1){
      
      st.all <- mergeSequenceTables(tables = seqtab)
      
      
      
      cat('\n# mergeSequenceTables done\n')
      
    }else{
      st.all = seqtab
      
      cat('\n# only one SequenceTable, no merging to do\n')
      
    }
  }
  
  ## ------------------------------------------------------------------------
  cat('\n# removeBimeraDenovo start\n')
  
  seqtab.raw <- removeBimeraDenovo(st.all, method = chimera_method,
                                   multithread = nthreads, 
                                   verbose = FALSE)
  
  cat('\n# removeBimeraDenovo done\n')
  
  num.chimera.removed <- ncol(st.all) - ncol(seqtab.raw)
  perc.num.chimera.removed <- round(100*num.chimera.removed/ncol(st.all) %>% round(2),2)
  reads.chimera.removed <- sum(colSums(st.all)) - sum(colSums(seqtab.raw))
  perc.reads.chimera.removed <- round(100*reads.chimera.removed/sum(colSums(st.all)) %>% round(2),2)
  
  cat(paste0('# ',num.chimera.removed," chimera were found and removed\n"))
  cat(paste0('# These represent ',perc.num.chimera.removed,'% of total ASVs and ',perc.reads.chimera.removed,'% of total reads\n'))
  
  ## ------------------------------------------------------------------------
  # Distribution of variants
  cat("\n# The variants (ASVs) have the following length distribution:\n")
  table(nchar(getSequences(seqtab.raw)))
  
  cat(paste0('\n# Reads shorter than ',trim_length[1],'bp and longer than ',trim_length[2], 'bp are going to be removed.\n'))
  
  
  plotLengthDistro <- function(st) {
    tot.svs <- table(nchar(colnames(st)))
    tot.reads <- tapply(colSums(st), nchar(colnames(st)), sum)
    df <- data.frame(Length=as.integer(c(names(tot.svs), names(tot.reads))),
                     Count=c(tot.svs, tot.reads),
                     Type=rep(c("ASVs", "Reads"), times=c(length(tot.svs), length(tot.reads))))
    
    p <- ggplot(data=df, aes(x=Length, y=Count, color=Type)) + geom_point() + facet_wrap(~Type, scales="free_y") + theme_bw() + xlab("Amplicon Length")
    
    return(p)
    
  }
  
  plotLengthDistro(seqtab.raw) + scale_y_log10() + 
    ggtitle(str_c("Overall Sequence / ASV length distribution ")) +
    geom_vline(xintercept = trim_length[1], size = 0.1, colour = "red", alpha = 0.8, linetype = 2, show.legend = "min") +
    geom_vline(xintercept = trim_length[2], size = 0.1, colour = "red", alpha = 0.8, linetype = 2 ) -> plot
  
  
  if (export == TRUE){
    ggsave(str_c(merged_run_dir,"/","seqtab_distrib",".pdf"),plot=plot, width = 9, height = 8)
  }
  ## ------------------------------------------------------------------------
  
  
  # Trim the unespecific amplifications from our dataset
  seqtab <- seqtab.raw[,nchar(colnames(seqtab.raw)) %in% seq(trim_length[1],
                                                             trim_length[2])]
  
  
  cat(paste0('\n# Reads shorter than ',trim_length[1],'bp and longer than ',trim_length[2], 'bp were removed.\n'))
  cat("\n# The variants (ASVs) after length filtering have the following length distribution:\n")
  table(nchar(getSequences(seqtab)))
  
  
  cat(paste0("\n# A total of ", round((sum(colSums(seqtab)) * 100) / sum(colSums(seqtab.raw)), digits = 2), "% reads were kept after length filtering.\n\n"))
  cat(paste0("\n# A total of ", round((dim(seqtab)[2] * 100) / dim(seqtab.raw)[2], digits = 2), "% ASVs were kept after length filtering.\n\n"))
  
  
  ## ------------------------------------------------------------------------
  full_join(summary, 
            data.frame(sample = rownames(st.all),
                       tabled_joined = rowSums(st.all),
                       chimera_out = rowSums(seqtab.raw),
                       length_filtered = rowSums(seqtab)), by='sample') %>% 
    # mutate(tabled_pc = round(tabled_joined /merged, 2)) %>%
    mutate(chimera_out_pc = round(chimera_out/tabled, 2)) %>% 
    mutate(length_filtered_pc = round(length_filtered/chimera_out, 2)) -> track
  
  #write_tsv(data.frame(track),str_c(output,"/",name,"_track_analysis.tsv"))
  
  cat("\n# The distribution of chimera reads kept is the following:\n")
  summary(track$chimera_out_pc)
  
  
  ## ------------------------------------------------------------------------
  if(collapseNoMis==TRUE){
    cat(str_c('\n# Saving uncollapsed .rds and fasta files as well as summary .tsv \n'))
    
    if (export == TRUE){
      saveRDS(seqtab, str_c(merged_run_dir,"/uncollapsed_no-chim-seqtab.rds"))
      
      uniquesToFasta(seqtab,
                     str_c(merged_run_dir,"/uncollapsed_no-chim-seqtab.fasta"),
                     ids= str_c("asv",c(1:ncol(seqtab)), ";size=", colSums(seqtab)))
      
      write_tsv(track, str_c(merged_run_dir,"/uncollapsed_track_analysis.tsv"))
    }
    cat('\n# You have decided to run collapseNoMismatch on your dataset. Please note that it is only helpful IF you are working with several sequencing runs and it might take a long time to run. You might want to go further (taxonomy, ...) on your uncollapsed seqtable while it is running \n')
    
    collapsed_100 <- collapseNoMismatch(seqtab,  
                                        minOverlap = minOverlap,
                                        identicalOnly = FALSE)
    
    physeq <- merge_phyloseq(otu_table(t(collapsed_100), taxa_are_rows=TRUE))
    ASV_seq <- Biostrings::DNAStringSet(taxa_names(physeq))
    names(ASV_seq) <- taxa_names(physeq)
    
    physeq <- merge_phyloseq(physeq, 
                             ASV_seq)
    
    taxa_names(physeq) <- paste0("ASV", str_pad(seq(ntaxa(physeq)),
                                                nchar(ntaxa(physeq)),
                                                pad = "0"))
    
    cat(str_c('\n# Saving collapsed .rds and fasta files as well as summary .tsv  \n'))
    
    if (export == TRUE){
      saveRDS(collapsed_100, str_c(merged_run_dir,"/minOverlap_",minOverlap,"_collapse_no_mismatch_no-chim-seqtab.rds"))
      
      saveRDS(physeq, str_c(merged_run_dir,"/physeq.rds"))
      
      uniquesToFasta(collapsed_100,
                     str_c(merged_run_dir,"/minOverlap_",minOverlap,"_collapse_no_mismatch_no-chim-seqtab.fasta"),
                     ids= str_c("asv",c(1:ncol(collapsed_100)), ";size=", colSums(collapsed_100)))
    }
    track %>% 
      mutate(collapsed_100 = rowSums(collapsed_100)) %>%
      mutate(collapsed_100_pc = round(collapsed_100 / length_filtered, digits = 10)) -> track.final
    
    if (export == TRUE){
      
      write_tsv(track.final, str_c(merged_run_dir,"/track_analysis.tsv"))
      
      cat(str_c('# Your final 100% clustered ASV table can be found in "', paste0(merged_run_dir,"/seqtab.rds"),'"\n'))
      cat(str_c('# A FASTA file with your final ASVs was written in "',paste0(merged_run_dir,"/seqtab.fasta"), '"\n'))
      
      cat(str_c('# In "',paste0(merged_run_dir,"/track_analysis.tsv"),"\" you will find a table where you can check the loss of reads in each step. Check it out to see if everything's correct!",'\n'))
      # cat(str_c('# You have to copy them to your local computer using "scp [your.user.id]@euler.ethz.ch:',str_c(output,"/",name,"_track_analysis.tsv "),'." and go further.\n'))
    }
  }
  ## ------------------------------------------------------------------------
  if(collapseNoMis==FALSE){
    
    physeq <- merge_phyloseq(otu_table(t(seqtab), taxa_are_rows=TRUE))
    ASV_seq <- Biostrings::DNAStringSet(taxa_names(physeq))
    names(ASV_seq) <- taxa_names(physeq)
    
    physeq <- merge_phyloseq(physeq, 
                             ASV_seq)
    
    taxa_names(physeq) <- paste0("ASV", str_pad(seq(ntaxa(physeq)),
                                                nchar(ntaxa(physeq)),
                                                pad = "0"))    
    
    if (export == TRUE){
      cat(str_c('\n# Saving .rds and fasta files as well as summary .tsv \n'))
      
      saveRDS(seqtab, str_c(merged_run_dir,"/no-chim-seqtab.rds"))
      
      saveRDS(physeq, str_c(merged_run_dir,"/physeq.rds"))
      
      uniquesToFasta(seqtab,
                     str_c(merged_run_dir,"/no-chim-seqtab.fasta"),
                     ids= str_c("asv",c(1:ncol(seqtab)), ";size=", colSums(seqtab)))
      
      write_tsv(track, str_c(merged_run_dir,"/track_analysis.tsv"))
      
      cat(str_c('# Your final ASV table can be found in "', paste0(merged_run_dir,"/seqtab.rds"),'"\n'))
      cat(str_c('# A FASTA file with your final ASVs was written in "',paste0(merged_run_dir,"/seqtab.fasta"), '"\n'))
      
      cat(str_c('# In "',paste0(merged_run_dir,"/track_analysis.tsv"),"\" you will find a table where you can check the loss of reads in each step. Check it out to see if everything's correct!",'\n'))
      # cat(str_c('# You have to copy them to your local computer using "scp [your.user.id]@euler.ethz.ch:',str_c(output,"/",name,"_track_analysis.tsv "),'." and go further.\n'))
    }
    cat('\n# chimera removal step is done. You can go further into the analysis (taxonomy, phylogeny) or explore collapseNoMismatch IF you are dealing with multiple runs ... it might take very long \n\n')
  }
  ## ------------------------------------------------------------------------
  
  if (return == TRUE){
    if(collapseNoMis==TRUE){
      return(list("seqtab" = collapsed_100,
                  "track" = track,
                  "plot" = plot,
                  "physeq" = physeq))
    }else{
      return(list("seqtab" = seqtab,
                  "track" = track,
                  "plot" = plot,
                  "physeq" = physeq))
    }
  }
  
}

