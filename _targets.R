library("conflicted")
library("MiscMetabar")
library("targets")
library("tarchetypes")
library("here")
library("autometric")

if (tar_active()) {
  log_start(
    path = "data/data_final/autometric_log.txt",
    seconds = 1
  )
}

here::i_am("_targets.R")
source(here("R/functions.R"))
lapply(list.files("~/Nextcloud/IdEst/Projets/MiscMetabar/R/", full.names = TRUE), source)

seq_len_min <- 200
fw_primer_sequences <- "GCATCGATGAAGAACGCAGC" 
rev_primer_sequences <- "TCCTCCGCTTATTGATATGC" 
n_threads <- 4
refseq_file_name <- "sh_general_release_dynamic_s_04.04.2024.fasta"
sam_data_file_name <- "sam_data.csv"
sample_col_name <- "Sample_names"
set.seed(22)


tar_plan(
  #> Place for file input
  tar_target(
    name = file_sam_data_csv,
    command = here("data/data_raw/metadata", sam_data_file_name),
    format = "file"
  ),

  tar_target(
    name = file_refseq_taxo,
    command = here("data/data_raw/refseq/", refseq_file_name),
    format = "file"
  ),

  tar_target(
    name = fastq_files_folder,
    command = here("data/data_raw/rawseq/"),
    format = "file"
  ),

  #> Match samples names from fastq files and metadata sam_data
  #> ———————————————————
  tar_target(s_d,
    sam_data_matching_names(
      path_sam_data = file_sam_data_csv,
      path_raw_seq =  fastq_files_folder,
      sample_col_name = sample_col_name,
      pattern_remove_fastq_files = "_L001.*",
      prefix = "samp_"
    )
  ),

  #> Paired end analysis
  #> ———————————————————

  ##> Remove primers
  tar_target(
    cutadapt,
    cutadapt_remove_primers(
      path_to_fastq = fastq_files_folder,
      primer_fw = fw_primer_sequences,
      primer_rev = rev_primer_sequences,
      folder_output = here("data/data_intermediate/seq_wo_primers/"),
      nproc = n_threads,
      return_file_path = TRUE,
      args_before_cutadapt = "source ~/miniforge3/etc/profile.d/conda.sh && conda activate cutadaptenv && "
    ),
    format = "file"
  ),
  tar_target(data_raw, {
    cutadapt
    list_fastq_files(path = here::here("data/data_intermediate/seq_wo_primers/"))
  }),

  ##> Classical dada2 pipeline
  tar_target(data_fnfs, data_raw$fnfs),
  tar_target(data_fnrs, data_raw$fnrs),
  ### Pre-filtered data with low stringency
  tar_target(
    filtered,
    filter_trim(
      output_fw = paste(
        getwd(),
        here("/data/data_intermediate/filterAndTrim_fwd"),
        sep = ""
      ),
      output_rev = paste(
        getwd(),
        here("/data/data_intermediate/filterAndTrim_rev"),
        sep = ""
      ),
      fw = data_fnfs,
      rev = data_fnrs,
      multithread = n_threads,
      compress = TRUE
    )
  ),

  ### Dereplicate fastq files
  tar_target(derep_fs, derepFastq(filtered[[1]]), format = "qs"),
  tar_target(derep_rs, derepFastq(filtered[[2]]), format = "qs"),
  ### Learns the error rates
  tar_target(err_fs, learnErrors(derep_fs, multithread = n_threads), format = "qs"),
  tar_target(err_rs, learnErrors(derep_rs, multithread = n_threads), format = "qs"),
  ### Make amplicon sequence variants
  tar_target(ddF, dada(derep_fs, err_fs, multithread = n_threads), format = "qs"),
  tar_target(ddR, dada(derep_rs, err_rs, multithread = n_threads), format = "qs"),
  ### Merge paired sequences
  tar_target(
    merged_seq,
    mergePairs(
      dadaF = ddF,
      dadaR = ddR,
      derepF = derep_fs,
      derepR = derep_rs
    ),
    format = "qs"
  ),
  ### Build a a table of ASV x Samples
  tar_target(seq_tab_Pairs, makeSequenceTable(merged_seq)),

  #> end Paired-end analysis
  #> ———————————————————————

  ##> Filtering sequences
  ### Remove chimera
  tar_target(seqtab_wo_chimera, chimera_removal_vs(seq_tab_Pairs)),
  ### Remove sequences based on length
  tar_target(seqtab, seqtab_wo_chimera[, nchar(colnames(seqtab_wo_chimera)) >= seq_len_min]),

  ##> Load sample data and rename samples
  tar_target(
    sam_tab,
    rename_samples(sample_data(s_d$sam_data), names_of_samples = s_d$sam_data$samples_names_common)
  ),
  tar_target(samp_n_otu_table,
      s_d$sam_names_matching$common_names[match(rownames(seqtab), s_d$sam_names_matching$raw_fastq)]
    ),

  tar_target(asv_tab, otu_table(
    rename_samples(
      otu_table(seqtab[!(duplicated(samp_n_otu_table) | duplicated(samp_n_otu_table, fromLast=TRUE)),], taxa_are_rows = FALSE),
        names_of_samples = samp_n_otu_table[!(duplicated(samp_n_otu_table) | duplicated(samp_n_otu_table, fromLast=TRUE))]),
    taxa_are_rows = FALSE
  )),

  tar_target(
    tax_tab,
    assignTaxonomy(
      seqtab,
      refFasta = file_refseq_taxo,
      taxLevels
      = c(
        "Kingdom",
        "Phyla",
        "Class",
        "Order",
        "Family",
        "Genus",
        "Species"
      ),
      multithread = n_threads
    )
  ),

  ##> Create the phyloseq object 'data_phyloseq' with
  ###   (i) table of asv,
  ###   ii) taxonomic table,
  ###   (iii) sample data and
  ###   (iv) references sequences


  tar_target(data_phyloseq, add_dna_to_phyloseq(
    phyloseq(asv_tab, sam_tab, tax_table(
      as.matrix(tax_tab, dimnames = rownames(tax_tab))
    ))
  )),

  tar_target(d_asv, 
    add_new_taxonomy_pq(data_phyloseq, ref_fasta = "data/data_raw/refseq/DADA2_EUK_SSU_v1.9_Fungi.fasta", suffix = "Eukaryome")
  ),

  ##> Create post-clustering ASV into OTU using vsearch
  tar_target(d_vs, asv2otu(
    d_asv, method = "vsearch", tax_adjust = 0
  )),
  ##> Clean post-clustering OTU using mumu
  tar_target(d_vs_mumu, mumu_pq(d_vs)$new_physeq),
  ##> Make a rarefied dataset
  tar_target(d_vs_mumu_rarefy, rarefy_even_depth(d_vs_mumu, sample.size = 2000)),

  ##> Create the phyloseq object 'd_asv' with
  tar_target(track_sequences_samples_clusters, track_wkflow(
    list(
      "Raw Forward sequences" = unlist(list_fastq_files(fastq_files_folder, paired_end = FALSE)),
      "Forward wo primers" = unlist(list_fastq_files(here::here("data/data_intermediate/seq_wo_primers/"), paired_end = FALSE)),
      "Forward sequences" = ddF,
      "Paired sequences" = seq_tab_Pairs,
      "Paired sequences without chimera" = seqtab_wo_chimera,
      "Paired sequences without chimera and longer than 200bp" = seqtab,
      "ASV denoising" = d_asv,
      "OTU after vsearch reclustering at 97%" = d_vs,
      "OTU vs after mumu cleaning algorithm" = d_vs_mumu,
      "OTU vs + mumu + rarefaction by sequencing depth" = d_vs_mumu_rarefy
    )
  )),
  tar_target(track_by_samples, track_wkflow_samples(
    list(
      "ASV denoising" = d_asv,
      "OTU after vsearch reclustering at 97%" = d_vs,
      "OTU vs after mumu cleaning algorithm" = d_vs_mumu,
      "OTU vs + mumu + rarefaction by sequencing depth" = d_vs_mumu_rarefy
    )
  )),
 
  ##> Build fastq quality report across the pipeline
  ### With raw sequences
  tar_target(
    quality_raw_seq,
    fastqc_agg(fastq_files_folder, qc.dir = here("data/data_final/quality_fastqc/raw_seq/"), multiqc=TRUE)
  ),
  ### After cutadapt
  tar_target(
    quality_seq_wo_primers,{
      cutadapt
    fastqc_agg(here("data/data_intermediate/seq_wo_primers/"), qc.dir = here("data/data_final/quality_fastqc/seq_wo_primers/"), multiqc=TRUE)
  }),
  ### After filtering and trimming (separate report for forward and reverse)
  tar_target(
    quality_seq_filtered_trimmed_FW,
    fastqc_agg(here(filtered[[1]]), qc.dir = here("data/data_final/quality_fastqc/filterAndTrim_fwd/"), multiqc=TRUE)
  ),
  tar_target(
    quality_seq_filtered_trimmed_REV,
    fastqc_agg(here(filtered[[1]]), qc.dir = here("data/data_final/quality_fastqc/filterAndTrim_rev/"), multiqc=TRUE)
  ),
  
  ##>  Build bioinformatic quarto report
  tar_target(bioinfo_report, {
    track_sequences_samples_clusters
    quarto::quarto_render(here::here("analysis", "01_bioinformatics.qmd"))
  }
  )#,
  # tar_target(build_website, {
  #   track_sequences_samples_clusters
  #   quarto::quarto_render(here::here())
  # }
  # )
)
