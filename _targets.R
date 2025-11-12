library("conflicted")
library("MiscMetabar")
library("targets")
library("tarchetypes")
library("here")

here::i_am("_targets.R")
source(here("R/functions.R"))

seq_len_min <- 200
fw_primer_sequences <- ""
rev_primer_sequences <- ""
n_threads <- 4
refseq_file_name <- ""
sam_data_file_name <- "sam_data.csv"

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

  #> Paired end analysis
  #> ———————————————————

  ##> Remove primers
  tar_target(
    cutadapt,
    cutadapt_remove_primers(
      path_to_fastq = here("data/data_raw/rawseq/"),
      primer_fw = fw_primer_sequences,
      primer_rev = rev_primer_sequences,
      folder_output = here("data/data_intermediate/seq_wo_primers/"),
      args_before_cutadapt = "source ~/miniforge3/etc/profile.d/conda.sh && conda activate cutadaptenv && "
    )
  ),
  tar_target(data_raw, {
    cutadapt
    list_fastq_files(path = here::here("data/data_intermediate/seq_wo_primers/"),
                     paired_end = FALSE)
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
      derepR = derep_rs,
      minOverlap = 8,
      maxMismatch = 1
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
    sample_data_with_new_names(
      paste0(here::here(), "/", file_sam_data_csv),
      names_of_samples = id_sam,
      samples_order = na.omit(match(id_asv_table, id_sam)),
      dec = ","
    )
  ),
  tar_target(asv_tab, otu_table(
    rename_samples(otu_table(seqtab, taxa_are_rows = FALSE), names_of_samples = id_asv_table),
    taxa_are_rows = FALSE
  )),
  ##> Create the phyloseq object 'data_phyloseq' with
  ###   (i) table of asv,
  ###   ii) taxonomic table,
  ###   (iii) sample data and
  ###   (iv) references sequences
  tar_target(data_phyloseq, add_dna_to_phyloseq(
    phyloseq(asv_tab, sample_data(sam_tab), tax_table(
      as.matrix(tax_tab_maarjam_species, dimnames = rownames(tax_tab_maarjam_species))
    ))
  )),

  ##> Create post-clustering ASV into OTU using vsearch
  tar_target(d_vs, asv2otu(
    data_phyloseq, method = "vsearch", tax_adjust = 0
  )),
  ##> Clean post-clustering OTU using mumu
  tar_target(d_vs_mumu, mumu_pq(
    d_vs, method = "vsearch", tax_adjust = 0
  )),
  ##> Make a rarefied dataset
  tar_target(d_vs_mumu_rarefy, rarefy_even_depth(d_vs_mumu, rngseed = 22)),

  ##> Create the phyloseq object 'data_phyloseq' with
  tar_target(track_sequences_samples_clusters, track_wkflow(
    list(
      "Paired sequences" = seq_tab_Pairs,
      "Paired sequences without chimera" = seqtab_wo_chimera,
      "Paired sequences without chimera and longer than 200bp" = seqtab,
      "ASV denoising" = data_phyloseq,
      "OTU after vsearch reclustering at 97%" = d_vs,
      "OTU vs after mumu cleaning algorithm" = d_vs_mumu,
      "OTU vs + mumu + rarefaction by sequencing depth" = d_vs_mumu_rarefy
    )
  )),
  tar_target(track_by_samples, track_wkflow_samples(
    list(
      "ASV denoising" = data_phyloseq,
      "OTU after vsearch reclustering at 97%" = d_vs,
      "OTU vs after mumu cleaning algorithm" = d_vs_mumu,
      "OTU vs + mumu + rarefaction by sequencing depth" = d_vs_mumu_rarefy
    )
  )),

  ##> Build fastq quality report across the pipeline
  ### With raw sequences
  tar_target(
    quality_raw_seq,
    fastqc_agg("data/data_raw/rawseq/", qc.dir = "data/data_final/quality_fastqc/raw_seq/")
  ),
  tar_target(
    quality_raw_seq_plot,
    fastqc_plot("data/data_final/quality_fastqc/raw_seq/")
  ),
  ### After cutadapt
  tar_target(
    quality_seq_wo_primers,
    fastqc_agg("data/data_intermediate/seq_wo_primers/", qc.dir = "data/data_final/quality_fastqc/seq_wo_primers/")
  ),
  tar_target(
    quality_seq_wo_primers_plot,
    fastqc_plot("data/data_final/quality_fastqc/seq_wo_primers/")
  ),
  ### After filtering and trimming (separate report for forward and reverse)
  tar_target(
    quality_seq_filtered_trimmed_FW,
    fastqc_agg("data/data_intermediate/filterAndTrim_fwd", qc.dir = "data/data_final/quality_fastqc/filterAndTrim_fwd/")
  ),
  tar_target(
    quality_seq_filtered_trimmed_FW_plot,
    fastqc_plot("data/data_final/quality_fastqc/filterAndTrim_fwd")
  ),
  tar_target(
    quality_seq_filtered_trimmed_REV,
    fastqc_agg("data/data_intermediate/filterAndTrim_rev", qc.dir = "data/data_final/quality_fastqc/filterAndTrim_rev/")
  ),
  tar_target(
    quality_seq_filtered_trimmed_REV_plot,
    fastqc_plot("data/data_final/quality_fastqc/filterAndTrim_rev")
  )
)

