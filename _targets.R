source("R/functions.R")
library("conflicted")
library("MiscMetabar")
library("targets")
library("tarchetypes")

here::i_am("_targets.R")

tar_plan(
  tar_target(
    name = file_sam_data_csv,
    command = "data/data_raw/metadata/sam_data.csv",
    format = "file"
  ),

  tar_target(
    name = file_refseq_taxo,
    command = "data/data_raw/refseq/XXX",
    format = "file"
  ),

  tar_target(
    name = file_refseq_taxo,
    command = "data/data_raw/refseq/XXX",
    format = "file"
  ),


)
