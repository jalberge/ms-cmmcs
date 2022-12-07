###### INIT

suppressPackageStartupMessages(library(circlize))

library(tidyverse)
library(data.table)
library(rstatix)
library(kableExtra)
library(RColorBrewer)
library(maftools)
library(cowplot)
library(readxl)
library(ggpubr)
library(scales)

setwd("~/Dropbox (Partners HealthCare)/Projects/ms-cmmcs/2_wgs/")

source("../0_utils/utils.R")
source("0_annotate_samples.R")

paired.participants <- paste0("CTF", 
                              c("001", "002", "003", 
                                "012", "013", "015", "016", "017", "018", "019", 
                                "021", "023", "025", 
                                "034", "036", 
                                "046", "047", "048", 
                                "050", "052", "053", "054", "055", "058"))

######## ANALYSIS

# MIXCR

full.files <- list.files("../data/mixcr/", pattern = ".clonotypes.ALL.txt", full.names = TRUE)
ids <- str_remove(basename(full.files), ".clonotypes.ALL.txt")

mixcr_clones <- lapply(full.files, read_mixcr_clones)
names(mixcr_clones) <- ids
mixr_clones_df <-
  rbindlist(mixcr_clones, use.names = TRUE, idcol = "sample_id")

summary.paired.clones <- mixr_clones_df |>
  inner_join(subset(clinicaldata, tissue != "WBCs"),
            by = c("sample_id" = "entity:sample_id")) |>
  dplyr::select(
    participant,
    tissue,
    sample_id,
    isotype,
    cloneId,
    cloneCount,
    cloneFraction,
    allVHitsWithScore,
    allDHitsWithScore,
    allJHitsWithScore,
    minQualCDR3,
    aaSeqCDR3
  ) |>
  rowwise() |>
  mutate(
    V_seg = extract_first_vdj(allVHitsWithScore),
    D_seg = extract_first_vdj(allDHitsWithScore),
    J_seg = extract_first_vdj(allJHitsWithScore)
  ) |>
  mutate_at(vars(cloneCount, cloneFraction), as.numeric) |>
  filter(participant %in% paired.participants) |>
  pivot_wider(
    id_cols = c(participant, aaSeqCDR3, V_seg, D_seg, J_seg, isotype),
    names_from = tissue,
    values_from = c(cloneCount, cloneFraction),
    values_fill = 0
  )

summary.paired.clones |>
  filter(isotype == "IGH") |>
  select(participant, V_seg, D_seg, J_seg, aaSeqCDR3, cloneFraction_BMPCs, cloneFraction_CMMCs) %>%
  mutate(cloneFraction_BMPCs=paste0(floor(100*cloneFraction_BMPCs), "%"), cloneFraction_CMMCs=paste0(floor(100*cloneFraction_CMMCs), "%")) |>
  write_tsv("data/mixcr_out_matched.txt")


