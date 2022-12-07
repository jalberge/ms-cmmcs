# Start here and load sample sheets mimicking the terra workspace
# TODO download directly from terra

setwd('~/Dropbox (Partners HealthCare)/Projects/ms-cmmcs/2_wgs/')

source("../0_utils/utils.R")

# Load sample details -----------------------------------------------------

samples <- read_tsv("annot/Sep6_paper2022_v3_samples.txt")
participants <- read_tsv("annot/Sep6_paper2022_v3_participants.txt")
pairs <- read_tsv("annot/Sep6_paper2022_v3_pairs.txt")
wgs.coverage <- read_csv("data/coverage_cmmcs.csv")

clinicaldata <- samples |>
  mutate(Tumor_Sample_Barcode=paste0(`entity:sample_id`, "_tumor")) |>
  left_join(participants, by=c("participant"="entity:participant_id")) |>
  left_join(pairs, by=c("entity:sample_id"="case_sample", "participant")) |>
  mutate(N_cells=as.numeric(str_replace(str_extract(`entity:sample_id`, "[:alnum:]{1,}$"), "(k|K)", "000")))|>
  left_join(wgs.coverage, by=c("entity:sample_id"="ID"))

# updated values from ankit's records
clinicaldata[ clinicaldata[,'entity:sample_id']=="CTF054_BMPCs", "N_cells"] <- 17243
clinicaldata[ clinicaldata[,'entity:sample_id']=="CTF055_BMPCs", "N_cells"] <- 20783
