
setwd('~/Dropbox (Partners HealthCare)/Projects/ms-cmmcs/2_wgs/')

source("../0_utils/utils.R")

# Load sample details -----------------------------------------------------

samples <- read_tsv("annot/Feb28_paper2022_v2_samples.txt")
participants <- read_tsv("annot/Feb28_paper2022_v2_participants.tsv")
pairs <- read_tsv("annot/Feb28_paper2022_v2_pairs.tsv")
# purity <- read_xlsx("data/purity_and_contam.xlsx", sheet="Sheet2")
wgs.coverage <- read_csv("data/coverage_cmmcs.csv")
# wgs.coverage <- read_csv("data/wgs_metrics_df.csv")
wgs.coverage <- wgs.coverage %>% select(-condition, -tissue)

clinicaldata <- samples |>
  mutate(Tumor_Sample_Barcode=paste0(`entity:sample_id`, "_tumor")) |>
  left_join(participants, by=c("participant"="entity:participant_id")) |>
  left_join(pairs, by=c("entity:sample_id"="case_sample", "participant")) |>
  mutate(N_cells=as.numeric(str_replace(str_extract(`entity:sample_id`, "[:alnum:]{1,}$"), "(k|K)", "000")))|>
  left_join(wgs.coverage, by=c("entity:sample_id"="ID"))
  # left_join(purity, by=c("entity:sample_id"="entity:pair_id"))
