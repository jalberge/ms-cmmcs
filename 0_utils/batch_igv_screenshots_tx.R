# This script takes bams and rebc sv output files as arguments and generate igv screenshots

# Params ------------------------------------------------------------------

setwd('~/Dropbox (Partners HealthCare)/Projects/ms-cmmcs/0_utils/')

source("../2_wgs/0_annotate_samples.R")

window.size=2000
max.window=25000
Plot.Normal <- FALSE

# Participant data --------------------------------------------------------

Pairs <- clinicaldata %>% 
  filter(tissue %in% c("CMMCs", "BMPCs")) %>% 
  pull(`entity:sample_id`)
sapply(paste0("../igvsessions/", Pairs), dir.create)

# Preprocess for structural variants --------------------------------------

sv <- read_tsv("../data/Paper_v3_revisions_SVs.txt")

positions.and.filenames <- sv %>%
  rowwise() %>%
  mutate(position1=paste0(chr1,":", pos1-window.size/2, "-", pos1+window.size/2), 
         position2=paste0(chr2,":", pos2-window.size/2, "-", pos2+window.size/2),
         ID_1=paste(individual, chr1, pos1, chr2, pos2, "left", sep="_"),
         ID_2=paste(individual, chr1, pos1, chr2, pos2, "right", sep="_"),
         filename.1=paste0(class, "_", ID_1, ".svg"),
         filename.2=paste0(class, "_", ID_2, ".svg")) %>%
  as.data.frame()

Pairs <- c("CTF048_BMPCs_9906", "CTF048_CMMCs_266", "CTF049_CMMCs_1513")

# Write IGV bash script for cohort structural variant ---------------------

script.file <- "../igvsessions/Oct17_igv_svs.script"

file.remove(script.file)

# the trick (thanks to Andrew) is to create a single script and call the instruction new for each new bam
# see https://github.com/igvteam/igv/wiki/Batch-commands

for(TumorID in Pairs) {
  
  Participant <- clinicaldata %>% filter(`entity:sample_id`==TumorID)  %>% pull(participant)
  # get basic info for igv
  # TumorID <- clinicaldata %>% filter(participant==Participant & tissue==Tissue) %>% pull(`entity:sample_id`)
  TumorBAM <- clinicaldata %>% filter(`entity:sample_id`==TumorID) %>% pull(analysis_ready_bam)
  TumorBAI <- clinicaldata %>% filter(`entity:sample_id`==TumorID) %>% pull(analysis_ready_bam_index)
  if(Plot.Normal) {
    NormalID <- clinicaldata %>% filter(participant==Participant & tissue=="WBCs") %>% pull(`entity:sample_id`)
    NormalBAM <- clinicaldata %>% filter(participant==Participant & tissue=="WBCs") %>% pull(analysis_ready_bam)
    NormalBAI <- clinicaldata %>% filter(participant==Participant & tissue=="WBCs") %>% pull(analysis_ready_bam_index)
  }
  
  write("new", script.file, append = TRUE)
  # init script
  write("genome hg19", script.file, append = TRUE)
  write("maxPanelHeight 3000", script.file, append = TRUE)
  
  write(paste0("load ", TumorBAM, " index=", TumorBAI, " name=", TumorID, "_tumor"), script.file, append = TRUE)
  write("colorBy UNEXPECTED_PAIR", script.file, append = TRUE)
  write("group MATE_CHROMOSOME", script.file, append = TRUE)
  write("collapse", script.file, append = TRUE)
  
  if(Plot.Normal) {
    write(paste0("load ", NormalBAM, " index=", NormalBAI), script.file, append = TRUE)
    write("colorBy UNEXPECTED_PAIR", script.file, append = TRUE)
    write("group MATE_CHROMOSOME", script.file, append = TRUE)
    write("collapse", script.file, append = TRUE)
  }
  
  write(paste0("snapshotDirectory ", " ../igvsessions/", TumorID, "/snapshots"), script.file, append = TRUE)
  
  # Translocations will require to group by chromosome, and two screenshots per SV
  sv.df <- positions.and.filenames %>% filter(individual == TumorID & class == "inter_chr" )
  
  if(nrow(sv.df>0)) {
    for(i in 1:nrow(sv.df)) {
      write(paste0("goto ", sv.df[i, "position1"]), script.file, append = TRUE)
      write(paste0("snapshot ", sv.df[i, "filename.1"]), script.file, append = TRUE)
      
      write(paste0("goto ", sv.df[i, "position2"]), script.file, append = TRUE)
      write(paste0("snapshot ", sv.df[i, "filename.2"]), script.file, append = TRUE)
    }
  }
  
  # other SVs not requiring mate chromosomes
  write("group PAIR_ORIENTATION", script.file, append = TRUE)
  write("colorBy UNEXPECTED_PAIR", script.file, append = TRUE)
  write("viewaspairs", script.file, append = TRUE)
  write("sort INSERTSIZE", script.file, append = TRUE)
  
  # local structural variants
  dels.df <- positions.and.filenames %>% filter(individual == TumorID  & class %in% c("deletion", "inversion", "long_range", "tandem_dup"))
  
  if(nrow(dels.df>0)) {
    for(i in 1:nrow(dels.df)) {
      
      if( dels.df[i, "span"] < max.window) {
        
        # quick hack
        write(paste0("goto ", dels.df[i, "chr1"], ":", dels.df[i, "pos1"]-500, "-", dels.df[i, "pos2"]+500 ), script.file, append = TRUE)
        write(paste0("snapshot ", dels.df[i, "filename.1"]), script.file, append = TRUE)
        
      } else {
        write(paste0("goto ", dels.df[i, "position1"]), script.file, append = TRUE)
        write(paste0("snapshot ", dels.df[i, "filename.1"]), script.file, append = TRUE)
        
        write(paste0("goto ", dels.df[i, "position2"]), script.file, append = TRUE)
        write(paste0("snapshot ", dels.df[i, "filename.2"]), script.file, append = TRUE)
      }
      
    }
  }
  
  write(paste0("saveSession ", "../igvsessions/", TumorID, "/", TumorID, ".xml"), script.file, append = TRUE)
  
}
write("exit", script.file, append = TRUE)


# Run script --------------------------------------------------------------

system(paste0("zsh ~/IGV_2.12.0/igv.sh -b ", script.file))

