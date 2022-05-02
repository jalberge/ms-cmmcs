---
title: "Heatmap generation for CMMC only samples"
author: "Ankit Dutta and JB Alberge"
date: February 16, 2022
output: 
  html_document:
    toc: true
    toc_depth: 2
    df_print: kable
    keep_md: true
---

Saving  CNV heatmap and paired tx data here for cmmc only cohort

# Init session


```r
source("0_annotate_samples.R")
source("../0_utils/circos_utils.R")

dir.create("figures")

# CNV.df <- read_tsv("data/Feb7_reprocessed.cnv.seg")
CNV.df <- read_tsv("data/March9_reprocessed.cnv.prod5.seg")
# maf <- read_maf("../data/Feb14.CMMCsPaper2022.aggregated.maf")
maf <- read_maf("../data/Paper2022_v2_filtered_crux.maf")
ctf.tx <- read_tsv("../data/FullPaper2022.aggregated.ctf.whitelist.feb28.txt")
svs <- read_xlsx("../data/CMMCsPaper2022_ConsensusSV_filtered_results_whitelist_Mar3.xlsx") 

contigs <- factor(1:22, levels=c(1:22))
sectors <- paste0("chr", contigs) # circlize uses UCSC chr notation. TODO update for hg38
```

# Post proc CNV


```r
contigs <- factor(1:22, levels=c(1:22))

CNV.df <- CNV.df %>% 
  inner_join(clinicaldata, by=c("entity:sample_id")) %>% 
  # filter(isRef==TRUE)
  filter(Chromosome %in% contigs) %>%
  mutate(Chromosome=factor(Chromosome, levels=contigs),
         Norm_Segment_Mean = pmax(0, pmin(Segment_Mean, 4)),
         YMIN=ifelse(tissue=="BMPCs", 0, 1), 
         YMAX =ifelse(tissue=="BMPCs", 1, 2))
```

# Post proc SNV indels


```r
maf <- maf |>
  mutate(sample_id=str_remove(Tumor_Sample_Barcode, "_tumor")) |>
  left_join(clinicaldata, by=c("sample_id"="entity:sample_id"))
```

# Post proc catch de fish


```r
ctf.tx <- ctf.tx %>% 
  filter(Keep=="TRUE") %>%
  left_join(clinicaldata, by=c("ID"="entity:sample_id")) %>%
  relocate(CHR.IG, IG_START, CHR.ONCO, PARTNER_START) %>% 
  mutate(chr1 = paste0("chr", CHR.IG), 
         chr2 = paste0("chr", CHR.ONCO),
         start1 = IG_START,
         end1 = IG_START,
         start2 = PARTNER_START,
         end2 = PARTNER_START) %>% # TODO update for hg38 chr
  mutate(chr1 = case_when(chr1=="chr23" ~ "ChrX",
                          chr1=="chr24" ~ "ChrY",
                          TRUE ~ chr1)) %>%
    group_by(ID, IG, ONCO) %>%
    slice_head(n=1) %>%
    filter(chr1 %in% sectors & chr2 %in% sectors)
```

# Post process SVs


```r
svs <- svs %>% filter(keep==1) %>% left_join(clinicaldata, by=c("individual"="entity:sample_id"))
```

# Re-format all for circlize

## CTF025 and CTF031 examples


```r
p.id="CTF025_CMMCs_198"

plot.single.circos(p.id = p.id, sectors = sectors, ctf = ctf.tx, cn=CNV.df, svs=svs, maf=maf, non.synonymous = non.synonymous, genes.of.interest = genes.of.interest, filename = "figures/" %+% p.id)
plot.single.circos(p.id = p.id, sectors = sectors, ctf = ctf.tx, cn=CNV.df, svs=svs, maf=maf, non.synonymous = non.synonymous, genes.of.interest = genes.of.interest)
```

![](6_circos_plots_snv_sv_cnv_files/figure-html/unnamed-chunk-6-1.png)<!-- -->




```r
p.id="CTF031_CMMCs_41500"

plot.single.circos(p.id = p.id, sectors = sectors, ctf = ctf.tx, cn=CNV.df, svs=svs, maf=maf, non.synonymous = non.synonymous, genes.of.interest = genes.of.interest, filename = "figures/" %+% p.id)
plot.single.circos(p.id = p.id, sectors = sectors, ctf = ctf.tx, cn=CNV.df, svs=svs, maf=maf, non.synonymous = non.synonymous, genes.of.interest = genes.of.interest)
```

![](6_circos_plots_snv_sv_cnv_files/figure-html/unnamed-chunk-7-1.png)<!-- -->

## Chromothripsis of chr 3 in CTF033


```r
p.id="CTF033_CMMCs_755"

p1 <- plot.single.circos(p.id = p.id, sectors = sectors, ctf = ctf.tx, cn=CNV.df, svs=svs, maf=maf, non.synonymous = non.synonymous, genes.of.interest = genes.of.interest)
```

![](6_circos_plots_snv_sv_cnv_files/figure-html/unnamed-chunk-8-1.png)<!-- -->

```r
p2 <- plot.single.circos(p.id = p.id, sectors = "chr3", ctf = ctf.tx, cn=CNV.df, svs=svs, maf=maf, non.synonymous = non.synonymous, genes.of.interest = genes.of.interest, axis.labels.cex = .5, plotType = c('axis', 'ideogram', 'labels'))
```

![](6_circos_plots_snv_sv_cnv_files/figure-html/unnamed-chunk-8-2.png)<!-- -->

```r
plot_grid(p1, p2)
```

![](6_circos_plots_snv_sv_cnv_files/figure-html/unnamed-chunk-8-3.png)<!-- -->

```r
plot.single.circos(p.id = p.id, sectors = sectors, ctf = ctf.tx, cn=CNV.df, svs=svs, maf=maf, non.synonymous = non.synonymous, genes.of.interest = genes.of.interest, filename = "figures/" %+% p.id %+% ".pdf")
```

```
## quartz_off_screen 
##                 2
```

```r
plot.single.circos(p.id = p.id, sectors = "chr3", ctf = ctf.tx, cn=CNV.df, svs=svs, maf=maf, non.synonymous = non.synonymous, genes.of.interest = genes.of.interest, filename = "figures/" %+% p.id %+% "_zoom.pdf", axis.labels.cex = .5, plotType = c('axis', 'ideogram', 'labels'))
```

```
## quartz_off_screen 
##                 2
```
## Chromoplexy of chr 7, 8, 18 in CTF034

### In CMMCs:


```r
p.id="CTF034_CMMCs_654"
sectors.ctf034 <- "chr" %+% c(7, 8, 18)

plot.single.circos(p.id = p.id, sectors = sectors.ctf034, ctf = ctf.tx, cn=CNV.df, svs=svs, maf=maf, non.synonymous = non.synonymous, genes.of.interest = genes.of.interest, axis.labels.cex = .5, plotType = c('axis', 'ideogram', 'labels'))
```

![](6_circos_plots_snv_sv_cnv_files/figure-html/unnamed-chunk-9-1.png)<!-- -->

```r
plot.single.circos(p.id = p.id, sectors = sectors.ctf034, ctf = ctf.tx, cn=CNV.df, svs=svs, maf=maf, non.synonymous = non.synonymous, genes.of.interest = genes.of.interest, axis.labels.cex = .5, plotType = c('axis', 'ideogram', 'labels'), filename = "figures/" %+% p.id %+% "_zoom.pdf")
```

```
## quartz_off_screen 
##                 2
```

### In BMPCs:


```r
p.id="CTF034_BMPCs_16800"
sectors.ctf034 <- "chr" %+% c(7, 8, 18)

plot.single.circos(p.id = p.id, sectors = sectors.ctf034, ctf = ctf.tx, cn=CNV.df, svs=svs, maf=maf, non.synonymous = non.synonymous, genes.of.interest = genes.of.interest, axis.labels.cex = .5, plotType = c('axis', 'ideogram', 'labels'))
```

![](6_circos_plots_snv_sv_cnv_files/figure-html/unnamed-chunk-10-1.png)<!-- -->

```r
plot.single.circos(p.id = p.id, sectors = sectors.ctf034, ctf = ctf.tx, cn=CNV.df, svs=svs, maf=maf, non.synonymous = non.synonymous, genes.of.interest = genes.of.interest, axis.labels.cex = .5, plotType = c('axis', 'ideogram', 'labels'), filename = "figures/" %+% p.id %+% "_zoom.pdf")
```

```
## quartz_off_screen 
##                 2
```


## CTF047 example


```r
p.id="CTF047_BMPCs"

plot.single.circos(p.id = p.id, sectors = sectors, ctf = ctf.tx, cn=CNV.df, svs=svs, maf=maf, non.synonymous = non.synonymous, genes.of.interest = genes.of.interest, filename = "figures/CTF047_PBMCs.pdf")
```

```
## quartz_off_screen 
##                 2
```

```r
plot.single.circos(p.id = p.id, sectors = sectors, ctf = ctf.tx, cn=CNV.df, svs=svs, maf=maf, non.synonymous = non.synonymous, genes.of.interest = genes.of.interest)
```

![](6_circos_plots_snv_sv_cnv_files/figure-html/unnamed-chunk-11-1.png)<!-- -->

```r
p.id="CTF047_CMMCs_75k"

plot.single.circos(p.id = p.id, sectors = sectors, ctf = ctf.tx, cn=CNV.df, svs=svs, maf=maf, non.synonymous = non.synonymous, genes.of.interest = genes.of.interest, filename = "figures/" %+% p.id %+% ".pdf")
```

```
## quartz_off_screen 
##                 2
```

```r
plot.single.circos(p.id = p.id, sectors = sectors, ctf = ctf.tx, cn=CNV.df, svs=svs, maf=maf, non.synonymous = non.synonymous, genes.of.interest = genes.of.interest)
```

![](6_circos_plots_snv_sv_cnv_files/figure-html/unnamed-chunk-11-2.png)<!-- -->


# Session info


```r
sessionInfo()
```

```
R version 4.1.1 (2021-08-10)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS Big Sur 10.16

Matrix products: default
BLAS:   /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRblas.0.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] circlize_0.4.13    scales_1.1.1       ggpubr_0.4.0       readxl_1.3.1      
 [5] cowplot_1.1.1      maftools_2.8.05    RColorBrewer_1.1-2 kableExtra_1.3.4  
 [9] rstatix_0.7.0      data.table_1.14.0  forcats_0.5.1      stringr_1.4.0     
[13] dplyr_1.0.7        purrr_0.3.4        readr_2.0.1        tidyr_1.1.3       
[17] tibble_3.1.3       ggplot2_3.3.5      tidyverse_1.3.1   

loaded via a namespace (and not attached):
 [1] fs_1.5.0            bit64_4.0.5         lubridate_1.7.10   
 [4] webshot_0.5.2       httr_1.4.2          tools_4.1.1        
 [7] backports_1.2.1     bslib_0.2.5.1       utf8_1.2.2         
[10] R6_2.5.0            DBI_1.1.1           colorspace_2.0-2   
[13] withr_2.4.2         tidyselect_1.1.1    bit_4.0.4          
[16] curl_4.3.2          compiler_4.1.1      cli_3.1.0          
[19] rvest_1.0.1         xml2_1.3.2          sass_0.4.0         
[22] systemfonts_1.0.3   digest_0.6.27       foreign_0.8-81     
[25] rmarkdown_2.10      svglite_2.0.0       rio_0.5.27         
[28] pkgconfig_2.0.3     htmltools_0.5.1.1   highr_0.9          
[31] dbplyr_2.1.1        rlang_0.4.11        GlobalOptions_0.1.2
[34] rstudioapi_0.13     shape_1.4.6         jquerylib_0.1.4    
[37] generics_0.1.0      jsonlite_1.7.2      vroom_1.5.4        
[40] zip_2.2.0           car_3.0-11          magrittr_2.0.1     
[43] Matrix_1.3-4        Rcpp_1.0.7          munsell_0.5.0      
[46] fansi_0.5.0         abind_1.4-5         lifecycle_1.0.0    
[49] stringi_1.7.3       yaml_2.2.1          carData_3.0-4      
[52] grid_4.1.1          parallel_4.1.1      crayon_1.4.1       
[55] lattice_0.20-44     haven_2.4.3         splines_4.1.1      
[58] hms_1.1.0           knitr_1.33          pillar_1.6.2       
[61] ggsignif_0.6.2      reprex_2.0.1        glue_1.4.2         
[64] evaluate_0.14       modelr_0.1.8        vctrs_0.3.8        
[67] tzdb_0.1.2          cellranger_1.1.0    gtable_0.3.0       
[70] assertthat_0.2.1    xfun_0.25           openxlsx_4.2.4     
[73] broom_0.7.9         survival_3.2-12     viridisLite_0.4.0  
[76] ellipsis_0.3.2     
```
