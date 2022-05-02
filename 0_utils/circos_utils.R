plot.single.circos <- function(p.id, sectors, cn, svs, ctf, maf, non.synonymous, genes.of.interest, filename=NULL, major.by=5E7, axis.labels.cex = 1E-4, plotType = c('axis', 'labels')) {
  
  # extract one patient
  tx <- ctf %>% filter(ID==p.id)
  
  cn <- cn |> filter(`entity:sample_id` == p.id)
  cn <- cnv.acs.to.circos(cn, sectors = sectors)
  
  bulk.snv <- maf |> filter(sample_id==p.id)
  bulk.snv <- snv.to.circos(bulk.snv, sectors = sectors)
  # bulk.snv$col <- colorize.snv.per.classification(bulk.snv$)
  maf.coding.interest <- maf |> 
    filter(Hugo_Symbol %in% genes.of.interest & 
             Variant_Classification %in% non.synonymous)
  snv <- maf.coding.interest |> filter(sample_id==p.id)
  snv <- snv.to.circos(snv, sectors = sectors)
  
  svs.circos.ready <- rebc_to_circos_svs(svs, sectors = sectors)
  sv.rect <- svs.circos.ready %>% 
    filter(individual == p.id & class %in% c("deletion", "tandem_dup"))
  sv.link <- svs.circos.ready %>% 
    filter(individual == p.id & class %in% c("long_range", "inter_chr", "inversion"))
  
  if (!is.null(filename)) {
    if (endsWith(filename, ".pdf")) {
      pdf(filename, width = 3, height = 3)
    } else if (endsWith(filename, ".png")) {
      png(filename, width = 3, height = 3, res = 300)
    }
  }
  
  circos.clear()
  
  par(mar = c(0.5, 0.5, 0.5, 0.5))
  # par(mar = c(0, 0, 0, 0))
  circos.par("cell.padding" = c(0, 0, 0, 0), "start.degree" = 90, "track.height" = 0.05)
  # 
  # circos.initializeWithIdeogram(species = "hg19", 
  #                               major.by = 1E8, axis.labels.cex = 1E-4, 
  #                               chromosome.index = sectors, plotType = c('axis', 'labels'))
  # 
  circos.initializeWithIdeogram(species = "hg19", 
                                major.by = major.by, axis.labels.cex = axis.labels.cex, 
                                # chromosome.index = sectors, plotType = c('axis', 'ideogram', 'labels'))
                                chromosome.index = sectors, plotType =  plotType)
  
  # cn <- CTF.cnv.pairs.flat.annotated |> filter(participant == pair & tissue==itissue)
  # cn <- cnv.acs.to.circos(cn)
  
  CN.bin <- !is.null(cn) && is.data.frame(cn) && nrow(cn)>=1
  SV.rect <- !is.null(sv.rect) && is.data.frame(sv.rect) && nrow(sv.rect)>=1
  
  if(CN.bin & SV.rect) {
    bed_list = list(cn, sv.rect[,c("chr1", "pos1", "pos2", "VAF", "col")])
  } else if (CN.bin & !SV.rect) {
    bed_list = list(cn)
  } else if (!CN.bin & SV.rect){
    bed_list = list(sv.rect[,c("chr1", "pos1", "pos2", "VAF", "col")])
  } else {
    bed_list=NA
  }
  
  if(length(bed_list)>0) {
    circos.genomicTrackPlotRegion(bed_list[[1]], 
                                  # numeric.column = c(4, 4),
                                  ylim = c(0, 10), 
                                  panel.fun = function(region, value, ...) {
                                    # i = getI(...)
                                    circos.genomicRect(region, value, col=value$col, border = NA, ...)
                                    # circos.genomicRect(region, value, col=value$col, border = NA, ...)
                                  }, bg.border = NA)
  }
  
  # 
  #   if(!is.null(sv.rect) && is.data.frame(sv.rect) && nrow(sv.rect)>=1) {
  #     # circos.genomicTrackPlotRegion(cn, ylim = c(0, 10), panel.fun = function(region, value, ...) {circos.genomicRect(region, value, col=value$col_l2r, border = NA, ...)
  #     # }, bg.border = NA)
  #     circos.genomicTrackPlotRegion(sv.rect[,c("chr1", "pos1", "pos2", "col")], ylim = c(0, 10),
  #                                   panel.fun = function(region, value, ...) {
  #                                     circos.genomicRect(region, value, col=value$col, border = NA, ...)
  #                                   }, bg.border = NA)
  #   }
  #   
  
  if(!is.null(bulk.snv) && is.data.frame(bulk.snv) && nrow(bulk.snv)>=1) {
    circos.genomicTrack(bulk.snv, numeric.column = 4, ylim = c(0,1),
                        panel.fun = function(region, value, ...) {
                          # numeric.column is automatically passed to `circos.genomicPoints()`
                          circos.genomicPoints(region, value, cex = .3, pch = 16, col = value$Color, ...)
                        },bg.border = NA)
  }
  
  if(!is.null(snv) && is.data.frame(snv) && nrow(snv)>=1) {
    # circos.genomicTrack(snv, numeric.column = 4, ylim = c(0,1),
    #                     panel.fun = function(region, value, ...) {
    #                       # numeric.column is automatically passed to `circos.genomicPoints()`
    #                       circos.genomicPoints(region, value, ...)
    #                     },bg.border = NA)
    circos.genomicLabels(snv, labels.column = 5, side = "inside", cex = 0.5, labels_height = 0.5, track.margin = c(0,0))
  }
  
  
  if(!is.null(sv.link) && is.data.frame(sv.link) && nrow(sv.link)>=1 ) {
    circos.genomicLink(sv.link[, c("chr1", "pos1", "pos1")], 
                       sv.link[, c("chr2", "pos2", "pos2")],
                       col = sv.link$col,
                       rou = .7, lwd = .5)
  }
  
  
  if(!is.null(tx) && is.data.frame(tx) && nrow(tx)>=1 ) {
    circos.genomicLink(tx[, c("chr1", "start1", "end1")], tx[, c("chr2", "start2", "end2")], col =
                         "darkred", rou = .75, lwd = 2, h.ratio=.3)
  }
  
  if(!is.null(filename)) {
    if(endsWith(filename, ".pdf") | endsWith(filename, ".png")) dev.off()
  } 
  
}
