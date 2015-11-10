works_with_R("3.2.2",
             data.table="1.9.6",
             "tdhock/namedCapture@a31a38be12d4dad4aec6257a47c0e2307679589c")

argv <- "~/PeakSegJoint-predictions-plots/chr5:156569171-156570800/peaksPlot.RData"

peaksPlot.RData <- normalizePath(argv[1])
out.dir <- dirname(peaksPlot.RData)

load(peaksPlot.RData)
load("immunoseq.RData")

peaksPlot$coverage[, donor := sub("_.*", "", sample.id)]
peaksPlot$coverage[, cell.type := sub("[0-9]*$", "", sample.group)]

peaksPlot$peaks[, donor := sub("_.*", "", sample.id)]
peaksPlot$peaks[, cell.type := sub("[0-9]*$", "", sample.group)]

plot.donor.vec <- unique(peaksPlot$coverage$donor)

chunk <- peaksPlot$peaks[1,]
setkey(chunk, chrom, zoomStart, zoomEnd)

genotype.code <- c(
  "0/0"="homozygousRef",
  "0/1"="heterozygous",
  "1/1"="homozygousAlt")

pattern <- paste0(
      "(?<GT>.*?",
      "(?<allele1>[0-9])",
      "/",
      "(?<allele2>[0-9])",
      ")",
      ":",
      "(?<AD>.*?)",
      ":",
      "(?<DP>.*?)",
      ":",
      "(?<GQ>.*?)",
      ":",
      "(?<PL>.*?)",
      "$")

allele.code <- c("0"="REF", "1"="ALT")

variants.by.pos <- list()
for(variant.type in c("rare")){
  type.variants <- immunoseq[[variant.type]]
  type.variants[, POS0 := POS]
  type.variants[, chrom := paste0("chr", CHROM)]
  setkey(type.variants, chrom, POS, POS0)
  over.dt <- foverlaps(type.variants, chunk, nomatch=0L)
  sample.col.vec <- grep("WB", names(over.dt), value=TRUE)
  names(sample.col.vec) <- sub("_.*", "", sample.col.vec)
  donor.vec <- intersect(names(sample.col.vec), plot.donor.vec)
  for(variant.i in 1:nrow(over.dt)){
    one.variant <- over.dt[variant.i, ]
    one.variant.calls <- one.variant[, sample.col.vec[donor.vec], with=FALSE]
    match.mat <- str_match_named(paste(one.variant.calls), pattern)
    variant.code <- one.variant[, c(REF=REF, ALT=ALT)]
    one.variant.dt <- data.table(
      donor=donor.vec,
      allele1code=allele.code[match.mat[, "allele1"]],
      allele2code=allele.code[match.mat[, "allele2"]]
      )
    one.variant.dt[, `:=`(
      allele1=variant.code[allele1code],
      allele2=variant.code[allele2code]
      )]
    variants.by.pos[[one.variant$chr_pos]] <-
      data.table(chrom=one.variant$chrom,
                 position=one.variant$POS,
                 one.variant.dt)
  }
}
variants <- do.call(rbind, variants.by.pos)

variants.RData <- sub("peaksPlot", "variants", peaksPlot.RData)
save(variants, file=variants.RData)

type.vec <- unique(peaksPlot$coverage$cell.type)
setkey(peaksPlot$coverage, cell.type)
setkey(peaksPlot$peaks, cell.type)
for(cell.type in type.vec[type.vec != "Input"]){
  type.coverage <- peaksPlot$coverage[cell.type]
  type.peaks <- peaksPlot$peaks[cell.type]
  donors.with.coverage <- unique(type.coverage$donor)
  type.variants <- variants[donor %in% donors.with.coverage,]
  n.profiles <- length(donors.with.coverage)
  alt.counts <- variants[, {
    list(nALT=sum(c(allele1code, allele2code)=="ALT"))
  }, by=donor]
  setkey(alt.counts, donor)
  peak.counts <- type.peaks[, {
    list(nPeaks=.N)
  }, by=donor]
  setkey(peak.counts, donor)
  type.donors <- data.table(donor=donors.with.coverage)
  donor.counts <- peak.counts[alt.counts[type.donors]]
  stopifnot(nrow(donor.counts) == nrow(type.donors))
  donor.counts.sorted <- donor.counts[order(nALT, nPeaks, decreasing=TRUE), ]
  dfac <- function(donor){
    factor(donor, donor.counts.sorted$donor)
  }
  type.coverage[, donor.fac := dfac(donor)]
  type.peaks[, donor.fac := dfac(donor)]
  type.variants[, donor.fac := dfac(donor)]
  
  gg <- 
    ggplot()+
      ggtitle(paste("cell type =", cell.type))+
      scale_y_continuous("aligned read coverage",
                         breaks=function(limits){
                           floor(limits[2])
                         })+
      scale_x_continuous(paste("position on", chunk$chrom,
                               "(kilo bases = kb)"))+
      coord_cartesian(xlim=chunk[, c(zoomStart, zoomEnd)/1e3])+
      theme_bw()+
      theme(panel.margin=grid::unit(0, "cm"))+
      facet_grid(donor.fac ~ ., scales="free")+
      geom_rect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
                    ymin=0, ymax=count),
                color="grey50",
                data=type.coverage)+
      geom_segment(aes(chromStart/1e3, 0,
                       xend=chromEnd/1e3, yend=0),
                   data=type.peaks,
                   size=4,
                   color="deepskyblue")+
      geom_text(aes(position/1e3, 0, label=allele1, color=allele1code),
                data=type.variants[!is.na(allele1code),],
                hjust=1,
                vjust=0)+
      geom_text(aes(position/1e3, 0, label=allele2, color=allele2code),
                data=type.variants[!is.na(allele2code),],
                hjust=0,
                vjust=0)+
      scale_color_manual("allele", values=c(REF="black", ALT="red"))
  height.pixels <- (n.profiles+1)*35
  png.name <- file.path(out.dir, paste0(cell.type, ".png"))
  print(png.name)
  png(png.name, width=1000, h=height.pixels, units="px")
  print(gg)
  dev.off()

}

