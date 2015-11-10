works_with_R("3.2.2", data.table="1.9.6")

library(namedCapture)
library(xtable)

getenv.or <- function(env.var, default){
  env.value <- Sys.getenv(env.var)
  if(env.value == ""){
    default
  }else{
    env.value
  }
}
qsub <- getenv.or("QSUB", "qsub")

argv <- "../results/H3K27ac/PeakSegJoint.predictions.RData"

argv <- commandArgs(trailingOnly=TRUE)

predictions.RData <- normalizePath(argv[1])

load(predictions.RData)

load("~/immunoseq.calls.RData")
iseq.mat <- immunoseq.calls$rare
chr.pos.pattern <- paste0(
  "(?<chrom>.*?)",
  "_",
  "(?<pos>.*)")
iseq.meta <- str_match_named(rownames(iseq.mat), chr.pos.pattern, list(
    chrom=function(x)paste0("chr", x),
    pos=as.integer))
iseq.dt <- with(iseq.meta, {
  data.table(chrom, pos, pos2=pos, variantID=rownames(iseq.mat))
})
setkey(iseq.dt, chrom, pos, pos2)

sample.pattern <- paste0(
  "(?<cellType>[a-zA-Z]+)",
  "(?<hub>[0-9]*)",
  "/",
  "(?<person>.*?)",
  "_")
start.end.pattern <- paste0(
  "(?<chrom>chr.*?)",
  ":",
  "(?<chromStart>.*?)",
  "-",
  "(?<chromEnd>.*)")
str(all.peaks.mat)

peak.meta <- str_match_named(colnames(all.peaks.mat), start.end.pattern, list(
  chromStart=as.integer, chromEnd=as.integer))
peak.meta.dt <- data.table(peak.meta, peakID=colnames(all.peaks.mat))
setkey(peak.meta.dt, chrom, chromStart, chromEnd)

peak.sample.meta <- str_match_named(
  rownames(all.peaks.mat), sample.pattern, list(hub=as.integer))
peak.sample.meta$sample.path <- rownames(all.peaks.mat)
match.by.cellType <- split(peak.sample.meta, peak.sample.meta$cellType)

non.Input.vec <- names(match.by.cellType)[names(match.by.cellType) != "Input"]
top100.by.type <- list()
for(cellType in non.Input.vec){
  cellType.peak.meta <- match.by.cellType[[cellType]]
  ## First consider only samples with immunoseq data.
  person.vec <- intersect(colnames(iseq.mat), paste(cellType.peak.meta$person))
  subset.meta <- subset(cellType.peak.meta, person %in% person.vec)
  person.count <- table(subset.meta$person)
  subset.meta$count <- person.count[paste(subset.meta$person)]
  unique.people <-
    subset(subset.meta, !(count == 2 & grepl("McGill", sample.path)))
  stopifnot(table(unique.people$person) == 1)
  cellType.iseq <- iseq.mat[, paste(unique.people$person)]
  cellType.peaks <- all.peaks.mat[paste(unique.people$sample.path), ]
  over.dt <- foverlaps(peak.meta.dt, iseq.dt, nomatch=0L)
  ## There are some peaks that contain several rare variants.
  over.dt[, list(variants=.N), by=peakID][, table(variants)]
  over.dt[, list(peaks=.N), by=variantID][, table(peaks)]
  pval.dt.list <- list()
  counts.by.variant <- list()
  for(variant.i in 1:nrow(over.dt)){
    over.row <- over.dt[variant.i, ]
    cat(sprintf("cellType=%s Fisher test on variant %4d / %4d\n",
                cellType,
                variant.i, nrow(over.dt)))
    peak.vec <- cellType.peaks[, over.row$peakID ]
    variant.vec <- cellType.iseq[over.row$variantID, ]
    count.tab <- table(variant.vec, peak.vec)
    counts.by.variant[[over.row$variantID]] <- count.tab
    m <- as.matrix(count.tab)
    p.value <- tryCatch({
      fisher.test(m)$p.value
    }, error=function(e){
      NA_real_
    })
    is.peak <- peak.vec == 1
    pval.dt.list[[variant.i]] <- 
      data.table(over.row,
                 up=sum(is.peak),
                 down=sum(!is.peak),
                 p.value,
                 up.alt=mean(variant.vec[is.peak], na.rm=TRUE),
                 down.alt=mean(variant.vec[!is.peak], na.rm=TRUE))
  }
  pval.dt <- do.call(rbind, pval.dt.list)
  pval.dt[, diff.alt := up.alt-down.alt]
  peakID.vec <- unique(c(
    pval.dt[order(-abs(diff.alt)), peakID][1:100],
    pval.dt[order(p.value), peakID][1:100]
    ))
  setkey(pval.dt, peakID)
  peak.dt <- pval.dt[peakID.vec, {
    data.table(variants=.N,
               min.p.value=min(p.value),
               max.diff.alt=max(abs(diff.alt)))
  }, by=.(peakID, up, down)]
  prefix <- paste0(sub(":", "-", peak.dt$peakID), "/", cellType)
  peak.dt$thumb <- sprintf('
<a href="%s.png">
  <img src="%s-thumb.png" />
</a>
', prefix, prefix)
  top100.by.type[[cellType]] <- peak.dt
}
top100 <- do.call(rbind, top100.by.type)

top100.sorted <- top100[order(min.p.value),]
xt <- xtable(top100.sorted, digits=3)
data.dir <- dirname(predictions.RData)
out.dir <- file.path(data.dir, "PeakSegJoint-predictions-plots")
dir.create(out.dir, showWarnings=FALSE)
print(xt, "html",
      file=file.path(out.dir, "top100.html"),
      sanitize.text.function=identity)
save(top100, file=file.path(out.dir, "top100.RData"))

PlotPeak.R <- system.file(
  "exec", "PlotPeak.R",
  package="PeakSegJoint",
  mustWork=TRUE)

peakID.vec <- unique(top100$peakID)
for(peak.i in seq_along(peakID.vec)){
  peakID <- peakID.vec[[peak.i]]
  cat(sprintf("%4d / %4d peaks %s\n", peak.i, length(peakID.vec), peakID))
  peak.dir <- file.path(out.dir, sub(":", "-", peakID))
  peaksPlot.RData <- file.path(peak.dir, "peaksPlot.RData")
  script.txt <-
    paste0("#!/bin/bash
#PBS -l nodes=1:ppn=4
#PBS -l walltime=00:30:00
#PBS -A bws-221-ae
#PBS -m ae
#PBS -M tdhock5@gmail.com
#PBS -o ", peak.dir, ".out
#PBS -e ", peak.dir, ".err
#PBS -V                                        
#PBS -N ", peakID, "
Rscript ", PlotPeak.R, " ", predictions.RData, " ", peakID, "
Rscript ~/bin/PlotISeq.R ", peaksPlot.RData, "\n")
  script.file <- paste0(peak.dir, ".sh")
  script.dir <- dirname(script.file)
  dir.create(script.dir, showWarnings=FALSE, recursive=TRUE)
  cat(script.txt, file=script.file)
  qsub.cmd <- paste(qsub, script.file)
  system(qsub.cmd)
}
  
