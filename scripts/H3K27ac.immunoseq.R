works_with_R("3.2.2", data.table="1.9.6",
             "tdhock/namedCapture@a31a38be12d4dad4aec6257a47c0e2307679589c")

load("immunoseq.calls.RData")
load("../results/H3K27ac/PeakSegJoint.predictions.RData")

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

H3K27ac.immunoseq <- list()
for(cellType in c("mono", "neutro")){
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
  str(iseq.dt)
  str(cellType.iseq)
  cellType.peaks <- all.peaks.mat[paste(unique.people$sample.path), ]
  str(cellType.peaks)
  str(peak.meta.dt)
  over.dt <- foverlaps(peak.meta.dt, iseq.dt, nomatch=0L)
  ## There are some peaks that contain several rare variants.
  over.dt[, list(variants=.N), by=peakID][, table(variants)]
  over.dt[, list(peaks=.N), by=variantID][, table(peaks)]
  pval.dt.list <- list()
  counts.by.variant <- list()
  for(variant.i in 1:nrow(over.dt)){
    over.row <- over.dt[variant.i, ]
    cat(sprintf("%4d / %4d variants\n", variant.i, nrow(over.dt)))
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
  pval.dt[order(p.value), ][1:20,]
  pval.dt[down==1,][order(p.value),][1:10,]
  pval.dt[down==2,][order(p.value),][1:10,]
  H3K27ac.immunoseq[[cellType]] <-
    list(counts.by.variant=counts.by.variant,
         p.values=pval.dt)
}

save(H3K27ac.immunoseq, file="H3K27ac.immunoseq.RData")
