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
iseq.dt <- with(iseq.meta, data.table(chrom, pos, pos2=pos))
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
peak.meta.dt <- data.table(peak.meta)
setkey(peak.meta.dt, chrom, chromStart, chromEnd)

peak.sample.meta <- str_match_named(
  rownames(all.peaks.mat), sample.pattern, list(hub=as.integer))
peak.sample.meta$sample.path <- rownames(all.peaks.mat)
match.by.cellType <- split(peak.sample.meta, peak.sample.meta$cellType)

for(cellType in c("mono", "neutro")){
  cellType.peak.meta <- match.by.cellType[[cellType]]
  ## First consider only samples with immunoseq data.
  person.vec <- intersect(iseq.person.vec, paste(cellType.peak.meta$person))
  subset.meta <- subset(cellType.peak.meta, person %in% person.vec)
  person.count <- table(subset.meta$person)
  subset.meta$count <- person.count[paste(subset.meta$person)]
  unique.people <-
    subset(subset.meta, !(count == 2 & grepl("McGill", sample.path)))
  stopifnot(table(unique.people$person) == 1)
  cellType.iseq <- iseq.mat[, paste(unique.people$person)]
  cellType.peaks <- all.peaks.mat[paste(unique.people$sample.path), ]
  stop("TODO compute which iseq SNPs fall in peaks")
}

