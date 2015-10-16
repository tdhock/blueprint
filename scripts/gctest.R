works_with_R("3.2.2", data.table="1.9.6",
             "tdhock/PeakSegJoint@2b1f266808136631bb7380c04ff257b9f3c24ef9",
             "tdhock/namedCapture@5178fcff3ab3fa22723a6ed274b12e5c02da366f")

## From http://genome.ucsc.edu/goldenpath/help/wiggle.html

## For example, this variableStep specification:

## variableStep chrom=chr2
## 300701 12.5
## 300702 12.5
## 300703 12.5
## 300704 12.5
## 300705 12.5

## is equivalent to:

## variableStep chrom=chr2 span=5
## 300701 12.5

## Our data:

## chr1:32,400,001-32,400,005 AGCTT 40%

## chr1:32,400,006-32,400,010 CATTC 40%

## chr1:32,400,016-32,400,020 AGACC 60%

pattern <- paste0(
  "(?<person>[^/]+)",
  "_",
  "(?<cellType>.*?)",
  "_",
  "(?<experiment>.*?)",
  "_",
  "(?<center>[^.]+)")
bedGraph.file.vec <- Sys.glob("../gctest/*.bedGraph")
match.mat <- str_match_named(bedGraph.file.vec, pattern)
match.df <- data.frame(match.mat)

coverage.list <- list()
coverageBins.list <- list()
for(bedGraph.file.i in seq_along(bedGraph.file.vec)){
  bedGraph.file <- bedGraph.file.vec[[bedGraph.file.i]]
  cat(sprintf("%4d / %4d %s\n", bedGraph.file.i, length(bedGraph.file.vec),
              bedGraph.file))
  file.coverage <- fread(bedGraph.file)
  setnames(file.coverage, c("chrom", "chromStart", "chromEnd", "count"))
  max.count <- max(file.coverage$count)
  file.coverage[, norm := count/max.count]
  file.bins <-
    binSum(file.coverage,
           bin.chromStart=32400000L,
           bin.size=100L,
           n.bins=255L)
  file.bins$norm <- file.bins$mean / max.count
  match.row <- data.frame(match.df[bedGraph.file.i, ])
  coverage.list[[bedGraph.file]] <- data.table(match.row, file.coverage)
  coverageBins.list[[bedGraph.file]] <- data.table(match.row, file.bins)
}
gctest <- do.call(rbind, coverage.list)
coverageBins <- do.call(rbind, coverageBins.list)

if(!file.exists("hg19.gc5Base.txt.gz"))download.file("http://hgdownload.cse.ucsc.edu/goldenPath/hg19/gc5Base/hg19.gc5Base.txt.gz", "hg19.gc5Base.txt.gz")

message("All data in hg19.gc5Base.txt")

gc5Base <- fread("chr1_32400000_32425500.wig")
setnames(gc5Base, c("position", "percent"))
gc5Base[, chromStart := position-1L]
gc5Base[, chromEnd := position+4L]
gc5Base[, count := as.integer(percent/100 * 5) ]

gcBins <-
  binSum(gc5Base,
         bin.chromStart=32400000L,
         bin.size=100L,
         n.bins=255L)

save(gctest, gc5Base,
     coverageBins, gcBins,
     file="gctest.RData")
