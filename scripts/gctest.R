works_with_R("3.2.2", data.table="1.9.6",
             "tdhock/namedCapture@5178fcff3ab3fa22723a6ed274b12e5c02da366f")

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
for(bedGraph.file.i in seq_along(bedGraph.file.vec)){
  bedGraph.file <- bedGraph.file.vec[[bedGraph.file.i]]
  cat(sprintf("%4d / %4d %s\n", bedGraph.file.i, length(bedGraph.file.vec),
              bedGraph.file))
  file.coverage <- fread(bedGraph.file)
  setnames(file.coverage, c("chrom", "chromStart", "chromEnd", "count"))
  match.row <- data.frame(match.df[bedGraph.file.i, ])
  coverage.list[[bedGraph.file]] <- data.table(match.row, file.coverage)
}
gctest <- do.call(rbind, coverage.list)

stop("TODO add GC content track, compute binned data")

save(gctest, file="gctest.RData")
