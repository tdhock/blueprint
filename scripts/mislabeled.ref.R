works_with_R("3.2.2",
             "tdhock/namedCapture@a31a38be12d4dad4aec6257a47c0e2307679589c",
             data.table="1.9.6")

load("hubs.nodup.RData")

pattern <- paste0(
  "(?<private>.*?)",
  " [|] *",
  "(?<long_type>.*?)",
  " [|] ",
  "(?<where>.*?)",
  " [|] ",
  "(?<experiment>[^ ]+)",
  " ",
  "(?<problem>.*)")

long2short <- c(
  mature_neutrophil="GR",
  Tcell="nTC",
  monocyte="Mono",
  "CD14-positive_CD16-negative_classical_monocyte"="Mono"
  )

samples.by.hub <- split(hubs.nodup, hubs.nodup$hub)
bad.txt.vec <- Sys.glob("../mislabeled_ref/*.txt")
mislabeled.ref.list <- list()
for(bad.txt in bad.txt.vec){
  bad.lines <- readLines(bad.txt)
  hub <- sub(".txt", "", basename(bad.txt))
  has.problems <- any(grepl("[|]", bad.lines))
  mislabeled.ref.list[[hub]] <- if(has.problems){
    hub.dt <- samples.by.hub[[hub]]
    bad.mat <- str_match_named(bad.lines[bad.lines!=""], pattern)
    bad.dt <- data.table(bad.mat)
    stopifnot(bad.dt$long_type %in% names(long2short))
    bad.dt[, cell.type := long2short[paste(long_type)]]
    setkey(hub.dt, private, cell.type, where, experiment)
    setkey(bad.dt, private, cell.type, where, experiment)
    merge.dt <- hub.dt[bad.dt]
    stopifnot(nrow(merge.dt) == nrow(bad.dt))
    stopifnot(!is.na(merge.dt))
    merge.dt
  }else{
    data.table()
  }
}

is.missing <- ! names(samples.by.hub) %in% names(mislabeled.ref.list)
if(any(is.missing)){
  print(names(samples.by.hub)[is.missing])
}

mislabeled <- data.frame(do.call(rbind, mislabeled.ref.list))

with(mislabeled, table(cell.type, where, experiment))

mislabeled$looks.like <- with(mislabeled, {
  ifelse(grepl("looks like", problem),
         sub("looks like ", "", problem),
         "Input")
})

with(mislabeled, table(problem, looks.like))
rownames(mislabeled) <- mislabeled$bigwig

hubs.nodup[, looks.like := mislabeled[paste(bigwig), "looks.like"]]
hubs.nodup[is.na(looks.like), looks.like := "OK"]
hubs.nodup[, url := sub(".*blueprint",
                   "http://hubs.hpc.mcgill.ca/~thocking", bigwig)]

hub.prefix <- "http://hubs.hpc.mcgill.ca/~thocking/hubs_nodup"
labeled.bigwigs <- hubs.nodup[, {
  data.frame(donor=private, cell.type, center=where,
             experiment, looks.like,
             hub=sprintf("%s/%s/hub.txt", hub.prefix, hub),
             url)
}]
table(labeled.bigwigs$looks.like)
with(labeled.bigwigs, table(experiment, looks.like))
stopifnot(sum(labeled.bigwigs$looks.like != "OK") == nrow(mislabeled))
write.csv(labeled.bigwigs, "labeled.bigwigs.csv")
save(labeled.bigwigs, file="labeled.bigwigs.RData")

write.csv(mislabeled, "mislabeled.csv")
save(mislabeled, file="mislabeled.RData")
