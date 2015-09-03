works_with_R("3.2.2", data.table="1.9.5")

public.private <- fread("public_to_private_names.csv")
setnames(public.private, c("public", "sample", "n", "private", "project"))
public.private[, cell.type := sub(".*_", "", sample)]

has.immunoseq <- fread("immunoseq_samples.txt", header=FALSE)
setnames(has.immunoseq, "immunoseqID")
has.immunoseq[, private := sub("_.*", "", immunoseqID)]
sum(!has.immunoseq$private %in% public.private$private)

public.with.immunoseq <- public.private[private %in% has.immunoseq$private, ]
table(public.with.immunoseq$cell.type)

setkey(has.immunoseq, private)
setkey(public.private, private)

public.private.immunoseq <-
  has.immunoseq[public.private, ][order(immunoseqID), ]

write.csv(public.private.immunoseq, "public.private.immunoseq.csv")
