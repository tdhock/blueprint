works_with_R("3.2.2", data.table="1.9.5")

public.private <- fread("../secret/public_to_private_names.csv")
setnames(public.private, c("public", "sample", "n", "private", "project"))
public.private[, cell.type := sub(".*_", "", sample)]

has.immunoseq <- fread("../secret/immunoseq_samples.txt", header=FALSE)
setnames(has.immunoseq, "immunoseqID")
has.immunoseq[, private := sub("_.*", "", immunoseqID)]
sum(!has.immunoseq$private %in% public.private$private)

public.with.immunoseq <- public.private[private %in% has.immunoseq$private, ]
table(public.with.immunoseq$cell.type)

setkey(has.immunoseq, private)
setkey(public.private, private)

public.private.immunoseq <-
  has.immunoseq[public.private, ][order(immunoseqID), ]

write.csv(public.private.immunoseq, "../secret/public.private.immunoseq.csv")

ihec <- fread("ihec-all-histone.txt")
setnames(ihec, sub(" *$", "", names(ihec)))
##No matches:
ihec[Donor %in% public.with.immunoseq$public, ]
ihec[Sample %in% has.immunoseq$private, ]

## For BluePrint samples that were analyzed at McGill, they were given
## different Sample IDs on the portal. For example, MS029801 is the
## same as the BluePrint sample with private ID S002EV.
(ihec.bp <- ihec[Sample %in% public.with.immunoseq$public, ])
setkey(public.private, public)
ihec.bp$private <- public.private[ihec.bp$Sample, private]
ihec.bp

## For BluePrint samples that were processed at Cambridge, they kept
## their "private" IDs on the portal as Donor. But they do have
## "public" IDs according to David.
(ihec.mcgill <- ihec[Donor %in% has.immunoseq$private &
     grepl("H3K27ac", `Track Name`), ])
setkey(public.private, private)
ihec.mcgill$public <- public.private[ihec.mcgill$Donor, public]
ihec.mcgill

## Seems like there are no samples that have been analyzed by both
## BluePrint/Cambridge and by McGill.
ihec[Sample %in% ihec.mcgill$public, ]
ihec[Donor %in% ihec.mcgill$public, ]
ihec[Sample %in% ihec.bp$private, ]
ihec[Donor %in% ihec.bp$private, ]

## This is a subset of the IHEC samples:
edcc <- fread("edcc-CD4positive-alpha-beta-T-cell.txt")
setnames(edcc, sub(" *$", "", names(edcc)))
edcc[Sample %in% public.with.immunoseq$public, ][order(Donor), ]
