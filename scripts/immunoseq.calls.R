works_with_R("3.2.2",
             data.table="1.9.6",
             "tdhock/namedCapture@a31a38be12d4dad4aec6257a47c0e2307679589c")

load("immunoseq.RData") #list of data.frames with columns:
sapply(immunoseq, dim)

## chr_pos: chromosome and position
## CHROM= chromosome number
## POS= position of the variant
## ID= rs ID (when unknown: “.”)
## REF= reference allele
## ALT= alternative allele
## QUAL= quality of the variant call
## S00*= sample variant call
## gq= number of samples that have good quality for this variant ( read depth>=10 and genotyping call >=40, mapping quality >=50)
## hom_ref= number of samples that are homozygous ref
## het= number of samples that are heterozygous
## hom_alt= number of samples that are homozygous alt
## het2= if there are samples that are heterozygous for a different alt allele.

## The sample variant call is (according to Andréanne):

## GT:AD:DP:GQ:PL
## GT= genotype (./.= no call, 0/0= hom ref, 0/1= heterozygous, 1/1= homozygous alt)
## DP= read depth
## GQ= genotyping quality

genotype.code <- c(
  "0/0"="homozygousRef",
  "0/1"="heterozygous",
  "1/1"="homozygousAlt")
alt.code <- c(
  homozygousRef=0,
  heterozygous=1,
  homozygousAlt=2)

## AD and PL I don’t know by heart (I don’t use them frequently) but
## you can find the meaning if you search for vcf file format.

##                            S0137X_WB_gDNA_1
## 1:                 0/1:14,8:22:99:219,0,418
## 2:                0/1:25,31:56:99:997,0,728
## 3:                0/1:12,15:27:99:439,0,325
## 4:                  0/0:23,0:23:63:0,63,861
## 5:                  0/0:20,0:20:60:0,60,679
## 6: 0/1:25,18,0:43:99:406,0,784,481,838,1319

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

immunoseq.calls <- list()
for(variant.type in names(immunoseq)){
  variant.dt <- immunoseq[[variant.type]]
  sample.vec <- grep("^S", names(variant.dt), value=TRUE)
  person.vec <- sub("_.*", "", sample.vec)
  stopifnot(table(person.vec) == 1) # no duplicate people.
  nAlt.mat <- matrix(NA, nrow(variant.dt), length(person.vec),
                     dimnames=list(variant.dt$chr_pos, person.vec))
  for(sample.i in seq_along(sample.vec)){
    s <- sample.vec[[sample.i]]
    cat(sprintf("%s %4d / %4d %s\n",
                variant.type, sample.i, length(sample.vec), s))
    variant.str <- variant.dt[[s]]
    match.mat <- str_match_named(variant.str, pattern)
    stopifnot(sum(is.na(match.mat[, 1])) == sum(is.na(variant.str)))
    nRef <- rowSums(match.mat[, c("allele1", "allele2")] == "0")
    nAlt <- 2-nRef
    person <- person.vec[[sample.i]]
    nAlt.mat[, person] <- nAlt
  }
  immunoseq.calls[[variant.type]] <- nAlt.mat
}

save(immunoseq.calls, file="immunoseq.calls.RData")
