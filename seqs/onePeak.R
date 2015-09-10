works_with_R("3.2.2", data.table="1.9.5", ggplot2="1.0.1",
             testthat="0.10.0")

getDNA <- function(chrom, from, to, genome="hg19"){
  u <-
    sprintf("http://genome.ucsc.edu/cgi-bin/das/%s/dna?segment=%s:%d,%d",
            genome, chrom, from, to)
  xml.lines <- readLines(u)
  dna.lines <- grep("[<]", xml.lines, invert=TRUE, value=TRUE)
  lower <- unlist(strsplit(dna.lines, split=""))
  ref.base <- toupper(lower)
  data.table(chrom, ref.pos=from:to,
             ref.base)
}

ref <- getDNA("chr1", 32398211, 32406079)
setkey(ref, chrom, ref.pos)

cigar.pattern <-
  paste0("(?<bases>[0-9]+)",
         "(?<code>[^0-9]+)")
cigar.code <-
  c(M="alignment match or mismatch",
    I="insertion to the reference",
    D="deletion from the reference",
    N="skipped region from the reference",
    S="soft clipping (clipped sequences present in SEQ)",
    H="hard clipping (clipped sequences NOT present in SEQ)",
    P="padding (silent deletion from padded reference)",
    "="="sequence match",
    X="sequence mismatch")
align <- function(ref.dt, read, start, cigar){
  stopifnot(is.character(read))
  stopifnot(length(read) == 1)
  stopifnot(is.integer(start))
  stopifnot(length(start) == 1)
  stopifnot(is.character(cigar))
  stopifnot(length(cigar) == 1)
  read.base <- strsplit(read, split="")[[1]]
  cigar.mat <- str_match_all_perl(cigar, cigar.pattern)[[1]]
  cigar.dt <-
    data.table(bases=as.integer(cigar.mat[, "bases"]),
               code=cigar.mat[, "code"])
  bad.cigar <- "[^ISDM]"
  if(any(grepl(bad.cigar, cigar.dt$code))){
    print(cigar.dt)
    stop("cigar matched ", bad.cigar)
  }
  stopifnot(sum(cigar.dt[grepl("[MIS=X]", code), bases]) == length(read.base))
  aligned <- 
    data.table(cigar=character(sum(cigar.dt$bases)))
  first.row <- 1
  for(cigar.i in 1:nrow(cigar.dt)){
    cigar <- cigar.dt[cigar.i, ]
    last.row <- first.row + cigar$bases - 1
    aligned$cigar[first.row:last.row] <- cigar$code
    first.row <- last.row + 1
  }
  ref.diff <- ifelse(grepl("[IS]", aligned$cigar), 0, 1)
  aligned$ref.pos <- ifelse(ref.diff==0, NA, start + cumsum(ref.diff) - 1L)
  read.diff <- ifelse(aligned$cigar == "D", 0L, 1L)
  aligned$read.pos <- ifelse(read.diff==0, NA, cumsum(read.diff))
  aligned$read.base <- read.base[aligned$read.pos]
  setkey(ref.dt, ref.pos)
  aligned$ref.base <- ref.dt[J(aligned$ref.pos)]$ref.base
  aligned
}  

## From http://genome.sph.umich.edu/wiki/SAM#What_is_a_CIGAR.3F
ref.txt <- "C  C  A  T  A  C  T     G  A  A  C  T  G  A  C  T  A  A  C"
read.no.spaces <- gsub(" ", "", "A  C  T  A  G  A  A     T  G  G  C  T")
ref.no.spaces <- gsub(" ", "", ref.txt)
test.ref <-
  data.table(ref.pos=1:19,
             ref.base=strsplit(ref.no.spaces, split="")[[1]])
test.cigar <- "3M1I3M1D5M"
test.start <- 5L
(test.computed <- align(test.ref, read.no.spaces, test.start, test.cigar))
ref.base  <- c("A", "C", "T",  NA, "G", "A", "A", "C", "T", "G", "A", "C", "T")
read.base <- c("A", "C", "T", "A", "G", "A", "A",  NA, "T", "G", "G", "C", "T")
test.expected <-
  data.table(read.pos=c(1:7, NA, 8:12),
             ref.pos=c(5:7, NA, 8:16),
             read.base, ref.base)
for(col.name in names(test.expected)){
  expect_equal(test.computed[[col.name]], test.expected[[col.name]])
}

## From the 26th line of onePeak.sam.
line2 <- getDNA("chr1", 32398237, 32398278)
cigar26 <- "2S40M"
seq26 <- "AAACTGACACCTCTATACGTAGTTAATAGTTCAAAAGATGAA"
"ACTGACACCTCTATACGTAGTTAATAGTTCAAAAGATGAA"
start26 <- 32398549L
test.computed <- align(ref, seq26, start26, cigar26)
test.expected <-
  data.table(ref.pos=c(NA, NA, start26:(start26+39)),
             read.pos=1:42)
for(col.name in names(test.expected)){
  expect_equal(test.computed[[col.name]], test.expected[[col.name]])
}

onePeak <- fread("onePeak.sam", drop=c(1:2, 5, 7:9, 11:21))
setnames(onePeak, c("chrom", "first", "cigar", "seq"))
table(nchar(onePeak$seq))
join.list <- list()
for(seq.i in 1:nrow(onePeak)){
  cat(sprintf("%4d / %4d reads\n", seq.i, nrow(onePeak)))
  one.seq <- onePeak[seq.i, ]
  align.dt <- with(one.seq, align(ref, seq, first, cigar))

  ggplot()+
    geom_text(aes(ref.pos, "seq", label=exp.base),
              data=align.dt)+
    geom_text(aes(ref.pos, "ref", label=ref.base),
              data=ref)

  join.list[[seq.i]] <- data.table(seq.i, align.dt)
}
join <- do.call(rbind, join.list)
setkey(join, ref.pos)
mutated <- join[ref.base != read.base, ]
others <- join[ref.pos %in% mutated$ref.pos, ]
by.pos <- split(others, others$ref.pos)
mutated.props <- sort(sapply(by.pos, with, mean(read.base != ref.base)))
by.pos[names(mutated.props)]
## TODO: plot these SNPs in context of the peak.
