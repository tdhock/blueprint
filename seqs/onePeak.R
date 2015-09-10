works_with_R("3.2.2", data.table="1.9.5", ggplot2="1.0")

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
  ref.diff <- ifelse(aligned$cigar == "I", 0, 1)
  aligned$ref.pos <- ifelse(ref.diff==0, NA, start + cumsum(ref.diff) - 1)
  read.diff <- ifelse(aligned$cigar == "D", 0, 1)
  aligned$read.pos <- ifelse(read.diff==0, NA, cumsum(read.diff))
  aligned$read.base <- read.base[aligned$read.pos]
  setkey(ref.dt, ref.pos)
  aligned$ref.base <- ref.dt[aligned$ref.pos]$ref.base
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
  all.equal(test.expected[[col.name]], test.computed[[col.name]])
}

line2 <- getDNA("chr1", 32398237, 32398278)
onePeak <- fread("onePeak.sam", drop=c(1:2, 5, 7:9, 11:21))
setnames(onePeak, c("chrom", "first", "cigar", "seq"))
table(nchar(onePeak$seq))
join.list <- list()
for(seq.i in 1:nrow(onePeak)){
  one.seq <- onePeak[seq.i, ]
  exp.base <- unlist(strsplit(one.seq$seq, split=""))
  cigar.mat <- cigar.list[[seq.i]]

  last <- one.seq$first + length(exp.base)-1
  setkey(seq.dt, chrom, position)
  
  ggplot()+
    geom_text(aes(position, "seq", label=exp.base),
              data=seq.dt)+
    geom_text(aes(position, "ref", label=ref.base),
              data=ref)

  one.join <- ref[seq.dt]
  join.list[[seq.i]] <- data.table(seq.i, one.join)
}
join <- do.call(rbind, join.list)
