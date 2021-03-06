works_with_R("3.2.2",
             Rsamtools="1.20.4",
             GenomicAlignments="1.2.1",
             IRanges="2.2.7",
             Gviz="1.13.7",
             data.table="1.9.6",
             ggplot2="1.0.1",
             testthat="0.10.0",
             "tdhock/animint@4257e8cf76eb5021a98010b6629b779a4f383b24",
             "tdhock/PeakSegJoint@e51c877e3bb0b5d90efe9a2fb5201e3c102d2543")

st <- c(2000000, 2070000, 2100000, 2160000)
ed <- c(2050000, 2130000, 2150000, 2170000)
a <- AnnotationTrack(start=st, end=ed)
plotTracks(a)

st <- c(2000000, 2070000, 2140000, 2160000)
ed <- c(2050000, 2130000, 2150000, 2170000)
a <- AnnotationTrack(start=st, end=ed)
plotTracks(a)

st <- c(2000000, 2020000, 2060000, 2090000, 2020000)
ed <- c(2050000, 2070000, 2100000, 2130000, 2030000)
a <- AnnotationTrack(start=st, end=ed)
plotTracks(a)

y <- disjointBins(IRanges(st, ed))
df <- data.frame(st, ed, y)
ggplot()+
  geom_segment(aes(st, y, xend=ed, yend=y),
               data=df)

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

ref <- getDNA("chr1", 32398180, 32406111)
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
    data.table(cigar=cigar.dt[, rep(code, bases)])
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
onePeak$seq.i <- 1:nrow(onePeak)
table(nchar(onePeak$seq))
join.list <- list()
onePeak$last <- NA
for(seq.i in 1:nrow(onePeak)){
  cat(sprintf("%4d / %4d reads\n", seq.i, nrow(onePeak)))
  one.seq <- onePeak[seq.i, ]
  align.dt <- with(one.seq, align(ref, seq, first, cigar))
  is.aligned <- which(align.dt$cigar == "M")
  last.aligned <- rev(is.aligned)[1]
  onePeak$last[seq.i] <- align.dt$ref.pos[last.aligned]
  ggplot()+
    geom_text(aes(ref.pos, "seq", label=exp.base),
              data=align.dt)+
    geom_text(aes(ref.pos, "ref", label=ref.base),
              data=ref)
  join.list[[seq.i]] <- data.table(seq.i, align.dt)
}
onePeak[, y := disjointBins(IRanges(first, last))]
join <- do.call(rbind, join.list)
setkey(join, ref.pos)
mutated <- join[ref.base != read.base, ]
others <- join[ref.pos %in% mutated$ref.pos, ]
by.pos <- split(others, others$ref.pos)

prop.mutated <- sapply(by.pos, with, mean(read.base != ref.base))
n.mutated <- sapply(by.pos, with, sum(read.base != ref.base))
n.total <- sapply(by.pos, nrow)
props <-
  data.table(ref.pos=as.integer(names(prop.mutated)),
             prop.mutated,
             n.mutated,
             n.total,
             fake.pos=seq_along(prop.mutated))

mutated.others.list <- list()
show.reads.list <- list()
for(fake.pos in seq_along(prop.mutated)){
  ref.pos <- names(prop.mutated)[fake.pos]
  one.pos <- by.pos[[ref.pos]]
  one.pos[, pos.seq.i := seq_along(seq.i)]
  mutated.others.list[[ref.pos]] <-
    data.table(fake.pos, one.pos)
  some.seqs <- onePeak[one.pos$seq.i, ]
  some.seqs$pos.seq.i <- one.pos$pos.seq.i
  show.reads.list[[ref.pos]] <-
    data.table(ref.pos=as.integer(ref.pos), fake.pos, some.seqs)
}
mutated.others <- do.call(rbind, mutated.others.list)
show.reads <- do.call(rbind, show.reads.list)

ggplot()+
  geom_tallrect(aes(xmin=ref.pos-10.5, xmax=ref.pos+10.5,
                    fill=prop.mutated),
                data=props)+
  scale_fill_gradient(low="black", high="red", limits=c(0, 1))+
  scale_color_manual("mutation", values=c(hg19="black", SNP="red"))+
  geom_text(aes(ref.pos, pos.seq.i, label=read.base,
                color=ifelse(read.base==ref.base, "hg19", "SNP")),
            data=mutated.others)

ggplot()+
  scale_color_manual("mutation", values=c(hg19="black", SNP="red"))+
  geom_text(aes(factor(ref.pos), pos.seq.i, label=read.base,
                color=ifelse(read.base==ref.base, "hg19", "SNP")),
            data=mutated.others)

ggplot()+
  geom_tallrect(aes(xmin=fake.pos-0.5, xmax=fake.pos+0.5,
                    fill=prop.mutated),
                alpha=0.5,
                data=props)+
  scale_fill_gradient(low="black", high="red", limits=c(0, 1))+
  scale_color_manual("mutation", values=c(hg19="black", SNP="red"))+
  geom_text(aes(fake.pos, pos.seq.i, label=read.base,
                color=ifelse(read.base==ref.base, "hg19", "SNP")),
            data=mutated.others)

ggplot()+
  scale_color_manual("mutation", values=c(hg19="black", SNP="red"))+
  ggtitle("this plot is misleading since it shows some reads twice")+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))+
  facet_grid(. ~ ref.pos, scales="free", space="free")+
  geom_segment(aes(first-0.5, pos.seq.i,
                   xend=last+0.5, yend=pos.seq.i),
               data=show.reads)+
  geom_text(aes(ref.pos, pos.seq.i,
                color=ifelse(read.base==ref.base, "hg19", "SNP"),
                label=read.base),
            data=mutated.others)

show.reads[, chromStart := first-1]
show.reads[, chromEnd := last]
clustered <- clusterPeaks(show.reads)
by.cluster <- split(clustered, clustered$cluster)
cluster.name <- "7"
show.seqs.by.cluster <- list()
setkey(props, ref.pos)
setkey(mutated.others, seq.i)
props.by.cluster <- list()
cluster.fac.levs <- paste0("cluster", names(by.cluster))
cluster.info.list <- list()
mo.clustered.list <- list()
for(cluster.i in seq_along(by.cluster)){
  cluster.name <- names(by.cluster)[[cluster.i]]
  cluster.fac <- factor(cluster.fac.levs[[cluster.i]], cluster.fac.levs)
  one.cluster <- by.cluster[[cluster.name]]
  ##cat(cluster.name, "\n");print(table(one.cluster$ref.pos))
  cluster.seq.vec <- unique(one.cluster$seq.i)
  cluster.seqs <- onePeak[cluster.seq.vec, ]
  cluster.info.list[[cluster.name]] <- info <- cluster.seqs[, {
    data.table(cluster.fac,
               clusterStart=min(first),
               clusterEnd=max(last))
  }]
  cluster.seqs$y <- cluster.seqs[, disjointBins(IRanges(first, last))]
  normalize <- function(x){
    info[, (x-clusterStart)/(clusterEnd-clusterStart)]
  }
  cluster.props <- props[J(unique(one.cluster$ref.pos))]
  cluster.props[, `:=`(before.norm=normalize(ref.pos-0.5),
                       after.norm=normalize(ref.pos+0.5))
                ]
  props.by.cluster[[cluster.name]] <-
    data.table(cluster.fac, cluster.props)
  cluster.mo <- mutated.others[J(cluster.seq.vec)]
  setkey(cluster.seqs, seq.i)
  cluster.mo$y <- cluster.seqs[J(cluster.mo$seq.i)]$y
  cluster.mo[, ref.pos.norm := normalize(ref.pos)]
  mo.clustered.list[[cluster.name]] <- data.table(cluster.fac, cluster.mo)
  cluster.seqs[, `:=`(first.norm=normalize(first-0.5),
                      last.norm=normalize(last+0.5))
               ]
  show.seqs.by.cluster[[cluster.name]] <-
    data.table(cluster.fac, cluster.seqs)
}
show.seqs <- do.call(rbind, show.seqs.by.cluster)
clustered.props <- do.call(rbind, props.by.cluster)
cluster.info <- do.call(rbind, cluster.info.list)
mo.clustered <- do.call(rbind, mo.clustered.list)

## All reads, even those that do not overlap SNPs.
mutated.others$global.y <- onePeak[mutated.others$seq.i]$y
allReads <- 
ggplot()+
  geom_tallrect(aes(xmin=(ref.pos-10.5)/1e3, xmax=(ref.pos+10.5)/1e3,
                    fill=prop.mutated),
                alpha=0.5,
                data=props)+
  scale_fill_gradient(low="black", high="red", limits=c(0, 1))+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))+
  ylab("number of aligned reads")+
  xlab("position on chr1 (kilo bases = kb)")+
  scale_color_manual("mutation", values=c(hg19="black", SNP="red"))+
  geom_segment(aes((first-0.5)/1e3, y,
                   xend=(last+0.5)/1e3, yend=y),
               data=onePeak)+
  geom_text(aes(ref.pos/1e3, global.y,
                color=ifelse(read.base==ref.base, "hg19", "SNP"),
                label=read.base),
            data=mutated.others)
pdf("figure-onePeak-allReads.pdf", w=20)
print(allReads)
dev.off()

## All reads that overlap at least one SNP:
mo.clustered <- do.call(rbind, mo.clustered.list)
someReads <- 
ggplot()+
  geom_tallrect(aes(xmin=(ref.pos-1.5)/1e3, xmax=(ref.pos+1.5)/1e3,
                    fill=prop.mutated),
                alpha=0.5,
                data=clustered.props)+
  scale_fill_gradient(low="black", high="red", limits=c(0, 1))+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))+
  facet_grid(. ~ cluster.fac, space="free", scales="free")+
  ylab("number of aligned reads")+
  scale_x_continuous("position on chr1 (kilo bases = kb)",
                     breaks=seq(32398, 32406, by=1))+
  scale_color_manual("mutation", values=c(hg19="black", SNP="red"))+
  geom_segment(aes((first-0.5)/1e3, y,
                   xend=(last+0.5)/1e3, yend=y),
               data=show.seqs)+
  geom_text(aes(ref.pos/1e3, y,
                color=ifelse(read.base==ref.base, "hg19", "SNP"),
                label=read.base),
            data=mo.clustered)

w <- 1200
viz <-
  list(allReads=ggplot()+
         theme_bw()+
         ggtitle("click to zoom")+
         theme_animint(width=w)+
         geom_text(aes(ref.pos/1e3, y,
                       color=ifelse(read.base==ref.base, "hg19", "SNP"),
                       label=read.base),
                   data=mo.clustered)+
         ylab("number of aligned reads")+
         xlab("position on chr1 (kilo bases = kb)")+
         scale_color_manual("mutation", values=c(hg19="black", SNP="red"))+
         geom_segment(aes((first-0.5)/1e3, y,
                          xend=(last+0.5)/1e3, yend=y),
                      data=onePeak)+
         geom_tallrect(aes(xmin=clusterStart/1e3, xmax=clusterEnd/1e3,
                           clickSelects=cluster.fac),
                       alpha=0.5,
                       data=cluster.info),

       zoom=ggplot()+
         theme_bw()+
         ggtitle("zoomed plot")+
         theme_animint(width=w)+
         scale_fill_gradient(low="black", high="red", limits=c(0, 1))+
         ylab("number of aligned reads")+
         scale_x_continuous("",
                            breaks=c())+
         scale_color_manual("mutation", values=c(hg19="black", SNP="red"))+
         geom_segment(aes(first.norm, y,
                          showSelected=cluster.fac,
                          xend=last.norm, yend=y),
                      data=show.seqs)+
         geom_text(aes(ref.pos.norm, y,
                       showSelected=cluster.fac,
                       color=ifelse(read.base==ref.base, "hg19", "SNP"),
                       label=read.base),
                   data=mo.clustered)+
         geom_tallrect(aes(xmin=before.norm, xmax=after.norm,
                           tooltip=sprintf("%d/%d=%.1f%% mutated",
                             n.mutated, n.total, prop.mutated*100),
                           fill=prop.mutated,
                           showSelected=cluster.fac),
                       alpha=0.5,
                       color=NA,
                       data=clustered.props),

       title="SNPs in one sample of ChIP-seq data")

animint2dir(viz, "figure-onePeak")

pdf("figure-onePeak.pdf", w=20)
print(someReads)
dev.off()
