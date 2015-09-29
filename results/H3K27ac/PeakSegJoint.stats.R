works_with_R("3.2.2", data.table="1.9.6")

load("PeakSegJoint.counts.RData")
load("PeakSegJoint.predictions.RData")

cell.type <- sub("[0-9]*/.*", "", rownames(all.peaks.mat))
type.tab <- table(cell.type)

wide <- PeakSegJoint.counts[, list(
  mono=sum(peaks[cell.type=="mono"]),
  Input=sum(peaks[cell.type=="Input"]),
  neutro=sum(peaks[cell.type=="neutro"])
  ), by=peak.name]
wide[, all := mono + neutro + Input]
wide[, max := ifelse(mono < neutro, neutro, mono)]
wide[, specificity := max/all]

posPattern <-
  paste0("(?<chrom>chr.*?)",
         ":",
         "(?<chromStart>.*?)",
         "-",
         "(?<chromEnd>.*?)",
         "$")
pos.mat <- str_match_perl(wide$peak.name, posPattern)

getStr <- function(cell.type){
  count <- get(cell.type)
  ifelse(0 < count,
         sprintf("%d/%d %s", count, type.tab[[cell.type]], cell.type),
         NA)
}

count.str <- wide[, {
  vec.list <- list()
  for(cell.type in c("mono", "neutro", "Input")){
    count <- get(cell.type)
    vec.list[[cell.type]] <- 
      ifelse(0 < count,
             sprintf("%d/%d%s", count, type.tab[[cell.type]], cell.type),
             NA)
  }
  type.mat <- do.call(cbind, vec.list)
  apply(type.mat, 1, function(type.vec){
    not.na <- type.vec[!is.na(type.vec)]
    paste(not.na, collapse=",")
  })
}]

bed.dt <- data.table(
  chrom=pos.mat[, "chrom"],
  chromStart=pos.mat[, "chromStart"],
  chromEnd=pos.mat[, "chromEnd"],
  name=count.str
  )

con <- gzfile("PeakSegJoint.types.bed.gz", "w")
header <- 
  paste('track',
        'visibility=pack',
        'name=PeakSegJoint',
        'description="PeakSegJoint at least',
        sprintf("1/%d", sum(type.tab)),
        'sample"')
writeLines(header, con)
write.table(bed.dt, file=con,
            sep="\t", quote=FALSE,
            row.names=FALSE,
            col.names=FALSE)
close(con)

ggplot()+
  geom_point(aes(all, specificity), data=wide)

ggplot()+
  geom_point(aes(max, specificity), data=wide)

dot.dt <- wide[, list(count=.N), by=.(mono, neutro)][order(count),]
dot.dt[, dotID := sprintf("mono%d.neutro%d", mono, neutro)]
dot.dt[, tooltip := paste(count, "peaks in", dotID)]

chromfactor <- function(x){
  factor(x, paste0("chr", c(1:22, "X", "Y", "M")))
}
genome.dt <- data.table(
  wide,
  chrom=chromfactor(pos.mat[,"chrom"]),
  peakStart=as.integer(pos.mat[,"chromStart"]),
  peakEnd=as.integer(pos.mat[,"chromEnd"]))
genome.dt[, dotID := sprintf("mono%d.neutro%d", mono, neutro)]
genome.dt[, bases := peakEnd-peakStart]
genome.dt[, middle := as.integer((peakStart+peakEnd)/2)]
genome.dt[, zoom.bases := bases * 10]
genome.dt[, zoomStart := as.integer(middle-zoom.bases)]
genome.dt[, zoomEnd := as.integer(middle+zoom.bases)]

hgTracks <- "http://genome.ucsc.edu/cgi-bin/hgTracks"
library(animint)
viz <- list(
  genome=ggplot()+
    ggtitle("Number of Input samples with peak")+
    theme_bw()+
    scale_y_discrete("hg19 chromosome", drop=FALSE)+
    scale_x_continuous("position (mega bases)")+
    theme_animint(width=1000)+
    geom_text(aes(125, chromfactor("chrY"),
                  showSelected=dotID,
                  label=tooltip),
              data=dot.dt)+
    geom_text(aes(middle/1e6, chrom,
                   href=sprintf("%s?db=hg19&position=%s:%d-%d",
                     hgTracks, chrom, zoomStart, zoomEnd),
                   label=Input,
                   showSelected=dotID),
               pch=1,
               data=genome.dt),
  first=list(dotID="mono184.neutro193"),
  counts=ggplot()+
    theme_bw()+
    geom_vline(xintercept=type.tab[["mono"]], color="grey")+
    geom_hline(yintercept=type.tab[["neutro"]], color="grey")+
    coord_equal()+
    theme_animint(width=800, height=800)+
    geom_point(aes(mono, neutro, color=log10(count),
                   tooltip=tooltip, 
                   clickSelects=dotID),
               alpha=0.7,
               data=dot.dt),
  title="PeakSegJoint on 393 H3K27ac+Input BluePrint samples")

animint2dir(viz, "PeakSegJoint-genome")

sorted <- wide[order(specificity, max, decreasing=TRUE), ]
head(sorted, 100)
