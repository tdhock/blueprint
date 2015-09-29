works_with_R("3.2.2", data.table="1.9.6")

load("PeakSegJoint.predictions.RData")

str(all.peaks.mat)

cell.type <- factor(sub("[0-9]*/.*", "", rownames(all.peaks.mat)))
count.total <- table(cell.type)
counts.by.peak <- list()
for(peak.i in 1:ncol(all.peaks.mat)){
  peak.name <- colnames(all.peaks.mat)[peak.i]
  cat(sprintf("%4d / %4d peaks\n", peak.i, ncol(all.peaks.mat)))
  peak.vec <- all.peaks.mat[, peak.name]
  count.tab <- table(cell.type[peak.vec==1])
  counts.by.peak[[peak.name]] <- 
    data.table(peak.name,
               cell.type=names(count.tab),
               peaks=as.integer(count.tab),
               samples=as.integer(count.total))
}
PeakSegJoint.counts <- do.call(rbind, counts.by.peak)

save(PeakSegJoint.counts, file="PeakSegJoint.counts.RData")

