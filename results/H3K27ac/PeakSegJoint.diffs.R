works_with_R("3.2.2", data.table="1.9.6")

objs <- load("PeakSegJoint.predictions.RData")

PeakSegJoint.diffs <- dist(all.peaks.mat, "manhattan")

save(PeakSegJoint.diffs, file="PeakSegJoint.diffs.RData")
