works_with_R("3.2.2",
             data.table="1.9.6",
             "tdhock/ggplot2@a8b06ddb680acdcdbd927773b1011c562134e4d2")

load("gctest.RData")
load("gctest.pred.RData")

gctest[, list(max=max(count)), by=.(person, cellType, experiment, center)]

table(gc5Base$percent)
gcTrack <- data.table(
  res.bases=5,
  person="GC", cellType="GC", experiment="GC", center="GC",
  gc5Base)
gcBinsTrack <- data.table(
  res.bases=100,
  person="GC", cellType="GC", experiment="GC", center="GC",
  gcBins)
gctest$res.bases <- 1
coverageBins$res.bases <- 100

##fig.norm <- 
  ggplot()+
    ggtitle("base pair resolution (black) and 100bp mean (red)")+
    geom_line(aes(position/1e3, percent/100), data=gcTrack)+
    geom_point(aes(chromStart/1e3, mean/5),
               color=bin.color,
               shape=1,
               data=gcBinsTrack)+
    scale_y_continuous("", breaks=function(limits){
      c(0.5, 1)
    })+
    geom_rect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
                  ymin=0, ymax=norm),
              data=gctest)+
    geom_point(aes(chromStart/1e3, norm),
               color=bin.color,
               shape=1,
               data=coverageBins)+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "cm"))+
    facet_grid(person + cellType + experiment + center ~ .)

pred.fig <- 
ggplot()+
  ggtitle("Random forest predictions for McGill samples")+
    geom_line(aes(position/1e3, percent/100), data=gcTrack)+
    scale_y_continuous("", breaks=function(limits){
      c(0.5, 1)
    })+
    geom_rect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
                  ymin=0, ymax=norm),
              data=gctest)+
    geom_point(aes((chromStart+chromEnd)/2e3, pred, color=set.name),
               pch=1,
               data=gctest.pred)+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "cm"))+
    facet_grid(person + cellType + experiment + center ~ .)

png("figure-gcpred-coverage.png", width=14, height=10, units="in", res=200)
print(pred.fig)
dev.off()

load("../results/H3K27ac/PeakSegJoint.predictions.RData")
all.peaks.dt <- data.table(all.peaks.df)
setkey(all.peaks.dt, sample.id)
type.code <- c(
  Mono="CD14-positive_CD16-negative_classical_monocyte",
  GR="mature_neutrophil")
train.sample.ids <- unique(gctest[, {
  paste(person, center, type.code[paste(cellType)], sep="_")
}])
not.input <- all.peaks.dt[train.sample.ids][sample.group != "Input", ]
more.NCMLS <- not.input[, {
  list(McGill=sum(grepl("McGill", sample.id)),
       NCMLS=sum(grepl("NCMLS", sample.id)))
}, by=peak.name][McGill < NCMLS, ]
stop("TODO train model using regions around more.NCMLS")

bin.color <- "#E41A1C"
fig.cov <- 
  ggplot()+
    ggtitle("base pair resolution (black) and 100bp mean (red)")+
    geom_line(aes(position/1e3, percent), data=gcTrack)+
    geom_point(aes(chromStart/1e3, mean/5 * 100),
               color=bin.color,
               shape=1,
               data=gcBinsTrack)+
    scale_y_continuous("", breaks=function(limits){
      floor(limits[2])
    })+
    geom_rect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
                  ymin=0, ymax=count),
              data=gctest)+
    geom_point(aes(chromStart/1e3, mean),
               color=bin.color,
               shape=1,
               data=coverageBins)+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "cm"))+
    facet_grid(person + cellType + experiment + center ~ ., scales="free")

png("figure-gctest-coverage.png", width=14, height=10, units="in", res=200)
print(fig.cov)
dev.off()
