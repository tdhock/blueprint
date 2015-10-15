works_with_R("3.2.2",
             data.table="1.9.6",
             "tdhock/ggplot2@a8b06ddb680acdcdbd927773b1011c562134e4d2")

load("gctest.RData")

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
