works_with_R("3.2.2",
             data.table="1.9.6",
             "tdhock/ggplot2@a8b06ddb680acdcdbd927773b1011c562134e4d2")

load("gctest.RData")

gctest[, list(max=max(count)), by=.(person, cellType, experiment, center)]

fig.cov <- 
  ggplot()+
    scale_y_continuous("", breaks=function(limits){
      floor(limits[2])
    })+
    geom_rect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
                  ymin=0, ymax=count),
              data=gctest)+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "cm"))+
    facet_grid(person + cellType + experiment + center ~ ., scales="free")

stop("TODO add GC content track, plot bins")

png("figure-gctest-coverage.png", width=14, height=10, units="in", res=200)
print(fig.cov)
dev.off()
