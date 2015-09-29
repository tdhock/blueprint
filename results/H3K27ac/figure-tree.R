works_with_R("3.2.2", ggdendro="0.1.14", data.table="1.9.6",
             flexclust="1.3.4",
             "tdhock/animint@d1c3deb15d52f613f7d3a35a57d1c94a079f0352")

objs <- load("PeakSegJoint.diffs.RData")
load("PeakSegJoint.predictions.RData")

method.vec <- c("ward.D", "ward.D2", "single", "complete", "average",
                "mcquitty", "median", "centroid")
get.type <- function(x)sub("[^a-zA-Z].*", "", x)
segs.list <- list()
labs.list <- list()
best.list <- list()
for(method in method.vec){
  hc <- hclust(PeakSegJoint.diffs, method=method)
  best.list[[method]] <- tryCatch({
    guess.mat <- cutree(hc, h=hc$height)
    true.str <- factor(get.type(rownames(guess.mat)))
    true.int <- as.integer(true.str)
    ARI.list <- list()
    for(guess.col in 1:ncol(guess.mat)){
      ARI.list[[guess.col]] <- randIndex(guess.mat[,guess.col], true.int)
    }
    ARI.vec <- do.call(c, ARI.list)
    best.ARI <- which.max(ARI.vec)
    data.table(method, ARI=ARI.vec, peaks=hc$height)[best.ARI,]
  }, error=function(e){
    data.table(method, ARI=NA, peaks=NA)
  })
  d.data <- dendro_data(hc)
  segs.list[[method]] <- data.table(method, d.data$segments)
  labs.list[[method]] <- data.table(method, d.data$labels)
}
best <- do.call(rbind, best.list)
segs <- do.call(rbind, segs.list)
labs <- do.call(rbind, labs.list)
labs[, type := get.type(label)]

(all.trees <-
  ggplot()+
    ggtitle(paste0("Cluster trees with different linkage methods ",
                  "(", nrow(all.peaks.mat), " samples, ",
                   ncol(all.peaks.mat), " peaks, L1 distance)"))+
  geom_text(aes(0, peaks, label=sprintf("ARI=%.2f", ARI)),
            color="black",
             data=best[!is.na(ARI),], hjust=1, vjust=-0.5)+
  geom_hline(aes(yintercept=peaks),
             color="black",
             data=best[!is.na(ARI),])+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))+
  facet_grid(method ~ ., scales="free")+
  geom_segment(aes(x, y,
                   xend=xend, yend=yend),
               color="grey50",
               size=1,
               data=segs)+
  geom_point(aes(x, y, label=label,
                 clickSelects=label,
                 color=type),
             pch=1,
             angle=90,
             hjust=1,
             alpha=0.75,
             data=labs)
 )

viz <-
  list(trees=all.trees+
         theme_animint(width=1500, height=1000),

       title="H3K27ac BluePrint data samples",

       first=list(label=character()),

       selector.types=list(label="multiple"))

animint2dir(viz, "figure-tree", open.browser=FALSE)
