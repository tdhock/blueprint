works_with_R("3.2.2", ggdendro="0.1.14", data.table="1.9.6",
             flexclust="1.3.4",
             "tdhock/animint@4257e8cf76eb5021a98010b6629b779a4f383b24")

objs <- load("PeakSegJoint.diffs.RData")
load("PeakSegJoint.predictions.RData")
load("../../scripts/mislabeled.RData")
bad.H3K27ac <- subset(mislabeled, experiment=="H3K27ac")
type.code <-
  c(Mono="CD14-positive_CD16-negative_classical_monocyte",
    GR="mature_neutrophil")
rownames(bad.H3K27ac) <- with(bad.H3K27ac, {
  paste(private, where, type.code[paste(cell.type)], sep="_")
})

method.vec <- c("ward.D", "ward.D2", "single", "complete", "average",
                "mcquitty", "median", "centroid")
get.type <- function(x)sub("[^a-zA-Z].*", "", x)
segs.list <- list()
labs.list <- list()
best.list <- list()
for(method in method.vec){
  hc <- hclust(PeakSegJoint.diffs, method=method)
  d.data <- dendro_data(hc)
  method.labs <- d.data$labels
  method.labs$type <- get.type(method.labs$label)
  label.vec <- sub(".*/", "", method.labs$label)
  prob.or.NA <-
    ifelse(method.labs$type=="Input", NA, bad.H3K27ac[label.vec, "problem"])
  method.labs$problem <- ifelse(is.na(prob.or.NA), "OK", "abnormal")
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
  segs.list[[method]] <- data.table(method, d.data$segments)
  labs.list[[method]] <- data.table(method, method.labs)
}
best <- do.call(rbind, best.list)
segs <- do.call(rbind, segs.list)
labs <- do.call(rbind, labs.list)

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
  scale_color_manual(values=c(OK="white", abnormal="black"))+
  geom_point(aes(x, y, label=label,
                 clickSelects=label,
                 color=problem,
                 fill=type),
             pch=21,
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
