works_with_R("3.2.2", data.table="1.9.6", 
             caret="6.0.41")

load("gctest.RData")

GC.vec <- gcBins$count
w.size <- 15
only.consider <- w.size:(length(GC.vec)-w.size)

H3K27ac <- coverageBins[experiment == "H3K27ac", ]
coverage.by.person <- split(H3K27ac, H3K27ac$person)
data.by.person <- list()
for(person in names(coverage.by.person)){
  person.coverage <- coverage.by.person[[person]]
  coverage.by.center <- split(person.coverage, person.coverage$center)
  output.vec <- coverage.by.center$NCMLS$norm
  McGill.vec <- coverage.by.center$McGill$norm
  change.mat <- sapply(only.consider, function(i){
    first.i <- max(1, i-w.size)
    last.i <- min(length(McGill.vec), i+w.size)
    before.vec <- McGill.vec[first.i:i]
    after.vec <- McGill.vec[i:last.i]
    before.diff <- diff(before.vec)
    after.diff <- diff(after.vec)
    c(before.pos=sum(before.diff[0 < before.diff]),
      before.neg=sum(before.diff[before.diff < 0]),
      after.pos=sum(after.diff[0 < after.diff]),
      after.neg=sum(after.diff[after.diff < 0]))
  })
  data.by.person[[person]] <- list(
    features=cbind(GC=GC.vec[only.consider], t(change.mat)),
    output=output.vec[only.consider])
}

test.person.vec <- c("S007VE", "S002KJ")
is.test <- names(data.by.person) %in% test.person.vec
set.list <- list(train.validation=!is.test, test=is.test)
data.by.set <- list()
for(set.name in names(set.list)){
  is.set <- set.list[[set.name]]
  set.by.person <- data.by.person[is.set]
  feature.mat.list <- list()
  output.vec.list <- list()
  for(person in names(set.by.person)){
    person.list <- set.by.person[[person]]
    feature.mat.list[[person]] <- person.list$features
    output.vec.list[[person]] <- person.list$output
  }
  data.by.set[[set.name]] <- 
    list(output.vec=do.call(c, output.vec.list),
         feature.mat=do.call(rbind, feature.mat.list))
}

fit <- with(data.by.set$train.validation, {
  train(feature.mat, output.vec, method="rf", preProcess=c("center", "scale"),
        tuneLength=10)
})

pred.vec <- predict(fit, data.by.set$test$feature.mat)
plot(pred.vec, data.by.set$test$output.vec - pred.vec)

pred.by.person <- list()
for(person in names(coverage.by.person)){
  person.coverage <- coverage.by.person[[person]]
  McGill.coverage <- person.coverage[center == "McGill", ]
  McGill.coverage$set.name <- ifelse(
    person %in% test.person.vec, "test", "train")
  McGill.coverage$pred <- NA
  person.data <- data.by.person[[person]]
  McGill.coverage$pred[only.consider] <- predict(fit, person.data$features)
  pred.by.person[[person]] <-  McGill.coverage
}

gctest.pred <- do.call(rbind, pred.by.person)

save(gctest.pred, file="gctest.pred.RData")
