topographicError <- function (sommap) {
  norm.data <- preprocessData(sommap$data, sommap$parameters$scaling)
  norm.proto <- preprocessProto(sommap$prototypes, sommap$parameters$scaling,
                                sommap$data)
  if (sommap$parameters$type=="numeric") {
    all.dist <- apply(norm.data, 1, function(x) {
      apply(norm.proto,1,function(y) sum((x-y)^2) )
    })
    ind.winner2 <- apply(all.dist,2,function(x) order(x)[2])
  } else if (sommap$parameters$type=="korresp") {
    nr <- nrow(sommap$data)
    nc <- ncol(sommap$data)
    all.dist.row <- apply(norm.data[1:nr,1:nc], 1, function(x) {
      apply(sommap$prototypes[,1:nc], 1, function(y) sum((x-y)^2) )
    })
    all.dist.col <- apply(norm.data[(nr+1):(nr+nc),(nc+1):(nr+nc)],
                          1, function(x) {
      apply(sommap$prototypes[,(nc+1):(nr+nc)], 1, function(y) sum((x-y)^2) )
    })
    ind.winner2 <- c(apply(all.dist.col,2,function(x) order(x)[2]),
                     apply(all.dist.row,2,function(x) order(x)[2]))
  } else if (sommap$parameters$type=="relational") {
    all.dist <- sapply(1:ncol(norm.proto), function(ind) {
      norm.proto%*%norm.data[ind,]-
        0.5*diag(norm.proto%*%norm.data%*%t(norm.proto))
    })
    ind.winner2 <- apply(all.dist,2,function(x) order(x)[2])
  }
  res.error <- mean(!sapply(1:nrow(sommap$data), function(x) {
      is.element(ind.winner2[x], selectNei(sommap$clustering[x],
                                           sommap$parameters$the.grid, 1,
                                           radius.type= "letremy"))
  }))
  return(res.error)
}

quantizationError <- function(sommap) {
  norm.data <- preprocessData(sommap$data, sommap$parameters$scaling)
  norm.proto <- preprocessProto(sommap$prototypes, sommap$parameters$scaling,
                                sommap$data)
  if (sommap$parameters$type=="numeric") {
    quantization.error <- sum(apply((norm.data-
                                       norm.proto[sommap$clustering,])^2,
                                    1,sum))/nrow(norm.data)
  } else if (sommap$parameters$type=="korresp") {
    nr <- nrow(sommap$data)
    nc <- ncol(sommap$data)
    quantization.error <- sum(apply((norm.data[1:nr,1:nc]-
                                       sommap$prototypes[sommap$clustering[
                                         (nc+1):(nc+nr)],1:nc])^2,1,sum))
    quantization.error <- quantization.error +
      sum(apply((norm.data[(nr+1):(nr+nc),(nc+1):(nr+nc)]-
                   sommap$prototypes[sommap$clustering[1:nc],(nc+1):(nr+nc)])^2,
                1,sum))
    quantization.error <- quantization.error / (nc+nr)
  } else if (sommap$parameters$type=="relational") {
    clust.proto <- norm.proto[sommap$clustering,]
    quantization.error <- clust.proto%*%norm.data - 0.5*
      tcrossprod(diag(clust.proto%*%norm.data%*%t(clust.proto)),
                 rep(1,ncol(norm.data)))
    quantization.error <- sum(diag(quantization.error))/nrow(norm.data)
  }
    
  quantization.error
}

kaskiLagusError <- function(sommap) {
  norm.data <- preprocessData(sommap$data, sommap$parameters$scaling)
  norm.proto <- preprocessProto(sommap$prototypes, sommap$parameters$scaling,
                                sommap$data)
  # Quantization error and computation of first and second closest prototypes
  if (sommap$parameters$type=="numeric") {
    obs.proto.dist <- t(apply(norm.data, 1, 
                              function(x) sqrt(colSums((t(norm.proto)-x)^2))))
    winners <- t(apply(obs.proto.dist, 1, function(x) order(x)[1:2]))
    quantization.error <- mean(obs.proto.dist[cbind(1:nrow(norm.data), 
                                                    winners[,1])])
  } else if (sommap$parameters$type=="korresp") {
    nr <- nrow(sommap$data)
    nc <- ncol(sommap$data)
    all.dist.row <- apply(norm.data[1:nr,1:nc], 1, function(x) {
      colSums((t(norm.proto[,1:nc])-x)^2)
    })
    all.dist.col <- apply(norm.data[(nr+1):(nr+nc),(nc+1):(nr+nc)], 1, 
                          function(x) 
                            sqrt(colSums((t(norm.proto[,(nc+1):(nr+nc)])-x)^2)))
    winners <- t(cbind(apply(all.dist.row,2,function(x) order(x)[1:2]),
                       apply(all.dist.col,2,function(x) order(x)[1:2])))
    quantization.error <- mean(c(all.dist.row[cbind(winners[1:nr,1],1:nr)],
                                 all.dist.col[cbind(winners[(nr+1):(nr+nc),1],
                                                    1:nc)]))
  } else if (sommap$parameters$type=="relational") {
    obs.proto.dist <- t(sapply(1:ncol(norm.proto), function(ind) {
      sqrt(norm.proto%*%norm.data[ind,]-
             0.5*diag(norm.proto%*%norm.data%*%t(norm.proto)))
    }))
    winners <- t(apply(obs.proto.dist,1,function(x) order(x)[1:2]))
    quantization.error <- mean(obs.proto.dist[cbind(1:nrow(norm.data), 
                                                    winners[,1])])
  }
  
  # Compute shortest path for all pairs of prototypes
  proto.dist <- protoDist(sommap, "neighbors")
  paths <- matrix(NA, ncol= length(proto.dist), nrow= length(proto.dist))
  for (i.prot in 1:nrow(paths))
    paths[i.prot, as.numeric(names(proto.dist[[i.prot]]))] <- 
    proto.dist[[i.prot]]
  if (sommap$parameters$type=="relational")
    paths <- sqrt(paths)
  paths <- allShortestPaths(paths)$length
  
  quantization.error + mean(paths[winners])
}

# main function
quality.somRes <- function(sommap, quality.type=c("all", "quantization",
                                                  "topographic", 
                                                  "kaski.lagus"), ...) {
  quality.type <- match.arg(quality.type)
  switch(quality.type,
         "all"=list("topographic"=topographicError(sommap),
                    "quantization"=quantizationError(sommap),
                    "kaski.lagus"=kaskiLagusError(sommap)),
         "topographic"=topographicError(sommap),
         "quantization"=quantizationError(sommap),
         "kaski.lagus"=kaskiLagusError(sommap)
         )
}

quality <- function(sommap, quality.type,...) {
  UseMethod("quality")
}
