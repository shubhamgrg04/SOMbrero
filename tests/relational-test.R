## Check that relational SOM gives positive coordinates that sum to 1

library(SOMbrero)

iris.dist <- dist(iris[1:30,1:4], method="minkowski", diag=TRUE, upper=TRUE, 
                  p=4)

rsom <- trainSOM(x.data=iris.dist, type="relational", maxit= 10, 
                 scaling= "none")
stopifnot(all.equal(as.vector(rowSums(rsom$prototypes)), 
                    rep(1, prod(rsom$parameters$the.grid$dim))))
stopifnot(!sum(rsom$prototypes<0))

rsom <- trainSOM(x.data=iris.dist, type="relational", maxit= 10, 
                 scaling= "sd")
stopifnot(all.equal(as.vector(rowSums(rsom$prototypes)), 
                    rep(1, prod(rsom$parameters$the.grid$dim))))
stopifnot(!sum(rsom$prototypes<0))

rsom <- trainSOM(x.data=iris.dist, type="relational", maxit= 10, 
                 scaling= "max")
stopifnot(all.equal(as.vector(rowSums(rsom$prototypes)), 
                    rep(1, prod(rsom$parameters$the.grid$dim))))
stopifnot(!sum(rsom$prototypes<0))

rsom <- trainSOM(x.data=iris.dist, type="relational", maxit= 10, 
                 scaling= "frobenius")
stopifnot(all.equal(as.vector(rowSums(rsom$prototypes)), 
                    rep(1, prod(rsom$parameters$the.grid$dim))))
stopifnot(!sum(rsom$prototypes<0))

rsom <- trainSOM(x.data=iris.dist, type="relational", maxit= 10, 
                 scaling= "cosine")
stopifnot(all.equal(as.vector(rowSums(rsom$prototypes)), 
                    rep(1, prod(rsom$parameters$the.grid$dim))))
stopifnot(!sum(rsom$prototypes<0))
