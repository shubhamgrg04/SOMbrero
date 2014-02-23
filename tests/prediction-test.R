## Check that predict.somRes works fine ("numeric" and "relational" cases)

library(SOMbrero)
data(lesmis)
data(presidentielles2002)

nsom <- trainSOM(iris[1:30,1:4], maxit=10, scaling= "none")
stopifnot(identical(predict(nsom, iris[1:30,1:4]), nsom$clustering))
stopifnot(predict(nsom, iris[1,1:4])==nsom$clustering[1])
stopifnot(identical(predict(nsom), nsom$clustering))

nsom <- trainSOM(iris[1:30,1:4], maxit=10, scaling= "center")
stopifnot(identical(predict(nsom, iris[1:30,1:4]), nsom$clustering))
stopifnot(predict(nsom, iris[1,1:4])==nsom$clustering[1])
stopifnot(identical(predict(nsom), nsom$clustering))

nsom <- trainSOM(iris[1:30,1:4], maxit=10, scaling= "unitvar")
stopifnot(identical(predict(nsom, iris[1:30,1:4]), nsom$clustering))
stopifnot(predict(nsom, iris[1,1:4])==nsom$clustering[1])
stopifnot(identical(predict(nsom), nsom$clustering))

rsom <- trainSOM(dissim.lesmis, type="relational", maxit=10, scaling= "none")
stopifnot(identical(predict(rsom, dissim.lesmis), rsom$clustering))
stopifnot(predict(rsom, dissim.lesmis[1,])==rsom$clustering[1])
stopifnot(identical(predict(nsom), nsom$clustering))

rsom <- trainSOM(dissim.lesmis, type="relational", maxit=10, scaling= "sd")
stopifnot(identical(predict(rsom, dissim.lesmis), rsom$clustering))
stopifnot(predict(rsom, dissim.lesmis[1,])==rsom$clustering[1])
stopifnot(identical(predict(nsom), nsom$clustering))

rsom <- trainSOM(dissim.lesmis, type="relational", maxit=10, scaling= "max")
stopifnot(identical(predict(rsom, dissim.lesmis), rsom$clustering))
stopifnot(predict(rsom, dissim.lesmis[1,])==rsom$clustering[1])
stopifnot(identical(predict(nsom), nsom$clustering))

rsom <- trainSOM(dissim.lesmis, type="relational", maxit=10,
                 scaling= "frobenius")
stopifnot(identical(predict(rsom, dissim.lesmis), rsom$clustering))
stopifnot(predict(rsom, dissim.lesmis[1,])==rsom$clustering[1])
stopifnot(identical(predict(nsom), nsom$clustering))

rsom <- trainSOM(dissim.lesmis, type="relational", maxit=10, scaling= "cosine")
sim.matrix <- -.5* (diag(1, nrow(rsom$data))- 1/nrow(rsom$data)) %*% 
  rsom$data %*% (diag(1, nrow(rsom$data))-1/nrow(rsom$data))
sim.matrix <- diag(round(sim.matrix, 10))
stopifnot(identical(predict(rsom, dissim.lesmis, auto.sim=sim.matrix),
                    rsom$clustering))
stopifnot(predict(rsom, dissim.lesmis[1,], auto.sim=sim.matrix[1])==
            rsom$clustering[1])
stopifnot(identical(predict(nsom), nsom$clustering))

korr <- trainSOM(presidentielles2002, type= "korresp", maxit= 10)
stopifnot(identical(predict(nsom), nsom$clustering))
