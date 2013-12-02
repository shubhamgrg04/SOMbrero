## Check that predict.somRes works fine ("numeric" and "relational" cases)

library(SOMbrero)

nsom <- trainSOM(iris[1:30,1:4], maxit=10)
stopifnot(identical(predict(nsom, iris[1:30,1:4]), nsom$clustering))
stopifnot(predict(nsom, iris[1,1:4])==nsom$clustering[1])

data(lesmis)
rsom <- trainSOM(dissim.lesmis, type="relational", maxit=10)
stopifnot(identical(predict(rsom, dissim.lesmis), rsom$clustering))
stopifnot(predict(rsom, dissim.lesmis[1,])==rsom$clustering[1])