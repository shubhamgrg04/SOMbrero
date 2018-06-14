context("Check the equivalence between results given by fast variant
        and slow variant of trainSOM function")

test_that("'fast' and 'slow' variants give identical results on USArrests data", {
  data(USArrests)
  set.seed(557)
  fsom <- my.som <- trainSOM(x.data=USArrests, variant="fast", dimension=c(5,8), nb.save=10, maxit=200, scaling="none", radius.type="gaussian", verbose="false")

  set.seed(557)
  ssom <- my.som <- trainSOM(x.data=USArrests, variant="slow", dimension=c(5,8), nb.save=10, maxit=200, scaling="none", radius.type="gaussian", verbose="false")

  expect_equal(fsom$clustering, ssom$clustering)
})

test_that("'fast' and 'slow' variants give identical results on presidentielles2002 data", {
  data(presidentielles2002)
  set.seed(557)
  fsom <-  trainSOM(x.data=presidentielles2002, variant="fast", dimension=c(8,8), type="korresp", scaling="chi2", nb.save=10, radius.type="letremy", verbose="true")
  
  set.seed(557)
  ssom <- trainSOM(x.data=presidentielles2002, variant="slow", dimension=c(8,8), type="korresp", scaling="chi2", nb.save=10, radius.type="letremy", verbose="true")
  
  expect_equal(fsom$clustering, ssom$clustering)
})