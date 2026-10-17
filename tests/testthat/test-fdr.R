test_that("Benjamini-Yekutieli rejections match p.adjust", {
  set.seed(7)
  cases = list(c(0.001, 0.013, 0.0175, 0.5), #step-up rejects p = .013 despite it missing its own cutoff
               c(0.03, 0.001, 0.001, 0.5),
               c(runif(20, 0, 0.01), runif(30)))
  for(p in cases){
    by = getBenjaminiYekutieli(p)
    expect_equal(by$reject, p.adjust(p, "BY") <= .05)
  }
})

test_that("Benjamini-Yekutieli critical values line up with the original p-value order", {
  p = c(0.5, 0.001, 0.02)
  cm = sum(1 / 1:3)
  expect_equal(getBenjaminiYekutieli(p)$cv, c(3, 1, 2) / (3 * cm) * .05)
})
