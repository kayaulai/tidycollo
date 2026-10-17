test_that("type frequencies, their logs and hapax proportions are consistent", {
  res = runMeasures(smallCorpus(), list("prod"))
  for(lemma in unique(res$x)){
    rows = res[res$x == lemma, ]
    expect_equal(rows$tf_x[1], nrow(rows))
    expect_equal(rows$log_tf_x[1], log2(nrow(rows)))
    expect_equal(rows$p_hapax_cond_x[1], sum(rows$n == 1) / sum(rows$n))
  }
  for(lemma in unique(res$y)){
    rows = res[res$y == lemma, ]
    expect_equal(rows$tf_y[1], nrow(rows))
    expect_equal(rows$log_tf_y[1], log2(nrow(rows)))
    expect_equal(rows$p_hapax_cond_y[1], sum(rows$n == 1) / sum(rows$n))
  }
})

test_that("entropy matches the reference Shannon entropy", {
  res = runMeasures(smallCorpus(), list("prod"))
  rows = res[res$x == "p", ]
  p = rows$n / sum(rows$n)
  expect_equal(rows$h_x_cond_y[1], -sum(p * log2(p)))
})
