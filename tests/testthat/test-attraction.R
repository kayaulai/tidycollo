test_that("contingency-table attraction measures match reference values", {
  res = runMeasures(smallCorpus())
  N = sum(res$n)
  for(i in seq_len(nrow(res))){
    row = res[i, ]
    tab = refTable(row$n, row$f_x, row$f_y, N)
    o = tab$obs
    e = tab$exp
    expect_equal(row$ef, e[1])
    expect_equal(row$pmi, log2(o[1] / e[1]))
    expect_equal(row$dice, 2 * o[1] / (row$f_x + row$f_y))
    expect_equal(row$res_pearson, (o[1] - e[1]) / sqrt(e[1]))
    expect_equal(row$t, (o[1] - e[1]) / sqrt(o[1]))
    expect_equal(row$chisq, sum((o - e)^2 / e))
    expect_equal(row$gsq, 2 * sum(refXLogXOverY(o, e)))
    expect_equal(row$lor, log((o[1] + .5) * (o[4] + .5) / ((o[2] + .5) * (o[3] + .5))))
  }
})

test_that("G^2 agrees with a Poisson GLM deviance and stays finite at large counts", {
  big = dplyr::tibble(doc = "A",
                      x = c(rep("a", 3000), rep("b", 3000)),
                      y = c(rep("u", 1800), rep("v", 1200), rep("u", 200), rep("v", 2800)))
  res = runMeasures(big, list("bam_test_no_fisher"))
  expect_true(all(is.finite(res$gsq)))

  row = res[res$x == "a" & res$y == "u", ]
  tab = matrix(c(1800, 1200, 200, 2800), 2, byrow = TRUE)
  cells = data.frame(count = as.vector(tab), r = factor(c(1, 2, 1, 2)), c = factor(c(1, 1, 2, 2)))
  dev = glm(count ~ r + c, family = poisson, data = cells)$deviance
  expect_equal(row$gsq, dev)
})

test_that("Fisher test p-values match fisher.test on the 2x2 table", {
  res = runMeasures(smallCorpus(), list("bam_test"))
  row = res[res$x == "q" & res$y == "z", ]
  o = refTable(row$n, row$f_x, row$f_y, sum(res$n))$obs
  tab = matrix(o, 2, byrow = TRUE)
  expect_equal(row$p_fisher_yates, fisher.test(tab)$p.value)
  expect_equal(row$p_fisher_yates_attract, fisher.test(tab, alternative = "greater")$p.value)
})

test_that("conditional probabilities, delta P and surprisal follow the slot_1 | slot_2 convention", {
  res = runMeasures(smallCorpus())
  row = res[res$x == "q" & res$y == "z", ]
  N = sum(res$n)
  expect_equal(row$p_x_cond_y, row$n / row$f_y)
  expect_equal(row$p_y_cond_x, row$n / row$f_x)
  expect_equal(row$Dp_x_on_y, row$n / row$f_y - (row$f_x - row$n) / (N - row$f_y))
  expect_equal(row$Dp_y_on_x, row$n / row$f_x - (row$f_y - row$n) / (N - row$f_x))
  expect_equal(row$surp_x, -log2(row$p_x_cond_y))
  expect_equal(row$surp_y, -log2(row$p_y_cond_x))
})

test_that("attraction KLDs match the Bernoulli KL divergence and handle p = 1", {
  refKLD = function(p, q) sum(refXLogXOverY(c(p, 1 - p), c(q, 1 - q), 2))
  for(data in list(smallCorpus(), repeatCorpus())){
    res = runMeasures(data, list("uam_str"))
    N = sum(res$n)
    for(i in seq_len(nrow(res))){
      row = res[i, ]
      expect_equal(row$kld_x_cond_y, refKLD(row$n / row$f_y, row$f_x / N))
      expect_equal(row$kld_y_cond_x, refKLD(row$n / row$f_x, row$f_y / N))
    }
    expect_false(any(is.na(res$kld_x_cond_y) | is.na(res$kld_y_cond_x)))
    expect_true(all(res$kld_x_cond_y >= 0 & res$kld_y_cond_x >= 0))
  }
})
