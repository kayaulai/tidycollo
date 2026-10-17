test_that("dispersion measures are numeric", {
  res = runMeasures(smallCorpus(), list("disp"))
  for(col in c("range", "DP", "DP_norm", "DP_nofreq", "kld_disp", "rosengren_s", "DC"))
    expect_type(res[[col]], "double")
})

test_that("dispersion uses token counts, not combination counts, as doc sizes", {
  res = runMeasures(repeatCorpus(), list("disp"))
  row = res[res$x == "p" & res$y == "z", ]
  #(p, z): 4 tokens in doc A and 1 in doc B; both docs have 5 tokens
  expect_equal(row$DP, refDP(c(4, 1), c(5, 5)))
})

test_that("dispersion measures match reference formulas with unequal doc sizes", {
  data = smallCorpus()
  res = runMeasures(data, list("disp"))
  docs = c("A", "B", "C", "D")
  sizes = as.numeric(table(factor(data$doc, levels = docs)))
  s = sizes / sum(sizes)
  nd = length(docs)

  for(i in seq_len(nrow(res))){
    row = res[i, ]
    v = as.numeric(table(factor(data$doc[data$x == row$x & data$y == row$y], levels = docs)))
    f = sum(v)
    p = v / f
    dp = sum(abs(p - s)) / 2

    expect_equal(row$range, mean(v > 0))
    expect_equal(row$DP, dp)
    expect_equal(row$DP_norm, dp / (1 - min(s)))
    expect_equal(row$kld_disp, sum(refXLogXOverY(p, s, 2)))
    expect_equal(row$kld_disp_norm, 1 - exp(-sum(refXLogXOverY(p, s, 2))))
    expect_equal(row$chisq_disp, sum((v - f * s)^2 / (f * s)))
    expect_equal(row$rosengren_s, sum(sqrt(s * v))^2 / f)
    expect_equal(row$rosengren_af, sum(sqrt(s * v))^2)
    expect_equal(row$DC, (sum(sqrt(v)) / nd)^2 / (f / nd))
    expect_equal(row$engvall_af, f * mean(v > 0))
    expect_equal(row$juilland_u, row$juilland_d * f)
  }
  expect_true(all(is.finite(res$kld_disp)))
})

test_that("DP floor is the exact minimum over all allocations", {
  set.seed(42)
  for(trial in 1:40){
    sizes = sample(1:8, sample(2:4, 1), replace = TRUE)
    f = sample(1:min(sum(sizes), 7), 1)
    allocs = refAllocations(f, sizes)
    dps = apply(allocs, 1, refDP, sizes = sizes)
    expect_equal(getDPFloor(f, sizes / sum(sizes)), min(dps))
  }
})

test_that("DP ceiling fills the smallest docs first with exactly f tokens", {
  #Note: this is not always the true maximum DP under doc-size constraints,
  #e.g. sizes (1, 2, 10), f = 2: (0, 2, 0) beats (1, 1, 0).
  refSmallestFirst = function(f, sizes){
    x = numeric(length(sizes))
    left = f
    for(i in order(sizes)){
      x[i] = min(sizes[i], left)
      left = left - x[i]
    }
    refDP(x, sizes)
  }
  sizes = c(2, 50, 50)
  expect_equal(getDPCeiling(3, sizes / sum(sizes), sizes), refDP(c(2, 1, 0), sizes))
  set.seed(42)
  for(trial in 1:40){
    sizes = sample(1:20, sample(2:6, 1), replace = TRUE)
    f = sample(1:sum(sizes), 1)
    expect_equal(getDPCeiling(f, sizes / sum(sizes), sizes), refSmallestFirst(f, sizes))
  }
})

test_that("DP_nofreq rescales DP between the floor and ceiling for its frequency", {
  data = smallCorpus()
  res = runMeasures(data, list("disp"))
  sizes = as.numeric(table(factor(data$doc, levels = c("A", "B", "C", "D"))))
  s = sizes / sum(sizes)
  for(i in seq_len(nrow(res))){
    row = res[i, ]
    lo = getDPFloor(row$n, s)
    hi = getDPCeiling(row$n, s, sizes)
    expect_equal(row$DP_nofreq, (row$DP - lo) / (hi - lo))
  }
})
