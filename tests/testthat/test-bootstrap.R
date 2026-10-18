bootCorpus = function(){
  #Combination (a, u) is ten times as frequent as (b, v) in every doc
  docs = paste0("t", 1:6)
  dplyr::bind_rows(lapply(docs, function(d){
    dplyr::tibble(text = d,
                  x = c(rep("a", 10), "b"),
                  y = c(rep("u", 10), "v"))
  }))
}

#Docs of different sizes and makeup; (s, w) only occurs in doc d1
mixedCorpus = function(){
  set.seed(2)
  docs = dplyr::bind_rows(lapply(paste0("d", 1:12), function(d){
    k = sample(5:15, 1)
    dplyr::tibble(text = d,
                  x = sample(c("p", "q", "r"), k, TRUE, prob = c(.5, .3, .2)),
                  y = sample(c("u", "v"), k, TRUE))
  }))
  dplyr::bind_rows(docs, dplyr::tibble(text = "d1", x = "s", y = "w"))
}

runBoot = function(data, ...){
  suppressMessages(capture.output(
    res <- getMeasuresBoot(data, c("x", "y"), c("x", "y"), "text", ...)))
  res
}

test_that("bootstrap quantiles are stored per combination and work with any doc_id", {
  set.seed(1)
  res = runBoot(bootCorpus(), measures = list("bam_str"), R = 5, type = "perc")
  med = res$quantiles[["50%"]]
  expect_equal(med$n[med$x == "a"], 60)
  expect_equal(med$n[med$x == "b"], 6)
})

test_that("a doc drawn twice in a resample counts as two docs", {
  set.seed(1)
  data = bootCorpus()
  rows_bydoc = split(seq_len(nrow(data)), factor(data$text, levels = unique(data$text)))
  for(i in 1:3){
    drawn = drawDocs(rows_bydoc, sample(6, 6, replace = TRUE))
    expect_equal(dplyr::n_distinct(drawn$docs), 6)
    expect_equal(length(drawn$rows), 66)
  }
})

test_that("combinations missing from a resample are dropped, not counted as 0", {
  set.seed(1)
  res = runBoot(mixedCorpus(), measures = list("bam_str"), R = 200, type = "perc",
                retainResamples = TRUE)
  rare = res$combos$x == "s"
  present = res$replicates$present[rare, ]
  #(s, w) is present exactly when d1 is drawn
  expect_true(any(!present) && any(present))
  expect_equal(res$presence$prop_present[res$presence$x == "s"], mean(present))
  expect_true(all(is.na(res$replicates$t[rare, "n", !present])))
  expect_equal(res$quantiles[["0%"]]$n[res$quantiles[["0%"]]$x == "s"],
               min(res$replicates$t[rare, "n", present]))
  expect_gte(res$quantiles[["0%"]]$n[res$quantiles[["0%"]]$x == "s"], 1)
  ci = dplyr::filter(res$intervals, x == "s", measure == "n")
  expect_equal(ci$n_resamples, sum(present))
})

test_that("percentile intervals and quantiles match quantile()", {
  set.seed(1)
  res = runBoot(mixedCorpus(), measures = list("bam_str"), R = 100, type = "perc",
                conf = c(.9, .95), retainResamples = TRUE)
  i = which(res$combos$x == "p" & res$combos$y == "u")
  t = res$replicates$t[i, "pmi", ]
  ci = dplyr::filter(res$intervals, x == "p", y == "u", measure == "pmi")
  expect_equal(ci$lower, unname(c(quantile(t, .05), quantile(t, .025))))
  expect_equal(ci$upper, unname(c(quantile(t, .95), quantile(t, .975))))
  expect_equal(ci$estimate, rep(res$original$results$pmi[i], 2))
  expect_equal(res$quantiles[["50%"]]$pmi[i], unname(median(t)))
})

test_that("BCa intervals agree with boot::boot.ci", {
  skip_if_not_installed("boot")
  set.seed(3)
  res = runBoot(mixedCorpus(), measures = list("bam_str"), R = 999,
                statistics = list(dice_pu = function(r) r$dice[r$x == "p" & r$y == "u"]),
                retainResamples = TRUE)
  reps = res$replicates
  bootCI = function(t0, t, jk){
    n = length(jk)
    out = list(t0 = t0, t = matrix(t), R = length(t), sim = "ordinary", stype = "i",
               strata = rep(1, n), weights = rep(1 / n, n), data = seq_len(n), call = quote(boot()))
    class(out) = "boot"
    boot::boot.ci(out, conf = .95, type = "bca", L = (n - 1) * (mean(jk) - jk))$bca[4:5]
  }
  ours = dplyr::filter(res$statistics, type == "bca")
  expected = bootCI(reps$stats_t0[1], reps$stats_t[1, ], reps$stats_jk[1, ])
  #boot interpolates quantiles on the normal scale, so allow a small difference
  expect_equal(c(ours$lower, ours$upper), expected, tolerance = .01)

  i = which(res$combos$x == "p" & res$combos$y == "u")
  ours = dplyr::filter(res$intervals, x == "p", y == "u", measure == "pmi", type == "bca")
  expected = bootCI(reps$t0[i, "pmi"], reps$t[i, "pmi", ], reps$jk[i, "pmi", ])
  expect_equal(c(ours$lower, ours$upper), expected, tolerance = .01)
})

test_that("intervals can be recalculated from retained replicates", {
  set.seed(1)
  res = runBoot(mixedCorpus(), measures = list("bam_str"), R = 100, conf = .95,
                statistics = list(n_pos = function(r) sum(r$pmi > 0)),
                retainResamples = TRUE)
  again = getBootIntervals(res, conf = .95)
  expect_equal(again$intervals, res$intervals)
  expect_equal(again$statistics, res$statistics)
  narrower = getBootIntervals(res, type = "perc", conf = .5)$statistics
  wider = dplyr::filter(res$statistics, type == "perc")
  expect_lte(wider$lower, narrower$lower)
  expect_gte(wider$upper, narrower$upper)
})

test_that("BCa is NA when the original value is outside every resample", {
  t = matrix(c(2, 3, 4, 5), nrow = 1)
  ci = getIntervalsFromMatrix(1, t, matrix(c(1, 2, 3), nrow = 1), "bca", .95)
  expect_true(is.infinite(ci$z0))
  expect_true(is.na(ci$lower) && is.na(ci$upper))
})

test_that("parallel resampling gives the same results as sequential", {
  skip_if_not_installed("future.apply")
  data = mixedCorpus()
  set.seed(1)
  sequential = runBoot(data, measures = list("bam_str"), R = 20, type = "perc")
  oplan = future::plan("sequential")
  on.exit(future::plan(oplan))
  set.seed(1)
  parallel = runBoot(data, measures = list("bam_str"), R = 20, type = "perc", parallel = TRUE)
  expect_equal(parallel$intervals, sequential$intervals)
  expect_equal(parallel$presence, sequential$presence)
})

test_that("resamples built from counts give the same measures as resampled tokens", {
  data = mixedCorpus()
  suppressMessages(capture.output(orig <- getMeasures(data, c("x", "y"), c("x", "y"), "text")))
  counts = orig$counts_by_doc_combo
  docs = unique(data$text)
  drawn = c(1, 1, 3, 5, 5, 5, 12)
  rows_bydoc = split(seq_len(nrow(counts)), factor(counts$text, levels = docs))
  fromCounts = drawDocs(rows_bydoc, drawn)
  combos = orig$results[c("x", "y")]
  rep_counts = getReplicateFromCounts(counts, fromCounts$rows, fromCounts$docs, combos, c("x", "y"),
                                      c("x", "y"), "text", list("bam_str", "prod", "disp"), c("n", "pmi", "DP"), list())

  tokens = dplyr::bind_rows(lapply(seq_along(drawn), function(k){
    dplyr::filter(data, text == docs[drawn[k]]) %>% dplyr::mutate(text = paste0(text, "_", k))
  }))
  suppressMessages(capture.output(direct <- getMeasures(tokens, c("x", "y"), c("x", "y"), "text",
                                                       list("bam_str", "prod", "disp"))$results))
  expect_equal(rep_counts$combo, matchCombos(direct, combos, c("x", "y")))
  expect_equal(unname(rep_counts$values), unname(as.matrix(direct[c("n", "pmi", "DP")])))
})

test_that("jackknifeSE matches the usual formula, with weights standing for copies", {
  x = matrix(c(1, 4, 2, 8, 5, 7), nrow = 2)
  refSE = function(v) sqrt((length(v) - 1) / length(v) * sum((v - mean(v))^2))
  expect_equal(jackknifeSE(x), c(refSE(c(1, 2, 5)), refSE(c(4, 8, 7))))
  expect_equal(jackknifeSE(x, c(2, 1, 1)), c(refSE(c(1, 1, 2, 5)), refSE(c(4, 4, 8, 7))))
  x[1, 2] = NA
  expect_equal(jackknifeSE(x)[1], refSE(c(1, 5)))
  arr = array(c(x, x), c(2, 1, 3), dimnames = list(NULL, "m", NULL))
  expect_equal(dim(jackknifeSE(arr)), c(2, 1))
})

test_that("the inner jackknife of a resample equals leaving out each drawn doc in turn", {
  data = mixedCorpus()
  suppressMessages(capture.output(orig <- getMeasures(data, c("x", "y"), c("x", "y"), "text")))
  counts = orig$counts_by_doc_combo
  rows_bydoc = split(seq_len(nrow(counts)), factor(counts$text, levels = unique(data$text)))
  combos = orig$results[c("x", "y")]
  getReplicate = function(rows, docs_new){
    getReplicateFromCounts(counts, rows, docs_new, combos, c("x", "y"), c("x", "y"), "text",
                           list("bam_str"), c("n", "pmi"), list(total = function(r) sum(r$n)))
  }
  drawn = c(2, 2, 2, 4, 7, 9, 9, 11, 1, 3, 6, 6)
  se = getInnerSE(drawn, rows_bydoc, getReplicate, nrow(combos), c("n", "pmi"),
                  list(total = sum), "jackknife")

  #Brute force: leave out each of the 12 positions
  reps = collectReplicates(length(drawn), nrow(combos), c("n", "pmi"), list(total = sum), function(k){
    d = drawDocs(rows_bydoc, drawn[-k])
    getReplicate(d$rows, d$docs)
  }, progress = FALSE)
  expect_equal(se$se, jackknifeSE(reps[[1]]))
  expect_equal(se$stats_se, jackknifeSE(reps[[3]]))
})

test_that("studentized intervals follow the bootstrap-t formula", {
  set.seed(1)
  res = runBoot(mixedCorpus(), measures = list("bam_str"), R = 100, type = c("perc", "stud"),
                statistics = list(dice_pu = function(r) r$dice[r$x == "p" & r$y == "u"]),
                retainResamples = TRUE)
  reps = res$replicates
  i = which(res$combos$x == "p" & res$combos$y == "u")
  t0 = reps$t0[i, "pmi"]
  t_stud = (reps$t[i, "pmi", ] - t0) / reps$se_t[i, "pmi", ]
  se_jack = jackknifeSE(matrix(reps$jk[i, "pmi", ], nrow = 1))
  ci = dplyr::filter(res$intervals, x == "p", y == "u", measure == "pmi", type == "stud")
  expect_equal(ci$lower, unname(t0 - quantile(t_stud, .975, na.rm = TRUE) * se_jack))
  expect_equal(ci$upper, unname(t0 - quantile(t_stud, .025, na.rm = TRUE) * se_jack))
  expect_equal(ci$se_jack, se_jack)
  expect_equal(ci$n_stud, sum(is.finite(t_stud)))
  expect_lt(ci$lower, ci$estimate)
  expect_gt(ci$upper, ci$estimate)

  stat = dplyr::filter(res$statistics, type == "stud")
  expect_true(is.finite(stat$lower) && is.finite(stat$upper))
  expect_equal(getBootIntervals(res, c("perc", "stud"))$intervals, res$intervals)
  #The jackknife of the original data was run for the scale, so BCa is available too
  expect_no_error(getBootIntervals(res, "bca"))
})

test_that("studentized intervals can use an inner bootstrap", {
  set.seed(1)
  res = runBoot(mixedCorpus(), measures = list("bam_str"), R = 20, type = "stud",
                inner = "bootstrap", R_inner = 10, retainResamples = TRUE)
  reps = res$replicates
  expect_null(reps$jk)
  i = which(res$combos$x == "p" & res$combos$y == "u")
  t0 = reps$t0[i, "pmi"]
  t_stud = (reps$t[i, "pmi", ] - t0) / reps$se_t[i, "pmi", ]
  ci = dplyr::filter(res$intervals, x == "p", y == "u", measure == "pmi")
  expect_equal(ci$lower, unname(t0 - quantile(t_stud, .975, na.rm = TRUE) * sd(reps$t[i, "pmi", ])))
  expect_error(getBootIntervals(res, "bca"), "jackknife")
})
