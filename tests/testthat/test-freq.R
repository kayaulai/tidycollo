#Tokens with an NA lemma in some rows, as with an empty slot
freqCorpus = function(){
  set.seed(7)
  data = dplyr::tibble(doc = sample(paste0("d", 1:8), 300, replace = TRUE),
                       x = sample(c("p", "q", "r", "s"), 300, TRUE, prob = c(.4, .3, .2, .1)),
                       y = sample(c("u", "v", "w"), 300, TRUE))
  data$x[1:15] = NA
  data
}

runFreq = function(data, ...){
  suppressMessages(capture.output(
    res <- getMeasures(data, c("x", "y"), c("x", "y"), "doc", ...)))
  res$results
}

test_that("a frequency table gives the same results as the tokens", {
  data = freqCorpus()
  expected = runFreq(data)
  freqs = dplyr::count(data, doc, x, y, name = "freq")
  expect_equal(runFreq(freqs[sample(nrow(freqs)), ], freq_col = "freq"), expected)

  #One combination and doc split over several rows
  split = dplyr::bind_rows(dplyr::mutate(freqs, freq = freq - 1L), dplyr::mutate(freqs, freq = 1L))
  expect_equal(runFreq(split, freq_col = "freq"), expected)
})

test_that("rows with frequency 0 are dropped with a message", {
  data = freqCorpus()
  freqs = dplyr::count(data, doc, x, y, name = "freq")
  zeros = dplyr::bind_rows(freqs, dplyr::tibble(doc = c("d1", "d2"), x = c("p", "zz"), y = c("zz", "zz"), freq = 0L))
  msgs = capture_messages(capture.output(getMeasures(zeros, c("x", "y"), c("x", "y"), "doc", freq_col = "freq")))
  expect_true(any(grepl("Dropped 2 combination-document rows with frequency 0", msgs)))
  expect_equal(runFreq(zeros, freq_col = "freq"), runFreq(data))
})

test_that("countsFromWide gives the same results as the tokens", {
  data = freqCorpus()
  freqs = dplyr::count(data, doc, x, y, name = "freq")
  docs = sort(unique(freqs$doc))
  wide = dplyr::distinct(freqs, x, y)
  for(d in docs){
    f = freqs[freqs$doc == d, ]
    wide[[d]] = f$freq[match(paste(wide$x, wide$y), paste(f$x, f$y))]
    wide[[d]][is.na(wide[[d]])] = 0L
  }
  long = countsFromWide(wide, c("x", "y"))
  expect_equal(nrow(long), nrow(wide) * length(docs))
  expect_equal(runFreq(long, freq_col = "freq"), runFreq(data))
})

test_that("frequencies must be whole numbers, at least 0, with no NAs", {
  freqs = dplyr::count(freqCorpus(), doc, x, y, name = "freq")
  for(bad in list(.5, -1, NA)){
    wrong = freqs
    wrong$freq[1] = bad
    expect_error(runFreq(wrong, freq_col = "freq"), "whole numbers")
  }
  expect_error(runFreq(freqs, freq_col = "nope"), "not in data")
})

test_that("keep_zero_rows keeps unattested combinations without changing attested ones", {
  data = freqCorpus()
  freqs = dplyr::count(data, doc, x, y, name = "freq")
  #Every combination in every doc, plus a new lemma t that is never attested
  grid = expand.grid(doc = unique(freqs$doc), x = c(unique(freqs$x), "t"), y = unique(freqs$y),
                     stringsAsFactors = FALSE)
  full = dplyr::left_join(dplyr::as_tibble(grid), freqs, by = c("doc", "x", "y")) %>%
    dplyr::mutate(freq = dplyr::coalesce(freq, 0L))

  expected = runFreq(data)
  kept = runFreq(full, freq_col = "freq", keep_zero_rows = TRUE)
  expect_equal(nrow(kept), nrow(grid) / length(unique(freqs$doc)))
  attested = dplyr::semi_join(kept, expected, by = c("x", "y"))
  #Ranks are among all listed combinations, so they can differ
  cols = setdiff(names(expected), grep("^rank", names(expected), value = TRUE))
  expect_equal(attested[cols], expected[cols])

  unattested = dplyr::filter(kept, x == "t")
  expect_true(all(unattested$n == 0))
  expect_true(all(is.na(unattested$DP) & is.na(unattested$range)))
  expect_true(all(is.na(unattested$h_x_cond_y)))
  expect_true(all(unattested$tf_x == 0))
})

test_that("keep_zero_rows counts a document listed only with 0 as a document", {
  data = freqCorpus()
  freqs = dplyr::count(data, doc, x, y, name = "freq")
  empty = dplyr::bind_rows(freqs, dplyr::tibble(doc = "empty", x = "p", y = "u", freq = 0L))
  kept = runFreq(empty, freq_col = "freq", keep_zero_rows = TRUE)
  dropped = runFreq(empty, freq_col = "freq")
  nd = length(unique(freqs$doc))
  expect_equal(kept$range, dropped$range * nd / (nd + 1))
  expect_equal(kept$pmi, dropped$pmi)
})

test_that("the bootstrap gives the same results from frequencies as from tokens", {
  data = freqCorpus()
  freqs = dplyr::count(data, doc, x, y, name = "freq")
  boot = function(...){
    set.seed(1)
    suppressMessages(capture.output(
      res <- getMeasuresBoot(..., c("x", "y"), c("x", "y"), "doc", list("bam_str"), R = 20, type = "perc")))
    res
  }
  fromTokens = boot(data)
  fromFreqs = boot(freqs, freq_col = "freq")
  expect_equal(fromFreqs$intervals, fromTokens$intervals)
  expect_equal(fromFreqs$presence, fromTokens$presence)
})
