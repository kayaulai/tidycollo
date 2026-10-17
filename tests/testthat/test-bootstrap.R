bootCorpus = function(){
  #Combination (a, u) is ten times as frequent as (b, v) in every doc
  docs = paste0("t", 1:6)
  dplyr::bind_rows(lapply(docs, function(d){
    dplyr::tibble(text = d,
                  x = c(rep("a", 10), "b"),
                  y = c(rep("u", 10), "v"))
  }))
}

runBoot = function(data, ...){
  suppressMessages(capture.output(
    res <- getMeasuresBoot(data, c("x", "y"), c("x", "y"), "text", list("bam_str"), ...)))
  res
}

test_that("bootstrap quantiles are stored per combination and work with any doc_id", {
  set.seed(1)
  res = runBoot(bootCorpus(), R = 5)
  med = res[[3]]
  expect_equal(med$n[med$x == "a"], 60)
  expect_equal(med$n[med$x == "b"], 6)
})

test_that("a doc drawn twice in a resample counts as two docs", {
  set.seed(1)
  data = bootCorpus() %>% dplyr::mutate(text = factor(text))
  resamples = suppressMessages(capture.output(
    r <- stratBoot(data, "text", 3, function(data_new, ...) list(results = data_new),
                   1, c("x", "y"), c("x", "y"), list())))
  for(data_new in r){
    expect_equal(dplyr::n_distinct(data_new$text), 6)
    expect_equal(nrow(data_new), 66)
  }
})
