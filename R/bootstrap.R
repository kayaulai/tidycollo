#Bootstrap intervals for co-occurrence measures, resampling documents.

#' Bootstrap co-occurrence measures by resampling documents.
#'
#' Documents are drawn with replacement, and [getMeasures()] is run again on
#' each resample. A document drawn twice counts as two documents.
#'
#' A combination that does not occur in a resample is left out of that
#' resample, instead of being counted as zero. This is intended: a real sample
#' can also miss rare combinations, so the interval of each combination
#' describes its measures given that it was observed. `presence` gives the
#' proportion of resamples in which each combination was observed.
#' Replicates where a measure is `NA` or `NaN` are also left out, so
#' `n_resamples` in `intervals` can be lower than `n_present`.
#'
#' Interval types:
#' * `"perc"`: percentile interval.
#' * `"bca"`: bias-corrected and accelerated interval. The acceleration comes
#'   from a leave-one-document-out jackknife on the original data, which costs
#'   one more [getMeasures()] call per document. The bias correction counts ties
#'   with the original value as half, which matters for discrete measures and
#'   for counts in `statistics`. The interval is `NA` when the original value is
#'   below or above every resample (the bias correction is infinite then). When
#'   all jackknife values are equal, the acceleration is taken as 0.
#' * `"stud"`: studentized (bootstrap-t) interval,
#'   \eqn{[\hat\theta - t^*_{1-\alpha/2} \hat{se}, \hat\theta - t^*_{\alpha/2} \hat{se}]},
#'   where \eqn{t^* = (\theta^* - \hat\theta) / se^*} for each resample. With
#'   `inner = "jackknife"`, \eqn{se^*} is a leave-one-document-out jackknife on
#'   the resample and \eqn{\hat{se}} the same jackknife on the original data
#'   (`se_jack`). This costs one [getMeasures()] call per distinct document
#'   in each resample: leaving out either copy of a document drawn twice gives
#'   the same result, so it is computed once and weighted twice. With
#'   `inner = "bootstrap"`, \eqn{se^*} is the SD over `R_inner` resamples of the
#'   resample and \eqn{\hat{se}} the SD over the outer resamples (`se`). This costs
#'   `R_inner` calls per resample. Resamples where \eqn{se^*} is 0 or missing give no
#'   \eqn{t^*} and are left out (`n_stud` counts the rest). This happens often for
#'   discrete measures.
#'
#' In both jackknifes, a combination that disappears when a document is left
#' out is left out of that jackknife sample, and the jackknife uses the
#' remaining samples.
#'
#' With `freq_col`, `data` can be a table of frequencies, which gives the same
#' results as the tokens it was counted from (under the same seed). With
#' `keep_zero_rows = TRUE`, a combination listed with frequency 0 in a document
#' stays in every resample that draws that document, with `n = 0`, instead of
#' being dropped as unobserved.
#'
#' @inheritParams getMeasures
#' @param R Number of resamples.
#' @param type Interval types, any of `"perc"`, `"bca"` and `"stud"`. See
#'   Details.
#' @param conf Confidence levels, e.g. `c(.9, .95)`.
#' @param quantiles Quantiles of the bootstrap distribution of each measure,
#'   returned in the same wide format as `getMeasures()$results`.
#' @param boot_measures Names of the columns of `getMeasures()$results` to
#'   bootstrap. `NULL` (the default) takes every numeric column.
#' @param statistics A named list of functions for statistics of the whole
#'   corpus. Each takes the `results` table of [getMeasures()] and returns a
#'   single number, e.g. `function(res) sum(res$pmi > 3)`.
#' @param inner How the studentized interval estimates the SE of each resample:
#'   `"jackknife"` or `"bootstrap"`. See Details.
#' @param R_inner Number of inner resamples when `inner = "bootstrap"`.
#' @param sampleRatio Number of documents drawn in each resample, as a
#'   proportion of the number of documents. The intervals assume 1.
#' @param retainResamples Whether to return the bootstrap and jackknife
#'   replicates. They are needed by [getBootIntervals()].
#' @param parallel Whether to process the resamples in parallel with
#'   [future.apply::future_lapply()], following the plan set with
#'   [future::plan()], e.g. `future::plan("multisession")`. For `"stud"`, each
#'   worker also runs the inner loop of its resamples. All documents are drawn
#'   beforehand, so results are the same as with `parallel = FALSE` under the
#'   same seed.
#'
#' @return A list:
#'   * `intervals`: a tibble with one row per combination, measure, interval
#'     type and confidence level, i.e.
#'     (combinations × measures × `length(type)` × `length(conf)`) rows. Columns:
#'     the two lemma columns, `measure`, `type`, `conf`, the original
#'     `estimate`, the `lower` and `upper` bounds, the number of resamples
#'     with a value (`n_resamples`), the bootstrap `se` (SD of the resamples)
#'     and `bias` (mean of the resamples minus `estimate`), `se_jack` (jackknife
#'     SE on the original data, `NA` if no jackknife was run), and, for their own
#'     types only, `z0` and `acc` (BCa bias correction and acceleration) and
#'     `n_stud` (number of studentized values used).
#'   * `presence`: a tibble with one row per combination: the lemma columns,
#'     the number (`n_present`) and proportion (`prop_present`) of resamples in
#'     which it was observed.
#'   * `quantiles`: a named list (`"2.5%"`, `"50%"`, ...) with one tibble per
#'     quantile. Each has one row per combination, the lemma columns and one
#'     column per measure in `boot_measures`.
#'   * `statistics`: as `intervals`, with a `statistic` column instead of the
#'     lemma and `measure` columns, or `NULL` if `statistics` is empty.
#'   * `original`: the output of [getMeasures()] on `data`.
#'   * `combos`, `replicates`: if `retainResamples`, the combinations (a tibble
#'     of the lemma columns, in the row order of `original$results`) and the
#'     replicates (see [getBootIntervals()]).
#' @export
#'
#' @examples
#' data = dplyr::tibble(doc = rep(c("A", "B", "C", "D"), each = 4),
#'                      x_lemma = rep(c("p", "p", "q", "r"), 4),
#'                      y_lemma = rep(c("z", "y", "z", "z"), 4))
#' set.seed(1)
#' res = getMeasuresBoot(data, lemma_names = c("x_lemma", "y_lemma"),
#'                       slot_names = c("x", "y"), doc_id = "doc",
#'                       measures = list("bam_str"), R = 50, type = "perc",
#'                       statistics = list(n_pmi_pos = function(res) sum(res$pmi > 0)))
#' res$intervals
#' res$statistics
getMeasuresBoot = function(data,
                           lemma_names = c("lemma_1", "lemma_2"),
                           slot_names = c("f", "s"),
                           doc_id = "doc",
                           measures = list("bam_str", "bam_test", "uam_str", "uam_rank", "prod", "disp"),
                           R = 1000,
                           type = c("perc", "bca"),
                           conf = .95,
                           quantiles = c(0, .025, .5, .975, 1),
                           boot_measures = NULL,
                           statistics = list(),
                           inner = c("jackknife", "bootstrap"),
                           R_inner = 50,
                           sampleRatio = 1,
                           retainResamples = F,
                           parallel = F,
                           freq_col = NULL,
                           keep_zero_rows = FALSE){
  type = match.arg(type, c("perc", "bca", "stud"), several.ok = TRUE)
  inner = match.arg(inner)
  if(length(statistics) > 0 && (is.null(names(statistics)) || any(names(statistics) == ""))){
    stop("All statistics must be named.")
  }
  #Resamples are built from the counts by combination and doc, not the tokens.
  #Counted outside quietly() so that the message on zero rows is shown.
  counts = getCountsByComboDoc(data, lemma_names, doc_id, freq_col, keep_zero_rows)
  counts_docs = as.character(counts[[doc_id]])
  docs_list = unique(counts_docs)
  nDocs = length(docs_list)
  if(nDocs < 2) stop("At least two documents are needed for resampling.")
  nDraw = round(nDocs * sampleRatio)
  if(nDraw < 1) stop("sampleRatio is too small to draw any documents.")
  if(sampleRatio != 1) warning("The intervals assume sampleRatio = 1.")

  message("Performing on original data ...")
  orig = quietly(getMeasuresFromCounts(counts, lemma_names, slot_names, doc_id, measures))
  combos = orig$results %>% select(all_of(lemma_names))
  if(is.null(boot_measures)){
    boot_measures = setdiff(names(orig$results)[sapply(orig$results, is.numeric)], lemma_names)
  } else if(!all(boot_measures %in% names(orig$results))){
    stop(glue("Not in the results of getMeasures: {paste(setdiff(boot_measures, names(orig$results)), collapse = ', ')}"))
  }

  rows_bydoc = split(seq_len(nrow(counts)), factor(counts_docs, levels = docs_list))
  getReplicate = function(rows, docs_new){
    getReplicateFromCounts(counts, rows, docs_new, combos, lemma_names, slot_names,
                           doc_id, measures, boot_measures, statistics)
  }

  replicates = list(t0 = as.matrix(orig$results[boot_measures]),
                    stats_t0 = vapply(statistics, function(f) as.numeric(f(orig$results)), numeric(1)),
                    inner = inner)

  #All random draws happen here, before any (possibly parallel) computation
  draws = lapply(seq_len(R), function(i) sample(nDocs, nDraw, replace = TRUE))
  stud = "stud" %in% type
  if(stud && inner == "bootstrap"){
    draws_inner = lapply(seq_len(R), function(i){
      lapply(seq_len(R_inner), function(j) sample(nDraw, nDraw, replace = TRUE))
    })
  }

  message("Now performing resampling and calculating measures ...")
  replicates[c("t", "present", "stats_t", "se_t", "stats_se_t")] = collectReplicates(R, nrow(combos), boot_measures, statistics, function(i){
    drawn = drawDocs(rows_bydoc, draws[[i]])
    rep_i = getReplicate(drawn$rows, drawn$docs)
    if(stud){
      se = getInnerSE(draws[[i]], rows_bydoc, getReplicate, nrow(combos), boot_measures, statistics,
                      inner, if(inner == "bootstrap") draws_inner[[i]])
      rep_i$se = se$se[rep_i$combo, , drop = FALSE]
      rep_i$stats_se = se$stats_se
    }
    rep_i
  }, parallel)

  if("bca" %in% type || (stud && inner == "jackknife")){
    message("Now performing the jackknife ...")
    replicates[c("jk", "jk_present", "stats_jk")] = collectReplicates(nDocs, nrow(combos), boot_measures, statistics, function(i){
      rows = unlist(rows_bydoc[-i], use.names = FALSE)
      getReplicate(rows, counts_docs[rows])
    }, parallel)[1:3]
  }

  message("Now extracting results ...")
  intervals = getBootIntervals(list(combos = combos, replicates = replicates), type, conf)

  quantileTables = lapply(quantiles, function(q){
    table = combos
    for(m in seq_along(boot_measures)){
      table[[boot_measures[m]]] = rowQuantiles(replicates$t[, m, , drop = FALSE], q)
    }
    table
  })
  names(quantileTables) = paste0(100 * quantiles, "%")

  n_present = rowSums(replicates$present)
  result = list(intervals = intervals$intervals,
                presence = combos %>% mutate(n_present = n_present, prop_present = n_present / R),
                quantiles = quantileTables,
                statistics = intervals$statistics,
                original = orig)
  if(retainResamples){
    result$combos = combos
    result$replicates = replicates
  }
  result
}


#' Recalculate bootstrap intervals.
#'
#' Calculates intervals from the replicates kept by [getMeasuresBoot()], e.g.
#' for other confidence levels, without resampling again.
#'
#' @param boot Output of [getMeasuresBoot()] with `retainResamples = TRUE`, or
#'   any list with:
#'   * `combos`: a tibble of the two lemma columns, one row per combination
#'     (`nCombos` rows).
#'   * `replicates`: a list with
#'     - `t0`: `nCombos × nMeasures` matrix of original values, with the
#'       measure names as column names;
#'     - `t`: `nCombos × nMeasures × R` array of resample values, `NA` where
#'       the combination is missing from the resample;
#'     - `present`: `nCombos × R` logical matrix, whether each combination is in
#'       each resample;
#'     - `stats_t0`: named vector of the original `statistics` (length
#'       `nStats`, possibly 0), and `stats_t`: `nStats × R` matrix;
#'     - `jk`, `jk_present`, `stats_jk`: as `t`, `present` and `stats_t`, for the
#'       `nDocs` leave-one-document-out jackknife samples (only if the jackknife
#'       was run);
#'     - `se_t`: as `t`, the inner SE of each resample value, and `stats_se_t`:
#'       as `stats_t` (only if `"stud"` was in `type`);
#'     - `inner`: `"jackknife"` or `"bootstrap"`.
#' @param type,conf As in [getMeasuresBoot()]. BCa needs the jackknife, and
#'   the studentized interval needs `se_t`, so they are only available if they
#'   were requested in `getMeasuresBoot()`.
#'
#' @return A list with `intervals` and `statistics`, as in [getMeasuresBoot()].
#' @export
getBootIntervals = function(boot, type = c("perc", "bca"), conf = .95){
  type = match.arg(type, c("perc", "bca", "stud"), several.ok = TRUE)
  reps = boot$replicates
  if(is.null(reps)) stop("No replicates found. Run getMeasuresBoot() with retainResamples = TRUE.")
  if("bca" %in% type && is.null(reps$jk)) stop("No jackknife found. Run getMeasuresBoot() with \"bca\" in type.")
  if("stud" %in% type && is.null(reps$se_t)) stop("No inner SEs found. Run getMeasuresBoot() with \"stud\" in type.")
  boot_measures = colnames(reps$t0)
  nCombos = nrow(boot$combos)
  inner = if(is.null(reps$inner)) "jackknife" else reps$inner

  #Flatten measures so that each row is one combination-measure pair
  flatten = function(x) if(is.null(x)) NULL else matrix(x, nrow = nCombos * length(boot_measures))
  keys = boot$combos[rep(seq_len(nCombos), length(boot_measures)), ] %>%
    mutate(measure = rep(boot_measures, each = nCombos))
  intervals = getIntervalsFromMatrix(as.vector(reps$t0), flatten(reps$t), flatten(reps$jk), type, conf,
                                     flatten(reps$se_t), inner)
  intervals = bind_cols(keys[intervals$.row, ], select(intervals, -".row"))

  statistics = NULL
  if(length(reps$stats_t0) > 0){
    statistics = getIntervalsFromMatrix(reps$stats_t0, reps$stats_t, reps$stats_jk, type, conf,
                                        reps$stats_se_t, inner)
    statistics = bind_cols(tibble(statistic = names(reps$stats_t0))[statistics$.row, ],
                           select(statistics, -".row"))
  }
  list(intervals = intervals, statistics = statistics)
}


#' Evaluate an expression without its messages and printed output.
#'
#' Used to silence the messages and progress bars of [getMeasures()].
#'
#' @param expr An expression.
#' @return The value of `expr`.
#' @keywords internal
#' @noRd
quietly = function(expr){
  utils::capture.output(res <- suppressMessages(expr))
  res
}

#' Rows and doc IDs of a resample.
#'
#' Each draw gets its own doc ID, so that a doc drawn twice counts as two docs.
#'
#' @param rows_bydoc A list with one integer vector per doc (`nDocs`
#'   elements), giving the rows of that doc in the counts table. Named by doc.
#' @param drawn Integer vector of drawn docs (indices into `rows_bydoc`), of
#'   any length.
#' @return A list with two vectors of the same length (the total number of
#'   rows of the drawn docs): `rows`, the rows of the counts table, and `docs`,
#'   the new doc ID of each row (`"<doc>_<position in drawn>"`).
#' @keywords internal
#' @noRd
drawDocs = function(rows_bydoc, drawn){
  rows = rows_bydoc[drawn]
  list(rows = unlist(rows, use.names = FALSE),
       docs = rep(paste0(names(rows_bydoc)[drawn], "_", seq_along(drawn)), lengths(rows)))
}

#' Row of each combination of a table in a table of combinations.
#'
#' @param res A data frame with the two lemma columns (any number of rows,
#'   combinations may repeat).
#' @param combos A data frame with the two lemma columns, one row per
#'   combination.
#' @param lemma_names The names of the two lemma columns.
#' @return An integer vector with one element per row of `res`: the row of
#'   `combos` with the same combination, or `NA` if there is none. `NA` lemmas
#'   match `NA`.
#' @keywords internal
#' @noRd
matchCombos = function(res, combos, lemma_names){
  res %>%
    select(all_of(lemma_names)) %>%
    left_join(combos %>% select(all_of(lemma_names)) %>% mutate(.row = row_number()),
              by = lemma_names) %>%
    pull(".row")
}

#' Measures and statistics on part of the counts table.
#'
#' @param counts The counts table (`getMeasures()$counts_by_doc_combo`): one
#'   row per combination per doc, with the lemma columns, the `doc_id` column
#'   and `n`.
#' @param rows Integer vector of rows of `counts` to use (rows may repeat).
#' @param docs_new Character vector of the same length as `rows`: the doc ID
#'   of each row.
#' @param combos The combinations of the original data, as in
#'   [matchCombos()].
#' @param lemma_names,slot_names,doc_id,measures As in [getMeasures()].
#' @param boot_measures Names of the measure columns to keep (`nMeasures`).
#' @param statistics Named list of `nStats` functions, as in
#'   [getMeasuresBoot()].
#' @return A list for the `k` combinations in the selected rows:
#'   * `combo`: integer vector (length `k`), the row of each in `combos`;
#'   * `values`: `k × nMeasures` numeric matrix of measures;
#'   * `stats`: numeric vector of the statistics (length `nStats`).
#' @keywords internal
#' @noRd
getReplicateFromCounts = function(counts, rows, docs_new, combos, lemma_names, slot_names,
                                  doc_id, measures, boot_measures, statistics){
  counts_new = counts[rows, ]
  counts_new[[doc_id]] = docs_new
  res = quietly(getMeasuresFromCounts(counts_new, lemma_names, slot_names, doc_id, measures))$results
  list(combo = matchCombos(res, combos, lemma_names),
       values = as.matrix(res[boot_measures]),
       stats = vapply(statistics, function(f) as.numeric(f(res)), numeric(1)))
}

#' Inner SEs of one resample, for the studentized interval.
#'
#' @param drawn Integer vector of the docs drawn for the resample (length
#'   `nDraw`), as in [drawDocs()].
#' @param rows_bydoc As in [drawDocs()].
#' @param getReplicate A function of `rows` and `docs_new` returning a list as
#'   in [getReplicateFromCounts()].
#' @param nCombos Number of combinations in the original data.
#' @param boot_measures,statistics As in [getReplicateFromCounts()].
#' @param inner `"jackknife"` (leave out one drawn doc at a time) or
#'   `"bootstrap"` (resample the drawn docs).
#' @param draws_inner For `"bootstrap"`: a list of `R_inner` integer vectors,
#'   each indexing positions in `drawn`.
#' @return A list with `se`, an `nCombos × nMeasures` matrix of SEs (`NA` for
#'   combinations not in the resample, or in fewer than two inner samples), and
#'   `stats_se`, a vector of SEs of the statistics (length `nStats`).
#' @keywords internal
#' @noRd
getInnerSE = function(drawn, rows_bydoc, getReplicate, nCombos, boot_measures, statistics,
                      inner, draws_inner = NULL){
  if(inner == "jackknife"){
    #Leaving out either copy of a doc drawn twice gives the same sample, so each
    #distinct doc is left out once and weighted by its number of copies
    distinct = unique(drawn)
    weights = tabulate(match(drawn, distinct))
    samples = lapply(distinct, function(d) drawn[-match(d, drawn)])
  } else {
    weights = NULL
    samples = lapply(draws_inner, function(positions) drawn[positions])
  }
  reps = collectReplicates(length(samples), nCombos, boot_measures, statistics, function(j){
    drawn_j = drawDocs(rows_bydoc, samples[[j]])
    getReplicate(drawn_j$rows, drawn_j$docs)
  }, progress = FALSE)
  if(inner == "jackknife"){
    list(se = jackknifeSE(reps[[1]], weights), stats_se = jackknifeSE(reps[[3]], weights))
  } else {
    list(se = rowSDs(reps[[1]]), stats_se = rowSDs(reps[[3]]))
  }
}

#' Run replicates and align them to the original combinations.
#'
#' @param n Number of replicates.
#' @param nCombos Number of combinations in the original data.
#' @param boot_measures Names of the measures (`nMeasures`).
#' @param statistics Named list of `nStats` functions.
#' @param getReplicate A function of the replicate number `i` returning a list
#'   as in [getReplicateFromCounts()], optionally with `se` (same shape as
#'   `values`) and `stats_se` (same length as `stats`).
#' @param parallel Whether to use [future.apply::future_lapply()].
#' @param progress Whether to show a progress bar (sequential only).
#' @return An unnamed list of five elements:
#'   1. `nCombos × nMeasures × n` array of values, `NA` where a combination is
#'      missing from a replicate;
#'   2. `nCombos × n` logical matrix, whether each combination is present;
#'   3. `nStats × n` matrix of statistics;
#'   4. as 1., for `se` (`NULL` if the replicates have no `se`);
#'   5. as 3., for `stats_se` (`NULL` likewise).
#' @keywords internal
#' @noRd
collectReplicates = function(n, nCombos, boot_measures, statistics, getReplicate,
                             parallel = F, progress = TRUE){
  if(parallel){
    if(!requireNamespace("future.apply", quietly = TRUE)) stop("parallel = TRUE needs the future.apply package.")
    reps = future.apply::future_lapply(seq_len(n), getReplicate, future.seed = TRUE)
  } else {
    if(progress) pb = txtProgressBar(min = 0, max = n, initial = 0)
    reps = lapply(seq_len(n), function(i){
      if(progress) setTxtProgressBar(pb, i)
      getReplicate(i)
    })
    if(progress) close(pb)
  }

  hasSE = n > 0 && !is.null(reps[[1]]$se)
  emptyValues = function() array(NA_real_, c(nCombos, length(boot_measures), n),
                                 dimnames = list(NULL, boot_measures, NULL))
  emptyStats = function() matrix(NA_real_, length(statistics), n, dimnames = list(names(statistics), NULL))
  t = emptyValues()
  present = matrix(FALSE, nCombos, n)
  stats_t = emptyStats()
  se_t = if(hasSE) emptyValues()
  stats_se_t = if(hasSE) emptyStats()
  for(i in seq_len(n)){
    t[reps[[i]]$combo, , i] = reps[[i]]$values
    present[reps[[i]]$combo, i] = TRUE
    stats_t[, i] = reps[[i]]$stats
    if(hasSE){
      se_t[reps[[i]]$combo, , i] = reps[[i]]$se
      stats_se_t[, i] = reps[[i]]$stats_se
    }
  }
  list(t, present, stats_t, se_t, stats_se_t)
}

#' Weighted jackknife SEs, ignoring missing samples.
#'
#' For each statistic, \eqn{\sqrt{(n - 1) / n \sum_j w_j (x_j - \bar x)^2}},
#' where \eqn{\bar x} is the weighted mean and \eqn{n = \sum_j w_j}, over
#' the jackknife samples where it is not `NA`.
#'
#' @param x An array whose last dimension is the `J` jackknife samples, e.g.
#'   `nCombos × nMeasures × J` or `nStats × J`.
#' @param weights Numeric vector of length `J`: how many leave-one-out samples
#'   each sample stands for. `NULL` means all 1.
#' @return An array of SEs with the dimensions of `x` without the last one (a
#'   vector if `x` is a matrix). `NA` where fewer than two samples have a
#'   value.
#' @keywords internal
#' @noRd
jackknifeSE = function(x, weights = NULL){
  dims = dim(x)
  dimNames = dimnames(x)
  J = dims[length(dims)]
  x = matrix(x, ncol = J)
  if(is.null(weights)) weights = rep(1, J)
  w = matrix(rep(weights, each = nrow(x)), nrow(x), J)
  w[is.na(x)] = 0
  n = rowSums(w)
  xbar = rowSums(w * x, na.rm = TRUE) / n
  se = sqrt((n - 1) / n * rowSums(w * (x - xbar)^2, na.rm = TRUE))
  se[n < 2] = NA
  if(length(dims) > 2) array(se, dims[-length(dims)], dimNames[-length(dims)]) else se
}

#' SDs over the last dimension of an array, ignoring `NA`.
#'
#' @param x An array whose last dimension is the samples, e.g.
#'   `nCombos × nMeasures × R` or `nStats × R`.
#' @return An array of SDs with the dimensions of `x` without the last one (a
#'   vector if `x` is a matrix). `NA` where fewer than two samples have a
#'   value.
#' @keywords internal
#' @noRd
rowSDs = function(x){
  dims = dim(x)
  x = matrix(x, ncol = dims[length(dims)])
  n = rowSums(!is.na(x))
  xbar = rowMeans(x, na.rm = TRUE)
  sds = sqrt(rowSums((x - xbar)^2, na.rm = TRUE) / (n - 1))
  sds[n < 2] = NA
  if(length(dims) > 2) array(sds, dims[-length(dims)]) else sds
}

#' Quantiles of each row, ignoring `NA`.
#'
#' Type 7 quantiles, as in [stats::quantile()].
#'
#' @param x An array whose last dimension is the samples, e.g. `k × R` or
#'   `k × 1 × R`; it is treated as a `k × R` matrix.
#' @param p Probabilities: one per row (length `k`), or one for all rows.
#'   `NA` gives `NA`.
#' @return A numeric vector of length `k`. `NA` for rows without values.
#' @keywords internal
#' @noRd
rowQuantiles = function(x, p){
  x = matrix(x, ncol = dim(x)[length(dim(x))])
  p = rep_len(p, nrow(x))
  n = rowSums(!is.na(x))
  sorted = t(apply(x, 1, sort, na.last = TRUE))
  if(ncol(x) == 1) sorted = t(sorted)
  h = (n - 1) * p + 1
  lo = floor(h)
  hi = ceiling(h)
  ok = n > 0 & !is.na(p)
  result = rep(NA_real_, nrow(x))
  rows = which(ok)
  x_lo = sorted[cbind(rows, lo[ok])]
  x_hi = sorted[cbind(rows, hi[ok])]
  #Avoid Inf - Inf when both ends are the same
  result[ok] = ifelse(x_lo == x_hi, x_lo, x_lo + (h[ok] - lo[ok]) * (x_hi - x_lo))
  result
}

#' Bootstrap intervals for each row of a matrix of replicates.
#'
#' @param t0 Numeric vector of original values (length `k`).
#' @param t `k × R` matrix of resample values, `NA` where missing.
#' @param jk `k × nDocs` matrix of jackknife values, `NA` where missing, or
#'   `NULL`. Needed for `"bca"`, and for `"stud"` with `inner = "jackknife"`.
#' @param type Interval types, any of `"perc"`, `"bca"` and `"stud"`.
#' @param conf Confidence levels.
#' @param se_t `k × R` matrix of the inner SE of each resample value, or
#'   `NULL`. Needed for `"stud"`.
#' @param inner `"jackknife"` or `"bootstrap"`: whether `"stud"` scales by
#'   the jackknife SE of the original data or the SD of the resamples.
#' @return A tibble with `k × length(type) × length(conf)` rows, in the order
#'   type, then conf, then row. Columns: `.row` (the row of `t`), `type`,
#'   `conf`, `estimate`, `lower`, `upper`, `n_resamples`, `se`, `bias`,
#'   `se_jack`, `z0`, `acc`, `n_stud`.
#' @keywords internal
#' @noRd
getIntervalsFromMatrix = function(t0, t, jk, type, conf, se_t = NULL, inner = "jackknife"){
  n_resamples = rowSums(!is.na(t))
  se = rowSDs(t)
  se_jack = if(is.null(jk)) rep(NA_real_, length(t0)) else jackknifeSE(jk)
  base = tibble(.row = seq_along(t0), estimate = t0, n_resamples = n_resamples,
                se = se, bias = rowMeans(t, na.rm = TRUE) - t0, se_jack = se_jack)

  if("bca" %in% type){
    #Ties with the original value count as half
    below = rowSums(t < t0, na.rm = TRUE) + rowSums(t == t0, na.rm = TRUE) / 2
    z0 = qnorm(below / n_resamples)
    d = rowMeans(jk, na.rm = TRUE) - jk
    ss = rowSums(d^2, na.rm = TRUE)
    acc = ifelse(ss == 0, 0, rowSums(d^3, na.rm = TRUE) / (6 * ss^1.5))
  }
  if("stud" %in% type){
    t_stud = (t - t0) / se_t
    t_stud[!is.finite(t_stud)] = NA
    n_stud = rowSums(!is.na(t_stud))
    scale = if(inner == "jackknife") se_jack else se
  }

  lapply(type, function(curr_type){
    lapply(conf, function(curr_conf){
      alpha = (1 + c(-curr_conf, curr_conf)) / 2
      z0_out = acc_out = n_stud_out = NA_real_
      if(curr_type == "perc"){
        bounds = lapply(alpha, function(a) rowQuantiles(t, a))
      } else if(curr_type == "bca"){
        bounds = lapply(qnorm(alpha), function(z){
          p = pnorm(z0 + (z0 + z) / (1 - acc * (z0 + z)))
          rowQuantiles(t, ifelse(is.finite(p), p, NA))
        })
        z0_out = z0
        acc_out = acc
      } else {
        #The upper quantile of t* gives the lower bound
        bounds = lapply(rev(alpha), function(a) t0 - rowQuantiles(t_stud, a) * scale)
        n_stud_out = n_stud
      }
      base %>%
        mutate(type = curr_type, conf = curr_conf, lower = bounds[[1]], upper = bounds[[2]],
               z0 = z0_out, acc = acc_out, n_stud = n_stud_out) %>%
        relocate(".row", "type", "conf", "estimate", "lower", "upper")
    }) %>% bind_rows
  }) %>% bind_rows
}
