# TODO

## Naming and downstream compatibility

- [ ] Decide how to handle names that break the `slot_1_cond_slot_2` convention without breaking
      downstream repos. One option: keep the current behaviour on a `classic` branch, and fix the names on `main`.
  - [ ] `h_{slot_1}_cond_{slot_2}` / `h_{slot_2}_cond_{slot_1}` (and the `_norm` versions): they are
        computed per slot-1 lemma (or per slot-2 lemma), so by the convention they are the other way round.
- [ ] Check downstream repos against output changes that are already on `main`, where the names stayed the same but the values changed:
  - [ ] `surp_{slot_1}` / `surp_{slot_2}` have swapped values (now `surp_x = -log2 p(x|y)`).
  - [ ] `log_tf_{slot_1}` / `log_tf_{slot_2}` have swapped values (now `log_tf_x = log2(tf_x)`).
  - [ ] `p_hapax_cond_*` now divides by the lemma's frequency instead of by corpus N.
  - [ ] `gsq`, `kld_*_cond_*`, all dispersion measures and `getBenjaminiYekutieli()$reject` were wrong
        before, so any results that used them need to be recomputed.

## Measures to review

- [ ] Normalised entropy (`getNormEntropyFromCounts`) divides by log2(tokens). The usual choice is
      log2(types). Decide which one is intended, and handle f = 1 (currently `NaN`).
- [ ] `MFDR` in `getBenjaminiYekutieli`: find a source for the formula, or remove it.
- [ ] `getUAMRank` ranks unsigned statistics (chisq, gsq), so strongly repelled pairs rank as high as
      attracted ones. Consider signed versions.
- [ ] `fisher.test(..., simulate.p.value = TRUE)` has no effect on 2x2 tables, so drop the argument.

## Bootstrap (complete revamp)

- [x] Rewrite `getMeasuresBoot` / `stratBoot`:
  - Dropping a combination that is missing from a resample (instead of counting it as 0) is intended,
    because real samples also miss rare combinations. That is documented now, and `presence` reports how often each combination is observed.
  - Resamples are joined to the original combinations once each, instead of filtering per row.
  - `existPerc` is returned as `presence`.
  - Exported and documented, along with `getBootIntervals`.
  - Added BCa (doc-level jackknife), corpus-level `statistics`, and `parallel` (future.apply).
- [x] Sped up `getMeasures`: dispersion is computed for all combinations at once (it was one tibble per
      combination), Fisher p-values come straight from the hypergeometric (no odds-ratio CI), and the
      bootstrap resamples count rows (`getMeasuresFromCounts`) instead of tokens. About 40× faster on an 80-doc corpus,
      with output identical to 1e-12. What's left is mostly dplyr overhead in the measure functions.
- [x] Studentized bootstrap (`type = "stud"`), with an inner jackknife (default) or inner bootstrap.
- [ ] Run the studentized bootstrap on HPC with more resamples (`parallel = TRUE` plus a future plan).
- [ ] Move `R/C04a_bootstrap.qmd` out of `R/`. R CMD check complains about non-.R files there.

## Input from tidytext / UD

Main features: adjacent bigrams, and construction-based combinations (POS patterns or dependencies).
Plan the tidytext and UD input paths together. Both should output the same tidy input for `getMeasures`.

- [ ] Adjacent bigrams: wrapper around `unnest_tokens(token = "ngrams", n = 2)` + `tidyr::separate`
      (works today by hand). Add `tidytext`/`tidyr` to Suggests.
- [ ] Construction-based input from POS-tagged / parsed tokens (e.g. `udpipe` output): extract slot pairs
      by POS pattern or dependency relation. This is where the ud.R stubs come in (`getCombosFromUD`, `getParent`,
      `getChildren`, etc.).
- [ ] Lemmatisation isn't in tidytext. `udpipe` gives lemmas and POS together, so it could feed both paths.

Not planned for now: window-based and sentence/document co-occurrence. If they are added, note
that `getComponentf`/`getBasicPs` take `f_x`, `f_y` and `N` from the pair table. That fits bigrams and
constructions, but surface/textual co-occurrence (Evert) should take them from the corpus instead.

## Package

- [ ] Run `devtools::check()` once there is free disk space (the C: drive was full).
- [ ] Fill in the `measures` documentation in `getMeasures` (the list of abbreviations).
- [ ] Decide whether `getMeasuresBoot`, `getBenjaminiYekutieli` and the ud.R functions should be exported.
- [ ] Finish the ud.R stubs (`getCombosFromUD`, `getParent`, `getChildren`, etc.).
