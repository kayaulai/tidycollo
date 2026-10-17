#Independent reference implementations and fixtures for checking the measures.

#x * log(x / y) with 0 * log(0) = 0
refXLogXOverY = function(x, y, base = exp(1)) ifelse(x == 0, 0, x * log(x / y, base))

#Four cells of the 2x2 table for one combination, plus expected frequencies
refTable = function(a, f_1, f_2, N){
  obs = c(a, f_1 - a, f_2 - a, N - f_1 - f_2 + a)
  exp = c(f_1 * f_2, f_1 * (N - f_2), (N - f_1) * f_2, (N - f_1) * (N - f_2)) / N
  list(obs = obs, exp = exp)
}

#All ways of putting f tokens into docs without exceeding each doc's size
refAllocations = function(f, sizes){
  if(length(sizes) == 1) return(if(f <= sizes) matrix(f) else NULL)
  do.call(rbind, lapply(0:min(f, sizes[1]), function(k){
    rest = refAllocations(f - k, sizes[-1])
    if(is.null(rest)) NULL else cbind(k, rest)
  }))
}

refDP = function(freqs, sizes) sum(abs(freqs / sum(freqs) - sizes / sum(sizes))) / 2

#Ten tokens over four docs of sizes 1-4
smallCorpus = function(){
  dplyr::tibble(doc = c(rep("A", 1), rep("B", 2), rep("C", 3), rep("D", 4)),
                x = c("p", "p", "q", "p", "q", "r", "p", "q", "r", "q"),
                y = c("z", "y", "z", "y", "z", "z", "y", "z", "z", "y"))
}

#Two docs of five tokens each; doc A repeats one combination four times
repeatCorpus = function(){
  dplyr::tibble(doc = c(rep("A", 5), rep("B", 5)),
                x = c("p", "p", "p", "p", "q", "p", "q", "r", "s", "t"),
                y = c("z", "z", "z", "z", "y", "z", "y", "y", "y", "y"))
}

runMeasures = function(data, measures = list("bam_str", "bam_test", "uam_str", "uam_rank", "prod", "disp")){
  suppressMessages(capture.output(
    res <- getMeasures(data, c("x", "y"), c("x", "y"), "doc", measures)))
  res$results
}
