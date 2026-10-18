#Core functions for getting co-occurrence measures.


#' Get common co-occurrence measures.
#'
#' @param data Tidy-formatted data, with one token of a certain combination in each line. It should include at least columns for lemmas of the two slots, plus one for the document in which the combination appears.
#' @param lemma_names Column names for the lemmas of the two slots.
#' @param slot_names Abbreviations for the two slots.
#' @param doc_id Column name containing documents.
#' @param measures A list of abbreviations. See Notes below. If you just want the basic measures of each type (frequency, attraction, productivity, dispersion) and don't plan to use fancier meaures (e.g. based on LNRE models) or composite (e.g. lexical gravity), the default should work fine. #TODO: actually write the list,.
#'
#' @return A list. The first value, results, is usually all you want. You are also provided with several other views of the original data, including a list of combinations with frequencies by document, a term-document matrix, and a contingency table.
#' @export
#'
#' @examples
#' data = dplyr::tibble(doc = c("A", "B", "B", "C", "C", "C"),
#'                      x_lemma = c("p", "p", "q", "p", "q", "r"),
#'                      y_lemma = c("z", "y", "z", "y", "z", "z"))
#' res = getMeasures(data, lemma_names = c("x_lemma", "y_lemma"),
#'                   slot_names = c("x", "y"), doc_id = "doc")
#' res$results
getMeasures = function(data,
                       lemma_names = c("lemma_1", "lemma_2"),
                       slot_names = c("f", "s"), doc_id = "doc",
                       measures = list("bam_str", "bam_test", "uam_str", "uam_rank", "prod", "disp")){

  message("Some basic operations ...")
  message("Data views")
  #Data views
  data_by_combo_doc = data %>%
    group_by(!!parse_expr(lemma_names[1]), !!parse_expr(lemma_names[2]), !!parse_expr(doc_id)) %>%
    count %>%
    ungroup
  getMeasuresFromCounts(data_by_combo_doc, lemma_names, slot_names, doc_id, measures)
}

#getMeasures from counts instead of tokens. data_by_combo_doc has one row per
#combination per doc: the two lemma columns, the doc_id column and the count n.
#Returns the same list as getMeasures.
getMeasuresFromCounts = function(data_by_combo_doc,
                                 lemma_names = c("lemma_1", "lemma_2"),
                                 slot_names = c("f", "s"), doc_id = "doc",
                                 measures = list("bam_str", "bam_test", "uam_str", "uam_rank", "prod", "disp")){
  #Basic info from input
  lemma_1 = lemma_names[1]
  lemma_2 = lemma_names[2]
  slot_1 = slot_names[1]
  slot_2 = slot_names[2]

  data_by_combo = data_by_combo_doc %>%
    group_by(!!parse_expr(lemma_1), !!parse_expr(lemma_2)) %>%
    summarise(n = sum(n), .groups = "drop")

  data_by_combo = data_by_combo %>% mutate(log_f = log(n))

  #Basic measures
  data_by_combo = getComponentf(data_by_combo, parse_expr(lemma_1), parse_expr(lemma_2), slot_1, slot_2)
  data_by_combo = getBasicPs(data_by_combo, slot_1, slot_2)
  data_by_combo = getTableComps(data_by_combo, slot_1, slot_2)

  message("Computing attraction measures ...")
  #Attraction
  if("bam_str" %in% measures) data_by_combo = getBAMStr(data_by_combo, slot_1, slot_2)
  if("bam_test" %in% measures) data_by_combo = getBAMTest(data_by_combo, slot_1, slot_2)
  else if("bam_test_no_fisher" %in% measures) data_by_combo = getBAMTest(data_by_combo, slot_1, slot_2, skip.fisher = T)
  if("uam_str" %in% measures) data_by_combo = getUAMStr(data_by_combo, slot_1, slot_2)
  if("uam_rank" %in% measures) data_by_combo = getUAMRank(data_by_combo, parse_expr(lemma_1), parse_expr(lemma_2), slot_1, slot_2)

  message("Computing productivity measures ...")
  #Productivity
  if("prod" %in% measures) data_by_combo = getProdSimple(data_by_combo, parse_expr(lemma_1), parse_expr(lemma_2), slot_1, slot_2)

  message("Computing dispersion measures ...")
  #Dispersion
  if("disp" %in% measures) data_by_combo = getDispSimple(data_by_combo, data_by_combo_doc, parse_expr(doc_id), parse_expr(lemma_1), parse_expr(lemma_2))

  message("Done!")
  list(results = data_by_combo,
       counts_by_doc_combo = data_by_combo_doc,
       tdm = NA,
       ct = NA)
  #TODO: Contingency table, TDM

}

getComponentf = function(df, lemma_1, lemma_2, slot_1, slot_2){
  df %>%
    group_by(!!lemma_1) %>%
    mutate("f_{slot_1}" := sum(n)) %>%
    ungroup %>%
    group_by(!!lemma_2) %>%
    mutate("f_{slot_2}" := sum(n)) %>%
    ungroup
}

getBasicPs = function(df, slot_1, slot_2){
  N = sum(df$n)
  f_1 = parse_expr(glue("f_{slot_1}"))
  f_2 = parse_expr(glue("f_{slot_2}"))
  p_1_cond_2 = parse_expr(glue("p_{slot_1}_cond_{slot_2}"))
  p_2_cond_1 = parse_expr(glue("p_{slot_2}_cond_{slot_1}"))
  p_1_cond_not_2 = parse_expr(glue("p_{slot_1}_cond_not_{slot_2}"))
  p_2_cond_not_1 = parse_expr(glue("p_{slot_2}_cond_not_{slot_1}"))

  df %>%
    mutate("p_{slot_1}_cond_{slot_2}" := (.data$n / !!f_2), #p(s1|s2)
           "p_{slot_2}_cond_{slot_1}" := (.data$n / !!f_1)) %>% #p(s2|s1)
    mutate("p_{slot_1}" := !!f_1 / N, #p(s1)
           "p_{slot_2}" := !!f_2 / N) %>% #p(s2)
    mutate(p = n / N) %>% #p(s1, s2)
    mutate("p_{slot_1}_cond_not_{slot_2}" := ((!!f_1 - .data$n) / (N - !!f_2)), #p(s1|~s2)
           "p_{slot_2}_cond_not_{slot_1}" := ((!!f_2 - .data$n) / (N - !!f_1))) %>% #p(s2|~s1)
    mutate("p_not_{slot_1}_cond_{slot_2}" := 1 - !!p_1_cond_2, #p(~s1|s2)
           "p_not_{slot_2}_cond_{slot_1}" := 1 - !!p_2_cond_1) %>% #p(~s2|s1)
    mutate("p_not_{slot_1}_cond_not_{slot_2}" := 1 - !!p_1_cond_not_2, #p(~s1|~s2)
           "p_not_{slot_2}_cond_not_{slot_1}" := 1 - !!p_2_cond_not_1) #p(~s2|~s1)
}


getTableComps = function(df, slot_1, slot_2){
  N = sum(df$n)
  f_1 = parse_expr(glue("f_{slot_1}"))
  f_2 = parse_expr(glue("f_{slot_2}"))
  p_1 = parse_expr(glue("p_{slot_1}"))
  p_2 = parse_expr(glue("p_{slot_2}"))

  df %>%
    mutate("f_not_{slot_1}" := N - !!f_1,
           "f_not_{slot_2}" := N - !!f_2) %>%
    mutate("f_{slot_1}_not_{slot_2}" := !!f_1 - .data$n,
           "f_{slot_2}_not_{slot_1}" := !!f_2 - .data$n) %>%
    mutate("f_not_{slot_1}_or_{slot_2}" := N - !!f_1 - !!f_2 + .data$n) %>%
    mutate(ef := N * !!p_1 * !!p_2,
           "ef_{slot_1}_not_{slot_2}" := N * !!p_1 * (1 - !!p_2),
           "ef_{slot_2}_not_{slot_1}" := N * !!p_2 * (1 - !!p_1)) %>%
    mutate("ef_not_{slot_1}_or_{slot_2}" := N * (1 - !!p_1) * (1 - !!p_2))
}

getBAMStr = function(df, slot_1, slot_2){
  N = sum(df$n)
  p_1 = parse_expr(glue("p_{slot_1}"))
  p_2 = parse_expr(glue("p_{slot_2}"))
  f_1_not_2 = parse_expr(glue("f_{slot_1}_not_{slot_2}"))
  f_2_not_1 = parse_expr(glue("f_{slot_2}_not_{slot_1}"))
  f_not_1_or_2 = parse_expr(glue("f_not_{slot_1}_or_{slot_2}"))
  f_1 = parse_expr(glue("f_{slot_1}"))
  f_2 = parse_expr(glue("f_{slot_2}"))

  df %>%
    mutate(pmi = log2(.data$p / (!!p_1 * !!p_2)),
           or = (as.numeric(n) / !!f_1_not_2) * (as.numeric(!!f_not_1_or_2) /
                   !!f_2_not_1),
           lor = log((n + .5) * (!!f_not_1_or_2 + .5) /
                          ((!!f_1_not_2 + .5) * (!!f_2_not_1 + .5))),
           dice = 2 * .data$n / (!!f_1 + !!f_2))
}


getBAMTest = function(df, slot_1, slot_2, skip.fisher = F){
  N = sum(df$n)
  p_1 = parse_expr(glue("p_{slot_1}"))
  p_2 = parse_expr(glue("p_{slot_2}"))
  f_1_not_2 = parse_expr(glue("f_{slot_1}_not_{slot_2}"))
  f_2_not_1 = parse_expr(glue("f_{slot_2}_not_{slot_1}"))
  ef_1_not_2 = parse_expr(glue("ef_{slot_1}_not_{slot_2}"))
  ef_2_not_1 = parse_expr(glue("ef_{slot_2}_not_{slot_1}"))
  f_not_1_or_2 = parse_expr(glue("f_not_{slot_1}_or_{slot_2}"))
  ef_not_1_or_2 = parse_expr(glue("ef_not_{slot_1}_or_{slot_2}"))
  f_1 = parse_expr(glue("f_{slot_1}"))
  f_2 = parse_expr(glue("f_{slot_2}"))
  p_1_cond_2 = parse_expr(glue("p_{slot_1}_cond_{slot_2}"))
  p_2_cond_1 = parse_expr(glue("p_{slot_2}_cond_{slot_1}"))
  p_1_cond_not_2 = parse_expr(glue("p_{slot_1}_cond_not_{slot_2}"))
  p_2_cond_not_1 = parse_expr(glue("p_{slot_2}_cond_not_{slot_1}"))

  df = df %>%
    mutate(res_pearson = (.data$n - .data$ef) / sqrt(.data$ef)) %>%
    mutate(chisq = .data$res_pearson^2 +
             (!!f_1_not_2 - !!ef_1_not_2)^2 / !!ef_1_not_2 +
             (!!f_2_not_1 - !!ef_2_not_1)^2 / !!ef_2_not_1 +
             (!!f_not_1_or_2 - !!ef_not_1_or_2)^2 / !!ef_not_1_or_2)  %>%
    mutate(t = (.data$n - .data$ef) / sqrt(.data$n))%>%
    #G^2 = 2 * sum(O * log(O / E)) over the four cells, with 0 * log(0) = 0
    mutate(gsq = 2 * (xLogXOverY(.data$n, .data$ef) +
                        xLogXOverY(!!f_1_not_2, !!ef_1_not_2) +
                        xLogXOverY(!!f_2_not_1, !!ef_2_not_1) +
                        xLogXOverY(!!f_not_1_or_2, !!ef_not_1_or_2)))

  if(!skip.fisher){
    df = df %>%
      mutate(p_fisher_yates = getFisherP(.data$n, !!f_1, !!f_2, N),
             p_fisher_yates_attract = getFisherP(.data$n, !!f_1, !!f_2, N, "greater"))
  }

  df
}

#p-values of Fisher's exact test on the 2x2 table of each combination, as in
#fisher.test(), without the odds ratio estimate and its CI. a, f_1 and f_2 are
#vectors of combination and slot frequencies; N is the number of tokens.
getFisherP = function(a, f_1, f_2, N, alternative = "two.sided"){
  if(alternative == "greater") return(phyper(a - 1, f_1, N - f_1, f_2, lower.tail = FALSE))
  vapply(seq_along(a), function(i){
    lo = max(0, f_2[i] - (N - f_1[i]))
    hi = min(f_2[i], f_1[i])
    d = dhyper(lo:hi, f_1[i], N - f_1[i], f_2[i], log = TRUE)
    d = exp(d - max(d))
    d = d / sum(d)
    #Tables at most as likely as the observed one, with fisher.test's tolerance
    sum(d[d <= d[a[i] - lo + 1] * (1 + 10^(-7))])
  }, numeric(1))
}



#x * log(x / y), taking 0 * log(0) = 0
xLogXOverY = function(x, y, base = exp(1)){
  ifelse(x == 0, 0, x * log(x / y, base))
}

getUAMStr = function(df, slot_1, slot_2){
  N = sum(df$n)
  f_1 = parse_expr(glue("f_{slot_1}"))
  f_2 = parse_expr(glue("f_{slot_2}"))
  p_1 = parse_expr(glue("p_{slot_1}"))
  p_2 = parse_expr(glue("p_{slot_2}"))
  p_1_cond_2 = parse_expr(glue("p_{slot_1}_cond_{slot_2}"))
  p_2_cond_1 = parse_expr(glue("p_{slot_2}_cond_{slot_1}"))
  p_1_cond_not_2 = parse_expr(glue("p_{slot_1}_cond_not_{slot_2}"))
  p_2_cond_not_1 = parse_expr(glue("p_{slot_2}_cond_not_{slot_1}"))
  kld_1_cond_2 = parse_expr(glue("kld_{slot_1}_cond_{slot_2}"))
  kld_2_cond_1 = parse_expr(glue("kld_{slot_2}_cond_{slot_1}"))


  df %>%
    mutate("surp_{slot_1}" := -log2(n / !!f_2), #-log2 p(s1|s2)
           "surp_{slot_2}" := -log2(n / !!f_1)) %>% #-log2 p(s2|s1)
    mutate("Dp_{slot_1}_on_{slot_2}" := !!p_1_cond_2 - !!p_1_cond_not_2,
           "Dp_{slot_2}_on_{slot_1}" := !!p_2_cond_1 - !!p_2_cond_not_1) %>%
    mutate("kld_{slot_1}_cond_{slot_2}" := xLogXOverY(!!p_1_cond_2, !!p_1, 2) +
             xLogXOverY(1 - !!p_1_cond_2, 1 - !!p_1, 2),
           "kld_{slot_2}_cond_{slot_1}" := xLogXOverY(!!p_2_cond_1, !!p_2, 2) +
             xLogXOverY(1 - !!p_2_cond_1, 1 - !!p_2, 2)) %>%
    mutate("kld_norm_{slot_1}_cond_{slot_2}" := 1 - exp(-!!kld_1_cond_2),
           "kld_norm_{slot_2}_cond_{slot_1}" := 1 - exp(-!!kld_2_cond_1))
}



getUAMRank = function(df, lemma_1, lemma_2, slot_1, slot_2, stats = c("chisq", "t", "gsq")){
  N = sum(df$n)

  for(stat in stats){
    df = df  %>%
      group_by(!!lemma_2) %>%
      mutate("rank_{stat}_{slot_1}_on_{slot_2}" := rank(!!parse_expr(stat)),
             "rank_inv_{stat}_{slot_1}_on_{slot_2}" := rank(-!!parse_expr(stat))) %>%
      ungroup %>%
      group_by(!!lemma_1) %>%
      mutate("rank_{stat}_{slot_2}_on_{slot_1}" := rank(!!parse_expr(stat)),
             "rank_inv_{stat}_{slot_2}_on_{slot_1}" := rank(-!!parse_expr(stat))) %>%
      ungroup
  }

  df
}



getEntropyFromCounts = function(counts){
  p = counts / sum(counts)
  -sum(p * log2(p))
}

getNormEntropyFromCounts = function(counts){
  p = counts / sum(counts)
  -sum(p * log2(p)) / log2(sum(counts))
}

getProdSimple = function(df, lemma_1, lemma_2, slot_1, slot_2){
  N = sum(df$n)

  df %>%
    group_by(!!lemma_1) %>%
    mutate("h_{slot_1}_cond_{slot_2}" := getEntropyFromCounts(n),
           "h_{slot_1}_cond_{slot_2}_norm" := getNormEntropyFromCounts(n),
           "tf_{slot_1}" := length(unique(!!lemma_2)),
           "log_tf_{slot_1}" := log2(length(unique(!!lemma_2))),
           "p_hapax_cond_{slot_1}" := sum(n == 1)/ sum(n)) %>% #Baayen's P
    ungroup() %>%
    group_by(!!lemma_2) %>%
    mutate("h_{slot_2}_cond_{slot_1}" := getEntropyFromCounts(n),
           "h_{slot_2}_cond_{slot_1}_norm" := getNormEntropyFromCounts(n),
           "tf_{slot_2}" := length(unique(!!lemma_1)),
           "log_tf_{slot_2}" := log2(length(unique(!!lemma_1))),
           "p_hapax_cond_{slot_2}" := sum(n == 1)/ sum(n)) %>% #Baayen's P
    ungroup()

}

getDPFloor_s2 = function(f, perc_docsize){
  draws = sapply(1:1000, function(x) rmultinom(f, 1, perc_docsize) %>% rowSums / f)
  apply(draws, 2, function(x) sum(abs(x - perc_docsize)) / 2) %>% min
}

#Strategy 3
#Largest-remainder (Hamilton) apportionment of f tokens to docs in proportion to
#doc size. This minimises sum(abs(freqs - f * perc_docsize)), and hence DP, over
#all integer allocations.
getDPFloor = function(f, perc_docsize){
  quotas = f * perc_docsize
  freqs_perdoc = floor(quotas)
  remainder = round(f - sum(freqs_perdoc))
  if(remainder > 0){
    top = order(quotas - freqs_perdoc, decreasing = TRUE)[1:remainder]
    freqs_perdoc[top] = freqs_perdoc[top] + 1
  }
  percs_perdoc = freqs_perdoc / f
  sum(abs(percs_perdoc - perc_docsize)) / 2
}

getDPCeiling = function(f, perc_docsize, abs_docsize){
  #This is slightly modified from Gries' method to account for very small docs.
  #We 'fill' the smallest doc first. Once that's filled, we go to the next smallest
  #doc, and so on.
  smallestFirst = order(perc_docsize)
  sizes = abs_docsize[smallestFirst]
  filledBefore = c(0, cumsum(sizes)[-length(sizes)])
  abs_bydoc_currcombo = numeric(length(sizes))
  abs_bydoc_currcombo[smallestFirst] = pmax(0, pmin(sizes, f - filledBefore))
  perc_bydoc_currcombo = abs_bydoc_currcombo / sum(abs_bydoc_currcombo)
  sum(abs(perc_docsize - perc_bydoc_currcombo)) / 2
}

#Taken directly from Gries (2022)
zero.to.one = function (x) { (y = x - min(x))/max(y) }

#Unused: the original one-combination-at-a-time version of the dispersion measures,
#kept as an easier-to-read reference for the vectorised getDispSimple below.
#sub_df has the rows of df_bydoc for one combination; floor and ceil are its DP
#floor and ceiling. Returns a one-row tibble with the lemmas and the measures.
getDispFromDFPart = function(sub_df, df_bydoc_totals, docs_list, lemma_1, lemma_2, doc_id,
                             floor, ceil){
  curr_lemma_1 = pull(sub_df, !!lemma_1)[1]
  curr_lemma_2 = pull(sub_df, !!lemma_2)[1]
  # print(glue("{curr_lemma_1} {curr_lemma_2}"))

  nd = length(docs_list) #no. of docs
  freqs_bydoc_no0 = sub_df$n
  absent_docs = nd - length(freqs_bydoc_no0)
  freqs_bydoc = c(freqs_bydoc_no0, rep(0, absent_docs))
  names(freqs_bydoc) = c(pull(sub_df, !!doc_id), setdiff(docs_list, pull(sub_df, !!doc_id)))
  f = sum(sub_df$n)

  range = length(freqs_bydoc_no0) / length(docs_list)
  sd_count = sd(freqs_bydoc) * sqrt((nd - 1) / nd)
  var_coef_count = sd_count / mean(freqs_bydoc)
  idf = log2(1 / range)

  #Percentage of the combo within each doc
  total_bydoc = df_bydoc_totals$n
  names(total_bydoc) = pull(df_bydoc_totals, !!doc_id)
  total_bydoc = total_bydoc[names(freqs_bydoc)]
  perc_bydoc_allcombos = freqs_bydoc / total_bydoc
  perc_bydoc_currcombo = freqs_bydoc / sum(freqs_bydoc)

  #Percentage of combo within each doc
  perc_docsize = total_bydoc / sum(total_bydoc)
  efreq_bydoc = perc_docsize * f

  #Names are skipped for any var names that are more than just a letter
  sd_perc = sd(perc_bydoc_allcombos) * sqrt((nd - 1) / nd)
  var_coef_perc = sd_perc / mean(perc_bydoc_allcombos)
  juilland_d = 1 - var_coef_count / sqrt(nd - 1)
  juilland_d_perc = 1 - var_coef_perc / sqrt(nd - 1)
  rosengren_s = (sum(sqrt(perc_docsize * freqs_bydoc))) ^ 2 / f #Unequal-parts version
  D2 = (log2(f) - sum(freqs_bydoc * log2(freqs_bydoc), na.rm = T) / f) / log2(nd)
  chisq = sum((freqs_bydoc - efreq_bydoc)^2 / efreq_bydoc)
  D3 = 1 - chisq / (4 * f)
  DC = ((sum(sqrt(freqs_bydoc)))/nd)^2 / (f / nd)
  DP = sum(abs(perc_docsize - perc_bydoc_currcombo)) / 2
  DP_norm = DP / (1 - min(perc_docsize))
  DP_nofreq = zero.to.one(c(floor, ceil, DP))[3]
  kld = sum(xLogXOverY(perc_bydoc_currcombo, perc_docsize, 2))
  kld_norm = 1 - exp(-kld)

  #Adj. freqs
  juilland_u = juilland_d * f
  rosengren_af = rosengren_s * f
  Um = (sum(freqs_bydoc)) * D2 + (1 - D2) * f / nd
  engvall_af = f * range

  tibble(!!deparse(lemma_1) := curr_lemma_1, !!deparse(lemma_2) := curr_lemma_2,
         range = range, sd = sd_count, var_coef = var_coef_count,
         sd_perc = sd_perc, var_coef_perc = var_coef_perc,
         idf = idf,
         juilland_d = juilland_d, juilland_d_perc = juilland_d_perc, juilland_u = juilland_u,
         D2 = D2, Um = Um, D3 = D3,
         rosengren_s = rosengren_s, rosengren_af = rosengren_af,
         chisq_disp = chisq, DC = DC,
         kld_disp = kld, kld_disp_norm = kld_norm,
         DP = DP, DP_norm = DP_norm, DP_nofreq = DP_nofreq,
         engvall_af = engvall_af)
}


#Dispersion measures for every combination at once. Each sum over docs is a sum
#over the docs where the combination occurs (the rows of df_bydoc), plus a closed
#form for the docs where it doesn't, so no combination-by-doc matrix is built.
getDispSimple = function(df, df_bydoc, doc_id, lemma_1, lemma_2){
  message(">>Starting dispersion calculation...")
  docs = as.character(pull(df_bydoc, !!doc_id))
  docs_list = unique(docs)
  nd = length(docs_list) #no. of docs
  doc_row = match(docs, docs_list)
  combo_row = matchCombos(df_bydoc, df, c(deparse(lemma_1), deparse(lemma_2)))
  #Sum of x over the docs where each combination occurs
  sumByCombo = function(x) as.vector(rowsum(rep_len(as.numeric(x), length(combo_row)), combo_row))

  #Token counts per doc (df_bydoc has one row per combo per doc)
  total_bydoc = as.vector(rowsum(as.numeric(df_bydoc$n), doc_row))
  perc_docsize = total_bydoc / sum(total_bydoc)

  v = as.numeric(df_bydoc$n) #Freq of the combo in each doc where it occurs
  s = perc_docsize[doc_row] #Size of those docs, as proportions
  f = sumByCombo(v)
  f_v = f[combo_row]
  k = sumByCombo(1) #No. of docs where the combo occurs
  absent = nd - k
  s_present = sumByCombo(s) #Proportion of the corpus in those docs

  #Population SDs, adding (0 - mean)^2 for each absent doc
  mean_count = f / nd
  sd_count = sqrt((sumByCombo((v - mean_count[combo_row])^2) + absent * mean_count^2) / nd)
  perc_bydoc_allcombos = v / total_bydoc[doc_row]
  mean_perc = sumByCombo(perc_bydoc_allcombos) / nd
  sd_perc = sqrt((sumByCombo((perc_bydoc_allcombos - mean_perc[combo_row])^2) + absent * mean_perc^2) / nd)
  perc_bydoc_currcombo = v / f_v

  range = k / nd
  var_coef_count = sd_count / mean_count
  idf = log2(1 / range)
  var_coef_perc = sd_perc / mean_perc
  juilland_d = 1 - var_coef_count / sqrt(nd - 1)
  juilland_d_perc = 1 - var_coef_perc / sqrt(nd - 1)
  rosengren_s = sumByCombo(sqrt(s * v))^2 / f #Unequal-parts version
  D2 = (log2(f) - sumByCombo(v * log2(v)) / f) / log2(nd)
  #An absent doc adds (0 - f * s)^2 / (f * s) = f * s to chisq
  chisq = sumByCombo((v - f_v * s)^2 / (f_v * s)) + f * (1 - s_present)
  D3 = 1 - chisq / (4 * f)
  DC = (sumByCombo(sqrt(v)) / nd)^2 / (f / nd)
  #An absent doc adds its size to the sum
  DP = (sumByCombo(abs(s - perc_bydoc_currcombo)) + (1 - s_present)) / 2
  DP_norm = DP / (1 - min(perc_docsize))
  kld = sumByCombo(perc_bydoc_currcombo * log2(perc_bydoc_currcombo / s))
  kld_norm = 1 - exp(-kld)

  #DP rescaled between its floor and ceiling for the combo's frequency
  freq_types = unique(f)
  dpCeils = sapply(freq_types, function(x) getDPCeiling(x, perc_docsize, total_bydoc))[match(f, freq_types)]
  dpFloors = sapply(freq_types, function(x) getDPFloor(x, perc_docsize))[match(f, freq_types)]
  dp_min = pmin(dpFloors, dpCeils, DP)
  DP_nofreq = (DP - dp_min) / (pmax(dpFloors, dpCeils, DP) - dp_min)

  #Not mutate(), where columns of df (e.g. chisq) would mask these variables
  bind_cols(df, tibble(range = range, sd = sd_count, var_coef = var_coef_count,
                       sd_perc = sd_perc, var_coef_perc = var_coef_perc,
                       idf = idf,
                       juilland_d = juilland_d, juilland_d_perc = juilland_d_perc,
                       juilland_u = juilland_d * f, #Adj. freqs
                       D2 = D2, Um = f * D2 + (1 - D2) * f / nd, D3 = D3,
                       rosengren_s = rosengren_s, rosengren_af = rosengren_s * f,
                       chisq_disp = chisq, DC = DC,
                       kld_disp = kld, kld_disp_norm = kld_norm,
                       DP = DP, DP_norm = DP_norm, DP_nofreq = DP_nofreq,
                       engvall_af = f * range))
}


getCompMeasures = function(){
  #TODO:
  # LFMD (log-frequency biased mutual dependency): MD + log2(P(w1w2))
  # Lexical gravity
}

getBenjaminiYekutieli = function(pval, alpha = .05){
  pval_order = order(pval)
  pval_arr = pval[pval_order]

  getHarmonic = function(m){
    result = 0
    for(i in 1:m) result = result + 1/i
    result
  }
  m = length(pval_arr)
  cm = getHarmonic(m)
  cv = 1:length(pval_arr) / (length(pval_arr) * cm) * alpha
  #Step-up: reject the k smallest p-values, where k is the largest i with p_(i) <= cv_i
  k = max(c(0, which(pval_arr <= cv)))
  reject = seq_len(m) <= k
  MFDR = alpha * (m + 1) /  (2 * m * (log(m) + -digamma(1))+ 1)

  cv_orig = numeric(m)
  reject_orig = logical(m)
  cv_orig[pval_order] = cv
  reject_orig[pval_order] = reject
  result = list(cv = cv_orig, reject = reject_orig, MFDR = MFDR)
  result
}


#TODO:
# From udpipe:
# MD (mutual dependency): log2(P(w1w2)^2 / P(w1) P(w2))
