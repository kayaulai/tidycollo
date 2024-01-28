# This file is part of the standard setup for testthat.
# It is recommended that you do not modify it.
#
# Where should you do additional test configuration?
# Learn more about the roles of various files in:
# * https://r-pkgs.org/testing-design.html#sec-tests-files-overview
# * https://testthat.r-lib.org/articles/special-files.html

library(testthat)
library(tidycollo)

test_check("tidycollo")

test_that("Extraction from CONLLU works", {
  #In the final thing I probably shouldn't import something here
  library(udpipe)
  conllu_all = udpipe::udpipe_read_conllu("inst/extdata/lzh_kyoto-ud-test.conllu")
  conllu_sents = conllu_all %>% group_by(doc_id, paragraph_id, sentence_id, sentence) %>% group_split
  conllu_sent = conllu_sents[[106]]
  tree = getTreeFromConllu(conllu_sent)
  conllu_sent_wtree = getTreeFromConllu(conllu_sent)


  getCombos = function(root, conllu_sent){
    print(glue("{conllu_sent$doc_id[1]} -  {conllu_sent$paragraph_id[1]} - {conllu_sent$sentence[1]}"))
    root_kids = as.numeric(names(root$children))
    rootID = as.numeric(root$name)
    sfp = conllu_sent %>%
      filter(.data$head_token_id == rootID,
             .data$xpos == "p,助詞,句末,*",
             .data$token_id > rootID
             ) %>% pull(lemma)
    adv_conj = conllu_sent %>%
      filter(
        .data$token_id < rootID,
        .data$upos %in% c("ADV", "SCONJ", "CCONJ")
        ) %>% pull(lemma)

    list(sfp = sfp, dc = adv_conj)
  }

  comboList = lapply(conllu_sents, function(conllu_sent){
    conllu_sent_wtree = getTreeFromConllu(conllu_sent)
    combos = getCombos(conllu_sent_wtree$root, conllu_sent_wtree$df)
    c(list(docID = conllu_sent$doc_id[1],
           paragraphID = conllu_sent$paragraph_id[1],
           sentenceID = conllu_sent$sentence_id[1],
           sentence = conllu_sent$sentence[1]), combos)
  })

  getCombos(conllu_sent_wtree$root, conllu_sent_wtree$df)

})
