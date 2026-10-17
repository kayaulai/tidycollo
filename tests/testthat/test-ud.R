test_that("getTreeFromConllu builds one node per token, in token order, with the right lemmas", {
  skip_if_not_installed("udpipe")
  conllu_all = udpipe::udpipe_read_conllu(
    system.file("extdata", "lzh_kyoto-ud-test.conllu", package = "tidycollo"))
  conllu_sents = conllu_all %>%
    dplyr::group_by(doc_id, paragraph_id, sentence_id, sentence) %>%
    dplyr::group_split()
  conllu_sent = conllu_sents[[106]]

  tree = getTreeFromConllu(conllu_sent)
  expect_length(tree$nodesList, nrow(conllu_sent))
  expect_equal(sapply(tree$nodesList, function(x) as.numeric(x$name)),
               as.numeric(conllu_sent$token_id))
  expect_equal(sapply(tree$nodesList, function(x) x$lemma), conllu_sent$lemma)
  expect_equal(as.numeric(tree$root$name),
               as.numeric(conllu_sent$token_id[conllu_sent$head_token_id == "0"]))
})
