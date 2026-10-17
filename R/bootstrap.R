# data = allPairs
# lemma_names = c("dc", "sfp")
# slot_names = c("dc", "sfp")
# doc_id = "doc"
# measures = list("bam_str", "bam_test", "uam_str", "prod")
# a = getMeasuresBoot(data, lemma_names, slot_names, doc_id, measures, R = 5)
getMeasuresBoot = function(data,
                           lemma_names = c("lemma_1", "lemma_2"),
                           slot_names = c("f", "s"),
                           doc_id = "doc",
                           measures = list("bam_str", "bam_test", "uam_str", "uam_rank", "prod", "disp"),
                           R = 100,
                           resampleFunc = stratBoot,
                           quantiles = c(0, .025, .5, .975, 1),
                           retainResamples = F,
                           sampleRatio = 1){
  data = data %>%
    mutate("{doc_id}" := as.factor(!!parse_expr(doc_id)))
  message("Performing on original data ...")
  orig = getMeasures(data, lemma_names, slot_names, doc_id, measures)
  resultTemplate = matrix(nrow = nrow(orig$result), ncol = ncol(orig$result)) %>% data.frame
  colnames(resultTemplate) = colnames(orig$result)
  resultTemplate[,1:2] = orig$result[,1:2]
  results = lapply(1:length(quantiles), function(x) resultTemplate)

  message("Now performing resampling and calculating measures ...")
  measuresResampled = resampleFunc(data, doc_id, R, getMeasures, sampleRatio, lemma_names, slot_names, measures)

  percentiles_all = list()
  existPerc = numeric(nrow(orig$results))
  pb = txtProgressBar(min = 0, max = nrow(orig$results), initial = 0)

  message("Now extracting results ...")
  for(i in 1:nrow(orig$results)){
    setTxtProgressBar(pb,i)

    currOrigRow = orig$results[i,]
    exists = sapply(measuresResampled, function(x){
      x %>% filter(!!parse_expr(lemma_names[1]) == currOrigRow[[lemma_names[1]]],
                   !!parse_expr(lemma_names[2]) == currOrigRow[[lemma_names[2]]]) %>%
        nrow > 0
    })
    existPerc[i] = mean(exists)
    if(sum(exists) > 0){
      currMeasures = sapply(which(exists), function(x){
        measuresResampled[[x]] %>% filter(!!parse_expr(lemma_names[1]) == currOrigRow[[lemma_names[1]]],
                                        !!parse_expr(lemma_names[2]) == currOrigRow[[lemma_names[2]]]) %>%
          select(-!!parse_expr(lemma_names[1]),
                 -!!parse_expr(lemma_names[2]))
      })

      percentiles = sapply(1:nrow(currMeasures), function(k){
        quantile(currMeasures[k,] %>% unlist, quantiles, na.rm = T)
      }) %>% matrix(nrow = length(quantiles)) %>% t
      rownames(percentiles) = rownames(currMeasures)
      percentiles_all[[i]] = data.frame(percentiles)
      for(j in 1:length(quantiles)){
        results[[j]][i, 3:ncol(orig$results)] = as.list(percentiles[,j])
      }
    }
  }
  close(pb)
  if(!retainResamples){
    results
  } else {
    list(resampleResults = measuresResampled,
         quantiles = results)
  }
}




stratBoot = function(data, doc_id, R, getMeasures, sampleRatio,
                     lemma_names, slot_names, measures){
  data_split = data %>%
    group_split(!!parse_expr(doc_id))
  nDocs = length(data_split)

  measuresResampled = lapply(1:R, function(i){
    drawn = sample(nDocs, round(nDocs * sampleRatio), replace = TRUE)
    #Give each draw its own doc ID so that a doc drawn twice counts as two docs
    data_new = lapply(seq_along(drawn), function(k){
      data_split[[drawn[k]]] %>%
        mutate("{doc_id}" := paste0(!!parse_expr(doc_id), "_", k))
    }) %>% bind_rows
    getMeasures(data_new, lemma_names, slot_names, doc_id, measures)$results
  })

  # for(i in 1:R){
  #   data_new = data_split[ceiling(runif(length(docs), 0, length(docs)))] %>% bind_rows
  #   measuresResampled[[i]] = getMeasures(data_new, lemma_names, slot_names, doc_id, measures)
  # }
  measuresResampled
}
