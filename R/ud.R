

getCombosFromUD = function(conllu, relation = expr(node() %any_down% node())){

}


tracePathInTree = function(nodeList, from, to){
  nodeList[from]
}



# getRowsFromCombos = function(combos){
#   n_sfp = length(combos$sfp)
#   n_dc = length(combos$dc)
#   if(n_sfp == 0) combos$sfp = NA
#   if(n_dc == 0) combos$dc = NA
#   expand_grid(sfp = combos$sfp, dc = combos$dc) %>%
#     mutate(docID = combos$docID,
#            paragraphID = combos$paragraphID,
#            sentenceID = combos$sentenceID,
#            sentence = combos$sentence)
# }

getTreeFromConllu = function(conllu){
  conllu = conllu %>% mutate(token_id = as.numeric(token_id), head_token_id = as.numeric(head_token_id))
  #Numeric because some corpora have decimal values

  #Start with the root
  currTokenID = (conllu %>% filter(.data$head_token_id == 0) %>% pull(.data$token_id))[1]
  currTokenLemma = (conllu %>% filter(.data$head_token_id == 0) %>% pull(.data$lemma))[1]
  root = Node$new(currTokenID, lemma = currTokenLemma)
  getDescendentsInConllu(conllu, currTokenID, root)
  nodesList = Traverse(root)
  nodesList = nodesList[sapply(nodesList, function(x) as.numeric(x$name))]
  list(df = conllu, root = root, nodesList = nodesList)
}

getTidyInputFromComboList = function(comboList){

}

getDescendentsInConllu = function(conllu, currTokenID, currNode){
  children = conllu %>% filter(.data$head_token_id == currTokenID) %>% pull(.data$token_id)
  for(currTokenID in children){
    currTokenLemma = (conllu %>% filter(head_token_id == 0) %>% pull(.data$lemma))[1]
    nextNode = currNode$AddChild(currTokenID, lemma = currTokenLemma)
    getDescendentsInConllu(conllu, currTokenID, nextNode)
  }
}

getParent = function(tree, deprel, order = "either"){

}

getChildren = function(tree, deprel, order = "either"){

}

getDescendents = function(trees, deprel, order = "either"){

}

getLineage = function(tree, deprel, order = "either"){

}
