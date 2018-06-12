# word is what we want from the go database
# returns IDs associated with given word
goTerms <- function(word) {
    go_terms <- Term(GO.db::GOTERM)
    pattern <- paste(".*", word, ".*", sep = "")
    return(names(go_terms[grepl(pattern, go_terms)]))
}

goIDs <- function(root_dir, ensembl_gene_ids, word, force = FALSE) {
    if(file_test("-f", file.path(root_dir, paste(word, "go.rds", sep = "_"))) && !force) {
        res <- readRDS(file.path(root_dir, paste(word, "go.rds", sep = "_")))
        return(res)
    }
    # remove version from gene id
    ensembl_gene_ids <- lapply(ensembl_gene_ids, 
                               function(x) 
                                   return(sub("\\..*", "", x)))
    
    # get the go ids we want to keep
    keep_ids <- goTerms(word)
    
    # incProgress(1/3, message = "Querying ENSEMBL")
    mart <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
    res <- biomaRt::getBM(mart = mart, 
                 attributes = c("ensembl_gene_id", "go_id"), 
                 filters = "ensembl_gene_id", 
                 values = ensembl_gene_ids)
    # incProgress(1/3, message = "Filtering GO Terms")
    res <- dplyr::filter(res, res$go_id != "")
    res <- dplyr::filter(res, res$go_id %in% keep_ids)
    res$go_term <- lapply(res$go_id, 
                          function(x)
                              return(Term(GOTERM[x])))
    saveRDS(res, file = file.path(root_dir, paste(word, "go.rds", sep = "_")))
    return(res)
}

# goDesc <- function(ensembl_gene_ids) {
#     ensembl_gene_ids <- lapply(ensembl_gene_ids, 
#                                function(x) 
#                                    return(sub("\\..*", "", x)))
#     mart <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
#     res <- biomaRt::getBM(mart = mart, 
#                           attributes = c("ensembl_gene_id", "go_id", ), 
#                           filters = "ensembl_gene_id", 
#                           values = ensembl_gene_ids)
# }