#'@title Count the total number of counting metric in a gene set for each sample
#'@description
#'This function count the total number of counting metric (i.e. alleles, locus or carrier-status from function count_vars_per_gene_per_sample) in a gene set for each sample.
#' @param table  A tibble which has first column representing gene name and the rest of columns representing the counting metrics of samples. An example is the output of the function count_vars_per_gene_per_sample.
#' @param gene.set  A vector of gene names specifying a gene set. Usually it is a subset of genes from the first column of input table.
#' @param gene.set.name  A string specifying the name of the gene set which will be used as the column name in the output
#'
#' @return
#' A tibble with two columns. The first column "sampleID" represents sample names and the second column named with *gene.set.name* represents the total number of metric in the input *gene.set*.
#' @importFrom rlang .data
#' @examples
#' \dontrun{
#' data("1000g")
#'
#' # generate allele count table for each gene and each sample
#' gene.ac <- count_vars_per_gene_per_sample(table = variants, gene.col = 8, gt.col = 12, type = "allele")
#'
#' # here just pick some genes and use as an example as a gene set input
#' gene.list <- gene.ac$GENE[10:14]
#'
#' # count the total number of alleles in the gene set for each sample
#' geneset.ac <- count_vars_in_gene_set_per_sample(table = gene.ac, gene.set = gene.list, gene.set.name = "genesetA")
#' }


count_vars_in_gene_set_per_sample <- function(table, gene.set, gene.set.name) {

  if(tibble::is_tibble(table) == FALSE){
    stop("The input table needs to be a tibble format")
  }

  table %>%
    dplyr::filter(dplyr::if_all(1, ~ . %in% gene.set)) %>%
    dplyr::summarise(dplyr::across(2:ncol(table),base::sum)) %>%
    tidyr::pivot_longer(dplyr::everything(), names_to = "sampleID", values_to = gene.set.name)
}
