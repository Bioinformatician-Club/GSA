#'@title count genotype of specified genes
#'@description the function is used to count the genotype of specified genes inputed by the user. A reference of which geneset has a higher impact on the disease.
#' @param Table  a dataframe contains a colume of gene names and columes of population genetype which have been processed into numeric form.
#' @param GeneSet  The name of gene(can be more than one) you are interested in.
#' @param Name  the name of geneset you come up with.
#'
#' @return The return of the function is a table containing two columes named by "Sample ID" and user's input "Name". The "Name" colume dislpays the sum of genotypes of the genes offered by the user.
#' @importFrom rlang .data
#' @examples
count_vars_in_gene_set_per_sample <- function(table, gene.set, gene.set.name) {

  if(tibble::is_tibble(table) == FALSE){
    stop("The input table needs to be a tibble format")
  }

  table %>%
    dplyr::filter(dplyr::if_all(1, ~ . %in% gene.set)) %>%
    dplyr::summarise(dplyr::across(2:ncol(table),base::sum)) %>%
    tidyr::pivot_longer(dplyr::everything(), names_to = "sampleID", values_to = gene.set.name)
}
