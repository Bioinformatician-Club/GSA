#' Count gene mutations
#'
#' @param table A table containing the information of each and every gene the samples have been tested on.
#' @param gene.col The number of the gene column when counting from the left.
#' @param gt.col The number of the first sample column when counting from the left.
#' @param type "Allele", "Snps", or "Carrier", indicating the table to be generated
#'
#' @return
#' type = Alleles, A table containing the number of mutated alleles for each gene in each sample.
#' type = Snps, A table containing the number of mutated locus in each gene for each sample.
#' type = Carrier, A table indicating the presence of mutation in each gene for each sample.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' count_vars_per_gene_per_sample(table = variants, gene.col = 8, gt.col = 12, type = "allele")
#' }

count_vars_per_gene_per_sample <- function (table, gene.col, gt.col, type = c("allele", "locus", "carrier")) {

  if(!(is.numeric(gene.col) && gene.col == round(gene.col))){
    stop("gene.col needs to be an integer.")
  }

  if(!(is.numeric(gt.col) && gt.col == round(gt.col))){
    stop("gene.col needs to be an integer.")
  }

  # the default type is "allele"
  type <- match.arg(type)

  # vectorizing functions is much faster
  gt2count <- function(genotypes){
    stringr::str_extract_all(genotypes, "\\d") %>% purrr::map_dbl(., ~as.numeric(.) %>% sum)
  }

  # convert the genotype table to allele count table
  table.allelecount <- table %>%
    dplyr::mutate(across(names(table)[gt.col]:tail(names(table),1),gt2count))

  # convert allele count table to carrier status table
  table.carrier <- table.allelecount %>%
    dplyr::mutate(across(names(table.allelecount)[gt.col]:tail(names(table.allelecount),1),~dplyr::if_else(. > 1, 1, .)))


  if(type == "allele") {

    gene.allelcount <- table.allelecount %>%
      dplyr::group_by_at(.vars = gene.col) %>%
      dplyr::summarise(dplyr::across(names(table.allelecount)[gt.col]:tail(names(table.allelecount),1), sum))

    return(gene.allelcount)

  }else {

    gene.locus <- table.carrier %>%
      dplyr::group_by_at(.vars = gene.col) %>%
      dplyr::summarise(dplyr::across(names(table.carrier)[gt.col]:tail(names(table.carrier),1), sum))

    if(type == "locus") {

      return(gene.locus)

    }else if(type == "carrier") {

      gene.carrier <- gene.locus %>%
        dplyr::mutate(dplyr::across(names(gene.locus)[2]:tail(names(gene.locus),1), ~dplyr::if_else(. > 1, 1, .)))

    return(gene.carrier)
    }
  }
}
