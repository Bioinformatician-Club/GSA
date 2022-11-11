#'@title Generate the metric table per gene per sample
#'@description
#'This function count the total number of alleles, the total number of locus, or the status of carrying any mutations per gene per sample depending on which type of counting metric you choose.
#' @param table A data frame which contains information of variants first and then followed by genotypes of all samples. The information of variants, annotated gene and sample names are required; Chrom, pos, ID, ref, alt, effect, impact, etc are optional.
#' @param gene.col The position of column representing annotated gene.
#' @param gt.col The position of column representing the genotype of the first sample. By default, all following columns after gt.col are genotypes of the rest samples.
#' @param type Control which type of counting metric to use.
#' * "allele" count the total number of allele count per gene per sample
#' * "locus" count the total number of genetic positions with mutation per gene per sample
#' * "carrier", return the status of the sample carrying at least one mutation for each gene and each sample. 1 for carrier, 0 for non-carrier.
#'
#' @return
#' A data frame which contains the count of allele/locus/carrier-status per gene per sample. The first column represents annotated gene, and the rest of columns represent the count of allele/locus/carrier-status of all samples, depending on which the parameter of *type* you choose
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data("1000g")
#'
#' # generate allele count table for each gene and each sample
#' gene.ac <- count_vars_per_gene_per_sample(table = variants, gene.col = 8, gt.col = 12, type = "allele")
#'
#' # generate locus count table for each gene and each sample
#' gene.locus <- count_vars_per_gene_per_sample(table = variants, gene.col = 8, gt.col = 12, type = "locus")
#'
#' # generate carrier status table for each gene and each sample
#' gene.carrier <- count_vars_per_gene_per_sample(table = variants, gene.col = 8, gt.col = 12, type = "carrier")
#'
#' }

count_vars_per_gene_per_sample <- function (table, gene.col, gt.col, type = c("allele", "locus", "carrier")) {

  if(!(is.numeric(gene.col) && gene.col == round(gene.col))){
    stop("gene.col needs to be an integer.")
  }

  if(!(is.numeric(gt.col) && gt.col == round(gt.col))){
    stop("gt.col needs to be an integer.")
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
