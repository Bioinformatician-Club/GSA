#' Run logistic regressions for multiple independent variables respectively at once
#'
#' @param ped
#' A data frame with a PED format supported by PLINK software to represent phenotype/outcome and covariates for the logistic regression. Running \code{data("1000g")} to load an example of ped file.
#' The input ped must be assigned with correct class/type for each variable as we normally do when running logistic regression in R. e.g.
#' * the outcome/dependent variable must be a binary variable, either 0 or 1.
#' * according to the natural property of some independent variables, they should be assigned with a proper class/type. e.g. usually the variable representing sex assigned as a factor.
#'
#' @param outcome
#' A string specifying the column name representing outcome in ped file.
#' @param interested.variables
#' A vector specifying column names of all your interested variables in \code{ped} for running logistic regressions in parallel with the same \code{covariates}.
#' @param covariates
#' A vector specifying column names of all covariates in \code{ped} which will be used in logistic regression.
#' @param tidy.output
#' Logical indicating whether or not to print the output in the tidied format. Defaults to FALSE. We recommend user to choose FALSE if user wanna obtain all parameters from logistic regression model.
#' @param conf.int
#' Logical indicating whether or not to include a confidence interval in the tidied output. Defaults to FALSE.
#' @param conf.level
#' The confidence level to use for the confidence interval if conf.int = TRUE. Must be strictly greater than 0 and less than 1. Defaults to 0.95, which corresponds to a 95 percent confidence interval.
#' @param exponentiate
#' Logical indicating whether or not to exponentiate the the coefficient estimates. This is typical for logistic regression. Defaults to FALSE.
#'
#' @return
#' If tidy.output is TRUE, this function return a tibble with columns:
#'
#' \code{interested.variables}
#'
#' \code{estimate}    The estimated value of the regression term.
#'
#' \code{std.error}   The standard error of the regression term.
#'
#' \code{statistic}   The value of a T-statistic to use in a hypothesis that the regression term is non-zero.
#'
#' \code{p.value}     The two-sided p-value associated with the observed statistic.
#'
#'
#' If tidy.output is TRUE, this function return a tibble with columns:
#'
#' \code{interested.variables}
#'
#' \code{data}        The raw data used in logistic regression
#'
#' \code{model}       The original output from the logistic regression.
#'
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
#'
#' library(tidyverse)
#'
#' # add one column representing the allele count in the gene set into ped file
#' ped.with.geneset.count <- ped %>%
#' as_tibble() %>%
#' left_join(geneset.ac, by = c("IND_ID" = "sampleID")) %>%
#' replace_na(list(genesetA = 0)) %>%
#' mutate(DISEASE = DISEASE - 1,
#'       SEX = factor(SEX))
#'
#' glm.results.tidy <- glm_for_multiple_tests(ped = ped.with.geneset.count, outcome = "DISEASE", interested.variables = c("genesetA", "AGE"), covariates= c("QT","SEX"), tidy.output = TRUE, conf.int = TRUE, exponentiate = F)
#'
#' glm.results.tidy
#'
#' glm.results <- glm_for_multiple_tests(ped = ped.with.geneset.count, outcome = "DISEASE", interested.variables = c("genesetA", "AGE"), covariates= c("QT","SEX"))
#'
#' glm.results
#' }
#'
glm_for_multiple_tests <- function(ped, outcome, interested.variables, covariates, tidy.output = FALSE, conf.int = FALSE, conf.level = 0.95, exponentiate = FALSE){

  interested.variables.tmp <- interested.variables
  ped.slim <- dplyr::select(ped, all_of(outcome), all_of(interested.variables.tmp), all_of(covariates))

  nest.data <- ped.slim %>%
    tidyr::pivot_longer(col = all_of(interested.variables.tmp), names_to = "interested.variables", values_to = "value") %>%
    dplyr::group_by(interested.variables) %>%
    tidyr::nest()

  nest.data.with.model <- nest.data %>%
    dplyr::mutate(model = purrr::map(data, ~stats::glm(as.formula(base::paste0(outcome,"~ .")), data = ., family = "binomial")))

  if(tidy.output) {

    nest.data.with.model.tidy <- nest.data.with.model %>%
      dplyr::mutate(summary = purrr::map(model, ~broom::tidy(.,conf.int = conf.int, conf.level = conf.level, exponentiate = exponentiate) %>% dplyr::filter(term == "value"))) %>%
      dplyr::select(interested.variables, summary) %>%
      tidyr::unnest(summary) %>%
      dplyr::select(-term) %>%
      dplyr::ungroup()

    return(nest.data.with.model.tidy)

  }else{
    return(nest.data.with.model)
  }
}
