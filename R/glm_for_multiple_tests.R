#'@title Use modeling to calculate correlation
#'
#'@description With covariants and genesets, this function would calculate the necessary values that uncovers the correlation between genesets and outcomes.
#' @param ped A raw dataframe from the bioinformation database.
#' @param outcome The type of outcome to be tested.
#' @param geneset The genesets(can be more than one)that users want to test the correlation with the outcome.
#' @param covariates The covariates(e.g. age, gender)help testing the correlation.
#'
#' @return The return of this function is a table with p-value, conf.high and conf.high that give users direction on which geneset has greatest correlation to the disease.
#' @importFrom rlang .data
#' @examples
glm_for_multiple_tests <- function(ped, outcome, interested.variables, covariates, tidy.output = FALSE, conf.int = FALSE, conf.level = 0.95, exponentiate = FALSE){

  ped.slim <- dplyr::select(ped, all_of(outcome), all_of(interested.variables), all_of(covariates))

  nest.data <- ped.slim %>%
    tidyr::pivot_longer(col = all_of(interested.variables), names_to = "target.variable", values_to = "value") %>%
    dplyr::group_by(target.variable) %>%
    tidyr::nest()

  nest.data.with.model <- nest.data %>%
    dplyr::mutate(model = purrr::map(data, ~stats::glm(as.formula(base::paste0(outcome,"~ .")), data = ., family = "binomial")))

  if(tidy.output) {

    nest.data.with.model.tidy <- nest.data.with.model %>%
      dplyr::mutate(summary = purrr::map(model, ~broom::tidy(.,conf.int = conf.int, conf.level = conf.level, exponentiate = exponentiate) %>% dplyr::filter(term == "value"))) %>%
      dplyr::select(target.variable, summary) %>%
      tidyr::unnest(summary) %>%
      dplyr::select(-term) %>%
      dplyr::ungroup()

    return(nest.data.with.model.tidy)

  }else{
    return(nest.data.with.model)
  }
}
