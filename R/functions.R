#' Descriptive stats for Lipidomic project data
#'
#' @param data
#'
#' @returns table with mean and SD
create_table_descriptive_stats <- function(data) {
  data |>
    dplyr::group_by(metabolite) |>
    dplyr::summarise(
      dplyr::across(value, list(mean = mean, sd = sd))
    ) |>
    dplyr::mutate(
      dplyr::across(where(is.numeric), \(x) round(x, digits = 1))
    ) |>
    dplyr::mutate(
      MeanSD = glue::glue("{value_mean} ({value_sd})")
    ) |>
    dplyr::select(
      Metabolite = metabolite,
      `Mean SD` = MeanSD
    )
}

#' Descriptive stats for Lipidomic project data
#'
#' @param data
#'
#' @returns ggplot hist
create_plot_distributions <- function(data) {
  data |>
    ggplot2::ggplot(ggplot2::aes(x = value)) +
    ggplot2::geom_histogram() +
    ggplot2::facet_wrap(~metabolite, scales = "free") +
    ggplot2::theme_minimal()
}


#' Do some cleaning to fix issues in the data.
#'
#' @param data The lipidomics data frame.
#'
#' @returns A data frame.
#'
clean_Cholesterol_dublicates <- function(data) {
  data |>
    dplyr::group_by(dplyr::pick(-value)) |>
    dplyr::summarise(value = mean(value), .groups = "keep") |>
    dplyr::ungroup()
}


#' Fix data to process it for model fitting.
#'
#' @param data The lipidomics data.
#'
#' @returns A data frame.
#'
preprocess <- function(data) {
  data |>
    dplyr::mutate(
      class = as.factor(class),
      value = scale(value)
    )
}


#' Do model fitting.
#'
#' @param my_data The lipidomics data.
#' @param my_formula The model formula
#'
#' @returns A data frame.
#'
fit_model <- function(my_data, my_formula) {
  glm(
    formula = my_formula,
    data = my_data,
    family = binomial
  ) |>
    broom::tidy(exponentiate = TRUE) |>
    dplyr::mutate(
      metabolite = unique(my_data$metabolite),
      model = format(my_formula),
      .before = everything()
    )
}

#' Get tidy model results
#'
#' @param data The lipidomics data.
#'
#' @returns A data frame.
#'
create_model_results <- function(data) {
  data |>
    dplyr::filter(metabolite == "Cholesterol") |>
    preprocess() |>
    fit_model(class ~ value)
}
