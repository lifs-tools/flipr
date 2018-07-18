a4 <- list(width = 8.27, height = 11.69)
a4r <- list(width = 11.69, height = 8.27)

dlnormMode <- function(meanlog, sdlog, shift=0) {
  exp(meanlog - (sdlog ^ 2)) - shift
}

#parameterized version of log-normal density function with additional "norm" expression and "scale" factor.
#both are multiplied with the actual log-normal density kernel.
dlnormPar <-
  function(x,
           meanlog = 0,
           sdlog = 1,
           scale = 1,
           shift = 0) {
    scale * dlnormMode(meanlog = meanlog, sdlog = sdlog, shift = shift) * dlnorm(x + shift,
                                                meanlog = meanlog,
                                                sdlog = sdlog,
                                                log = FALSE)
}
#' @importFrom magrittr %>%
writeGroup <- function(tibble, outputPrefix, skip = FALSE) {
  if (skip == FALSE) {
    outData <-
      dplyr::select(tibble,
                    combinationId,
                    precursorCollisionEnergy,
                    scanRelativeIntensity)
    readr::write_csv(outData, path = file.path(paste(
      make.names(unique(tibble$combinationId)), outputPrefix, ".csv", sep = "-"
    )))
  }
  return(tibble)
}

#' @importFrom magrittr %>%
#' @export
fits <- function(tibble, outputPrefix, group, instrumentId, skipGroupOutput = TRUE, start_lower, start_upper, lower, upper) {
  message(paste(
    "Fitting scanRelativeIntensities of fragments for:",
    instrumentId,
    outputPrefix,
    group,
    sep = " "
  ))
  message("Creating nls.tibbleId")
  nls.tibbleId <- tibble %>%
    dplyr::mutate(polarity = replace(polarity, polarity == "POSITIVE", "(+)")) %>%
    dplyr::mutate(polarity = replace(polarity, polarity == "NEGATIVE", "(-)")) %>%
    tidyr::unite(
      combinationId,
      species,
      fragment,
      adduct,
      polarity,
      `foundMassRange[ppm]`,
      group,
      sep = "|",
      remove = FALSE
    )

  message("Creating nls.tibble")
  nls.tibble <-
    nls.tibbleId %>% dplyr::group_by(., combinationId) %>% dplyr::mutate(samplesPerCombinationId = dplyr::n()) %>% dplyr::do(writeGroup(., outputPrefix, skip =
                                                                    skipGroupOutput))
  message("Writing data for fit")
  readr::write_tsv(nls.tibble, path = file.path(paste0(
    outputPrefix, "-data-for-fit.tsv"
  )))

  message("Creating nls.tibble.nested")
  nls.tibble.nested <- nls.tibble %>%
    tidyr::nest()
  # run nls fits with automatics model selection based on AIC
  message("Running nls.multstart on nls.tibble.nested")
  fits <- nls.tibble.nested %>%
    dplyr::mutate(fit = purrr::map(
      data,
      ~ nls.multstart::nls_multstart(
        scanRelativeIntensity ~ flipr::dlnormPar(precursorCollisionEnergy, meanlog, sdlog, scale, shift),
        data = .x,
        iter = 500,
        start_lower = start_lower,
        start_upper = start_upper,
        supp_errors = 'Y',
        na.action = na.omit,
        lower = lower,
        upper = upper
      )
    ,safely(x, otherwise = NA)),
    fitFun = "dlnormPar")
  message(paste("# of fits:", nrow(fits), sep = " "))
  stopifnot(length(fits) > 0)

  message("Extracting fit info")
  print(fits)
  # get fit information / statistics
  info <- fits %>%
    tidyr::unnest(fit %>% purrr::map(broom::glance))

  message("Extracting fit parameters")
  # get fit parameters
  params <- fits %>%
    tidyr::unnest(fit %>% purrr::map(broom::tidy))

  message("Calculating confidence intervals")
  # calculate confidence intervals for parameters
  CI <- fits %>%
    tidyr::unnest(fit %>% purrr::map(
      ~ nlstools::confint2(.x) %>%
        data.frame() %>%
        dplyr::rename(., conf.low = X2.5.., conf.high = X97.5..)
    )) %>%
    dplyr::group_by(., combinationId) %>%
    dplyr::mutate(., term = c('meanlog', 'sdlog', 'scale', 'shift')) %>%
    dplyr::ungroup()

  message("Creating LipidCreator parameters")
  # create output for lipidcreator import
  lipidCreatorParams <- params %>%
    tidyr::separate(
      combinationId,
      c("species", "fragment", "adduct", "polarity", "ppmMassRange", "group"),
      sep = "\\|",
      remove = TRUE
    ) %>% dplyr::rename(
      class = species,
      frag = fragment,
      model = fitFun,
      ParKey = term,
      ParValue = estimate
    ) %>% dplyr::mutate(instrument = instrumentId,
                        ppmMassRange = as.numeric(ppmMassRange)) %>%
    dplyr::select(
      "instrument",
      "class",
      "frag",
      "adduct",
      "ppmMassRange",
      "group",
      "model",
      "ParKey",
      "ParValue"
    ) %>% dplyr::filter(ppmMassRange == min(ppmMassRange))

  # TODO additional renormalization after grouping for each "species", adduct, polarity and ppmMassRange -> scale needs to fall into the range of 0-1

  message("Writing LipidCreator parameters")
  readr::write_tsv(lipidCreatorParams, path = file.path(paste0(
    outputPrefix, "-lipidcreator-parameters.tsv"
  )))

  # merge parameters and CI estimates
  message("Merging fit parameters and CI estimates")
  params <- merge(params, CI, by = dplyr::intersect(names(params), names(CI)))
  # params <- params %>% mutate(combinationId = paste(fragment, adduct, polarity, sep = "-", collapse=T) )

  message("Writing merged parameters")
  readr::write_tsv(params, path = file.path(paste0(outputPrefix, "-parameters.tsv")))

  message("Calculating predictions from data")
  preds_from_data <- fits %>%  tidyr::unnest(fit %>% purrr::map(broom::augment))

  message("Creating new x-values for predictions")
  # use tibble with combinationId for merging and create x-values for each fit (covering the complete x-axis ranges)
  new_preds <- nls.tibbleId %>%
    dplyr::do(.,
       data.frame(
         precursorCollisionEnergy = seq(min(0), max(.$precursorCollisionEnergy), by =
                                          1),
         stringsAsFactors = FALSE
       ))

  message("Calculating min/max x-value range for each fit")
  # calculate the minimum and maximum range of x-values for each unique combination (fit) and add some slack
  max_min <- dplyr::group_by(nls.tibbleId, combinationId) %>%
    dplyr::summarise(
      .,
      min_pce = plyr::round_any(min(precursorCollisionEnergy), 10, f = floor) - 1,
      max_pce = plyr::round_any(max(precursorCollisionEnergy), 10, f = ceiling) + 1
    ) %>%
    dplyr::ungroup()

  message("Recalculating predictions based on fit parameters")
  # recalculate predictions based on fit parameters
  preds <- fits %>%
    tidyr::unnest(fit %>% purrr::map(augment, newdata = new_preds)) %>%
    merge(., max_min, by = "combinationId") %>%
    dplyr::group_by(., combinationId) %>%
    dplyr::filter(
      .,
      precursorCollisionEnergy > unique(min_pce) &
        precursorCollisionEnergy < unique(max_pce)
    ) %>%
    dplyr::rename(., scanRelativeIntensity = .fitted) %>%
    dplyr::ungroup()

  message("Extracting fit details")
  # get details of fits
  fitinfo <-
    info %>% dplyr::select(combinationId, fitFun, sigma, isConv, finTol, logLik, AIC, BIC)
  message("Writing fit predictions")
  readr::write_tsv(preds, path = file.path(paste(outputPrefix, "-predictions.tsv", sep =
                                            "")))
  message("Writing fit info")
  readr::write_tsv(fitinfo, path = file.path(paste(outputPrefix, "-fit-info.tsv", sep =
                                              "")))
  list(
    fits = fits,
    params = params,
    CI = CI,
    preds = preds,
    nls.tibble = nls.tibble,
    preds_from_data = preds_from_data
  )
}
