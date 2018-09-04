#' #huber loss function
#' #' @export
#' huberLossPw <- function(x, hgamma = 0.05) {
#'   stopifnot(is.double(x))
#'   if (abs(x) <= hgamma) {
#'     0.5 * (x) ^ 2
#'   } else {
#'     hgamma * abs(x) - ((hgamma ^ 2) / 2)
#'   }
#' }
#'
#' #' @export
#' huberLoss <- function(x, hgamma = 0.05) {
#'   sapply(x, huberLossPw, hgamma = hgamma)
#' }
#'
#' huberLoss2 <- function(x1, x2, hgamma = 0.05) {
#'   huberLoss(x1 - x2, hgamma = hgamma)
#' }

#' Mode function of shifted \code{dlnorm}.
#' @param meanlog the \code{meanlog} argument for \code{dlnorm}.
#' @param sdlog the \code{sdlog} argument for \code{dlnorm}
#' @param shift the \code{shift} parameter modulates the shift of the log-normal density function or 0.
#' @return the value of the mode function for the given parameters.
#' @export
dlnormMode <- function(meanlog, sdlog, shift = 0) {
  exp(meanlog - (sdlog ^ 2)) - shift
}

#' Parameterized version of the log-normal density function \code{dlnorm} with additional "norm" expression and "scale" factor.
#' Both are multiplied with the actual log-normal density kernel.
#'
#' @param x the \code{x} argument for the log-normal density function
#' @param meanlog the \code{meanlog} argument for \code{dlnorm} or 0.
#' @param sdlog the \code{sdlog} argument for \code{dlnorm} or 1.
#' @param scale the \code{scale} parameter modulates the height of the log-normal density function or 1.
#' @param shift the \code{shift} parameter modulates the shift of the log-normal density function on the x-axis or 0.
#' @return The value of the log-normal density function for \code{x} and all parameters.
#' @export
dlnormPar <-
  function(x,
           meanlog = 0,
           sdlog = 1,
           scale = 1,
           shift = 0) {
    scale * dlnormMode(meanlog = meanlog,
                       sdlog = sdlog,
                       shift = shift) * dlnorm(x + shift,
                                               meanlog = meanlog,
                                               sdlog = sdlog,
                                               log = FALSE)
  }

#' Selects \code{combinationId}, \code{precursorCollisionEnergy} and \code{scanRelativeIntenty} to write to file.
#' @param tibble the data frame / tibble to write.
#' @param outputPrefix the output prefix for the file.
#' @param skip if true, skip writing.
#' @return the input tibble.
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

#' Add parameter for instrumentId (cvTerm)  and output format suitable
#' for LipidCreator import.
#' @param params the fit parameters to modify.
#' @param instrumentId the instrumentId cvTerm to use.
#' @return a tibble in the correct format for LipidCreator with ParKey and ParValue.
#' @importFrom magrittr %>%
createLipidCreatorParameters <-
  function(params, instrumentId = "MS:XXXXXX") {
    message("Creating LipidCreator parameters")
    lcModel <- unique(params[["fitFun"]])
    lcParams <- params %>% dplyr::select("combinationId",
                                         "term",
                                         "estimate")
    lcParams <-
      lcParams %>% dplyr::group_by(combinationId) %>% tidyr::unite(.,
                                                                   paramComb,
                                                                   term,
                                                                   estimate,
                                                                   sep = "|",
                                                                   remove = TRUE)
    lcParams$combinationId <- as.character(lcParams$combinationId)
    lcParams$paramComb <- as.character(lcParams$paramComb)
    modelFrame <-
      tibble::as.tibble(data.frame(paramComb = rep(
        paste0("model", "|", lcModel), length.out = length(unique(lcParams$combinationId))
      )))
    modelFrame$combinationId <-
      as.character(unique(lcParams$combinationId))
    modelFrame$paramComb <- as.character(modelFrame$paramComb)
    lcParams <- dplyr::bind_rows(lcParams, modelFrame)

    lipidCreatorParams <- lcParams %>%
      tidyr::separate(
        combinationId,
        c(
          "species",
          "fragment",
          "adduct",
          "polarity",
          "calculatedMass",
          "ppmMassRange",
          "group"
        ),
        sep = "\\|",
        remove = TRUE
      ) %>%
      tidyr::separate(paramComb,
                      c("ParKey", "ParValue"),
                      sep = "\\|",
                      remove = TRUE) %>%
      dplyr::rename(class = species, frag = fragment) %>%
      dplyr::mutate(instrument = instrumentId,
                    ppmMassRange = as.numeric(ppmMassRange)) %>%
      dplyr::select(
        "instrument",
        "class",
        "frag",
        "adduct",
        #"calculatedMass",
        "ppmMassRange",
        "group",
        #      "model",
        "ParKey",
        "ParValue"
      ) %>% dplyr::filter(ppmMassRange == min(ppmMassRange))

    #print(lipidCreatorParams)

    lipidCreatorParams <-
      lipidCreatorParams %>%
      dplyr::group_by(., instrument, class, frag, adduct, group) %>%
      dplyr::arrange(., instrument, class, frag, adduct, group, ParKey, ParValue)
    return(lipidCreatorParams)
  }

#' Calculates the nonlinear fits for the given tibble. Data will be grouped by
#' combinationId, species, fragment, adduct, polarity, calculatedMass, foundMassRange[ppm], and group columns.
#' Calculates non linear regression models by performing an interative grid search within the coordinate bounds provided
#' by \code{lower} and \code{upper} vectors, starting at \code{start_lower} and \code{start_upper}.
#' Intermediate models are scored by the AIC value until convergence has been achieved as to the defaults of nls.multstart.
#'
#' @param tibble the data tibble containing fragment data.
#' @param outputPrefix the output prefix for the file.
#' @param group the group that this data belongs to.
#' @param instrumentId the instrumentId of this data.
#' @param skipGroupOutput whether to skip writing of grouping information (for debugging).
#' @param start_lower the lower bound to start the parameter grid search.
#' @param start_upper the upper bound to start the parameter grid search.
#' @param lower the lower bound for the parameter grid search.
#' @param upper the upper bound for the parameter grid search.
#' @param minSamplesPerCombinationId the minimum number of data points required per fragment / adduct / ppm combination to be considered for model calculation.
#' @param max_iter the maximum number of iterations of the model to calculate.
#' @importFrom magrittr %>%
#' @export
fits <-
  function(tibble,
           outputPrefix,
           group,
           instrumentId,
           skipGroupOutput = TRUE,
           start_lower,
           start_upper,
           lower,
           upper,
           minSamplesPerCombinationId = 50,
           max_iter = 500) {
    message(
      paste(
        "Fitting scanRelativeIntensities of fragments for:",
        instrumentId,
        outputPrefix,
        group,
        sep = " "
      )
    )
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
        calculatedMass,
        `foundMassRange[ppm]`,
        group,
        sep = "|",
        remove = FALSE
      )

    message("Creating nls.tibble")
    nls.tibble <-
      nls.tibbleId %>% dplyr::group_by(., combinationId) %>% dplyr::mutate(samplesPerCombinationId = dplyr::n()) %>% dplyr::do(writeGroup(., outputPrefix, skip =

                                                                                                                                            skipGroupOutput))
    message("Writing unfiltered data before fit")
    nls.tibble.unfiltered <- nls.tibble
    combinations.unfiltered <- length(unique(nls.tibble.unfiltered$combinationId))
    readr::write_tsv(nls.tibble.unfiltered, path = file.path(paste0(
      outputPrefix, "-data-for-fit-unfiltered.tsv"
    )))

    nls.tibble <-
      nls.tibble %>% dplyr::filter(samplesPerCombinationId >= minSamplesPerCombinationId)
    nrow.removed <- combinations.unfiltered - length(unique(nls.tibble$combinationId))
    message(paste("Filtered", nrow.removed, "cases from data where samplesPerCombinationId <",minSamplesPerCombinationId))
    message(paste("Remaining data rows:", nrow(nls.tibble)))
    if(nrow(nls.tibble)==0) {
      stop("No rows remaining for model calculation after filtering!")
    }
    message("Writing data for fit")
    #nls.tibble$weights <- 1/nls.tibble$sriVarPerCombinationId
    readr::write_tsv(nls.tibble, path = file.path(paste0(outputPrefix, "-data-for-fit.tsv")))

    message("Creating nls.tibble.nested")
    nls.tibble.nested <- nls.tibble %>%
      tidyr::nest()
    # run nls fits with automatics model selection based on AIC
    message(paste("Running nls.multstart on nls.tibble.nested with at most", max_iter, "iterations"))
    #ls( environment() )
    fits <- nls.tibble.nested %>%
      dplyr::mutate(fit = purrr::map(
        data,
        ~ nls.multstart::nls_multstart(
          scanRelativeIntensity ~ flipr::dlnormPar(precursorCollisionEnergy, meanlog, sdlog, scale, shift),
          data = .x,
          iter = max_iter,
          start_lower = start_lower,
          start_upper = start_upper,
          supp_errors = 'Y',
          na.action = na.omit,
          lower = lower,
          upper = upper
        )
        ,
        safely(x, otherwise = NA)
      ),
      fitFun = "dlnormPar")
    stopifnot(length(fits) > 0)

    message(paste("Extracting fit info for", nrow(fits), "fits", sep = " "))
    print(fits)
    # get fit information / statistics
    #options(error = dump.frames(to.file=TRUE))
    info <- fits %>%
      tidyr::unnest(fit %>% purrr::map(broom::glance))
    message("Extracting fit parameters")
    # get fit parameters
    params <- fits %>%
      tidyr::unnest(fit %>% purrr::map(broom::tidy))

    message("Writing raw parameters")
    readr::write_tsv(params, path = file.path(paste0(outputPrefix, "-raw-parameters.tsv")))

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

    # check that we only have a single model
    stopifnot(length(unique(params$fitFun)) == 1)

    lipidCreatorParams <-
      createLipidCreatorParameters(params, instrumentId)
    # TODO additional renormalization after grouping for each "species", adduct, polarity and ppmMassRange -> scale needs to fall into the range of 0-1
    message("Writing LipidCreator parameters")
    readr::write_csv(lipidCreatorParams, path = file.path(paste0(
      outputPrefix, "-lipidcreator-parameters.csv"
    )))

    # merge parameters and CI estimates
    message("Merging fit parameters and CI estimates")
    params <-
      merge(params, CI, by = dplyr::intersect(names(params), names(CI)))
    # params <- params %>% mutate(combinationId = paste(fragment, adduct, polarity, sep = "-", collapse=T) )

    message("Writing merged parameters")
    readr::write_tsv(params, path = file.path(paste0(outputPrefix, "-parameters.tsv")))

    message("Calculating predictions from data")
    preds_from_data <-
      fits %>%  tidyr::unnest(fit %>% purrr::map(broom::augment))

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
      tidyr::unnest(fit %>% purrr::map(broom::augment, newdata = new_preds)) %>%
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
      nls.tibble.unfiltered = nls.tibble.unfiltered,
      nls.tibble = nls.tibble,
      preds_from_data = preds_from_data
    )
  }
