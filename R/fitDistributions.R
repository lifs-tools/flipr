library(fitdistrplus)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(survival)
library(MASS)
library(nls.multstart)
library(broom)
library(purrr)
library(tidyr)
library(nlstools)

a4 <- list(width = 8.27, height = 11.69)
a4r <- list(width = 11.69, height = 8.27)

dlnormMode <- function(meanlog, sdlog, shift) {
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
    scale * dlnormMode(meanlog, sdlog) * dlnorm(x + shift,
                                                meanlog = meanlog,
                                                sdlog = sdlog,
                                                log = FALSE)
}

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

fits <- function(tibble, outputPrefix, skipGroupOutput = TRUE) {
  print(paste(
    "Fitting scanRelativeIntensities of fragments for:",
    outputPrefix,
    sep = " "
  ))
  nls.tibbleId <- tibble %>%
    dplyr::mutate(polarity = dplyr::replace(polarity, polarity == "POSITIVE", "(+)")) %>%
    dplyr::mutate(polarity = dplyr::replace(polarity, polarity == "NEGATIVE", "(-)")) %>%
    tidyr::unite(
      combinationId,
      species,
      fragment,
      adduct,
      polarity,
      `foundMassRange[ppm]`,
      sep = "|",
      remove = FALSE
    )

  nls.tibble <-
    nls.tibbleId %>% dplyr::group_by(., combinationId) %>% dplyr::mutate(samplesPerCombinationId = dplyr::n()) %>% do(writeGroup(., outputPrefix, skip =
                                                                    skipGroupOutput))
  readr::write_tsv(nls.tibble, path = file.path(paste0(
    outputPrefix, "-data-for-fit.tsv"
  )))


  nls.tibble.nested <- nls.tibble %>%
    tidyr::nest()
  # run nls fits with automatics model selection based on AIC
  fits <- nls.tibble.nested %>%
    dplyr::mutate(fit = purrr::map(
      data,
      ~ nls.multstart::nls_multstart(
        scanRelativeIntensity ~ flipr::dlnormPar(precursorCollisionEnergy, meanlog, sdlog, scale, shift),
        data = .x,
        iter = 500,
        start_lower = c(
          meanlog = -10,
          sdlog = 0.01,
          scale = min(.x$scanRelativeIntensity),
          shift = -200
        ),
        start_upper = c(
          meanlog = 10,
          sdlog = 10,
          scale = max(.x$scanRelativeIntensity),
          shift = 200
        ),
        supp_errors = 'Y',
        na.action = na.omit,
        lower = c(
          meanlog = -20,
          sdlog = 0.0001,
          scale = 0.000001,
          shift = -1000
        ),
        upper = c(
          meanlog = 20,
          sdlog = 20,
          scale = 5,
          shift = 1000
        )
      )
    ),
    fitFun = "dlnormPar")
  print(paste("Length of fits:", length(fits), sep = " "))
  stopifnot(length(fits) > 0)

  # get fit information / statistics
  info <- fits %>%
    tidyr::unnest(fit %>% purrr::map(glance))

  # get fit parameters
  params <- fits %>%
    tidyr::unnest(fit %>% purrr::map(tidy))

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

  # create output for lipidcreator import
  lipidCreatorParams <- params %>%
    tidyr::separate(
      combinationId,
      c("species", "fragment", "adduct", "polarity", "ppmMassRange"),
      sep = "\\|",
      remove = TRUE
    ) %>% dplyr::rename(
      class = species,
      frag = fragment,
      model = fitFun,
      ParKey = term,
      ParValue = estimate
    ) %>% dplyr::mutate(instrument = NA,
                        ppmMassRange = as.numeric(ppmMassRange)) %>%
    dplyr::select(
      "instrument",
      "class",
      "frag",
      "adduct",
      "ppmMassRange",
      "model",
      "ParKey",
      "ParValue"
    ) %>% dplyr::filter(ppmMassRange == min(ppmMassRange))

  print(lipidCreatorParams)
  readr::write_tsv(lipidCreatorParams, path = file.path(paste0(
    outputPrefix, "-lipidcreator-parameters.tsv"
  )))

  # merge parameters and CI estimates
  params <- merge(params, CI, by = dplyr::intersect(names(params), names(CI)))
  # params <- params %>% mutate(combinationId = paste(fragment, adduct, polarity, sep = "-", collapse=T) )
  readr::write_tsv(params, path = file.path(paste(outputPrefix, "-parameters.tsv", sep =
                                             "")))
  preds_from_data <- fits %>%  tidyr::unnest(fit %>% purrr::map(augment))

  # use tibble with combinationId for merging and create x-values for each fit (covering the complete x-axis ranges)
  new_preds <- nls.tibbleId %>%
    do(.,
       data.frame(
         precursorCollisionEnergy = seq(min(0), max(.$precursorCollisionEnergy), by =
                                          1),
         stringsAsFactors = FALSE
       ))

  # calculate the minimum and maximum range of x-values for each unique combination (fit) and add some slack
  max_min <- dplyr::group_by(nls.tibbleId, combinationId) %>%
    dplyr::summarise(
      .,
      min_pce = plyr::round_any(min(precursorCollisionEnergy), 10, f = floor) - 1,
      max_pce = plyr::round_any(max(precursorCollisionEnergy), 10, f = ceiling) + 1
    ) %>%
    dplyr::ungroup()

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

  # get details of fits
  fitinfo <-
    info %>% dplyr::select(combinationId, fitFun, sigma, isConv, finTol, logLik, AIC, BIC)
  readr::write_tsv(preds, path = file.path(paste(outputPrefix, "-predictions.tsv", sep =
                                            "")))
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

plotFits <-
  function(nlsFitPlotsOutputList,
           outputPrefix,
           format = "pdf") {
    #Residualplot
    preds_from_data <- nlsFitPlotsOutputList$preds_from_data
    preds_from_data <-
      preds_from_data %>% tidyr::separate(
        combinationId,
        c(
          "species",
          "fragment",
          "adduct",
          "polarity",
          "foundMassRange[ppm]"
        ),
        sep = "\\|",
        remove = FALSE
      )
    resplot <- ggplot2::ggplot() +
      ggplot2::geom_point(ggplot2::aes(precursorCollisionEnergy, .resid, colour = fragment),
                 size = 2,
                 preds_from_data) +
      ggplot2::geom_rug(
        ggplot2::aes(precursorCollisionEnergy, .resid, colour = fragment),
        alpha = 0.5,
        sides = "b",
        preds_from_data
      ) +
      ggplot2::geom_hline(yintercept = 0,
                 col = "red",
                 linetype = "dashed") +
      ggplot2::facet_wrap(
        fragment + adduct ~ polarity,
        labeller = ggplot2::labeller(.multi_line = TRUE),
        ncol = 6
      ) +
      ggplot2::theme_bw(base_size = 12, base_family = 'Helvetica') +
      ggplot2::labs(title = paste0(unique(preds_from_data$species)), colour = 'Fragment') +
      ggplot2::ylab(expression(paste("Residuals (", Delta,"(",y,",",yhat,")", ")",sep = " "))) +
      ggplot2::xlab('HCD Collision Energy [eV]')# +
    #theme(legend.position = c(0.9, 0.15))
    ggplot2::ggsave(
      resplot,
      filename = paste0(outputPrefix, "-residuals.", format),
      width = a4r$width,
      height = a4r$height
    )

    preds <- nlsFitPlotsOutputList$preds
    preds <-
      preds %>% tidyr::separate(
        combinationId,
        c(
          "species",
          "fragment",
          "adduct",
          "polarity",
          "foundMassRange[ppm]"
        ),
        sep = "\\|",
        remove = FALSE
      )

    params <- nlsFitPlotsOutputList$params

    CI <- nlsFitPlotsOutputList$CI
    nls.tibble <- nlsFitPlotsOutputList$nls.tibble
    nls.tibble.mean <-
      nls.tibble %>% dplyr::group_by(fragment,
                              adduct,
                              polarity,
                              `foundMassRange[ppm]`,
                              precursorCollisionEnergy) %>% dplyr::summarize(m = mean(scanRelativeIntensity))
    print(nls.tibble.mean)

    fitplot <- ggplot2::ggplot() +
      ggplot2::geom_point(
        ggplot2::aes(precursorCollisionEnergy, scanRelativeIntensity, colour = fragment),
        size = 2,
        nls.tibble
      ) +
      ggplot2::geom_rug(
        ggplot2::aes(precursorCollisionEnergy, scanRelativeIntensity, colour = fragment),
        alpha = 0.5,
        sides = "b",
        nls.tibble
      ) +
      ggplot2::geom_point(
        ggplot2::aes(precursorCollisionEnergy, m),
        data = nls.tibble.mean,
        colour = "grey80",
        size = 0.5,
        alpha = 0.75
      ) +
      #    geom_ribbon(aes(precursorCollisionEnergy, ymin = CI_low, ymax = CI_high), fill = 'lightgray', alpha = .2, preds) +
      ggplot2::geom_line(
        ggplot2::aes(precursorCollisionEnergy, scanRelativeIntensity, group = adduct),
        colour = "blue",
        alpha = 0.5,
        preds
      ) +
      ggplot2::facet_wrap(
        fragment + adduct + `foundMassRange[ppm]` ~ polarity,
        labeller = ggplot2::labeller(.multi_line = TRUE),
        ncol = 6
      ) +
      ggplot2::theme_bw(base_size = 12, base_family = 'Helvetica') +
      ggplot2::labs(title = paste0(unique(preds_from_data$species)), colour = 'Fragment') +
      ggplot2::ylab('Relative Intensity') +
      ggplot2::xlab('HCD Collision Energy [eV]')# +
    #theme(legend.position = c(0.9, 0.15))
    ggplot2::ggsave(
      fitplot,
      filename = paste0(outputPrefix, "-fit.", format),
      width = a4r$width,
      height = a4r$height
    )

    #separate out columns that are united in combinationId
    params <-
      params %>% tidyr::separate(
        combinationId,
        c(
          "species",
          "fragment",
          "adduct",
          "polarity",
          "foundMassRange[ppm]"
        ),
        sep = "\\|",
        remove = FALSE
      )

    ciplot <- ggplot2::ggplot(params, ggplot2::aes(col = fragment)) +
      ggplot2::geom_point(ggplot2::aes(combinationId, estimate)) +
      ggplot2::facet_wrap(adduct ~ term, scale = 'free_x', ncol = 2) +
      ggplot2::geom_errorbar(ggplot2::aes(combinationId, ymin = conf.low, ymax = conf.high)) +
      ggplot2::coord_flip() +
      ggplot2::theme_bw(base_size = 12, base_family = 'Helvetica') +
      ggplot2::labs(title = paste0(unique(preds_from_data$species)),
           caption = "Conf.Int. between [2.5%, 97.5%]",
           colour = 'Fragment') +
      ggplot2::xlab('Fragment') +
      ggplot2::ylab('Log-Normal Parameter Estimate')
    ggplot2::ggsave(
      ciplot,
      filename = paste0(outputPrefix, "-confint.", format),
      width = a4r$width,
      height = a4r$height
    )
}
