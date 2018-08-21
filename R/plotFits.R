#' @importFrom magrittr %>%
#' @export
plotParameterConfidenceIntervals <- function(nlsFitPlotsOutputList, combinationId, outputPrefix, plotFormat="pdf") {

  #separate out columns that are united in combinationId
  params <-
    nlsFitPlotsOutputList$params %>% tidyr::separate(
      combinationId,
      c(
        "species",
        "fragment",
        "adduct",
        "polarity",
        "calculatedMass",
        "foundMassRange[ppm]",
        "group"
      ),
      sep = "\\|",
      remove = FALSE
    )

  params$fragadd <- paste(params$fragment, params$adduct, sep = " ")
  params$fragadd <-
    factor(params$fragadd, levels = unique(params[order(params$calculatedMass),]$fragadd))
  params <- params %>% dplyr::arrange(calculatedMass, combinationId)
  ciplot <- ggplot2::ggplot(params, ggplot2::aes(col = fragadd)) +
    ggplot2::geom_point(ggplot2::aes(combinationId, estimate)) +
    ggplot2::facet_wrap(adduct ~ term, scale = 'free_x',labeller = ggplot2::label_wrap_gen(multi_line=FALSE), ncol = 2) +
    ggplot2::geom_errorbar(ggplot2::aes(combinationId, ymin = conf.low, ymax = conf.high)) +
    ggplot2::coord_flip() +
    ggplot2::theme_bw(base_size = 12, base_family = 'Helvetica') +
    ggplot2::labs(title = paste0(unique(params$species)," ",unique(params$group)),
                  caption = "Conf.Int. between [2.5%, 97.5%]",
                  colour = 'Fragment') +
    ggplot2::xlab('Fragment') +
    ggplot2::ylab('Log-Normal Parameter Estimate')
  ggplot2::ggsave(
    ciplot,
    filename = paste0(outputPrefix, "-confint.", plotFormat),
    width = a4r$width,
    height = a4r$height
  )
}

#' @importFrom magrittr %>%
#' @export
plotPredictedFits <- function(nlsFitPlotsOutputList, outputPrefix, plotFormat="pdf") {

  #predictions plot
  preds <-
    nlsFitPlotsOutputList$preds %>% tidyr::separate(
      combinationId,
      c(
        "species",
        "fragment",
        "adduct",
        "polarity",
        "calculatedMass",
        "foundMassRange[ppm]",
        "group"
      ),
      sep = "\\|",
      remove = FALSE
    )

  preds$fragadd <- paste(preds$fragment, preds$adduct, sep = " ")
  preds$fragadd <-
    factor(preds$fragadd, levels = unique(preds[order(preds$calculatedMass),]$fragadd))

  nls.tibble <- nlsFitPlotsOutputList$nls.tibble
  nls.tibble$fragadd <- paste(nls.tibble$fragment, nls.tibble$adduct, sep = " ")
  nls.tibble$fragadd <-
    factor(nls.tibble$fragadd, levels = unique(nls.tibble[order(nls.tibble$calculatedMass),]$fragadd))
  nls.tibble.mean <-
    nls.tibble %>% dplyr::group_by(fragment,
                                   adduct,
                                   polarity,
                                   `foundMassRange[ppm]`,
                                   precursorCollisionEnergy, group) %>% dplyr::summarize(m = mean(scanRelativeIntensity))
  nls.tibble.mean$fragadd <- paste(nls.tibble.mean$fragment, nls.tibble.mean$adduct, sep = " ")
  nls.tibble.mean$fragadd <-
    factor(nls.tibble.mean$fragadd, levels = levels(nls.tibble$fragadd))
  fitplot <- ggplot2::ggplot() +
    ggplot2::geom_point(
      ggplot2::aes(precursorCollisionEnergy, scanRelativeIntensity, colour = fragadd),
      size = 2,
      nls.tibble
    ) +
    ggplot2::geom_rug(
      ggplot2::aes(precursorCollisionEnergy, scanRelativeIntensity, colour = fragadd),
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
      fragadd + `foundMassRange[ppm]` ~ polarity,
      labeller = ggplot2::label_wrap_gen(multi_line=FALSE),
      ncol = 6
    ) +
    ggplot2::theme_bw(base_size = 12, base_family = 'Helvetica') +
    ggplot2::labs(title = paste0(unique(preds$species)," ",unique(preds$group)), colour = 'Fragment') +
    ggplot2::ylab('Relative Intensity') +
    ggplot2::xlab('Collision Energy [eV]')
  #theme(legend.position = c(0.9, 0.15))
  ggplot2::ggsave(
    fitplot,
    filename = paste0(outputPrefix, "-fit.", plotFormat),
    width = a4r$width,
    height = a4r$height
  )
}

#' @importFrom magrittr %>%
#' @export
plotResiduals <- function(nlsFitPlotsOutputList, outputPrefix, plotFormat = "pdf") {
  #Residualplot
  preds_from_data <-
    nlsFitPlotsOutputList$preds_from_data %>% tidyr::separate(
      combinationId,
      c(
        "species",
        "fragment",
        "adduct",
        "polarity",
        "calculatedMass",
        "foundMassRange[ppm]",
        "group"
      ),
      sep = "\\|",
      remove = FALSE
    )
  preds_from_data$fragadd <- paste(preds_from_data$fragment, preds_from_data$adduct, sep = " ")
  preds_from_data$fragadd <-
    factor(preds_from_data$fragadd, levels = unique(preds_from_data[order(preds_from_data$calculatedMass),]$fragadd))
  resplot <- ggplot2::ggplot(preds_from_data, ggplot2::aes(col = fragadd)) +
    ggplot2::geom_point(ggplot2::aes(precursorCollisionEnergy, .resid),
                        size = 2,
                        preds_from_data) +
    ggplot2::geom_rug(
      ggplot2::aes(precursorCollisionEnergy, .resid),
      #colour = preds_from_data$fragadd,
      alpha = 0.5,
      sides = "b",
      preds_from_data
    ) +
    ggplot2::geom_hline(yintercept = 0,
                        col = "red",
                        linetype = "dashed") +
    ggplot2::facet_wrap(
      fragadd ~ `foundMassRange[ppm]`,
      labeller = ggplot2::label_wrap_gen(multi_line=FALSE),
      ncol = 6
    ) +
    ggplot2::theme_bw(base_size = 12, base_family = 'Helvetica') +
    ggplot2::labs(title = paste0(unique(preds_from_data$species)," ",unique(preds_from_data$group)), colour = 'Fragment') +
    ggplot2::ylab(expression(paste("Residuals (", Delta,"(",y,",",yhat,")", ")",sep = " "))) +
    ggplot2::xlab('Collision Energy [eV]') +
    ggplot2::scale_y_continuous(limits = c(-0.1, 0.1))
  # +
  #theme(legend.position = c(0.9, 0.15))
  ggplot2::ggsave(
    resplot,
    filename = paste0(outputPrefix, "-residuals.", plotFormat),
    width = a4r$width,
    height = a4r$height
  )
}

#' @importFrom magrittr %>%
#' @export
plotFits <-
  function(nlsFitPlotsOutputList,
           outputPrefix,
           plotFormat = "pdf") {
    plotResiduals(nlsFitPlotsOutputList=nlsFitPlotsOutputList, outputPrefix=outputPrefix, plotFormat=plotFormat)
    plotPredictedFits(nlsFitPlotsOutputList=nlsFitPlotsOutputList, outputPrefix=outputPrefix, plotFormat=plotFormat)
    plotParameterConfidenceIntervals(nlsFitPlotsOutputList=nlsFitPlotsOutputList, combinationId=combinationId, outputPrefix=outputPrefix, plotFormat=plotFormat)
  }