#' @export
plotRawTicVsTotalIonCurrent <- function(data, basename, plotFormat="png", plotDimensions=list(width=11.69, height=8.27)) {
  rawTicVsTotalIonCurrentPlot <-
    ggplot2::qplot(
      x = rawTic,
      y = totalIonCurrent,
      data = data,
      color = fragadd,
      geom = "point",
      main = paste(unique(data$species), unique(data$polarity), sep = " ")
    )  + ggplot2::labs(colour = 'Fragment') +
    ggplot2::facet_wrap(fragadd ~ data$`foundMassRange[ppm]`, ncol =6, labeller = ggplot2::label_context(multi_line=FALSE)) +
    ggplot2::ylab("Total Ion Current [a.u.]") + ggplot2::xlab("Raw Total Ion Current [a.u.]")
  ggplot2::ggsave(
    rawTicVsTotalIonCurrentPlot,
    filename = paste0(basename, "-rawTic-vs-Tic.", plotFormat),
    width = plotDimensions$width,
    height = plotDimensions$height
  )
}
#' @export
plotFragmentPpmBoxplot <- function(data, basename, plotFormat="png", plotDimensions=list(width=11.69, height=8.27)) {
  fpbPlot <-
    ggplot2::qplot(
      x = fragadd,
      y = (calculatedMass - foundMass) ,
      data = data,
      color = fragadd,
      geom = c("violin"),
      main = paste(unique(data$species), unique(data$polarity), sep = " ")
    ) + ggplot2::labs(colour = 'Fragment') +
    ggplot2::facet_wrap( ~ data$`foundMassRange[ppm]`, ncol =
                           6, labeller = ggplot2::label_context(multi_line=FALSE)) + ggplot2::ylab(expression(paste(
                             Delta, " of calculated and found m/z", sep = " "
                           ))) + ggplot2::xlab("Fragment") + ggplot2::geom_jitter(ggplot2::aes(
                             x = fragadd,
                             y = (calculatedMass - foundMass),
                             color = fragadd
                           ),
                           data = data,
                           alpha = 0.1)
  ggplot2::ggsave(
    fpbPlot,
    filename = paste0(basename, "-frag-ppm-box.", plotFormat),
    width = plotDimensions$width,
    height = plotDimensions$height
  )
  fpbPlot
}
#' @export
plotPrecCollEnergyVsFoundIntensity <- function(data, basename, plotFormat="png", plotDimensions=list(width=11.69, height=8.27)) {
  message("precursorCollisionEnergy-vs-foundIntensity")

  pceFiPlot <-
    ggplot2::qplot(
      x = precursorCollisionEnergy,
      y = foundIntensity,
      data = data,
      color = fragadd,
      geom = c("point"),
      main = paste(
        unique(data$species),
        unique(data$polarity),
        sep = " "
      )
    ) + ggplot2::geom_smooth(
      method = "loess",
      colour = "blue",
      se = TRUE,
      level = 0.95
    ) + ggplot2::scale_x_continuous(breaks = seq(from=-10, to=plyr::round_any(max(data$precursorCollisionEnergy), 10, f = ceiling)+10, by=10)) +
    ggplot2::ylab("Absolute Intensity [a.u.]") +
    ggplot2::xlab("Precursor Collision Energy [eV]") +
    ggplot2::labs(colour = 'Fragment') +
    ggplot2::facet_wrap(fragadd ~ data$`foundMassRange[ppm]`, ncol = 6, labeller = ggplot2::label_context(multi_line=FALSE))
  ggplot2::ggsave(
    pceFiPlot,
    filename = paste0(basename, "-precCE-vs-I.", plotFormat),
    width = plotDimensions$width,
    height = plotDimensions$height
  )
}
#' @export
plotPrecCollEnergyVsScanRelativeIntensityNormalized <- function(data, basename, plotFormat="png", plotDimensions=list(width=11.69, height=8.27)) {
  message("precursorCollisionEnergy-vs-foundIntensity-scan-relative-normalized")
  precCeVsIsrnPlot <-
    ggplot2::qplot(
      x = precursorCollisionEnergy,
      y = scanRelativeIntensity,
      data = data,
      color = fragadd,
      geom = c("point"),
      main = paste(
        unique(data$species),
        unique(data$polarity),
        sep = " "
      )
    ) + ggplot2::geom_smooth(
      method = "loess",
      colour = "blue",
      se = TRUE,
      level = 0.95
    ) + ggplot2::scale_x_continuous(breaks = seq(from=-10, to=plyr::round_any(max(data$precursorCollisionEnergy), 10, f = ceiling)+10, by=10)) +
    ggplot2::ylab("Scan Relative Intensity") +
    ggplot2::xlab("Precursor Collision Energy [eV]") +
    ggplot2::labs(colour = 'Fragment') +
    ggplot2::facet_wrap(fragadd ~ data$`foundMassRange[ppm]`, ncol = 6, labeller = ggplot2::label_context(multi_line=FALSE))
  ggplot2::ggsave(
    precCeVsIsrnPlot,
    filename = paste0(basename, "-precCE-vs-I-srn.", plotFormat),
    width = plotDimensions$width,
    height = plotDimensions$height
  )
  precCeVsIsrnPlot
}
#' @export
plotPrecCollEnergyVsScanRelativeIntensityOverlay <- function(data, basename, plotFormat="png", plotDimensions=list(width=11.69, height=8.27)) {

  message(
    "precursorCollisionEnergy-vs-foundIntensity-scan-relative-normalized-overlay"
  )
  precCeVsIsrnOverlayPlot <-
    ggplot2::qplot(
      x = precursorCollisionEnergy,
      y = scanRelativeIntensity,
      data = data,
      color = fragadd,
      fill = fragadd,
      # geom = c("boxplot"),
      main = paste(
        unique(data$species),
        unique(data$polarity),
        sep = " "
      )
    ) +
    ggplot2::geom_area(position = "identity", alpha = 0.3) +
    ggplot2::scale_x_continuous(breaks = seq(from=-10, to=plyr::round_any(max(data$precursorCollisionEnergy), 10, f = ceiling)+10, by=10)) +
    ggplot2::ylab("Scan Relative Intensity") +
    ggplot2::xlab("Precursor Collision Energy [eV]") +
    ggplot2::labs(colour = 'Fragment', fill = 'Fragment') +
    ggplot2::facet_wrap(~ data$`foundMassRange[ppm]`, ncol = 2, labeller = ggplot2::label_context(multi_line=FALSE))
  ggplot2::ggsave(
    precCeVsIsrnOverlayPlot,
    filename = paste0(basename, "-precCE-vs-I-srn-overlay.", plotFormat),
    width = plotDimensions$width,
    height = plotDimensions$height
  )
  precCeVsIsrnOverlayPlot
}
#' @export
plotPrecCollEnergyVsMassErrorPpm <- function(data, basename, plotFormat="png", plotDimensions=list(width=11.69, height=8.27)) {
  message("precursorCollisionEnergy-vs-mass-error-ppm")
  precCeVsMerrPpmmPlot <-
    ggplot2::qplot(
      x = precursorCollisionEnergy,
      y = `foundMassError[ppm]`,
      data = data,
      color = fragadd,
      geom = c("point"),
      main = paste(
        unique(data$species),
        unique(data$polarity),
        sep = " "
      )
    ) + ggplot2::geom_smooth(
      method = "loess",
      colour = "blue",
      se = TRUE,
      level = 0.95
    ) + ggplot2::scale_x_continuous(breaks = seq(from=-10, to=plyr::round_any(max(data$precursorCollisionEnergy), 10, f = ceiling)+10, by=10)) +
    ggplot2::ylab("Mass Error [ppm]") +
    ggplot2::xlab("Precursor Collision Energy [eV]") +
    ggplot2::labs(colour = 'Fragment') +
    ggplot2::facet_wrap(fragadd ~ data$`foundMassRange[ppm]`, ncol = 6, labeller = ggplot2::label_context(multi_line=FALSE))
  ggplot2::ggsave(
    precCeVsMerrPpmmPlot,
    filename = paste0(basename, "-precCE-vs-merr-ppm.", plotFormat),
    width = plotDimensions$width,
    height = plotDimensions$height
  )
  precCeVsMerrPpmmPlot
}
#' @export
plotMassDensityDistribution <- function(data, basename, plotFormat="png", nppmLevels, plotDimensions=list(width=11.69, height=8.27)) {
  message("m/z density distribution")
  stopifnot(!missing(data))
  stopifnot(!missing(basename))
  nppms <- length(unique(data$`foundMassRange[ppm]`))
  #plot m/z density distributions
  mddPlot <-
    ggplot2::ggplot(ggplot2::aes(
      x = foundMass,
      color = fragadd,
      fill = fragadd
    ), data = data) +
    ggplot2::ggtitle(paste(
      unique(data$species),
      unique(data$polarity),
      sep = " "
    )) + ggplot2::geom_density(ggplot2::aes(y = ..count.. / sum(..count..)), data = data) +
    ggplot2::geom_vline(ggplot2::aes(xintercept = calculatedMass),
                        data = data,
                        linetype = 3) +
    ggplot2::geom_text(
      x = data$calculatedMass,
      y = 0,
      label = round(data$calculatedMass, digits = 4),
      check_overlap = TRUE,
      colour = "blue",
      hjust = "left",
      vjust = "top",
      size = 3
    ) +
    ggplot2::geom_vline(
      ggplot2::aes(xintercept = foundMassLowerBound),
      colour = "darkblue",
      data = data,
      linetype = 2
    ) +
    ggplot2::geom_vline(
      ggplot2::aes(xintercept = foundMassUpperBound),
      colour = "darkblue",
      data = data,
      linetype = 2
    ) +
    ggplot2::facet_wrap(~ `foundMassRange[ppm]`, nrow = nppms, labeller = ggplot2::label_context(multi_line=FALSE)) +
    ggplot2::xlab("Fragment m/z") +
    ggplot2::ylab("Normalized Count") +
    ggplot2::labs(colour = 'Fragment', fill = 'Fragment')
  ggplot2::ggsave(
    mddPlot,
    filename = paste0(basename, "-mass-dens-distr.", plotFormat),
    width = plotDimensions$width,
    height = plotDimensions$height
  )
  mddPlot
}
#' @export
plotMzVsMerrPpm <- function(data, basename, plotFormat="png", plotDimensions=list(width=11.69, height=8.27)) {
  # plot m/z vs mass error
  message("m/z-vs-mass-error-ppm")
  mzVsMerrPpmPlot <- ggplot2::qplot(
    x = foundMass,
    y = `foundMassError[ppm]`,
    data = data,
    color = fragadd,
    geom = c("violin"),
    main = paste(
      unique(data$species),
      unique(data$polarity),
      sep = " "
    )
  ) + ggplot2::ylab("Mass Error [ppm]") +
    ggplot2::xlab("Fragment m/z") +
    ggplot2::labs(colour = 'Fragment') +
    ggplot2::facet_wrap(~ data$`foundMassRange[ppm]`, ncol = 6, labeller = ggplot2::label_context(multi_line=FALSE)) +
    ggplot2::geom_jitter(
      ggplot2::aes(x = foundMass, y = `foundMassError[ppm]`, color = fragadd),
      data = data,
      alpha = 0.1
    )
  ggplot2::ggsave(
    mzVsMerrPpmPlot,
    filename = paste0(basename, "-mz-vs-merr-ppm.", plotFormat),
    width = plotDimensions$width,
    height = plotDimensions$height
  )
  plotMzVsMerrPpm
}