#' Plots the raw, non-normalized TIC values.
#' @param data the FIP data to plot.
#' @param basename the basename (identifying a case) for the plot.
#' @param plotFormat the format, passed to \code{ggplot2::ggsave}.
#' @param plotDimensions the dimensions of the plot.
#' @param color_scale the shared color scale to identify fragment and adduct pairs.
#' @return the ggplot object.
#' @export
plotRawTicVsTotalIonCurrent <- function(data, basename, plotFormat="png", plotDimensions=list(width=11.69, height=8.27), color_scale = ggplot2::scale_colour_hue()) {
  rawTicVsTotalIonCurrentPlot <-
    ggplot2::qplot(
      x = rawTic,
      y = totalIonCurrent,
      data = data,
      color = fragadd,
      geom = "point",
      main = paste(unique(data$species), unique(data$polarity), unique(data$group), sep = " ")
    )  + ggplot2::labs(colour = 'Fragment') +
    ggplot2::facet_wrap(fragadd ~ data$`foundMassRange[ppm]`, ncol =6, labeller = ggplot2::label_wrap_gen(multi_line=FALSE)) +
    ggplot2::theme_bw(base_size = 12, base_family = 'Helvetica') +
    ggplot2::theme(axis.text.x = ggplot2::element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1
    )) +
    ggplot2::ylab("Total Ion Current [a.u.]") + ggplot2::xlab("Raw Total Ion Current [a.u.]") +
    color_scale
  ggplot2::ggsave(
    rawTicVsTotalIonCurrentPlot,
    filename = paste0(basename, "-rawTic-vs-Tic.", plotFormat),
    width = plotDimensions$width,
    height = plotDimensions$height
  )
  return(rawTicVsTotalIonCurrentPlot)
}

#' Plots a PPM deviation boxplot for all fragments.
#' @param data the FIP data to plot.
#' @param basename the basename (identifying a case) for the plot.
#' @param plotFormat the format, passed to \code{ggplot2::ggsave}.
#' @param plotDimensions the dimensions of the plot.
#' @param color_scale the shared color scale to identify fragment and adduct pairs.
#' @return the ggplot object.
#' @export
plotFragmentPpmBoxplot <- function(data, basename, plotFormat="png", plotDimensions=list(width=11.69, height=8.27), color_scale = ggplot2::scale_colour_hue()) {
  fpbPlot <-
    ggplot2::qplot(
      x = fragadd,
      y = (calculatedMass - foundMass) ,
      data = data,
      color = fragadd,
      geom = c("violin"),
      main = paste(unique(data$species), unique(data$polarity), unique(data$group), sep = " ")
    ) + ggplot2::labs(colour = 'Fragment') +
    ggplot2::theme_bw(base_size = 12, base_family = 'Helvetica') +
    ggplot2::facet_wrap( ~ data$`foundMassRange[ppm]`, ncol =
                           6, labeller = ggplot2::label_wrap_gen(multi_line=FALSE)) + ggplot2::ylab(expression(paste(
                             Delta, " of calculated and found m/z", sep = " "
                           ))) + ggplot2::xlab("Fragment") + ggplot2::geom_jitter(ggplot2::aes(
                             x = fragadd,
                             y = (calculatedMass - foundMass),
                             color = fragadd
                           ),
                           data = data,
                           alpha = 0.1) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
    color_scale
  ggplot2::ggsave(
    fpbPlot,
    filename = paste0(basename, "-frag-ppm-box.", plotFormat),
    width = plotDimensions$width,
    height = plotDimensions$height
  )
  return(fpbPlot)
}

#' Plots precursor collision energy vs found raw intensity.
#' @param data the FIP data to plot.
#' @param basename the basename (identifying a case) for the plot.
#' @param plotFormat the format, passed to \code{ggplot2::ggsave}.
#' @param plotDimensions the dimensions of the plot.
#' @param color_scale the shared color scale to identify fragment and adduct pairs.
#' @return the ggplot object.
#' @export
plotPrecCollEnergyVsFoundIntensity <- function(data, basename, plotFormat="png", plotDimensions=list(width=11.69, height=8.27), color_scale = ggplot2::scale_colour_hue()) {
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
        unique(data$group),
        sep = " "
      )
    ) + ggplot2::geom_smooth(
      method = "loess",
      colour = "blue",
      se = TRUE,
      level = 0.95
    ) + ggplot2::scale_x_continuous(breaks = seq(from=-10, to=plyr::round_any(max(data$precursorCollisionEnergy), 10, f = ceiling)+10, by=10)) +
    ggplot2::theme_bw(base_size = 12, base_family = 'Helvetica') +
    ggplot2::ylab("Absolute Intensity [a.u.]") +
    ggplot2::xlab(flipr::collisionEnergyLabel(data)) +
    ggplot2::labs(colour = 'Fragment') +
    ggplot2::facet_wrap(fragadd ~ data$`foundMassRange[ppm]`, ncol = 6, labeller = ggplot2::label_wrap_gen(multi_line=FALSE)) +
    color_scale
  ggplot2::ggsave(
    pceFiPlot,
    filename = paste0(basename, "-precCE-vs-I.", plotFormat),
    width = plotDimensions$width,
    height = plotDimensions$height
  )
  return(pceFiPlot)
}

#' Plots precursor collision energy vs. total scan intensity normalized relative fragment intensity in separate panels.
#' @param data the FIP data to plot.
#' @param basename the basename (identifying a case) for the plot.
#' @param plotFormat the format, passed to \code{ggplot2::ggsave}.
#' @param plotDimensions the dimensions of the plot.
#' @param color_scale the shared color scale to identify fragment and adduct pairs.
#' @return the ggplot object.
#' @export
plotPrecCollEnergyVsScanRelativeIntensityNormalized <- function(data, basename, plotFormat="png", plotDimensions=list(width=11.69, height=8.27), color_scale = ggplot2::scale_colour_hue()) {
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
        unique(data$group),
        sep = " "
      )
    ) + ggplot2::geom_point(
      ggplot2::aes(x = precursorCollisionEnergy, y = foundIntensity, color="darkgray", alpha=0.7), show.legend = FALSE
    ) + ggplot2::geom_smooth(
      method = "loess",
      colour = "blue",
      se = TRUE,
      level = 0.95
    ) + ggplot2::scale_x_continuous(breaks = seq(from=-10, to=plyr::round_any(max(data$precursorCollisionEnergy), 10, f = ceiling)+10, by=10)) +
    ggplot2::theme_bw(base_size = 12, base_family = 'Helvetica') +
    ggplot2::ylab("Scan Relative Intensity") +
    ggplot2::xlab(flipr::collisionEnergyLabel(data)) +
    ggplot2::labs(colour = 'Fragment') +
    ggplot2::facet_wrap(fragadd ~ data$`foundMassRange[ppm]`, ncol = 6, labeller = ggplot2::label_wrap_gen(multi_line=FALSE)) + ggplot2::ylim(0,1) +
    color_scale
  ggplot2::ggsave(
    precCeVsIsrnPlot,
    filename = paste0(basename, "-precCE-vs-I-srn.", plotFormat),
    width = plotDimensions$width,
    height = plotDimensions$height
  )
  return(precCeVsIsrnPlot)
}

#' Plots precursor collision energy vs. total scan intensity normalized relative fragment intensity overlaid on one panel.
#' @param data the FIP data to plot.
#' @param basename the basename (identifying a case) for the plot.
#' @param plotFormat the format, passed to \code{ggplot2::ggsave}.
#' @param plotDimensions the dimensions of the plot.
#' @param color_scale the shared color scale to identify fragment and adduct pairs.
#' @return the ggplot object.
#' @export
plotPrecCollEnergyVsScanRelativeIntensityOverlay <- function(data, basename, plotFormat="png", plotDimensions=list(width=11.69, height=8.27), color_scale = ggplot2::scale_colour_hue()) {

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
        unique(data$group),
        sep = " "
      )
    ) +
    ggplot2::geom_area(position = "identity", alpha = 0.3) +
    ggplot2::scale_x_continuous(breaks = seq(from=-10, to=plyr::round_any(max(data$precursorCollisionEnergy), 10, f = ceiling)+10, by=10)) +
    ggplot2::ylab("Scan Relative Intensity") +
    ggplot2::xlab(flipr::collisionEnergyLabel(data)) +
    ggplot2::labs(colour = 'Fragment', fill = 'Fragment') +
    ggplot2::facet_wrap(~ data$`foundMassRange[ppm]`, ncol = 2, labeller = ggplot2::label_wrap_gen(multi_line=FALSE)) + ggplot2::ylim(0,1) +
    ggplot2::theme_bw(base_size = 12, base_family = 'Helvetica') +
    color_scale
  ggplot2::ggsave(
    precCeVsIsrnOverlayPlot,
    filename = paste0(basename, "-precCE-vs-I-srn-overlay.", plotFormat),
    width = plotDimensions$width,
    height = plotDimensions$height
  )
  return(precCeVsIsrnOverlayPlot)
}

#' Plots the precursor collision energy vs. the fragment mass error in ppm.
#' @param data the FIP data to plot.
#' @param basename the basename (identifying a case) for the plot.
#' @param plotFormat the format, passed to \code{ggplot2::ggsave}.
#' @param plotDimensions the dimensions of the plot.
#' @param color_scale the shared color scale to identify fragment and adduct pairs.
#' @return the ggplot object.
#' @export
plotPrecCollEnergyVsMassErrorPpm <- function(data, basename, plotFormat="png", plotDimensions=list(width=11.69, height=8.27), color_scale = ggplot2::scale_colour_hue()) {
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
        unique(data$group),
        sep = " "
      )
    ) + ggplot2::geom_smooth(
      method = "loess",
      colour = "blue",
      se = TRUE,
      level = 0.95
    ) + ggplot2::scale_x_continuous(breaks = seq(from=-10, to=plyr::round_any(max(data$precursorCollisionEnergy), 10, f = ceiling)+10, by=10)) +
    ggplot2::ylab("Mass Error [ppm]") +
    ggplot2::xlab(flipr::collisionEnergyLabel(data)) +
    ggplot2::labs(colour = 'Fragment') +
    ggplot2::facet_wrap(fragadd ~ data$`foundMassRange[ppm]`, ncol = 6, labeller = ggplot2::label_wrap_gen(multi_line=FALSE)) +
    ggplot2::theme_bw(base_size = 12, base_family = 'Helvetica') +
    color_scale
  ggplot2::ggsave(
    precCeVsMerrPpmmPlot,
    filename = paste0(basename, "-precCE-vs-merr-ppm.", plotFormat),
    width = plotDimensions$width,
    height = plotDimensions$height
  )
  return(precCeVsMerrPpmmPlot)
}

#' Plots the mass density distribution.
#' @param data the FIP data to plot.
#' @param basename the basename (identifying a case) for the plot.
#' @param plotFormat the format, passed to \code{ggplot2::ggsave}.
#' @param plotDimensions the dimensions of the plot.
#' @param color_scale the shared color scale to identify fragment and adduct pairs.
#' @return the ggplot object.
#' @export
plotMassDensityDistribution <- function(data, basename, plotFormat="png", plotDimensions=list(width=11.69, height=8.27), color_scale = ggplot2::scale_colour_hue()) {
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
      unique(data$group),
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
    ggplot2::facet_wrap(~ `foundMassRange[ppm]`, nrow = nppms, labeller = ggplot2::label_wrap_gen(multi_line=FALSE)) +
    ggplot2::theme_bw(base_size = 12, base_family = 'Helvetica') +
    ggplot2::xlab("Fragment m/z") +
    ggplot2::ylab("Normalized Count") +
    ggplot2::labs(colour = 'Fragment', fill = 'Fragment') +
    color_scale
  ggplot2::ggsave(
    mddPlot,
    filename = paste0(basename, "-mass-dens-distr.", plotFormat),
    width = plotDimensions$width,
    height = plotDimensions$height
  )
  return(mddPlot)
}
#' Plots m/z vs. the mass error in ppm.
#' @param data the FIP data to plot.
#' @param basename the basename (identifying a case) for the plot.
#' @param plotFormat the format, passed to \code{ggplot2::ggsave}.
#' @param plotDimensions the dimensions of the plot.
#' @param color_scale the shared color scale to identify fragment and adduct pairs.
#' @return the ggplot object.
#' @export
plotMzVsMerrPpm <- function(data, basename, plotFormat="png", plotDimensions=list(width=11.69, height=8.27), color_scale = ggplot2::scale_colour_hue()) {
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
      unique(data$group),
      sep = " "
    )
  ) + ggplot2::ylab("Mass Error [ppm]") +
    ggplot2::xlab("Fragment m/z") +
    ggplot2::labs(colour = 'Fragment') +
    ggplot2::facet_wrap(~ data$`foundMassRange[ppm]`, ncol = 6, labeller = ggplot2::label_wrap_gen(multi_line=FALSE)) +
    ggplot2::theme_bw(base_size = 12, base_family = 'Helvetica') +
    ggplot2::geom_jitter(
      ggplot2::aes(x = foundMass, y = `foundMassError[ppm]`, color = fragadd),
      data = data,
      alpha = 0.1
    ) + color_scale
  ggplot2::ggsave(
    mzVsMerrPpmPlot,
    filename = paste0(basename, "-mz-vs-merr-ppm.", plotFormat),
    width = plotDimensions$width,
    height = plotDimensions$height
  )
  return(plotMzVsMerrPpm)
}
#' Plots the scan relative intensity histogram.
#' @param data the FIP data to plot.
#' @param basename the basename (identifying a case) for the plot.
#' @param plotFormat the format, passed to \code{ggplot2::ggsave}.
#' @param plotDimensions the dimensions of the plot.
#' @param color_scale the shared color scale to identify fragment and adduct pairs.
#' @return the ggplot object.
#' @export
plotScanRelativeIntensityHistogram <- function(data, basename, plotFormat="png", plotDimensions=list(width=11.69, height=8.27), color_scale = ggplot2::scale_colour_hue()) {

  message(
    "foundIntensity-scan-relative-normalized-histogram"
  )
  plot <-
    ggplot2::ggplot(ggplot2::aes(
      x = scanRelativeIntensity,
      color = fragadd,
      fill = fragadd
    ), data = data) +
    ggplot2::ggtitle(paste(
      unique(data$species),
      unique(data$polarity),
      unique(data$group),
      sep = " "
    )) +
    ggplot2::geom_histogram(ggplot2::aes(y =..ndensity..), alpha=0.1, binwidth = 0.01, show.legend = FALSE) +
    ggplot2::geom_density(ggplot2::aes(y =..scaled..), alpha=0.05, color="darkgray", fill="darkgray", show.legend = FALSE) +
    ggplot2::geom_rug(ggplot2::aes(x = scanRelativeIntensity, y = 0, color=fragadd), position = ggplot2::position_jitter(height = 0)) +
    #ggplot2::scale_x_continuous(breaks = seq(from=-10, to=plyr::round_any(max(data$precursorCollisionEnergy), 10, f = ceiling)+10, by=10)) +
    ggplot2::xlab("Scan Relative Intensity") +
    ggplot2::ylab("Scaled Density") +
    ggplot2::labs(color = 'Fragment', fill = 'Fragment') +
    color_scale +
    ggplot2::xlim(0,1) +
    ggplot2::facet_wrap(fragadd ~ data$`foundMassRange[ppm]`, ncol = 2, labeller = ggplot2::label_wrap_gen(multi_line=FALSE)) +
    ggplot2::theme_bw(base_size = 12, base_family = 'Helvetica')
  ggplot2::ggsave(
    plot,
    filename = paste0(basename, "-I-srn-histo.", plotFormat),
    width = plotDimensions$width,
    height = plotDimensions$height
  )
  return(plot)
}