#' Plots the number of data samples for each combinationId.
#' @param nlsFitPlotsOutputList the FIP fits output list.
#' @param outputPrefix the basename (identifying a case) for the plot.
#' @param plotFormat the format, passed to \code{ggplot2::ggsave}.
#' @param plotDimensions the dimensions of the plot.
#' @param color_scale the shared color scale to identify fragment and adduct pairs.
#' @return the ggplot object.
#' @importFrom magrittr %>%
#' @export
plotNumberOfSamplesPerCombinationId <-
  function(nlsFitPlotsOutputList,
           outputPrefix,
           plotFormat = "png",
           plotDimensions = list(width = 11.69, height = 8.27),
           color_scale = ggplot2::scale_colour_hue()) {
    data <- nlsFitPlotsOutputList$nls.tibble.unfiltered
    data$fragadd <- paste(data$fragment, data$adduct, sep = " ")
    data$fragadd <-
      factor(data$fragadd, levels = unique(data[order(data$calculatedMass), ]$fragadd))
    plot <- ggplot2::ggplot(data) +
      ggplot2::geom_boxplot(ggplot2::aes(x = combinationId, y = samplesPerCombinationId, colour =
                                           combinationId)) +
      ggplot2::theme_bw(base_size = 12, base_family = 'Helvetica') +
      ggplot2::labs(
        title = paste0(unique(data$species), " ", unique(data$group)),
        caption = "Number of samples for model training",
        colour = 'Fragment'
      ) +
      ggplot2::coord_flip() +
      ggplot2::xlab('Fragment') +
      ggplot2::ylab('Number of Samples') +
      color_scale
    ggplot2::ggsave(
      plot,
      filename = paste0(outputPrefix, "-samples-per-combinationId.", plotFormat),
      width = plotDimensions$width,
      height = plotDimensions$height
    )
    return(plot)
  }
#' Plots the confidence intervals for the parameters used for model optimization.
#' @param nlsFitPlotsOutputList the FIP fits output list.
#' @param combinationId the combinationId for the case to plot.
#' @param outputPrefix the basename (identifying a case) for the plot.
#' @param plotFormat the format, passed to \code{ggplot2::ggsave}.
#' @param plotDimensions the dimensions of the plot.
#' @param color_scale the shared color scale to identify fragment and adduct pairs.
#' @return the ggplot object.
#' @importFrom magrittr %>%
#' @export
plotParameterConfidenceIntervals <-
  function(nlsFitPlotsOutputList,
           combinationId,
           outputPrefix,
           plotFormat = "png",
           plotDimensions = list(width = 11.69, height = 8.27),
           color_scale = ggplot2::scale_colour_hue()) {
    #separate out columns that are united in combinationId
    params <-
      nlsFitPlotsOutputList$params %>% tidyr::separate(
        combinationId,
        c(
          "species",
          "precursorAdduct",
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
      factor(params$fragadd, levels = unique(params[order(params$calculatedMass), ]$fragadd))
    params <- params %>% dplyr::arrange(calculatedMass, combinationId)
    ciplot <- ggplot2::ggplot(params, ggplot2::aes(col = fragadd)) +
      ggplot2::geom_point(ggplot2::aes(combinationId, estimate)) +
      ggplot2::facet_wrap(
        adduct ~ term,
        scale = 'free_x',
        labeller = ggplot2::label_wrap_gen(multi_line = FALSE),
        ncol = 2
      ) +
      ggplot2::geom_errorbar(ggplot2::aes(combinationId, ymin = conf.low, ymax = conf.high)) +
      ggplot2::coord_trans(y = "log10") +
      ggplot2::coord_flip() +
      ggplot2::theme_bw(base_size = 12, base_family = 'Helvetica') +
      ggplot2::labs(
        title = paste0(unique(params$species), " ", unique(params$group)),
        caption = "Conf.Int. between [2.5%, 97.5%]",
        colour = 'Fragment'
      ) +
      ggplot2::xlab('Fragment') +
      ggplot2::ylab('Log-Normal Parameter Estimate') +
      color_scale
    ggplot2::ggsave(
      ciplot,
      filename = paste0(outputPrefix, "-confint.", plotFormat),
      width = plotDimensions$width,
      height = plotDimensions$height
    )
    return(ciplot)
  }

#' Plots the predicted fits.
#' @param nlsFitPlotsOutputList the FIP fits output list.
#' @param outputPrefix the basename (identifying a case) for the plot.
#' @param plotFormat the format, passed to \code{ggplot2::ggsave}.
#' @param plotDimensions the dimensions of the plot.
#' @param color_scale the shared color scale to identify fragment and adduct pairs.
#' @return the ggplot object.
#' @importFrom magrittr %>%
#' @export
plotPredictedFits <-
  function(nlsFitPlotsOutputList,
           outputPrefix,
           plotFormat = "png",
           plotDimensions = list(width = 11.69, height = 8.27),
           color_scale = ggplot2::scale_colour_hue()) {
    #predictions plot
    preds <-
      nlsFitPlotsOutputList$preds %>% tidyr::separate(
        combinationId,
        c(
          "species",
          "precursorAdduct",
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
      factor(preds$fragadd, levels = unique(preds[order(preds$calculatedMass), ]$fragadd))
    preds$`foundMassRange[ppm]` <-
      factor(preds$`foundMassRange[ppm]`,
             levels = unique(preds$`foundMassRange[ppm]`))

    nls.tibble <- nlsFitPlotsOutputList$nls.tibble
    nls.tibble$fragadd <-
      paste(nls.tibble$fragment, nls.tibble$adduct, sep = " ")
    nls.tibble$fragadd <-
      factor(nls.tibble$fragadd, levels = unique(nls.tibble[order(nls.tibble$calculatedMass), ]$fragadd))
    nls.tibble$`foundMassRange[ppm]` <-
      factor(nls.tibble$`foundMassRange[ppm]`,
             levels = sort(as.numeric(unique(nls.tibble$`foundMassRange[ppm]`))))
    nls.tibble.mean <-
      nls.tibble %>% dplyr::group_by(
        fragment,
        adduct,
        polarity,
        `foundMassRange[ppm]`,
        precursorCollisionEnergy,
        group
      ) %>% dplyr::summarize(m = mean(scanRelativeIntensity))
    nls.tibble.mean$fragadd <-
      paste(nls.tibble.mean$fragment, nls.tibble.mean$adduct, sep = " ")
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
      ggplot2::geom_line(
        ggplot2::aes(precursorCollisionEnergy, scanRelativeIntensity, group = adduct),
        colour = "blue",
        alpha = 0.5,
        preds
      ) +
      ggplot2::facet_wrap(
        fragadd + `foundMassRange[ppm]` ~ polarity,
        labeller = ggplot2::label_wrap_gen(multi_line = FALSE),
        ncol = 6
      ) +
      ggplot2::theme_bw(base_size = 12, base_family = 'Helvetica') +
      ggplot2::labs(title = paste0(unique(preds$species), " ", unique(preds$group)), colour = 'Fragment') +
      ggplot2::ylab('Relative Intensity') +
      ggplot2::xlab(flipr::collisionEnergyLabel(nlsFitPlotsOutputList$nls.tibble)) +
      color_scale
    ggplot2::ggsave(
      fitplot,
      filename = paste0(outputPrefix, "-fit.", plotFormat),
      width = plotDimensions$width,
      height = plotDimensions$height
    )
    return(fitplot)
  }

#' Plots the mean of the squared sum of residuals.
#' @param nlsFitPlotsOutputList the FIP fits output list.
#' @param outputPrefix the basename (identifying a case) for the plot.
#' @param plotFormat the format, passed to \code{ggplot2::ggsave}.
#' @param plotDimensions the dimensions of the plot.
#' @param color_scale the shared color scale to identify fragment and adduct pairs.
#' @return the ggplot object.
#' @importFrom magrittr %>%
#' @export
plotResidualsMeanSumSq <-
  function(nlsFitPlotsOutputList,
           outputPrefix,
           plotFormat = "png",
           plotDimensions = list(width = 11.69, height = 8.27),
           color_scale = ggplot2::scale_colour_hue()) {
    #residuals plot of mean sum of squared residuals
    data <- nlsFitPlotsOutputList$res_normality %>% tidyr::separate(
      combinationId,
      c(
        "species",
        "precursorAdduct",
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
    data$fragadd <- paste(data$fragment, data$adduct, sep = " ")
    data$fragadd <-
      factor(data$fragadd, levels = unique(data[order(data$calculatedMass), ]$fragadd))
    data$`foundMassRange[ppm]` <-
      factor(data$`foundMassRange[ppm]`,
             levels = sort(as.numeric(unique(data$`foundMassRange[ppm]`))))
    resplot <- ggplot2::ggplot(data = data) +
      ggplot2::geom_point(
        ggplot2::aes(
          x = meanResSSq,
          y = fragadd,
          shape = `foundMassRange[ppm]`,
          col = fragadd
        ),
        size = 2,
        data
      ) +
      ggplot2::theme_bw(base_size = 12, base_family = 'Helvetica') +
      ggplot2::labs(
        title = paste0(unique(data$species), " ", unique(data$group)),
        colour = 'Fragment',
        shape = 'Mass Range [ppm]'
      ) +
      ggplot2::ylab("Fragment + Adduct") +
      ggplot2::xlab(expression(paste(log[10], bar(mu), " Res. SSq", sep = " "))) +
      # ggplot2::scale_shape_manual(values=c(4, 1), breaks=c("TRUE", "FALSE"), labels=c("Normal","Non-Normal")) +
      ggplot2::coord_trans(x = "log10") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(
        angle = 90,
        vjust = 0.5,
        hjust = 1
      )) +
      color_scale
    ggplot2::ggsave(
      resplot,
      filename = paste0(outputPrefix, "-residuals-mean-ssq.", plotFormat),
      width = plotDimensions$width,
      height = plotDimensions$height
    )
    return(resplot)
  }

#' Plots the residuals quantile-quantile plot between standardized residuals and an assumed normal distribution.
#' @param nlsFitPlotsOutputList the FIP fits output list.
#' @param outputPrefix the basename (identifying a case) for the plot.
#' @param plotFormat the format, passed to \code{ggplot2::ggsave}.
#' @param plotDimensions the dimensions of the plot.
#' @param color_scale the shared color scale to identify fragment and adduct pairs.
#' @return the ggplot object.
#' @importFrom magrittr %>%
#' @export
plotResidualsQQ <-
  function(nlsFitPlotsOutputList,
           outputPrefix,
           plotFormat = "png",
           plotDimensions = list(width = 11.69, height = 8.27),
           color_scale = ggplot2::scale_colour_hue()) {
    #residuals quantile-quantile plot
    preds_from_data <-
      dplyr::left_join(
        nlsFitPlotsOutputList$preds_from_data,
        nlsFitPlotsOutputList$res_normality,
        by = c("combinationId")
      ) %>% tidyr::separate(
        combinationId,
        c(
          "species",
          "precursorAdduct",
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
    preds_from_data$`foundMassRange[ppm]` <-
      factor(preds_from_data$`foundMassRange[ppm]`,
             levels = sort(as.numeric(unique(preds_from_data$`foundMassRange[ppm]`))))
    preds_from_data$fragadd <-
      paste(preds_from_data$fragment, preds_from_data$adduct, sep = " ")
    preds_from_data$fragadd <-
      factor(preds_from_data$fragadd,
             levels = unique(preds_from_data[order(as.numeric(preds_from_data$calculatedMass)), ]$fragadd))
    resplot <- ggplot2::ggplot(preds_from_data) +
      ggplot2::aes(shape = isNormal,
                   sample = .std.resid,
                   col = fragadd) + ggplot2::stat_qq() + ggplot2::stat_qq_line(col = "red", linetype = "dashed") +
      ggplot2::facet_wrap(
        fragadd ~ `foundMassRange[ppm]`,
        labeller = ggplot2::label_wrap_gen(multi_line = FALSE),
        ncol = 6
      ) +
      ggplot2::theme_bw(base_size = 12, base_family = 'Helvetica') +
      ggplot2::labs(
        title = paste0(
          unique(preds_from_data$species),
          " ",
          unique(preds_from_data$group)
        ),
        colour = 'Fragment',
        shape = 'Res. S-W-Normality'
      ) +
      ggplot2::ylab(expression(paste(
        "Standardized Residuals (", Delta, "(", y, ",", yhat, ")", ")", sep = " "
      ))) +
      ggplot2::xlab('Theoretical Normal Distribution') +
      ggplot2::scale_shape_manual(
        values = c(4, 1),
        breaks = c("TRUE", "FALSE"),
        labels = c("Normal", "Non-Normal")
      ) +
      color_scale
    ggplot2::ggsave(
      resplot,
      filename = paste0(outputPrefix, "-residuals-qq.", plotFormat),
      width = plotDimensions$width,
      height = plotDimensions$height
    )
    return(resplot)
  }

#' Plots the residuals between predicted and measured values.
#' @param nlsFitPlotsOutputList the FIP fits output list.
#' @param outputPrefix the basename (identifying a case) for the plot.
#' @param plotFormat the format, passed to \code{ggplot2::ggsave}.
#' @param plotDimensions the dimensions of the plot.
#' @param color_scale the shared color scale to identify fragment and adduct pairs.
#' @return the ggplot object.
#' @importFrom magrittr %>%
#' @export
plotResiduals <-
  function(nlsFitPlotsOutputList,
           outputPrefix,
           plotFormat = "png",
           plotDimensions = list(width = 11.69, height = 8.27),
           color_scale = ggplot2::scale_colour_hue()) {
    #Residualplot
    preds_from_data <-
      nlsFitPlotsOutputList$preds_from_data %>% tidyr::separate(
        combinationId,
        c(
          "species",
          "precursorAdduct",
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
    preds_from_data$`foundMassRange[ppm]` <-
      factor(preds_from_data$`foundMassRange[ppm]`,
             levels = sort(as.numeric(unique(preds_from_data$`foundMassRange[ppm]`))))
    preds_from_data$fragadd <-
      paste(preds_from_data$fragment, preds_from_data$adduct, sep = " ")
    preds_from_data$fragadd <-
      factor(preds_from_data$fragadd,
             levels = unique(preds_from_data[order(as.numeric(preds_from_data$calculatedMass)), ]$fragadd))
    resplot <-
      ggplot2::ggplot(preds_from_data, ggplot2::aes(col = fragadd)) +
      ggplot2::geom_point(ggplot2::aes(precursorCollisionEnergy, .resid),
                          size = 2,
                          preds_from_data) +
      ggplot2::geom_rug(
        ggplot2::aes(precursorCollisionEnergy, .resid),
        alpha = 0.5,
        sides = "b",
        preds_from_data
      ) +
      ggplot2::geom_hline(yintercept = 0,
                          col = "red",
                          linetype = "dashed") +
      ggplot2::facet_wrap(
        fragadd ~ `foundMassRange[ppm]`,
        labeller = ggplot2::label_wrap_gen(multi_line = FALSE),
        ncol = 6
      ) +
      ggplot2::theme_bw(base_size = 12, base_family = 'Helvetica') +
      ggplot2::labs(title = paste0(
        unique(preds_from_data$species),
        " ",
        unique(preds_from_data$group)
      ), colour = 'Fragment') +
      ggplot2::ylab(expression(paste(
        "Residuals (", Delta, "(", y, ",", yhat, ")", ")", sep = " "
      ))) +
      ggplot2::xlab(flipr::collisionEnergyLabel(nlsFitPlotsOutputList$nls.tibble)) +
      ggplot2::scale_y_continuous(limits = c(-0.1, 0.1)) +
      color_scale
    # +
    #theme(legend.position = c(0.9, 0.15))
    ggplot2::ggsave(
      resplot,
      filename = paste0(outputPrefix, "-residuals.", plotFormat),
      width = plotDimensions$width,
      height = plotDimensions$height
    )
    return(resplot)
  }

#' Creates the x-axis plot label depending on precursorActivationType and precursorCollisionEnergyUnit as new columns precursorCollisionEnergyUnitLabel.
#' @param data the FIP input data.
#' @importFrom magrittr %>%
#' @export
collisionEnergyLabel <- function(data) {
  paste0(
    unique(data$precursorCollisionEnergyUnitLabel),
    " ",
    unique(data$precursorActivationType)
  )
}

#' Plots information about the fits (residuals, predictions, confidence intervals).
#' @param nlsFitPlotsOutputList the FIP fits output list.
#' @param outputPrefix the basename (identifying a case) for the plot.
#' @param plotFormat the format, passed to \code{ggplot2::ggsave}.
#' @param plotDimensions the dimensions of the plot.
#' @param color_scale the shared color scale to identify fragment and adduct pairs.
#' @importFrom magrittr %>%
#' @export
plotFits <-
  function(nlsFitPlotsOutputList,
           outputPrefix,
           plotFormat = "png",
           plotDimensions = list(width = 11.69, height = 8.27),
           color_scale = ggplot2::scale_colour_hue()) {
    plotNumberOfSamplesPerCombinationId(
      nlsFitPlotsOutputList = nlsFitPlotsOutputList,
      outputPrefix = outputPrefix,
      plotFormat = plotFormat,
      plotDimensions = plotDimensions,
      color_scale = color_scale
    )
    plotResiduals(
      nlsFitPlotsOutputList = nlsFitPlotsOutputList,
      outputPrefix = outputPrefix,
      plotFormat = plotFormat,
      plotDimensions = plotDimensions,
      color_scale = color_scale
    )
    plotResidualsQQ(
      nlsFitPlotsOutputList = nlsFitPlotsOutputList,
      outputPrefix = outputPrefix,
      plotFormat = plotFormat,
      plotDimensions = plotDimensions,
      color_scale = color_scale
    )
    plotResidualsMeanSumSq(
      nlsFitPlotsOutputList = nlsFitPlotsOutputList,
      outputPrefix = outputPrefix,
      plotFormat = plotFormat,
      plotDimensions = plotDimensions,
      color_scale = color_scale
    )
    plotPredictedFits(
      nlsFitPlotsOutputList = nlsFitPlotsOutputList,
      outputPrefix = outputPrefix,
      plotFormat = plotFormat,
      plotDimensions = plotDimensions,
      color_scale = color_scale
    )
    plotParameterConfidenceIntervals(
      nlsFitPlotsOutputList = nlsFitPlotsOutputList,
      combinationId = combinationId,
      outputPrefix = outputPrefix,
      plotFormat = plotFormat,
      plotDimensions = plotDimensions,
      color_scale = color_scale
    )
  }