library(ggplot2)
library(dplyr)

varyNormParm <- function() {
  xrange = c(-1, 15)
  n = 300
  dist.fun = "dnorm"
  mean = seq(from = 0, to = 10, by = 0.5)
  sd = seq(from = 0, to = 10, by = 1)
  grid <- expand.grid(mean = mean, sd = sd)

  distFrame <- data.frame()

  for (i in 1:nrow(grid)) {
    grid[i, ]
    df <- data.frame(
      x = seq(
        from = xrange[[1]],
        to = xrange[[2]],
        length.out = n
      ),
      mean = as.factor(grid[i, ]$mean),
      sd = as.factor(grid[i, ]$sd)
    )
    df$y <-
      do.call(dist.fun, list(
        x = df$x,
        mean = grid[i, ]$mean,
        sd = grid[i, ]$sd
      ))
    distFrame <- rbind(distFrame, df)
  }
  colnames(distFrame) <- c("x", "mean", "sd", "y")
  #View(distFrame)
  ggplot2::qplot(
    x = x,
    y = y,
    color = sd,
    data = distFrame,
    geom = c("line")
  ) + ggplot2::theme(legend.position = 'right') + ggplot2::facet_wrap( ~ mean)
}

varyWeibullParm <- function() {
  xrange = c(-1, 15)
  n = 300
  dist.fun = "dweibull"
  shape = seq(from = 0, to = 10, by = 0.5)
  scale = seq(from = 0, to = 10, by = 1)
  grid <- expand.grid(shape = shape, scale = scale)
  #head(grid)
  #View(grid)
  distFrame <- data.frame()

  #View(distFrame)
  for (i in 1:nrow(grid)) {
    grid[i, ]
    df <- data.frame(
      x = seq(
        from = xrange[[1]],
        to = xrange[[2]],
        length.out = n
      ),
      shape = as.factor(grid[i, ]$shape),
      scale = as.factor(grid[i, ]$scale)
    )
    df$y <-
      do.call(dist.fun,
              list(
                x = df$x,
                shape = grid[i, ]$shape,
                scale = grid[i, ]$scale
              ))
    distFrame <- rbind(distFrame, df)
  }
  colnames(distFrame) <- c("x", "shape", "scale", "y")
  #View(distFrame)
  ggplot2::qplot(
    x = x,
    y = y,
    color = shape,
    data = distFrame,
    geom = c("line")
  ) + ggplot2::theme(legend.position = 'right') + ggplot2::facet_wrap( ~ scale)
  #    plotDist(xrange=xrange, n = n, dist.fun=dist.fun, sdlog=sdlogvalue, meanlog = meanlog)

}

varyDlnormParm <- function() {
  xrange = c(-1, 50)
  n = 300
  dist.fun = "dlnormPar"

  meanlog = seq(from = 0, to = 10, by = 0.5)
  sdlog = seq(from = 1, to = 5, by = 0.5)
  grid <- expand.grid(meanlog = meanlog, sdlog = sdlog)
  #head(grid)
  #View(grid)
  distFrame <- data.frame()

  #View(distFrame)
  for (i in 1:nrow(grid)) {
    grid[i, ]
    df <- data.frame(
      x = seq(
        from = xrange[[1]],
        to = xrange[[2]],
        length.out = n
      ),
      sdlog = as.factor(grid[i, ]$sdlog),
      meanlog = as.factor(grid[i, ]$meanlog)
    )
    df$y <-
      do.call(dist.fun,
              list(
                x = df$x,
                sdlog = grid[i, ]$sdlog,
                meanlog = grid[i, ]$meanlog
              ))
    distFrame <- rbind(distFrame, df)
  }
  colnames(distFrame) <- c("x", "meanlog", "sdlog", "y")
  #View(distFrame)
  ggplot2::qplot(
    x = x,
    y = y,
    color = meanlog,
    data = distFrame,
    geom = c("line")
  ) + ggplot2::theme(legend.position = 'right') + ggplot2::facet_wrap( ~ sdlog)

}

plotDist <- function(xrange = c(-10, 10),
                     n = 500,
                     dist.fun,
                     ...) {
  df <-
    data.frame(x = seq(
      from = xrange[[1]],
      to = xrange[[2]],
      length.out = n
    ))
  #df$y <-
  #print(df)
  #integral <- integrate(f=eval(dist.fun), lower = xrange[[1]], upper = xrange[[2]])
  #print(integral)
  plot <- ggplot2::ggplot(df, ggplot2::aes(x = x)) +
    ggplot2::stat_function(fun = eval(dist.fun), args = list(...))
  plot
}


fit <- function(tibble, dist.fun, ...) {
  tibble$fragment <- as.factor(tibble$fragment)
  tibble$`foundMassRange[ppm]` <-
    as.factor(tibble$`foundMassRange[ppm]`)
  lvl <- levels(tibble$fragment)
  #print(paste("Levels of fragment:", lvl, sep = " "))
  # plots_neg_frags <- lapply(dfl_neg_frags, function(dt) {
  #   plot1 <- qplot(x=precursorCollisionEnergy, y=scanRelativeIntensity, data=dt)
  #   plot1
  # })
  #grid.arrange(grobs=plots_neg_frags, ncol=1)
  # plots_frags <- lapply(dfl_polarities, function(dt) {
  #
  splits <- split(tibble, tibble$`foundMassRange[ppm]`)
  plots <- lapply(splits, function(splitset) {
    polFraAdd <-
      splitset %>% dplyr::group_by(polarity, `foundMassRange[ppm]`, fragment, adduct)
    #print(paste("Groups:", groups(polFraAdd), sep = " "))
    #calculate fit based on estimated distribution coefficients and selected distribution
    polFraAddSummary <-
      polFraAdd %>% dplyr::summarise(
        MeanEpc = mean(precursorCollisionEnergy),
        sdEpc = sd(precursorCollisionEnergy)
      )
    #print(polFraAddSummary)
    plot <- ggplot2::qplot(
      x = precursorCollisionEnergy,
      y = scanRelativeIntensity,
      data = splitset,
      color = fragment
    ) + ggplot2::facet_wrap(fragment + adduct ~ polarity) +
      ggplot2::ggtitle(paste0(unique(splitset$`foundMassRange[ppm]`), " [ppm]")) +
      ggplot2::stat_function(fun = eval(dist.fun),
                    args = list(...),
                    colour = "blue")# + geom_quantile(colour="red")
    plot
  })
  # plot1
  # })
  return(plots)
}
