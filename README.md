# FlipR Package #

## Installation ##

Install devtools:

  install.packages("devtools")
  
Run

  install_git("https://gitlab.isas.de/hoffmann/flipr.git")

Done!

## Usage ##

  library("flipr")

  prm_file <- system.file("extdata", "12-HETE-d8_fip.tsv", package = "flipr", mustWork = TRUE)

  fits <- flip(projectDir = dirname(prm_file), plotFormat = "pdf", filePattern = basename(prm_file), dataPlots = FALSE)
