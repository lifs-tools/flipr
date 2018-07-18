# FlipR Package #

## Installation ##

Install devtools:

```R
  install.packages("devtools")
```
  
Run

```R
  install_git("https://gitlab.isas.de/hoffmann/flipr.git")
```

Done!

## Usage ##

```R
  library("flipr")
  prm_file <- system.file("extdata", "12-HETE-d8_fip.tsv", package = "flipr", mustWork = TRUE)
  fits <- flip(projectDir = dirname(prm_file), plotFormat = "pdf", filePattern = basename(prm_file), dataPlots = FALSE)
```

## Installation for development ##
flipr uses packrat to manage its dependencies in a clean state during development. This does mean, that you need to install flipr into your local R installation, if you want to use your local development version outside of the flipr directory.

```R
   library(devtools)
   install_local("/path/to/flipr/packrat/lib/x86_64-pc-linux-gnu/3.4.4/flipr/")
```