#!/bin/bash
R -e "library('devtools')" -e "library('pkgdown')" -e "devtools::document()" -e "pkgdown::build_site()"
