# The suggested name for this image is: do1-aps-feris.isas.de:5000/isas/suprg-base

FROM rocker/tidyverse:3.6.2

MAINTAINER nils.hoffmann@isas.de

ADD install.R /tmp/

# invalidates cache every 24 hours
ADD http://master.bioconductor.org/todays-date /tmp/

#need to install the packages globally
RUN R -f /tmp/install.R
LABEL application.version="1.0.6-testing"
LABEL application.build.timestamp="2020-01-22T17:15:11+01:00"
LABEL application.name="flipr"
LABEL application.authors="Nils Hoffmann"
LABEL application.maintainer=""
LABEL application.license="MIT"
LABEL application.from.image=""
LABEL application.scm=""
