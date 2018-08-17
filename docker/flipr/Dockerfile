# The suggested name for this image is: do1-aps-feris.isas.de:5000/isas/suprg-base

FROM rocker/tidyverse

MAINTAINER nils.hoffmann@isas.de

ADD docker/flipr/install.R /tmp/

# invalidates cache every 24 hours
ADD http://master.bioconductor.org/todays-date /tmp/

RUN apt-get update && apt-get install -y -t libssl-dev libxml2-dev

#need to install the packages globally
RUN R -f /tmp/install.R
