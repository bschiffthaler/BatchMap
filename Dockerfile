FROM rocker/r-base

RUN apt update && apt install -y libssl-dev libxml2-dev libcurl4-openssl-dev

RUN Rscript -e 'install.packages(c("tidyverse", "devtools", "roxygen2", "Rcpp", "RcppArmadillo"))'

RUN Rscript -e 'devtools::install_github("bschiffthaler/BatchMap")'

ENTRYPOINT ["R"]
