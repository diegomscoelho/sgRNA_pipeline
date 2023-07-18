FROM r-base:4.2.1

# Install bioconductor
RUN R -e 'install.packages("remotes");' 
RUN R -e 'remotes::install_version("dplyr", version = "1.1.0")'
RUN R -e 'remotes::install_version("tidyr", version = "1.3.0")'

WORKDIR /home

# Default command to execute at startup of the container
CMD R --no-save
