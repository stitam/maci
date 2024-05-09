FROM debian:bookworm-20230904

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update

# Install R and R packages from apt
RUN apt-get install -y \
  r-base \
  r-bioc-biostrings \
  r-cran-ape \
  r-cran-biocmanager \
  r-cran-data.table \
  r-cran-devtools \
  r-cran-ggrepel \
  r-cran-hmisc \
  r-cran-htmlwidgets \
  r-cran-igraph \
  r-cran-lubridate \
  r-cran-magrittr \
  r-cran-reshape2 \
  r-cran-rjava \
  r-cran-seqinr \
  r-cran-stringdist \
  r-cran-tibble \
  r-cran-tidyr \
  r-cran-tidyverse

# Install ggmap
RUN R -e 'options(warn = 2); install.packages("cli", dependencies = TRUE)'
RUN R -e 'options(warn = 2); install.packages("ggmap", dependencies = TRUE)'

# Install leaflet
RUN apt-get install -y \
  libgdal-dev \
  libudunits2-dev
RUN R -e 'options(warn = 2); install.packages("leaflet", dependencies = TRUE)'

#Install R packages from CRAN
RUN R -e 'options(warn = 2); install.packages("leaflet.minicharts", dependencies = TRUE)'
RUN R -e 'options(warn = 2); install.packages("tictoc", dependencies = TRUE)'
RUN R -e 'options(warn = 2); install.packages("tidytree", dependencies = TRUE)'
RUN R -e 'options(warn = 2); install.packages("treedater", dependencies = TRUE)'

# Install vegan
RUN apt-get install -y r-cran-vegan

# Install scatterpie
RUN R -e 'options(warn = 2); install.packages("scatterpie", dependencies = TRUE)'

# Install optparse
RUN apt-get install -y r-cran-optparse

# Install ggpubr, scales
RUN apt-get install -y r-cran-ggpubr r-cran-scales

# Install ggplotify
RUN R -e 'options(warn = 2); install.packages("ggplotify", dependencies = TRUE)'

# Install circlize, complexheatmap
RUN apt-get install -y r-cran-circlize r-bioc-complexheatmap

# Install ggpubr
RUN apt-get install -y cmake r-cran-ggpubr

# Install countrycode
RUN R -e 'options(warn = 2); install.packages("countrycode", dependencies = TRUE)'

# Install abdiv
RUN R -e 'options(warn = 2); install.packages("abdiv", dependencies = TRUE)'

# Install extrafont
RUN R -e 'options(warn = 2); install.packages("extrafont", dependencies = TRUE)'

# Install new fonts
RUN apt-get install -y cabextract wget
RUN wget http://ftp.br.debian.org/debian/pool/contrib/m/msttcorefonts/ttf-mscorefonts-installer_3.8_all.deb
RUN dpkg -i ttf-mscorefonts-installer_3.8_all.deb
RUN R -e "options(warn = 2); extrafont::font_import(path = '/usr/share/fonts/truetype/msttcorefonts', prompt = FALSE)"

# Install ggh4x
RUN R -e 'options(warn = 2); install.packages("ggh4x", dependencies = TRUE)'

# Install qualpalr from GitHub instead of CRAN
RUN R -e 'options(warn = 2); devtools::install_github("jolars/qualpalr@v.0.4.3")'

# Install ggtree
RUN R -e 'options(warn = 2); devtools::install_github("YuLab-SMU/ggtree@b0c16e561989d7c503eb942ce653b6b75cd65db8")'

# Install ggimage, ggnewscale, ggtreeExtra
RUN R -e 'options(warn = 2); install.packages("ggimage", dependencies = TRUE)'
RUN R -e 'options(warn = 2); install.packages("ggnewscale", dependencies = TRUE)'
RUN R -e 'options(warn = 2); BiocManager::install("ggtreeExtra", update = FALSE, ask = FALSE)'

# Install specific verion of ggplot2
RUN R -e 'options(warn = 2); devtools::install_github("tidyverse/ggplot2@v3.4.1")'

# Install specific version of dplyr
RUN R -e 'options(warn = 2); devtools::install_github("tidyverse/dplyr@v1.1.3")'

# Install specific version of ggtext
RUN R -e 'options(warn = 2); devtools::install_github("wilkelab/ggtext@0.1.0")'

# Install naniar
RUN R -e 'options(warn = 2); install.packages("naniar", dependencies = TRUE)'

# Set R library path to avoid binding from host system
ENV R_LIBS_USER /usr/local/lib/R/site-library