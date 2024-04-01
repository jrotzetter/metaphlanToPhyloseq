
<!-- README.md is generated from README.Rmd. Please edit that file -->

# metaphlanToPhyloseq

<!-- badges: start -->

[![GitHub
release](https://img.shields.io/github/release/jrotzetter/metaphlanToPhyloseq?include_prereleases=&sort=semver&color=blue)](https://github.com/jrotzetter/metaphlanToPhyloseq/releases/ "View releases")
[![License](https://img.shields.io/badge/License-MIT-blue)](#license "View license summary")
[![issues -
metaphlanToPhyloseq](https://img.shields.io/github/issues/jrotzetter/metaphlanToPhyloseq)](https://github.com/jrotzetter/metaphlanToPhyloseq/issues "View open issues")
[![Made with
R](https://img.shields.io/badge/R-4.3.3-blue?logo=r&logoColor=white)](https://cran.r-project.org/ "Go to CRAN homepage")
[![Made with
R](https://img.shields.io/badge/RStudio-2023.12.1_Build_402-blue?logo=rstudio&logoColor=white)](https://posit.co/products/open-source/rstudio/ "Go to RSTUDIO IDE homepage")
<!-- badges: end -->

## Overview

metaphlanToPhyloseq is a simple R package to transform *MetaPhlAn 4*
taxonomic microbiome abundance profiles into the right format for easy
creation of a phyloseq object.

## Installation

You can install the development version of metaphlanToPhyloseq from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("jrotzetter/metaphlanToPhyloseq", build_vignettes = TRUE)
```

Alternatively you can also use:

``` r
# install.packages("pak")
pak::pak("jrotzetter/metaphlanToPhyloseq")
```

## Dependencies

- dplyr (\>= 1.1.4)
- phyloseq (\>= 1.46.0)
- utils (\>= 4.3.3)

## Usage

``` r
library(metaphlanToPhyloseq)

# Filter data to only include the specified taxonomic rank
single_profile <- filter_taxa_lvl(df = single_abundance_profile, taxa_lvl = "Genus")
merged_profiles <- filter_taxa_lvl(df = merged_abundance_profiles, taxa_lvl = "s")

# Keep only the columns of interest (clade_name and relative_abundance)
single_profile <- single_profile[, c(1, 3)]

# Create a phyloseq object
physeq_single <- metaphlan_to_phyloseq(
  mtphlan_profile = single_profile,
  taxa_lvl = "genus",
  use_taxa_names = TRUE
)

physeq_merged <- metaphlan_to_phyloseq(
  merged_profiles,
  taxa_lvl = "Species"
)
```

## To-do

- Deploy pkgdown/GitHub Pages
- Add valid, package-specific maintainer e-mail address
- Rewrite some functions (e.g., shorten_clade_names()) to also work with
  data that wasn’t pre-filtered to a chosen taxonomic rank when not
  directly loading a file from a path
- Allow metadata to be loaded directly from path in
  metaphlan_to_phyloseq()
- Add tests

## Getting help

If you encounter a bug, please file an issue with a minimal reproducible
example on
[GitHub](https://github.com/jrotzetter/metaphlanToPhyloseq/issues). For
questions or help with MetaPhlAn, please visit the corresponding
[bioBakery](https://forum.biobakery.org/c/microbial-community-profiling/metaphlan/7)
forum. For help with phyloseq, helpful tutorials and articles can be
found on [GitHub](https://joey711.github.io/phyloseq/index.html)

## License

Released under [MIT](https://choosealicense.com/licenses/mit/) by
[@jrotzetter](https://github.com/jrotzetter).

This license means:

- You can freely copy, modify, distribute and reuse this software.
- The *original license* must be included with copies of this software.
- Please *link back* to this repo if you use a significant portion of
  the source code.
- The software is provided “as is”, without warranty of any kind.
