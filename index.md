# metaphlanToPhyloseq

## Overview

metaphlanToPhyloseq is a simple R package to transform **MetaPhlAn 4**
taxonomic microbiome abundance profiles into the right format for easy
creation of a phyloseq object.

***NOTE**: The package is stable and fully functional but is no longer
under active development. Updates and support will be provided on a
best-effort basis.*

***Recommendation**: Users looking for more advanced features and
expanded capabilities are encouraged to explore the
[microEDA](https://github.com/jrotzetter/microEDA) package.*

## Installation

You can install the development version of metaphlanToPhyloseq from
[GitHub](https://github.com/jrotzetter/metaphlanToPhyloseq) with:

``` r

# install.packages("pak")
pak::pak("jrotzetter/metaphlanToPhyloseq")
```

Or the `remotes` package:

``` r

# install.packages("remotes")
remotes::install_github("jrotzetter/metaphlanToPhyloseq", build_vignettes = TRUE)
```

## Dependencies

- dplyr (\>= 1.1.4)
- phyloseq (\>= 1.46.0)
- utils

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
    use_taxa_names = TRUE,
    merged_profiles = FALSE
)

physeq_merged <- metaphlan_to_phyloseq(
  merged_profiles,
  taxa_lvl = "Species"
)

# Above steps are optional as the filtering to a specific taxonomic rank can be
# done directly from within the function or not at all
physeq <- metaphlan_to_phyloseq(
  merged_abundance_profiles,
  taxa_lvl = NULL,
  use_taxa_names = TRUE
)
```

For more details please see
[`vignette("metaphlanToPhyloseq")`](https://jrotzetter.github.io/metaphlanToPhyloseq/articles/metaphlanToPhyloseq.md)
or the help pages in the documentation. Both are also available online
at <https://jrotzetter.github.io/metaphlanToPhyloseq/>.

## Getting help

If you encounter a bug, please file an issue with a minimal reproducible
example on
[GitHub](https://github.com/jrotzetter/metaphlanToPhyloseq/issues). For
questions or help with MetaPhlAn, please visit the corresponding
[bioBakery](https://forum.biobakery.org/c/microbial-community-profiling/metaphlan/7)
forum. For help with phyloseq, helpful tutorials and articles can be
found on [GitHub](https://joey711.github.io/phyloseq/index.html).

## License

Released under [MIT](https://choosealicense.com/licenses/mit/) by
[@jrotzetter](https://github.com/jrotzetter).

This license means:

- You can freely copy, modify, distribute and reuse this software.
- The *original license* must be included with copies of this software.
- Please *link back* to this repo if you use a significant portion of
  the source code.
- The software is provided “as is”, without warranty of any kind.
