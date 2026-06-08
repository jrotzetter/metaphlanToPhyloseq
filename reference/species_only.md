# Merged MetaPhlAn profiles data to species level

Example data of multiple MetaPhlAn profiles which were merged into one
file, and then filtered with `filter_taxa_lvl(df, taxa_lvl = "s")` to
only include species-level.

## Usage

``` r
species_only
```

## Format

### `species_only`

A data frame with 12 rows and 7 columns:

- clade_name:

  The taxonomic lineage of the taxon reported on this line. Taxon names
  are prefixed with one-letter codes to help indicate their rank.

- SRS0144...:

  The sample names for which sample-specific abundance profiles are
  found in a column.

## Source

<https://github.com/biobakery/biobakery/wiki/MetaPhlAn-4.1>
