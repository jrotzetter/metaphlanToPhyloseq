# Shorten clade names in a dataset to chosen taxonomic level

This function shortens the taxonomic names of clades in a given dataset
based on a specified taxonomic level, using the first letter of the
taxonomic rank + "\_\_" as rank identifiers.

## Usage

``` r
shorten_clade_names(
  data,
  taxa_lvl,
  apply_to_colnames = TRUE,
  selected_cols = NULL
)
```

## Arguments

- data:

  The input dataset. Assumes clade names are of the structure
  ***k\_\_Bacteria\|p\_\_Firmicutes\|c\_\_Bacilli\|o\_\_Lactobacillales\|f\_\_Lactobacillaceae***
  and found in the 'clade_name' column or the column names themselves.

- taxa_lvl:

  The taxonomic level at which the clade names should be shortened.
  Valid options include 'kingdom', 'phylum', 'class', 'order', 'family',
  'genus', 'species' or 't' (SGB). First letter abbreviations (e.g.,
  's') are also accepted.

- apply_to_colnames:

  `Logical` indicating whether the shortening should be applied to
  column names or row values. Default is `TRUE`.

- selected_cols:

  A `character` vector specifying the columns to which the shortening
  should be applied. If `NULL` (the default), the shortening is applied
  to all columns.

## Value

The dataset with the clade names shortened based on the specified
taxonomic level.

In the case where there are entries not matching the chosen taxonomic
rank, these are either returned 'as is', or if they follow the same
structure, the name will be shortened to the last taxonomic entry (see
rows 2, 5 and 6 of the example)

## Note

This function is not intended to be used with the workflow for the
creation of phyloseq objects as the full sequence of taxonomic names is
needed for the creation of the taxonomy table in
[`get_taxa_table()`](https://jrotzetter.github.io/metaphlanToPhyloseq/reference/get_taxa_table.md).
It may however be useful for analyses or plots created directly
with/from the dataframes.

This function uses the dplyr package for data manipulation.

## Author

Jérémy Rotzetter

## Examples

``` r
head(merged_abundance_profiles$clade_name)
#> [1] "UNCLASSIFIED"                                        
#> [2] "k__Bacteria"                                         
#> [3] "k__Bacteria|p__Proteobacteria"                       
#> [4] "k__Bacteria|p__Actinobacteria"                       
#> [5] "k__Bacteria|p__Proteobacteria|c__Gammaproteobacteria"
#> [6] "k__Bacteria|p__Actinobacteria|c__Actinomycetia"      

taxa_shortened <- shorten_clade_names(
  merged_abundance_profiles,
  "Phylum",
  apply_to_colnames = FALSE,
  selected_cols = "clade_name"
)

head(taxa_shortened$clade_name)
#> [1] "UNCLASSIFIED"           "k__Bacteria"            "Proteobacteria"        
#> [4] "Actinobacteria"         "c__Gammaproteobacteria" "c__Actinomycetia"      
```
