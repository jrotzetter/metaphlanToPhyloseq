# Filter MetaPhlAn profile by taxonomic rank

Function to filter rows of a MetaPhlAn profile based on the selected
taxonomic level. Used if profile not yet separated into taxonomic ranks.

## Usage

``` r
filter_taxa_lvl(df, taxa_lvl)
```

## Arguments

- df:

  A data frame. Assumes clade names are of the structure
  ***k\_\_Bacteria\|p\_\_Firmicutes\|c\_\_Bacilli\|o\_\_Lactobacillales\|f\_\_Lactobacillaceae***
  and found in the 'clade_name' column.

- taxa_lvl:

  A `character` string. The taxonomic level to extract ('kingdom',
  'phylum', 'class', 'order', 'family', 'genus', 'species' or 't'
  (SGB)). First letter abbreviations (e.g., 's') are also accepted.

## Value

A filtered data frame.

## Author

Jérémy Rotzetter

## Examples

``` r
filtered_profile <- filter_taxa_lvl(merged_abundance_profiles, "Class")
head(filtered_profile)
#>                                             clade_name SRS014459.Stool
#> 1 k__Bacteria|p__Proteobacteria|c__Gammaproteobacteria               0
#> 2       k__Bacteria|p__Actinobacteria|c__Actinomycetia               0
#> 3           k__Bacteria|p__Firmicutes|c__Negativicutes               0
#> 4           k__Bacteria|p__Bacteroidota|c__Bacteroidia               0
#> 5                 k__Bacteria|p__Firmicutes|c__Bacilli               0
#>   SRS014464.Anterior_nares SRS014470.Tongue_dorsum SRS014472.Buccal_mucosa
#> 1                 81.20128                 0.00000                 14.7449
#> 2                 18.79872                 0.00000                  0.0000
#> 3                  0.00000                74.02531                  0.0000
#> 4                  0.00000                25.97469                  0.0000
#> 5                  0.00000                 0.00000                 85.2551
#>   SRS014476.Supragingival_plaque SRS014494.Posterior_fornix
#> 1                        0.00000                          0
#> 2                       85.43048                          0
#> 3                        0.00000                          0
#> 4                        0.00000                          0
#> 5                       14.56952                        100
```
