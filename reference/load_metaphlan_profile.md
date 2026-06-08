# Load MetaPhlAn profile(s)

Function to load a MetaPhlAn profile from a file path.

## Usage

``` r
load_metaphlan_profile(path, merged_profiles = TRUE)
```

## Arguments

- path:

  A `character` string containing the path to the file to read.

- merged_profiles:

  `Logical`; if `TRUE` (the default) the file to be loaded is assumed to
  be multiple merged MetaPhlAn profiles.

## Value

A data frame containing the MetaPhlAn results.

## Details

If a single MetaPhlAn profile is to be loaded, `merged_profiles` should
be set to `FALSE`. THe single profile is assumed to contain the column
names `clade_name`, `NCBI_tax_id`, `relative_abundance` and
`additional_species`.

## Author

Jérémy Rotzetter

## Examples

``` r
path <- list_example("merged_abundance_table.txt")
data <- load_metaphlan_profile(path = path)
head(data)
#>                                             clade_name SRS014459.Stool
#> 1                                         UNCLASSIFIED             100
#> 2                                          k__Bacteria               0
#> 3                        k__Bacteria|p__Proteobacteria               0
#> 4                        k__Bacteria|p__Actinobacteria               0
#> 5 k__Bacteria|p__Proteobacteria|c__Gammaproteobacteria               0
#> 6       k__Bacteria|p__Actinobacteria|c__Actinomycetia               0
#>   SRS014464.Anterior_nares SRS014470.Tongue_dorsum SRS014472.Buccal_mucosa
#> 1                  0.00000                       0                  0.0000
#> 2                100.00000                     100                100.0000
#> 3                 81.20128                       0                 14.7449
#> 4                 18.79872                       0                  0.0000
#> 5                 81.20128                       0                 14.7449
#> 6                 18.79872                       0                  0.0000
#>   SRS014476.Supragingival_plaque SRS014494.Posterior_fornix
#> 1                        0.00000                          0
#> 2                      100.00000                        100
#> 3                        0.00000                          0
#> 4                       85.43048                          0
#> 5                        0.00000                          0
#> 6                       85.43048                          0
```
