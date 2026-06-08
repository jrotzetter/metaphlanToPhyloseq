# Remove pattern from column names

Function to remove a pattern from column names.

## Usage

``` r
clean_colnames(df, pattern)
```

## Arguments

- df:

  A data frame or matrix.

- pattern:

  A `character` string with the pattern to match.

## Value

Returns the input table with the specified pattern removed from the
column names.

## Author

Jérémy Rotzetter

## Examples

``` r
colnames(merged_abundance_profiles)
#> [1] "clade_name"                     "SRS014459.Stool"               
#> [3] "SRS014464.Anterior_nares"       "SRS014470.Tongue_dorsum"       
#> [5] "SRS014472.Buccal_mucosa"        "SRS014476.Supragingival_plaque"
#> [7] "SRS014494.Posterior_fornix"    

names_cleaned <- clean_colnames(merged_abundance_profiles, "SRS0144\\d{2}.")

colnames(names_cleaned)
#> [1] "clade_name"           "Stool"                "Anterior_nares"      
#> [4] "Tongue_dorsum"        "Buccal_mucosa"        "Supragingival_plaque"
#> [7] "Posterior_fornix"    
```
