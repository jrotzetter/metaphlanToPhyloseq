# Only keep shared samples between metadata and MetaPhlAn profile

The `clean_metadata()` function is designed to clean and prepare
metadata and a MetaPhlAn profile for further analysis by only keeping
shared samples.

## Usage

``` r
clean_metadata(metadata, mtphlan_profile, sample_column, remove_spaces = TRUE)
```

## Arguments

- metadata:

  A data frame containing metadata information.

- mtphlan_profile:

  A data frame containing metaphlan profile information. Samples are
  columns.

- sample_column:

  A `character` string. The column in the metadata containing the sample
  names. Should match column names of the MetaPhlAn profile.

- remove_spaces:

  Should spaces (" ") be replaced with "\_" in the column names?

## Value

The cleaned metadata as a data frame.

## Author

Jérémy Rotzetter

## Examples

``` r
metadata <- data.frame(
  Sample_name = c("A", "B", "C", "D", "F", "G"),
  Sex = c(rep(c("F", "M")))
)

mtphlan_profile <- data.frame(
  clade_name = c("x", "y", "z"),
  A = c(1, 2, 3),
  B = c(4, 5, 6),
  D = c(7, 8, 9),
  Z = c(10, 11, 12)
)

cleaned_metadata <- clean_metadata(metadata, mtphlan_profile, "Sample_name")
#> Metadata and profile Sample_names match.

print(cleaned_metadata)
#>   Sample_name Sex
#> 1           A   F
#> 2           B   M
#> 4           D   M
```
