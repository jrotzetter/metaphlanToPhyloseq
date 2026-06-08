# Apply threshold filtering to a data frame

This function applies threshold filtering to a data frame, summing
values below the threshold for each column, storing the sum in the new
'Other' row, but only if all values within a row (taxon) across all
columns are below the threshold.

## Usage

``` r
filter_threshold(df, threshold, merged_profiles = TRUE)
```

## Arguments

- df:

  The input data frame of merged MetaPhlAn profiles.

- threshold:

  The threshold value for filtering.

- merged_profiles:

  `Logical`; if `TRUE` (the default) the file to be loaded is assumed to
  be multiple merged MetaPhlAn profiles.

## Value

The filtered data frame.

## Details

Function can be applied to a single MetaPhlAn profile by setting
`merged_profiles = FALSE`. Likewise also works with merged profiles by
slicing the data frame to only the clade_name and the numeric sample
column of interest. In both cases a `clade_name` column is required.

## Author

Jérémy Rotzetter

## Examples

``` r
df <- data.frame(
  clade_name = c("A", "B", "C", "D"),
  col1 = c(10, 20, 5, 65),
  col2 = c(8, 15, 3, 38),
  col3 = c(5, 35, 4, 6)
)
print(df)
#>   clade_name col1 col2 col3
#> 1          A   10    8    5
#> 2          B   20   15   35
#> 3          C    5    3    4
#> 4          D   65   38    6
df_thresh <- filter_threshold(df, 11)
print(df_thresh)
#>   clade_name col1 col2 col3
#> 1          B   20   15   35
#> 2          D   65   38    6
#> 3      Other   15   11    9
```
