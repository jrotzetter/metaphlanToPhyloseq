# Set a column as row names

This function sets a specified column as the row names of a data frame.

## Usage

``` r
col_to_rownames(df, column)
```

## Arguments

- df:

  The input data frame.

- column:

  A `character` string. The name of the column to be set as row names.

## Value

The input data frame with the specified column set as row names.

## Examples

``` r
df <- data.frame(
  clade_name = c("A", "B", "C"),
  col1 = c(10, 20, 5),
  col2 = c(8, 15, 3),
  col3 = c(64, 89, 32)
)
col_to_rownames(df, "clade_name")
#>   col1 col2 col3
#> A   10    8   64
#> B   20   15   89
#> C    5    3   32
```
