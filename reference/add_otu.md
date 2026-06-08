# Add OTUs as row names

This function adds numbered Operational Taxonomic Units (OTUs) to the
row names of a data frame. Optionally can also add a OTU column to the
data frame.

## Usage

``` r
add_otu(df, add_column = FALSE)
```

## Arguments

- df:

  The input MetaPhlAn data frame.

- add_column:

  `Logical` value indicating whether to add an OTU column to the data
  frame. Default is `FALSE`.

## Value

The data frame with OTU row names added.

## Author

Jérémy Rotzetter

## Examples

``` r
df <- data.frame(x = c(1, 2, 3), y = c(4, 5, 6))
add_otu(df) # Adds OTU row names to the data frame
#>      x y
#> Otu1 1 4
#> Otu2 2 5
#> Otu3 3 6
add_otu(df, add_column = TRUE) # Adds OTU column and row names to the data frame
#>      x y  Otu
#> Otu1 1 4 Otu1
#> Otu2 2 5 Otu2
#> Otu3 3 6 Otu3
```
