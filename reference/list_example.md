# List example files

This function lists available example files, found in the `inst/extdata`
directory.

## Usage

``` r
list_example(path = NULL)
```

## Source

This function was adapted from `readxl::readxl_example()` and originally
created by *Hadley Wickham & Jennifer Bryan*.

## Arguments

- path:

  The name of an example file. If not provided, the function lists the
  available example files.

## Value

A character vector.

## Examples

``` r
list_example() # Lists example files contained in the package
#> [1] "SRS014470-Tongue_dorsum_profile.txt" "merged_abundance_table.txt"         
list_example("merged_abundance_table.txt") # Lists the path of the specified file
#> [1] "/home/runner/work/_temp/Library/metaphlanToPhyloseq/extdata/merged_abundance_table.txt"
```
