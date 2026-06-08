# Modify the case of the first letter in each word of a character string

These functions take a character string as input and return a new
character string with the first letter of each word modified to either
uppercase or lowercase.

## Usage

``` r
capitalize_str(character_string)

lowercase_str(character_string)
```

## Arguments

- character_string:

  A `character` vector containing the strings to be modified.

## Value

A character vector with the first letter of each word modified to either
uppercase or lowercase.

## Examples

``` r
capitalize_str("hello world") # Returns "Hello World"
#> [1] "Hello world"
capitalize_str(c("hello", "world")) # Returns c("Hello", "World")
#> [1] "Hello" "World"
lowercase_str("Hello World") # Returns "hello world"
#> [1] "hello World"
lowercase_str(c("Hello", "World")) # Returns c("hello", "world")
#> [1] "hello" "world"
```
