# Convert MetaPhlAn profile to phyloseq object

This function converts a MetaPhlAn profile to a phyloseq object when
given the file path or a pre-loaded table.

## Usage

``` r
metaphlan_to_phyloseq(
  mtphlan_profile,
  taxa_lvl = NULL,
  metadata = NULL,
  sample_column = NULL,
  use_taxa_names = FALSE,
  merged_profiles = TRUE
)
```

## Arguments

- mtphlan_profile:

  The MetaPhlAn profile to be converted. It can be either a file path or
  a data frame of MetaPhlAn profile(s).

- taxa_lvl:

  Optional taxonomic level to filter the profile to. Valid options are
  'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species' or
  't' (SGB). First letter abbreviations are also accepted.

- metadata:

  Optional metadata for the samples. If provided, it should be a data
  frame.

- sample_column:

  `Character` string. The column in the metadata containing the sample
  names. Should match column names of the MetaPhlAn profile.

- use_taxa_names:

  `Logical` indicating whether to use taxonomic names instead of OTUs in
  the resulting phyloseq object. Default is `FALSE`.

- merged_profiles:

  `Logical`; if `TRUE` (the default) the file to be loaded is assumed to
  be multiple merged MetaPhlAn profiles.

## Value

A phyloseq object representing the MetaPhlAn profile.

## Note

If `mtphlan_profile` is an already loaded object, the clade names should
not be shortened, otherwise the taxonomic table cannot be created.
Should you wish to use the taxonomic names of the specified rank as the
row names in the phyloseq object, set `use_taxa_names = TRUE`.

## Author

Jérémy Rotzetter

## Examples

``` r
metaphlan_to_phyloseq(species_only, taxa_lvl = "s")
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 12 taxa and 6 samples ]
#> tax_table()   Taxonomy Table:    [ 12 taxa by 7 taxonomic ranks ]
```
