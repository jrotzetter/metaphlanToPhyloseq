# SRS014470-Tongue_dorsum_profile data

Example data of a single MetaPhlAn profile output file, containing the
computed taxon abundances.

## Usage

``` r
single_abundance_profile
```

## Format

### `single_abundance_profile`

A data frame with 15 rows and 4 columns:

- clade_name:

  The taxonomic lineage of the taxon reported on this line. Taxon names
  are prefixed with one-letter codes to help indicate their rank.

- NCBI_tax_id:

  The NCBI-equivalent taxon IDs of the named taxa from clade_name.

- relative_abundance:

  The taxon's relative abundance in %.

- additional_species:

  Additional species names for cases where the metagenome profile
  contains clades that represent multiple species. The species listed in
  column 1 is the representative species in such cases.

## Source

<https://github.com/biobakery/biobakery/wiki/MetaPhlAn-4.1>
