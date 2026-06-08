# Introduction to metaphlanToPhyloseq

metaphlanToPhyloseq is a simple R package to transform *MetaPhlAn 4*
taxonomic microbiome abundance profiles into the right format for easy
creation of a phyloseq object. It also comes bundled with a few useful
functions for easy pattern removal from column names or row values, as
well as filtering for taxa below a specified threshold and bundling them
together in their own entry.

Principally there are two ways to use this package: 1. Automatic:
Loading the profile(s) directly from a file path and convert directly to
a phyloseq object. 2. Manual: Import data separately into R before
phyloseq creation. This may require additional steps but potentially
allows for greater control of which data is included.

Both ways will now shortly be illustrated with example data included in
the package. But first let’s load the package.

``` r

library(metaphlanToPhyloseq)
```

## 1. Automatic: Directly from path

This part is relatively straight forward. Simply specify the file path
to the MetaPhlAn profile(s) to be converted with `mtphlan_profile`. If
multiple profiles (`merged_profiles = TRUE`) are used and metadata for
the samples is available, add it under `metadata` and point to the
column containing the sample names in `sample_column`. The function will
then compare the profile sample names with those found in the metadata
and only keep those that are shared between the two. `taxa_lvl` is used
to specify the taxonomic level to filter the profile to (optional).
Valid options are ‘kingdom’, ‘phylum’, ‘class’, ‘order’, ‘family’,
‘genus’, ‘species’ or ‘t’ for Species Genome Bin (SGB). First letter
abbreviations are also accepted (except for SGB). With
`use_taxa_names = TRUE` it is possible to use the taxonomic names
instead of numbered OTUs in the resulting phyloseq object. Compare the
row names between the otu table and taxonomy table below.

``` r

path <- system.file("extdata/SRS014470-Tongue_dorsum_profile.txt", package = "metaphlanToPhyloseq")
single_profile <- metaphlan_to_phyloseq(
  mtphlan_profile = path,
  metadata = NULL,
  taxa_lvl = "Species",
  sample_column = NULL,
  use_taxa_names = TRUE,
  merged_profiles = FALSE
)
head(single_profile@otu_table)
#> OTU Table:          [2 taxa and 1 samples]
#>                      taxa are rows
#>                      relative_abundance
#> Veillonella_dispar             74.02531
#> Prevotella_histicola           25.97469

path <- system.file("extdata/merged_abundance_table.txt", package = "metaphlanToPhyloseq")
merged_profiles <- metaphlan_to_phyloseq(
  path,
  taxa_lvl = "o",
  use_taxa_names = FALSE
)
head(merged_profiles@tax_table)
#> Taxonomy Table:     [6 taxa by 4 taxonomic ranks]:
#>      Kingdom    Phylum           Class                 Order              
#> Otu1 "Bacteria" "Proteobacteria" "Gammaproteobacteria" "Moraxellales"     
#> Otu2 "Bacteria" "Actinobacteria" "Actinomycetia"       "Corynebacteriales"
#> Otu3 "Bacteria" "Firmicutes"     "Negativicutes"       "Veillonellales"   
#> Otu4 "Bacteria" "Bacteroidota"   "Bacteroidia"         "Bacteroidales"    
#> Otu5 "Bacteria" "Firmicutes"     "Bacilli"             "Lactobacillales"  
#> Otu6 "Bacteria" "Proteobacteria" "Gammaproteobacteria" "Pasteurellales"
```

## 2. Manual: Load and prepare data with the help of utility functions

### 2.1 With metaphlan_to_phyloseq()

If the MetaPhlAn profiles are not loaded directly through
[`metaphlan_to_phyloseq()`](https://jrotzetter.github.io/metaphlanToPhyloseq/reference/metaphlan_to_phyloseq.md),
whether for editing names or subsetting columns or rows first,
[`metaphlan_to_phyloseq()`](https://jrotzetter.github.io/metaphlanToPhyloseq/reference/metaphlan_to_phyloseq.md)
can be used just as easily.

``` r

# Structure of a single MetaPhlAn microbiome profile
head(single_abundance_profile, 3)
#>                    clade_name NCBI_tax_id relative_abundance additional_species
#> 1                 k__Bacteria           2          100.00000                   
#> 2   k__Bacteria|p__Firmicutes      2|1239           74.02531                   
#> 3 k__Bacteria|p__Bacteroidota       2|976           25.97469

# Create a phyloseq object from a single pre-loaded profile...
physeq_single <- metaphlan_to_phyloseq(
    mtphlan_profile = single_abundance_profile,
    taxa_lvl = "genus",
    merged_profiles = FALSE
)

# ... or from multiple merged profiles
class(merged_abundance_profiles)
#> [1] "data.frame"

physeq_merged <- metaphlan_to_phyloseq(
  merged_abundance_profiles,
  taxa_lvl = "Species"
)
```

### 2.2 Manually without metaphlan_to_phyloseq()

Below is an example workflow for manually loaded MetaPhlAn profiles
without the
[`metaphlan_to_phyloseq()`](https://jrotzetter.github.io/metaphlanToPhyloseq/reference/metaphlan_to_phyloseq.md)
function and with examples for removing patterns from column names and
filtering low abundance entries.

``` r

merged_profiles <- merged_abundance_profiles

colnames(merged_profiles)
#> [1] "clade_name"                     "SRS014459.Stool"               
#> [3] "SRS014464.Anterior_nares"       "SRS014470.Tongue_dorsum"       
#> [5] "SRS014472.Buccal_mucosa"        "SRS014476.Supragingival_plaque"
#> [7] "SRS014494.Posterior_fornix"

# Remove the SRS0144... from the columnames, keeping only their sample type/origin
merged_profiles <- clean_colnames(merged_profiles, pattern = "SRS0144\\d{2}.")
colnames(merged_profiles)
#> [1] "clade_name"           "Stool"                "Anterior_nares"      
#> [4] "Tongue_dorsum"        "Buccal_mucosa"        "Supragingival_plaque"
#> [7] "Posterior_fornix"

# Filter to only include the chosen taxonomic rank
merged_profiles_species <- filter_taxa_lvl(merged_profiles, "Species")

# Merge all rows, where all values per sample are below the specified threshold, into 'Other' row 
merged_profiles_filtered <- filter_threshold(merged_profiles_species, threshold = 15)

# Create the taxonomic table
taxmat <- get_taxa_table(merged_profiles_filtered, taxa_lvl = "Species")

# Add the row names from the taxonomic table so they are the same in both
rownames(merged_profiles_filtered) <- rownames(taxmat)

# Remove clade_name column, otherwise phyloseq cannot create the otu_table due to it being character
merged_profiles_filtered$clade_name <- NULL

# Create the otu_table- and taxonomyTable-class objects
otutab <- phyloseq::otu_table(merged_profiles_filtered, taxa_are_rows = TRUE)
taxtab <- phyloseq::tax_table(taxmat)

# Create the phyloseq object
physeq <- phyloseq::phyloseq(otutab, taxtab)
```
