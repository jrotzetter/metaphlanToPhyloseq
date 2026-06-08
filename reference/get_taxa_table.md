# Get taxonomic table

Function to extract taxonomic information from a data frame and prepare
it for phyloseq.

## Usage

``` r
get_taxa_table(
  df,
  taxa_lvl = NULL,
  taxa_are_rows = FALSE,
  use_taxa_names = FALSE
)
```

## Arguments

- df:

  The input data frame.

- taxa_lvl:

  The taxonomic level to extract ('kingdom', 'phylum', 'class', 'order',
  'family', 'genus', 'species' or 't' (SGB)). First letter abbreviations
  (e.g., 's') are also accepted. Is deactivated with `NULL` by default.

- taxa_are_rows:

  `Logical`; if `TRUE`, taxa information is in the row names of the data
  frame, otherwise it is assumed to be in the 'clade_name' column.

- use_taxa_names:

  `Logical`; if `TRUE`, assign taxa names as row names, otherwise assign
  Otu as row names.

## Value

A taxonomic matrix.

## Note

Will work with unfiltered as well as filtered to a specific taxonomic
level input data frames.

## Author

Jérémy Rotzetter

## Examples

``` r
taxmat <- get_taxa_table(species_only, taxa_lvl = "s")
head(taxmat)
#>      Kingdom    Phylum           Class                 Order              
#> Otu1 "Bacteria" "Proteobacteria" "Gammaproteobacteria" "Moraxellales"     
#> Otu2 "Bacteria" "Actinobacteria" "Actinomycetia"       "Corynebacteriales"
#> Otu3 "Bacteria" "Firmicutes"     "Negativicutes"       "Veillonellales"   
#> Otu4 "Bacteria" "Bacteroidota"   "Bacteroidia"         "Bacteroidales"    
#> Otu5 "Bacteria" "Firmicutes"     "Bacilli"             "Lactobacillales"  
#> Otu6 "Bacteria" "Proteobacteria" "Gammaproteobacteria" "Pasteurellales"   
#>      Family               Genus            
#> Otu1 "Moraxellaceae"      "Moraxella"      
#> Otu2 "Corynebacteriaceae" "Corynebacterium"
#> Otu3 "Veillonellaceae"    "Veillonella"    
#> Otu4 "Prevotellaceae"     "Prevotella"     
#> Otu5 "Streptococcaceae"   "Streptococcus"  
#> Otu6 "Pasteurellaceae"    "Haemophilus"    
#>      Species                               
#> Otu1 "Moraxella_nonliquefaciens"           
#> Otu2 "Corynebacterium_pseudodiphtheriticum"
#> Otu3 "Veillonella_dispar"                  
#> Otu4 "Prevotella_histicola"                
#> Otu5 "Streptococcus_mitis"                 
#> Otu6 "Haemophilus_haemolyticus"            
```
