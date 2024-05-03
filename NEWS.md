# metaphlanToPhyloseq 0.2.0

* `metaphlan_to_phyloseq()` now no longer requires data to be pre-filtered to
a chosen taxonomic rank when not directly loading from a file.
* Taxonomic rank is now optional for `metaphlan_to_phyloseq()` and deactivated
by default, i.e. `taxa_lvl = NULL`, suppressing the filtering. Simply specify
a valid taxonomic rank with `taxa_lvl` to filter to this rank.
* Added support for Species Genome Bin (SGB) as an additional taxonomic rank
(set `taxa_lvl = 't'`).

# metaphlanToPhyloseq 0.1.0

* `shorten_clade_names()` now also works with data that wasn't pre-filtered to a
chosen taxonomic rank, returning the taxon without rank prefix for the chosen
rank as before, and the last taxon with rank prefix for all others for clarity
and as a reminder that data was not filtered to a specific taxonomic rank.
* The planned feature, which would have made it possible to load metadata
directly from a path in `metaphlan_to_phyloseq()`, was discarded due to the
large number of file formats in which metadata can be available.
* Fixed valid taxonomic rank check for cases in which 'order' was written out.
