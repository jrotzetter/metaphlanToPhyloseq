# Changelog

## metaphlanToPhyloseq 0.2.1

#### Maintenance Status Update

- **Status Change**: The package is now in **maintenance-only mode**. It
  remains stable and fully functional for converting MetaPhlAn4 profiles
  to `phyloseq` objects.
- **Development Scope**: No new feature development is planned. Future
  updates will be restricted to critical bug fixes and compatibility
  adjustments.
- **Recommendation**: Users seeking more advanced capabilities are
  encouraged to explore the
  [microEDA](https://github.com/jrotzetter/microEDA) package.

## metaphlanToPhyloseq 0.2.0

- [`metaphlan_to_phyloseq()`](https://jrotzetter.github.io/metaphlanToPhyloseq/reference/metaphlan_to_phyloseq.md)
  now no longer requires data to be pre-filtered to a chosen taxonomic
  rank when not directly loading from a file.
- Taxonomic rank is now optional for
  [`metaphlan_to_phyloseq()`](https://jrotzetter.github.io/metaphlanToPhyloseq/reference/metaphlan_to_phyloseq.md)
  and deactivated by default, i.e. `taxa_lvl = NULL`, suppressing the
  filtering. Simply specify a valid taxonomic rank with `taxa_lvl` to
  filter to this rank.
- Added support for Species Genome Bin (SGB) as an additional taxonomic
  rank (set `taxa_lvl = 't'`).

## metaphlanToPhyloseq 0.1.0

- [`shorten_clade_names()`](https://jrotzetter.github.io/metaphlanToPhyloseq/reference/shorten_clade_names.md)
  now also works with data that wasn’t pre-filtered to a chosen
  taxonomic rank, returning the taxon without rank prefix for the chosen
  rank as before, and the last taxon with rank prefix for all others for
  clarity and as a reminder that data was not filtered to a specific
  taxonomic rank.
- The planned feature, which would have made it possible to load
  metadata directly from a path in
  [`metaphlan_to_phyloseq()`](https://jrotzetter.github.io/metaphlanToPhyloseq/reference/metaphlan_to_phyloseq.md),
  was discarded due to the large number of file formats in which
  metadata can be available.
- Fixed valid taxonomic rank check for cases in which ‘order’ was
  written out.
