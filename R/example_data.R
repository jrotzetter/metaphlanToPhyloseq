#' SRS014470-Tongue_dorsum_profile data
#'
#' Example data of a single MetaPhlAn profile output file, containing the
#'  computed taxon abundances.
#'
#' @format ## `single_abundance_profile`
#' A data frame with 15 rows and 4 columns:
#' \describe{
#'   \item{clade_name}{The taxonomic lineage of the taxon reported on this line.
#'   Taxon names are prefixed with one-letter codes to help indicate their rank.}
#'   \item{NCBI_tax_id}{The NCBI-equivalent taxon IDs of the named taxa from
#'   clade_name.}
#'   \item{relative_abundance}{The taxon's relative abundance in %.}
#'   \item{additional_species}{Additional species names for cases where the
#'   metagenome profile contains clades that represent multiple species. The
#'   species listed in column 1 is the representative species in such cases.}
#' }
#' @source <https://github.com/biobakery/biobakery/wiki/MetaPhlAn-4.1>
"single_abundance_profile"


#' Merged MetaPhlAn profiles data
#'
#' Example data of multiple MetaPhlAn profiles merged into one file.
#'
#' @format ## `merged_abundance_profiles`
#' A data frame with 59 rows and 7 columns:
#' \describe{
#'   \item{clade_name}{The taxonomic lineage of the taxon reported on this line.
#'   Taxon names are prefixed with one-letter codes to help indicate their rank.}
#'   \item{SRS0144...}{The sample names for which sample-specific abundance
#'    profiles are found in a column.}
#' }
#' @source <https://github.com/biobakery/biobakery/wiki/MetaPhlAn-4.1>
"merged_abundance_profiles"

#' Merged MetaPhlAn profiles data to species level
#'
#' Example data of multiple MetaPhlAn profiles which were merged into one file,
#' and then filtered with `filter_taxa_lvl(df, taxa_lvl = "s")` to only include
#' species-level.
#'
#' @format ## `species_only`
#' A data frame with 12 rows and 7 columns:
#' \describe{
#'   \item{clade_name}{The taxonomic lineage of the taxon reported on this line.
#'   Taxon names are prefixed with one-letter codes to help indicate their rank.}
#'   \item{SRS0144...}{The sample names for which sample-specific abundance
#'    profiles are found in a column.}
#' }
#' @source <https://github.com/biobakery/biobakery/wiki/MetaPhlAn-4.1>
"species_only"
