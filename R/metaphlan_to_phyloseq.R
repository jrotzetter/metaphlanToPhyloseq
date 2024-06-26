#' Convert MetaPhlAn profile to phyloseq object
#'
#' @description This function converts a MetaPhlAn profile to a phyloseq object
#'  when given the file path or a pre-loaded table.
#'
#' @param mtphlan_profile The MetaPhlAn profile to be converted. It can be
#'  either a file path or a data frame of MetaPhlAn profile(s).
#' @param taxa_lvl Optional taxonomic level to filter the profile to. Valid
#'  options are 'kingdom', 'phylum', 'class', 'order', 'family', 'genus',
#'  'species' or 't' (SGB). First letter abbreviations are also accepted.
#' @param metadata Optional metadata for the samples. If provided, it should be
#'  a data frame.
#' @param sample_column `Character` string. The column in the metadata
#'  containing the sample names. Should match column names of the MetaPhlAn
#'  profile.
#' @param use_taxa_names `Logical` indicating whether to use taxonomic names
#'  instead of OTUs in the resulting phyloseq object. Default is `FALSE`.
#' @param merged_profiles `Logical`; if `TRUE` (the default) the file to be
#' loaded is assumed to be multiple merged MetaPhlAn profiles.
#'
#' @returns A phyloseq object representing the MetaPhlAn profile.
#'
#' @note
#' If `mtphlan_profile` is an already loaded object, the clade names should not
#' be shortened, otherwise the taxonomic table cannot be created. Should you
#' wish to use the taxonomic names of the specified rank as the row names in
#' the phyloseq object, set `use_taxa_names = TRUE`.
#'
#' @export
#'
#' @examples
#' metaphlan_to_phyloseq(species_only, taxa_lvl = "s")
#'
#' @author Jérémy Rotzetter
metaphlan_to_phyloseq <- function(
    mtphlan_profile,
    taxa_lvl = NULL,
    metadata = NULL,
    sample_column = NULL,
    use_taxa_names = FALSE,
    merged_profiles = TRUE) {
  stopifnot(is.null(taxa_lvl) | is.character(taxa_lvl))

  # Load metaphlan profile from file
  if (is.character(mtphlan_profile)) {
    if (merged_profiles) {
      mtphlan_profile <- load_metaphlan_profile(mtphlan_profile)
    } else {
      mtphlan_profile <- load_metaphlan_profile(
        mtphlan_profile,
        merged_profiles = FALSE
      )
      mtphlan_profile <- mtphlan_profile[, c(1, 3)]
    }
  }

  if (merged_profiles && any(c("NCBI_tax_id", "additional_species")
  %in% names(mtphlan_profile))) {
    stop(paste(
      "Loaded file seems to be a single MetaPhlAn profile!",
      "Please set merged_profiles to FALSE."
    ))
  }

  if (!merged_profiles) {
    if (!"relative_abundance" %in% names(mtphlan_profile)) {
      stop(paste(
        "Loaded profile does not appear to be a single MetaPhlAn profile!",
        "Please set merged_profiles to TRUE."
      ))
    }

    index <- which(names(mtphlan_profile) %in% c("clade_name", "relative_abundance"))
    mtphlan_profile <- mtphlan_profile[index]
  }

  if (is.character(taxa_lvl)) {
    # Convert taxa_lvl to lowercase
    taxa_lvl <- lowercase_str(taxa_lvl)

    # Check if taxa_lvl is valid
    valid_taxa_lvls <- c(
      "k", "p", "c", "o", "f", "g", "s", "t", "kingdom", "phylum",
      "class", "order", "family", "genus", "species"
    )
    if (!(taxa_lvl %in% valid_taxa_lvls)) {
      stop(paste(
        "Invalid taxa_lvl. Please choose one of 'kingdom', 'phylum',",
        "'class', 'order', 'family', 'genus', 'species', or 't (SGB)'."
      ))
    }

    # # Check if selected taxon level matches metaphlan profile
    # if (!check_taxa_lvl(mtphlan_profile, taxa_lvl)) {
    #   stop("Selected taxon level does not match with filtered MetaPhlAn profile!")
    # }

    # Filter to selected taxon level
    mtphlan_profile <- filter_taxa_lvl(mtphlan_profile, taxa_lvl)

    if (nrow(mtphlan_profile) == 0) {
      stop("Selected taxon level does not match with already filtered MetaPhlAn profile!")
    }
  }

  # Process metadata if provided
  if (!is.null(metadata) && merged_profiles == TRUE) {
    # Get the index of the sample_column
    sample_column_index <- which(names(metadata) == sample_column)
    # Replace spaces with "_" in colnames
    colnames(metadata) <- gsub(
      x = colnames(metadata),
      pattern = " ",
      replacement = "_"
    )
    sample_column <- names(metadata[sample_column_index])

    # Keep only sample names shared between metadata and metaphlan profile
    shared_names <- base::intersect(
      colnames(mtphlan_profile),
      metadata[[sample_column]]
    )
    mtphlan_profile_cleaned <- mtphlan_profile[, c("clade_name", shared_names)]
    metadata_cleaned <- metadata[match(shared_names, metadata[[sample_column]]), ]

    # Check if all sample names match between metadata and profile
    message(
      ifelse(
        all(
          metadata_cleaned[[sample_column]] %in% names(
            mtphlan_profile_cleaned
          )[-which(names(mtphlan_profile_cleaned) == "clade_name")]
        ),
        "Metadata and profile Sample_names match.",
        "Check Sample_names in 'metadata' and 'mtphlan_profile' for chosen taxa!"
      )
    )

    metadata_cleaned <- as.data.frame(metadata_cleaned)
    rownames(metadata_cleaned) <- metadata_cleaned[[sample_column]]
  } else {
    # If no metadata provided, just use metaphlan profile as is
    mtphlan_profile_cleaned <- mtphlan_profile
    metadata_cleaned <- metadata
  }
  # Get taxa table from profile
  taxa_table <- get_taxa_table(
    mtphlan_profile_cleaned,
    taxa_lvl = taxa_lvl,
    taxa_are_rows = FALSE,
    use_taxa_names = use_taxa_names
  )
  rownames(mtphlan_profile_cleaned) <- rownames(taxa_table)
  # Remove clade_name column, otherwise phyloseq cannot create the otu_table
  mtphlan_profile_cleaned$clade_name <- NULL

  # Create phyloseq object
  physeq <- phyloseq::phyloseq(
    phyloseq::otu_table(as.matrix(mtphlan_profile_cleaned), taxa_are_rows = TRUE),
    phyloseq::tax_table(taxa_table),
    phyloseq::sample_data(metadata_cleaned, errorIfNULL = FALSE)
  )

  return(physeq)
}
