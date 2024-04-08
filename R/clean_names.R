#' Remove pattern from column names
#'
#' @description Function to remove a pattern from column names.
#'
#' @param df A data frame or matrix.
#' @param pattern A `character` string with the pattern to match.
#'
#' @returns Returns the input table with the specified pattern removed from the
#'  column names.
#'
#' @export
#'
#' @examples
#' colnames(merged_abundance_profiles)
#'
#' names_cleaned <- clean_colnames(merged_abundance_profiles, "SRS0144\\d{2}.")
#'
#' colnames(names_cleaned)
#'
#' @author Jérémy Rotzetter
clean_colnames <- function(df, pattern) {
  # Check that df is a data frame or matrix
  stopifnot(is.data.frame(df) | is.matrix(df))
  # Replace the pattern in the column names with blank
  names(df) <- gsub(pattern, "", names(df))
  return(df)
}


#' Only keep shared samples between metadata and MetaPhlAn profile
#'
#' @description The `clean_metadata()` function is designed to clean and prepare
#' metadata and a MetaPhlAn profile for further analysis by only keeping shared
#' samples.
#'
#' @param metadata A data frame containing metadata information.
#' @param mtphlan_profile A data frame containing metaphlan profile information.
#'  Samples are columns.
#' @param sample_column A `character` string. The column in the metadata
#'  containing the sample names. Should match column names of the MetaPhlAn
#'   profile.
#' @param remove_spaces Should spaces (" ") be replaced with "_" in the column
#'  names?
#'
#' @returns The cleaned metadata as a data frame.
#'
#' @export
#'
#' @examples
#' metadata <- data.frame(
#'   Sample_name = c("A", "B", "C", "D", "F", "G"),
#'   Sex = c(rep(c("F", "M")))
#' )
#'
#' mtphlan_profile <- data.frame(
#'   clade_name = c("x", "y", "z"),
#'   A = c(1, 2, 3),
#'   B = c(4, 5, 6),
#'   D = c(7, 8, 9),
#'   Z = c(10, 11, 12)
#' )
#'
#' cleaned_metadata <- clean_metadata(metadata, mtphlan_profile, "Sample_name")
#'
#' print(cleaned_metadata)
#'
#' @author Jérémy Rotzetter
clean_metadata <- function(
    metadata,
    mtphlan_profile,
    sample_column,
    remove_spaces = TRUE) {
  if (remove_spaces) {
    # Get the index of the sample_column
    sample_column_index <- which(names(metadata) == sample_column)
    # Replace spaces with "_" in colnames
    colnames(metadata) <- gsub(
      x = colnames(metadata),
      pattern = " ",
      replacement = "_"
    )

    sample_column <- names(metadata[sample_column_index])
  }
  # Get shared column names between mtphlan_profile and metadata
  shared_names <- base::intersect(
    colnames(mtphlan_profile),
    metadata[[sample_column]]
  )
  # Create a new metaphlan_profile with only the "clade_name" column and the
  # shared column names
  mtphlan_profile_cleaned <- mtphlan_profile[, c("clade_name", shared_names)]
  # Create a new metadata with only the rows that have shared column names
  metadata_cleaned <- metadata[match(shared_names, metadata[[sample_column]]), ]

  # Check if all sample names are present in both metadata and mtphlan_profile
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

  return(metadata_cleaned)
}


#' Shorten clade names in a dataset to chosen taxonomic level
#'
#' @description This function shortens the taxonomic names of clades in a given
#'  dataset based on a specified taxonomic level, using the first letter of the
#'  taxonomic rank + "__" as rank identifiers.
#'
#' @param data The input dataset. Assumes clade names are of the structure
#'  ***k__Bacteria|p__Firmicutes|c__Bacilli|o__Lactobacillales|f__Lactobacillaceae***
#'  and found in the 'clade_name' column or the column names themselves.
#' @param taxa_lvl The taxonomic level at which the clade names should be
#'  shortened. Valid options include 'kingdom', 'phylum', 'class', 'order',
#'   'family', 'genus', or 'species'. First letter abbreviations (e.g., 's')
#'    are also accepted.
#' @param apply_to_colnames `Logical` indicating whether the shortening should
#'  be applied to column names or row values. Default is `TRUE`.
#' @param selected_cols A `character` vector specifying the columns to which
#'  the shortening should be applied. If `NULL` (the default), the shortening
#'  is applied to all columns.
#'
#' @returns The dataset with the clade names shortened based on the specified
#'  taxonomic level.
#'
#'  In the case where there are entries not matching the chosen taxonomic rank,
#'  these are either returned 'as is', or if they follow the same structure,
#'  the name will be shortened to the last taxonomic entry
#'  (see rows 2, 5 and 6 of the example)
#'
#' @export
#'
#' @note
#' This function is not intended to be used with the workflow for the creation
#' of phyloseq objects as the full sequence of taxonomic names is needed for
#' the creation of the taxonomy table in [get_taxa_table()]. It may however be
#' useful for analyses or plots created directly with/from the dataframes.
#'
#' This function uses the dplyr package for data manipulation.
#'
#' @examples
#' head(merged_abundance_profiles$clade_name)
#'
#' taxa_shortened <- shorten_clade_names(
#'   merged_abundance_profiles,
#'   "Phylum",
#'   apply_to_colnames = FALSE,
#'   selected_cols = "clade_name"
#' )
#'
#' head(taxa_shortened$clade_name)
#'
#' @author Jérémy Rotzetter
shorten_clade_names <- function(
    data,
    taxa_lvl,
    apply_to_colnames = TRUE,
    selected_cols = NULL) {
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop(
      "Package \"dplyr\" must be installed to use this function.",
      call. = FALSE
    )
  }

  # Check if taxa_lvl is missing
  if (missing(taxa_lvl)) {
    stop(
      paste0(
        "The 'taxa_lvl' parameter is missing. Please choose one of ",
        "'kingdom', 'phylum', 'class', 'order', 'family', 'genus' or 'species'."
      )
    )
  }

  taxa_lvl <- lowercase_str(taxa_lvl)

  # Check if taxa_lvl is valid
  valid_taxa_lvls <- c(
    "k", "p", "c", "o", "f", "g", "s", "kingdom", "phylum",
    "class", "order", "family", "genus", "species"
  )
  if (!(taxa_lvl %in% valid_taxa_lvls)) {
    stop(paste0(
      "Invalid taxa_lvl. Please choose one of 'kingdom', 'phylum', ",
      "'class', 'order', 'family', 'genus', or 'species'."
    ))
  }

  first_letter <- substr(taxa_lvl, 1, 1)

  # Function to get the last taxonomic entry
  extract_name <- function(name, first_letter) {
    pattern <- "\\b[kpcofgs]__\\w+\\b$"
    matches <- regmatches(name, gregexpr(pattern, name))[[1]]
    last_match <- matches[length(matches)]

    if (length(last_match) > 0) {
      if (grepl(paste0("\\b", first_letter, "__\\w+\\b$"), last_match)) {
        # If the first_letter pattern is found, extract the short name
        short_name <- sub(
          paste0("\\b", first_letter, "__(\\w+\\b$)"),
          "\\1",
          last_match
        )
      } else { # If the first_letter pattern is not found, use the entire last
        # match as the short name
        short_name <- last_match
      }
    } else { # If there are no matches, use original 'name' as the short name
      short_name <- name
    }
    return(short_name)
  }

  # Apply shortening to column names
  if (apply_to_colnames) {
    # Shorten all columns if none specified
    if (is.null(selected_cols)) {
      colnames(data) <- sapply(
        colnames(data),
        extract_name,
        first_letter = first_letter
      )
    } else { # Otherwise only shorten selected cols
      colnames(data)[colnames(data) %in% selected_cols] <- sapply(
        colnames(data)[colnames(data) %in% selected_cols],
        extract_name,
        first_letter = first_letter
      )
    }
    # Apply shortening to row values
  } else {
    if (is.null(selected_cols)) {
      # Shorten all columns if none specified
      if (any(sapply(data, is.numeric))) {
        stop("There are numeric entries! Please select only non-numeric columns.")
        # Shorten non-numeric columns
      } else {
        data <- data |>
          dplyr::rowwise() |>
          dplyr::mutate(dplyr::across(
            dplyr::everything(),
            ~ extract_name(., first_letter = first_letter)
          ))
      }
      # Otherwise only shorten specified columns
    } else {
      # Check selected columns for numerics
      if (any(sapply(data[, selected_cols], is.numeric))) {
        stop("One or more selected column(s) is numeric!")
        # Shorten non-numeric columns
      } else {
        data <- data |>
          dplyr::rowwise() |>
          dplyr::mutate(dplyr::across(
            dplyr::all_of(selected_cols),
            ~ extract_name(., first_letter = first_letter)
          ))
      }
    }
  }
  return(data)
}
