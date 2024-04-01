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
#'  dataset based on a specified taxonomic level, using '|' as a separator.
#'
#' @param data The input dataset. Assumes clade names are of the structure
#'  ***k__Bacteria|p__Firmicutes|c__Bacilli|o__Lactobacillales|f__Lactobacillaceae***
#'  and found in the 'clade_name' column.
#' @param taxa_lvl The taxonomic level at which the clade names should be
#'  shortened. Valid options include 'kingdom', 'phylum', 'class', 'order',
#'   'family', 'genus', or 'species'. First letter abbreviations (e.g., 's')
#'    are also accepted.
#' @param apply_to_colnames `Logical` indicating whether the shortening should
#'  be applied to column names or row values. Default is TRUE.
#' @param selected_cols A `character` vector specifying the columns to which
#'  the shortening should be applied. If `NULL`, the shortening is applied to all
#'   columns. Default is `NULL`.
#'
#' @returns The dataset with the clade names shortened based on the specified
#'  taxonomic level.
#'
#' @export
#'
#' @details This function uses the dplyr package for data manipulation.
#'
#' @examples
#' head(species_only$clade_name)
#'
#' taxa_shortened <- shorten_clade_names(
#'   species_only,
#'   "species",
#'   apply_to_colnames = FALSE,
#'   selected_cols = "clade_name")
#'
#'   head(taxa_shortened$clade_name)
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
    "class", "family", "genus", "species"
  )
  if (!(taxa_lvl %in% valid_taxa_lvls)) {
    stop(paste0(
      "Invalid taxa_lvl. Please choose one of 'kingdom', 'phylum', ",
      "'class', 'order', 'family', 'genus', or 'species'."
    ))
  }

  first_letter <- substr(taxa_lvl, 1, 1)

  # Apply shortening to column names
  if (apply_to_colnames) {
    # Shorten all columns if none specified
    if (is.null(selected_cols)) {
      colnames(data) <- sub(
        paste0(".*\\|", first_letter, "__"), "", colnames(data)
      )
      # colnames(data) <- sub(paste0("(.*)\\|", first_letter, "__(.*)\\|.*"), "\\2", colnames(data))
    } else { # Otherwise only shorten selected cols
      colnames(data)[colnames(data) %in% selected_cols] <- sub(
        paste0(
          ".*\\|", first_letter, "__"
        ),
        "",
        colnames(data)[colnames(data) %in% selected_cols]
      )
    }
    # Apply shortening to row values
  } else {
    if (is.null(selected_cols)) {
      # Shorten all columns if none specified
      if (any(sapply(data, is.numeric))) {
        stop("There are numeric entries!")
        # Shorten non-numeric columns
      } else {
        data <- data |>
          dplyr::rowwise() |>
          dplyr::mutate(dplyr::across(dplyr::everything(), ~ sub(
            paste0(".*\\|", first_letter, "__"), "", .
          )))
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
          dplyr::mutate(dplyr::across(dplyr::all_of(selected_cols), ~ sub(
            paste0(".*\\|", first_letter, "__"), "", .
          )))
      }
    }
  }
  return(data)
}
