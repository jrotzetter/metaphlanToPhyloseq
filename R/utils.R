#' Add OTUs as row names
#'
#' @description This function adds numbered Operational Taxonomic Units (OTUs)
#'  to the row names of a data frame. Optionally can also add a OTU column to
#'  the data frame.
#'
#' @param df The input MetaPhlAn data frame.
#' @param add_column `Logical` value indicating whether to add an OTU column
#'  to the data frame. Default is `FALSE`.
#'
#' @returns The data frame with OTU row names added.
#'
#' @export
#'
#' @examples
#' df <- data.frame(x = c(1, 2, 3), y = c(4, 5, 6))
#' add_otu(df) # Adds OTU row names to the data frame
#' add_otu(df, add_column = TRUE) # Adds OTU column and row names to the data frame
#'
#' @author Jérémy Rotzetter
add_otu <- function(df, add_column = FALSE) {
  otu_index <- paste0("Otu", 1:nrow(df))
  if (add_column) {
    df$Otu <- otu_index
  }
  rownames(df) <- otu_index
  return(df)
}


#' Set a column as row names
#'
#' @description This function sets a specified column as the row names of a
#'  data frame.
#'
#' @param df The input data frame.
#' @param column A `character` string. The name of the column to be set as
#'  row names.
#'
#' @returns The input data frame with the specified column set as row names.
#'
#' @importFrom utils hasName
#' @export
#'
#' @examples
#' df <- data.frame(
#'   clade_name = c("A", "B", "C"),
#'   col1 = c(10, 20, 5),
#'   col2 = c(8, 15, 3),
#'   col3 = c(64, 89, 32)
#' )
#' col_to_rownames(df, "clade_name")
#'
col_to_rownames <- function(df, column) {
  if (!hasName(df, column)) {
    stop("The selected column does not exist!")
  }
  df <- as.data.frame(df)
  rownames(df) <- df[[column]]
  df[[column]] <- NULL
  return(df)
}


#' Modify the case of the first letter in each word of a character string
#'
#' @description These functions take a character string as input and return a
#'  new character string with the first letter of each word modified to either
#'  uppercase or lowercase.
#'
#' @param character_string A `character` vector containing the strings to be
#'  modified.
#'
#' @returns A character vector with the first letter of each word modified to
#'  either uppercase or lowercase.
#'
#' @examples
#' capitalize_str("hello world") # Returns "Hello World"
#' capitalize_str(c("hello", "world")) # Returns c("Hello", "World")
#' lowercase_str("Hello World") # Returns "hello world"
#' lowercase_str(c("Hello", "World")) # Returns c("hello", "world")
#'
#' @export
#' @rdname first_letter_case
capitalize_str <- function(character_string) {
  sapply(character_string, function(x) {
    paste0(toupper(substring(x, 1, 1)), substring(x, 2))
  }, USE.NAMES = FALSE)
}

#' @export
#' @rdname first_letter_case
lowercase_str <- function(character_string) {
  sapply(character_string, function(x) {
    paste0(tolower(substring(x, 1, 1)), substring(x, 2))
  }, USE.NAMES = FALSE)
}


#' Get full taxonomic rank name
#'
#' @description Function to get the full taxonomic rank name based on the
#'  abbreviated taxonomic level.
#'
#' @param taxa_lvl A `character` string with the shortened taxonomic level
#'  (e.g., "p", "c", ..., "g", "s").
#'
#' @returns A `character` vector of the extended taxonomic level.
#'
#' @author Jérémy Rotzetter
#'
#' @noRd
get_full_name <- function(taxa_lvl) {
  stopifnot(is.character(taxa_lvl))
  # Dictionary mapping abbreviated taxonomic levels to their full names
  taxa_dict <- c(
    "k" = "Kingdom",
    "p" = "Phylum",
    "c" = "Class",
    "o" = "Order",
    "f" = "Family",
    "g" = "Genus",
    "s" = "Species"
  )
  # Check if the taxonomic level is in the dictionary
  if (taxa_lvl %in% names(taxa_dict)) {
    return(unname(taxa_dict[taxa_lvl])) # Return the full name
  } else {
    return(taxa_lvl) # Return the input taxonomic level if not found in the dictionary
  }
}


#' Get numbered sample columns
#'
#' @description This function retrieves the column names from a data frame that
#' contain a number anywhere in the name. It is assumed only sample names in the
#' MetaPhlAn profile contain a number.
#'
#' @param df The input data frame. Intended to be a MetaPhlAn profile.
#'
#' @returns A `character` vector of column names containing a number.
#'
#' @examples
#' df <- data.frame(
#'   a = c("Taxon 1", "Taxon 2", "Taxon 3"),
#'   Sample_1 = 1:3,
#'   Sample_2 = 4:6,
#'   Sample_3 = 7:9
#' )
#' get_sample_cols(df)
#'
#' @author Jérémy Rotzetter
#'
#' @noRd
get_sample_cols <- function(df) {
  r <- ".*[0-9].*" # match column names that contain a number anywhere
  sample_cols <- names(df)[grepl(r, names(df))]
  return(sample_cols)
}


#' Check the validity of the selected taxonomic level
#'
#' @description This function compares the first letter of the provided
#'  taxa_lvl to the last matched taxa rank in mtphlan_profile, to check if they
#'   are the same.
#'
#' @param mtphlan_profile The MetaPhlAn profile data containing clade names in
#'  the 'clade_name' column.
#' @param taxa_lvl A `character` string. The taxa level to be checked. Valid
#'  options include 'kingdom', 'phylum', 'class', 'order', 'family', 'genus',
#'  'species' or 't' (SGB). First letter abbreviations (e.g., 's') are also accepted.
#'
#' @returns This function does not return any value. It stops the execution with
#'  an error message if the taxa level is missing or invalid.
#'
#' @author Jérémy Rotzetter
#'
#' @noRd
check_taxa_lvl <- function(mtphlan_profile, taxa_lvl) {
  # Check if taxa_lvl is missing
  if (missing(taxa_lvl)) {
    stop(paste0(
      "The 'taxa_lvl' parameter is missing. Please choose one of ",
      "'kingdom', 'phylum', 'class', 'order', 'family', 'genus' or 'species'."
    ))
  }

  taxa_lvl <- lowercase_str(taxa_lvl)

  # Check if taxa_lvl is valid
  valid_taxa_lvls <- c(
    "k", "p", "c", "o", "f", "g", "s", "t", "kingdom", "phylum",
    "class", "order", "family", "genus", "species"
  )
  if (!(taxa_lvl %in% valid_taxa_lvls)) {
    stop(paste0(
      "Invalid taxa_lvl. Please choose one of 'kingdom', 'phylum',",
      "'class', 'order', 'family', 'genus', 'species', or 't (SGB)'."
    ))
  }

  first_letter <- substr(taxa_lvl, 1, 1)
  first_letter <- lowercase_str(first_letter)

  matches <- regmatches(
    mtphlan_profile$clade_name,
    gregexpr("([kpcofgst])(?=__)",
      mtphlan_profile$clade_name,
      perl = TRUE
    )
  )
  letter <- utils::tail(matches[[1]], 1)

  all(first_letter == letter)

  # # stringr 1.5.1 solution
  # # Extract the last occurrence of the pattern [kpcofgs]__
  # last_occurrence <- stringr::str_extract(clade_names, "(?<=\\|)[kpcofgs]__(?!.*\\|)")
  #
  # # Remove the "__" from the extracted strings
  # last_occurrence <- stringr::str_remove(last_occurrence, "__")
  #
  # all(first_letter == last_occurrence)
}
