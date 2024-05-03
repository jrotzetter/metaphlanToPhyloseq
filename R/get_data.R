#' Load MetaPhlAn profile(s)
#'
#' @description Function to load a MetaPhlAn profile from a file path.
#'
#' @param path A `character` string containing the path to the file to read.
#' @param merged_profiles `Logical`; if `TRUE` (the default) the file to be
#' loaded is assumed to be multiple merged MetaPhlAn profiles.
#'
#' @returns A data frame containing the MetaPhlAn results.
#'
#' @details
#' If a single MetaPhlAn profile is to be loaded, `merged_profiles` should be
#' set to `FALSE`. THe single profile is assumed to contain the column names
#' `clade_name`, `NCBI_tax_id`, `relative_abundance` and `additional_species`.
#'
#'
#' @importFrom utils read.table
#' @export
#'
#' @examples
#' path <- list_example("merged_abundance_table.txt")
#' data <- load_metaphlan_profile(path = path)
#' head(data)
#'
#' @author Jérémy Rotzetter
load_metaphlan_profile <- function(path, merged_profiles = TRUE) {
  stopifnot(is.character(path))

  if (merged_profiles) {
    data <- read.table(path, header = TRUE, sep = "\t")
    if (!any(colnames(data) == "clade_name")) {
      stop(paste(
        "Loaded file seems to be missing the \"clade_name\" column!"
      ))
    }
    if (any(c("NCBI_tax_id", "additional_species") %in% names(data))) {
      stop(paste(
        "Loaded file seems to be a single MetaPhlAn profile!",
        "Please set merged_profiles to FALSE."
      ))
    }
  } else {
    data <- read.table(path, header = FALSE, sep = "\t")
    if (length(data) != 4) {
      stop(paste0(
        "Loaded file has wrong number of columns to be a single ",
        "MetaPhlAn profile!\n  Only the columns \"clade_name\", \"NCBI_tax_id\", ",
        "\"relative_abundance\", \"additional_species\" should be present."
      ))
    }
    colnames(data) <- c(
      "clade_name",
      "NCBI_tax_id",
      "relative_abundance",
      "additional_species"
    )
  }
  return(data)
}


#' Get taxonomic table
#'
#' @description Function to extract taxonomic information from a data frame and
#'  prepare it for phyloseq.
#'
#' @param df The input data frame.
#' @param taxa_lvl The taxonomic level to extract ('kingdom', 'phylum', 'class',
#'  'order', 'family', 'genus', 'species' or 't' (SGB)). First letter abbreviations
#'   (e.g., 's') are also accepted. Is deactivated with `NULL` by default.
#' @param taxa_are_rows `Logical`; if `TRUE`, taxa information is in the row
#'  names of the data frame, otherwise it is assumed to be in the 'clade_name'
#'  column.
#' @param use_taxa_names `Logical`; if `TRUE`, assign taxa names as row names,
#'  otherwise assign Otu as row names.
#'
#' @returns A taxonomic matrix.
#'
#' @note Will work with unfiltered as well as filtered to a specific taxonomic
#' level input data frames.
#'
#' @importFrom utils hasName
#' @export
#'
#' @examples
#' taxmat <- get_taxa_table(species_only, taxa_lvl = "s")
#' head(taxmat)
#'
#' @author Jérémy Rotzetter
get_taxa_table <- function(
    df,
    taxa_lvl = NULL,
    taxa_are_rows = FALSE,
    use_taxa_names = FALSE) {
  # Function to trim taxon names
  trim_taxa_names <- function(x) {
    match <- gsub("^[kpcofgst]__", "", as.character(x))
    return(match)
  }

  stopifnot(is.null(taxa_lvl) | is.character(taxa_lvl))

  taxa_cols <- c(
    "Kingdom",
    "Phylum",
    "Class",
    "Order",
    "Family",
    "Genus",
    "Species",
    "SGB"
  )

  # Check if 'clade_name' column exists when taxa_are_rows is FALSE
  if (!taxa_are_rows & !hasName(df, "clade_name")) {
    stop("The 'clade_name' column does not exist!")
  }

  # Check if rownames are numbers
  if (taxa_are_rows &&
    (all(grepl("^\\d+$", rownames(df))))) {
    warning(paste0(
      "Rownames are numbers!'taxa_are_rows' set to FALSE."
    ))
    taxa_are_rows <- FALSE
  }

  # Split clade_name or rownames by pipe to separate taxon levels
  if (taxa_are_rows) {
    split_parts <- strsplit(rownames(df), "\\|")
  } else {
    # Split the clade_name column
    split_parts <- strsplit(df$clade_name, "\\|")
  }

  # Determine the maximum number of parts
  max_parts <- max(lengths(split_parts))

  # Create a new dataframe with NA values
  df_taxa <- data.frame(matrix(NA, nrow = length(split_parts), ncol = max_parts))

  # Assign the split parts to the new dataframe
  for (i in 1:length(split_parts)) {
    df_taxa[i, 1:length(split_parts[[i]])] <- split_parts[[i]]
  }
  names(df_taxa) <- taxa_cols[1:ncol(df_taxa)]

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

    if (any(df_taxa[1] == "Other")) {
      df_taxa[df_taxa[1] == "Other", ] <- "Other"
    }

    # Check if all rows have the same number of non-NA entries, to verify if
    # metaphlan profile was prefiltered to specifc taxonomic rank or not
    if (!all(rowSums(!is.na(df_taxa)) == rowSums(!is.na(df_taxa))[1])) {
      stop("Please filter the metaphlan profile to a specific taxonomic rank first.")
    }

    # Check if selected taxon level matches metaphlan profile
    if (!check_taxa_lvl(df, taxa_lvl)) {
      stop("Selected taxon level does not match with filtered metaphlan profile!")
    }

    # Trim taxon names
    for (col in names(df_taxa)) {
      df_taxa[[col]] <- trim_taxa_names(df_taxa[[col]])
    }
  }

  # Set new rownames
  if (use_taxa_names) {
    # rownames(df_taxa) <- df_taxa[, ncol(df_taxa)]
    rownames(df_taxa) <- apply(df_taxa, 1, function(x) utils::tail(stats::na.omit(x), 1))
    # Trim taxon names
    for (col in names(df_taxa)) {
      df_taxa[[col]] <- trim_taxa_names(df_taxa[[col]])
    }
  } else {
    # Trim taxon names
    for (col in names(df_taxa)) {
      df_taxa[[col]] <- trim_taxa_names(df_taxa[[col]])
    }
    otu_index <- paste0("Otu", 1:nrow(df))

    rownames(df_taxa) <- otu_index
  }

  return(as.matrix(df_taxa))
}


#' List example files
#'
#' This function lists available example files, found in the `inst/extdata`
#'  directory.
#'
#' @param path The name of an example file. If not provided, the function lists
#'  the available example files.
#'
#' @returns A character vector.
#'
#' @examples
#' list_example() # Lists example files contained in the package
#' list_example("merged_abundance_table.txt") # Lists the path of the specified file
#'
#' @export
#'
#' @source This function was adapted from [readxl::readxl_example()] and
#'  originally created by *Hadley Wickham & Jennifer Bryan*.
list_example <- function(path = NULL) {
  if (is.null(path)) {
    dir(system.file("extdata", package = "metaphlanToPhyloseq"))
  } else {
    system.file("extdata", path, package = "metaphlanToPhyloseq", mustWork = TRUE)
  }
}
