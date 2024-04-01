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
    data <- read.table(path, skip = 1, header = TRUE, sep = "\t")
    if (!any(colnames(data) == "clade_name")) {
      stop(paste0(
        "Loaded file seems to not be merged MetaPhlAn profiles! ",
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
#'  'order', 'family', 'genus', or 'species'). First letter abbreviations
#'   (e.g., 's') are also accepted.
#' @param taxa_are_rows `Logical`; if `TRUE`, taxa information is in the row
#'  names of the data frame, otherwise it is assumed to be in the 'clade_name'
#'  column.
#' @param use_taxa_names `Logical`; if `TRUE`, assign taxa names as row names,
#'  otherwise assign Otu as row names.
#'
#' @returns A taxonomic matrix.
#'
#' @note Assumes the input data frame has already been filtered to a specific
#' taxonomic level.
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
    taxa_lvl,
    taxa_are_rows = FALSE,
    use_taxa_names = FALSE) {
  # Function to trim taxon names
  trim_taxa_names <- function(x) {
    match <- gsub("^[kpcofgs]__", "", as.character(x))
    return(match)
  }
  # Check if 'clade_name' column exists when taxa_are_rows is FALSE
  if (!taxa_are_rows & !hasName(df, "clade_name")) {
    stop("The 'clade_name' column does not exist!")
  }

  # Check if rownames are numbers or missing the taxonomic level separator
  if (taxa_are_rows &&
    (all(grepl("^\\d+$", rownames(df))) ||
      !any(grepl("\\|", rownames(df))))) {
    warning(paste0(
      "Rownames are numbers, or do not contain the '|' separator. ",
      "'taxa_are_rows' set to FALSE."
    ))
    taxa_are_rows <- FALSE
  }

  # Check if taxa_lvl is missing
  if (missing(taxa_lvl)) {
    stop(paste0(
      "The 'taxa_lvl' parameter is missing. Please choose one of ",
      "'kingdom', 'phylum', 'class', 'order', 'family', 'genus' or 'species'."
    ))
  }
  # Convert taxa_lvl to lowercase
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
  # Check if selected taxon level matches metaphlan profile
  if (!check_taxa_lvl(df, taxa_lvl)) {
    stop("Selected taxon level does not match metaphlan profile!")
  }
  # Get full name and capitalize taxa_lvl
  taxa_lvl <- get_full_name(taxa_lvl)
  taxa_lvl <- capitalize_str(taxa_lvl)

  # Split clade_name or rownames by pipe to separate taxon levels
  if (taxa_are_rows) {
    df_taxa <- as.data.frame(strsplit(rownames(df), "|", fixed = TRUE))
  } else {
    df_taxa <- as.data.frame(strsplit(df$clade_name, "|", fixed = TRUE))
  }

  df_taxa <- as.data.frame(t(df_taxa))
  # Set column names and value to extract
  taxa_cols <- c(
    "Kingdom",
    "Phylum",
    "Class",
    "Order",
    "Family",
    "Genus",
    "Species"
  )
  taxa_dict <- c(
    "Kingdom" = 1,
    "Phylum" = 2,
    "Class" = 3,
    "Order" = 4,
    "Family" = 5,
    "Genus" = 6,
    "Species" = 7
  )
  value <- taxa_dict[[taxa_lvl]]
  taxa_cols <- taxa_cols[1:value]
  names(df_taxa) <- taxa_cols
  # Trim taxon names
  for (col in names(df_taxa)) {
    df_taxa[[col]] <- trim_taxa_names(df_taxa[[col]])
  }
  # Set new rownames
  if (use_taxa_names) {
    rownames(df_taxa) <- df_taxa[, ncol(df_taxa)]
  } else {
    otu_index <- paste0("Otu", 1:nrow(df))
    # df_taxa$Otu <- otu_index # in case you would like an Otu column
    # taxa_cols <- names(df_taxa)[!names(df_taxa) %in% "Otu"]
    # taxa_cols <- taxa_cols[!grepl('Otu', taxa_cols)]

    # for (col in taxa_cols) {
    #   # df_taxa[nrow(df_taxa), col] <- 'Other'
    #   df_taxa[df_taxa == 'Other', col] <- 'Other'
    # }
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
