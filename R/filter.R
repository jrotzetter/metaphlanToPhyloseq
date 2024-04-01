#' Filter MetaPhlAn profile by taxonomic rank
#'
#' @description Function to filter rows of a MetaPhlAn profile based on the
#' selected taxonomic level. Used if profile not yet separated into taxonomic
#' ranks.
#'
#' @param df A data frame. Assumes clade names are of the structure
#'  ***k__Bacteria|p__Firmicutes|c__Bacilli|o__Lactobacillales|f__Lactobacillaceae***
#'  and found in the 'clade_name' column.
#' @param taxa_lvl A `character` string. The taxonomic level to extract
#'  ('kingdom', 'phylum', 'class', 'order', 'family', 'genus', or 'species').
#'   First letter abbreviations (e.g., 's') are also accepted.
#'
#' @returns A filtered data frame.
#'
#' @export
#'
#' @examples
#' filtered_profile <- filter_taxa_lvl(merged_abundance_profiles, "Class")
#' head(filtered_profile)
#'
#' @author Jérémy Rotzetter
filter_taxa_lvl <- function(df, taxa_lvl) {
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
  # Extract the first letter of taxa_lvl
  first_letter <- substr(taxa_lvl, 1, 1)

  # Filter the dataframe based on taxa_lvl
  # df_taxa_lvl <- df[grepl(paste0("\\|", first_letter,"__[^|]*$"), df$clade_name), ]
  df_taxa_lvl <- df[grepl(paste0(first_letter, "__[^|]*$"), df$clade_name), ]
  rownames(df_taxa_lvl) <- NULL

  return(df_taxa_lvl)
}


#' Apply threshold filtering to a data frame
#'
#' @description This function applies threshold filtering to a data frame,
#'  summing values below the threshold for each column, storing the sum in the
#'  new 'Other' row, but only if all values within a row (taxon) across all
#'  columns are below the threshold.
#'
#' @param df The input data frame of merged MetaPhlAn profiles.
#' @param threshold The threshold value for filtering.
#' @param merged_profiles `Logical`; if `TRUE` (the default) the file to be
#' loaded is assumed to be multiple merged MetaPhlAn profiles.
#'
#' @details
#' Function can be applied to a single MetaPhlAn profile by setting
#' `merged_profiles = FALSE`. Likewise also works with merged profiles
#' by slicing the data frame to only the clade_name and the numeric sample
#' column of interest. In both cases a `clade_name` column is required.
#'
#' @returns The filtered data frame.
#'
#' @examples
#' df <- data.frame(
#'   clade_name = c("A", "B", "C", "D"),
#'   col1 = c(10, 20, 5, 65),
#'   col2 = c(8, 15, 3, 38),
#'   col3 = c(5, 35, 4, 6)
#' )
#' print(df)
#' df_thresh <- filter_threshold(df, 11)
#' print(df_thresh)
#'
#' @export
#'
#' @author Jérémy Rotzetter
filter_threshold <- function(df, threshold, merged_profiles = TRUE) {
  # Add a new row to hold aggregated counts below threshold
  index_to_append <- nrow(df) + 1
  df[index_to_append, "clade_name"] <- "Other"
  # Get column names of only numeric columns
  sample_cols <- colnames(df[sapply(df, is.numeric)])
  # For each column except clade_name
  for (col in sample_cols) {
    if (col != "clade_name") {
      # Sum values below threshold for current column, excluding Other clade
      # df[df$clade_name == 'Other', col] <- sum(df[df[col] < threshold & df$clade_name != "Other", col])
      all_below_thresh <- rowSums(df[sample_cols] < threshold) == length(sample_cols)
      df[df$clade_name == "Other", col] <- sum(
        df[all_below_thresh & df$clade_name != "Other", col]
      )
    }
  }
  if (merged_profiles) {
    # sample_cols <- get_sample_cols(df)
    sample_cols <- colnames(df[sapply(df, is.numeric)])
    ids <- df$clade_name
    # Create thresholded version of data
    df_thresh_removed <- df[sample_cols]
    df_thresh_removed <- apply(
      df_thresh_removed, 2, function(x) ifelse(x < threshold, 0, x)
    )
    df_thresh_removed <- as.data.frame(df_thresh_removed)
    df_thresh_removed$clade_name <- ids

    # Drop rows with rowsum of 0, except if clade_name == 'Other'
    df <- df[rowSums(
      df_thresh_removed[, sample_cols]
    ) != 0 | df$clade_name == "Other", ]
    # Reset row names
    rownames(df) <- NULL
  } else {
    df <- subset(df, df[[sample_cols]] >= threshold | df$clade_name == "Other")
  }

  return(df)
}
