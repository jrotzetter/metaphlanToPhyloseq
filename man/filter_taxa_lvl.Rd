% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filter.R
\name{filter_taxa_lvl}
\alias{filter_taxa_lvl}
\title{Filter MetaPhlAn profile by taxonomic rank}
\usage{
filter_taxa_lvl(df, taxa_lvl)
}
\arguments{
\item{df}{A data frame. Assumes clade names are of the structure
\emph{\strong{k__Bacteria|p__Firmicutes|c__Bacilli|o__Lactobacillales|f__Lactobacillaceae}}
and found in the 'clade_name' column.}

\item{taxa_lvl}{A \code{character} string. The taxonomic level to extract
('kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species' or 't'
(SGB)). First letter abbreviations (e.g., 's') are also accepted.}
}
\value{
A filtered data frame.
}
\description{
Function to filter rows of a MetaPhlAn profile based on the
selected taxonomic level. Used if profile not yet separated into taxonomic
ranks.
}
\examples{
filtered_profile <- filter_taxa_lvl(merged_abundance_profiles, "Class")
head(filtered_profile)

}
\author{
Jérémy Rotzetter
}
