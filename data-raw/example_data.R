## code to prepare the example datasets goes here

path_single_profile <- list_example("SRS014470-Tongue_dorsum_profile.txt")
single_abundance_profile <- load_metaphlan_profile(
  path_single_profile,
  merged_profiles = FALSE)

usethis::use_data(single_abundance_profile, overwrite = TRUE)

path_merged_profiles <- list_example("merged_abundance_table.txt")
merged_abundance_profiles <- load_metaphlan_profile(path_merged_profiles)

usethis::use_data(merged_abundance_profiles, overwrite = TRUE)

species_only <- filter_taxa_lvl(merged_abundance_profiles, "s")
usethis::use_data(species_only, overwrite = TRUE)
