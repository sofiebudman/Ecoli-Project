load_seq_from_gen_bank <- function(IDs, seq.names = IDs, species.names = TRUE, as.character = TRUE) {
  base_url <- "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&id="
  sequences <- list()  # Empty list to store results

  for (i in seq_along(IDs)) {
    # URL to get info
    url <- paste0(base_url, IDs[i], "&report=fasta&retmode=text")

    # Read data from NCBI
    raw_data <- readLines(url, warn = FALSE)

    # Extract and clean
    sequence <- paste(raw_data[-1], collapse = "")  # Remove > header

    # if statement for species name
    species <- NA
    if (species.names) {
      species_url <- paste0("https://www.ncbi.nlm.nih.gov/nuccore/", IDs[i])
      species_page <- readLines(species_url, warn = FALSE)


      match <- grep("title=\".*?\\(", species_page, value = TRUE)
      if (length(match) > 0) {
        species <- sub(".*title=\"(.*?) \\(.*", "\\1", match[1])
      }
    }

    # Assign sequence name
    name <- if (!is.na(species) && species.names) species else seq.names[i]

    # Store sequence
    sequences[[name]] <- if (as.character) sequence else strsplit(sequence, "")[[1]]
  }

  return(sequences)
}
