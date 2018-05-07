### Generate the citation info for the test data, to be put into the new spreadsheet format for the compadrinos

library(ProjectTemplate)
load.project()

# The subsetting is needed to convert the data frame back to a character vector
doi_list <- read.csv("reports/doi_list.csv", stringsAsFactors = FALSE, 
                     colClasses = "character", header = FALSE)[[1]]

focal_entries <- match(doi_list, comadre$metadata$DOI.ISBN)
focal_db <- comadre$metadata[focal_entries, ]

output_fields <- c("Authors", "Journal", "YearPublication", "DOI.ISBN", "AdditionalSource")
write.csv(focal_db[output_fields], "reports/test_src.csv", quote = FALSE, 
          fileEncoding = "UTF-8", row.names = FALSE)
