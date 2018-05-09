### Generate the sample of studies to examine

library(ProjectTemplate)
load.project()

# Pull out all the non-human studies in Comadre
non_human <- subsetDB(comadre, GenusAccepted != "Homo")
study_list <- unique(non_human$metadata[c("Authors", "Journal", "YearPublication", 
                                        "DOI.ISBN")])

# Restrict to peer-reviewed (approximately; there are a few journal articles w/out 
#  DOIs that will be dropped)
study_list <- subset(study_list, DOI.ISBN != "None")

# Augment to full metadata
study_list <- comadre$metadata[row.names(study_list), ]

# Make numeric years (have to drop trailing letters)
study_list$YearPub <- as.numeric(sub("[a-z]", "", study_list$YearPublication))

# Look at frequencies through time
hist(study_list$YearPub)

# Look at studies with developmental stages
study_list_dev <- subset(study_list, MatrixCriteriaOntogeny %in% c("Based on age", "Yes"))
hist(study_list_dev$YearPub)

# Strategy: take all the 20th c. studies (65) and add 60 from 21st c.
study_list_20 <- subset(study_list, YearPub < 2000)
study_list_21 <- subset(study_list, YearPub >= 2000)

sample_21 <- sample(row.names(study_list_21), size = 60, replace = FALSE)

full_sample <- c(row.names(study_list_20), sample_21)

# Now augment to ensure good dev coverage in developmental stages
