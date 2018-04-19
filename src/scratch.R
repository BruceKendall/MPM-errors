library(ProjectTemplate)
load.project()

## experiment with downloading pdfs of papers
library(fulltext)

# Pick out recent papers with ontognetic stages
ssData <- subsetDB(comadre, MatrixCriteriaOntogeny %in% c("Based on age", "Yes") & as.numeric(YearPublication) > 2010)
dois <- unique(ssData$metadata$DOI.ISBN)

ft_search(dois[1])
test_pdfs <- ft_get(dois, type = "pdf")

doitest <- "10.1093/jme/tjv112"
tempdata <- subsetDB(comadre, DOI.ISBN == doitest)
tempdata
convert2flat(tempdata)
