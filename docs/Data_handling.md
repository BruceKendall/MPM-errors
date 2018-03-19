# Data handling
## Source data
The current versions are:

* Comadre: 2.0.1
* Compadre: 4.0.1

The data will need to be downladed from http://www.compadre-db.org. Please ensure that the raw data files (e.g., `data/COMADRE_v.2.0.1.Rdata`), the cached data (e.g., `cache/comadre.Rdata`) and the flat data tables (e.g., `cache/flat_comadre.Rdata`) are excluded from the git repository.

## Usage tips
Be sure to set `munging: TRUE` in `config/global.dcf` the first time you run the project with new com[p]adre data. Then set it to `FALSE` to speed future loading.

## General strategy
The general strategy is to use the `convert2flat()` function from the **Mage** package (to be installed with `devtools::install_github("jonesor/compadreDB/mage")`) to link the three components of the database together, creating `flat_comadre` and `flat_compadre`. This makes it easier to draw upon information in the `matrixClass` array, but requires extra effort to reference elements of the vectors and matrices. The **Mage** package provides a way to recreate the matrices; the `get_element()` function in `lib/helpers.R` allows access to elements of the character vectors.

Working with the flat datafiles, one constructs new columns containing the new information about the populations. If these can't be calculated directly from existing information in the database, but require qualitative assessment of variables (e.g. classnames) or reference to the original paper, please bring the new information in via a lookup table that is added to the `data` directory.

Eventually, a copy of the metadata table will be created that has all of the original metadata as well as extra columns with the additional classification.
