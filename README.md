# MPM-errors

This is the repository for analysis and writing associated with the paper "Persistent problems in the construction of matrix population models," to be submitted to the special issue of *Ecological Modelling* on matrix population models.

To reproduce the analyses, you will need to:
* Clone or download the repository
* Download the COMPADRE and COMADRE databases from http://www.compadre-db.org. Put the files in the `data` subdirectory. Note that reproducibility is only guaranteed if you use v. 2.0.1 of COMADRE and v. 4.0.1 of COMPADRE
* Install the **ProjectTemplate** library, and run `ProjectTemplate::load.project()` at the beginning of your session. More details are below.

**Collaborators:** Please use the issue tracker on github (https://github.com/BruceKendall/MPM-errors/issues) to make comments and to-do items.

## Required libraries
- **Mage** (`devtools::install_github("jonesor/compadreDB/Mage")`)



## ProjectTemplate

The compadre-plus project was set up with the **ProjectTemplate** package, as described below:

ProjectTemplate is an R package that helps you organize your statistical
analysis projects. Since you're reading this file, we'll assume that you've
already called `create.project()` to set up this project and all of its
contents.

To load your new project, you'll first need to `setwd()` into the directory
where this README file is located. Then you need to run the following two
lines of R code:

	library('ProjectTemplate')
	load.project()

After you enter the second line of code, you'll see a series of automated
messages as ProjectTemplate goes about doing its work. This work involves:
* Reading in the global configuration file contained in `config`.
* Loading any R packages you listed in he configuration file.
* Reading in any datasets stored in `data` or `cache`.
* Preprocessing your data using the files in the `munge` directory.

Once that's done, you can execute any code you'd like. For every analysis
you create, we'd recommend putting a separate file in the `src` directory.
If the files start with the two lines mentioned above:

	library('ProjectTemplate')
	load.project()

You'll have access to all of your data, already fully preprocessed, and
all of the libraries you want to use.

For more details about ProjectTemplate, see http://projecttemplate.net
