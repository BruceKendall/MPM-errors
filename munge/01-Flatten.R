# Flatten the databases
flat_compadre <- convert2flat(compadre, Aonly = FALSE)
flat_comadre <- convert2flat(comadre, Aonly = FALSE)

cache("flat_comadre")
cache("flat_compadre")
