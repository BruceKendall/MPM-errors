# Get the i'th element of a string-form stage information vector produced by
# convert2flat(). If str_vec is a vector of such objects then a vector giving the i'th
# element of each is returned. By default, i = 1.
get_element <- function(str_vec, i = 1) {
  plyr::laply(strsplit(str_vec, split=" [|] "), function(x) x[i])
}
