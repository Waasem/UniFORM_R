Package: UniFORM
Title: Universal ImmunoFluorescence Normalization for Multiplex Tissue Imaging
Version: 1.0.0
Authors@R: 
    person(given = "Asim Waqas",
           family = "Waqas", 
           role = c("aut", "cre"),
           email = "wangmar@ohsu.edu; asim.waqas@moffitt.org")
Description: R implementation of UniFORM for normalizing multiplex immunofluorescence data.
License: MIT
Encoding: UTF-8
LazyData: true
RoxygenNote: 7.2.0
Imports: 
    R6,
    ggplot2,
    dplyr,
    tidyr,
    mixtools,
    signal
Suggests: 
    testthat (>= 3.0.0),
    knitr,
    rmarkdown
Config/testthat/edition: 3