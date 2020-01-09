-   [preparation](#preparation)
-   [input files](#input-files)

<h1 align="center">
Analysis correlation of coverage of duplex group’s two segments
</h1>
<p align="center">
Zhuoer Dong
</p>
<p align="center">
2018-12-19 (rendered on `r Sys.Date()`)
</p>

------------------------------------------------------------------------

``` r
# install needed R packages
remotes::update_packages(c("magrittr", "ggplot2", "dplyr", "stringr", "purrr", "tidyr", "testthat", "rlang", "cowplot"), upgrade = TRUE)
remotes::install_github(c("dongzhuoer/paristools"), upgrade = TRUE)
```

preparation
===========

``` r
library(magrittr)
library(dplyr)
library(ggplot2)
library(stringr)
library(purrr)
```

Set global theme for ggplot

``` r
theme_set(theme_classic(base_size = 24))
```

helper function

``` r
as_integer <- . %>% as.integer() %T>% {.[is.na(.)] = 0L}
testthat::expect_identical(as_integer(c("12", "")), c(12L, 0L))
```

input files
===========
