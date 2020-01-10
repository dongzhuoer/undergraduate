-   [preparation](#preparation)
-   [Cladistics 2015](#cladistics-2015)
-   [Scientific Reports 2016](#scientific-reports-2016)
-   [Cladistics 2017](#cladistics-2017)

<h1 align="center">
Reproduce Chronostratigraphic Background in Three Papers by Prof. Xie
</h1>
<p align="center">
Zhuoer Dong
</p>
<p align="center">
2018-12-13 (rendered on `r Sys.Date()`)
</p>

------------------------------------------------------------------------

preparation
===========

install needed R packages

``` r
remotes::install_github(c("dongzhuoer/ggcsgb"), upgrade = TRUE)
```

``` r
library(magrittr)
```

Set global theme for ggplot

``` r
ggplot2::theme_set(
    ggplot2::theme(
        axis.line.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.line.x = ggplot2::element_line()
    ) + ggplot2::theme_grey()
)
```

Cladistics 2015
===============

![Original Figure 3 in
<a href="https://doi.org/10.1111/cla.12137" class="uri">https://doi.org/10.1111/cla.12137</a>](ggcsgb/2015-Fig3.jpg)

``` r
"ggcsgb/2015.tree" %>% treeio::read.beast() %>% ggtree::ggtree() +
    ggtree::geom_tiplab(size = 3) + ggplot2::geom_blank(ggplot2::aes(x = 370)) + 
    ggcsgb::chrono_strati_arg() + 
    ggcsgb::chrono_strati_bar() + 
    ggcsgb::chrono_strati_axis() + 
    ggcsgb::chrono_strati_label()
```

<img src="ggcsgb/2015rep-1.png" width="1152" />

Scientific Reports 2016
=======================

![Original Figure 3 in
<a href="https://doi.org/10.1038/srep38939" class="uri">https://doi.org/10.1038/srep38939</a>](ggcsgb/2016-Fig3.png)

``` r
"ggcsgb/2016.tree" %>% treeio::read.beast() %>% ggtree::ggtree() +
    ggtree::geom_tiplab(size = 3) + ggplot2::geom_blank(ggplot2::aes(x = 570)) + 
    ggcsgb::chrono_strati_arg() + 
    ggcsgb::chrono_strati_bar() + 
    ggcsgb::chrono_strati_axis() + 
    ggcsgb::chrono_strati_label()
```

<img src="ggcsgb/2016rep-1.png" width="1152" />

Cladistics 2017
===============

![Original Figure 6 in
<a href="https://doi.org/10.1111/cla.12232" class="uri">https://doi.org/10.1111/cla.12232</a>](ggcsgb/2017-Fig6.jpg)

``` r
"ggcsgb/2017.tree" %>% treeio::read.beast() %>% ggtree::ggtree() +
    ggtree::geom_tiplab(size = 3) + ggplot2::geom_blank(ggplot2::aes(x = 390)) + 
    ggcsgb::chrono_strati_arg(offset_first = 1) + 
    ggcsgb::chrono_strati_bar() + 
    ggcsgb::chrono_strati_axis() + 
    ggcsgb::chrono_strati_label()
```

<img src="ggcsgb/2017rep-1.png" width="1152" />
