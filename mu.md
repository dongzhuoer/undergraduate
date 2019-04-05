# KEGG enrichment of DE genes

朱老师介绍的一个学生，穆成龙，就简单的 enrich 和 dotplot。不过那是我第一次遇到选取合适的 cutoff 和要分析的 factor 的难题，最后取了一个生物学上还能说过去的结果。

这个项目算是可重复科研的一个早期尝试吧，当时还想着让 Travis 自动运行代码，并把结果更新到 GitHub，而且为了显示所有的因素组合，都开始 programmatically generate knitr child document。



## read data

First, let's look at our data. The first column Entrez gene id, others are normalized expression value (measured by FPKM) of 20 samples. In sample name, after `BF`, `H` means high-fat feed, `K` (vs `WC`) means SDHB knock-out.

```r
# I remove some columns from original `.xlsx` to save space
fpkm <- readxl::read_excel('data-raw/穆成龙+SDHB+KO+RNA+seq+result.FPKM.xlsx') %>% 
    dplyr::select(1, dplyr::starts_with('BF')) %T>% print;
#> # A tibble: 18,458 x 21                                                       
#>    gene_id BFKC106_FPKM BFKC128_FPKM BFKC144_FPKM BFKC61_FPKM BFKC63_FPKM
#>      <dbl>        <dbl>        <dbl>        <dbl>       <dbl>       <dbl>
#>  1   11657       47449.       30442.       34459.      39868.      23542.
#>  2   11816       17594.       16645.       16740.      20437.      15255.
#> ...
```



## differential expression analysis

Thanks for the first answer of this [quesiton](https://support.bioconductor.org/p/56275/):

> If FPKM is really all you have, then convert the values to a log2 scale and do an ordinary limma analysis as you would for microarray data, using eBayes() with trend=TRUE.

First, let us calculate differentially expressed gene with **limma**. We will examine the difference caused by feed, i.e, high-fat feed or normal feed

```r
# the first column is entrez id
mat <- fpkm %>% dplyr::select(-1) %T>% print;
#> # A tibble: 18,458 x 20
#>    BFKC106_FPKM BFKC128_FPKM BFKC144_FPKM BFKC61_FPKM BFKC63_FPKM BFKH24_FPKM
#>           <dbl>        <dbl>        <dbl>       <dbl>       <dbl>       <dbl>
#>  1       47449.       30442.       34459.      39868.      23542.      43868.
#>  2       17594.       16645.       16740.      20437.      15255.      22071.
#> ...

# high-fat feed or normal feed
level <- stringr::str_detect(colnames(mat), 'H') %T>% print;
#>  [1] FALSE FALSE FALSE FALSE FALSE  TRUE  TRUE  TRUE  TRUE  TRUE FALSE FALSE
#> [13] FALSE FALSE FALSE  TRUE  TRUE  TRUE  TRUE  TRUE

# limma analysis with FPKM
DE <- mat %>% log2 %>% 
    limma::lmFit(design = cbind(Intercept = 1, level)) %>% 
    limma::eBayes(trend = T)

# subset differentially expressed genes	using a p value cutoff of 0.05
DE_gene <- dplyr::filter(fpkm, !!DE$p.value[ , 2] < 0.05)[[1]] %T>% {print(head(.))}
#> [1] 11807 11806 13627 14862 16644 12266
```



## KEGG enrichment

Now we can perform KEGG enrichment and plot the result

```r
ekegg <- clusterProfiler::enrichKEGG(DE_gene, 'mmu') %T>% 
    readr::write_rds('data/kegg-enrich.rds')

clusterProfiler::dotplot(ekegg, showCategory = Inf) + 
    ggplot2::labs(title = 'KEGG enrichment of DE (p < 0.05) genes casued by food') +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 12))
```
![](image/mu/kegg-enrich.png)



## browse kegg pathway 

```r
ekegg <- readr::read_rds('data/kegg-enrich.rds')
clusterProfiler::browseKEGG(ekegg, ekegg@result$ID[5])
```

![](image/mu/mmu01200.png)

I saved the result on 2018-05-31, when I rerun the code on 2019-04-05, it seems quite different.



## epilogue

I tried a variety of combinations (the result is not shown), i.e, using 0.01 vs 0.01 for DE cutoff +

- the factor of food vs genotype
- the factor of food in knocked out vs wild-type mouse
- the factor of genotype given high-fat cs normal feed 

When I look back on 2019-04-05, I find a interesting point. If you looks at the first row and then the second row (I use p<0.01 for DE), it seems like knock-out of SDHB has similar effect as high hat feed.

![](image/mu/various-combinations.png)
