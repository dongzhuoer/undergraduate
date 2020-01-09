# Overview

孙小雅师姐人挺好的，可惜我那时能力不足，尤其是不会教别人。最后混了个三作也是不错的。

```r
library(tidyverse)
```






# convert annotation, 每行一基因 to 每行一GO

输入文件是每个基因有多个 GO，输入要求每个 GO 一行，对应的基因用 `;` 串联起来。

最开始小雅师姐要求用 ontoloty 和 term 作为 key，对应的多个 GO ID 也用 `;` 串联起来。但我后来修改代码后，和原来的结果比对，发现这些 GO ID 好像都是一样的，然后我突然意识到  ontoloty 和 term 其实也就确定了 GO ID，验证之后我就把三者一起作为 key 了。

```r
# I remove extra columns from original `.xlsx` to save space
annotation_gene <- readxl::read_excel('xiaoya/annotation-gene.xlsx')

annotation_go <- annotation_gene %>% 
    filter(!is.na(go)) %>% select(gene = gene_id, go) %>% 
    mutate(go = str_split(go, '`')) %>% unnest() %>%
    mutate(
        ontoloty = str_extract(go, '(?<=\\^)\\w+'), 
        term = str_replace(go, '^GO:\\d+\\^\\w+\\^', ''),
        id = str_extract(go, '^GO:\\d+')
    ) %>% select(-go) %>% 
    group_by(ontoloty, term, id) %>%
    summarize(n = n(), genes = paste0(gene, collapse = ';')) %T>% print
#> # A tibble: 10,478 x 5
#> # Groups:   ontoloty, term [10,478]
#>    ontoloty     term                  id         n genes                        
#>    <chr>        <chr>                 <chr>  <int> <chr>                        
#>  1 biological_… 'de novo' AMP biosyn… GO:00…     3 CL1815.Contig1_All;CL1815.Co…
#>  2 biological_… 'de novo' CTP biosyn… GO:00…     1 Unigene15161_All 

write_tsv(annotation_go, 'annotation_go.tsv')
```






# merge reads count files

就是把几个 `.csv` 文件合并一下，代码亲测能用，就不优化了。

```r
genes <- readr::read_tsv('xiaoya/Gene.list', 'gene', 'c')

datas <- plyr::llply(
    dir('xiaoya', 'AlignReadsperGene', full.names = T),
    function(file) {
        data <- readr::read_tsv(file, c('gene', 'reads'), 'ci') %>% 
            merge(genes, ., by = 'gene', all.x = T, sort = T) %>% 
            dplyr::mutate(reads = ifelse(is.na(reads), 0, reads))
        
        sample <- stringr::str_extract(file, '\\w+(?=\\.AlignReadsperGene)')
        colnames(data)[2] = sample

        dplyr::as_tibble(data)
    }
)

# datas are already sorted when merge, so I can bind them directly
result <- plyr::llply(datas, . %>% dplyr::select(2)) %>% 
    dplyr::bind_cols() %>% dplyr::bind_cols(dplyr::select(datas[[1]], 1), .)

readr::write_csv(result, 'reads-per-gene.csv')
```






# differential expression enrichment

`trinotate_annotation_report.tsv` : `Trinotate Trinotate.sqlite report > ...`



## prepare data

```r
higher <- readr::read_csv('xiaoya/higher.csv')
lower <- readr::read_csv('xiaoya/lower.csv')
annotation <- readr::read_tsv('xiaoya/trinotate.tsv.gz')

df <- annotation %>% dplyr::rename(gene_id = '#gene_id') %>% 
    dplyr::select(gene_id, eggnog: gene_ontology_pfam) %>% 
    dplyr::mutate(
        eggnog = stringr::str_extract(eggnog, '\\w+(?=\\^)'), 
        GO = paste(gene_ontology_blast, gene_ontology_pfam)
    ) %>% dplyr::mutate(
        GO = stringr::str_extract_all(GO, 'GO:\\d+') %>% 
            plyr::laply(. %>% unique %>% paste0(collapse = ','))
    )
```



## GO

```r
split_GO <- . %>% {stringr::str_split(., ',')[[1]]}
GO2gene <- df %>% dplyr::select(1, GO) %>% dplyr::filter(nchar(GO) > 0) %>% 
    plyr::ddply(
        'gene_id', 
        . %>% {data.frame(GO = split_GO(.$GO), gene_id = .$gene_id)}
    ) %>% dplyr::as_tibble()
enricher_go <- function(change) {
    gene <- merge(change, df, by = 1, all.x = T) %>% 
        dplyr::filter(nchar(GO) > 0) %>% {.$ID}
    clusterProfiler::enricher(gene, TERM2GENE = GO2gene)
}
```

```r
egoh <- enricher_go(higher) %T>% write_rds('xiaoya/enrich-go-higher.rds')
enrichplot::dotplot(egoh)
```

![](xiaoya/enrich-go-higher.png)

```r
egol <- enricher_go(lower) %T>% write_rds('xiaoya/enrich-go-lower.rds')
enrichplot::dotplot(egol)
```

![](xiaoya/enrich-go-lower.png)

Aboving is [Fig 4](https://doi.org/10.1093/jisesa/iey114#F4) of the paper, though I forget how to get GO term (title).



## KEGG

```r
split_KEGG <-  . %>% {stringr::str_split(., '`')[[1]]} %>% 
    stringr::str_subset('KEGG')
KEGG2gene <- df %>% dplyr::select(1, Kegg) %>% dplyr::filter(Kegg != '.') %>% 
    plyr::ddply(
        'Kegg', 
        . %>% {data.frame(Kegg = split_KEGG(.$Kegg), gene_id = .$gene_id)}
    ) %>% dplyr::as_tibble()
enricher_kegg <- function(change) {
    gene <- merge(change, df, by = 1, all.x = T) %>% 
        dplyr::filter(Kegg != '.') %>% {.$ID}
    clusterProfiler::enricher(gene, TERM2GENE = KEGG2gene)
}
```

```r
ekeggh <- enricher_kegg(higher)
enrichplot::dotplot(ekeggh)
```

![](xiaoya/enrich-kegg-higher.png)

```r
ekeggl <- enricher_kegg(lower)
enrichplot::dotplot(ekeggl)
```

![](xiaoya/enrich-kegg-lower.png)



## KO

```r
KO2gene <- df %>% dplyr::filter(stringr::str_detect(Kegg, 'KO'))
KO2gene %<>% dplyr::mutate(KO = stringr::str_extract(Kegg, 'KO:K\\d+'))
KO2gene %<>% dplyr::select(KO, gene_id)
enricher_ko <- function(change) {
    gene <- merge(change, df, by = 1, all.x = T) %>% 
        dplyr::filter(stringr::str_detect(Kegg, 'KO')) %>% {.$ID}
    clusterProfiler::enricher(gene, TERM2GENE = KO2gene)
}
```

```r
ekoh <- enricher_ko(higher)
enrichplot::dotplot(ekoh)
```

![](xiaoya/enrich-ko-higher.png)

```r
ekol <- enricher_ko(lower)
enrichplot::dotplot(ekol)
```

![](xiaoya/enrich-ko-lower.png)
