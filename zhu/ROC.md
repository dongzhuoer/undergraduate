```r
library(tidyverse)
```



# 05/09 

> 做一个用LEM4 表达数据预测内分泌治疗预测的ROC分析给我吧，GEO GSE2034,GSE4922, GSE22219
> 
> 最多两天时间
> 
> --- Pro. Zhu, 18:20


- [ANKLE2](https://www.ncbi.nlm.nih.gov/gene/23141)
- [AURKA](https://www.ncbi.nlm.nih.gov/gene/6790)



## GSE22219

There is neither ANKLE2 nor AURKA in GSE22219, I even searched all NM_*（ GI_* is not provided in NCBI).



## GSE2034

```r
matrix_raw <- readr::read_lines('data-raw/GSE/GSE2034_series_matrix.txt.gz')
info <- matrix_raw %>% {
    tibble(
	    sample = str_subset(., '^!Sample_geo_accession') %>% 
            str_extract_all('GSM\\d+') %>% {.[[1]]},
	    relapse = str_subset(., 'bone relapses') %>% 
            str_extract_all('(?<=: )[01]') %>% {.[[1]] == '1'}
    )
} %T>% print

matrix <- qGSEA::read_gse_matrix('data-raw/GSE/GSE2034_series_matrix.txt.gz')
chip <- qGSEA::read_gse_soft('data-raw/GSE/GSE2034_family.soft.gz')
expression_df <- chip %>% filter(symbol %in% c('ANKLE2', 'AURKA')) %>%
	inner_join(matrix) %T>% print
```

```r
expression_df_collapsed <- expression_df %>% group_by(symbol) %>% 
	summarise_at(-1:-2, . %>% max(na.rm = T)) %T>% print

value_df <- expression_df_collapsed %>% {
    tibble(
        sample = colnames(.)[-1], 
        ANKLE2 = filter(., symbol == 'ANKLE2') %>% select(-1) %>% {as.matrix(.)[1, ]},
        AURKA  = filter(., symbol == 'AURKA')  %>% select(-1) %>% {as.matrix(.)[1, ]}
    )
} %T>% print

gene2relapse <- inner_join(value_df, info) %T>% print
```

```r
genes <- c('AURKA', 'ANKLE2')
rocs <- lapply(
    genes, 
    . %>% {pROC::roc(gene2relapse[['relapse']], gene2relapse[[.]])}
)

cowplot::plot_grid(
	plotlist = lapply(
        rocs, 
        . %>% {
            pROC::ggroc(., color = 'blue') + 
                geom_abline(intercept = 1, color = 'green')
        }
    ),
	labels = sapply(rocs, . %>% pROC::auc() %>% formatC(3)) %>% 
        paste0(genes, '    AUC=', .)
)
```



## GSE4922

```
ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE4nnn/GSE4922/soft/GSE4922_family.soft.gz
ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE4nnn/GSE4922/matrix/GSE4922-GPL96_series_matrix.txt.gz
ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE4nnn/GSE4922/matrix/GSE4922-GPL97_series_matrix.txt.gz
```

from Ivshina2006

- Singapore cohort, no information on recurrence or cause of death was available  

- Uppsala cohort, Information pertaining to breast cancer therapy, clinical follow-up, and sample processing are described elsewhere 

```bash
cat data-raw/GSE/GSE4922_family.soft | grep -P '(?:ANKLE2|AURKA)\s' > data-raw/GSE/GSE4922-ANKLE2_AURKA.soft
cat data-raw/GSE/GSE4922-GPL96_series_matrix.txt | grep -P '^!Sample_' > data-raw/GSE/GSE4922-GPL96-sample.info
```

```r
soft2 <- read_lines('data-raw/GSE/GSE4922-ANKLE2_AURKA.soft') %>% 
	{tibble(probe = str_extract(., '\\w+'), symbol = str_extract(., 'ANKLE2|AURKA'))}
matrix_raw2 <- read_lines('data-raw/GSE/GSE4922-GPL96_series_matrix.txt.gz')

sample_df2 <- matrix_raw2 %>% str_subset('^!Sample_') %>% str_split('\t') %>% 
    do.call(rbind, .)%>% plyr::aaply(2, . %>% paste0(collapse = '\t')) %>% 
    paste0(collapse = '\n')  %>% read_tsv(F)

# "endocrine therapy only", "No systemic therapy"
tmp <- sample_df2 %>% select(sample = 2, 18:19) %>% slice(-1) %>% 
    filter_at(-1, all_vars(!is.na(.))) %>% 
    mutate_at(-1, . %>% str_sub(-1, -1) %>% {. == '1'}) %>% filter(X18)

info2 <- bind_cols(
	sample_df2[1,] %>% str_which('^!Sample_geo_accession') %>% 
        select(sample_df2, sample = .) %>% slice(-1),
	sample_df2[2,] %>% str_which('ER status:') %>% 
        select(sample_df2, ER = .) %>% slice(-1) %>% 
        mutate_at(1, . %>% str_sub(-1, -1) %>% {. == '+'}),
	sample_df2[2,] %>% str_which('DFS EVENT') %>% 
        select(sample_df2, relapse = .) %>% slice(-1) %>% 
        mutate_at(1, . %>% str_extract('(?<=: )\\d$') %>% {. == '1'})
) %T>% print

matrix <- qGSEA::read_gse_matrix('data-raw/GSE/GSE4922-GPL96_series_matrix.txt.gz')
expression_df2 <- matrix %>% rename(probe = ID_REF) %>% inner_join(soft2, .)
```

```r
expression_df_collapsed2 <- expression_df2 %>% group_by(symbol) %>% 
	summarise_at(-1:-2, . %>% max(na.rm = T)) %T>% print

value_df2   <- expression_df_collapsed2 %>% {
    tibble(
        sample = colnames(.)[-1], 
        ANKLE2 = filter(., symbol == 'ANKLE2') %>% select(-1) %>% {as.matrix(.)[1, ]},
        AURKA =  filter(., symbol == 'AURKA')  %>% select(-1) %>% {as.matrix(.)[1, ]}
    )
} %T>% print

gene2relapse2 <- inner_join(value_df2, info2) %>% inner_join(tmp) %T>% print %>% filter(ER)

genes <- c('AURKA', 'ANKLE2')
rocs <- lapply(
    genes, 
    . %>% {pROC::roc(gene2relapse2[['relapse']], gene2relapse2[[.]])}
)

cowplot::plot_grid(
	plotlist = lapply(
        rocs, 
        . %>% {
            pROC::ggroc(., color = 'blue') + 
                geom_abline(intercept = 1, color = 'green')
        }
    ),
	labels = sapply(rocs, . %>% pROC::auc() %>% formatC(3)) %>% 
        paste0(genes, '    AUC=', .)
)
```


