# Overview

extra work during thesis under Prof. Zhu's supervision

```r
library(tidyverse)
```






# 04/19 Lem4 and mki67 co-expression

Prof. Zhu want to know the relationship between lem4 and mki67, he appoints several dataset, and finally choose the best result, GSE2990.

那时发现，虽然 RNA-seq 的 BioProject 也会在 GEO 注册一个数据集，但是下载下来的文件都是空的，没有任何表达数据，还得靠自己 assembly。

```r 
lem4_mki67s <- mclapply(
    c('GSE2990', 'GSE2034', 'GSE16446', 'GSE20685', 'GSE17705', 'GSE11121'),
    . %>% {
        chip <- qGSEA::read_chip(paste0('MKI67/', ., '.chip'))
        matrix <- qGSEA::read_txt(paste0('MKI67/', ., '.txt'))
        lem4_mki67 <- qGSEA::make_phenotype(matrix, chip, c('ANKLE2', 'MKI67')) 
        add_column(lem4_mki67, accession = .)
    }
) %>% bind_rows()
```

```r
ggplot(lem4_mki67s, aes(log10(ANKLE2), log10(MKI67))) + 
    geom_point() + geom_smooth(method = 'lm') + 
    facet_wrap(~accession, scales = "free") +
    ggpmisc::stat_poly_eq(
        formula = y ~ x, 
        aes(label = paste(..rr.label.., sep = "~~~")), 
        parse = TRUE
    ) 
```

```r
ggplot(filter(lem4_mki67s, accession == 'GSE2990'), aes(log10(ANKLE2), log10(MKI67))) + 
    geom_point() + geom_smooth(method = 'lm', se = F, color = 'black')  +
    ggpmisc::stat_poly_eq(
        formula = y ~ x, 
        aes(label = paste0('atop(', ..eq.label.., ' ,', ..rr.label.., ')'), size = 20), 
        parse = TRUE
    ) + cowplot::theme_cowplot(12) + 
    labs(x = 'ANKLE2 (Log10)', y = 'MKI67 (Log10)')
```
![](zhu/lem4-mki67-GSE2990.png)






# Lem4 expression in PRJNA390636 (GSE100075)



## 1. download raw data

```r
read_tsv('zhu/zhu/PRJNA390636.txt') %>% 
    .$fastq_ftp %>% str_split(';') %>% unlist %>% 
    paste0('ftp://', .) %T>% print %>% write_lines('sra.md')
#>  [1] "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR568/000/SRR5687500/SRR5687500_1.fastq.gz"
#>  [2] "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR568/000/SRR5687500/SRR5687500_2.fastq.gz"
#> ...
#> [30] "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR568/004/SRR5687514/SRR5687514_2.fastq.gz"
```

then download `.fastq.gz` in `sra.md` into `GSE100075/fastq/`.

```r
# check that the following R code includes all input files 
identical(
    dir('GSE100075/fastq/') %>% stringr::str_remove('_\\d.fastq.gz') %>% unique,
    formatC(0:14, width = 2, flag = '0') %>% paste0('SRR56875', .)
)
```

- download ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/grch38.tar.gz
- ftp://ftp.ensembl.org/pub/release-92/gtf/homo_sapiens/Homo_sapiens.GRCh38.84.gtf



## 2. assembly RNA-seq

although all support specifying threads number, only hisat can _really_ run in parallel

```r
# intermediate .sam file is quite large, hisat is IO intensive while samtools runs quickly
hisat_command <- 'hisat2 -p 24 --dta -x GSE100075/grch38/genome -1 GSE100075/fastq/SRR5687500_1.fastq.gz -2 GSE100075/fastq/SRR5687500_2.fastq.gz -S GSE100075/assembly/SRR5687500.sam; samtools sort -@ 24 -o GSE100075/assembly/SRR5687500.bam GSE100075/assembly/SRR5687500.sam; rm GSE100075/assembly/SRR5687500.sam'
formatC(0:14, width = 2, flag = '0') %>% paste0('SRR56875', .) %>% sapply(. %>% stringr::str_replace_all(hisat_command, 'SRR5687500', .))      %>% paste0(collapse = '; ')   %>% paste0('nohup bash -c "', ., '" &> /dev/null &') %>% write_lines('temp.sh')

# we use mclapply for paralleling, but IO limits CPU use efficientcy
stringtie_command <- 'stringtie -p 2 -G GSE100075/Homo_sapiens.GRCh38.84.gtf -o GSE100075/assembly/SRR5687500.gtf -l SRR5687500 GSE100075/assembly/SRR5687500.bam'
formatC(0:14, width = 2, flag = '0') %>% paste0('SRR56875', .) %>% sapply(. %>% stringr::str_replace_all(stringtie_command, 'SRR5687500', .))  %>% paste0('\'', ., '\'', collapse = ', ') %>% paste0('nohup R -e "parallel::mclapply(c(', ., '), system, mc.preschedule = F)" &> /dev/null &') %>% write_lines('temp.sh')
```

```bash
ls GSE100075/assembly/SRR*.gtf > GSE100075/assembly/mergelist.txt

stringtie --merge -p 24 -G GSE100075/Homo_sapiens.GRCh38.84.gtf -o GSE100075/assembly/merged.gtf GSE100075/assembly/mergelist.txt
```



## 3. calculate expression

```r
# we don't use parallel duo to IO bottleneck
stringtie_command2 <- 'stringtie -e -B -p 24 -G GSE100075/assembly/merged.gtf -o GSE100075/assembly/SRR5687500-ballgown.gtf -l SRR5687500 GSE100075/assembly/SRR5687500.bam'
formatC(0:14, width = 2, flag = '0') %>% paste0('SRR56875', .) %>% sapply(. %>% stringr::str_replace_all(stringtie_command2, 'SRR5687500', .)) %>% paste0(collapse = '; ')   %>% paste0('nohup bash -c "', ., '" &> /dev/null &') %>% write_lines('temp.sh')
```

```bash
for id in `ls *ballgown* | sed 's/-ballgown.gtf//'`; do cat "$HOME/GSE100075/assembly/$id-ballgown.gtf" | grep -P 'ANKLE2'  | grep -P 'transcript\t' | sed -r "s/StringTie/$id/"; done > GSE100075/GSE100075-ANKLE2.gtf
```


## 4. plot in R

`GSE100075-biosample_result.xml` comes from https://www.ncbi.nlm.nih.gov/biosample?Db=biosample&DbFrom=bioproject&Cmd=Link&LinkName=bioproject_biosample&LinkReadableName=BioSample&ordinalpos=1&IdsFromResult=390636

```r
sample_table <- inner_join(
	x <- read_tsv('zhu/zhu/PRJNA390636.txt') %>% select(sample_accession, run_accession),
	y <- xml2::read_xml('zhu/GSE100075-biosample_result.xml') %>% 
        xml2::as_list() %>% {.[[1]]} %>% 
        sapply(. %>% {c(id = .$Ids$Id[[1]], description = .$Description$Title[[1]])}) %>% 
		t %>% as_tibble(),
	by = c('sample_accession' = 'id')
) %>% select(sample = run_accession, description) %T>% print;
#> # A tibble: 15 x 2
#>    sample     description               
#>    <chr>      <chr>                     
#>  1 SRR5687500 MCF7_WT_repl1             
#> ...       
#> 15 SRR5687514 SUM44_LTED_repl3   

df <- read_lines('zhu/GSE100075-ANKLE2.gtf') %>% {tibble(
		sample = str_extract(., 'SRR\\d+'), 
		transcript = str_extract(., 'transcript_id "\\.\\w+'), 
		FPKM = str_extract(., '(?<=FPKM ")[\\.\\d]+') %>% as.numeric, 
		TPM = str_extract(., '(?<=TPM ")[\\.\\d]+') %>% as.numeric
)} %>% inner_join(sample_table) %T>% print;
#> # A tibble: 180 x 5
#>    sample     transcript     FPKM      TPM description  
#>    <chr>      <chr>         <dbl>    <dbl> <chr>        
#>  1 SRR5687500 NA         0        0        MCF7_WT_repl1
#>  2 SRR5687500 NA         0.998    1.47     MCF7_WT_repl1
#>  ...
```

```r
# 总体来看，FPKM 与 TPM 有很强的线性关系
ggplot(df, aes(FPKM, TPM)) + geom_point() + geom_smooth(method = "lm")
```
![](zhu/FPKM-TPM-GSE100075.png)

```r
# In every sample, FPKM 与 TPM 完美线性相关
ggplot(df, aes(FPKM, TPM)) + geom_point() + geom_smooth(method = "lm") + facet_wrap(~sample, scales = 'free')
ggsave('FPKM-TPM-per-sample-GSE100075.png', width = 6, height = 6)
```
![](zhu/FPKM-TPM-per-sample-GSE100075.png)

```r
df %>% group_by(sample) %>% 
    dplyr::summarise(TPM = sum(TPM), description = description[1]) %>%
	mutate(genotype = str_remove(description, '_repl\\d')) %>% 
    group_by(genotype) %>% summarise(TPM_mean = mean(TPM), TPM_std = sd(TPM)) %>%
	mutate(cell_line = str_extract(genotype, '\\w+?(?=_)'))  %T>% print %>%
	ggplot(aes(genotype, TPM_mean, fill = cell_line)) + 
        geom_bar(position=position_dodge(), stat = "identity") +
        geom_errorbar(aes(ymin = TPM_mean - TPM_std, ymax = TPM_mean + TPM_std), width = .2, position = position_dodge(.9)) + 
        labs(y = 'TPM') +
        coord_flip()
        theme(axis.text.x = element_text(angle = 60))
#> # A tibble: 5 x 4
#>   genotype             TPM_mean TPM_std cell_line
#>   <chr>                   <dbl>   <dbl> <chr>    
#> 1 MCF7_LTED_ESR1_WT       23.3    4.44  MCF7     
#> ... 
#> 5 SUM44_WT                15.8    2.99  SUM44    
```
![](zhu/ANKLE-GSE100075.png)






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
matrix_raw <- readr::read_lines('GSE/GSE2034_series_matrix.txt.gz')
info <- matrix_raw %>% {
    tibble(
	    sample = str_subset(., '^!Sample_geo_accession') %>% 
            str_extract_all('GSM\\d+') %>% {.[[1]]},
	    relapse = str_subset(., 'bone relapses') %>% 
            str_extract_all('(?<=: )[01]') %>% {.[[1]] == '1'}
    )
} %T>% print

matrix <- qGSEA::read_gse_matrix('GSE/GSE2034_series_matrix.txt.gz')
chip <- qGSEA::read_gse_soft('GSE/GSE2034_family.soft.gz')
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
cat GSE/GSE4922_family.soft | grep -P '(?:ANKLE2|AURKA)\s' > GSE/GSE4922-ANKLE2_AURKA.soft
cat GSE/GSE4922-GPL96_series_matrix.txt | grep -P '^!Sample_' > GSE/GSE4922-GPL96-sample.info
```

```r
soft2 <- read_lines('GSE/GSE4922-ANKLE2_AURKA.soft') %>% 
	{tibble(probe = str_extract(., '\\w+'), symbol = str_extract(., 'ANKLE2|AURKA'))}
matrix_raw2 <- read_lines('GSE/GSE4922-GPL96_series_matrix.txt.gz')

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

matrix <- qGSEA::read_gse_matrix('GSE/GSE4922-GPL96_series_matrix.txt.gz')
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


