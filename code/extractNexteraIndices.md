Generate barcode sheets
================

``` r
library(tidyverse)
```

    ## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.1 ──

    ## ✓ ggplot2 3.3.5     ✓ purrr   0.3.4
    ## ✓ tibble  3.1.6     ✓ dplyr   1.0.7
    ## ✓ tidyr   1.1.4     ✓ stringr 1.4.0
    ## ✓ readr   2.1.1     ✓ forcats 0.5.1

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

``` r
library(bioseq)
library(DNABarcodes)
```

    ## Loading required package: Matrix

    ## 
    ## Attaching package: 'Matrix'

    ## The following objects are masked from 'package:tidyr':
    ## 
    ##     expand, pack, unpack

    ## Loading required package: parallel

``` r
Nextera_I7_Primers <- readxl::read_excel("../data/41586_2015_BFnature14590_MOESM36_ESM.xlsx") %>% 
  rename(Index_ID = 1, FullPrimer = 2) %>%
  filter(str_detect(Index_ID, "v2_Ad2\\.")) %>%
  separate(Index_ID, into = c("id1", "id2", "I7_Primer_Sequence"), sep = "_", remove = FALSE) %>%
  unite(col = Index_ID, id1:id2, sep = "_") %>%
  mutate(I7_Index = I7_Primer_Sequence) %>%
  mutate(I7_Index = dna(I7_Index)) %>%
  mutate(I7_Index = seq_complement(I7_Index)) %>%
  mutate(I7_Index = seq_reverse((I7_Index))) %>%
  select(Index_ID, I7_Index, I7_Primer_Sequence, FullPrimer) %>%
  write_csv(file = "../barcodes/Nextera_I7_Primers.csv")

Nextera_I5_Primers <- readxl::read_excel("../data/41586_2015_BFnature14590_MOESM36_ESM.xlsx") %>% 
  rename(Index_ID = 1, FullPrimer = 2) %>%
  filter(str_detect(Index_ID, "v2_Ad1\\.")) %>%
  separate(Index_ID, into = c("id1", "id2", "I5_Primer_Sequence"), sep = "_", remove = FALSE) %>%
  unite(col = Index_ID, id1:id2, sep = "_") %>%
  mutate(I5_Reverse_Complement = I5_Primer_Sequence) %>%
  mutate(I5_Reverse_Complement = dna(I5_Reverse_Complement)) %>%
  mutate(I5_Reverse_Complement = seq_complement(I5_Reverse_Complement)) %>%
  mutate(I5_Reverse_Complement = seq_reverse((I5_Reverse_Complement))) %>%
  select(Index_ID, I5_Primer_Sequence, I5_Reverse_Complement, FullPrimer) %>%
  write_csv(file = "../barcodes/Nextera_I5_Primers.csv")
```
