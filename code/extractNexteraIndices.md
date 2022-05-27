Generate barcode sheets
================

``` r
library(tidyverse)
```

    ## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.1 ──

    ## ✔ ggplot2 3.3.6     ✔ purrr   0.3.4
    ## ✔ tibble  3.1.7     ✔ dplyr   1.0.9
    ## ✔ tidyr   1.2.0     ✔ stringr 1.4.0
    ## ✔ readr   2.1.2     ✔ forcats 0.5.1

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()

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
Nextera_i7_Primers <- readxl::read_excel("../data/41586_2015_BFnature14590_MOESM36_ESM.xlsx") %>% 
  rename(Nextera_i7 = 1, i7_Primer = 2) %>%
  filter(str_detect(Nextera_i7, "v2_Ad2\\.")) %>%
  separate(Nextera_i7, into = c("id1", "id2", "i7_Index"), sep = "_", remove = FALSE) %>%
  unite(col = Nextera_i7, id1:id2, sep = "_") %>%
  mutate(i7_Adapter_Sequence = i7_Index) %>%
  mutate(i7_Adapter_Sequence = dna(i7_Adapter_Sequence)) %>%
  mutate(i7_Adapter_Sequence = seq_complement(i7_Adapter_Sequence)) %>%
  mutate(i7_Adapter_Sequence = seq_reverse((i7_Adapter_Sequence))) %>%
  select(Nextera_i7, i7_Index, i7_Adapter_Sequence, i7_Primer) %>%
  write_csv(file = "../barcodes/Nextera_i7_Adapters.csv")

Nextera_i5_Primers <- readxl::read_excel("../data/41586_2015_BFnature14590_MOESM36_ESM.xlsx") %>% 
  rename(Nextera_i5 = 1, i5_Primer = 2) %>%
  filter(str_detect(Nextera_i5, "v2_Ad1\\.")) %>%
  separate(Nextera_i5, into = c("id1", "id2", "i5_Adapter_Sequence"), sep = "_", remove = FALSE) %>%
  unite(col = Nextera_i5, id1:id2, sep = "_") %>%
  mutate(i5_Index = i5_Adapter_Sequence) %>%
  mutate(i5_Index = dna(i5_Index)) %>%
  mutate(i5_Index = seq_complement(i5_Index)) %>%
  mutate(i5_Index = seq_reverse((i5_Index))) %>%
  select(Nextera_i5, i5_Index, i5_Adapter_Sequence, i5_Primer)  %>%
  write_csv(file = "../barcodes/Nextera_i5_Adapters.csv")
```
