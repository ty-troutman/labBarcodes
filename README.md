README
================

Generate barcode sheet for Nextera PCR primers from **[Supplementary
Table
1](https://static-content.springer.com/esm/art%3A10.1038%2Fnature14590/MediaObjects/41586_2015_BFnature14590_MOESM36_ESM.xlsx)**:  
[Buenrostro JD, Wu B, Litzenburger UM, Ruff D, Gonzales ML, Snyder MP,
Chang HY, Greenleaf WJ. Single-cell chromatin accessibility reveals
principles of regulatory variation. Nature. 2015 Jul
23;523(7561):486-90. doi: 10.1038/nature14590. Epub 2015 Jun 17. PMID:
26083756; PMCID: PMC4685948.](https://pubmed.ncbi.nlm.nih.gov/26083756/)

``` r
Nextera_I7_Primers <- readxl::read_excel("data/41586_2015_BFnature14590_MOESM36_ESM.xlsx") %>% 
  rename(Index_ID = 1, FullPrimer = 2) %>%
  filter(str_detect(Index_ID, "v2_Ad2\\.")) %>%
  separate(Index_ID, into = c("id1", "id2", "I7_Primer_Sequence"), sep = "_", remove = FALSE) %>%
  unite(col = Index_ID, id1:id2, sep = "_") %>%
  mutate(I7_Index = I7_Primer_Sequence) %>%
  mutate(I7_Index = dna(I7_Index)) %>%
  mutate(I7_Index = seq_complement(I7_Index)) %>%
  mutate(I7_Index = seq_reverse((I7_Index))) %>%
  select(Index_ID, I7_Index, I7_Primer_Sequence, FullPrimer) %>%
  write_csv(file = "barcodes/Nextera_I7_Primers.csv")

Nextera_I5_Primers <- readxl::read_excel("data/41586_2015_BFnature14590_MOESM36_ESM.xlsx") %>% 
  rename(Index_ID = 1, FullPrimer = 2) %>%
  filter(str_detect(Index_ID, "v2_Ad1\\.")) %>%
  separate(Index_ID, into = c("id1", "id2", "I5_Primer_Sequence"), sep = "_", remove = FALSE) %>%
  unite(col = Index_ID, id1:id2, sep = "_") %>%
  mutate(I5_Reverse_Complement = I5_Primer_Sequence) %>%
  mutate(I5_Reverse_Complement = dna(I5_Reverse_Complement)) %>%
  mutate(I5_Reverse_Complement = seq_complement(I5_Reverse_Complement)) %>%
  mutate(I5_Reverse_Complement = seq_reverse((I5_Reverse_Complement))) %>%
  select(Index_ID, I5_Primer_Sequence, I5_Reverse_Complement, FullPrimer) %>%
  write_csv(file = "barcodes/Nextera_I5_Primers.csv")
```

### Compare Hamming Distance of I7 Nextera indexes to I7 NextFlex 289:384 indexes

![I7 hamming distance &gt; 2](plots/i7_hamming_greater2.png)

![I7 hamming distance &gt; 2](plots/i7_hamming_greater2_unclustered.png)

### Compare Hamming Distance of I5 Nextera indexes to I5 NextFlex 289:384 indexes

![I5 hamming distance &gt; 2](plots/i5_hamming_greater2.png)

![I5 hamming distance &gt; 2](plots/i5_hamming_greater2_unclustered.png)
