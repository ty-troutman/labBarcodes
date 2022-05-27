Hamming Distance Calculations
================

Setup environment

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
library(pheatmap)
library(viridis)
```

    ## Loading required package: viridisLite

Hamming Distance of I7 Nextera indexes to I7 NextFlex 289:384 indexes

``` r
Nextera_I7_Primers <- read_csv("../barcodes/Nextera_i7_Adapters.csv") %>%
  rename(Index_ID = 1, I7_Index = 2)
```

    ## Rows: 96 Columns: 4
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr (4): Nextera_i7, i7_Index, i7_Adapter_Sequence, i7_Primer
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
Nextera_I5_Primers <- read_csv("../barcodes/Nextera_i5_Adapters.csv") %>%
  rename(Index_ID = 1, I5_Index = 2)
```

    ## Rows: 92 Columns: 4
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr (4): Nextera_i5, i5_Index, i5_Adapter_Sequence, i5_Primer
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
NextFlex <- read_csv("../barcodes/NEXTFLEX_Adapters.csv") %>%
  rename(Index_ID = 1, I7_Index = 2, I5_Index = 3)
```

    ## Rows: 384 Columns: 4
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr (4): NEXTFLEX_ID, i7_Index, i5_Index, i5_Adapter_Sequence
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
I7_Merged <-
  bind_rows(
  Nextera_I7_Primers %>% select(Index_ID, I7_Index),
  NextFlex %>% slice(289:384) %>% select(Index_ID, I7_Index)
)

tmp <- barcode.set.distances(pull(I7_Merged, I7_Index), metric = c("hamming"))
rownames(tmp) <- I7_Merged$Index_ID
colnames(tmp) <- I7_Merged$Index_ID
x1 <- tmp@x
nameX <- tmp@Dimnames
nameY <- tmp@Dimnames
tmp <- as.matrix(tmp, row.names = tmp@Dimnames)
tmp <- as.data.frame(tmp)
tmp <- tmp %>% rownames_to_column(., var = "rowname") %>% as_tibble()
tmp2 <- tmp %>% pivot_longer(cols = -1, names_to = "barcode", values_to = "hamming")

tmp <- barcode.set.distances(pull(I7_Merged, I7_Index), metric = c("seqlev"))
rownames(tmp) <- I7_Merged$Index_ID
colnames(tmp) <- I7_Merged$Index_ID
x1 <- tmp@x
nameX <- tmp@Dimnames
nameY <- tmp@Dimnames
tmp <- as.matrix(tmp, row.names = tmp@Dimnames)
tmp <- as.data.frame(tmp)
tmp <- tmp %>% rownames_to_column(., var = "rowname") %>% as_tibble()
tmp <- tmp %>% pivot_longer(cols = -1, names_to = "barcode", values_to = "seqlev")

i7 <- full_join(tmp, tmp2)
```

    ## Joining, by = c("rowname", "barcode")

``` r
i7 %>% sample_frac(size = 1, replace = FALSE) %>% 
  filter(str_detect(rowname, "UDI")) %>%
  filter(str_detect(barcode, "v2")) %>%
  rename(NEXTFLEX_ID = 1, Nextera_i7 = 2) %>%
  write_csv("../indexDistances/i7_distances.csv")

i7_tmp <-
  i7 %>% sample_frac(size = 1, replace = FALSE) %>% 
  filter(str_detect(rowname, "UDI")) %>%
  filter(str_detect(barcode, "v2")) %>%
  select(-seqlev) %>%
  filter(hamming > 2) %>%
  pivot_wider(names_from = barcode, values_from = hamming) %>%
  arrange(rowname) %>%
  column_to_rownames(var = "rowname")

pheatmap(
  i7_tmp,
  cluster_rows = T,
  cluster_cols = T,
  fontsize = 10,
  fontsize_row = 6,
  fontsize_col = 6,
  border_color = NA,
  color = viridis(
    n = 6,
    direction = -1,
    option = "rocket",
    end = 0.8
  ),
  clustering_method = "ward.D",
  na_col = "white",
  main = "I7 barcode Hamming distances, Hamming > 2",
  filename = "../plots/i7_hamming_greater2.pdf",
  width = 10,
  height = 7.5
)

pheatmap(
  i7_tmp,
  cluster_rows = F,
  cluster_cols = F,
  fontsize = 10,
  fontsize_row = 6,
  fontsize_col = 6,
  border_color = NA,
  color = viridis(
    n = 6,
    direction = -1,
    option = "rocket",
    end = 0.8
  ),
  na_col = "white",
  main = "I7 barcode Hamming distances, Hamming > 2",
  filename = "../plots/i7_hamming_greater2_unclustered.pdf",
  width = 10,
  height = 7.5
)

pheatmap(
  i7_tmp,
  cluster_rows = T,
  cluster_cols = T,
  fontsize = 10,
  fontsize_row = 6,
  fontsize_col = 6,
  border_color = NA,
  color = viridis(
    n = 6,
    direction = -1,
    option = "rocket",
    end = 0.8
  ),
  clustering_method = "ward.D",
  na_col = "white",
  main = "I7 barcode Hamming distances, Hamming > 2",
  filename = "../plots/i7_hamming_greater2.png",
  width = 10,
  height = 7.5
)

pheatmap(
  i7_tmp,
  cluster_rows = F,
  cluster_cols = F,
  fontsize = 10,
  fontsize_row = 6,
  fontsize_col = 6,
  border_color = NA,
  color = viridis(
    n = 6,
    direction = -1,
    option = "rocket",
    end = 0.8
  ),
  na_col = "white",
  main = "I7 barcode Hamming distances, Hamming > 2",
  filename = "../plots/i7_hamming_greater2_unclustered.png",
  width = 10,
  height = 7.5
)
```

Compare Hamming Distance of I5 Nextera indexes to I5 NextFlex 289:384
indexes

``` r
I5_Merged <-
  bind_rows(
  Nextera_I5_Primers %>% select(Index_ID, I5_Index),
  NextFlex %>% slice(289:384) %>% select(Index_ID, I5_Index)
)

tmp <- barcode.set.distances(pull(I5_Merged, I5_Index), metric = c("hamming"))
rownames(tmp) <- I5_Merged$Index_ID
colnames(tmp) <- I5_Merged$Index_ID
x1 <- tmp@x
nameX <- tmp@Dimnames
nameY <- tmp@Dimnames
tmp <- as.matrix(tmp, row.names = tmp@Dimnames)
tmp <- as.data.frame(tmp)
tmp <- tmp %>% rownames_to_column(., var = "rowname") %>% as_tibble()
tmp2 <- tmp %>% pivot_longer(cols = -1, names_to = "barcode", values_to = "hamming")

tmp <- barcode.set.distances(pull(I5_Merged, I5_Index), metric = c("seqlev"))
rownames(tmp) <- I5_Merged$Index_ID
colnames(tmp) <- I5_Merged$Index_ID
x1 <- tmp@x
nameX <- tmp@Dimnames
nameY <- tmp@Dimnames
tmp <- as.matrix(tmp, row.names = tmp@Dimnames)
tmp <- as.data.frame(tmp)
tmp <- tmp %>% rownames_to_column(., var = "rowname") %>% as_tibble()
tmp <- tmp %>% pivot_longer(cols = -1, names_to = "barcode", values_to = "seqlev")

i5 <- full_join(tmp, tmp2)
```

    ## Joining, by = c("rowname", "barcode")

``` r
i5 %>%
  filter(str_detect(rowname, "UDI")) %>%
  filter(str_detect(barcode, "v2")) %>%
  rename(NEXTFLEX_ID = 1, Nextera_i5 = 2) %>%
  write_csv("../indexDistances/i5_distances.csv")

i5_tmp <-
  i5 %>%
  filter(str_detect(rowname, "UDI")) %>%
  filter(str_detect(barcode, "v2")) %>%
  select(-seqlev) %>%
  filter(hamming > 2) %>%
  pivot_wider(names_from = barcode, values_from = hamming) %>%
  arrange(rowname) %>%
  column_to_rownames(var = "rowname")

pheatmap(
  i5_tmp,
  cluster_rows = T,
  cluster_cols = T,
  fontsize = 10,
  fontsize_row = 6,
  fontsize_col = 6,
  border_color = NA,
  color = viridis(
    n = 6,
    direction = -1,
    option = "rocket",
    end = 0.8
  ),
  clustering_method = "ward.D",
  na_col = "white",
  main = "I5 barcode Hamming distances, Hamming > 2",
  filename = "../plots/i5_hamming_greater2.pdf",
  width = 10,
  height = 7.5
)

pheatmap(
  i5_tmp,
  cluster_rows = F,
  cluster_cols = F,
  fontsize = 10,
  fontsize_row = 6,
  fontsize_col = 6,
  border_color = NA,
  color = viridis(
    n = 6,
    direction = -1,
    option = "rocket",
    end = 0.8
  ),
  na_col = "white",
  main = "I5 barcode Hamming distances, Hamming > 2",
  filename = "../plots/i5_hamming_greater2_unclustered.pdf",
  width = 10,
  height = 7.5
)

pheatmap(
  i5_tmp,
  cluster_rows = T,
  cluster_cols = T,
  fontsize = 10,
  fontsize_row = 6,
  fontsize_col = 6,
  border_color = NA,
  color = viridis(
    n = 6,
    direction = -1,
    option = "rocket",
    end = 0.8
  ),
  clustering_method = "ward.D",
  na_col = "white",
  main = "I5 barcode Hamming distances, Hamming > 2",
  filename = "../plots/i5_hamming_greater2.png",
  width = 10,
  height = 7.5
)

pheatmap(
  i5_tmp,
  cluster_rows = F,
  cluster_cols = F,
  fontsize = 10,
  fontsize_row = 6,
  fontsize_col = 6,
  border_color = NA,
  color = viridis(
    n = 6,
    direction = -1,
    option = "rocket",
    end = 0.8
  ),
  na_col = "white",
  main = "I5 barcode Hamming distances, Hamming > 2",
  filename = "../plots/i5_hamming_greater2_unclustered.png",
  width = 10,
  height = 7.5
)
```

## Finished!
