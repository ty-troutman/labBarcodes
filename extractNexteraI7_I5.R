library(DNABarcodes)
library(tidyverse)
library(pheatmap)
library(viridis)

nextflex <- readxl::read_excel("UDI_Index_Sequences_v21.06_v1-New.xlsx") %>% 
  rename(barcode = 1, i7 = 2, i5 = 3, i5RevComp = 4) 

buenrostroi7 <- readxl::read_excel("Buenrostro_ATACSeq_2015/41586_2015_BFnature14590_MOESM36_ESM.xlsx") %>% 
  rename(barcode = 1, primer = 2) %>%
  filter(str_detect(barcode, "v2_Ad2")) %>%
  separate(barcode, into = c("drop", "drop2", "id", "i7"), remove = F) %>%
  filter(!is.na(i7))

buenrostroi5 <- readxl::read_excel("Buenrostro_ATACSeq_2015/41586_2015_BFnature14590_MOESM36_ESM.xlsx") %>% 
  rename(barcode = 1, primer = 2) %>%
  filter(str_detect(barcode, "v2_Ad1")) %>%
  separate(barcode, into = c("drop", "drop2", "id", "i5"), remove = F) %>%
  filter(!is.na(i5))

i7 <-
  bind_rows(
    nextflex %>% select(barcode, i7) %>% slice(97:192),
    buenrostroi7 %>% select(barcode, i7)
  )

tmp <- barcode.set.distances(pull(i7, i7), metric = c("hamming"))
rownames(tmp) <- i7$barcode
colnames(tmp) <- i7$barcode
x1 <- tmp@x
nameX <- tmp@Dimnames
nameY <- tmp@Dimnames
tmp <- as.matrix(tmp, row.names = tmp@Dimnames)
tmp <- as.data.frame(tmp)
tmp <- tmp %>% rownames_to_column(., var = "rowname") %>% as_tibble()
tmp2 <- tmp %>% pivot_longer(cols = -1, names_to = "barcode", values_to = "hamming")

tmp <- barcode.set.distances(pull(i7, i7), metric = c("seqlev"))
rownames(tmp) <- i7$barcode
colnames(tmp) <- i7$barcode
x1 <- tmp@x
nameX <- tmp@Dimnames
nameY <- tmp@Dimnames
tmp <- as.matrix(tmp, row.names = tmp@Dimnames)
tmp <- as.data.frame(tmp)
tmp <- tmp %>% rownames_to_column(., var = "rowname") %>% as_tibble()
tmp <- tmp %>% pivot_longer(cols = -1, names_to = "barcode", values_to = "seqlev")

i7 <- full_join(tmp, tmp2)

i7_tmp <-
  i7 %>% sample_frac(size = 1, replace = FALSE) %>% 
  filter(str_detect(rowname, "UDI")) %>%
  filter(str_detect(barcode, "v2")) %>%
  select(-hamming) %>%
  filter(seqlev > 2) %>% 
  pivot_wider(names_from = barcode, values_from = seqlev) %>%
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
    n = 17,
    direction = -1,
    option = "rocket",
    end = 0.8
  ),
  clustering_method = "ward.D",
  na_col = "white",
  main = "I7 barcode distances, seqlev > 2",
  # filename = "i7_hammingGreater2.pdf",
  # width = 10,
  # height = 7.5
)

i5 <-
  bind_rows(
  nextflex %>% select(barcode, i5) %>% slice(97:192),
  buenrostroi5 %>% select(barcode, i5)
)

tmp <- barcode.set.distances(pull(i5, i5), metric = c("hamming"))
rownames(tmp) <- i5$barcode
colnames(tmp) <- i5$barcode
x1 <- tmp@x
nameX <- tmp@Dimnames
nameY <- tmp@Dimnames
tmp <- as.matrix(tmp, row.names = tmp@Dimnames)
tmp <- as.data.frame(tmp)
tmp <- tmp %>% rownames_to_column(., var = "rowname") %>% as_tibble()
tmp2 <- tmp %>% pivot_longer(cols = -1, names_to = "barcode", values_to = "hamming")

tmp <- barcode.set.distances(pull(i5, i5), metric = c("seqlev"))
rownames(tmp) <- i5$barcode
colnames(tmp) <- i5$barcode
x1 <- tmp@x
nameX <- tmp@Dimnames
nameY <- tmp@Dimnames
tmp <- as.matrix(tmp, row.names = tmp@Dimnames)
tmp <- as.data.frame(tmp)
tmp <- tmp %>% rownames_to_column(., var = "rowname") %>% as_tibble()
tmp <- tmp %>% pivot_longer(cols = -1, names_to = "barcode", values_to = "seqlev")

i5 <- full_join(tmp, tmp2)

tmp_i5 <-
  i5 %>% sample_frac(size = 1, replace = FALSE) %>% 
  filter(str_detect(rowname, "UDI")) %>%
  filter(str_detect(barcode, "v2")) %>%
  # filter(rowname != barcode) %>%
  select(-hamming) %>%
  filter(seqlev > 2) %>% 
  pivot_wider(names_from = barcode, values_from = seqlev) %>%
  arrange(rowname) %>%
  column_to_rownames(var = "rowname")

pheatmap(
  tmp_i5,
  cluster_rows = T,
  cluster_cols = T,
  fontsize = 10,
  fontsize_row = 6,
  fontsize_col = 6,
  border_color = NA,
  color = viridis(
    n = 17,
    direction = -1,
    option = "rocket",
    end = 0.8
  ),
  clustering_method = "ward.D",
  na_col = "white",
  main = "I5 barcode distances, seqlev > 2",
  # filename = "i7_hammingGreater2.pdf",
  # width = 10,
  # height = 7.5
)


i7 %>% 
  pivot_longer(cols = seqlev:hamming, names_to = "metric", values_to = "distance") %>%
  filter(str_detect(rowname, "UDI")) %>%
  filter(str_detect(barcode, "v2")) %>%
  ggplot(aes(x = rowname, y = barcode, fill = distance)) +
  geom_raster() +
  scale_fill_viridis(option = "rocket", direction  = -1, end = 0.8) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = 6),
    axis.text.y = element_text(size = 6)
  ) +
  facet_wrap(~ metric)

i7 %>% 
  filter(str_detect(rowname, "UDI")) %>%
  filter(str_detect(barcode, "v2")) %>%
  mutate(diff = hamming/seqlev) %>%
  ggplot(aes(x = rowname, y = barcode, fill = diff)) +
  geom_raster() +
  scale_fill_viridis(option = "rocket", direction  = -1, end = 0.8) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = 6),
    axis.text.y = element_text(size = 6)
  ) 
