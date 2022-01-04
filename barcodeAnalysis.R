library(DNABarcodes)
library(tidyverse)
library(bioseq)
library(pheatmap)
library(viridis)
library(cowplot)

buenrostro_i5 <- readxl::read_excel("41586_2015_BFnature14590_MOESM36_ESM.xlsx") %>% 
  rename(barcode = 1, primer = 2) %>%
  filter(str_detect(barcode, "v2_Ad1")) %>%
  separate(barcode, into = c("drop", "drop2", "id", "i5"), remove = FALSE) %>%
  filter(!is.na(i5))

buenrostro_i7 <- readxl::read_excel("41586_2015_BFnature14590_MOESM36_ESM.xlsx") %>% 
  rename(barcode = 1, primer = 2) %>%
  filter(str_detect(barcode, "v2_Ad2")) %>%
  separate(barcode, into = c("drop", "drop2", "id", "i7"), remove = FALSE) %>%
  filter(!is.na(i7))

nextflex <- readxl::read_excel("../UDI_Index_Sequences_v21.06_v1-New.xlsx") %>% 
  rename(barcode = 1, i7 = 2, i5 = 3, i5RevComp = 4) 

merge <-
  bind_rows(
    nextflex %>% select(barcode, i7) %>% slice(1:192),
    buenrostro_i7 %>% select(barcode, i7)
  )

analyse.barcodes(pull(merge, i7))

setDistance <- pull(merge, i7) %>%
  barcode.set.distances(barcodes = ., metric = "hamming")

rownames(setDistance) <- merge$barcode
colnames(setDistance) <- merge$barcode

x1 <- setDistance@x
nameX <- setDistance@Dimnames
nameY <- setDistance@Dimnames

tmp_ <- as.matrix(setDistance, row.names = setDistance@Dimnames)
tmp_ <- as.data.frame(tmp_)
i7_hamming <- tmp_ %>% rownames_to_column(., var = "rowname") %>% as_tibble()

i7Plot <- i7_hamming %>% pivot_longer(cols = -1, names_to = "barcode", values_to = "hamming") %>%
  filter(str_detect(rowname, "UDI")) %>% filter(str_detect(barcode, "v2")) %>%
  ggplot(aes(x = rowname, y = barcode, fill = hamming)) + geom_tile() +
  scale_fill_viridis(direction = 1, option = "rocket") + 
  coord_fixed() +
  xlab("Nextflex i7") + ylab("Buenrostro i7") + labs(fill = "Hamming Distance") +
  theme(axis.text.y = element_text(size = 2), 
        axis.text.x = element_text(size = 2, angle = 90, hjust = 1))

merge <-
  bind_rows(
    nextflex %>% select(barcode, i5) %>% slice(1:192),
    buenrostro_i5 %>% select(barcode, i5)
  )

analyse.barcodes(pull(merge, i5))

setDistance <- pull(merge, i5) %>%
  barcode.set.distances(barcodes = ., metric = "hamming")

rownames(setDistance) <- merge$barcode
colnames(setDistance) <- merge$barcode

x1 <- setDistance@x
nameX <- setDistance@Dimnames
nameY <- setDistance@Dimnames

tmp_ <- as.matrix(setDistance, row.names = setDistance@Dimnames)
tmp_ <- as.data.frame(tmp_)
i5_hamming <- tmp_ %>% rownames_to_column(., var = "rowname") %>% as_tibble()

i5Plot <- i5_hamming %>% pivot_longer(cols = -1, names_to = "barcode", values_to = "hamming") %>%
  filter(str_detect(rowname, "UDI")) %>% filter(str_detect(barcode, "v2")) %>%
  ggplot(aes(x = rowname, y = barcode, fill = hamming)) + geom_tile() +
  scale_fill_viridis(direction = 1, option = "rocket") + 
  coord_fixed() +
  xlab("Nextflex i5") + ylab("Buenrostro i5") + labs(fill = "Hamming Distance") +
  theme(axis.text.y = element_text(size = 2), 
        axis.text.x = element_text(size = 2, angle = 90, hjust = 1))

ggsave(filename = "hammingDistances_Nextflex_BuenrostroATAC.pdf", 
       plot = plot_grid(i7Plot, i5Plot, nrow = 1, axis = l),
       height = 8.5, width = 11, device = cairo_pdf
)

i5_hamming %>% pivot_longer(cols = -1, names_to = "barcode", values_to = "hamming") %>%
  filter(str_detect(rowname, "UDI")) %>% filter(str_detect(barcode, "v2")) %>%
  filter(hamming < 3) %>% distinct(barcode) %>% print(n = 100)

i7_hamming %>% pivot_longer(cols = -1, names_to = "barcode", values_to = "hamming") %>%
  filter(str_detect(rowname, "UDI")) %>% filter(str_detect(barcode, "v2")) %>%
  filter(hamming < 3) %>% distinct(barcode) %>% print(n = 100)
