#!/home/users/kjyi/tools/R/R-3.4.0/bin/Rscript
library(tidyverse)
library(doMC)
registerDoMC(20)

# samtools view -c
chr <- c(1:22, "X", "Y", "MT")
sample <- dir("bamfiles", "bam$", full.names = TRUE)
read_count_by_chromosome <- 
    foreach(l = chr, .combine = cbind) %:%
    foreach(s = sample, .combine = rbind) %dopar% {
    cmd <- paste("samtools view -c", s, l)
    x <- system(cmd, intern = TRUE)
    return(as.integer(x))
    # return(cmd)
} %>% as.data.frame()
colnames(read_count_by_chromosome) <- chr
rownames(read_count_by_chromosome) <- stringr::str_extract(sample, "chRCC\\d+")

rownames(read_count_by_chromosome) <- rownames(read_count_by_chromosome) %>% stringr::str_extract("\\d+") %>% stringr::str_pad(2, pad = "0") %>% paste0("chRCC", .)

read_count_by_chromosome$autosome <- read_count_by_chromosome %>%
    select(-X, -Y, -MT) %>%
    rowSums()

read_count_by_chromosome <- read_count_by_chromosome %>% 
    rownames_to_column("case") %>%
    arrange(desc(case)) %>%
    mutate(mt_copy = read_count_by_chromosome$MT/read_count_by_chromosome$autosome/16969*2860000000*2)
read_count_by_chromosome %>% write_tsv("read_count_by_chromosome.tsv")
read_count_by_chromosome %>%
    select(case, mt_copy) %>% column_to_rownames("case") %>% as.matrix() %>% t() %>%
    barplot(horiz = T, las = 2)
