suppressMessages(library(tidyverse))
library(doMC)
registerDoMC(10)

# --- Collect Files ---

snp.f <- dir("varscan", ".snp$", full.names = TRUE)
name <- snp.f %>% stringr::str_extract("chRCC\\d+")
snp <-
    foreach(s = 1:length(snp.f), .combine = "bind_rows") %do% {
        o <- read_tsv(snp.f[s], col_names=TRUE, col_types=cols(.default="c"), comment="##")
        o$case = name[s]
        return(o)
    }

indel.f <- dir("varscan", ".indel$", full.names = TRUE)
name <- indel.f %>% stringr::str_extract("chRCC\\d+")
indel <-
    foreach(s = 1:length(indel.f), .combine = "bind_rows") %do% {
        o <- read_tsv(indel.f[s], col_names=TRUE, col_types=cols(.default="c"), comment="##")
        o$case = name[s]
        return(o)
    }
indel

rm(o, snp.f, indel.f, name, s)

# --- cleanup header ---
# --- snp ---
snp$VAF <- snp$Sample1 %>% stringr::str_split_fixed(":", 8) %>% .[,7] %>% stringr::str_replace("%", "") %>% as.numeric() %>% "/"(100) %>% round(4)
x <- snp$Annot %>% stringr::str_extract("[^;]+(?=\\;)") %>% stringr::str_split_fixed(",", 8) %>% .[ , c(1:3, 6)] %>% as.data.frame(stringsAsFactors = FALSE)
colnames(x) <- c("category", "gene", "strand", "aa")
snp <- snp %>% bind_cols(x)
rm(x)
snp <- select(snp, `#CHROM`:QUAL, case:strand, aa, Nearest_5p, Nearest_3p, everything())
snp %>% write_tsv("varscan/all.snps")
# --- indel ---
indel$VAF <- indel$Sample1 %>% stringr::str_split_fixed(":", 8) %>% .[,7] %>% stringr::str_replace("%", "") %>% as.numeric() %>% "/"(100) %>% round(4)
indel <- indel %>% select(`#CHROM`:QUAL, case, VAF, promoter:other_info, everything())
indel %>% write_tsv("varscan/all.indels")


# 
snp %>% filter(VAF < .97, category=="nsSNP") %>% View
indel %>% filter(cds!="c")
save(snp, indel, file = "varscan/all.Rdata")
