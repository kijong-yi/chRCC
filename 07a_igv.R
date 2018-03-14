suppressMessages(library(tidyverse))
library(stringr)
library(doMC); registerDoMC(20)

csv <- system("ls mtoolbox/*/*annotation.csv", intern = TRUE)
annotated <- foreach(c = csv, .combine = bind_rows) %do% {
    read_tsv(c) %>%
        filter(HF > 0.03) %>%
        filter(is.na(Haplogroup))  %>%
        select(-RSRS, -MHCS, -rCRS, -Haplogroup)
}

getRange <- function(v) {
    x <- stringr::str_extract(v, "\\d+") %>%
        as.numeric()
    y <- paste0((x - 342), "-", (x + 341))
    y[which(x < 343)] <- "1-684"
    y[which(x > 16288)] <- "15886-16569"
    return(y)
    }

batchlist <- annotated %>%
    filter(HF > .10) %>%
    filter(RSRS == "yes", MHCS == "yes", rCRS == "yes") %>%
    filter(!`Aa Change` %in% c("syn", NA)) %>%
    transmute(s0 = ifelse(Sample != lag(Sample, default = "a"), "new\n", ""),
              s1 = ifelse(Sample != lag(Sample, default = "a"), 
                          paste0("load /home/users/kjyi/Projects/chromophobe/mtoolbox/OUT_", 
                                 Sample,
                                 "/OUT2-sorted.bam\n"), 
                          ""),
              s2 = ifelse(Sample != lag(Sample, default = "a"), 
                          "collapse OUT2-sorted.bam\n", ""),
              s3 = paste0("goto chrRSRS:", getRange(`Variant Allele`), "\n"),
              s4 = paste0("snapshot ", Sample, "_", 
                          str_pad(round(as.numeric(HF) * 100, 0), 3, "left", "0"), "_",
                          str_pad(Locus, 8, "right", "_"), "_",
                          `Variant Allele`, "_", `Aa Change`, ".png\n"))
x <- batchlist %>%
    mutate_all(funs(replace(., is.na(.), ""))) %>%
    as.matrix %>% t %>% c()
h <- "# Running IGV with batch file (screen shot saver)
# https://software.broadinstitute.org/software/igv/batch
new
genome /home/users/kjyi/igv/genomes/chRSRS.genome
expand
snapshotDirectory /home/users/kjyi/Projects/chromophobe/igv/RSRS/coding_mutations
"
cat(h, x, file = "07b_igv.bat", sep = "")
