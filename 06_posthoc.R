suppressMessages(library(tidyverse))
library(stringr)
library(doMC); registerDoMC(20)
load("varscan/all.Rdata")

vcf <- read_tsv("mtoolbox/VCF_file.vcf", comment = "##")
summary <- read_tsv("mtoolbox/summary_20170804_192716.txt", skip = 12, n_max = 17)
# MT dna depth
summary %>% arrange(Sample) %>% .$`Per base depth` %>% barplot(horiz = T)

snp
indel
snp %>% group_by(case) %>% summarise(n = n()) %>% summary(n)
snp %>% filter(VAF > .5) %>% group_by(case) %>% summarise(n = n()) %>% summary(n)
indel %>% group_by(case) %>% summarise(n = n()) %>% summary(n)

# filter
readxl::excel_sheets("database/mitomap.xlsx")
readxl::excel_sheets("database/mtdb.xlsx")
mitomap1 <- readxl::read_xlsx("database/mitomap.xlsx", 1)
mitomap2 <- readxl::read_xlsx("database/mitomap.xlsx", 2)
mitomap <- bind_rows(mitomap1, mitomap2)
mitomap0 <- with(mitomap, paste0(Position, Nucleotide_Change))
snp$in_mitomap <- (with(snp, paste0(POS, REF, "-", ALT)) %in% mitomap0)
mtdb <- readxl::read_xlsx("database/mtdb.xlsx")
mtdb$A[is.na(mtdb$A)] <- 0
mtdb$G[is.na(mtdb$G)] <- 0
mtdb$T[is.na(mtdb$T)] <- 0
mtdb$C[is.na(mtdb$C)] <- 0
mtdb[ , 3:6] %>% apply(1, function(x) sort(x)[3])
mtdb[ , 3:6] %>% apply(1, function(x) sort(x)[3]) %>% density() %>% plot(xlim = c(0,200))
abline(h = 0)
abline(v = c(0,10,20,50,100), col = rainbow(5)) #20
mtdb0 <- c(paste0(mtdb$Posn.[mtdb$A > 20], "A"),
           paste0(mtdb$Posn.[mtdb$C > 20], "C"),
           paste0(mtdb$Posn.[mtdb$T > 20], "T"),
           paste0(mtdb$Posn.[mtdb$G > 20], "G"))
snp$in_mtdb <- with(snp, paste0(POS, ALT)) %in% mtdb0
haplogroups <- read_tsv("database/haplogroups.txt")
snp$in_pylotree16 <- with(snp, paste0(POS, ALT)) %in% haplogroups$POSITIONnucleotidic_change

table(snp$VAF > .50)

snp %>% filter(VAF > .50) %>%
    transmute(not_in_mitomap = !in_mitomap, not_in_mtdb = !in_mtdb, not_in_phylotree = !in_pylotree16) %>%
    colSums()
snp %>% filter(!in_mitomap, !in_mtdb, !in_pylotree16, VAF > .50)
table(snp$in_mitomap)
table(snp$in_mtdb)
table(snp$in_pylotree16)
table(snp$in_mitomap|snp$in_mtdb|snp$in_pylotree16)



snp %>% filter(!in_mitomap, !in_mtdb, !in_pylotree16, VAF > .50) %>%
    mutate(Category = ifelse(category == "LocInfo"|category == "ncUTR", "nc", category)) %>%
    transmute(case, POS, REF, ALT, HF = round(VAF * 100, 1), gene, category)
snp %>% filter(!in_mtdb, !in_pylotree16, VAF > .50) %>%
    transmute(case, POS, REF, ALT, HF = round(VAF * 100, 1), gene, category)
snp %>% filter(!in_mtdb, !in_pylotree16, category == "nsSNP", VAF > .50) %>%
    transmute(case, POS, REF, ALT, HF = round(VAF * 100, 1), gene, category) %>%
    group_by(gene) %>%
    summarise(n = n())
snp %>% filter(!in_mtdb, !in_pylotree16, category == "nsSNP", VAF > .50) %>%
    transmute(case, POS, REF, ALT, HF = round(VAF * 100, 1), gene, category) %>%
    group_by(case, gene) %>%
    summarise(n = n()) %>%
    ungroup() %>% filter(gene != "MT-CYB")
snp %>% filter(!in_mtdb, !in_pylotree16, category == "nsSNP", VAF > .50, POS != 14766) %>%
    transmute(case, POS, REF, ALT, HF = round(VAF * 100, 1), gene, aa) %>%
    print(n = 22)



#mtoolbox
csv <- system("ls mtoolbox/*/*annotation.csv", intern = TRUE)
annotated <- foreach(c = csv, .combine = bind_rows) %do% {
    read_tsv(c)
}
annotated
summary <- read_tsv("mtoolbox/summary_20170804_192716.txt", skip = 12, n_max = 17)
summary %>% 
    transmute(Sample, `Predicted haplogroups` = `Best predicted haplogroup(s)`, `N. of variants` = `N. of homoplasmic variants` +
                  `N. of heteroplasmic variants HF>=0.5`) %>%
    arrange(Sample)
annotated %>%
    filter(is.na(Haplogroup),
           is.na(`Other Haplogroups`),
           RSRS == "yes",
           MHCS == "yes",
           rCRS == "yes",
           (HF > .50 | is.na(HF)),
           `Aa Change` != "syn") %>%
    select(Sample:HF, Locus, `Aa Change`)


annotated %>% filter(`Aa Change` == "Stop-gain", HF > 0.2)
select(Sample:HF, Locus, `Aa Change`)

vcf <- read_tsv("mtoolbox/VCF_file.vcf", comment = "##")
vcf %>% select(POS, REF, ALT, chRCC04, chRCC15) %>%
    filter(between(POS, 11860, 11870) | between(POS, 12380, 12390))

`%::%` <- function(x, n) {
    s <- str_split(x, ":")
    idx <- 1:5
    names(idx) <- c("GT", "DP", "HF", "CILOW", "CIUP")
    sapply(s,function(x) ifelse(is.na(x[idx[n]]), 0, x[idx[n]]))}
`%,%` <- function(x, n) {
    sapply(str_split(x, ","),
           function(x) ifelse(is.na(x[n]), 0, x[n]))}
`%+%` <- function(x, o) sapply(str_split(x, ","), function(y) sum(as.numeric(y)))
`%m%` <- function(x, o) sapply(str_split(x, ","), function(y) min(as.numeric(y)))
`%M%` <- function(x, o) sapply(str_split(x, ","), function(y) max(as.numeric(y)))
vcf %>% select(POS, REF, ALT, chRCC04) %>%
    filter(between(POS, 11860, 11872)) %>% .$chRCC04 %::% "HF" %,% 1:2 %>% .[1:2,1:2] %>% as.numeric %>% sum
vcf %>% select(POS, REF, ALT, chRCC15) %>%
    filter(between(POS, 12379, 12390)) %>% .$chRCC15 %::% "HF" %,% 1 %>% .[1:3] %>% as.numeric %>% sum

library(extrafont)
extrafont::font_import()
extrafont::fonttable()
# Functional annotation
if (F) {
    annotated %>%
        filter(is.na(Haplogroup),
               is.na(`Other Haplogroups`),
               RSRS == "yes",
               MHCS == "yes",
               rCRS == "yes",
               (HF > .50),
               `Aa Change` != "syn") %>% 
        select(Sample:HF, Locus, `Aa Change`, `Codon Position`, `Disease Score`,
               `MutPred pred`:`Mitomap Associated Disease(s)`,
               `1000 Genomes Heteroplasmy`) %>%
        as.data.frame() %>%
        column_to_rownames("Sample")
}


#################
# violin plot       --------------
#################

# data prep
dt_tmp <- annotated %>%
    mutate(selected = is.na(Haplogroup) &
               is.na(`Other Haplogroups`) &
               RSRS == "yes" &
               MHCS == "yes" &
               rCRS == "yes" &
               (HF > .50) &
               `Aa Change` != "syn") %>% 
    mutate(dot = ".") %>%
    select(Sample, `Variant Allele`, selected, Locus, `Aa Change`,
           `Disease Score`, dot, `MutPred pred`:`SNPs&GO prob`) %>% View
colnames(dt_tmp) %>% as.matrix(ncol = 1)
colnames(dt_tmp) <- NULL

dt <- rbind(cbind(dt_tmp[ , 1:5], pridictor = "Disease Score", dt_tmp[ , c(7, 6)]),
                            cbind(dt_tmp[ , 1:5], pridictor = "MutPhred", dt_tmp[ , c(8, 9)]),
            cbind(dt_tmp[ , 1:5], pridictor = "PolyPhen-2 HumDiv", dt_tmp[ , c(10, 11)]),
            cbind(dt_tmp[ , 1:5], pridictor = "PolyPhen-2 HumVar", dt_tmp[ , c(12, 13)]),
            cbind(dt_tmp[ , 1:5], pridictor = "PANTHER", dt_tmp[ , c(14, 15)]),
            cbind(dt_tmp[ , 1:5], pridictor = "PhD-SNP", dt_tmp[ , c(16, 17)]),
            cbind(dt_tmp[ , 1:5], pridictor = "SNPs&GO", dt_tmp[ , c(18, 19)]))
colnames(dt) <- c("Sample", "Variant Allele", "selected", "Locus", "Aa Change", "predictor", "decision", "prob")
suppressWarnings(dt$prob <- as.numeric(dt$prob))
df <- dt[!is.na(dt$prob), ]
rm(dt_tmp)

# http://shinyapps.org/apps/RGraphCompendium/index.php
source('../codes/violinplot.R')

mytheme <- theme_bw(base_size = 16, base_family = "DejaVu Sans") +
    theme(axis.text.x     = element_text(size = 14),
          axis.title.y    = element_text(vjust = +1.5),
          panel.grid.major  = element_blank(),
          panel.grid.minor  = element_blank(),
          legend.background = element_blank(),
          legend.key = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.line  = element_line(colour = "black"))
# creat variable for plotting
df$predictor <- as.factor(df$predictor)
df$grpN <- as.numeric(df$predictor)
df$sv <- paste(df$Locus, df$`Aa Change`, paste0("(", df$Sample, 
                                                ")"))

probs <- c(0, 0.25, 0.5, 0.75, 1)
# get quantiles
qtiles <- ldply(unique(df$grpN),
                function(gr) quantile(round(df$prob[df$grpN == gr],digits = 4), 
                                      probs, na.rm = T, type = 3))
freqs  <- ldply(unique(df$grpN),
                function(gr) table(cut(df$prob[df$grpN == gr],
                                       breaks = qtiles[gr, ],
                                       na.rm = T,
                                       include.lowest = T,
                                       right = T)))
labels <- sapply(unique(df$grpN),
                 function(gr)levels(cut(round(df$prob[df$grpN == gr],
                                              digits = 4),
                                        breaks = qtiles[gr, ],
                                        na.rm = T,
                                        include.lowest = T,
                                        right = T)))

# Get regular violinplot using package ggplot2
g.pv <- ggplot(df,aes(x = predictor,y = prob)) + 
    geom_violin(aes(group = predictor),
                scale = "width", color = "grey30", fill = "grey30",
                trim = T, adjust = .7)

# Cut at quantiles using vioQtile() in C-3PR
g.pv0 <- vioQtile(g.pv, qtiles, probs)

fivecol <- c("#45B3FF", "#49E83E", "#FFD432", "#E84B30", "#B243FF")

g.pv1 <- g.pv0 + 
    geom_point(data = subset(df, selected == TRUE),
               color = "black", size = 5, shape = 16) +
    geom_point(data = subset(df, selected == TRUE),
               aes(color = sv), size = 4.5, shape = 16) +
    scale_color_manual(name = "Variants", values = fivecol)
tt <- "Distribution of functional annotation scores of \nunfiltered & filtered variants"
g.pv2 <- g.pv1 +
    xlab("") + ylab("") + 
    ggtitle(tt) +
    mytheme + theme(axis.text.x = element_text(angle = -45, hjust = 0),
                    title = element_text(size = 13),
                    plot.title = element_text(hjust = 0.5),
                    legend.title = element_text(size = 15))
g.pv2
svg("figure/unfiltered_violine.svg", 980/72, 549/72 )
g.pv2
dev.off()
##

# other cancer somatic ns mutations
readxl::excel_sheets("database/elife-02935-supp2-v2.xlsx")
readxl::excel_sheets("database/elife_case.xlsx")
readxl::excel_sheets("database/elife-02935-supp3-v2.xlsx")
patho_table <- read_tsv("database/patho_table.txt", col_types = cols(.default = "c"))
patho_table[1000:1003, ] %>% as.data.frame
cg_somatic <- readxl::read_xlsx("database/elife-02935-supp2-v2.xlsx", 1) %>%
    select(sample_Index:`Confirmed_disease_association*`)
cg_germline <- readxl::read_xlsx("database/elife-02935-supp2-v2.xlsx", 2)
case_info <- readxl::read_xlsx("database/elife_case.xlsx")
types <- case_info$TumourTypes
names(types) <- case_info$index
cg <- full_join(cg_somatic, cg_germline)
cg$type <- types[cg$sample_Index]
cg %>% filter(`VAF(T)` > .50) %>%
    filter(!wt_allele %in% c("T", "C", "G", "A") | !var_allele %in% c("T", "C", "G", "A")) %>%
    nrow()
# 59/2158
cg_snp <- cg %>% filter(`VAF(T)` > .50) %>%
    filter(wt_allele %in% c("T", "C", "G", "A") | var_allele %in% c("T", "C", "G", "A")) %>%
    mutate(variant = paste0(mtDNA_pos, var_allele))
cg_snp$variant
cg_ready <- cg_snp %>% select(variant, type)

patho_table_ready <- patho_table %>% 
    select(Variant, `Disease Score`, `MutPred Prob`, `PolyPhen-2 HumDiv Prob`, 
           `PolyPhen-2 HumVar Prob`, `PANTHER Prob`, `PhD-SNP Prob`,
           `SNPs&GO Prob`) %>% na.omit %>%
    gather(-Variant, key = "predictor", value = "prob")


cg_scored <- foreach(i = 1:nrow(cg_ready), .combine = bind_rows) %do% {
    t <- filter(patho_table_ready, Variant == cg_ready[[i, 1]])
    if (nrow(t) != 0) mutate(t, type = cg_ready[[i, 2]]) else NULL
}
suppressWarnings(cg_scored$prob <- as.numeric(cg_scored$prob))
cg_scored <- cg_scored[!is.na(cg_scored$prob), ]
df2 <- cg_scored
# is.na(df2) %>% colSums # 893/2269 na
# violin plot
arr <- unique(df$predictor)
names(arr) <- unique(df2$predictor)[c(1,2,3,4,7,5,6)]
# arr
df2$predictor <- arr[df2$predictor]
df2$predictor <- factor(df2$predictor, levels = unique(df$predictor))
df2$grpN <- as.numeric(df2$predictor)
# df2$sv <- paste(df2$Locus, df2$`Aa Change`, paste0("(", df2$Sample,
#                                                 ")"))

# get quantiles
source('../codes/violinplot.R')

# Get regular violinplot using package ggplot2
g2.pv <- ggplot(df2, aes(x = predictor,y = prob)) + 
    geom_violin(aes(group = predictor),
                scale = "width", color = "grey30", fill = "grey30",
                trim = T, adjust = .7)
# Cut at quantiles using vioQtile() in C-3PR
g2.pv0 <- vioQtile(g2.pv)

fivecol <- c("#45B3FF", "#49E83E", "#FFD432", "#E84B30", "#B243FF")

g2.pv1 <- g2.pv0 + 
    geom_point(data = subset(df, selected == TRUE),
               color = "black", size = 5, shape = 16) +
    geom_point(data = subset(df, selected == TRUE),
               aes(color = sv), size = 4.5, shape = 16) +
    scale_color_manual(name = "Variants", values = fivecol)
tt2 <- "Distribution of functional annotation scores of \nVariants in other study"
g2.pv2 <- g2.pv1 +
    xlab("") + ylab("") + 
    ggtitle(tt2) +
    mytheme + theme(axis.text.x = element_text(angle = -45, hjust = 0),
                    title = element_text(size = 13),
                    plot.title = element_text(hjust = 0.5),
                    legend.title = element_text(size = 15))
g2.pv2
svg("figure/cgp_violine.svg", 980/72, 549/72 )
g2.pv2
dev.off()


##
`%;%` <- function(x, n) {
    s <- str_split(x, ";")
    idx <- 1:5
    names(idx) <- c("GT", "DP", "HF", "CILOW", "CIUP")
    sapply(s,function(x) ifelse(is.na(x[idx[n]]), 0, x[idx[n]]))} 
annotated$tmp <- annotated$`CI_lower;CI_upper` %;% 2
annotated %>% filter(tmp != "1.0") %>% View
    
    
## mitochondrial read count
# samtools view -c
# library(tidyverse)
# rm(`%:%`)
# library(doMC)
# registerDoMC(20)
# 
# chr <- c(1:22, "X", "Y", "MT")
# sample <- dir(".bamfiles", "bam$", full.names = TRUE)
# if(F){ # takes long time
#     read_count_by_chromosome <- 
#         foreach(l = chr, .combine = cbind) %:%
#         foreach(s = sample, .combine = rbind) %dopar% {
#             cmd <- paste("samtools view -c", s, l)
#             x <- system(cmd, intern = TRUE)
#             return(as.integer(x))
#             # return(cmd)
#         } %>% as.data.frame()
# }
# colnames(read_count_by_chromosome) <- chr
# rownames(read_count_by_chromosome) <- stringr::str_extract(sample, "chRCC\\d+")
# 
# 
# 
# read_count_by_chromosome$autosome <- read_count_by_chromosome %>%
#     select(-X, -Y, -MT) %>%
#     rowSums()
# 
# read_count_by_chromosome <- read_count_by_chromosome %>% 
#     mutate(mt_copy = read_count_by_chromosome$MT/read_count_by_chromosome$autosome/16969*2860000000*2)
# read_count_by_chromosome %>% write_tsv("read_count_by_chromosome.tsv")
read_count_by_chromosome <- read_tsv("read_count_by_chromosome.tsv")
read_count_by_chromosome %>% as.data.frame() %>%
    select(case, mt_copy) %>% column_to_rownames("case") %>% as.matrix() %>% t() %>%
    barplot(horiz = T, las = 2)
ann = c(".", ".", "yello", "frameshift", ".", ".", "green", "frameshift", "red", ".", "purple", ".", ".", "frameshift", ".", ".", "blue")
read_count_by_chromosome %>% arrange(case) %>% transmute(case, mt_copy, ann = factor(ann, levels = c(".", "yello", "green", "red", "purple", "blue", "frameshift"))) -> read_count_by_chromosome2
colours <- c("grey60","#45B3FF", "#49E83E", "#FFD432", "#E84B30", "#B243FF", "grey20")
mc <- read_count_by_chromosome2 %>% 
    ggplot(aes(x = reorder(stringr::str_extract(case, "\\d+"), mt_copy), y = mt_copy, fill = ann)) + geom_bar(stat = "identity") +
    xlab("Sample") + ylab("Mitochondrial copy number") +
    scale_fill_manual(name = "Coding Variant", 
        labels = c("", "MT-ND5 N207S", "MT-CYB E271K", "MT-NDS V243I", "MT-ND6 A142V", "MT-ATP6 M140T", "Framshift"), 
                      values = colours )
mc
svg("figure/mitochondrial_copy.svg", 980/72, 549/72 )
mc
dev.off()

### filter cg also
rm(list = ls())
library(tidyverse)
library(doMC);registerDoMC(23)

readxl::excel_sheets("database/elife-02935-supp2-v2.xlsx")

# cg_somatic
cg_somatic <- readxl::read_xlsx("database/elife-02935-supp2-v2.xlsx", 1) %>%
    select(sample_Index:`Confirmed_disease_association*`)
case_info <- readxl::read_xlsx("database/elife_case.xlsx")
types <- case_info$TumourTypes
names(types) <- case_info$index
cg_somatic$type <- types[cg_somatic$sample_Index]
rm(case_info, types)
cg_somatic <- cg_somatic %>% filter(`VAF(T)` > .50) %>%
    filter(wt_allele %in% c("T", "C", "G", "A") | var_allele %in% c("T", "C", "G", "A")) %>%
    mutate(variant = paste0(mtDNA_pos, var_allele))

# filters -> cg_ready
mtdb <- readxl::read_xlsx("database/mtdb.xlsx")
mtdb$A[is.na(mtdb$A)] <- 0
mtdb$G[is.na(mtdb$G)] <- 0
mtdb$T[is.na(mtdb$T)] <- 0
mtdb$C[is.na(mtdb$C)] <- 0
mtdb0 <- c(paste0(mtdb$Posn.[mtdb$A > 20], "A"),
           paste0(mtdb$Posn.[mtdb$C > 20], "C"),
           paste0(mtdb$Posn.[mtdb$T > 20], "T"),
           paste0(mtdb$Posn.[mtdb$G > 20], "G"))
haplogroups <- read_tsv("database/haplogroups.txt")

cg_somatic$in_mtdb <- with(cg_somatic, paste0(mtDNA_pos, var_allele)) %in% mtdb0
cg_somatic$in_pylotree16 <- with(cg_somatic, paste0(mtDNA_pos, var_allele)) %in% haplogroups$POSITIONnucleotidic_change

cg_ready <- cg_somatic %>% filter(`VAF(T)` > .5, !in_mtdb, !in_pylotree16) %>% select(variant, type, sample_Index)
rm(cg_somatic, haplogroups, mtdb, mtdb0)

# prepare annotation table : patho_table_ready
patho_table <- read_tsv("database/patho_table.txt", col_types = cols(.default = "c"))

patho_table_ready <- patho_table %>% 
    select(Variant, `Disease Score`, `MutPred Prob`, `PolyPhen-2 HumDiv Prob`, 
           `PolyPhen-2 HumVar Prob`, `PANTHER Prob`, `PhD-SNP Prob`,
           `SNPs&GO Prob`) %>% na.omit %>%
    gather(-Variant, key = "predictor", value = "prob")
rm(patho_table)
#annotation : df3
cg_scored <- foreach(i = 1:nrow(cg_ready), .combine = bind_rows) %do% {
    t <- filter(patho_table_ready, Variant == cg_ready[[i, 1]])
    if (nrow(t) != 0) mutate(t, type = cg_ready[[i, 2]], sample_Index = cg_ready[[i, 3]]) else NULL
}
suppressWarnings(cg_scored$prob <- as.numeric(cg_scored$prob))
cg_scored <- cg_scored[!is.na(cg_scored$prob), ]
df3 <- cg_scored
df3$predictor <- stringr::str_replace(df3$predictor, " Prob", "")
rm(cg_ready, cg_scored, patho_table_ready, t, i)
df3 %>% group_by(predictor) %>% dplyr::summarize(n = n())
df3$sample_Index %>% unique %>% length


# filtered 5 mutations
csv <- system("ls mtoolbox/*/*annotation.csv", intern = TRUE)
annotated <- foreach(c = csv, .combine = bind_rows) %do% {
    read_tsv(c)
}; rm(csv, c)
dt_selected <- annotated %>%
    filter(is.na(Haplogroup), is.na(`Other Haplogroups`), RSRS == "yes",
           MHCS == "yes", rCRS == "yes", (HF > .50), `Aa Change` != "syn") %>% 
    mutate(dot = ".") %>%
    mutate(dot2 = ".") %>%
    select(Sample, `Variant Allele`, dot, Locus, `Aa Change`,
           `Disease Score`, dot2, `MutPred pred`:`SNPs&GO prob`)
colnames(dt_selected) <- NULL
dt_selected.t <- rbind(cbind(dt_selected[ , 1:5], pridictor = "Disease Score", dt_selected[ , c(7, 6)]),
            cbind(dt_selected[ , 1:5], pridictor = "MutPred", dt_selected[ , c(8, 9)]),
            cbind(dt_selected[ , 1:5], pridictor = "PolyPhen-2 HumDiv", dt_selected[ , c(10, 11)]),
            cbind(dt_selected[ , 1:5], pridictor = "PolyPhen-2 HumVar", dt_selected[ , c(12, 13)]),
            cbind(dt_selected[ , 1:5], pridictor = "PANTHER", dt_selected[ , c(14, 15)]),
            cbind(dt_selected[ , 1:5], pridictor = "PhD-SNP", dt_selected[ , c(16, 17)]),
            cbind(dt_selected[ , 1:5], pridictor = "SNPs&GO", dt_selected[ , c(18, 19)]))
colnames(dt_selected.t) <- c("Sample", "Variant Allele", "selected", "Locus", "Aa Change", "predictor", "decision", "prob")
suppressWarnings(dt_selected.t$prob <- as.numeric(dt_selected.t$prob))
dt_selected.to <- dt_selected.t[!is.na(dt_selected.t$prob), ]
rm(dt_selected, dt_selected.t, annotated)
dt_selected.to$sv <- paste(dt_selected.to$Locus, dt_selected.to$`Aa Change`, paste0("(", dt_selected.to$Sample, ")"))
ds <- dt_selected.to %>% select(`Variant Allele`, predictor, prob, sv)
levels(ds$predictor) == unique(df3$predictor)
# creat variable for plotting
df3$predictor <- factor(df3$predictor, levels = levels(ds$predictor))
df3$grpN <- as.numeric(df3$predictor)

# get quantiles
source('../codes/violinplot.R')

# Get regular violinplot using package ggplot2
g3.pv <- ggplot(df3, aes(x = predictor, y = prob)) + 
    geom_violin(aes(group = predictor),
                scale = "width", color = "grey30", fill = "grey30",
                trim = T, adjust = .7)
# Cut at quantiles using vioQtile() in C-3PR
g3.pv0 <- vioQtile(g3.pv)

fivecol <- c("#45B3FF", "#49E83E", "#FFD432", "#E84B30", "#B243FF")

g3.pv1 <- g3.pv0 + 
    geom_point(data = dt_selected.to ,
               color = "black", size = 5, shape = 16) +
    geom_point(data = dt_selected.to ,
               aes(color = sv), size = 4.5, shape = 16) +
    scale_color_manual(name = "Variants", values = fivecol)
tt3 <- "Distribution of functional annotation scores of \nVariants in other study"
mytheme <- theme_bw(base_size = 16, base_family = "DejaVu Sans") +
    theme(axis.text.x     = element_text(size = 14),
          axis.title.y    = element_text(vjust = +1.5),
          panel.grid.major  = element_blank(),
          panel.grid.minor  = element_blank(),
          legend.background = element_blank(),
          legend.key = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.line  = element_line(colour = "black"))
g3.pv2 <- g3.pv1 +
    xlab("") + ylab("") + 
    ggtitle(tt3) +
    mytheme + theme(axis.text.x = element_text(angle = -45, hjust = 0),
                    title = element_text(size = 13),
                    plot.title = element_text(hjust = 0.5),
                    legend.title = element_text(size = 15))
g3.pv2
svg("figure/cgp_violine.svg", 980/72, 549/72 )
g3.pv2
dev.off()
