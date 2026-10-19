suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(optparse))

parser <- OptionParser(description = "SparCC corrleation pvalue")
parser <- add_option(parser, c("-c", "--corr"), help="SparCC corrletion file")
parser <- add_option(parser, c("-p", "--pvalue"), help="SparCC pvalue file")
parser <- add_option(parser, c("-s", "--species"), help="id species map")
parser <- add_option(parser, c("--corr_num"), type="double", default=0.1, help="correlation")
parser <- add_option(parser, c("--pvalue_num"), type="double", default=0.0001, help="pvalue")
parser <- add_option(parser, c("-o", "--out"),  help="out file name prefix")

parse <- parse_args(parser)
corr_fname <- parse$corr
species_f <- parse$species
#tool <- parse$tools
pvalue_f <- parse$pvalue
corr_num <- parse$corr_num
pvalue_num <- parse$pvalue_num
out_f <- parse$out

suffix <- corr_fname %>% str_sub(-3, -1)
#
#corr_fname <- "select_110_sample/false_positive/fastsparcc_res/kssd_correlation.tsv"
#pvalue_f <- "select_110_sample/false_positive/fastsparcc_res/kssd_pvalues.tsv"
#species_f <- "select_110_sample/false_positive/comb_profile/kssd_species.txt"
if(suffix == "csv") {
  corr_f <- read.csv(corr_fname, row.names = 1, header = T)
  pvalue_df <- read.csv(pvalue_f, row.names = 1, header = T)
} else {
  corr_f <- read.table(corr_fname, row.names = 1, header = F)
  pvalue_df <- read.table(pvalue_f, row.names = 1, header = F)
}

sp <- read_tsv(species_f, col_names = c("num", "species"))

species <- sp$species

mat_df <- function(mat, col_name) {
  row.names(mat) <- col_name
  colnames(mat) <- col_name
  mat[lower.tri(mat, diag = T)]=NA
  ind <- which(is.na(mat)==F,arr.ind = T)
  df <- data.frame(from = colnames(mat)[ind[,1]],
                   to = colnames(mat)[ind[,2]], corr = mat[ind],stringsAsFactors = F)
  return(df)
}
corr_df <- mat_df(corr_f, species)

corr_df$pvalue <- mat_df(pvalue_df, species)[,3]

result <- corr_df %>% filter(abs(corr) >= corr_num, pvalue <= pvalue_num)
#result <- corr_df
write.csv(result, paste0(out_f, ".csv"), quote = F, row.names = F)


#用法举例
#Rscript sparcc_corr_pvalue.R -c PRJNA823736_correlation.tsv -p PRJNA823736_pvalues.tsv -s PRJNA823736_species.txt -o PRJNA823736_corr_pvalue
