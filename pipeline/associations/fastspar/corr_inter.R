
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(VennDiagram))

parser <- OptionParser(description = "Corrleation intersection")
parser <- add_option(parser, c("-f", "--file_list"), help="corr pvalue df, sep=, ")
parser <- add_option(parser, c("-n", "--name_list"), help="name, sep=, ")
parse <- parse_args(parser)
df_li <- parse$file_list
n_li <- parse$name_list

df_n_li <- str_split(df_li, ",") %>% unlist()
n_li <- str_split(n_li, ",") %>% unlist()
#df_n_li <- paste0("design_200_sample/new_genome2/test_corr/",
#                  c("golden_corr_pvalue.csv", "gtdb214_L3K11_thr50_corr_pvalue.csv", "metaphlan4_corr_pvalue.csv"))
#n_li <- c("golden", "kssd", "metaphlan4")

corr_pvalue_li <- lapply(df_n_li, function(x) {read_csv(x, col_names = T)})

link_map_unique <- function(link) {
  sort_link <- sort(link)
  return(paste0(sort_link[1], "\t", sort_link[2]))
}

get_links <- function(corr_pvalue_df) {
  n <- nrow(corr_pvalue_df)
  links_v = vector(length = n)
  for(i in 1:n) {
    links_v[i] <- link_map_unique(c(corr_pvalue_df$from[i], corr_pvalue_df$to[i]))
  }
  return(links_v)
}

links_li <- lapply(corr_pvalue_li, function(x) {get_links(x)})

for(i in seq_along(links_li)) {
  write.table(links_li[[i]], paste0(n_li[i], "_links.txt"), row.names = F, quote = F, col.names = F)
}

#get_same_links <- function(links_v1, links_v2) {
#  v_index1 <- vector()
#  v_index2 <- vector()
#  for(i in seq_along(links_v1)) {
#    for(j in seq_along(links_v2)) {
#      if(links_v1[i]==links_v2[j]) {
#        v_index1 <- c(v_index1, i)
#        v_index2 <- c(v_index2, j)
#      }
#    }
#  }
#  return(list(v_index1, v_index2))
#}
#
#all_comb <- combn(1:length(n_li), 2)
#same_links <- vector("list", length = ncol(all_comb))
#for(i in 1:ncol(all_comb)) {
#  comb1 <- all_comb[,i][1]
#  comb2 <- all_comb[,i][2]
#  same_links[[i]] <- get_same_links(links_li[[comb1]], links_li[[comb2]])
#}
#
#comb_links_tbl <- function(tbl1, tbl2, index1, index2, n1, n2) {
#  sub_tbl1 <- tbl1[index1, ]
#  names(sub_tbl1) <- paste(n1, names(tbl1), sep = "_")
#  sub_tbl2 <- tbl2[index2, ]
#  names(sub_tbl2) <- paste(n2, names(tbl2), sep = "_")
#  comb_tbl <- cbind(sub_tbl1, sub_tbl2)
#  return(comb_tbl)
#}
#
#for(i in seq_along(same_links)) {
#  index1 <- same_links[[i]][[1]]
#  index2 <- same_links[[i]][[2]]
#  comb1 <- all_comb[,i][1]
#  comb2 <- all_comb[,i][2]
#  df <- comb_links_tbl(corr_pvalue_li[[comb1]], corr_pvalue_li[[comb2]],
#                        index1, index2, n_li[comb1], n_li[comb2])
#  out_f <- paste(n_li[comb1], "inter", n_li[comb2], "corr", "pvalue", sep = "_")
#  write.csv(df, file = paste0(out_f, ".csv"), row.names = F, quote = F)
#}


#用法举例
#Rscript corr_inter.R -f PRJEB51353_corr_pvalue,PRJNA63661_corr_pvalue.csv,PRJNA799247_corr_pvalue.csv,PRJNA695570_corr_pvalue.csv -n PRJEB51353,PRJNA63661,PRJNA799247,PRJNA695570

