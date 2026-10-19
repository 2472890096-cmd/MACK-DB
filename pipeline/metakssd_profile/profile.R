#! /usr/bin/R

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(optparse))

parser <- OptionParser(description = "calculate kssd species abundance:")
parser <- add_option(parser, c("-l", "--list"), help="sample list")
parser <- add_option(parser, c("-i", "--input"), default="./", help="sample input dir, default=./")
parser <- add_option(parser, c("-s", "--sharekmer"), type="integer", default=12, help="share kmer, default=12")
parser <- add_option(parser, c("-o", "--out"),  help="sample out dir")
parse <- parse_args(parser)

li <- parse$l
in_dir <- parse$input
sk <- parse$sharekmer
out_dir <- parse$out

poissonN <- c(0,0.02,0.08,0.69,1.57,2.49,3.41,4.31,5.20,6.08,6.94)
poissonN_thr <- 5
avgpct9899_offset <- 3

change_avgpct9899 <- function(kssd_df, thr=poissonN_thr, offset=avgpct9899_offset, possion_v=poissonN) {
    depth <- vector(length=nrow(kssd_df))
    for(i in 1:nrow(kssd_df)) {
        if(kssd_df[i,]$X5 > thr) {
            depth[i] <- kssd_df[i,]$X5 - offset
        } else {
              depth[i] <- possion_v[kssd_df[i,]$X5 + 1]
        }
    }
    kssd_df$depth <- depth
    return(kssd_df)
}

add_ratio <- function(kssd_df) {
  kssd_df$ratio <- round(kssd_df$depth/sum(kssd_df$depth)*100, 4)
  kssd_df <- kssd_df %>% select(X2, ratio)
  names(kssd_df) <- c("gtdb_sp", "ratio")
  return(kssd_df)
}

# Run list
run_li <- read_tsv(li, col_names="run")
# profile
kssd_profile_li <- lapply(paste0(in_dir, "/", run_li$run, "_all.txt"),
                      function(x) {read_tsv(x, col_names = F)})
names(kssd_profile_li) <- run_li$run
# filter by share kmer
kssd_profile_li <- lapply(kssd_profile_li, function(x) {x %>% filter(X3 >= sk)})
# delete run
nrow_li <- lapply(kssd_profile_li, function(x) {nrow(x)}) %>% unlist()
del_index <- which(nrow_li != 0)
#print(del_index)
kssd_profile_li <- kssd_profile_li[del_index]
#for(i in seq_along(kssd_profile_li)) {
#  kssd_profile_li[[i]]$X1 <- kssd_profile_li[[i]]$X1 %>% str_extract("[A-Z].+")
#}

kssd_profile_li <- lapply(kssd_profile_li, function(x) {change_avgpct9899(x)})

kssd_profile_li <- lapply(kssd_profile_li, function(x) {add_ratio(x)})

if(!dir.exists(out_dir)) {
    dir.create(out_dir)
}

for(i in seq_along(kssd_profile_li)) {
    write.table(kssd_profile_li[[i]] %>% select(gtdb_sp, ratio),
               file = paste0(out_dir, "/", names(kssd_profile_li)[i], "_profile.txt"),
               row.names = F, col.names = F, quote = F, sep = "\t")
}
