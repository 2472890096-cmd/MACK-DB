#! /usr/bin/R

suppressPackageStartupMessages(library(optparse))

parser <- OptionParser(description = "ratio profile tbl shape into count tbl")
parser <- add_option(parser, c("-i", "--input"), help="input ratio table")
parser <- add_option(parser, c("-n", "--num"),type="integer", default=1000000, help="input ratio table, default 1000000")
parser <- add_option(parser, c("-o", "--out"), help="out count table")

parse <- parse_args(parser)
input_f <- parse$input
num <- parse$num
out_f <- parse$out

ratio_profile <- read.table(input_f, row.names = 1, header = T)

count_profile <- ratio_profile * num
count_profile <- rbind(0:(ncol(count_profile)-1), count_profile)
row.names(count_profile) <- c("OTU_id", 0:(nrow(count_profile)-2))

write.table(count_profile, out_f, row.names = T, col.names = F, quote = F, sep = "\t")
