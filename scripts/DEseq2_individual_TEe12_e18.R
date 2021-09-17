#!/usr/bin/env Rscript
library(DESeq2)
library(optparse)
library(tidyverse)
library(readr)
library(stringr)
library(tidyr)

# arguments to provide
option_list = list(
  make_option(c("-c", "--counts"), type="character", default="results/counts.txt", help="counts tabulated file from Feature Counts", metavar="character"),
  make_option(c("-s", "--samplefile"), type="character", default="data/samples2.tsv", help="sample files used to get conditions for DESEq2 model fit", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default="results/deseqNew.csv", help="where to place differential expression files", metavar="character"),
  make_option(c("-f", "--helperfile"), type="character", default="results/helperfile.csv", help="helper file needed for the creation of the clustering and the heatmaps", metavar="character"),
  make_option(c("-m", "--maxfraction"), type="double", default=1.0, help="maximum fraction of the total number of genes to be allowed to be differential between two conditions to be included (number between 0 and 1)", metavar="double")
) 

# parse the command-line arguments and pass them to a list called 'opt'
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# Import the htseq count table

count_table <- "test2.txt"
count <- read_delim(count_table,
  "\t",
  escape_double = FALSE,
  col_names = FALSE,
  trim_ws = TRUE)

colnames(count) <- c("name", "chr", "start", "end", "strand", "width", "e12a", "e12b", "e12c", "e18a", "e18b")
fileExperiment = "samples.tsv"
designExp= read_delim(fileExperiment, "\t", escape_double = FALSE, trim_ws = TRUE)
row.names(designExp) <- designExp$sample
countTable=count[,7:ncol(count)]
row.names(countTable) <- count$name
designExp=designExp[names(countTable),]
designExp$condition=factor(designExp$condition)
dds=DESeqDataSetFromMatrix(countData=countTable,colData=designExp,design = ~condition)
dds$condition=relevel(dds$condition, ref='e12')
dds=DESeq(dds)
res2=results(dds,contrast=c("condition","e18","e12"))
rld=rlog(dds,blind=TRUE)

gtf <- "/exports/humgen/jihed/TE_Transcript/temp/TE_repeat_masker.gtf"
test <- read_delim(gtf,
  "\t",
  escape_double = FALSE,
  col_names = FALSE,
  trim_ws = TRUE)

colnames(test) <- c("chr", "source", "type", "start", "end", "width", "stand", "nothing", "long")
a <- separate(
    data=test,
    col = long,
    sep = "\ ",
    into = c("gene_id", "TE_name", "transcript_id", "TE_dup", "family_id", "Family", "class_id", "Class")
)

b <- str_sub(a$TE_dup, 2,-3)
c<-noquote(b)
d <- cbind (test, c)

data <- "DESeq2_TE_locus_Thymus_WT_e18vse12.txt"

df <- read_delim(data,
  "\t",
  escape_double = FALSE,
  col_names = TRUE,
  trim_ws = TRUE)

colnames(df) <- c("name_TE", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")
colnames(d) <- c("chr", "source", "type", "start", "end", "width", "stand", "nothing", "long", "name_TE")

# Export the tables

df_res2 <- as.data.frame(res2)
write.table(df_res2, file = "DESeq2_TE_locus_Thymus_WT_e18vse12.txt", sep="\t", quote = FALSE, col.names=NA)
df_res3 <- df_res2 %>%
  drop_na()
write.table(df_res3, file = "DESeq2_TE_locus_Thymus_WT_e18vse12_woNA.txt", sep="\t", quote = FALSE, col.names=NA)

df_res4 <- as.data.frame(merge)
df_res5 <- df_res4 %>%
    drop_na()

write.table(df_res5, file = "DESeq2_TE_locus_Thymus_WT_e18vse12_woNA_coordinate.txt", sep="\t", quote = FALSE, col.names=NA)
