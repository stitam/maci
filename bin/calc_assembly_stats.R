# NOTE
# Here, we define contig length as contig length without any "n" characters.
# While fasta files are naturally sorted by contig length, subtracting the
# number of "n" characters may impact the order. This is currently not handled
# in this script. The script always returns the "adjusted" length of the first
# contig.

rm(list=ls())
library(dplyr)
args <- commandArgs(trailingOnly = TRUE)

sample_id <- args[1]
f <- seqinr::read.fasta(args[2], forceDNAtolower = FALSE)

contig <- data.frame(
    name = names(f),
    gc = unname(sapply(f, function(x) {seqinr::GC(x)})),
    n_count = unname(sapply(f, function(x) {
        nuc_count <- table(x)
        nonamb <- c("a", "A", "t", "T", "g", "G", "c", "C")
        if (any(names(nuc_count) %in% nonamb == FALSE)) {
            out <- sum(nuc_count[which(names(nuc_count) %in% nonamb == FALSE)])
            return(out)
        } else {
            return(0)
        }
    })),
    length = unname(sapply(f, length))
)
contig <- contig[order(contig$length, decreasing = TRUE),]
contig$cumsum_quantile <- cumsum(contig$length)/sum(contig$length, na.rm = TRUE)

df <- data.frame(
    assembly = sample_id,
    genome_size = sum(contig$length, na.rm = TRUE),
    gc_content = round(sum(contig$gc*contig$length)/sum(contig$length), 3),
    n_count = sum(contig$n_count),
    contig_count = nrow(contig),
    longest_contig = (contig$length - contig$n_count)[1],
    N50 = contig$length[min(which(contig$cumsum_quantile > 0.5))],
    N95 = contig$length[min(which(contig$cumsum_quantile > 0.95))]
)

write.table(
    df,
    file = paste0(sample_id, ".tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
)
