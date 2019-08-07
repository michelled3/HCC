library("ggplot2")
library("RColorBrewer")

library('DESeq2')

# Load expression data
human_raw = read.csv('/cellar/users/mdow/Projects/HCC/RNA_analysis/tcga_iClust_rawcounts.csv',row.names=1,check.names=FALSE)

# Load grouping data

myClust = read.csv('human/prolif_cond.csv',row.names=1)
colnames(myClust) = c('PID','condition')
# Sort tables

human_srt <- human_raw[,order(names(human_raw))]
human_srt_n = as.data.frame(sapply(human_srt, as.integer))
rownames(human_srt_n) <- rownames(human_srt)

myClust_srt <- myClust[order(myClust$PID),]

run_deseq <- function(exp_df,cond_df,fname){
    print(fname)
    dds <- DESeqDataSetFromMatrix(countData = exp_df, colData = cond_df, design = ~ condition)
    # Eliminate rows with low expression and columns with all 0s
    x = 10
    dds <- dds[rowSums(counts(dds)) > x, colSums(counts(dds)) > 0]
    
    cts <- counts(dds)
    geoMeans <- apply(cts, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
    dds <- estimateSizeFactors(dds, geoMeans=geoMeans)

    dds <- DESeq(dds)
    
    res <- results(dds,contrast=c('condition','high_prolif','low_prolif'))
    write.csv(as.data.frame(res), file=fname)

    return(dds)
}

dd <- run_deseq(human_srt_n,myClust_srt,'prolif_high_vs_low.csv')
print('Done Prolif')


