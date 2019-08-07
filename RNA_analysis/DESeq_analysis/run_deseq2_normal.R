library("ggplot2")
library("RColorBrewer")

library('DESeq2')


args = commandArgs(trailingOnly=TRUE)

print(args[1])
print(args[2])


# Load expression data
indir='/cellar/users/mdow/Projects/HCC/RNA_analysis/'
exp_indir = paste(indir,'tcga_all_rawcounts.csv',sep='')
exp_df = read.csv(exp_indir,check.names=FALSE,row.names=1)
print(dim(exp_df))

#exp_df_n = as.data.frame(sapply(exp_df, as.integer))

exp_df_n = as.data.frame(sapply(exp_df, as.integer))
rownames(exp_df_n) <- rownames(exp_df)

print(dim(exp_df_n))

run_deseq <- function(exp_df,cond_df,fname){
    dds <- DESeqDataSetFromMatrix(countData = exp_df, colData = cond_df, design = ~ condition)
    # Eliminate rows with low expression and columns with all 0s
    x = 10
    dds <- dds[rowSums(counts(dds)) > x, colSums(counts(dds)) > 0]
    
    cts <- counts(dds)
    geoMeans <- apply(cts, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
    dds <- estimateSizeFactors(dds, geoMeans=geoMeans)

    dds <- DESeq(dds)

    groups = unique(cond_df$condition)
    test_group= as.character(groups[groups != 'normal'])
    #print(groups)
    print(test_group)
    res <- results(dds,contrast=c('condition',test_group,'normal'))
    write.csv(as.data.frame(res), file=fname)

    return(dds)
}

# Change the condition to be argument
print('start')
indir_cond = '/cellar/users/mdow/Projects/HCC/RNA_analysis/DESeq_analysis/human/'

infname=args[1]
this_infname = paste(indir_cond,infname,sep='')
print(this_infname)
condDf = read.csv(this_infname,check.names=FALSE,row.names=1)
print(dim(condDf)) 
print(head(condDf)) 
                      
outfname=args[2]

condDf_srt = condDf[order(rownames(condDf)),]
exp_srt = exp_df_n[,order(colnames(exp_df_n))]
                      
dd_out <- run_deseq(exp_srt,condDf_srt,outfname)

print(outfname)