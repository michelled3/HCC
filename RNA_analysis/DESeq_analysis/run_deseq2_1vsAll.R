library("ggplot2")
library("RColorBrewer")

library('DESeq2')

# Load expression data
human_raw = read.csv('tcga_lihc_rawcounts.csv',row.names=1,check.names=FALSE)

# Load grouping data

myClust = read.csv('human_mouseClust.csv',check.names=FALSE)
colnames(myClust) = c('PID','condition')

# Sort tables

human_srt <- human_raw[,order(names(human_raw))]
human_srt_n = as.data.frame(sapply(human_srt, as.integer))
rownames(human_srt_n) <- rownames(human_srt)

myClust_srt <- myClust[order(myClust$PID),]
myClust_srt_n <- as.data.frame(sapply(myClust, as.numeric))
rownames(myClust_srt_n) <- myClust_srt$PID

# Compare paired conditions

#cond12 = myClust_srt[which(myClust_srt$condition == 0|myClust_srt$condition == 1),]
#cond13 = myClust_srt[which(myClust_srt$condition == 0|myClust_srt$condition == 2),]
#cond23 = myClust_srt[which(myClust_srt$condition == 1|myClust_srt$condition == 2),]
#print(c(dim(cond12),dim(cond13),dim(cond23)))

#df12 = human_srt_n[cond12$PID]
#df13 = human_srt_n[cond13$PID]
#df23 = human_srt_n[cond23$PID]
#print(c(dim(df12),dim(df13),dim(df23)))

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
    
    res <- results(dds)
    write.csv(as.data.frame(res), file=fname)

    return(dds)
}

# Create 1 vs all table
test = myClust_srt
test$condition[test$condition != 0] <- 1 # change 1 and 2 to 1
                      
dd_12 <- run_deseq(human_srt_n,test,'group1_vs_2and3.csv')
print('Done 1_vs_2and3')
                      
test = myClust_srt
test$condition[test$condition != 1] <- 0 # change 2 and 0 to 0
                      
dd_13 <- run_deseq(human_srt_n,test,'group2_vs_1and3.csv')
print('Done 2_vs_1and3')
                      
test = myClust_srt
test$condition[test$condition != 2] <- 0 # change 1 and 0 to 0
                      
dd_23 <- run_deseq(human_srt_n,test,'group3_vs_1and2.csv')
print('Done 3_vs_1and2')


