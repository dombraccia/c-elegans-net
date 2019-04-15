# coexpression net generation script

# LOADING DATA #
library(tidyverse)
miRNA_bulk <- read.delim("../data/spatiotemporal-data/index.html?acc=GSE115884&format=file&file=GSE115884_miRNA_bulk.tsv.gz")
miRNA_slices <- read.delim("../data/spatiotemporal-data/index.html?acc=GSE115884&format=file&file=GSE115884_miRNA_slices.tsv.gz")
mRNA_bulk <- read.delim("../data/spatiotemporal-data/index.html?acc=GSE115884&format=file&file=GSE115884_mRNA_bulk.tsv.gz")
mRNA_slices <- read.delim("../data/spatiotemporal-data/index.html?acc=GSE115884&format=file&file=GSE115884_mRNA_slices.tsv.gz")
################

# RAW TABLE -> TRANSCRIPT LEVEL EXPRESSION TABLE #
# extracting relevant columns
mRNA_slices_N2 <- filter(mRNA_slices, genotype == "N2")
mRNA_slices_N2 <- mRNA_slices_N2[!is.na(mRNA_slices_N2$slice_index),] # removing rows with NA's for slice index
sample_names <- c("N2_mRNA_A1", "N2_mRNA_A2", "N2_mRNA_A3", 
                  "N2_mRNA_P1", "N2_mRNA_P2", "N2_mRNA_P3") #manually setting sample names
#mRNA_slices_N2$sample_name_num <- as.numeric((mRNA_slices_N2$sample_name))
n_slices <- 13 # only using the slices with all viable samples
n_samples <- length(sample_names)
tmp <- list() # initializing list
for (i in 1:n_slices) {
  for (j in 1:n_samples) {
    k <- 6 * (i - 1) + j
    tmp[[k]] <- filter(mRNA_slices_N2, slice_index == i, sample_name == sample_names[j])
  }
}

# naive approach before using tximport #
#gene_ids <- as.character(unique(mRNA_slices_N2$gene_id))
gene_exp <- tmp[[1]]
gene_exp <- gene_exp %>% group_by(gene_id) %>%
  summarise(gene_cpm = mean(cpm))
colnames(gene_exp)[2]<-paste("slice", "_", 1,"_", "sample", "_", 1, sep = "")

for (k in 2:78) {
  gene_exp <- cbind(gene_exp, (tmp[[k]] %>% group_by(gene_id) %>% summarise(gene_cpm = mean(cpm)))[, 2])
  colnames(gene_exp)[k+1]<- paste("slice", "_", round((ceiling(k/6))), "_", "sample", "_", round((6*(1+(k/6)-ceiling(k/6)))), sep = "")
}

# removing ERCC spike-in counts
gene_exp <- gene_exp[-(1:86), ]

# making gene_id column into rownames and deleting first column
rownames(gene_exp) <- gene_exp$gene_id
gene_exp <- gene_exp[, -1]

# removing any non-zero rows from matrix
gene_exp_nonzero <- gene_exp %>% 
  mutate(row_sum_nonzero = rowSums(gene_exp) != 0) %>% 
  filter(row_sum_nonzero == TRUE) %>%
  select(-row_sum_nonzero)
##################################################

# PREPPING SLICES LIST BEFORE NETWORK GENERATION #
slices <- list()
for (i in 1:13) {
  slices[[i]] <- gene_exp_nonzero[, colnames(gene_exp_nonzero)[(6*(i-1)+1):(6*i)]]
}
##################################################

# GENERATING amat, zmat and pmat AND GRAPH GENERATION #

# making matricies of interest
amat <- matrix(0, ncol = nrow(slices[[1]]), nrow = nrow(slices[[1]]))
zmat <- matrix(0, ncol = nrow(slices[[1]]), nrow = nrow(slices[[1]]))
pmat <- matrix(0, ncol = nrow(slices[[1]]), nrow = nrow(slices[[1]]))
n <- ncol(slices[[1]])
for (i in 1:nrow(slices[[1]])) {
  for (j in i:nrow(slices[[1]])) {
    #calculate current pearson & store
    amat[i, j] <- cor(as.numeric(slices[[1]][i,]), as.numeric(slices[[1]][j,]))
    amat[j, i] <- amat[i, j]
    zmat[i, j] <- atanh(amat[i, j])
    zmat[j, i] <- zmat[i, j]
    pmat[i, j] <- 2 * pnorm(abs(zmat[i, j]), 0, sqrt(1/(n-3)), lower.tail=FALSE)
    pmat[j, i] <- pmat[i, j]
  }
}

#calculating adjusted p-values and making graphs based on p-vals
n <- length(pmat)
pcorrmat <- matrix(0, dim(amat)[1], dim(amat)[2])
for(i in 1:nrow(amat)){ 
  for(j in 1:nrow(amat)){
    rowi <- amat[i, -c(i, j)]
    rowj <- amat[j, -c(i, j)] 
    tmp <- (amat[i, j] - rowi * rowj) / sqrt((1 - rowi^2) * (1 - rowj^2)) 
    tmp.zvals <- (0.5) * log((1 + tmp) / (1 - tmp)) 
    tmp.s.zvals <- sqrt(n - 1 - 3) * tmp.zvals 
    tmp.pvals <- 2 * pnorm(abs(tmp.s.zvals), 0, 1, lower.tail = FALSE) 
    pcorrmat[i, j] <- max(tmp.pvals)
  } 
}
pcorradjmat_bon <- matrix(p.adjust(pcorrmat, method = "bonferroni", n = length(pcorrmat)), nrow = nrow(slices[[1]]), ncol = nrow(slices[[1]]))
pcorradjmat_BH <- matrix(p.adjust(pcorrmat, method = "BH", n = length(pcorrmat)), nrow = nrow(slices[[1]]), ncol = nrow(slices[[1]]))

pcorr_graph_BH <- matrix(as.numeric(pcorradjmat_BH < 0.001), nrow = nrow(slices[[1]]), ncol = nrow(slices[[1]]))
pcorr_graph_bon <- matrix(as.numeric(pcorradjmat_bon < 0.001), nrow = nrow(slices[[1]]), ncol = nrow(slices[[1]]))
sum(pcorr_graph_BH)/sum(pcorr_graph_bon)
sum(pcorr_graph_BH)/n
##########################################################