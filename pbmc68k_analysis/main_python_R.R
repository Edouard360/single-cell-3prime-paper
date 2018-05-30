# Title     : main_python_R.R
# Objective : To save as python .svmlight datasets the necessary resources
# Created by: edouardm
# Created on: 29/05/2018

# Take scripts from:
# main_process_68k_pbmc.R
# main_process_pure_pbmc.R

rm(list = ls()) # clear workspace
# ----------------------------
# load relevant libraries
# ----------------------------
library(svd)
library(ggplot2)
library(dplyr)
library(plyr)
library(Rtsne)
library(data.table)
library(Matrix)
library(sparsity)
# -------------------------------------
# specify paths and load functions
# -------------------------------------
FILE_DIR <- "./"
DATA_DIR <- "data/68k"
PROG_DIR <- "single-cell-3prime-paper/pbmc68k_analysis/"
RES_DIR <- "data/68k"

source(file.path(PROG_DIR, "util.R"))
source(file.path(PROG_DIR, "select_pure_pbmc.R"))

# -----------------------
# 1 - load purified data
# -----------------------
print("Loading purified data")
pure_pbmcs <- readRDS(file.path(DATA_DIR, "all_pure_pbmc_data.rds"))
all_data_pure_pbmcs <- pure_pbmcs$all_data
all_json <- pure_pbmcs$all_json
all_metrics <- pure_pbmcs$all_metrics
all_mol_info <- pure_pbmcs$all_mol_info

symbols_pure <- all_data_pure_pbmcs[[1]]$hg19$gene_symbols
genes_pure <- all_data_pure_pbmcs[[1]]$hg19$genes

# -------------------------------------------------------------------------
# downsample mapped reads/cell
# so that all samples have the same # of confidently mapped reads/cell
# -------------------------------------------------------------------------
set.seed(1)
rpc <- all_metrics %>%
    mutate(conf_mapped_rpc = raw_rpc *
        conf_mapped_frac *
        good_bc_frac *
        good_umi_frac) %>%
    select(sample_id, description, conf_mapped_rpc)
tgt_rpc <- floor(min(rpc$conf_mapped_rpc)) # 13995
subsampled_purified_mats <- lapply(1 : length(all_data_pure_pbmcs), function(i) { # subsample the matrix to match tgt_rpc
    cat(sprintf("%d...\n", i))
    .downsample_gene_bc_mtx(all_json[[i]], all_data_pure_pbmcs[[i]], all_mol_info[[i]], tgt_rpc, "conf_mapped_reads")[[1]]
})

# -------------------------------------------
# run PCA on down-sampled data / ignore TSNE
# the steps take a long time to complete
# -------------------------------------------
# using 10 PCs provides enough information to filter out undesired populations from purified populations
print("Computing PCA")
set.seed(1)
all_pure_pca <- lapply(1 : length(subsampled_purified_mats), function(i) pure_pca_i <- .do_propack(subsampled_purified_mats[[i]], 10))
all_pure_tsne <- lapply(1 : length(all_pure_pca), function(i) NULL)

# ------------------------------------
# curate the purified populations
# ------------------------------------
pure_id <- c("CD34+", "CD56+ NK", "CD4+/CD45RA+/CD25- Naive T", "CD4+/CD25 T Reg", "CD8+/CD45RA+ Naive Cytotoxic",
"CD4+/CD45RO+ Memory", "CD8+ Cytotoxic T", "CD19+ B", "CD4+ T Helper2", "CD14+ Monocyte", "Dendritic")
FIG_DIR = NULL
sub_idx <- list(data.frame(sample = 1, use = (get_pure_pop_idx(genes, pure_id[1], all_pure_pca[[1]], all_pure_tsne[[1]], FIG_DIR))),
data.frame(sample = 2, use = (get_pure_pop_idx(genes, pure_id[2], all_pure_pca[[2]], all_pure_tsne[[2]], FIG_DIR))),
data.frame(sample = 3, use = (get_pure_pop_idx(genes, pure_id[3], all_pure_pca[[3]], all_pure_tsne[[3]], FIG_DIR))),
data.frame(sample = 4, use = (get_pure_pop_idx(genes, pure_id[4], all_pure_pca[[4]], all_pure_tsne[[4]], FIG_DIR))),
data.frame(sample = 5, use = (get_pure_pop_idx(genes, pure_id[5], all_pure_pca[[5]], all_pure_tsne[[5]], FIG_DIR))),
data.frame(sample = 6, use = (get_pure_pop_idx(genes, pure_id[6], all_pure_pca[[6]], all_pure_tsne[[6]], FIG_DIR))),
data.frame(sample = 7, use = (get_pure_pop_idx(genes, pure_id[7], all_pure_pca[[7]], all_pure_tsne[[7]], FIG_DIR))),
data.frame(sample = 8, use = (get_pure_pop_idx(genes, pure_id[8], all_pure_pca[[8]], all_pure_tsne[[8]], FIG_DIR))),
data.frame(sample = 9, use = (get_pure_pop_idx(genes, pure_id[9], all_pure_pca[[9]], all_pure_tsne[[9]], FIG_DIR))),
data.frame(sample = 10, use = (get_pure_pop_idx(genes, pure_id[10], all_pure_pca[[10]], all_pure_tsne[[10]], FIG_DIR))),
data.frame(sample = 10, use = (get_pure_pop_idx(genes, pure_id[11], all_pure_pca[[10]], all_pure_tsne[[10]], FIG_DIR))))
pure_select_11 <- lapply(1 : length(sub_idx), function(i) {subsampled_purified_mats[[sub_idx[[i]]$sample[1]]][sub_idx[[i]]$use,]})
pure_use_genes <- which(colSums(do.call(rBind, lapply(pure_select_11, function(x) x))) > 1)
pure_select_use_genes <- lapply(1 : length(pure_select_11), function(i) pure_select_11[[i]][, pure_use_genes])
pure_avg <- do.call(rbind, lapply(pure_select_use_genes, .train_multinomial))

pure_11 <- list(pure_id = pure_id, pure_avg = pure_avg, pure_use_genes = pure_use_genes)

write.csv(pure_11$pure_use_genes, file = file.path(RES_DIR, "pure_use_genes.csv"))

saveRDS(pure_11, file = file.path(RES_DIR, "my_pure_select_11types.rds"))

rm(all_pure_pca)
# dim(all_data_pure_pbmcs[[1]]$hg19$mat) # out : 9232 32738
# dim(subsampled_purified_mats[[1]]) # out : 9232 32738 (unchanged)
# dim(pure_select_11[[1]]) # [1]  6412 32738 - after pca remove cells too distant
# dim(pure_select_use_genes[[1]]) # out : 6412 19841 - remove dimensions where nothing happens

# ------------------------------------------
# fork (my python end)
# ------------------------------------------
m <- pure_select_11[[1]]
labels <- rep(1, dim(pure_select_11[[1]])[[1]])
for (i in 2 : 11) {
    m <- rbind(m, pure_select_11[[i]])
    labels <- c(labels, rep(i, dim(pure_select_11[[i]])[[1]]))
}
write.svmlight(m, labels, file.path(RES_DIR, "pure_full.svmlight"))
rm(m)
rm(labels)
# ------------------------------------------------------------
# 2 - load 68k PBMC data, 11 purified PBMC data and meta-data
# ------------------------------------------------------------
pbmc_68k <- readRDS(file.path(DATA_DIR, "pbmc68k_data.rds")) # pbmc_68k$all_data[[1]]$hg19$genes
all_data_68k <- pbmc_68k$all_data

symbols_68k <- all_data_68k$hg19$gene_symbols
genes_68k <- pbmc_68k$all_data[[1]]$hg19$genes
# ------------------------------------------------------------

for (i in 2 : 10) {
    print(sprintf("Gene Coincidence %i", i))
    print(all(all_data_pure_pbmcs[[1]]$hg19$genes == all_data_pure_pbmcs[[i]]$hg19$genes))
}
print("Checking coincidence pure and 68k")
print(all(genes_68k == genes_pure))

# all_data_68k_select_use_genes <- all_data_68k$hg19$mat[,pure_use_genes]
# all_data_68k_select_use_genes <- as(all_data_68k_select_use_genes, "dgCMatrix")
pure_11 <- readRDS(file.path(RES_DIR, "my_pure_select_11types.rds"))
purified_ref_11 <- load_purified_pbmc_types(pure_11, pbmc_68k$ens_genes)

# ------------------------------------------
# 68k prepro
# normalize by RNA content (umi counts) and select the top 1000 most variable genes
# --------------------------------------------------------------------------------------
# all_data_pure_pbmcs <- read.svmlight(file.path(RES_DIR,"pure.svmlight"))
# pure_11 <- readRDS(file.path(RES_DIR,"my_pure_select_11types.rds"))
# pbmc_68k <- readRDS(file.path(DATA_DIR,"pbmc68k_data.rds")) TODO remove dependence... purified_ref_11 is problematic
print("Normalizing 68k (just for getting labels)")
m <- all_data_68k[[1]]$hg19$mat
l <- .normalize_by_umi(m)
m_n <- l$m
df <- .get_variable_gene(m_n)
disp_cut_off <- sort(df$dispersion_norm, decreasing = T)[1000]
df$used <- df$dispersion_norm >= disp_cut_off

set.seed(0)
m_n_1000 <- m_n[, head(order(- df$dispersion_norm), 1000)]
# ---------------------------------------------------------------------------------------------------------------------------
# assign IDs by comparing the transcriptome profile of each cell to the reference profile from purified PBMC populations
# this produces Fig. 3j in the manuscript
# ---------------------------------------------------------------------------------------------------------------------------
print("Filtering genes (just for getting labels)")
m_filt <- m_n_1000
use_genes_n <- order(- df$dispersion_norm)
m_filt <- m_n[, head(use_genes_n, 1000)]
use_genes_n_id <- symbols_68k[l$use_genes][order(- df$dispersion_norm)]
use_genes_n_ens <- genes_68k[l$use_genes][order(- df$dispersion_norm)]
# first_1000_genes_to_use <- l$use_genes[head(use_genes_n,1000)] #TODO this line is better
print("Compare by correlation")
z_1000_11 <- .compare_by_cor(m_filt, use_genes_n_ens[1 : 1000], purified_ref_11)  # instead of purified_ref_11 # sum(purified_ref_11[1,1:11])
# reassign IDs, as therere some overlaps in the purified pbmc populations
test <- .reassign_pbmc_11(z_1000_11)
cls_id <- factor(colnames(z_1000_11)[test])

labels <- rep(0, length(cls_id))
for (i in 1 : 11) {
    labels[cls_id == pure_id[i]] = i
}

print("Writting 68k_assignments file")
write.svmlight(as(m,"dgCMatrix"), labels, file.path(RES_DIR, "68k_assignments.svmlight"))

# set.seed(0)
# # Only do pca + tsne + plot on the first :4000
# m_n_1000<-m_n[1:4000,head(order(-df$dispersion_norm),1000)]
# pca_n_1000<-.do_propack(m_n_1000,50)
# tsne_n_1000<-Rtsne(pca_n_1000$pca,pca=F)
# tdf_n_1000<-data.frame(tsne_n_1000$Y)
# tdf_n_1000$cls_id<-cls_id[1:4000]
# # adjust ordering of cells for plotting aesthetics
# tdf_mod <- tdf_n_1000[tdf_n_1000$cls_id!='CD4+/CD45RA+/CD25- Naive T',]
# tdf_mod <- rbind(tdf_mod,tdf_n_1000[tdf_n_1000$cls_id=='CD4+/CD45RA+/CD25- Naive T',])
# tdf_mod_2 <- tdf_mod[tdf_mod$cls_id!='CD56+ NK',]
# tdf_mod_2 <- rbind(tdf_mod_2,tdf_mod[tdf_mod$cls_id=='CD56+ NK',])
# out <- ggplot(tdf_mod_2,aes(X1,X2,col=cls_id))+geom_point(size=0,alpha=1)+theme_classic()+.set_pbmc_color_11()
#
# .save_figure(out,'./','pca','tsne',width=8.0,height=4)
