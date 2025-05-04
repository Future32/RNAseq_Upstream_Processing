# 1. 加载包
library(dplyr)
library(org.Hs.eg.db)
library(readr)

# 2. 设置路径和样本
rsem_dir <- "/home/future/rnaseq_batch_250424/08_rsem"
samples <- list.files(rsem_dir, pattern = "*.genes.results") %>%
  gsub(".genes.results", "", .)

files <- file.path(rsem_dir, paste0(samples, ".genes.results"))
names(files) <- samples

# 3. 批量读取所有 RSEM 结果，并整合
rsem_list <- lapply(files, read_tsv)

# 4. 统一提取 Ensembl ID
gene_ids <- rsem_list[[1]]$gene_id
gene_ids_clean <- gsub("\\..*", "", gene_ids)  # 去掉版本号

# 5. Ensembl ID → Gene Symbol 注释
gene_annot <- AnnotationDbi::select(org.Hs.eg.db,
                                    keys = gene_ids_clean,
                                    columns = "SYMBOL",
                                    keytype = "ENSEMBL")

# 6. 合并表达量数据
count_matrix <- do.call(cbind, lapply(rsem_list, function(x) x$expected_count))
colnames(count_matrix) <- samples

tpm_matrix <- do.call(cbind, lapply(rsem_list, function(x) x$TPM))
colnames(tpm_matrix) <- samples

fpkm_matrix <- do.call(cbind, lapply(rsem_list, function(x) x$FPKM))
colnames(fpkm_matrix) <- samples

# 加入 Ensembl ID
count_matrix <- data.frame(ENSEMBL = gene_ids_clean, count_matrix)
tpm_matrix <- data.frame(ENSEMBL = gene_ids_clean, tpm_matrix)
fpkm_matrix <- data.frame(ENSEMBL = gene_ids_clean, fpkm_matrix)

# 合并注释
count_annot <- left_join(gene_annot, count_matrix, by = "ENSEMBL")
tpm_annot <- left_join(gene_annot, tpm_matrix, by = "ENSEMBL")
fpkm_annot <- left_join(gene_annot, fpkm_matrix, by = "ENSEMBL")

# 7. 合并重复Gene Symbol（取平均）
## 先去掉缺失的Symbol
count_annot <- count_annot[!is.na(count_annot$SYMBOL), ]
tpm_annot <- tpm_annot[!is.na(tpm_annot$SYMBOL), ]
fpkm_annot <- fpkm_annot[!is.na(fpkm_annot$SYMBOL), ]

## 使用 dplyr::group_by + summarise_all(mean)
count_final <- count_annot %>%
  group_by(SYMBOL) %>%
  summarise(across(where(is.numeric), mean)) %>%
  as.data.frame()

tpm_final <- tpm_annot %>%
  group_by(SYMBOL) %>%
  summarise(across(where(is.numeric), mean)) %>%
  as.data.frame()

fpkm_final <- fpkm_annot %>%
  group_by(SYMBOL) %>%
  summarise(across(where(is.numeric), mean)) %>%
  as.data.frame()

## 设置基因名为行名
rownames(count_final) <- count_final$SYMBOL
rownames(tpm_final) <- tpm_final$SYMBOL
rownames(fpkm_final) <- fpkm_final$SYMBOL

count_final <- count_final[, -1]
tpm_final <- tpm_final[, -1]
fpkm_final <- fpkm_final[, -1]

# 8. 保存为CSV
write.csv(count_final, "rsem_counts_geneName.csv")
write.csv(tpm_final, "rsem_tpm_geneName.csv")
write.csv(fpkm_final, "rsem_fpkm_geneName.csv")
