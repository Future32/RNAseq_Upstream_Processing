# 1. 加载需要的包
library(readr)
library(dplyr)

# 2. 设置路径和样本名
rsem_dir <- "/home/future/rnaseq_batch_250424/08_rsem"
samples <- list.files(rsem_dir, pattern = ".genes.results$", full.names = FALSE)
samples <- gsub(".genes.results", "", samples)

# 3. 读取所有结果
count_list <- list()
tpm_list <- list()
fpkm_list <- list()

for (sample in samples) {
  file <- file.path(rsem_dir, paste0(sample, ".genes.results"))
  dat <- read_tsv(file, show_col_types = FALSE)

  # 提取各类型
  count_list[[sample]] <- dat$expected_count
  tpm_list[[sample]] <- dat$TPM
  fpkm_list[[sample]] <- dat$FPKM
}

# 4. 组合成矩阵
gene_ids <- read_tsv(file.path(rsem_dir, paste0(samples[1], ".genes.results")), show_col_types = FALSE)$gene_id

count_matrix <- do.call(cbind, count_list)
rownames(count_matrix) <- gene_ids
colnames(count_matrix) <- samples

tpm_matrix <- do.call(cbind, tpm_list)
rownames(tpm_matrix) <- gene_ids
colnames(tpm_matrix) <- samples

fpkm_matrix <- do.call(cbind, fpkm_list)
rownames(fpkm_matrix) <- gene_ids
colnames(fpkm_matrix) <- samples

# 5. 保存输出
write.csv(count_matrix, file = "rsem_count_matrix.csv") # RSEM的count不是整数，非DESeq2标准输入
write.csv(tpm_matrix, file = "rsem_tpm_matrix.csv")
write.csv(fpkm_matrix, file = "rsem_fpkm_matrix.csv")
