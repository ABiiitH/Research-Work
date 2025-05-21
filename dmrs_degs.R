# =============================
# 📦 Load Required Packages
# =============================
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

install.packages(c("dplyr", "data.table", "readxl"))
BiocManager::install("impute")
BiocManager::install("limma")

library(tidyverse)
library(data.table)
library(dplyr)
library(readxl)
library(readr)
library(impute)
library(limma)
library(genefilter)

# =============================
# 📂 Load Metadata
# =============================
metadata <- read_excel("data/S1 - Sample Information.xlsx", sheet = "Sheet1")
head(metadata)
table(metadata$type, metadata$ER)

tumor_samples <- metadata %>% filter(type == "TUMOUR") %>% pull(samp)
normal_samples <- metadata %>% filter(type == "ADJNORMAL") %>% pull(samp)

# =============================
# 📂 Load & Clean Methylation Data
# =============================
prom_meth <- read_csv("data/promoter_avg_meth_filt.csv") %>%
  select(-name3.chr)

prom_meth_clean <- prom_meth %>%
  select(-start, -end) %>%
  group_by(name) %>%
  summarise(across(everything(), mean, na.rm = TRUE)) %>%
  ungroup()

meth_matrix <- prom_meth_clean %>%
  column_to_rownames("name") %>%
  as.matrix()

# Filter and Impute Missing Values
meth_matrix_filtered <- meth_matrix[, colMeans(is.na(meth_matrix)) <= 0.8]
imputed <- impute.knn(meth_matrix_filtered, k = 5)
imputed_df <- as.data.frame(imputed$data) %>% rownames_to_column("name")

# =============================
# 🔁 Long Format + Group Labels (Methylation)
# =============================
meth_long_labeled <- imputed_df %>%
  pivot_longer(-name, names_to = "Sample_ID", values_to = "Beta_Methylation") %>%
  rename(Gene = name) %>%
  mutate(Group = case_when(
    Sample_ID %in% tumor_samples ~ "Tumor",
    Sample_ID %in% normal_samples ~ "Normal"
  )) %>%
  filter(!is.na(Group))

# =============================
# 📊 Compute Delta Beta and Identify DMRs
# =============================
beta_avg <- meth_long_labeled %>%
  group_by(Gene, Group) %>%
  summarise(mean_beta = mean(Beta_Methylation), .groups = "drop")

delta_beta_df <- beta_avg %>%
  pivot_wider(names_from = Group, values_from = mean_beta, names_prefix = "mean_beta_") %>%
  mutate(delta_beta = mean_beta_Tumor - mean_beta_Normal)

dmr_genes <- delta_beta_df %>%
  filter(abs(delta_beta) > 0.2) %>%
  arrange(desc(abs(delta_beta)))

# =============================
# 📈 Differential Methylation (limma)
# =============================
group <- ifelse(colnames(meth_matrix) %in% tumor_samples, "Tumor", "Normal") %>% factor()
design <- model.matrix(~ group)

fit <- lmFit(meth_matrix, design)
fit <- eBayes(fit)
topTable <- topTable(fit, coef = "groupTumor", adjust = "fdr", number = Inf) %>%
  rownames_to_column("name")

dmrs_limma <- topTable %>%
  filter(abs(logFC) > 0.2, adj.P.Val < 0.05)

# Compare Limma with delta_beta
comparison_df <- inner_join(delta_beta_df, dmrs_limma, by = "name")

ggplot(comparison_df, aes(x = delta_beta, y = logFC)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  theme_minimal() +
  labs(title = "Δβ vs logFC from limma",
       x = "Δβ (Tumor - Normal)",
       y = "logFC (from limma)")

cor.test(comparison_df$delta_beta, comparison_df$logFC)

# =============================
# 📂 Load & Filter Expression Data
# =============================
expr_mat <- read_csv("data/expression_matrix.csv")
expr_clean <- expr_mat %>% select(name, starts_with("MB_"))

expr_matrix <- expr_clean %>%
  column_to_rownames("name") %>%
  as.matrix()

expr_unlogged <- 2^expr_matrix  # Convert log2 values to original scale

# Filter with genefilter
ffun <- filterfun(pOverA(0.20, 100), cv(0.7, 10))
filter_result <- genefilter(expr_unlogged, ffun)
expr_filtered <- expr_unlogged[filter_result, ]
expr_filtered_log2 <- log2(expr_filtered)

# =============================
# 🔁 Long Format + Group Labels (Expression)
# =============================
expr_long_filtered <- expr_filtered_log2 %>%
  as.data.frame() %>%
  rownames_to_column("Gene") %>%
  pivot_longer(-Gene, names_to = "Sample_ID", values_to = "Expression")

expr_long_labeled <- expr_long_filtered %>%
  mutate(Group = case_when(
    Sample_ID %in% tumor_samples ~ "Tumor",
    Sample_ID %in% normal_samples ~ "Normal"
  )) %>%
  filter(!is.na(Group))

# =============================
# 📊 Compute Delta Expression
# =============================
expr_avg <- expr_long_labeled %>%
  group_by(Gene, Group) %>%
  summarise(mean_expr = mean(Expression, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Group, values_from = mean_expr, names_prefix = "mean_expr_") %>%
  mutate(delta_expr = mean_expr_Tumor - mean_expr_Normal)

# =============================
# 🔁 Merge DMRs and DEGs
# =============================
dmr_expr_merged <- inner_join(dmr_genes, expr_avg, by = "Gene")
cor.test(dmr_expr_merged$delta_beta, dmr_expr_merged$delta_expr, method = "spearman")

# =============================
# 📈 Δβ vs ΔExpr Plot
# =============================
ggplot(dmr_expr_merged, aes(x = delta_beta, y = delta_expr)) +
  geom_point(color = "steelblue", alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "darkred", linetype = "dashed") +
  theme_minimal(base_size = 14) +
  labs(
    title = expression("ΔExpression vs ΔMethylation (Δβ)"),
    x = expression(Delta*"β (Methylation Tumor - Normal)"),
    y = expression(Delta*"Expr (Expression Tumor - Normal)")
  )