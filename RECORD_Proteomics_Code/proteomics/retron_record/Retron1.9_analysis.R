project_root <- normalizePath(Sys.getenv("RECORD_CODE_ROOT", unset = "."), winslash = "/", mustWork = TRUE)
results_dir <- file.path(project_root, "results", "proteomics", "retron_record")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

#1.导入数据
library(readr)
Retron<- read_tsv(file.path(project_root, "data", "proteomics", "retron_record", "DG4S_scaling.tsv"))
str(Retron)

library(dplyr)

# 1️⃣ 重命名列（WT → cMyc，MUT → Ctrl，TELO 保留）
# Normalize DIA-NN sample headers to filenames on both Windows and Linux.
sample_headers <- basename(gsub("\\\\", "/", names(Retron)))
if (anyDuplicated(sample_headers)) {
  stop("Sample-header normalization produced duplicate column names; check the input table.")
}
names(Retron) <- sample_headers

Retron <- Retron %>%
  rename(
    `cMyc-1` = `WT-1.raw`,
    `cMyc-2` = `WT-2.raw`,
    `cMyc-3` = `WT-3.raw`,
    
    `TELO-1` = `TELO-1.raw`,
    `TELO-2` = `TELO-2.raw`,
    `TELO-3` = `TELO-3.raw`,
    
    `Ctrl-1` = `MUT-1.raw`,
    `Ctrl-2` = `MUT-2.raw`,
    `Ctrl-3` = `MUT-3.raw`
  )

# 2️⃣ 构建 Ctrl、cMyc 和 TELO 子集
Ctrl <- Retron %>%
  select(Protein.Group, Protein.Names, Genes, First.Protein.Description,
         `Ctrl-1`, `Ctrl-2`, `Ctrl-3`)

cMyc <- Retron %>%
  select(Protein.Group, Protein.Names, Genes, First.Protein.Description,
         `cMyc-1`, `cMyc-2`, `cMyc-3`)

TELO <- Retron %>%
  select(Protein.Group, Protein.Names, Genes, First.Protein.Description,
         `TELO-1`, `TELO-2`, `TELO-3`)

# 3️⃣ 合并为 VennData
VennData <- Ctrl %>%
  left_join(cMyc,
            by = c("Protein.Group","Protein.Names","Genes","First.Protein.Description")) %>%
  left_join(TELO,
            by = c("Protein.Group","Protein.Names","Genes","First.Protein.Description"))

# 4️⃣ 查看结果
head(VennData)

# 5️⃣ 保存结果
write_tsv(VennData, file.path(results_dir, "Retron_raw.tsv"))

#筛选特异性蛋白
library(dplyr)
library(readr)

# 1️⃣ cMyc 特有蛋白：Ctrl 全为0，TELO 全为0，cMyc 任意>0
cMyc_unique <- VennData %>%
  filter(`Ctrl-1` == 0 & `Ctrl-2` == 0 & `Ctrl-3` == 0 &
           `TELO-1` == 0 & `TELO-2` == 0 & `TELO-3` == 0 &
           (`cMyc-1` > 0 | `cMyc-2` > 0 | `cMyc-3` > 0))

# 2️⃣ TELO 特有蛋白：Ctrl 全为0，cMyc 全为0，TELO 任意>0
TELO_unique <- VennData %>%
  filter(`Ctrl-1` == 0 & `Ctrl-2` == 0 & `Ctrl-3` == 0 &
           `cMyc-1` == 0 & `cMyc-2` == 0 & `cMyc-3` == 0 &
           (`TELO-1` > 0 | `TELO-2` > 0 | `TELO-3` > 0))

# 3️⃣ 保存结果
write_tsv(cMyc_unique, file.path(results_dir, "cMyc_unique.tsv"))
write_tsv(TELO_unique, file.path(results_dir, "TELO_unique.tsv"))

# 4️⃣ 查看前几行
head(cMyc_unique)
head(TELO_unique)

library(dplyr)
library(readr)

# cMyc & TELO 都有表达，Ctrl 全为0
cMyc_TELO_intersect <- VennData %>%
  filter(`Ctrl-1` == 0 & `Ctrl-2` == 0 & `Ctrl-3` == 0 &
           (`cMyc-1` > 0 | `cMyc-2` > 0 | `cMyc-3` > 0) &
           (`TELO-1` > 0 | `TELO-2` > 0 | `TELO-3` > 0))

# 保存结果
write_tsv(cMyc_TELO_intersect, file.path(results_dir, "cMyc_TELO_intersect.tsv"))

# 查看前几行
head(cMyc_TELO_intersect)

library(dplyr)
library(ggvenn)

# ===============================
# 1️⃣ 准备蛋白集合
# ===============================
cMyc_detected <- VennData %>%
  filter(`cMyc-1` > 0 | `cMyc-2` > 0 | `cMyc-3` > 0) %>%
  pull(Genes)

TELO_detected <- VennData %>%
  filter(`TELO-1` > 0 | `TELO-2` > 0 | `TELO-3` > 0) %>%
  pull(Genes)

Ctrl_detected <- VennData %>%
  filter(`Ctrl-1` > 0 | `Ctrl-2` > 0 | `Ctrl-3` > 0) %>%
  pull(Genes)

# ===============================
# 2️⃣ 创建列表，用于 ggvenn
# ===============================
venn_list <- list(
  cMyc = cMyc_detected,
  TELO = TELO_detected,
  Ctrl = Ctrl_detected
)

# ===============================
# 3️⃣ 绘制韦恩图
# ===============================
ggvenn(
  venn_list,
  fill_color = c("#F8766D", "#00BA38", "#619CFF"), # 可以换颜色
  stroke_size = 0.5,
  set_name_size = 6,
  text_size = 5,
  show_elements = FALSE # 不显示元素名字
)

#韦恩图2
library(dplyr)
library(ggvenn)

# ===============================
# 1️⃣ 准备蛋白集合
# ===============================
cMyc_detected <- VennData %>%
  filter(`cMyc-1` > 0 | `cMyc-2` > 0 | `cMyc-3` > 0) %>%
  pull(Genes)

TELO_detected <- VennData %>%
  filter(`TELO-1` > 0 | `TELO-2` > 0 | `TELO-3` > 0) %>%
  pull(Genes)

Ctrl_detected <- VennData %>%
  filter(`Ctrl-1` > 0 | `Ctrl-2` > 0 | `Ctrl-3` > 0) %>%
  pull(Genes)

# ===============================
# 2️⃣ 创建列表
# ===============================
venn_list <- list(
  cMyc = cMyc_detected,
  TELO = TELO_detected,
  Ctrl = Ctrl_detected
)

# ===============================
# 3️⃣ 绘制美化版韦恩图
# ===============================
ggvenn(
  venn_list,
  fill_color = c("#FF9999", "#99CC99", "#9999FF"), # 颜色柔和
  stroke_size = 1,
  set_name_size = 6,
  text_size = 5,
  show_elements = FALSE,        # 不显示元素名
  stroke_color = "black"        # 圆边颜色
) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    legend.position = "none"
  )

#列最小值填充

library(dplyr)
library(readr)

Retron_filled <- VennData %>%
  mutate(across(c(`Ctrl-1`, `Ctrl-2`, `Ctrl-3`,
                  `cMyc-1`, `cMyc-2`, `cMyc-3`,
                  `TELO-1`, `TELO-2`, `TELO-3`),
                ~ ifelse(. == 0, min(.[. > 0], na.rm = TRUE), .)))

head(Retron_filled)

write_tsv(Retron_filled, file.path(results_dir, "Retron_filled_min.tsv"))

#取log2 ，FC，p和adjP
library(dplyr)
library(readr)

# 1️⃣ 计算 log2 转换
Retron_log2 <- Retron_filled %>%
  mutate(across(c(`Ctrl-1`, `Ctrl-2`, `Ctrl-3`,
                  `cMyc-1`, `cMyc-2`, `cMyc-3`,
                  `TELO-1`, `TELO-2`, `TELO-3`),
                ~ log2(.)))

# 2️⃣ 计算每组平均值
Retron_stats <- Retron_log2 %>%
  rowwise() %>%
  mutate(
    Ctrl_mean  = mean(c(`Ctrl-1`, `Ctrl-2`, `Ctrl-3`), na.rm = TRUE),
    cMyc_mean  = mean(c(`cMyc-1`, `cMyc-2`, `cMyc-3`), na.rm = TRUE),
    TELO_mean  = mean(c(`TELO-1`, `TELO-2`, `TELO-3`), na.rm = TRUE),
    
    # log2FC（以 Ctrl 为参考）
    log2FC_cMyc  = cMyc_mean - Ctrl_mean,
    log2FC_TELO  = TELO_mean - Ctrl_mean
  ) %>%
  ungroup()

# 3️⃣ 计算 t 检验 P 值
Retron_stats <- Retron_stats %>%
  rowwise() %>%
  mutate(
    P_cMyc = tryCatch(t.test(c(`cMyc-1`,`cMyc-2`,`cMyc-3`),
                             c(`Ctrl-1`,`Ctrl-2`,`Ctrl-3`))$p.value,
                      error = function(e) NA),
    P_TELO = tryCatch(t.test(c(`TELO-1`,`TELO-2`,`TELO-3`),
                             c(`Ctrl-1`,`Ctrl-2`,`Ctrl-3`))$p.value,
                      error = function(e) NA)
  ) %>%
  ungroup()

# 4️⃣ 计算 padj（Benjamini-Hochberg）
Retron_stats <- Retron_stats %>%
  mutate(
    padj_cMyc = p.adjust(P_cMyc, method = "BH"),
    padj_TELO = p.adjust(P_TELO, method = "BH")
  )

# 5️⃣ 保存到新表
write_tsv(Retron_stats, file.path(results_dir, "Retron_stats.tsv"))
library(ggplot2)
library(ggrepel)
library(dplyr)

# 阈值
fc_cutoff <- 1
p_cutoff  <- 0.05

# 去掉 NA
Retron_stats <- Retron_stats %>% 
  filter(!is.na(log2FC_cMyc) & !is.na(P_cMyc))

# 添加显著性标记（只使用原始 p 值）
Retron_stats <- Retron_stats %>%
  mutate(Significant = case_when(
    log2FC_cMyc > fc_cutoff & P_cMyc < p_cutoff ~ "Up",
    TRUE                                        ~ "Not Significant"
  ))

# 选前20个上调蛋白作为标签
top_up <- Retron_stats %>%
  filter(Significant == "Up") %>%
  arrange(P_cMyc) %>%
  slice_head(n = 20)

# 绘制火山图
volcano_cMyc <- ggplot(Retron_stats, aes(x = log2FC_cMyc, y = -log10(P_cMyc), color = Significant)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c("Up" = "#E64B35", "Not Significant" = "grey80")) +
  geom_vline(xintercept = fc_cutoff, linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed", color = "black") +
  geom_text_repel(
    data = top_up,
    aes(label = Genes),
    size = 3.2,
    box.padding = 0.6,
    point.padding = 0.6,
    nudge_y = 0.5,
    min.segment.length = 0.3,
    segment.size = 0.4,
    force = 2,
    max.overlaps = 30
  ) +
  theme_bw(base_size = 14) +
  labs(
    title = "cMyc DG4S possible binding protein",
    x = "log2 Fold Change (cMyc / Ctrl)",
    y = "-log10(p-value)"
  ) +
  theme(
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )

print(volcano_cMyc)

# 提取 cMyc 上调蛋白
cMyc_up <- Retron_stats %>%
  filter(Significant == "Up") %>%
  arrange(P_cMyc)  # 按 p 值排序，可选

# 查看前几行
head(cMyc_up)

# 如果需要保存到本地 CSV
write.csv(cMyc_up, "cMyc_up_proteins.csv", row.names = FALSE)

# 提取所有 cMyc 上调蛋白（FC > 0）
cMyc_FCpos <- Retron_stats %>%
  filter(log2FC_cMyc > 0) %>%
  arrange(desc(log2FC_cMyc))  # 按 log2FC 降序排列，可选

# 查看前几行
head(cMyc_FCpos)

# 保存到 CSV 文件（可选）
write.csv(cMyc_FCpos, "cMyc_FC_positive_proteins.csv", row.names = FALSE)

#cMyc上调火山图
library(ggplot2)
library(ggrepel)
library(dplyr)
library(readr)

# 阈值
fc_cutoff <- 1
p_cutoff  <- 0.05

# 去掉 NA
Retron_stats <- Retron_stats %>% 
  filter(!is.na(log2FC_cMyc) & !is.na(P_cMyc))

# 筛选上调蛋白（log2FC > 0）
cMyc_Up <- Retron_stats %>% filter(log2FC_cMyc > 0)

# 添加显著性标记
cMyc_Up <- cMyc_Up %>%
  mutate(Significant = ifelse(log2FC_cMyc > fc_cutoff & P_cMyc < p_cutoff, "Up", "Not Significant"))

# 保存所有上调蛋白到新表
write.csv(cMyc_Up, "cMyc_UpProteins_all.csv", row.names = FALSE)

# 标注显著上调前20蛋白
top_up <- cMyc_Up %>%
  filter(Significant == "Up") %>%
  arrange(desc(log2FC_cMyc)) %>%
  slice_head(n = 20)

cMyc_Up$label <- ifelse(cMyc_Up$Genes %in% top_up$Genes, cMyc_Up$Genes, NA)

# 火山图绘制
volcano_cMyc <- ggplot(cMyc_Up, aes(x = log2FC_cMyc, y = -log10(P_cMyc), color = Significant)) +
  geom_point(alpha = 0.7, size = 2.5) +
  scale_color_manual(values = c("Up"="#E64B35", "Not Significant"="grey80")) +
  geom_vline(xintercept = fc_cutoff, linetype="dashed", color="black", size=0.8, alpha=0.7) +
  geom_hline(yintercept = -log10(p_cutoff), linetype="dashed", color="black", size=0.8, alpha=0.7) +
  geom_text_repel(aes(label = label),
                  size = 3.5,
                  box.padding = 0.8,
                  point.padding = 0.5,
                  nudge_y = 0.5,
                  min.segment.length = 0.3,
                  segment.size = 0.4,
                  max.overlaps = 30,
                  force = 2) +
  theme_bw(base_size = 14) +
  labs(title="cMyc DG4S possible binding protein",
       x="log2 Fold Change (cMyc / Ctrl)",
       y="-log10(p-value)") +
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"))

# 打印火山图
print(volcano_cMyc)



#Telo火山图
library(ggplot2)
library(ggrepel)
library(dplyr)

# 阈值
fc_cutoff <- 1
p_cutoff  <- 0.05

# 去掉 NA
Retron_stats <- Retron_stats %>% 
  filter(!is.na(log2FC_TELO) & !is.na(P_TELO))

# 添加显著性标记（全局）
Retron_stats <- Retron_stats %>%
  mutate(Significant = ifelse(log2FC_TELO > fc_cutoff & P_TELO < p_cutoff, "Up", "Not Significant"))

# 选前20个显著上调蛋白作为标签
top_up <- Retron_stats %>%
  filter(Significant == "Up") %>%
  arrange(P_TELO) %>%
  slice_head(n = 20)

# 火山图绘制（只标上调标签）
volcano_Telo_global <- ggplot(Retron_stats, aes(x = log2FC_TELO, y = -log10(P_TELO), color = Significant)) +
  geom_point(alpha = 0.7, size = 2.5) +
  scale_color_manual(values = c("Up"="#E64B35", "Not Significant"="grey80")) +
  geom_vline(xintercept = fc_cutoff, linetype="dashed", color="black", size=0.8, alpha=0.7) +
  geom_hline(yintercept = -log10(p_cutoff), linetype="dashed", color="black", size=0.8, alpha=0.7) +
  geom_text_repel(data = top_up,
                  aes(label = Genes),
                  size = 3.5,
                  box.padding = 0.8,
                  point.padding = 0.5,
                  nudge_y = 0.5,
                  min.segment.length = 0.3,
                  segment.size = 0.4,
                  max.overlaps = 30,
                  force = 2) +
  theme_bw(base_size = 14) +
  labs(title="Telo DG4S global volcano plot",
       x="log2 Fold Change (Telo / Ctrl)",
       y="-log10(p-value)") +
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"))

# 打印火山图
print(volcano_Telo_global)

#Telo上调火山图
library(ggplot2)
library(ggrepel)
library(dplyr)

# 阈值
fc_cutoff <- 1
p_cutoff  <- 0.05

# 去掉 NA
Retron_stats <- Retron_stats %>% 
  filter(!is.na(log2FC_TELO) & !is.na(P_TELO))

# 筛选上调蛋白（log2FC > 0）
Telo_Up <- Retron_stats %>% filter(log2FC_TELO > 0)

# 添加显著性标记
Telo_Up <- Telo_Up %>%
  mutate(Significant = ifelse(log2FC_TELO > fc_cutoff & P_TELO < p_cutoff, "Up", "Not Significant"))

# 保存所有上调蛋白到新表（可选）
write.csv(Telo_Up, "Telo_UpProteins_all.csv", row.names = FALSE)

# 标注显著上调前20蛋白
top_up <- Telo_Up %>%
  filter(Significant == "Up") %>%
  arrange(desc(log2FC_TELO)) %>%
  slice_head(n = 20)

Telo_Up$label <- ifelse(Telo_Up$Genes %in% top_up$Genes, Telo_Up$Genes, NA)

# 火山图绘制
volcano_Telo <- ggplot(Telo_Up, aes(x = log2FC_TELO, y = -log10(P_TELO), color = Significant)) +
  geom_point(alpha = 0.7, size = 2.5) +
  scale_color_manual(values = c("Up"="#E64B35", "Not Significant"="grey80")) +
  geom_vline(xintercept = fc_cutoff, linetype="dashed", color="black", size=0.8, alpha=0.7) +
  geom_hline(yintercept = -log10(p_cutoff), linetype="dashed", color="black", size=0.8, alpha=0.7) +
  geom_text_repel(aes(label = label),
                  size = 3.5,
                  box.padding = 0.8,
                  point.padding = 0.5,
                  nudge_y = 0.5,
                  min.segment.length = 0.3,
                  segment.size = 0.4,
                  max.overlaps = 30,
                  force = 2) +
  theme_bw(base_size = 14) +
  labs(title="Telo DG4S possible binding protein",
       x="log2 Fold Change (Telo / Ctrl)",
       y="-log10(p-value)") +
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"))

# 打印火山图
print(volcano_Telo)




