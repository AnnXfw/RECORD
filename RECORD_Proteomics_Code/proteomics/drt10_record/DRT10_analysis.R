project_root <- normalizePath(Sys.getenv("RECORD_CODE_ROOT", unset = "."), winslash = "/", mustWork = TRUE)
results_dir <- file.path(project_root, "results", "proteomics", "drt10_record")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

#1.导入数据
library(readr)
DRT10 <- read_tsv(
  file.path(project_root, "data", "proteomics", "drt10_record", "DRT10-RE.txt"),
  locale = locale(encoding = "UTF-16LE")
)

str(DRT10)

#2调整下
library(dplyr)
library(readr)

# 1️⃣ Ctrl 组（重命名）
Ctrl <- DRT10 %>%
  select(Protein.Group, Protein.Names, Genes, First.Protein.Description,
         `dDRT10-1`, `dDRT10-2`, `dDRT10-3`) %>%
  rename(`Ctrl-1` = `dDRT10-1`,
         `Ctrl-2` = `dDRT10-2`,
         `Ctrl-3` = `dDRT10-3`)

# 2️⃣ 实验组（DRT10）——保持原名
Exp <- DRT10 %>%
  select(Protein.Group, Protein.Names, Genes, First.Protein.Description,
         `DRT10-1`, `DRT10-2`, `DRT10-3`)

# 3️⃣ 合并（按蛋白信息对齐）
VennData <- Ctrl %>%
  left_join(Exp, 
            by = c("Protein.Group", "Protein.Names", "Genes", "First.Protein.Description"))

# 4️⃣ 查看整理结果
head(VennData)

# 5️⃣ 保存文件
write_tsv(VennData, file.path(results_dir, "VennData.tsv"))


#3 DRT10特有的蛋白
library(dplyr)
library(readr)

# 提取 DRT10 特有蛋白，保留所有列
DRT10_unique <- VennData %>%
  filter(
    `Ctrl-1` == 0 & `Ctrl-2` == 0 & `Ctrl-3` == 0 &      # Ctrl 全为 0
      (`DRT10-1` > 0 | `DRT10-2` > 0 | `DRT10-3` > 0)     # DRT10 至少一个 > 0
  )

# 导出为新表，保留所有信息
write_tsv(DRT10_unique, file.path(results_dir, "DRT10_unique_proteins_full.tsv"))



#韦恩图
library(dplyr)
library(ggvenn)

# 读取数据
library(readr)

DRT10 <- read_tsv(
  file.path(project_root, "data", "proteomics", "drt10_record", "DRT10-RE.txt"),
  locale = locale(encoding = "UTF-16LE")
)

# 定义两组
Ctrl_cols  <- c("dDRT10-1", "dDRT10-2", "dDRT10-3")
Exp_cols   <- c("DRT10-1",  "DRT10-2",  "DRT10-3")

# 提取两组的蛋白（非零即检测到）
Ctrl_proteins <- DRT10 %>%
  filter(rowSums(select(., all_of(Ctrl_cols)) > 0) > 0) %>%
  pull(Protein.Group) %>%
  unique()

DRT10_proteins <- DRT10 %>%
  filter(rowSums(select(., all_of(Exp_cols)) > 0) > 0) %>%
  pull(Protein.Group) %>%
  unique()

# ---- 计算特有 & 交集 ----
DRT10_unique <- setdiff(DRT10_proteins, Ctrl_proteins)  # 实验组特有
Ctrl_unique  <- setdiff(Ctrl_proteins, DRT10_proteins)  # Ctrl 特有
Intersection <- intersect(DRT10_proteins, Ctrl_proteins) # 交集

# ---- 导出 ----
write_tsv(data.frame(Protein = DRT10_unique),
          file.path(results_dir, "DRT10_unique.tsv"))

write_tsv(data.frame(Protein = Intersection),
          file.path(results_dir, "Intersection.tsv"))

#绘制韦恩图
venn_list <- list(
  DRT10 = DRT10_proteins,
  CTRL = Ctrl_proteins
  
)

ggvenn(venn_list,
       show_elements = FALSE,
       fill_color = c("lightblue", "pink"),
       stroke_size = 0.5,
       set_name_size = 8,
       text_size = 8)


library(ggplot2)

# 构建数据框，所有类别英文
df_venn_bar <- data.frame(
  Category = c("Ctrl Unique", "DRT10 Unique", "Intersection"),
  Count = c(length(Ctrl_unique),
            length(DRT10_unique),
            length(Intersection))
)

# 按数量从高到低排序
df_venn_bar$Category <- factor(df_venn_bar$Category, 
                               levels = df_venn_bar$Category[order(-df_venn_bar$Count)])

# 绘制柱状图
ggplot(df_venn_bar, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 0.7, color = "black", size = 0.5) +  # 宽柱子让间距更近
  geom_text(aes(label = Count), vjust = -0.5, size = 5) +                 # 显示数量
  scale_fill_manual(values = c("Ctrl Unique" = "lightblue",
                               "DRT10 Unique" = "pink",
                               "Intersection" = "lightgreen")) +
  theme_minimal() +
  labs(title = "Protein Count",
       x = "",
       y = "Protein Count") +
  theme(
    text = element_text(size = 14),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)  # 标题居中
  )

library(readr)

# 导出 DRT10 特有蛋白到新表
write_tsv(data.frame(Protein = DRT10_unique),
          "DRT10_unique_proteins.tsv")


#火山图数据处理

library(dplyr)
library(readr)

# 假设 DRT10 已经读入
# DRT10 <- read_tsv(file.path(project_root, "data", "proteomics", "drt10_record", "DRT10-RE.tsv"))

# 选择数值列（比如 dDRT10-1 ~ DRT10-3 对应第7~12列，根据你的表格调整）
num_cols <- 7:12

# 替换 0 值为该列最小非零值
DRT10_filled <- DRT10
for (col in num_cols) {
  non_zero_min <- min(DRT10[[col]][DRT10[[col]] > 0], na.rm = TRUE)
  DRT10_filled[[col]][DRT10_filled[[col]] == 0] <- non_zero_min
}

# 导出为新表
write_tsv(DRT10_filled, file.path(results_dir, "DRT10_zero_filled.tsv"))

library(dplyr)
library(readr)

# 读入已填充零值的表格
# DRT10_filled <- read_tsv(file.path(results_dir, "DRT10_zero_filled.tsv"))

# 定义组列
ctrl_cols <- c("dDRT10-1", "dDRT10-2", "dDRT10-3")
drt_cols  <- c("DRT10-1", "DRT10-2", "DRT10-3")

# log2 转换
DRT10_log2 <- DRT10_filled %>%
  mutate(across(all_of(c(ctrl_cols, drt_cols)), ~ log2(.)))

# 计算均值和 FC
DRT10_log2 <- DRT10_log2 %>%
  rowwise() %>%
  mutate(
    Mean_Ctrl = mean(c_across(all_of(ctrl_cols))),
    Mean_DRT  = mean(c_across(all_of(drt_cols))),
    FC = Mean_DRT - Mean_Ctrl  # log2 fold change
  ) %>%
  ungroup()

# 计算 p 值，恒量行返回 NA
p_values <- apply(DRT10_log2[, c(ctrl_cols, drt_cols)], 1, function(x) {
  ctrl <- x[1:3]
  drt  <- x[4:6]
  if (sd(ctrl) == 0 & sd(drt) == 0) {
    return(NA)  # 恒量行
  } else {
    return(t.test(ctrl, drt)$p.value)
  }
})

# 多重检验校正
adjP <- p.adjust(p_values, method = "BH")  # Benjamini-Hochberg

# 整合结果
DRT10_results <- DRT10_log2 %>%
  select(Protein.Group, Protein.Names, Genes, First.Protein.Description, 
         all_of(ctrl_cols), all_of(drt_cols), Mean_Ctrl, Mean_DRT, FC) %>%
  mutate(p = p_values, adjP = adjP)

# 导出新表
write_tsv(DRT10_results, file.path(results_dir, "DRT10_log2_FC_p_adjP.tsv"))



library(ggplot2)

# 假设 DRT10_results 已经生成
# DRT10_results <- read_tsv(file.path(results_dir, "DRT10_log2_FC_p_adjP.tsv"))

# 创建显著性标签（可根据阈值调整）
DRT10_results <- DRT10_results %>%
  mutate(
    negLog10P = -log10(p),
    Significance = case_when(
      !is.na(p) & FC >= 1 & p < 0.05 ~ "Up",
      !is.na(p) & FC <= -1 & p < 0.05 ~ "Down",
      TRUE ~ "NotSignificant"
    )
  )

# 火山图
library(ggplot2)
library(ggrepel)
library(dplyr)
library(readr)

# 读入结果表
DiffProt <- read_tsv(file.path(results_dir, "DRT10_log2_FC_p_adjP.tsv"))

# 替换 NA p 值，避免 -log10(NA) 报错
DiffProt <- DiffProt %>%
  mutate(pvalue = ifelse(is.na(p), 1, p),
         log2FC = FC)  # 用 FC 列作为 log2FC

# 阈值
fc_cutoff <- 1
p_cutoff  <- 0.05

# 显著性标记
DiffProt$Significant <- "Not Significant"
DiffProt$Significant[DiffProt$log2FC > fc_cutoff & DiffProt$pvalue < p_cutoff] <- "Up"
DiffProt$Significant[DiffProt$log2FC < -fc_cutoff & DiffProt$pvalue < p_cutoff] <- "Down"

# 筛选显著上调前20
top_up <- DiffProt[DiffProt$Significant == "Up", ]
top_up <- top_up[order(-top_up$log2FC), ]
top_up <- head(top_up, 20)

# 标注列
DiffProt$label <- ifelse(DiffProt$Genes %in% top_up$Genes, DiffProt$Genes, NA)

# 火山图
volcano <- ggplot(DiffProt, aes(x = log2FC, y = -log10(pvalue), color = Significant)) +
  geom_point(alpha = 0.7, size = 2.5) +
  scale_color_manual(values = c("Up"="#E64B35", "Down"="#4DBBD5", "Not Significant"="grey80")) +
  geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype="dashed", color="black") +
  geom_hline(yintercept = -log10(p_cutoff), linetype="dashed", color="black") +
  geom_text_repel(aes(label = label),
                  size = 3.5,
                  box.padding = 0.5,
                  point.padding = 0.5,
                  nudge_y = 0.5,
                  max.overlaps = Inf,
                  force = 2) +
  theme_bw(base_size = 14) +
  labs(title="Telo dG4s possible binding Protein",
       x="log2 Fold Change (DRT10 / Ctrl)",
       y="-log10(p-value)") +
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))

print(volcano)



library(dplyr)
library(readr)

# 假设 DiffProt 已经存在
# DiffProt <- read_tsv(file.path(results_dir, "DRT10_log2_FC_p_adjP.tsv"))

# 数值列
num_cols <- c("dDRT10-1","dDRT10-2","dDRT10-3","DRT10-1","DRT10-2","DRT10-3")

# 提取所有上调基因（FC > 0）
all_up_genes <- DiffProt %>%
  filter(FC > 0) %>%
  select(all_of(c("Protein.Group","Protein.Names","Genes","First.Protein.Description",
                  num_cols,"Mean_Ctrl","Mean_DRT","FC","p","adjP")))

# 导出到新表
write_tsv(all_up_genes, file.path(results_dir, "DRT10_All_Upregulated_Genes_log2.tsv"))

#上调基因火山图
library(ggplot2)
library(ggrepel)
library(dplyr)
library(readr)

# 读入结果表
DiffProt <- read_tsv(file.path(results_dir, "DRT10_log2_FC_p_adjP.tsv"))

# 如果 log2FC 列不存在，使用 FC
DiffProt <- DiffProt %>%
  mutate(
    log2FC = FC,
    pvalue = ifelse(is.na(p), 1, p)  # 避免 NA 导致 -log10 报错
  )

# 阈值
fc_cutoff <- 1
p_cutoff  <- 0.05

# 筛选上调蛋白
UpProteins <- DiffProt %>% filter(log2FC > 0)

# 添加显著性标记
UpProteins$Significant <- "Not Significant"
UpProteins$Significant[UpProteins$log2FC > fc_cutoff & UpProteins$pvalue < p_cutoff] <- "Up"

# 标注显著上调前20蛋白
top_up <- UpProteins %>%
  filter(Significant == "Up") %>%
  arrange(desc(log2FC)) %>%
  slice_head(n = 20)

UpProteins$label <- ifelse(UpProteins$Genes %in% top_up$Genes, UpProteins$Genes, NA)

# 火山图绘制
volcano_up <- ggplot(UpProteins, aes(x = log2FC, y = -log10(pvalue), color = Significant)) +
  geom_point(alpha = 0.7, size = 2.5) +
  scale_color_manual(values = c("Up"="#E64B35", "Not Significant"="grey80")) +
  geom_vline(xintercept = fc_cutoff, linetype="dashed", color="black") +
  geom_hline(yintercept = -log10(p_cutoff), linetype="dashed", color="black") +
  geom_text_repel(aes(label = label),
                  size = 3.5,
                  box.padding = 1,
                  point.padding = 0.8,
                  nudge_y = 0.8,
                  min.segment.length = 0.5,
                  segment.size = 0.5,
                  max.overlaps = 20,
                  force = 3) +
  theme_bw(base_size = 14) +
  labs(title="Telo dG4S possible binding protein",
       x="log2 Fold Change (DRT10 / Ctrl)",
       y="-log10(p-value)") +
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))

print(volcano_up)





