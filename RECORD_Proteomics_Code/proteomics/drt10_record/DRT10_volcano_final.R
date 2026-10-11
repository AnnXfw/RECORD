project_root <- normalizePath(Sys.getenv("RECORD_CODE_ROOT", unset = "."), winslash = "/", mustWork = TRUE)
results_dir <- file.path(project_root, "results", "proteomics", "drt10_record")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

#1.瀵煎叆鏁版嵁
library(readr)
DRT10 <- read_tsv(
  file.path(project_root, "data", "proteomics", "drt10_record", "DRT10-RE.txt"),
  locale = locale(encoding = "UTF-16LE")
)

str(DRT10)

#2璋冩暣涓?
library(dplyr)
library(readr)

# 1锔忊儯 Ctrl 缁勶紙閲嶅懡鍚嶏級
Ctrl <- DRT10 %>%
  select(Protein.Group, Protein.Names, Genes, First.Protein.Description,
         `dDRT10-1`, `dDRT10-2`, `dDRT10-3`) %>%
  rename(`Ctrl-1` = `dDRT10-1`,
         `Ctrl-2` = `dDRT10-2`,
         `Ctrl-3` = `dDRT10-3`)

# 2锔忊儯 瀹為獙缁勶紙DRT10锛夆€斺€斾繚鎸佸師鍚?
Exp <- DRT10 %>%
  select(Protein.Group, Protein.Names, Genes, First.Protein.Description,
         `DRT10-1`, `DRT10-2`, `DRT10-3`)

# 3锔忊儯 鍚堝苟锛堟寜铔嬬櫧淇℃伅瀵归綈锛?
VennData <- Ctrl %>%
  left_join(Exp, 
            by = c("Protein.Group", "Protein.Names", "Genes", "First.Protein.Description"))

# 4锔忊儯 鏌ョ湅鏁寸悊缁撴灉
head(VennData)

# 5锔忊儯 淇濆瓨鏂囦欢
write_tsv(VennData, file.path(results_dir, "VennData.tsv"))


#3 DRT10鐗规湁鐨勮泲鐧?
library(dplyr)
library(readr)

# 鎻愬彇 DRT10 鐗规湁铔嬬櫧锛屼繚鐣欐墍鏈夊垪
DRT10_unique <- VennData %>%
  filter(
    `Ctrl-1` == 0 & `Ctrl-2` == 0 & `Ctrl-3` == 0 &      # Ctrl 鍏ㄤ负 0
      (`DRT10-1` > 0 | `DRT10-2` > 0 | `DRT10-3` > 0)     # DRT10 鑷冲皯涓€涓?> 0
  )

# 瀵煎嚭涓烘柊琛紝淇濈暀鎵€鏈変俊鎭?
write_tsv(DRT10_unique, file.path(results_dir, "DRT10_unique_proteins_full.tsv"))



#闊︽仼鍥?
library(dplyr)
library(ggvenn)

# 璇诲彇鏁版嵁
library(readr)

DRT10 <- read_tsv(
  file.path(project_root, "data", "proteomics", "drt10_record", "DRT10-RE.txt"),
  locale = locale(encoding = "UTF-16LE")
)

# 瀹氫箟涓ょ粍
Ctrl_cols  <- c("dDRT10-1", "dDRT10-2", "dDRT10-3")
Exp_cols   <- c("DRT10-1",  "DRT10-2",  "DRT10-3")

# 鎻愬彇涓ょ粍鐨勮泲鐧斤紙闈為浂鍗虫娴嬪埌锛?
Ctrl_proteins <- DRT10 %>%
  filter(rowSums(select(., all_of(Ctrl_cols)) > 0) > 0) %>%
  pull(Protein.Group) %>%
  unique()

DRT10_proteins <- DRT10 %>%
  filter(rowSums(select(., all_of(Exp_cols)) > 0) > 0) %>%
  pull(Protein.Group) %>%
  unique()

# ---- 璁＄畻鐗规湁 & 浜ら泦 ----
DRT10_unique <- setdiff(DRT10_proteins, Ctrl_proteins)  # 瀹為獙缁勭壒鏈?
Ctrl_unique  <- setdiff(Ctrl_proteins, DRT10_proteins)  # Ctrl 鐗规湁
Intersection <- intersect(DRT10_proteins, Ctrl_proteins) # 浜ら泦

# ---- 瀵煎嚭 ----
write_tsv(data.frame(Protein = DRT10_unique),
          file.path(results_dir, "DRT10_unique.tsv"))

write_tsv(data.frame(Protein = Intersection),
          file.path(results_dir, "Intersection.tsv"))

#缁樺埗闊︽仼鍥?
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

# 鏋勫缓鏁版嵁妗嗭紝鎵€鏈夌被鍒嫳鏂?
df_venn_bar <- data.frame(
  Category = c("Ctrl Unique", "DRT10 Unique", "Intersection"),
  Count = c(length(Ctrl_unique),
            length(DRT10_unique),
            length(Intersection))
)

# 鎸夋暟閲忎粠楂樺埌浣庢帓搴?
df_venn_bar$Category <- factor(df_venn_bar$Category, 
                               levels = df_venn_bar$Category[order(-df_venn_bar$Count)])

# 缁樺埗鏌辩姸鍥?
ggplot(df_venn_bar, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 0.7, color = "black", size = 0.5) +  # 瀹芥煴瀛愯闂磋窛鏇磋繎
  geom_text(aes(label = Count), vjust = -0.5, size = 5) +                 # 鏄剧ず鏁伴噺
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
    plot.title = element_text(hjust = 0.5)  # 鏍囬灞呬腑
  )

library(readr)

# 瀵煎嚭 DRT10 鐗规湁铔嬬櫧鍒版柊琛?
write_tsv(data.frame(Protein = DRT10_unique),
          "DRT10_unique_proteins.tsv")


#鐏北鍥炬暟鎹鐞?

library(dplyr)
library(readr)

# 鍋囪 DRT10 宸茬粡璇诲叆
# DRT10 <- read_tsv(file.path(project_root, "data", "proteomics", "drt10_record", "DRT10-RE.tsv"))

# 閫夋嫨鏁板€煎垪锛堟瘮濡?dDRT10-1 ~ DRT10-3 瀵瑰簲绗?~12鍒楋紝鏍规嵁浣犵殑琛ㄦ牸璋冩暣锛?
num_cols <- 7:12

# 鏇挎崲 0 鍊间负璇ュ垪鏈€灏忛潪闆跺€?
DRT10_filled <- DRT10
for (col in num_cols) {
  non_zero_min <- min(DRT10[[col]][DRT10[[col]] > 0], na.rm = TRUE)
  DRT10_filled[[col]][DRT10_filled[[col]] == 0] <- non_zero_min
}

# 瀵煎嚭涓烘柊琛?
write_tsv(DRT10_filled, file.path(results_dir, "DRT10_zero_filled.tsv"))

library(dplyr)
library(readr)

# 璇诲叆宸插～鍏呴浂鍊肩殑琛ㄦ牸
# DRT10_filled <- read_tsv(file.path(results_dir, "DRT10_zero_filled.tsv"))

# 瀹氫箟缁勫垪
ctrl_cols <- c("dDRT10-1", "dDRT10-2", "dDRT10-3")
drt_cols  <- c("DRT10-1", "DRT10-2", "DRT10-3")

# log2 杞崲
DRT10_log2 <- DRT10_filled %>%
  mutate(across(all_of(c(ctrl_cols, drt_cols)), ~ log2(.)))

# 璁＄畻鍧囧€煎拰 FC
DRT10_log2 <- DRT10_log2 %>%
  rowwise() %>%
  mutate(
    Mean_Ctrl = mean(c_across(all_of(ctrl_cols))),
    Mean_DRT  = mean(c_across(all_of(drt_cols))),
    FC = Mean_DRT - Mean_Ctrl  # log2 fold change
  ) %>%
  ungroup()

# 璁＄畻 p 鍊硷紝鎭掗噺琛岃繑鍥?NA
p_values <- apply(DRT10_log2[, c(ctrl_cols, drt_cols)], 1, function(x) {
  ctrl <- x[1:3]
  drt  <- x[4:6]
  if (sd(ctrl) == 0 & sd(drt) == 0) {
    return(NA)  # 鎭掗噺琛?
  } else {
    return(t.test(ctrl, drt)$p.value)
  }
})

# 澶氶噸妫€楠屾牎姝?
adjP <- p.adjust(p_values, method = "BH")  # Benjamini-Hochberg

# 鏁村悎缁撴灉
DRT10_results <- DRT10_log2 %>%
  select(Protein.Group, Protein.Names, Genes, First.Protein.Description, 
         all_of(ctrl_cols), all_of(drt_cols), Mean_Ctrl, Mean_DRT, FC) %>%
  mutate(p = p_values, adjP = adjP)

# 瀵煎嚭鏂拌〃
write_tsv(DRT10_results, file.path(results_dir, "DRT10_log2_FC_p_adjP.tsv"))

# ===============================
# Volcano plot, rewritten cleanly
# ===============================
library(ggplot2)
library(ggrepel)
library(dplyr)
library(readr)

# ---- Easy-to-adjust parameters ----
output_dir <- results_dir
fc_cutoff <- 1
p_cutoff <- 0.05
plot_width <- 6.4
plot_height <- 5.2
plot_dpi <- 600
label_all_up <- TRUE

# ---- Read statistical result table ----
volcano_data <- read_tsv(
  file.path(results_dir, "DRT10_log2_FC_p_adjP.tsv"),
  show_col_types = FALSE
) %>%
  mutate(
    log2FC = FC,
    pvalue = ifelse(is.na(p) | p <= 0, 1, p),
    neg_log10_p = -log10(pvalue),
    Significant = case_when(
      log2FC >= fc_cutoff & pvalue < p_cutoff ~ "Up",
      log2FC <= -fc_cutoff & pvalue < p_cutoff ~ "Down",
      TRUE ~ "Not significant"
    )
  )

# ---- Export all proteins with positive fold change ----
num_cols <- c("dDRT10-1", "dDRT10-2", "dDRT10-3", "DRT10-1", "DRT10-2", "DRT10-3")
all_up_genes <- volcano_data %>%
  filter(FC > 0) %>%
  select(any_of(c("Protein.Group", "Protein.Names", "Genes", "First.Protein.Description",
                  num_cols, "Mean_Ctrl", "Mean_DRT", "FC", "p", "adjP")))
write_tsv(all_up_genes, file.path(results_dir, "DRT10_All_Upregulated_Genes_log2.tsv"))

# ---- Label all significant upregulated genes ----
get_label_data <- function(data) {
  data %>%
    filter(Significant == "Up", !is.na(Genes), Genes != "") %>%
    arrange(pvalue, desc(log2FC)) %>%
    distinct(Genes, .keep_all = TRUE)
}

# ---- Shared plotting function ----
make_volcano <- function(data, title, up_only = FALSE) {
  plot_data <- if (up_only) {
    data %>% filter(log2FC > 0)
  } else {
    data
  }

  label_data <- get_label_data(plot_data)
  max_abs_fc <- max(abs(plot_data$log2FC), na.rm = TRUE)
  x_limit <- max(fc_cutoff * 1.8, ceiling(max_abs_fc)) * 1.08

  base_plot <- ggplot(plot_data, aes(x = log2FC, y = neg_log10_p)) +
    geom_point(aes(color = Significant), alpha = 0.78, size = 1.55, stroke = 0) +
    geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed",
               color = "grey55", linewidth = 0.35) +
    geom_text_repel(
      data = label_data,
      aes(label = Genes),
      size = 3.5,
      color = "grey12",
      box.padding = 0.6,
      point.padding = 0.65,
      min.segment.length = 0,
      force = 6,
      force_pull = 0.12,
      max.time = 3,
      max.iter = 30000,
      segment.color = "grey60",
      segment.linewidth = 0.25,
      max.overlaps = Inf,
      seed = 123
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.04, 0.16))) +
    labs(
      title = title,
      x = "log2 Fold Change (Telo/ CTRL)",
      y = "-log10(p-value)",
      color = NULL
    ) +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
      axis.title = element_text(size = 13),
      axis.text = element_text(color = "grey25", size = 11),
      legend.position = "top",
      legend.justification = "center",
      legend.text = element_text(size = 11),
      legend.key.width = unit(1.0, "lines"),
      aspect.ratio = 3 / 4,
      plot.margin = margin(8, 12, 8, 8)
    )

  if (up_only) {
    base_plot +
      geom_vline(xintercept = fc_cutoff, linetype = "dashed", color = "grey55", linewidth = 0.35) +
      scale_color_manual(values = c("Up" = "#D84A3A", "Not significant" = "#C8CCD0"),
                         breaks = c("Up", "Not significant")) +
      coord_cartesian(xlim = c(0, x_limit), clip = "off")
  } else {
    base_plot +
      geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed",
                 color = "grey55", linewidth = 0.35) +
      scale_color_manual(values = c("Up" = "#D84A3A", "Down" = "#2F80B7", "Not significant" = "#C8CCD0"),
                         breaks = c("Up", "Down", "Not significant")) +
      coord_cartesian(xlim = c(-x_limit, x_limit), clip = "off")
  }
}

# ---- Draw and save plots ----
volcano_all <- make_volcano(volcano_data, "Telo dG4 enriched protein", up_only = FALSE)
print(volcano_all)
ggsave(file.path(output_dir, "DRT10_volcano_all_up_labeled.TIF"), volcano_all,
       width = plot_width, height = plot_height, dpi = 1200)

volcano_up <- make_volcano(volcano_data, "Telo dG4 enriched protein", up_only = TRUE)
print(volcano_up)
ggsave(file.path(output_dir, "DRT10_up_volcano_all_labeled.TIF"), volcano_up,
       width = plot_width, height = plot_height, dpi = 1200)
