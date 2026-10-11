project_root <- normalizePath(Sys.getenv("RECORD_CODE_ROOT", unset = "."), winslash = "/", mustWork = TRUE)
results_dir <- file.path(project_root, "results", "proteomics", "retron_record")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

#1.瀵煎叆鏁版嵁
library(readr)
Retron<- read_tsv(file.path(project_root, "data", "proteomics", "retron_record", "DG4S_scaling.tsv"))
str(Retron)

library(dplyr)

# 1锔忊儯 閲嶅懡鍚嶅垪锛圵T 鈫?cMyc锛孧UT 鈫?Ctrl锛孴ELO 淇濈暀锛?
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

# 2锔忊儯 鏋勫缓 Ctrl銆乧Myc 鍜?TELO 瀛愰泦
Ctrl <- Retron %>%
  select(Protein.Group, Protein.Names, Genes, First.Protein.Description,
         `Ctrl-1`, `Ctrl-2`, `Ctrl-3`)

cMyc <- Retron %>%
  select(Protein.Group, Protein.Names, Genes, First.Protein.Description,
         `cMyc-1`, `cMyc-2`, `cMyc-3`)

TELO <- Retron %>%
  select(Protein.Group, Protein.Names, Genes, First.Protein.Description,
         `TELO-1`, `TELO-2`, `TELO-3`)

# 3锔忊儯 鍚堝苟涓?VennData
VennData <- Ctrl %>%
  left_join(cMyc,
            by = c("Protein.Group","Protein.Names","Genes","First.Protein.Description")) %>%
  left_join(TELO,
            by = c("Protein.Group","Protein.Names","Genes","First.Protein.Description"))

intensity_cols <- c("Ctrl-1", "Ctrl-2", "Ctrl-3",
                    "cMyc-1", "cMyc-2", "cMyc-3",
                    "TELO-1", "TELO-2", "TELO-3")

VennData_raw <- VennData

# Venn plot / detected protein filtering: missing values mean not detected, so fill NA as 0.
VennData <- VennData_raw %>%
  mutate(across(all_of(intensity_cols), ~ ifelse(is.na(.), 0, .)))

# 4锔忊儯 鏌ョ湅缁撴灉
head(VennData)

# 5锔忊儯 淇濆瓨缁撴灉
write_tsv(VennData, file.path(results_dir, "Retron_raw.tsv"))

#绛涢€夌壒寮傛€ц泲鐧?
library(dplyr)
library(readr)

# 1锔忊儯 cMyc 鐗规湁铔嬬櫧锛欳trl 鍏ㄤ负0锛孴ELO 鍏ㄤ负0锛宑Myc 浠绘剰>0
cMyc_unique <- VennData %>%
  filter(`Ctrl-1` == 0 & `Ctrl-2` == 0 & `Ctrl-3` == 0 &
           `TELO-1` == 0 & `TELO-2` == 0 & `TELO-3` == 0 &
           (`cMyc-1` > 0 | `cMyc-2` > 0 | `cMyc-3` > 0))

# 2锔忊儯 TELO 鐗规湁铔嬬櫧锛欳trl 鍏ㄤ负0锛宑Myc 鍏ㄤ负0锛孴ELO 浠绘剰>0
TELO_unique <- VennData %>%
  filter(`Ctrl-1` == 0 & `Ctrl-2` == 0 & `Ctrl-3` == 0 &
           `cMyc-1` == 0 & `cMyc-2` == 0 & `cMyc-3` == 0 &
           (`TELO-1` > 0 | `TELO-2` > 0 | `TELO-3` > 0))

# 3锔忊儯 淇濆瓨缁撴灉
write_tsv(cMyc_unique, file.path(results_dir, "cMyc_unique.tsv"))
write_tsv(TELO_unique, file.path(results_dir, "TELO_unique.tsv"))

# 4锔忊儯 鏌ョ湅鍓嶅嚑琛?
head(cMyc_unique)
head(TELO_unique)

library(dplyr)
library(readr)

# cMyc & TELO 閮芥湁琛ㄨ揪锛孋trl 鍏ㄤ负0
cMyc_TELO_intersect <- VennData %>%
  filter(`Ctrl-1` == 0 & `Ctrl-2` == 0 & `Ctrl-3` == 0 &
           (`cMyc-1` > 0 | `cMyc-2` > 0 | `cMyc-3` > 0) &
           (`TELO-1` > 0 | `TELO-2` > 0 | `TELO-3` > 0))

# 淇濆瓨缁撴灉
write_tsv(cMyc_TELO_intersect, file.path(results_dir, "cMyc_TELO_intersect.tsv"))

# 鏌ョ湅鍓嶅嚑琛?
head(cMyc_TELO_intersect)

library(dplyr)
library(ggvenn)

# ===============================
# 1锔忊儯 鍑嗗铔嬬櫧闆嗗悎
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
# 2锔忊儯 鍒涘缓鍒楄〃锛岀敤浜?ggvenn
# ===============================
venn_list <- list(
  cMyc = cMyc_detected,
  TELO = TELO_detected,
  Ctrl = Ctrl_detected
)

# ===============================
# 3锔忊儯 缁樺埗闊︽仼鍥?
# ===============================
ggvenn(
  venn_list,
  fill_color = c("#F8766D", "#00BA38", "#619CFF"), # 鍙互鎹㈤鑹?
  stroke_size = 0.5,
  set_name_size = 6,
  text_size = 5,
  show_elements = FALSE # 涓嶆樉绀哄厓绱犲悕瀛?
)

#闊︽仼鍥?
library(dplyr)
library(ggvenn)

# ===============================
# 1锔忊儯 鍑嗗铔嬬櫧闆嗗悎
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
# 2锔忊儯 鍒涘缓鍒楄〃
# ===============================
venn_list <- list(
  cMyc = cMyc_detected,
  TELO = TELO_detected,
  Ctrl = Ctrl_detected
)

# ===============================
# 3锔忊儯 缁樺埗缇庡寲鐗堥煢鎭╁浘
# ===============================
ggvenn(
  venn_list,
  fill_color = c("#FF9999", "#99CC99", "#9999FF"), # 棰滆壊鏌斿拰
  stroke_size = 1,
  set_name_size = 6,
  text_size = 5,
  show_elements = FALSE,        # 涓嶆樉绀哄厓绱犲悕
  stroke_color = "black"        # 鍦嗚竟棰滆壊
) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    legend.position = "none"
  )

#鍒楁渶灏忓€煎～鍏?

library(dplyr)
library(readr)

fill_missing_with_column_min <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  min_positive <- min(x[is.finite(x) & x > 0], na.rm = TRUE)
  if (!is.finite(min_positive)) {
    min_positive <- 1
  }
  ifelse(is.na(x) | !is.finite(x) | x <= 0, min_positive, x)
}

# Volcano plot statistics: fill NA / 0 / non-finite intensity values with each column's minimum positive value.
Retron_filled <- VennData_raw %>%
  mutate(across(all_of(intensity_cols), fill_missing_with_column_min))

head(Retron_filled)

write_tsv(Retron_filled, file.path(results_dir, "Retron_filled_min.tsv"))

#鍙杔og2 锛孎C锛宲鍜宎djP
library(dplyr)
library(readr)

# 1锔忊儯 璁＄畻 log2 杞崲
Retron_log2 <- Retron_filled %>%
  mutate(across(c(`Ctrl-1`, `Ctrl-2`, `Ctrl-3`,
                  `cMyc-1`, `cMyc-2`, `cMyc-3`,
                  `TELO-1`, `TELO-2`, `TELO-3`),
                ~ log2(.)))

# 2锔忊儯 璁＄畻姣忕粍骞冲潎鍊?
Retron_stats <- Retron_log2 %>%
  rowwise() %>%
  mutate(
    Ctrl_mean  = mean(c(`Ctrl-1`, `Ctrl-2`, `Ctrl-3`), na.rm = TRUE),
    cMyc_mean  = mean(c(`cMyc-1`, `cMyc-2`, `cMyc-3`), na.rm = TRUE),
    TELO_mean  = mean(c(`TELO-1`, `TELO-2`, `TELO-3`), na.rm = TRUE),
    
    # log2FC锛堜互 Ctrl 涓哄弬鑰冿級
    log2FC_cMyc  = cMyc_mean - Ctrl_mean,
    log2FC_TELO  = TELO_mean - Ctrl_mean
  ) %>%
  ungroup()

# 3锔忊儯 璁＄畻 t 妫€楠?P 鍊?
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

# 4锔忊儯 璁＄畻 padj锛圔enjamini-Hochberg锛?
Retron_stats <- Retron_stats %>%
  mutate(
    padj_cMyc = p.adjust(P_cMyc, method = "BH"),
    padj_TELO = p.adjust(P_TELO, method = "BH")
  )

# 5锔忊儯 淇濆瓨鍒版柊琛?
write_tsv(Retron_stats, file.path(results_dir, "Retron_stats.tsv"))
# Volcano plots: compact publication-style output
library(ggplot2)
library(ggrepel)
library(dplyr)

fc_cutoff <- 1
p_cutoff  <- 0.05
label_n   <- 6
output_dir <- results_dir

make_volcano_data <- function(data, fc_col, p_col) {
  data %>%
    filter(!is.na(.data[[fc_col]]), !is.na(.data[[p_col]]), .data[[p_col]] > 0) %>%
    mutate(
      neg_log10_p = -log10(.data[[p_col]]),
      Significant = case_when(
        .data[[fc_col]] >= fc_cutoff & .data[[p_col]] < p_cutoff ~ "Up",
        .data[[fc_col]] <= -fc_cutoff & .data[[p_col]] < p_cutoff ~ "Down",
        TRUE ~ "Not significant"
      )
    )
}

pick_volcano_labels <- function(data, fc_col, p_col, n = label_n) {
  data %>%
    filter(Significant == "Up", !is.na(Genes), Genes != "") %>%
    arrange(.data[[p_col]], desc(abs(.data[[fc_col]]))) %>%
    distinct(Genes, .keep_all = TRUE) %>%
    slice_head(n = n)
}

plot_volcano <- function(data, fc_col, p_col, title, x_label) {
  plot_data <- make_volcano_data(data, fc_col, p_col)
  label_data <- pick_volcano_labels(plot_data, fc_col, p_col)
  max_abs_fc <- max(abs(plot_data[[fc_col]]), na.rm = TRUE)
  x_limit <- max(fc_cutoff * 1.8, ceiling(max_abs_fc))

  ggplot(plot_data, aes(x = .data[[fc_col]], y = neg_log10_p)) +
    geom_point(aes(color = Significant), alpha = 0.80, size = 1.65, stroke = 0) +
    geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed",
               color = "grey55", linewidth = 0.35) +
    geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed",
               color = "grey55", linewidth = 0.35) +
    geom_text_repel(
      data = label_data,
      aes(label = Genes),
      size = 3.8,
      color = "grey12",
      box.padding = 0.35,
      point.padding = 0.25,
      min.segment.length = 0,
      segment.color = "grey60",
      segment.linewidth = 0.25,
      max.overlaps = Inf,
      seed = 123
    ) +
    scale_color_manual(
      values = c("Up" = "#D84A3A", "Down" = "#2F80B7", "Not significant" = "#C8CCD0"),
      breaks = c("Up", "Down", "Not significant")
    ) +
    coord_cartesian(xlim = c(-x_limit, x_limit), clip = "off") +
    labs(title = title, x = x_label, y = "-log10(p-value)", color = NULL) +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
      axis.title = element_text(size = 13),
      axis.text = element_text(color = "grey25", size = 11),
      legend.position = "top",
      legend.justification = "center",
      legend.text = element_text(size = 11),
      legend.key.width = unit(1.0, "lines"),
      plot.margin = margin(8, 12, 8, 8)
    )
}

cMyc_volcano_data <- make_volcano_data(Retron_stats, "log2FC_cMyc", "P_cMyc")
cMyc_up <- cMyc_volcano_data %>%
  filter(Significant == "Up") %>%
  arrange(P_cMyc)
cMyc_FCpos <- cMyc_volcano_data %>%
  filter(log2FC_cMyc > 0) %>%
  arrange(desc(log2FC_cMyc))

write.csv(cMyc_up, file.path(output_dir, "cMyc_up_proteins.csv"), row.names = FALSE)
write.csv(cMyc_FCpos, file.path(output_dir, "cMyc_FC_positive_proteins.csv"), row.names = FALSE)
write.csv(cMyc_FCpos, file.path(output_dir, "cMyc_UpProteins_all.csv"), row.names = FALSE)

volcano_cMyc <- plot_volcano(
  Retron_stats,
  fc_col = "log2FC_cMyc",
  p_col = "P_cMyc",
  title = "cMyc DG4S possible binding protein",
  x_label = "log2 Fold Change (cMyc / Ctrl)"
)
print(volcano_cMyc)
ggsave(file.path(output_dir, "volcano_cMyc_4x3_large_text.png"), volcano_cMyc,
       width = 6.4, height = 4.8, dpi = 600)

Telo_volcano_data <- make_volcano_data(Retron_stats, "log2FC_TELO", "P_TELO")
Telo_up <- Telo_volcano_data %>%
  filter(Significant == "Up") %>%
  arrange(P_TELO)
Telo_FCpos <- Telo_volcano_data %>%
  filter(log2FC_TELO > 0) %>%
  arrange(desc(log2FC_TELO))

write.csv(Telo_up, file.path(output_dir, "Telo_up_proteins.csv"), row.names = FALSE)
write.csv(Telo_FCpos, file.path(output_dir, "Telo_UpProteins_all.csv"), row.names = FALSE)

volcano_Telo <- plot_volcano(
  Retron_stats,
  fc_col = "log2FC_TELO",
  p_col = "P_TELO",
  title = "Telo DG4S possible binding protein",
  x_label = "log2 Fold Change (Telo / Ctrl)"
)
print(volcano_Telo)
ggsave(file.path(output_dir, "volcano_Telo_4x3_large_text.png"), volcano_Telo,
       width = 6.4, height = 4.8, dpi = 600)







