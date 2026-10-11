project_root <- normalizePath(Sys.getenv("RECORD_CODE_ROOT", unset = "."), winslash = "/", mustWork = TRUE)
results_dir <- file.path(project_root, "results", "proteomics", "retron_record")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

#1.鐎电厧鍙嗛弫鐗堝祦
library(readr)
Retron<- read_tsv(file.path(project_root, "data", "proteomics", "retron_record", "DG4S_scaling.tsv"))
str(Retron)

library(dplyr)

# 1閿斿繆鍎?闁插秴鎳￠崥宥呭灙閿涘湹T 閳?cMyc閿涘UT 閳?Ctrl閿涘ELO 娣囨繄鏆€閿?
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

# 2閿斿繆鍎?閺嬪嫬缂?Ctrl閵嗕恭Myc 閸?TELO 鐎涙劙娉?
Ctrl <- Retron %>%
  select(Protein.Group, Protein.Names, Genes, First.Protein.Description,
         `Ctrl-1`, `Ctrl-2`, `Ctrl-3`)

cMyc <- Retron %>%
  select(Protein.Group, Protein.Names, Genes, First.Protein.Description,
         `cMyc-1`, `cMyc-2`, `cMyc-3`)

TELO <- Retron %>%
  select(Protein.Group, Protein.Names, Genes, First.Protein.Description,
         `TELO-1`, `TELO-2`, `TELO-3`)

# 3閿斿繆鍎?閸氬牆鑻熸稉?VennData
VennData <- Ctrl %>%
  left_join(cMyc,
            by = c("Protein.Group","Protein.Names","Genes","First.Protein.Description")) %>%
  left_join(TELO,
            by = c("Protein.Group","Protein.Names","Genes","First.Protein.Description"))

# 4閿斿繆鍎?閺屻儳婀呯紒鎾寸亯
head(VennData)

# 5閿斿繆鍎?娣囨繂鐡ㄧ紒鎾寸亯
write_tsv(VennData, file.path(results_dir, "Retron_raw.tsv"))

#缁涙盯鈧澹掑鍌涒偓褑娉查惂?
library(dplyr)
library(readr)

# 1閿斿繆鍎?cMyc 閻楄婀侀摂瀣閿涙trl 閸忋劋璐?閿涘ELO 閸忋劋璐?閿涘畱Myc 娴犵粯鍓?0
cMyc_unique <- VennData %>%
  filter(`Ctrl-1` == 0 & `Ctrl-2` == 0 & `Ctrl-3` == 0 &
           `TELO-1` == 0 & `TELO-2` == 0 & `TELO-3` == 0 &
           (`cMyc-1` > 0 | `cMyc-2` > 0 | `cMyc-3` > 0))

# 2閿斿繆鍎?TELO 閻楄婀侀摂瀣閿涙trl 閸忋劋璐?閿涘畱Myc 閸忋劋璐?閿涘ELO 娴犵粯鍓?0
TELO_unique <- VennData %>%
  filter(`Ctrl-1` == 0 & `Ctrl-2` == 0 & `Ctrl-3` == 0 &
           `cMyc-1` == 0 & `cMyc-2` == 0 & `cMyc-3` == 0 &
           (`TELO-1` > 0 | `TELO-2` > 0 | `TELO-3` > 0))

# 3閿斿繆鍎?娣囨繂鐡ㄧ紒鎾寸亯
write_tsv(cMyc_unique, file.path(results_dir, "cMyc_unique.tsv"))
write_tsv(TELO_unique, file.path(results_dir, "TELO_unique.tsv"))

# 4閿斿繆鍎?閺屻儳婀呴崜宥呭殤鐞?
head(cMyc_unique)
head(TELO_unique)

library(dplyr)
library(readr)

# cMyc & TELO 闁姤婀佺悰銊ㄦ彧閿涘瓔trl 閸忋劋璐?
cMyc_TELO_intersect <- VennData %>%
  filter(`Ctrl-1` == 0 & `Ctrl-2` == 0 & `Ctrl-3` == 0 &
           (`cMyc-1` > 0 | `cMyc-2` > 0 | `cMyc-3` > 0) &
           (`TELO-1` > 0 | `TELO-2` > 0 | `TELO-3` > 0))

# 娣囨繂鐡ㄧ紒鎾寸亯
write_tsv(cMyc_TELO_intersect, file.path(results_dir, "cMyc_TELO_intersect.tsv"))

# 閺屻儳婀呴崜宥呭殤鐞?
head(cMyc_TELO_intersect)

library(dplyr)
library(ggvenn)

# ===============================
# 1閿斿繆鍎?閸戝棗顦摂瀣闂嗗棗鎮?
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
# 2閿斿繆鍎?閸掓稑缂撻崚妤勩€冮敍宀€鏁ゆ禍?ggvenn
# ===============================
venn_list <- list(
  cMyc = cMyc_detected,
  TELO = TELO_detected,
  Ctrl = Ctrl_detected
)

# ===============================
# 3閿斿繆鍎?缂佹ê鍩楅棅锔戒患閸?
# ===============================
ggvenn(
  venn_list,
  fill_color = c("#F8766D", "#00BA38", "#619CFF"), # 閸欘垯浜掗幑銏ゎ杹閼?
  stroke_size = 0.5,
  set_name_size = 6,
  text_size = 5,
  show_elements = FALSE # 娑撳秵妯夌粈鍝勫帗缁辩姴鎮曠€?
)

#闂婏附浠奸崶?
library(dplyr)
library(ggvenn)

# ===============================
# 1閿斿繆鍎?閸戝棗顦摂瀣闂嗗棗鎮?
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
# 2閿斿繆鍎?閸掓稑缂撻崚妤勩€?
# ===============================
venn_list <- list(
  cMyc = cMyc_detected,
  TELO = TELO_detected,
  Ctrl = Ctrl_detected
)

# ===============================
# 3閿斿繆鍎?缂佹ê鍩楃紘搴″閻楀牓鐓㈤幁鈺佹禈
# ===============================
ggvenn(
  venn_list,
  fill_color = c("#FF9999", "#99CC99", "#9999FF"), # 妫版粏澹婇弻鏂挎嫲
  stroke_size = 1,
  set_name_size = 6,
  text_size = 5,
  show_elements = FALSE,        # 娑撳秵妯夌粈鍝勫帗缁辩姴鎮?
  stroke_color = "black"        # 閸﹀棜绔熸０婊嗗
) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    legend.position = "none"
  )

#閸掓娓剁亸蹇撯偓鐓庯綖閸?

library(dplyr)
library(readr)

Retron_filled <- VennData %>%
  mutate(across(c(`Ctrl-1`, `Ctrl-2`, `Ctrl-3`,
                  `cMyc-1`, `cMyc-2`, `cMyc-3`,
                  `TELO-1`, `TELO-2`, `TELO-3`),
                ~ ifelse(. == 0, min(.[. > 0], na.rm = TRUE), .)))

head(Retron_filled)

write_tsv(Retron_filled, file.path(results_dir, "Retron_filled_min.tsv"))

#閸欐潝og2 閿涘瓗C閿涘閸滃畮djP
library(dplyr)
library(readr)

# 1閿斿繆鍎?鐠侊紕鐣?log2 鏉烆剚宕?
Retron_log2 <- Retron_filled %>%
  mutate(across(c(`Ctrl-1`, `Ctrl-2`, `Ctrl-3`,
                  `cMyc-1`, `cMyc-2`, `cMyc-3`,
                  `TELO-1`, `TELO-2`, `TELO-3`),
                ~ log2(.)))

# 2閿斿繆鍎?鐠侊紕鐣诲В蹇曠矋楠炲啿娼庨崐?
Retron_stats <- Retron_log2 %>%
  rowwise() %>%
  mutate(
    Ctrl_mean  = mean(c(`Ctrl-1`, `Ctrl-2`, `Ctrl-3`), na.rm = TRUE),
    cMyc_mean  = mean(c(`cMyc-1`, `cMyc-2`, `cMyc-3`), na.rm = TRUE),
    TELO_mean  = mean(c(`TELO-1`, `TELO-2`, `TELO-3`), na.rm = TRUE),
    
    # log2FC閿涘牅浜?Ctrl 娑撳搫寮懓鍐跨礆
    log2FC_cMyc  = cMyc_mean - Ctrl_mean,
    log2FC_TELO  = TELO_mean - Ctrl_mean
  ) %>%
  ungroup()

# 3閿斿繆鍎?鐠侊紕鐣?t 濡偓妤?P 閸?
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

# 4閿斿繆鍎?鐠侊紕鐣?padj閿涘湐enjamini-Hochberg閿?
Retron_stats <- Retron_stats %>%
  mutate(
    padj_cMyc = p.adjust(P_cMyc, method = "BH"),
    padj_TELO = p.adjust(P_TELO, method = "BH")
  )

# 5閿斿繆鍎?娣囨繂鐡ㄩ崚鐗堟煀鐞?
write_tsv(Retron_stats, file.path(results_dir, "Retron_stats.tsv"))
library(ggplot2)
library(ggrepel)
library(dplyr)

# Volcano plot beautification only. Upstream data processing is unchanged.
fc_cutoff <- 1
p_cutoff  <- 0.05
label_n   <- 6
force_label_genes_cMyc <- c("ACD", "ZNF420", "ZNF616", "GLI4", "DDX28", "ZNF808", "UP11", "ZNF564", "ZNF57", "MRM1")
force_label_genes_Telo <- c("ACTL6B", "ZNF320", "TRUB2", "GLI4", "HEY1", "KLF16", "ING4", "ZNF616", "NCOA4", "TEAD4", "DDX28", "ASH1L", "ZBTB12", "ZBTB420")
output_dir <- results_dir

prepare_volcano_data <- function(data, fc_col, p_col) {
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

make_volcano_plot <- function(data, fc_col, p_col, title, x_label, force_label_genes = character()) {
  plot_data <- prepare_volcano_data(data, fc_col, p_col)
  label_data <- plot_data %>%
    filter(!is.na(Genes), Genes %in% force_label_genes) %>%
    distinct(Genes, .keep_all = TRUE)

  max_abs_fc <- max(abs(plot_data[[fc_col]]), na.rm = TRUE)
  x_limit <- max(fc_cutoff * 1.8, ceiling(max_abs_fc)) * 1.08

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
      box.padding = 0.55,
      point.padding = 0.55,
      min.segment.length = 0,
      force = 5,
      force_pull = 0.15,
      max.time = 2,
      max.iter = 20000,
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
    scale_y_continuous(expand = expansion(mult = c(0.04, 0.14))) +
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
      aspect.ratio = 3 / 4,
      plot.margin = margin(8, 12, 8, 8)
    )
}

# cMyc volcano plot
volcano_cMyc <- make_volcano_plot(
  Retron_stats,
  fc_col = "log2FC_cMyc",
  p_col = "P_cMyc",
  title = "cMyc dG4 enriched protein",
  x_label = "log2 Fold Change (cMyc dG4/CTRL)",
  force_label_genes = force_label_genes_cMyc
)
print(volcano_cMyc)
ggsave(file.path(output_dir, "volcano_cMyc_4x3_beautified.png"), volcano_cMyc,
       width = 6.4, height = 5.2, dpi = 600)

# Telo volcano plot
volcano_Telo <- make_volcano_plot(
  Retron_stats,
  fc_col = "log2FC_TELO",
  p_col = "P_TELO",
  title = "Telo dG4 enriched protein",
  x_label = "log2 Fold Change (Telo dG4/CTRL)",
  force_label_genes = force_label_genes_Telo
)
print(volcano_Telo)
ggsave(file.path(output_dir, "volcano_Telo_4x3_beautified.png"), volcano_Telo,
       width = 6.4, height = 5.2, dpi = 600)













