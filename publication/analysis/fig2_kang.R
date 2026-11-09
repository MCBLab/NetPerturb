#!/usr/bin/env Rscript
# Figure 2 -- Kang 2018 PBMC (control cells) benchmark, from the NetPerturb run
# `determined_lamarck` (--network hdwgcna --binding agonist --sctknk, seed 1).
#
# Inputs  : ../data/kang2018/*.txt.gz                      (pipeline outputs)
#           target_expression_from_report.tsv              (transcribed from the
#           target_detection_pct_from_report.tsv            report's section 8.2;
#                                                           replace with the run's
#                                                           downsample/target_expression.tsv)
# Outputs : ../Fig/fig2_kang.pdf, kang_stats.txt (numbers quoted in the text)
#
# Run from publication/analysis:  Rscript fig2_kang.R

suppressMessages({
  library(data.table); library(ggplot2); library(patchwork)
})
D <- "../data/kang2018/"
out <- file("kang_stats.txt", "w"); say <- function(...) { cat(..., "\n", sep = ""); cat(..., "\n", sep = "", file = out) }

# ---- palette (dataviz reference palette; slots 1-2 validated all-pairs) ----
col_mono  <- "#2a78d6"   # monocytes
col_other <- "#eb6834"   # other identities
ink2 <- "#52514e"; grid <- "#e6e5e1"
seq_ramp <- c("#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b")

theme_np <- theme_minimal(base_size = 7.5) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = grid, linewidth = .3),
        axis.text = element_text(colour = ink2), plot.title = element_text(face = "bold", size = 8),
        plot.title.position = "plot", legend.key.size = unit(8, "pt"))

tiers <- list(
  cascade = c("IFNAR1", "IFNAR2", "IFNAR1;IFNAR2", "JAK1", "TYK2", "JAK1;TYK2", "STAT1", "STAT2", "STAT1;STAT2", "IRF9"),
  ISG     = c("IRF1", "IRF7", "STAT3", "IFIT1", "IFIT3", "IFI6", "ISG15", "MX1", "OAS1", "CXCL10"),
  lineage = c("CD19", "MS4A1", "IL4R", "TLR4", "CSF1R"))
tier_of <- setNames(rep(names(tiers), lengths(tiers)), unlist(tiers))
tier_lab <- c(cascade = "IFN cascade", ISG = "ISGs", lineage = "Lineage controls")
ct_lab <- c(CD14__Monocytes = "CD14+ Mono", FCGR3A__Monocytes = "FCGR3A+ Mono", B_cells = "B",
            CD4_T_cells = "CD4 T", CD8_T_cells = "CD8 T", NK_cells = "NK")
ct_order <- names(ct_lab)

# ---- perturbation scores + expression -------------------------------------
ps <- fread(paste0(D, "perbscore_all_targets.txt.gz"))
ps[, `:=`(lg = log10(perb_score), floor = perb_score < 1e-10, tier = tier_of[target],
          mono = grepl("Monocytes", cell_type))]
ex <- merge(melt(fread("target_expression_from_report.tsv"), id.vars = "cell_type", variable.name = "gene", value.name = "expr"),
            melt(fread("target_detection_pct_from_report.tsv"), id.vars = "cell_type", variable.name = "gene", value.name = "pct"))
ps[, pct := mapply(function(t, ct) max(ex[cell_type == ct & gene %in% strsplit(t, ";")[[1]], pct]), target, cell_type)]
ps[, expr := mapply(function(t, ct) max(ex[cell_type == ct & gene %in% strsplit(t, ";")[[1]], expr]), target, cell_type)]
ps[, rk := frank(-perb_score, ties.method = "min"), by = target]

say("== Scores: ", nrow(ps), " pairs, ", uniqueN(ps$target), " targets, ", uniqueN(ps$cell_type), " identities")
say("floor scores (<1e-10): ", sum(ps$floor), " -- all with target detection ", paste(unique(ps[floor == TRUE, pct]), collapse = ","), "%")
top1 <- ps[rk == 1, .N, by = .(tier, mono)]
say("targets whose top identity is a monocyte: ", sum(ps[rk == 1 & mono == TRUE, .N]), "/", uniqueN(ps$target))
for (t in names(tiers)) say("  ", t, ": ", ps[rk == 1 & mono & tier == t, .N], "/", length(tiers[[t]]))

nf <- ps[floor == FALSE]
fit <- lm(lg ~ target + cell_type, data = nf)
cte <- coef(fit)[grep("^cell_type", names(coef(fit)))]
say("cell-type main effects vs B cells (log10): ", paste(sprintf("%s=%.2f", sub("cell_type", "", names(cte)), cte), collapse = "; "))
nf[, resid := resid(fit)]
rt <- nf[, .(mono = mean(resid[mono]), other = mean(resid[!mono])), by = tier]
for (i in seq_len(nrow(rt))) say(sprintf("residual %s: monocytes %+.2f, other %+.2f (diff %.2f log10)", rt$tier[i], rt$mono[i], rt$other[i], rt$mono[i] - rt$other[i]))
w <- wilcox.test(nf[tier == "ISG" & mono, resid], nf[tier == "ISG" & !mono, resid])
say(sprintf("ISG residual monocyte vs other: Wilcoxon p = %.2g", w$p.value))
say(sprintf("Spearman(log score, target expression), non-floor pairs: rho = %.2f (n = %d)",
            cor(nf$lg, nf$expr, method = "spearman"), nrow(nf)))
fit2 <- lm(lg ~ target + expr + mono, data = nf)
say(sprintf("lm(log score ~ target + expr + monocyte): monocyte term %.2f log10 (p = %.1g)",
            coef(fit2)["monoTRUE"], summary(fit2)$coefficients["monoTRUE", 4]))
art <- nf[pct < 1]
say("non-floor scores with target detected in <1% of cells: ", nrow(art), " pairs; ranked first in ", art[rk == 1, .N],
    " (", paste(art[rk == 1, paste(target, cell_type)], collapse = ", "), ")")

# Panel A: heatmap
hm <- copy(ps)
hm[, target_f := factor(target, levels = unlist(tiers))]
hm[, ct_f := factor(ct_lab[cell_type], levels = rev(ct_lab[ct_order]))]
hm[, tier_f := factor(tier_lab[tier], levels = tier_lab)]
hm[, fill := fifelse(floor, NA_real_, lg)]
pA <- ggplot(hm, aes(target_f, ct_f)) +
  geom_tile(aes(fill = fill), colour = "white", linewidth = .4) +
  geom_point(data = hm[pct < 1 & !floor], shape = 4, size = 1.1, colour = "#e34948", stroke = .5) +
  geom_text(data = hm[rk == 1], label = "•", colour = "white", size = 2.6, family = "sans") +
  facet_grid(. ~ tier_f, scales = "free_x", space = "free_x") +
  scale_fill_gradientn(colours = seq_ramp, na.value = "#f0efec", name = expression(log[10]~score),
                       limits = c(-6.2, -2.5), oob = scales::squish) +
  labs(x = NULL, y = NULL, title = "A  Perturbation score (hdWGCNA, agonist)") + theme_np +
  theme(axis.text.x = element_text(angle = 50, hjust = 1, size = 6), panel.grid = element_blank(),
        strip.text = element_text(face = "bold", colour = ink2))

# Panel B: residual after removing target and cell-type main effects
nf[, tier_f := factor(tier_lab[tier], levels = tier_lab)]
nf[, grp := factor(ifelse(mono, "Monocytes", "Other identities"), levels = c("Monocytes", "Other identities"))]
pB <- ggplot(nf, aes(tier_f, resid, colour = grp)) +
  geom_hline(yintercept = 0, colour = ink2, linewidth = .3) +
  geom_point(position = position_jitterdodge(jitter.width = .15, dodge.width = .6, seed = 1), size = .9, alpha = .7) +
  stat_summary(fun = mean, geom = "crossbar", position = position_dodge(.6), width = .45, linewidth = .35) +
  scale_colour_manual(values = c(col_mono, col_other), name = NULL) +
  labs(x = NULL, y = "Residual log10 score", title = "B  Target-specific signal") + theme_np +
  theme(legend.position = "bottom")

# ---- Panel C: ISG partners among the strongest network edges -------------
isg <- c("IFIT1","IFIT2","IFIT3","IFIT5","IFI6","IFI44","IFI44L","IFI35","IFI16","IFITM1","IFITM2","IFITM3","ISG15","ISG20",
         "MX1","MX2","OAS1","OAS2","OAS3","OASL","RSAD2","CXCL10","IRF7","IRF1","IRF9","STAT1","STAT2","HERC5","HERC6","XAF1",
         "EPSTI1","PLSCR1","SAMD9L","SAMD9","BST2","LY6E","DDX58","IFIH1","EIF2AK2","USP18","TRIM22","SP110","PARP9","DTX3L",
         "LGALS3BP","GBP1","GBP2","GBP4","GBP5","CMPK2","APOBEC3A","TNFSF10","DDX60","HELZ2","PARP14","RTP4","SIGLEC1","CCL2","CCL8")
tc <- fread(paste0(D, "top_connections_all_targets.txt.gz"))
tc[, `:=`(isg = partner %in% isg, tier = tier_of[target], mono = grepl("Monocytes", cell_type))]
tcs <- tc[, .(frac = mean(isg), n = .N), by = .(tier, mono)]
for (i in seq_len(nrow(tcs))) say(sprintf("top-15 partners that are ISGs: %s, %s: %.1f%% (n = %d edges)", tcs$tier[i], ifelse(tcs$mono[i], "monocytes", "other"), 100 * tcs$frac[i], tcs$n[i]))
tcs[, tier_f := factor(tier_lab[tier], levels = tier_lab)]
tcs[, grp := factor(ifelse(mono, "Monocytes", "Other identities"), levels = c("Monocytes", "Other identities"))]
pC <- ggplot(tcs, aes(tier_f, frac, fill = grp)) +
  geom_col(position = position_dodge(.75), width = .7) +
  geom_text(aes(label = sprintf("%.0f%%", 100 * frac)), position = position_dodge(.75), vjust = -.4, size = 2.2, colour = ink2) +
  scale_fill_manual(values = c(col_mono, col_other), name = NULL, guide = "none") +
  scale_y_continuous(labels = scales::percent, expand = expansion(c(0, .15))) +
  labs(x = NULL, y = "ISGs among top-15 partners", title = "C  Network neighbourhood") + theme_np

# ---- Panel D: knockout target-dependence ---------------------------------
ko <- fread(paste0(D, "sctenifoldknk_all_targets.txt.gz"))
sig <- ko[p.adj < 0.05]
ribo <- "^RP[SL][0-9]|^RPLP|^RPSA"
say(sprintf("KO: %d combinations, %d significant rows (FDR<0.05); ribosomal genes %.0f%% of hits vs %.1f%% of tested genes",
            uniqueN(ko[, .(cell_type, target)]), nrow(sig), 100 * mean(grepl(ribo, sig$gene)), 100 * mean(grepl(ribo, ko$gene))))
jac <- sig[, {
  s <- split(gene, target); p <- combn(names(s), 2)
  .(jaccard = apply(p, 2, function(x) length(intersect(s[[x[1]]], s[[x[2]]])) / length(union(s[[x[1]]], s[[x[2]]]))),
    n_sig_median = as.numeric(median(lengths(s))))
}, by = cell_type]
jsum <- jac[, .(median_j = median(jaccard), n_sig = n_sig_median[1]), by = cell_type]
rho <- ko[, {
  m <- as.matrix(dcast(.SD, gene ~ target, value.var = "FC")[, -1]); cm <- suppressWarnings(cor(m, method = "spearman", use = "pairwise"))
  .(rho = median(cm[upper.tri(cm)], na.rm = TRUE))
}, by = cell_type]
jsum <- merge(jsum, rho)
for (i in seq_len(nrow(jsum))) say(sprintf("KO %s: median %g hits/target, median pairwise Jaccard %.2f, median FC Spearman %.2f",
                                           jsum$cell_type[i], jsum$n_sig[i], jsum$median_j[i], jsum$rho[i]))
say("CD14+ monocyte hits, all targets: ", paste(sig[cell_type == "CD14__Monocytes", unique(gene)], collapse = ", "))
jac[, ct_f := factor(ct_lab[cell_type], levels = ct_lab[ct_order])]
jac[, grp := factor(ifelse(grepl("Mono", cell_type), "Monocytes", "Other identities"), levels = c("Monocytes", "Other identities"))]
pD <- ggplot(jac, aes(ct_f, jaccard, colour = grp)) +
  geom_boxplot(outlier.shape = NA, width = .55, linewidth = .35) +
  geom_point(position = position_jitter(width = .15, seed = 1), size = .35, alpha = .35) +
  geom_text(data = jsum[, .(ct_f = factor(ct_lab[cell_type], levels = ct_lab[ct_order]), n_sig)],
            aes(ct_f, 1.09, label = paste0("n=", n_sig)), inherit.aes = FALSE, size = 2.1, colour = ink2) +
  scale_colour_manual(values = c(col_mono, col_other), guide = "none") +
  scale_y_continuous(limits = c(0, 1.12), breaks = seq(0, 1, .25)) +
  labs(x = NULL, y = "Jaccard, DR genes between targets", title = "D  Knockout: overlap across targets") + theme_np +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

# ---- Panel E: running time per step (report section 9.1, Nextflow trace) --
rt <- data.table(
  step  = c("DOWNSAMPLE", "HDWGCNA", "RANK_SCORE", "SCTENIFOLDKNK_BUILD", "SCTENIFOLDKNK_KO", "GSEA_SCTENIFOLDKNK"),
  tasks = c(1, 6, 25, 6, 150, 1),
  mean_s = c(33.0, 61, 31.3, 43 * 60 + 14, 5.5, 5 * 60 + 39),
  total_s = c(33.0, 6 * 60 + 9, 13 * 60 + 2, 4 * 3600 + 19 * 60, 13 * 60 + 50, 5 * 60 + 39))
naive_h <- 150 * rt[step == "SCTENIFOLDKNK_BUILD", mean_s] / 3600
split_h <- (rt[step == "SCTENIFOLDKNK_BUILD", total_s] + rt[step == "SCTENIFOLDKNK_KO", total_s]) / 3600
say(sprintf("compute: %.1f h total over all tasks; wall-clock 8m10s (resumed run; see report)", sum(rt$total_s) / 3600))
say(sprintf("knockout: build-once = %.1f CPU-task h; one scTenifoldKnk() call per pair (150 builds) ~ %.0f h -> %.0fx less",
            split_h, naive_h, naive_h / split_h))
rt[, step_f := factor(step, levels = rev(step))]
pE <- ggplot(rt, aes(y = step_f)) +
  geom_segment(aes(x = 1, xend = mean_s, yend = step_f), colour = seq_ramp[4], linewidth = 3.2) +
  geom_text(aes(x = mean_s, label = paste0("\u00d7", tasks, ifelse(tasks == 1, " task", " tasks"))), hjust = -.15, size = 2.2, colour = ink2) +
  scale_x_log10(breaks = c(1, 10, 60, 600, 3600), labels = c("1 s", "10 s", "1 min", "10 min", "1 h"),
                limits = c(1, 3600 * 8), expand = expansion(c(0, 0))) +
  labs(x = "Mean time per task (log scale)", y = NULL, title = "E  Running time per step") + theme_np +
  theme(panel.grid.major.y = element_blank())

fig <- pA / (pB | pC) / (pD | pE) + plot_layout(heights = c(1.15, 1, 1))
ggsave("../Fig/fig2_kang.pdf", fig, width = 7.1, height = 7.4, device = cairo_pdf)
ggsave("fig2_kang_preview.png", fig, width = 7.1, height = 7.4, dpi = 150)
close(out)
