#!/usr/bin/env Rscript
# Figure 2 -- the two rank-score benchmarks (seed 1):
#   A-C  Kang 2018 control PBMCs, hdWGCNA networks, agonist mode, 25 IFN-axis / ISG / lineage targets
#   D-F  van Galen 2019 AML bone marrow, GENIE3 networks, antagonist mode, 55 AML drug targets
#
# Inputs  : ../data/kang2018/{perbscore,top_connections}_all_targets.txt.gz
#           target_expression_from_report.tsv, target_detection_pct_from_report.tsv
#             (Kang target expression, transcribed from that report's section 8.2)
#           ../data/vangalen2021/genie3/perbscore_all_targets.txt
#           vangalen/genie3/cell_counts.tsv (from that report; identity_qc / network_qc
#             too when the run's report carries a Data quality section)
# Outputs : ../Fig/fig2_benchmarks.pdf, benchmark_stats.txt (numbers quoted in the text)
#
# Run from publication/analysis:  Rscript fig2_benchmarks.R

suppressMessages({ library(data.table); library(ggplot2); library(patchwork) })
out <- file("benchmark_stats.txt", "w")
say <- function(...) { cat(..., "\n", sep = ""); cat(..., "\n", sep = "", file = out) }

# ---- palette (dataviz reference palette) --------------------------------
col_a <- "#2a78d6"; col_b <- "#eb6834"; col_c <- "#1baf7a"; ink2 <- "#52514e"; grid <- "#e6e5e1"
seq_ramp <- c("#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b")
theme_np <- theme_minimal(base_size = 7.5) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = grid, linewidth = .3),
        axis.text = element_text(colour = ink2), plot.title = element_text(face = "bold", size = 8),
        plot.title.position = "plot", legend.key.size = unit(8, "pt"))
# "No score": an exact zero, or a value below 1e-6 x the run's median. The
# absolute scale of a score depends on the network method (GENIE3 importances
# sit three orders of magnitude below hdWGCNA's rank-normalised weights), so
# the floor is relative to each run: on Kang it catches the nine
# floating-point residues at 1e-15, on van Galen the 296 exact zeros, and
# nothing else in either.
no_score <- function(x) x == 0 | x < 1e-6 * median(x[x > 0])
lims <- function(lg) quantile(lg, c(.02, .98), na.rm = TRUE)

# =========================================================================
# Kang 2018 (A-C)
# =========================================================================
say("################ Kang 2018, control PBMCs, agonist")
tiers <- list(
  cascade = c("IFNAR1", "IFNAR2", "IFNAR1;IFNAR2", "JAK1", "TYK2", "JAK1;TYK2", "STAT1", "STAT2", "STAT1;STAT2", "IRF9"),
  ISG     = c("IRF1", "IRF7", "STAT3", "IFIT1", "IFIT3", "IFI6", "ISG15", "MX1", "OAS1", "CXCL10"),
  lineage = c("CD19", "MS4A1", "IL4R", "TLR4", "CSF1R"))
tier_of <- setNames(rep(names(tiers), lengths(tiers)), unlist(tiers))
tier_lab <- c(cascade = "IFN cascade", ISG = "ISGs", lineage = "Controls")
ct_lab <- c(CD14__Monocytes = "CD14+ Mono", FCGR3A__Monocytes = "FCGR3A+ Mono", B_cells = "B",
            CD4_T_cells = "CD4 T", CD8_T_cells = "CD8 T", NK_cells = "NK")

ps <- fread("../data/kang2018/perbscore_all_targets.txt.gz")
ps[, `:=`(lg = log10(perb_score), floor = no_score(perb_score), tier = tier_of[target], mono = grepl("Monocytes", cell_type))]
ex <- merge(melt(fread("target_expression_from_report.tsv"), id.vars = "cell_type", variable.name = "gene", value.name = "expr"),
            melt(fread("target_detection_pct_from_report.tsv"), id.vars = "cell_type", variable.name = "gene", value.name = "pct"))
ps[, pct  := mapply(function(t, ct) max(ex[cell_type == ct & gene %in% strsplit(t, ";")[[1]], pct]),  target, cell_type)]
ps[, expr := mapply(function(t, ct) max(ex[cell_type == ct & gene %in% strsplit(t, ";")[[1]], expr]), target, cell_type)]
ps[, rk := frank(-perb_score, ties.method = "min"), by = target]

say("scores: ", nrow(ps), " pairs, ", uniqueN(ps$target), " targets, ", uniqueN(ps$cell_type), " identities; no score (zero or < 1e-6 x median): ", sum(ps$floor),
    "; non-floor range ", paste(signif(range(ps[floor == FALSE, perb_score]), 2), collapse = "-"))
say("targets whose top identity is a monocyte: ", ps[rk == 1 & mono, .N], "/", uniqueN(ps$target),
    " (cascade ", ps[rk == 1 & mono & tier == "cascade", .N], "/10, ISG ", ps[rk == 1 & mono & tier == "ISG", .N], "/10, lineage ", ps[rk == 1 & mono & tier == "lineage", .N], "/5)")
nf <- ps[floor == FALSE]
fit <- lm(lg ~ target + cell_type, data = nf); cte <- coef(fit)[grep("^cell_type", names(coef(fit)))]
say("identity main effects vs B (log10): ", paste(sprintf("%s=%.2f", sub("cell_type", "", names(cte)), cte), collapse = "; "))
nf[, resid := resid(fit)]
rt <- nf[, .(mono = mean(resid[mono]), other = mean(resid[!mono])), by = tier]
for (i in seq_len(nrow(rt))) say(sprintf("residual %s: monocytes %+.2f, other %+.2f (diff %+.2f log10)", rt$tier[i], rt$mono[i], rt$other[i], rt$mono[i] - rt$other[i]))
say(sprintf("ISG residual monocyte vs other: Wilcoxon p = %.2g", wilcox.test(nf[tier == "ISG" & mono, resid], nf[tier == "ISG" & !mono, resid])$p.value))
say(sprintf("Spearman(log score, target expression), non-floor: rho = %.2f (n = %d)", cor(nf$lg, nf$expr, method = "spearman"), nrow(nf)))
fit2 <- lm(lg ~ target + expr + mono, data = nf)
say(sprintf("lm(log score ~ target + expr + monocyte): monocyte term %.2f log10 (p = %.1g)", coef(fit2)["monoTRUE"], summary(fit2)$coefficients["monoTRUE", 4]))
art <- nf[pct < 1]
say("non-floor scores with the target detected in <1% of cells: ", nrow(art), " pairs; ranked first in ", art[rk == 1, .N], " (", paste(art[rk == 1, paste(target, cell_type)], collapse = ", "), ")")

hm <- copy(ps)
hm[, `:=`(target_f = factor(target, levels = unlist(tiers)), ct_f = factor(ct_lab[cell_type], levels = rev(ct_lab)),
          tier_f = factor(tier_lab[tier], levels = tier_lab), fill = fifelse(floor, NA_real_, lg))]
pA <- ggplot(hm, aes(target_f, ct_f)) +
  geom_tile(aes(fill = fill), colour = "white", linewidth = .4) +
  geom_point(data = hm[pct < 1 & !floor], shape = 4, size = 1.1, colour = "#e34948", stroke = .5) +
  geom_text(data = hm[rk == 1], label = "•", colour = "white", size = 2.6) +
  facet_grid(. ~ tier_f, scales = "free_x", space = "free_x") +
  scale_fill_gradientn(colours = seq_ramp, na.value = "#f0efec", name = expression(log[10]~score), limits = lims(hm$fill), oob = scales::squish) +
  labs(x = NULL, y = NULL, title = "A  Kang 2018 PBMCs (hdWGCNA): perturbation score (agonist)") + theme_np +
  theme(axis.text.x = element_text(angle = 50, hjust = 1, size = 6), panel.grid = element_blank(), strip.text = element_text(face = "bold", colour = ink2))

nf[, `:=`(tier_f = factor(tier_lab[tier], levels = tier_lab), grp = factor(ifelse(mono, "Monocytes", "Other identities"), levels = c("Monocytes", "Other identities")))]
pB <- ggplot(nf, aes(tier_f, resid, colour = grp)) +
  geom_hline(yintercept = 0, colour = ink2, linewidth = .3) +
  geom_point(position = position_jitterdodge(jitter.width = .15, dodge.width = .6, seed = 1), size = .9, alpha = .7) +
  stat_summary(fun = mean, geom = "crossbar", position = position_dodge(.6), width = .45, linewidth = .35) +
  scale_colour_manual(values = c(col_a, col_b), name = NULL) +
  labs(x = NULL, y = "Residual log10 score", title = "B  Target-specific signal") + theme_np + theme(legend.position = "bottom")

isg <- c("IFIT1","IFIT2","IFIT3","IFIT5","IFI6","IFI44","IFI44L","IFI35","IFI16","IFITM1","IFITM2","IFITM3","ISG15","ISG20","MX1","MX2","OAS1","OAS2","OAS3","OASL","RSAD2","CXCL10","IRF7","IRF1","IRF9","STAT1","STAT2","HERC5","HERC6","XAF1","EPSTI1","PLSCR1","SAMD9L","SAMD9","BST2","LY6E","DDX58","IFIH1","EIF2AK2","USP18","TRIM22","SP110","PARP9","DTX3L","LGALS3BP","GBP1","GBP2","GBP4","GBP5","CMPK2","APOBEC3A","TNFSF10","DDX60","HELZ2","PARP14","RTP4","SIGLEC1","CCL2","CCL8")
tc <- fread("../data/kang2018/top_connections_all_targets.txt.gz")
tc[, `:=`(isg = partner %in% isg, tier = tier_of[target], mono = grepl("Monocytes", cell_type))]
tcs <- tc[, .(frac = mean(isg), n = .N), by = .(tier, mono)]
for (i in seq_len(nrow(tcs))) say(sprintf("top-15 partners that are ISGs, %s, %s: %.1f%% (n = %d edges)", tcs$tier[i], ifelse(tcs$mono[i], "monocytes", "other"), 100 * tcs$frac[i], tcs$n[i]))
tcs[, `:=`(tier_f = factor(tier_lab[tier], levels = tier_lab), grp = factor(ifelse(mono, "Monocytes", "Other identities"), levels = c("Monocytes", "Other identities")))]
pC <- ggplot(tcs, aes(tier_f, frac, fill = grp)) +
  geom_col(position = position_dodge(.75), width = .7) +
  geom_text(aes(label = sprintf("%.0f%%", 100 * frac)), position = position_dodge(.75), vjust = -.4, size = 2.2, colour = ink2) +
  scale_fill_manual(values = c(col_a, col_b), guide = "none") +
  scale_y_continuous(labels = scales::percent, expand = expansion(c(0, .15))) +
  labs(x = NULL, y = "ISGs among top partners", title = "C  Network neighbourhood") + theme_np

# =========================================================================
# van Galen 2019 (D-F)
# =========================================================================
say("\n################ van Galen 2019, AML bone marrow, GENIE3, antagonist")
vg <- fread("../data/vangalen2021/genie3/perbscore_all_targets.txt")
vg[, `:=`(lg = log10(perb_score), floor = no_score(perb_score))]
vg[, rk := frank(-perb_score, ties.method = "min"), by = target]
pairs <- data.table(normal = c("HSC", "Prog", "GMP", "ProMono", "Mono", "cDC"),
                    like   = c("HSC-like", "Prog-like", "GMP-like", "ProMono-like", "Mono-like", "cDC-like"))
id_order <- c("B", "Plasma", "T", "CTL", "NK", "earlyEry", "lateEry", as.vector(rbind(pairs$normal, pairs$like)))
cc <- fread("vangalen/genie3/cell_counts.tsv")
say("identities kept: ", uniqueN(vg$cell_type), " (", sum(cc$`Cells used`), " cells); dropped: ",
    paste(cc[Status != "kept", paste0(Identity, " (", `Cells in input`, ")")], collapse = ", "))
say("scores: ", nrow(vg), " pairs, ", uniqueN(vg$target), " targets (", sum(!grepl(";", unique(vg$target))), " single, ",
    sum(grepl(";", unique(vg$target))), " joint); no score (zero or < 1e-6 x median): ", sum(vg$floor), " (exact zeros: ", sum(vg$perb_score == 0), ")",
    "; non-floor range ", paste(signif(range(vg[floor == FALSE, perb_score]), 2), collapse = "-"))
say("identities with most zero scores: ", paste(vg[floor == TRUE, .N, by = cell_type][order(-N)][1:4, paste0(cell_type, " ", N)], collapse = ", "),
    "; targets with most: ", paste(vg[floor == TRUE, .N, by = target][order(-N)][1:4, paste0(target, " ", N)], collapse = ", "))
top <- vg[rk == 1, .N, by = cell_type][order(-N)]
say("top identity per target: ", paste(sprintf("%s %d", top$cell_type, top$N), collapse = ", "))
vnf <- vg[floor == FALSE]
vfit <- lm(lg ~ target + cell_type, data = vnf); vce <- coef(vfit)[grep("^cell_type", names(coef(vfit)))]
off <- data.table(identity = c("B", sub("cell_type", "", names(vce))), offset = c(0, unname(vce)))
say("identity main effects vs B (log10), extremes: ", paste(sprintf("%s=%+.2f", off[order(-offset)][c(1:3, .N - 1, .N), identity], off[order(-offset)][c(1:3, .N - 1, .N), offset]), collapse = "; "))
if (file.exists("vangalen/network_qc.tsv") && file.exists("vangalen/identity_qc.tsv")) {
  nq <- fread("vangalen/network_qc.tsv"); iq <- fread("vangalen/identity_qc.tsv")
  m <- merge(merge(off, nq[, .(identity, density)], by = "identity"), iq[, .(identity, median_counts, sparsity_gene4use)], by = "identity")
  say(sprintf("Spearman(identity offset, network density) = %.2f; (offset, median depth) = %.2f; (offset, sparsity) = %.2f",
              cor(m$offset, m$density, method = "spearman"), cor(m$offset, m$median_counts, method = "spearman"), cor(m$offset, m$sparsity_gene4use, method = "spearman")))
} else say("(no identity_qc / network_qc tables for this run: its report predates the Data quality section)")

# malignant-like minus matched normal, per target and pair, both non-floor
w <- dcast(vg, target ~ cell_type, value.var = "lg"); wf <- dcast(vg, target ~ cell_type, value.var = "floor")
dd <- rbindlist(lapply(seq_len(nrow(pairs)), function(i) data.table(
  target = w$target, pair = pairs$like[i],
  diff = w[[pairs$like[i]]] - w[[pairs$normal[i]]],
  ok = !wf[[pairs$like[i]]] & !wf[[pairs$normal[i]]])))[ok == TRUE]
pd <- dd[, .(median = median(diff), share_like_higher = mean(diff > 0), n = .N), by = pair]
for (i in seq_len(nrow(pd))) say(sprintf("%-13s minus normal: median %+.2f log10, malignant higher for %d/%d targets", pd$pair[i], pd$median[i], round(pd$share_like_higher[i] * pd$n[i]), pd$n[i]))
say(sprintf("all pairs pooled: median %+.2f, malignant higher for %.0f%% of %d target x pair comparisons", median(dd$diff), 100 * mean(dd$diff > 0), nrow(dd)))
blast <- c("FLT3", "KIT",  "IL3RA", "CD33")
for (t in c(blast, "BCL2")) say(sprintf("%-6s top-3 identities: %s | malignant minus normal: %s", t, paste(vg[target == t][order(rk)][1:3, cell_type], collapse = ", "),
                                        paste(dd[target == t, sprintf("%s %+.2f", sub("-like", "", pair), diff)], collapse = ", ")))
for (t in grep("FLT3;|;FLT3", unique(vg$target), value = TRUE)) say(sprintf("%-11s top-3: %s", t, paste(vg[target == t][order(rk)][1:3, cell_type], collapse = ", ")))
say("blast markers: median malignant-minus-normal over the progenitor pairs (HSC, Prog, GMP): ",
    sprintf("%.2f", median(dd[target %in% blast & pair %in% c("HSC-like", "Prog-like", "GMP-like"), diff])),
    " vs ", sprintf("%.2f", median(dd[!target %in% blast & pair %in% c("HSC-like", "Prog-like", "GMP-like"), diff])), " for the other targets")

# Joint targets against their members: does A;B score like A, like B, or
# like neither? Spearman over the identities with a score for both, and which
# member's top identity the pair's top identity matches.
joint <- unique(vg$target)[grepl(";", unique(vg$target))]
score_of <- function(tg) { x <- vg[target == tg][match(id_order, cell_type)]; ifelse(x$floor, NA_real_, x$lg) }
jt <- rbindlist(lapply(joint, function(j) {
  g <- strsplit(j, ";")[[1]]
  a <- score_of(j)
  rho <- sapply(g, function(m) { b <- score_of(m); ok <- !is.na(a) & !is.na(b); if (sum(ok) < 5) NA_real_ else cor(a[ok], b[ok], method = "spearman") })
  top_pair <- vg[target == j & rk == 1, cell_type]; tops <- sapply(g, function(m) vg[target == m & rk == 1, cell_type])
  data.table(pair = j, member = g, rho = rho, top_matches = tops == top_pair,
             member_median = sapply(g, function(m) median(vg[target == m & floor == FALSE, lg])))
}))
jt[, stronger := member_median == max(member_median), by = pair]
say(sprintf("joint targets: %d; Spearman with the stronger member median %.2f, with the weaker %.2f; the pair's top identity is the stronger member's in %d/%d pairs, the weaker's in %d",
            length(joint), median(jt[stronger == TRUE, rho], na.rm = TRUE), median(jt[stronger == FALSE, rho], na.rm = TRUE),
            sum(jt[stronger == TRUE, top_matches]), length(joint), sum(jt[stronger == FALSE, top_matches])))
say(sprintf("joint targets score higher than single genes: median log10 score %.2f vs %.2f (non-zero pairs), %.1fx", median(vg[floor == FALSE & grepl(";", target), lg]), median(vg[floor == FALSE & !grepl(";", target), lg]),
            10^(median(vg[floor == FALSE & grepl(";", target), lg]) - median(vg[floor == FALSE & !grepl(";", target), lg]))))
say("pairs whose top identity is neither member's: ", paste(jt[, .(none = !any(top_matches)), by = pair][none == TRUE, pair], collapse = ", "))

# D: heatmap of every target, single genes then joint targets, pairs adjacent
single <- setdiff(unique(vg$target), joint)
hv <- copy(vg)
hv[, `:=`(target_f = factor(target, levels = c(sort(single), sort(joint))), ct_f = factor(cell_type, levels = rev(id_order)), fill = fifelse(floor, NA_real_, lg),
          group = factor(fifelse(grepl(";", target), "Joint targets", "Single-gene targets"), levels = c("Single-gene targets", "Joint targets")))]
pD <- ggplot(hv, aes(target_f, ct_f)) +
  geom_tile(aes(fill = fill), colour = "white", linewidth = .3) +
  geom_text(data = hv[rk == 1], label = "\u2022", colour = "white", size = 2.2) +
  facet_grid(. ~ group, scales = "free_x", space = "free_x") +
  scale_fill_gradientn(colours = seq_ramp, na.value = "#f0efec", name = expression(log[10]~score), limits = lims(hv$fill), oob = scales::squish) +
  labs(x = NULL, y = NULL, title = "D  van Galen 2019 AML (GENIE3): perturbation score (antagonist)") + theme_np +
  theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 4.8), panel.grid = element_blank(),
        axis.text.y = element_text(size = 6), strip.text = element_text(face = "bold", colour = ink2), panel.spacing.x = unit(4, "pt"))

# E: malignant-like minus matched normal
dd[, pair_f := factor(sub("-like", "", pair), levels = pairs$normal)]
dd[, grp := fifelse(target %in% blast, "FLT3 / KIT / IL3RA / CD33", fifelse(target == "BCL2", "BCL2", "other targets"))]
dd[, grp := factor(grp, levels = c("other targets", "FLT3 / KIT / IL3RA / CD33", "BCL2"))]
pE <- ggplot(dd, aes(pair_f, diff)) +
  geom_hline(yintercept = 0, colour = ink2, linewidth = .3) +
  geom_point(data = dd[grp == "other targets"], position = position_jitter(width = .18, seed = 1), size = .6, alpha = .35, colour = "#9a9892") +
  stat_summary(data = dd, fun = median, geom = "crossbar", width = .5, linewidth = .35, colour = ink2) +
  geom_point(data = dd[grp != "other targets"], aes(colour = grp), size = 1.6) +
  scale_colour_manual(values = c(`FLT3 / KIT / IL3RA / CD33` = col_a, BCL2 = col_b), name = NULL) +
  labs(x = "Malignant-like vs. matched normal identity", y = "Malignant minus normal, log10 score", title = "E  Malignant vs. normal counterpart") +
  theme_np + theme(legend.position = "bottom")

# F: identity offsets
off[, `:=`(identity_f = factor(identity, levels = id_order),
           kind = fifelse(grepl("-like$", identity), "malignant", fifelse(identity %in% pairs$normal, "normal counterpart", "other")))]
pF <- ggplot(off, aes(identity_f, offset, colour = kind)) +
  geom_hline(yintercept = 0, colour = ink2, linewidth = .3) +
  geom_segment(aes(xend = identity_f, y = 0, yend = offset), linewidth = .6) +
  geom_point(size = 2) +
  scale_colour_manual(values = c(malignant = col_a, `normal counterpart` = col_b, other = "#9a9892"), name = NULL) +
  labs(x = NULL, y = "Identity effect, log10 (vs. B)", title = "F  Identity-wide offset") +
  theme_np + theme(axis.text.x = element_text(angle = 50, hjust = 1, size = 6), legend.position = "bottom")

fig <- (pA | (pB / pC)) + plot_layout(widths = c(1.6, 1))
fig <- fig / pD / (pE | pF) + plot_layout(heights = c(1.1, 1.15, 0.95))
# Panel D on its own, sized for a 16:9 slide, with larger type.
pD_slide <- pD +
  geom_point(data = hv[rk == 1], colour = "white", size = 2.6) +
  labs(title = "van Galen 2019 AML (GENIE3 networks): perturbation score, antagonist mode",
       subtitle = paste0("Rows: 19 cell identities, malignant (-like) under their normal counterpart.\n",
                         "White dot: the identity ranked first for that target. Grey: the target has no edge in that identity's network.")) +
  theme(text = element_text(size = 14), plot.title = element_text(size = 16, face = "bold"), plot.subtitle = element_text(size = 11.5, colour = ink2),
        axis.text.x = element_text(size = 9.5, angle = 55, hjust = 1), axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 13, face = "bold", colour = ink2),
        legend.title = element_text(size = 12), legend.text = element_text(size = 11), legend.key.size = unit(16, "pt"))
ggsave("../Fig/fig2D_vangalen_heatmap.png", pD_slide, width = 15.3, height = 6.5, dpi = 300, bg = "white")
ggsave("../Fig/fig2D_vangalen_heatmap.pdf", pD_slide, width = 15.3, height = 6.5, device = cairo_pdf, bg = "white")

ggsave("../Fig/fig2_benchmarks.pdf", fig, width = 7.1, height = 7.6, device = cairo_pdf)
ggsave("fig2_benchmarks_preview.png", fig, width = 7.1, height = 7.6, dpi = 150)
close(out)
