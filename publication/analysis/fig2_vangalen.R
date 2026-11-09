#!/usr/bin/env Rscript
# Figure 2 -- the van Galen 2019 AML bone marrow benchmark of the rank-score
# track (seed 1, antagonist mode, 55 AML drug targets, 19 identities).
#   A  log10 perturbation score, every identity x target pair
#   B  malignant-like minus matched normal identity, per target
#   C  identity-wide offset from a target + identity decomposition
#
# Inputs  : ../data/vangalen2021/<method>/perbscore_all_targets.txt
#           vangalen/<method>/{cell_counts,identity_qc,network_qc,dq_flags,runtime}.tsv
#             (pulled from that run's report by vangalen/extract_report_tables.py;
#              identity_qc / network_qc / dq_flags only exist when the report
#              has a Data quality section)
# Outputs : ../Fig/fig2_vangalen_<method>.{pdf,png}, ../Fig/fig2A_vangalen_heatmap_<method>.{pdf,png} (slide),
#           vangalen_stats_<method>.txt (numbers quoted in the text)
#
# Run from publication/analysis:  Rscript fig2_vangalen.R [hdwgcna|genie3]   (default hdwgcna)

suppressMessages({ library(data.table); library(ggplot2); library(patchwork) })
method <- commandArgs(trailingOnly = TRUE)[1]; if (is.na(method)) method <- "hdwgcna"
other_method <- setdiff(c("hdwgcna", "genie3"), method)
method_lab <- c(hdwgcna = "hdWGCNA", genie3 = "GENIE3")
qc_dir <- file.path("vangalen", method)
out <- file(sprintf("vangalen_stats_%s.txt", method), "w")
say <- function(...) { cat(..., "\n", sep = ""); cat(..., "\n", sep = "", file = out) }

# ---- palette (dataviz reference palette) --------------------------------
col_a <- "#2a78d6"; col_b <- "#eb6834"; ink2 <- "#52514e"; grid <- "#e6e5e1"
seq_ramp <- c("#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b")
theme_np <- theme_minimal(base_size = 7.5) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = grid, linewidth = .3),
        axis.text = element_text(colour = ink2), plot.title = element_text(face = "bold", size = 8),
        plot.title.position = "plot", legend.key.size = unit(8, "pt"))
# "No score": an exact zero, or a value below 1e-6 x the run's median. The
# absolute scale of a score depends on the network method, so the floor is
# relative to each run. GENIE3 leaves a target with no edge at exactly zero
# (296 pairs); hdWGCNA's signed, rescaled TOM leaves it at 1e-15 to 1e-18, and
# every target the report flags as isolated lands there.
no_score <- function(x) x == 0 | x < 1e-6 * median(x[x > 0])
lims <- function(lg) quantile(lg, c(.02, .98), na.rm = TRUE)
read_scores <- function(m) {
  x <- fread(sprintf("../data/vangalen2021/%s/perbscore_all_targets.txt", m))
  x[, `:=`(lg = log10(perb_score), floor = no_score(perb_score))]
  x[, rk := frank(-perb_score, ties.method = "min"), by = target][]
}

say("################ van Galen 2019, AML bone marrow, ", method_lab[method], ", ", paste(unique(fread(sprintf("../data/vangalen2021/%s/perbscore_all_targets.txt", method))$binding), collapse = "/"))
vg <- read_scores(method)
pairs <- data.table(normal = c("HSC", "Prog", "GMP", "ProMono", "Mono", "cDC"),
                    like   = c("HSC-like", "Prog-like", "GMP-like", "ProMono-like", "Mono-like", "cDC-like"))
id_order <- c("B", "Plasma", "T", "CTL", "NK", "earlyEry", "lateEry", as.vector(rbind(pairs$normal, pairs$like)))
cc <- fread(file.path(qc_dir, "cell_counts.tsv"))
say("identities kept: ", uniqueN(vg$cell_type), " (", sum(cc$`Cells used`), " cells); dropped: ",
    paste(cc[Status != "kept", paste0(Identity, " (", `Cells in input`, ")")], collapse = ", "))
say("genes: ", cc[Status == "kept", `gene4use`][1], " in the shared universe (gene4use); per identity ", paste(range(cc[Status == "kept", `Genes in network`]), collapse = "-"), " expressed")
say("scores: ", nrow(vg), " pairs, ", uniqueN(vg$target), " targets (", sum(!grepl(";", unique(vg$target))), " single, ",
    sum(grepl(";", unique(vg$target))), " joint); no score (zero or < 1e-6 x median): ", sum(vg$floor), " (exact zeros: ", sum(vg$perb_score == 0), ")",
    "; non-floor range ", paste(signif(range(vg[floor == FALSE, perb_score]), 2), collapse = "-"),
    "; median ", signif(median(vg$perb_score), 2))
if (sum(vg$floor)) say("identities with most no-score pairs: ", paste(vg[floor == TRUE, .N, by = cell_type][order(-N)][1:min(5, .N), paste0(cell_type, " ", N)], collapse = ", "),
    "; targets with most: ", paste(vg[floor == TRUE, .N, by = target][order(-N)][1:min(5, .N), paste0(target, " ", N)], collapse = ", "),
    "; malignant (-like) identities hold ", vg[floor == TRUE & grepl("-like$", cell_type), .N], " of ", sum(vg$floor))

# Data-quality flags from the report: a target "barely detected" (< 5% of the
# identity's cells), "at noise level" or "isolated" in an identity. A joint
# target is flagged only when every member is.
dq <- if (file.exists(file.path(qc_dir, "dq_flags.tsv"))) fread(file.path(qc_dir, "dq_flags.tsv")) else NULL
flag_of <- function(kind) {
  if (is.null(dq)) return(rep(FALSE, nrow(vg)))
  f <- dq[grepl(paste0(" ", kind, "$"), Check), paste(Identity, sub(paste0(" ", kind, "$"), "", Check))]
  mapply(function(t, ct) all(paste(ct, strsplit(t, ";")[[1]]) %in% f), vg$target, vg$cell_type)
}
vg[, `:=`(barely = flag_of("barely detected"), noise = flag_of("at noise level"), isolated = flag_of("isolated"))]
if (!is.null(dq)) {
  say("report flags: barely detected (<5% of cells) ", sum(vg$barely), " pairs, at noise level ", sum(vg$noise), ", isolated ", sum(vg$isolated),
      "; isolated pairs that are no-score: ", vg[isolated == TRUE, sum(floor)], "/", sum(vg$isolated))
  say("barely-detected pairs that still score: ", vg[barely & !floor, .N], "; ranked first in ", vg[barely & !floor & rk == 1, .N],
      if (vg[barely & !floor & rk == 1, .N]) paste0(" (", paste(vg[barely & !floor & rk == 1, paste(target, cell_type)], collapse = ", "), ")") else "")
  say(sprintf("median log10 score, barely detected vs. not (both scored): %.2f vs %.2f", median(vg[barely & !floor, lg]), median(vg[!barely & !floor, lg])))
  nflag <- dq[Area != "Target", .N, by = Check]; say("identity-level flags: ", paste(sprintf("%s %d", nflag$Check, nflag$N), collapse = ", "))
}

top <- vg[rk == 1, .N, by = cell_type][order(-N)]
say("top identity per target: ", paste(sprintf("%s %d", top$cell_type, top$N), collapse = ", "))
say("targets whose top identity is malignant (-like): ", vg[rk == 1 & grepl("-like$", cell_type), uniqueN(target)], "/", uniqueN(vg$target))
vnf <- vg[floor == FALSE]
vfit <- lm(lg ~ target + cell_type, data = vnf); vce <- coef(vfit)[grep("^cell_type", names(coef(vfit)))]
off <- data.table(identity = c("B", sub("cell_type", "", names(vce))), offset = c(0, unname(vce)))
say(sprintf("variance of log10 score explained: target %.0f%%, identity %.0f%% (sequential SS, target first)",
            100 * anova(vfit)["target", "Sum Sq"] / sum(anova(vfit)[, "Sum Sq"]), 100 * anova(vfit)["cell_type", "Sum Sq"] / sum(anova(vfit)[, "Sum Sq"])))
say("identity main effects vs B (log10), extremes: ", paste(sprintf("%s=%+.2f", off[order(-offset)][c(1:3, .N - 2, .N - 1, .N), identity], off[order(-offset)][c(1:3, .N - 2, .N - 1, .N), offset]), collapse = "; "),
    sprintf("; span %.2f", diff(range(off$offset))))
say(sprintf("mean offset, malignant identities %+.2f vs their normal counterparts %+.2f vs others %+.2f",
            off[identity %in% pairs$like, mean(offset)], off[identity %in% pairs$normal, mean(offset)], off[!identity %in% c(pairs$like, pairs$normal), mean(offset)]))
m <- NULL
if (file.exists(file.path(qc_dir, "network_qc.tsv")) && file.exists(file.path(qc_dir, "identity_qc.tsv"))) {
  nq <- fread(file.path(qc_dir, "network_qc.tsv")); iq <- fread(file.path(qc_dir, "identity_qc.tsv"))
  m <- merge(merge(off, nq[, .(identity, density, isolated_frac, strength_gini, sft_r2)], by = "identity"),
             iq[, .(identity, n_cells, median_counts, sparsity_gene4use)], by = "identity")
  for (v in c("density", "isolated_frac", "strength_gini", "median_counts", "sparsity_gene4use", "n_cells"))
    say(sprintf("Spearman(identity offset, %s) = %+.2f", v, cor(m$offset, m[[v]], method = "spearman")))
  say(sprintf("network density %.3f-%.3f; malignant identities %.4f vs normal counterparts %.4f (mean)", min(nq$density), max(nq$density),
              nq[identity %in% pairs$like, mean(density)], nq[identity %in% pairs$normal, mean(density)]))
} else say("(no identity_qc / network_qc tables for this run: its report predates the Data quality section)")

# malignant-like minus matched normal, per target and pair, both scored
w <- dcast(vg, target ~ cell_type, value.var = "lg"); wf <- dcast(vg, target ~ cell_type, value.var = "floor")
dd <- rbindlist(lapply(seq_len(nrow(pairs)), function(i) data.table(
  target = w$target, pair = pairs$like[i],
  diff = w[[pairs$like[i]]] - w[[pairs$normal[i]]],
  ok = !wf[[pairs$like[i]]] & !wf[[pairs$normal[i]]])))[ok == TRUE]
pd <- dd[, .(median = median(diff), share_like_higher = mean(diff > 0), n = .N,
             p = wilcox.test(diff, mu = 0, exact = FALSE)$p.value), by = pair]
for (i in seq_len(nrow(pd))) say(sprintf("%-13s minus normal: median %+.2f log10, malignant higher for %d/%d targets (Wilcoxon signed-rank p = %.2g)",
                                         pd$pair[i], pd$median[i], round(pd$share_like_higher[i] * pd$n[i]), pd$n[i], pd$p[i]))
say(sprintf("all pairs pooled: median %+.2f, malignant higher for %.0f%% of %d target x pair comparisons", median(dd$diff), 100 * mean(dd$diff > 0), nrow(dd)))
dd[, resid_diff := diff - (off[match(pair, identity), offset] - off[match(sub("-like", "", pair), identity), offset])]
say(sprintf("after removing the two identities' offsets: median %+.2f, malignant higher for %.0f%%", median(dd$resid_diff), 100 * mean(dd$resid_diff > 0)))
blast <- c("FLT3", "KIT", "IL3RA", "CD33")
for (t in c(blast, "BCL2")) say(sprintf("%-6s top-3 identities: %s | malignant minus normal: %s", t, paste(vg[target == t][order(rk)][1:3, cell_type], collapse = ", "),
                                        paste(dd[target == t, sprintf("%s %+.2f", sub("-like", "", pair), diff)], collapse = ", ")))
for (t in grep("FLT3;|;FLT3", unique(vg$target), value = TRUE)) say(sprintf("%-11s top-3: %s", t, paste(vg[target == t][order(rk)][1:3, cell_type], collapse = ", ")))
prog_pairs <- c("HSC-like", "Prog-like", "GMP-like")
say("blast markers: median malignant-minus-normal over the progenitor pairs (HSC, Prog, GMP): ",
    sprintf("%.2f (n = %d)", median(dd[target %in% blast & pair %in% prog_pairs, diff]), dd[target %in% blast & pair %in% prog_pairs, .N]),
    " vs ", sprintf("%.2f (n = %d)", median(dd[!target %in% blast & pair %in% prog_pairs, diff]), dd[!target %in% blast & pair %in% prog_pairs, .N]), " for the other targets",
    sprintf("; Wilcoxon p = %.2g", wilcox.test(dd[target %in% blast & pair %in% prog_pairs, diff], dd[!target %in% blast & pair %in% prog_pairs, diff])$p.value))
tg_prog <- dd[pair %in% prog_pairs, .(med = median(diff), n = .N, higher = sum(diff > 0)), by = target][n >= 2][order(-med)]
say("progenitor pairs, targets scored in >= 2 of 3, highest median malignant excess: ", paste(tg_prog[1:10, sprintf("%s %+.2f (%d/%d)", target, med, higher, n)], collapse = ", "))
say("  ... and lowest: ", paste(tg_prog[.N:(.N - 5), sprintf("%s %+.2f (%d/%d)", target, med, higher, n)], collapse = ", "))
say("  rank of the blast markers among ", nrow(tg_prog), ": ", paste(tg_prog[, .(target, r = .I)][target %in% blast, sprintf("%s %d", target, r)], collapse = ", "))

# Joint targets against their members: does A;B score like A, like B, or
# like neither?
joint <- unique(vg$target)[grepl(";", unique(vg$target))]
score_of <- function(tg) { x <- vg[target == tg][match(id_order, cell_type)]; ifelse(x$floor, NA_real_, x$lg) }
jt <- rbindlist(lapply(joint, function(j) {
  g <- strsplit(j, ";")[[1]]
  a <- score_of(j)
  rho <- sapply(g, function(mm) { b <- score_of(mm); ok <- !is.na(a) & !is.na(b); if (sum(ok) < 5) NA_real_ else cor(a[ok], b[ok], method = "spearman") })
  top_pair <- vg[target == j & rk == 1, cell_type]; tops <- sapply(g, function(mm) vg[target == mm & rk == 1, cell_type])
  data.table(pair = j, member = g, rho = rho, top_matches = tops == top_pair,
             member_median = sapply(g, function(mm) median(vg[target == mm & floor == FALSE, lg])))
}))
jt[, stronger := member_median == max(member_median), by = pair]
say(sprintf("joint targets: %d; Spearman with the stronger member median %.2f, with the weaker %.2f; the pair's top identity is the stronger member's in %d/%d pairs, the weaker's in %d",
            length(joint), median(jt[stronger == TRUE, rho], na.rm = TRUE), median(jt[stronger == FALSE, rho], na.rm = TRUE),
            sum(jt[stronger == TRUE, top_matches]), length(joint), sum(jt[stronger == FALSE, top_matches])))
say(sprintf("joint targets vs single genes: median log10 score %.2f vs %.2f (scored pairs), %.1fx", median(vg[floor == FALSE & grepl(";", target), lg]), median(vg[floor == FALSE & !grepl(";", target), lg]),
            10^(median(vg[floor == FALSE & grepl(";", target), lg]) - median(vg[floor == FALSE & !grepl(";", target), lg]))))
say("joint targets whose top identity is a member's: ", paste(jt[, .(any = any(top_matches)), by = pair][any == TRUE, pair], collapse = ", "))

# ---- the same benchmark under the other network method -------------------
# Both runs share identities and targets; how far do they agree?
if (file.exists(sprintf("../data/vangalen2021/%s/perbscore_all_targets.txt", other_method))) {
  say("\n---- agreement with the ", method_lab[other_method], " run (same identities and targets)")
  ov <- read_scores(other_method)
  mm <- merge(vg[, .(cell_type, target, lg, floor, rk)], ov[, .(cell_type, target, lg_o = lg, floor_o = floor, rk_o = rk)], by = c("cell_type", "target"))
  both <- mm[!floor & !floor_o]
  say(sprintf("pairs scored in both: %d of %d; Spearman of log10 scores %.2f", nrow(both), nrow(mm), cor(both$lg, both$lg_o, method = "spearman")))
  pt <- both[, .(rho = if (.N >= 5) cor(lg, lg_o, method = "spearman") else NA_real_), by = target]
  say(sprintf("per target, Spearman across identities: median %.2f (IQR %.2f to %.2f, %d targets)", median(pt$rho, na.rm = TRUE),
              quantile(pt$rho, .25, na.rm = TRUE), quantile(pt$rho, .75, na.rm = TRUE), sum(!is.na(pt$rho))))
  say("same top identity in both runs: ", mm[rk == 1 & rk_o == 1, .N], "/", uniqueN(mm$target), " targets (", paste(mm[rk == 1 & rk_o == 1, paste(target, cell_type)], collapse = ", "), ")")
  ofit <- lm(lg ~ target + cell_type, data = ov[floor == FALSE]); oce <- coef(ofit)[grep("^cell_type", names(coef(ofit)))]
  oo <- merge(off, data.table(identity = c("B", sub("cell_type", "", names(oce))), offset_o = c(0, unname(oce))), by = "identity")
  say(sprintf("identity offsets: span %.2f here vs %.2f in %s; Spearman between the two runs' offsets %.2f", diff(range(off$offset)), diff(range(oo$offset_o)),
              method_lab[other_method], cor(oo$offset, oo$offset_o, method = "spearman")))
  wo <- dcast(ov, target ~ cell_type, value.var = "lg"); wfo <- dcast(ov, target ~ cell_type, value.var = "floor")
  ddo <- rbindlist(lapply(seq_len(nrow(pairs)), function(i) data.table(
    target = wo$target, pair = pairs$like[i], diff_o = wo[[pairs$like[i]]] - wo[[pairs$normal[i]]],
    ok = !wfo[[pairs$like[i]]] & !wfo[[pairs$normal[i]]])))[ok == TRUE]
  dm <- merge(dd[, .(target, pair, diff)], ddo[, .(target, pair, diff_o)], by = c("target", "pair"))
  say(sprintf("malignant-minus-normal differences scored in both: %d; same sign in %.0f%%; Spearman %.2f", nrow(dm), 100 * mean(sign(dm$diff) == sign(dm$diff_o)),
              cor(dm$diff, dm$diff_o, method = "spearman")))
  pdo <- ddo[, .(median_o = median(diff_o)), by = pair]
  say("per-pair medians, here vs ", method_lab[other_method], ": ", paste(merge(pd, pdo, by = "pair")[match(pairs$like, pair), sprintf("%s %+.2f vs %+.2f", pair, median, median_o)], collapse = "; "))
  say(sprintf("no-score pairs: %d here, %d in %s, %d in both", sum(mm$floor), sum(mm$floor_o), method_lab[other_method], mm[floor & floor_o, .N]))
}

rt_file <- file.path(qc_dir, "runtime.tsv")
if (file.exists(rt_file)) {
  rt <- fread(rt_file); say("\n---- running time (report, Nextflow trace)")
  for (i in seq_len(nrow(rt))) say(sprintf("%-12s %3d tasks, mean %s (min %s, max %s), total %s", rt$Step[i], rt$Tasks[i], rt$Mean[i], rt$Min[i], rt$Max[i], rt$Total[i]))
}

# =========================================================================
# Figure
# =========================================================================
single <- setdiff(unique(vg$target), joint)
hv <- copy(vg)
hv[, `:=`(target_f = factor(target, levels = c(sort(single), sort(joint))), ct_f = factor(cell_type, levels = rev(id_order)), fill = fifelse(floor, NA_real_, lg),
          group = factor(fifelse(grepl(";", target), "Joint targets", "Single-gene targets"), levels = c("Single-gene targets", "Joint targets")))]
ttl <- sprintf("van Galen 2019 AML (%s networks): perturbation score, antagonist mode", method_lab[method])
pA <- ggplot(hv, aes(target_f, ct_f)) +
  geom_tile(aes(fill = fill), colour = "white", linewidth = .3) +
  geom_point(data = hv[barely & !floor], shape = 4, size = .6, colour = "#e34948", stroke = .35) +
  geom_text(data = hv[rk == 1], label = "•", colour = "white", size = 2.2) +
  facet_grid(. ~ group, scales = "free_x", space = "free_x") +
  scale_fill_gradientn(colours = seq_ramp, na.value = "#f0efec", name = expression(log[10]~score), limits = lims(hv$fill), oob = scales::squish) +
  labs(x = NULL, y = NULL, title = paste("A ", ttl)) + theme_np +
  theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 4.8), panel.grid = element_blank(),
        axis.text.y = element_text(size = 6), strip.text = element_text(face = "bold", colour = ink2), panel.spacing.x = unit(4, "pt"))

dd[, pair_f := factor(sub("-like", "", pair), levels = pairs$normal)]
dd[, grp := fifelse(target %in% blast, "FLT3 / KIT / IL3RA / CD33", fifelse(target == "BCL2", "BCL2", "other targets"))]
dd[, grp := factor(grp, levels = c("other targets", "FLT3 / KIT / IL3RA / CD33", "BCL2"))]
pB <- ggplot(dd, aes(pair_f, diff)) +
  geom_hline(yintercept = 0, colour = ink2, linewidth = .3) +
  geom_point(data = dd[grp == "other targets"], position = position_jitter(width = .18, seed = 1), size = .6, alpha = .35, colour = "#9a9892") +
  stat_summary(data = dd, fun = median, geom = "crossbar", width = .5, linewidth = .35, colour = ink2) +
  geom_point(data = dd[grp != "other targets"], aes(colour = grp), size = 1.6) +
  scale_colour_manual(values = c(`FLT3 / KIT / IL3RA / CD33` = col_a, BCL2 = col_b), name = NULL) +
  labs(x = "Malignant-like vs. matched normal identity", y = "log10 (malignant / normal score)", title = "B  Malignant vs. normal counterpart") +
  theme_np + theme(legend.position = "bottom")

off[, `:=`(identity_f = factor(identity, levels = id_order),
           kind = factor(fifelse(grepl("-like$", identity), "malignant", fifelse(identity %in% pairs$normal, "normal counterpart", "other")),
                         levels = c("malignant", "normal counterpart", "other")))]
kind_col <- c(malignant = col_a, `normal counterpart` = col_b, other = "#9a9892")
pC <- ggplot(off, aes(identity_f, offset, colour = kind)) +
  geom_hline(yintercept = 0, colour = ink2, linewidth = .3) +
  geom_segment(aes(xend = identity_f, y = 0, yend = offset), linewidth = .6) +
  geom_point(size = 2) +
  scale_colour_manual(values = kind_col, name = NULL) +
  labs(x = NULL, y = "Identity effect, log10 (vs. B)", title = "C  Identity-wide offset") +
  theme_np + theme(axis.text.x = element_text(angle = 50, hjust = 1, size = 6), legend.position = "bottom")

fig <- pA / ((pB | pC) + plot_layout(widths = c(1, 1.15))) + plot_layout(heights = c(1.15, 1))

ggsave(sprintf("../Fig/fig2_vangalen_%s.pdf", method), fig, width = 7.1, height = 6.0, device = cairo_pdf)
ggsave(sprintf("../Fig/fig2_vangalen_%s.png", method), fig, width = 7.1, height = 6.0, dpi = 300, bg = "white")
ggsave(sprintf("fig2_vangalen_%s_preview.png", method), fig, width = 7.1, height = 6.0, dpi = 150, bg = "white")

# Panel A on its own, sized for a 16:9 slide, with larger type.
pA_slide <- pA +
  geom_point(data = hv[rk == 1], colour = "white", size = 2.6) +
  labs(title = ttl,
       subtitle = paste0("Rows: 19 cell identities, malignant (-like) under their normal counterpart.\n",
                         "White dot: the identity ranked first for that target. Grey: the target has no edge in that identity's network.",
                         if (any(hv$barely)) " Red x: target detected in <5% of the identity's cells." else "")) +
  theme(text = element_text(size = 14), plot.title = element_text(size = 16, face = "bold"), plot.subtitle = element_text(size = 11.5, colour = ink2),
        axis.text.x = element_text(size = 9.5, angle = 55, hjust = 1), axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 13, face = "bold", colour = ink2),
        legend.title = element_text(size = 12), legend.text = element_text(size = 11), legend.key.size = unit(16, "pt"))
pA_slide$layers[[2]]$aes_params$size <- 2.2
ggsave(sprintf("../Fig/fig2A_vangalen_heatmap_%s.png", method), pA_slide, width = 15.3, height = 6.5, dpi = 300, bg = "white")
ggsave(sprintf("../Fig/fig2A_vangalen_heatmap_%s.pdf", method), pA_slide, width = 15.3, height = 6.5, device = cairo_pdf, bg = "white")
close(out)
