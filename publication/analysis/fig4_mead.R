#!/usr/bin/env Rscript
# Figure 4 -- intestinal organoid stem cells (Mead et al. 2022), GENIE3
# networks, antagonist mode, 206 drug targets over 9 identities. Does the
# in-silico score single out the cycling stem II/III pool for Xpo1, as the
# paper's KPT-330 (selinexor) phenotype -- Paneth enrichment through stem
# II/III-specific cell-cycle and stress changes -- would suggest?
#   A  Xpo1 across the differentiation axis: its percentile among the 206
#      targets within each identity (free of the identity offset), raw score
#      printed
#   B  every target's residual score in stem III (x) against stem II (y):
#      where Xpo1 sits next to the paper's tested targets and the cell-cycle
#      targets
# Run from publication/analysis:  Rscript fig4_mead.R
# Writes ../Fig/fig4_mead.{pdf,png} and mead_stats.txt

suppressMessages({ library(data.table); library(ggplot2); library(patchwork); library(ggrepel) })
out <- file("mead_stats_tmp.txt", "w"); say <- function(...) { cat(..., "\n", sep = ""); cat(..., "\n", sep = "", file = out) }

col_a <- "#2a78d6"; col_b <- "#eb6834"; col_c <- "#1baf7a"; col_d <- "#eda100"; ink <- "#0b0b0b"; ink2 <- "#52514e"; grid <- "#e6e5e1"
seq_ramp <- c("#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b")
theme_np <- theme_minimal(base_size = 7.5) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = grid, linewidth = .3),
        axis.text = element_text(colour = ink2), plot.title = element_text(face = "bold", size = 8),
        plot.title.position = "plot", legend.key.size = unit(8, "pt"))

# Which network method produced the scores in data/mead2022, read from the
# report that sits beside them, so a rerun after swapping the run's files in
# labels its own outputs instead of overwriting the other method's.
method <- local({
  h <- paste(readLines("../data/mead2022/netperturb_report.html", warn = FALSE, n = 60000), collapse = " ")
  m <- regmatches(h, regexpr("<code>--network</code>: [a-z0-9]+", h))
  if (length(m)) sub(".*: ", "", m[1]) else "unknown"
})
stats_file <- paste0(sub("_tmp$", "", tools::file_path_sans_ext(summary(out)$description)), "_", method, ".txt")
say("network method: ", method)

ids <- c("stem_I", "stem_II", "stem_III", "early_secretory", "early_Paneth", "Paneth", "enteroendocrine", "early_enterocyte", "enterocyte")
id_lab <- c(stem_I = "stem I", stem_II = "stem II", stem_III = "stem III", early_secretory = "early secretory", early_Paneth = "early Paneth",
            Paneth = "Paneth", enteroendocrine = "enteroendocrine", early_enterocyte = "early enterocyte", enterocyte = "enterocyte")
stem23 <- c("stem_II", "stem_III")
paper <- c("Xpo1", "Aurka", "Aurkb", "Cdk4", "Cdk6", "Map2k1", "Map2k2", "Trp53")
cycle <- c("Cdk1", "Cdk2", "Plk1", "Wee1", "Chek1", "Top2a", "Kif11", "Birc5")

d <- fread("../data/mead2022/perbscore_all_targets.txt")
d[, lg := log10(perb_score)]
d[, floor := perb_score == 0 | perb_score < 1e-6 * median(perb_score[perb_score > 0])]
say("pairs ", nrow(d), "; targets ", uniqueN(d$target), "; identities ", uniqueN(d$cell_type), "; no score ", sum(d$floor))
nf <- d[floor == FALSE]
fit <- lm(lg ~ target + cell_type, data = nf)
nf[, resid := resid(fit)]
nf[, pct := 100 * (frank(lg) - 0.5) / .N, by = cell_type]
ce <- coef(fit)[grep("^cell_type", names(coef(fit)))]
off <- data.table(identity = c(sort(unique(nf$cell_type))[1], sub("cell_type", "", names(ce))), offset = c(0, unname(ce)))
off[, offset := offset - mean(offset)]
say("identity offsets (log10, centred): ", paste(off[order(-offset), sprintf("%s %+.2f", id_lab[identity], offset)], collapse = ", "))

# ---- A: Xpo1 along the differentiation axis ---------------------------------
xa <- nf[target == "Xpo1"][match(ids, cell_type)]
xa[, `:=`(identity = factor(id_lab[ids], levels = id_lab[ids]), group = fifelse(ids %in% stem23, "stem II / III", fifelse(ids == "stem_I", "stem I", "differentiating")))]
say("Xpo1 percentile among the 206 targets, by identity: ", paste(xa[!is.na(pct), sprintf("%s %.0f", id_lab[cell_type], pct)], collapse = ", "),
    "; no score in ", paste(xa[is.na(pct), id_lab[ids[is.na(xa$pct)]]], collapse = ", "))
say("Xpo1 residual (log10, target and identity effects removed): ", paste(xa[!is.na(resid), sprintf("%s %+.2f", id_lab[cell_type], resid)], collapse = ", "))
pA <- ggplot(xa[!is.na(pct)], aes(identity, pct, colour = group)) +
  geom_segment(aes(xend = identity, y = 0, yend = pct), linewidth = .6) +
  geom_point(size = 2.6) +
  geom_text(aes(label = sprintf("%.1f", lg), y = pct + 6), size = 2.1, colour = ink2) +
  scale_colour_manual(values = c(`stem II / III` = col_a, `stem I` = col_b, differentiating = "#9a9892"), name = NULL) +
  scale_y_continuous(limits = c(0, 108), breaks = c(0, 25, 50, 75, 100)) +
  labs(x = NULL, y = "Xpo1 percentile among 206 targets", title = "A  Xpo1 across the differentiation axis",
       subtitle = paste0(c(genie3 = "GENIE3", scrank = "scRank", hdwgcna = "hdWGCNA")[method], " networks; labels give the raw log10 score")) +
  theme_np + theme(axis.text.x = element_text(angle = 40, hjust = 1), legend.position = "bottom",
                   plot.subtitle = element_text(size = 7, colour = ink2))

# ---- B: every target's effect in stem III (x) against stem II (y) -----------
# The residual is the log10 score after removing the target's and the
# identity's main effects: how much more (or less) this identity's network
# relies on the target than the target's average and the identity's average
# predict. Upper right = relatively more important in both cycling stem pools.
sb <- dcast(nf[cell_type %in% stem23], target ~ cell_type, value.var = "resid")[!is.na(stem_II) & !is.na(stem_III)]
sb[, class := fifelse(target == "Xpo1", "Xpo1", fifelse(target %in% paper, "tested in Mead et al.", fifelse(target %in% cycle, "cell-cycle target", "other")))]
sb[, class := factor(class, levels = c("other", "cell-cycle target", "tested in Mead et al.", "Xpo1"))]
sb[, both := (stem_II + stem_III) / 2]
sb[, rank_both := frank(-both)]
say("stem II vs stem III residuals: ", nrow(sb), " targets scored in both; Spearman ", sprintf("%.2f", cor(sb$stem_II, sb$stem_III, method = "spearman")),
    "; Xpo1 stem II ", sprintf("%+.2f", sb[target == "Xpo1", stem_II]), ", stem III ", sprintf("%+.2f", sb[target == "Xpo1", stem_III]),
    ", rank ", sb[target == "Xpo1", rank_both], " of ", nrow(sb), " on their mean")
say("upper-right quadrant (both > 0): ", sb[stem_II > 0 & stem_III > 0, .N], " targets; top 10 by mean: ",
    paste(sb[order(rank_both)][1:10, sprintf("%s (%+.2f / %+.2f)", target, stem_III, stem_II)], collapse = ", "))
say("paper and cell-cycle targets (stem III / stem II): ", paste(sb[class != "other"][order(rank_both), sprintf("%s %+.2f/%+.2f", target, stem_III, stem_II)], collapse = ", "))
# Only the targets of the compounds the paper itself tested: the XPO1
# inhibitor KPT-330 and its Aurora, CDK4/6, MEK and p53 arms. The axis range
# and the ranks come from all targets scored in both stem states, so dropping
# the rest of the library from the picture does not move these points or hide
# where they sit in it.
pb <- sb[class %in% c("Xpo1", "tested in Mead et al.")]
pb[, rank_II := frank(-sb$stem_II)[match(target, sb$target)]]
pb[, rank_III := frank(-sb$stem_III)[match(target, sb$target)]]
pb[, lab := sprintf("%s  %d/%d", target, rank_III, rank_II)]
say("panel B shows ", nrow(pb), " paper-tested targets; ranks out of ", nrow(sb), " scored in both stem states")
for (i in order(pb$rank_both)) say(sprintf("  %-7s stem III %+.2f (rank %3d), stem II %+.2f (rank %3d)", pb$target[i], pb$stem_III[i], pb$rank_III[i], pb$stem_II[i], pb$rank_II[i]))
pB <- ggplot(pb, aes(stem_III, stem_II)) +
  annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = Inf, fill = col_a, alpha = .06) +
  geom_hline(yintercept = 0, colour = ink2, linewidth = .3) + geom_vline(xintercept = 0, colour = ink2, linewidth = .3) +
  geom_abline(slope = 1, intercept = 0, colour = "#c3c2b7", linewidth = .3, linetype = "dashed") +
  geom_point(aes(colour = class), size = 2.6) +
  geom_text_repel(aes(label = lab, colour = class), size = 2.2, max.overlaps = 40, min.segment.length = .1,
                  segment.colour = "#c3c2b7", seed = 3, show.legend = FALSE) +
  scale_colour_manual(values = c(`tested in Mead et al.` = col_b, Xpo1 = col_a), breaks = c("Xpo1", "tested in Mead et al."), name = NULL) +
  labs(x = "Stem III: residual log10 score", y = "Stem II: residual log10 score",
       title = "B  The targets the paper tested",
       subtitle = paste0("labels: rank in stem III / stem II, of ", nrow(sb), " targets scored in both")) +
  coord_equal() + theme_np + theme(legend.position = "bottom", plot.subtitle = element_text(size = 7, colour = ink2))

fig <- pA | pB
fig <- fig + plot_layout(widths = c(1, 1.1))
ggsave(paste0("../Fig/fig4_mead_", method, ".pdf"), fig, width = 7.1, height = 3.9, device = cairo_pdf)
ggsave(paste0("../Fig/fig4_mead_", method, ".png"), fig, width = 7.1, height = 3.9, dpi = 300, bg = "white")
ggsave(paste0("fig4_mead_", method, "_preview.png"), fig, width = 7.1, height = 3.9, dpi = 150)
close(out)
invisible(file.rename("mead_stats_tmp.txt", stats_file))
