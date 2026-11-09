#!/usr/bin/env Rscript
# Figure 3 -- which targets single out the malignant progenitors (HSC-like,
# Prog-like, GMP-like) in the van Galen 2019 AML run (antagonist), for the
# network method given on the command line (default hdwgcna).
#   A  per target, log10(-like / normal) in each of the three progenitor pairs,
#      ranked by the median; targets with a score on both sides of >= 2 pairs
#   B  that matched-pair selectivity against specificity for the malignant
#      progenitor compartment (mean over the three -like identities minus the
#      mean over the other 16), single genes and joint targets told apart
#   C  the scores themselves, for the targets at the top of A, in the six
#      progenitor identities
# Run from publication/analysis:  Rscript fig3_progenitors.R [hdwgcna|genie3]
# Writes ../Fig/fig3_progenitors_<method>{,_slide}.{pdf,png} and fig3_stats_<method>.txt

suppressMessages({ library(data.table); library(ggplot2); library(patchwork); library(ggrepel) })
method <- commandArgs(trailingOnly = TRUE)[1]; if (is.na(method)) method <- "hdwgcna"
out <- file(sprintf("fig3_stats_%s.txt", method), "w"); say <- function(...) { cat(..., "\n", sep = ""); cat(..., "\n", sep = "", file = out) }

col_pair <- c(HSC = "#2a78d6", Prog = "#eb6834", GMP = "#1baf7a")   # categorical slots 1-3
red <- "#d03b3b"; ink <- "#0b0b0b"; ink2 <- "#52514e"; grid <- "#e6e5e1"; seq_ramp <- c("#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b")
theme_np <- theme_minimal(base_size = 7.5) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = grid, linewidth = .3),
        axis.text = element_text(colour = ink2), plot.title = element_text(face = "bold", size = 8),
        plot.title.position = "plot", legend.key.size = unit(8, "pt"))

vg <- fread(sprintf("../data/vangalen2021/%s/perbscore_all_targets.txt", method))
say("network method: ", method)
vg[, lg := log10(perb_score)]
vg[, floor := perb_score == 0 | perb_score < 1e-6 * median(perb_score[perb_score > 0])]
prog <- c("HSC", "Prog", "GMP"); like <- paste0(prog, "-like")
like6 <- grep("-like$", unique(vg$cell_type), value = TRUE)

# A zero means the target has no edge in that identity's network: a fact about
# the network, not a missing measurement. Dropping those pairs threw away most
# of the comparisons, and the ones it threw away hardest were the informative
# ones, where the target is wired in the malignant cell and absent from its
# normal counterpart. So a zero side is replaced by the run's smallest non-zero
# score and the ratio is reported as a bound (open symbols), not an exact
# value. A pair with zeros on both sides carries nothing and is still dropped.
floor_lg <- log10(min(vg[floor == FALSE, perb_score]))
say("scores: ", nrow(vg), " identity x target pairs, ", sum(vg$floor), " with no score (no edge to the target)")
say("a zero side is read as the run's smallest non-zero score, 10^", sprintf("%.1f", floor_lg))
vg[, lg_f := fifelse(floor, floor_lg, lg)]
w  <- dcast(vg, target ~ cell_type, value.var = "lg_f")
wf <- dcast(vg, target ~ cell_type, value.var = "floor")

# ---- A: matched progenitor pairs, zero sides kept as bounds ----------------
dd <- rbindlist(lapply(seq_along(prog), function(i) data.table(
  target = w$target, pair = prog[i], diff = w[[like[i]]] - w[[prog[i]]],
  z_like = wf[[like[i]]], z_norm = wf[[prog[i]]])))
dd[, kind := fifelse(!z_like & !z_norm, "both scored",
             fifelse(z_like & z_norm, "neither scored",
             fifelse(z_norm, "only the malignant side is wired", "only the normal side is wired")))]
say("progenitor pairs: ", dd[kind == "both scored", .N], " both scored, ",
    dd[kind == "only the malignant side is wired", .N], " malignant-only (lower bound), ",
    dd[kind == "only the normal side is wired", .N], " normal-only (upper bound), ",
    dd[kind == "neither scored", .N], " dropped")
dd <- dd[kind != "neither scored"][, bound := kind != "both scored"]
# A bound's size is set by how low the detection floor sits, not by biology,
# so bounded pairs are capped at 10x either way: they say which side is wired,
# not by how much.
cap <- 1
say("bounded pairs are capped at ", 10^cap, "x; ", dd[bound == TRUE & abs(diff) > cap, .N], " of ",
    dd[bound == TRUE, .N], " bounds hit the cap")
dd[bound == TRUE, diff := pmax(pmin(diff, cap), -cap)]
a <- dd[, .(pairs_scored = .N, exact = sum(!bound), higher = sum(diff > 0), median_diff = median(diff)), by = target]
say("targets with >= 2 usable pairs: ", a[pairs_scored >= 2, .N], " of ", uniqueN(vg$target),
    " (", dd[bound == FALSE, .N, by = target][N >= 2, .N], " when zero-sided pairs were dropped)")
a_all <- copy(a)[, joint := grepl(";", target)]
a <- a[pairs_scored >= 2]
a[, `:=`(joint = grepl(";", target), consistent = higher == pairs_scored | higher == 0,
         any_bound = exact < pairs_scored)]
a <- merge(a, dd[bound == FALSE, .(exact_median = median(diff)), by = target], by = "target", all.x = TRUE)
setorder(a, median_diff, exact_median, na.last = FALSE)
dd_a <- dd[target %in% a$target]
dd_a[, `:=`(target_f = factor(target, levels = a$target), pair = factor(pair, levels = prog))]
a[, target_f := factor(target, levels = a$target)]
say("Targets with a score on both sides of at least two progenitor pairs: ", nrow(a), " of ", uniqueN(vg$target))
say("Malignant progenitors higher in every scored pair (>= 2): ", paste(a[higher == pairs_scored][order(-median_diff), sprintf("%s (%.1fx, %d pairs)", target, 10^median_diff, pairs_scored)], collapse = ", "))
say("Normal progenitors higher in every scored pair: ", paste(a[higher == 0][order(median_diff), sprintf("%s (%.1fx, %d)", target, 10^median_diff, pairs_scored)], collapse = ", "))

pA <- ggplot(dd_a, aes(diff, target_f)) +
  geom_vline(xintercept = 0, colour = ink2, linewidth = .3) +
  geom_segment(data = a, aes(x = 0, xend = median_diff, y = target_f, yend = target_f), colour = "#c3c2b7", linewidth = .5) +
  geom_point(aes(colour = pair, shape = bound), size = 1.6, alpha = .9, stroke = .7) +
  geom_point(data = a, aes(median_diff, target_f), shape = 124, size = 3, colour = ink) +
  scale_colour_manual(values = col_pair, name = "-like vs. normal:", labels = c("HSC", "Prog", "GMP")) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1), labels = c("both scored", "one side unwired (bound)"), name = NULL) +
  scale_x_continuous(breaks = function(l) seq(ceiling(l[1]), floor(l[2]), by = if (diff(l) > 5) 2 else 1), labels = function(x) ifelse(x == 0, "1x", ifelse(x > 0, paste0(10^abs(x), "x"), paste0("1/", 10^abs(x), "x")))) +
  labs(x = "Malignant / normal score (bar: median)", y = NULL,
       title = "A  Malignant vs. normal progenitor, per pair") +
  theme_np + theme(axis.text.y = element_text(size = 6, face = ifelse(a$joint, "plain", "bold")), legend.position = "bottom",
                   legend.box = "vertical", panel.grid.major.y = element_blank())

# ---- B: matched-pair selectivity vs compartment specificity ----------------
# The compartment axis stays on the identities where the target is actually
# wired: "among the networks that carry it, is it stronger in the malignant
# progenitors?". Folding the zeros in here as well would make both axes of
# this panel measure the same thing -- whether the target is wired at all --
# and their correlation rose from 0.6 to 0.8 when we tried it.
b <- vg[floor == FALSE, .(n_like = sum(cell_type %in% like),
                          spec = mean(lg[cell_type %in% like]) - mean(lg[!cell_type %in% like]),
                          top3_like = sum(cell_type[order(-lg)][1:3] %in% like)), by = target][n_like >= 2]
ab <- merge(a[, .(target, median_diff, pairs_scored, higher, joint)], b, by = "target")
# BCL2 is called out in red: venetoclax's target, and the paper's
# counter-example, since its score follows expression rather than the
# population the drug acts on.
# BCL2 is called out in red. It is scored in all three malignant progenitors,
# so its compartment specificity is solid, but it has no score in normal HSC
# or Prog -- no edge in those networks -- so its pair selectivity rests on the
# GMP pair alone and it misses the two-pair rule the rest of the panel obeys.
# It is added back explicitly, with the pair count on its label.
if (!"BCL2" %in% ab$target && "BCL2" %in% a_all$target && "BCL2" %in% b$target) {
  ab <- rbind(ab, merge(a_all[target == "BCL2", .(target, median_diff, pairs_scored, higher, joint)], b, by = "target"), fill = TRUE)
}
ab[, class := fifelse(target == "BCL2", "BCL2", fifelse(joint, "joint target", "single gene"))]
ab[, class := factor(class, levels = c("single gene", "joint target", "BCL2"))]
ab[, label := fifelse(median_diff > 0.35 | spec > 0.35 | median_diff < -0.35, target, "")]
ab[target == "BCL2", label := sprintf("BCL2 (%d pair%s)", pairs_scored, fifelse(pairs_scored == 1, "", "s"))]
say(sprintf("Spearman(matched-pair selectivity, compartment specificity) over %d targets: %.2f", nrow(ab), cor(ab$median_diff, ab$spec, method = "spearman")))
say("BCL2: pair selectivity ", sprintf("%+.2f", ab[target == "BCL2", median_diff]), " log10 on ", ab[target == "BCL2", pairs_scored],
    " pair(s), compartment specificity ", sprintf("%+.2f", ab[target == "BCL2", spec]), " over ", ab[target == "BCL2", n_like],
    " of the 3 malignant progenitors")
say("Both above 0.3 log10 (2x): ", paste(ab[median_diff > 0.3 & spec > 0.3][order(-spec), sprintf("%s (pair %.1fx, compartment %.1fx)", target, 10^median_diff, 10^spec)], collapse = "; "))
pB <- ggplot(ab, aes(median_diff, spec)) +
  annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = Inf, fill = "#2a78d6", alpha = .06) +
  geom_hline(yintercept = 0, colour = ink2, linewidth = .3) + geom_vline(xintercept = 0, colour = ink2, linewidth = .3) +
  geom_point(aes(shape = class, colour = class, size = class)) +
  geom_text_repel(aes(label = label, colour = class), size = 2.1, max.overlaps = 30, min.segment.length = .2,
                  segment.colour = "#c3c2b7", seed = 1, show.legend = FALSE) +
  scale_shape_manual(values = c(`single gene` = 16, `joint target` = 17, BCL2 = 16), name = NULL) +
  scale_colour_manual(values = c(`single gene` = "#2a78d6", `joint target` = "#eb6834", BCL2 = "#d03b3b"), name = NULL) +
  scale_size_manual(values = c(`single gene` = 1.8, `joint target` = 1.8, BCL2 = 3), guide = "none") +
  labs(x = "Median log10 (malignant / normal), progenitor pairs", y = "Log10 (malignant progenitors / other identities)",
       title = "B  Pair selectivity vs. compartment specificity") +
  theme_np + theme(legend.position = "bottom")

# ---- C: the scores behind the top targets ---------------------------------
top <- a[higher == pairs_scored & median_diff > 0][order(-median_diff)][1:min(10, .N), target]
cc <- vg[target %in% top & cell_type %in% c(prog, like)]
cc[, `:=`(target_f = factor(target, levels = rev(top)),
          ct_f = factor(cell_type, levels = as.vector(rbind(prog, like))),
          fill = fifelse(floor, NA_real_, lg))]
pC <- ggplot(cc, aes(ct_f, target_f, fill = fill)) +
  geom_tile(colour = "white", linewidth = .5) +
  geom_text(aes(label = ifelse(is.na(fill), "", sprintf("%.1f", fill)),
                colour = !is.na(fill) & fill > quantile(fill, .6, na.rm = TRUE)), size = 2.2, show.legend = FALSE) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = ink)) +
  scale_fill_gradientn(colours = seq_ramp, na.value = "#f0efec", name = expression(log[10]~score), limits = quantile(cc$fill, c(.02, .98), na.rm = TRUE), oob = scales::squish) +
  labs(x = NULL, y = NULL, title = "C  The most selective targets, scores in the progenitor identities") +
  theme_np + theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 40, hjust = 1, face = ifelse(grepl("-like", levels(cc$ct_f)), "bold", "plain")),
                   axis.text.y = element_text(size = 6.5))
say("Panel C targets: ", paste(top, collapse = ", "))

# ---- D: malignant identities only, overall perturbation strength ----------
# Not selectivity. How hard each target perturbs the malignant cells
# themselves, as the mean raw log10 score over the six -like identities it is
# wired in; a target has to be wired in at least four of the six to be ranked.
like6 <- grep("-like$", unique(vg$cell_type), value = TRUE)
md <- vg[cell_type %in% like6, .(n_scored = sum(!floor), mean_lg = mean(lg_f)), by = target][n_scored >= 4]
setorder(md, -mean_lg); md[, rank := .I]
say("malignant-only strength: ", nrow(md), " targets wired in >= 4 of the 6 -like identities")
say("strongest 12: ", paste(md[1:12, sprintf("%s (%.1f, %d/6)", target, mean_lg, n_scored)], collapse = ", "))
say("where the selective targets land: ", paste(md[target %in% c("KIT","HDAC1","IDH1","MTOR","PSMB5","BCL2","FLT3","CD33","DNMT1")][order(rank), sprintf("%s %d/%d", target, rank, nrow(md))], collapse = ", "))
md_top <- md[1:20][, `:=`(target_f = factor(target, levels = rev(target)),
                          class = fifelse(target == "BCL2", "BCL2", fifelse(grepl(";", target), "joint target", "single gene")))]
pD <- ggplot(md_top, aes(mean_lg, target_f)) +
  geom_segment(aes(x = min(md_top$mean_lg) - .25, xend = mean_lg, yend = target_f), colour = "#c3c2b7", linewidth = .5) +
  geom_point(aes(colour = class), size = 2.2) +
  geom_text(aes(label = sprintf("%d/6", n_scored), x = mean_lg + .08), size = 2, colour = ink2, hjust = 0) +
  scale_colour_manual(values = c(`single gene` = "#2a78d6", `joint target` = "#eb6834", BCL2 = red), name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.12))) +
  labs(x = "Mean log10 score across the malignant identities", y = NULL,
       title = "D  Malignant cells only: overall perturbation strength",
       subtitle = "strongest 20 targets wired in >= 4 of the 6 -like identities") +
  theme_np + theme(legend.position = "bottom", panel.grid.major.y = element_blank(),
                   axis.text.y = element_text(size = 6))

fig <- (pA | pB) / (pC | pD) + plot_layout(heights = c(1.3, 1))
ggsave(sprintf("../Fig/fig3_progenitors_%s.pdf", method), fig, width = 7.4, height = 8.4, device = cairo_pdf)
ggsave(sprintf("../Fig/fig3_progenitors_%s.png", method), fig, width = 7.4, height = 8.4, dpi = 300, bg = "white")
ggsave(sprintf("fig3_progenitors_%s_preview.png", method), fig, width = 7.4, height = 8.4, dpi = 150)

# ---- slide version: 16:9, larger type, panel A trimmed to the extremes ----
# Fifty rows of target labels do not survive slide-sized type, so A keeps the
# twelve most malignant-selective and the six most normal-selective targets
# and says so; B, C and D are unchanged.
keep <- c(head(a[order(-median_diff), target], 12), tail(a[order(-median_diff), target], 6))
a_s <- a[target %in% keep]; setorder(a_s, median_diff)
dd_s <- dd[target %in% keep][, `:=`(target_f = factor(target, levels = a_s$target), pair = factor(pair, levels = prog))]
a_s[, target_f := factor(target, levels = a_s$target)]
pA_s <- pA %+% dd_s +
  labs(subtitle = "the 12 most malignant-selective and 6 most normal-selective targets") +
  theme(axis.text.y = element_text(size = 9, face = ifelse(a_s$joint, "plain", "bold")))
pA_s$layers[[2]]$data <- a_s   # the median bars and their lollipop stems
pA_s$layers[[4]]$data <- a_s

pB_s <- pB + labs(y = "Malignant progenitors vs. rest (log10)", x = "Median log10 (malignant / normal)")
fig_s <- (pA_s | pB_s) / (pC | pD) + plot_layout(heights = c(1.15, 1))
fig_s <- fig_s & theme(text = element_text(size = 11), plot.title = element_text(face = "bold", size = 13),
                       plot.subtitle = element_text(size = 9.5, colour = ink2),
                       axis.text = element_text(size = 9.5, colour = ink2),
                       legend.text = element_text(size = 10), legend.key.size = unit(12, "pt"))
ggsave(sprintf("../Fig/fig3_progenitors_%s_slide.pdf", method), fig_s, width = 13.3, height = 7.5, device = cairo_pdf)
ggsave(sprintf("../Fig/fig3_progenitors_%s_slide.png", method), fig_s, width = 13.3, height = 7.5, dpi = 300, bg = "white")
close(out)
