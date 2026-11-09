#!/usr/bin/env Rscript
# Supplementary figures -- hdWGCNA rank-score runs on four AML/bone-marrow
# cohorts re-annotated by Zeng et al. 2025 (Blood Cancer Discov 6(4):307).
#   figS_zeng_method      van Galen (Zeng object, hdWGCNA) against the
#                         van Galen 2021 GENIE3 and hdWGCNA runs
#     A  log10 score scatter, shared targets x identities
#     B  per-target Spearman rho across identities
#     C  malignant (-like) minus normal, per target, both comparisons
#   figS_zeng_datasets    consistency across beneyto / lasry / naldini / vangalen
#     A  pairwise Spearman at four levels (target mean, identity effect,
#        all pairs, target x identity residual)
#     B  identity main effect per dataset
#     C  no-score rate per identity x dataset
#     D  top identity per target, targets agreeing in >= 3 datasets
#     E  top-identity agreement against a label-permutation null
#     F  BCL2: within-identity rank among all targets, per dataset
# BCL2 (venetoclax) is highlighted throughout.
#
# Inputs  : ../data/Zeng2025/<ds>/<ds>.txt, ../data/vangalen2021/<method>/perbscore_all_targets.txt
# Outputs : ../Fig/figS_zeng_{method,datasets}.{pdf,png}, zeng_stats.txt
#
# Run from publication/analysis:  Rscript figS_zeng.R

suppressMessages({ library(data.table); library(ggplot2); library(patchwork) })
set.seed(1)
out <- file("zeng_stats.txt", "w")
say <- function(...) { cat(..., "\n", sep = ""); cat(..., "\n", sep = "", file = out) }

# ---- palette (dataviz reference palette, as fig2_vangalen.R) -------------
col_a <- "#2a78d6"; col_b <- "#eb6834"; ink2 <- "#52514e"; grid <- "#e6e5e1"
seq_ramp <- c("#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b")
div_ramp <- c("#b5411b", "#eb6834", "#f6c3a8", "#f7f6f2", "#a9cdf5", "#2a78d6", "#0d366b")
id_cols <- c(HSC_MPP = "#0d366b", GMP = "#2a78d6", ProMono = "#86b6ef", Monocyte = "#eb6834",
             cDC = "#b5411b", Early_Erythroid = "#3f9b5a", Late_Erythroid = "#a3d4a8")
theme_np <- theme_minimal(base_size = 7.5) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = grid, linewidth = .3),
        axis.text = element_text(colour = ink2), plot.title = element_text(face = "bold", size = 8),
        plot.title.position = "plot", legend.key.size = unit(8, "pt"))
no_score <- function(x) x == 0 | x < 1e-6 * median(x[x > 0])
hl <- "BCL2"; col_hl <- "#d6336c"
rho <- function(a, b) cor(a, b, method = "spearman", use = "complete.obs")
read_scores <- function(f) { x <- fread(f); x[, `:=`(lg = log10(perb_score), floor = no_score(perb_score))][] }

############################################################################
# 1. Method comparison on van Galen
############################################################################
runs <- list(genie3  = read_scores("../data/vangalen2021/genie3/perbscore_all_targets.txt"),
             hdw_old = read_scores("../data/vangalen2021/hdwgcna/perbscore_all_targets.txt"),
             zeng    = read_scores("../data/Zeng2025/vangalen/vangalen.txt"))
ref_lab <- c(genie3 = "GENIE3 (van Galen 2021 run)", hdw_old = "hdWGCNA (van Galen 2021 run)")
pairs <- data.table(normal = c("HSC", "Prog", "GMP", "ProMono", "Mono", "cDC"),
                    like   = c("HSC-like", "Prog-like", "GMP-like", "ProMono-like", "Mono-like", "cDC-like"))

mm <- rbindlist(lapply(names(ref_lab), function(r) {
  m <- merge(runs[[r]], runs$zeng, by = c("cell_type", "target"), suffixes = c(".ref", ".z"))
  m[, ref := ref_lab[r]]
}))
mm[, both := !floor.ref & !floor.z]
say("################ van Galen: Zeng object (hdWGCNA) vs earlier runs")
for (r in ref_lab) say(r, ": ", uniqueN(mm[ref == r, target]), " targets x ", uniqueN(mm[ref == r, cell_type]), " identities; both scored ",
                       mm[ref == r, sum(both)], "; Spearman ", sprintf("%.2f", mm[ref == r & both, rho(lg.ref, lg.z)]))

sc_lab <- mm[both == TRUE, .(lab = sprintf("rho = %.2f\nn = %d", rho(lg.ref, lg.z), .N)), by = ref]
pA <- ggplot(mm[both == TRUE], aes(lg.ref, lg.z)) +
  geom_point(size = .8, alpha = .6, colour = col_a, stroke = 0) +
  geom_point(data = mm[both == TRUE & target == hl], size = 1.6, colour = col_hl) +
  ggrepel::geom_text_repel(data = mm[both == TRUE & target == hl], aes(label = cell_type), size = 1.8, colour = col_hl,
                           min.segment.length = 0, segment.size = .2, max.overlaps = Inf, seed = 1) +
  geom_text(data = sc_lab, aes(x = -Inf, y = Inf, label = lab), hjust = -.1, vjust = 1.2, size = 2.3, colour = ink2) +
  facet_wrap(~ref, scales = "free_x") +
  labs(x = expression(log[10]~score~"(earlier run)"), y = expression(log[10]~score~"(Zeng, hdWGCNA)"),
       title = "A  Same identity x target pairs, van Galen (BCL2 in pink)") + theme_np

pt <- mm[both == TRUE, if (.N >= 5) .(rho = rho(lg.ref, lg.z), n = .N), by = .(ref, target)]
ord <- pt[ref == ref_lab["genie3"]][order(rho), target]
pt[, target := factor(target, levels = unique(c(ord, setdiff(target, ord))))]
say("per-target rho across identities, median: ", paste(pt[, sprintf("%s %.2f", ref, median(rho)), by = ref]$V1, collapse = "; "))
pB <- ggplot(pt, aes(rho, target, colour = ref)) +
  geom_vline(xintercept = 0, colour = ink2, linewidth = .3) +
  geom_point(size = 1.4) +
  scale_colour_manual(values = setNames(c(col_b, col_a), ref_lab), name = NULL) +
  labs(x = "Spearman rho across identities", y = NULL, title = "B  Per-target identity profile") +
  theme_np + theme(legend.position = "bottom", legend.direction = "vertical",
                   axis.text.y = element_text(colour = ifelse(levels(pt$target) == hl, col_hl, ink2),
                                              face = ifelse(levels(pt$target) == hl, "bold", "plain")))

delta <- function(x) rbindlist(lapply(seq_len(nrow(pairs)), function(i) {
  w <- dcast(x[cell_type %in% unlist(pairs[i]) & floor == FALSE], target ~ cell_type, value.var = "lg")
  if (all(unlist(pairs[i]) %in% names(w))) w[, .(pair = pairs$like[i], target, d = get(pairs$like[i]) - get(pairs$normal[i]))][!is.na(d)]
}))
dl <- rbindlist(lapply(names(ref_lab), function(r) {
  m <- merge(delta(runs[[r]]), delta(runs$zeng[target %in% runs[[r]]$target]), by = c("pair", "target"), suffixes = c(".ref", ".z"))
  m[, ref := ref_lab[r]]
}))
dl_lab <- dl[, .(lab = sprintf("rho = %.2f\nsign agrees %d/%d", rho(d.ref, d.z), sum(sign(d.ref) == sign(d.z)), .N)), by = ref]
say("malignant - normal delta, Zeng vs earlier: ", paste(dl_lab[, paste(ref, gsub("\n", ", ", lab))], collapse = "; "))
lim <- max(abs(c(dl$d.ref, dl$d.z)))
pC <- ggplot(dl, aes(d.ref, d.z)) +
  geom_hline(yintercept = 0, colour = grid) + geom_vline(xintercept = 0, colour = grid) +
  geom_abline(linetype = 2, colour = ink2, linewidth = .3) +
  geom_point(aes(colour = pair), size = 1) +
  geom_point(data = dl[target == hl], aes(colour = pair), size = 2.2, shape = 21, stroke = .7, fill = NA, colour = col_hl) +
  ggrepel::geom_text_repel(data = dl[target == hl], aes(label = pair), size = 1.8, colour = col_hl,
                           min.segment.length = 0, segment.size = .2, max.overlaps = Inf, box.padding = .4, seed = 1) +
  geom_text(data = dl_lab, aes(x = -Inf, y = Inf, label = lab), hjust = -.1, vjust = 1.2, size = 2.3, colour = ink2) +
  facet_wrap(~ref) + coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
  scale_colour_manual(values = c("#0d366b", "#2a78d6", "#86b6ef", "#eb6834", "#b5411b", "#3f9b5a"), name = NULL) +
  labs(x = expression(Delta~log[10]~"(earlier run)"), y = expression(Delta~log[10]~"(Zeng)"),
       title = "C  Malignant (-like) minus normal identity, per target (BCL2 ringed)") + theme_np

fig1 <- (pA / pC) | pB
fig1 <- fig1 + plot_layout(widths = c(2.3, 1))
ggsave("../Fig/figS_zeng_method.pdf", fig1, width = 7.4, height = 6.4, device = cairo_pdf)
ggsave("../Fig/figS_zeng_method.png", fig1, width = 7.4, height = 6.4, dpi = 300, bg = "white")

############################################################################
# 2. Consistency across the four Zeng datasets
############################################################################
dsn <- c("beneyto", "lasry", "naldini", "vangalen")
# van Galen keeps its own labels in the Zeng object; map the normal ones onto
# the Zeng annotation. Prog has no single counterpart and is left out.
map_vg <- c(HSC = "HSC_MPP", GMP = "GMP", ProMono = "ProMono", Mono = "Monocyte", cDC = "cDC",
            pDC = "pDC", earlyEry = "Early_Erythroid", lateEry = "Late_Erythroid")
d <- rbindlist(lapply(dsn, function(k) {
  x <- read_scores(sprintf("../data/Zeng2025/%s/%s.txt", k, k))[, ds := k]
  if (k == "vangalen") x <- x[cell_type %in% names(map_vg)][, cell_type := map_vg[cell_type]]
  x
}))
cts <- Reduce(intersect, split(d$cell_type, d$ds)); tgs <- Reduce(intersect, split(d$target, d$ds))
dall <- d[target %in% tgs]   # every identity, shared targets (panel F)
d <- d[cell_type %in% cts & target %in% tgs]
id_order <- c("HSC_MPP", "GMP", "ProMono", "Monocyte", "cDC", "Early_Erythroid", "Late_Erythroid")
say("\n################ Zeng 2025 cross-dataset (hdWGCNA)")
say("shared identities ", length(cts), ", shared targets ", length(tgs), "; no-score rate ",
    paste(d[, sprintf("%s %.2f", ds, mean(floor)), by = ds]$V1, collapse = ", "))

dd <- d[floor == FALSE]
dd[, res := resid(lm(lg ~ target + cell_type)), by = ds]
ide <- dd[, { f <- lm(lg ~ 0 + cell_type + target); co <- coef(f)[grep("^cell_type", names(coef(f)))]
              .(cell_type = sub("cell_type", "", names(co)), eff = co - mean(co)) }, by = ds]
pw <- function(M, lev) rbindlist(lapply(dsn, function(a) lapply(dsn, function(b)
  data.table(a = a, b = b, level = lev, r = rho(M[[a]], M[[b]])))) |> unlist(recursive = FALSE))
levs <- c("Target mean", "Identity effect", "All scored pairs", "Target x identity residual")
cm <- rbind(pw(dcast(dd[, .(m = mean(lg)), by = .(ds, target)], target ~ ds, value.var = "m"), levs[1]),
            pw(dcast(ide, cell_type ~ ds, value.var = "eff"), levs[2]),
            pw(dcast(dd, cell_type + target ~ ds, value.var = "lg"), levs[3]),
            pw(dcast(dd, cell_type + target ~ ds, value.var = "res"), levs[4]))
cm[, `:=`(level = factor(level, levels = levs), a = factor(a, dsn), b = factor(b, rev(dsn)))]
for (l in levs) say("pairwise rho, ", l, ": ", paste(sprintf("%.2f", cm[level == l & as.integer(a) < (5 - as.integer(b)), r]), collapse = " "))
cm2 <- cm[as.integer(a) <= 5 - as.integer(b)]
pD <- ggplot(cm2, aes(a, b, fill = r)) +
  geom_tile(colour = "white", linewidth = .4) +
  geom_text(aes(label = sprintf("%.2f", r), colour = abs(r) > .6), size = 2) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black"), guide = "none") +
  scale_fill_gradientn(colours = div_ramp, limits = c(-1, 1), name = "Spearman") +
  facet_wrap(~level, nrow = 1) + labs(x = NULL, y = NULL, title = "A  Agreement between datasets, by level") +
  theme_np + theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 40, hjust = 1))

ide[, cell_type := factor(cell_type, id_order)]
pE <- ggplot(ide, aes(cell_type, eff, colour = ds, group = ds)) +
  geom_hline(yintercept = 0, colour = ink2, linewidth = .3) +
  geom_line(linewidth = .4) + geom_point(size = 1.3) +
  scale_colour_manual(values = c(beneyto = col_b, lasry = "#3f9b5a", naldini = "#86b6ef", vangalen = "#0d366b"), name = NULL) +
  labs(x = NULL, y = expression(identity~effect~(log[10])), title = "B  Identity-wide offset") +
  theme_np + theme(axis.text.x = element_text(angle = 40, hjust = 1), legend.position = "top")

fr <- d[, .(f = mean(floor)), by = .(ds, cell_type)][, cell_type := factor(cell_type, rev(id_order))]
pF <- ggplot(fr, aes(ds, cell_type, fill = f)) +
  geom_tile(colour = "white", linewidth = .4) +
  geom_text(aes(label = sprintf("%.0f%%", 100 * f), colour = f > .5), size = 2) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black"), guide = "none") +
  scale_fill_gradientn(colours = seq_ramp, limits = c(0, 1), labels = scales::percent, name = "no score") +
  labs(x = NULL, y = NULL, title = "C  Pairs without a score") +
  theme_np + theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 40, hjust = 1))

top <- d[floor == FALSE][order(-perb_score)][, .SD[1], by = .(ds, target)][, .(ds, target, top = cell_type)]
# Only targets scored in at least one shared identity of every dataset
TT <- dcast(top, target ~ ds, value.var = "top"); TT <- TT[complete.cases(TT)]
M <- as.matrix(TT[, ..dsn])
agree <- apply(M, 1, function(x) max(table(x)))
stat <- function(M) { a <- apply(M, 1, function(x) max(table(x))); c(ge3 = sum(a >= 3), all4 = sum(a == 4)) }
obs <- stat(M); null <- t(replicate(2000, stat(apply(M, 2, sample))))
say("top identity agreement over ", nrow(M), " targets: >=3 datasets ", obs["ge3"], " (null ", round(mean(null[, "ge3"]), 1),
    ", p ", signif(mean(null[, "ge3"] >= obs["ge3"]), 2), "); all 4 ", obs["all4"], " (null ", round(mean(null[, "all4"]), 1),
    ", p ", signif(mean(null[, "all4"] >= obs["all4"]), 2), ")")
TT[, `:=`(agree = agree, mode = apply(M, 1, function(x) names(which.max(table(x)))))]
say("same top in all 4: ", paste(TT[agree == 4, paste0(target, " (", mode, ")")], collapse = ", "))
tl <- melt(TT[agree >= 3], id.vars = c("target", "agree", "mode"), variable.name = "ds", value.name = "top")
tl[, target := factor(target, levels = TT[agree >= 3][order(factor(mode, id_order), -agree, target), rev(target)])]
tl[, top := factor(top, id_order)]
pG <- ggplot(tl, aes(ds, target, fill = top)) +
  geom_tile(colour = "white", linewidth = .4) +
  geom_text(data = tl[, .(ds = "vangalen", lab = ifelse(agree[1] == 4, "4/4", "")), by = target], aes(ds, target, label = lab),
            inherit.aes = FALSE, nudge_x = .75, size = 1.9, colour = ink2) +
  scale_fill_manual(values = id_cols, name = "top identity", drop = FALSE) +
  scale_x_discrete(expand = expansion(add = c(.5, 1))) +
  labs(x = NULL, y = NULL, title = sprintf("D  Top identity, %d targets agreeing in >= 3/4", nrow(TT[agree >= 3]))) +
  theme_np + theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 40, hjust = 1), axis.text.y = element_text(size = 5.5))

nl <- data.table(ge3 = null[, "ge3"])
pH <- ggplot(nl, aes(ge3)) +
  geom_histogram(binwidth = 1, fill = "#cde2fb", colour = "white") +
  geom_vline(xintercept = obs["ge3"], colour = col_b, linewidth = .7) +
  annotate("text", x = obs["ge3"], y = Inf, label = sprintf("observed %d", obs["ge3"]), hjust = 1.05, vjust = 1.5, size = 2.3, colour = col_b) +
  labs(x = "targets with the same top identity in >= 3/4 datasets", y = "permutations",
       title = "E  Against shuffled top identities") + theme_np

# BCL2 across datasets. Raw scores are not comparable between runs, so show
# its rank among the shared targets within each identity (1 = top target);
# no-score pairs rank last and are drawn grey.
br <- dall[, .(target, floor, pct = frank(-perb_score, ties.method = "min") / .N, rk = frank(-perb_score, ties.method = "min"), n = .N), by = .(ds, cell_type)][target == hl]
br[, top := rk == min(rk[!floor]) & !floor, by = ds]
say("\n################ BCL2")
for (k in dsn) say(k, ": ", paste(br[ds == k][order(rk), sprintf("%s %s", cell_type, ifelse(floor, "no-score", sprintf("%d/%d", rk, n)))], collapse = ", "))
bt <- top[target == hl]; say("BCL2 top identity among the shared 7: ", paste(bt$ds, bt$top, sep = " ", collapse = ", "))
all_ids <- c(id_order[1:4], "LMPP", "MEP_MkP", "EoBasoMast", id_order[5], "pDC", "Early_Lymphoid", "ProB_PreB", id_order[6:7])
br[, cell_type := factor(cell_type, all_ids)][, ds := factor(ds, rev(dsn))]
pI <- ggplot(br, aes(cell_type, ds)) +
  geom_tile(aes(fill = ifelse(floor, NA, pct)), colour = "white", linewidth = .4) +
  geom_text(aes(label = ifelse(floor, "-", rk), colour = !floor & pct < .25), size = 2) +
  geom_tile(data = br[top == TRUE], fill = NA, colour = col_hl, linewidth = .7) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black"), guide = "none") +
  scale_fill_gradientn(colours = rev(seq_ramp), limits = c(0, 1), labels = scales::percent, na.value = "#f0efec",
                       name = "rank among\ntargets (top %)") +
  scale_x_discrete(drop = FALSE) +
  labs(x = NULL, y = NULL, title = sprintf("F  %s: rank among %d targets within each identity (pink box = best rank per dataset; - = no score)", hl, length(tgs))) +
  theme_np + theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 40, hjust = 1))

fig2 <- pD / ((pE / pF / pH) | pG) / pI + plot_layout(heights = c(1, 3.2, .75))
ggsave("../Fig/figS_zeng_datasets.pdf", fig2, width = 7.4, height = 10.6, device = cairo_pdf)
ggsave("../Fig/figS_zeng_datasets.png", fig2, width = 7.4, height = 10.6, dpi = 300, bg = "white")
close(out)
