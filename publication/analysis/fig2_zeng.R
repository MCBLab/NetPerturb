#!/usr/bin/env Rscript
# Figures 2-4 of the BMC Bioinformatics manuscript -- multi-cohort benchmark
# of the rank-score track on four AML cohorts from the Zeng et al. 2025 atlas
# (BoneMarrowMap broad identity x scMalignantFinder malignant/normal call;
# GENIE3, antagonist, 55 targets).
#   fig2_zeng_overview   design, reproducibility, network-size offsets
#     A identities scored per cohort; B target mean against the other cohorts;
#     C variance decomposition; D identity offset against cells used;
#     E agreement between cohorts at four levels
#   fig3_zeng_stem       targets that act on stem identities
#     A target x compartment map; B stem specificity per cohort;
#     C malignant against normal stem compartment; D sign consistency
#     against a label-permutation null; E robustness of the ranking
#   fig4_zeng_malignant  malignant identities against their normal counterparts
#     A per target over the 30 matched pairs; B blast markers minus other
#     targets per pair; C selected targets per pair
#   fig2_zeng_cohorts    the six-panel summary used by the Bioinformatics
#                        Advances draft (archive/)
# The target-specific effect throughout is the residual of log10 score ~
# target + identity, fitted per cohort.
#
# Inputs  : ../data/Zeng2025/<ds>/<ds>[_]perbscore_all_targets.txt
#           zeng2025/<ds>_{cells,runtime}.tsv (tables pulled from each run's report)
#           ../data/vangalen2021/genie3/perbscore_all_targets.txt (earlier van Galen run)
# Outputs : ../Fig/fig{2_zeng_overview,3_zeng_stem,4_zeng_malignant,2_zeng_cohorts}.{pdf,png},
#           zeng2025/target_effects.tsv, zeng_cohort_stats.txt (numbers quoted in the text)
#
# Run from publication/analysis:  Rscript fig2_zeng.R

suppressMessages({ library(data.table); library(ggplot2); library(patchwork) })
set.seed(1)
out <- file("zeng_cohort_stats.txt", "w")
say <- function(...) { cat(..., "\n", sep = ""); cat(..., "\n", sep = "", file = out) }

# ---- palette (dataviz reference palette, as the other figures) ----------
ink2 <- "#52514e"; grid <- "#e6e5e1"; na_fill <- "#f0efec"
seq_ramp <- c("#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b")
div_ramp <- c("#b5411b", "#eb6834", "#f6c3a8", "#f7f6f2", "#a9cdf5", "#2a78d6", "#0d366b")
dsn <- c("beneyto", "lasry", "naldini", "vangalen")
ds_lab <- c(beneyto = "Beneyto", lasry = "Lasry", naldini = "Naldini", vangalen = "van Galen")
ds_col <- c(beneyto = "#2a78d6", lasry = "#eb6834", naldini = "#1baf7a", vangalen = "#eda100")   # categorical slots 1-4
theme_np <- theme_minimal(base_size = 7.5) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = grid, linewidth = .3),
        axis.text = element_text(colour = ink2), plot.title = element_text(face = "bold", size = 8),
        plot.title.position = "plot", legend.key.size = unit(8, "pt"))
rho <- function(a, b) cor(a, b, method = "spearman", use = "pairwise.complete.obs")
no_score <- function(x) x == 0 | x < 1e-6 * median(x[x > 0])

# ---- hierarchy: BoneMarrowMap broad identities, three compartments --------
stem <- c("HSC_MPP", "LMPP")
prog <- c("Cycling_Progenitor", "Early_Lymphoid", "Early_GMP", "Late_GMP", "MEP", "Megakaryocyte_Precursor",
          "EoBasoMast_Precursor", "Early_Erythroid")
mature <- c("Late_Erythroid", "Pro-Monocyte", "Monocyte", "cDC", "pDC", "B", "Plasma_Cell", "Naive_T",
            "CD4_Memory_T", "CD8_Memory_T", "NK")
base_order <- c(stem, prog, mature)
compart_of <- function(b) factor(fifelse(b %in% stem, "Stem", fifelse(b %in% prog, "Progenitor", "Mature")),
                                 c("Stem", "Progenitor", "Mature"))

f <- c(beneyto = "beneyto/beneytoperbscore_all_targets.txt", lasry = "lasry/lasryperbscore_all_targets.txt",
       naldini = "naldini/naldiniperbscore_all_targets.txt", vangalen = "vangalen/vangalen_perbscore_all_targets.txt")
d <- rbindlist(lapply(dsn, function(k) fread(file.path("../data/Zeng2025", f[k]))[, ds := k]))
d[, `:=`(lg = log10(perb_score), floor = no_score(perb_score)), by = ds]
d[, `:=`(malig = grepl("-like$", cell_type), base = sub("-like$", "", cell_type))][, compart := compart_of(base)]
stopifnot(all(d$base %in% base_order))
cc <- rbindlist(lapply(dsn, function(k) fread(sprintf("zeng2025/%s_cells.tsv", k))[, ds := k]))[, cell_type := gsub(" ", "_", Identity)]
d <- merge(d, cc[, .(ds, cell_type, cells = `Cells used`, genes = `Genes in network`)], by = c("ds", "cell_type"))
single <- sort(unique(d[!grepl(";", target), target])); joint <- sort(unique(d[grepl(";", target), target]))

say("################ Zeng 2025 atlas, four cohorts, GENIE3, ", paste(unique(d$binding), collapse = "/"))
for (k in dsn) {
  x <- cc[ds == k]; kept <- x[Status == "kept"]
  rt <- fread(sprintf("zeng2025/%s_runtime.tsv", k))
  say(sprintf("%s: %d identities kept (%d normal, %d malignant), %s cells; %d dropped (< --min_cells); %d genes offered, %s per network; GENIE3 per identity mean %s, max %s",
              ds_lab[k], nrow(kept), sum(!grepl("-like$", kept$cell_type)), sum(grepl("-like$", kept$cell_type)),
              format(sum(kept$`Cells used`), big.mark = ","), sum(x$Status != "kept"), kept$gene4use[1],
              paste(range(kept$`Genes in network`), collapse = "-"), rt[Step == "GENIE3", Mean], rt[Step == "GENIE3", Max]))
}
say("total: ", uniqueN(d[, paste(ds, cell_type)]), " identity networks, ", format(sum(cc$`Cells used`), big.mark = ","), " cells, ",
    nrow(d), " scores; 55 targets (", length(single), " single, ", length(joint), " joint) in every cohort")
say("no score (< 1e-6 x median): ", sum(d$floor), " of ", nrow(d), " (", paste(d[floor == TRUE, .N, by = target][order(-N), paste(target, N)], collapse = ", "), ")")

# ---- decomposition: log10 score = target + identity + residual ------------
dd <- d[floor == FALSE]
dd[, res := resid(lm(lg ~ target + cell_type)), by = ds]
r2 <- dd[, .(t = summary(lm(lg ~ target))$r.squared, i = summary(lm(lg ~ cell_type))$r.squared), by = ds]
say("variance explained by target / identity: ", paste(r2[, sprintf("%s %.0f%%/%.0f%%", ds, 100 * t, 100 * i)], collapse = ", "))
ide <- dd[, { fit <- lm(lg ~ 0 + cell_type + target); co <- coef(fit)[grep("^cell_type", names(coef(fit)))]
              .(cell_type = sub("cell_type", "", names(co)), eff = co - mean(co)) }, by = ds]
ide <- merge(ide, unique(d[, .(ds, cell_type, malig, base, compart, cells, genes)]), by = c("ds", "cell_type"))
say("identity effect span (log10): ", paste(ide[, sprintf("%s %.2f", ds, diff(range(eff))), by = ds]$V1, collapse = ", "))
say("identity effect vs cells used, Spearman: ", paste(ide[, sprintf("%s %.2f", ds, rho(eff, cells)), by = ds]$V1, collapse = ", "),
    "; pooled ", sprintf("%.2f", ide[, rho(eff, cells)]))
say("identity effect vs genes in network, Spearman: ", paste(ide[, sprintf("%s %.2f", ds, rho(eff, genes)), by = ds]$V1, collapse = ", "))
say("cells used, median malignant vs normal identity: ", ide[malig == TRUE, median(cells)], " vs ", ide[malig == FALSE, median(cells)])
fmt_coef <- function(fit, nm) { s <- summary(fit)$coef; sprintf("%s %+.3f (P=%.2g)", nm, s[nm, 1], s[nm, 4]) }
f0 <- lm(eff ~ ds + compart + malig, ide); f1 <- lm(eff ~ ds + log10(cells) + compart + malig, ide)
say("identity effect ~ cohort + compartment + malignant: ", paste(sapply(c("compartProgenitor", "compartMature", "maligTRUE"), fmt_coef, fit = f0), collapse = "; "))
say("  + log10(cells): ", paste(sapply(c("log10(cells)", "compartProgenitor", "compartMature", "maligTRUE"), fmt_coef, fit = f1), collapse = "; "),
    sprintf("; R2 %.2f -> %.2f", summary(f0)$r.squared, summary(f1)$r.squared))

# ---- agreement between cohorts -------------------------------------------
shared <- Reduce(intersect, split(dd$cell_type, dd$ds))
say("identities scored in all four cohorts (", length(shared), "): ", paste(shared, collapse = ", "))
pw <- function(W, lev) rbindlist(combn(dsn, 2, function(p) data.table(a = p[1], b = p[2], level = lev, r = rho(W[[p[1]]], W[[p[2]]])), simplify = FALSE))
stem_spec <- function(x) x[, .(s = mean(res[compart == "Stem"]) - mean(res[compart != "Stem"])), by = .(ds, target)]
S <- dcast(stem_spec(dd), target ~ ds, value.var = "s")
levs <- c("Target mean", "Identity offset", "Target x identity", "Stem specificity")
agr <- rbind(pw(dcast(dd[, .(m = mean(lg)), by = .(ds, target)], target ~ ds, value.var = "m"), levs[1]),
             pw(dcast(ide[cell_type %in% shared], cell_type ~ ds, value.var = "eff"), levs[2]),
             pw(dcast(dd[cell_type %in% shared], cell_type + target ~ ds, value.var = "res"), levs[3]),
             pw(S, levs[4]))
agr[, level := factor(level, levs)]
for (l in levs) say("pairwise Spearman, ", l, ": ", paste(sprintf("%.2f", range(agr[level == l, r])), collapse = "-"), " (median ", sprintf("%.2f", median(agr[level == l, r])), ")")
old <- fread("../data/vangalen2021/genie3/perbscore_all_targets.txt")[, floor := no_score(perb_score)][floor == FALSE]
vo <- merge(old[, .(o = mean(log10(perb_score))), by = target], dd[ds == "vangalen", .(n = mean(lg)), by = target], by = "target")
say("target mean, this van Galen run vs the earlier run on the original labels (Supplementary Fig.): Spearman ", sprintf("%.2f", rho(vo$o, vo$n)), " over ", nrow(vo), " targets")

# ---- stem specificity: permutation of identity labels within cohort -------
R <- dcast(dd, ds + cell_type ~ target, value.var = "res")
lab0 <- compart_of(sub("-like$", "", R$cell_type)) == "Stem"
say("stem identities per cohort: ", paste(R[lab0, .(x = paste(cell_type, collapse = "+")), by = ds][, paste0(ds, " ", x)], collapse = "; "))
spec_of <- function(lab) sapply(dsn, function(k) { i <- R$ds == k; X <- as.matrix(R[i, -(1:2)])
  colMeans(X[lab[i], , drop = FALSE], na.rm = TRUE) - colMeans(X[!lab[i], , drop = FALSE], na.rm = TRUE) })
S0 <- spec_of(lab0); m0 <- rowMeans(S0)
B <- 10000
perm <- replicate(B, { lab <- lab0; for (k in dsn) { i <- which(R$ds == k); lab[i] <- sample(lab0[i]) }; rowMeans(spec_of(lab)) })
# a permutation that puts only no-score pairs on the stem side gives NA; it is left out
nB <- rowSums(!is.na(perm))
te <- data.table(target = names(m0), stem = m0, stem_pos = rowSums(S0 > 0),
                 p_stem_up = (rowSums(perm >= m0, na.rm = TRUE) + 1) / (nB + 1), p_stem_dn = (rowSums(perm <= m0, na.rm = TRUE) + 1) / (nB + 1))
te[, `:=`(q_stem_up = p.adjust(p_stem_up, "BH"), q_stem_dn = p.adjust(p_stem_dn, "BH"))]
te[, stem_class := fifelse(stem_pos == 4 & p_stem_up < .05, "stem", fifelse(stem_pos == 0 & p_stem_dn < .05, "not stem", ""))]
say("\n################ stem specificity (mean residual, stem minus other identities; ", B, " within-cohort label permutations)")
say("FDR < 0.05, towards stem: ", paste(te[q_stem_up < .05][order(-stem), sprintf("%s %+.2f (q=%.2g)", target, stem, q_stem_up)], collapse = ", "),
    "; away from stem: ", paste(te[q_stem_dn < .05][order(stem), sprintf("%s %+.2f (q=%.2g)", target, stem, q_stem_dn)], collapse = ", "))
say("positive in 4/4 cohorts and permutation P < 0.05: ", paste(te[stem_class == "stem"][order(-stem), sprintf("%s %+.2f (P=%.3f)", target, stem, p_stem_up)], collapse = ", "))
say("positive in 4/4, P >= 0.05: ", paste(te[stem_pos == 4 & stem_class == ""][order(-stem), sprintf("%s %+.2f (P=%.2f)", target, stem, p_stem_up)], collapse = ", "))
say("negative in 4/4 cohorts: ", paste(te[stem_pos == 0][order(stem), sprintf("%s %+.2f (P=%.3f)", target, stem, p_stem_dn)], collapse = ", "))
set_test <- function(g, nm) { obs <- mean(m0[g]); nul <- colMeans(perm[g, , drop = FALSE], na.rm = TRUE)
  say(sprintf("set %s (%s): mean stem specificity %+.2f, permutation P (two-sided) = %.3g", nm, paste(g, collapse = ","), obs,
              min(1, 2 * min(mean(nul >= obs), mean(nul <= obs)) + 1 / length(nul)))) }
set_test(c("AURKA", "AURKB", "PLK1", "CHEK1", "WEE1"), "mitotic/checkpoint kinases")
set_test(c("KIT", "FLT3", "CD33"), "blast markers (pre-specified in the earlier van Galen run)")
hs <- dd[, .(s = mean(res[base == "HSC_MPP"]) - mean(res[compart != "Stem"])), by = .(ds, target)][, .(s = mean(s)), by = target]
say("robustness: stem = HSC/MPP only, Spearman with HSC/MPP+LMPP ", sprintf("%.2f", rho(hs$s, m0[hs$target])))
dd[, pct := frank(-lg) / .N, by = .(ds, cell_type)]
rk <- dd[, .(s = mean(pct[compart != "Stem"]) - mean(pct[compart == "Stem"])), by = .(ds, target)][, .(s = mean(s)), by = target]
say("robustness: within-identity rank instead of residual, Spearman ", sprintf("%.2f", rho(rk$s, m0[rk$target])))

# ---- malignant minus matched normal identity -----------------------------
mp <- dcast(dd[, .(ds, target, base, malig, res)], ds + target + base ~ malig, value.var = "res")[!is.na(`TRUE`) & !is.na(`FALSE`)]
mp[, dl := `TRUE` - `FALSE`]
say("\n################ malignant (-like) minus matched normal identity (residual)")
say("matched pairs: ", uniqueN(mp[, paste(ds, base)]), " (", paste(mp[, .(n = uniqueN(base)), by = ds][, paste(ds, n)], collapse = ", "), ")")
mw <- mp[, .(mal = mean(dl), mal_pos = mean(dl > 0), n_pairs = .N, p_mal = wilcox.test(dl)$p.value), by = target][, q_mal := p.adjust(p_mal, "BH")]
say("top 6: ", paste(mw[order(-mal)][1:6, sprintf("%s %+.2f (%d%% of pairs, P=%.2g, q=%.2g)", target, mal, round(100 * mal_pos), p_mal, q_mal)], collapse = ", "))
say("bottom 4: ", paste(mw[order(mal)][1:4, sprintf("%s %+.2f (P=%.2g, q=%.2g)", target, mal, p_mal, q_mal)], collapse = ", "))
bm <- c("KIT", "FLT3", "CD33")
pp <- mp[, .(d = mean(dl[target %in% bm]) - mean(dl[!target %in% bm])), by = .(ds, base)]
say(sprintf("blast markers (KIT,FLT3,CD33) minus other targets, per pair: median %+.2f, higher in %d/%d pairs, signed-rank P=%.2g",
            median(pp$d), sum(pp$d > 0), nrow(pp), wilcox.test(pp$d)$p.value))
sm <- mp[base %in% stem][, .(d = mean(dl)), by = .(ds, target)]
sw <- dcast(sm, target ~ ds, value.var = "d"); sc <- setdiff(names(sw), "target")
sw[, `:=`(stem_mal = rowMeans(.SD), stem_mal_pos = rowSums(.SD > 0), stem_mal_n = rowSums(!is.na(.SD))), .SDcols = sc]
say("malignant minus normal stem (", paste(mp[base %in% stem, unique(paste(ds, base))], collapse = ", "), "), higher in every cohort: ",
    paste(sw[stem_mal_pos == stem_mal_n][order(-stem_mal), sprintf("%s %+.2f", target, stem_mal)], collapse = ", "))
te <- merge(merge(te, mw, by = "target"), sw[, .(target, stem_mal, stem_mal_pos, stem_mal_n)], by = "target", all.x = TRUE)

# ---- per-compartment map --------------------------------------------------
cm <- dd[, .(r = mean(res)), by = .(ds, target, compart, malig)][, .(r = mean(r), n = .N), by = .(target, compart, malig)]
cm[, col := factor(paste(compart, ifelse(malig, "M", "N")), paste(rep(c("Stem", "Progenitor", "Mature"), each = 2), c("M", "N")))]
for (cmp in levels(cm$col)) te[, (paste0("map_", gsub(" ", "_", cmp))) := cm[col == cmp][match(te$target, target), r]]
say("\n################ stem targets: malignant vs normal stem compartment (map, mean over cohorts; normal stem absent in Beneyto)")
for (k in te[stem_class != "" | stem_pos == 4][order(-stem), target])
  say(sprintf("%-12s stem-spec %+.2f | stem malig. %+.2f, stem normal %+.2f | malig.-normal stem %+.2f (%d/%d cohorts > 0)",
              k, te[target == k, stem], te[target == k, map_Stem_M], te[target == k, map_Stem_N],
              te[target == k, stem_mal], te[target == k, stem_mal_pos], te[target == k, stem_mal_n]))
fwrite(te[order(-stem)], "zeng2025/target_effects.tsv", sep = "\t")

# ---- further numbers for the text ---------------------------------------
say("\n################ further detail")
vd <- dd[, { ss <- anova(lm(lg ~ target + cell_type))[["Sum Sq"]]; .(Target = ss[1] / sum(ss), Identity = ss[2] / sum(ss), Residual = ss[3] / sum(ss)) }, by = ds]
say("variance shares target/identity/residual: ", paste(vd[, sprintf("%s %.0f/%.0f/%.0f%%", ds, 100 * Target, 100 * Identity, 100 * Residual)], collapse = ", "))
sign_cons <- function(S) sum(rowSums(S > 0) == 4 | rowSums(S < 0) == 4, na.rm = TRUE)
B2 <- 2000
cons_null <- replicate(B2, { lab <- lab0; for (k in dsn) { i <- which(R$ds == k); lab[i] <- sample(lab0[i]) }; sign_cons(spec_of(lab)) })
cons_obs <- sign_cons(S0)
say(sprintf("stem specificity with the same sign in all four cohorts: %d targets (%d towards, %d away); %d permutations: mean %.1f, P = %.3g",
            cons_obs, sum(rowSums(S0 > 0) == 4), sum(rowSums(S0 < 0) == 4), B2, mean(cons_null), (sum(cons_null >= cons_obs) + 1) / (B2 + 1)))
jt <- rbindlist(lapply(dsn, function(k) rbindlist(lapply(joint, function(t) {
  g <- strsplit(t, ";")[[1]]; w <- dcast(dd[ds == k & target %in% c(t, g)], cell_type ~ target, value.var = "lg")
  if (!all(c(t, g) %in% names(w))) return(NULL)
  m <- sapply(g, function(z) mean(w[[z]], na.rm = TRUE)); st <- g[which.max(m)]; wk <- g[which.min(m)]
  data.table(ds = k, target = t, r_strong = rho(w[[t]], w[[st]]), r_weak = rho(w[[t]], w[[wk]]),
             fold = 10^median(w[[t]] - pmax(w[[g[1]]], w[[g[2]]]), na.rm = TRUE)) }))))
say(sprintf("joint targets: Spearman across identities with the stronger member %.2f, the weaker %.2f (medians over %d joint x cohort); score %.1f-fold the higher member",
            median(jt$r_strong), median(jt$r_weak), nrow(jt), median(jt$fold)))
pp <- merge(pp, unique(dd[, .(base, compart)]), by = "base")
say("blast-marker excess by compartment: ", paste(pp[, sprintf("%s %d/%d pairs, median %+.2f", compart[1], sum(d > 0), .N, median(d)), by = compart]$V1, collapse = "; "))
cl <- mp[, .(d = mean(dl)), by = .(ds, target)][, .(pos = sum(d > 0), n = .N), by = target]
say("malignant > normal in the cohort mean of matched pairs, all four cohorts: ", paste(cl[pos == 4, target], collapse = ", "))

############################################################################
# Figure
############################################################################
ord <- function(tg) te[target %in% tg][order(stem), target]
lv <- c(ord(single), ord(joint))
grp <- function(t) factor(ifelse(grepl(";", t), "Joint", "Single gene"), c("Single gene", "Joint"))

# A coverage
cv <- merge(CJ(ds = dsn, base = base_order, status = c("N", "M")), cc[Status == "kept", .(ds, base = sub("-like$", "", cell_type),
            status = ifelse(grepl("-like$", cell_type), "M", "N"), cells = `Cells used`)], all.x = TRUE, by = c("ds", "base", "status"))
cv <- cv[base %in% cv[!is.na(cells), unique(base)]]
cv[, `:=`(base = factor(base, rev(base_order)), status = factor(status, c("N", "M")), dsl = factor(ds_lab[ds], ds_lab))]
pA <- ggplot(cv, aes(status, base, fill = cells)) +
  geom_tile(colour = "white", linewidth = .5) +
  scale_fill_gradientn(colours = seq_ramp, trans = "log10", na.value = na_fill, breaks = c(200, 1000, 3000), name = "cells") +
  facet_grid(~dsl) +
  scale_y_discrete(labels = function(x) gsub("_", " ", x)) +
  labs(x = "N normal, M malignant", y = NULL, title = "A  Identities scored") +
  theme_np + theme(panel.grid = element_blank(), strip.text = element_text(size = 6.5), panel.spacing = unit(2, "pt"),
                   legend.position = "bottom", legend.key.width = unit(14, "pt"), legend.key.height = unit(5, "pt"))

# B offset vs cells
ide[, `:=`(dsl = factor(ds_lab[ds], ds_lab), st = ifelse(malig, "malignant", "normal"))]
rb <- ide[, .(lab = sprintf("%s (rho %.2f)", ds_lab[ds[1]], rho(eff, cells))), by = ds]
pB <- ggplot(ide, aes(cells, eff)) +
  geom_hline(yintercept = 0, colour = ink2, linewidth = .3) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = ink2, linewidth = .4, linetype = 2) +
  geom_point(aes(colour = dsl, shape = st), size = 1.6, stroke = .6) +
  scale_x_log10(breaks = c(200, 500, 1000, 3000)) +
  scale_colour_manual(values = setNames(ds_col, ds_lab), labels = setNames(rb$lab, ds_lab[rb$ds]), name = NULL) +
  scale_shape_manual(values = c(malignant = 16, normal = 1), name = NULL) +
  guides(colour = guide_legend(nrow = 2, order = 1), shape = guide_legend(order = 2)) +
  labs(x = "cells in the identity's network", y = expression(identity~offset~(log[10])), title = "B  Offsets follow network size") +
  theme_np + theme(legend.position = "bottom", legend.box = "vertical", legend.spacing.y = unit(0, "pt"), legend.margin = margin())

# C agreement
pC <- ggplot(agr, aes(r, level)) +
  geom_vline(xintercept = 0, colour = ink2, linewidth = .3) +
  geom_point(size = 1.4, colour = "#2a78d6", alpha = .8) +
  stat_summary(fun = median, geom = "point", shape = 124, size = 4, colour = "#0d366b") +
  scale_y_discrete(limits = rev(levs)) + scale_x_continuous(limits = c(-.1, 1), breaks = c(0, .5, 1)) +
  labs(x = "Spearman, 6 cohort pairs", y = NULL, title = "C  Agreement between cohorts") + theme_np

# D compartment map, E stem specificity, F malignant - normal (shared target axis)
cm[, `:=`(target = factor(target, lv), g = grp(target))]
lim <- quantile(abs(cm$r), .98)
pD <- ggplot(cm, aes(col, target, fill = pmax(pmin(r, lim), -lim))) +
  geom_tile(colour = "white", linewidth = .3) +
  scale_fill_gradientn(colours = rev(div_ramp), limits = c(-lim, lim), name = "target-specific\neffect (log10)") +
  facet_grid(g ~ ., scales = "free_y", space = "free_y") +
  scale_x_discrete(labels = function(x) sub(" M$", "\nmalig.", sub(" N$", "\nnormal", x))) +
  labs(x = NULL, y = NULL, title = "D  Where each target acts") +
  theme_np + theme(panel.grid = element_blank(), strip.text = element_blank(), axis.text.x = element_text(size = 5.6, angle = 90, hjust = 1, vjust = .5),
                   axis.text.y = element_text(size = 5.6, colour = ink2),
                   legend.position = "bottom", legend.title.position = "top", legend.key.width = unit(12, "pt"), legend.key.height = unit(5, "pt"))
se <- melt(S, id.vars = "target", variable.name = "ds", value.name = "s")[, `:=`(target = factor(target, lv), g = grp(target), dsl = factor(ds_lab[ds], ds_lab))]
sm_ <- te[, .(target = factor(target, lv), g = grp(target), stem, stem_class = fifelse(stem_class == "", ".", stem_class))]
pE <- ggplot(se, aes(s, target)) +
  geom_vline(xintercept = 0, colour = ink2, linewidth = .3) +
  geom_point(aes(colour = dsl), size = .9, alpha = .85) +
  geom_point(data = sm_, aes(stem, target, fill = stem_class), shape = 23, size = 1.6, colour = "black", stroke = .3, inherit.aes = FALSE) +
  scale_fill_manual(values = c(stem = "#0d366b", `not stem` = "#b5411b", `.` = "white"), breaks = c("stem", "not stem"),
                    labels = c("4/4 cohorts towards stem, P<0.05", "4/4 away from stem, P<0.05"), name = NULL) +
  scale_colour_manual(values = setNames(ds_col, ds_lab), guide = "none") +
  facet_grid(g ~ ., scales = "free_y", space = "free_y") +
  labs(x = "stem minus other identities", y = NULL, title = "E  Stem specificity") +
  theme_np + theme(axis.text.y = element_blank(), strip.text = element_blank(), legend.position = "bottom", legend.direction = "vertical")
mpl <- mp[, .(d = mean(dl)), by = .(ds, target)][, `:=`(target = factor(target, lv), g = grp(target), dsl = factor(ds_lab[ds], ds_lab))]
mwl <- mw[, .(target = factor(target, lv), g = grp(target), mal, sig = p_mal < .05)]
pF <- ggplot(mpl, aes(d, target)) +
  geom_vline(xintercept = 0, colour = ink2, linewidth = .3) +
  geom_point(aes(colour = dsl), size = .9, alpha = .85) +
  geom_point(data = mwl, aes(mal, target, fill = sig), shape = 23, size = 1.6, colour = "black", stroke = .3, inherit.aes = FALSE) +
  scale_fill_manual(values = c(`TRUE` = "#0d366b", `FALSE` = "white"), labels = c(`TRUE` = "P<0.05 over 30 pairs", `FALSE` = "n.s."), name = NULL) +
  scale_colour_manual(values = setNames(ds_col, ds_lab), name = NULL, guide = "none") +
  facet_grid(g ~ ., scales = "free_y", space = "free_y") +
  labs(x = "malignant minus normal", y = NULL, title = "F  Malignant vs normal") +
  theme_np + theme(axis.text.y = element_blank(), strip.text = element_blank(), legend.position = "bottom", legend.direction = "vertical")

top <- pA + pB + pC + plot_layout(widths = c(1.25, 1, 1))
bot <- pD + pE + pF + plot_layout(widths = c(.55, 1, 1))
fig <- top / bot + plot_layout(heights = c(1, 2.6))
ggsave("../Fig/fig2_zeng_cohorts.pdf", fig, width = 7.2, height = 8.6, device = cairo_pdf)
ggsave("../Fig/fig2_zeng_cohorts.png", fig, width = 7.2, height = 8.6, dpi = 300, bg = "white")

save_fig <- function(p, name, w, h) {
  ggsave(sprintf("../Fig/%s.pdf", name), p, width = w, height = h, device = cairo_pdf)
  ggsave(sprintf("../Fig/%s.png", name), p, width = w, height = h, dpi = 300, bg = "white")
}
retitle <- function(p, t) p + labs(title = t)

############################################################################
# Figure 2 -- design, reproducibility, network-size offsets
############################################################################
tm <- dd[, .(m = mean(lg)), by = .(ds, target)]
tm[, other := sapply(seq_len(.N), function(i) tm[ds != ds[i] & target == target[i], mean(m)])]
tm[, dsl := factor(ds_lab[ds], ds_lab)]
rt <- tm[, .(lab = sprintf("%s %.2f", ds_lab[ds[1]], rho(m, other))), by = ds]
p2B <- ggplot(tm, aes(other, m, colour = dsl)) +
  geom_abline(linetype = 2, colour = ink2, linewidth = .3) +
  geom_point(size = 1.1, alpha = .8) +
  scale_colour_manual(values = setNames(ds_col, ds_lab), guide = "none") +
  annotate("text", x = -Inf, y = Inf, hjust = -.05, vjust = 1.15, size = 2.1, colour = ink2,
           label = paste0("Spearman\n", paste(rt$lab, collapse = "\n"))) +
  labs(x = expression(mean~log[10]~score*","~other~three~cohorts), y = expression(mean~log[10]~score*","~this~cohort),
       title = "B  Target means replicate") + theme_np
vdl <- melt(vd, id.vars = "ds", variable.name = "part", value.name = "f")[, dsl := factor(ds_lab[ds], rev(ds_lab))]
p2C <- ggplot(vdl, aes(f, dsl, fill = part)) +
  geom_col(width = .7, colour = "white", linewidth = .5, position = position_stack(reverse = TRUE)) +
  geom_text(aes(label = ifelse(f > .15, sprintf("%.0f%%", 100 * f), "")), position = position_stack(vjust = .5, reverse = TRUE), size = 2.1,
            colour = ifelse(vdl$part == "Residual", "black", "white")) +
  scale_fill_manual(values = c(Target = "#2a78d6", Identity = "#eb6834", Residual = "#d6d4cf"), name = NULL,
                    labels = c(Target = "Target", Identity = sprintf("Identity (%.0f-%.0f%%)", 100 * min(vd$Identity), 100 * max(vd$Identity)), Residual = "Residual")) +
  scale_x_continuous(labels = scales::percent, expand = c(0, 0)) +
  labs(x = expression(share~of~variance~of~log[10]~scores), y = NULL, title = "C  What the scores vary with") +
  theme_np + theme(legend.position = "bottom", panel.grid.major.y = element_blank())
fig2 <- (retitle(pA, "A  Identities scored") | ((p2B | p2C) / (retitle(pB, "D  Offsets follow network size") | retitle(pC, "E  Agreement between cohorts")))) +
  plot_layout(widths = c(1, 2.1))
save_fig(fig2, "fig2_zeng_overview", 7.2, 6.2)

############################################################################
# Figure 3 -- targets that act on stem identities
############################################################################
cls_col <- c(stem = "#0d366b", `not stem` = "#b5411b", other = "#b0aea8")
sc3 <- te[, .(target, x = map_Stem_N, y = map_Stem_M, cls = fifelse(stem_class == "", "other", stem_class), joint = grepl(";", target))]
lab3 <- sc3[cls != "other" | target %in% c("FLT3", "HMGCR", "TP53", "BCL2", "CHEK1")]
lim3 <- max(abs(c(sc3$x, sc3$y)), na.rm = TRUE) * 1.05
p3C <- ggplot(sc3, aes(x, y)) +
  geom_hline(yintercept = 0, colour = grid) + geom_vline(xintercept = 0, colour = grid) +
  geom_abline(linetype = 2, colour = ink2, linewidth = .3) +
  geom_point(aes(colour = cls, shape = joint), size = 1.6, stroke = .6) +
  ggrepel::geom_text_repel(data = lab3, aes(label = target), size = 2, colour = "black", min.segment.length = 0,
                           segment.size = .2, max.overlaps = Inf, box.padding = .3, seed = 1) +
  annotate("text", x = -lim3, y = lim3, hjust = 0, vjust = 1, size = 2.1, colour = ink2, label = "malignant stem") +
  annotate("text", x = lim3, y = -lim3, hjust = 1, vjust = 0, size = 2.1, colour = ink2, label = "normal stem") +
  scale_colour_manual(values = cls_col, breaks = c("stem", "not stem"), labels = c("towards stem, 4/4", "away from stem, 4/4"), name = NULL) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1), labels = c("single gene", "joint"), name = NULL) +
  coord_equal(xlim = c(-lim3, lim3), ylim = c(-lim3, lim3)) +
  labs(x = "effect in normal stem identities", y = "effect in malignant stem identities", title = "C  Malignant or normal stem cells") +
  theme_np + theme(legend.position = "bottom", legend.box = "vertical", legend.spacing.y = unit(0, "pt"), legend.margin = margin())
p3D <- ggplot(data.table(n = cons_null), aes(n)) +
  geom_histogram(binwidth = 1, fill = "#a9cdf5", colour = "white") +
  geom_vline(xintercept = cons_obs, colour = "#0d366b", linewidth = .7) +
  annotate("text", x = cons_obs, y = Inf, hjust = 1.05, vjust = 1.4, size = 2.1, colour = "#0d366b",
           label = sprintf("observed %d", cons_obs)) +
  labs(x = "targets with the same sign in all 4 cohorts", y = "permutations", title = "D  Consistency beyond chance") + theme_np
rob <- rbind(data.table(target = hs$target, alt = hs$s, main = m0[hs$target], what = sprintf("HSC/MPP only, rho %.2f", rho(hs$s, m0[hs$target]))),
             data.table(target = rk$target, alt = rk$s, main = m0[rk$target], what = sprintf("ranks, rho %.2f", rho(rk$s, m0[rk$target]))))
p3E <- ggplot(rob, aes(main, alt)) +
  geom_hline(yintercept = 0, colour = grid) + geom_vline(xintercept = 0, colour = grid) +
  geom_point(size = 1, colour = "#2a78d6", alpha = .8) +
  geom_point(data = rob[target %in% c("IL3RA", "KDM1A", "KIT")], size = 1.4, colour = "#0d366b") +
  facet_wrap(~what, scales = "free_y") +
  labs(x = "stem specificity (HSC/MPP + LMPP, residual)", y = "alternative definition", title = "E  Robustness of the ranking") +
  theme_np + theme(strip.text = element_text(size = 6.5))
lg3 <- theme(legend.position = "bottom", legend.direction = "vertical")
pE3 <- retitle(pE, "B  Stem specificity per cohort") +
  scale_colour_manual(values = setNames(ds_col, ds_lab), name = NULL) +
  guides(colour = guide_legend(nrow = 2, override.aes = list(size = 1.6)), fill = guide_legend(ncol = 1)) +
  theme(legend.position = "bottom", legend.box = "vertical", legend.spacing.y = unit(0, "pt"), legend.margin = margin())
fig3 <- (retitle(pD, "A  Where each target acts") + pE3 +
           (p3C / p3D / p3E + plot_layout(heights = c(1.5, .7, .8))) + plot_layout(widths = c(.55, .9, 1.25)))
save_fig(fig3, "fig3_zeng_stem", 7.2, 8.6)

############################################################################
# Figure 4 -- malignant identities against their normal counterparts
############################################################################
lv4 <- c(mw[target %in% single][order(mal), target], mw[target %in% joint][order(mal), target])
m4 <- mp[, .(ds, target, base, dl)][, `:=`(target = factor(target, lv4), g = grp(target), dsl = factor(ds_lab[ds], ds_lab))]
w4 <- mw[, .(target = factor(target, lv4), g = grp(target), mal, sig = p_mal < .05)]
p4A <- ggplot(m4, aes(dl, target)) +
  geom_vline(xintercept = 0, colour = ink2, linewidth = .3) +
  geom_point(aes(colour = dsl), size = .55, alpha = .55, position = position_jitter(height = .15, width = 0, seed = 1)) +
  geom_point(data = w4, aes(mal, target, fill = sig), shape = 23, size = 1.7, colour = "black", stroke = .3, inherit.aes = FALSE) +
  scale_fill_manual(values = c(`TRUE` = "#0d366b", `FALSE` = "white"), labels = c(`TRUE` = "mean, P<0.05", `FALSE` = "mean, n.s."), name = NULL) +
  scale_colour_manual(values = setNames(ds_col, ds_lab), name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 1.6, alpha = 1), nrow = 2), fill = guide_legend(nrow = 2)) +
  facet_grid(g ~ ., scales = "free_y", space = "free_y") +
  coord_cartesian(xlim = quantile(m4$dl, c(.005, .995))) +
  labs(x = "malignant minus normal (target-specific effect)", y = NULL, title = "A  Per target, 30 matched pairs") +
  theme_np + theme(axis.text.y = element_text(size = 5.6), strip.text = element_blank(), legend.position = "bottom", legend.box = "vertical",
                   legend.spacing.y = unit(0, "pt"), legend.margin = margin())
pp[, `:=`(pair = paste0(ds_lab[ds], " · ", gsub("_", " ", base)), compart = factor(compart, c("Stem", "Progenitor", "Mature")))]
pp[, pair := factor(pair, pp[order(compart, d), rev(pair)])]
p4B <- ggplot(pp, aes(d, pair)) +
  geom_vline(xintercept = 0, colour = ink2, linewidth = .3) +
  geom_segment(aes(x = 0, xend = d, yend = pair), colour = grid, linewidth = .6) +
  geom_point(aes(colour = factor(ds_lab[ds], ds_lab)), size = 1.6) +
  scale_colour_manual(values = setNames(ds_col, ds_lab), guide = "none") +
  facet_grid(compart ~ ., scales = "free_y", space = "free_y") +
  labs(x = "KIT, FLT3, CD33 minus other targets", y = NULL, title = sprintf("B  Blast markers: %d/%d pairs", sum(pp$d > 0), nrow(pp))) +
  theme_np + theme(axis.text.y = element_text(size = 5.6), strip.text.y = element_blank())
sel4 <- c("FLT3", "BCL2", "CD33", "KIT", "IL3RA", "IDH2", "AURKA", "CDK19", "CHEK1")
h4 <- mp[target %in% sel4, .(ds, base, target, dl)][pp[, .(ds, base, pair, compart)], on = .(ds, base)]
h4[, target := factor(target, sel4)]
lim4 <- quantile(abs(h4$dl), .97)
p4C <- ggplot(h4, aes(target, pair, fill = pmax(pmin(dl, lim4), -lim4))) +
  geom_tile(colour = "white", linewidth = .3) +
  scale_fill_gradientn(colours = rev(div_ramp), limits = c(-lim4, lim4), breaks = c(-1, 0, 1), name = "malignant\nminus normal") +
  facet_grid(compart ~ ., scales = "free_y", space = "free_y") +
  labs(x = NULL, y = NULL, title = "C  Per pair") +
  theme_np + theme(panel.grid = element_blank(), axis.text.y = element_blank(), strip.text.y = element_text(angle = -90, size = 6.5),
                   axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5), legend.position = "bottom",
                   legend.title.position = "top", legend.key.width = unit(12, "pt"), legend.key.height = unit(5, "pt"))
fig4 <- p4A + p4B + p4C + plot_layout(widths = c(1, 1, .75))
save_fig(fig4, "fig4_zeng_malignant", 7.2, 7.6)
close(out)
