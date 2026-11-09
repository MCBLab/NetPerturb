#!/usr/bin/env Rscript
# Figure 5 -- does the in-silico score agree with the Mead et al. 2022 screen?
#
# The screen (supplementary FIG 1b) scores 434 compounds at four doses on
# Paneth/LYZ abundance: repSSMD.LYZ.NS and .S, positive meaning more Paneth
# cells, which is the hit direction (KPT-330 is the paper's best validated
# hit, +2.69 in FIG 1d). A compound's screen score here is its highest
# repSSMD across doses, since the strong hits act at low dose and turn
# negative at toxic ones (KPT-330: +1.02 at 0.08 uM, -1.71 at 10 uM with
# ATP -2.99).
#
# Compounds are mapped to target genes through ChEMBL mechanism-of-action
# records (mead/chembl_map.py -> mead/compound_targets_chembl.tsv), human
# symbols taken to mouse by capitalisation, then matched against the 206
# scored targets. A compound with several targets takes the mean of their
# scores, and a target hit by several compounds keeps the best screen score.
#
# Panels
#   A  screen score vs. stem II residual, and vs. stem III residual
#   B  the same for every identity: Spearman of the screen score against the
#      identity's residual, so any agreement can be placed on the lineage
#   C  the compounds that map, ranked by screen score, with their targets'
#      stem II and stem III residuals
#
# Run from publication/analysis:  Rscript fig5_screen.R
# Writes ../Fig/fig5_screen.{pdf,png} and screen_stats.txt

suppressMessages({ library(data.table); library(ggplot2); library(patchwork); library(ggrepel); library(readxl) })
out <- file("screen_stats_tmp.txt", "w"); say <- function(...) { cat(..., "\n", sep = ""); cat(..., "\n", sep = "", file = out) }

col_a <- "#2a78d6"; col_b <- "#eb6834"; ink <- "#0b0b0b"; ink2 <- "#52514e"; grid <- "#e6e5e1"
theme_np <- theme_minimal(base_size = 7.5) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = grid, linewidth = .3),
        axis.text = element_text(colour = ink2), plot.title = element_text(face = "bold", size = 8),
        plot.title.position = "plot", legend.key.size = unit(8, "pt"),
        plot.subtitle = element_text(size = 7, colour = ink2))

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

# ---- the pipeline's scores ------------------------------------------------
d <- fread("../data/mead2022/perbscore_all_targets.txt")
d[, lg := log10(perb_score)]
d[, floor := perb_score == 0 | perb_score < 1e-6 * median(perb_score[perb_score > 0])]
nf <- d[floor == FALSE]
nf[, resid := resid(lm(lg ~ target + cell_type, data = nf))]
res <- dcast(nf, target ~ cell_type, value.var = "resid")

# ---- the screen -----------------------------------------------------------
xl <- "../data/mead2022/41551_2022_863_MOESM3_ESM.xlsx"
b <- as.data.table(read_excel(xl, sheet = "FIG 1b"))
setnames(b, c("Treatment", "Dose", "FC.ATP", "FC.LYZ.NS", "FC.LYZ.S", "Z.ATP", "Z.LYZ.NS", "Z.LYZ.S", "SSMD.ATP", "SSMD.LYZ.NS", "SSMD.LYZ.S"))
for (cc in names(b)[-1]) set(b, j = cc, value = suppressWarnings(as.numeric(b[[cc]])))
screen <- b[!is.na(SSMD.LYZ.NS), .(ssmd_ns = max(SSMD.LYZ.NS), ssmd_s = max(SSMD.LYZ.S, na.rm = TRUE),
                                   dose_at_best = Dose[which.max(SSMD.LYZ.NS)]), by = Treatment]
say("screen: ", nrow(screen), " compounds; KPT-330 best repSSMD.LYZ.NS ", sprintf("%.2f", screen[Treatment == "KPT-330", ssmd_ns]),
    " at ", screen[Treatment == "KPT-330", dose_at_best], " uM, rank ", screen[order(-ssmd_ns), which(Treatment == "KPT-330")], " of ", nrow(screen))

# ---- compound -> target ---------------------------------------------------
ct <- fread("mead/compound_targets_chembl.tsv")
ct <- ct[organism == "Homo sapiens" & human_genes != ""]
ct <- ct[, .(gene = unlist(strsplit(human_genes, ";"))), by = .(compound, mechanism_of_action)]
ct[, mouse := paste0(substr(gene, 1, 1), tolower(substr(gene, 2, nchar(gene))))]
scored <- unique(d$target)
single <- ct[mouse %in% scored, .(compound, target = mouse)]
# a joint target counts for a compound when the compound hits any of its
# genes; this run scored single genes only, so the branch is usually empty
joint <- scored[grepl(";", scored)]
map <- if (length(joint) > 0) {
  jm <- rbindlist(lapply(joint, function(j) data.table(target = j, mouse = strsplit(j, ";")[[1]])))
  unique(rbind(single, merge(ct[, .(compound, mouse)], jm, by = "mouse", allow.cartesian = TRUE)[, .(compound, target)]))
} else unique(single)
say("ChEMBL: ", uniqueN(ct$compound), " compounds with a human mechanism target; ", uniqueN(map$compound),
    " of them hit at least one of the ", length(scored), " scored targets (", uniqueN(map$target), " targets covered)")

cm <- merge(map, screen[, .(compound = Treatment, ssmd_ns, ssmd_s)], by = "compound")
cm <- merge(cm, res, by = "target")
agg <- cm[, c(.(n_targets = .N, targets = paste(sort(unique(target)), collapse = ";"), ssmd_ns = ssmd_ns[1], ssmd_s = ssmd_s[1]),
              lapply(.SD, mean, na.rm = TRUE)), by = compound, .SDcols = ids]
say("compounds in the correlation: ", nrow(agg), " (those with a screen score and >= 1 scored target)")

# ---- correlations ---------------------------------------------------------
cor_tab <- rbindlist(lapply(ids, function(id) {
  x <- agg[[id]]; ok <- !is.na(x)
  data.table(identity = id,
             rho_ns = cor(agg$ssmd_ns[ok], x[ok], method = "spearman"),
             p_ns = cor.test(agg$ssmd_ns[ok], x[ok], method = "spearman", exact = FALSE)$p.value,
             rho_s = cor(agg$ssmd_s[ok], x[ok], method = "spearman"),
             n = sum(ok))
}))
setorder(cor_tab, -rho_ns)
for (i in seq_len(nrow(cor_tab))) say(sprintf("Spearman(screen LYZ.NS, %s residual) = %+.2f (p = %.3f, n = %d); LYZ.S %+.2f",
                                              id_lab[cor_tab$identity[i]], cor_tab$rho_ns[i], cor_tab$p_ns[i], cor_tab$n[i], cor_tab$rho_s[i]))
say("KPT-330 in the mapped set: targets ", agg[compound == "KPT-330", targets], "; screen ", sprintf("%.2f", agg[compound == "KPT-330", ssmd_ns]),
    " (rank ", agg[order(-ssmd_ns), which(compound == "KPT-330")], " of ", nrow(agg), "); stem II ", sprintf("%+.2f", agg[compound == "KPT-330", stem_II]),
    " (rank ", agg[order(-stem_II), which(compound == "KPT-330")], "), stem III ", sprintf("%+.2f", agg[compound == "KPT-330", stem_III]),
    " (rank ", agg[order(-stem_III), which(compound == "KPT-330")], ")")
top5 <- agg[order(-ssmd_ns)][1:5]
say("top 5 by screen score: ", paste(sprintf("%s (%.2f; stem II %+.2f, III %+.2f)", top5$compound, top5$ssmd_ns, top5$stem_II, top5$stem_III), collapse = "; "))

# ---- A: screen vs stem II and stem III ------------------------------------
la <- melt(agg, id.vars = c("compound", "ssmd_ns", "targets"), measure.vars = c("stem_II", "stem_III"), variable.name = "identity", value.name = "resid")[!is.na(resid)]
la[, identity := factor(id_lab[as.character(identity)], levels = id_lab[c("stem_II", "stem_III")])]
la[, hit := compound == "KPT-330"]
lab_a <- la[hit | compound %in% agg[order(-ssmd_ns)][1:5, compound]]
ann <- cor_tab[identity %in% c("stem_II", "stem_III")]
ann[, identity := factor(id_lab[identity], levels = id_lab[c("stem_II", "stem_III")])]
pA <- ggplot(la, aes(resid, ssmd_ns)) +
  geom_hline(yintercept = 0, colour = grid, linewidth = .3) + geom_vline(xintercept = 0, colour = grid, linewidth = .3) +
  geom_smooth(method = "lm", se = FALSE, colour = "#c3c2b7", linewidth = .5, formula = y ~ x) +
  geom_point(aes(colour = hit, size = hit)) +
  geom_text_repel(data = lab_a, aes(label = compound, colour = hit), size = 2, max.overlaps = 20, seed = 1,
                  segment.colour = "#c3c2b7", show.legend = FALSE) +
  geom_text(data = ann, aes(x = Inf, y = -Inf, label = sprintf("rho = %+.2f (p = %.2f, n = %d)", rho_ns, p_ns, n)),
            hjust = 1.05, vjust = -0.8, size = 2.2, colour = ink2, inherit.aes = FALSE) +
  facet_wrap(~ identity, nrow = 1) +
  scale_colour_manual(values = c(`FALSE` = "#8f8e8a", `TRUE` = col_a), guide = "none") +
  scale_size_manual(values = c(`FALSE` = 1.1, `TRUE` = 2.6), guide = "none") +
  labs(x = "Target's residual log10 score in that identity", y = "Screen: best repSSMD (LYZ.NS)",
       title = "A  Screen effect on Paneth cells vs. in-silico score in the cycling stem pools",
       subtitle = "one point per compound; KPT-330 (XPO1) in blue; positive repSSMD = more Paneth cells") +
  theme_np + theme(strip.text = element_text(face = "bold", colour = ink2))

# ---- B: correlation across identities -------------------------------------
cb <- copy(cor_tab)
cb[, identity_f := factor(id_lab[identity], levels = id_lab[ids])]
pB <- ggplot(cb, aes(rho_ns, identity_f)) +
  geom_vline(xintercept = 0, colour = ink2, linewidth = .3) +
  geom_segment(aes(x = 0, xend = rho_ns, yend = identity_f), colour = "#c3c2b7", linewidth = .5) +
  geom_point(aes(colour = p_ns < 0.05), size = 2.4) +
  geom_text(aes(label = sprintf("%+.2f", rho_ns), x = rho_ns + sign(rho_ns) * .035), size = 2.1, colour = ink2) +
  scale_colour_manual(values = c(`FALSE` = "#8f8e8a", `TRUE` = col_b), labels = c("n.s.", "p < 0.05"), name = NULL) +
  scale_x_continuous(limits = c(-.45, .45)) +
  labs(x = "Spearman: screen score vs. residual score", y = NULL,
       title = "B  Agreement by identity",
       subtitle = paste0("n = ", min(cor_tab$n), "-", max(cor_tab$n), " compounds per identity")) +
  theme_np + theme(legend.position = "bottom", panel.grid.major.y = element_blank())

# ---- C: the mapped compounds, ranked --------------------------------------
cc <- melt(agg[order(-ssmd_ns)][1:20], id.vars = c("compound", "ssmd_ns"), measure.vars = c("stem_II", "stem_III"),
           variable.name = "identity", value.name = "resid")[!is.na(resid)]
cc[, `:=`(compound_f = factor(compound, levels = rev(agg[order(-ssmd_ns)][1:20, compound])),
          identity = factor(id_lab[as.character(identity)], levels = id_lab[c("stem_II", "stem_III")]))]
pC <- ggplot(cc, aes(resid, compound_f)) +
  geom_vline(xintercept = 0, colour = ink2, linewidth = .3) +
  geom_line(aes(group = compound_f), colour = "#c3c2b7", linewidth = .4) +
  geom_point(aes(colour = identity), size = 1.9) +
  scale_colour_manual(values = c(col_a, col_b), name = NULL) +
  labs(x = "Residual log10 score", y = NULL,
       title = "C  The 20 strongest screen hits that map to a scored target",
       subtitle = "ordered by screen score, strongest at the top") +
  theme_np + theme(legend.position = "bottom", panel.grid.major.y = element_blank(), axis.text.y = element_text(size = 6))

fig <- pA / (pB | pC) + plot_layout(heights = c(1, 1.15))
ggsave(paste0("../Fig/fig5_screen_", method, ".pdf"), fig, width = 7.1, height = 6.6, device = cairo_pdf)
ggsave(paste0("../Fig/fig5_screen_", method, ".png"), fig, width = 7.1, height = 6.6, dpi = 300, bg = "white")
ggsave(paste0("fig5_screen_", method, "_preview.png"), fig, width = 7.1, height = 6.6, dpi = 150)
fwrite(agg, "mead/compound_scores_joined.tsv", sep = "\t")
close(out)
invisible(file.rename("screen_stats_tmp.txt", stats_file))
