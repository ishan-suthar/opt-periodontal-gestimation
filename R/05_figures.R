#==============================================================================
# 05_figures.R
# Every figure is drawn from the CSVs written by scripts 02-04, so the figures
# always match the tables.
#
#   Figure 1   Adjusted LS-mean trajectories by visit and arm (95% CI)
#   Figure 2   Visit 5 adjusted control-minus-treatment differences (95% CI)
#   Figure 3   Mixed-model change contrast vs IPCW g-estimate (-psi), same scale
#   Figure S1  Participant flow and follow-up status
#   Figure S2  Bootstrap distributions of the IPCW psi, point estimate and BC-a CI
#
# Each figure is saved as PNG (300 dpi) and PDF in output/figures/.
#==============================================================================

COL_CONTROL <- "#7F7F7F"
COL_TREAT   <- "#1B7F6B"
COL_IPCW    <- "#1B7F6B"
COL_PANEL   <- "#404040"

theme_ms <- theme_classic(base_size = 11) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 10),
        legend.position = "bottom", legend.title = element_blank(),
        panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3))

save_fig <- function(p, name, w, h) {
  ggsave(file.path(FIG_DIR, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
  pdf_dev <- if (capabilities("cairo")) cairo_pdf else "pdf"
  ggsave(file.path(FIG_DIR, paste0(name, ".pdf")), p, width = w, height = h, device = pdf_dev)
  cat("Saved", name, "\n")
}

# Site-based outcomes and PLI are on different scales, so they get separate panels.
PANEL_LEVELS <- c("Site-based outcomes (percentage points)", "Plaque index (Silness-L\u00f6e units)")
scale_group  <- function(o) factor(ifelse(o == "PlI", PANEL_LEVELS[2], PANEL_LEVELS[1]),
                                   levels = PANEL_LEVELS)
short_label  <- function(o) factor(o, rev(OUTCOMES), sub(" \\(.*", "", OUTCOME_LABELS[rev(OUTCOMES)]))

panel <- read.csv(out_path("Panel_AllEstimates.csv"), stringsAsFactors = FALSE) %>% filter(Covariates == "Primary")
lsm   <- read.csv(out_path("Panel_LSMeans.csv"),      stringsAsFactors = FALSE) %>% filter(Covariates == "Primary")
gest  <- read.csv(out_path("Gest_Bootstrap.csv"),     stringsAsFactors = FALSE)
# BC-a intervals need many replicates; in TEST_MODE (25 replicates) they can be
# NA, so fall back to percentile intervals there. Full runs always have BC-a.
gest  <- gest %>% mutate(bca_lcl = coalesce(bca_lcl, perc_lcl), bca_ucl = coalesce(bca_ucl, perc_ucl))
reps  <- read.csv(out_path("Gest_BootReplicates.csv"), stringsAsFactors = FALSE)
fu    <- read.csv(out_path("FollowUp_Counts.csv"),    stringsAsFactors = FALSE)

# ---- Figure 1 ----------------------------------------------------------------
f1 <- lsm %>% mutate(
  Outcome = factor(Outcome, OUTCOMES, OUTCOME_LABELS[OUTCOMES]),
  Visit   = factor(Visit, c("Baseline", "Visit3", "Visit5"),
                   c("Baseline\n(13-16 wk)", "Visit 3\n(21-24 wk)", "Visit 5\n(29-32 wk)")),
  Group   = factor(Group, c("C", "T"), c("Control", "Treatment")))
pd <- position_dodge(width = 0.25)
p1 <- ggplot(f1, aes(Visit, emmean, colour = Group, group = Group)) +
  geom_line(position = pd, linewidth = 0.6) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.12, position = pd, linewidth = 0.4) +
  geom_point(position = pd, size = 2) +
  facet_wrap(~ Outcome, scales = "free_y", ncol = 3) +
  scale_colour_manual(values = c(Control = COL_CONTROL, Treatment = COL_TREAT)) +
  labs(x = NULL, y = "Adjusted least-squares mean (95% CI)") +
  theme_ms + theme(panel.grid.major.x = element_blank())
save_fig(p1, "Figure1_LSMean_Trajectories", 9, 6)

# ---- Figure 2 ----------------------------------------------------------------
f2 <- panel %>% transmute(Outcome, est = V5_diff, lcl = V5_lcl, ucl = V5_ucl) %>%
  mutate(Panel = scale_group(Outcome), Outcome = short_label(Outcome))
p2 <- ggplot(f2, aes(est, Outcome)) +
  geom_vline(xintercept = 0, linewidth = 0.4) +
  geom_errorbarh(aes(xmin = lcl, xmax = ucl), height = 0.15, colour = COL_TREAT) +
  geom_point(size = 2.4, colour = COL_TREAT) +
  facet_wrap(~ Panel, scales = "free", ncol = 2) +
  labs(x = "Adjusted Visit 5 difference, control minus treatment (95% CI)", y = NULL) +
  theme_ms
save_fig(p2, "Figure2_V5_Differences", 9, 3.8)

# ---- Figure 3 ----------------------------------------------------------------
m_lab <- c("Mixed model, change Baseline to Visit 5", "IPCW g-estimate (-psi)")
f3 <- bind_rows(
  panel %>% transmute(Outcome, Method = m_lab[1], est = Change_diff, lcl = Change_lcl, ucl = Change_ucl),
  gest %>% filter(Spec == "base") %>%
    transmute(Outcome, Method = m_lab[2], est = -psi, lcl = -bca_ucl, ucl = -bca_lcl)) %>%
  mutate(Panel = scale_group(Outcome), Outcome = short_label(Outcome),
         Method = factor(Method, m_lab))
leg <- c(m_lab[1], expression(paste("IPCW g-estimate (", -Psi, ")")))
pd3 <- position_dodge(width = 0.5)
p3 <- ggplot(f3, aes(est, Outcome, colour = Method, shape = Method)) +
  geom_vline(xintercept = 0, linewidth = 0.4) +
  geom_errorbarh(aes(xmin = lcl, xmax = ucl), height = 0.2, position = pd3) +
  geom_point(size = 2.4, position = pd3) +
  facet_wrap(~ Panel, scales = "free", ncol = 2) +
  scale_colour_manual(values = c(COL_PANEL, COL_IPCW), labels = leg) +
  scale_shape_manual(values = c(16, 17), labels = leg) +
  labs(x = "Reduction with treatment, change from Baseline to Visit 5 (95% CI)", y = NULL) +
  theme_ms
save_fig(p3, "Figure3_Panel_vs_IPCW", 9, 4.2)

# ---- Figure S1: participant flow ---------------------------------------------
fc <- fu[fu$Group == "C", ]; ft <- fu[fu$Group == "T", ]
pct <- function(k, n) sprintf("N = %d (%.1f%%)", k, 100 * k / n)
box <- function(x, y, w, h, title, body, bold = FALSE, fill = "grey96") {
  x <- unname(x)
  list(rect = data.frame(xmin = x - w / 2, xmax = x + w / 2, ymin = y - h / 2, ymax = y + h / 2, fill = fill),
       text = data.frame(x = x, y = c(y + h * 0.2, y - h * 0.2), label = c(title, body),
                         face = c(if (bold) "bold" else "plain", "plain")))
}
xs <- c(C = 3, T = 9); W <- 4.6; H <- 0.62
bx <- list(box(6, 10, 3.2, H, "Randomized", sprintf("N = %d", fc$n + ft$n), TRUE),
           box(xs["C"], 8.6, 2.8, H, "Control",   sprintf("N = %d", fc$n), TRUE),
           box(xs["T"], 8.6, 2.8, H, "Treatment", sprintf("N = %d", ft$n), TRUE))
ys <- c(7.1, 6.1, 5.1, 4.1)
for (a in c("C", "T")) {
  f <- if (a == "C") fc else ft
  lab <- list(c("Attended both follow-up visits",       pct(f$both, f$n)),
              c("Missed Visit 3, returned at Visit 5",  pct(f$missV3_retV5, f$n)),
              c("Attended Visit 3, missed Visit 5",     pct(f$V3_missV5, f$n)),
              c("Missed both follow-up visits",         pct(f$neither, f$n)))
  for (k in 1:4) bx[[length(bx) + 1]] <- box(xs[a] + 0.45, ys[k], W, H, lab[[k]][1], lab[[k]][2])
  bx[[length(bx) + 1]] <- box(xs[a] + 0.45, 2.9, W, H, "Missed at least one follow-up visit",
                              pct(f$incomplete, f$n), fill = "white")
  bx[[length(bx) + 1]] <- box(xs[a] + 0.45, 2.05, W, H, "Examined at Visit 5",
                              sprintf("N = %d", f$examined_V5), fill = "white")
}
rects <- do.call(rbind, lapply(bx, `[[`, "rect")); texts <- do.call(rbind, lapply(bx, `[[`, "text"))
segs <- rbind(
  data.frame(x = 6, xend = 6, y = 10 - H / 2, yend = 9.3),
  data.frame(x = xs["C"], xend = xs["T"], y = 9.3, yend = 9.3),
  data.frame(x = xs["C"], xend = xs["C"], y = 9.3, yend = 8.6 + H / 2),
  data.frame(x = xs["T"], xend = xs["T"], y = 9.3, yend = 8.6 + H / 2),
  do.call(rbind, lapply(xs, function(x) rbind(
    data.frame(x = x, xend = x, y = 8.6 - H / 2, yend = 7.85),
    data.frame(x = x, xend = x - 2.3, y = 7.85, yend = 7.85),
    data.frame(x = x - 2.3, xend = x - 2.3, y = 7.85, yend = ys[4])))))
arrows <- do.call(rbind, lapply(xs, function(x)
  data.frame(x = x - 2.3, xend = x + 0.45 - W / 2, y = ys, yend = ys)))
chi <- fu$chisq[1]; chip <- fu$chisq_p[1]
pS1 <- ggplot() +
  geom_segment(data = segs, aes(x, y, xend = xend, yend = yend), linewidth = 0.4, colour = "grey25") +
  geom_segment(data = arrows, aes(x, y, xend = xend, yend = yend), linewidth = 0.4, colour = "grey25",
               arrow = arrow(length = unit(0.08, "in"), type = "closed")) +
  geom_rect(data = rects, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = I(fill)),
            colour = "grey25", linewidth = 0.4) +
  geom_text(data = texts, aes(x, y, label = label, fontface = face), size = 3.3) +
  geom_hline(yintercept = 3.5, linetype = "dashed", colour = "grey60", linewidth = 0.3) +
  annotate("text", x = 6, y = 1.3, size = 3.3, parse = TRUE,
           label = sprintf('italic("Incomplete follow-up differed by arm: ")*italic(chi)^2*italic("(1) = %.2f, p = %.3f")', chi, chip)) +
  coord_cartesian(xlim = c(0.4, 11.9), ylim = c(1, 10.5), expand = FALSE) +
  theme_void()
save_fig(pS1, "FigureS1_Participant_Flow", 8, 6.09)

# ---- Figure S2 ----------------------------------------------------------------
fs2     <- reps %>% filter(Spec == "base", is.finite(psi)) %>%
  mutate(Outcome = factor(Outcome, OUTCOMES, OUTCOME_LABELS[OUTCOMES]))
fs2_ref <- gest %>% filter(Spec == "base") %>%
  mutate(Outcome = factor(Outcome, OUTCOMES, OUTCOME_LABELS[OUTCOMES]))
pS2 <- ggplot(fs2, aes(psi)) +
  geom_histogram(bins = 40, fill = "grey75", colour = "white", linewidth = 0.2) +
  geom_vline(data = fs2_ref, aes(xintercept = psi), colour = "#B2182B", linewidth = 0.7) +
  geom_vline(data = fs2_ref, aes(xintercept = bca_lcl), colour = COL_TREAT, linetype = "dashed", linewidth = 0.6) +
  geom_vline(data = fs2_ref, aes(xintercept = bca_ucl), colour = COL_TREAT, linetype = "dashed", linewidth = 0.6) +
  facet_wrap(~ Outcome, scales = "free", ncol = 3) +
  labs(x = expression(paste("Bootstrap estimate of ", Psi, " (IPCW)")), y = "Number of resamples") +
  theme_ms
save_fig(pS2, "FigureS2_Bootstrap_Psi", 9, 5.5)
cat("\n")
