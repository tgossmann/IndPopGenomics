###############################################
###  Folded Four-Lineage SFS Demonstration  ###
###  + Effect of Mosaic Haploid Reference   ###
###  + Loss of hidden-ref 1/4 configurations ###
###############################################

pdf("SFS_mosaic_effect.pdf", width = 10, height = 6)

# total haplotypes = 4 (four lineages)
n <- 4

# --- SFS setup (allele counts 1,2,3 among 4 lineages) ---
# approximate neutral SFS shape ~ 1/i
weights <- 1 / (1:(n - 1))                 # (1, 1/2, 1/3)
weights <- weights / sum(weights)

S_initial <- 4000
sfs_neutral <- weights * S_initial         # counts for i = 1,2,3

# negative Tajima's D: excess rare (i=1 and i=3 after unfolding)
sfs_neg <- sfs_neutral
sfs_neg[c(1, 3)] <- sfs_neg[c(1, 3)] * 2.5
sfs_neg <- sfs_neg * (sum(sfs_neutral) / sum(sfs_neg))

# positive Tajima’s D: excess intermediate (i=2)
sfs_pos <- sfs_neutral
sfs_pos[2] <- sfs_pos[2] * 2.5
sfs_pos <- sfs_pos * (sum(sfs_neutral) / sum(sfs_pos))

# --- fold into 2 bins (MAF = 1/4 and 2/4) ---
fold_sfs <- function(sfs) {
  maf1 <- sfs[1] + sfs[3]  # minor allele count = 1 (or 3) out of 4
  maf2 <- sfs[2]           # minor allele count = 2 out of 4
  c(`1/4` = maf1, `2/4` = maf2)
}

fold_neutral <- fold_sfs(sfs_neutral)
fold_neg     <- fold_sfs(sfs_neg)
fold_pos     <- fold_sfs(sfs_pos)

###############################################
###  Mosaic haploid effects                  ###
###############################################
# Two conceptual distortions when using a mosaic haploid reference:
# (A) "Collapse": a fraction of true 2/4 sites (shared polymorphism / ref het)
#     is reassigned to the 1/4 class because the reference is treated as fixed.
#
# (B) "Loss": a fraction of true 1/4 sites is entirely missed because the minor
#     allele is on the *hidden* reference haplotype and becomes indistinguishable
#     from the mosaic reference base (appears monomorphic).
#
# We represent (B) explicitly as a third bar "lost".

mosaic_effect <- function(folded_sfs, f_collapse = 0.5, f_loss1 = 0.2) {
  stopifnot(all(c("1/4", "2/4") %in% names(folded_sfs)))
  
  # (B) loss from the original 1/4 bin (hidden reference singleton configs)
  lost1 <- folded_sfs["1/4"] * f_loss1
  one_obs_after_loss <- folded_sfs["1/4"] - lost1
  
  # (A) collapse: shift some of 2/4 into observed 1/4
  shift <- folded_sfs["2/4"] * f_collapse
  
  c(`1/4` = one_obs_after_loss + shift,
    `2/4` = folded_sfs["2/4"] - shift,
    `lost` = lost1)
}

# Tunable parameters
f_collapse <- 0.5  # fraction of 2/4 reassigned to 1/4 in mosaic mapping
f_loss1    <- 0.2  # fraction of 1/4 that is completely missed (hidden ref haplotype)

fold_neutral_mosaic <- mosaic_effect(fold_neutral, f_collapse, f_loss1)
fold_neg_mosaic     <- mosaic_effect(fold_neg,     f_collapse, f_loss1)
fold_pos_mosaic     <- mosaic_effect(fold_pos,     f_collapse, f_loss1)

###############################################
###  Plotting: 2 x 3 grid                   ###
###############################################

# Colors (use the same two colors across all panels; add a faint "lost" color)
col_12   <- c("grey80", "grey40")  # 1/4, 2/4
col_lost <- "grey95"              # lost

# Global y-limit so panels are comparable (include "lost" bars)
ymax <- max(
  fold_neutral["1/4"] + fold_neutral["2/4"],
  fold_neg["1/4"]     + fold_neg["2/4"],
  fold_pos["1/4"]     + fold_pos["2/4"],
  sum(fold_neutral_mosaic),
  sum(fold_neg_mosaic),
  sum(fold_pos_mosaic)
) * 1.10

par(mfrow = c(2, 3),
    mar = c(4.2, 4.2, 4.8, 1.0),
    cex.lab = 1.2,
    cex.axis = 1.1)

plot_panel <- function(vals, panel_title, show_lost = FALSE) {
  if (!show_lost) {
    barplot(vals[c("1/4", "2/4")],
            col = col_12,
            ylim = c(0, ymax),
            xlab = "Minor allele frequency",
            ylab = "Number of sites",
            names.arg = c("1/4", "2/4"),
            las = 1)
  } else {
    barplot(vals[c("1/4", "2/4", "lost")],
            col = c(col_12, col_lost),
            ylim = c(0, ymax),
            xlab = "Minor allele frequency",
            ylab = "Number of sites",
            names.arg = c("1/4", "2/4", "lost"),
            las = 1)
  }
  
  # Force stable top labels (avoids barplot main rendering glitches)
  mtext(panel_title, side = 3, line = 1.2, cex = 1.2)
  
  # Optional: annotate retained vs lost totals for mosaic panels
  if (show_lost) {
    total_true <- sum(vals)
    total_obs  <- vals["1/4"] + vals["2/4"]
    lost_pct   <- 100 * (vals["lost"] / total_true)
    txt <- sprintf("observed = %.0f, lost = %.0f (%.1f%%)",
                   total_obs, vals["lost"], lost_pct)
    mtext(txt, side = 3, line = 0.1, cex = 0.9)
  }
}

# Row 1: true four-lineage folded SFS (2 bars)
plot_panel(fold_neutral, "(a) Neutral (D = 0)", show_lost = FALSE)
plot_panel(fold_neg,     "(b) Negative (D < 0)", show_lost = FALSE)
plot_panel(fold_pos,     "(c) Positive (D > 0)", show_lost = FALSE)

# Row 2: mosaic mapping outcome (2 bars + explicit "lost")
plot_panel(fold_neutral_mosaic, "(d) Neutral with mosaic reference", show_lost = TRUE)
plot_panel(fold_neg_mosaic,     "(e) Negative (D < 0) with mosaic",  show_lost = TRUE)
plot_panel(fold_pos_mosaic,     "(f) Positive (D > 0) with mosaic",  show_lost = TRUE)

# Legend once (in last panel area; tweak position if needed)
legend("topright",
       legend = c("1/4", "2/4", "lost (hidden ref 1/4)"),
       fill   = c(col_12, col_lost),
       bty    = "n",
       cex    = 0.95)

dev.off()