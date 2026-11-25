
set.seed(2026)
genes3 <- c('PARP1', 'OGG1', 'APEX1')
claimed <- c('PARP2', 'MCM2', 'PBRM1', 'SMARCC1', 'RBM14', 'TARDBP', 'HNRNPR', 'YEATS2')

for (p in c('survival', 'dplyr', 'boot', 'data.table'))
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
suppressPackageStartupMessages({library(survival); library(dplyr); library(boot)})
dir.create('outputs_ber', showWarnings = FALSE); dir.create('outputs_ber/figures', showWarnings = FALSE)
writeLines(capture.output(sessionInfo()), 'outputs_ber/sessionInfo.txt')
cat('Working folder:', getwd(), '\n')
if (!any(grepl('^data_clinical_patient', list.files())))
  stop('data_clinical_patient is not in this folder. Use setwd() to your data folder first.')

find_file <- function(base) {
  all <- list.files()
  hit <- all[all == base | tools::file_path_sans_ext(all) == base]
  if (length(hit) != 1) stop('Need exactly one file for: ', base, ' | found: ', paste(hit, collapse = ', '))
  hit
}
fmt_p <- function(p) ifelse(p < 0.001, 'P < 0.001', sprintf('P = %.3f', p))

# 1. DATA
clin <- read.delim(find_file('data_clinical_patient'), comment.char = '#', stringsAsFactors = FALSE)
clin2 <- data.frame(patient_id = clin$PATIENT_ID,
                    os_time = suppressWarnings(as.numeric(clin$OS_MONTHS)),
                    os_event = as.integer(grepl('^1', clin$OS_STATUS)),
                    age = if ('AGE' %in% names(clin)) suppressWarnings(as.numeric(clin$AGE)) else NA_real_,
                    stringsAsFactors = FALSE)
clin2$os_event[is.na(clin2$os_time)] <- NA

ex <- read.delim(find_file('data_mrna_seq_v2_rsem'), check.names = FALSE, stringsAsFactors = FALSE)
ex <- ex[!is.na(ex$Hugo_Symbol) & ex$Hugo_Symbol != '', ]
ex <- ex[!duplicated(ex$Hugo_Symbol), ]
mat <- as.matrix(ex[, !(names(ex) %in% c('Hugo_Symbol', 'Entrez_Gene_Id'))])
rownames(mat) <- ex$Hugo_Symbol
mode(mat) <- 'numeric'
mat <- log2(mat + 1)
mat <- mat[rowSums(is.na(mat)) == 0, ]
mat <- mat[, grepl('-01$', substr(colnames(mat), 1, 15))]
colnames(mat) <- substr(colnames(mat), 1, 12)
mat <- mat[, !duplicated(colnames(mat))]
stopifnot(all(genes3 %in% rownames(mat)))
cat('Tumours with RNA-seq:', ncol(mat), '\n')

hf <- find_file('TCGA.HRD_withSampleID')
h <- if (grepl('\\.xlsx?$', hf)) as.data.frame(readxl::read_excel(hf)) else data.table::fread(hf, data.table = FALSE)
hrd_col <- if ('HRD' %in% names(h)) 'HRD' else if ('HRD-sum' %in% names(h)) 'HRD-sum' else stop('No HRD column found')
is_bc <- sapply(h, function(x) mean(grepl('^TCGA-', as.character(x))) > 0.5)
bc <- as.character(h[[names(h)[is_bc][1]]])
keep <- nchar(bc) < 15 | substr(bc, 14, 15) == '01'
hrd <- data.frame(patient_id = substr(bc, 1, 12), hrd = as.numeric(h[[hrd_col]]))[keep, ]
hrd <- hrd[!duplicated(hrd$patient_id), ]

d0 <- data.frame(patient_id = colnames(mat), t(mat[genes3, ]), stringsAsFactors = FALSE)
d0 <- d0 %>% left_join(clin2, by = 'patient_id') %>% left_join(hrd, by = 'patient_id')
d  <- d0[!is.na(d0$os_time) & !is.na(d0$os_event), ]
d_age <- d[!is.na(d$age), ]
flow <- data.frame(step = c('Tumours with RNA-seq', 'With survival data', 'With survival data and age', 'With genomic HRD score'),
                   n = c(nrow(d0), nrow(d), nrow(d_age), sum(!is.na(d0$hrd))),
                   deaths = c(NA, sum(d$os_event), sum(d_age$os_event), NA))
cat('\n=== B1. COHORT ===\n'); print(flow)
write.csv(flow, 'outputs_ber/B1_cohort_flow.csv', row.names = FALSE)

# 2. EXPRESSION SUMMARY
summ <- do.call(rbind, lapply(genes3, function(g) {
  x <- d0[[g]]; q <- quantile(x, c(.25, .5, .75))
  data.frame(gene = g, n = length(x), median = q[2], q1 = q[1], q3 = q[3], min = min(x), max = max(x))
}))
rownames(summ) <- NULL
cat('\n=== B2. EXPRESSION, log2(RSEM + 1), TUMOURS ===\n'); print(round(summ[, -1], 2))
write.csv(summ, 'outputs_ber/B2_expression_summary.csv', row.names = FALSE)
png('outputs_ber/figures/FigB1_expression_distribution.png', width = 2400, height = 900, res = 300)
par(mfrow = c(1, 3), mar = c(4.5, 4.5, 3, 1))
for (g in genes3) hist(d0[[g]], breaks = 25, col = '#1b9e77', border = 'white', main = g,
                       xlab = 'log2(RSEM + 1)', ylab = 'Tumours', cex.main = 1)
dev.off()

# 3. PAIRWISE CO-EXPRESSION
rho_ci <- function(x, y) {
  ct <- suppressWarnings(cor.test(x, y, method = 'spearman', exact = FALSE))
  b <- boot(data.frame(x, y), function(z, i) cor(z$x[i], z$y[i], method = 'spearman'), R = 2000)
  ci <- boot.ci(b, type = 'perc')$percent[4:5]
  c(rho = unname(ct$estimate), p = ct$p.value, lo = ci[1], hi = ci[2])
}
pairs_df <- data.frame(t(sapply(list(c('PARP1', 'OGG1'), c('PARP1', 'APEX1'), c('OGG1', 'APEX1')),
                                function(pr) rho_ci(d0[[pr[1]]], d0[[pr[2]]]))))
pairs_df <- cbind(pair = c('PARP1-OGG1', 'PARP1-APEX1', 'OGG1-APEX1'), n = nrow(d0), pairs_df)
cat('\n=== B3. PAIRWISE SPEARMAN CORRELATION (tumours) ===\n'); print(pairs_df, digits = 3)
write.csv(pairs_df, 'outputs_ber/B3_pairwise_correlation.csv', row.names = FALSE)
png('outputs_ber/figures/FigB2_pairwise_scatter.png', width = 2400, height = 2400, res = 300)
pairs(d0[, genes3], pch = 19, col = adjustcolor('grey30', 0.5), cex = 0.6,
      panel = function(x, y, ...) { points(x, y, ...); lines(lowess(x, y), col = '#d95f02', lwd = 2) },
      labels = paste0(genes3, '\nlog2(RSEM + 1)'))
dev.off()

# 4. CORRELATION WITH GENOMIC HRD
dh <- d0[!is.na(d0$hrd), ]
hrd_df <- data.frame(t(sapply(genes3, function(g) rho_ci(dh[[g]], dh$hrd))))
hrd_df <- cbind(gene = genes3, n = nrow(dh), hrd_df); rownames(hrd_df) <- NULL
hrd_df$p_BH <- p.adjust(hrd_df$p, method = 'BH')
cat('\n=== B4. SPEARMAN CORRELATION WITH GENOMIC HRD SCORE ===\n'); print(hrd_df, digits = 3)
write.csv(hrd_df, 'outputs_ber/B4_correlation_with_HRD.csv', row.names = FALSE)

# 5. SURVIVAL
for (g in genes3) { d[[paste0('z_', g)]] <- as.numeric(scale(d[[g]])) }
d_age <- d[!is.na(d$age), ]
surv_rows <- lapply(genes3, function(g) {
  f1 <- coxph(as.formula(paste0('Surv(os_time, os_event) ~ z_', g)), data = d)
  f2 <- coxph(as.formula(paste0('Surv(os_time, os_event) ~ z_', g, ' + age')), data = d_age)
  s1 <- summary(f1); s2 <- summary(f2)
  data.frame(gene = g,
             n_unadj = f1$n, HR_unadj = s1$conf.int[1, 1], lo_unadj = s1$conf.int[1, 3], hi_unadj = s1$conf.int[1, 4], p_unadj = s1$coefficients[1, 5],
             n_age = f2$n, HR_age = s2$conf.int[1, 1], lo_age = s2$conf.int[1, 3], hi_age = s2$conf.int[1, 4], p_age = s2$coefficients[1, 5],
             PH_global_p = cox.zph(f2)$table['GLOBAL', 'p'])
})
surv_df <- do.call(rbind, surv_rows)
surv_df$p_unadj_BH <- p.adjust(surv_df$p_unadj, method = 'BH')
surv_df$p_age_BH <- p.adjust(surv_df$p_age, method = 'BH')
cat('\n=== B5. COX MODELS, HR PER 1 SD HIGHER log2 EXPRESSION ===\n'); print(surv_df, digits = 3)
write.csv(surv_df, 'outputs_ber/B5_survival_cox.csv', row.names = FALSE)

joint <- coxph(Surv(os_time, os_event) ~ z_PARP1 + z_OGG1 + z_APEX1 + age, data = d_age)
sj <- summary(joint)
joint_df <- data.frame(term = rownames(sj$conf.int), HR = sj$conf.int[, 1], lo = sj$conf.int[, 3], hi = sj$conf.int[, 4], p = sj$coefficients[, 5]); rownames(joint_df) <- NULL
cat('\n=== B6. JOINT MODEL (all three genes + age; n =', joint$n, ', events =', joint$nevent, ') ===\n'); print(joint_df, digits = 3)
write.csv(joint_df, 'outputs_ber/B6_joint_model.csv', row.names = FALSE)

km_plot <- function(file, data, grp, title, cols = c('#d95f02', '#1b9e77'), tmax = 120, step = 20) {
  data$grp <- factor(grp); lv <- levels(data$grp)
  fit <- survfit(Surv(os_time, os_event) ~ grp, data = data)
  sdf <- survdiff(Surv(os_time, os_event) ~ grp, data = data)
  p <- 1 - pchisq(sdf$chisq, df = length(sdf$n) - 1)
  tt <- seq(0, tmax, by = step)
  nr <- matrix(summary(fit, times = tt, extend = TRUE)$n.risk, nrow = length(lv), byrow = TRUE)
  png(file, width = 2000, height = 1900, res = 300)
  par(mar = c(9, 4.5, 3, 1))
  plot(fit, col = cols[seq_along(lv)], lwd = 2.5, mark.time = FALSE, xlim = c(0, tmax), xlab = '',
       ylab = 'Overall survival probability', main = title, cex.main = 0.95)
  mtext('Months from diagnosis', side = 1, line = 2.4)
  legend('topright', bty = 'n', lwd = 2.5, col = cols[seq_along(lv)], legend = paste0(lv, ' (n = ', table(data$grp)[lv], ')'), cex = 0.8)
  legend('bottomleft', bty = 'n', legend = paste('Log-rank', fmt_p(p)), cex = 0.9)
  mtext('Number at risk', side = 1, line = 3.8, at = 0, adj = 0, cex = 0.7, font = 2)
  for (i in seq_along(lv)) for (j in seq_along(tt)) mtext(nr[i, j], side = 1, line = 4.6 + 0.8 * (i - 1), at = tt[j], cex = 0.7, col = cols[i])
  dev.off()
  p
}
ms <- data.frame(gene = genes3, n = nrow(d), HR_high_vs_low = NA, lo = NA, hi = NA, p_cox = NA, p_logrank = NA)
for (k in seq_along(genes3)) {
  g <- genes3[k]
  grp <- factor(ifelse(d[[g]] > median(d[[g]]), 'High expression', 'Low expression'), levels = c('Low expression', 'High expression'))
  dd <- d; dd$grp <- grp
  f <- coxph(Surv(os_time, os_event) ~ grp, data = dd); s <- summary(f)
  ms[k, c('HR_high_vs_low', 'lo', 'hi', 'p_cox')] <- c(s$conf.int[1, 1], s$conf.int[1, 3], s$conf.int[1, 4], s$coefficients[1, 5])
  ms$p_logrank[k] <- km_plot(sprintf('outputs_ber/figures/FigB%d_KM_%s.png', 2 + k, g), d, grp, paste('Overall survival by', g, 'expression (median split)'))
}
cat('\n=== B7. MEDIAN-SPLIT SURVIVAL (high vs low) ===\n'); print(ms, digits = 3)
write.csv(ms, 'outputs_ber/B7_median_split.csv', row.names = FALSE)

rows <- rbind(
  data.frame(label = sprintf('%s, per 1 SD (age-adjusted, n=%d)', surv_df$gene, surv_df$n_age), HR = surv_df$HR_age, lo = surv_df$lo_age, hi = surv_df$hi_age),
  data.frame(label = sprintf('%s, joint model (n=%d)', genes3, joint$n), HR = joint_df$HR[1:3], lo = joint_df$lo[1:3], hi = joint_df$hi[1:3]))
k <- nrow(rows); xl <- c(min(0.5, min(rows$lo) * 0.9), max(2, max(rows$hi) * 1.1))
png('outputs_ber/figures/FigB6_forest_plot.png', width = 2900, height = 1600, res = 300)
par(mar = c(5.5, 17, 2.5, 9))
plot(NA, xlim = xl, ylim = c(0.5, k + 0.5), log = 'x', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', main = 'Overall survival: hazard ratios per 1 SD higher expression', cex.main = 0.95)
at <- c(0.5, 0.75, 1, 1.5, 2); axis(1, at = at[at >= xl[1] & at <= xl[2]], labels = at[at >= xl[1] & at <= xl[2]])
mtext('Hazard ratio (95% CI), log scale', side = 1, line = 2.6)
mtext('<-- lower hazard of death  |  higher hazard of death -->', side = 1, line = 4, cex = 0.7)
abline(v = 1, lty = 2, col = 'grey40'); y <- k:1
segments(rows$lo, y, rows$hi, y, lwd = 2.5, col = '#1b9e77'); points(rows$HR, y, pch = 15, cex = 1.3, col = '#1b9e77')
axis(2, at = y, labels = rows$label, las = 1, tick = FALSE, cex.axis = 0.75)
xr <- 10^par('usr')[2]
text(xr, y, sprintf('%.2f (%.2f-%.2f)', rows$HR, rows$lo, rows$hi), adj = 0, xpd = NA, cex = 0.75)
text(xr, k + 0.5, 'HR (95% CI)', adj = 0, xpd = NA, cex = 0.75, font = 2)
dev.off()
write.csv(rows, 'outputs_ber/B8_forest_values.csv', row.names = FALSE)

# 6. GENOME-WIDE CO-EXPRESSION
expressed <- apply(mat, 1, function(x) sd(x) > 0 && mean(x > 0) >= 0.5)
mat2 <- mat[expressed, ]
cat('\nGenes with variable expression in >= 50% of tumours:', nrow(mat2), '\n')
writeLines(rownames(mat2), 'outputs_ber/background_genes_expressed.txt')
cm <- cor(t(mat2), t(mat2[genes3, , drop = FALSE]), method = 'spearman')
top100 <- list()
for (g in genes3) {
  r <- cm[, g]; r <- r[names(r) != g]; r <- sort(r, decreasing = TRUE)
  top100[[g]] <- r[1:100]
  write.csv(data.frame(rank = 1:100, gene = names(top100[[g]]), spearman_rho = as.numeric(top100[[g]])),
            sprintf('outputs_ber/B9_top100_%s.csv', g), row.names = FALSE)
}
ov <- data.frame(comparison = c('PARP1 & OGG1', 'PARP1 & APEX1', 'OGG1 & APEX1', 'all three'),
                 n_shared = c(length(intersect(names(top100$PARP1), names(top100$OGG1))),
                              length(intersect(names(top100$PARP1), names(top100$APEX1))),
                              length(intersect(names(top100$OGG1), names(top100$APEX1))),
                              length(Reduce(intersect, lapply(top100, names)))))
ov$genes <- c(paste(intersect(names(top100$PARP1), names(top100$OGG1)), collapse = ' '),
              paste(intersect(names(top100$PARP1), names(top100$APEX1)), collapse = ' '),
              paste(intersect(names(top100$OGG1), names(top100$APEX1)), collapse = ' '),
              paste(Reduce(intersect, lapply(top100, names)), collapse = ' '))
cat('\n=== B9. OVERLAP OF THE TOP-100 CO-EXPRESSED GENES (Spearman) ===\n'); print(ov, right = FALSE)
cat('\nRange of Spearman rho within each top-100 list:\n')
print(sapply(top100, function(x) round(range(x), 2)))
write.csv(ov, 'outputs_ber/B9_top100_overlap.csv', row.names = FALSE)

cl <- intersect(claimed, rownames(mat2))
cat('\nClaimed genes missing from the expressed set:', paste(setdiff(claimed, cl), collapse = ' '), '\n')
if (length(cl)) {
  pear <- cor(t(mat2[cl, , drop = FALSE]), t(mat2[genes3, , drop = FALSE]))
  chk <- do.call(rbind, lapply(genes3, function(g) {
    r <- cm[, g]; r <- r[names(r) != g]
    rk <- rank(-r)
    data.frame(seed = g, claimed_gene = cl, spearman_rho = round(r[cl], 3), rank_of_all = as.integer(rk[cl]),
               in_top100 = rk[cl] <= 100, pearson_r = round(pear[cl, g], 3))
  }))
  rownames(chk) <- NULL
  cat('\n=== B10. CHECK OF THE EIGHT GENES NAMED IN THE OLD REPORT ===\n'); print(chk)
  write.csv(chk, 'outputs_ber/B10_claimed_gene_check.csv', row.names = FALSE)
}
cat('\nBER SCRIPT DONE. Results are in the outputs_ber folder.\n')

