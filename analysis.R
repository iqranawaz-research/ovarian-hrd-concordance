set.seed(2026)

genes <- c('PARP1','OGG1','APEX1','PARP2','BRCA1','BRCA2','RAD51','RAD51B',
           'RAD51C','RAD51D','PALB2','FANCA','FANCD2','ATM','ATR','CHEK1','CHEK2')

need <- c('survival','dplyr','boot','data.table','readxl','rstudioapi')
for (p in need) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')
if (!requireNamespace('GSVA', quietly = TRUE)) BiocManager::install('GSVA', ask = FALSE, update = FALSE)
suppressPackageStartupMessages({library(survival); library(dplyr); library(boot); library(GSVA)})

if (requireNamespace('rstudioapi', quietly = TRUE) && rstudioapi::isAvailable()) {
  p <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) '')
  if (!is.null(p) && nzchar(p)) setwd(dirname(p))
}
cat('Working folder:', getwd(), '\n')
if (!any(grepl('^data_clinical_patient', list.files())))
  stop('data_clinical_patient is not in this folder. Save analysis.R in the same folder as your data files, then run again.')

dir.create('outputs', showWarnings = FALSE)
writeLines(capture.output(sessionInfo()), 'outputs/sessionInfo.txt')

find_file <- function(base) {
  all <- list.files()
  hit <- all[all == base | tools::file_path_sans_ext(all) == base]
  if (length(hit) != 1) stop('Need exactly one file for: ', base, ' | found: ', paste(hit, collapse = ', '))
  hit
}

# ---- 1. CLINICAL ----
clin <- read.delim(find_file('data_clinical_patient'), comment.char = '#', stringsAsFactors = FALSE)
if (!all(c('PATIENT_ID','OS_MONTHS','OS_STATUS') %in% names(clin)))
  stop('Clinical columns not as expected. Found: ', paste(names(clin), collapse = ', '))
age_col   <- if ('AGE' %in% names(clin)) 'AGE' else NA
stage_col <- grep('PATHOLOGIC.*STAGE', names(clin), value = TRUE)[1]
cat('Age column:', age_col, '| Stage column:', stage_col, '\n')
st <- if (is.na(stage_col)) rep(NA_character_, nrow(clin)) else toupper(as.character(clin[[stage_col]]))
clin2 <- data.frame(
  patient_id = clin$PATIENT_ID,
  os_time  = as.numeric(clin$OS_MONTHS),
  os_event = as.integer(grepl('^1', clin$OS_STATUS)),
  age = if (is.na(age_col)) rep(NA_real_, nrow(clin)) else as.numeric(clin[[age_col]]),
  stage_grp = ifelse(grepl('^STAGE IV', st), 'IV',
                     ifelse(grepl('^STAGE III', st), 'III',
                            ifelse(grepl('^STAGE I', st), 'I-II', NA))),
  stringsAsFactors = FALSE)

# ---- 2. EXPRESSION + SIGNATURE SCORE ----
ex <- read.delim(find_file('data_mrna_seq_v2_rsem'), check.names = FALSE, stringsAsFactors = FALSE)
ex <- ex[!is.na(ex$Hugo_Symbol) & ex$Hugo_Symbol != '', ]
ex <- ex[!duplicated(ex$Hugo_Symbol), ]
mat <- as.matrix(ex[, !(names(ex) %in% c('Hugo_Symbol','Entrez_Gene_Id'))])
rownames(mat) <- ex$Hugo_Symbol
mode(mat) <- 'numeric'
mat <- log2(mat + 1)
mat <- mat[rowSums(is.na(mat)) == 0, ]
mat <- mat[, grepl('-01$', substr(colnames(mat), 1, 15))]
colnames(mat) <- substr(colnames(mat), 1, 12)
mat <- mat[, !duplicated(colnames(mat))]
cat('Expression samples:', ncol(mat), '\n')
miss <- setdiff(genes, rownames(mat))
if (length(miss)) cat('NOT FOUND in expression file:', paste(miss, collapse = ', '), '\n')
genes_used <- intersect(genes, rownames(mat))
cat('Genes used:', length(genes_used), 'of', length(genes), '(your report says 15 of 17)\n')
if (length(genes_used) < 2) stop('Too few signature genes found in the expression file.')
sc <- if (packageVersion('GSVA') >= '1.50.0') {
  gsva(ssgseaParam(mat, list(HRBER = genes_used)))
} else {
  gsva(mat, list(HRBER = genes_used), method = 'ssgsea')
}
hrber <- data.frame(patient_id = colnames(sc), hrber = as.numeric(sc['HRBER', ]))

# ---- 3. HRD SCORES ----
hrd_score_col <- Sys.getenv('HRD_COL', unset = NA)
hf <- find_file('TCGA.HRD_withSampleID')
h <- if (grepl('\\.xlsx?$', hf)) as.data.frame(readxl::read_excel(hf)) else data.table::fread(hf, data.table = FALSE)
cat('\nColumns in the HRD file:\n'); print(names(h)); print(head(h, 3))
is_bc <- sapply(h, function(x) mean(grepl('^TCGA-', as.character(x))) > 0.5)
if (!any(is_bc)) stop('No column with TCGA barcodes found in the HRD file.')
bc_col <- names(h)[is_bc][1]
num_cols <- names(h)[sapply(h, is.numeric)]
if (is.na(hrd_score_col)) {
  if (length(num_cols) == 1) {
    hrd_score_col <- num_cols
  } else if (interactive()) {
    cat('\nSummary of the numeric columns:\n'); print(summary(h[num_cols]))
    pick <- utils::menu(num_cols, title = 'Type the NUMBER of the column that is the combined HRD score used in your report (then press Enter):')
    if (pick == 0) stop('No column chosen.')
    hrd_score_col <- num_cols[pick]
  } else {
    stop('Several numeric columns: ', paste(num_cols, collapse = ', '))
  }
}
cat('HRD barcode column:', bc_col, '| HRD score column:', hrd_score_col, '\n')
bc <- as.character(h[[bc_col]])
keep <- nchar(bc) < 15 | substr(bc, 14, 15) == '01'
hrd <- data.frame(patient_id = substr(bc, 1, 12), hrd = as.numeric(h[[hrd_score_col]]))[keep, ]
if (any(duplicated(hrd$patient_id))) cat('Note: duplicate patients in HRD file; keeping the first.\n')
hrd <- hrd[!duplicated(hrd$patient_id), ]

# ---- 4. MERGE + PATIENT FLOW ----
df <- clin2 %>% left_join(hrber, by = 'patient_id') %>% left_join(hrd, by = 'patient_id')
flow <- data.frame(
  step = c('Patients in clinical file','With HRD score','With expression score','With both'),
  n = c(nrow(df), sum(!is.na(df$hrd)), sum(!is.na(df$hrber)), sum(!is.na(df$hrd) & !is.na(df$hrber))))
cat('\n=== PATIENT FLOW ===\n'); print(flow)
write.csv(flow, 'outputs/01_patient_flow.csv', row.names = FALSE)
d <- df %>% filter(!is.na(hrd), !is.na(hrber), !is.na(os_time), !is.na(os_event))
cat('ANALYSIS n =', nrow(d), '| deaths =', sum(d$os_event), '\n')

# ---- 5. CORRELATION + 95% CI ----
ct <- cor.test(d$hrd, d$hrber, method = 'spearman', exact = FALSE)
bt <- boot(d, function(x, i) cor(x$hrd[i], x$hrber[i], method = 'spearman'), R = 2000)
ci <- boot.ci(bt, type = 'perc')$percent[4:5]
res_corr <- data.frame(n = nrow(d), rho = unname(ct$estimate), p = ct$p.value, ci_low = ci[1], ci_high = ci[2])
cat('\n=== CORRELATION ===\n'); print(res_corr)
write.csv(res_corr, 'outputs/02_correlation.csv', row.names = FALSE)

# ---- 6. GROUPS AT 2, 3, 4 CUT-OFFS ----
make_groups <- function(d, k) {
  hh <- ntile(d$hrd, k); ss <- ntile(d$hrber, k)
  factor(ifelse(hh == 1 & ss == k, 'LowHRD_HighExpr',
                ifelse(hh == k & ss == 1, 'HighHRD_LowExpr', 'Other')),
         levels = c('Other','LowHRD_HighExpr','HighHRD_LowExpr'))
}
out <- list()
for (k in c(2, 3, 4)) {
  d$g <- make_groups(d, k)
  sm <- summary(coxph(Surv(os_time, os_event) ~ g, data = d))
  nm <- sub('^g', '', rownames(sm$conf.int))
  out[[as.character(k)]] <- data.frame(cutoff_groups = k, group = nm, n = as.integer(table(d$g)[nm]),
                                       HR = sm$conf.int[, 'exp(coef)'], lo = sm$conf.int[, 'lower .95'], hi = sm$conf.int[, 'upper .95'],
                                       p = sm$coefficients[, 'Pr(>|z|)'])
}
res_groups <- do.call(rbind, out); rownames(res_groups) <- NULL
cat('\n=== GROUPS BY CUT-OFF (reference group = Other) ===\n'); print(res_groups)
write.csv(res_groups, 'outputs/03_groups_by_cutoff.csv', row.names = FALSE)

# ---- 7. CONTINUOUS MODEL (NO CUT-OFFS) ----
d$hrd_z  <- as.numeric(scale(d$hrd))
d$expr_z <- as.numeric(scale(d$hrber))
m_cont <- coxph(Surv(os_time, os_event) ~ hrd_z * expr_z, data = d)
cat('\n=== CONTINUOUS MODEL ===\n'); print(summary(m_cont)); print(cox.zph(m_cont))
capture.output(summary(m_cont), cox.zph(m_cont), file = 'outputs/04_continuous_model.txt')

# ---- 8. ADJUSTED FOR AGE AND STAGE ----
d$stage_grp <- factor(d$stage_grp)
covs <- c('age','stage_grp')[sapply(c('age','stage_grp'), function(v) any(!is.na(d[[v]])))]
rhs  <- paste(c('hrd_z * expr_z', covs), collapse = ' + ')
m_adj <- coxph(as.formula(paste('Surv(os_time, os_event) ~', rhs)), data = d)
cat('\n=== ADJUSTED CONTINUOUS MODEL ===\n'); print(summary(m_adj))
cat('Adjusted model n =', m_adj$n, '| events =', m_adj$nevent, '\n')
d$g3 <- make_groups(d, 3)
m_g3 <- coxph(as.formula(paste('Surv(os_time, os_event) ~', paste(c('g3', covs), collapse = ' + '))), data = d)
cat('\n=== TERTILE GROUPS, ADJUSTED ===\n'); print(summary(m_g3))
capture.output(summary(m_adj), summary(m_g3), file = 'outputs/05_adjusted_models.txt')
saveRDS(d, 'outputs/analysis_table.rds')

cat('\nDONE. Open the "outputs" folder next to this script. Send me the lines printed above.\n')
# ======================= 9. FOLLOW-UP CHECKS =========================
both <- df[!is.na(df$hrd) & !is.na(df$hrber), ]
ct2 <- cor.test(both$hrd, both$hrber, method = 'spearman', exact = FALSE)
bt2 <- boot(both, function(x, i) cor(x$hrd[i], x$hrber[i], method = 'spearman'), R = 2000)
ci2 <- boot.ci(bt2, type = 'perc')$percent[4:5]
res_corr_all <- data.frame(n = nrow(both), rho = unname(ct2$estimate), p = ct2$p.value, ci_low = ci2[1], ci_high = ci2[2])
cat('\n=== A. CORRELATION, ALL PATIENTS WITH BOTH SCORES ===\n'); print(res_corr_all)
write.csv(res_corr_all, 'outputs/06_correlation_all_with_both_scores.csv', row.names = FALSE)

cat('\n=== B. PATIENT WITH BOTH SCORES BUT NO SURVIVAL DATA ===\n')
print(both[is.na(both$os_time) | is.na(both$os_event), c('patient_id','hrd','hrber','os_time','os_event')])
cat('\nOS_STATUS values:\n'); print(table(clin$OS_STATUS, useNA = 'ifany'))
cat('\nStage values (all NA means stage was NOT adjusted for):\n')
if (!is.na(stage_col)) print(table(clin[[stage_col]], useNA = 'ifany'))

cat('\n=== C. SAME PATIENTS FOR EVERY MODEL (age-complete) ===\n')
d_age <- d[!is.na(d$age), ]
cat('n =', nrow(d_age), '| deaths =', sum(d_age$os_event), '\n')
m_g_un  <- coxph(Surv(os_time, os_event) ~ g3, data = d_age)
m_g_age <- coxph(Surv(os_time, os_event) ~ g3 + age, data = d_age)
m_c_age <- coxph(Surv(os_time, os_event) ~ hrd_z + expr_z + age, data = d_age)
m_i_age <- coxph(Surv(os_time, os_event) ~ hrd_z * expr_z + age, data = d_age)
cat('\nGroup model, unadjusted:\n'); print(round(summary(m_g_un)$conf.int, 3)); print(round(summary(m_g_un)$coefficients[, 'Pr(>|z|)'], 4))
cat('\nGroup model, adjusted for age:\n'); print(round(summary(m_g_age)$conf.int, 3)); print(round(summary(m_g_age)$coefficients[, 'Pr(>|z|)'], 4))
cat('\nPrimary additive model (HRD + expression + age):\n'); print(round(summary(m_c_age)$conf.int, 3)); print(round(summary(m_c_age)$coefficients[, 'Pr(>|z|)'], 4))
cat('\nProportional hazards, group model:\n'); print(cox.zph(m_g_age))
cat('\nProportional hazards, continuous model + age:\n'); print(cox.zph(m_c_age))
cat('\nInteraction test (likelihood ratio):\n'); print(anova(m_c_age, m_i_age))

cat('\n=== D. DOES EXPRESSION ADD ANYTHING WITHIN LOW-HRD PATIENTS? ===\n')
d$hrd_t <- ntile(d$hrd, 3); d$expr_t <- ntile(d$hrber, 3)
lowH <- d[d$hrd_t == 1, ]
cat('Low-HRD tertile: n =', nrow(lowH), '| deaths =', sum(lowH$os_event), '\n')
print(table(expression_tertile = lowH$expr_t))
fitL <- coxph(Surv(os_time, os_event) ~ expr_z, data = lowH)
cat('Per 1 SD expression, within low-HRD patients:\n'); print(round(summary(fitL)$conf.int, 3)); print(round(summary(fitL)$coefficients[, 'Pr(>|z|)'], 4))
lowH2 <- lowH[lowH$expr_t %in% c(1, 3), ]; lowH2$high_expr <- as.integer(lowH2$expr_t == 3)
fitL2 <- coxph(Surv(os_time, os_event) ~ high_expr, data = lowH2)
cat('High vs low expression tertile, within low-HRD patients (n =', nrow(lowH2), '):\n')
print(round(summary(fitL2)$conf.int, 3)); print(round(summary(fitL2)$coefficients[, 'Pr(>|z|)'], 4))
cat('\nLowest HRD tertile versus everyone else (n =', nrow(d), '):\n')
print(round(summary(coxph(Surv(os_time, os_event) ~ I(hrd_t == 1), data = d))$conf.int, 3))

cat('\n=== E. MULTIPLE-TESTING ADJUSTMENT FOR THE SIX GROUP COMPARISONS ===\n')
print(cbind(res_groups[, c('cutoff_groups','group','HR','p')], BH_adjusted_p = round(p.adjust(res_groups$p, method = 'BH'), 3)))

cat('\n=== F. GRADE SENSITIVITY ===\n')
smp <- read.delim('data_clinical_sample.txt', comment.char = '#', stringsAsFactors = FALSE)
gr <- smp[, c('PATIENT_ID', 'GRADE')]; names(gr) <- c('patient_id', 'grade_raw')
gr <- gr[!duplicated(gr$patient_id), ]
d2 <- merge(d, gr, by = 'patient_id', all.x = TRUE, sort = FALSE)
d2$grade_grp <- ifelse(d2$grade_raw %in% c('G3','G4'), 'G3-4', ifelse(d2$grade_raw %in% c('G1','G2'), 'G1-2', NA))
d2$grade_grp <- factor(d2$grade_grp, levels = c('G3-4','G1-2'))
d2 <- d2[!is.na(d2$age) & !is.na(d2$grade_grp), ]
cat('Grade-complete cohort: n =', nrow(d2), '| deaths =', sum(d2$os_event), '\n'); print(table(d2$grade_grp))
m_ref <- coxph(Surv(os_time, os_event) ~ hrd_z + expr_z + age, data = d2)
m_grd <- coxph(Surv(os_time, os_event) ~ hrd_z + expr_z + age + grade_grp, data = d2)
m_ggr <- coxph(Surv(os_time, os_event) ~ g3 + age + grade_grp, data = d2)
cat('\nContinuous model, same patients, WITHOUT grade:\n'); print(round(summary(m_ref)$conf.int, 3))
cat('\nContinuous model, WITH grade:\n'); print(round(summary(m_grd)$conf.int, 3)); print(round(summary(m_grd)$coefficients[, 'Pr(>|z|)'], 4))
cat('\nGroup model, WITH age and grade:\n'); print(round(summary(m_ggr)$conf.int, 3)); print(round(summary(m_ggr)$coefficients[, 'Pr(>|z|)'], 4))

cat('\nALL DONE. Now save the Console text as outputs/console_log.txt\n')
# ======================= 10. FIGURES (saved as PNG in outputs/figures) ======
dir.create('outputs/figures', showWarnings = FALSE)
pal <- c('#1b9e77', '#d95f02', '#7570b3')
xmax <- 120   # x-axis of survival curves shown up to 120 months

km_plot <- function(file, data, grp, title, cols = pal) {
  data$grp <- factor(grp)
  fit <- survfit(Surv(os_time, os_event) ~ grp, data = data)
  sd  <- survdiff(Surv(os_time, os_event) ~ grp, data = data)
  p   <- 1 - pchisq(sd$chisq, df = length(sd$n) - 1)
  lv  <- levels(data$grp)
  png(file, width = 2000, height = 1600, res = 300)
  par(mar = c(4.5, 4.5, 3, 1))
  plot(fit, col = cols[seq_along(lv)], lwd = 2.5, mark.time = FALSE, xlim = c(0, xmax),
       xlab = 'Months from diagnosis', ylab = 'Overall survival probability', main = title, cex.main = 0.95)
  legend('topright', bty = 'n', lwd = 2.5, col = cols[seq_along(lv)],
         legend = paste0(lv, ' (n = ', table(data$grp)[lv], ')'), cex = 0.8)
  legend('bottomleft', bty = 'n', legend = sprintf('Log-rank P = %.3f', p), cex = 0.9)
  dev.off()
}

# Figure 1: HR/BER signature score vs genomic HRD score
png('outputs/figures/Fig1_signature_vs_HRD.png', width = 2000, height = 1700, res = 300)
par(mar = c(4.5, 4.5, 3, 1))
plot(both$hrd, both$hrber, pch = 19, col = adjustcolor('grey30', 0.6),
     xlab = 'Genomic HRD score', ylab = 'HR/BER signature score (ssGSEA)',
     main = 'HR/BER expression signature vs genomic HRD', cex.main = 0.95)
lines(lowess(both$hrd, both$hrber), col = '#d95f02', lwd = 2.5)
legend('topright', bty = 'n', cex = 0.85,
       legend = sprintf('Spearman rho = %.3f (95%% CI %.2f to %.2f), n = %d',
                        res_corr_all$rho, res_corr_all$ci_low, res_corr_all$ci_high, res_corr_all$n))
dev.off()

# Figure 2: lowest HRD tertile versus all other patients
km_plot('outputs/figures/Fig2_KM_lowest_HRD_tertile.png', d,
        ifelse(d$hrd_t == 1, 'Lowest HRD tertile', 'Other patients'),
        'Overall survival by genomic HRD', cols = c('#d95f02', '#1b9e77'))

# Figure 3: concordance groups (tertiles)
km_plot('outputs/figures/Fig3_KM_concordance_groups.png', d, d$g3,
        'Overall survival by HRD/expression concordance group (tertiles)')

# Figure 4: within low-HRD patients, by expression tertile
lowH$expr_grp <- ifelse(lowH$expr_t == 1, '1 Low expression', ifelse(lowH$expr_t == 2, '2 Middle expression', '3 High expression'))
km_plot('outputs/figures/Fig4_KM_lowHRD_by_expression.png', lowH, lowH$expr_grp,
        'Overall survival by expression tertile within low-HRD patients')

# Figure 5: forest plot of the main hazard ratios
fr <- function(label, m, term) {
  s <- summary(m)$conf.int
  data.frame(label = label, HR = s[term, 'exp(coef)'], lo = s[term, 'lower .95'], hi = s[term, 'upper .95'])
}
m_low <- coxph(Surv(os_time, os_event) ~ I(hrd_t == 1), data = d)
g3row <- res_groups[res_groups$cutoff_groups == 3 & res_groups$group == 'LowHRD_HighExpr', ]
rows <- rbind(
  fr(sprintf('HRD, per 1 SD (age-adjusted, n=%d)', nrow(d_age)), m_c_age, 'hrd_z'),
  fr(sprintf('Expression, per 1 SD (age-adjusted, n=%d)', nrow(d_age)), m_c_age, 'expr_z'),
  fr(sprintf('Lowest HRD tertile vs rest (n=%d)', nrow(d)), m_low, 'I(hrd_t == 1)TRUE'),
  data.frame(label = sprintf('Low-HRD/high-expr vs other (n=%d, unadjusted)', g3row$n), HR = g3row$HR, lo = g3row$lo, hi = g3row$hi),
  fr('Low-HRD/high-expr vs other (age-adjusted)', m_g_age, 'g3LowHRD_HighExpr'),
  fr(sprintf('Expression, per 1 SD, within low-HRD (n=%d)', nrow(lowH)), fitL, 'expr_z'))
k <- nrow(rows)
png('outputs/figures/Fig5_forest_plot.png', width = 2600, height = 1500, res = 300)
par(mar = c(4.5, 19, 2.5, 1))
plot(NA, xlim = c(0.4, 3.2), ylim = c(0.5, k + 0.5), log = 'x', yaxt = 'n',
     xlab = 'Hazard ratio (95% CI), log scale', ylab = '', main = 'Overall survival: main hazard ratios', cex.main = 0.95)
abline(v = 1, lty = 2, col = 'grey40')
y <- k:1
segments(rows$lo, y, rows$hi, y, lwd = 2.5, col = '#1b9e77')
points(rows$HR, y, pch = 15, cex = 1.3, col = '#1b9e77')
axis(2, at = y, labels = rows$label, las = 1, tick = FALSE, cex.axis = 0.75)
text(3.2, y, sprintf('%.2f (%.2f-%.2f)', rows$HR, rows$lo, rows$hi), adj = 1, cex = 0.7, pos = 2, offset = 0.2)
dev.off()
write.csv(rows, 'outputs/08_forest_plot_values.csv', row.names = FALSE)
cat('\nFIGURES DONE. See the outputs/figures folder.\n')