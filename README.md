# Ovarian cancer (TCGA-OV): HRD concordance and BER-gene expression
Two exploratory, reproducible R analyses of the TCGA ovarian serous cystadenocarcinoma cohort. TCGA data are not included.
## 1. HRD concordance (this folder)
`analysis.R` tests whether a 15-gene HR/BER expression signature agrees with genomic HRD score and how both relate to overall survival.
Data: cBioPortal study "Ovarian Serous Cystadenocarcinoma (TCGA, PanCancer Atlas)", downloaded 1 September 2026; HRD scores from file TCGA.HRD_withSampleID.txt on the NCI GDC page "The Immune Landscape of Cancer" (https://gdc.cancer.gov/about-data/publications/panimmune; file ID 66dd07d7-6366-4774-83c3-5ad1e22b177e; Thorsson et al., Immunity 2018, doi:10.1016/j.immuni.2018.03.023), downloaded 15 November 2026. Columns ai1, lst1 and hrd-loh are scar components, and HRD is their sum.
Run: place the data files beside the script, then `Sys.setenv(HRD_COL='HRD'); source('analysis.R')`. Seed 2026; R 4.6.1, GSVA 2.6.3 (outputs/sessionInfo.txt).
Notes: use 06_correlation_all_with_both_scores.csv for the correlation (n = 290). 11 patients with blank survival status and 3 with a status but no survival time are excluded. Stage and residual disease are unavailable; models adjust for age, with grade as a sensitivity analysis.
## 2. BER genes: PARP1, OGG1, APEX1 (folder ber_expression)
`ber_script.R` analyses co-expression, correlation with genomic HRD score and overall survival of the three genes in 300 TCGA-OV primary tumours. Same data as above. Run: `source('ber_script.R')` from the data folder; results go to `outputs_ber/`. Tumour-versus-normal plots were made by hand in GEPIA2 (accessed 21 November 2026; OV; |log2FC| cutoff 1, P cutoff 0.01, "Match TCGA normal and GTEx data" = yes, log scale = yes) and are in the folder `gepia2/`.
All findings are exploratory.
