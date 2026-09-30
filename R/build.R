## =============================================================================
## R/build.R -- LOCAL-ONLY PHI pipeline. Run from the project root:
##
##     Rscript R/build.R
##
## Reads the raw RPDR extract and the ICI data pull, derives the cohort, and
## writes three things:
##   cache/patient_level.rds               EMPI-keyed, gitignored, never rendered
##   cache/summary.rds                     aggregate-only list read by the .qmd pages
##   output/sIL2R_ICI_patient_level.xlsx   EMPI-keyed, gitignored, for chart review
##
## Console output is limited to aggregate counts.
## =============================================================================

source("R/functions.R")
if (!file.exists("R/paths_local.R"))
  stop("Create R/paths_local.R from R/paths_template.R (local file paths).")
source("R/paths_local.R")   # defines RPDR_DIR, RPDR_PREFIX, ICI_FILE

dir.create("cache",  showWarnings = FALSE)
dir.create("output", showWarnings = FALSE)
rp <- function(x) file.path(RPDR_DIR, paste0(RPDR_PREFIX, "_", x, ".txt"))
S  <- list()   # aggregate-only summary

## ---- 1. RPDR patient set ------------------------------------------------------
mrn <- read_rpdr(rp("Mrn"), select = c("Enterprise_Master_Patient_Index", "Status"))
setnames(mrn, "Enterprise_Master_Patient_Index", "EMPI")
rpdr_set <- unique(mrn$EMPI)

## ---- 2. RPDR IL2R labs --------------------------------------------------------
lab <- read_rpdr(rp("Lab"), select = c("EMPI", "Seq_Date_Time", "Group_Id", "Test_Id",
                                       "Loinc_Code", "Test_Description", "Result",
                                       "Result_Text", "Reference_Units", "Abnormal_Flag"))
S$lab_rows     <- nrow(lab)
S$lab_patients <- uniqueN(lab$EMPI)

il2 <- lab[Group_Id %in% IL2R_GROUPS]
lab_pts <- unique(lab$EMPI)                       # anyone with any Lab row (IL2R or CRE)
S$lab_groups <- lab[, .(rows = .N, patients = uniqueN(EMPI)), by = Group_Id]
rm(lab); invisible(gc())
il2[, draw_date := parse_rpdr_date(Seq_Date_Time)]
stopifnot(!anyNA(il2$draw_date))
il2[, status := classify_il2r_result(Result, Result_Text)]
il2[, value  := fifelse(status == STATUS_LEVELS[1], extract_num(Result),
                fifelse(status == STATUS_LEVELS[2], extract_num(Result_Text), NA_real_))]
il2[, unit := fifelse(Test_Id %chin% c("5200007493", "954.1945"), "pg/mL", "U/mL")]

S$il2_tests <- il2[, .(records = .N, patients = uniqueN(EMPI),
                       first_year = year(min(draw_date)), last_year = year(max(draw_date))),
                   by = .(Group_Id, Test_Id, Loinc_Code, Test_Description, unit)][order(-records)]
S$il2_status <- il2[, .(records = .N, patients = uniqueN(EMPI)), by = status][order(status)]
S$il2_status_by_test <- dcast(il2[, .N, by = .(Test_Id, status)], Test_Id ~ status,
                              value.var = "N", fill = 0L)
## Short, digit-free status phrases behind the non-numeric classes (lab
## vocabulary, not patient data); anything rarer/longer is only counted.
nn <- il2[status %in% STATUS_LEVELS[3:4]]
nn[, phrase := fifelse(!grepl("[0-9]", Result) & nchar(Result) <= 40,
                       fifelse(nzchar(trimws(Result)), toupper(trimws(Result)), "(blank Result)"),
                       "(other / free text)")]
S$il2_nonnum_phrases <- nn[, .N, by = .(status, phrase)][order(status, -N)]
## numeric values: shape check (anything that is not "number [unit] [flag]")
num <- il2[status == STATUS_LEVELS[1]]
S$il2_numeric_irregular <- sum(!grepl("^\\s*[0-9]+(\\.[0-9]+)?(\\s+[A-Za-z/]+)?(\\s+[HL])?\\s*$", num$Result))
S$il2_value_summary <- il2[!is.na(value), .(n = .N, median = median(value),
                                            p25 = quantile(value, .25), p75 = quantile(value, .75)),
                           by = unit]
S$il2_by_year <- il2[, .(records = .N, patients = uniqueN(EMPI),
                         valid_records = sum(status %in% VALID_STATUS)),
                     by = .(year = year(draw_date))][order(year)]
S$il2_same_day_multi <- il2[status %in% VALID_STATUS, .N, by = .(EMPI, draw_date)][N > 1, .N]

## ---- 3. ICI data pull ---------------------------------------------------------
ici <- read_ici(ICI_FILE)                       # EMPI MRN -> EMPI; date-only start
ici[, ici_date := SCHEDULEDSTARTDTS]
ici[, blinded  := is_blinded(MEDICATIONDSC)]
ici[, generic  := ici_generic(MEDICATIONDSC)]

S$ici_rows          <- nrow(ici)
S$ici_patients      <- uniqueN(ici$EMPI)
S$ici_cols          <- names(ici)
S$ici_empi_9digit   <- mean(grepl("^[0-9]{9}$", ici$EMPI))
S$rpdr_empi_9digit  <- mean(grepl("^[0-9]{9}$", il2$EMPI))
S$ici_date_na       <- sum(is.na(ici$ici_date))
S$ici_date_range    <- range(ici$ici_date, na.rm = TRUE)
S$ici_hour          <- ici[, .N, by = .(hour = SCHEDULEDSTARTDTS_raw_hour)][order(hour)]
S$ici_by_year       <- ici[, .(records = .N, patients = uniqueN(EMPI)),
                           by = .(year = year(ici_date))][order(year)]
S$ici_mar_action    <- ici[, .N, by = MARACTIONDSC][order(-N)]
S$ici_order_status  <- ici[, .N, by = ORDERSTATUSDSC][order(-N)]
S$ici_route         <- ici[, .N, by = ROUTEDSC][order(-N)]
S$ici_drug_desc_n   <- uniqueN(ici$MEDICATIONDSC)
S$ici_unmatched_desc <- ici[is.na(generic), .N]
S$ici_generic       <- ici[, .(records = .N, patients = uniqueN(EMPI)), by = generic][order(-records)]
S$ici_blinded       <- ici[, .(records = .N, patients = uniqueN(EMPI),
                               descriptions = uniqueN(MEDICATIONDSC)), by = blinded]
S$ici_blinded_only_patients <- ici[, .(all_blinded = all(blinded)), by = EMPI][all_blinded == TRUE, .N]

## ---- 4. Linkage + cohort -----------------------------------------------------
ici_empi <- unique(ici$EMPI)
il2_pts  <- unique(il2$EMPI)
valid_pts <- unique(il2[status %in% VALID_STATUS, EMPI])
post_pts  <- unique(il2[status %in% VALID_STATUS & draw_date >= CUTOFF_DATE, EMPI])

S$flow <- data.table(
  step = c("RPDR patient set (IL2R lab OR Epic IL2R procedure code)",
           "  with >=1 IL2R lab record in the Lab file",
           "  with >=1 performed IL2R result (not cancelled/credited)",
           "  with >=1 performed IL2R result on/after 2015-01-01"),
  n = c(length(rpdr_set), length(il2_pts), length(valid_pts), length(post_pts)))
S$excluded <- data.table(
  reason = c("In RPDR set, but no IL2R row in the Lab file (no draw date; see Section 3)",
             "IL2R records all cancelled / credited",
             "All performed IL2R values before 2015-01-01"),
  n = c(sum(!rpdr_set %chin% il2_pts),
        length(setdiff(il2_pts, valid_pts)),
        length(setdiff(valid_pts, post_pts))))
S$no_lab_in_ici <- sum(!rpdr_set %chin% il2_pts & rpdr_set %chin% ici_empi)
S$il2_not_in_set <- sum(!il2_pts %chin% rpdr_set)
S$overlap_any_il2 <- sum(il2_pts %chin% ici_empi)

## primary definition
pts <- derive_cohort(il2, ici)
S$answers <- answer_counts(pts)
S$groups  <- pts[, .N, by = group][order(group)]
S$groups_by_first_year <- dcast(pts[, .N, by = .(year = year(first_post_draw), group)],
                                year ~ group, value.var = "N", fill = 0L)
S$same_day_A <- pts[group == GROUP_LEVELS[1] & same_day_draw == TRUE, .N]
S$same_day_B <- pts[group == GROUP_LEVELS[2] & same_day_draw == TRUE, .N]

A <- pts[group == GROUP_LEVELS[1]]
S$A_n_after <- A[, .N, by = .(n_after = fifelse(n_after_dates >= 5L, "5+", as.character(n_after_dates)))][order(n_after)]
S$A_days_to_first <- A[, .(n = .N, median = median(days_ici_to_first_after),
                           p25 = quantile(days_ici_to_first_after, .25),
                           p75 = quantile(days_ici_to_first_after, .75))]
S$A_time_bins <- A[, .N, by = .(bin = cut(days_ici_to_first_after,
                                          c(0, 90, 180, 365, 730, Inf),
                                          labels = c("1-90 days", "91-180 days", "181-365 days",
                                                     "1-2 years", "> 2 years")))][order(bin)]
S$A_first_regimen <- A[, .N, by = ici_first_regimen][order(-N)]
S$A_ici_start_year <- A[, .N, by = .(year = year(ici_start))][order(year)]
S$A_measurements_total <- A[, sum(n_after_dates)]

## ---- 5. Sensitivity grid -------------------------------------------------------
sens <- list(
  "Primary"                                        = list(),
  "Numeric results only"                           = list(valid_status = NUMERIC_ONLY),
  "Include blinded or-placebo trial drug"          = list(include_blinded = TRUE),
  "Same-day ICI counts as prior"                   = list(same_day_is_prior = TRUE),
  "Cutoff 2016-01-01 (after calendar 2015)"        = list(cutoff = as.Date("2016-01-01")),
  "Every IL2R record (incl. cancelled)"            = list(valid_status = STATUS_LEVELS))
S$sensitivity <- rbindlist(lapply(names(sens), function(nm) {
  p <- do.call(derive_cohort, c(list(il2 = il2, ici = ici), sens[[nm]]))
  as.data.table(c(list(analysis = nm), as.list(answer_counts(p))))
}))

## ---- 6. RPDR medication cross-check (ICI outside the ICI pull) -----------------
med <- read_rpdr(rp("Med"), select = c("EMPI", "Medication_Date", "Medication"))
S$med_rows <- nrow(med)
med <- med[is_ici_med(Medication)]
med[, med_date := parse_rpdr_date(Medication_Date)]
med[, generic  := ici_generic(Medication)]
S$med_ici_rows <- nrow(med)
S$med_ici_patients <- uniqueN(med$EMPI)
rpdr_ici <- med[!is.na(med_date), .(rpdr_ici_first = min(med_date),
                                    rpdr_ici_n     = .N,
                                    rpdr_ici_drugs = paste(sort(unique(generic)), collapse = "; ")),
                by = EMPI]
pts <- merge(pts, rpdr_ici, by = "EMPI", all.x = TRUE)
pts[, rpdr_ici_any := !is.na(rpdr_ici_first)]
pts[, rpdr_ici_before_last_draw := rpdr_ici_any & rpdr_ici_first < last_post_draw]
pts[, rpdr_ici_pre2015 := rpdr_ici_any & rpdr_ici_first < CUTOFF_DATE]
S$med_by_group <- pts[, .(patients = .N,
                          rpdr_ici_any = sum(rpdr_ici_any),
                          rpdr_ici_before_last_draw = sum(rpdr_ici_before_last_draw),
                          rpdr_ici_first_pre2015 = sum(rpdr_ici_pre2015)), by = group][order(group)]
S$no_lab_rpdr_ici <- sum(!rpdr_set %chin% il2_pts & rpdr_set %chin% rpdr_ici$EMPI)
rm(med); invisible(gc())

## ---- 6b. RPDR diagnosis cross-check: malignancy codes (chart-review triage) ----
dia <- read_rpdr(rp("Dia"), select = c("EMPI", "Date", "Code_Type", "Code"))
dia[, ca := (Code_Type == "ICD10" & grepl(CANCER_ICD10, Code)) |
            (Code_Type == "ICD9"  & grepl(CANCER_ICD9,  Code))]
dia[, hx := (Code_Type == "ICD10" & grepl(CANCER_HX_ICD10, Code)) |
            (Code_Type == "ICD9"  & grepl(CANCER_HX_ICD9,  Code))]
dia <- dia[ca | hx]
dia[, dx_date := parse_rpdr_date(Date)]
ca <- dia[ca == TRUE & !is.na(dx_date), .(cancer_dx_first = min(dx_date)), by = EMPI]
hx <- dia[hx == TRUE, .(cancer_hx_any = TRUE), by = EMPI]
rm(dia); invisible(gc())
pts <- merge(pts, ca, by = "EMPI", all.x = TRUE)
pts <- merge(pts, hx, by = "EMPI", all.x = TRUE)
pts[, cancer_dx_any := !is.na(cancer_dx_first)]
pts[, cancer_dx_before_last_draw := cancer_dx_any & cancer_dx_first <= last_post_draw]
pts[is.na(cancer_hx_any), cancer_hx_any := FALSE]
S$dx_by_group <- pts[, .(patients = .N,
                         cancer_dx_any = sum(cancer_dx_any),
                         cancer_dx_before_last_draw = sum(cancer_dx_before_last_draw)),
                     by = group][order(group)]

## chart-review triage for group C (hierarchical: first tier that applies)
pts[, review_tier := NA_character_]
pts[group == GROUP_LEVELS[4], review_tier := fcase(
  rpdr_ici_any,               "1. ICI in RPDR Med file",
  cancer_dx_before_last_draw, "2. Malignancy code on/before last IL2R draw",
  cancer_dx_any,              "3. Malignancy code only after last IL2R draw",
  cancer_hx_any,              "4. Cancer-history / chemotherapy code only",
  default =                   "5. No malignancy, history or chemotherapy code")]
S$C_triage <- pts[group == GROUP_LEVELS[4], .N, by = .(tier = review_tier)][order(tier)]

## Group B with an RPDR-Med ICI before an IL2R draw: possibly Q1-positive
S$B_rpdr_earlier <- pts[group == GROUP_LEVELS[2] & rpdr_ici_before_last_draw == TRUE, .N]

## ---- 6c. The RPDR patients with an IL2R order but no result row ----------------
nl <- data.table(EMPI = setdiff(rpdr_set, il2_pts))
nl[, in_ici_file := EMPI %chin% ici_empi]
nl <- merge(nl, ici[blinded == FALSE & !is.na(ici_date), .(ici_start = min(ici_date)), by = EMPI],
            by = "EMPI", all.x = TRUE)
nl <- merge(nl, rpdr_ici, by = "EMPI", all.x = TRUE)
nl <- merge(nl, ca, by = "EMPI", all.x = TRUE)
nl <- merge(nl, hx, by = "EMPI", all.x = TRUE)
nl[, `:=`(rpdr_ici_any = !is.na(rpdr_ici_first), cancer_dx_any = !is.na(cancer_dx_first),
          cancer_hx_any = !is.na(cancer_hx_any) & cancer_hx_any)]
S$nolab <- nl[, .(patients = .N, rpdr_ici_any = sum(rpdr_ici_any),
                  cancer_dx_any = sum(cancer_dx_any),
                  cancer_or_hx = sum(cancer_dx_any | cancer_hx_any)),
              by = in_ici_file][order(-in_ici_file)]
S$no_lab_cancer <- nl[, sum(cancer_dx_any)]
S$nolab_other_lab <- nl[, sum(EMPI %chin% lab_pts)]   # Lab rows, but creatinine only
S$nolab_no_lab    <- nl[, sum(!EMPI %chin% lab_pts)]  # no Lab rows at all

S$il2_last_month <- format(max(il2$draw_date), "%Y-%m")
S$il2_first_month_epic <- format(min(il2[Test_Id == "5200007493", draw_date]), "%Y-%m")

## ---- 7. Save -----------------------------------------------------------------
S$built <- Sys.time()
S$cutoff <- CUTOFF_DATE
S$files <- list(rpdr = paste0(RPDR_PREFIX, "_{Mrn,Lab,Med}.txt"), ici = basename(ICI_FILE))
saveRDS(S, "cache/summary.rds")
saveRDS(list(pts = pts, il2 = il2[, .(EMPI, draw_date, Test_Id, status, value, unit)]),
        "cache/patient_level.rds")

## local chart-review workbook (EMPI-keyed; gitignored; never published)
fmt_cols <- c("EMPI", "group", "review_tier", "n_post_dates", "first_post_draw",
              "last_post_draw", "in_ici_file", "ici_start", "ici_first_regimen",
              "n_after_dates", "first_after", "days_ici_to_first_after", "same_day_draw",
              "rpdr_ici_any", "rpdr_ici_first", "rpdr_ici_drugs", "rpdr_ici_before_last_draw",
              "cancer_dx_any", "cancer_dx_first", "cancer_dx_before_last_draw", "cancer_hx_any")
out <- pts[, ..fmt_cols][order(group, review_tier, first_post_draw)]
out[, group := as.character(group)]
nl_out <- nl[, .(EMPI, in_ici_file, ici_start, rpdr_ici_any, rpdr_ici_first, rpdr_ici_drugs,
                 cancer_dx_any, cancer_dx_first, cancer_hx_any)][
  order(-in_ici_file, -rpdr_ici_any, -cancer_dx_any)]
writexl::write_xlsx(list(
  README = data.frame(note = c(
    "LOCAL ONLY - contains EMPI. Do not email, upload or commit.",
    paste("Built", format(S$built, "%Y-%m-%d %H:%M"), "by R/build.R"),
    "Cohort: >=1 performed soluble IL2R result on/after 2015-01-01 (RPDR Lab).",
    "ICI start = earliest SCHEDULEDSTARTDTS (date) in ICI_data_pull, excluding blinded or-placebo products.",
    "ici_first_regimen = every ICI given on the ICI start date.",
    "rpdr_ici_* = ICI mentions in the RPDR Med file (independent of the ICI pull).",
    "cancer_dx_* = ICD-10 C00-C96 (incl. C4A/C7A/C7B) / ICD-9 140-209.3 in the RPDR Dia file.",
    "cancer_hx_any = ICD-10 Z85, Z51.1x, Z92.21 / ICD-9 V10, V58.1x.",
    "review_tier (Group C only): 1 = ICI in RPDR Med ... 5 = no cancer-related code.",
    "no_IL2R_result_row = RPDR patients who qualified via the Epic IL2R order code but have no Lab row (no draw date).")),
  not_in_ICI_cohort          = out[group == GROUP_LEVELS[4]],
  ICI_cohort_RPDR_ICI_earlier = out[group == GROUP_LEVELS[2] & rpdr_ici_before_last_draw == TRUE],
  ICI_before_IL2R            = out[group == GROUP_LEVELS[1]],
  all_post2015_IL2R          = out,
  no_IL2R_result_row         = nl_out),
  "output/sIL2R_ICI_patient_level.xlsx")

cat("\n== aggregate results ==\n")
print(S$flow); print(S$excluded)
print(S$answers)
print(S$groups)
print(S$sensitivity)
print(S$med_by_group)
print(S$dx_by_group)
print(S$C_triage)
print(S$nolab)
cat("Group B with RPDR-Med ICI before last draw:", S$B_rpdr_earlier, "\n")
print(S$A_first_regimen)
cat("Built OK:", format(S$built), "\n")
