## =============================================================================
## R/functions.R -- helpers for the soluble IL-2 receptor (sIL2R) x ICI pipeline
##
## Pure functions. Nothing in this file prints row-level data; identifier
## columns (EMPI etc.) are used only as join keys inside data.table operations.
## =============================================================================

suppressPackageStartupMessages(library(data.table))

## ---- constants ----------------------------------------------------------------

## RPDR Lab Group_Id values that are soluble IL-2 receptor (CD25) tests.
## (Test_Id: el:5200007493, wc954.1945, mcsq-il2ra, dcsqe-il2ra)
IL2R_GROUPS <- c("INTERL2-RS", "INTERL2", "IL2R-ALPHA")

## "After 2015" = drawn on or after this date. Patients whose IL2R values are
## ALL before this date are excluded.
CUTOFF_DATE <- as.Date("2015-01-01")

## Small-cell threshold for anything rendered to the website (cells 1-10 hidden).
SMALL_CELL <- 11L

## Identifier / quasi-identifier columns: never printed, never rendered.
RESTRICTED_COLS <- c(
  "EMPI", "EMPI MRN", "Enterprise_Master_Patient_Index", "EPIC_PMRN", "PMRN",
  "MRN", "MGH MRN", "DFCI MRN", "MGH_MRN", "BWH_MRN", "FH_MRN", "SRH_MRN",
  "NWH_MRN", "NSMC_MRN", "MCL_MRN", "MEE_MRN", "DFC_MRN", "WDH_MRN",
  "IncomingId", "PATIENTID", "PATIENTENCOUNTERID", "Encounter_number",
  "Accession", "Date_of_Birth", "BIRTHDTS", "Zip_code", "Date_Of_Death"
)

## Immune checkpoint inhibitors: generic names, sponsor codes, US brand names.
ICI_REGEX <- paste(c(
  "pembrolizumab", "nivolumab", "ipilimumab", "atezolizumab", "durvalumab",
  "avelumab", "cemiplimab", "dostarlimab", "tremelimumab", "relatlimab",
  "retifanlimab", "toripalimab", "tislelizumab", "cosibelimab", "ivonescimab",
  "fianlimab", "favezelimab", "vibostolimab",
  "mk-?3475", "bms-?936558", "medi-?4736", "mpdl-?3280a", "tsr-?042",
  "regn-?2810", "msb-?0010718c", "bms-?986213", "mk-?4280a", "mk-?7684a",
  "keytruda", "opdivo", "yervoy", "tecentriq", "imfinzi", "bavencio",
  "libtayo", "jemperli", "imjudo", "opdualag", "zynyz", "loqtorzi",
  "tevimbra", "unloxcyt"), collapse = "|")

## Code / brand -> generic, used only for aggregate drug tables.
ICI_GENERIC_MAP <- c(
  "mk-?3475" = "pembrolizumab", "keytruda" = "pembrolizumab",
  "bms-?936558" = "nivolumab", "opdivo" = "nivolumab",
  "yervoy" = "ipilimumab",
  "mpdl-?3280a" = "atezolizumab", "tecentriq" = "atezolizumab",
  "medi-?4736" = "durvalumab", "imfinzi" = "durvalumab",
  "msb-?0010718c" = "avelumab", "bavencio" = "avelumab",
  "regn-?2810" = "cemiplimab", "libtayo" = "cemiplimab",
  "tsr-?042" = "dostarlimab", "jemperli" = "dostarlimab",
  "imjudo" = "tremelimumab",
  ## fixed-dose co-formulations
  "opdualag" = "nivolumab/relatlimab", "bms-?986213" = "nivolumab/relatlimab",
  "mk-?4280a" = "favezelimab/pembrolizumab", "mk-?7684a" = "pembrolizumab/vibostolimab",
  "zynyz" = "retifanlimab", "loqtorzi" = "toripalimab", "tevimbra" = "tislelizumab",
  "unloxcyt" = "cosibelimab")

## IL2R result strings that mean the test was NOT performed.
CANCEL_REGEX <- paste(c(
  "credit", "refus", "cancel", "no specimen", "no spec\\.", "duplicate",
  "wrong order", "not performed", "not done", "\\bqns\\b",
  "quantity not sufficient", "insufficient", "not received", "reject",
  "unacceptable", "unsatisfactory", "clotted", "hemoly", "specimen lost",
  "lost specimen", "problem tube", "no longer needed", "reordered",
  "discontinued", "incorrect(ly)? (order|collect|label)"), collapse = "|")

## Malignant neoplasm codes (RPDR Dia): ICD-10 C00-C96 incl. C4A/C7A/C7B;
## ICD-9 140-208 plus 209.0-209.3 (malignant neuroendocrine tumours).
CANCER_ICD10 <- "^C(0[0-9]|[1-8][0-9A-Z]|9[0-6])"
CANCER_ICD9  <- "^((1[4-9][0-9]|20[0-8])(\\.|$)|209\\.?[0-3])"
## Cancer history / antineoplastic treatment codes: ICD-10 Z85, Z51.1x (incl.
## Z51.12 immunotherapy encounter), Z92.21; ICD-9 V10, V58.1x.
CANCER_HX_ICD10 <- "^(Z85|Z51\\.?1|Z92\\.?21)"
CANCER_HX_ICD9  <- "^(V10|V58\\.?1)"

## Leading number, optionally with thousands separators ("12,345").
NUM_REGEX <- "^\\s*[<>]?\\s*[0-9]+(,[0-9]{3})*(\\.[0-9]+)?"

STATUS_LEVELS <- c("numeric (Result)", "numeric (Result_Text)",
                   "reported, non-discrete", "cancelled / credited")
VALID_STATUS  <- STATUS_LEVELS[1:3]   # a "measurement" = any non-cancelled result
NUMERIC_ONLY  <- STATUS_LEVELS[1:2]   # sensitivity analysis

GROUP_LEVELS <- c(
  "A. ICI before >=1 post-2015 IL2R draw",
  "B. In ICI cohort, ICI start on/after every post-2015 draw",
  "B2. In ICI cohort, blinded (or-placebo) trial drug only",
  "C. Not in ICI cohort")

## ---- readers ------------------------------------------------------------------

## RPDR text files: pipe-delimited, no quoting, read everything as character.
read_rpdr <- function(path, select = NULL) {
  stopifnot(file.exists(path))
  fread(path, sep = "|", quote = "", colClasses = "character", fill = TRUE,
        select = select, showProgress = FALSE, na.strings = NULL)
}

## ICI data pull. Two changes are applied on read, so the rest of the pipeline
## (and every downstream file) matches the RPDR conventions:
##   1. `EMPI MRN` (column G) is renamed to `EMPI`
##   2. SCHEDULEDSTARTDTS "YYYY-MM-DDThh:mm:ssZ" becomes a Date (YYYY-MM-DD)
read_ici <- function(path) {
  stopifnot(file.exists(path))
  ici <- fread(path, colClasses = "character", showProgress = FALSE)
  stopifnot("EMPI MRN" %in% names(ici), "SCHEDULEDSTARTDTS" %in% names(ici))
  setnames(ici, "EMPI MRN", "EMPI")
  ici[, SCHEDULEDSTARTDTS_raw_hour := as.integer(substr(SCHEDULEDSTARTDTS, 12, 13))]
  ici[, SCHEDULEDSTARTDTS := parse_ici_date(SCHEDULEDSTARTDTS)]
  ici[]
}

## ---- parsers ------------------------------------------------------------------

## RPDR dates: Lab "MM/DD/YYYY HH:MM", Med "M/D/YYYY". as.Date() ignores the
## trailing time component.
parse_rpdr_date <- function(x) as.Date(x, format = "%m/%d/%Y")

## ICI SCHEDULEDSTARTDTS: ISO-8601 "YYYY-MM-DDThh:mm:ssZ" -> calendar Date.
## The date is taken exactly as written (no time-zone shift): the trailing "Z"
## is an export artefact, not a real UTC offset -- scheduled hours cluster at
## 08:00-17:00, i.e. local clinic hours (see data_management.qmd, Table 3).
parse_ici_date <- function(x) {
  ok <- is.na(x) | grepl("^\\d{4}-\\d{2}-\\d{2}T\\d{2}:\\d{2}:\\d{2}Z?$", x)
  if (!all(ok)) stop(sum(!ok), " SCHEDULEDSTARTDTS values are not ISO-8601")
  as.Date(substr(x, 1, 10))
}

extract_num <- function(x) {
  out <- rep(NA_real_, length(x))
  hit <- grepl(NUM_REGEX, x)
  out[hit] <- as.numeric(gsub(",", "", sub(
    "^\\s*[<>]?\\s*([0-9]+(,[0-9]{3})*(\\.[0-9]+)?).*$", "\\1", x[hit])))
  out
}

## Classify every IL2R lab row into one of STATUS_LEVELS.
classify_il2r_result <- function(result, result_text) {
  result      <- fifelse(is.na(result), "", result)
  result_text <- fifelse(is.na(result_text), "", result_text)
  res_num <- grepl(NUM_REGEX, result)
  txt_num <- grepl(NUM_REGEX, result_text)
  cancel  <- !res_num & (grepl(CANCEL_REGEX, result, ignore.case = TRUE) |
                         (!txt_num & grepl(CANCEL_REGEX, result_text, ignore.case = TRUE)))
  factor(fcase(res_num,            STATUS_LEVELS[1],
               cancel,             STATUS_LEVELS[4],
               txt_num,            STATUS_LEVELS[2],
               default =           STATUS_LEVELS[3]),
         levels = STATUS_LEVELS)
}

is_ici_med <- function(x) grepl(ICI_REGEX, x, ignore.case = TRUE)

## Blinded placebo-controlled trial product: the patient may have received
## placebo. "RELATLIMAB OR PLACEBO-NIVOLUMAB" is NOT blinded for ICI exposure
## (nivolumab is given in both arms), nor are "X OR <other ICI>" products.
is_blinded <- function(x) {
  grepl("or placebo", x, ignore.case = TRUE) &
    !grepl("placebo-nivolumab", x, ignore.case = TRUE)
}

## Product label for aggregate tables: the first ICI (or co-formulation) named
## in the description. In "X OR PLACEBO-Y" only Y is certain, so Y is used.
ici_generic <- function(x) {
  x <- sub("^.*placebo-", "", x, ignore.case = TRUE)
  m <- tolower(regmatches(x, regexpr(ICI_REGEX, x, ignore.case = TRUE)))
  out <- rep(NA_character_, length(x))
  out[grepl(ICI_REGEX, x, ignore.case = TRUE)] <- m
  for (k in names(ICI_GENERIC_MAP)) out <- sub(paste0("^", k, "$"), ICI_GENERIC_MAP[[k]], out)
  out
}

## ---- cohort derivation --------------------------------------------------------

## Patient-level classification of the post-cutoff IL2R cohort.
##   il2  : IL2R rows with EMPI, draw_date, status
##   ici  : ICI rows with EMPI, ici_date, blinded
## Returns one row per patient with >=1 valid IL2R draw on/after `cutoff`.
derive_cohort <- function(il2, ici,
                          cutoff            = CUTOFF_DATE,
                          valid_status      = VALID_STATUS,
                          include_blinded   = FALSE,
                          same_day_is_prior = FALSE) {
  post <- il2[status %in% valid_status & draw_date >= cutoff]

  pts <- post[, .(n_post_draws    = .N,
                  n_post_dates    = uniqueN(draw_date),
                  first_post_draw = min(draw_date),
                  last_post_draw  = max(draw_date)), by = EMPI]

  ici_use <- if (include_blinded) ici else ici[blinded == FALSE]
  ici_use <- ici_use[!is.na(ici_date)]
  start   <- ici_use[, .(ici_start = min(ici_date)), by = EMPI]
  ## regimen = every ICI given on the start date (e.g. "ipilimumab + nivolumab")
  reg <- ici_use[start, on = .(EMPI, ici_date = ici_start), nomatch = 0L][
    , .(ici_first_regimen = paste(sort(unique(ici_generic(MEDICATIONDSC))), collapse = " + ")),
    by = EMPI]
  start <- merge(start, reg, by = "EMPI", all.x = TRUE)

  pts[, in_ici_file := EMPI %chin% unique(ici$EMPI)]
  pts <- merge(pts, start, by = "EMPI", all.x = TRUE)

  after <- merge(post[, .(EMPI, draw_date)], start[, .(EMPI, ici_start)], by = "EMPI")
  after <- if (same_day_is_prior) after[draw_date >= ici_start] else after[draw_date > ici_start]
  a <- after[, .(n_after_draws  = .N,
                 n_after_dates  = uniqueN(draw_date),
                 first_after    = min(draw_date)), by = EMPI]
  pts <- merge(pts, a, by = "EMPI", all.x = TRUE)
  pts[is.na(n_after_draws), `:=`(n_after_draws = 0L, n_after_dates = 0L)]

  ## same-day: an IL2R draw on the ICI start date (reported separately)
  sd <- merge(post[, .(EMPI, draw_date)], start[, .(EMPI, ici_start)], by = "EMPI")[
    draw_date == ici_start, unique(EMPI)]
  pts[, same_day_draw := EMPI %chin% sd]

  pts[, group := factor(fcase(
    n_after_dates >= 1L,                GROUP_LEVELS[1],
    in_ici_file & !is.na(ici_start),    GROUP_LEVELS[2],
    in_ici_file &  is.na(ici_start),    GROUP_LEVELS[3],
    default =                           GROUP_LEVELS[4]), levels = GROUP_LEVELS)]
  pts[, days_ici_to_first_after := as.integer(first_after - ici_start)]
  pts[]
}

## The three requested counts from a derive_cohort() frame.
answer_counts <- function(pts) {
  c(post_cutoff_patients    = nrow(pts),
    q1_ici_before_draw      = sum(pts$n_after_dates >= 1L),
    q2_ge2_after_ici        = sum(pts$n_after_dates >= 2L),
    q2_ge2_after_ici_draws  = sum(pts$n_after_draws >= 2L),
    q3_not_in_ici_cohort    = sum(!pts$in_ici_file),
    in_ici_not_before       = sum(pts$in_ici_file & pts$n_after_dates == 0L))
}

## ---- display helpers (website) ------------------------------------------------

## Format a count; 1-10 shown as "<11". Zero stays 0.
fmt_n <- function(n, threshold = SMALL_CELL) {
  n <- as.numeric(n)
  ifelse(is.na(n), "--",
         ifelse(n > 0 & n < threshold, paste0("<", threshold),
                formatC(n, format = "d", big.mark = ",")))
}

## Percent that is hidden whenever its numerator is a small cell.
fmt_pct <- function(n, d, digits = 1, threshold = SMALL_CELL) {
  p <- 100 * as.numeric(n) / as.numeric(d)
  ifelse(is.na(p), "--",
         ifelse(n > 0 & n < threshold, "--", formatC(p, format = "f", digits = digits)))
}

## Percent shown only when the formatted count itself is shown (a percent next
## to a hidden count would reveal it).
pct_shown <- function(n_fmt, n, d, digits = 1) {
  ifelse(grepl("^[0-9]", n_fmt), fmt_pct(n, d, digits), "--")
}

## Change from a reference count; changes of 1-10 in either direction are
## shown as "+<11" / "-<11" so a small difference cannot be read off.
fmt_delta <- function(x, ref, threshold = SMALL_CELL) {
  d <- as.numeric(x) - as.numeric(ref)
  ifelse(d == 0, "same",
         ifelse(abs(d) < threshold, ifelse(d > 0, "+<11", "-<11"),
                paste0(ifelse(d > 0, "+", "-"),
                       formatC(abs(d), format = "d", big.mark = ","))))
}

## Suppress a set of mutually exclusive parts of a SHOWN total. If exactly one
## part is small it could be recovered by subtraction, so the smallest other
## non-zero part is also hidden ("suppressed").
fmt_parts <- function(n, threshold = SMALL_CELL) {
  n <- as.numeric(n)
  small <- n > 0 & n < threshold
  comp  <- rep(FALSE, length(n))
  if (sum(small) == 1L) {
    cand <- which(!small & n > 0)
    if (length(cand)) comp[cand[which.min(n[cand])]] <- TRUE
  }
  out <- fmt_n(n, threshold)
  out[comp] <- "suppressed"
  out
}

## "n of d": hidden when n OR its complement (d - n) is a small cell, since
## the group size d is shown elsewhere.
fmt_n_of <- function(n, d, threshold = SMALL_CELL) {
  n <- as.numeric(n); d <- as.numeric(d)
  hide_n    <- n > 0 & n < threshold
  hide_comp <- (d - n) > 0 & (d - n) < threshold
  out <- fmt_n(n, threshold)
  out[!hide_n & hide_comp] <- "suppressed"
  out
}

## Collapse rows whose count is a small cell into one "Other" row.
collapse_small <- function(dt, label_col, n_col = "N", other = "Other",
                           threshold = SMALL_CELL) {
  dt <- copy(as.data.table(dt))
  small <- dt[[n_col]] < threshold
  if (sum(small) < 2L) return(dt)
  keep <- dt[!small]
  oth  <- data.table(other, sum(dt[[n_col]][small]))
  setnames(oth, c(label_col, n_col))
  rbind(keep[, c(label_col, n_col), with = FALSE], oth)
}

## Source lines of named top-level objects in a file (with their leading
## "##" comment block) -- used to show the executed definitions on the site.
src_of <- function(names, file = "R/functions.R") {
  ex    <- parse(file, keep.source = TRUE)
  refs  <- attr(ex, "srcref")
  lines <- readLines(file)
  out   <- character()
  for (nm in names) {
    idx <- which(vapply(ex, function(e) is.call(e) &&
                          identical(e[[1]], as.name("<-")) &&
                          identical(e[[2]], as.name(nm)), logical(1)))
    stopifnot(length(idx) == 1L)
    s <- refs[[idx]][1]; e <- refs[[idx]][3]
    while (s > 1 && grepl("^##", lines[s - 1])) s <- s - 1
    out <- c(out, lines[s:e], "")
  }
  out
}
