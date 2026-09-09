source(if (file.exists("paths.R")) "paths.R" else "scripts/paths.R")
source("scripts/pipeline_lib.R")
suppressMessages({library(openxlsx); library(readxl)})

# the icd-9-cm codes the validation data has no icda-8 match for, checked
# against the 858 three character icda-8 rubrics in ICD_Codes_Labels.xlsx.
# six are marked wrong, where the icda-8 rubric is unambiguous and there is no
# competing candidate in the surrounding range. everything else is left as
# validation data correct. a plausible candidate is not a confirmed mapping and
# icd-8 inclusion notes are not in the repo, so the rest is not second guessed.

ORIG <- "data/original"
VAL  <- file.path(ORIG, "ICD_Codes_Files_and_Validation_Data/Validation_Data .xlsx")
LAB  <- file.path(ORIG, "ICD_Codes_Files_and_Validation_Data/ICD_Codes_Labels.xlsx")
OUT  <- out_path("52_codes_no_icda8.xlsx")

W <- "validation data wrong"; C <- "validation data correct"

# only the six where the icda-8 rubric is unambiguous are marked wrong. every
# other code is left as validation data correct, because a plausible candidate
# is not the same as a confirmed mapping and guessing would be worse than
# leaving the validation alone.
chk <- list(
  c("165", W, "163", "ICD-8 stops at 163 other and unspecified respiratory organs, there is no other home for ill-defined respiratory and intrathoracic sites"),
  c("175", W, "174", "ICD-8 did not split breast cancer by sex. 174 malignant neoplasm of breast is the only breast rubric and it covers both"),
  c("179", W, "182", "ICD-8 uterus is 180 cervix, 181 chorionepithelioma, 182 other. Part unspecified has only one home, 182"),
  c("555", W, "563", "563 chronic enteritis and ulcerative colitis is the only chronic inflammatory bowel rubric in ICD-8, regional enteritis sits there"),
  c("576", W, "576", "576 other diseases of gallbladder and biliary ducts is the same rubric under the same number"),
  c("745", W, "746", "746 congenital anomalies of heart is the only heart rubric. Note ICD-8 745 is ear face and neck, the numbers do not line up between the two classifications"),
  c("176", C, "", "no Kaposi sarcoma rubric in ICD-8, it entered ICD-9-CM later"),
  c("230", C, "", "ICD-8 230 to 239 is neoplasm of unspecified nature. There is no carcinoma in situ block in ICD-8"),
  c("231", C, "", "no carcinoma in situ block in ICD-8"),
  c("233", C, "", "no carcinoma in situ block in ICD-8"),
  c("234", C, "", "no carcinoma in situ block in ICD-8"),
  c("249", C, "", "ICD-8 has only 250 diabetes mellitus. Secondary diabetes is a later ICD-9-CM concept"),
  c("305", C, "", "ICD-8 has 304 drug dependence only. Nondependent abuse is a later concept"),
  c("315", C, "", "ICD-8 315 is unspecified mental retardation, a different concept"),
  c("327", C, "", "no sleep disorder rubric in ICD-8, organic sleep disorders entered ICD-9-CM in 2005"),
  c("338", C, "", "no pain NEC rubric in ICD-8, it entered ICD-9-CM in 2006"),
  c("495", C, "", "extrinsic allergic alveolitis has no ICD-8 rubric"),
  c("496", C, "", "chronic airway obstruction NEC is an ICD-9 concept")
)

excl <- read_excel(VAL, sheet = "Validation_ICD9_ICD8_Excld")
lab  <- read_excel(LAB, sheet = "CCS ICD-9-CM-3Level")
a8   <- read_excel(LAB, sheet = "ICDA-8-3Level")
a8lab <- setNames(as.character(a8$ICDA_8_LABEL), as.character(a8$ICDA_8))

m <- as.data.frame(do.call(rbind, chk), stringsAsFactors = FALSE)
names(m) <- c("code", "verdict", "cand", "note")

d <- data.frame(`ICD-9-CM` = as.character(excl$`ICD-9-CM`), check.names = FALSE,
                stringsAsFactors = FALSE)
i <- match(d$`ICD-9-CM`, as.character(lab$ICD_9_CM))
d$`ICD-9-CM label` <- as.character(lab$ICD_9_CM_LABEL)[i]
j <- match(d$`ICD-9-CM`, m$code)
d$`Validation correct?`   <- ifelse(is.na(j), C, m$verdict[j])
d$`Closest ICDA-8 code`   <- ifelse(is.na(j) | !nzchar(m$cand[j]), "none",
                                    paste(m$cand[j], a8lab[m$cand[j]]))
d$Notes <- ifelse(is.na(j), "", m$note[j])
d <- d[order(as.numeric(d$`ICD-9-CM`)), ]

## chapter the icd-9 code sits in, and the icda-8 codes it actually travels
## with in the co-occurrence data. she asked for both, the point being that the
## validation says no valid match but the data may still show a companion
CH9 <- c("1 infectious and parasitic", "2 neoplasms",
         "3 endocrine, nutritional, metabolic", "4 blood and blood-forming organs",
         "5 mental disorders", "6 nervous system and sense organs",
         "7 circulatory system", "8 respiratory system", "9 digestive system",
         "10 genitourinary system", "11 pregnancy and childbirth",
         "12 skin and subcutaneous tissue", "13 musculoskeletal and connective tissue",
         "14 congenital anomalies", "15 perinatal conditions",
         "16 symptoms, signs and ill-defined conditions", "17 injury and poisoning")

d$`ICD-9 chapter` <- CH9[vapply(d$`ICD-9-CM`, find_icd9cm_chapter, integer(1))]

co <- as.data.frame(load_cooccurrence_df(
  file.path(ORIG, "Co_occurrence/icd_8_9_co_occurrence_3d.xlsx")))
co$ICD_9_CM_Code <- as.character(co$ICD_9_CM_Code)
co$Co_Occurrence_Frequency <- as.numeric(co$Co_Occurrence_Frequency)

top3 <- function(k) {
  s <- co[co$ICD_9_CM_Code == k, ]
  if (!nrow(s)) return(c("none in the co-occurrence data", "", ""))
  s <- s[order(-s$Co_Occurrence_Frequency), ]
  out <- sprintf("%s %s (n=%d)", s$ICDA_8_Code, s$ICDA_8_LABEL, s$Co_Occurrence_Frequency)
  c(out, "", "", "")[1:3]
}
t3 <- t(vapply(d$`ICD-9-CM`, top3, character(3)))
d$`Top co-occurring ICDA-8`   <- t3[, 1]
d$`2nd co-occurring ICDA-8`   <- t3[, 2]
d$`3rd co-occurring ICDA-8`   <- t3[, 3]
d$`Co-occurrence rows` <- vapply(d$`ICD-9-CM`, function(k) sum(co$ICD_9_CM_Code == k), integer(1))

d <- d[, c("ICD-9-CM", "ICD-9-CM label", "ICD-9 chapter", "Validation correct?",
           "Closest ICDA-8 code", "Top co-occurring ICDA-8", "2nd co-occurring ICDA-8",
           "3rd co-occurring ICDA-8", "Co-occurrence rows", "Notes")]

wb <- createWorkbook()
addWorksheet(wb, "no ICDA-8 match")
writeData(wb, 1, d, headerStyle = createStyle(textDecoration = "bold", valign = "bottom"))
addStyle(wb, 1, createStyle(valign = "top", wrapText = TRUE),
         rows = 2:(nrow(d) + 1), cols = 1:ncol(d), gridExpand = TRUE, stack = TRUE)
setColWidths(wb, 1, cols = 1:10, widths = c(10, 40, 30, 22, 40, 46, 46, 46, 12, 60))
freezePane(wb, 1, firstActiveRow = 2)

addWorksheet(wb, "chapters")
ch <- as.data.frame(table(d$`ICD-9 chapter`), stringsAsFactors = FALSE)
names(ch) <- c("ICD-9 chapter", "codes with no ICDA-8 match")
ch <- ch[order(-ch$`codes with no ICDA-8 match`), ]
writeData(wb, 2, ch, headerStyle = createStyle(textDecoration = "bold"))
setColWidths(wb, 2, cols = 1:2, widths = c(44, 28))
freezePane(wb, 2, firstActiveRow = 2)

## none of the 52 appear in the icd-9 to icda-8 co-occurrence file at all, and
## that is the whole reason the pipeline stays silent on them. worth putting
## beside the icd-10 track, where all 9 excluded codes do have co-occurrence
## data and all 9 get mapped
c10 <- as.character(read_excel(file.path(ORIG, "Co_occurrence/icd_10_9_co_occurrence_3c.xlsx"),
                               sheet = 1)$ICD_9_CM_Code3)
e10 <- as.character(read_excel(VAL, sheet = "Validation_ICD9_ICD10_Excld")$`ICD-9-CM`)
m10 <- unique(as.character(read_excel(VAL, sheet = "Validation_ICD9_ICD10")$`ICD-9-CM`))
m8  <- unique(as.character(read_excel(VAL, sheet = "Validaion_ICD9_ICD8")$`ICD-9-CM`))
c8  <- unique(co$ICD_9_CM_Code)

cov <- data.frame(
  Track = c("ICD-9-CM to ICD-10-CA", "ICD-9-CM to ICD-10-CA",
            "ICD-9-CM to ICDA-8", "ICD-9-CM to ICDA-8"),
  Group = c("no valid match in the validation data", "has a valid match",
            "no valid match in the validation data", "has a valid match"),
  Codes = c(length(e10), length(m10), nrow(d), length(m8)),
  `With co-occurrence data` = c(sum(e10 %in% c10), sum(m10 %in% c10),
                                sum(d$`ICD-9-CM` %in% c8), sum(m8 %in% c8)),
  check.names = FALSE)
cov$Share <- sprintf("%.0f%%", 100 * cov$`With co-occurrence data` / cov$Codes)

addWorksheet(wb, "co-occurrence coverage")
writeData(wb, 3, cov, headerStyle = createStyle(textDecoration = "bold"))
setColWidths(wb, 3, cols = 1:5, widths = c(26, 40, 10, 24, 10))
freezePane(wb, 3, firstActiveRow = 2)

saveWorkbook(wb, OUT, overwrite = TRUE)
cat("wrote", OUT, "\n\n")
print(cov, row.names = FALSE)
cat("\n")
print(table(d$`Validation correct?`))
cat("\nby chapter:\n"); print(ch, row.names = FALSE)
cat("\ncodes with no co-occurrence data at all:",
    sum(d$`Co-occurrence rows` == 0), "\n")
