# ============================================================
# 脚本：NPS_R03_AnalyticDataset.R
# 功能：构建分析数据集 —— LOD/√2替代、变量标准化、多模块合并
# 项目：NHANES邻苯二甲酸酯与卒中关联研究 V5
# 创建日期：2026-05-08
#
# V5对齐V4策略：
#   - 8种代谢物（与V4一致：MBP, MBzP, MECPP, MEHHP, MEOHP, MCPP, MEP, MiBP）
#   - LOD替代：Below LOD → LOD/√2
#   - 排除：年龄<20、孕妇、卒中缺失、URXUCR无效、代谢物不完整、协变量不完整
#   - URXUCR：2015-2018从ALB_CR模块获取
#   - 高血压：BPQ020 | BPQ040A
#   - 糖尿病：Borderline → NA
#   - 饮酒：ALQ101+ALQ120Q（V4逻辑）
#   - 高脂血症：BPQ080（V4逻辑）
# ============================================================

library(dplyr)
library(tidyr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript 02_build_analytic_dataset.R <reproduction_workspace> <reproduction_outputs>")
}
workspace_root <- gsub("\\\\", "/", args[[1]])
output_root <- gsub("\\\\", "/", args[[2]])
raw_root <- file.path(workspace_root, "data_raw")
derived_root <- file.path(workspace_root, "data_derived")
dir.create(derived_root, recursive = TRUE, showWarnings = FALSE)
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

cat("══════════════════════════════════════════════\n")
cat(" NPS R03：构建分析数据集（对齐V4）\n")
cat("══════════════════════════════════════════════\n\n")

# ── 1. 元数据定义 ──────────────────────────────────────────

cycles <- data.frame(
  suffix  = c("C",    "D",    "E",    "F",    "G",    "H",    "I",    "J"),
  years   = c("2003_2004","2005_2006","2007_2008","2009_2010",
              "2011_2012","2013_2014","2015_2016","2017_2018"),
  pth_src = c("L24PH","PHTHTE","PHTHTE","PHTHTE",
              "PHTHTE","PHTHTE","PHTHTE","PHTHTE"),
  stringsAsFactors = FALSE
)

# V4的8种代谢物（含URXMC1=MCPP，不含MCNP/MCOP/MNBP）
v4_metab_vars <- c(
  "URXMBP"  = "MBP",
  "URXMZP"  = "MBzP",
  "URXECP"  = "MECPP",
  "URXMHH"  = "MEHHP",
  "URXMOH"  = "MEOHP",
  "URXMC1"  = "MCPP",
  "URXMEP"  = "MEP",
  "URXMIB"  = "MiBP"
)

# V5额外代谢物（MiNP检出率仅26%，移出主队列）
v5_extra_metab <- c(
  "URXCNP"  = "MCNP",
  "URXCOP"  = "MCOP"
  # URXMNP = MiNP（非MnBP），检出率26%，低于纳入门槛，暂不纳入
)

all_metab_vars <- c(v4_metab_vars, v5_extra_metab)

# LOD标记列
lod_flag_vars <- c(
  "URXMEP"  = "URDMEPLC",  "URXMBP"  = "URDMBPLC",
  "URXMZP"  = "URDMZPLC",  "URXMHP"  = "URDMHPLC",
  "URXMOH"  = "URDMOHLC",  "URXMHH"  = "URDMHHLC",
  "URXECP"  = "URDECPLC",  "URXMIB"  = "URDMIBLC",
  "URXCNP"  = "URDCNPLC",  "URXCOP"  = "URDCOPLC",
  "URXMC1"  = "URDMC1LC"
)

# 周期特定 LOD 值（ng/mL）
# 来源：NHANES Codebook L24PH_C 至 PHTHTE_J
# Values are fixed cycle-specific constants used by the registered reconstruction chain.
# URXCNP/URXCOP：2003-04 周期列缺失（P-04 已确认），设 NA 后 !is.na() 自动跳过
# URXMC1：2003-04 LOD=0.16（填充值≈0.113，与实测 0.1 吻合）；2013-18 LOD=0.40
lod_by_cycle <- list(
  "2003_2004" = c(URXMEP=0.264, URXMBP=0.40,  URXMZP=0.072, URXMHP=0.90,
                  URXMOH=0.45,  URXMHH=0.32,  URXECP=0.25,  URXMIB=0.26,
                  URXCNP=NA,    URXCOP=NA,    URXMC1=0.16),
  "2005_2006" = c(URXMEP=0.528, URXMBP=0.6,   URXMZP=0.216, URXMHP=1.2,
                  URXMOH=0.7,   URXMHH=0.7,   URXECP=0.6,   URXMIB=0.3,
                  URXCNP=0.6,   URXCOP=0.7,   URXMC1=0.20),
  "2007_2008" = c(URXMEP=0.462, URXMBP=0.6,   URXMZP=0.216, URXMHP=1.1,
                  URXMOH=0.6,   URXMHH=0.7,   URXECP=0.5,   URXMIB=0.3,
                  URXCNP=0.5,   URXCOP=0.7,   URXMC1=0.20),
  "2009_2010" = c(URXMEP=0.462, URXMBP=0.4,   URXMZP=0.216, URXMHP=0.5,
                  URXMOH=0.2,   URXMHH=0.2,   URXECP=0.2,   URXMIB=0.2,
                  URXCNP=0.2,   URXCOP=0.2,   URXMC1=0.20),
  "2011_2012" = c(URXMEP=0.6,   URXMBP=0.4,   URXMZP=0.3,   URXMHP=0.5,
                  URXMOH=0.2,   URXMHH=0.2,   URXECP=0.2,   URXMIB=0.2,
                  URXCNP=0.2,   URXCOP=0.2,   URXMC1=0.20),
  "2013_2014" = c(URXMEP=1.2,   URXMBP=0.4,   URXMZP=0.3,   URXMHP=0.8,
                  URXMOH=0.2,   URXMHH=0.4,   URXECP=0.4,   URXMIB=0.8,
                  URXCNP=0.2,   URXCOP=0.3,   URXMC1=0.40),
  "2015_2016" = c(URXMEP=1.2,   URXMBP=0.4,   URXMZP=0.3,   URXMHP=0.8,
                  URXMOH=0.2,   URXMHH=0.4,   URXECP=0.4,   URXMIB=0.8,
                  URXCNP=0.2,   URXCOP=0.3,   URXMC1=0.40),
  "2017_2018" = c(URXMEP=1.2,   URXMBP=0.4,   URXMZP=0.3,   URXMHP=0.8,
                  URXMOH=0.2,   URXMHH=0.4,   URXECP=0.4,   URXMIB=0.8,
                  URXCNP=0.2,   URXCOP=0.3,   URXMC1=0.40)
)

# ── NHANES变量转换 ─────────────────────────────────────────

nh_num <- function(x) {
  if (is.factor(x)) return(as.numeric(x))
  as.numeric(x)
}

# ── 2. 加载数据（精确复现V4 R02逻辑）────────────────────────

cat("── 加载数据（V4 R02风格）──\n\n")

load_cycle <- function(suffix, years, pth_src) {
  # DEMO（含RIDEXPRG）
  d <- readRDS(file.path(raw_root, years, paste0("DEMO_", suffix, ".rds")))
  keep_d <- intersect(c("SEQN","RIAGENDR","RIDAGEYR","RIDRETH1","DMDEDUC2",
                         "DMDMARTL","INDFMPIR","SDMVPSU","SDMVSTRA","WTMEC2YR",
                         "RIDEXPRG"), names(d))
  d <- d[, keep_d]
  d$cycle <- years

  # MCQ
  m <- readRDS(file.path(raw_root, years, paste0("MCQ_", suffix, ".rds")))
  m <- m[, intersect(c("SEQN","MCQ160F"), names(m))]

  # BPQ（含BPQ040A、BPQ080）
  bp <- readRDS(file.path(raw_root, years, paste0("BPQ_", suffix, ".rds")))
  bp <- bp[, intersect(c("SEQN","BPQ020","BPQ040A","BPQ080"), names(bp))]

  # DIQ
  di <- readRDS(file.path(raw_root, years, paste0("DIQ_", suffix, ".rds")))
  di <- di[, intersect(c("SEQN","DIQ010"), names(di))]

  # SMQ
  sm <- readRDS(file.path(raw_root, years, paste0("SMQ_", suffix, ".rds")))
  sm <- sm[, intersect(c("SEQN","SMQ020"), names(sm))]

  # ALQ（V4逻辑：ALQ111→ALQ101重命名）
  al <- readRDS(file.path(raw_root, years, paste0("ALQ_", suffix, ".rds")))
  alq_vars <- intersect(c("SEQN","ALQ101","ALQ111","ALQ120Q","ALQ130"), names(al))
  al <- al[, alq_vars, drop=FALSE]
  if ("ALQ111" %in% names(al) & !"ALQ101" %in% names(al)) {
    al$ALQ101 <- al$ALQ111
    cat("  ", suffix, ": ALQ111 → ALQ101 重命名\n")
  }

  # BMX
  bx <- readRDS(file.path(raw_root, years, paste0("BMX_", suffix, ".rds")))
  bx <- bx[, intersect(c("SEQN","BMXBMI"), names(bx))]

  # PHTHTE（代谢物）
  f <- readRDS(file.path(raw_root, years, paste0(pth_src, "_", suffix, ".rds")))
  if (is.null(f)) f <- data.frame(SEQN = integer(0))

  # Weight-Check-02 fix: Replace WTMEC2YR with phthalate subsample weight
  # 2003-2010, 2013-2018: WTSB2YR; 2011-2012: WTSA2YR
  ph_wt_var <- ifelse(years == "2011_2012", "WTSA2YR", "WTSB2YR")
  if (ph_wt_var %in% names(f)) {
    d$WTMEC2YR <- f[[ph_wt_var]][match(d$SEQN, f$SEQN)]
    cat("  ", suffix, ": WTMEC2YR replaced with", ph_wt_var, "from phthalate file\n")
  } else {
    cat("  [WARN] ", suffix, ":", ph_wt_var, "not found in phthalate file, keeping DEMO WTMEC2YR\n")
  }

  # URXUCR来源（V4逻辑：2015-2018从ALB_CR取）
  if (years %in% c("2015_2016","2017_2018")) {
    alb_name <- if (years == "2015_2016") "ALB_CR_I" else "ALB_CR_J"
    alb <- readRDS(file.path(raw_root, years, paste0(alb_name, ".rds")))
    ucr <- alb[, intersect(c("SEQN","URXUCR"), names(alb)), drop=FALSE]
    cat("  ", suffix, ": URXUCR from", alb_name, "\n")
  } else {
    ucr <- f[, intersect(c("SEQN","URXUCR"), names(f)), drop=FALSE]
  }

  # LOD替代（周期特定 LOD）
  cycle_lod <- lod_by_cycle[[years]]
  for (mv in names(all_metab_vars)) {
    if (!mv %in% names(f)) {
      f[[mv]] <- NA
      next
    }
    flag_col <- lod_flag_vars[mv]
    if (flag_col %in% names(f)) {
      flag_vals <- as.character(f[[flag_col]])
      below <- flag_vals == "Below lower detection limit"
      below[is.na(below)] <- FALSE
      lod_val <- cycle_lod[mv]
      if (sum(below) > 0 && !is.na(lod_val)) {
        f[[mv]][below] <- lod_val / sqrt(2)
      }
    }
  }

  # 选择代谢物列
  metab_keep <- c("SEQN", names(all_metab_vars))
  f <- f[, intersect(metab_keep, names(f)), drop=FALSE]
  for (mv in names(all_metab_vars)) {
    if (!mv %in% names(f)) f[[mv]] <- NA
  }

  # 合并
  df <- d
  df <- left_join(df, m,  by="SEQN")
  df <- left_join(df, bp, by="SEQN")
  df <- left_join(df, di, by="SEQN")
  df <- left_join(df, sm, by="SEQN")
  df <- left_join(df, al, by="SEQN")
  df <- left_join(df, bx, by="SEQN")
  df <- left_join(df, f,  by="SEQN")
  df <- left_join(df, ucr, by="SEQN")
  return(df)
}

raw_list <- lapply(1:nrow(cycles), function(i) {
  load_cycle(cycles$suffix[i], cycles$years[i], cycles$pth_src[i])
})
raw <- bind_rows(raw_list)
cat("\n原始合并完成：n =", nrow(raw), "\n\n")

# ── 3. 纳入排除（精确对齐V4 R03）──────────────────────────

cat("── 纳入排除 ──\n")

n0 <- nrow(raw)
n_initial <- n0

# Step 1: Age >= 20
df <- raw[!is.na(raw$RIDAGEYR) & raw$RIDAGEYR >= 20, ]
n_age_excluded <- n0 - nrow(df)
cat("1. Age >= 20：排除", n0 - nrow(df), "→", nrow(df), "\n")

# Step 2: 排除孕妇
n0 <- nrow(df)
pregnant <- as.character(df$RIDEXPRG) == "Yes, positive lab pregnancy test or self-reported pregnant at exam"
pregnant[is.na(pregnant)] <- FALSE
df <- df[!pregnant, ]
n_pregnancy_excluded <- n0 - nrow(df)
cat("2. 排除孕妇：排除", n0 - nrow(df), "→", nrow(df), "\n")

# Step 3: MCQ160F Yes/No
n0 <- nrow(df)
mcq_vals <- as.character(df$MCQ160F)
df <- df[mcq_vals %in% c("Yes","No"), ]
n_stroke_missing_excluded <- n0 - nrow(df)
cat("3. MCQ160F Yes/No：排除", n0 - nrow(df), "→", nrow(df), "\n")

# Step 4: URXUCR有效
# BEGIN S3 SNAPSHOT: isolated in-memory branch; never serialized.
s3_excluded <- df[is.na(df$URXUCR) | df$URXUCR <= 0,
                  c("SEQN", "BPQ020", "BPQ040A", "DIQ010", "BPQ080", "SMQ020"),
                  drop = FALSE]
# END S3 SNAPSHOT
n0 <- nrow(df)
df <- df[!is.na(df$URXUCR) & df$URXUCR > 0, ]
n_creatinine_unavailable_excluded <- n0 - nrow(df)
cat("4. URXUCR有效：排除", n0 - nrow(df), "→", nrow(df), "\n")

# Step 5: V4八代谢物ln有效（非NA、非Infinite、非负值）
n0 <- nrow(df)
for (mv in names(v4_metab_vars)) {
  clean <- df[[mv]]
  clean[clean <= 0] <- NA
  ln_name <- paste0(v4_metab_vars[mv], "_ln")
  df[[ln_name]] <- log(clean / df$URXUCR * 1000)
}
v4_ln_vars <- paste0(v4_metab_vars, "_ln")
mat <- as.matrix(df[, v4_ln_vars])
keep_metab <- rowSums(!is.na(mat) & !is.infinite(mat)) == length(v4_metab_vars)
df <- df[keep_metab, ]
n_incomplete_phthalates_excluded <- n0 - nrow(df)
cat("5. V4八代谢物ln有效：排除", n0 - nrow(df), "→", nrow(df), "\n")

# ── 4. 协变量编码（精确对齐V4 R03）──────────────────────────

cat("\n── 协变量编码 ──\n")

# 卒中
df$stroke <- ifelse(as.character(df$MCQ160F) == "Yes", 1L, 0L)

# 性别
df$sex <- factor(
  ifelse(as.character(df$RIAGENDR) == "Female", "Female",
  ifelse(as.character(df$RIAGENDR) == "Male", "Male", NA)),
  levels = c("Female", "Male"))

# 种族（V4编码）
df$race_eth <- as.character(df$RIDRETH1)
df$race_eth_f <- factor(
  ifelse(df$race_eth == "Mexican American", "Mexican American",
  ifelse(df$race_eth == "Other Hispanic", "Other Hispanic",
  ifelse(df$race_eth == "Non-Hispanic White", "Non-Hispanic White",
  ifelse(df$race_eth == "Non-Hispanic Black", "Non-Hispanic Black",
  ifelse(df$race_eth %in% c("Non-Hispanic Asian",
         "Other Race - Including Multi-Racial","Other/Multi-racial"),
         "Other/Multiracial", NA))))),
  levels = c("Non-Hispanic White","Mexican American","Other Hispanic",
             "Non-Hispanic Black","Other/Multiracial"))

# 教育（V4编码：college or above / HS / less than HS）
edu <- as.character(df$DMDEDUC2)
df$edu_f <- factor(
  ifelse(edu %in% c("Less than 9th grade","Less Than 9th Grade",
                     "9-11th grade (Includes 12th grade with no diploma)",
                     "9-11th Grade (Includes 12th grade with no diploma)"),
         "Less than high school",
  ifelse(edu %in% c("High school graduate/GED or equivalent",
                     "High School Grad/GED or Equivalent"),
         "High school graduate",
  ifelse(edu %in% c("Some college or AA degree","Some College or AA degree",
                     "College graduate or above","College Graduate or above"),
         "College or above", NA))),
  levels = c("College or above","High school graduate","Less than high school"))

# 婚姻（V4编码：3类）
mar <- as.character(df$DMDMARTL)
df$marital_f <- factor(
  ifelse(mar %in% c("Married","Living with partner"), "Married/partnered",
  ifelse(mar == "Never married", "Never married",
  ifelse(mar %in% c("Widowed","Divorced","Separated"),
         "Separated/divorced/widowed", NA))),
  levels = c("Married/partnered","Never married","Separated/divorced/widowed"))

# PIR
df$pir <- df$INDFMPIR

# BMI
df$bmi <- df$BMXBMI

# 吸烟（V4编码：二分类 Ever/Never）
df$smoking_f <- factor(
  ifelse(as.character(df$SMQ020) == "Yes", "Ever",
  ifelse(as.character(df$SMQ020) == "No", "Never", NA)),
  levels = c("Never","Ever"))

# 饮酒（V4编码：ALQ101+ALQ120Q复合逻辑）
alq101 <- as.character(df$ALQ101)
alq120q <- nh_num(df$ALQ120Q)
df$alcohol_f <- factor(
  ifelse(alq101 == "No", "No",
  ifelse(alq101 == "Yes" & !is.na(alq120q) & alq120q == 0, "No",
  ifelse(alq101 == "Yes", "Yes", NA))),
  levels = c("No","Yes"))

# 高血压（V4编码：BPQ020 | BPQ040A）
bpq020 <- as.character(df$BPQ020)
bpq040a <- as.character(df$BPQ040A)
df$hypertension_f <- factor(
  ifelse(bpq020 == "Yes", "Yes",
  ifelse(!is.na(bpq040a) & bpq040a == "Yes", "Yes",
  ifelse(bpq020 == "No", "No", NA))),
  levels = c("No","Yes"))

# 糖尿病（V4编码：Borderline → NA）
df$diabetes_f <- factor(
  ifelse(as.character(df$DIQ010) == "Yes", "Yes",
  ifelse(as.character(df$DIQ010) == "No", "No", NA)),
  levels = c("No","Yes"))

# 高脂血症（V4编码：BPQ080）
df$hyperlipidemia_f <- factor(
  ifelse(as.character(df$BPQ080) == "Yes", "Yes",
  ifelse(as.character(df$BPQ080) == "No", "No", NA)),
  levels = c("No","Yes"))

cat("编码完成。\n\n")

# ── 5. 完整病例筛选（对齐V4 R03 line 279-286）──────────────

cat("── 完整病例筛选 ──\n")

# V4的covariate_vars（R03 line 273-275）
v4_covariate_vars <- c("sex","RIDAGEYR","race_eth_f","edu_f","marital_f",
                        "pir","bmi","smoking_f","alcohol_f",
                        "hypertension_f","diabetes_f","hyperlipidemia_f")
v4_design_vars <- c("SDMVPSU","SDMVSTRA","WTMEC2YR")

# V4的条件：stroke非NA + 所有代谢物ln有效 + 所有协变量非NA + 权重>0
# stroke已保证非NA（Step3已排除）
# 代谢物ln已保证有效（Step5已排除）
# 检查协变量完整
covar_complete <- rowSums(!is.na(df[, v4_covariate_vars])) == length(v4_covariate_vars)
design_ok <- rowSums(!is.na(df[, v4_design_vars])) == length(v4_design_vars) & df$WTMEC2YR > 0

analytic_main <- df[covar_complete & design_ok, ]
n_covariates_weights_excluded <- nrow(df) - nrow(analytic_main)

# 生成V4风格的权重
n_cycles <- length(unique(analytic_main$cycle))
analytic_main$WTMEC_comb <- analytic_main$WTMEC2YR / n_cycles

cat("  主队列（对齐V4）：n =", nrow(analytic_main), "\n")
cat("  卒中：", sum(analytic_main$stroke == 1), "\n")
cat("  Female：", sum(analytic_main$sex == "Female"), "\n")
cat("  Male：", sum(analytic_main$sex == "Male"), "\n")
cat("  权重除数：", n_cycles, "周期\n\n")

# ── 6. 生成V5完整代谢物版本（含MCNP/MCOP）──

cat("── V5完整代谢物版本 ──\n")

# V5额外代谢物也做LOD替代 + 肌酐校正
for (mv in names(v5_extra_metab)) {
  if (mv %in% names(analytic_main)) {
    clean <- analytic_main[[mv]]
    clean[clean <= 0] <- NA
    ln_name <- paste0(v5_extra_metab[mv], "_ln")
    analytic_main[[ln_name]] <- log(clean / analytic_main$URXUCR * 1000)
  }
}

# 主分析集的10种代谢物ln变量（V4的8种 + 额外2种：MCNP/MCOP）
v5_metab_cols <- c(v4_ln_vars, paste0(v5_extra_metab, "_ln"))
for (mv in v5_metab_cols) {
  if (mv %in% names(analytic_main)) {
    n_na <- sum(is.na(analytic_main[[mv]]) | is.infinite(analytic_main[[mv]]))
    cat(sprintf("  %-10s: 缺失=%d / %d\n", mv, n_na, nrow(analytic_main)))
  }
}

# ── 7. 保存 ────────────────────────────────────────────────

cat("\n── 保存数据集 ──\n")

eligible_cycles <- c(
  "2005_2006", "2007_2008", "2009_2010", "2011_2012",
  "2013_2014", "2015_2016", "2017_2018"
)
mixture_vars <- c(
  "MBP_ln", "MBzP_ln", "MECPP_ln", "MEHHP_ln", "MEOHP_ln",
  "MCPP_ln", "MEP_ln", "MiBP_ln", "MCNP_ln", "MCOP_ln"
)
mixture_covars <- c(
  "RIDAGEYR", "race_eth_f", "edu_f", "marital_f", "pir", "bmi",
  "smoking_f", "alcohol_f", "hypertension_f", "diabetes_f", "hyperlipidemia_f"
)
eligible <- analytic_main[
  analytic_main$cycle %in% eligible_cycles,
  c("SEQN", "cycle", mixture_vars, "stroke", "sex", mixture_covars),
  drop = FALSE
]
mixture_complete <- complete.cases(
  eligible[, c(mixture_vars, "stroke", "sex", mixture_covars), drop = FALSE]
)
mixture_sample <- eligible[mixture_complete, , drop = FALSE]

anchors <- data.frame(
  metric = c(
    "initial_pooled_sample", "age_exclusion", "pregnancy_exclusion",
    "missing_stroke_exclusion", "creatinine_data_unavailable_exclusion",
    "incomplete_phthalates_exclusion", "missing_covariates_or_nonpositive_weights_exclusion",
    "final_analytic_n", "stroke_events", "female_n", "male_n",
    "common_cycle_n", "common_cycle_stroke_events",
    "mixture_female_n", "mixture_female_events", "mixture_male_n", "mixture_male_events"
  ),
  actual = c(
    n_initial, n_age_excluded, n_pregnancy_excluded,
    n_stroke_missing_excluded, n_creatinine_unavailable_excluded,
    n_incomplete_phthalates_excluded, n_covariates_weights_excluded,
    nrow(analytic_main), sum(analytic_main$stroke == 1),
    sum(analytic_main$sex == "Female"), sum(analytic_main$sex == "Male"),
    nrow(eligible), sum(eligible$stroke == 1),
    sum(mixture_sample$sex == "Female"),
    sum(mixture_sample$stroke[mixture_sample$sex == "Female"] == 1),
    sum(mixture_sample$sex == "Male"),
    sum(mixture_sample$stroke[mixture_sample$sex == "Male"] == 1)
  ),
  expected = c(
    80312, 35522, 941, 68, 23403, 7049, 4200,
    9129, 375, 4611, 4518, 8311, 339, 4179, 175, 4132, 164
  ),
  stringsAsFactors = FALSE
)
anchors$pass <- anchors$actual == anchors$expected
write.csv(
  anchors, file.path(output_root, "phase_a_anchor_results.csv"),
  row.names = FALSE, na = "", fileEncoding = "UTF-8"
)

participant_flow <- data.frame(
  step = c(
    "initial pooled sample", "after age exclusion", "after pregnancy exclusion",
    "after stroke-status exclusion", "after creatinine-data-unavailable exclusion",
    "after incomplete-phthalates exclusion", "final after covariates/weights exclusion"
  ),
  retained_n = c(
    n_initial,
    n_initial - n_age_excluded,
    n_initial - n_age_excluded - n_pregnancy_excluded,
    n_initial - n_age_excluded - n_pregnancy_excluded - n_stroke_missing_excluded,
    n_initial - n_age_excluded - n_pregnancy_excluded - n_stroke_missing_excluded - n_creatinine_unavailable_excluded,
    nrow(df),
    nrow(analytic_main)
  ),
  stringsAsFactors = FALSE
)
write.csv(
  participant_flow, file.path(output_root, "phase_a_participant_flow.csv"),
  row.names = FALSE, na = "", fileEncoding = "UTF-8"
)

if (!all(anchors$pass)) {
  print(anchors[!anchors$pass, , drop = FALSE])
  stop("PHASE_A_ANCHOR_MISMATCH")
}

saveRDS(analytic_main, file.path(derived_root, "NPS_DATA_MAIN_analytic_V5.rds"))
cat("已保存：workspace/data_derived/NPS_DATA_MAIN_analytic_V5.rds\n")

df_female <- analytic_main[analytic_main$sex == "Female", ]
df_male   <- analytic_main[analytic_main$sex == "Male", ]
saveRDS(df_female, file.path(derived_root, "NPS_DATA_FEMALE_analytic_V5.rds"))
saveRDS(df_male,   file.path(derived_root, "NPS_DATA_MALE_analytic_V5.rds"))
cat("已保存：FEMALE (n=", nrow(df_female), "), MALE (n=", nrow(df_male), ")\n")

# ── 8. 验证报告 ────────────────────────────────────────────

cat("\n══════════════════════════════════════════════\n")
cat(" 验证报告（V4对齐）\n")
cat("══════════════════════════════════════════════\n\n")

cat(sprintf("  %-20s  %8s  %8s\n", "指标", "V5实际", "V4目标"))
cat(paste(rep("-", 50), collapse=""), "\n")
cat(sprintf("  %-20s  %8d  %8s\n", "总样本量", nrow(analytic_main), "9,129"))
cat(sprintf("  %-20s  %8d  %8s\n", "卒中", sum(analytic_main$stroke==1), "375"))
cat(sprintf("  %-20s  %8d  %8s\n", "Female", sum(analytic_main$sex=="Female"), "4,611"))
cat(sprintf("  %-20s  %8d  %8s\n", "Male", sum(analytic_main$sex=="Male"), "4,518"))
cat(sprintf("  %-20s  %8d  %8s\n", "差异", nrow(analytic_main)-9129, ""))
cat(sprintf("  %-20s  %8s  %8s\n", "URXSG in data",
            "URXSG" %in% names(analytic_main), "---"))

cat("\n── 周期分布 ──\n")
print(table(analytic_main$cycle))

cat("\n── 排除流程汇总 ──\n")
cat(sprintf("  原始总人数：%d\n", nrow(raw)))
cat(sprintf("  Age < 20：%d\n", nrow(raw[is.na(raw$RIDAGEYR) | raw$RIDAGEYR < 20,])))
cat(sprintf("  孕妇：%d\n", n_pregnancy_excluded))
cat(sprintf("  MCQ160F缺失/无效：%d\n", n_stroke_missing_excluded))
cat(sprintf("  URXUCR无效：%d\n", n_creatinine_unavailable_excluded))
cat(sprintf("  V4代谢物不完整：%d\n", n_incomplete_phthalates_excluded))
cat(sprintf("  协变量不完整/权重无效：%d\n", n_covariates_weights_excluded))
cat(sprintf("  最终：%d\n", nrow(analytic_main)))

cat("\n✓ R03 分析数据集构建完成（对齐V4）。\n")

# BEGIN S3 SUMMARY: four corrected variables only; aggregate output only.
stopifnot(nrow(analytic_main) == 9129L, nrow(s3_excluded) == 23403L,
          !anyDuplicated(analytic_main$SEQN), !anyDuplicated(s3_excluded$SEQN),
          length(intersect(analytic_main$SEQN, s3_excluded$SEQN)) == 0L)
s3_bpq020 <- as.character(s3_excluded$BPQ020)
s3_bpq040a <- as.character(s3_excluded$BPQ040A)
# Preserve the existing nested ifelse, including its NA propagation.
s3_ex <- list(
  hypertension = ifelse(s3_bpq020 == "Yes", 1L,
    ifelse(!is.na(s3_bpq040a) & s3_bpq040a == "Yes", 1L,
    ifelse(s3_bpq020 == "No", 0L, NA_integer_))),
  diabetes = ifelse(as.character(s3_excluded$DIQ010) == "Yes", 1L,
    ifelse(as.character(s3_excluded$DIQ010) == "No", 0L, NA_integer_)),
  hyperlipidemia = ifelse(as.character(s3_excluded$BPQ080) == "Yes", 1L,
    ifelse(as.character(s3_excluded$BPQ080) == "No", 0L, NA_integer_)),
  smoking = ifelse(as.character(s3_excluded$SMQ020) == "Yes", 1L,
    ifelse(as.character(s3_excluded$SMQ020) == "No", 0L, NA_integer_))
)
s3_in <- list(
  hypertension = as.integer(analytic_main$hypertension_f == "Yes"),
  diabetes = as.integer(analytic_main$diabetes_f == "Yes"),
  hyperlipidemia = as.integer(analytic_main$hyperlipidemia_f == "Yes"),
  smoking = as.integer(analytic_main$smoking_f == "Ever")
)
s3_rows <- lapply(names(s3_in), function(v) {
  p1 <- mean(s3_in[[v]], na.rm = TRUE)
  p2 <- mean(s3_ex[[v]], na.rm = TRUE)
  smd <- abs(p1 - p2) / sqrt((p1 * (1 - p1) + p2 * (1 - p2)) / 2)
  do.call(rbind, lapply(c("analytic", "excluded_creatinine_unavailable"), function(g) {
    x <- if (g == "analytic") s3_in[[v]] else s3_ex[[v]]
    data.frame(variable = v, group = g, group_n = length(x),
      positive_n = sum(x == 1L, na.rm = TRUE), nonmissing_n = sum(!is.na(x)),
      missing_n = sum(is.na(x)), prevalence = mean(x, na.rm = TRUE), SMD = smd,
      stringsAsFactors = FALSE)
  }))
})
s3_summary <- do.call(rbind, s3_rows)
stopifnot(nrow(s3_summary) == 8L, all(is.finite(s3_summary$SMD)))
write.csv(s3_summary, file.path(output_root, "s3_four_variable_summary.csv"),
          row.names = FALSE, na = "", fileEncoding = "UTF-8")
cat("S3 aggregate summary written; no participant-level S3 output.\n")
# END S3 SUMMARY
