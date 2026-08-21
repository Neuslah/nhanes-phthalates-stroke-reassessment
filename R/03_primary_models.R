# ============================================================
# 脚本：NPS_R20_MainRegression.R
# 功能：主分析——加权Logistic回归
#   - Model 1：年龄 + 种族（sex-stratified）
#   - Model 2：+ 教育、婚姻、PIR
#   - Model 3：+ BMI、吸烟、饮酒、高血压、糖尿病、高脂血症
#   - 10种代谢物（V4的8种 + MCNP/MCOP），逐个单独入模
#   - 按性别分层（Female / Male）
#   - 代谢物直接使用 _ln 变量（已在R03中完成 log(metab/UCR*1000)）
# 输出：NPS_OUT_TBL_MainRegression_Female_V5.csv
#        NPS_OUT_TBL_MainRegression_Male_V5.csv
# ============================================================

library(dplyr)
library(survey)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript 03_primary_models.R <analytic_rds> <output_directory>")
}
analytic_rds <- gsub("\\\\", "/", args[[1]])
output_dir <- gsub("\\\\", "/", args[[2]])
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("══════════════════════════════════════════════\n")
cat(" NPS R20：主分析回归\n")
cat("══════════════════════════════════════════════\n\n")

# ── 1. 读入分析集 ──────────────────────────────────────────

analytic_main <- readRDS(analytic_rds)
cat("读入分析集：n =", nrow(analytic_main), "\n")

# V4的8种代谢物（ln转换后）
v4_metabs <- c("MBP_ln", "MBzP_ln", "MECPP_ln", "MEHHP_ln",
               "MEOHP_ln", "MCPP_ln", "MEP_ln",  "MiBP_ln")

# V5额外代谢物
v5_extra  <- c("MCNP_ln", "MCOP_ln")

# 全部10种
all_metabs <- c(v4_metabs, v5_extra)

# 确认全部10种变量存在；缺失时 fail-fast，避免静默改变分析范围。
missing_metabs <- setdiff(all_metabs, names(analytic_main))
if (length(missing_metabs)) {
  stop("Missing required metabolites: ", paste(missing_metabs, collapse = ", "))
}
avail_metabs <- all_metabs
cat("可用代谢物：", length(avail_metabs), "种\n\n")

# ── 2. 构建调查设计对象 ────────────────────────────────────

svy_all <- svydesign(
  ids     = ~SDMVPSU,
  strata  = ~SDMVSTRA,
  weights = ~WTMEC_comb,
  data    = analytic_main,
  nest    = TRUE
)
svy_female <- subset(svy_all, sex == "Female")
svy_male   <- subset(svy_all, sex == "Male")

cat("全样本 n =", nrow(analytic_main), "\n")
cat("女性   n =", sum(analytic_main$sex == "Female"), "\n")
cat("男性   n =", sum(analytic_main$sex == "Male"), "\n\n")

# ── 3. 模型公式（精确对齐V4 R20）───────────────────────────

# Model 1: 年龄 + 种族
covars_m1_strat <- "RIDAGEYR + race_eth_f"

# Model 2: + 教育、婚姻、PIR
covars_m2_strat <- paste(covars_m1_strat,
                           "+ edu_f + marital_f + pir")

# Model 3: + BMI、吸烟、饮酒、高血压、糖尿病、高脂血症
covars_m3_strat <- paste(covars_m2_strat,
                           "+ bmi + smoking_f + alcohol_f",
                           "+ hypertension_f + diabetes_f + hyperlipidemia_f")

# ── 4. 回归函数 ────────────────────────────────────────────

run_svyglm <- function(design, outcome, exposure, covars) {
  formula_str <- paste(outcome, "~", exposure, "+", covars)
  f <- as.formula(formula_str)
  tryCatch({
    fit <- svyglm(f, design = design,
                  family = quasibinomial(link = "logit"))
    coef_val <- coef(fit)[exposure]
    se_val   <- sqrt(vcov(fit)[exposure, exposure])
    OR    <- exp(coef_val)
    OR_lo <- exp(coef_val - 1.96 * se_val)
    OR_hi <- exp(coef_val + 1.96 * se_val)
    p_val <- 2 * pnorm(-abs(coef_val / se_val))
    list(OR=OR, OR_lo=OR_lo, OR_hi=OR_hi, p=p_val, converged=TRUE)
  }, error = function(e) {
    cat("  [ERROR]", formula_str, ":", e$message, "\n")
    list(OR=NA, OR_lo=NA, OR_hi=NA, p=NA, converged=FALSE)
  })
}

fmt_or <- function(OR, lo, hi) {
  sprintf("%.3f (%.3f-%.3f)", OR, lo, hi)
}

# ── 5. 主循环 ──────────────────────────────────────────────

cat("── 开始回归分析 ──\n\n")

results <- list()

for (metab in avail_metabs) {
  cat("  代谢物：", metab, "\n")

  # 确定协变量（V4的8种用全集，V5额外的如果有NA则用complete case）
  n_na <- sum(is.na(analytic_main[[metab]]))

  for (model_id in 1:3) {

    # 选择协变量
    cv_strat <- switch(model_id,
                       covars_m1_strat, covars_m2_strat, covars_m3_strat)

    # 女性
    res_f   <- run_svyglm(svy_female, "stroke", metab, cv_strat)
    # 男性
    res_m   <- run_svyglm(svy_male,   "stroke", metab, cv_strat)

    results[[length(results) + 1]] <- data.frame(
      Metabolite  = metab,
      Model       = paste0("Model", model_id),
      OR_Female   = res_f$OR,
      CI_Female   = fmt_or(res_f$OR, res_f$OR_lo, res_f$OR_hi),
      P_Female    = round(res_f$p, 4),
      OR_Male     = res_m$OR,
      CI_Male     = fmt_or(res_m$OR, res_m$OR_lo, res_m$OR_hi),
      P_Male      = round(res_m$p, 4),
      stringsAsFactors = FALSE
    )

    cat(sprintf("    Model%d | F: %s | M: %s\n",
                model_id,
                fmt_or(res_f$OR,   res_f$OR_lo,   res_f$OR_hi),
                fmt_or(res_m$OR,   res_m$OR_lo,   res_m$OR_hi)))
  }
  cat("\n")
}

results_df <- do.call(rbind, results)
rownames(results_df) <- NULL

# ── 6. 分性别输出 ──────────────────────────────────────────

cat("── 保存结果 ──\n")

# Female 表
female_out <- results_df[, c("Metabolite", "Model",
                              "OR_Female", "CI_Female", "P_Female")]
names(female_out) <- c("Metabolite", "Model", "OR", "95% CI", "P value")
write.csv(female_out, file.path(output_dir, "NPS_OUT_TBL_MainRegression_Female_V5.csv"),
          row.names=FALSE)
cat("已保存：output_tables/NPS_OUT_TBL_MainRegression_Female_V5.csv\n")

# Male 表
male_out <- results_df[, c("Metabolite", "Model",
                            "OR_Male", "CI_Male", "P_Male")]
names(male_out) <- c("Metabolite", "Model", "OR", "95% CI", "P value")
write.csv(male_out, file.path(output_dir, "NPS_OUT_TBL_MainRegression_Male_V5.csv"),
          row.names=FALSE)
cat("已保存：output_tables/NPS_OUT_TBL_MainRegression_Male_V5.csv\n")

# ── 7. 核查报告 ─────────────────────────────────────────────

cat("\n══════════════════════════════════════════════\n")
cat(" 核查报告\n")
cat("══════════════════════════════════════════════\n\n")

# V4关键指标
mbzp_f_m3 <- results_df[results_df$Metabolite == "MBzP_ln" &
                          results_df$Model == "Model3", ]
cat("── V4关键指标对照 ──\n")
cat(sprintf("  Female MBzP Model3 OR = %.3f  95%%CI: %s  P = %.4f\n",
            mbzp_f_m3$OR_Female, mbzp_f_m3$CI_Female, mbzp_f_m3$P_Female))
cat("  V4目标: OR = 1.135, 95%CI: 0.968\u20131.332, P = 0.118\n\n")

# Model 3 全结果（Female）
cat("── Model 3 全结果（Female）──\n")
m3_f <- results_df[results_df$Model == "Model3",
                    c("Metabolite", "OR_Female", "CI_Female", "P_Female")]
print(m3_f, row.names=FALSE)

# Model 3 全结果（Male）
cat("\n── Model 3 全结果（Male）──\n")
m3_m <- results_df[results_df$Model == "Model3",
                    c("Metabolite", "OR_Male", "CI_Male", "P_Male")]
print(m3_m, row.names=FALSE)

# MBzP 三模型趋势（Female）
cat("\n── MBzP 三模型趋势（Female）──\n")
mbzp_trend <- results_df[results_df$Metabolite == "MBzP_ln",
                          c("Model", "OR_Female", "CI_Female", "P_Female")]
print(mbzp_trend, row.names=FALSE)

cat("\n✓ R20 主分析回归完成。\n")
