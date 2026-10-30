#!/usr/bin/env Rscript
# plexz 多轮模拟测试
# 测试覆盖: 基本拟合、检测上下限、多因子、边界情况、噪声鲁棒性、4PL/5PL

script_dir <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) getwd())
if (identical(Sys.getenv("R_PACKAGE_NAME"), "plexz")) {
  library(plexz)
} else {
  pkg_dir <- file.path(dirname(script_dir), "R")
  if (!dir.exists(pkg_dir)) pkg_dir <- file.path(getwd(), "plexz", "R")
  if (!dir.exists(pkg_dir)) pkg_dir <- file.path(dirname(getwd()), "plexz", "R")
  for (f in sort(list.files(pkg_dir, full.names = TRUE))) source(f)
}

pass_count <- 0
fail_count <- 0

assert <- function(cond, msg) {
  if (isTRUE(cond)) {
    pass_count <<- pass_count + 1
    cat("  [PASS]", msg, "\n")
  } else {
    fail_count <<- fail_count + 1
    cat("  [FAIL]", msg, "\n")
  }
}

cat("========================================\n")
cat(" msdconvert 多轮模拟测试\n")
cat("========================================\n\n")

# ─────────────────────────────────────────────
# 第 1 轮: 基本线性拟合 + 检测限
# ─────────────────────────────────────────────
cat("--- 第 1 轮: 基本线性拟合 + 检测限 ---\n")
set.seed(100)
sig1 <- seq(100, 10000, length.out = 8)
conc1 <- seq(1, 100, length.out = 8) + rnorm(8, 0, 0.5)
unk1 <- c(50, 3000, 5000, 12000)

res1 <- msd_to_concentration(sig1, conc1, unk1, method = "linear")
assert(length(res1$concentrations) == 4, "线性拟合返回 4 个浓度值")
assert(all(!is.na(res1$concentrations)), "所有浓度值非 NA")
assert(res1$status[1] == "below_LLOQ", "信号 50 低于 LLOQ")
assert(res1$status[4] == "above_ULOQ", "信号 12000 高于 ULOQ")
assert(res1$status[2] == "in_range", "信号 3000 在范围内")
assert(abs(res1$LLOQ - min(conc1)) < 1e-10, "LLOQ 等于最低标准品浓度")
assert(abs(res1$ULOQ - max(conc1)) < 1e-10, "ULOQ 等于最高标准品浓度")
cat("\n")

# ─────────────────────────────────────────────
# 第 2 轮: 对数拟合 + 正值校验
# ─────────────────────────────────────────────
cat("--- 第 2 轮: 对数拟合 ---\n")
set.seed(200)
sig2 <- exp(seq(3, 8, length.out = 7))
conc2 <- seq(0.5, 50, length.out = 7) + rnorm(7, 0, 0.3)
unk2 <- c(exp(5), exp(6.5))

res2 <- msd_to_concentration(sig2, conc2, unk2, method = "log")
assert(length(res2$concentrations) == 2, "对数拟合返回 2 个浓度值")
assert(all(res2$concentrations > 0), "对数拟合浓度为正")
assert(res2$status[1] == "in_range", "exp(5) 在范围内")

err2 <- tryCatch(msd_to_concentration(c(-1, 2, 3), c(1, 2, 3), c(1), method = "log"),
                 error = function(e) e$message)
assert(grepl("signal values > 0", err2), "负信号值触发错误")
cat("\n")

# ─────────────────────────────────────────────
# 第 3 轮: 幂函数拟合
# ─────────────────────────────────────────────
cat("--- 第 3 轮: 幂函数拟合 ---\n")
set.seed(300)
sig3 <- seq(10, 1000, length.out = 6)
conc3 <- 2 * sig3^0.8 + rnorm(6, 0, 1)
unk3 <- c(100, 500)

res3 <- msd_to_concentration(sig3, conc3, unk3, method = "power")
assert(length(res3$concentrations) == 2, "幂函数拟合返回 2 个浓度值")
assert(all(res3$concentrations > 0), "幂函数拟合浓度为正")
cat("\n")

# ─────────────────────────────────────────────
# 第 4 轮: 4PL 拟合 (S 型曲线)
# ─────────────────────────────────────────────
cat("--- 第 4 轮: 4PL 拟合 ---\n")
set.seed(400)
true_a <- 0.5; true_b <- 1.2; true_c <- 5000; true_d <- 200
sig4 <- seq(100, 50000, length.out = 12)
conc4 <- true_d + (true_a - true_d) / (1 + (sig4 / true_c)^true_b) + rnorm(12, 0, 0.5)
conc4 <- pmax(conc4, 0.01)
unk4 <- c(50, 5000, 30000)

res4 <- msd_to_concentration(sig4, conc4, unk4, method = "4pl")
assert(length(res4$concentrations) == 3, "4PL 返回 3 个浓度值")
assert(all(is.finite(res4$concentrations)), "4PL 浓度值有限")
assert(res4$status[1] == "below_LLOQ", "4PL: 信号 50 低于 LLOQ (min signal=100)")
assert(res4$status[2] == "in_range", "4PL: 信号 5000 在范围内")
cat("\n")

# ─────────────────────────────────────────────
# 第 5 轮: 5PL 拟合
# ─────────────────────────────────────────────
cat("--- 第 5 轮: 5PL 拟合 ---\n")
set.seed(500)
sig5 <- seq(100, 80000, length.out = 10)
conc5 <- 300 + (0.2 - 300) / (1 + (sig5 / 10000)^1.5)^0.8 + rnorm(10, 0, 0.3)
conc5 <- pmax(conc5, 0.01)
unk5 <- c(1000, 20000)

res5 <- msd_to_concentration(sig5, conc5, unk5, method = "5pl")
assert(length(res5$concentrations) == 2, "5PL 返回 2 个浓度值")
assert(all(is.finite(res5$concentrations)), "5PL 浓度值有限")
cat("\n")

# ─────────────────────────────────────────────
# 第 6 轮: 多因子批量处理
# ─────────────────────────────────────────────
cat("--- 第 6 轮: 多因子批量处理 ---\n")
set.seed(600)
il6_signal <- pmax(seq(100, 50000, length.out = 8) + rnorm(8, 0, 200), 50)
tnf_signal <- pmax(seq(200, 80000, length.out = 8) + rnorm(8, 0, 300), 100)
ifn_signal <- pmax(seq(50, 30000, length.out = 8) + rnorm(8, 0, 100), 20)

calib6 <- data.frame(
  factor = rep(c("IL-6", "TNF-a", "IFN-g"), each = 8),
  signal = c(il6_signal, tnf_signal, ifn_signal),
  concentration = c(
    c(0.5, 1, 2, 5, 10, 25, 50, 100),
    c(1, 2, 5, 10, 20, 50, 100, 200),
    c(0.2, 0.5, 1, 2, 5, 10, 20, 50)
  )
)
samples6 <- data.frame(
  factor = c("IL-6", "IL-6", "IL-6", "TNF-a", "TNF-a", "IFN-g", "IFN-g"),
  signal = c(30, 12000, 60000, 500, 50000, 30, 15000)
)

result6 <- msd_multi_factor(calib6, samples6,
  method = list("IL-6" = "log", "TNF-a" = "linear", "IFN-g" = "quadratic"))

assert(nrow(result6) == 7, "多因子返回 7 行")
assert(all(c("factor", "signal", "concentration", "status", "LLOQ", "ULOQ", "method") %in% names(result6)),
       "输出包含所有必要列")
assert(result6$status[1] == "below_LLOQ", "IL-6 信号 30 低于 LLOQ")
assert(result6$status[3] == "above_ULOQ", "IL-6 信号 60000 高于 ULOQ")
assert(all(abs(result6$LLOQ[1:3] - 0.5) < 1e-10), "IL-6 LLOQ = 0.5")
assert(all(abs(result6$ULOQ[1:3] - 100) < 1e-10), "IL-6 ULOQ = 100")
assert(all(result6$method[1:3] == "log"), "IL-6 使用 log 方法")

# 检测限汇总
lim_summary <- summarize_limits(result6)
assert(nrow(lim_summary) == 3, "汇总返回 3 个因子")
assert("LLOQ" %in% names(lim_summary), "汇总包含 LLOQ 列")
assert("ULOQ" %in% names(lim_summary), "汇总包含 ULOQ 列")
cat("  检测限汇总:\n")
print(lim_summary)
cat("\n")

# ─────────────────────────────────────────────
# 第 7 轮: 边界情况 - 不足 3 个标准品
# ─────────────────────────────────────────────
cat("--- 第 7 轮: 边界情况 ---\n")
calib7 <- data.frame(
  factor = c("X", "X"),
  signal = c(100, 200),
  concentration = c(1, 2)
)
samples7 <- data.frame(factor = "X", signal = 150)

res7 <- tryCatch(
  msd_to_concentration(c(100, 200), c(1, 2), c(150), method = "linear"),
  error = function(e) e$message
)
assert(grepl("At least 3", res7), "少于 3 个标准品触发错误")

res7b <- suppressWarnings(msd_multi_factor(calib7, samples7, method = "linear"))
assert(res7b$status == "error", "多因子中不足 3 点标记为 error")
cat("\n")

# ─────────────────────────────────────────────
# 第 8 轮: 高噪声鲁棒性测试
# ─────────────────────────────────────────────
cat("--- 第 8 轮: 高噪声鲁棒性 ---\n")
set.seed(800)
sig8 <- seq(100, 10000, length.out = 10)
conc8 <- seq(1, 100, length.out = 10) + rnorm(10, 0, 3)
unk8 <- c(2000, 5000, 8000)

res8 <- msd_to_concentration(sig8, conc8, unk8, method = "linear")
assert(all(is.finite(res8$concentrations)), "高噪声下线性拟合仍产出有限值")

cmp8 <- compare_methods(sig8, conc8, methods = c("linear", "quadratic", "cubic"))
assert(nrow(cmp8) == 3, "方法对比返回 3 行")
assert(all(!is.na(cmp8$R2)), "高噪声下所有方法均有 R2")
cat("\n")

# ─────────────────────────────────────────────
# 第 9 轮: 4PL 多因子 + 检测限展示
# ─────────────────────────────────────────────
cat("--- 第 9 轮: 4PL 多因子 ---\n")
set.seed(900)
n_pts <- 12
vegf_sig <- seq(100, 100000, length.out = n_pts)
vegf_conc <- 100 + (0.1 - 100) / (1 + (vegf_sig / 10000)^1.5) + rnorm(n_pts, 0, 0.5)
vegf_conc <- pmax(vegf_conc, 0.01)

mcp_sig <- seq(50, 60000, length.out = n_pts)
mcp_conc <- 50 + (0.05 - 50) / (1 + (mcp_sig / 8000)^1.2) + rnorm(n_pts, 0, 0.3)
mcp_conc <- pmax(mcp_conc, 0.01)

calib9 <- data.frame(
  factor = c(rep("VEGF", n_pts), rep("MCP-1", n_pts)),
  signal = c(vegf_sig, mcp_sig),
  concentration = c(vegf_conc, mcp_conc)
)
samples9 <- data.frame(
  factor = c("VEGF", "VEGF", "MCP-1", "MCP-1"),
  signal = c(50, 20000, 100, 40000)
)

result9 <- msd_multi_factor(calib9, samples9,
  method = list("VEGF" = "4pl", "MCP-1" = "4pl"))

assert(nrow(result9) == 4, "4PL 多因子返回 4 行")
assert(all(result9$status != "error"), "4PL 多因子无 error 状态")
assert(result9$status[1] == "below_LLOQ", "VEGF 信号 50 低于 LLOQ")
assert(all(result9$method == "4pl"), "方法列正确记录 4pl")

lim9 <- summarize_limits(result9)
assert(nrow(lim9) == 2, "4PL 汇总返回 2 个因子")
cat("  4PL 检测限汇总:\n")
print(lim9)
cat("\n")

# ─────────────────────────────────────────────
# 第 10 轮: 方法对比 + 自动选择最优
# ─────────────────────────────────────────────
cat("--- 第 10 轮: 方法对比与自动选择 ---\n")
set.seed(1000)
sig10 <- seq(50, 50000, length.out = 12)
conc10 <- 200 + (0.5 - 200) / (1 + (sig10 / 8000)^1.3) + rnorm(12, 0, 0.5)
conc10 <- pmax(conc10, 0.01)

cmp10 <- compare_methods(sig10, conc10,
  methods = c("linear", "quadratic", "log", "power", "4pl"))
assert(nrow(cmp10) == 5, "对比 5 种方法")
best <- cmp10$method[1]
assert(!is.na(cmp10$R2[1]), "最优方法有 R2 值")
cat("  最优方法:", best, " R2 =", round(cmp10$R2[1], 4), "\n")

res10 <- msd_to_concentration(sig10, conc10, c(2000, 10000), method = best)
assert(all(is.finite(res10$concentrations)), "最优方法预测正常")
cat("\n")

# ─────────────────────────────────────────────
# 第 11 轮: 重复性测试 - 多次随机数据验证
# ─────────────────────────────────────────────
cat("--- 第 11 轮: 重复性测试 (20 次) ---\n")
repeat_pass <- 0
repeat_fail <- 0
for (trial in 1:20) {
  set.seed(trial * 77)
  n <- sample(6:15, 1)
  sig_t <- sort(runif(n, 50, 50000))
  conc_t <- sort(runif(n, 0.1, 500))
  unk_t <- runif(3, 10, 60000)

  res_t <- tryCatch(
    msd_to_concentration(sig_t, conc_t, unk_t, method = "linear"),
    error = function(e) NULL
  )
  if (!is.null(res_t) && all(is.finite(res_t$concentrations))) {
    repeat_pass <- repeat_pass + 1
  } else {
    repeat_fail <- repeat_fail + 1
  }
}
assert(repeat_pass == 20, sprintf("20 次重复线性拟合全部成功 (%d/20)", repeat_pass))
cat("\n")

# ─────────────────────────────────────────────
# 第 12 轮: 指数拟合
# ─────────────────────────────────────────────
cat("--- 第 12 轮: 指数拟合 ---\n")
set.seed(1200)
sig12 <- seq(1, 10, length.out = 8)
conc12 <- 2 * exp(0.5 * sig12) + rnorm(8, 0, 1)
unk12 <- c(3, 7)

res12 <- msd_to_concentration(sig12, conc12, unk12, method = "exponential")
assert(length(res12$concentrations) == 2, "指数拟合返回 2 个浓度值")
assert(all(is.finite(res12$concentrations)), "指数拟合浓度值有限")
assert(all(res12$concentrations > 0), "指数拟合浓度为正")
cat("\n")

# ─────────────────────────────────────────────
# 第 13 轮: NA / 非有限值鲁棒性
# ─────────────────────────────────────────────
cat("--- 第 13 轮: NA 与非有限值处理 ---\n")
set.seed(1300)
sig13 <- c(100, 500, NA, 2000, Inf, 10000, 20000)
conc13 <- c(1, 5, 10, NA, 100, 200, 400)
unk13 <- c(800, NA, 5000)

res13 <- suppressWarnings(msd_to_concentration(sig13, conc13, unk13, method = "linear"))
assert(length(res13$concentrations) == 3, "含 NA 校准数据仍可拟合")
assert(all(is.finite(res13$concentrations[c(1, 3)])), "有效样本浓度有限")
assert(is.na(res13$concentrations[2]), "NA 信号预测浓度为 NA")
assert(res13$status[2] == "unknown", "NA 信号标记为 unknown")
cat("\n")

res13b <- msd_to_concentration(c(1, 2, 3, 4), c(1, 2, 3, 4),
                               c(NA, Inf, 2), method = "log")
assert(identical(res13b$status, c("unknown", "unknown", "in_range")),
       "log 方法的 NA/Inf 未知信号标记为 unknown")
res13c <- msd_to_concentration(c(1, 2, 3, 4), c(1, 2, 3, 4),
                               c(NA, Inf, 2), method = "power")
assert(identical(res13c$status, c("unknown", "unknown", "in_range")),
       "power 方法的 NA/Inf 未知信号标记为 unknown")

# ─────────────────────────────────────────────
# 第 14 轮: 绘图与保存
# ─────────────────────────────────────────────
cat("--- 第 14 轮: 绘图与保存 ---\n")
set.seed(1400)
sig14 <- seq(100, 50000, length.out = 10)
conc14 <- 200 + (0.5 - 200) / (1 + (sig14 / 8000)^1.3) + rnorm(10, 0, 0.5)
conc14 <- pmax(conc14, 0.01)
unk14 <- c(500, 5000, 30000)

tmpdir <- tempdir()
f_png <- file.path(tmpdir, "test_calib.png")
f_pdf <- file.path(tmpdir, "test_calib.pdf")
f_panel <- file.path(tmpdir, "test_panel.pdf")

ev14 <- save_plot(f_png,
  plot_calibration(sig14, conc14, unk14, method = "linear"),
  width = 7, height = 5, dpi = 150)
assert(file.exists(ev14), "单曲线 PNG 导出成功")
assert(file.size(ev14) > 1000, "PNG 文件大小合理")

ev14b <- save_plot(f_pdf,
  plot_calibration(sig14, conc14, unk14, method = "linear", log_x = TRUE),
  width = 7, height = 5)
assert(file.exists(ev14b), "对数轴 PDF 导出成功")

calib14 <- data.frame(
  factor = rep(c("A", "B", "C"), each = 8),
  signal = c(sig14[1:8], sig14[1:8] * 2, sig14[1:8] * 0.5),
  concentration = rep(seq(1, 100, length.out = 8), 3)
)
samples14 <- data.frame(factor = c("A", "B", "C"), signal = c(5000, 8000, 2000))

ev14c <- save_plot(f_panel,
  plot_multi_factor(calib14, samples14,
    method = list("A" = "linear", "B" = "quadratic", "C" = "log"), ncols = 2),
  width = 9, height = 7)
assert(file.exists(ev14c), "多因子面板图导出成功")

err14 <- tryCatch(save_plot(file.path(tmpdir, "test.xyz"), plot(1)),
                  error = function(e) e$message)
assert(grepl("Unsupported", err14), "不支持的格式触发错误")
cat("\n")

# ─────────────────────────────────────────────
# 第 15 轮: 非法输入与错误信息质量
# ─────────────────────────────────────────────
cat("--- 第 15 轮: 错误处理 ---\n")
err15a <- tryCatch(msd_to_concentration(c(1,2,3), c(1,2), c(1), method = "linear"),
                   error = function(e) e$message)
assert(grepl("same length", err15a), "长度不匹配报错明确")

err15b <- tryCatch(msd_to_concentration(c(1,2,3), c(1,2,3), c(1), method = "nope"),
                   error = function(e) e$message)
assert(grepl("Unknown method", err15b), "未知方法报错并提示可用方法")

err15c <- tryCatch(
  msd_multi_factor(data.frame(a = 1), data.frame(b = 1)),
  error = function(e) e$message)
assert(grepl("must contain columns", err15c), "缺列报错明确")

res15d <- suppressWarnings(msd_multi_factor(
  data.frame(factor = "X", signal = 100, concentration = 1),
  data.frame(factor = "X", signal = 50)))
assert(res15d$status == "error", "单点校准标记为 error 且不中断")
cat("\n")

# ─────────────────────────────────────────────
# 第 16 轮: 线性轴模式绘图阴影坐标一致性
# ─────────────────────────────────────────────
cat("--- 第 16 轮: 线性轴绘图 (log_x=FALSE) ---\n")
set.seed(1600)
sig16 <- seq(100, 50000, length.out = 10)
conc16 <- seq(1, 200, length.out = 10) + rnorm(10, 0, 1)
unk16 <- c(5000, 30000)

f16 <- file.path(tempdir(), "test_linear_axis.png")
err16 <- tryCatch({
  save_plot(f16,
    plot_calibration(sig16, conc16, unk16, method = "linear",
                     log_x = FALSE, show_limits = TRUE),
    width = 7, height = 5, dpi = 100)
  NULL
}, error = function(e) e$message)
assert(is.null(err16) && file.exists(f16), "线性轴 + 阴影绘制不报错")
cat("\n")

# ─────────────────────────────────────────────
# 第 17 轮: 4PL 难收敛数据 — 未收敛模型不被静默接受
# ─────────────────────────────────────────────
cat("--- 第 17 轮: 4PL 收敛校验 ---\n")
set.seed(1700)
sig17 <- c(100, 200, 300, 400, 500)
conc17 <- c(1, 100, 2, 99, 3)

res17 <- tryCatch(
  msd_to_concentration(sig17, conc17, c(250), method = "4pl"),
  error = function(e) e$message)
assert(is.character(res17) && grepl("converge", res17),
       "无法收敛的数据触发明确错误而非静默输出")
cat("\n")

# ─────────────────────────────────────────────
# 第 18 轮: eval_fit 含 NA 时的行为
# ─────────────────────────────────────────────
cat("--- 第 18 轮: eval_fit NA 处理 ---\n")
sig18 <- c(100, 500, NA, 2000, 10000)
conc18 <- c(1, 5, 10, NA, 100)
ev18 <- eval_fit(sig18, conc18, method = "linear")
assert(!is.na(ev18$R2), "eval_fit 含 NA 仍返回有效 R2")
assert(ev18$n_points == 3, "eval_fit 正确报告有效点数")
cat("\n")

# ─────────────────────────────────────────────
# 第 19 轮: 反向曲线（信号随浓度递减）
# ─────────────────────────────────────────────
cat("--- 第 19 轮: 反向曲线判定 ---\n")
sig19 <- c(1000, 800, 600, 400, 200, 100)
conc19 <- c(1, 5, 10, 50, 100, 500)
unk19 <- c(900, 300, 50)

res19 <- msd_to_concentration(sig19, conc19, unk19, method = "linear")
lim19 <- compute_limits(sig19, conc19)
assert(lim19$direction == "decreasing", "反向曲线方向为 decreasing")
assert(res19$status[1] == "in_range", "反向: 信号 900 在 100-1000 范围内")
assert(res19$status[2] == "in_range", "反向: 中间信号 300 在范围内")
assert(res19$status[3] == "above_ULOQ", "反向: 低信号 50 < signal_lo=100 对应高浓度 -> above_ULOQ")

unk19b <- c(1200, 80)
res19b <- msd_to_concentration(sig19, conc19, unk19b, method = "linear")
assert(res19b$status[1] == "below_LLOQ", "反向: 信号 1200 > signal_hi=1000 对应低浓度 -> below_LLOQ")
assert(res19b$status[2] == "above_ULOQ", "反向: 信号 80 < signal_lo=100 对应高浓度 -> above_ULOQ")
cat("\n")

# ─────────────────────────────────────────────
# 第 20 轮: Inf/-Inf 样本标记为 unknown
# ─────────────────────────────────────────────
cat("--- 第 20 轮: Inf/-Inf 样本处理 ---\n")
sig20 <- seq(100, 10000, length.out = 6)
conc20 <- seq(1, 100, length.out = 6)
unk20 <- c(500, Inf, -Inf, NA)

res20 <- msd_to_concentration(sig20, conc20, unk20, method = "linear")
assert(res20$status[1] == "in_range", "正常信号 in_range")
assert(res20$status[2] == "unknown", "Inf 标记为 unknown")
assert(res20$status[3] == "unknown", "-Inf 标记为 unknown")
assert(res20$status[4] == "unknown", "NA 标记为 unknown")
cat("\n")

# ─────────────────────────────────────────────
# 第 21 轮: summarize_limits 含 NA 状态
# ─────────────────────────────────────────────
cat("--- 第 21 轮: summarize_limits NA 状态 ---\n")
fake_result <- data.frame(
  factor = rep("X", 4),
  signal = c(1, 2, 3, 4),
  concentration = c(1, 2, 3, 4),
  status = c("in_range", NA, "below_LLOQ", "in_range"),
  LLOQ = rep(1, 4), ULOQ = rep(100, 4),
  method = rep("linear", 4),
  stringsAsFactors = FALSE
)
lim21 <- summarize_limits(fake_result)
assert(!is.na(lim21$n_in_range), "NA 状态不影响 n_in_range")
assert(lim21$n_in_range == 2, "n_in_range = 2")
assert(lim21$n_below == 1, "n_below = 1")
assert(lim21$n_unknown == 1, "NA 归入 n_unknown")
cat("\n")

# ─────────────────────────────────────────────
# 第 22 轮: 外推样本不撑爆 y 轴
# ─────────────────────────────────────────────
cat("--- 第 22 轮: 外推样本裁剪 ---\n")
set.seed(2200)
sig22 <- seq(100, 500, length.out = 6)
conc22 <- seq(1, 50, length.out = 6)
unk22 <- c(300, 1e9)

f22 <- file.path(tempdir(), "test_clip.png")
err22 <- tryCatch({
  save_plot(f22,
    plot_calibration(sig22, conc22, unk22, method = "linear", log_x = FALSE),
    width = 6, height = 4, dpi = 72)
  NULL
}, error = function(e) e$message)
assert(is.null(err22) && file.exists(f22), "外推样本不导致绘图错误")
cat("\n")

# ─────────────────────────────────────────────
# 第 23 轮: 对数轴绘图与拟合数据一致
# ─────────────────────────────────────────────
cat("--- 第 23 轮: 对数轴数据一致性 ---\n")
sig23 <- c(0, 100, 200, 300, 400, 500)
conc23 <- c(1, 5, 10, 20, 30, 50)

w23 <- tryCatch({
  f23 <- file.path(tempdir(), "test_logx_consistency.png")
  save_plot(f23,
    plot_calibration(sig23, conc23, method = "linear", log_x = TRUE),
    width = 6, height = 4, dpi = 72)
  NULL
}, warning = function(w) conditionMessage(w), error = function(e) e$message)
assert(!is.null(w23) && grepl("non-positive", w23), "含零信号时给出警告")
cat("\n")

# ─────────────────────────────────────────────
# 第 24 轮: compare_methods 输出含 error 列
# ─────────────────────────────────────────────
cat("--- 第 24 轮: compare_methods error 列 ---\n")
sig24 <- c(100, 200, 300)
conc24 <- c(1, 2, 3)
cmp24 <- compare_methods(sig24, conc24, methods = c("linear", "4pl"))
assert("error" %in% names(cmp24), "输出包含 error 列")
assert(cmp24$error[cmp24$method == "linear"] == "", "成功方法 error 为空")
assert(grepl("at least 5", cmp24$error[cmp24$method == "4pl"], ignore.case = TRUE),
       "4PL 点数不足时记录错误原因")
cat("\n")

# ─────────────────────────────────────────────
# 第 25 轮: plot_multi_factor 跳过无校准因子时发出警告
# ─────────────────────────────────────────────
cat("--- 第 25 轮: 跳过无校准因子警告 ---\n")
calib25 <- data.frame(factor = rep("A", 5), signal = seq(100, 500, length.out = 5),
                      concentration = 1:5)
samples25 <- data.frame(factor = c("A", "ZZZ"), signal = c(200, 300))
w25 <- tryCatch({
  save_plot(file.path(tempdir(), "test_skip.pdf"),
    plot_multi_factor(calib25, samples25, ncols = 1),
    width = 5, height = 4)
  NULL
}, warning = function(w) conditionMessage(w), error = function(e) e$message)
assert(!is.null(w25) && grepl("ZZZ", w25), "跳过无校准因子时警告包含因子名")
cat("\n")

# ─────────────────────────────────────────────
# 第 26 轮: 平台注册表
# ─────────────────────────────────────────────
cat("--- 第 26 轮: 平台注册表 ---\n")
plats <- list_platforms()
assert(length(plats) >= 6, "至少支持 6 个平台")
assert("luminex" %in% plats, "Luminex 平台已注册")
assert("elisa" %in% plats, "ELISA 平台已注册")
assert("simoa" %in% plats, "Simoa 平台已注册")
assert("cba" %in% plats, "CBA 平台已注册")

lp <- get_platform("luminex")
assert(lp$signal_label == "MFI", "Luminex 信号标签为 MFI")
assert("5pl" %in% lp$recommended_methods, "Luminex 推荐方法包含 5pl")

ep <- get_platform("elisa")
assert(ep$signal_label == "OD", "ELISA 信号标签为 OD")
assert(ep$typical_range[2] <= 5, "ELISA 典型信号范围合理")

err_plat <- tryCatch(get_platform("nonexist"), error = function(e) e$message)
assert(grepl("Unknown platform", err_plat), "未知平台报错明确")
cat("\n")

# ─────────────────────────────────────────────
# 第 27 轮: 自定义平台注册
# ─────────────────────────────────────────────
cat("--- 第 27 轮: 自定义平台注册 ---\n")
register_platform("myassay", "My Custom Assay", "RFU",
                  signal_unit = "Relative Fluorescence Units",
                  typical_range = c(50, 100000),
                  recommended_methods = c("4pl", "log"))
mp <- get_platform("myassay")
assert(mp$signal_label == "RFU", "自定义平台信号标签正确")
assert("myassay" %in% list_platforms(), "自定义平台出现在列表中")
cat("\n")

# ─────────────────────────────────────────────
# 第 28 轮: Luminex 多因子转换
# ─────────────────────────────────────────────
cat("--- 第 28 轮: Luminex 多因子转换 ---\n")
set.seed(2800)
n_lx <- 10
lx_sig <- seq(20, 20000, length.out = n_lx)
lx_conc <- 500 + (0.5 - 500) / (1 + (lx_sig / 3000)^1.2) + rnorm(n_lx, 0, 1)
lx_conc <- pmax(lx_conc, 0.01)

calib_lx <- data.frame(
  factor = rep(c("IL-2", "IL-10"), each = n_lx),
  signal = c(lx_sig, lx_sig * 1.5),
  concentration = rep(lx_conc, 2)
)
samples_lx <- data.frame(
  factor = c("IL-2", "IL-2", "IL-10"),
  signal = c(100, 5000, 25000)
)

result_lx <- msd_multi_factor(calib_lx, samples_lx,
  method = list("IL-2" = "5pl", "IL-10" = "4pl"), platform = "luminex")

assert(nrow(result_lx) == 3, "Luminex 多因子返回 3 行")
assert(all(result_lx$platform == "Luminex"), "平台列正确标记为 Luminex")
assert(all(c("LLOQ", "ULOQ", "status", "method") %in% names(result_lx)),
       "Luminex 输出包含必要列")
cat("\n")

# ─────────────────────────────────────────────
# 第 29 轮: ELISA 单因子转换
# ─────────────────────────────────────────────
cat("--- 第 29 轮: ELISA 单因子转换 ---\n")
set.seed(2900)
od_sig <- seq(0.05, 3.5, length.out = 8)
od_conc <- 100 + (0.1 - 100) / (1 + (od_sig / 1.5)^2) + rnorm(8, 0, 0.5)
od_conc <- pmax(od_conc, 0.01)

res_elisa <- convert_signal(od_sig, od_conc, unknown = c(0.5, 2.0),
                            platform = "elisa")
assert(res_elisa$platform == "ELISA", "ELISA 平台标识正确")
assert(length(res_elisa$concentrations) == 2, "ELISA 返回 2 个浓度")
assert(all(is.finite(res_elisa$concentrations)), "ELISA 浓度有限")
cat("\n")

# ─────────────────────────────────────────────
# 第 30 轮: 多平台批量转换
# ─────────────────────────────────────────────
cat("--- 第 30 轮: 多平台批量转换 ---\n")
set.seed(3000)
calib_msd <- data.frame(
  factor = rep("TNF-a", 8),
  signal = seq(100, 50000, length.out = 8),
  concentration = c(1, 5, 10, 50, 100, 250, 500, 1000)
)
samples_msd <- data.frame(factor = "TNF-a", signal = c(5000, 30000))

calib_lum <- data.frame(
  factor = rep("IFN-g", 8),
  signal = seq(30, 18000, length.out = 8),
  concentration = c(0.5, 2, 5, 20, 50, 100, 200, 500)
)
samples_lum <- data.frame(factor = "IFN-g", signal = c(2000, 15000))

batch_result <- multi_platform_convert(list(
  list(calib_data = calib_msd, sample_data = samples_msd,
       platform = "msd", method = "linear"),
  list(calib_data = calib_lum, sample_data = samples_lum,
       platform = "luminex", method = "log")
))

assert(nrow(batch_result) == 4, "多平台批量返回 4 行")
assert(all(c("msd", "luminex") %in% tolower(batch_result$platform)),
       "多平台结果包含两种平台")
assert(all(!is.na(batch_result$concentration)), "多平台浓度均有效")
cat("\n")

# ─────────────────────────────────────────────
# 第 31 轮: Simoa 平台转换
# ─────────────────────────────────────────────
cat("--- 第 31 轮: Simoa 平台转换 ---\n")
set.seed(3100)
simoa_sig <- seq(0.01, 50, length.out = 8)
simoa_conc <- c(0.01, 0.05, 0.1, 0.5, 1, 5, 10, 50)

res_simoa <- msd_to_concentration(simoa_sig, simoa_conc,
                                  unknown = c(0.5, 10), method = "linear",
                                  platform = "simoa")
assert(res_simoa$platform == "Simoa", "Simoa 平台标识正确")
assert(all(is.finite(res_simoa$concentrations)), "Simoa 浓度有限")
cat("\n")

# ─────────────────────────────────────────────
# 第 32 轮: convert_signal 自动方法选择
# ─────────────────────────────────────────────
cat("--- 第 32 轮: convert_signal 自动方法选择 ---\n")
set.seed(3200)
auto_sig <- seq(20, 20000, length.out = 10)
auto_conc <- 300 + (0.5 - 300) / (1 + (auto_sig / 3000)^1.3) + rnorm(10, 0, 0.5)
auto_conc <- pmax(auto_conc, 0.01)

res_auto <- suppressMessages(
  convert_signal(auto_sig, auto_conc, unknown = c(1000), platform = "luminex"))
assert(res_auto$method == "5pl", "Luminex 自动选择 5pl 方法")
assert(is.finite(res_auto$concentrations), "自动选择方法预测正常")
cat("\n")

# ─────────────────────────────────────────────
# 第 33 轮: 外推样本 y 轴几何断言（问题 B 根因验证）
# ─────────────────────────────────────────────
cat("--- 第 33 轮: 外推样本不撑爆 y 轴（几何断言） ---\n")
set.seed(3300)
sig33 <- seq(100, 500, length.out = 6)
conc33 <- seq(1, 5, length.out = 6)
unk33 <- c(300, 1e9)

f33 <- file.path(tempdir(), "test_ylim.png")
grDevices::png(f33, width = 6, height = 4, units = "in", res = 72)
ev33 <- plot_calibration(sig33, conc33, unk33, method = "linear")
invisible(grDevices::dev.off())

assert(!is.null(ev33$ylim), "plot_calibration 返回 ylim")
assert(diff(ev33$ylim) < 10,
       sprintf("ylim 跨度与校准浓度同量级 (diff=%.2f)", diff(ev33$ylim)))
assert(ev33$n_clipped >= 1, sprintf("1e9 外推样本被裁剪 (n_clipped=%d)", ev33$n_clipped))
cat("\n")

# ─────────────────────────────────────────────
# 第 34 轮: convert_signal 点数不足时方法回退
# ─────────────────────────────────────────────
cat("--- 第 34 轮: convert_signal 方法回退 ---\n")
set.seed(3400)
sig34 <- seq(30, 18000, length.out = 5)
conc34 <- 200 + (0.5 - 200) / (1 + (sig34 / 3000)^1.3) + rnorm(5, 0, 0.3)
conc34 <- pmax(conc34, 0.01)

res34 <- suppressMessages(
  convert_signal(sig34, conc34, unknown = c(2000), platform = "luminex"))
assert(res34$method != "5pl", "5 点数据回退避开 5pl（需 6 点）")
assert(is.finite(res34$concentrations), "回退后方法预测正常")

sig34b <- seq(0.05, 3.0, length.out = 4)
conc34b <- c(1, 10, 50, 100) + rnorm(4, 0, 0.5)
res34b <- suppressMessages(
  tryCatch(convert_signal(sig34b, conc34b, unknown = c(1.0), platform = "elisa"),
           error = function(e) NULL))
assert(!is.null(res34b) && is.finite(res34b$concentrations),
       "ELISA 4 点数据回退到可用方法")
cat("\n")

# ─────────────────────────────────────────────
# 第 35 轮: multi_platform_convert 失败批次容错
# ─────────────────────────────────────────────
cat("--- 第 35 轮: 多平台批处理失败容错 ---\n")
set.seed(3500)
good_calib <- data.frame(factor = rep("A", 6),
                          signal = seq(100, 10000, length.out = 6),
                          concentration = 1:6)
good_samples <- data.frame(factor = "A", signal = 3000)

bad_batch <- list(calib_data = good_calib, sample_data = good_samples,
                  platform = "no_such_platform")
good_batch <- list(calib_data = good_calib, sample_data = good_samples,
                   platform = "msd", method = "linear")

res35 <- suppressWarnings(multi_platform_convert(list(bad = bad_batch, good = good_batch)))
assert("batch" %in% names(res35), "输出包含 batch 列")
assert(any(res35$batch == "bad" & res35$status == "error"), "失败批次标记为 error 且可识别")
assert(any(res35$batch == "good" & res35$status != "error"), "正常批次仍返回结果")

res35b <- suppressWarnings(multi_platform_convert(list(x = list(platform = "msd"))))
assert(any(res35b$status == "error"), "缺字段批次标记为 error 而非中断")
cat("\n")

# ─────────────────────────────────────────────
# 第 36 轮: register_platform 输入校验
# ─────────────────────────────────────────────
cat("--- 第 36 轮: register_platform 输入校验 ---\n")
err36a <- tryCatch(register_platform("", "X", "S"), error = function(e) e$message)
assert(grepl("name", err36a, ignore.case = TRUE), "空名称被拒绝")

err36b <- tryCatch(register_platform("p2", "X", "S", typical_range = c(1, 2, 3)),
                   error = function(e) e$message)
assert(grepl("typical_range", err36b), "typical_range 长度非 2 被拒绝")

err36c <- tryCatch(register_platform("p3", "X", "S", typical_range = c(10, 1)),
                   error = function(e) e$message)
assert(grepl("typical_range", err36c), "typical_range min>=max 被拒绝")

err36d <- tryCatch(register_platform("p4", "X", ""),
                   error = function(e) e$message)
assert(grepl("signal_label", err36d), "空 signal_label 被拒绝")
cat("\n")

# ─────────────────────────────────────────────
# 第 37 轮: typical_range 平台错配警告
# ─────────────────────────────────────────────
cat("--- 第 37 轮: typical_range 校验 ---\n")
w37 <- tryCatch({
  msd_to_concentration(c(1e9, 2e9, 3e9), c(1, 2, 3), 1.5e9,
                       method = "linear", platform = "elisa")
  NULL
}, warning = function(w) conditionMessage(w), error = function(e) e$message)
assert(!is.null(w37) && grepl("typical range", w37), "ELISA 传入 1e9 量级信号触发警告")
cat("\n")

# ─────────────────────────────────────────────
# 第 38 轮: 非单调校准曲线使用完整观测范围
# ─────────────────────────────────────────────
cat("--- 第 38 轮: 非单调校准信号范围 ---\n")
w38 <- tryCatch(compute_limits(c(100, 5000, 300, 1000), 1:4),
                warning = function(w) conditionMessage(w))
lim38 <- suppressWarnings(compute_limits(c(100, 5000, 300, 1000), 1:4))
assert(grepl("not monotonic", w38), "非单调校准信号给出警告")
assert(lim38$signal_lo == 100 && lim38$signal_hi == 5000,
       "检测范围包含全部校准信号")
err38 <- tryCatch(compute_limits(numeric(), numeric()), error = function(e) e$message)
assert(grepl("At least 2", err38), "空检测限输入明确报错")
err38b <- tryCatch(compute_limits(c(1, 2), c(5, 5)), error = function(e) e$message)
assert(grepl("distinct concentration", err38b), "相同浓度无法定义检测限并明确报错")
cat("\n")

# ─────────────────────────────────────────────
# 第 39 轮: 自定义平台注册方法校验与自动回退
# ─────────────────────────────────────────────
cat("--- 第 39 轮: 平台方法校验与自动回退 ---\n")
err39 <- tryCatch(register_platform("badmethod", "Bad", "RFU",
                                    recommended_methods = "missing"),
                  error = function(e) e$message)
assert(grepl("Unknown recommended method", err39), "未知推荐方法在注册时拒绝")
register_platform("fallback", "Fallback", "RFU",
                  recommended_methods = c("5pl", "linear"))
res39 <- suppressMessages(convert_signal(1:5, 1:5, 3, platform = "fallback"))
assert(res39$method == "linear", "拟合失败后自动尝试次选推荐方法")
cat("\n")

# ─────────────────────────────────────────────
# 第 40 轮: 退化数据、NA 因子和非法方法输入
# ─────────────────────────────────────────────
cat("--- 第 40 轮: 输入边界校验 ---\n")
err40a <- tryCatch(msd_multi_factor(
  data.frame(factor = "A", signal = 1:3, concentration = 1:3),
  data.frame(factor = NA_character_, signal = 2)), error = function(e) e$message)
assert(grepl("factor must not contain NA", err40a), "NA 因子值明确拒绝")
err40b <- tryCatch(msd_multi_factor(
  data.frame(factor = "A", signal = 1:3, concentration = 1:3),
  data.frame(factor = "A", signal = 2), method = c("linear", "log")),
  error = function(e) e$message)
assert(grepl("single method", err40b), "多元素方法向量明确拒绝")
err40c <- tryCatch(fit_4pl(c(-1, 1, 2, 3, 4), 1:5), error = function(e) e$message)
assert(grepl("> 0", err40c), "4PL 非正信号明确拒绝")
cat("\n")

# ─────────────────────────────────────────────
# 第 41 轮: 多平台自动回退和绘图非法未知值
# ─────────────────────────────────────────────
cat("--- 第 41 轮: 多平台回退与绘图边界 ---\n")
auto_batch <- multi_platform_convert(list(list(
  calib_data = data.frame(factor = rep("A", 5), signal = 1:5, concentration = 1:5),
  sample_data = data.frame(factor = "A", signal = 3), platform = "luminex")))
assert(nrow(auto_batch) == 1 && auto_batch$status != "error" &&
       auto_batch$method %in% c("4pl", "cubic"),
       "多平台自动方法流程尝试可用后备模型")
w41 <- tryCatch({
  save_plot(file.path(tempdir(), "test_plot_invalid_unknown.pdf"),
    plot_calibration(1:5, 1:5, c(-1, 3), method = "log"))
  NULL
}, warning = function(w) conditionMessage(w), error = function(e) e$message)
assert(!is.null(w41) && grepl("non-positive unknown", w41),
       "对数模型绘图剔除非正未知信号并警告")
cat("\n")

# ─────────────────────────────────────────────
# 汇总
# ─────────────────────────────────────────────
cat("========================================\n")
cat(sprintf(" 测试完成: %d 通过, %d 失败, 共 %d 项\n",
            pass_count, fail_count, pass_count + fail_count))
cat("========================================\n")

if (fail_count > 0) {
  cat("存在失败项，请检查!\n")
  quit(status = 1)
} else {
  cat("全部通过!\n")
}
