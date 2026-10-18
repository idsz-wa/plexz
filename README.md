# plexz

多平台免疫分析信号值 → 浓度值转换工具包。支持 MSD、Luminex、ELISA、Simoa、CBA、ELLA 等检测平台，提供多种标准曲线拟合方法、多因子批量处理、检测上下限（LLOQ/ULOQ）判定，以及美化后的标准曲线绘图与导出。

## 安装

```r
# 方式一：直接 source（无需安装，适合脚本调用）
pkg_dir <- "plexz/R"
for (f in sort(list.files(pkg_dir, full.names = TRUE))) source(f)

# 方式二：作为 R 包安装
# 在终端执行：
# R CMD INSTALL plexz/
# 然后在 R 中：
# library(plexz)
```

仅依赖 R 基础包（stats / graphics / grDevices / tools），无需安装第三方依赖。

## 支持的检测平台

| 平台 | 全称 | 信号类型 | 典型范围 | 推荐方法 |
|------|------|----------|----------|----------|
| `msd` | Meso Scale Discovery | Signal (ECL) | 10 – 200,000 | 4pl, 5pl, log |
| `luminex` | Luminex xMAP | MFI | 10 – 25,000 | 5pl, 4pl, cubic |
| `elisa` | ELISA | OD | 0.01 – 4.0 | 4pl, log, quadratic |
| `simoa` | Quanterix Simoa | AEB | 0.001 – 100 | 4pl, linear, log |
| `cba` | BD Cytometric Bead Array | MFI | 5 – 50,000 | 4pl, 5pl, log |
| `ella` | ProteinSimple ELLA | Signal | 100 – 500,000 | 4pl, 5pl, log |

查看全部平台：`list_platforms()`

注册自定义平台：

```r
register_platform("myassay", "My Custom Assay", "RFU",
                  typical_range = c(50, 100000),
                  recommended_methods = c("4pl", "log"))
```

## 支持的拟合方法

| 方法名 | 模型 | 适用场景 |
|--------|------|----------|
| `linear` | y = a + b·x | 信号与浓度线性相关 |
| `quadratic` | y = a + b·x + c·x² | 轻度弯曲的曲线 |
| `cubic` | 三次多项式 | 弯曲较明显的曲线 |
| `log` | y = a + b·ln(x) | 对数关系（信号需为正） |
| `power` | y = a·x^b | 幂函数关系（信号和浓度均需为正） |
| `exponential` | y = a·exp(b·x) | 指数增长关系（浓度需为正） |
| `4pl` | 四参数逻辑斯蒂 | 典型免疫分析 S 型曲线（至少 5 个点） |
| `5pl` | 五参数逻辑斯蒂 | 非对称 S 型曲线（至少 6 个点） |

查看全部方法：`list_methods()`

## 快速开始

### 1. 单因子转换（指定平台）

```r
# MSD 平台示例
signal <- c(100, 500, 1000, 5000, 10000, 25000, 50000)
conc   <- c(1, 5, 10, 50, 100, 250, 500)
unknown <- c(2000, 30000, 80000)

res <- msd_to_concentration(signal, conc, unknown, method = "4pl", platform = "msd")
res$concentrations   # 预测浓度
res$status           # "in_range" / "below_LLOQ" / "above_ULOQ" / "unknown"
res$LLOQ             # 定量下限
res$ULOQ             # 定量上限
res$platform         # "MSD"
```

### 2. 自动方法选择（按平台推荐）

```r
# Luminex 平台：自动选择 5pl
res <- convert_signal(signal, conc, unknown, platform = "luminex")
# 或手动指定方法
res <- convert_signal(signal, conc, unknown, platform = "elisa", method = "4pl")
```

### 3. 多因子批量转换

```r
calib_data <- data.frame(
  factor = rep(c("IL-6", "TNF-a"), each = 7),
  signal = c(signal, signal * 1.2),
  concentration = rep(conc, 2)
)
sample_data <- data.frame(
  factor = c("IL-6", "IL-6", "TNF-a"),
  signal = c(2000, 30000, 60000)
)

result <- msd_multi_factor(calib_data, sample_data,
  method = list("IL-6" = "4pl", "TNF-a" = "linear"),
  platform = "luminex")

print(result)
summarize_limits(result)
```

### 4. 多平台批量转换

```r
batch_result <- multi_platform_convert(list(
  list(calib_data = calib_msd, sample_data = samples_msd,
       platform = "msd", method = "4pl"),
  list(calib_data = calib_lum, sample_data = samples_lum,
       platform = "luminex", method = "5pl"),
  list(calib_data = calib_elisa, sample_data = samples_elisa,
       platform = "elisa")  # method 为 NULL 时自动选择平台推荐方法
))
```

### 5. 方法对比与自动选择

```r
cmp <- compare_methods(signal, conc,
  methods = c("linear", "quadratic", "log", "power", "4pl"))
print(cmp)  # 按 R2 降序，含 error 列记录失败原因

best <- cmp$method[1]
res <- msd_to_concentration(signal, conc, unknown, method = best)
```

### 6. 绘图

```r
# 单条标准曲线（轴标签自动使用平台信号名称）
plot_calibration(signal, conc, unknown, method = "4pl", platform = "luminex")
# X 轴显示为 "MFI (Luminex)"

# 对数轴
plot_calibration(signal, conc, unknown, method = "4pl",
                 platform = "msd", log_x = TRUE)

# 多因子面板图
plot_multi_factor(calib_data, sample_data,
  method = list("IL-6" = "4pl", "TNF-a" = "linear"),
  platform = "luminex", ncols = 2)

# 导出图片
save_plot("curve_luminex.pdf",
  plot_calibration(signal, conc, unknown, method = "5pl", platform = "luminex"),
  width = 7, height = 5)
```

## 输出状态说明

| status | 含义 |
|--------|------|
| `in_range` | 样本信号在标准曲线范围内，浓度可靠 |
| `below_LLOQ` | 信号低于最低标准品，浓度为外推值，建议标注 < LLOQ |
| `above_ULOQ` | 信号高于最高标准品，浓度为外推值，建议稀释后重测 |
| `unknown` | 样本信号为 NA 或非有限值（Inf/-Inf），无法判定 |
| `error` | 该因子拟合失败或标准品不足 3 个 |

## 健壮性说明

- 校准数据中的 `NA` 或非有限值会被自动剔除并给出警告
- 4PL/5PL 使用多组初始值策略并严格检查收敛状态（`convInfo$isConv`），未收敛的拟合不会被静默接受
- 拟合失败时多因子流程不会中断，该因子标记为 `error`，其余因子正常输出
- 对数 / 幂函数拟合会自动校验正值要求，非法输入给出明确报错
- 支持信号随浓度递减的反向曲线（如竞争法），自动检测方向并正确判定
- 超出绘图范围的样本点被裁剪，不会撑爆坐标轴
- 对数轴模式下非正信号点被剔除并发出警告，绘图与拟合数据一致
- `convert_signal()` 按平台推荐方法逐级回退：数据点数不足或拟合失败时自动尝试次选方法
- 校准信号超出平台 `typical_range` 过半时发出警告，提示检查平台是否选错
- `multi_platform_convert()` 单批失败不中断整体，失败批次标记 `status = "error"` 并带 `batch` 列，便于对账

## 关于检测限的说明

本包中 `LLOQ` / `ULOQ` 定义为标准曲线的最低 / 最高标准品浓度（即标准曲线范围）。这是常见的简化定义，适用于日常实验中的快速判定。若结果需进入正式报告或法规文件，建议按照行业规范（如 CV ≤ 20%、回收率 80–120%）重新计算真正的定量限。

## 运行测试

```bash
Rscript plexz/tests/test_plexz.R
```

测试包含 37 轮模拟：基本拟合、检测限判定、反向曲线、多因子批处理、边界情况、高噪声鲁棒性、4PL/5PL 收敛、重复性、绘图导出、NA/Inf 处理、外推裁剪几何断言、多平台转换、自定义平台注册与输入校验、自动方法选择与回退、批处理失败容错、平台范围校验等。

## 函数速查

| 函数 | 用途 |
|------|------|
| `msd_to_concentration()` | 单因子信号→浓度转换（指定平台） |
| `convert_signal()` | 自动选择平台推荐方法的转换 |
| `msd_multi_factor()` | 多因子批量转换（指定平台） |
| `multi_platform_convert()` | 多平台批量转换 |
| `compare_methods()` | 多方法 R2/RMSE 对比 |
| `eval_fit()` | 单一方法拟合质量评估 |
| `plot_calibration()` | 单条标准曲线绘图 |
| `plot_multi_factor()` | 多因子面板绘图 |
| `save_plot()` | 导出图片文件 |
| `summarize_limits()` | 检测限与样本状态汇总 |
| `list_methods()` | 列出全部拟合方法 |
| `list_platforms()` | 列出全部支持平台 |
| `get_platform()` | 获取平台配置 |
| `register_platform()` | 注册自定义平台 |