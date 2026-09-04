# c=24.9 自包含包验证报告

验证日期：2026-09-05（Asia/Hong_Kong）

## 结论

**PASSED。** 本包已从目标目录和最终 ZIP 的全新解压目录分别执行。缺失依赖、
区间认证、结构认证、总账、TeX 编译和 PDF 渲染检查全部通过。

## 自包含性

- `c26_interval_certificate.py`：存在于包根目录，并被
  `c249_interval_certificate.py` 直接导入；
- `bjs_source/Explicit_Chen_-_New.tex`：存在并被总验证器哈希跟踪；
- `wu_source/ChenDoubleSieve1_Paper.tex`：存在并被总验证器哈希跟踪；
- 原审计输入和被修订原稿：均在 `audit_inputs/`；
- Python 锁定依赖：三个 wheel 均在 `wheelhouse/`，已在全新 Python 3.12
  虚拟环境中用 `--no-index` 成功安装。

## 从头重算结果

- 34 个 Arb/ACB 区间积分：PASSED；
- 136 项结构断言：PASSED；
- TeX 数值表与 11 行有符号总账逐项一致：PASSED；
- 7 项有限预算总和：PASSED；
- 主项：`[0.023536843646471438769041841946393 +/- 1.22e-31]`；
- 总预算：`[0.012005000000000000000000000000000 +/- 3e-38]`；
- 净余量：`[0.011531843646471438769041841946393 +/- 1.22e-31]`；
- 严格正下界：`0.0108`。

## 文档认证

- pdfTeX / TeX Live 2025，固定执行三遍：PASSED；
- 最终 PDF：27 页；
- 未解析引用、空白或异常稀疏页面：未发现；
- 以 5 张 contact sheet 覆盖第 1--27 页进行人工视觉检查：未发现文字裁切、
  重叠、破损表格、黑块或不可读字符；
- 完整运行日志：`validation/full_reproduction.log`。

## ZIP 独立复测

最终 ZIP 被解压到包外的全新目录后，执行：

```powershell
python verify_bundle.py --require-generated
python c249_verify.py --check-tex
reproduce.ps1 -CorePython <离线安装得到的 Python> `
              -PdfPython <同一 Python>
```

三项均返回成功；解压目录未使用原工作区文件。ZIP 的外置 SHA-256 校验文件与
ZIP 并列交付，以避免自包含文件对自身哈希的循环依赖。

