# c=24.9 自包含证明与认证包

作者：Cheng Huang  
单位：Shanghai Jiaotong University

本目录给出以下无条件结论的完整复现材料：对每个满足
`log log N >= 24.9` 的偶数 `N`，有 `D_{1,2}(N) > 0`。

这不是仅含最终 PDF 的归档。包内包括证明源码、34 个 Arb/ACB 区间积分的
计算核心、136 项结构断言、唯一正向总账、两篇外部定理原始 TeX、原审计输入、
最终 PDF 及其重新编译和渲染检查程序。

## 关键依赖已经内置

最容易漏掉的传递依赖在包内的固定位置如下：

- `c249_interval_certificate.py` 直接导入根目录的
  `c26_interval_certificate.py`；后者是 34 个 Arb/ACB 区间积分的计算核心。
- `c249_verify.py` 审查
  `bjs_source/Explicit_Chen_-_New.tex`。
- `c249_verify.py` 审查
  `wu_source/ChenDoubleSieve1_Paper.tex`。

完整机器可读依赖图见 `dependency_graph.json`。`source_manifest.sha256` 固定所有
不可变输入；`verify_bundle.py` 只用 Python 标准库即可先检查包是否缺件或被改动。

## 环境

- Windows PowerShell 5.1 或 PowerShell 7；
- Python 3.12；
- `pdflatex`（本包认证时使用 TeX Live 2025）；
- `pdftoppm`（Poppler）；
- `c249_requirements.txt` 中锁定的 Python 包。

包内 `wheelhouse` 已带有 Windows / CPython 3.12 对应的锁定版 wheel。建议在解压
目录内建立独立虚拟环境，并离线安装：

```powershell
py -3.12 -m venv .venv
.\install_offline.ps1 -Python .\.venv\Scripts\python.exe
```

离线安装脚本等价于：

```powershell
.\.venv\Scripts\python.exe -m pip install --no-index `
  --find-links .\wheelhouse -r .\c249_requirements.txt
```

## 一键完整复现

在本目录执行：

```powershell
powershell -ExecutionPolicy Bypass -File .\reproduce.ps1 `
  -CorePython .\.venv\Scripts\python.exe
```

脚本依次完成：

1. 核对包内必需文件和 SHA-256 清单；
2. 从源码重新计算 34 个 Arb/ACB 区间积分；
3. 重新执行 136 项结构断言；
4. 重新生成三份 JSON 总账并逐项核对 TeX 表格；
5. 三遍编译 `explicit_chen_c249_final.tex`；
6. 渲染全部 27 页，检查缺页、空页和未解析引用；
7. 再次核对生成结果与源文件哈希。

若核心证明和 PDF 检查使用两个不同的现有 Python，可分别指定：

```powershell
.\reproduce.ps1 -CorePython C:\path\core-python.exe `
                 -PdfPython C:\path\pdf-python.exe
```

`-SkipPdfQa` 只用于机器没有 Poppler 时的核心算术复现；它不构成完整复现。

## 可分别运行的命令

```powershell
python .\verify_bundle.py --require-generated
python .\c249_verify.py --check-tex
pdflatex -interaction=nonstopmode -halt-on-error `
  -output-directory=output/pdf explicit_chen_c249_final.tex
python .\c249_pdf_qa.py
```

## 认证时应得到的关键输出

- Arb/ACB 区间积分数：`34`；
- 结构断言数：`136`；
- 有符号总账行数：`11`；
- 有限预算数：`7`；
- 主项：`[0.023536843646471438769041841946393 +/- 1.22e-31]`；
- 总预算：`[0.012005000000000000000000000000000 +/- 3e-38]`；
- 净余量：`[0.011531843646471438769041841946393 +/- 1.22e-31]`；
- 严格正下界：`0.0108`；
- PDF 页数：`27`。

最终证明声明及每项审计关闭状态见 `final_audit.md`。原始审计标准与被修订论文
分别保存在 `audit_inputs/reply.md` 和
`audit_inputs/explicit_chen_c2336_final.tex`，因此复现不依赖包外的 `D:` 盘路径。
本次打包、离线安装、目标目录运行和 ZIP 解压复测的记录见
`VALIDATION_REPORT.md`；完整运行输出保存在 `validation/full_reproduction.log`。

