# 七分之八 Chen 结果：唯一正式复现入口

正文：`qrh_chen_51046.tex`。结果和边界：`qrh_final_audit.md`。

## 从零运行

使用 Python 3.12.8 和 python-flint 0.9.0（本次核验环境）。进入文件所在目录：

```text
python qrh_chen_reproduce.py
python qrh_chen_interface_checks.py
python qrh_chen_positivity_search.py
```

第一条从头生成全部有限输入和起点，写 `qrh_chen_reproduction.json`。
第二条输出主起点51046的独立常数/量纲门限，写 `qrh_chen_interface_checks.json`。
第三条先核对第一条的输入哈希，再穷举已声明范围，写 `qrh_chen_positivity_search.json`。
执行第一条不需要任何预先生成的 JSON、旧 GRH 论文、旧 c249 程序或联网。
其余两条不能替代第一条。

## 必需源文件

```text
qrh_chen_reproduce.py
qrh_chen_certificate.py
qrh_chen_interface_checks.py
qrh_chen_positivity_search.py
conditional_chen_prefix_block_certificate.py
conditional_chen_rankin_certificate.py
conditional_chen_uniform_product_certificate.py
qrh_chen_51046.tex
qrh_chen_research_20261007.md
bjs_source/Explicit_Chen_-_New.tex
qrh_source/paper.tex
qrh_source/Nonvanishing.lean
```

`conditional_chen_rankin_certificate.py` 在正式路径中只提供完整整数素数筛；
文件名前缀不是 GRH 输入。中等素数质量和产品预算都在正式路径中重新生成。
三个引用源码/记录文件只作出处和哈希留档，不从它们提取旧数值预算。
`Nonvanishing.lean` 不是完整 Lean 工程，不能单靠这个文件验证外部零点定理。

本轮已在 `tmp/qrh_clean_final_20261007/` 中仅复制上述12个源文件，
从无 JSON 的目录运行全流程。Python 库由环境提供，不包含在这些源文件中。

## 数学认证边界

所有预算和积分由有理参数重新求值。完整数值入口使用192-bit Arb/ACB，
三段积分由 `acb.integral` 验证；见正文附录的分段、解析分支和无奇点说明。
有限 PASS 不证明外部零点定理、一般解析引理、形式化正确性或全局最优。
这些逻辑不能通过哈希或数值门限代替；正文给出新接口的书面证明及精确引用。

`qrh_chen_certificate.py` 中的 `inputs()` 和直接执行入口仍是历史两分支探索。
不要把 `python qrh_chen_certificate.py` 当作正式复现命令；其历史 JSON 不在新证书链中。
`qrh_chen_research_20261007.md` 和早期接口审计是时间顺序的研究记录，
其中较旧的51200目标已被当前正文51046取代。权威总账是正文式(15)及
`qrh_chen_reproduction.json` 的 `positivity_endpoint`，不是历史记录里的搜索结果。

## PDF

可用既有 TeX 环境对 `qrh_chen_51046.tex` 编译两遍。正文是独立文件，
不需要输入 BJS/Wu 的整份源码；这些文件只用于出处留档和复现哈希。
PDF 位于 `output/pdf/qrh_chen_51046.pdf`，个人署名和机构为空。
