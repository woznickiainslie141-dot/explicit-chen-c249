# 七分之八零点界：最终结果与核查边界

日期：2026-10-07。当前正文为 `qrh_chen_51046.tex`，无个人署名、机构或 PDF 作者元数据。
本记录针对新零点界下的新证明，不替代旧 c=24.9 稿的逐项审计结论。

## 可以声称的结果

设 (Z) 为：所有 Dirichlet L 函数对所有模数和高度统一地在 Re(s)>7/8 无零点，
允许主字符在 1 的极点。正文证明了

\[
 (Z)\Longrightarrow
 D_{1,2}(N)>9\cdot10^{-7}\frac{U_NN}{(\log N)^2}>0
 \quad\text{对每个偶数 }N\ge e^{51046},
\]

其中

\[
 U_N=2e^{-\gamma}\prod_{p>2}\left(1-\frac1{(p-1)^2}\right)
       \prod_{\substack{p>2\\p\mid N}}\frac{p-1}{p-2}.
\]

这对应 `c=log(51046)=10.84048246596857468749...`，不是 `N>=exp(exp(51046))`。
采用外部新论文 Theorem 1.1 作为已证定理输入，就得到无条件的上述 Chen 阈值。
不需要 GRH，亦不把只有 zeta 的零点界误当作全部 Dirichlet 字符的界。

外部定理出处固定为 OpenAI 仓库提交
`adc7f1241b42e322a6451854ab7e4b4c146bf78a` 的
`preprints/The-Quasi-Riemann-Hypothesis-September-30-2026/build/paper.tex`。
本轮从原始网址取回字节核对：本地源码 SHA-256 为
`42a5ee0febca59fd1def55cfd6c6808c322ef4237d7726303ea52f711deac6a3`。
本项目没有重新证明这项外部重大定理，也没有在本机编译其整个 Lean 依赖链。
下载的 `Nonvanishing.lean` 只是出处留档，不是本机正式验证证明。

## 各接口的处理

下表中 CLOSED 指本文新增接口的书面证明及所需有限运算已补齐、核对；
不指该接口已形式化，也不表示已经有独立同行审查。

| 接口 | 证据与处理 | 状态 |
|---|---|---|
| 7/8 界导出有效字符和 | Lemma 2.1：主字符去极点、无零圆盘、BC/Cauchy 常数、平滑 Perron 移线、差商和部分求和全部写出；指数 61/64 | CLOSED under (Z) |
| 不暗用 RH 的素数质量、倒素数和、Mertens 产品 | Lemma 2.2：由新字符和界推导，包括常数匹配、Euler 尾项和 N 的大素因子删除 | CLOSED under (Z) |
| 产品条件不仅检查一个上端点 | §3.1：完整有限峰/谷有序合并；远端采用 BJS Lemmas 16--17、JS Lemma 3.2；含所有中间实端点 | CLOSED with cited tails |
| 有限预筛与删除均匀性 | Lemmas 3.1--3.2：first-failure 恒等式、支持缩减和正块算子证明；348510 个中等素数全部生成 | CLOSED |
| 合成筛系数不会放大余项 | §3.4：逐点乘法不等式，分离素数范围，系数绝对值 <=1；支持 <kappa Q_* D0 D | CLOSED |
| BJS 显式筛适用性 | §3.4：固定密度 1/(p-1)，所有中间端点的产品条件，D>=x 或 x^2，s 范围；Table 1 常数106、108 | CLOSED with cited sieve theorem |
| 大模数中诱导的小原始导数 | Lemma 4.1：先完整重排 d=mr，全部 r<=Q 成本进入低导数预算，不能仅删除西格尔零项 | CLOSED under (Z) |
| 高导数平均输入 | 精确保留 BJS Lemma 30 证明末式中的 (r,N)=1 和 log Y；应用时 r|d 自动满足；不借用整条引理的例外模数前提 | CLOSED with cited mean theorem |
| AP 从 psi 到 pi 和序列中心化 | Lemma 4.1：prime-power、初段、部分求和、删除 N 的素因子均收费；得到 c4(N)L^2<2.21 | CLOSED |
| 固定截断的矩形估计 | Lemma 5.1：低导数、诱导删除、dyadic 最后一段和几何和逐项证明；五常数39、108、26、88、106 | CLOSED under (Z) and cited large sieve |
| Chen 三项接入 | §6：引用 BJS Lemma 37/Nathanson Theorem 10.2；保留全部有限损失和额外1；不复用 Wu 十五项或旧四素数 switching | CLOSED with cited combinatorial inequality |
| 素因子行的模型和积分 | §6.2：固定 q 的筛、s_q、变量替换、Stieltjes 误差，以及 (q,d)->qd 的单射余项汇总 | CLOSED |
| 三素数行唯一性与主质量 | §6.3：唯一的 <y 素因子及排序因子恢复；短盒扩大、可能的负差平移、两个单调核与末盒费用 | CLOSED |
| 矩形模型删除费用量纲 | 式(14)是 L^4 e^{-L/8} B(L)，已修正原先少计的一个 L，并独立由原始费用重建 | CLOSED |
| 唯一数值总账 | 式(15)包含 M0,M1,M2,E_AP,E_R,E_del,E_f；没有事后拟合费用 | CLOSED |
| 全半直线 | §7.1：固定有理 alpha，逐项解析导数与正包络；不把离散素数/盒子当连续函数，不以端点采样替代证明 | CLOSED |
| 起点区间与印刷小数 | 主复现程序直接检查论文列出的全部粗区间、支持余量、质量门限和完整枚举计数 | CLOSED by finite certification |
| 独立复现路径 | 干净目录只含声明的12个源文件，开始时不含任何旧 JSON；从头运行完整素数生成和 Arb/ACB 积分 | CLOSED by recorded run |
| 外部零点定理的独立证明/Lean 构建 | 不属于本地数值认证；正文明确作为外部定理引用 | NOT CLAIMED |
| 所有方法/连续参数的全局最优 | 有限五分支搜索不能证明这个命题 | NOT CLAIMED |

## 数值与参数搜索

固定参数为 alpha=2158047/500000000、delta=1/1000000、T=5000000、
u=11/2、w0=5、kappa=11/169、epsilon=1/41219。
`qrh_chen_reproduction.json` 中起点完整净余量被严格包在
`(9.8468e-7,9.8469e-7)`。

穷举范围是整数 `50000<=log N<=51200`；每个分支找到首个正端点后停止，
未找到的分支检查完整范围。共5461个严格比较。

| delta | 该配方范围内首个正整数 log N |
|---|---:|
| 1e-4 | 51095 |
| 1e-5 | 51049 |
| 1e-6 | 51046 |
| 1e-7 | 51066 |
| 1e-8 | 范围内未通过 |

51045 的对应搜索配方净余量小于 -1.5399e-6；它不是对其他参数或方法的不可能性证明。
论文在51046之后固定 alpha，不沿用随 L 取整的搜索配方。

## 复现与排版

环境：Python 3.12.8、python-flint 0.9.0、192-bit Arb/ACB。
完整命令、所需文件和历史记录边界见 `qrh_REPRODUCE.md`。
正文已用现有 TeX Live 2025 编译为16页 PDF；无未解析引用或溢出警告，逐页视觉核对。
内置编辑器保留打开；其编译器报平台目录错误，因此使用本机已有 TeX，未安装新工具。
旧 c=24.9 稿和旧条件稿未被新结果覆盖，也没有修改 GitHub。
