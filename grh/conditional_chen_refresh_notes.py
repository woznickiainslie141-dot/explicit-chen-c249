"""Refresh the current numerical notes from the endpoint certificate.

This writes human-readable ledgers, preserving the assumption discussion and
the explicitly historical part of the original research memorandum.
"""
import json
import math
import hashlib
from fractions import Fraction
from pathlib import Path
from flint import arb, ctx

ROOT = Path(__file__).resolve().parent
ctx.prec = 192


def decimal_bound(value, places=13, upper=False):
    x = arb(value)
    scale = 10**places
    n = math.ceil(float(x)*scale) if upper else math.floor(float(x)*scale)
    if upper:
        while not x < arb(n)/scale:
            n += 1
    else:
        while not x > arb(n)/scale:
            n -= 1
    return f"{n//scale}.{n%scale:0{places}d}"


def compilation_status():
    """Report success only when the saved build and review match current bytes."""
    try:
        build = json.loads((ROOT/"conditional_chen_compile_certificate.json").read_text(encoding="utf-8"))
        review = json.loads((ROOT/"conditional_chen_pdf_review.json").read_text(encoding="utf-8"))
        source_hash = hashlib.sha256((ROOT/build["source"]).read_bytes()).hexdigest()
        pdf_hash = hashlib.sha256((ROOT/build["pdf"]).read_bytes()).hexdigest()
        if (build["status"] == "PASS" and source_hash == build["source_sha256"]
                == review["source_sha256"] and pdf_hash == build["pdf_sha256"]
                == review["pdf_sha256"] and review["all_pages_rendered"]):
            return (
                f"本机 TeX Live 已成功编译现行稿件为 {review['pages']} 页 PDF，"
                "路径 output/pdf/conditional_chen_grh_bridge.pdf。交叉引用已稳定，"
                "没有编译错误、排版溢出或重复链接；全部页面已渲染并作视觉检查。"
                "复现命令为 python conditional_chen_build_pdf.py；编译和页面审查记录分别为 "
                "conditional_chen_compile_certificate.json、conditional_chen_pdf_review.json。"
                "这些检查不认证数学证明。内置 LaTeX 沙箱仍报 "
                "Unable to find standard directories for platform，应用内自动预览尚未恢复。")
    except (OSError, ValueError, KeyError):
        pass
    return "现行稿件的编译与 PDF 页面审查尚未验证，或其记录已不匹配当前源文件。"


def weighted_distribution_note():
    report = json.loads((ROOT/"conditional_chen_weighted_distribution_certificate.json")
                        .read_text(encoding="utf-8"))
    assert report["status"].startswith("PASS")
    for filename, expected in report["input_sha256"].items():
        assert hashlib.sha256((ROOT/filename).read_bytes()).hexdigest() == expected
    rows = ["| 分布水平 | C≤1 | C≤10³ | C≤10⁴ | C≤10⁶ |",
            "|---|---:|---:|---:|---:|"]
    for row in report["endpoint_rows"]:
        cases = [case for case in row["endpoint_cases"] if case["C_upper_budget"]]
        assert [case["C_upper_budget"] for case in cases] == [1,1000,10000,1000000]
        rows.append("| "+row["theta"]+" | "+" | ".join(str(case["log_start"]) for case in cases)+" |")
    return (
        "此前固定四组权的下界只需实际筛权的单向有符号余项 R_N^-≥−C N/log⁴N，"
        "对筛法解析域中的所有偶数 N≥X 成立；解析域为 L>3log(10^12)、T<z、2<s<3。"
        "原完整全模数绝对误差输入以同一 C、X 推出它；"
        "这项特定量化输入本身不是已有的标准命名猜想，也未由 GRH 证明。\n\n"
        "先对全部素数的补数施加下界权，再删除整除 N 的素数，得到 "
        "S(A,z)≥pi(N)G_N^-+R_N^-−omega(N)。删除费用缩为 "
        "L³ exp(−L)/log 2；排除补数一另扣 L² exp(−L)。"
        "使用实际正主项分配余项预算，得到下表解析起点 A：\n\n"
        + "\n".join(rows) + "\n\n"
        "每格结论都是：若相应量化输入的常数不超过列示上界，"
        "则所有偶数 N≥max{X,exp(A)} 满足 pi₂(N)>10^-5 N/L²。"
        "列示数字是允许预算，不是已经证明的分布常数。X 在每项结论中保留。\n\n"
        "例如已提供 theta=99/100、C=10^6 的同一全模数量化输入时，"
        "原分配要求 L≥sqrt(10^7)≈3162.27766，新分配只需 L≥829；"
        "两者均须 N≥X。这改进了常数到阈值的归约，不提供普通 EH 下的数值纪录。\n\n"
        "完整全模数输入仍受整数性限制；这些限制不能直接移植到单个有符号和。"
        "即使新余项费用为零，现行固定参数表达式的首个正整数起点仍为 "
        "722、200、158、136，前一个整数处的主筛因子严格为负，不能据此拼接 B。\n\n"
        "实际权重保留 lambda_N^-(p)=−1 的大素数系数。沿 N=2^(3000k) 的四组水平，"
        "Bertrand 定理给出 p>N^(theta/3)，故无法满足平衡三段分解所需的零系数条件。"
        "这只排除了直接套用三重可因子化定理，没有判定所有二重可因子化变体。"
        "当前支持也不能直接缩为 MPZ 的平滑模数而忽略剩余项。\n\n"
        "新增证书以 192 位 Arb 认证端点和标量恒等式，并用精确有理数检查 "
        "10748 个复合权情形、99158 个逐点模式、21496 个有符号恒等式。"
        "有限大素数产品与矩输入复用并校验哈希，未在该脚本重算；"
        "没有认证分布上界、完整拼接或全局最优。文献适用条件详见 "
        "conditional_chen_weighted_distribution_20261004.md。")


def dense_grh_note():
    report = json.loads((ROOT/"conditional_chen_dense_band_grh_certificate.json")
                        .read_text(encoding="utf-8"))
    assert report["status"].startswith("PASS")
    for name, expected in report["input_sha256"].items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest() == expected, name
    low, strong = report["GRH_positivity"], report["GRH_at_31000"]
    rows = []
    for row in [low,strong]:
        rows.append("| "+str(row["log_N0"])+" | "+row["alpha"]+" | "+row["beta"]+" | >"
                    +decimal_bound(row["margin_ball"])+" |")
    pre = low["presieve"]
    return (
        f"前一组 GRH 阈值为 exp({low['log_N0']})，声明 pi₂(N)>10^-5 U_N N/L²。"
        "在 exp(31000) 之后，计数下界系数提高到 0.001。前一组 exp(30824) 和 "
        "exp(31000) 下系数 0.0006 的证书保留为对照。此次阈值比值为 exp(77)。\n\n"
        "定义权重、T=5000000、w₀=5、Q_*=15、u=49/8、M=12、κ=11/169 "
        "与普通指数集合均不变。扩大有限 majorant 集合为 "
        "Gamma_G_new={a/25:2≤a≤9}；完整整数筛和 192 位正项矩递推生成的新有限文件给出 "
        f"奇偶有限质量差分别 <{decimal_bound(pre['upper_relative_Rosser_gap_ball'],12,True)}、"
        f"<{decimal_bound(pre['lower_relative_Rosser_gap_ball'],12,True)}。"
        f"取 epsilon={low['epsilon']}、eta+={low['eta_plus']}、eta−={low['eta_minus']}。\n\n"
        "令 H=sqrt(N)/L^beta，并保留原 GRH AP 常数的正上包络 "
        "P(L)=q_+(L)+1/(16pi)+2 exp(-L/2)/log 2。"
        "实际方自由模数计数给出 R_A≤0.608 P(L) N/L^(beta−1)。"
        "Chen 不等式中此余项的总系数为 3/2，因此只把原 AP 首项改为 "
        "0.912 P(L)/(U_* L^(beta−3))；两个其余 AP 费用、矩形及有限费用全部保留。\n\n"
        "| log N₀ | alpha | beta | 认证最终裕量 |\n"
        "|---|---|---|---:|\n"+"\n".join(rows)+"\n\n"
        f"较低起点的支持裕量 >{decimal_bound(low['support_slack_ball'])}；"
        "矩形余项 <10^-248，有限损失 <10^-1659。较强计数界的相应界为 "
        "<10^-250、<10^-1673。固定 alpha、beta 后，支持裕量和 H 递增，"
        "q_+ 及全部减项递减，覆盖整个半直线。仍只假设 GRH，没有新增分布假设。\n\n"
        "dense_band_grh_certificate 记录输入哈希；当前端点复算复用新生成的 "
        "dense_band_grh_finite.json，原统一产品文件也按哈希复用。"
        "它复现原 L=30824 的五项标量至 10^-45；修改后的 AP、支持和最终裕量又以 "
        "直接公式独立核对至 10^-45。完整重算命令为 "
        "python conditional_chen_dense_band_grh_certificate.py；只复算端点时加 --reuse-finite。\n\n"
        "只扩大 majorant 集合时起点为 exp(30779)；再优化 AP 系数和 beta 得到上述结果。"
        "beta=3.50 到 3.57 的八个候选不是全局优化证明。当前仍未完成两端拼接，"
        "外部解析引理应用与稿件仍需独立数学审查。")


def prefix_block_note():
    sparse = json.loads((ROOT/'conditional_chen_sparse_modulus_certificate.json').read_text(encoding='utf-8'))
    refined = json.loads((ROOT/'conditional_chen_prefix_block_refinement.json').read_text(encoding='utf-8'))
    checks = json.loads((ROOT/'conditional_chen_prefix_block_endpoint_checks.json').read_text(encoding='utf-8'))
    exact = json.loads((ROOT/'conditional_chen_prefix_block_checks.json').read_text(encoding='utf-8'))
    for result in [sparse,refined,checks,exact]:
        assert result['status'].startswith('PASS')
        for name,expected in result['input_sha256'].items():
            assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==expected,name
    low,strong=sparse['GRH_positivity'],sparse['GRH_at_31000']
    rows=['| 分布水平 | C≤1 | C≤10³ | C≤10⁴ | C≤10⁶ |','|---|---:|---:|---:|---:|']
    for theta in ['3/4','9/10','19/20','99/100']:
        cells=[r for r in refined['distribution_endpoints'] if r['theta']==theta]
        rows.append('| '+theta+' | '+' | '.join(str(r['log_start']) for r in cells)+' |')
    parameters=['| 层级、C 预算 | T、w₀、u、h | epsilon、eta+、eta− | C1、C2 |',
                '|---|---|---|---|']
    for row in refined['distribution_endpoints']:
        case=row['case']; p=case['presieve']
        parameters.append(f"| {row['theta']}、{row['C_upper_budget']} | {p['T']}、{p['exact_cutoff']}、{p['u']}、{p['log_bin_width']} | "
                          f"{case['epsilon']}、{case['eta_plus']}、{case['eta_minus']} | {case['C1']}、{case['C2']} |")
    text=(
        f"前一组仅假设 GRH 的已认证阈值为 exp({low['log_N0']})，声明 pi₂(N)>10^-5 U_N N/L²。"
        "exp(31000) 后系数提高到 0.0179。此前 exp(30747)、0.001 的证书与证明保留为历史对照。"
        "这是未经独立研究审稿的稿件结论，不是已发表纪录。\n\n"
        "新有限递推保留所有先前已通过的 Rosser 检查。按 floor(log p/h) 合并素数块，"
        "对每一步使用必要的可能通过及可能失败条件；取整含混区同时计入两个正上界。"
        "正系数块多项式 D(a)=product(1+t_p+t_p a)、C(a)=(D(a)−1)/(1+a)，t_p=1/(p−2)，"
        "保证删素数占优。截去块内次数大于 K 的项不会移除真实路径。它只改变质量差的上界。\n\n"
        "该组 GRH 参数为 T=5000000、w₀=5、Q_*=15、u=87/16、h=1/200，"
        "K=22、J=16775、1810 个非空块，完整处理 348510 个中间素数。"
        f"epsilon={low['epsilon']}、eta+={low['eta_plus']}、eta−={low['eta_minus']}；κ=11/169。"
        f"支持成本 c≈{float(arb(low['effective_log_support_cost_ball'])):.12f}。\n\n"
        "新增支持数量引理给出 #非零模数≤2^m D H₀^sigma product(1+p^(−sigma))。"
        "对 sigma=4/5，完整 Euler 乘积给出相对于全支持区间的比例 "
        f"<{decimal_bound(low['sparse_modulus_count']['selected']['count_fraction_ball'],12,True)}<1/1653。"
        "该界适用于删素数及所有实际 A、A_q 筛权，A_q 的级别为 D/q；"
        "对 A_q 总和，把 q 吸收到大因子中得到 qd₁<D；q 是 qd 中唯一不小于 z 的素因子，"
        "故 (q,d)→qd 对整个求和单射，整体像集仍受同一数量界支配。"
        "双线性余项仍按最大模数 H 收费。令 beta=3、H=sqrt(N)/L³，"
        "AP 首项为 3P(L)/(2×1653 U_*)，两个其余 AP 项、矩形及有限费用全部保留。\n\n"
        "| log N₀ | alpha | beta | 支持裕量 | 认证最终裕量 |\n"
        "|---|---|---|---:|---:|\n"
        + '\n'.join('| '+str(r['log_N0'])+' | '+r['alpha']+' | '+r['beta']+' | >'
              +decimal_bound(r['support_slack_ball'])+' | >'+decimal_bound(r['margin_ball'])+' |'
              for r in [low,strong])+'\n\n'
        "低起点的矩形及有限费用分别 <10^-212、<10^-1446；31000 的对应界为 "
        "<10^-250、<10^-1673。固定 alpha、beta 后支持和 H 递增，所有扣项递减，"
        "结论覆盖全部半直线。逐项直接公式独立复现两组各八项量至 10^-43。\n\n"
        "只改进前缀递推和取整时起点为 exp(28279)，31000 后系数 0.0116；"
        "再用实际模数计数得到当前结果。先前粗分块给出 exp(28334)。"
        "一千万、两千万的所测 GRH 参数均较差；有限搜索没有证明全局最优。\n\n"
        "分布路线以同一 signedcriterion 对新的指定筛权族认证如下解析起点 A：\n\n"
        +'\n'.join(rows)+'\n\n'
        "每格严格保留结论 N≥max{X,exp(A)}，声明 pi₂(N)>10^-5 N/L²。"
        "若仅给有符号输入，它必须用于该格下表指定的权族；完整全模数输入以相同 C、X "
        "蕴含全部所选有符号输入。C 是允许预算，X 未指定；普通 EH 不给这些数值。\n\n"
        +'\n'.join(parameters)+'\n\n'
        "新的五组分布权中最低支持成本约 31.052421204，仍超过 B=4×10^18 所允许的 "
        "log B/3≈14.277608679，因此即使有符号误差为零也不能直接拼接。"
        "全绝对误差输入的整数性约束继续有效；它不能直接用于单个有符号和。\n\n"
        f"精确对照覆盖 {exact['deleted_prime_subsets']} 个删素数子集、"
        f"{exact['exact_first_failure_identities']} 个首次失败恒等式、"
        f"{exact['pointwise_divisor_patterns']} 个逐点模式、{exact['exact_block_coefficient_groups']} 个系数组。"
        "支持计数另有 18 个精确小例；16 个分布格的主项、裕量及解析域通过直接公式复核。"
        "有限输入、缓存复用和完整产品重算均记录 SHA-256；这些检查不证明一般引理、"
        "外部解析输入、独立研究审查、全局最优或完整拼接。\n\n"
        "复现顺序：prefix_block_checks.py、prefix_block_endpoint_certificate.py、"
        "prefix_block_refinement.py、sparse_modulus_certificate.py、prefix_block_endpoint_checks.py，"
        "各文件名均以 conditional_chen_ 开头。一般性证明在稿件的 Positive prefix blocks "
        "及 Counting the actual modulus support 引理中。")
    (ROOT/'conditional_chen_prefix_block_20261004.md').write_text(text+'\n',encoding='utf-8',newline='\n')
    return text,sparse


def rough_interval_note():
    report=json.loads((ROOT/'conditional_chen_rough_interval_certificate.json').read_text(encoding='utf-8'))
    checks=json.loads((ROOT/'conditional_chen_rough_interval_checks.json').read_text(encoding='utf-8'))
    small=json.loads((ROOT/'conditional_chen_binned_support_checks.json').read_text(encoding='utf-8'))
    for result in [report,checks,small]:
        assert result['status'].startswith('PASS')
        for name,h in result['input_sha256'].items():
            assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==h,name
    r=report['GRH_positivity'];count=r['binned_support_count'];start,anchor=r['valid_log_interval']
    assert checks['analytic_join']['covered_log_interval']==[start,anchor]
    text=(
        f"前一组仅假设 GRH 的阈值为 exp({start})，声明所有偶数 N≥exp({start}) 满足 "
        "pi₂(N)>10^-5 U_N N/log²N；exp(31000) 后系数仍为 0.0179。"
        f"此前 exp({anchor}) 的半直线证书保留，并用作这一轮解析区间的接点。"
        "新结论有一般计数引理和严格数值证书，仍未经独立研究审稿，不是已发表纪录。\n\n"
        "这一轮没有增加 AP 分布假设。第一项改进是正项对数直方图："
        "b(p)=floor(log p/h)，J=ceil(log H₀/h)，K 为小素数最小乘积所决定的真实次数上限。"
        "每块的 product(1+p^(−sigma)a) 只保留次数≤K，再按标签相乘并截去标签≥J。"
        "真实 d<H₀ 的标签仍在其中；以 exp(sigma min(log H₀,h(B+K))) 放大，"
        "它贡献至少 1。多出的正项不影响上界，删除素数只会降低块系数。\n\n"
        f"现行计数取 sigma={count['sigma']}、h={count['log_bin_width']}、K={count['maximum_actual_degree']}、"
        f"J={count['label_limit']}，{count['nonempty_blocks']} 个非空块；"
        f"完整处理 {count['medium_prime_count']} 个中等素数。A₀<5.488142×10³⁰，"
        f"nu₀=2^omega(Q_*) A₀/(Q_*H₀)<{r['binned_modulus_count_budget']}。"
        "质量差仍用 h=1/200 的先前正项分块结果；两个网格负责不同的上界。\n\n"
        "第二项改进是大因子 d₁ 不含任何 p≤T 的素因子。对整数序列使用密度 1/p，"
        "在 2、3、5 完整排除，并对中等素数用同一有限 Rosser 上界权。"
        "其归一化正项递推 t=1/(p−1)≤1/(p−2)，故质量差受现有 eta+ 支配。"
        "得到 R_T(D)≤a_T D+b_T，其中 a_T=(1+eta+) product_(p≤T)(1−1/p)，"
        "b_T=2Q_*H₀；因子 2 来自整数筛所需的偶素数，不能遗漏。"
        "每个除数的 floor 余项绝对值≤1，有限边界费用严格保留。\n\n"
        f"另一份完整整数素数筛重新生成 {checks['complete_integer_primes_independently_generated']} 个素数，"
        f"确认 V_T<0.036397999733 和 a_T<{r['rough_density_budget']}。"
        "实际模数数量因此≤nu₀(a_T+b_T/D)H。对 A_q 总和，把 q 吸收到大因子中，"
        "qd₁<D 仍为 T-rough；q 是 qd 中唯一不小于 z 的素因子，故整个 (q,d)→qd 映射单射。"
        "该总和的整体像集受同一粗整数数量界支配，不逐个 q 重复支付 floor 费用。"
        "双线性费用继续按最大模数 H 收费。\n\n"
        f"新解析区间为 {start}≤L=log N≤{anchor}，固定 alpha={r['alpha']}、beta={r['beta']}，"
        "D=exp((1/2−alpha)L)、H=sqrt(N)/L^beta。beta<3 时 AP 首项可能随 L 增大，"
        "不能直接声明同一参数覆盖半直线。以 "
        "3 nu*/(2U*) (a*+b_T exp(−tL₀)) P(L₀) L*^(3−beta) 作常数费用上界，"
        "明确覆盖整个闭区间。支持裕量和 H 递增，其他正扣项递减；"
        f"在 L*={anchor} 接入旧 GRH 半直线，两个范围用同一计数系数和 U_N 归一化，"
        "接点包含在两段中。因此结论对所有 L≥L₀ 成立，而不是只检查五个点。\n\n"
        f"起点 s≈{float(arb(r['s_ball'])):.6f}；支持裕量>{decimal_bound(r['support_slack_ball'])}；"
        f"AP 常数费用<{decimal_bound(r['AP_budget_ball'],15,True)}；"
        f"最终裕量>{decimal_bound(r['margin_ball'],15)}。"
        "矩形费用<10^-194，有限费用<10^-1337，整数筛边界比例<10^-5302，均保留正收费。\n\n"
        "独立公式复现 12 项端点量至相对误差 10^-43，并重新检查旧尾部证书。"
        "48 个有限直方图情形枚举 9216 个子集，核对 17172 个有理数系数；"
        "整数粗筛另检查 16384 个除数 floor 项与 24316 个整数模式。"
        "跨 q 聚合另有 8608 个模数模式的精确单射与整体计数检查。"
        "135 个端点候选是有限参数搜索，没有证明全局最优。\n\n"
        "复现顺序：python conditional_chen_binned_support_checks.py，"
        "python conditional_chen_rough_interval_certificate.py，"
        "python conditional_chen_rough_interval_checks.py。"
        "支持直方图按生成器 SHA-256 校验缓存，整数产品本轮完整重算；"
        "独立标量复核没有调用生产端点评估器。\n\n"
        "这里接上的是两个解析估计的有效范围，尚未拼接至有限验证 B=4×10^18。"
        "现行直接分布筛的最小支持成本仍为 31.052421204…，而 B 只容许小于 14.278；"
        "需要进一步改变筛权或解析接口。\n\n"+compilation_status())
    (ROOT/'conditional_chen_rough_interval_20261004.md').write_text(text+'\n',encoding='utf-8',newline='\n')
    return text,report


def accepted_interval_note():
    report=json.loads((ROOT/'conditional_chen_accepted_interval_certificate.json').read_text(encoding='utf-8'))
    checks=json.loads((ROOT/'conditional_chen_accepted_interval_checks.json').read_text(encoding='utf-8'))
    small=json.loads((ROOT/'conditional_chen_accepted_support_checks.json').read_text(encoding='utf-8'))
    for result in [report,checks,small]:
        assert result['status'].startswith('PASS')
        for name,h in result['input_sha256'].items():
            assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==h,name
    weak,strong=report['GRH_positivity'],report['GRH_at_31000']
    start,anchor=weak['valid_log_interval']
    assert checks['analytic_join']['covered_log_interval']==[start,anchor]
    assert checks['strong_half_line']['claim']==strong['claimed_margin']
    coefficient=str(float(Fraction(strong['claimed_margin'])))
    text=(f"现行仅假设 GRH 的阈值为 exp({start})，所有偶数 N≥exp({start}) 满足 "
        f"pi₂(N)>10^-5 U_N N/log²N；N≥exp(31000) 时系数提高到 {coefficient}。"
        f"较低起点比前一组 exp({anchor}) 降低 {anchor-start}，新增分布假设为零。"
        "这是本稿未经独立研究审稿的结果，不是已发表纪录。\n\n"
        "新计数直接保留有限 Rosser 支持的逐步接受条件。上权在奇数长度检查，"
        "下权在偶数长度检查；每个子集记录仅上权、仅下权、两权均存活和长度奇偶，"
        "形成六个互斥状态。完整对数块中 b=floor(log p/h)，"
        "c_b=ceil((log D₀−3log p_min)/h)。真实检查通过蕴含前缀标签 B<c_b；"
        "低标签保留所检查的符号，高标签只删除它，未检查的符号保留。"
        "初始为空乘积、两权均存活、偶数长度；两符号都死亡时永久删除。\n\n"
        "每个相同标签块选择 j 个素数的次数恰为 binom(m_b,j)，"
        "因此用 sum_(j≤K) binom(m_b,j) A_b^j 更新。每次只截去标签≥J；"
        "真实接受乘积 d<H₀ 的次数≤K、最终及所有前缀标签<J，均不会被遗漏。"
        "每个保留子集只有一个存活状态，两权均接受的乘积不会重复计数。"
        "删除素数后保留完整 K、J 和 c_b，仅降低块素数个数 m_b；"
        "所有算子为正，故逐状态、逐系数受完整集合支配。"
        "一般证明见稿件六状态并集引理。\n\n"
        "弱端点 T=5000000、w₀=5、u=87/16、κ=11/169。支持计数用 h=1/100，"
        "质量差仍用 h=1/200，两项分别认证。完整 348510 个中等素数形成 992 块，"
        "K=22、J=8115；A_acc=655257348020585286337311701158，"
        "是前一未过滤正项上界的约 1/8.375。"
        f"nu_acc<1.007455×10^-6<{weak['binned_modulus_count_budget']}。"
        "用既有整数粗筛 R_T(D)≤a_T D+b_T 支配大因子。"
        "对整个 A_q 总和，将 q 吸收到 T-rough 大因子，利用 q 为唯一≥z 的素因子，"
        "保持 (q,d)→qd 单射；该整体像集也受新并集数支配。"
        "floor 边界、prime-modulus centering 和双线性最大模数费用均保留。\n\n"
        f"弱参数 alpha={weak['alpha']}、beta={weak['beta']}，新闭区间 "
        f"{start}≤L=log N≤{anchor} 使用常数 AP 上界。"
        "其后接上 [24799,26820] 和 [26820,∞)，都采用同一 10^-5 声明和 U_N 归一化。"
        "支持裕量递增、其他扣项递减、接点包含在两段中，故覆盖所有 L≥L₀。"
        f"起点支持裕量>{decimal_bound(weak['support_slack_ball'])}，"
        f"AP 费用<{decimal_bound(weak['AP_budget_ball'],15,True)}，"
        f"最终裕量>{decimal_bound(weak['margin_ball'],15)}。"
        "矩形费用<10^-188，有限费用<10^-1300，正的粗筛边界比例<10^-5189。\n\n"
        f"较强端点另用 u=11/2、alpha={strong['alpha']}、beta=3，"
        f"nu_acc<{strong['binned_modulus_count_budget']}，eta+=1/3231、eta−=1/3242。"
        "L^(3−beta)=1，P(L)、粗筛边界和所有其余费用在整条半直线上递减，"
        f"L=31000 的严格裕量>{decimal_bound(strong['margin_ball'],15)}，"
        f"因此可声明系数 {coefficient} 对所有 L≥31000 成立。"
        "较大系数由这一新半直线证明，不经过旧 0.0179 结论作数值提升。\n\n"
        "独立标量程序对两个端点各复现 12 项量至相对误差 10^-43，"
        "重新生成全部 348513 个整数素数，并检查六状态整数计数之和。"
        "精确有限检查枚举 4608 个完整子集、114876 个整数直方图系数、"
        "104976 个删素数子集模式和 384 个删素数生产直方图；"
        "既有跨 q 聚合另检查 8608 个模数模式。"
        f"{len(report['candidate_rows'])} 个弱端点候选和 9 个较强端点候选属于有限参数搜索，"
        "没有证明全局最优。样本点的复算只检查实现；连续范围依据一般引理和单调性证明。\n\n"
        "复现顺序：python conditional_chen_accepted_support_checks.py，"
        "python conditional_chen_accepted_interval_certificate.py，"
        "python conditional_chen_accepted_interval_checks.py。计数缓存核对生成器哈希，"
        "独立端点检查没有调用生产端点评估器。\n\n"
        "尚未拼接至有限验证 B=4×10^18，也没有完成独立研究审稿或证明现有理论全局边界。"
        "既有两阶段分布筛的支持成本障碍仍在，下一步仍需优化筛权及验证端。\n\n"
        +compilation_status())
    (ROOT/'conditional_chen_accepted_interval_20261004.md').write_text(text+'\n',encoding='utf-8',newline='\n')
    return text,report


def full_rosser_local_note():
    checks=json.loads((ROOT/'conditional_chen_full_rosser_bridge_checks.json').read_text(encoding='utf-8'))
    assert checks['status'].startswith('PASS')
    for name,h in checks['input_sha256'].items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==h,name
    best=checks['local_scalar_reproductions'][-1]
    assert best=={'w':2,'u':'61/20','C_budget':3489,'kappa':'3/25'}
    for row in checks['centered_pi_input_counterexamples']:
        assert row['status'].startswith('REFUTED')
        assert arb(row['necessary_A_epsilon_1_1000_lower_ball'])>arb(row['candidate_A_capacity_ball'])
    return (
        "新增全素数有限 Rosser 构造 T=2000000、w=2、u=61/20，支持低于 N^(99/100)，"
        "下权主质量下界约 0.0444000801948。若实际带符号余项 R_W(N)≥−3489N/log⁴N "
        "对每个偶数 4×10^18≤N≤8×10^18 成立，则整段有 pi₂(N)>10^-5 N/log²N。"
        "这一余项条件未证明，不能声称仅 GRH、标准 EH 或标准 Montgomery 已给出此拼接。"
        "连续范围由固定权与单调性证明，尚未覆盖 8×10^18 以上。\n\n"
        "保留全部接受检查的加权支持计数已完成。三组 epsilon=1/1000 的逐模数中心化 pi "
        "输入所容许 A 约 0.255103、0.255829、0.196596，实际筛权中的三个反例分别强制 "
        "A>0.643081、0.781228、0.559409，因此这些具体数值输入已被否定。"
        "反例使用实际 a=N，素数都有纯整数验证的 Lucas 证书；没有反驳标准渐近猜想或 "
        "整个带符号余项条件。\n\n"
        "新独立检查重算 148932 个奇素数，检查 9216 个完整子集、152316 个加权系数、"
        "768 个删素数生产直方图及 26244 个点态模式。一般证明和三个精确反例详见 "
        "[局部接口研究记录](conditional_chen_full_rosser_bridge_20261004.md)，"
        "数值及素性证书见 conditional_chen_full_rosser_bridge_checks.json。"
        "旧两阶段支持成本障碍只约束那些固定权；新的全素数构造绕过了该成本，"
        "带符号误差控制和后续连续覆盖仍是未完成项。"
    )


def main():
    c = json.loads((ROOT/"conditional_chen_sharp_endpoint_certificate.json").read_text())
    base = json.loads((ROOT/"conditional_chen_composite_certificate.json").read_text())
    weighted_text = weighted_distribution_note()
    dense_text = dense_grh_note()
    prefix_text,_ = prefix_block_note()
    rough_text,_ = rough_interval_note()
    current_text,current = accepted_interval_note()
    current_text+='\n\n## 全素数局部接口与有限数值输入反例\n\n'+full_rosser_local_note()
    current_coefficient=str(float(Fraction(current['GRH_at_31000']['claimed_margin'])))
    current_text+='\n\n## 前一组未过滤支持计数与粗整数界\n\n'+rough_text
    current_text+='\n\n## 更早 GRH 支持数量界与现行分布归约\n\n'+prefix_text
    dense = json.loads((ROOT/"conditional_chen_dense_band_grh_certificate.json")
                       .read_text(encoding="utf-8"))
    g = c["GRH"]
    positive = c["GRH_positivity"]
    levels = [("EH_3_4", "3/4", 20), ("EH_9_10", "9/10", 10),
              ("EH_19_20", "19/20", 10), ("EH_99_100", "99/100", 10)]
    old = (ROOT/"conditional_chen_results_20261003.md").read_text(encoding="utf-8")
    assumptions = old.split("## 原文与新增假设的关系\n", 1)[1].split("\n## ", 1)[0].strip()
    table = ["| 分布层级 | T、u | epsilon、eta+、eta− | 条件性阈值 |",
             "|---|---|---|---|"]
    fees = []
    for key, theta, budget in levels:
        e, p = c[key], c[key]["presieve"]
        table.append(f"| {theta} | {p['T']}、{p['u']} | {e['epsilon']}、{e['eta_plus']}、{e['eta_minus']} | "
                     f"max{{X,exp(max{{{e['log_start']},sqrt({budget}C)}})}} |")
        fees.append(f"- 层级 {theta}：支持成本约 {float(arb(e['effective_log_support_cost_ball'])):.12f}；"
                    f"奇偶质量差分别 <{decimal_bound(p['upper_relative_Rosser_gap_ball'],12,True)}、"
                    f"<{decimal_bound(p['lower_relative_Rosser_gap_ball'],12,True)}；"
                    f"最终裕量 >{decimal_bound(e['surplus_ball'])}。")
    eh_table = "\n".join(table)
    eh_fees = "\n".join(fees)
    tradeoff = json.loads((ROOT/"conditional_chen_distribution_tradeoff_certificate.json").read_text(encoding="utf-8"))
    assert tradeoff["status"].startswith("PASS")
    tradeoff_rows = "\n".join(
        f"| {r['theta']} | >{r['displayed_C_lower']} | >{r['C_equals_one_onset_log_witness']} | "
        f"{r['budget']} | >{r['current_budget_endpoint_log_witness']} |"
        for r in tradeoff["rows"])
    tradeoff_text = (
        "全模数绝对误差和有整数性下界。令 y=log x，若 exp(theta*y)≥72，输入必须满足\n\n"
        "C>I_theta(y)=y³/(2theta)[exp(−(1−theta)y)−exp(−y)]−y⁴ exp(−y/2)，对所有 y≥log X。\n\n"
        "因此 C 与 X 不能独立指定为小数值。192 位 Arb 认证以下必要限制：\n\n"
        "| 层级 | X≤B 时的 C 下界 | C=1 时的 log X 下界 | 原固定预算 b | 该预算下 R 的下界 |\n"
        "|---|---:|---:|---:|---:|\n" + tradeoff_rows + "\n\n"
        "这里 B=4×10^18，R=max{log X,A,sqrt(bC)}。若 R≤M，输入在 x=exp(M) 强制 "
        "C>I_theta(M)>M²/b≥R²/b≥C，矛盾。所列是严格见证下界，未声称最优，"
        "通过下界也不意味着分布估计成立。\n\n"
        "特别是在原 sqrt(10C) 固定预算下，99/100 行的实际阈值必定大于 exp(834)，"
        "解析起点 141、142 并不能成为该预算下的实际起点。"
        "若取 C=1，则 X>exp(2246)。其余三行的固定解析起点仍大于所列 R 下界。\n\n"
        "进一步，在现行直接线性筛表达式及 U_*=1.32 exp(−gamma) 归一化下，"
        "主系数至多 K_theta=7.92 log(3theta−1)/(3theta)，所以正端点要求 C<K_theta R²。"
        "在 theta=99/100，I_theta(492)>K_theta*492²；故即使放宽误差预算和降低支持成本，"
        "此表达式仍需 R>492。这不约束所有 Chen 方法或其他余项模型。\n\n"
        "新增 distribution_tradeoff_certificate.py 除认证严格数值见证外，还直接用精确有理数"
        "核算 19701 个整数 (x,H) 的素数模数余类计数。有限对照不替代一般证明；"
        "证书没有提供任何分布上界或可行 C、X。")
    positivity_rows = "\n".join(
        f"| {theta} | max{{X,exp(max{{{c['distribution_positivity'][key]['log_start']},sqrt({budget}C)}})}} | "
        f">{decimal_bound(c['distribution_positivity'][key]['surplus_ball'])} |"
        for key, theta, budget in levels)
    smaller_grh = (
        "前一组较弱计数下界的 GRH 起点为 exp(30824)，"
        "声明 pi₂(N)>10^-5 U_N N/L²。此行保持相同有限权和 gap 输入，"
        "另取 alpha=4238/10^6、s=3.966096；"
        f"认证余量 >{decimal_bound(positive['margin_ball'])}，支持裕量 >0.002084。"
        "所有支持、参数、积分和费用另行检查，矩形余项 <10^-248、有限损失 <10^-1660；"
        "固定参数的单调性覆盖 L≥30824。exp(31000) 的较强下界继续保留。")
    lq_positive = arb(positive["log_N0"])/arb(10).log()
    digits_positive = math.floor(float(lq_positive))+1
    assert digits_positive-1 < lq_positive < digits_positive
    smaller_grh += (
        f" 此较低阈值的 log log N₀≈{float(arb(positive['log_log_N0_ball'])):.12f}，"
        f"约 {digits_positive} 位十进制数字。")
    terms = "\n".join(
        f"| {label} | {float(arb(g[key])):.12f} |"
        for label, key in [("M₀", "lower_ball"), ("M₁", "half_upper_ball"),
                           ("M₂", "half_switch_ball"), ("AP 费用", "AP_budget_ball")])
    lq = arb(g["log_N0"])/arb(10).log()
    digits = math.floor(float(lq))+1
    assert digits-1 < lq < digits
    grh_parameters = (
        f"w₀=5、Q_*=15、T=5,000,000、D₀=T^(49/8)、M=12、alpha=4215/10^6、"
        f"epsilon={g['epsilon']}、eta+={g['eta_plus']}、eta−={g['eta_minus']}、"
        "(C1,C2)=(106,108)、delta=10^-8。普通指数集合 "
        "Sigma_G={0,1/10,7/50,9/50,11/50,7/25,9/25}；带约束平方上界使用 "
        "Gamma_G={2/25,3/25,4/25,1/5,6/25}。")
    grh_bounds = (
        "全部 348510 个中间素数的有限递推给出奇偶界 "
        f"<{decimal_bound(g['presieve']['upper_relative_Rosser_gap_ball'],12,True)}、"
        f"<{decimal_bound(g['presieve']['lower_relative_Rosser_gap_ball'],12,True)}，"
        "总界 <0.000720336096。统一产品误差 <2.426054×10^-5。"
        f"实际支持成本约 {float(arb(g['effective_log_support_cost_ball'])):.12f}，"
        f"端点支持裕量 >{decimal_bound(g['support_slack_ball'],6)}；"
        f"认证最终裕量 >{decimal_bound(g['margin_ball'])}，声明保守下界 0.0006。")
    evidence = """- composite_certificate：本轮完整整数筛与 192 位有限输入重算通过；粗支持基线为 exp(31650)，分布起点 823、225、176、150。
- sharp_endpoint_certificate：复用上述有限输入，记录输入 SHA-256；重新扫描完整中间素数集求支持最大值，认证所有新端点费用。八项 GRH 标量与粗支持基线复算一致至 10^-45；没有声称本脚本重算 10^8 产品和所有矩。
- uniform_product_certificate/checks：完整生成 5761454 个奇素数、检查 5760226 个跳跃；独立检查 4096 个整数序列、1954 个精确素数产品端点对及 9 个尾产品情形。
- prefix_band_certificate/checks：实际生产使用完整素数集；独立精确代数对照使用非素数的二幂形式原子，共 12 组配置、2688 个首次失败/删除子集检查、25 个产品带实例，包括一个两原子首次失败。它们不是对实际素数分布的额外假设。
- prefix_band_eh_search：96 个严格标量候选；所选有限输入由主证书重新生成，支持缩减后的起点由端点证书重新认证。未证明全局最优。
- sharp_support_certificate：完整相邻素数扫描给出 κ=11/169、5/49；52488 个精确系数/删除子集情形，42226 个接受系数通过支持检查。
- rosser_checks：602368 个逐点模式、16 个系数情形、3088 个质量恒等式/删除子集检查，覆盖 5、9/2、35/8、25/4、17/4、49/8、39/8、33/8 八组有理水平。
- degreewise_rankin_checks、quadratic_rankin_checks：各保留九组配置及 2304 个精确子集对照，作为历史方法的检查。
- weight_checks：此前 2347008 个有限逐点模式检查保留。
- distribution_tradeoff_certificate：192 位 Arb 认证 C、X 及实际端点的必要限制，独立直接计数核算 19701 个整数 (x,H)；未提供分布上界。
- weighted_distribution_certificate：保留实际下界权的单向余项，重新分配已给常数 C 的误差预算；10748 个精确复合权情形和 21496 个有符号恒等式通过。校验复用有限输入哈希，未提供 C、X。
- source_audit：检查源结构、引用、现行数字和证书依赖哈希；不能证明外部解析定理或验证 PDF。"""
    obstruction = (
        "前一组固定分布权在 B=4×10^18 要求 c≤(theta−2/3)log B。层级 3/4、9/10 "
        "允许成本仅为 3.569402170、9.994326075。即使 theta 趋近 1，允许成本最多 "
        "log B/3≈14.277608679，仍小于这一组最高层级成本 37.960176484。"
        "因此即使把 AP 误差压到零，这些固定筛权仍不能拼接。这个障碍仅针对现行构造，"
        "不是所有条件性 Chen 方法的不可能性结论。"
        "本轮全素数有限 Rosser 权已通过在 B 的支持检查；这些旧成本不适用于该新构造，"
        "但新构造的带符号余项仍未证明，8×10^18 以上也尚未接到解析范围。")
    limits = (
        "稿件仍需独立研究审稿，外部筛定理和 AP/双线性输入的应用是明确依赖。"
        "未知的分布常数 C、X 没有被赋予便利数值。旧 exp(8500)、exp(15000) "
        "固定上端接口未闭合；旧 QEH_B 与 theta=3/4 有限预算被整数性下界反驳并撤回。"
        "这些问题没有反驳已发表的 GRH Chen 定理。\n\n"
        + compilation_status()
        + "完整拼接和整体独立审查尚未完成，研究目标继续进行。")
    results = f"""# 条件性显式 Chen：当前稿件、假设与审查状态

更新：2026 年 10 月 4 日。英文源稿为 conditional_chen_grh_bridge.tex。

当前稿件仅假设 GRH 的已认证阈值为 exp({current['GRH_positivity']['log_N0']})，声明 pi₂(N)>10^-5 U_N N/log²N。
exp(31000) 之后的较强下界为 pi₂(N)>{current_coefficient} U_N N/log²N。
这是未经独立研究审稿的稿件结果，不能称为已发表纪录。带未知 C、X 的四条分布归约
原固定误差预算下，较弱计数下界的解析起点为 774、210、165、141；
较强计数下界的起点为 786、212、166、142。完整拼接仍未完成。

## 原文与新增假设的关系

{assumptions}

## 当前 GRH 改进与认证

{current_text}

## 前一组有限网格改进（保留对照）

{dense_text}

## 前一组 GRH 参数与费用（保留对照）

{grh_parameters}

{grh_bounds}

实际支持为 κ Q_*D₀D，其中 κ=11/169。对 L=log N≥31000，
κ Q_*D₀D≤H=sqrt(N)/L^(7/2)。水平 D₀ 及有限质量差保持原定义，
支持节省 log(169/11)≈2.732003442，没有改变筛权或增加猜想。

较强计数界起点 exp(31000) 的 log log N₀≈{float(arb(g['log_log_N0_ball'])):.12f}，
阈值约 {digits} 位十进制数字。
F(s)≈{float(arb(g['F_s_ball'])):.13f}，J_h≈{float(arb(g['h_integral_ball'])):.13f}，
s=3.966280。积分采用 Acb；原加权引理的端点费用保留。

| 费用项 | 认证数值，展示时四舍五入 |
|---|---:|
{terms}
| 矩形余项 | <10^-248 |
| 有限损失 | <10^-1670 |

固定参数下支持裕量递增，所有减项在 L≥31000 递减，故覆盖整个半直线。

{smaller_grh}

## 首次失败的带约束平方上界

首次失败前缀 d 的最小素数为 r。D₀=T^u、u>3 时，短前缀的直接乘积界
与较长前缀中先前已通过的检查共同强制
d<D₀，而失败本身要求 dr²≥D₀。因此 δ=log(dr²/D₀) 满足 0≤δ<2log r。
在这个区间，exp(σδ)[1+cδ+bδ(2log r−δ)]²≥1（σ,c,b≥0），
且整个多项式平方在实轴上非负。先用完整素数集的非负期望占优，
再优化 Gram 二次式，所以删除素数后仍有效；不假设优化商本身单调。

证书对原始对数矩采用正项指数生成函数递推，再以 Arb 包含有符号的平移和 Gram 组合。
只接受符号严格认证的优化器。原生多项式截至四次计算所有所需矩，未丢掉子集质量；
浮点截零不参与证明。普通分次数上界、旧二次上界保留为可选基线。

支持缩减另有独立引理：末前缀若被检查，则 d<D₀/r²；若末素数未被检查，
前一个已通过前缀的最小素数 s>r 给出 d<D₀ r/s²。完整中间素数集上的
κ=max{{1/p_min²,max_(r<s) r/s²}} 统一控制所有删除子集。
固定 s 时只需检查它的前一个完整素数。w₀=5 时 κ=11/169，w₀=3 时 κ=5/49。

## 带 C、X 的分布归约

输入为 E_theta(x)=sum_(d≤x^theta) μ²(d) max_(a,d)=1 |π(x;d,a)−π(x)/φ(d)|
≤C x/log⁴x，对所有 x≥X。此路线不使用 RH、GRH 或 GEH，也不指定 C、X。

{eh_table}

四行均用 w₀=3、Q_*=3、κ=5/49、M=10、(C1,C2)=(113,114)，
Sigma_E={{0}}∪{{a/20:2≤a≤12}}、Gamma_E={{0,3/25,6/25,9/25}}。
实际费用 c=log(κ Q_*)+u log T；有限 gap 仍按 D₀=T^u 计算。

{eh_fees}

第一行声明 pi₂>0.01 N/L²，后三行声明 pi₂>0.02 N/L²。
无条件产品界、删除整除 N 的素数及排除整数 1 的费用均保留。
筛因子在 2<s<3 递增，负费用递减，故覆盖所列阈值之后的全部偶数。
标准 EH 不提供这里的数值 C、X，不能将这些起点当作普通 EH 下的有效数值纪录。

## 分布常数与适用起点的必要约束

{tradeoff_text}

如果仅要求 pi₂(N)>10^-5 N/L²，保持同一有限权、同一误差预算和同一 C、X 输入，
可另行认证以下更低起点：

| 分布层级 | 较弱计数下界对应的条件性阈值 | 认证最终裕量 |
|---|---|---:|
{positivity_rows}

## 可复现性与审查边界

按顺序运行 python conditional_chen_composite_certificate.py 和
python conditional_chen_sharp_endpoint_certificate.py；前者完整生成有限输入，
后者复用带哈希的输入重新核算缩减支持后的全部端点费用。
conditional_chen_refresh_notes.py 从现行证书刷新本说明及审查总账。

{evidence}

## 实际筛权余项与常数预算改进

{weighted_text}

## 完整拼接的剩余障碍

{obstruction}

{limits}
"""
    audit = f"""# 条件性 Chen：当前证明审查总账

2026-10-04。现行稿件 conditional_chen_grh_bridge.tex；现行 GRH 端点证书
conditional_chen_accepted_interval_certificate.json，分布端点为 prefix_block_refinement.json。
前一组 sharp_endpoint_certificate.json 和 composite_certificate.json 保留为对照。

| 内容 | 当前状态 |
|---|---|
| 原发表版 exp(exp(14)) | 已发表，仅假设 GRH |
| 本稿 exp({current['GRH_positivity']['log_N0']})；exp(31000) 下系数 {current_coefficient} | 仅假设 GRH，未经独立研究审稿 |
| 较弱分布计数界起点 774、210、165、141 | 参数化归约，C、X 未知，不要求 RH/GRH |
| 较强分布计数界起点 786、212、166、142 | 同一有限输入，较强声明 |
| 完整拼接至 4×10^18 | 未完成 |

GRH 包括这里的 RH；EH/GEH 是另一类常用分布猜想，GRH 不已知推出它们。
本稿分布输入的未知常数版本可由标准 EH 转换得到；指定数值 C、X 是额外定量输入。
四条分布路线没有同时要求 GRH，没有使用 GEH。

## 当前 GRH 改进与认证

{current_text}

## 前一组有限网格改进（保留对照）

{dense_text}

## 前一组 GRH 参数和费用（保留对照）

{grh_parameters}

{grh_bounds}

s=3.966280；F(s)≈{float(arb(g['F_s_ball'])):.13f}；
J_h≈{float(arb(g['h_integral_ball'])):.13f}。κ=11/169 的支持费用已实际计入。
最终端点费用由 192 位 Arb/Acb 认证，矩形余项 <10^-248、有限损失 <10^-1670，
均以正数收费。固定参数下覆盖 L≥31000。

{smaller_grh}

AP 估计排除偶数模数并取 r(1)=0；Aq 的误差为 r(qd)−r(q)/φ(d)，
(q,d)→qd 单射，费用 R_A+0.55 L R_prime 与替换主质量均保留。
矩形互素删除、双线性常数 m<3、末矩形及所有端点计入。
三素数质量采用计数上包络；归一化 U_N=2 exp(−gamma) C_N。

## 分布参数和费用

{eh_table}

四行用 w₀=3、Q_*=3、κ=5/49、M=10、(C1,C2)=(113,114)。
Sigma_E={{0}}∪{{a/20:2≤a≤12}}，Gamma_E={{0,3/25,6/25,9/25}}。
先前 96 个严格标量候选决定有限权参数；当前起点使用缩减支持另行认证。

{eh_fees}

较弱计数下界 pi₂>10^-5 N/L² 的四行单独认证如下，预算和 C、X 输入保持不变：

| 分布层级 | 较弱计数界阈值 | 最终裕量 |
|---|---|---:|
{positivity_rows}

## 解析与数值依赖

{tradeoff_text}

首次失败恒等式、删素数占优、次数限制、带约束平方 majorant、
Gram 组合、支持缩减及复合权系数界在稿件中分别证明。
所有有限次数最小值都控制同一贡献；先完成非负期望的完整素数占优再优化。
定义水平 D₀ 保持不变，不能将 κD₀ 误代入原 gap 公式。
原始矩用正项递推，平移及优化用 Arb 包含，未用浮点截零。

{evidence}

这些检查不能替代一般引理和外部已证筛定理的应用审查。
端点证书记录 finite_inputs_regenerated_in_this_run=false，
support_maxima_regenerated_in_this_run=true，并保存输入 SHA-256。
复用与重算的范围明确区分；粗支持基线 exp(31650) 保留可复现。

## 实际筛权余项与常数预算改进

{weighted_text}

## 拼接障碍与未完成项

{obstruction}

{limits}
"""
    for filename, value in [("conditional_chen_results_20261003.md", results),
                            ("conditional_chen_audit_20261003.md", audit)]:
        (ROOT/filename).write_text(value, encoding="utf-8", newline="\n")
    research_path = ROOT/"conditional_chen_research_20261003.md"
    historical = "**追加更正"+research_path.read_text(encoding="utf-8").split("**追加更正", 1)[1]
    research = f"""# 显式 Chen 两端拼接：条件性方案、数值筛查与待证接口

更新日期：2026 年 10 月 4 日。最新推导以 conditional_chen_results_20261003.md
和 conditional_chen_grh_bridge.tex 为准；下方追加更正之后保留较早探索记录。

最新稿件仅假设 GRH 的已认证阈值为 exp({current['GRH_positivity']['log_N0']})，exp(31000) 下系数提高到 {current_coefficient}，
使用保持所有先前检查的正项分块、六状态实际支持并集计数、整数粗筛及有限解析区间衔接
和 κ=11/169 的实际支持界；有限输入完整重算的粗支持基线为 exp(31650)，
当前费用复用带 SHA-256 的有限输入，另行扫描支持最大值和认证所有端点条件。
带未知 C、X 的分布路线在 3/4、9/10、19/20、99/100 的起点为
原固定误差预算的较弱计数界起点 774、210、165、141，较强计数界起点
786、212、166、142，κ=5/49 已计入费用。
96 个严格标量搜索候选用于有限参数选择；没有证明全局最优。
原发表版和本稿 GRH 路线只假设 GRH；EH 型路线不同时要求 GRH，
标准 EH 不提供具体数值 C、X。稿件仍需独立研究审稿，完整拼接仍未完成。

"""
    research_path.write_text(research+"## 当前 GRH 改进与认证\n\n"+current_text
                             +"\n\n## 前一组有限网格改进（保留对照）\n\n"+dense_text
                             +"\n\n## 分布常数与起点约束\n\n"+tradeoff_text
                             +"\n\n## 实际筛权余项与常数预算改进\n\n"+weighted_text
                             +"\n\n"+historical, encoding="utf-8", newline="\n")
    assert base["GRH"]["log_N0"] == 31650
    print("PASS: three current ledgers and the prefix-block note refreshed from checked certificates")


if __name__ == "__main__":
    main()
