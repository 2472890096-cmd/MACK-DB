#!/usr/bin/env python3
"""Validate completed summary arithmetic and render the Stage 3C report."""
import csv
import hashlib
import json
import math
import statistics
from pathlib import Path

OUT = Path(__file__).resolve().parents[1]
BASE = OUT.parent
REPORT = BASE / 'reports/03c_null_model_validation.md'


def read(path):
    with Path(path).open(newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()


def write(path,rows):
    with Path(path).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n')
        w.writeheader();w.writerows(rows)


def quantile(values,q):
    values=sorted(values);index=(len(values)-1)*q;lo=math.floor(index);hi=math.ceil(index)
    return values[lo]+(values[hi]-values[lo])*(index-lo)


def near(x,y):
    return math.isclose(float(x),float(y),rel_tol=1e-12,abs_tol=1e-12)


def main():
    status=json.loads((OUT/'qc/03c_run_status.json').read_text())
    assert status['overall_status']=='PASS' and status['random_evaluations_completed']==18000
    rows=read(OUT/'summary/03c_null_model_summary.tsv')
    assert len(rows)==18
    protocol=read(OUT/'qc/03c_chen_wallen_protocol_consistency.tsv')
    assert all(r['status']=='MATCH' for r in protocol)
    independent=[]
    swap_audit=[]
    for r in rows:
        t=r['null_type'][-1];n=r['network'];E=int(r['edge_count'])
        path=OUT/('null_a' if t=='A' else 'null_b')/(n+'_null'+t+'_iterations.tsv')
        iterations=read(path)
        assert len(iterations)==1000
        counts=[int(x['recovered_count']) for x in iterations]
        rates=[x/E for x in counts]
        observed=int(r['observed_recovered']);mean=statistics.mean(rates);sd=statistics.stdev(rates)
        checks={
            'count_mean':near(r['null_mean_recovered'],statistics.mean(counts)),
            'count_sample_sd':near(r['null_sd_recovered'],statistics.stdev(counts)),
            'rate_mean':near(r['null_mean_recovery_rate'],mean),
            'rate_sample_sd':near(r['null_sd_recovery_rate'],sd),
            'rate_median':near(r['null_median_recovery_rate'],statistics.median(rates)),
            'linear_q025':near(r['null_q025_recovery_rate'],quantile(rates,.025)),
            'linear_q975':near(r['null_q975_recovery_rate'],quantile(rates,.975)),
            'minimum_maximum':near(r['null_min_recovery_rate'],min(rates)) and near(r['null_max_recovery_rate'],max(rates)),
            'enrichment':near(r['recovery_enrichment'],(observed/E)/mean) if mean else math.isinf(float(r['recovery_enrichment'])),
            'z_score':near(r['z_score'],((observed/E)-mean)/sd) if sd else math.isinf(float(r['z_score'])),
            'zero_sd_flag':r['null_sd_zero']==('Yes' if sd==0 else 'No'),
            'n_ge_observed':int(r['n_ge_observed'])==sum(x>=observed for x in counts),
            'corrected_P':near(r['empirical_p'],(1+sum(x>=observed for x in counts))/1001) and float(r['empirical_p'])>0,
        }
        for name,ok in checks.items():
            independent.append(dict(network=n,null_type=r['null_type'],check=name,status='PASS' if ok else 'FAIL'))
        if t=='B':
            successful=[int(x['successful_swaps']) for x in iterations]
            attempts=[int(x['swap_attempts']) for x in iterations]
            assert set(successful)=={50*E}
            swap_audit.append(dict(network=n,iterations=1000,E=E,required_successful_swaps_per_iteration=50*E,
                minimum_successful_swaps=min(successful),maximum_successful_swaps=max(successful),
                minimum_attempts=min(attempts),maximum_attempts=max(attempts),mean_attempts=statistics.mean(attempts),
                all_node_degrees_preserved='Yes' if all(x['per_node_degree_preserved']=='Yes' for x in iterations) else 'No',
                all_50E_successful_swaps='Yes',status='PASS'))
    write(OUT/'qc/03c_independent_summary_validation.tsv',independent)
    assert all(r['status']=='PASS' for r in independent)
    write(OUT/'qc/03c_null_b_swap_audit.tsv',swap_audit)
    integrity=read(OUT/'qc/03c_input_integrity.tsv')
    assert all(r['unchanged']=='True' for r in integrity)
    qc=read(OUT/'qc/03c_null_model_qc.tsv')
    assert all(r['status']=='PASS' for r in qc)
    # Report includes all required inputs and results, with explicit units and interval interpretation.
    graphs=read(OUT/'qc/03c_external_graph_audit.tsv')
    config=json.loads((OUT/'qc/03c_run_config.json').read_text())
    lines=['# External Validation Stage 3C', '', 'overall_status = PASS', '',
        '**18,000 / 18,000 randomized-network evaluations 已完成；Chen/Wallen 核心协议全部 MATCH。**', '',
        '## 正式协议及零标准差规则', '',
        '正式方法来自 `/home/wangmengyang/new_co_abundance_consistency/chen_wallen_null_models_final_v1.py`。'
        '执行脚本通过 AST 直接提取 `encode_edge`、`init_b_worker`、`null_b_worker` 与 `summarize`，未运行原脚本 main。'
        'Null B 的 RNG、交换提议、拒绝条件、成功次数循环和 recovery 计算保留原样，仅在最终返回处添加 graph QC。'
        '旧脚本与陈/Wallen 正式结果均未改动。函数来源、行号和 hash 见 `qc/03c_legacy_function_reuse.tsv`。', '',
        '前次停止原因已按用户新指令解决：null SD = 0 时，正式规则为 `z_score = inf`，并记 `null_sd_zero = Yes`；'
        '否则 `null_sd_zero = No`。`inf` 是旧 protocol 的定义，不作计算异常。前次停止审计保存在 `history/protocol_stop_before_zero_sd_resolution/`，不作为当前结果。', '',
        '## 参考数据库与 Observed', '',
        '每张网络使用 Stage 3B manifest 指定并已通过 QC 的 study-specific source-excluded MACK-DB，检验真实网络恢复是否依赖其自身数据来源之外的证据。'
        'Stage 3C 只读取派生数据库，没有重新执行 source exclusion，也没有把所有论文来源合并删除。Lacruz_MHO 无来源交集，按 manifest 复用 baseline 数据库。', '',
        'Observed 直接读取 Stage 3B summary 与逐 pair overlap，再与相应参考数据库进行独立 exact overlap 核对，全部一致。'
        'External pairs 从 Stage 3A manifest 读取；继续使用冻结 GTDB R226 species 名称，保留大小写和 suffix，不调用旧脚本的标签清洗函数，不重新进行 taxonomy mapping。', '',
        '| Network | N | E | Reference pairs | Observed recovered / E | Observed rate | Role |',
        '|---|---:|---:|---:|---:|---:|---|']
    for g in graphs:
        lines.append('| {} | {} | {} | {} | {}/{} | {:.4f}% | {} |'.format(g['network'],g['node_count'],g['edge_count'],
            g['source_excluded_database_pair_count'],g['observed_recovered'],g['edge_count'],100*float(g['observed_recovery_rate']),g['primary_or_exploratory']))
    lines += ['', '## Null A / Null B 定义', '',
        '**Null A**：固定外部图的 GTDB species 集合 V 与 E，在 V 的全部 N(N−1)/2 个 unordered pairs 中无放回均匀抽取 E 条；'
        '不排除真实边，不只从数据库已有 pair 抽样。图的顶点集合固定为 V，允许抽样后部分顶点成为 isolated nodes；不以“全部顶点必须在边中出现”为条件重新抽样。', '',
        '**Null B**：每个 iteration 独立从真实外部图开始，执行标准 degree-preserving double-edge swap。'
        '随机选择两条边和两种替代配对之一，拒绝会产生 self-loop 或重复边的提议。每次必须完成 **50 × E successful swaps**。'
        'attempt limit 沿用旧脚本 `max(target × 1000, target + 1,000,000)`，达到上限而未完成目标则 FAIL，不降低目标。', '',
        '两种 null 均每张网络 **1,000 iterations**，合计 9 × 2 × 1,000 = **18,000**；'
        'primary statistic 都是 recovered pair count / E，不使用 sign concordance，不按 positive/negative 分开随机。', '',
        '随机数固定为 `MASTER_SEED = 20260917`。使用旧方法的分组偏移方式扩展到九张网络：'
        '`seed = MASTER_SEED + 100000 × (2 × network_index + null_index) + iteration`，'
        '其中 network_index 按 run_config 的顺序从 0 开始，null_index 为 A=1、B=2，iteration=1…1000。18,000 个 seed 无重复。'
        '随机边序与候选池使用确定性排序；并行任务完成顺序不影响 seed 或结果，最终 iteration 文件按 iteration 排序。', '',
        '运行环境：Python {}；NumPy {}；pandas {}；{} 个计算进程。运行配置与 prepared input hash 保存在 `qc/03c_run_config.json`。'.format(
            config['python_version'].split()[0],config['numpy_version'],config['pandas_version'],config['workers']), '',
        '## 统计规则', '',
        '- 复用旧脚本的样本标准差（ddof=1）和 NumPy 默认线性插值的 2.5% / 97.5% quantiles。',
        '- 95% interval 是 **null distribution 的经验百分位区间**，不是 null mean 的置信区间；paper table 的 `95CI` 字段沿用任务要求命名。',
        '- enrichment = observed recovery rate / null mean recovery rate。',
        '- SD > 0 时 z = (observed rate − null mean rate) / null SD rate；旧脚本用等价 count 形式计算，已独立核对。SD = 0 时 z = inf。',
        '- `n_ge_observed` 为 null recovered count ≥ observed count 的次数（E 固定，与 rate 比较等价）。',
        '- 单侧 empirical P = (n_ge_observed + 1) / 1001，最小为 1/1001 = 0.000999000999…，不报告 P=0。',
        '- Null A 与 Null B 分别报告 P；未合并 P，未引入 BH-FDR 作为筛选标准，未按显著性删除结果。', '',
        '## Null A 结果', '',
        '| Network | Mean recovery | 95% interval | Enrichment | z | null_sd_zero | n_ge | Empirical P |',
        '|---|---:|---|---:|---:|---|---:|---:|']
    def render(r):
        return '| {} | {:.4f}% | {:.4f}%–{:.4f}% | {:.6f} | {:.6f} | {} | {} | {:.12g} |'.format(
            r['network'],100*float(r['null_mean_recovery_rate']),100*float(r['null_q025_recovery_rate']),100*float(r['null_q975_recovery_rate']),
            float(r['recovery_enrichment']),float(r['z_score']),r['null_sd_zero'],r['n_ge_observed'],float(r['empirical_p']))
    lines += [render(r) for r in rows if r['null_type']=='Null A']
    lines += ['', '## Null B 结果', '', '| Network | Mean recovery | 95% interval | Enrichment | z | null_sd_zero | n_ge | Empirical P |',
        '|---|---:|---|---:|---:|---|---:|---:|']
    lines += [render(r) for r in rows if r['null_type']=='Null B']
    nonsig=[r for r in rows if float(r['empirical_p'])>0.05]
    lines += ['', '**Empirical P > 0.05：** ' + ('；'.join('{} {} (P={:.12g})'.format(r['network'],r['null_type'],float(r['empirical_p'])) for r in nonsig) if nonsig else '无。'), '']
    zero=[r for r in rows if r['null_sd_zero']=='Yes']
    lines += ['**Null SD = 0：** '+ ('；'.join(r['network']+' '+r['null_type'] for r in zero)+'。这些 z=inf 来自正式 Chen/Wallen protocol 的零标准差边界规则。' if zero else '无；18 项结果的 null_sd_zero 均为 No。'), '',
        'Lacruz_MHO 保留 exploratory 标记；其余网络角色保持不变。', '',
        '## Null B 交换完成情况', '',
        '| Network | Iterations | 每次要求与实际成功 swaps | Attempts min–max | Individual degrees preserved |',
        '|---|---:|---:|---|---|']
    for r in swap_audit:
        lines.append('| {} | {} | {} | {}–{} | {} |'.format(r['network'],r['iterations'],r['required_successful_swaps_per_iteration'],r['minimum_attempts'],r['maximum_attempts'],r['all_node_degrees_preserved']))
    lines += ['', '9,000 个 Null B iterations 全部达到精确 50E 次成功交换，总成功交换次数 {:,}。'.format(status['total_successful_swaps']), '',
        '## QC 与冻结完整性', '',
        '- 核心 protocol consistency：{} 项全部 MATCH。'.format(len(protocol)),
        '- 执行 QC：{} 项全部 PASS。'.format(len(qc)),
        '- 独立 summary 算术检查：{} 项全部 PASS（使用标准库样本标准差、独立线性 quantile 与尾部计数核对）。'.format(len(independent)),
        '- 每个随机图检查 E、self-loop、duplicate edge、node pool 与 exact matching；每个 Null B 图额外逐节点验证 degree vector，并验证排序 degree sequence 与完整 node set。',
        '- 保留每个 network × null_type 的 iteration 1、500、1000 完整边表，共 54 份；其 edge hash、overlap 和相应 degree 检查全部通过。未保存全部 18,000 份完整图。',
        '- 全部 {} 个输入文件 hash 前后一致，包括 Stage1/2/2B/2C/3A/3A.1/3B 正式文件和旧 Chen/Wallen 方法证据；Stage3B 所有 reference hashes 与正式 manifest 一致。'.format(len(integrity)),
        '- 原始 human_gut.info SHA256 保持 `4b67935940a765d28518a2dd612587f347c0e47e985323910fcdde51950a7834`。',
        '- 无新的 taxonomy、source exclusion、FastSpar、cross-paper recurrence 或其他后续分析。', '',
        '## 正式输出路径', '', '- 报告：`'+str(REPORT)+'`',
        '- 全部绝对路径、大小和 SHA256：`'+str(OUT/'qc/03c_formal_output_manifest.tsv')+'`', '',
        '以下路径均相对于 `'+str(OUT)+'`：', '']
    for p in sorted(OUT.rglob('*')):
        if p.is_file() and not any(x in p.relative_to(OUT).parts for x in ['history','prepared','snapshots']) and p.name!='03c_formal_output_manifest.tsv':
            lines.append('- `'+str(p.relative_to(OUT))+'`')
    lines += ['', 'Prepared 图/参考节点池文件在 `prepared/`；54 份代表图在 `snapshots/`；前次停止记录在 `history/`。这些文件均逐项列入输出 manifest，历史记录单独标记。', '',
        '全部输出 TSV 中的 recovery rate 使用 0–1 单位；本报告表格转换为百分比。完成 Stage 3C 后停止。', '']
    REPORT.write_text('\n'.join(lines))
    mf=OUT/'qc/03c_formal_output_manifest.tsv'
    paths=[p for p in sorted(OUT.rglob('*')) if p.is_file() and p!=mf and p.suffix!='.tmp']+[REPORT]
    write(mf,[dict(absolute_path=str(p),size_bytes=p.stat().st_size,sha256=sha(p),
        role='historical_stop_audit' if '/history/' in str(p) else 'current_Stage3C_output') for p in paths])
    print('Report:',REPORT)
    print('Output manifest:',mf,'entries:',len(paths))
    print('Independent summary checks:',len(independent),'PASS')


if __name__=='__main__':main()
