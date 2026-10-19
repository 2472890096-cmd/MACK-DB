"""Read back delivered matrices against immutable originals, then report completion."""
from prepare_inputs import *
import re

def matrix_rows(path):
    with Path(path).open(newline='') as f:return list(csv.reader(f,delimiter='\t'))

def main():
    extraction=readtsv(QC/'stage5B1A_sample_extraction_qc.tsv')
    conversion=readtsv(QC/'stage5B1A_count_conversion_qc.tsv')
    sources=readtsv(QC/'stage5B1A_resolved_sample_list_sources.tsv')
    sourceby={r['network']:r for r in sources}
    results=[]
    def check(network,test,ok,detail):
        results.append(dict(network=network,check=test,status='PASS' if ok else 'FAIL',details=str(detail)))
    for network,n in NETWORKS.items():
        src=sourceby[network];identities=readtsv(src['resolved_list'])
        runorder=[r['run_accession'] for r in identities]
        abundance=matrix_rows(SUB/(network+'_abundance.matrix'))
        counts=matrix_rows(COUNT/(network+'_count.matrix'))
        original=matrix_rows(src['source_matrix'])
        expected_cols=[r.get('column_index',r.get('mackdb_column_index')) for r in identities]
        positions=[original[0].index(col) for col in expected_cols]
        selected=[[row[0]]+[row[j] for j in positions] for row in original[1:]]
        maps=readtsv(QC/(network+'_sample_column_map.tsv'))
        check(network,'exact_sample_identity_and_order',abundance[0]==['OTU_id']+runorder and len(runorder)==len(set(runorder))==n,'Abundance headers equal exact frozen sample-list run IDs in the same order')
        check(network,'all_original_species_and_numeric_tokens_retained',abundance[1:]==selected,f'{len(selected)} original taxa; every selected numeric string identical; no filtering/normalization')
        check(network,'count_header_and_dimensions',counts[0]==['OTU_id']+[str(i) for i in range(n)] and len(counts)==len(abundance) and all(len(row)==n+1 for row in counts),'Original converter 0-based sample labeling and OTU_id header verified')
        check(network,'count_species_identity_and_order',[r[0] for r in counts[1:]]==[r[0] for r in original[1:]],'Input row IDs already sequential; original converter label reset preserves all original species IDs and order here')
        integer=finite=nonnegative=exact=True;checked=0
        for a,c in zip(abundance[1:],counts[1:]):
            for x,y in zip(a[1:],c[1:]):
                val=Decimal(y);target=Decimal(x)*Decimal('1000000');checked+=1
                finite=finite and val.is_finite()
                integer=integer and val==val.to_integral_value()
                nonnegative=nonnegative and val>=0
                exact=exact and val==target
        check(network,'integer_finite_nonnegative_exact_scaling',integer and finite and nonnegative and exact,f'{checked} cells checked with exact decimal arithmetic; no tolerance/rounding applied')
        check(network,'sample_column_map_readback',[r['run_accession'] for r in maps]==runorder and [r['count_sample_label'] for r in maps]==counts[0][1:] and [r['input_column_label'] for r in maps]==expected_cols,'Frozen run identity -> original column -> abundance header -> count column is reversible')
    # Original-file integrity is compared to the initial preflight baseline, including the failed first dependency attempt.
    baseline=json.loads((QC/'stage5B1A_input_snapshot.json').read_text())
    unchanged=all(Path(p).is_file() and sha(p)==r['sha256'] for p,r in baseline.items())
    check('ALL','original_inputs_and_scripts_unchanged',unchanged,f'{len(baseline)} original/frozen files match initial Stage5B-1A SHA256 baseline')
    check('ALL','five_extraction_and_conversion_QC_rows_pass',len(extraction)==len(conversion)==5 and all(r['sample_match_status']=='PASS' and r['missing_samples']=='0' and r['requested_samples']==r['matched_samples']==r['output_samples'] for r in extraction) and all(r['status']=='PASS' for r in conversion),'Requested/matched sample counts and all numeric conversion checks passed for all five networks')
    membership=readtsv(QC/'stage5B1A_group_membership_qc.tsv')
    check('ALL','sample_group_partitions_pass',all(r['status']=='PASS' for r in membership),'Control and ME/CFS are disjoint and their exact union is Common; Wallen PD/NHC disjoint')
    inferred=[p for p in ROOT.rglob('*.tsv') if p.name.endswith(('_correlation.tsv','_covariance.tsv','_pvalues.tsv'))]
    check('ALL','no_inference_outputs_created',not inferred,'No FastSpar correlation/covariance/pvalue files exist in this stage; command template is text only')
    write(QC/'stage5B1A_independent_readback_qc.tsv',results)
    event('Independent delivery readback against original matrices and frozen exact lists',STAGE,QC/'stage5B1A_independent_readback_qc.tsv',int(any(r['status']=='FAIL' for r in results)))
    if any(r['status']=='FAIL' for r in results):raise RuntimeError('Independent readback QC failed; do not mark complete')

    table=[]
    for r in extraction:
        n=r['network'];project=sourceby[n]['bioproject']
        table.append(f'| {n} | {project} | {r["requested_samples"]} / {r["matched_samples"]} / {r["missing_samples"]} | {r["taxa_count"]} | [{n}_abundance.matrix](abundance_subset/{n}_abundance.matrix) | [{n}_count.matrix](count_matrix/{n}_count.matrix) | PASS |')
    cells=sum(int(r['numeric_cells_checked']) for r in conversion)
    matrixpaths=[p for n in NETWORKS for p in [SUB/(n+'_abundance.matrix'),COUNT/(n+'_count.matrix')]]
    filelist='\n'.join('- `'+str(p)+'`' for p in matrixpaths)
    zerorows=sum(int(r['all_zero_taxa_retained']) for r in conversion)
    zerosamples=sum(int(r['all_zero_sample_columns']) for r in conversion)
    report=ROOT/'STAGE5B1A_REPORT.md'
    report.write_text(f'''# Stage 5B-1A — Exact sample subset abundance and count matrices

**状态：COMPLETE，2026-09-23。** 五组abundance subset与count matrix均已生成并通过QC。共检查{cells:,}个输出数值单元格；原始丰度字符串逐项不变，count严格等于abundance × 1,000,000，全部为有限非负整数。没有执行FastSpar、network overlap、Jaccard、sign concordance或任何taxa过滤。

## 输入与结果

| Network | Source BioProject | Requested / Matched / Missing samples | Taxa（原样保留） | Abundance subset | Count matrix | QC |
|---|---|---:|---:|---|---|---|
{chr(10).join(table)}

原始输入矩阵：

- Guo三组：[{AB}/PRJNA751448_abundance.matrix]({AB}/PRJNA751448_abundance.matrix)，197个原始sample columns、377个species rows。
- Wallen两组：[{AB}/PRJNA834801_abundance.matrix]({AB}/PRJNA834801_abundance.matrix)，725个原始sample columns、598个species rows。

Guo清单读取自Stage5B-0.5的`01_sample_rescue/sample_lists/Guo_{{Control,MECFS,Common}}_exact_samples.tsv`，逐项确认PASS_EXACT资格及精确本地列映射。Wallen按已冻结阶段的约定保留在`00_audit/Wallen_PD_exact_sample_list.tsv`与`00_audit/Wallen_NHC_exact_sample_list.tsv`；Stage5B-0.5没有重复复制这两份清单。本轮按清单中的BioProject和abundance路径自动定位PRJNA834801，未重新做样本rescue，也未向旧目录写文件。[清单解析来源](qc/stage5B1A_resolved_sample_list_sources.tsv)

## 样本身份与矩阵格式

样本匹配依据是**exact run accession + 已冻结的原始profile逐列验证索引 + 原ordered run list**，不是依据数量相等。每个run都核对原矩阵列号、原run list位置及原始矩阵路径。输入原列头为数字编号；输出abundance列头使用经验证的真实SRR run accession，列顺序严格沿用exact sample list的记录顺序。

所有species ID及其行顺序保持不变，数值字符串直接从原文件选列复制，未归一化、未log transform、未round、未删低丰度species。转换前的5组样本匹配均已通过；若有任何缺失，程序会在调用R前终止。

原始`count_sparcc.R`会将count表的列名重编号为0..n−1，也会重编号行ID。因此每组都保存独立的身份映射：

- `qc/<network>_sample_column_map.tsv`：参与者、BioSample、run、原始列、abundance列头、count列编号及清单源行。
- `qc/<network>_taxa_row_map.tsv`：原species ID、abundance行ID、count行ID和行位置。本次两份原矩阵的species ID本来就是连续编号，故count中的species ID及行顺序也完全不变；没有重新计算taxonomy。
- [合并样本映射](qc/stage5B1A_resolved_exact_samples.tsv)。不同count文件的数字列标签是各自的局部编号，不能跨文件直接当成同一个样本。

Guo Control与ME/CFS严格互斥，91+106人的精确run集合等于Common197人。Wallen PD/NHC严格互斥，共724人；原725列中未被两份冻结清单选中的run为`SRR19064765`，原列`449`，未补入任一组。[样本组QC](qc/stage5B1A_group_membership_qc.tsv)、[未选入run](qc/stage5B1A_source_runs_not_selected.tsv)

## 原始count conversion执行

五组均实际调用未修改的原脚本：

```bash
R_LIBS_USER="{LOG}/R_library" \\
TMPDIR="{LOG}/r_tmp" \\
/usr/bin/Rscript /home/wangmengyang/project/script/count_sparcc.R \\
  -i "{SUB}/<network>_abundance.matrix" \\
  -n 1000000 \\
  -o "{COUNT}/<network>_count.matrix"
```

初次转换发现系统R缺少`optparse`，错误原样保存在`logs/Guo_Control_count_conversion.initial_missing_optparse.stderr.log`。随后仅在本阶段`logs/R_library`安装[CRAN optparse 1.8.2](https://cran.r-project.org/web/packages/optparse/index.html)，使用R 4.5.1重新执行，五次最终转换退出码均为0。没有修改全局R库、原R脚本或原abundance文件。依赖源码URL、SHA256、安装输出、R版本与库路径分别记录在`logs/dependency_source_manifest.tsv`、安装日志和`logs/R_environment.txt`。

`-n 1000000`是原始方法的乘数，不是把每个sample列和重新归一化到1,000,000；输出列和随原始丰度而变化。本轮未改变这个含义，也未把这些缩放值声明为原始测序read counts。

## QC结果

- [样本提取QC](qc/stage5B1A_sample_extraction_qc.tsv)：5组PASS；requested=matched=output；missing=0；duplicate sample columns=0；原species顺序、所有丰度数值字符串完全一致。
- [count转换QC](qc/stage5B1A_count_conversion_qc.tsv)：5组PASS；species及sample数一致；首行`OTU_id`存在；非整数=0；精确乘数关系不符=0；最大绝对差=0；NA/非有限值=0；负值=0。共检查{cells:,}个数值单元格。
- [独立回读QC](qc/stage5B1A_independent_readback_qc.tsv)：{len(results)}项PASS，0 FAIL。直接重读冻结样本清单、原矩阵和交付文件，对身份、行序、数值及转换比例交叉核验。
- 本次子集中全零taxa行合计{zerorows}、全零sample列合计{zerosamples}；程序仅记录，不执行任何自动删除。
- [原始输入完整性](qc/stage5B1A_input_integrity.tsv)：{len(baseline)}个输入/原脚本SHA256在本阶段执行前后全部一致；[上一阶段冻结哈希核对](qc/stage5B1A_prior_frozen_input_checks.tsv)同时通过所有有历史基线的输入。
- [执行日志](logs/execution_log.tsv)记录command、timestamp、input、output、exit code；保留首次依赖失败及后续成功，不将失败历史覆盖。最终产物哈希见[artifact manifest](qc/stage5B1A_artifact_manifest.tsv)。

## FastSpar pipeline只读审计与下一步模板

已读取`/home/wangmengyang/project/script/run_fastspar.sh`，未执行。其三个FastSpar二进制的绝对路径均存在，本阶段只检查路径，没有执行二进制或`--help`。

| 项目 | 原pipeline设置 |
|---|---|
| 主拟合 | iterations=50；exclude_iterations=20；threshold=0.2 |
| Bootstrap表 | number=10000 |
| 每个bootstrap拟合 | 仅显式传`-i 5`；exclude_iterations、threshold、threads未显式设置，未臆测默认值 |
| P值 | permutations=10000 |
| 线程与并行 | 主流程T=CPU processor行数−3；bootstrap拟合通过GNU parallel，原脚本未指定job数 |
| Seed | 原脚本未显式设置 |
| 输入 | species×sample的tab分隔count OTU表，首行`OTU_id` |
| 输出 | `<PREFIX>_correlation.tsv`、`<PREFIX>_covariance.tsv`、`<PREFIX>_pvalues.tsv` |
| 清理 | 原脚本最后删除两个临时bootstrap目录；完整模板明确记录该行为 |

[完整命令模板](qc/stage5B1A_fastspar_command_template.txt)包含逐步命令以及五个已生成输入的绝对路径。将来手动运行时，替换输入/输出前缀即可，例如：

```bash
# 仅作为未来人工运行示例；本阶段未执行。未来输出父目录须预先存在。
bash /home/wangmengyang/project/script/run_fastspar.sh \\
  {COUNT}/Guo_Control_count.matrix \\
  /ABSOLUTE/PATH/FUTURE_OUTPUT_DIR/Guo_Control
```

本阶段只准备same-sample输入，不改变前阶段对Guo中间网络/Agg对象或Wallen方法比较范围的区分；没有计算任何网络或网络比较指标。

## 交付绝对路径

{filelist}
- `{report}`

本阶段到此停止。FastSpar运行、taxa筛选、taxonomy重新计算及网络比较均未启动。
''',encoding='utf-8')
    # Validate report links within this stage; all actual generated output paths remain inside ROOT.
    for link in re.findall(r'\]\(([^)]+)\)',report.read_text()):
        if link.startswith(('http://','https://')):continue
        p=Path(link) if link.startswith('/') else ROOT/link
        if p.name=='stage5B1A_artifact_manifest.tsv':continue # Written immediately below.
        assert p.exists(),'Broken local report link: '+str(p)
    event('Render STAGE5B1A_REPORT.md with final readback/QC and pipeline template',QC,report)
    paths=matrixpaths+[report]+sorted(ROOT.glob('*.py'))+sorted(QC.glob('*'))+[LOG/'dependency_source_manifest.tsv',LOG/'R_environment.txt']
    paths=[p for p in paths if p.is_file() and p.name!='stage5B1A_artifact_manifest.tsv']
    write(QC/'stage5B1A_artifact_manifest.tsv',[dict(path=str(p),size_bytes=p.stat().st_size,sha256=sha(p)) for p in paths])
    completed=dict(stage='Stage 5B-1A',status='COMPLETE',timestamp=timestamp(),networks=NETWORKS,
        extraction_pass=5,conversion_pass=5,independent_readback_pass=len(results),qc_failures=0,numeric_cells_verified=cells,
        immutable_input_files=len(baseline),fastspar_executed=False,network_comparisons_executed=False,report=str(report))
    (LOG/'completion.json').write_text(json.dumps(completed,indent=2)+'\n')
    (LOG/'preparation_state.json').write_text(json.dumps(completed,indent=2)+'\n')
    event('Stage 5B-1A COMPLETE',ROOT,report)
    terminal='Stage 5B-1A COMPLETE\n'+'\n'.join(str(p) for p in matrixpaths+[report])+'\n'
    (LOG/'completion_terminal.txt').write_text(terminal)
    print(terminal)

if __name__=='__main__':main()
