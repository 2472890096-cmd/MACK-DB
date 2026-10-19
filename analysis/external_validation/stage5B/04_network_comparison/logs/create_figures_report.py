from pathlib import Path
from collections import Counter
import csv,datetime,hashlib,json,math,os,platform
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
ROOT=Path(__file__).resolve().parents[1]
def read(p):
 with Path(p).open() as f:return list(csv.DictReader(f,delimiter='\t'))
def write(p,rs):
 with Path(p).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rs[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rs)
def sha(p):
 with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def fmt(v,d=3):return 'NA' if str(v)=='NA' else f'{float(v):.{d}f}'
def pct(v):return 'NA' if str(v)=='NA' else f'{100*float(v):.2f}%'
def table(headers,rows):return '\n'.join(['| '+' | '.join(headers)+' |','| '+' | '.join('---' for _ in headers)+' |']+['| '+' | '.join(str(v) for v in r)+' |' for r in rows])
def main():
 primary=read(ROOT/'07_statistics/stage5B2_primary_summary.tsv');nulls={r['network']:r for r in read(ROOT/'06_null_model/stage5B2_degree_null_summary.tsv')};signs={(r['network'],r['analysis']):r for r in read(ROOT/'05_sign_concordance/stage5B2_sign_concordance_summary.tsv')};sens=read(ROOT/'07_statistics/stage5B2_wallen_threshold_sensitivity.tsv')
 assert len(primary)==5 and all(r['status']=='PASS' for r in nulls.values())
 labels=['Guo\nControl','Guo\nME/CFS','Guo\nCommon','Wallen\nPD','Wallen\nNHC'];colors=['#3568A8']*3+['#13877B']*2
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.titlesize':12,'axes.labelsize':11,'pdf.fonttype':42,'ps.fonttype':42,'savefig.dpi':400,'axes.spines.top':False,'axes.spines.right':False})
 def save(fig,name):
  for ext in ['png','pdf']:fig.savefig(ROOT/'08_figures'/f'{name}.{ext}',bbox_inches='tight',facecolor='white')
  plt.close(fig)
 fig,ax=plt.subplots(figsize=(8,4.4),layout='constrained');x=np.arange(5);vals=[100*float(r['published_edge_recovery']) for r in primary]
 ax.bar(x,vals,color=colors,width=.62)
 for i,r in enumerate(primary):ax.text(i,vals[i]+1.2,f"{r['shared_edges']} / {r['published_eligible_edges']}\n{vals[i]:.1f}%",ha='center',va='bottom',fontsize=10)
 ax.set(xticks=x,xticklabels=labels,ylabel='Published-edge recovery (%)',ylim=(0,66),title='A   Recovery of published associations');ax.set_axisbelow(True);ax.grid(axis='y',color='#E6E6E6');ax.axvline(2.5,color='#AAAAAA',ls=':',lw=.8)
 ax.legend(handles=[Patch(color=colors[0],label='Cross-framework'),Patch(color=colors[3],label='Cross-pipeline FastSpar')],frameon=False,loc='upper left',fontsize=9)
 fig.supxlabel('Denominator: eligible published edges in common mapped analysis taxa',fontsize=9)
 save(fig,'Figure_A_published_edge_recovery')
 fig,axes=plt.subplots(1,2,figsize=(9,4.3),layout='constrained',gridspec_kw={'width_ratios':[3,2]})
 for ax,inds,title in [(axes[0],[0,1,2],'Guo: cross-framework'),(axes[1],[3,4],'Wallen: cross-pipeline FastSpar')]:
  for j,i in enumerate(inds):
   r=primary[i];nr=nulls[r['network']];mean=float(nr['null_mean_shared']);lo=float(nr['null_2_5_percentile']);hi=float(nr['null_97_5_percentile']);obs=int(r['shared_edges'])
   ax.errorbar(j-.09,mean,yerr=[[mean-lo],[hi-mean]],fmt='o',color='#777777',capsize=5,label='Null mean + 95% interval' if j==0 else None)
   ax.scatter(j+.09,obs,s=65,color=colors[i],marker='D',zorder=3,label='Observed' if j==0 else None)
   ax.annotate(str(obs),(j+.09,obs),xytext=(0,8),textcoords='offset points',ha='center')
  ax.set(xticks=np.arange(len(inds)),xticklabels=[labels[i] for i in inds],title=title,ylabel='Shared edge count',xlim=(-.5,len(inds)-.5));ax.set_ylim(0,max(int(primary[i]['shared_edges']) for i in inds)*1.22);ax.grid(axis='y',color='#E6E6E6');ax.set_axisbelow(True)
 axes[0].legend(frameon=False,fontsize=8,loc='upper left');fig.suptitle('B   Observed overlap versus degree-preserving null',fontsize=13)
 fig.supxlabel('1,000 replicates per network; 50E successful swaps per replicate; empirical P = 0.001 for each',fontsize=9)
 save(fig,'Figure_B_observed_vs_degree_null')
 fig,ax=plt.subplots(figsize=(8,4.3),layout='constrained');same=[];opp=[];unknown=[]
 for r in primary:
  z=signs[(r['network'],'primary')];den=int(z['shared_edges']);same.append(100*int(z['same_sign'])/den);opp.append(100*int(z['opposite_sign'])/den);unknown.append(100*int(z['unknown_or_conflicting_sign'])/den)
 ax.bar(x,same,color='#277F8E',width=.62,label='Same sign');ax.bar(x,opp,bottom=same,color='#D87149',width=.62,label='Opposite sign')
 if any(unknown):ax.bar(x,unknown,bottom=np.array(same)+opp,color='#BDBDBD',width=.62,label='Unknown/conflicting')
 for i,r in enumerate(primary):ax.text(i,102,f"{signs[(r['network'],'primary')]['same_sign']} / {r['shared_edges']}\n{same[i]:.2f}%",ha='center',va='bottom',fontsize=9)
 ax.set(xticks=x,xticklabels=labels,ylabel='Shared edges (%)',ylim=(0,123),title='C   Sign concordance of shared associations');ax.set_yticks([0,25,50,75,100]);ax.legend(frameon=False,loc='upper right',ncol=2,fontsize=9)
 save(fig,'Figure_C_sign_concordance')
 fig,axes=plt.subplots(1,2,figsize=(8.4,4.2),layout='constrained');spearman={}
 for ax,n in zip(axes,['Wallen_PD','Wallen_NHC']):
  z=read(ROOT/'07_statistics'/f'{n}_shared_r_spearman.tsv')[0];spearman[n]=z
  shared=read(next(r['shared_pair_file'] for r in primary if r['network']==n));valid=[r for r in shared if r['published_r']!='NA'];vx=[float(r['published_r']) for r in valid];vy=[float(r['fastspar_r']) for r in valid]
  ax.scatter(vx,vy,s=15,alpha=.45,c='#13877B',edgecolors='none',rasterized=True);ax.plot([-1,1],[-1,1],ls='--',lw=1,color='#555555');ax.set(xlim=(-1,1),ylim=(-1,1),xlabel='Published FastSpar r',ylabel='MACK-DB rerun FastSpar r',title=n.replace('_',' '));ax.set_aspect('equal',adjustable='box')
  ax.text(.04,.97,f"n = {len(valid)}\nSpearman ρ = {float(z['spearman_rho']):.3f}\nP = {float(z['spearman_p']):.2e}",transform=ax.transAxes,va='top',fontsize=9,bbox={'facecolor':'white','alpha':.85,'edgecolor':'none'})
 fig.suptitle('D   Shared-edge correlation concordance',fontsize=13);fig.supxlabel('Same samples; different upstream abundance pipelines. Dashed line: y = x.',fontsize=9)
 save(fig,'Figure_D_Wallen_shared_r_concordance')
 # Record cause counts at two explicitly different denominators.
 causes=[];mapping=[];universe=[]
 for r in primary:
  n=r['network'];mapping.append(read(ROOT/'01_published_networks'/f'{n}_mapping_summary.tsv')[0]);universe.append(read(ROOT/'03_common_universe'/f'{n}_universe_audit.tsv')[0])
  raw=read(ROOT/'04_overlap'/f'{n}_edge_fate.tsv');label='primary' if n.startswith('Guo') else 'paper_matched';uniq=read(ROOT/'04_overlap'/f'{n}_{label}_unique_pair_fate.tsv')
  for scope,data in [('raw_published_records',raw),('unique_mapped_GTDB_pairs',uniq)]:
   counts=Counter(x['reason_if_not_recovered'] or 'recovered' for x in data)
   for reason,count in sorted(counts.items()):causes.append(dict(network=n,analysis='primary',counting_unit=scope,reason=reason,count=count,total=len(data)))
 write(ROOT/'07_statistics/stage5B2_edge_fate_reason_summary.tsv',causes)
 # User-requested sensitivity aliases make denominator and metric explicit.
 for r in sens:
  r['shared_published_edges']=r['shared_edges'];r['jaccard_or_endpoint_restricted_jaccard']=r['jaccard'];r['jaccard_scope']=r['jaccard_type']
 write(ROOT/'07_statistics/stage5B2_wallen_threshold_sensitivity.tsv',sens)
 # Frozen-input integrity: compare initial read hashes; check Stage 1A artifacts as well.
 inventory=read(ROOT/'00_inputs/input_manifest.tsv');integrity=[]
 priorpath=ROOT.parent/'02_prepare_fastspar_input/qc/stage5B1A_artifact_manifest.tsv';prior={r['path']:r['sha256'] for r in read(priorpath)}
 for r in inventory:
  digest=sha(r['path']);ok=digest==r['sha256_before'];older=prior.get(r['path'],'not_in_Stage1A_manifest');priorok=older=='not_in_Stage1A_manifest' or digest==older
  integrity.append(dict(path=r['path'],sha256_before=r['sha256_before'],sha256_after=digest,stage1A_prior_sha256=older,status='PASS' if ok and priorok else 'FAIL'))
 write(ROOT/'qc/stage5B2_input_integrity.tsv',integrity);assert all(r['status']=='PASS' for r in integrity)
 # Save reproducibility settings and scope, including intentionally different upstream inference parameters.
 config=dict(stage='Stage5B-2',python=platform.python_version(),numpy=np.__version__,matplotlib=matplotlib.__version__,networks=[r['network'] for r in primary],Guo_threshold='abs(r)>=0.3 and P<=0.0001',Wallen_primary_threshold='abs(r)>0.2 and P<0.05',Wallen_sensitivity_threshold='abs(r)>=0.3 and P<=0.0001',null_replicates=1000,successful_swaps_per_replicate='50E',max_attempts='max(10000,100*50E)',null_base_seed=20260930,null_sd_ddof=1,null_node_set='all common mapped analysis taxa including isolates',replicate_restart='original eligible published graph each replicate',empirical_p='(1 + count(null_shared >= observed))/(1000+1)',covariance_used_in_edge_analysis=False)
 (ROOT/'00_inputs/analysis_config.json').write_text(json.dumps(config,indent=2)+'\n')
 commands='''All writes restricted to 04_network_comparison. No FastSpar executable called.\nUse /home/wangmengyang/anaconda3/bin/python -s -B and set MPLCONFIGDIR=<output>/logs/matplotlib and NUMBA_CACHE_DIR=<output>/logs/numba.\nExecution order:\n1. logs/build_reference.py\n2. logs/stage5b2_analysis.py\n3. logs/run_degree_null.py\n4. logs/independent_recalculation.py\n5. logs/create_figures_report.py\nNull run initial JSON logging error was fixed by converting numpy.bool_ to Python bool; identical fixed seeds reused. No statistical settings changed. Initial error retained in logs/degree_null_initial_logging_error.log.\n'''
 (ROOT/'logs/commands.txt').write_text(commands)
 qc=[r for r in read(ROOT/'qc/stage5B2_qc_summary.tsv') if r['check'] not in ['frozen_input_integrity_and_Stage1A_manifest','output_nonempty_and_figure_delivery']]
 qc.append(dict(network='ALL',check='frozen_input_integrity_and_Stage1A_manifest',status='PASS',details=f'{len(integrity)} input files unchanged; all applicable Stage1A artifact hashes match',timestamp=datetime.datetime.now(datetime.timezone.utc).isoformat()))
 files=list(ROOT.rglob('*.tsv'))+list((ROOT/'08_figures').glob('*'))
 assert all(p.stat().st_size>0 for p in files)
 qc.append(dict(network='ALL',check='output_nonempty_and_figure_delivery',status='PASS',details=f'{len(files)} nonempty TSV/figure artifacts; four figures each PNG + PDF',timestamp=datetime.datetime.now(datetime.timezone.utc).isoformat()))
 write(ROOT/'qc/stage5B2_qc_summary.tsv',qc);assert not any(r['status']=='FAIL' for r in qc)
 # Manuscript text reflects eligible-denominator estimates, and never equates published edges with truth.
 manuscript=("Using identical biological sample sets, we compared published associations with MACK-DB MetaKSSD/GTDB R226 FastSpar reruns within common, harmonized analysis taxa. For Guo, cross-framework comparisons of published FlashWeave-S networks recovered 8/102 (7.84%), 17/87 (19.54%), and 31/130 (23.85%) eligible associations in controls, ME/CFS, and the combined cohort, respectively; all shared edges had concordant signs. For Wallen, paper-matched FastSpar thresholds recovered 431/865 (49.83%) PD and 596/1,175 (50.72%) NHC associations, with sign concordance of 99.30% and 98.49%. Shared-edge FastSpar coefficients were positively correlated across pipelines (Spearman rho = 0.671, n = 431, and 0.830, n = 595, respectively). Observed overlap exceeded the 97.5th percentile of degree-preserving null distributions for all five networks (1,000 independent-start randomizations each; empirical one-sided P = 0.001 each). Under the stricter MACK-DB threshold, Wallen recovery decreased to 19.54% and 18.21%, while sign concordance remained 100% and 99.53%. These results indicate nonrandom but incomplete association reproducibility across profiling/inference frameworks; they do not establish direct ecological interactions or causality.")
 (ROOT/'reports/manuscript_results.txt').write_text(manuscript+'\n')
 lines=['# Stage 5B-2 — Same-sample published-network versus MACK-DB FastSpar reproducibility analysis','',f'状态：**Stage 5B-2 COMPLETE**。UTC：{datetime.datetime.now(datetime.timezone.utc).isoformat()}。','',
 '## 比较范围与固定输入','',
 '本阶段未重跑FastSpar。直接读取Stage 5B-1B correlation/pvalue，保留原始矩阵、已冻结taxonomy和Stage 5B-0/0.5/1A/1B文件。样本数依次为91、106、197、490、234。所有结果仅写入本目录。',
 'Guo为**same-sample cross-framework / cross-pipeline reproducibility**：作者FlashWeave-S，与MACK-DB MetaKSSD/GTDB R226 + FastSpar相比，同时存在上游profiling、taxonomy、协变量处理和网络推断差异，不能归因于单一算法。Wallen为**same-sample cross-pipeline FastSpar reproducibility**：两侧均为FastSpar，但上游profiling/abundance/taxonomy不同，原作者与MACK-DB推断参数也并非完全一致（作者主拟合100 iterations、threshold 0.1、1000 permutations；本次原生产脚本50 iterations、exclude_iterations 20、threshold 0.2、10000 permutations）。这里只匹配Wallen的edge-selection criterion，不称为相同abundance矩阵的技术重复。','',
 '## 数字taxa ID恢复与published来源','',
 'Stage 1A的row-map中original/abundance/count ID均为数字，不能直接当物种名。本阶段逐行验证correlation和pvalue行列ID、count ID、abundance行ID及原project行ID完全一致且唯一；再按`GDTB226_species_abundance_result/species/<project>_species.matrix`查回production species label，去掉有证据的数字前缀，并通过冻结R226全部143614个物种名的唯一序列化逆映射恢复正式拼写。保留`_A`等GTDB后缀，不按词形猜物种。复核sample column map及Stage1A artifact hashes。详见`00_inputs/*_restored_taxa.tsv`与QC。',
 'Guo原始Table S7第10行区块标签、第26行Taxa1/Taxa2/Weight字段直接定位三个单次拟合对象，第27行起提取；Common=297、Control=210、ME/CFS=223，正负边数亦与原表一致。不使用AggControl/AggMECFS边集。作者边无附加weight阈值，FlashWeave weight仅保留描述与sign；不计算Guo weight-r Spearman。',
 'Wallen直接从原作者Supplementary Data 5/6恢复带方向的r和P，严格abs(r)>0.2且P<0.05后PD=1817、NHC=2509。逐边及r/P与Wallen_etal.xlsx和独立PD/NHC component文件完全一致，证实并非union或source-excluded网络。仅复用v44冻结taxon-level决策；映射后component union与冻结2268-pair reference一致只用于来源核验，不用于定义任何单组网络。','',
 '## Taxonomy harmonization','',
 'Guo使用Stage2B冻结决策；Wallen使用现有v44冻结rescue决策。完整tested universe中不在Wallen既有端点表的物种，只复用其他冻结表已有、名称确定性格式归一化后唯一的决策；没有既有决策、歧义或既有拒绝者保留unresolved，不重新做fuzzy matching或one-to-many expansion。每条决策都记录文件和记录ID。所有映射后self-loop删除；重复无向pair合并，保留来源列表；不同r不平均，存在多个不同published r时唯一pair的r=NA；sign不一致则sign=NA。','',
 table(['Network','Raw edges','Mapped/raw endpoints','Unresolved edge records','Self-loops','Duplicate records collapsed','Final unique GTDB pairs'],[[r['network'],r['raw_edges'],r['mapped_endpoints']+'/'+r['unique_raw_endpoints'],r['unresolved_edge_records'],r['self_loops_after_mapping'],r['duplicate_undirected_pairs_collapsed'],r['final_unique_gtdb_pairs']] for r in mapping]),'',
 '## 完整analysis universe及共同可比范围','',
 '五组均有完整作者analysis taxon身份集合的证据，而不是将final edge endpoints默认当tested universe。Guo Methods明确同一组100种，并明确Common网络包含全部100种；从Common原始单次拟合列表取回这些身份，Control保留其无边节点。Wallen补充表legend明确697/499种全部tested；在|r|筛选前读取完整P<0.05表，端点恰覆盖全部697/499种，因此完整身份可恢复。证据保存在`00_inputs/published_universe_evidence.tsv`及methods/legend文本。',
 '**完整raw tested身份集合不等于完整GTDB映射覆盖。** 下表明确列出unresolved物种；未解决的taxonomy与未进入FastSpar input的端点在计算recovery前排除。这里的Jaccard是**common-mapped-analysis-universe Jaccard**，只涉及共同可表示的分析物种，不是两个原始全taxa网络的unrestricted Jaccard。','',
 table(['Network','Published tested raw taxa','Mapped raw taxa','Unresolved raw taxa','Unique mapped GTDB taxa','FastSpar input taxa','Common taxa'],[[r['network'],r['published_analysis_raw_taxa'],r['published_analysis_mapped_raw_taxa'],r['published_analysis_unresolved_raw_taxa'],r['published_analysis_unique_gtdb_taxa'],r['fastspar_input_taxa'],r['common_taxa']] for r in universe]),'',
 '## Primary结果','',
 'Guo FastSpar：abs(r)>=0.3且P<=0.0001；Wallen primary：abs(r)>0.2且P<0.05，严格不等式保持原样。仅从i<j上三角提取，按GTDB物种字符串排序形成canonical无向pair。Published-edge recovery=shared/published eligible；FastSpar directional overlap=shared/FastSpar eligible，不称precision，作者网络不是gold standard。','',
 table(['Network','Raw pub','Final GTDB pub','Eligible pub','FS eligible','Shared','Recovery','FS→published overlap','Jaccard','Sign concordance'],[[r['network'],r['published_raw_edges'],r['published_final_gtdb_edges'],r['published_eligible_edges'],r['fastspar_edges'],r['shared_edges'],pct(r['published_edge_recovery']),pct(r['fastspar_to_published_overlap']),fmt(r['jaccard']),pct(r['sign_concordance'])] for r in primary]),'',
 'Final GTDB pub列包含映射成功但端点不在FastSpar input的pairs；它不是recovery分母。FS eligible仅计共同universe内边；全FastSpar输入taxa上的阈值边数另外保存为`fastspar_edges_all_input_taxa`，不混入Jaccard。','',
 '## Sign与Wallen r-value concordance','',
 table(['Network','PP','NN','PN','NP','Unknown','Same/shared'],[[r['network'],signs[(r['network'],'primary')]['positive_positive'],signs[(r['network'],'primary')]['negative_negative'],signs[(r['network'],'primary')]['positive_negative'],signs[(r['network'],'primary')]['negative_positive'],signs[(r['network'],'primary')]['unknown_or_conflicting_sign'],signs[(r['network'],'primary')]['same_sign']+'/'+r['shared_edges']] for r in primary]),'',
 table(['Network','Shared edges','Reliable unique-r pairs','Spearman rho','Spearman P','Nonunique r excluded'],[[n,z['shared_edges'],z['pairs_with_unique_reliable_published_r'],fmt(z['spearman_rho']),f"{float(z['spearman_p']):.4e}",z['excluded_nonunique_or_unavailable_r']] for n,z in spearman.items()]),'',
 'Wallen NHC一个shared GTDB pair合并了不同published r，因此只从Spearman/scatter排除，保留在edge presence和sign比较；原始r列表及来源仍在published edge文件中。相关性只描述已通过双方边选择的shared subset，不能推广为全部矩阵的一致性。','',
 '## Degree-preserving null','',
 '对共同universe中的eligible published graph做无向double-edge swaps，包含共同节点中的isolates。每次replicate从原图重新开始，保持节点集合、边数、每节点精确度数；禁止self-loop与duplicate；不额外限制连通性或正负度数。每组1000次，每次成功交换50E，maximum attempts=max(10000,100×50E)。固定base seed=20260930，SeedSequence生成每次独立seed。随机图档案、seed、成功交换数、attempts、success rate和逐次shared count均保留；5000/5000 replicate达到目标并通过度数与边数复核。','',
 table(['Network','Observed','Null mean','Null SD','Null 95% interval','Enrichment','z','Empirical P'],[[r['network'],r['observed_shared_count'],fmt(r['null_mean_shared']),fmt(r['null_sd_shared']),f"[{fmt(r['null_2_5_percentile'],2)}, {fmt(r['null_97_5_percentile'],2)}]",fmt(r['null_enrichment']),fmt(r['null_z']),r['empirical_p_display']] for r in nulls.values()]),'',
 '五组observed均高于null 97.5th percentile。Empirical one-sided P=(1+count(null>=observed))/1001；本次均为1/1001，报告P=0.001，不能写P<0.001。这里为每个comparison的未作多重检验调整的P。SD使用ddof=1；若SD=0则z=NA，若mean=0则enrichment=NA。随机化结果说明在此degree-preserving null下overlap高于随机期望，不证明生态直接作用。','',
 '## Wallen threshold sensitivity','',
 table(['Network','Analysis','FS eligible','Shared/eligible pub','Recovery','Jaccard','Sign concordance'],[[r['network'],'paper-matched' if 'primary' in r['analysis_label'] else 'MACKDB-standard',r['fastspar_edges'],r['shared_edges']+'/'+r['published_eligible_edges'],pct(r['published_edge_recovery']),fmt(r['jaccard']),pct(r['sign_concordance'])] for r in sens]),'',
 '作者published网络始终不变。严格数据库阈值下，PD/NHC recovery降至19.54%/18.21%，而sign一致性仍为100%/99.53%。因此部分跨pipeline边及方向具有稳定性，但恢复水平明显依赖边选择阈值；不能声称阈值改变后恢复率不变。本阶段null仅针对primary，未对sensitivity声称独立null显著性。','',
 '## 未恢复边的原因','',
 '以下按去重后的mapped GTDB pair计数；taxonomy未解决和self-loop损失见上方raw mapping表，不当作FastSpar biological failure。缺失input taxa者不进入eligible分母。',
 table(['Network','Mapped pair not in FS input','Eligible: r fails','Eligible: P fails','Eligible: both fail','Recovered'],[[n]+[next((r['count'] for r in causes if r['network']==n and r['counting_unit']=='unique_mapped_GTDB_pairs' and r['reason']==reason),0) for reason in ['taxon_not_in_fastspar_input','abs_r_below_threshold','p_above_threshold','both_r_and_p_fail','recovered']] for n in [r['network'] for r in primary]]),'',
 '逐条原因、原始物种、mapping状态、FastSpar r/P、是否过阈值、shared/sign均见`04_overlap/<network>_edge_fate.tsv`；该文件每个raw published edge一行，故涉及duplicate的记录不能直接加总为unique-pair recovery分母。另提供去重pair fate和`07_statistics/stage5B2_edge_fate_reason_summary.tsv`。','',
 '## 可以用于manuscript Results的摘要','',manuscript,'',
 '## 图与交付QC','',
 '- Figure A：eligible published-edge recovery，标注shared/published eligible。\n- Figure B：observed与degree-preserving null mean及中央95%分位区间；Guo/Wallen使用各自纵轴尺度。\n- Figure C：仅shared edges的same/opposite sign。\n- Figure D：Wallen primary shared-edge r散点与y=x；NHC排除1个非唯一published r pair。',
 '四图均提供400 dpi PNG和嵌入字体PDF，位于`08_figures/`。',
 f'独立程序未导入主分析或rewiring函数，重新读取最终edge TSV、原矩阵与全部随机图，完成134项数值回算，全部PASS。{len(integrity)}个登记输入SHA256前后不变；适用的Stage1A原manifest亦匹配。taxonomy不完整性以WARN明确记录，无关键FAIL。',
 '程序与执行顺序保存在`logs/`；analysis_config及输入SHA256在`00_inputs/`。初次null仅在输出JSON日志时遇到NumPy bool序列化错误，已修复并使用相同固定seed执行；没有更改方法、阈值或挑选结果。该初始错误日志保留。','',
 '**Stage 5B-2 COMPLETE**','']
 report=ROOT/'STAGE5B2_REPORT.md';report.write_text('\n'.join(lines))
 captions='\n\n'.join(['Figure A. Published-edge recovery within shared mapped analysis taxa; denominators exclude unresolvable or unavailable endpoints. Guo uses published FlashWeave-S edges versus MACK-DB-standard FastSpar; Wallen uses paper-matched FastSpar selection.','Figure B. Observed shared counts versus means and 2.5–97.5 percentiles of 1,000 degree-preserving randomizations (50E successful swaps each; independent original-graph starts). Panel scales differ.','Figure C. Same- and opposite-sign fractions among shared primary associations.','Figure D. Published versus rerun Wallen FastSpar r among shared primary pairs, with y=x. NHC omits one mapped pair with nonunique published coefficients from this panel only. All panels describe association reproducibility, not causality.'])
 (ROOT/'reports/figure_captions.txt').write_text(captions+'\n')
 required=[report,ROOT/'07_statistics/stage5B2_primary_summary.tsv',ROOT/'07_statistics/stage5B2_wallen_threshold_sensitivity.tsv',ROOT/'05_sign_concordance/stage5B2_sign_concordance_summary.tsv',ROOT/'06_null_model/stage5B2_degree_null_summary.tsv',ROOT/'qc/stage5B2_qc_summary.tsv',ROOT/'qc/stage5B2_independent_recalculation.tsv']
 required += [ROOT/'08_figures'/f'{stem}.{ext}' for stem in ['Figure_A_published_edge_recovery','Figure_B_observed_vs_degree_null','Figure_C_sign_concordance','Figure_D_Wallen_shared_r_concordance'] for ext in ['png','pdf']]
 assert all(p.is_file() and p.stat().st_size>0 for p in required)
 result=dict(status='COMPLETE',stage='Stage5B-2',timestamp=datetime.datetime.now(datetime.timezone.utc).isoformat(),networks=5,valid_null_replicates=5000,critical_QC_failures=0,independent_recalculation='PASS',input_integrity='PASS',report=str(report))
 (ROOT/'logs/completion.json').write_text(json.dumps(result,indent=2)+'\n');(ROOT/'logs/state.json').write_text(json.dumps(result,indent=2)+'\n')
 artifacts=[p for p in ROOT.rglob('*') if p.is_file() and 'matplotlib' not in p.parts and p.name!='stage5B2_delivery_manifest.tsv' and not p.name.endswith('.log')]
 write(ROOT/'qc/stage5B2_delivery_manifest.tsv',[dict(path=str(p),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(artifacts)])
 print('Stage 5B-2 COMPLETE')
 print(table(['Network','Published final GTDB edges','Published eligible edges','FastSpar edges','Shared edges','Published-edge recovery','Jaccard','Sign concordance','Null enrichment','Empirical P'],[[r['network'],r['published_final_gtdb_edges'],r['published_eligible_edges'],r['fastspar_edges'],r['shared_edges'],pct(r['published_edge_recovery']),fmt(r['jaccard'],4),pct(r['sign_concordance']),fmt(r['null_enrichment'],2),'0.001'] for r in primary]))
 for n,z in spearman.items():print(n,'Spearman rho',z['spearman_rho'])
 for p in required[:7]:print(p)
if __name__=='__main__':main()
