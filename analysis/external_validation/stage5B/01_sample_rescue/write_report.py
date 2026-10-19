from rescue_io import *
import collections

def main():
    summary=rows(R/'stage5B_05_sample_rescue_summary.tsv');by={r['network']:r for r in summary}
    qc=rows(Q/'stage5B_05_sample_rescue_qc.tsv');counts=collections.Counter(r['status'] for r in qc)
    integrity=rows(Q/'frozen_input_integrity.tsv');manifest=rows(E/'stage5B_05_web_evidence_manifest.tsv')
    ready=[r for r in summary if r['next_stage_ready']=='YES']
    readyrows=[]
    for r in ready:
        readyrows.append(dict(network=r['network'],classification=r['classification'],participants=r['exact_samples_recovered'],runs=r['exact_runs_recovered'],
            exact_sample_list=str(S/(r['network']+'_exact_samples.tsv')),public_network_definition='Table S7 '+{'Guo_Control':'E:G (Control Network)','Guo_MECFS':'H:J (ME/CFS Network)','Guo_Common':'B:D (Common Network)'}[r['network']],
            published_edge_source=str(OLD/'raw/guo2023/1-s2.0-S193131282300029X-mmc8.xlsx'),eligibility_scope='sample provenance and local availability only',
            next_stage_constraint='Use published intermediate edge object, not frozen AggControl/AggMECFS; covariate/conditional-weight comparison protocol still needs specification; no fitting performed in Stage5B-0.5'))
    write(R/'stage5B_05_next_stage_ready.tsv',readyrows)
    recommended=[]
    order=['Chen','Chen_LLD','Chen_500FG','Chen_300OB','Chen_IBD','Zhou_MS','Zhou_HHC','Guo_Control','Guo_MECFS','Guo_Common','Guo_AggControl','Guo_AggMECFS','Lacruz_MHNO','Lacruz_MHO','Lacruz_MUNO','Lacruz_MUO','Li_IBS']
    reasons={
        'Chen':'四队列独立拟合后的 union；最终筛选身份未落实',
        'Chen_LLD':'公开 EGA 映射1135；没有独立证据证明最终拟合集合等于整个数据集',
        'Chen_500FG':'原项目471 BioSamples/875 runs；最终450选择未知；另517为后来重测',
        'Chen_300OB':'公开 EGA297；论文网络298，缺最终筛选及本地输入',
        'Chen_IBD':'公开 EGA355；论文网络496，缺最终筛选及本地输入',
        'Zhou_MS':'648 runs仅本地29；209人仅3人所需runs齐全',
        'Zhou_HHC':'684 runs仅本地37；209人仅4人所需runs齐全',
        'Guo_Control':'Table S7独立Control拟合；91份首访样本全部存在',
        'Guo_MECFS':'Table S7独立ME/CFS拟合；106份首访样本全部存在',
        'Guo_Common':'Table S7全197份样本拟合；全部存在',
        'Guo_AggControl':'197名父级参与者齐全；Common+验证后组特异边，非单次拟合，具体子抽样未恢复',
        'Guo_AggMECFS':'197名父级参与者齐全；Common+验证后组特异边，非单次拟合，具体子抽样未恢复',
        'Lacruz_MHNO':'完整103人；本地5人', 'Lacruz_MHO':'完整19人；本地0人',
        'Lacruz_MUNO':'完整377人；本地25人','Lacruz_MUO':'完整432人；本地20人',
        'Li_IBS':'13个discovery项目；最终逐样本输入未知，当前仅其中1个项目有本地矩阵'}
    for n in order:
        r=by[n];sample=r['exact_samples_recovered']
        if sample=='0':sample='0 certified（候选另列）'
        if n.startswith('Guo_Agg'):sample='197 parent participants'
        recommended.append(f'| {n} | {r["classification"]} | {sample} | {r["samples_found_in_local_abundance"]} | '+('Yes，下一阶段' if r['next_stage_ready']=='YES' else 'No')+f' | {reasons[n]} |')
    text=f'''# Stage 5B-0.5 — Published exact sample-set rescue

完成日期：2026-09-23。审计范围为 Chen2020、Zhou2022、Guo2023、Lacruz2026、Li2024。Stage 5B-0 保持冻结，Wallen 仅引用已完成结果。本轮新增可用的 exact sample inputs 为 **Guo_Control（91）、Guo_MECFS（106）、Guo_Common（197）**，对应作者已经发表的三个中间网络。

**Zhou 的209 MS和209 HHC公开身份链已完整恢复，但本地数据不足；Lacruz完整931人已确认，本地仅50人。Chen、Li仍无法严格重建最终网络样本集。** 本轮没有执行 FastSpar、其他网络拟合、网络 overlap、Jaccard、randomization、taxonomy 修改或冻结结果修改。

## 1. 输出与判定口径

- [逐网络汇总](01_sample_rescue/stage5B_05_sample_rescue_summary.tsv)：17行，含10个本轮涉及的冻结网络输入、4个Chen队列组成对象、3个Guo已发表中间网络。`audit_scope` 明确区分。
- [逐样本/run主表](01_sample_rescue/stage5B_05_sample_master.tsv)：13,948行；同一run在不同网络可分别出现，同一网络内不重复。
- [独立样本清单目录](01_sample_rescue/sample_lists/)；`*_exact_samples.tsv` 中只有证据充分的样本，或一个明确的 `status=unresolved` 空身份占位行。未确认的公开候选另存 `*_unresolved_candidates.tsv`。
- [网络资源证据索引](01_sample_rescue/evidence/stage5B_05_web_evidence_manifest.tsv)：{len(manifest)}条去重资源记录，保留下载失败、挑战页面、替代来源及SHA256。
- [样本救援QC](01_sample_rescue/qc/stage5B_05_sample_rescue_qc.tsv)：PASS={counts['PASS']}、INFO={counts['INFO']}、WARN={counts['WARN']}、FAIL={counts['FAIL']}。
- [下一阶段可用输入](01_sample_rescue/stage5B_05_next_stage_ready.tsv)、[参与者完整输入覆盖](01_sample_rescue/evidence/participant_complete_input_coverage.tsv)、[所有本地额外run去向](01_sample_rescue/evidence/local_extra_run_disposition.tsv)。

`exact_samples_recovered` 以参与者/生物学输入为单位；`exact_runs_recovered` 是测序run，不能互换。`samples_found_in_local_abundance` 只计**所有已知必需runs均在本地**的完整参与者输入；另列 `participants_with_any_local_run` 和 `local_runs_recovered`。`recovery_fraction` 为公开参与者身份恢复比例，`local_complete_fraction` 才是完整本地覆盖比例；公开恢复100%不等于可重建。

U 表示**最终网络纳入关系未解决**，不是公开accession本身不存在。Chen/Li 的0 exact表示本轮没有认证最终拟合身份，并非其公开队列没有样本，也不代表实际重叠为0。无法判断本地最终样本数时写 `unresolved`。E5逐行解释唯一证据链，不按疾病模糊匹配或按期望数量凑组。

Guo Agg两行有完整父级集合，但没有单一拟合集合，且作者验证子抽样的具体成员未恢复。本报告按 **PARTIAL** 保守标注其复合拟合来源；这不表示197份首访样本丢失，也不把同一BioProject错误描述为多个cohort。

## 2. Chen2020

论文：[Nature Communications](https://www.nature.com/articles/s41467-020-17840-y)。作者在四个队列分别建网，再比较/汇总；现有冻结Chen边集是跨队列union，不能视作PRJNA319574或PRJNA942468上的单次拟合。

本轮成功取得[作者网络代码](01_sample_rescue/cache/Chen_author_network.R)，路径 `Projects/Microbial co-abundance network/co-abundance network main codes.R`，master页面提交为 `723045925f0ef87a39587d27837e1a3f00429cf4`。代码明确SPIEC-EASI `glasso`，与SparCC结果联合筛选；脚本读通用输入路径，没有实际样本清单。脚本使用自定义经验edge FDR，不能将其自动写成BH校正；本轮只纠正记录，没有执行代码。

| Cohort | 原始/公开来源规模 | 论文网络人数 | 公开来源与本轮映射 | exact最终拟合集合 | 本地情况 |
|---|---|---:|---|---|---|
| LLD | 招募1500；EGA描述1200 | 1135 | EGAD00001001991 / EGAS00001001704；1135个EGAN与1135个EGAR、subject ID | 未认证 | 没有已验证的本地EGA run列映射 |
| 500FG | 招募534；原项目471 BioSamples、875 runs | 450 | PRJNA319574；另PRJNA942468有517份后来重测 | 未认证最终450选择 | 原项目466 runs；后来重测517 runs，均不能按数量代替最终450 |
| 300OB | 招募302 | 298 | 论文DAC EGAC00001001143进一步定位EGAD00001005083 / EGAS00001003508；公开297个样本/run | 未认证 | 无已验证本地输入 |
| IBD | EGA文字描述461、当前公开355 | 496 | EGAD00001004194 / EGAS00001002702；355个样本/run及subject ID | 未认证 | 无已验证本地输入 |

EGA的**公开元数据已实际获取**，不是因受控数据访问就停止搜索。使用官方公开API保存完整dataset→sample→subject和dataset→run→sample映射；EGAN是EGA样本accession，不能冒充NCBI BioSample；EGAR也不能当SRR/ERR。LLD数目一致仍缺“最终网络筛选等于当前公开dataset全集”的独立证据，因此仍保留U。见[队列选择审计](01_sample_rescue/evidence/Chen_cohort_selection_audit.tsv)及各 `Chen_*_unresolved_candidates.tsv`。[EGA公开API说明](https://ega-archive.org/discovery/metadata/public-metadata-api/)

PRJNA942468的[官方项目描述缓存](01_sample_rescue/cache/Chen_2024_ENA_study.xml)说明：同一500FG粪便样本再次测序，提取方案为QIAamp Fast DNA Stool Mini Kit；2020论文则明确使用AllPrep DNA/RNA Mini Kit，并将500FG数据指向PRJNA319574。这些直接来源差异支持不将PRJNA942468认证为原始建网runs或替代输入。其首次公开日期为2023年，仅作背景；公开日期本身不能证明测序日期。原PRJNA319574的404个BioSample有两个runs、67个有一个run，875 runs不能按875人解释。即使能认出重测的HV编号，也不能将其替换成原始测序输入。

新增检查了官方[Source Data](01_sample_rescue/cache/Chen_Source_Data.xlsx)和[补充文件说明](01_sample_rescue/cache/Chen_supplement_descriptions.pdf)：Source Data包含队列标签与图表数值，没有逐样本身份，不能按Excel行序配对。补充文件说明也纠正了旧审计中一个文件标号：MOESM4是Supplementary Data1的species网络，MOESM6才是Data3；旧文件未改动。

结果：四个组成队列均 **FAIL**（最终网络成员未唯一认证），原Chen union为 **MULTI_COHORT_NOT_RECONSTRUCTABLE**。不是因为875/517/1135等公开身份没有救回，而是没有足够证据将它们认证为最终拟合集合。

## 3. Zhou2022

论文：[Cell](https://doi.org/10.1016/j.cell.2022.08.021)。关键新增来源是**正式出版Table S1**，正确PII为S0092867422011151；[官方xlsx](https://ars.els-cdn.com/content/image/1-s2.0-S0092867422011151-mmc2.xlsx)已保存到[本地](01_sample_rescue/cache/Zhou_publisher_Table_S1.xlsx)。其中Sequenced samples有7801条记录，Sample phenotype有1152名参与者，较旧作者仓库样本表更完整。

选择规则完全使用作者字段：MS=`disease=MS AND treatment_status=Untreated`；HHC=`disease=Control AND treatment_control=Untreated_HHC`。精确连接 `iMSMS_ID` 后只保留 `target_gene=WGS AND paired_used=Yes`；再逐条核对public Run→BioSample及experiment accession。不是先要求209再凑数。

| Network | 公开参与者 | 公开sample/BioSample | 精确runs | 本地runs | 至少1个run的参与者 | 全部runs齐全的参与者 | 分类 |
|---|---:|---:|---:|---:|---:|---:|---|
| Zhou_MS | 209 | 507 | 648 | 29 | 25 | 3 | PARTIAL |
| Zhou_HHC | 209 | 518 | 684 | 37 | 34 | 4 | PARTIAL |

两组都完整恢复209名参与者；run与BioSample逐条唯一，MS/HHC的参与者、BioSample、run三层均互斥。HHC相比旧阶段36 runs/33人，新增了1个run及1名可见参与者，但仍远不足以重建完整209人输入。

论文WGS方法明确将每名参与者各份样本的sequencing counts相加。完整清单保留采样日期、wet/dry、Cohort、publication sample ID、BioSample和全部runs，不能选一个“代表run”，也不能把不同采样/制备/重复测序当作独立参与者。这里判断本地可重建时要求所有必需runs齐全，因此是3/209和4/209，而不是25/209、34/209。将来如数据补齐，应按作者的counts汇总语义处理，不能直接将不同run的相对丰度相加。

出版表将其中MS112和HHC116个runs标为WGS，但公共元数据的library_strategy为OTHER；它们的BioSample、experiment、run身份仍精确一致。本轮保留作者明确的WGS注释，并在每行标记差异，没有按通用library_strategy错误删除。

作者[GitHub仓库](https://github.com/BaranziniLab/iMSMS_study)已检查脚本、metadata及R对象。缓存seqmeta.rda只有1151个元数据对象，而taxa.list中的species对象有1152列；3个新缓存对象/工作流的Git blob哈希与已保存树 `1638c6a0c05966b5a25582190f97875fce2fad5a` 一致。[检查表](01_sample_rescue/qc/repository_content_identity.tsv)

实际网络脚本引用 `species.treatment.rda`、`seqmeta1152.rda`，在已检查树中没有这些最终对象；因此本轮恢复的是**正式出版定义的完整组及样本链**，不声称拿到最终执行时序列化的拟合对象。Dryad DOI10.7272/Q60C4T26的版本和文件列表已获取，文件下载遇到401/403后改用出版商正式Table S1；不声称两版本文件二进制相同。

两组都**不能**以当前本地矩阵做same-sample重建，不能PASS_EXACT。清单：[MS](01_sample_rescue/sample_lists/Zhou_MS_exact_samples.tsv)、[HHC](01_sample_rescue/sample_lists/Zhou_HHC_exact_samples.tsv)，文件内 `status=partial` 表示本地不完整。

## 4. Guo2023

论文：[Cell Host & Microbe](https://doi.org/10.1016/j.chom.2023.01.004)。本轮直接下载了PRJNA751448对应的197份[BioSample属性XML](01_sample_rescue/cache/Guo_direct_BioSample_attributes.xml)，并与独立SRA快照及ENA run映射交叉核对。

逐样本分组使用作者提交的 **SampleType=Control/Case**，同时检查 `mecfs_status=N/Y` 一致、`participant_id=Original_ID`、唯一 `STOOL_ID`。论文明确只分析第一次随访，并明确Control网络使用所有91份对照、ME/CFS网络使用所有106份病例、Common网络使用全部197份。因此证据链是：论文的全组/首访定义与该论文数据提交 → 明确作者分组标签和粪便/参与者ID → 唯一BioSample → SRA/ENA Run → 原始profile逐列核实的本地矩阵。数量91/106是检验结果，不是分组依据。

BioSample标签本身为E4；完整sample-to-network assignment跨论文、公开元数据和本地列证据，保守记E5，每行均写明唯一推断链。[逐样本属性核对](01_sample_rescue/evidence/Guo_explicit_group_attribute_audit.tsv)

| 已发表单次拟合对象 | 参与者=样本=BioSample=runs | 本地完整覆盖 | Table S7位置 | 分类 |
|---|---:|---:|---|---|
| Guo_Control | 91 | 91 | E:G，Control Network | PASS_EXACT |
| Guo_MECFS | 106 | 106 | H:J，ME/CFS Network | PASS_EXACT |
| Guo_Common | 197 | 197 | B:D，Common Network | PASS_EXACT |

197个参与者、STOOL_ID、BioSample、run分别唯一，所有runs属于PRJNA751448，Control/ME/CFS严格互斥，本地全部覆盖。可用清单：[Control](01_sample_rescue/sample_lists/Guo_Control_exact_samples.tsv)、[ME/CFS](01_sample_rescue/sample_lists/Guo_MECFS_exact_samples.tsv)、[Common](01_sample_rescue/sample_lists/Guo_Common_exact_samples.tsv)。本轮没有提取或重新处理网络边。

必须保留原冻结Agg对象与中间对象的区别。原始补充表 `{OLD}/raw/guo2023/1-s2.0-S193131282300029X-mmc8.xlsx` 第10行给出对象定义，第12行报告Common/Control/ME/CFS的边数为297/210/223。AggControl和AggMECFS则分别结合Common网络与经过验证的组特异边，作者报告388/358条；这些是读取作者表格的定义与计数，没有计算网络重叠。

因此，把91份对照跑一次FastSpar后直接对比冻结AggControl，会把含全197份样本拟合信息的边表误认作纯对照拟合。两张Agg表的197名父级参与者全部已恢复，但构造流程包括组特异验证子抽样，具体抽样成员/种子没有恢复。本轮对Agg对象标 **PARTIAL**，`sample_set_status` 明确写为complete parent memberships / composite not single fit / subsampling memberships unresolved。

三项PASS仅证明**样本来源和本地可提取性**。FlashWeave的条件关联、作者协变量及权重语义与边际SparCC相关不等价；下一阶段必须先明确对比协议，不能把“样本一致”升级为“仅推断算法不同”。本轮不会将作者中间边表静默替换成冻结GTDB Agg边表。

## 5. Lacruz-Pleguezuelos2026

论文：[Nature Communications](https://www.nature.com/articles/s41467-026-72588-1)。复核已冻结Supplementary Data3的931名参与者及原始sample-source manifest，并下载[作者v1.0存档](https://zenodo.org/records/19484218)，归档commit前缀fd8c517。2.13MB zip的MD5与Zenodo记录一致；版本内whole_cohort_data.csv的SHA256与已有作者丰度文件相同，931个指定输入样本都在作者列名中。

网络脚本读取 `physeq_adj.rds`，按MetObesity四个phenotype筛选样本，再筛taxa。batch处理保留合并metadata，无法因此给缺失的本地参与者补数。AI4Food多visit按作者metadata脚本的phenotype/visit规则选每人一份；本轮逐条沿用冻结manifest的subject、visit和run身份，未自行改为baseline或选任意run。MetaCardis一名参与者可能对应多个技术测序runs，这些runs不能当多个人。

| Cohort/project | 论文参与者 | 对应runs | 本地目标runs |
|---|---:|---:|---:|
| Karlsson / PRJEB1786（ERP002469） | 145 | 145 | 0 |
| MetaCardis / PRJEB37249（ERP120556） | 249 | 1360 | 0 |
| MetaCardis / PRJEB38742（ERP122188） | 378 | 1814 | 0 |
| Feng / PRJEB7774（ERP008729） | 61 | 61 | 50 |
| AI4Food / PRJEB87701（ERP170896） | 98 | 98 | 0 |

所有五个公开来源都检查了canonical BioProject及明确的ERP别名在本地abundance根目录的文件名，实际只有PRJEB7774存在。另发现相关MetaCardis项目PRJEB41311（ERP125059）本地41列，对238个公开runs检查后，与完整931人输入的run、BioSample、精确alias均无交集，因此不能仅凭MetaCardis名称补入。[五来源审计](01_sample_rescue/evidence/Lacruz_all_source_cohorts.tsv)、[项目及别名本地检查](01_sample_rescue/evidence/local_project_and_alias_audit.tsv)

| Phenotype | 完整参与者/BioSample | 完整runs | 本地完整参与者/runs | 缺失参与者 | 分类 |
|---|---:|---:|---:|---:|---|
| MHNO | 103 | 348 | 5 / 5 | 98 | MULTI_COHORT_NOT_RECONSTRUCTABLE |
| MHO | 19 | 19 | 0 / 0 | 19 | MULTI_COHORT_NOT_RECONSTRUCTABLE |
| MUNO | 377 | 1250 | 25 / 25 | 352 | MULTI_COHORT_NOT_RECONSTRUCTABLE |
| MUO | 432 | 1861 | 20 / 20 | 412 | MULTI_COHORT_NOT_RECONSTRUCTABLE |

四个phenotype互斥，完整身份已恢复，共931人、3478 runs；本地缺少大部分来源，四组都不能严格重建。独立清单记录E1出版manifest、E2作者列名和ENA映射的E5联合链，并明确 `status=partial`。公开仓库未包含最终 `physeq_adj.rds` 和部分中间metadata对象，这一限制已保留，不影响出版manifest本身的逐样本核验。

## 6. Li2024

论文：[Advanced Science](https://doi.org/10.1002/advs.202308313)。Methods将网络定义为来自discovery cohorts的60个signature species的SparCC关联；“Li_IBS”是研究网络名称，不足以证明仅含IBS病例，不能据此丢掉discovery对照。

Table S1明确区分14个discovery assay cohorts（13个16S、1个shotgun），对应13个唯一项目：

`PRJCA002555, PRJEB11252, PRJEB42304, PRJNA268708, PRJNA373876, PRJNA391149, PRJNA46339, PRJNA475187, PRJNA541075, PRJNA541572, PRJNA544721, PRJNA566284, PRJNA589202`。

其中PRJEB42304在Table S1两行出现：16S为63HC+78IBS=141，WGS为59HC+80IBS=139。不能把它当整个network输入，也不能把两个assay计数直接相加成独立参与者数。本地139个PRJEB42304 WGS runs完整存在，但最终network逐样本纳入关系仍未恢复，不能用数量139相等认证exact。

Validation项目包括PRJEB11419（AGP，16S与WGS）、PRJEB37924（baseline），以及2024年追加功能验证的PRJEB19857（baseline）和PRJNA812699（metatranscriptome）。本轮新增识别PRJEB11419本地467列，并验证原始profile值与矩阵列一致；其角色是validation，不能加入primary discovery network来增加“本地覆盖”。Treatment部分还包含PRJNA279031、PRJNA391915，以及部分discovery项目的治疗样本；项目角色重叠不等于治疗后样本用于discovery。[完整角色审计](01_sample_rescue/evidence/Li_discovery_validation_treatment_audit.tsv)

检查了supplement ZIP中的Excel、PDF、DOCX及作者[仓库](https://github.com/junhuili/IBS_co-occurrence_network)，main页面commit为 `7c1233ae7a57ee542c262b1be3849051c3cc507e`；仓库只有README与species-pair相关表压缩包，没有逐样本manifest或构建输入脚本。补充表为队列计数、差异/功能统计及模型结果，未提供最终网络sample→run清单。论文引用Figshare记录12525242是Kraken2参考数据库，不是样本表，未下载其大型参考payload。

已实际下载其公开项目RunInfo/可获得的SRA属性。PRJCA002555在NCBI/ENA不返回完整可用runs后，继续访问[NGDC官方项目](https://ngdc.cncb.ac.cn/bioproject/browse/PRJCA002555)，取得154个SAMC BioSample及真实archive链接CRA002557、CRA002559；archive页面返回浏览器安全检查，故未声称获得GSA完整runs。早期探索的CRA002235没有官方关联证据，明确排除，不列为发现的来源。

Li主表保留1838个公开run候选和154个无run的GSA样本候选，均标U；它们不是最终network exact清单。对公共SRA一个experiment包包含多BioSample的情况，解析器按run内部唯一Pool/Member精确选择，而非取包中第一个BioSample；NCBI/ENA冲突检查为0。

结果：最终network N、QC后选择和逐样本角色未唯一确定，且本地只存在13个discovery项目中的1个，**MULTI_COHORT_NOT_RECONSTRUCTABLE**。139与467是本地项目列数，不是已认证的network贡献数。

## 7. 可复核性、资源限制和完整性

本轮主动访问出版商、supplement、作者仓库、Zenodo、Dryad、Figshare、NCBI/ENA、NGDC和EGA；按标题/DOI/作者/项目以及metadata/sample/repository等词进行了检索。搜索记录保存在evidence目录；没有把搜索摘要用作样本分配证据。针对GitHub不只读首页：[文件清单](01_sample_rescue/evidence/repository_file_inventory.tsv)记录branch/commit，[关键词检索](01_sample_rescue/evidence/repository_keyword_search.tsv)记录具体脚本行及路径。Guo未发现论文对应的独立sample-manifest仓库；通用FlashWeave/ANCOM工具仓库不是该研究的样本证据。

出版商部分HTML、GitHub API、Dryad下载、GSA run页出现403/401或HTTP200挑战页面。原错误及替代来源都已记录；挑战页面显式标为 `blocked_challenge_response`。Nature旧完整全文及新官方补充文件、GitHub raw文件/页面及版本zip等可用来源继续推进了工作。公开资源不存在与本次访问失败严格区分。

共检查29个BioProject及其明确别名，规范化73,240个公开runs；大型项目用完整ENA RunInfo作为NCBI批量下载替代，PRJCA002555的完整run缺口单独披露。[最终metadata可用性](01_sample_rescue/evidence/project_metadata_final_availability.tsv)。EGA三个数据集的公开元数据另存JSON；未尝试下载受控测序文件或任何大型FASTQ。

本地只枚举指定abundance根目录的1577个矩阵文件名，并对已识别项目/别名作精确查找，没有全盘扫描。数字矩阵列不能当run ID：沿用Stage5B-0已验证的7项目2376列，并对本轮新增PRJEB11419、PRJEB41311共508列检查原始ordered run list及每个保留taxon的原profile数值，共2884列、9个矩阵。[可追踪列号索引](01_sample_rescue/evidence/verified_local_run_columns.tsv)

QC覆盖重复run、case/control参与者冲突、run→BioSample唯一性、BioSample歧义、公开/恢复样本计数、本地完整输入、重复列、缺失样本、额外本地runs、互斥组及证据完整性。WARN表示真实数据缺口或未解决拟合来源，未将其隐藏成PASS。Chen/Li部分本地run的最终纳入仍未解释，保留UNRESOLVED；Zhou/Lacruz非目标样本明确排除；Guo没有额外本地run。

已对{len(integrity)}个登记原始文件核对前后SHA256，全部不变；其中Stage5B-0最初149个冻结文件有执行前快照。[输入完整性](01_sample_rescue/qc/frozen_input_integrity.tsv)。所有下载、脚本、审计产物位于01_sample_rescue。仅按任务第十九节明确指定，在Stage5B根目录额外写入本报告，报告副本也保留在01_sample_rescue/reports。未执行sudo或对外发送消息。

审计脚本：`build_rescue.py`、`finalize_evidence.py`、`write_report.py`、`validate_delivery.py`。下载/提取和执行事件记录在[execution_log.tsv](01_sample_rescue/logs/execution_log.tsv)；最终产物SHA256另存qc/delivery_artifact_manifest.tsv。重新运行审计脚本不会运行网络推断。

## 8. 下一步推荐表

以下Local abundance match count仍以**完整生物学输入**为单位。Agg行明确是父级集合，不是完整子抽样记录。

| Network | Classification | Exact sample count | Local abundance match count | Can rerun FastSpar? | Reason |
|---|---|---:|---:|---|---|
{chr(10).join(recommended)}

**NEXT_STAGE_READY_NETWORKS（本轮新增，仅样本资格）**：

- Guo_Control — 91 exact first-visit samples/runs。
- Guo_MECFS — 106 exact first-visit samples/runs。
- Guo_Common — 197 exact first-visit samples/runs。

这三个对象均为 **PASS_EXACT**，允许下一阶段使用对应完整清单准备输入；本轮不启动拟合。新中间对象的物种映射和条件/边际关联对比协议仍应在下一阶段明确，不能默认为已解决。

**ALREADY_READY_FROM_STAGE5B_0（只读沿用，不计本轮新增）**：Wallen_PD=490、Wallen_NHC=234，原Class B及既有exact清单不变。Wallen原方法已经是FastSpar，应作为参数/流程复现对照，不能把它描述为不同推断算法验证。本轮不重复救援，也不修改原正式11输入的A/B/C/D分级表。
'''
    # Explicit report location requested in section19 is the sole output outside R.
    report=STAGE/'STAGE5B_05_REPORT.md';report.write_text(text,encoding='utf-8')
    (R/'reports/STAGE5B_05_REPORT.md').write_text(text,encoding='utf-8')
    event('write_report.py: render report and sample-ready input manifest',R,report)
    print('REPORT_WRITTEN',report)

if __name__=='__main__':main()
