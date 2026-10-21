# get_link_info5.py规范性和潜在科学风险

## 结论

**中文日志/历史注释本身是STYLE_ONLY；静默异常、支持向量错位、未实现tool模式及舍入策略涉及真实诊断或科学输出风险，不能归为单纯风格。** 本轮没有修改该脚本，也没有据合成反例宣布论文正式统计错误。

源完整SHA见FULL_SCRIPT_REFERENCE.tsv；逐问题行证据、分类和建议见GET_LINK_INFO5_RISK_REGISTER.tsv。`get_info`函数AST与已归档独立复现源完全相同，且现存生产共享脚本也具有相同逐值5位舍入和再平均逻辑。路径移植没有引入这一精度策略；没有全部历史生产版本快照，不能证明所有批次均使用同一函数。

## 实际微型诊断

四次运行均只读取合成2–3项目文件，单次约0.5秒，退出码均0，输出隔离在本审计目录，没有科研队列或正式数据写入：

1. 有效2项目数据但belonging声明999：脚本忽略声明并输出proj_count2，没有错误提示。
2. 第一个项目的CSV缺失：输出保留2个project ID和2个sample count，但corr/P列表仅1项。不能按位置给其分数指派原项目；均值也改为可找到的子集。
3. `--tool ccl_freq`：帮助声明该模式，实际没有赋值info_file分支，异常被吞，输出均值NA和空corr/P，仍退出0。
4. 3项目corr=0.3000049,0.3000049,0.3000059：先舍入再平均得到0.30000；原值均值再舍入是0.30001。1e-6的个体P在当前5位舍入后成为0。没有修改这套历史算法去提高精度。

证据：GET_LINK_INFO5_SYNTHETIC_DIAGNOSTICS.json、GET_LINK_INFO5_ROUNDING_DIAGNOSTIC.json及本地隔离fixture/log目录（未随GitHub打包）。初始2项目舍入案例只验证逐值舍入/向量，不把它冒称均值差异证明；第4例才直接证明均值差异。

## 支持数、名称、样本数和实际tool

get_dic忽略输入project_count，保留ID列表、不检查唯一性；get_info对每个支持ID无条件追加sample count，对成功查到的pair才追加corr/P。发生错误时没有占位或关联键保护。即使所有文件存在，重复项目ID/重复边及CSV解析异常仍需预先排除。sample_count只由header推导、不是读取/验证全部矩阵；失败写NA。NaN/Inf和corr/P值域未检查。

`--tool fastspar`是唯一有实际文件分支的模式。fastCCLasso应该使用现有专用reader；不以相似文件名猜测兼容性。`-s/--source`虽然解析，但当前main不调用get_environment，它不会校验环境映射。main依赖 `_core_links_filtered.belonging` basename并按下划线解析环境。PIPELINE已只修正文档的命名命令，不修改程序。

## 对论文的含义与修正边界

这是条件性SCIENTIFIC_OUTPUT_RISK：如果正式支持表存在缺失/错位/重复输入，可能改变environment聚合系数、证据归属及下游source-exclusion统计。该任务没有遍历正式全表，因此实际受影响记录数和论文数值变化均UNKNOWN，不能写成0，也不能说已经证实错误。已验收论文结果/归档本轮保持原状。

当前建议先把该入口标为条件历史流程，运行前检查所有输入及完整对齐；将来若授权改代码，应增加明确异常、精确tool choices、支持ID/edge/有限值/重复及environment检查，配套这些小fixture。舍入策略的改变需单独科学影响评估，不能作为“清理代码”偷偷改变。

## 其他77脚本范围

复用此前93源码审计、Phase1调用/移除审计、Stage2实际运行证据。重点保留EPIcore兼容副本、同内容workflow依赖和DKI R定义，不重做无界筛选。其他现有风险：原始prevalence实现round(M*f)与数学ceil的边界差异；EPI原core有效值兼容行为依赖入口检查；DKI/R批处理exit0可能含项目FAILED，必须看status；部分外部main忽略CLI且需要14held辅助源；旧Stage5B main不等于终稿S9。已在模块README/索引记录，不修改公式或结果。
