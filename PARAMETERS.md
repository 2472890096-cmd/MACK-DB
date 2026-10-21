# Scientific parameters and evidence levels

`FROZEN_SOURCE` means the reviewed actual source/retained caller. `ACCEPTED_ARCHIVE` means prior accepted provenance copied in the sole input directory. `EXECUTED_SYNTHETIC` means this task actually exercised the implementation on tiny synthetic input. None means full-cohort re-execution.

| Stage | Required scientific configuration | Evidence and caveat |
|---|---|---|
| Sketch / Query | Fixed realERR14036375: MetaKSSD2.21, pinnedL3K11/ID1481063132/22mer/nominal1/4096, `dist -A -p1 -m2`, `composite -r markerDB -q sampleSketch` | ACTUAL_SINGLE_RUN_20261021; exact key/count and Query matching; full original historical argv remains limited |
| Profile | `sharekmer=12`, keep X3>=12; X5<=5 uses `[0,.02,.08,.69,1.57,2.49,...]` lookup; X5>5 subtract3; percentage round4 | unchanged099 + EXECUTED_SYNTHETIC; 12 is shared-marker cutoff, not k-mer length |
| Species matrix | f starts .1, +.1 until <=2000 taxa or f>=1; worker cutoff `round(n*f)` | 053+101, algorithm unchanged; dictionary/ID ordering across historical Python versions is not newly validated |
| Counts | input abundance percentages ×1,000,000, no new rounding/read-count interpretation | unchanged092 |
| FastSpar | main iterations50/exclude20/threshold.2; 10000 resamples; bootstrap fit5; 10000 permutations | 100 source; thread count is now explicit (default1), scientific counts unchanged; bootstrap not run |
| Final network filter | abs(r)>=.3 AND P<=1e-4, unique upper triangle | EXECUTED_SYNTHETIC102; always pass `--corr_num 0.3`; retained generic R default .1 is historical, not final analysis setting |
| fastCCLasso | nboot100,kmax20; species seed20260804, genus seed20260806; original runner correlation filter .1 followed by final .3 filter | retained061/064 and060; official upstream core not bundled. Actual official boot_fastCCLasso uses the mean correlation over nboot+1 estimates (100 bootstraps plus original), t=rbar*sqrt((n-2)/(1-rbar^2)), u=pt(t,n-2), P=min(u,1-u), **without ×2**; not a sign-frequency bootstrap P. Source SHA and formula evidence: EXTERNAL_DEPENDENCIES.md. Documentation corrected in v0.2; no core execution/change |
| Environment support | >=2 distinct BioProjects (f=2), within environment | retained096/059; do not count repeated rows as projects |
| Zi-Pi modules | signed spin-glass `implementation="neg"`, per connected component; spins25,gamma1,gamma.minus1,seed12345; full network for scientific roles | unchanged074 function expressions; no optimization/plot execution in small suite |
| Zi/Pi roles | nonweighted within-module degree; sample SD; Pi=1−sum((kis/ki)^2); z<=2.5,p<=.62 are Peripheral, z<=2.5,p>.62 Connectors, z>2.5,p<=.62 Module hubs, both greater Network hubs | actual microeco2.0.0 EXECUTED_SYNTHETIC; valid undefined z/Pi yields NA rather than forced zero; NA roles never candidates |
| EPI v3 | mode=prevalence, inclusive .25–.75; remove focal and normalize remaining community; D1 cross-group mean BC; D2 BC of group mean compositions; Q unique off-diagonal-pair ECDF p=.1 with ties and first-jump fallback | exact frozen core SHA98d832…; actual original whole-program execution + independent scalar oracle <=1e-15 |
| EPI candidates | each metric strictly score > mean+2 sampleSD; >=2 finite scores needed; final >=2 of3 includes triple support; empty candidates can be successful | actual workers, scalar threshold and postprocess tests |
| DKI model | max-epochs1000, patience80, lr.005,batch32,relu,auto bottleneck,normalize-z true,time-mode endpoint,rk4 step.05,CPU,seed20260629 | retained104/105 override defaults; core defaults3000/200 must not replace production caller; no training in task |
| DKI R score | focal removed from qtrn and remaining predicted baseline normalized (`q_null`), BC(q_null,qtst)×(1−focal original abundance); median over valid cases | version MACKDB_DKI_R_20261012_v1; all090 function expressions unchanged; actual90+080 execution; old full-qtrn counterexample .32 vs accepted .16 |
| DKI candidates | full valid taxon background; sampleSD; K=min(N,max(3,ceil(.01N))); descending score then name; score >= mean+2SD among TopK | unchanged076 selection functions; zero-valued valid scores follow >= semantics, unlike EPI `>` |
| Figure3 | species/genus independent; seven mutually exclusive Broad combinations; Strict same project; >=2 includes All3; unique names across contexts | unchanged concordance_release.calc; synthetic sets + independent set operations tested |

Figure 3 accepted frozen main scope is1571, sensitivity1570/1568; PRJEB51353 is outside these scopes. The 22-value reference SHA is listed in the original supporting manifest. Accepted main values:

- Species A: 4575,3604,4786,879,1141,356,434 in the seven documented combination order.
- Genus B: 2108,1260,1826,330,621,161,273.
- C Broad species >=2/All3:2491/416; genus:982/243.
- D Strict species >=2/All3:970/25; genus:399/19.

These 22 numbers are preserved accepted references; this task did **not** recompute full-scope values because the 12 scientific inputs are not bundled. Assertions against the original mixed baseline and accepted all-four-R impact remain in the portable Figure3 release entry.

## v0.2 frozen-input validation

Figure 3 final main22 and both sensitivity22 sets were actually recomputed from frozen candidate records and independently cross-checked. Strict All3 remains25 species/19 genera. This extends numerical input validation only; all scientific parameters remain unchanged; original84 bytes are preserved in provenance, with the current retained set separately manifested. That full-data v0.2 evidence is held locally; this candidate includes scope/aggregate references and the unchanged statistics code, not row-level inputs.

The [real single-Run pins](examples/real_fastq/ERR14036375/reference_pins.json) add observed input/reference/binary hashes; they do not change the Profile core, EPI/DKI versions, association cutoffs or any manuscript result. The optional historicalSRAstream plan is not an executed alternative. All other parameter evidence retains its originally labelled scope.
