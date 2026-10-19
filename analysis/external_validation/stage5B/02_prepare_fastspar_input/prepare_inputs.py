"""Stage 5B-1A: exact column selection and unmodified MACK-DB count conversion only.

All mutations are restricted to this stage's directory. No inference executable is invoked.
"""

# v0.1: repository-local IO configuration; scientific functions below are unchanged.
import sys as _path_sys
from pathlib import Path as _RepoPath
_repo_root = next(p for p in _RepoPath(__file__).resolve().parents if (p / "lib/mackdb_paths.py").is_file())
_path_sys.path.insert(0, str(_repo_root / "lib"))
from mackdb_paths import resolve_path
from pathlib import Path
from decimal import Decimal
import csv, hashlib, json, datetime, subprocess, shlex, os, sys, collections

ROOT=Path(__file__).resolve().parent
STAGE=ROOT.parent
RESCUE=STAGE/'01_sample_rescue'
AB=Path(str(resolve_path('species_input', 'abundance')))
CONVERTER=Path(str(resolve_path('count_script')))
PIPELINE=Path(str(resolve_path('fastspar_script')))
RSCRIPT=Path('/usr/bin/Rscript')
NETWORKS={'Guo_Control':91,'Guo_MECFS':106,'Guo_Common':197,'Wallen_PD':490,'Wallen_NHC':234}
SUB=ROOT/'abundance_subset';COUNT=ROOT/'count_matrix';QC=ROOT/'qc';LOG=ROOT/'logs'
SCALE=Decimal(1000000)

def timestamp():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()
def readtsv(p):
    with Path(p).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def write(p,rs,fields=None):
    p=Path(p);assert ROOT in p.resolve().parents
    rs=list(rs)
    with p.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields or list(rs[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rs)
def event(command,inp,out,code=0,start=None):
    p=LOG/'execution_log.tsv';new=not p.exists()
    with p.open('a',encoding='utf-8',newline='') as f:
        w=csv.writer(f,delimiter='\t',lineterminator='\n')
        if new:w.writerow(['command','timestamp','input','output','exit_code','completed_at'])
        w.writerow([command,start or timestamp(),str(inp),str(out),code,timestamp()])
def readmatrix(p):
    with p.open() as f:
        header=f.readline().rstrip('\r\n').split('\t')
        data=[line.rstrip('\r\n').split('\t') for line in f]
    assert header[0]=='OTU_id',str(p)
    assert len(set(header[1:]))==len(header)-1,str(p)
    assert data and all(len(r)==len(header) for r in data),str(p)
    assert len({r[0] for r in data})==len(data),str(p)
    assert all(Decimal(x).is_finite() and Decimal(x)>=0 for r in data for x in r[1:]),str(p)
    return header,data
def source_list(network):
    requested=RESCUE/'sample_lists'/(network+'_exact_samples.tsv')
    if requested.is_file():return requested,'Stage5B-0.5 exact list'
    assert network.startswith('Wallen_'),'Missing current rescue list: '+network
    fallback=STAGE/'00_audit'/(network+'_exact_sample_list.tsv')
    assert fallback.is_file(),'Missing frozen Wallen exact list'
    return fallback,'Stage5B-0 frozen exact list; Stage5B-0.5 deliberately reused Wallen without copying/rescuing'

def template():
    original=PIPELINE.read_text()
    for token in ['--iterations 50','--exclude_iterations 20','--threshold 0.2','--number 10000','-i 5','--permutations 10000']:
        assert token in original,'Pipeline parameters changed: '+token
    txt='''Stage 5B-1A FastSpar pipeline audit — READ ONLY; NOT EXECUTED

Source script: /home/wangmengyang/project/script/run_fastspar.sh
Source SHA256: {digest}

INPUT: tab-delimited OTU table; first line OTU_id followed by sample labels;
       species/taxa in rows, selected samples in columns, nonnegative integer values.
       Unmodified count_sparcc.R resets both sample and row labels to 0-based indices.
       Use qc/<network>_sample_column_map.tsv and qc/<network>_taxa_row_map.tsv to recover identity.

OUTPUTS: <PREFIX>_correlation.tsv, <PREFIX>_covariance.tsv, <PREFIX>_pvalues.tsv.
         Correlation/covariance/pvalue objects are taxa-by-taxa, not sample tables.
         No inference output is produced during Stage 5B-1A.

MAIN FIT: iterations=50; exclude_iterations=20; threshold=0.2.
BOOTSTRAP: 10000 resampled OTU tables; each bootstrap fit passes only -i 5.
           Bootstrap exclude_iterations, threshold and threads are NOT explicitly set
           by the source script; their effective defaults require the installed tool's
           documentation in the future. They must not be silently reported as main-fit values.
PVALUES: --permutations 10000.
RESOURCES: source script sets T=(processor lines in /proc/cpuinfo)-3;
           GNU parallel uses its default job count for the 10000 bootstrap fits.
RANDOM SEED: not explicitly specified by the original script.
BINARIES: /home/changqing/miniconda3/bin/fastspar, fastspar_bootstrap, fastspar_pvalues.
          Exact paths were checked for existence only; no FastSpar binary, including --help,
          was executed during this stage.
CLEANUP: original script ends with rm -rf of both bootstrap directories.
         The expanded reference below leaves this cleanup commented for inspection.

Simplest future invocation (NOT EXECUTED; replace only input/output prefix):

bash /home/wangmengyang/project/script/run_fastspar.sh \\
  /ABSOLUTE/PATH/NETWORK_count.matrix \\
  /ABSOLUTE/PATH/FUTURE_OUTPUT_DIR/NETWORK

The parent of FUTURE_OUTPUT_DIR/NETWORK must exist. It is not created in this stage.

Expanded original command structure (reference text, NOT an executed shell script):

PROFILE=/ABSOLUTE/PATH/NETWORK_count.matrix
PREFIX=/ABSOLUTE/PATH/FUTURE_OUTPUT_DIR/NETWORK
MY_PATH=/home/changqing/miniconda3/bin
cpu=$(cat /proc/cpuinfo | grep processor | wc -l)
T=$((cpu-3))
CORR=${{PREFIX}}_correlation.tsv
COVAR=${{PREFIX}}_covariance.tsv
BOOSTSTARP_COUNTS=${{PREFIX}}_bootstrap_counts
BOOSTSTARP_CORR=${{PREFIX}}_bootstrap_correlation
PVALUE=${{PREFIX}}_pvalues.tsv

"$MY_PATH/fastspar" --iterations 50 --exclude_iterations 20 --threshold 0.2 \\
  --otu_table "$PROFILE" --correlation "$CORR" --covariance "$COVAR" --threads "$T"
mkdir "$BOOSTSTARP_COUNTS"
"$MY_PATH/fastspar_bootstrap" --otu_table "$PROFILE" --number 10000 \\
  --prefix "$BOOSTSTARP_COUNTS/fake_data" --threads "$T"
mkdir "$BOOSTSTARP_CORR"
parallel "$MY_PATH/fastspar" --otu_table {{}} --correlation "$BOOSTSTARP_CORR/cor_{{/}}" \\
  --covariance "$BOOSTSTARP_CORR/cov_{{/}}" -i 5 ::: "$BOOSTSTARP_COUNTS"/*
"$MY_PATH/fastspar_pvalues" --otu_table "$PROFILE" --correlation "$CORR" \\
  --prefix "$BOOSTSTARP_CORR/cor_fake_data_" --permutations 10000 --outfile "$PVALUE" --threads "$T"
# Original automatic cleanup, shown only for audit:
# rm -rf "$BOOSTSTARP_COUNTS" "$BOOSTSTARP_CORR"

Concrete prepared inputs (replace FUTURE_OUTPUT_DIR before a future manual run):

{commands}

No seed, parameter, taxonomy, filtering, sign or network comparison decision is made here.
'''.format(digest=sha(PIPELINE),commands='\n'.join('bash '+str(PIPELINE)+' '+str(COUNT/(n+'_count.matrix'))+' /ABSOLUTE/PATH/FUTURE_OUTPUT_DIR/'+n for n in NETWORKS))
    (QC/'stage5B1A_fastspar_command_template.txt').write_text(txt)
    write(QC/'stage5B1A_pipeline_binary_path_audit.tsv',[dict(path=str(Path('/home/changqing/miniconda3/bin')/n),exists=(Path('/home/changqing/miniconda3/bin')/n).is_file(),executed='NO') for n in ['fastspar','fastspar_bootstrap','fastspar_pvalues']])
    event('Read run_fastspar.sh; write command template only (no execution)',PIPELINE,QC/'stage5B1A_fastspar_command_template.txt')

def main():
    assert str(ROOT)==str(resolve_path('stage5b', '02_prepare_fastspar_input'))
    for d in [SUB,COUNT,QC,LOG,LOG/'r_tmp']:d.mkdir(parents=True,exist_ok=True)
    event(shlex.join([sys.executable,'-s','-B',str(Path(__file__).resolve())]),RESCUE,ROOT)
    summarypath=RESCUE/'stage5B_05_sample_rescue_summary.tsv'
    summary={r['network']:r for r in readtsv(summarypath)}
    verifiedpath=RESCUE/'evidence/verified_local_run_columns.tsv'
    index=readtsv(verifiedpath);runs={r['run_accession']:r for r in index};assert len(runs)==len(index)
    frozenpath=RESCUE/'qc/input_snapshot.json';frozen=json.loads(frozenpath.read_text())
    deliveredpath=RESCUE/'qc/delivery_artifact_manifest.tsv';delivered={r['path']:r['sha256'] for r in readtsv(deliveredpath)}
    sourcefiles={summarypath,verifiedpath,frozenpath,deliveredpath,CONVERTER,PIPELINE,STAGE/'00_audit/stage5B_published_component_eligibility.tsv'}
    requested={};sources={};selection={};selection_qc=[];matched=[];source_audit=[];rowmaps=[];input_checks=[]
    # Preflight every group before writing matrices or performing any conversion.
    for network,expected in NETWORKS.items():
        lp,resolution=source_list(network);sourcefiles.add(lp);rs=readtsv(lp)
        assert len(rs)==expected and len({r['run_accession'] for r in rs})==expected,(network,'sample list mismatch')
        if network.startswith('Guo_'):
            assert summary[network]['classification']=='PASS_EXACT' and summary[network]['next_stage_ready']=='YES'
            assert all(r['exact_network_assignment']=='YES' and r['found_in_local_abundance']=='YES' and r['status']=='exact' for r in rs)
        else:assert all(r['certified_for_exact_biological_samples']=='Yes' for r in rs)
        normalized=[]
        for ordinal,r in enumerate(rs):
            run=r['run_accession'];v=runs.get(run);assert v,(network,run,'missing verified run')
            project=r.get('bioproject',r.get('source_bioproject'))
            matrix=AB/(project+'_abundance.matrix')
            assert matrix==Path(r.get('local_abundance_file',r.get('abundance_file')))==Path(v['abundance_file'])
            assert v['source_bioproject']==project and v['column_index']==r.get('column_index',r.get('mackdb_column_index'))
            assert v['profile_values_match_every_retained_taxon']=='Yes'
            sourcefiles.add(matrix);sourcefiles.add(Path(v['run_list_file']))
            if matrix not in sources:sources[matrix]=readmatrix(matrix)
            head,data=sources[matrix];ci=int(v['column_position_1based'])
            assert head[ci]==v['column_index'] and len(head)-1>0
            originalruns=Path(v['run_list_file']).read_text().splitlines()
            assert len(originalruns)==len(head)-1 and len(originalruns)==len(set(originalruns)) and originalruns[ci-1]==run
            row=dict(network=network,output_sample_position_1based=ordinal+1,count_sample_label=str(ordinal),abundance_sample_label=run,
                run_accession=run,biosample=r.get('biosample',r.get('biosample_accession')),participant_id=r.get('participant_id',r.get('publication_subject_id')),
                publication_sample_id=r.get('publication_sample_id',r.get('publication_subject_id')),group=r.get('group',r.get('publication_group')),
                bioproject=project,input_matrix=str(matrix),input_column_label=head[ci],input_sample_position_1based=ci,
                exact_sample_list=str(lp),sample_list_row_1based=ordinal+2,verified_column_index_file=str(verifiedpath),
                evidence_level=r.get('evidence_level','frozen_Stage5B-0_certified'),identity_rule=v['identity_rule'])
            normalized.append(row)
        assert len({r['input_matrix'] for r in normalized})==1,'Current five networks must each have one project'
        assert len({r['participant_id'] for r in normalized})==expected and len({r['biosample'] for r in normalized})==expected
        requested[network]=rs;selection[network]=normalized;matched.extend(normalized)
        p=Path(normalized[0]['input_matrix']);head,data=sources[p]
        selection_qc.append(dict(network=network,input_matrix=str(p),input_total_samples=len(head)-1,requested_samples=expected,matched_samples=len(normalized),missing_samples=0,
            output_samples=len(normalized),taxa_count=len(data),duplicate_sample_columns=len(normalized)-len({r['run_accession'] for r in normalized}),sample_match_status='PASS',
            species_order_status='pending',abundance_value_status='pending',exact_sample_list=str(lp),sample_list_resolution=resolution,output_sample_order='exact sample list row order',abundance_header_policy='verified exact run accessions',
            output_matrix=str(SUB/(network+'_abundance.matrix'))))
        source_audit.append(dict(network=network,requested_list=str(RESCUE/'sample_lists'/(network+'_exact_samples.tsv')),resolved_list=str(lp),resolution=resolution,bioproject=normalized[0]['bioproject'],source_matrix=str(p)))
    write(QC/'stage5B1A_resolved_sample_list_sources.tsv',source_audit)
    write(QC/'stage5B1A_sample_extraction_qc.tsv',selection_qc)
    # Protect and compare all source files to the previous frozen artifact hashes where available.
    baseline={str(p):dict(sha256=sha(p),size=p.stat().st_size) for p in sorted(sourcefiles)}
    for p,x in baseline.items():
        prior=delivered.get(p,frozen.get(p,{}).get('sha256'))
        if prior:assert x['sha256']==prior,'Input changed since frozen rescue: '+p
        input_checks.append(dict(path=p,current_sha256=x['sha256'],prior_sha256=prior or '',prior_freeze_check='PASS' if prior else 'CURRENT_STAGE_BASELINE'))
    snapshot=QC/'stage5B1A_input_snapshot.json'
    if snapshot.exists():
        initial=json.loads(snapshot.read_text())
        assert all(p in baseline and baseline[p]==x for p,x in initial.items()),'Input changed since initial Stage5B-1A preflight'
    else:snapshot.write_text(json.dumps(baseline,indent=2)+'\n')
    write(QC/'stage5B1A_prior_frozen_input_checks.tsv',input_checks)
    write(QC/'stage5B1A_resolved_exact_samples.tsv',matched)
    # Preserve every original numerical string and every species ID, in original row order.
    for q in selection_qc:
        network=q['network'];rs=selection[network];matrix=Path(q['input_matrix']);head,data=sources[matrix]
        chosen=[int(r['input_sample_position_1based']) for r in rs];out=SUB/(network+'_abundance.matrix')
        with out.open('w',encoding='utf-8',newline='') as f:
            f.write('\t'.join([head[0]]+[r['run_accession'] for r in rs])+'\n')
            for row in data:f.write('\t'.join([row[0]]+[row[j] for j in chosen])+'\n')
        oh,od=readmatrix(out)
        assert oh[1:]==[r['run_accession'] for r in rs]
        assert [r[0] for r in od]==[r[0] for r in data]
        assert all(row[1:]==[orig[j] for j in chosen] for row,orig in zip(od,data))
        q.update(species_order_status='PASS',abundance_value_status='PASS_BYTE_IDENTICAL_NUMERIC_TOKENS')
        write(QC/(network+'_sample_column_map.tsv'),rs)
        taxamap=[dict(network=network,row_position_1based=i+1,original_species_id=r[0],abundance_species_id=r[0],count_species_id=str(i)) for i,r in enumerate(data)]
        write(QC/(network+'_taxa_row_map.tsv'),taxamap);rowmaps.extend(taxamap)
        event('Exact accession/verified-column extraction; all species and abundance strings retained',matrix,out)
    write(QC/'stage5B1A_sample_extraction_qc.tsv',selection_qc)
    assert all(q['sample_match_status']=='PASS' and q['matched_samples']==q['requested_samples'] and q['missing_samples']==0 for q in selection_qc)
    conversion=[]
    env=dict(os.environ,TMPDIR=str(LOG/'r_tmp'),R_LIBS_USER=str(LOG/'R_library'))
    for network in NETWORKS:
        inp=SUB/(network+'_abundance.matrix');out=COUNT/(network+'_count.matrix')
        cmd=[str(RSCRIPT),str(CONVERTER),'-i',str(inp),'-n','1000000','-o',str(out)]
        start=timestamp()
        with (LOG/(network+'_count_conversion.stdout.log')).open('w') as stdout,(LOG/(network+'_count_conversion.stderr.log')).open('w') as stderr:
            proc=subprocess.run(cmd,cwd=ROOT,env=env,stdout=stdout,stderr=stderr)
        event('R_LIBS_USER='+shlex.quote(env['R_LIBS_USER'])+' TMPDIR='+shlex.quote(env['TMPDIR'])+' '+shlex.join(cmd),inp,out,proc.returncode,start)
        assert proc.returncode==0,'Count conversion failed: '+network
        ah,ad=readmatrix(inp);ch,cd=readmatrix(out)
        mismatch=noninteger=negative=nonfinite=0;maximum=Decimal(0)
        assert len(cd)==len(ad) and len(ch)==len(ah)
        for a,c in zip(ad,cd):
            for x,y in zip(a[1:],c[1:]):
                val=Decimal(y);target=Decimal(x)*SCALE
                if not val.is_finite():nonfinite+=1;continue
                difference=abs(val-target);maximum=max(maximum,difference)
                mismatch+=val!=target;noninteger+=val!=val.to_integral_value();negative+=val<0
        taxa_ok=[r[0] for r in cd]==[str(i) for i in range(len(ad))]
        col_ok=ch[1:]==[str(i) for i in range(len(ah)-1)]
        assert taxa_ok and col_ok
        # The unmodified script resets row IDs. In these source matrices they are already sequential.
        species_preserved=[r[0] for r in ad]==[r[0] for r in cd]
        totals=[sum(Decimal(r[j]) for r in cd) for j in range(1,len(ch))]
        zerotaxa=sum(all(Decimal(v)==0 for v in r[1:]) for r in cd)
        result=dict(network=network,abundance_input=str(inp),count_output=str(out),conversion_script=str(CONVERTER),scale=1000000,exit_code=proc.returncode,
            abundance_taxa_count=len(ad),count_taxa_count=len(cd),abundance_sample_count=len(ah)-1,count_sample_count=len(ch)-1,
            taxa_count_status='PASS',sample_count_status='PASS',OTU_id_header_status='PASS' if ch[0]=='OTU_id' else 'FAIL',
            noninteger_values=noninteger,integer_status='PASS' if noninteger==0 else 'FAIL',scaling_mismatches=mismatch,max_absolute_scaling_difference=str(maximum),
            abundance_times_1000000_status='PASS' if mismatch==0 else 'FAIL',NA_or_nonfinite_values=nonfinite,negative_values=negative,
            species_id_order_status='PASS' if species_preserved else 'FAIL',sample_label_reset_status='PASS',all_zero_taxa_retained=zerotaxa,all_zero_sample_columns=sum(v==0 for v in totals),
            min_sample_total=str(min(totals)),max_sample_total=str(max(totals)),numeric_cells_checked=(len(ah)-1)*len(ad),
            status='PASS' if not(noninteger or mismatch or negative or nonfinite) and species_preserved else 'FAIL')
        conversion.append(result);write(QC/'stage5B1A_count_conversion_qc.tsv',conversion)
        assert result['status']=='PASS',(network,result)
        print(network,'samples',len(ch)-1,'taxa',len(cd),'extraction PASS; original R conversion PASS',flush=True)
    gs={n:{r['run_accession'] for r in selection[n]} for n in NETWORKS}
    partition=[dict(check='Guo_Control_MECFS_sample_exclusivity',status='PASS' if not gs['Guo_Control']&gs['Guo_MECFS'] else 'FAIL',details='91 and106 disjoint exact run identities'),
        dict(check='Guo_Common_is_exact_case_control_union',status='PASS' if gs['Guo_Control']|gs['Guo_MECFS']==gs['Guo_Common'] else 'FAIL',details='197 distinct exact identities, verified by run strings'),
        dict(check='Wallen_PD_NHC_sample_exclusivity',status='PASS' if not gs['Wallen_PD']&gs['Wallen_NHC'] else 'FAIL',details='490 and234 disjoint exact run identities')]
    excluded=[]
    for project,nets in [('PRJNA751448',['Guo_Control','Guo_MECFS']),('PRJNA834801',['Wallen_PD','Wallen_NHC'])]:
        included=set().union(*(gs[n] for n in nets))
        for r in index:
            if r['source_bioproject']==project and r['run_accession'] not in included:
                excluded.append(dict(project=project,run_accession=r['run_accession'],column_index=r['column_index'],reason='Not in either frozen exact subgroup list; no replacement or supplementation'))
    write(QC/'stage5B1A_group_membership_qc.tsv',partition)
    write(QC/'stage5B1A_source_runs_not_selected.tsv',excluded,['project','run_accession','column_index','reason'])
    assert all(r['status']=='PASS' for r in partition)
    template()
    integrity=[dict(path=p,before_sha256=x['sha256'],after_sha256=sha(p),status='PASS' if sha(p)==x['sha256'] else 'FAIL') for p,x in baseline.items()]
    write(QC/'stage5B1A_input_integrity.tsv',integrity)
    assert all(r['status']=='PASS' for r in integrity),'Original input changed'
    event('Before/after source SHA256 verification',str(len(integrity))+' frozen/original files',QC/'stage5B1A_input_integrity.tsv')
    (LOG/'preparation_state.json').write_text(json.dumps(dict(state='MATRICES_GENERATED_PENDING_INDEPENDENT_READBACK',timestamp=timestamp(),source_files=len(integrity),networks=NETWORKS,fastspar_executed=False),indent=2)+'\n')

if __name__=='__main__':
    try:main()
    except Exception as ex:
        if LOG.exists():
            event('prepare_inputs.py ERROR',ROOT,str(ex),1)
            (LOG/'preparation_error.json').write_text(json.dumps(dict(error=repr(ex),timestamp=timestamp()),indent=2)+'\n')
        raise
