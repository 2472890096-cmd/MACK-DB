#!/usr/bin/env python3
"""Stage 3C: frozen input adapters + verbatim Chen/Wallen randomization kernel.

Run with python3 -s -B. Only Stage 3C outputs and its report are writable.
No taxonomy mapping, source exclusion, FastSpar, or downstream analyses.
"""

# v0.1: repository-local IO configuration; scientific functions below are unchanged.
import sys as _path_sys
from pathlib import Path as _RepoPath
_repo_root = next(p for p in _RepoPath(__file__).resolve().parents if (p / "lib/mackdb_paths.py").is_file())
_path_sys.path.insert(0, str(_repo_root / "lib"))
from mackdb_paths import resolve_path
import ast
import collections
import copy
import csv
import hashlib
import itertools
import json
import math
import multiprocessing as mp
import os
import random
import sys
import time
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd

OUT = Path(__file__).resolve().parents[1]
BASE = OUT.parent
STAGE3B = BASE / 'stage3b_source_exclusion'
LEGACY = resolve_path('shared_protocol', 'chen_wallen_null_models_final_v1.py')
HISTORY = OUT / 'history/protocol_stop_before_zero_sd_resolution'
AUTHORIZATION = Path(str(resolve_path('null_protocol_evidence')))
MASTER_SEED, N_ITER, SWAPS_PER_EDGE = 20260917, 1000, 50
WORKERS = min(16, len(os.sched_getaffinity(0)))
NETWORKS = ['Zhou_MS', 'Zhou_HHC', 'Guo_AggControl', 'Guo_AggMECFS', 'Lacruz_MHNO', 'Lacruz_MHO', 'Lacruz_MUNO', 'Lacruz_MUO', 'Li_IBS']
EXPECTED = dict(zip(NETWORKS, [(596,4818,1284), (585,4086,1201), (91,335,242), (91,312,210),
                             (235,649,556), (254,623,501), (192,486,436), (164,417,388), (54,599,375)]))
CHECKS = []
ROW_FIELDS = ['network', 'null_type', 'iteration', 'seed', 'random_edge_count', 'recovered_count', 'recovery_rate',
              'successful_swaps', 'swap_attempts', 'swap_acceptance_rate', 'degree_sequence_preserved',
              'per_node_degree_preserved', 'node_set_preserved', 'node_pool_count', 'sampled_endpoint_count',
              'all_endpoints_in_node_pool', 'no_self_loop', 'no_duplicate_edge', 'edge_count_preserved',
              'random_edge_sha256', 'reference_database_sha256', 'analysis_status']


def read(path):
    with Path(path).open(newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def write(path, rows, fields=None):
    rows = list(rows)
    path = Path(path)
    temp = path.with_suffix(path.suffix + '.tmp')
    with temp.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields or list(rows[0]), delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)
    temp.replace(path)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(8*1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def check(name, ok, detail=''):
    CHECKS.append(dict(check=name, status='PASS' if ok else 'FAIL', detail=str(detail)))
    if not ok:
        write(OUT / 'qc/03c_null_model_qc.tsv', CHECKS)
        raise RuntimeError((name, detail))


def seed_base(network, null_type):
    # Same disjoint study/null offset scheme as formal Chen/Wallen, expanded to 9 networks.
    return MASTER_SEED + 100000 * (2*NETWORKS.index(network) + (1 if null_type == 'A' else 2))


def edge_hash(keys):
    return hashlib.sha256(('\n'.join(map(str, sorted(keys))) + '\n').encode()).hexdigest()


def degree_vector(keys, n):
    d = [0]*n
    for k in keys:
        u,v = divmod(k,n)
        d[u] += 1
        d[v] += 1
    return d


def load_legacy():
    """Compile only pure functions. Original script/main/source exclusion never runs."""
    source = LEGACY.read_text()
    tree = ast.parse(source)
    fs = {n.name:n for n in tree.body if isinstance(n, ast.FunctionDef)}
    selected = [copy.deepcopy(fs[n]) for n in ['encode_edge', 'init_b_worker', 'null_b_worker', 'summarize']]
    worker = selected[2]
    check('legacy_worker_has_one_final_return', isinstance(worker.body[-1], ast.Return))
    # The original RNG, proposals, rejection logic, successful-swap loop and counting stay verbatim.
    original_body = ast.dump(ast.Module(body=worker.body[:-1], type_ignores=[]))
    old_return = worker.body[-1].value
    worker.body[-1].value = ast.Call(func=ast.Name(id='_audit_b_result', ctx=ast.Load()),
        args=[old_return, ast.Name(id='U', ctx=ast.Load()), ast.Name(id='V', ctx=ast.Load()), ast.Name(id='edge_set', ctx=ast.Load())], keywords=[])
    check('legacy_randomization_kernel_AST_unchanged',
          ast.dump(ast.Module(body=worker.body[:-1], type_ignores=[])) == original_body,
          'Only final return is wrapped in additional graph QC; no random draws added')
    module = ast.fix_missing_locations(ast.Module(body=selected, type_ignores=[]))
    exec(compile(module, str(LEGACY), 'exec'), globals())
    evidence = []
    for name in ['encode_edge','init_b_worker','null_b_worker','summarize']:
        fragment = ast.get_source_segment(source, fs[name])
        evidence.append(dict(function=name, source_file=str(LEGACY), source_sha256=sha(LEGACY),
            start_line=fs[name].lineno, end_line=fs[name].end_lineno,
            original_function_sha256=hashlib.sha256(fragment.encode()).hexdigest(),
            reuse='verbatim' if name != 'null_b_worker' else 'verbatim randomization core; return wrapped for per-iteration QC'))
    write(OUT / 'qc/03c_legacy_function_reuse.tsv', evidence)


def update_protocol():
    rows = read(HISTORY / 'qc/03c_chen_wallen_protocol_consistency.tsv')
    source_hash = sha(LEGACY)
    check('formal_legacy_script_frozen', all(r['evidence_sha256'] == source_hash for r in rows))
    for r in rows:
        if r['item'] == 'zero_null_SD_z_score':
            check('legacy_zero_SD_inf_branch_confirmed', 'if sd_count > 0 else float("inf")' in LEGACY.read_text())
            r.update(item='zero_sd_zscore_rule', stage3c_requirement='null SD=0 -> z=inf; null_sd_zero=Yes',
                     status='MATCH', consequence='Conflict resolved by explicit user instruction: ' + str(AUTHORIZATION))
        else:
            check('prior_protocol_check_' + r['item'], r['status'] == 'PASS')
            r['status'] = 'MATCH'
            if r['item'] == 'reproducibility_strategy':
                r['consequence'] = 'MASTER_SEED=20260917; disjoint 100000 offsets by network/null, plus iteration; all 18000 seeds unique'
    check('all_core_protocol_rows_MATCH', all(r['status'] == 'MATCH' for r in rows), len(rows))
    write(OUT / 'qc/03c_chen_wallen_protocol_consistency.tsv', rows)
    return rows


def prepare(integrity):
    hashes = {r['input_file']: r['sha256_before'] for r in integrity}
    previous_graphs = {r['network']:r for r in read(HISTORY / 'qc/03c_external_graph_audit.tsv')}
    summaries = {r['network']:r for r in read(STAGE3B / 'summary/03b_source_exclusion_summary.tsv')}
    references = {}
    for r in read(STAGE3B / 'source_excluded_networks/03b_database_manifest.tsv'):
        check('reference_frozen_' + r['database_group'], hashes[r['file']] == r['sha256'])
        for n in r['networks'].split(';'):
            check('unique_reference_' + n, n not in references)
            references[n] = r
    for f in ['03b_source_exclusion_qc.tsv', '03b_source_provenance_qc.tsv']:
        check('upstream_' + f, all(r['status']=='PASS' for r in read(STAGE3B / 'qc' / f)))
    pair_cache, studies, manifest = {}, {}, []
    for m in read(BASE / 'stage3a/qc/03a_input_manifest.tsv'):
        n,path = m['network'],m['input_file']
        check(n+'_pair_file_hash', hashes[path] == m['sha256'])
        if path not in pair_cache:
            pair_cache[path] = read(path)
        rows = [r for r in pair_cache[path] if r['network_id']==m['source_network_id']]
        pairs = sorted((r['species_a'],r['species_b']) for r in rows)
        check(n+'_exact_frozen_pair_keys', all(r['pair_key']==' || '.join(sorted((r['species_a'],r['species_b']))) for r in rows))
        nodes = sorted({x for p in pairs for x in p})
        N,E = len(nodes),len(pairs)
        check(n+'_simple_graph', len(set(pairs))==E and all(a<b for a,b in pairs))
        detail_path=STAGE3B/'overlap'/(n+'_source_excluded_overlap.tsv')
        details=read(detail_path)
        recovered=sum(r['source_excluded_recovered']=='Yes' for r in details)
        s=summaries[n]
        check(n+'_confirmed_N_E_observed', (N,E,recovered)==EXPECTED[n] and E==int(s['final_gtdb_pairs']) and recovered==int(s['source_excluded_recovered']))
        check(n+'_previous_graph_unchanged', N==int(previous_graphs[n]['node_count']) and E==int(previous_graphs[n]['edge_count']) and recovered==int(previous_graphs[n]['observed_recovered']))
        check(n+'_Stage3B_detail_keys', len(details)==E and {r['external_pair_key'] for r in details}=={' || '.join(p) for p in pairs})
        check(n+'_Stage3B_observed_rate', float(s['source_excluded_recovery_rate'])==recovered/E)
        primary='exploratory' if n=='Lacruz_MHO' else 'primary'
        check(n+'_primary_status', s['primary_or_exploratory']==m['primary_or_exploratory']==primary)
        idx={x:i for i,x in enumerate(nodes)}
        edges=[(idx[a],idx[b]) for a,b in pairs]
        obs={u*N+v for u,v in edges}
        ref=references[n]
        studies[n]=dict(network=n,n=N,E=E,nodes=nodes,edges=edges,original_keys=obs,
                        degrees=degree_vector(obs,N),ref_keys=set(),observed_recovered=recovered,observed_rate=recovered/E,
                        primary_or_exploratory=primary,reference_file=ref['file'],reference_sha256=ref['sha256'],
                        reference_pair_count=int(ref['retained_associations']),idx=idx)
        manifest.append(dict(network=n, external_pair_file=path,external_pair_sha256=m['sha256'],source_network_id=m['source_network_id'],
            edge_count=E,node_count=N,observed_input_file=str(detail_path),observed_input_sha256=hashes[str(detail_path)],
            source_excluded_database=ref['file'],reference_sha256=ref['sha256'],database_group=ref['database_group'],primary_or_exploratory=primary))
    decoder={r['database_label']:r['gtdb_r226_species'] for r in read(BASE/'stage3a/qc/03a_database_label_decoding.tsv')}
    groups=collections.defaultdict(list)
    for n,s in studies.items():groups[s['reference_file']].append(n)
    for path,names in groups.items():
        total=0
        with open(path) as f:
            r=csv.reader(f,delimiter='\t');header=next(r)
            check('reference_species_columns_'+names[0], header[:2]==['species_a','species_b'])
            for row in r:
                total+=1
                a,b=decoder[row[0]],decoder[row[1]]
                for n in names:
                    s=studies[n];idx=s['idx']
                    if a in idx and b in idx:
                        s['ref_keys'].add(encode_edge(idx[a],idx[b],s['n']))
        for n in names:
            s=studies[n]
            check(n+'_reference_pair_count',total==s['reference_pair_count'],total)
            check(n+'_independent_exact_observed_overlap',len(s['original_keys'] & s['ref_keys'])==s['observed_recovered'])
        print('Reference read and checked:',names,total,flush=True)
    graphs=[]
    for n in NETWORKS:
        s=studies[n];s.pop('idx')
        s['possible_keys']=[u*s['n']+v for u,v in itertools.combinations(range(s['n']),2)]
        for k in ['original_keys','ref_keys']:s[k]=sorted(s[k])
        serial=json.dumps(s,sort_keys=True,separators=(',',':'))
        (OUT/'prepared'/(n+'.json')).write_text(serial+'\n')
        graphs.append(dict(network=n,node_count=s['n'],edge_count=s['E'],possible_pair_count=len(s['possible_keys']),
            density=s['E']/len(s['possible_keys']),self_loop_count=0,duplicate_edge_count=0,
            primary_or_exploratory=s['primary_or_exploratory'],observed_recovered=s['observed_recovered'],observed_recovery_rate=s['observed_rate'],
            source_excluded_database_pair_count=s['reference_pair_count'],reference_pairs_within_external_node_pool=len(s['ref_keys']),graph_input_status='PASS'))
    write(OUT/'qc/03c_external_graph_audit.tsv',graphs)
    write(OUT/'qc/03c_input_manifest.tsv',manifest)
    return studies


def initialize_b(s):
    global _STUDY
    _STUDY=s
    init_b_worker(s['edges'],s['ref_keys'],s['n'],SWAPS_PER_EDGE*s['E'],seed_base(s['network'],'B'))


def validate_keys(keys,s):
    keys=list(keys);N=s['n'];unique=set(keys)
    assert len(keys)==len(unique)==s['E']
    assert all(0<=k//N<k%N<N for k in keys)
    endpoints={v for k in keys for v in divmod(k,N)}
    return endpoints


def result_row(s,null,rep,keys,recovered):
    endpoints=validate_keys(keys,s)
    row=dict(network=s['network'],null_type='Null '+null,iteration=rep,seed=seed_base(s['network'],null)+rep,
        random_edge_count=len(keys),recovered_count=recovered,recovery_rate=recovered/s['E'],
        successful_swaps='NA',swap_attempts='NA',swap_acceptance_rate='NA',degree_sequence_preserved='NA',
        per_node_degree_preserved='NA',node_set_preserved='Yes',node_pool_count=s['n'],sampled_endpoint_count=len(endpoints),
        all_endpoints_in_node_pool='Yes',no_self_loop='Yes',no_duplicate_edge='Yes',edge_count_preserved='Yes',
        random_edge_sha256=edge_hash(keys),reference_database_sha256=s['reference_sha256'],analysis_status='PASS')
    # Null A has vertex set V, including any sampled isolates. Never condition on all vertices having edges.
    if rep in (1,500,1000):row['_snapshot']=sorted(keys)
    return row


def _audit_b_result(legacy_row,U,V,edge_set):
    s=_STUDY
    assert len(U)==len(V)==s['E']
    assert len({u*s['n']+v for u,v in zip(U,V)})==s['E']
    assert {u*s['n']+v for u,v in zip(U,V)}==edge_set
    validate_keys(edge_set,s)
    degree=degree_vector(edge_set,s['n'])
    assert degree==s['degrees'] and sorted(degree)==sorted(s['degrees'])
    assert all(d>0 for d in degree)
    assert legacy_row['accepted_swaps']==SWAPS_PER_EDGE*s['E']
    assert legacy_row['attempts']>=legacy_row['accepted_swaps']
    assert len(edge_set.intersection(_B_REF))==legacy_row['recovered_count']
    row=result_row(s,'B',legacy_row['replicate'],edge_set,legacy_row['recovered_count'])
    assert row['seed']==legacy_row['seed']
    row.update(successful_swaps=legacy_row['accepted_swaps'],swap_attempts=legacy_row['attempts'],
        swap_acceptance_rate=legacy_row['swap_acceptance_rate'],degree_sequence_preserved='Yes',per_node_degree_preserved='Yes')
    return row


def one_a(s,rep):
    rng=random.Random(seed_base(s['network'],'A')+rep)
    keys=rng.sample(s['possible_keys'],s['E'])
    reference=set(s['ref_keys'])
    recovered=sum(k in reference for k in keys)
    assert recovered==len(set(keys)&reference)
    return result_row(s,'A',rep,keys,recovered)


def save_snapshot(s,null,row):
    keys=row.pop('_snapshot',None)
    if keys is None:return
    rows=[dict(network=s['network'],null_type='Null '+null,iteration=row['iteration'],seed=row['seed'],
        species_a=s['nodes'][k//s['n']],species_b=s['nodes'][k%s['n']],
        canonical_pair_key=s['nodes'][k//s['n']]+' || '+s['nodes'][k%s['n']]) for k in keys]
    write(OUT/'snapshots'/('{}_null{}_iteration_{}.tsv'.format(s['network'],null,row['iteration'])),rows)


def verify_rows(s,null,rows,complete=True):
    ids=[int(r['iteration']) for r in rows]
    assert len(ids)==len(set(ids)) and all(1<=i<=N_ITER for i in ids)
    if complete:assert sorted(ids)==list(range(1,N_ITER+1))
    for r in rows:
        assert r['network']==s['network'] and r['null_type']=='Null '+null
        assert int(r['seed'])==seed_base(s['network'],null)+int(r['iteration'])
        assert int(r['random_edge_count'])==s['E'] and int(r['node_pool_count'])==s['n']
        assert float(r['recovery_rate'])==int(r['recovered_count'])/s['E']
        assert r['reference_database_sha256']==s['reference_sha256'] and r['analysis_status']=='PASS'
        assert all(r[k]=='Yes' for k in ['no_self_loop','no_duplicate_edge','edge_count_preserved','node_set_preserved','all_endpoints_in_node_pool'])
        if null=='B':
            assert int(r['successful_swaps'])==SWAPS_PER_EDGE*s['E']
            assert int(r['swap_attempts'])>=int(r['successful_swaps'])
            assert r['degree_sequence_preserved']==r['per_node_degree_preserved']=='Yes'
            assert int(r['sampled_endpoint_count'])==s['n']


def run_model(s,null):
    path=OUT/('null_a' if null=='A' else 'null_b')/(s['network']+'_null'+null+'_iterations.tsv')
    old=read(path) if path.exists() else []
    verify_rows(s,null,old,False)
    done={int(r['iteration']) for r in old}
    missing=[i for i in range(1,N_ITER+1) if i not in done]
    print(s['network'],'Null '+null,'resume',len(done),'/',N_ITER,flush=True)
    if missing:
        with path.open('a',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=ROW_FIELDS,delimiter='\t',lineterminator='\n')
            if f.tell()==0:writer.writeheader()
            def consume(rows):
                count=len(done)
                for row in rows:
                    save_snapshot(s,null,row)
                    writer.writerow(row);f.flush();count+=1
                    if count%100==0:
                        print(s['network'],'Null '+null,count,'/',N_ITER,flush=True)
            if null=='A':consume(one_a(s,rep) for rep in missing)
            else:
                with mp.get_context('fork').Pool(WORKERS,initializer=initialize_b,initargs=(s,)) as pool:
                    consume(pool.imap_unordered(null_b_worker,missing,chunksize=1))
    rows=sorted(read(path),key=lambda r:int(r['iteration']))
    verify_rows(s,null,rows)
    write(path,rows,ROW_FIELDS)
    check(s['network']+'_Null'+null+'_all_1000_iterations_QC',True,'Every graph, seed, reference, E and applicable per-node degree / 50E checks passed')
    # Validate saved edge-list evidence against iteration statistics without extra randomizations.
    for rep in (1,500,1000):
        snap=read(OUT/'snapshots'/('{}_null{}_iteration_{}.tsv'.format(s['network'],null,rep)))
        idx={name:i for i,name in enumerate(s['nodes'])}
        keys=[encode_edge(idx[r['species_a']],idx[r['species_b']],s['n']) for r in snap]
        validate_keys(keys,s)
        r=rows[rep-1]
        check(s['network']+'_Null'+null+'_snapshot_'+str(rep),
            edge_hash(keys)==r['random_edge_sha256'] and len(set(keys)&set(s['ref_keys']))==int(r['recovered_count']) and
            (null=='A' or degree_vector(keys,s['n'])==s['degrees']))
    frame=pd.DataFrame(rows)
    study=SimpleNamespace(name=s['network'],nodes=s['nodes'],E=s['E'],observed_recovered=s['observed_recovered'],observed_rate=s['observed_rate'])
    # This is the verbatim legacy summarize() function compiled above.
    legacy=summarize(study,'Null '+null,frame,N_ITER)
    sd=legacy['null_sd_recovery_rate'];z=legacy['z_score']
    check(s['network']+'_Null'+null+'_zero_SD_rule',
        (sd==0 and math.isinf(z)) or (sd>0 and math.isclose(z,(s['observed_rate']-legacy['null_mean_recovery_rate'])/sd,rel_tol=1e-12,abs_tol=1e-10)))
    counts=np.asarray([int(r['recovered_count']) for r in rows],dtype=float)
    n_ge=sum(x>=s['observed_recovered'] for x in counts)
    check(s['network']+'_Null'+null+'_empirical_P', legacy['empirical_one_sided_P']==(n_ge+1)/1001 and legacy['empirical_one_sided_P']>0)
    return dict(network=s['network'],analysis_status='PASS',primary_or_exploratory=s['primary_or_exploratory'],
        node_count=s['n'],edge_count=s['E'],source_excluded_database_pair_count=s['reference_pair_count'],
        observed_recovered=s['observed_recovered'],observed_recovery_rate=s['observed_rate'],null_type='Null '+null,null_iterations=len(rows),
        null_mean_recovered=legacy['null_mean_recovered_count'],null_sd_recovered=legacy['null_sd_recovered_count'],
        null_mean_recovery_rate=legacy['null_mean_recovery_rate'],null_sd_recovery_rate=sd,
        null_median_recovery_rate=float(np.median(counts/s['E'])),null_q025_recovery_rate=legacy['null_2.5pct_recovery_rate'],
        null_q975_recovery_rate=legacy['null_97.5pct_recovery_rate'],null_min_recovery_rate=float(np.min(counts/s['E'])),
        null_max_recovery_rate=float(np.max(counts/s['E'])),recovery_enrichment=legacy['observed_to_null_mean_enrichment'],
        z_score=z,null_sd_zero='Yes' if sd==0 else 'No',n_ge_observed=legacy['null_replicates_ge_observed'],empirical_p=legacy['empirical_one_sided_P'],
        master_seed=MASTER_SEED,reference_database_sha256=s['reference_sha256'],reference_database=s['reference_file'])


def main():
    started=time.time()
    protocol=update_protocol()
    load_legacy()
    previous=read(HISTORY/'qc/03c_input_integrity.tsv')
    integrity=[]
    print('Checking frozen hashes before execution...',flush=True)
    for r in previous:
        h=sha(r['input_file'])
        check('frozen_prior_hash:'+r['input_file'],h==r['sha256_after'])
        integrity.append(dict(input_file=r['input_file'],scope='frozen_upstream_or_legacy_protocol',size_bytes=Path(r['input_file']).stat().st_size,
            sha256_before=h,prior_recorded_sha256=r['sha256_after'],matches_prior_recorded_hash='PASS'))
    for p in [Path(__file__).resolve(),AUTHORIZATION,HISTORY/'qc/03c_external_graph_audit.tsv',HISTORY/'qc/03c_input_manifest.tsv']:
        integrity.append(dict(input_file=str(p),scope='Stage3C_execution_or_confirmed_input',size_bytes=p.stat().st_size,
            sha256_before=sha(p),prior_recorded_sha256='NA',matches_prior_recorded_hash='NO_PRIOR_HASH'))
    studies=prepare(integrity)
    config=dict(master_seed=MASTER_SEED,iterations=N_ITER,successful_swaps_per_edge=SWAPS_PER_EDGE,
        workers=WORKERS,network_order=NETWORKS,seed_formula='master_seed + 100000*(2*zero_based_network_index + (1 for A, 2 for B)) + iteration',
        legacy_script=str(LEGACY),legacy_sha256=sha(LEGACY),execution_script_sha256=sha(Path(__file__)),
        python_version=sys.version,python_executable=sys.executable,numpy_version=np.__version__,pandas_version=pd.__version__,
        zero_sd_zscore_rule='inf',sd_ddof=1,quantile_method='numpy default linear',
        prepared_hashes={n:sha(OUT/'prepared'/(n+'.json')) for n in NETWORKS},
        null_a_isolated_nodes='Allowed; vertex set V fixed, no conditioning on sampled endpoint coverage',
        exact_pair_rule='Frozen GTDB R226 species names; unordered; preserve case and suffixes; no legacy label normalization',
        reference_policy='Read Stage3B study-specific derivative files; do not rerun source exclusion')
    cp=OUT/'qc/03c_run_config.json'
    if cp.exists():check('resume_configuration_identical',json.loads(cp.read_text())==config)
    else:cp.write_text(json.dumps(config,indent=2)+'\n')
    seeds=[dict(network=n,null_type='Null '+t,iteration=i,seed=seed_base(n,t)+i) for n in NETWORKS for t in 'AB' for i in range(1,N_ITER+1)]
    check('all_18000_seeds_unique',len({r['seed'] for r in seeds})==18000)
    write(OUT/'qc/03c_iteration_seed_manifest.tsv',seeds)
    write(OUT/'qc/03c_preflight_qc.tsv',[dict(check='protocol_consistency',status='PASS',detail='All core rows MATCH, including zero_sd_zscore_rule=inf'),
        dict(check='input_graph_and_observed',status='PASS',detail='All nine confirmed graphs, observed values, and study-specific reference hashes verified')])
    (OUT/'qc/03c_run_status.json').write_text(json.dumps(dict(overall_status='RUNNING',protocol_consistency='MATCH',random_evaluations_required=18000),indent=2)+'\n')
    results=[]
    for n in NETWORKS:
        for null in 'AB':
            results.append(run_model(studies[n],null))
    check('18000_randomized_networks_completed',sum(r['null_iterations'] for r in results)==18000)
    write(OUT/'summary/03c_null_model_summary.tsv',results)
    wide=[]
    for n in NETWORKS:
        rs={r['null_type']:r for r in results if r['network']==n};s=studies[n]
        row=dict(network=n,analysis_status='PASS',primary_or_exploratory=s['primary_or_exploratory'],node_count=s['n'],edge_count=s['E'],
            observed_recovered=s['observed_recovered'],observed_recovery_rate=s['observed_rate'])
        for t in 'AB':
            r=rs['Null '+t];prefix='Null'+t+'_'
            for short,long in [('mean','null_mean_recovery_rate'),('enrichment','recovery_enrichment'),('z','z_score'),('null_sd_zero','null_sd_zero'),
                               ('n_ge_observed','n_ge_observed'),('empirical_p','empirical_p')]:row[prefix+short]=r[long]
            row[prefix+'95CI']='{};{}'.format(r['null_q025_recovery_rate'],r['null_q975_recovery_rate'])
        wide.append(row)
    write(OUT/'summary/03c_null_model_paper_table.tsv',wide)
    print('Final frozen hash verification...',flush=True)
    for r in integrity:
        r['sha256_after']=sha(r['input_file']);r['unchanged']=str(r['sha256_after']==r['sha256_before'])
    check('all_frozen_hashes_unchanged',all(r['unchanged']=='True' for r in integrity),len(integrity))
    write(OUT/'qc/03c_input_integrity.tsv',integrity)
    CHECKS.append(dict(check='overall_status',status='PASS',detail='18,000 primary evaluations complete; all protocol, graph, swap and input-integrity checks passed'))
    write(OUT/'qc/03c_null_model_qc.tsv',CHECKS)
    zero=[r['network']+' '+r['null_type'] for r in results if r['null_sd_zero']=='Yes']
    nonsig=[r['network']+' '+r['null_type'] for r in results if r['empirical_p']>0.05]
    status=dict(overall_status='PASS',protocol_consistency='MATCH',random_evaluations_completed=18000,
        null_a_iterations=9000,null_b_iterations=9000,null_b_all_successful_swaps_50E=True,
        total_successful_swaps=1000*50*sum(studies[n]['E'] for n in NETWORKS),null_sd_zero_tests=zero,empirical_p_gt_005=nonsig,
        frozen_inputs_verified=len(integrity),qc_checks=len(CHECKS),elapsed_seconds=time.time()-started)
    (OUT/'qc/03c_run_status.json').write_text(json.dumps(status,indent=2)+'\n')
    print(json.dumps(status,indent=2),flush=True)


if __name__=='__main__':
    try:main()
    except Exception as exc:
        (OUT/'qc/03c_run_status.json').write_text(json.dumps(dict(overall_status='FAIL',error=str(exc),
            null_iterations_policy='Never lower N_ITER or 50E target; partial TSVs retained for audited resume'),indent=2)+'\n')
        raise
