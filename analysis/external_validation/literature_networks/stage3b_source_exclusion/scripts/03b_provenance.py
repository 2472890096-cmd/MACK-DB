#!/usr/bin/env python3
"""Read frozen metadata and downloaded primary evidence; write Stage 3B only."""
import collections
import csv
import hashlib
import json
import re
import xml.etree.ElementTree as ET
from pathlib import Path

import openpyxl

OUT = Path(__file__).resolve().parents[1]
BASE = OUT.parent
RAW = OUT / 'provenance/raw'


def read(path):
    with Path(path).open(newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def write(path, rows, fields=None):
    rows = list(rows)
    with Path(path).open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields or list(rows[0]), delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)


def main():
    qc = []

    def check(name, ok, detail):
        qc.append(dict(check=name, status='PASS' if ok else 'FAIL', detail=str(detail)))
        if not ok:
            write(OUT / 'qc/03b_source_provenance_qc.tsv', qc)
            raise AssertionError((name, detail))

    alias = {}

    def add_alias(paper, reported, canonical, evidence, kind=None):
        key = (paper, reported)
        if key in alias:
            assert alias[key]['canonical_project_accession'] == canonical
            return
        alias[key] = dict(paper=paper, reported_accession=reported,
                          accession_type=kind or ('BioProject' if reported.startswith('PRJ') else 'study'),
                          canonical_project_accession=canonical, evidence=evidence, status='PASS')

    studies = {'Zhou2022': ['PRJEB32762'], 'Guo2023': ['PRJNA751448'],
               'Lacruz2026': ['PRJEB87701', 'PRJEB7774', 'PRJEB1786', 'PRJEB41311', 'PRJEB38742', 'PRJEB37249']}
    runs = {}
    for paper, projects in studies.items():
        for project in projects:
            f = RAW / (project + '_ena_runs.tsv')
            data = read(f)
            check(project + '_ENA_identity', bool(data) and {r['study_accession'] for r in data} == {project}, len(data))
            add_alias(paper, project, project, str(f))
            for secondary in sorted({r['secondary_study_accession'] for r in data}):
                add_alias(paper, secondary, project, str(f))
            if paper == 'Lacruz2026':
                for r in data:
                    assert r['run_accession'] not in runs
                    runs[r['run_accession']] = r
    check('Zhou_ERP115476_alias', alias[('Zhou2022', 'ERP115476')]['canonical_project_accession'] == 'PRJEB32762', 'one source project')

    # Extract only Table S1's actual labelled cohort rows, retaining row-level provenance.
    workbook = RAW / 'ADVS-11-2308313-s002.xlsx'
    sheet = openpyxl.load_workbook(workbook, read_only=True, data_only=True)['Table_S1']
    cohorts = []
    for rowno, row in enumerate(sheet.values, 1):
        if row[2] not in ('Discovery', 'Validation', 'Rifaximin treatment'):
            continue
        accession = re.search(r'PRJ[A-Z]+\d+', row[1]).group()
        cohorts.append(dict(table='Table_S1', excel_row=rowno, sequencing=row[0], reported_cohort=row[1],
                            role=row[2], canonical_project_accession=accession, reference=row[14],
                            notes=row[16] or '', source_file=str(workbook)))
    discovery = [r for r in cohorts if r['role'] == 'Discovery']
    li_sources = sorted({r['canonical_project_accession'] for r in discovery})
    validation = sorted({r['canonical_project_accession'] for r in cohorts if r['role'] == 'Validation'})
    treatment = sorted({r['canonical_project_accession'] for r in cohorts if r['role'] == 'Rifaximin treatment'})
    check('Li_complete_14_discovery_cohorts', len(discovery) == 14, [r['reported_cohort'] for r in discovery])
    check('Li_13_16S_plus_one_shotgun', collections.Counter(r['sequencing'] for r in discovery) ==
          {'16S rRNA amplicon': 13, 'Shotgun metagenome': 1}, 'Table S1 rows 4-17')
    check('Li_14_cohorts_are_13_projects', len(li_sources) == 13 and
          sum(r['canonical_project_accession'] == 'PRJEB42304' for r in discovery) == 2,
          'PRJEB42304: 16S and WGS; never count as two independent projects')
    check('Li_validation_excluded', set(validation) == {'PRJEB37924', 'PRJEB11419', 'PRJEB19857', 'PRJNA812699'}
          and not set(validation) & set(li_sources), validation)
    check('Li_rifaximin_discovery_overlap_preserved', set(treatment) & set(li_sources) ==
          {'PRJEB11252', 'PRJNA475187', 'PRJCA002555'}, treatment)
    check('Li_rifaximin_only_excluded', set(treatment) - set(li_sources) == {'PRJNA279031', 'PRJNA391915'}, 'not network sources')
    for r in cohorts:
        r['include_in_network_source_set'] = 'Yes' if r['canonical_project_accession'] in li_sources else 'No'
        add_alias('Li2024', r['reported_cohort'], r['canonical_project_accession'],
                  str(workbook) + ':Table_S1:row' + str(r['excel_row']), 'labelled BioProject')
    for project in li_sources + ['PRJNA834801']:
        paper = 'Wallen2022' if project == 'PRJNA834801' else 'Li2024'
        if project == 'PRJCA002555':
            f = RAW / 'PRJCA002555_ngdc.html'
            check('Li_PRJCA_native_accession', 'PRJCA002555' in f.read_text() and 'CRA002557' in f.read_text(),
                  'Official NGDC project retained as PRJCA002555. No unproven INSDC equivalent; page search-box example is not an alias.')
            add_alias(paper, project, project, str(f), 'NGDC BioProject')
            continue
        f = RAW / (project + '_ena_study.xml')
        x = ET.parse(f).getroot().find('PROJECT')
        check(project + '_project_identity', x.attrib['accession'] == project and
              x.findtext('IDENTIFIERS/PRIMARY_ID') == project, str(f))
        add_alias(paper, project, project, str(f))
        for el in x.findall('IDENTIFIERS/SECONDARY_ID'):
            add_alias(paper, el.text, project, str(f))
        # Parent, umbrella, related project and page-example IDs are deliberately not aliases.
    write(OUT / 'provenance/03b_li_table_s1_cohort_audit.tsv', cohorts)

    # Verify the frozen Stage 1 sample manifest against the source workbook without rewriting it.
    frozen = read(BASE / 'parsed/lacruz2026_sample_source_manifest.tsv')
    sw = openpyxl.load_workbook(BASE / 'raw/lacruz2026/Lacruz2026_Supplementary_Data3.xlsx', read_only=True, data_only=True)['SuppData3']
    srows = list(sw.values)[2:]
    check('Lacruz_frozen_manifest_matches_Data3', len(srows) == len(frozen) == 931 and
          all(tuple(str(v) for v in row) == tuple(r[k] for k in ('subject_id', 'run_accession', 'study_name', 'phenotype'))
              for row, r in zip(srows, frozen)), '931 exact source metadata records; frozen files not regenerated')
    with (BASE / 'raw/lacruz2026/whole_cohort_data.csv').open() as f:
        headers = next(csv.reader(f))[1:]
    check('Lacruz_abundance_samples_unique', len(headers) == len(set(headers)) == 931, len(headers))
    hs = set(headers)
    evidence, counts, sets, matched = [], collections.Counter(), collections.defaultdict(set), []
    for r in frozen:
        sample_runs = r['run_accession'].split(';')
        check('Lacruz_run_resolution_' + r['subject_id'], all(x in runs for x in sample_runs), r['run_accession'])
        info = [runs[x] for x in sample_runs]
        projects = {x['study_accession'] for x in info}
        check('Lacruz_sample_single_project_' + r['subject_id'], len(projects) == 1, sorted(projects))
        project = next(iter(projects))
        subject = r['subject_id']
        if r['study_name'] == 'AI4Food':
            aliases = {x['sample_alias'] for x in info}
            assert len(aliases) == 1
            ena_alias = next(iter(aliases))
            match = re.fullmatch(r'A4F_(V\d+)_(\d+)', ena_alias)
            check('Lacruz_AI4Food_ENA_subject_' + subject, bool(match) and match[2] == subject, ena_alias)
            abundance_id = 'A4F_' + match[2] + '_' + match[1]
            rule = 'ENA sample_alias A4F_Vn_subject -> abundance A4F_subject_Vn; exact subject and visit'
        else:
            ena_alias = ';'.join(sorted({x['sample_alias'] for x in info}))
            abundance_id, rule = subject, 'exact subject_id'
        check('Lacruz_abundance_membership_' + subject, abundance_id in hs, abundance_id)
        matched.append(abundance_id)
        network = 'Lacruz_' + r['phenotype']
        counts[(network, project)] += 1
        sets[network].add(project)
        evidence.append(dict(network=network, subject_id=subject, phenotype=r['phenotype'], study_name=r['study_name'],
                             run_accessions=r['run_accession'], canonical_project_accession=project,
                             ena_sample_alias=ena_alias, abundance_sample_id=abundance_id,
                             abundance_membership='Yes', identity_rule=rule,
                             frozen_metadata_file=str(BASE / 'parsed/lacruz2026_sample_source_manifest.tsv'),
                             metadata_excel_row=r['source_row'], ena_evidence=str(RAW / (project + '_ena_runs.tsv'))))
    check('Lacruz_all_samples_bijective', len(set(matched)) == 931 and set(matched) == hs, '833 direct + 98 ENA-resolved AI4Food')
    check('Lacruz_phenotype_counts', collections.Counter(r['phenotype'] for r in frozen) ==
          {'MHNO': 103, 'MHO': 19, 'MUNO': 377, 'MUO': 432}, 'Data3 matches create_visualize_networks.R lines 18-21')
    code = (RAW / 'lacruz_create_networks.R').read_text()
    check('Lacruz_network_preparation_evidence', all(('MetObesity == "' + g + '"') in code for g in ['MHNO', 'MHO', 'MUNO', 'MUO'])
          and all(('edgeweights_' + g + '.csv') in code for g in ['mhno', 'mho', 'muno', 'muo']), str(RAW / 'lacruz_create_networks.R'))
    check('Lacruz_PRJEB41311_not_in_network_samples', not any('PRJEB41311' in v for v in sets.values()),
          'Mentioned by article but zero Data3 run IDs map to this project; excluded from all source sets')
    check('Lacruz_source_sets_are_different', len({tuple(sorted(v)) for v in sets.values()}) == 3, dict((k, sorted(v)) for k, v in sets.items()))
    write(OUT / 'provenance/03b_lacruz_sample_project_evidence.tsv', evidence)
    write(OUT / 'provenance/03b_lacruz_phenotype_project_counts.tsv',
          [dict(network=n, canonical_project_accession=p, construction_sample_count=counts[(n, p)],
                include_in_source_set='Yes' if counts[(n, p)] else 'No')
           for n in ['Lacruz_MHNO', 'Lacruz_MHO', 'Lacruz_MUNO', 'Lacruz_MUO'] for p in studies['Lacruz2026']])
    config = {}
    for n in ['Zhou_MS', 'Zhou_HHC']:
        config[n] = dict(paper='Zhou2022', sources=['PRJEB32762'], database_group='Zhou')
    for n in ['Guo_AggControl', 'Guo_AggMECFS']:
        config[n] = dict(paper='Guo2023', sources=['PRJNA751448'], database_group='Guo')
    for n in ['Lacruz_MHNO', 'Lacruz_MHO', 'Lacruz_MUNO', 'Lacruz_MUO']:
        config[n] = dict(paper='Lacruz2026', sources=sorted(sets[n]),
                         database_group='Lacruz_MUNO_MUO' if n in ('Lacruz_MUNO', 'Lacruz_MUO') else n)
    config['Li_IBS'] = dict(paper='Li2024', sources=li_sources, database_group='Li')
    # Match source statements in primary papers, not only user-provided accessions.
    check('Zhou2022_article_source_statement', 'ERP115476' in (RAW / 'Zhou_article.html').read_text(),
          'https://pmc.ncbi.nlm.nih.gov/articles/10143502/; local full HTML + complete ENA run metadata')
    for paper, fname, expected in [('Guo2023', 'Guo_fulltext.xml', 'PRJNA751448')]:
        xml = ET.parse(RAW / fname)
        paragraphs = [''.join(p.itertext()).strip() for p in xml.iter('p')]
        paragraphs = [p for p in paragraphs if expected in p]
        check(paper + '_article_source_statement', bool(paragraphs), expected + '; ' + str(RAW / fname))
    lit = ET.parse(RAW / 'Li_fulltext.xml')
    text = ' '.join(lit.getroot().itertext())
    check('Li_network_constructed_from_discovery', 'microbiomes of the discovery cohorts' in text and '742 edges' in text,
          'Li primary Methods: network construction; raw frozen 742-edge selection unchanged')
    (OUT / 'provenance/03b_network_source_config.json').write_text(json.dumps(config, indent=2) + '\n')
    write(OUT / 'provenance/03b_accession_alias_map.tsv', sorted(alias.values(), key=lambda r:(r['paper'], r['reported_accession'])))
    check('aliases_unambiguous', all(r['status'] == 'PASS' for r in alias.values()), len(alias))
    qc.append(dict(check='overall_status', status='PASS', detail='All source provenance checks passed.'))
    write(OUT / 'qc/03b_source_provenance_qc.tsv', qc)
    print(json.dumps(config, indent=2))
    print('provenance QC:', len(qc), 'PASS')


if __name__ == '__main__':
    main()
