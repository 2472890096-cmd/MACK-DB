"""Explicit data/code path registry for opt-in supplementary analyses.
No scientific parameters, credentials, network access, or automatic downloads.
"""
from pathlib import Path
import json,os
ROOT=Path(__file__).resolve().parents[1]
DEFAULTS={
 'external_root':'analysis/external_validation',
 'literature':'analysis/external_validation/literature_networks',
 'stage5b':'analysis/external_validation/stage5B',
 'shared_protocol':'analysis/external_validation/shared_protocol',
 'gtdb_r226':'data/reference/gtdb_r226',
 'project_data':'data/project_inputs',
 'species_input':'data/species',
 'genus_input':'data/genus',
 'environment':'data/env',
 'environment_network':'data/networks/human_gut.info',
 'legacy_external_inputs':'data/external/legacy_inputs',
 'null_protocol_evidence':'data/external/null_protocol_authorization.txt',
 'count_script':'pipeline/abundance/count_sparcc.R',
 'fastspar_script':'pipeline/associations/fastspar/run_fastspar.sh',
 'genus_capture':'pipeline/epi_v3/genus/scripts/capture_v3.R',
}
def resolve_path(key,*parts):
 config={}
 if os.environ.get('MACKDB_PATHS_CONFIG'):
  config=json.loads(Path(os.environ['MACKDB_PATHS_CONFIG']).read_text())
 value=config.get(key,DEFAULTS.get(key))
 if value is None:raise KeyError('Unconfigured MACK-DB path: '+key)
 p=Path(value).expanduser()
 if not p.is_absolute():p=ROOT/p
 return p.joinpath(*parts)
