#!/usr/bin/env python3
"""Small safety/IO regressions; no network or model training; not a real-run PASS."""
import gzip
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('single_run_example', ROOT/'pipeline/metakssd_profile/run_single_run.py')
M = importlib.util.module_from_spec(spec)
spec.loader.exec_module(M)


class Contracts(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='mackdb-single-test-')
        self.d = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def mate(self, name, identifier='test/1', sequence='ACGT'):
        f = self.d/name
        with gzip.open(f, 'wb') as h:
            h.write(('@'+identifier+'\n'+sequence+'\n+\n'+'I'*len(sequence)+'\n').encode())
        return f

    def test_pair_and_base_totals(self):
        a, b = self.mate('R1.gz'), self.mate('R2.gz', 'test/2')
        r = M.validate_pair(a, b, {'read_count': 2, 'base_count': 8})
        self.assertEqual(r['mate_records'], [1, 1])
        self.assertTrue(r['gzip_trailers_validated'])

    def test_wrong_pair_id_rejected(self):
        with self.assertRaisesRegex(ValueError, 'Paired record'):
            M.validate_pair(self.mate('a.gz'), self.mate('b.gz', 'other/2'), {'read_count': 2, 'base_count': 8})

    def test_truncated_gzip_rejected(self):
        a, b = self.mate('a.gz'), self.mate('b.gz', 'test/2')
        a.write_bytes(a.read_bytes()[:-5])
        with self.assertRaises((EOFError, OSError)):
            M.validate_pair(a, b, {'read_count': 2, 'base_count': 8})

    def test_wrong_read_count_rejected(self):
        with self.assertRaisesRegex(ValueError, 'totals'):
            M.validate_pair(self.mate('a.gz'), self.mate('b.gz', 'test/2'), {'read_count': 4, 'base_count': 8})

    def test_foreign_scientific_parameter_rejected(self):
        cfg = dict(run_accession='ERR14036375', fastq_r1='a.gz', fastq_r2='b.gz', metakssd_bin=sys.executable,
                   shuffle='ref.shuf', markerdb_dir='ref', rscript=sys.executable, scratch_dir='/tmp', threads=8)
        f = self.d/'config.json'; f.write_text(json.dumps(cfg))
        with self.assertRaisesRegex(ValueError, 'only documented'):
            M.settings(f)

    def test_relative_input_paths_and_path_spaces(self):
        cfg = dict(run_accession='ERR14036375', fastq_r1='input with spaces/a.gz', fastq_r2='b.gz',
                   metakssd_bin=sys.executable, shuffle='ref.shuf', markerdb_dir='ref', rscript=sys.executable, scratch_dir='/tmp')
        f = self.d/'config.json'; f.write_text(json.dumps(cfg))
        self.assertEqual(M.settings(f)['fastq_r1'], self.d/'input with spaces/a.gz')

    def test_wrong_input_sha_rejected_before_tool_execution(self):
        a, b = self.mate('a.gz'), self.mate('b.gz', 'test/2')
        cfg = dict(run_accession='ERR14036375', fastq_r1='a.gz', fastq_r2='b.gz',
                   metakssd_bin=sys.executable, shuffle='ref.shuf', markerdb_dir='ref', rscript=sys.executable, scratch_dir='/tmp')
        f = self.d/'config.json'; f.write_text(json.dumps(cfg))
        expected = dict(run_accession='ERR14036375', metakssd_sha256='0'*64,
                        fastq={'R1': {'bytes': a.stat().st_size, 'md5': M.digest(a, 'md5'), 'sha256': '0'*64}})
        with mock.patch.object(M, 'pins', return_value=expected), mock.patch.object(M, 'probe', side_effect=AssertionError('Tool must not execute')):
            with self.assertRaisesRegex(ValueError, 'FASTQ SHA256 mismatch'):
                M.preflight(f)

    def test_existing_output_refused_without_reading_config(self):
        with self.assertRaisesRegex(ValueError, 'already exists'):
            M.execute(self.d/'missing_config.json', self.d)

    def test_owned_cleanup_preserves_inputs_and_rejects_other_directory(self):
        original = self.d/'user_input'; original.write_text('keep')
        owned = self.d/'mackdb-upstream-test'; owned.mkdir(); (owned/'.owner').write_text('testtoken')
        (owned/'generated.fastq').write_text('generated')
        M.safe_cleanup(owned, self.d, 'testtoken')
        self.assertFalse(owned.exists()); self.assertEqual(original.read_text(), 'keep')
        with self.assertRaises(ValueError):
            M.safe_cleanup(self.d, self.d.parent, 'testtoken')
        self.assertTrue(original.exists())

    def test_stage_nonzero_exit_is_failure(self):
        with self.assertRaisesRegex(ValueError, 'failed'):
            M.run_command(self.d, 'nonzero', [sys.executable, '-c', 'raise SystemExit(7)'])
        row = json.loads((self.d/'nonzero.status.json').read_text())
        self.assertEqual(row['exit_code'], 7)

    def test_timeout_stops_only_owned_child_session(self):
        with self.assertRaisesRegex(RuntimeError, 'timed out'):
            M.run_command(self.d, 'timeout', [sys.executable, '-c', 'import time;time.sleep(5)'], timeout=.15)
        self.assertEqual(json.loads((self.d/'timeout.status.json').read_text())['status'], 'TIMEOUT_OR_INTERRUPTED')

    def test_sra_mode_is_nonexecuting_plan(self):
        q = subprocess.run([sys.executable, ROOT/'pipeline/metakssd_profile/run_single_run.py', 'plan-sra',
                            '--sra-file', self.d/'nonexistent.sra'], capture_output=True, text=True, check=True)
        r = json.loads(q.stdout)
        self.assertEqual(r['status'], 'NOT_VERIFIED_SRA_STREAM_PLAN_ONLY')
        self.assertFalse(r['execution_supported']); self.assertIsNone(r['aws_download_command'])
        self.assertEqual(list(self.d.iterdir()), [])


if __name__ == '__main__':
    unittest.main(verbosity=2)
