#!/usr/bin/env python3
"""Small adapter checks; no random simulation or cohort analysis."""
import ast
import importlib.util
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock
from types import SimpleNamespace

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('paper_support_replay', HERE / 'replay.py')
replay = importlib.util.module_from_spec(spec)
spec.loader.exec_module(replay)


class ReplayTests(unittest.TestCase):
    def test_original_null_summary_sample_sd_and_upper_tail(self):
        tree = ast.parse(replay.LEGACY.read_text())
        function = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == 'summarize')
        ns = dict(np=replay.np, pd=replay.pd)
        exec(compile(ast.Module(body=[function], type_ignores=[]), str(replay.LEGACY), 'exec'), ns)
        study = SimpleNamespace(name='synthetic', E=4, nodes=range(4), observed_recovered=3, observed_rate=.75)
        result = ns['summarize'](study, 'Null A', replay.pd.DataFrame({'recovered_count':[1, 3]}), 2)
        self.assertEqual(result['null_mean_recovered_count'], 2)
        self.assertEqual(result['null_sd_recovered_count'], math.sqrt(2))
        self.assertEqual(result['empirical_one_sided_P'], 2/3)
        zero = ns['summarize'](study, 'Null A', replay.pd.DataFrame({'recovered_count':[1, 1]}), 2)
        self.assertEqual(zero['z_score'], float('inf'))

    def test_duplicate_keys_rejected(self):
        with self.assertRaisesRegex(ValueError, 'unique_keys'):
            replay.compare_rows([{'k':'a'}, {'k':'a'}], [{'k':'a'}], ['k'], 'duplicate')

    def test_changed_scientific_source_rejected_before_compile(self):
        with mock.patch.object(replay, 'sha', return_value='0' * 64):
            with self.assertRaisesRegex(ValueError, 'scientific_source_SHA256'):
                replay.load_functions(replay.LEGACY, ['summarize'], {})

    def test_output_escape_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            replay.OUTPUT = Path(tmp).resolve()
            with self.assertRaisesRegex(ValueError, 'escapes'):
                replay.write(Path(tmp).parent / 'escape.tsv', [{'a':1}])

    def test_existing_output_rejected_before_input_access(self):
        with tempfile.TemporaryDirectory() as tmp:
            proc = subprocess.run([sys.executable, '-B', str(HERE/'replay.py'), 'nulls',
                '--literature-root', tmp, '--chen-wallen-root', tmp, '--output', tmp],
                capture_output=True, text=True, timeout=15)
            self.assertEqual(proc.returncode, 2)
            self.assertIn('new directory', proc.stderr)


if __name__ == '__main__':
    unittest.main()
