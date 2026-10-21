#!/usr/bin/env python3
"""Small IO/command-planning tests; no scientific inference or training."""
from pathlib import Path
import importlib.util
import json
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'examples/real_bioproject/PRJEB81106/run_workflow.py'
spec = importlib.util.spec_from_file_location('portable_project', SOURCE)
workflow = importlib.util.module_from_spec(spec)
spec.loader.exec_module(workflow)


class PortableWorkflowTests(unittest.TestCase):
    def test_epi_has_no_network_dependency(self):
        self.assertEqual(workflow.plan_steps(['epi']), ['prepare', 'epi'])

    def test_zipi_consumes_new_network_stage(self):
        self.assertEqual(workflow.plan_steps(['zipi']), ['prepare', 'fastspar', 'zipi'])

    def test_dki_score_requires_fresh_training(self):
        self.assertEqual(workflow.plan_steps(['dki-score']), ['prepare', 'dki-input', 'dki-train', 'dki-score'])

    def test_independent_branches_are_deduplicated(self):
        self.assertEqual(workflow.plan_steps(['zipi', 'epi', 'epi']), ['prepare', 'fastspar', 'zipi', 'epi'])

    def test_bundled_complete_input_hashes_and_schema(self):
        self.assertEqual(len(workflow.checked_bundle()), 64)

    def test_resume_rejects_changed_output(self):
        with tempfile.TemporaryDirectory() as name:
            folder = Path(name)
            (folder / 'result.tsv').write_text('a\tb\n1\t2\n')
            record = dict(status='COMPLETE_VALIDATED', artifacts=workflow.artifact_records(folder))
            workflow.atomic_json(folder / 'SUCCESS.json', record)
            self.assertEqual(workflow.verify_stage(folder)['status'], 'COMPLETE_VALIDATED')
            (folder / 'result.tsv').write_text('changed\n')
            with self.assertRaises(ValueError):
                workflow.verify_stage(folder)

    def test_malformed_prediction_cannot_be_success(self):
        with tempfile.TemporaryDirectory() as name:
            p = Path(name) / 'qtst.csv'
            p.write_text('0.5,0.5\n')
            workflow.finite_csv(p, 1, 2, normalized=True)
            p.write_text('nan,1\n')
            with self.assertRaises(ValueError):
                workflow.finite_csv(p, 1, 2, normalized=True)


if __name__ == '__main__':
    unittest.main()
