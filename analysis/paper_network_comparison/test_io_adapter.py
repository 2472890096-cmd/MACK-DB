#!/usr/bin/env python3
"""Small path/source-safety tests; no scientific data or analysis is executed."""
import hashlib
import tempfile
import unittest
from pathlib import Path

from run_original_statistics import install_io_rebase, load_functions
from validate_statistics import unique_index
from fixed_order_recovery import archived_steps


class AdapterTests(unittest.TestCase):
    def test_duplicate_archived_position_is_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            source = Path(folder) / "archive.tsv"
            source.write_text("analysis\torder_id\tposition\n" + "".join(
                f"primary_all_known_source_excluded\t22\t{i}\n" for i in [1, 2, 3, 4, 5, 6, 6]))
            with self.assertRaisesRegex(ValueError, "Missing/duplicate"):
                archived_steps(source)

    def test_duplicate_actual_statistic_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "Duplicate actual statistics key"):
            unique_index([{"network": "Guo_Common"}, {"network": "Guo_Common"}], ["network"], 2)

    def test_rebase_is_idempotent_with_nested_destination(self):
        ns = {"ROOT": Path("/archive/work/out"), "OUT": Path("/archive/work/out")}
        translate = install_io_rebase(ns, {"/archive": "/archive/relocated"})
        value = translate("/archive/data/file.tsv")
        self.assertEqual(value, "/archive/relocated/data/file.tsv")
        self.assertEqual(translate(value), value)
        self.assertEqual(translate("/archive/work/out/table.tsv"), "/archive/work/out/table.tsv")

    def test_prefix_boundary_and_hash_are_preserved(self):
        ns = {"ROOT": Path("/result"), "read": lambda path: [
            {"path": "/old/data/x.tsv", "sha256": "abc123", "name": "Genus_A species"}]}
        translate = install_io_rebase(ns, {"/old/data": "/new/data"})
        self.assertEqual(translate("/old/database/x.tsv"), "/old/database/x.tsv")
        self.assertEqual(ns["read"]("/old/data/file.tsv")[0],
                         {"path": "/new/data/x.tsv", "sha256": "abc123", "name": "Genus_A species"})

    def test_wrong_source_rejected_before_execution(self):
        with tempfile.TemporaryDirectory() as folder:
            source = Path(folder) / "producer.py"
            source.write_text("raise AssertionError('top-level must not run')\ndef value(): return 3\n")
            with self.assertRaises(ValueError):
                load_functions(source, "0" * 64, {})

    def test_top_level_actions_excluded_and_function_unchanged(self):
        with tempfile.TemporaryDirectory() as folder:
            source = Path(folder) / "producer.py"
            source.write_text("raise AssertionError('must not run')\ndef value(): return 3\n")
            ns = {}
            load_functions(source, hashlib.sha256(source.read_bytes()).hexdigest(), ns)
            self.assertEqual(ns["value"](), 3)


if __name__ == "__main__":
    unittest.main()
