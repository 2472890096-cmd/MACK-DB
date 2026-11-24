#!/usr/bin/env python3
"""Small path/source-safety tests; no scientific data or analysis is executed."""
import ast
import hashlib
import importlib.util
import json
import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from run_original_statistics import PINS, default_source, digest, install_io_rebase, load_functions
from validate_statistics import unique_index
from fixed_order_recovery import archived_steps


class AdapterTests(unittest.TestCase):
    def test_bundled_sources_and_helpers_have_accepted_hashes(self):
        root = Path(__file__).resolve().parents[2]
        for mode in ["b", "c", "d"]:
            source = default_source(mode)
            self.assertTrue(source.is_relative_to(root))
            self.assertEqual(digest(source), PINS[mode])
            if mode != "b":
                self.assertEqual(digest(source.with_name("common.py")), PINS[mode + "_common"])

    def test_matched_subgroup_data_root_survives_source_relocation(self):
        """Moving scientific source must not redirect data IO into its code folder."""
        root = Path(__file__).resolve().parents[2]
        spec = importlib.util.spec_from_file_location("matched_subgroup_paths", root / "lib/mackdb_paths.py")
        paths = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(paths)
        sources = default_source("b").parent
        self.assertEqual(sources, root / "analysis/paper_network_comparison/original/matched_subgroups")
        with tempfile.TemporaryDirectory() as folder:
            data = Path(folder) / "separate-inputs"
            config = Path(folder) / "paths.json"
            config.write_text(json.dumps({"stage5b": str(data)}))
            with patch.dict(os.environ, {"MACKDB_PATHS_CONFIG": str(config)}):
                for name in ["stage5b2_analysis.py", "build_reference.py",
                             "independent_recalculation.py", "run_degree_null.py"]:
                    source = sources / name
                    tree = ast.parse(source.read_text())
                    value = next(n.value for n in tree.body if isinstance(n, ast.Assign)
                                 and any(isinstance(t, ast.Name) and t.id == "ROOT" for t in n.targets))
                    expression = compile(ast.Expression(value), str(source), "eval")
                    for code_path in [source, Path(folder) / "relocated-code" / name]:
                        actual = eval(expression, {"resolve_path": paths.resolve_path,
                                                   "Path": Path, "__file__": str(code_path)})
                        self.assertEqual(actual, data / "04_network_comparison")

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

    def test_broad_project_rebase_preserves_exact_source_but_relocates_data(self):
        with tempfile.TemporaryDirectory() as folder:
            project = Path(folder) / "project"
            source = project / "checkout/analysis/source.py"
            source.parent.mkdir(parents=True)
            source.write_text("def value(): return 3\n")
            seen = []
            ns = {"ROOT": Path(folder) / "output", "__file__": str(source),
                  "register": lambda path: seen.append(Path(path))}
            local_data = Path(folder) / "relocated-data"
            translate = install_io_rebase(ns, {str(project): str(local_data)})
            ns["register"](source)
            self.assertEqual(seen, [source])
            self.assertTrue(seen[0].is_file())
            self.assertEqual(translate(project / "cohort/abundance.tsv"),
                             str(local_data / "cohort/abundance.tsv"))
            # Even a neighboring file in the checkout is not exempted merely
            # because it is under the repository directory.
            self.assertEqual(translate(source.with_name("saved_data.tsv")),
                             str(local_data / "checkout/analysis/saved_data.tsv"))

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
