"""Controller integration checks using a fake external wolframscript executable.

Run with: python3 -m unittest discover -s Scripts/tests -v
No Wolfram kernel or mathematical benchmark is launched.
"""
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import tempfile
import unittest


REPO = Path(__file__).resolve().parents[2]
RUNNER_FILES = (
    "Scripts/run_curated_exact.py", "Scripts/RunCuratedExact.wls",
    "Scripts/CuratedExactScenarios.wls", "Scripts/ExactArtifactIntegrity.wls",
    "MFGraphs/validationTools.wl",
)

FAKE_WORKER = r'''#!/usr/bin/env python3
import hashlib, json, os, sys
from pathlib import Path
worker, manifest_path, output, source = map(Path, sys.argv[2:6])
case, method = sys.argv[6:8]
live = Path(os.environ["FAKE_LIVE_ROOT"])
mode = os.environ.get("FAKE_MODE", "complete")
count_file = live / "calls.json"
count = json.loads(count_file.read_text()) if count_file.exists() else 0
count_file.write_text(json.dumps(count + 1))
runner = worker.parent.parent
paths = {"package": source / "MFGraphs" / "solver.wl",
         "worker": worker,
         "scenario": runner / "Scripts" / "CuratedExactScenarios.wls",
         "integrity": runner / "Scripts" / "ExactArtifactIntegrity.wls",
         "validator": runner / "MFGraphs" / "validationTools.wl"}
observed = {name: hashlib.sha256(path.read_bytes()).hexdigest() for name, path in paths.items()}
artifact = output / "result.wxf"
artifact.write_bytes(b"fake exact artifact, intentionally not a WXF expression")
summary = {"Case": case, "Method": method, "ResultType": "ExactDetermined",
           "Soundness": "Proved", "Completeness": "Proved",
           "ArtifactIntegrity": "ExactContentMatch", "ObservedSources": observed}
(output / "summary.json").write_text(json.dumps(summary))
done = {"Complete": True, "ArtifactSHA256": hashlib.sha256(artifact.read_bytes()).hexdigest()}
(output / "completed.json").write_text(json.dumps(done))
if mode == "mutate_live" and count == 0:
    for path in [live / "MFGraphs" / "solver.wl", live / "baseline" / "MFGraphs" / "solver.wl",
                 live / "Scripts" / "RunCuratedExact.wls", live / "Scripts" / "CuratedExactScenarios.wls",
                 live / "Scripts" / "ExactArtifactIntegrity.wls", live / "MFGraphs" / "validationTools.wl"]:
        path.write_text("changed after the first sample")
if count == 0:
    if mode == "solve_timeout":
        summary.update(ResultType="Timeout", Soundness="NotChecked", Completeness="NotChecked")
        (output / "summary.json").write_text(json.dumps(summary))
    elif mode == "bad_checksum":
        done["ArtifactSHA256"] = "0" * 64
        (output / "completed.json").write_text(json.dumps(done))
    elif mode == "false_completion":
        done["Complete"] = False
        (output / "completed.json").write_text(json.dumps(done))
    elif mode.startswith("missing_"):
        (output / {"completion": "completed.json", "summary": "summary.json",
                   "artifact": "result.wxf"}[mode.removeprefix("missing_")]).unlink()
    elif mode in ("invalid_completion", "invalid_summary"):
        (output / ("completed.json" if mode.endswith("completion") else "summary.json")).write_text("{")
    elif mode in ("completion_not_object", "summary_not_object"):
        (output / ("completed.json" if mode.startswith("completion") else "summary.json")).write_text("[]")
    elif mode == "corrupt_artifact":
        artifact.write_bytes(b"artifact bytes changed after completion")
(output / "expected-raw-files.json").write_text(json.dumps(
    {p.name: p.read_bytes().hex() for p in output.iterdir() if p.is_file()}))
if mode == "nonzero_exit" and count == 0:
    sys.exit(7)
'''


class CuratedRunnerTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="curated-runner-test-")
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name).resolve()
        for relative in RUNNER_FILES:
            target = self.root / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(REPO / relative, target)
        for package, value in ((self.root, "current solver"), (self.root / "baseline", "baseline solver")):
            (package / "MFGraphs" / "Kernel").mkdir(parents=True, exist_ok=True)
            (package / "MFGraphs" / "solver.wl").write_text(value)
            (package / "MFGraphs" / "MFGraphs.wl").write_text("fake package loader")
            (package / "MFGraphs" / "Kernel" / "init.m").write_text("fake kernel entrypoint")
        executable = self.root / "bin" / "wolframscript"
        executable.parent.mkdir()
        executable.write_text(FAKE_WORKER)
        executable.chmod(0o755)
        self.environment = dict(os.environ, PATH=str(executable.parent) + os.pathsep + os.environ["PATH"],
                                FAKE_LIVE_ROOT=str(self.root), PYTHONDONTWRITEBYTECODE="1")

    def run_controller(self, mode="complete", extra=()):
        results = self.root / "Results" / "exact-foundation"
        previous = set(results.glob("*-curated-*"))
        (self.root / "calls.json").write_text("0")
        result = subprocess.run(
            [sys.executable, "Scripts/run_curated_exact.py", "--cases", "diamond", *extra],
            cwd=self.root, env=dict(self.environment, FAKE_MODE=mode), text=True,
            capture_output=True, timeout=20)
        runs = list(set(results.glob("*-curated-*")) - previous)
        self.assertEqual(len(runs), 1, result.stdout + result.stderr)
        return result, runs[0]

    def test_live_edits_cannot_change_any_executed_source_in_later_samples(self):
        expected = {
            "worker": hashlib.sha256((self.root / "Scripts/RunCuratedExact.wls").read_bytes()).hexdigest(),
            "scenario": hashlib.sha256((self.root / "Scripts/CuratedExactScenarios.wls").read_bytes()).hexdigest(),
            "integrity": hashlib.sha256((self.root / "Scripts/ExactArtifactIntegrity.wls").read_bytes()).hexdigest(),
            "validator": hashlib.sha256((self.root / "MFGraphs/validationTools.wl").read_bytes()).hexdigest(),
        }
        packages = {"baseline": hashlib.sha256(b"baseline solver").hexdigest(),
                    "current": hashlib.sha256(b"current solver").hexdigest()}
        result, run = self.run_controller("mutate_live", ["--compare-root", str(self.root / "baseline"),
                                                         "--repetitions", "2"])
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        rows = json.loads((run / "results.json").read_text())
        self.assertEqual(len(rows), 4)
        for row in rows:
            self.assertEqual(row["ObservedSources"], dict(expected, package=packages[row["Source"]]))
            command = json.loads((run / row["Directory"] / "command.json").read_text())
            self.assertTrue(Path(command[2]).is_relative_to(run))
            self.assertTrue(Path(command[5]).is_relative_to(run))
        manifest = json.loads((run / "manifest.json").read_text())
        with tarfile.open(run / "sources.tar.gz") as archive:
            for label, files in dict(manifest["sources"], runner=manifest["runner_sources"]).items():
                for relative, expected_hash in files.items():
                    self.assertEqual(hashlib.sha256(archive.extractfile(f"{label}/{relative}").read()).hexdigest(),
                                     expected_hash)

    def test_every_artifact_failure_has_a_durable_row_and_final_failure_marker(self):
        modes = ("bad_checksum", "false_completion", "missing_completion", "missing_summary",
                 "missing_artifact", "invalid_completion", "invalid_summary",
                 "completion_not_object", "summary_not_object", "corrupt_artifact", "nonzero_exit")
        for mode in modes:
            with self.subTest(mode=mode):
                result, run = self.run_controller(mode, ["--repetitions", "2"])
                self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
                rows = json.loads((run / "results.json").read_text())
                self.assertEqual(len(rows), 2)
                self.assertNotEqual(rows[0]["RunStatus"], "COMPLETE")
                self.assertTrue(rows[0]["Error"])
                self.assertEqual(rows[1]["RunStatus"], "COMPLETE")
                marker = json.loads((run / "run-completed.json").read_text())
                self.assertEqual(marker, {"Complete": False, "ExpectedRows": 2, "ActualRows": 2})
                with (run / "results.csv").open() as stream:
                    self.assertEqual(len(list(csv.DictReader(stream))), 2)
                for row in rows:
                    sample = run / row["Directory"]
                    for name, raw_hex in json.loads((sample / "expected-raw-files.json").read_text()).items():
                        self.assertEqual((sample / name).read_bytes().hex(), raw_hex)

    def test_worker_launch_failure_is_recorded_for_every_requested_sample(self):
        (self.root / "bin/wolframscript").write_text("#!/nonexistent-curated-test-interpreter\n")
        result, run = self.run_controller(extra=["--repetitions", "2"])
        self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
        rows = json.loads((run / "results.json").read_text())
        self.assertEqual([row["RunStatus"] for row in rows], ["PROCESS_FAILURE", "PROCESS_FAILURE"])
        self.assertTrue(all("FileNotFoundError" in row["Error"] for row in rows))
        self.assertEqual(json.loads((run / "run-completed.json").read_text()),
                         {"Complete": False, "ExpectedRows": 2, "ActualRows": 2})

    def test_solve_timeout_preserves_the_existing_explicit_skip_policy(self):
        result, run = self.run_controller("solve_timeout", ["--repetitions", "3"])
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        rows = json.loads((run / "results.json").read_text())
        self.assertEqual([row["RunStatus"] for row in rows],
                         ["COMPLETE", "SKIPPED_TIMEOUT", "SKIPPED_TIMEOUT"])
        self.assertEqual(rows[0]["ResultType"], "Timeout")
        self.assertEqual(json.loads((self.root / "calls.json").read_text()), 1)
        self.assertEqual(json.loads((run / "run-completed.json").read_text()),
                         {"Complete": True, "ExpectedRows": 3, "ActualRows": 3})


if __name__ == "__main__":
    unittest.main()
