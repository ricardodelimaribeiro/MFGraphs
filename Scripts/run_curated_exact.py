#!/usr/bin/env python3
"""Sequential fresh-kernel exact benchmarks, with durable, verified artifacts.

Examples:
  python3 Scripts/run_curated_exact.py --cases diamond,triangle --timeout 5
  python3 Scripts/run_curated_exact.py --cases grid4 --methods dnf,lex --timeout 60
  python3 Scripts/run_curated_exact.py --compare-root /path/to/baseline/source --repetitions 2

No inventory sweep, numeric fallback, or solution-cache retrieval is performed.
The worker uses an identical small generic warm-up in every fresh kernel.
"""
import argparse
import csv
import datetime as dt
import fcntl
import hashlib
import json
import os
from pathlib import Path
import platform
import signal
import subprocess
import sys
import time
import uuid
import tarfile

ROOT = Path(__file__).resolve().parents[1]
CASES = ("competing-exits", "diamond", "triangle", "braess-split", "braess-congest",
         "camilli-simple", "grid3", "grid4", "jamarat", "case23")
METHODS = ("dnf", "reduce", "lex", "block-edge", "linear-net")
RUNNER_FILES = ("Scripts/run_curated_exact.py", "Scripts/RunCuratedExact.wls",
                "Scripts/CuratedExactScenarios.wls", "Scripts/ExactArtifactIntegrity.wls",
                "MFGraphs/validationTools.wl")


def command_output(args):
    result = subprocess.run(args, cwd=ROOT, text=True, capture_output=True)
    return result.stdout.strip() if result.returncode == 0 else result.stderr.strip()


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def freeze_files(root, paths, destination):
    """Hash and archive the same read-only bytes that workers will load."""
    hashes = {}
    for path in paths:
        relative = path.relative_to(root)
        target = destination / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(path.read_bytes())
        target.chmod(0o444)
        hashes[str(relative)] = digest(target)
    return hashes


class ArtifactIntegrityError(ValueError):
    pass


def read_sample_artifacts(sample, case, method):
    completion = json.loads((sample / "completed.json").read_text())
    summary = json.loads((sample / "summary.json").read_text())
    if not isinstance(completion, dict) or not isinstance(summary, dict):
        raise ValueError("Completion and summary must be JSON objects")
    if completion.get("Complete") is not True:
        raise ArtifactIntegrityError("Worker did not certify completed artifacts")
    checksum = digest(sample / "result.wxf")
    if completion.get("ArtifactSHA256") != checksum:
        raise ArtifactIntegrityError("Result checksum does not match the completion record")
    if summary.get("ArtifactIntegrity") != "ExactContentMatch":
        raise ArtifactIntegrityError("Worker did not certify exact content integrity")
    if summary.get("Case") != case or summary.get("Method") != method:
        raise ValueError("Summary identifies a different case or method")
    for field in ("ResultType", "Soundness", "Completeness"):
        if not isinstance(summary.get(field), str):
            raise ValueError(f"Summary has no valid {field}")
    return summary, checksum


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cases", default=",".join(CASES))
    parser.add_argument("--methods", default="dnf")
    parser.add_argument("--timeout", type=float, default=10)
    parser.add_argument("--validation-timeout", type=float, default=10)
    parser.add_argument("--repetitions", type=int, default=1)
    parser.add_argument("--source-root", type=Path, default=ROOT)
    parser.add_argument("--compare-root", type=Path)
    parser.add_argument("--baseline-methods", help="Methods available in the preserved source (defaults to --methods)")
    parser.add_argument("--tag", default="curated")
    args = parser.parse_args()
    cases, methods = args.cases.split(","), args.methods.split(",")
    if set(cases) - set(CASES) or set(methods) - set(METHODS):
        parser.error("Unknown case or method; inspect --help and CuratedExactScenarios.wls")
    if not 0 < args.timeout <= 300 or not 0 < args.validation_timeout <= 300 or args.repetitions < 1:
        parser.error("Stage limits must be positive and at most 300 seconds; repetitions >= 1")
    # Persistent advisory lock: a crashed process releases it automatically.
    lock_path = ROOT / "Results" / "exact-foundation" / ".runner.lock"
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    with lock_path.open("a") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            parser.error("Another curated run is active; timing experiments must run sequentially")
        stamp = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        run = lock_path.parent / f"{stamp}-{args.tag}-{uuid.uuid4().hex[:8]}"
        run.mkdir(exist_ok=False)
        sources = [("current", args.source_root.resolve())]
        if args.compare_root:
            sources.insert(0, ("baseline", args.compare_root.resolve()))
        baseline_methods = args.baseline_methods.split(",") if args.baseline_methods else methods
        if set(baseline_methods) - set(METHODS):
            parser.error("Unknown baseline method")
        variants = [(label, root, method) for label, root in sources
                    for method in (baseline_methods if label == "baseline" else methods)]
        source_hashes = {}
        execution_sources = {}
        for label, root in sources:
            destination = run / "source" / label
            paths = sorted((root / "MFGraphs").glob("*.wl"))
            paths += [root / "MFGraphs/Kernel/init.m"]
            if (root / "PacletInfo.m").exists():
                paths.append(root / "PacletInfo.m")
            source_hashes[label] = freeze_files(root, paths, destination)
            execution_sources[label] = str(destination)
        runner_root = run / "source" / "runner"
        runner_hashes = freeze_files(ROOT, [ROOT / name for name in RUNNER_FILES], runner_root)
        manifest = {
            "started_utc": stamp, "command": sys.argv,
            "commit": command_output(["git", "rev-parse", "HEAD"]),
            "status": command_output(["git", "status", "--short"]),
            "platform": platform.platform(), "machine": platform.machine(),
            "cpu": command_output(["sysctl", "-n", "machdep.cpu.brand_string"]),
            "memory_bytes": command_output(["sysctl", "-n", "hw.memsize"]),
            "warmup": "Fresh kernel per sample; one 2-vertex DNF warm-up, max 5 seconds",
            "solve_cache": "Bypassed; direct system solver calls",
            "skip_rule": "After a solve timeout, record SKIPPED_TIMEOUT for later repetitions of that exact case/source/method",
            "order": "Variants rotate by repetition within each case; fresh kernel per repetition",
            "sources": source_hashes, "runner_sources": runner_hashes,
            "execution_sources": execution_sources, "runner_root": str(runner_root),
            "options": {k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()},
        }
        (run / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
        (run / "worktree.patch").write_text(command_output(["git", "diff", "--binary"]) + "\n")
        with tarfile.open(run / "sources.tar.gz", "w:gz") as archive:
            for label, hashes in dict(source_hashes, runner=runner_hashes).items():
                for name in hashes:
                    archive.add(run / "source" / label / name, arcname=f"{label}/{name}")
        print(f"RUN_DIRECTORY={run}", flush=True)
        rows, skipped = [], set()
        for case in cases:
            for rep in range(1, args.repetitions + 1):
                offset = (rep - 1) % len(variants)
                for label, source, method in variants[offset:] + variants[:offset]:
                    key = (case, label, method)
                    row = {"Case": case, "Source": label, "Method": method, "Repetition": rep,
                           "Position": len(rows) + 1}
                    if key in skipped:
                        row["RunStatus"] = "SKIPPED_TIMEOUT"
                    else:
                        sample = run / f"{len(rows)+1:03d}-{case}-{label}-{method}-r{rep}"
                        sample.mkdir()
                        cmd = ["wolframscript", "-file", str(runner_root / "Scripts/RunCuratedExact.wls"),
                               str(run / "manifest.json"), str(sample), execution_sources[label], case, method,
                               str(args.timeout), str(args.validation_timeout)]
                        row["Directory"] = sample.name
                        (sample / "command.json").write_text(json.dumps(cmd, indent=2) + "\n")
                        started = time.monotonic()
                        with (sample / "stdout.log").open("w") as output:
                            try:
                                # Outside-kernel watchdog covers startup/export and interrupts
                                # that the kernel's TimeConstrained might fail to handle.
                                proc = subprocess.Popen(cmd, cwd=ROOT, stdout=output, stderr=subprocess.STDOUT,
                                                        start_new_session=True)
                                proc.wait(timeout=args.timeout + 2 * args.validation_timeout + 60)
                                row["ExitCode"] = proc.returncode
                                if proc.returncode != 0:
                                    row["RunStatus"] = "EXIT_FAILURE"
                                    row["Error"] = f"Worker exited with status {proc.returncode}"
                            except subprocess.TimeoutExpired:
                                os.killpg(proc.pid, signal.SIGKILL)
                                proc.wait()
                                row["ExitCode"] = None
                                row["RunStatus"] = "PROCESS_TIMEOUT"
                                row["Error"] = "Worker exceeded the outer process time limit"
                            except OSError as error:
                                row["ExitCode"] = None
                                row["RunStatus"] = "PROCESS_FAILURE"
                                row["Error"] = f"{type(error).__name__}: {error}"
                        row["ProcessSeconds"] = time.monotonic() - started
                        try:
                            summary, checksum = read_sample_artifacts(sample, case, method)
                            # Worker metadata must not overwrite controller status/identity.
                            row.update({key: value for key, value in summary.items() if key not in row})
                            row.setdefault("RunStatus", "COMPLETE")
                            row["ArtifactSHA256"] = checksum
                        except (OSError, ValueError, UnicodeError) as error:
                            if isinstance(error, ArtifactIntegrityError):
                                artifact_status = "INTEGRITY_FAILURE"
                            elif isinstance(error, FileNotFoundError):
                                artifact_status = "INCOMPLETE_ARTIFACT"
                            else:
                                artifact_status = "INVALID_ARTIFACT"
                            row.setdefault("RunStatus", artifact_status)
                            row["ArtifactStatus"] = artifact_status
                            row["Error"] = "; ".join(filter(None, [row.get("Error"),
                                f"{type(error).__name__}: {error}"]))
                        if row.get("ResultType") == "Timeout" or row["RunStatus"] == "PROCESS_TIMEOUT":
                            skipped.add(key)
                    rows.append(row)
                    # Durable after every sample, including skipped/failed ones.
                    (run / "results.json").write_text(json.dumps(rows, indent=2) + "\n")
                    headers = list(dict.fromkeys(k for r in rows for k in r))
                    with (run / "results.csv").open("w", newline="") as stream:
                        writer = csv.DictWriter(stream, fieldnames=headers)
                        writer.writeheader(); writer.writerows(rows)
                    print(f"{case} {label} {method} rep={rep}: {row['RunStatus']} "
                          f"{row.get('ResultType', '')} solve={row.get('SolveSeconds', '-')} "
                          f"sound={row.get('Soundness', '-')} complete={row.get('Completeness', '-')}", flush=True)
        integrity = {str(p.relative_to(run)): digest(p) for p in sorted(run.rglob("*")) if p.is_file()}
        (run / "sha256.json").write_text(json.dumps(integrity, indent=2) + "\n")
        successful_run = all(r["RunStatus"] in ("COMPLETE", "SKIPPED_TIMEOUT") for r in rows)
        (run / "run-completed.json").write_text(json.dumps({"Complete": successful_run,
            "ExpectedRows": len(cases) * len(variants) * args.repetitions, "ActualRows": len(rows)}, indent=2) + "\n")
        return 0 if successful_run else 1


if __name__ == "__main__":
    sys.exit(main())
