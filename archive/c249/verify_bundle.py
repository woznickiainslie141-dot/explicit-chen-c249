#!/usr/bin/env python3
"""Verify that the c=24.9 bundle is complete and internally consistent.

This checker uses only the Python standard library.  It verifies every immutable
input against source_manifest.sha256, checks all required transitive files, and
cross-checks the generated proof ledgers against the current source files.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent

REQUIRED = (
    "README.md",
    "README_REPRODUCE.md",
    "CITATION.cff",
    "NOTICE.md",
    ".gitignore",
    "VALIDATION_REPORT.md",
    "reproduce.ps1",
    "install_offline.ps1",
    "verify_bundle.py",
    "source_manifest.sha256",
    "dependency_graph.json",
    "c249_requirements.txt",
    "c249_verify.py",
    "c249_interval_certificate.py",
    "c249_structural_certificate.py",
    "c26_interval_certificate.py",
    "c249_diagnostic.py",
    "c249_pdf_qa.py",
    "explicit_chen_c249_final.tex",
    "final_audit.md",
    "c249_status.md",
    "bjs_source/Explicit_Chen_-_New.tex",
    "wu_source/ChenDoubleSieve1_Paper.tex",
    "audit_inputs/reply.md",
    "audit_inputs/explicit_chen_c2336_final.tex",
    "source_archives/bjs2207.09452v6.tar",
    "source_archives/wu0705.1652.tar",
)

GENERATED = (
    "c249_final_ledger.json",
    "c249_interval_ledger.json",
    "c249_structural_ledger.json",
    "output/pdf/explicit_chen_c249_final.pdf",
)

EXPECTED_MAIN = "[0.023536843646471438769041841946393 +/- 1.22e-31]"
EXPECTED_BUDGET = "[0.012005000000000000000000000000000 +/- 3e-38]"
EXPECTED_NET = "[0.011531843646471438769041841946393 +/- 1.22e-31]"


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def fail(message: str) -> None:
    raise SystemExit(f"BUNDLE CHECK FAILED: {message}")


def checked_path(relative: str) -> Path:
    path = (ROOT / relative).resolve()
    try:
        path.relative_to(ROOT.resolve())
    except ValueError:
        fail(f"manifest path escapes bundle root: {relative}")
    return path


def read_manifest() -> dict[str, str]:
    result: dict[str, str] = {}
    manifest = ROOT / "source_manifest.sha256"
    for number, raw in enumerate(manifest.read_text(encoding="utf-8").splitlines(), 1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        pieces = line.split(maxsplit=1)
        if len(pieces) != 2 or len(pieces[0]) != 64:
            fail(f"malformed manifest line {number}")
        relative = pieces[1].lstrip("*").replace("\\", "/")
        if relative in result:
            fail(f"duplicate manifest entry: {relative}")
        result[relative] = pieces[0].lower()
    return result


def check_required(require_generated: bool) -> None:
    names = REQUIRED + (GENERATED if require_generated else ())
    missing = [name for name in names if not checked_path(name).is_file()]
    if missing:
        fail("missing required files: " + ", ".join(missing))
    empty = [name for name in names if checked_path(name).stat().st_size == 0]
    if empty:
        fail("empty required files: " + ", ".join(empty))


def check_manifest() -> None:
    entries = read_manifest()
    for relative, expected in entries.items():
        path = checked_path(relative)
        if not path.is_file():
            fail(f"manifested file is missing: {relative}")
        actual = digest(path)
        if actual != expected:
            fail(f"SHA-256 mismatch for {relative}: {actual} != {expected}")
    for relative in REQUIRED:
        if relative != "source_manifest.sha256" and relative not in entries:
            fail(f"immutable required file is not manifested: {relative}")
    print(f"Immutable inputs: {len(entries)} SHA-256 entries verified")


def check_dependency_interfaces() -> None:
    verifier = (ROOT / "c249_verify.py").read_text(encoding="utf-8")
    interval = (ROOT / "c249_interval_certificate.py").read_text(encoding="utf-8")
    expected_fragments = {
        "c249_verify.py": (
            "import c249_interval_certificate",
            "import c249_structural_certificate",
            '"bjs_source/Explicit_Chen_-_New.tex"',
            '"wu_source/ChenDoubleSieve1_Paper.tex"',
        ),
        "c249_interval_certificate.py": (
            "import c26_interval_certificate as kernel",
        ),
    }
    bodies = {
        "c249_verify.py": verifier,
        "c249_interval_certificate.py": interval,
    }
    for name, fragments in expected_fragments.items():
        for fragment in fragments:
            if fragment not in bodies[name]:
                fail(f"expected dependency interface changed in {name}: {fragment}")
    print("Dependency interfaces: direct and transitive local imports resolved")


def check_wheelhouse() -> None:
    wheels = [path.name.lower() for path in (ROOT / "wheelhouse").glob("*.whl")]
    expected = {
        "python-flint==0.9.0": ("python_flint-0.9.0-",),
        "pypdf==6.10.0": ("pypdf-6.10.0-",),
        "Pillow==12.3.0": ("pillow-12.3.0-",),
    }
    for requirement, prefixes in expected.items():
        if not any(name.startswith(tuple(prefix.lower() for prefix in prefixes)) for name in wheels):
            fail(f"offline wheel missing for {requirement}")
    print(f"Offline Python dependencies: {len(wheels)} wheel files present")


def check_ledgers() -> None:
    final = json.loads((ROOT / "c249_final_ledger.json").read_text(encoding="utf-8"))
    interval = json.loads((ROOT / "c249_interval_ledger.json").read_text(encoding="utf-8"))
    structural = json.loads((ROOT / "c249_structural_ledger.json").read_text(encoding="utf-8"))

    expected = {
        "main_before_charges": EXPECTED_MAIN,
        "sum_of_finite_budgets": EXPECTED_BUDGET,
        "net_arithmetic_enclosure": EXPECTED_NET,
    }
    for key, value in expected.items():
        if final.get(key) != value:
            fail(f"unexpected {key}: {final.get(key)!r}")
    if final.get("target") != "all even N with log log N >= 24.9, positivity only":
        fail("unexpected theorem target in final ledger")
    if final.get("tex_table_checked") is not True:
        fail("TeX table was not checked")
    if final.get("rigorous_positive_floor") != "0.0108":
        fail("unexpected rigorous positive floor")
    if len(final.get("rows", [])) != 11:
        fail("final signed-row ledger does not contain 11 rows")
    if len(final.get("normalized_finite_budgets", {})) != 7:
        fail("final budget ledger does not contain 7 budgets")
    if len(interval.get("integrals", [])) != 34:
        fail("interval ledger does not contain 34 Arb/ACB integrals")
    if len(structural.get("passed_assertions", [])) != 136:
        fail("structural ledger does not contain 136 passed assertions")

    tracked = final.get("sha256", {})
    if not isinstance(tracked, dict) or not tracked:
        fail("final ledger has no tracked source hashes")
    for relative, expected_hash in tracked.items():
        actual_hash = digest(checked_path(relative))
        if actual_hash != expected_hash:
            fail(f"final-ledger source mismatch for {relative}")

    pdf = ROOT / "output/pdf/explicit_chen_c249_final.pdf"
    if pdf.is_file():
        data_start = pdf.read_bytes()[:8]
        if not data_start.startswith(b"%PDF-") or pdf.stat().st_size < 100_000:
            fail("generated PDF is missing, truncated, or malformed")

    print("Generated ledgers: 34 integrals, 136 structural assertions, 11 rows, 7 budgets")
    print(f"Main   = {EXPECTED_MAIN}")
    print(f"Budget = {EXPECTED_BUDGET}")
    print(f"Net    = {EXPECTED_NET}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--require-generated",
        action="store_true",
        help="also require and cross-check all generated ledgers and the final PDF",
    )
    args = parser.parse_args()
    check_required(args.require_generated)
    check_manifest()
    check_dependency_interfaces()
    check_wheelhouse()
    if args.require_generated:
        check_ledgers()
    print("BUNDLE CHECK PASSED")


if __name__ == "__main__":
    main()

