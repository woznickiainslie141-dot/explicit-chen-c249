"""Build the current manuscript using the installed TeX Live compiler.

Run: python conditional_chen_build_pdf.py
This project fallback does not repair Codex's sandboxed live preview.
"""

from pathlib import Path
import hashlib
import json
import re
import shutil
import subprocess


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "conditional_chen_grh_bridge.tex"
OUTPUT = ROOT / "output" / "pdf"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    compiler = shutil.which("pdflatex")
    if not compiler:
        raise SystemExit("pdflatex was not found; the existing TeX Live installation is required.")
    OUTPUT.mkdir(parents=True, exist_ok=True)
    version = subprocess.run([compiler, "--version"], check=True, capture_output=True,
                             text=True, errors="replace").stdout.splitlines()[0]
    source_hash = digest(SOURCE)
    command = [compiler, "-interaction=nonstopmode", "-halt-on-error",
               "-file-line-error", "-no-shell-escape", "-synctex=1",
               "-output-directory=" + str(OUTPUT), SOURCE.name]
    log_path = OUTPUT / (SOURCE.stem + ".log")
    for pass_number in range(1, 5):
        result = subprocess.run(command, cwd=ROOT, capture_output=True, timeout=120)
        transcript = result.stdout + result.stderr
        (OUTPUT / f"conditional_chen_build_pass{pass_number}.stdout.log").write_bytes(transcript)
        if result.returncode:
            raise SystemExit(f"TeX pass {pass_number} failed. See {log_path}")
        log = log_path.read_text(encoding="utf-8", errors="replace")
        rerun = any(term in log for term in (
            "There were undefined references", "Rerun to get cross-references right",
            "Rerun to get outlines right"))
        if not rerun:
            break
    else:
        raise SystemExit("References did not stabilize after four TeX passes.")
    checks = {
        "undefined_references_or_citations": len(re.findall(r"(?:Reference|Citation) `[^\n]+undefined", log)),
        "overfull_boxes": len(re.findall(r"Overfull \\[hv]box", log)),
        "duplicate_destinations": log.count("duplicate ignored"),
        "latex_errors": len(re.findall(r"^!|LaTeX Error:", log, re.MULTILINE)),
    }
    pdf = OUTPUT / (SOURCE.stem + ".pdf")
    assert digest(SOURCE) == source_hash, "Source changed while compilation was running"
    assert pdf.is_file() and pdf.stat().st_size > 0
    report = {
        "status": "PASS" if not any(checks.values()) else "WARNINGS",
        "scope": "TeX compilation and log checks; mathematical correctness is not certified",
        "source": SOURCE.name, "source_sha256": source_hash,
        "compiler": compiler, "compiler_version": version,
        "command": command, "passes": pass_number,
        "pdf": str(pdf.relative_to(ROOT)), "pdf_sha256": digest(pdf),
        "checks": checks,
        "codex_builtin_preview": "still fails: Unable to find standard directories for platform",
    }
    (ROOT / "conditional_chen_compile_certificate.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps({"status": report["status"], "passes": pass_number,
                      "pdf": report["pdf"], "checks": checks}, ensure_ascii=False))
    if any(checks.values()):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
