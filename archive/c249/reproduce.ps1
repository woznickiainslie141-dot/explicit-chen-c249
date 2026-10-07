param(
    [string]$CorePython = "python",
    [string]$PdfPython = "",
    [switch]$SkipPdfQa
)

$ErrorActionPreference = "Stop"
Set-Location -LiteralPath $PSScriptRoot
if ([string]::IsNullOrWhiteSpace($PdfPython)) {
    $PdfPython = $CorePython
}

function Invoke-Checked {
    param(
        [Parameter(Mandatory = $true)][string]$Program,
        [Parameter(ValueFromRemainingArguments = $true)][string[]]$Arguments
    )
    & $Program @Arguments
    if ($LASTEXITCODE -ne 0) {
        throw "Command failed with exit code ${LASTEXITCODE}: $Program $($Arguments -join ' ')"
    }
}

Write-Host "[1/6] Checking the unpacked bundle and immutable-input hashes"
Invoke-Checked $CorePython "verify_bundle.py" "--require-generated"

Write-Host "[2/6] Checking Python proof dependencies"
Invoke-Checked $CorePython "-c" "import flint; print('python-flint', flint.__version__)"

Write-Host "[3/6] Recomputing all interval and structural certificates"
Invoke-Checked $CorePython "c249_verify.py" "--check-tex"

if (-not (Get-Command "pdflatex" -ErrorAction SilentlyContinue)) {
    throw "pdflatex is required for full document reproduction but was not found on PATH"
}
Write-Host "[4/6] Rebuilding the proof PDF with three fixed TeX passes"
$latexArgs = @(
    "-interaction=nonstopmode",
    "-halt-on-error",
    "-file-line-error",
    "-output-directory=output/pdf",
    "explicit_chen_c249_final.tex"
)
1..3 | ForEach-Object { Invoke-Checked "pdflatex" @latexArgs }

$logPath = Join-Path $PSScriptRoot "output/pdf/explicit_chen_c249_final.log"
$logText = Get-Content -LiteralPath $logPath -Raw
$fatalWarnings = @(
    "There were undefined references",
    "LaTeX Warning: Reference .* undefined",
    "LaTeX Warning: Citation .* undefined",
    "Rerun to get cross-references right"
)
foreach ($pattern in $fatalWarnings) {
    if ($logText -match $pattern) {
        throw "Unresolved TeX cross-reference warning matched: $pattern"
    }
}

if ($SkipPdfQa) {
    Write-Host "[5/6] PDF render QA skipped by explicit switch"
} else {
    if (-not (Get-Command "pdftoppm" -ErrorAction SilentlyContinue)) {
        throw "pdftoppm is required for full PDF rendering QA but was not found on PATH"
    }
    Write-Host "[5/6] Rendering all PDF pages and running text/layout preflight"
    Invoke-Checked $PdfPython "-c" "import pypdf, PIL; print('pypdf', pypdf.__version__, 'Pillow', PIL.__version__)"
    Invoke-Checked $PdfPython "c249_pdf_qa.py"
}

Write-Host "[6/6] Rechecking the regenerated ledgers and final artifact"
Invoke-Checked $CorePython "verify_bundle.py" "--require-generated"
Write-Host "FULL C=24.9 REPRODUCTION PASSED"

