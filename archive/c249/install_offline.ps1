param(
    [string]$Python = "python"
)

$ErrorActionPreference = "Stop"
Set-Location -LiteralPath $PSScriptRoot
& $Python -m pip install --no-index --find-links (Join-Path $PSScriptRoot "wheelhouse") -r (Join-Path $PSScriptRoot "c249_requirements.txt")
if ($LASTEXITCODE -ne 0) {
    throw "Offline dependency installation failed with exit code $LASTEXITCODE"
}
Write-Host "OFFLINE PYTHON DEPENDENCIES INSTALLED"

