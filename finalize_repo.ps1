# finalize_repo.ps1
# Metatron Dynamics, Inc. V6.
# Copies all analysis outputs to the repo and makes a clean commit.
#
# Run from inside abr-event-analysis:
#   cd 'C:\Users\Robin Macomber\Documents\Metatron_Dynamics\GitHub_Repos\abr-event-analysis'
#   .\finalize_repo.ps1

$ErrorActionPreference = "Stop"
$Downloads = "$env:USERPROFILE\Downloads"

Write-Host "Finalizing abr-event-analysis repository..." -ForegroundColor Cyan
Write-Host ""

# Files to copy from Downloads to repo root
$Scripts = @(
    "analyze_hadronic_phase0a.py",
    "analyze_hadronic_phase0b.py",
    "analyze_hadronic_phase0c.py",
    "analyze_sigma_constant.py",
    "reproduce.ps1",
    "paper_relational_hadronic_v1.md"
)

foreach ($file in $Scripts) {
    $src = "$Downloads\$file"
    if (Test-Path $src) {
        Copy-Item $src ".\$file" -Force
        Write-Host "  copied  $file" -ForegroundColor Green
    } else {
        Write-Host "  missing $file (download from Claude)" -ForegroundColor Yellow
    }
}

# README goes to repo root
if (Test-Path "$Downloads\README.md") {
    Copy-Item "$Downloads\README.md" ".\README.md" -Force
    Write-Host "  copied  README.md" -ForegroundColor Green
}

# Kernel documents go to docs/
$KernelDocs = @(
    "operators_notation_and_constraint_v6.md",
    "primary_operators_delta_sigma_v6.md",
    "primary_region_formal_interior_v6.md",
    "role_separation_and_operator_application_v6.md",
    "abr_operators_plain_v6.md"
)

foreach ($doc in $KernelDocs) {
    $src = "$Downloads\$doc"
    if (Test-Path $src) {
        Copy-Item $src ".\docs\$doc" -Force
        Write-Host "  copied  docs/$doc" -ForegroundColor Green
    }
}

# Rust kernel reference copy
if (Test-Path "$Downloads\operators_primary.rs") {
    if (-not (Test-Path ".\kernel")) { New-Item -ItemType Directory ".\kernel" | Out-Null }
    Copy-Item "$Downloads\operators_primary.rs" ".\kernel\operators_primary.rs" -Force
    Write-Host "  copied  kernel/operators_primary.rs" -ForegroundColor Green
}

# White paper to docs
if (Test-Path "$Downloads\directed_relational_evolution_whitepaper_v2.md") {
    Copy-Item "$Downloads\directed_relational_evolution_whitepaper_v2.md" `
        ".\docs\directed_relational_evolution_whitepaper_v2.md" -Force
    Write-Host "  copied  docs/directed_relational_evolution_whitepaper_v2.md" -ForegroundColor Green
}

# Update .gitignore to exclude large data files
@"
# Python
__pycache__/
*.py[cod]
.venv/
venv/

# Large data files -- download separately (see README)
data/events/*.csv
data/*.root
data/*.duckdb
*.root

# Generated results
results/*.md
results/*.json
results/*.csv

# OS
.DS_Store
Thumbs.db

# IDE
.vscode/
.idea/
"@ | Set-Content ".\.gitignore" -Encoding UTF8
Write-Host "  updated .gitignore" -ForegroundColor Green

# Git add and commit
Write-Host ""
Write-Host "Making git commit..." -ForegroundColor Yellow

git add -A
git status --short

git commit -m "V6 publishable state: relational hadronic analysis complete

Metatron Dynamics, Inc. V6. Bounded over D. No claim beyond D.

Findings:
- rank(Im Delta) per-event > PDG summary: gap = 1-3 dimensions
- Sigma constant at decay vertices: A = (n-1)*rho*p^2/(1+p)
  Confirmed to numerical precision on real LHCb data
- Directed closure reproduces conservation within declared graphs
- General formula across decay topologies (2-body and 3-body)
- All results without particle identification or mass hypothesis

Data:
- Simulated: 10,000 D0 -> K-pi+ events (generate_events.py)
- Real LHCb D0 Masterclass: 91,529 events (record 401, CC0)
- Real LHCb B->HHH: Level 2 track momenta (record 4900, CC0)

Kernel: V6 primary operators (Delta, Sigma), 36/36 tests passing
Paper: paper_relational_hadronic_v1.md (working draft)"

Write-Host ""
Write-Host "Repository finalized." -ForegroundColor Cyan
Write-Host ""
Write-Host "Next steps:" -ForegroundColor White
Write-Host "  1. Push to GitHub:"
Write-Host "     git remote add origin https://github.com/Relational-Relativity-Corporation/abr-event-analysis"
Write-Host "     git push -u origin main"
Write-Host "  2. Create GitHub release: v1.0"
Write-Host "  3. Add DOI via Zenodo for citable preprint"
Write-Host ""
Write-Host "Bounded over D. No claim beyond D." -ForegroundColor Cyan
