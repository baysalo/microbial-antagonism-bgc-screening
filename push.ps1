param(
    [string]$RepoUrl = "https://github.com/baysalo/microbial-antagonism-bgc-screening.git",
    [string]$Branch = "main",
    [string]$CommitMessage = "Add microbial antagonism, SignalP 6.0, DeepTMHMM, and antiSMASH BGC pipeline"
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

if (-not (Get-Command git -ErrorAction SilentlyContinue)) {
    throw "Git bulunamadi. Lutfen Git for Windows kurun ve PATH'e ekleyin."
}

# Git kimligi tanimli degilse uyar / kontrol et
$gitName = git config user.name 2>$null
$gitEmail = git config user.email 2>$null
if ([string]::IsNullOrWhiteSpace($gitName) -or [string]::IsNullOrWhiteSpace($gitEmail)) {
    throw "Git user.name veya user.email ayarlanmamis. Once su komutlari calistirin:`ngit config --global user.name `"Omur Baysal`"`ngit config --global user.email `"eposta@ornek.com`" "
}

if (-not (Test-Path ".git")) {
    Write-Host "Yerel Git deposu baslatiliyor..." -ForegroundColor Cyan
    git init
}

# Mevcut dali main yap
git branch -M $Branch

# README.md yoksa olustur
if (-not (Test-Path "README.md")) {
    @"
# Microbial Antagonism & BGC Screening Pipeline

Genomic and antagonistic candidate screening workflow integrating heuristic scoring, **SignalP 6.0**, **DeepTMHMM** topology filtering, and **antiSMASH** biosynthetic gene cluster (BGC) mapping (including *Bacillus pumilus* vs. *Staphylococcus aureus* NCTC 8325).
"@ | Out-File -FilePath "README.md" -Encoding utf8
}

git add .

$staged = git status --porcelain
if (-not [string]::IsNullOrWhiteSpace($staged)) {
    git commit -m $CommitMessage
} else {
    Write-Host "Yeni veya degismis commit edilecek dosya yok." -ForegroundColor Yellow
}

$remotes = git remote
if ($remotes -contains "origin") {
    git remote set-url origin $RepoUrl
} else {
    git remote add origin $RepoUrl
}

# Uzak depoda onceden README/license varsa cakismayi onlemek icin senkronize et
git fetch origin $Branch 2>$null
$remoteExists = ($LASTEXITCODE -eq 0)

if ($remoteExists) {
    git pull origin $Branch --allow-unrelated-histories --no-rebase
}

Write-Host "$RepoUrl adresine ($Branch) gonderiliyor..." -ForegroundColor Cyan
git push -u origin $Branch

Write-Host "Push basarili: https://github.com/baysalo/microbial-antagonism-bgc-screening" -ForegroundColor Green