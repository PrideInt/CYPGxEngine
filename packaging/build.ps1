<#
    Builds both distributions of the PGxEngine.

        powershell -ExecutionPolicy Bypass -File packaging\build.ps1

    Run it from the project root, not from packaging\.

    packaging\build.sh does the same thing and also covers Linux and macOS.
    This one exists so Windows needs no bash.

    Produces, under dist\ :

      PGxEngine-portable\   a runnable jar plus the engine. Needs Java 17+ on
                            the target machine; runs on Windows, Linux and
                            Apple Silicon.

      PGxEngine\           a Windows application folder. Needs nothing
                            installed: double-click "PGxEngine.exe".

    The interpreter used to freeze the engine is chosen automatically: the
    first one found that can import PyInstaller. Override it if you need a
    specific one —

        $env:PYTHON = "py -3.11"; powershell -File packaging\build.ps1

    It does not have to be the interpreter you develop against. It only has
    to run the engine, and PyInstaller usually lags new CPython releases by
    a few months.
#>

$ErrorActionPreference = "Stop"

$AppName    = "PGxEngine"
$AppVersion = "0.1.0"
$MainClass  = "pgx.client.Launcher"
$MainJar    = "pgx-client.jar"

function Step($text) { Write-Host "`n== $text" -ForegroundColor Cyan }
function Info($text) { Write-Host "   $text" }
function Need($command, $hint) {
    if (-not (Get-Command $command -ErrorAction SilentlyContinue)) {
        throw "$command is not on PATH. $hint"
    }
}

# // PRECONDITIONS

Step "Checking the toolchain"
Need "mvn"      "Install Maven."
Need "java"     "Install a JDK 17 or newer."
Need "jpackage" "It ships with the JDK, not the JRE - if java works but jpackage does not, you have a JRE."

function Split-Command($text) {
    $parts = $text -split ' +'
    $rest = if ($parts.Length -gt 1) { $parts[1..($parts.Length - 1)] } else { @() }
    return @{ Exe = $parts[0]; Args = $rest }
}

function Choose-Python {
    $candidates = @()
    if ($env:PYTHON) { $candidates += $env:PYTHON }
    $candidates += @("py -3.12", "py -3.11", "py -3.10", "py -3.13", "python3", "python")

    foreach ($candidate in $candidates) {
        $split = Split-Command $candidate
        if (-not (Get-Command $split.Exe -ErrorAction SilentlyContinue)) { continue }
        & $split.Exe @($split.Args) -c "import PyInstaller" 2>$null | Out-Null
        if ($LASTEXITCODE -eq 0) { return $split }
    }
    return $null
}

$python = Choose-Python
if (-not $python) {
    Write-Host "`nNo interpreter with PyInstaller. Tried:" -ForegroundColor Red
    foreach ($candidate in @($env:PYTHON, "py -3.12", "py -3.11", "py -3.10", "py -3.13", "python3", "python")) {
        if (-not $candidate) { continue }
        $split = Split-Command $candidate
        if (Get-Command $split.Exe -ErrorAction SilentlyContinue) {
            $version = (& $split.Exe @($split.Args) -c "import sys; print(sys.version.split()[0])" 2>$null)
            Write-Host ("   {0,-14} {1} - no PyInstaller" -f $candidate, $version)
        } else {
            Write-Host ("   {0,-14} not installed" -f $candidate)
        }
    }
    Write-Host "`nInstall it into one of them:`n"
    Write-Host "    python -m pip install -r requirements-dev.txt`n"
    Write-Host "PyInstaller lags new CPython releases, so if your newest Python refuses,"
    Write-Host "use an older one - it only has to run the engine:`n"
    Write-Host '    $env:PYTHON = "py -3.11"; powershell -File packaging\build.ps1' "`n"
    throw "No usable Python."
}

$pyVersion = & $python.Exe @($python.Args) -c "import sys; print(sys.version.split()[0])"
$piVersion = & $python.Exe @($python.Args) -c "import PyInstaller; print(PyInstaller.__version__)"
Info "python:      $($python.Exe) $($python.Args) ($pyVersion)"
Info "pyinstaller: $piVersion"

if (-not (Test-Path "data\pgx.db")) {
    throw "data\pgx.db is missing. Build it first: python -m pgx.query.cpic --genes CYP2C19 CYP2D6 CYP2B6 --out data\pgx.db"
}

Remove-Item -Recurse -Force build, dist -ErrorAction SilentlyContinue
New-Item -ItemType Directory -Force -Path build | Out-Null

# // ENGINE

Step "Freezing the engine"
& $python.Exe @($python.Args) -m PyInstaller packaging\pgx-engine.spec --noconfirm `
    --distpath build\engine --workpath build\pyinstaller --log-level WARN
if ($LASTEXITCODE -ne 0) { throw "PyInstaller failed." }

$engine = "build\engine\pgx-engine"
if (-not (Test-Path "$engine\pgx-engine.exe")) { throw "pgx-engine.exe was not produced in $engine." }

# The knowledge base rides beside the executable rather than inside it, so a
# rebuild against a newer CPIC release does not mean refreezing the engine.
Copy-Item "data\pgx.db" "$engine\pgx.db"
Info ("engine: {0} MB" -f [math]::Round((Get-ChildItem $engine -Recurse | Measure-Object Length -Sum).Sum / 1MB, 1))

# // JARS

Step "Building the cross-platform jar"
mvn -q -f client\pom.xml -Pdist,cross clean package
if ($LASTEXITCODE -ne 0) { throw "Maven failed (cross-platform jar)." }
Copy-Item "client\target\$MainJar" "build\pgx-client-all.jar"

Step "Building the Windows-only jar"
mvn -q -f client\pom.xml -Pdist clean package
if ($LASTEXITCODE -ne 0) { throw "Maven failed (Windows jar)." }

# Everything in the jpackage input directory is copied into the application
# image, so the engine goes in here and lands beside the jar. No post-copy,
# and the layout comes out the same as the bash script's.
New-Item -ItemType Directory -Force -Path build\jpackage-input | Out-Null
Copy-Item "client\target\$MainJar" "build\jpackage-input\$MainJar"
Copy-Item $engine "build\jpackage-input\engine" -Recurse

# // PORTABLE

Step "Assembling dist\PGxEngine-portable"
$portable = "dist\PGxEngine-portable"
New-Item -ItemType Directory -Force -Path $portable | Out-Null
Copy-Item "build\pgx-client-all.jar" "$portable\$MainJar"
Copy-Item $engine "$portable\engine" -Recurse

@'
PGxEngine - portable
=====================

Needs Java 17 or newer. Check with:  java -version

Then, from this folder:

    java -jar pgx-client.jar

The application starts the engine in engine\ by itself, on a port it picks
at the time. Nothing listens outside this machine and nothing is sent
anywhere. Close the window and the engine stops with it.

Educational use only. Not clinical decision support. Not validated for
patient care.

Windows, Linux and Apple Silicon. On an Intel Mac the window will fail to
open: rebuild with the `mac` JavaFX classifier in place of `mac-aarch64`,
see client\pom.xml.

The engine folder here was frozen on Windows and only runs on Windows. On
another operating system, rebuild there with packaging/build.sh.

If it does not start, the engine's own output is in your temp directory as
pgx-engine.log.
'@ | Set-Content "$portable\README.txt" -Encoding UTF8

# // WINDOWS APPLICATION

Step "Building the Windows application"
$jpackageArgs = @(
    "--type", "app-image",
    "--name", $AppName,
    "--app-version", $AppVersion,
    "--input", "build\jpackage-input",
    "--main-jar", $MainJar,
    "--main-class", $MainClass,
    "--dest", "dist",
    "--vendor", "BMI 530, Stony Brook University",
    "--description", "CPIC guideline lookup - educational use only"
)
if (Test-Path "packaging\app.ico") { $jpackageArgs += @("--icon", "packaging\app.ico") }

& jpackage @jpackageArgs
if ($LASTEXITCODE -ne 0) { throw "jpackage failed." }

$appDir = "dist\$AppName"
if (-not (Test-Path $appDir)) { throw "jpackage reported success but $appDir is not there." }

@'
PGxEngine
==========

Double-click "PGxEngine.exe". Nothing needs installing - Java and Python
are both inside this folder.

The application starts its own engine, on the loopback interface only, and
stops it when you close the window.

Educational use only. Not clinical decision support. Not validated for
patient care.

Built for Windows. It will not run on another operating system; build it
there with packaging/build.sh instead.

If the window opens but says no engine was found, the engine folder is
missing from beside the jar. The engine's own output goes to
pgx-engine.log in your temp directory.
'@ | Set-Content "$appDir\README.txt" -Encoding UTF8

# // DONE

Step "Done"
foreach ($path in @($portable, $appDir)) {
    $mb = [math]::Round((Get-ChildItem $path -Recurse | Measure-Object Length -Sum).Sum / 1MB, 1)
    Write-Host ("   {0,-34} {1} MB" -f $path, $mb)
}
Write-Host "`n   Zip either folder to hand it over.`n"
