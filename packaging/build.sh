#!/usr/bin/env bash
#
# Builds both distributions of the CYPGxEngine.
#
#     bash packaging/build.sh
#
# Run it from the project root. Works under Git Bash on Windows, and on Linux
# and macOS.
#
# Produces, under dist/ :
#
#   CYPGxEngine-portable/   a runnable jar plus the engine. Needs Java 17+ on
#                         the target machine; runs on Windows, Linux and
#                         Apple Silicon.
#
#   CYPGxEngine/           a native application folder for THIS platform.
#   CYPGxEngine.app/       Needs nothing installed.
#
# Both carry their own engine and knowledge base, and both start the engine
# themselves. Neither needs Python.
#
# The interpreter used to freeze the engine is chosen automatically: the
# first one found that can import PyInstaller. Override it if you need a
# specific one —
#
#     PYTHON=python3.11 bash packaging/build.sh
#     PYTHON="py -3.10" bash packaging/build.sh
#
# It does not have to be the interpreter you develop against. It only has to
# run the engine, and PyInstaller usually lags new CPython releases by a few
# months.

set -euo pipefail

APP_NAME="CYPGxEngine"
APP_VERSION="0.1.0"
MAIN_CLASS="pgx.client.Launcher"
MAIN_JAR="cypgx-client.jar"

# // OUTPUT

bold=$'\033[1m'; cyan=$'\033[36m'; red=$'\033[31m'; dim=$'\033[2m'; off=$'\033[0m'
step() { printf '\n%s== %s%s\n' "$cyan$bold" "$*" "$off"; }
info() { printf '   %s\n' "$*"; }
note() { printf '   %s%s%s\n' "$dim" "$*" "$off"; }
die()  { printf '\n%sBuild failed:%s %s\n\n' "$red$bold" "$off" "$*" >&2; exit 1; }

# // PLATFORM

case "$(uname -s)" in
    MINGW*|MSYS*|CYGWIN*) PLATFORM=windows; EXE=".exe" ;;
    Darwin)               PLATFORM=mac;     EXE="" ;;
    Linux)                PLATFORM=linux;   EXE="" ;;
    *) die "Unrecognised platform: $(uname -s)" ;;
esac

# Git Bash hands POSIX paths to native Windows tools, which jpackage does not
# always understand. Convert them where it matters.
native_path() {
    if [[ "$PLATFORM" == windows ]]; then cygpath -w "$1"; else printf '%s' "$1"; fi
}

case "$PLATFORM" in
    windows) ICON="packaging/app.ico" ;;
    mac)     ICON="packaging/app.icns" ;;
    linux)   ICON="client/src/main/resources/brand/icon-256.png" ;;
esac

# // PRECONDITIONS

step "Checking the toolchain"
info "platform: $PLATFORM"

need() {
    command -v "$1" >/dev/null 2>&1 || die "$1 is not on PATH. $2"
}
need mvn      "Install Maven."
need java     "Install a JDK 17 or newer."
need jpackage "It ships with the JDK, not the JRE — if java works but jpackage does not, you have a JRE."

# Pick an interpreter that can actually freeze something.
choose_python() {
    local candidates=() c
    [[ -n "${PYTHON:-}" ]] && candidates+=("$PYTHON")
    candidates+=(python3.12 python3.11 python3.10 python3.13 python3 python)

    for c in "${candidates[@]}"; do
        command -v "${c%% *}" >/dev/null 2>&1 || continue
        if $c -c "import PyInstaller" >/dev/null 2>&1; then
            printf '%s' "$c"
            return 0
        fi
    done
    return 1
}

if ! PY="$(choose_python)"; then
    printf '\n%sNo interpreter with PyInstaller.%s Tried:\n\n' "$red$bold" "$off" >&2
    for c in "${PYTHON:-}" python3.12 python3.11 python3.10 python3.13 python3 python; do
        [[ -z "$c" ]] && continue
        if command -v "${c%% *}" >/dev/null 2>&1; then
            printf '   %-14s %s — no PyInstaller\n' \
                "$c" "$($c -c 'import sys; print(sys.version.split()[0])' 2>/dev/null || echo '?')" >&2
        else
            printf '   %-14s not installed\n' "$c" >&2
        fi
    done
    printf '\nInstall it into one of them:\n\n    python -m pip install -r requirements-dev.txt\n' >&2
    printf '\nPyInstaller lags new CPython releases, so if your newest Python refuses,\n' >&2
    printf 'use an older one — it only has to run the engine:\n\n' >&2
    printf '    PYTHON=python3.11 bash packaging/build.sh\n\n' >&2
    exit 1
fi

info "python:   $PY ($($PY -c 'import sys; print(sys.version.split()[0])'))"
info "pyinstaller: $($PY -c 'import PyInstaller; print(PyInstaller.__version__)')"
info "java:     $(java -version 2>&1 | head -1 | tr -d '\r')"

[[ -f data/cypgx.db ]] || die "data/cypgx.db is missing. Build it first:
    python -m pgx.query.cpic --genes CYP2C19 CYP2D6 CYP2B6 --out data/cypgx.db"

[[ -f "$ICON" ]] || { note "no icon at $ICON — building without one"; ICON=""; }

rm -rf build dist
mkdir -p build

# // ENGINE

step "Freezing the engine"
$PY -m PyInstaller packaging/cypgx-engine.spec --noconfirm \
    --distpath build/engine --workpath build/pyinstaller --log-level WARN \
    || die "PyInstaller failed."

ENGINE="build/engine/cypgx-engine"
[[ -x "$ENGINE/cypgx-engine$EXE" ]] || die "cypgx-engine$EXE was not produced in $ENGINE."

# The knowledge base rides beside the executable rather than inside it, so a
# rebuild against a newer CPIC release does not mean refreezing the engine.
cp data/cypgx.db "$ENGINE/cypgx.db"
info "engine: $(du -sh "$ENGINE" | cut -f1)"

# // JARS

step "Building the cross-platform jar"
mvn -q -f client/pom.xml -Pdist,cross clean package || die "Maven failed (cross-platform jar)."
cp "client/target/$MAIN_JAR" build/cypgx-client-all.jar

step "Building the $PLATFORM-only jar"
mvn -q -f client/pom.xml -Pdist clean package || die "Maven failed ($PLATFORM jar)."

# Everything in the jpackage input directory is copied into the application
# image, so the engine goes in here and lands beside the jar. No post-copy,
# and the layout comes out the same on every platform.
mkdir -p build/jpackage-input
cp "client/target/$MAIN_JAR" "build/jpackage-input/$MAIN_JAR"
cp -R "$ENGINE" build/jpackage-input/engine

# // PORTABLE

step "Assembling dist/CYPGxEngine-portable"
PORTABLE="dist/CYPGxEngine-portable"
mkdir -p "$PORTABLE"
cp build/cypgx-client-all.jar "$PORTABLE/$MAIN_JAR"
cp -R "$ENGINE" "$PORTABLE/engine"

cat > "$PORTABLE/README.txt" <<'EOF'
CYPGxEngine — portable
=====================

Needs Java 17 or newer. Check with:  java -version

Then, from this folder:

    java -jar cypgx-client.jar

The application starts the engine in engine/ by itself, on a port it picks
at the time. Nothing listens outside this machine and nothing is sent
anywhere. Close the window and the engine stops with it.

Educational use only. Not clinical decision support. Not validated for
patient care.

Windows, Linux and Apple Silicon. On an Intel Mac the window will fail to
open: rebuild with the `mac` JavaFX classifier in place of `mac-aarch64`,
see client/pom.xml.

The engine folder here was frozen on one platform and only runs on that
one. On a different operating system, rebuild there.

If it does not start, the engine's own output is in your temp directory as
cypgx-engine.log.
EOF

# // NATIVE APPLICATION

step "Building the native application"
JPACKAGE_ARGS=(
    --type app-image
    --name "$APP_NAME"
    --app-version "$APP_VERSION"
    --input "$(native_path build/jpackage-input)"
    --main-jar "$MAIN_JAR"
    --main-class "$MAIN_CLASS"
    --dest "$(native_path dist)"
    --vendor "BMI 530, Stony Brook University"
    --description "CPIC guideline lookup — educational use only"
)
[[ -n "$ICON" ]] && JPACKAGE_ARGS+=(--icon "$(native_path "$ICON")")

jpackage "${JPACKAGE_ARGS[@]}" || die "jpackage failed."

if [[ "$PLATFORM" == mac ]]; then APP_DIR="dist/$APP_NAME.app"; else APP_DIR="dist/$APP_NAME"; fi
[[ -d "$APP_DIR" ]] || die "jpackage reported success but $APP_DIR is not there."

cat > "$APP_DIR/README.txt" <<EOF
CYPGxEngine
==========

Nothing needs installing — Java and Python are both inside this folder.

The application starts its own engine, on the loopback interface only, and
stops it when you close the window.

Educational use only. Not clinical decision support. Not validated for
patient care.

Built for $PLATFORM. It will not run on another operating system; build it
there instead.

If the window opens but says no engine was found, the engine folder is
missing from beside the jar. The engine's own output goes to
cypgx-engine.log in your temp directory.
EOF

# // DONE

step "Done"
for path in "$PORTABLE" "$APP_DIR"; do
    printf '   %-34s %s\n' "$path" "$(du -sh "$path" | cut -f1)"
done
printf '\n   Zip either folder to hand it over.\n\n'
