# PyInstaller spec for the PGx engine.
#
#   pyinstaller packaging/cypgx-engine.spec --noconfirm
#
# Run from the project root. Produces dist/cypgx-engine/, a folder containing
# cypgx-engine.exe and its _internal/ directory.
#
# One-dir rather than one-file deliberately: a one-file build unpacks itself
# into a temporary directory on every launch, which costs seconds of startup
# and is what most antivirus heuristics object to. The client starts this
# process and waits for it, so startup time is time the user watches.

import os

# PyInstaller resolves paths in a spec relative to the SPEC FILE's directory,
# not the working directory, and injects SPECPATH to say where that is. A
# path written relative to the project root gets "packaging/" prepended and
# is not found.
try:
    spec_dir = os.path.abspath(SPECPATH)
except NameError:                                   # exec'd outside PyInstaller
    spec_dir = os.path.abspath(os.path.dirname(__file__))

project_root = os.path.dirname(spec_dir)
entry_point = os.path.join(spec_dir, "engine_main.py")

# uvicorn resolves its loop, protocol and lifespan implementations by importing
# them by name at runtime. PyInstaller's static analysis cannot see any of it,
# so every one has to be declared or the frozen engine dies on first request.
hiddenimports = [
    "uvicorn.logging",
    "uvicorn.loops",
    "uvicorn.loops.auto",
    "uvicorn.loops.asyncio",
    "uvicorn.protocols",
    "uvicorn.protocols.http",
    "uvicorn.protocols.http.auto",
    "uvicorn.protocols.http.h11_impl",
    "uvicorn.protocols.http.httptools_impl",
    "uvicorn.protocols.websockets",
    "uvicorn.protocols.websockets.auto",
    "uvicorn.protocols.websockets.websockets_impl",
    "uvicorn.protocols.websockets.wsproto_impl",
    "uvicorn.lifespan",
    "uvicorn.lifespan.on",
    "uvicorn.lifespan.off",
    "h11",
    "httptools",
    "websockets",
    "anyio",
    "pgx",
    "pgx.api.main",
    "pgx.core.engine",
    "pgx.models",
]

analysis = Analysis(
    [entry_point],
    # So that `import pgx` resolves to the checkout rather than anything
    # installed. project_root is the directory holding pgx/ and packaging/.
    pathex=[project_root],
    binaries=[],
    # The knowledge base is NOT frozen in. It is 16 MB of generated data with
    # its own version stamp, and it is copied beside the executable by the
    # build script so it can be replaced without rebuilding the engine.
    datas=[],
    hiddenimports=hiddenimports,
    hookspath=[],
    runtime_hooks=[],
    excludes=["tkinter", "matplotlib", "numpy", "pandas", "requests"],
    noarchive=False,
)

pyz = PYZ(analysis.pure)

exe = EXE(
    pyz,
    analysis.scripts,
    [],
    exclude_binaries=True,
    name="cypgx-engine",
    debug=False,
    strip=False,
    upx=False,
    console=True,
)

coll = COLLECT(
    exe,
    analysis.binaries,
    analysis.datas,
    strip=False,
    upx=False,
    name="cypgx-engine",
)
