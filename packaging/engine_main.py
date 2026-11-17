"""PyInstaller entry point.

A plain script rather than `python -m pgx.service`, because PyInstaller's
analysis starts from a file and follows its imports. Everything the frozen
engine needs is reachable from here.
"""

import sys

from pgx.service import main

if __name__ == "__main__":
    sys.exit(main())
