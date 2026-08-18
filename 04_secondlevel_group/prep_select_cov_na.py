#!/usr/bin/env python3
"""Compatibility entry point for the shared local covariate builder.

Both second-level models consume the same covariate columns, so the na entry
point now writes the same ignored ``select.cov.local.1D`` file.
"""

import os
from pathlib import Path
import sys

SCRIPT_DIR = Path(__file__).resolve().parent
SHARED = SCRIPT_DIR / "prep_select_cov_local.py"
PYTHON = Path(os.environ.get(
    "GROUP_ANALYSIS_PYTHON",
    SCRIPT_DIR / ".venv/bin/python",
))
if not PYTHON.exists():
    PYTHON = Path(sys.executable)
os.execv(
    str(PYTHON),
    [str(PYTHON), str(SHARED), "--output", str(SCRIPT_DIR / "select.cov.local.1D"),
     *sys.argv[1:]],
)
