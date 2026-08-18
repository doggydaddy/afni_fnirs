#!/usr/bin/env python3
"""Compatibility entry point that writes the shared data_aux covariate file.

Prefer ``prep_select_cov_local.py`` for routine use. This entry point preserves
the historical destination ``data_aux/select.cov.1D``.
"""

import os
from pathlib import Path
import sys

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent
SHARED = SCRIPT_DIR / "prep_select_cov_local.py"
PYTHON = Path(os.environ.get(
    "GROUP_ANALYSIS_PYTHON",
    REPO_ROOT / "04_secondlevel_group/.venv/bin/python",
))
if not PYTHON.exists():
    PYTHON = Path(sys.executable)
os.execv(
    str(PYTHON),
    [str(PYTHON), str(SHARED), "--output", str(REPO_ROOT.parent / "data_aux/select.cov.1D"),
     *sys.argv[1:]],
)
