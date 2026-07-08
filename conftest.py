# conftest.py
# Metatron Dynamics, Inc. V6. Bounded over D. No claim beyond D.
#
# Adds the repo root to sys.path so pytest can find the declaration,
# operators, and analysis packages from any working directory.

import sys
import os

sys.path.insert(0, os.path.dirname(__file__))
