"""Pytest path setup: make `import proteus...` work without installing the package.

The stages live under src/proteus; tests import them directly (the smoke test
shells out via PYTHONPATH=src, but unit tests import in-process). Prepend src/
to sys.path so `from proteus.utils import ...` resolves from a bare checkout.
"""
import os
import sys

_SRC = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src")
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)
