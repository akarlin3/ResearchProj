"""AST-based loader that imports the EXACT pure functions from the existing
analysis modules without running their module-level pipelines.

The paper's analysis.py scripts (absorption-recampaign, transient-tests) execute
top-to-bottom on import (no __main__ guard) and would regenerate/overwrite the
beta=0.05 outputs. To honour "reuse modules, no reimplementation" we parse each
module, keep only its top-level imports and function definitions, and exec that
into a fresh namespace. Module-level constants/data the functions close over
(BR, A_vals, graze, RNG, ...) are injected by the caller before invocation.
"""
import ast
from pathlib import Path


def load_funcs(path):
    """Return a namespace dict with all top-level imports + function defs from
    `path` exec'd into it. No module-level executable statements are run."""
    src = Path(path).read_text()
    tree = ast.parse(src)
    future, imports, funcs = [], [], []
    for node in tree.body:
        if isinstance(node, ast.ImportFrom) and node.module == "__future__":
            future.append(node)
        elif isinstance(node, (ast.Import, ast.ImportFrom)):
            imports.append(node)
        elif isinstance(node, ast.FunctionDef):
            funcs.append(node)
    ns = {}
    # Execute imports individually so an optional/missing dependency (e.g. an ML
    # import used only by functions we do not call) doesn't block loading.
    for node in future + imports:
        m = ast.Module(body=[node], type_ignores=[])
        ast.fix_missing_locations(m)
        try:
            exec(compile(m, filename=str(path), mode="exec"), ns)
        except ModuleNotFoundError:
            pass
    mod = ast.Module(body=funcs, type_ignores=[])
    ast.fix_missing_locations(mod)
    exec(compile(mod, filename=str(path), mode="exec"), ns)
    return ns
