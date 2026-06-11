import os
from pathlib import Path


def use_agg_backend():
    os.environ.setdefault("MPLCONFIGDIR", str(Path(".matplotlib-cache").resolve()))
    os.environ.setdefault("XDG_CACHE_HOME", str(Path(".cache").resolve()))

    import matplotlib

    matplotlib.use("Agg")
