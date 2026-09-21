"""KIRC-Hetionet pipeline package.

Submodules are imported lazily so that ``import kirc_hetionet`` stays cheap and
a single broken optional dependency cannot take the whole package down.
"""

__all__ = ["hetionet", "kirc", "gene_mapping", "context", "readme_log", "report", "drive"]


def __getattr__(name):
    if name in __all__:
        import importlib

        return importlib.import_module(f"{__name__}.{name}")
    raise AttributeError(name)
