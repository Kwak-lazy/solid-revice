"""Hetionet v1.0 download, loading and type inspection.

Keeps the behaviour of the original notebook cells (download -> load nodes /
edges -> inspect node kinds -> inspect edge types) behind small reusable
functions.  Loaded frames are cached in module state so the notebook can call
``load_nodes()`` repeatedly without re-parsing a 2.2M-row file.
"""

from __future__ import annotations

import hashlib
import urllib.request
from pathlib import Path

import pandas as pd

import config

_CACHE: dict[str, pd.DataFrame] = {}


# ---------------------------------------------------------------------------
# Download
# ---------------------------------------------------------------------------
def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _download(url: str, dest: Path, force: bool = False) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0 and not force:
        print(f"[skip] {dest.name} already present ({dest.stat().st_size:,} bytes)")
        return dest
    print(f"[get ] {url}")
    urllib.request.urlretrieve(url, dest)
    print(f"[ok  ] {dest} ({dest.stat().st_size:,} bytes)")
    return dest


def download_hetionet(force: bool = False) -> dict[str, Path]:
    """Download the Hetionet v1.0 node and edge TSVs into data/raw/.

    The edge file lives in git-lfs; ``config.HETIONET_EDGES_URL`` therefore
    points at the media endpoint.  A 133-byte result means an LFS pointer was
    served instead of the payload, which we detect via the checksum.
    """
    nodes = _download(config.HETIONET_NODES_URL, config.HETIONET_NODES_FILE, force)
    edges = _download(config.HETIONET_EDGES_URL, config.HETIONET_EDGES_FILE, force)

    digest = _sha256(edges)
    if digest != config.HETIONET_EDGES_SHA256:
        raise RuntimeError(
            f"Hetionet edge file checksum mismatch for {edges}.\n"
            f"  expected {config.HETIONET_EDGES_SHA256}\n"
            f"  got      {digest}\n"
            "A tiny file usually means a git-lfs pointer was downloaded "
            "instead of the payload."
        )
    print(f"[ok  ] edge checksum verified ({digest[:12]}...)")
    return {"nodes": nodes, "edges": edges}


def download_hgnc(force: bool = False) -> Path:
    """Download the HGNC complete set (Ensembl <-> HGNC <-> Entrez bridge)."""
    return _download(config.HGNC_URL, config.HGNC_FILE, force)


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------
def load_nodes(force_reload: bool = False) -> pd.DataFrame:
    """Hetionet nodes: columns ``id``, ``name``, ``kind``."""
    if "nodes" not in _CACHE or force_reload:
        if not config.HETIONET_NODES_FILE.exists():
            download_hetionet()
        _CACHE["nodes"] = pd.read_csv(
            config.HETIONET_NODES_FILE, sep="\t", dtype=str, keep_default_na=False
        )
    return _CACHE["nodes"]


def load_edges(force_reload: bool = False) -> pd.DataFrame:
    """Hetionet edges: columns ``source``, ``metaedge``, ``target``."""
    if "edges" not in _CACHE or force_reload:
        if not config.HETIONET_EDGES_FILE.exists():
            download_hetionet()
        _CACHE["edges"] = pd.read_csv(
            config.HETIONET_EDGES_FILE, sep="\t", dtype=str, keep_default_na=False
        )
    return _CACHE["edges"]


# ---------------------------------------------------------------------------
# Type inspection
# ---------------------------------------------------------------------------
def node_kind_counts(nodes: pd.DataFrame | None = None) -> pd.DataFrame:
    nodes = load_nodes() if nodes is None else nodes
    out = nodes["kind"].value_counts().rename_axis("kind").reset_index(name="n_nodes")
    return out


def metaedge_counts(edges: pd.DataFrame | None = None) -> pd.DataFrame:
    edges = load_edges() if edges is None else edges
    out = (
        edges["metaedge"]
        .value_counts()
        .rename_axis("metaedge")
        .reset_index(name="n_edges")
    )
    return out


def describe_hetionet(nodes=None, edges=None) -> None:
    """Print the node-type / edge-type overview the original notebook showed."""
    nodes = load_nodes() if nodes is None else nodes
    edges = load_edges() if edges is None else edges
    print(f"Hetionet nodes : {len(nodes):,}")
    print(f"Hetionet edges : {len(edges):,}")
    print("\n--- node kinds ---")
    print(node_kind_counts(nodes).to_string(index=False))
    print("\n--- metaedges ---")
    print(metaedge_counts(edges).to_string(index=False))


# ---------------------------------------------------------------------------
# Gene nodes
# ---------------------------------------------------------------------------
def gene_nodes(nodes: pd.DataFrame | None = None) -> pd.DataFrame:
    """Hetionet Gene nodes with the Entrez id parsed out of ``Gene::<entrez>``.

    Returns ``hetionet_gene_id``, ``hetionet_gene_symbol``, ``entrez_id``.
    """
    nodes = load_nodes() if nodes is None else nodes
    genes = nodes[nodes["kind"] == "Gene"].copy()
    genes = genes.rename(
        columns={"id": "hetionet_gene_id", "name": "hetionet_gene_symbol"}
    )
    genes["entrez_id"] = (
        genes["hetionet_gene_id"]
        .str.removeprefix(config.HETIONET_GENE_PREFIX)
        .astype(str)
    )
    return genes[["hetionet_gene_id", "hetionet_gene_symbol", "entrez_id"]].reset_index(
        drop=True
    )
