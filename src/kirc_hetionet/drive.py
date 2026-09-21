"""Google Drive synchronization.

Everything here reports honestly: a file that could not be written is never
printed as ``[OK]``.  ``save_project_to_drive()`` returns a structured report
and ``verify_drive_backup()`` re-checks the files on disk afterwards rather
than trusting the copy step.

Outside Colab the Drive path usually does not exist.  Set
``KIRC_HETIONET_DRIVE_DIR`` to mirror the same logic into a local directory
(useful for testing the backup code without Colab).
"""

from __future__ import annotations

import os
import shutil
from pathlib import Path

import config


def drive_project_dir() -> Path:
    return Path(os.environ.get("KIRC_HETIONET_DRIVE_DIR", config.DRIVE_PROJECT_DIR))


# ---------------------------------------------------------------------------
# Mount
# ---------------------------------------------------------------------------
def mount_drive(force_remount: bool = False) -> bool:
    """Mount Google Drive when running in Colab.  Returns True on success."""
    if not config.in_colab():
        target = drive_project_dir()
        if "KIRC_HETIONET_DRIVE_DIR" in os.environ:
            print(f"[info] not in Colab; using local mirror {target}")
            return True
        print("[skip] not running in Colab - Google Drive cannot be mounted here")
        return False
    try:
        from google.colab import drive as _drive

        _drive.mount(config.DRIVE_MOUNT_POINT, force_remount=force_remount)
        print(f"[ok  ] Drive mounted at {config.DRIVE_MOUNT_POINT}")
        return True
    except Exception as exc:  # pragma: no cover - Colab only
        print(f"[FAIL] Drive mount failed: {exc}")
        return False


def drive_available() -> bool:
    target = drive_project_dir()
    try:
        target.mkdir(parents=True, exist_ok=True)
        return target.is_dir()
    except OSError as exc:
        print(f"[FAIL] Drive project dir not writable ({target}): {exc}")
        return False


# ---------------------------------------------------------------------------
# Notebook snapshot
# ---------------------------------------------------------------------------
def save_current_notebook(dest: Path | None = None) -> dict:
    """Write the *currently running* notebook to ``dest``.

    Colab does not expose the notebook as a file, so the live JSON is requested
    from the frontend.  When that is unavailable (plain Jupyter, headless run,
    frontend not responding) the repository copy at ``config.NOTEBOOK_PATH`` is
    used instead, and the method actually used is reported back.

    Returns ``{"ok", "method", "path", "detail"}``.
    """
    dest = (drive_project_dir() / "notebooks" / config.NOTEBOOK_NAME) if dest is None else Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)

    if config.in_colab():
        try:  # pragma: no cover - Colab only
            import json

            from google.colab import _message

            payload = _message.blocking_request("get_ipynb", request="", timeout_sec=60)
            notebook = payload["ipynb"] if isinstance(payload, dict) else None
            if notebook:
                dest.write_text(json.dumps(notebook, ensure_ascii=False, indent=1),
                                encoding="utf-8")
                return {
                    "ok": True,
                    "method": "colab:get_ipynb (live notebook)",
                    "path": dest,
                    "detail": "",
                }
            detail = "colab get_ipynb returned no notebook payload"
        except Exception as exc:  # pragma: no cover - Colab only
            detail = f"colab get_ipynb unavailable: {exc}"
    else:
        detail = "not running in Colab"

    src = config.NOTEBOOK_PATH
    if src.exists():
        shutil.copy2(src, dest)
        return {
            "ok": True,
            "method": f"repository copy ({src})",
            "path": dest,
            "detail": detail + " - copied the repository notebook instead; "
                               "unsaved in-session edits are NOT included",
        }

    return {
        "ok": False,
        "method": None,
        "path": dest,
        "detail": f"{detail}; and no notebook at {src}",
    }


# ---------------------------------------------------------------------------
# Project sync
# ---------------------------------------------------------------------------
def _copy(src: Path, dst: Path) -> int:
    """Copy a file or a directory tree.  Returns the number of files written."""
    if src.is_dir():
        n = 0
        for item in sorted(src.rglob("*")):
            if item.is_dir() or "__pycache__" in item.parts or item.name.startswith("."):
                continue
            target = dst / item.relative_to(src)
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(item, target)
            n += 1
        return n
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)
    return 1


def save_project_to_drive(include_notebook: bool = True) -> dict:
    """Synchronize code, notebook, README, data/processed and results to Drive.

    Returns a report dict; nothing is claimed as saved unless the copy
    succeeded.
    """
    target = drive_project_dir()
    report = {"target": target, "copied": {}, "skipped": {}, "failed": {}, "notebook": None}

    if not drive_available():
        print(f"[FAIL] Drive project dir unavailable: {target}")
        print("       Nothing was saved.")
        report["failed"]["<mount>"] = f"{target} unavailable"
        return report

    print(f"--- syncing project to {target} ---")
    for rel in config.DRIVE_SYNC_PATHS:
        src = config.PROJECT_DIR / rel
        if not src.exists():
            report["skipped"][rel] = "not present locally"
            print(f"  [skip] {rel} (not present locally)")
            continue
        try:
            n = _copy(src, target / rel)
            report["copied"][rel] = n
            print(f"  [ok  ] {rel}  ({n} file{'s' if n != 1 else ''})")
        except OSError as exc:
            report["failed"][rel] = str(exc)
            print(f"  [FAIL] {rel}: {exc}")

    if include_notebook:
        nb = save_current_notebook()
        report["notebook"] = nb
        if nb["ok"]:
            print(f"  [ok  ] notebook via {nb['method']}")
            if nb["detail"]:
                print(f"         note: {nb['detail']}")
        else:
            print(f"  [FAIL] notebook not saved: {nb['detail']}")

    return report


# ---------------------------------------------------------------------------
# Verification
# ---------------------------------------------------------------------------
def verify_drive_backup(verbose: bool = True) -> dict:
    """Check that the expected files really exist in Drive after a backup."""
    target = drive_project_dir()
    status = {}
    for rel, label in config.DRIVE_VERIFY_PATHS:
        path = target / rel
        try:
            ok = path.exists() and path.stat().st_size > 0
        except OSError:
            ok = False
        status[label] = {"ok": ok, "path": path}
        if verbose:
            print(f"[{'OK  ' if ok else 'FAIL'}] {label}"
                  + ("" if ok else f"  (missing: {path})"))
    status["all_ok"] = all(v["ok"] for k, v in status.items() if isinstance(v, dict))
    return status
