"""README.md as the running project log.

The README is not a static description: it accumulates status, progress,
results, issues, decisions and a change log across runs.  These helpers edit it
section-by-section so a re-run refreshes the current state without deleting
anything that was recorded earlier.
"""

from __future__ import annotations

import datetime as _dt
import re
from pathlib import Path

import config

SECTIONS = [
    "Current Status",
    "1. Project Goal",
    "2. Current Pipeline",
    "3. Data",
    "4. Gene ID Standardization",
    "5. Graph Context",
    "6. Experiment Progress",
    "7. Results",
    "8. Problems / Issues",
    "9. Decisions",
    "10. Next Steps",
    "11. Change Log",
]

TITLE = "# KIRC-Hetionet Project"


def today() -> str:
    return _dt.date.today().isoformat()


# ---------------------------------------------------------------------------
# Low-level section editing
# ---------------------------------------------------------------------------
def _read(path: Path) -> str:
    if path.exists():
        return path.read_text(encoding="utf-8")
    return skeleton()


def skeleton() -> str:
    parts = [TITLE, ""]
    for name in SECTIONS:
        parts += [f"## {name}", "", "_(not yet recorded)_", ""]
    return "\n".join(parts)


def _split_sections(text: str):
    """-> (preamble, [(heading, body), ...]) for level-2 headings."""
    matches = list(re.finditer(r"^## (.+?)\s*$", text, flags=re.MULTILINE))
    if not matches:
        return text, []
    preamble = text[: matches[0].start()]
    out = []
    for i, m in enumerate(matches):
        end = matches[i + 1].start() if i + 1 < len(matches) else len(text)
        out.append((m.group(1).strip(), text[m.end(): end]))
    return preamble, out


def _join_sections(preamble: str, sections) -> str:
    chunks = [preamble.rstrip("\n"), ""]
    for name, body in sections:
        chunks.append(f"## {name}")
        chunks.append(body.strip("\n"))
        chunks.append("")
    return "\n".join(chunks).rstrip("\n") + "\n"


def set_section(name: str, body: str, path: Path | None = None) -> None:
    """Replace a section's body, creating the section if it is missing."""
    path = config.README_PATH if path is None else Path(path)
    text = _read(path)
    preamble, sections = _split_sections(text)
    body = body.strip("\n")

    names = [n for n, _ in sections]
    if name in names:
        sections = [(n, body if n == name else b) for n, b in sections]
    else:
        # insert at the canonical position
        order = {s: i for i, s in enumerate(SECTIONS)}
        sections.append((name, body))
        sections.sort(key=lambda kv: order.get(kv[0], len(order)))

    path.write_text(_join_sections(preamble, sections), encoding="utf-8")


def get_section(name: str, path: Path | None = None) -> str:
    path = config.README_PATH if path is None else Path(path)
    _, sections = _split_sections(_read(path))
    for n, body in sections:
        if n == name:
            return body.strip("\n")
    return ""


def append_dated_block(
    section: str, lines, date: str | None = None, path: Path | None = None
) -> None:
    """Append a ``### <date>`` block to a section.

    A block for the same date is replaced (so re-running today does not spam
    the log); every earlier date is preserved untouched.
    """
    date = today() if date is None else date
    body = get_section(section, path)
    if body.strip() in {"", "_(not yet recorded)_"}:
        body = ""

    block = "\n".join([f"### {date}", ""] + list(lines))

    pattern = re.compile(
        rf"^### {re.escape(date)}\s*$.*?(?=^### |\Z)", flags=re.MULTILINE | re.DOTALL
    )
    if pattern.search(body):
        body = pattern.sub(block + "\n\n", body)
    else:
        body = (body.rstrip("\n") + "\n\n" + block) if body.strip() else block

    set_section(section, body, path)


# ---------------------------------------------------------------------------
# High-level updates
# ---------------------------------------------------------------------------
def update_status(
    stage: str,
    status: str,
    completed,
    next_steps,
    path: Path | None = None,
) -> None:
    lines = [
        f"Stage: {stage}",
        "",
        f"Status: {status}",
        "",
        f"Last updated: {today()}",
        "",
        "Completed:",
    ]
    lines += [f"- {item}" for item in completed]
    lines += ["", "Next:"]
    lines += [f"- {item}" for item in next_steps]
    set_section("Current Status", "\n".join(lines), path)


def update_results(
    comparison,
    extra_lines=(),
    path: Path | None = None,
) -> None:
    """Write the measured comparison table into section 7.

    Only values that come out of an actual run are written - ``None`` is
    rendered as ``n/a`` rather than guessed.
    """
    def cell(value):
        if value is None:
            return "n/a"
        try:
            return f"{int(value):,}"
        except (TypeError, ValueError):
            return str(value)

    contexts = [c for c in comparison.columns if c != "metric"]
    pretty = {"pathway": "Pathway", "biological_process": "Biological Process"}

    header = "| Metric | " + " | ".join(pretty.get(c, c) for c in contexts) + " |"
    rule = "|---|" + "---:|" * len(contexts)
    rows = [
        "| " + str(r["metric"]) + " | "
        + " | ".join(cell(r[c]) for c in contexts) + " |"
        for _, r in comparison.iterrows()
    ]

    lines = [f"Measured on {today()} from an actual pipeline run.", "", header, rule]
    lines += rows
    lines += ["", "Result files:", "", "```text", "results/", "├── pathway/",
              "│   ├── context_nodes.tsv", "│   ├── gene_context_edges.tsv",
              "│   └── subgraph_nodes.tsv", "├── biological_process/",
              "│   ├── context_nodes.tsv", "│   ├── gene_context_edges.tsv",
              "│   └── subgraph_nodes.tsv", "└── comparison/",
              "    ├── context_comparison.tsv",
              "    └── context_comparison_detail.tsv", "```"]
    lines += list(extra_lines)
    set_section("7. Results", "\n".join(lines), path)


def log_progress(items, date: str | None = None, path: Path | None = None) -> None:
    """``items``: list of ``(done: bool, text: str)`` or plain strings."""
    lines = []
    for item in items:
        if isinstance(item, tuple):
            done, text = item
            lines.append(f"- [{'x' if done else ' '}] {text}")
        else:
            lines.append(f"- [x] {item}")
    append_dated_block("6. Experiment Progress", lines, date, path)


def log_change(items, date: str | None = None, path: Path | None = None) -> None:
    append_dated_block("11. Change Log", [f"- {i}" for i in items], date, path)


def log_issues(issues, date: str | None = None, path: Path | None = None) -> None:
    """``issues``: list of dicts with ``what/where/why/status`` keys."""
    lines = []
    for issue in issues:
        lines.append(f"- **{issue.get('what', '(unnamed)')}**")
        lines.append(f"  - Where: {issue.get('where', 'n/a')}")
        lines.append(f"  - Cause: {issue.get('cause', 'n/a')}")
        lines.append(f"  - Resolved: {issue.get('status', 'open')}")
    if not lines:
        lines = ["- none recorded in this run"]
    append_dated_block("8. Problems / Issues", lines, date, path)


def set_next_steps(items, path: Path | None = None) -> None:
    set_section("10. Next Steps", "\n".join(f"- [ ] {i}" for i in items), path)


def set_decisions(items, path: Path | None = None) -> None:
    set_section("9. Decisions", "\n".join(f"- {i}" for i in items), path)
