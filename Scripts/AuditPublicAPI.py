#!/usr/bin/env python3
from __future__ import annotations

import re
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = ROOT / 'MFGraphs'
TEST_DIR = ROOT / 'MFGraphs' / 'Tests'
OUT_DIR = ROOT / 'docs' / 'planning'
API_REFERENCE = ROOT / 'API_REFERENCE.md'

USAGE_RE = re.compile(r'(?m)^([A-Za-z$][A-Za-z0-9$]*)::usage\s*=')
HEADER_RE = re.compile(r'^##\s+(.+?)\s*$', re.M)
PRIVATE_BEGIN_RE = re.compile(r'Begin\["`Private`"\]')

STOP_NAMES = {
    'Module', 'Block', 'With', 'Function', 'If', 'Which', 'Switch', 'Do', 'Table',
    'Map', 'Select', 'Cases', 'DeleteCases', 'Join', 'Lookup', 'AssociationThread',
    'Return', 'Failure', 'True', 'False', 'Null', 'Missing', 'Length', 'AnyTrue',
    'Options', 'Get', 'EndPackage',
    'AllTrue', 'VectorQ', 'NumericQ', 'Insert', 'While', 'Quiet', 'Check', 'Print',
    'Echo', 'Set', 'SetDelayed', 'ReplaceAll', 'Rule', 'RuleDelayed', 'First', 'Last',
    'Rest', 'Part', 'Flatten', 'Sort', 'Union', 'Normal', 'N', 'And', 'Or', 'Not'
}

CATEGORY_OVERRIDES = {
    'IsFeasible': ('core', 'keep'),
    'makeScenario': ('core', 'keep'),
    'validateScenario': ('core', 'keep'),
    'completeScenario': ('core', 'keep'),
    'ScenarioData': ('advanced', 'keep'),
    'scenario': ('advanced', 'review'),
    'scenarioQ': ('advanced', 'keep'),
    'Data2Equations': ('compatibility', 'keep-with-deprecation'),
    'FinalStep': ('compatibility', 'keep-with-deprecation'),
    'RemoveDuplicates': ('compatibility', 'keep-with-deprecation'),
    'ReplaceSolution': ('compatibility', 'keep-with-deprecation'),
    'alpha$': ('accidental', 'fix-export'),
    'j$': ('accidental', 'fix-export'),
    'DataG': ('compatibility', 'review'),
}


def read_text(path: Path) -> str:
    return path.read_text(encoding='utf-8', errors='ignore')


def build_usage_map() -> dict[str, list[str]]:
    out: dict[str, list[str]] = defaultdict(list)
    for path in sorted(SRC_DIR.rglob('*.wl')):
        rel = str(path.relative_to(ROOT))
        for sym in USAGE_RE.findall(read_text(path)):
            out[sym].append(rel)
    return dict(out)


def exported_symbols() -> list[str]:
    if API_REFERENCE.exists():
        return HEADER_RE.findall(read_text(API_REFERENCE))
    return sorted(build_usage_map())


def classify_symbol(symbol: str) -> tuple[str, str, str]:
    if symbol in CATEGORY_OVERRIDES:
        category, decision = CATEGORY_OVERRIDES[symbol]
    elif symbol.startswith('$'):
        category, decision = 'runtime-global', 'review'
    elif symbol in {'alpha', 'Cost', 'j', 'u', 'z'}:
        category, decision = 'symbolic-head', 'keep'
    elif re.fullmatch(r'[IUS]\d+', symbol):
        category, decision = 'example-parameter', 'review'
    elif symbol.startswith(('Build', 'Compute', 'Encode', 'Decode', 'Extract', 'Lookup', 'Select', 'Use')):
        category, decision = 'advanced', 'review'
    else:
        category, decision = 'advanced', 'review'

    if category == 'accidental':
        stability = 'broken-export'
    elif category == 'compatibility':
        stability = 'compatibility'
    elif decision == 'keep':
        stability = 'stable-candidate'
    else:
        stability = 'review-needed'
    return category, decision, stability


def tests_for_symbol(symbol: str, test_files: list[Path]) -> list[str]:
    pattern = re.compile(r'(?<![A-Za-z0-9$`])' + re.escape(symbol) + r'(?![A-Za-z0-9$`])')
    hits = []
    for path in test_files:
        if pattern.search(read_text(path)):
            hits.append(str(path.relative_to(ROOT)))
    return hits


def private_candidates() -> dict[str, list[str]]:
    out: dict[str, list[str]] = {}
    usage_names = set(build_usage_map())
    for path in sorted(SRC_DIR.rglob('*.wl')):
        txt = read_text(path)
        m = PRIVATE_BEGIN_RE.search(txt)
        if not m:
            continue
        names = []
        for line in txt[m.end():].splitlines():
            if not line.strip() or line.lstrip().startswith('(*'):
                continue
            if line[:1].isspace():
                continue
            mfun = re.match(r'^([A-Za-z$][A-Za-z0-9$]*)\[', line)
            mvar = re.match(r'^([A-Za-z$][A-Za-z0-9$]*)\s*=', line)
            name = mfun.group(1) if mfun else (mvar.group(1) if mvar else None)
            if not name:
                continue
            if name in STOP_NAMES or name == 'End':
                continue
            if name in usage_names:
                continue
            if name.startswith('$'):
                continue
            if name not in names:
                names.append(name)
            if len(names) >= 12:
                break
        if names:
            out[str(path.relative_to(ROOT))] = names
    return out


def write_inventory(rows: list[dict], priv: dict[str, list[str]]) -> Path:
    out = OUT_DIR / 'public_api_symbol_inventory.md'
    lines = [
        '# MFGraphs Phase 1 public API symbol inventory',
        '',
        'This document inventories the currently exported symbols in the **MFGraphs public context** and classifies each symbol as part of the Phase 1 publicization review.',
        '',
        f'The current exported surface contains **{len(rows)}** symbols.',
        '',
        '## Inventory table',
        '',
        '| Symbol | Declared in | Category | Decision | Stability | Test hits | Notes |',
        '|---|---|---|---|---|---:|---|',
    ]
    for row in rows:
        notes = []
        if row['duplicate_declaration']:
            notes.append('duplicate usage declaration')
        if row['weird_export']:
            notes.append('appears accidental')
        if row['test_hits'] == 0:
            notes.append('no direct test reference found')
        if row['symbol'].startswith('$'):
            notes.append('runtime global')
        if re.fullmatch(r'[IUS]\d+', row['symbol']):
            notes.append('example parameter symbol')
        note_text = '; '.join(notes) if notes else '—'
        declared_in = ', '.join(f'`{p}`' for p in row['declared_in']) if row['declared_in'] else '—'
        lines.append(f"| `{row['symbol']}` | {declared_in} | {row['category']} | {row['decision']} | {row['stability']} | {row['test_hits']} | {note_text} |")
    lines += [
        '',
        '## Notable private implementation families',
        '',
        'The list below is intentionally approximate. It highlights representative helper definitions visible after entering the private section of each source module, which makes them useful candidates for later wrapper extraction or publicization review.',
        '',
        '| Module | Representative private helper names |',
        '|---|---|',
    ]
    for module, defs in priv.items():
        lines.append(f"| `{module}` | {', '.join(f'`{d}`' for d in defs)} |")
    lines += [
        '',
        '## Phase 1 observations',
        '',
        'The exported surface is broader than the currently documented user narrative. It mixes core workflow functions, advanced solver helpers, runtime globals, symbolic heads, compatibility aliases, and several example-parameter symbols. It also contains suspicious exports such as `alpha$` and `j$`, which appear to be packaging artefacts rather than intentional API symbols.',
    ]
    out.write_text('\n'.join(lines) + '\n', encoding='utf-8')
    return out


def write_gap_report(rows: list[dict]) -> Path:
    out = OUT_DIR / 'public_api_gap_report.md'
    no_tests = [r for r in rows if r['test_hits'] == 0]
    accidental = [r for r in rows if r['weird_export']]
    compat = [r for r in rows if r['category'] == 'compatibility']
    globals_ = [r for r in rows if r['category'] == 'runtime-global']
    params = [r for r in rows if r['category'] == 'example-parameter']
    low_trust = [r for r in rows if r['test_hits'] == 0 and r['category'] not in {'example-parameter', 'accidental'}]

    lines = [
        '# MFGraphs Phase 1 public API gap report',
        '',
        'This report summarizes the main gaps between the current exported symbol surface and the desired end state of a fully public, fully documented, and trusted MFGraphs API.',
        '',
        '## Summary',
        '',
        '| Metric | Count |',
        '|---|---:|',
        f'| Exported symbols in current surface | {len(rows)} |',
        f'| Symbols with no direct test reference | {len(no_tests)} |',
        f'| Compatibility or deprecation symbols | {len(compat)} |',
        f'| Runtime globals exported as public symbols | {len(globals_)} |',
        f'| Example-parameter symbols exported as public symbols | {len(params)} |',
        f'| Suspicious or accidental exports | {len(accidental)} |',
        '',
        '## Highest-priority gaps',
        '',
        '| Gap | Why it matters | Recommended Phase 2 action |',
        '|---|---|---|',
        '| Accidental exports such as `alpha$` and `j$` | These weaken confidence in the API boundary and pollute generated docs | Fix package export hygiene before broadening the public surface |',
        '| Broad advanced-helper surface with sparse direct tests | Publicization without verification would reduce trust instead of increasing it | Add a symbol-to-test coverage matrix and targeted low-level unit tests |',
        '| Runtime globals exposed without explicit stability policy | Users cannot tell whether these are supported controls or incidental internals | Decide which globals remain public and document the rest as internal controls |',
        '| Example-parameter symbols exported beside functional API symbols | These clutter the generated reference and blur the distinction between examples and supported functions | Decide whether parameters stay exported, move to examples-only docs, or get grouped specially |',
        '| Compatibility aliases mixed into the main public surface | This makes the API look larger and less coherent than it really is | Mark them consistently in docs and consider dedicated compatibility sections |',
        '',
        '## Symbols with no direct test reference',
        '',
        '| Symbol | Category | Decision |',
        '|---|---|---|',
    ]
    for row in low_trust:
        lines.append(f"| `{row['symbol']}` | {row['category']} | {row['decision']} |")
    lines += [
        '',
        '## Suspicious exports to fix first',
        '',
        '| Symbol | Reason |',
        '|---|---|',
    ]
    for row in accidental:
        lines.append(f"| `{row['symbol']}` | Exported in the generated API surface but not suitable as a stable public symbol |")
    lines += [
        '',
        '## Recommendation',
        '',
        'Before making additional helpers intentionally public, Phase 2 should first normalize the exported surface: remove accidental exports, label compatibility symbols clearly, and define a stable classification rule for globals, symbolic heads, and example parameters. Only then should the repository expand usage contracts and per-symbol verification.',
    ]
    out.write_text('\n'.join(lines) + '\n', encoding='utf-8')
    return out


def main() -> None:
    usage_map = build_usage_map()
    exported = exported_symbols()
    test_files = sorted(TEST_DIR.rglob('*.mt'))
    rows = []
    for symbol in exported:
        declared_in = usage_map.get(symbol, [])
        category, decision, stability = classify_symbol(symbol)
        hits = tests_for_symbol(symbol, test_files)
        rows.append({
            'symbol': symbol,
            'declared_in': declared_in,
            'category': category,
            'decision': decision,
            'stability': stability,
            'test_hits': len(hits),
            'duplicate_declaration': len(declared_in) > 1,
            'weird_export': category == 'accidental',
        })
    rows.sort(key=lambda r: r['symbol'].lower())
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    inv = write_inventory(rows, private_candidates())
    gap = write_gap_report(rows)
    print(inv)
    print(gap)


if __name__ == '__main__':
    main()
