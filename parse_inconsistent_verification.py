import json
from pathlib import Path

path = Path('/Users/ribeirrd/Documents/GitHub/MFGraphs/Results/inconsistent_examples_verification.json')
data = json.loads(path.read_text())
for item in data:
    alt = item['AltTransitionFlows']
    actual = alt.count('&&') + 1 if alt else 0
    standard = item['StandardClauseCount']
    extra = max(actual - standard, 0)
    print(json.dumps({
        'Name': item['Name'],
        'ConsistentQ': item['ConsistentQ'],
        'StandardClauseCount': standard,
        'ActualClauseCountDerived': actual,
        'DerivedExtraClauseCount': extra,
        'SampleCrossPairClauses': [
            clause for clause in item.get('ExtraCrossPairClauses', [])
            if ('j[1, 2, 4] == 0, MFGraphs`Private`j[3, 2, 1] == 0' in clause)
            or ('j[1, 2, 3] == 0, MFGraphs`Private`j[4, 2, 1] == 0' in clause)
            or ('j[2, 1, 3] == 0, MFGraphs`Private`j[3, 1, en1] == 0' in clause)
            or ('j[1, 3, 4] == 0, MFGraphs`Private`j[2, 3, 1] == 0' in clause)
        ][:4],
        'CriticalMessage': item['CriticalMessage']
    }))
