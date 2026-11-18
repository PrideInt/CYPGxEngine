import sys
from collections import Counter, defaultdict

import requests

API_URL = "https://api.cpicpgx.org/v1"
PAGE_SIZE = 1000

GENE = sys.argv[1].upper() if len(sys.argv) > 1 else "CYP2D6"

# Page through an endpoint, optionally with a stable sort order
def fetch_pages(endpoint, order=None, **args):
    rows = []
    offset = 0

    while True:
        query = dict(args, limit=str(PAGE_SIZE), offset=str(offset))

        if order:
            query["order"] = order

        response = requests.get(f"{API_URL}/{endpoint}", params=query, timeout=60)

        response.raise_for_status()
        batch = response.json()

        if not batch:
            return rows

        rows.extend(batch)
        offset += len(batch)

def summarize(rows):
    counts = Counter(row.get("diplotype") for row in rows)

    return {
        "fetched": len(rows),
        "distinct": len(counts),
        "repeated": sum(1 for total in counts.values() if total > 1),
        "counts": counts
    }

# Run the same fetch twice to see whether the row set is reproducible
def run_mode(label, order):
    print(f"\n=== {label} ===")

    first = summarize(fetch_pages("diplotype", order=order, genesymbol=f"eq.{GENE}"))
    second = summarize(fetch_pages("diplotype", order=order, genesymbol=f"eq.{GENE}"))

    for index, run in enumerate((first, second), start=1):
        print(
            f"  run {index}: {run['fetched']} fetched, "
            f"{run['distinct']} distinct, {run['repeated']} repeated"
        )

    stable = first["distinct"] == second["distinct"]
    print(f"  reproducible: {'YES' if stable else 'NO'}")

    if not stable:
        only_first = set(first["counts"]) - set(second["counts"])
        only_second = set(second["counts"]) - set(first["counts"])
        print(f"  diplotypes only in run 1: {len(only_first)}")
        print(f"  diplotypes only in run 2: {len(only_second)}")

    return first

# Report what distinguishes rows that share a diplotype string
def report_duplicates(rows):
    grouped = defaultdict(list)

    for row in rows:
        grouped[row.get("diplotype")].append(row)

    duplicated = {name: group for name, group in grouped.items() if len(group) > 1}

    print(f"\n=== Duplicate analysis ({len(duplicated)} repeated strings) ===")

    if not duplicated:
        print("  No duplicates. The primary key is safe as written.")
        return

    differing = Counter()

    for group in duplicated.values():
        for column in group[0]:
            values = {str(row.get(column)) for row in group}
            if len(values) > 1:
                differing[column] += 1

    if not differing:
        print("  Duplicate rows are byte-identical.")
        print("  Consistent with the same row being fetched twice across page")
        print("  boundaries, not with CPIC publishing a diplotype more than once.")
    else:
        print("  Columns that differ within a duplicated diplotype:")
        for column, total in differing.most_common():
            print(f"    {column:<24} differs in {total} of {len(duplicated)} groups")

print(f"Gene: {GENE}")

unordered = run_mode("Unordered paging, current cpic.py behaviour", None)
ordered = run_mode("Ordered paging, order=diplotype.asc", "diplotype.asc")

report_duplicates(unordered)

print("\n=== Verdict ===")

if ordered["repeated"] == 0 and unordered["repeated"] > 0:
    print("  Ordered paging removes the repeats. Add the order argument to")
    print(f"  cpic.py. Expect {ordered['distinct']} {GENE} diplotypes.")
elif ordered["repeated"] > 0:
    print("  Repeats persist even with ordering. CPIC genuinely returns the")
    print("  same diplotype string more than once; widen the primary key.")
else:
    print("  No repeats in either mode. Nothing to fix.")
