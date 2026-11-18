import argparse
import csv
import os
import sys
from collections import Counter

sys.path.insert(0, ".")

from pgx.activity import (
    band_activity_score,
    parse_activity_score,
    total_activity_score
)
from pgx.query.cpic import fetch_all

GENE = "CYP2D6"

def ensure_directory(path):
    directory = os.path.dirname(os.path.abspath(path))

    if directory:
        os.makedirs(directory, exist_ok=True)

# Three independent checks against the values CPIC publishes for each diplotype
def check(row, gene):
    first = row.get("activityvalue1")
    second = row.get("activityvalue2")
    published_total = row.get("totalactivityscore")
    published_phenotype = row.get("generesult")

    lookup = row.get("lookupkey") or {}
    published_key = lookup.get(gene) if isinstance(lookup, dict) else None

    computed_total = total_activity_score(first, second)
    key_score, key_is_bound = parse_activity_score(published_key)
    computed_phenotype = band_activity_score(key_score, key_is_bound)

    total_parsed, _ = parse_activity_score(published_total)

    if computed_total is None or total_parsed is None:
        sum_outcome = "not_comparable"
    elif computed_total == total_parsed:
        sum_outcome = "agree"
    else:
        sum_outcome = "disagree"

    if published_phenotype is None:
        phenotype_outcome = "not_comparable"
    elif computed_phenotype == published_phenotype:
        phenotype_outcome = "agree"
    else:
        phenotype_outcome = "disagree"

    return {
        "diplotype": row.get("diplotype"),
        "activity_value_1": first,
        "activity_value_2": second,
        "computed_total": str(computed_total) if computed_total is not None else "",
        "published_total": published_total if published_total is not None else "",
        "sum": sum_outcome,
        "lookup_key": published_key or "",
        "bucketed": "yes" if key_is_bound else "no",
        "computed_phenotype": computed_phenotype,
        "published_phenotype": published_phenotype or "",
        "phenotype": phenotype_outcome
    }

def run(gene, report_path):
    rows = fetch_all("diplotype", order="diplotype.asc", genesymbol=f"eq.{gene}")

    results = [check(row, gene) for row in rows]
    seen = {}

    for result in results:
        seen.setdefault(result["diplotype"], result)

    results = list(seen.values())

    sums = Counter(result["sum"] for result in results)
    phenotypes = Counter(result["phenotype"] for result in results)
    bucketed = Counter(
        result["phenotype"] for result in results if result["bucketed"] == "yes"
    )

    ensure_directory(report_path)

    with open(report_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(results[0]))
        writer.writeheader()
        writer.writerows(results)

    return results, sums, phenotypes, bucketed

# The CLI stuff
parser = argparse.ArgumentParser(
    description="Validate activity score arithmetic and banding against CPIC"
)
parser.add_argument("--gene", default=GENE)
parser.add_argument("--report", default="results/activity_score_validation.csv")
parser.add_argument("--show-disagreements", type=int, default=10)

if __name__ == "__main__":
    arguments = parser.parse_args()
    gene = arguments.gene.upper()

    results, sums, phenotypes, bucketed = run(gene, arguments.report)

    print(f"{gene}: {len(results)} distinct diplotypes\n", file=sys.stderr)

    comparable = sums["agree"] + sums["disagree"]
    rate = (100.0 * sums["agree"] / comparable) if comparable else 0.0
    print(
        f"  activity value 1 + 2 == published total: "
        f"{sums['agree']}/{comparable} ({rate:.2f}%)",
        file=sys.stderr
    )
    if sums["not_comparable"]:
        print(
            f"      {sums['not_comparable']} rows lack a per-allele value",
            file=sys.stderr
        )

    comparable = phenotypes["agree"] + phenotypes["disagree"]
    rate = (100.0 * phenotypes["agree"] / comparable) if comparable else 0.0
    print(
        f"  banded lookup key == published phenotype: "
        f"{phenotypes['agree']}/{comparable} ({rate:.2f}%)",
        file=sys.stderr
    )

    if bucketed:
        print("\n  bucketed (>=) lookup keys:", file=sys.stderr)

        for outcome, count in bucketed.most_common():
            print(f"      {outcome:18} {count}", file=sys.stderr)

    disagreements = [
        result for result in results
        if result["sum"] == "disagree" or result["phenotype"] == "disagree"
    ]

    if disagreements:
        print(f"\n  {len(disagreements)} disagreements, first few:", file=sys.stderr)
        
        for result in disagreements[: arguments.show_disagreements]:
            print(
                f"      {result['diplotype']:28} "
                f"values {result['activity_value_1']}+{result['activity_value_2']} "
                f"= {result['computed_total']} vs {result['published_total']}; "
                f"key {result['lookup_key']!r} -> {result['computed_phenotype']} "
                f"vs {result['published_phenotype']}",
                file=sys.stderr
            )

    print(f"\nWrote {arguments.report}", file=sys.stderr)
