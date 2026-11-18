import argparse
import io
import json
import pathlib
import sys
from collections import Counter, defaultdict

DRUGS = [
    "escitalopram",
    "sertraline",
    "paroxetine",
    "fluvoxamine",
    "venlafaxine",
    "vortioxetine"
]

ACTIONABLE = {"Strong", "Moderate", "Optional"}

REDIRECTION_HINT = (
    "This looks like UTF-16, which Windows shell redirection writes. "
    "Rerun with python -m pgx.cli --batch <csv> --out <file> so Python "
    "writes the file itself."
)

# Windows shell redirection may write UTF-16, so decode by the byte order mark
def decode_reports(path):
    data = pathlib.Path(path).read_bytes()

    for mark, encoding in (
        (b"\xff\xfe\x00\x00", "utf-32"),
        (b"\x00\x00\xfe\xff", "utf-32"),
        (b"\xff\xfe", "utf-16"),
        (b"\xfe\xff", "utf-16")
    ):
        if data.startswith(mark):
            raise ValueError(f"{path}: {REDIRECTION_HINT}")

    return data.decode("utf-8-sig")

def load_reports(path):
    reports = []

    with io.StringIO(decode_reports(path)) as handle:
        for number, line in enumerate(handle, start=1):
            line = line.strip()

            if not line:
                continue

            try:
                reports.append(json.loads(line))
            except json.JSONDecodeError as error:
                raise ValueError(
                    f"{path} line {number} is not JSON ({error})."
                ) from None

    if not reports:
        raise ValueError(f"{path} contains no reports.")

    return reports

# How many samples carried each phenotype, per gene
def phenotype_census(reports):
    census = defaultdict(Counter)

    for report in reports:
        for result in report["genes"]:
            if result["phenotype"]:
                census[result["gene"]][result["phenotype"]] += 1

    return census

# For each drug: how many samples got a recommendation, and of what strength
def drug_census(reports, drugs):
    reached = Counter()
    strengths = defaultdict(Counter)
    genes_used = defaultdict(Counter)

    for report in reports:
        seen = set()

        for recommendation in report["drugs"]:
            drug = recommendation["drug"].lower()

            if drugs and drug not in drugs:
                continue

            seen.add(drug)
            strengths[drug][recommendation["classification"] or "none"] += 1
            genes_used[drug]["+".join(recommendation["genes"]) or "none"] += 1

        for drug in seen:
            reached[drug] += 1

    return reached, strengths, genes_used

# The advice given to a normal metabolizer, per drug AND per set of genes that
# informed the match. A sample typed on one gene of a two-gene drug legitimately
# receives different wording from one typed on both, even when every phenotype
# is normal, so a single per-drug baseline would score that difference as a
# change in dosing.
def baseline_key(drug, recommendation):
    return (drug, tuple(recommendation["genes"]))

def standard_text(reports, drugs):
    baseline = {}

    for report in reports:
        phenotypes = {
            result["gene"]: result["phenotype"] for result in report["genes"]
        }

        for recommendation in report["drugs"]:
            drug = recommendation["drug"].lower()

            if drugs and drug not in drugs:
                continue

            informing = recommendation["genes"]

            if not informing:
                continue

            if any(phenotypes.get(g) != "Normal Metabolizer" for g in informing):
                continue

            baseline.setdefault(
                baseline_key(drug, recommendation),
                recommendation["recommendation"]
            )

    return baseline

# Recommendations that differ from the normal-metabolizer baseline for that drug
def differs_from_standard(reports, drugs, baseline):
    changed = Counter()
    offered = Counter()

    for report in reports:
        for recommendation in report["drugs"]:
            drug = recommendation["drug"].lower()

            if drugs and drug not in drugs:
                continue

            if recommendation["classification"] == "No Recommendation":
                continue

            offered[drug] += 1
            key = baseline_key(drug, recommendation)

            if key not in baseline:
                continue

            if recommendation["recommendation"] != baseline[key]:
                changed[drug] += 1

    return changed, offered

# Samples where at least one drug carried an actionable recommendation
def actionable_samples(reports, drugs):
    total = 0

    for report in reports:
        for recommendation in report["drugs"]:
            if drugs and recommendation["drug"].lower() not in drugs:
                continue

            if recommendation["classification"] in ACTIONABLE:
                total += 1
                break

    return total

# Samples where at least one drug differs from its normal-metabolizer baseline.
# This is the meaningful figure: receiving a recommendation is near-universal,
# receiving a changed one is not.
def samples_with_a_change(reports, drugs, baseline):
    total = 0
    per_sample = Counter()

    for report in reports:
        changed = 0

        for recommendation in report["drugs"]:
            drug = recommendation["drug"].lower()

            if drugs and drug not in drugs:
                continue

            if recommendation["classification"] == "No Recommendation":
                continue

            key = baseline_key(drug, recommendation)

            if key in baseline and recommendation["recommendation"] != baseline[key]:
                changed += 1

        per_sample[changed] += 1

        if changed:
            total += 1

    return total, per_sample

def main():
    parser = argparse.ArgumentParser(
        description="Summarize a batch of PGx reports"
    )
    parser.add_argument("reports", nargs="?", default="results/getrm_reports.jsonl")
    parser.add_argument("--drugs", nargs="+", default=DRUGS)
    parser.add_argument("--all-drugs", action="store_true")
    arguments = parser.parse_args()

    drugs = None if arguments.all_drugs else {d.lower() for d in arguments.drugs}
    reports = load_reports(arguments.reports)

    print(f"{len(reports)} samples from {arguments.reports}\n")
    print("=== Phenotype distribution ===")

    for gene, counts in sorted(phenotype_census(reports).items()):
        total = sum(counts.values())
        print(f"\n  {gene} ({total} typed)")

        for phenotype, count in counts.most_common():
            share = 100.0 * count / total
            print(f"      {phenotype:34} {count:4}  {share:5.1f}%")

    reached, strengths, genes_used = drug_census(reports, drugs)

    print("\n=== Drug coverage ===")
    for drug in sorted(reached, key=lambda d: -reached[d]):
        share = 100.0 * reached[drug] / len(reports) if reports else 0.0
        print(f"\n  {drug} ({reached[drug]}/{len(reports)} samples, {share:.1f}%)")

        for classification, count in strengths[drug].most_common():
            print(f"      {classification:24} {count:4}")
        if len(genes_used[drug]) > 1:
            print("      genes informing the match:")

            for combination, count in genes_used[drug].most_common():
                print(f"        {combination:30} {count:4}")

    baseline = standard_text(reports, drugs)
    changed, offered = differs_from_standard(reports, drugs, baseline)

    print("\n=== Recommendations differing from the normal-metabolizer baseline ===")
    covered = {drug for drug, genes in baseline}

    for drug in sorted(offered, key=lambda d: -changed[d]):
        if drug not in covered:
            print(f"  {drug:16} no normal-metabolizer baseline in this cohort")
            continue
        share = 100.0 * changed[drug] / offered[drug] if offered[drug] else 0.0
        print(f"  {drug:16} {changed[drug]:4}/{offered[drug]:<4} {share:5.1f}%")

    total = actionable_samples(reports, drugs)
    share = 100.0 * total / len(reports) if reports else 0.0
    print(
        f"\n{total}/{len(reports)} samples ({share:.1f}%) carry at least one "
        "recommendation of any strength"
    )

    changed_total, per_sample = samples_with_a_change(reports, drugs, baseline)
    changed_share = 100.0 * changed_total / len(reports) if reports else 0.0
    print(
        f"{changed_total}/{len(reports)} samples ({changed_share:.1f}%) differ "
        "from standard dosing for at least one drug"
    )

    print("\n  drugs changed per sample:")
    for count in sorted(per_sample):
        share = 100.0 * per_sample[count] / len(reports)
        print(f"      {count} drugs  {per_sample[count]:4}  {share:5.1f}%")

    warned = sum(1 for report in reports if report["warnings"])
    print(f"\n{warned}/{len(reports)} samples produced at least one warning")

if __name__ == "__main__":
    main()
