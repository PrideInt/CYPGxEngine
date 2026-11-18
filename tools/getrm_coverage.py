import argparse
import csv
import sys
from collections import Counter, defaultdict

sys.path.insert(0, ".")

from pgx.core.engine import normalize_diplotype, open_database
from pgx.notation import NOT_A_DIPLOTYPE, expand_unknown_copies, normalize_getrm

UNRESOLVED = "no_cpic_match"

# Try every candidate reading against the CPIC diplotype table, best first
def resolve(connection, gene, raw):
    interpretation = normalize_getrm(raw)
    attempted = []

    for candidate in interpretation["candidates"]:
        for expanded in expand_unknown_copies(candidate):
            forms = [expanded]

            try:
                canonical = normalize_diplotype(expanded)
                if canonical != expanded:
                    forms.append(canonical)
            except ValueError:
                pass

            for form in forms:
                attempted.append(form)

                row = connection.execute(
                    "SELECT diplotype, phenotype FROM diplotype "
                    "WHERE gene = ? AND diplotype = ?",
                    (gene, form)
                ).fetchone()

                if row is not None:
                    return {
                        "resolved": True,
                        "matched": row["diplotype"],
                        "phenotype": row["phenotype"],
                        "note": interpretation["note"],
                        "attempted": attempted,
                        "rank": interpretation["candidates"].index(candidate)
                    }

    return {
        "resolved": False,
        "matched": None,
        "phenotype": None,
        "note": interpretation["note"] if not interpretation["candidates"]
        else UNRESOLVED,
        "attempted": attempted,
        "rank": None
    }

def load_samples(csv_path, genes):
    with open(csv_path, newline="") as handle:
        reader = csv.DictReader(handle)

        if reader.fieldnames is None:
            raise ValueError(f"No header row found in {csv_path!r}")

        missing = [gene for gene in genes if gene not in reader.fieldnames]

        if missing:
            raise ValueError(f"Columns not found: {', '.join(missing)}")

        return list(reader)

def run(database_path, csv_path, genes, normalized_path, report_path):
    connection = open_database(database_path)
    samples = load_samples(csv_path, genes)

    audit = []
    normalized = []
    outcomes = defaultdict(Counter)
    unmatched = defaultdict(Counter)

    for row in samples:
        sample = row.get("sample", "")
        resolved_row = {"sample": sample}

        for gene in genes:
            raw = (row.get(gene) or "").strip()

            if not raw:
                outcomes[gene]["not_typed"] += 1
                resolved_row[gene] = ""
                continue

            result = resolve(connection, gene, raw)

            audit.append({
                "sample": sample,
                "gene": gene,
                "raw": raw,
                "matched": result["matched"] or "",
                "phenotype": result["phenotype"] or "",
                "note": result["note"],
                "resolved": "yes" if result["resolved"] else "no",
                "used_alternate": "yes" if result["rank"] else "no",
                "attempted": " | ".join(result["attempted"])
            })

            if result["resolved"]:
                outcomes[gene]["resolved"] += 1
                resolved_row[gene] = result["matched"]
            else:
                outcomes[gene][result["note"]] += 1
                unmatched[gene][raw] += 1
                resolved_row[gene] = ""

        normalized.append(resolved_row)

    with open(normalized_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["sample"] + genes)
        writer.writeheader()
        writer.writerows(normalized)

    with open(report_path, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "sample", "gene", "raw", "matched", "phenotype",
                "note", "resolved", "used_alternate", "attempted"
            ]
        )
        writer.writeheader()
        writer.writerows(audit)

    return samples, outcomes, unmatched

# The CLI stuff
parser = argparse.ArgumentParser(
    description="Resolve GeT-RM consensus calls against the CPIC diplotype table"
)
parser.add_argument("samples", nargs="?", default="data/getrm_diplotypes.csv")
parser.add_argument("--db", default="data/cypgx.db")
parser.add_argument("--genes", nargs="+", default=["CYP2C19", "CYP2D6", "CYP2B6"])
parser.add_argument("--normalized", default="data/getrm_normalized.csv")
parser.add_argument("--report", default="results/getrm_resolution.csv")
parser.add_argument("--show-unmatched", type=int, default=12)

if __name__ == "__main__":
    arguments = parser.parse_args()
    genes = [gene.upper() for gene in arguments.genes]

    samples, outcomes, unmatched = run(
        arguments.db,
        arguments.samples,
        genes,
        arguments.normalized,
        arguments.report
    )

    print(f"{len(samples)} samples from {arguments.samples}", file=sys.stderr)

    for gene in genes:
        counts = outcomes[gene]
        typed = sum(count for note, count in counts.items() if note != "not_typed")
        resolved = counts["resolved"]
        rate = (100.0 * resolved / typed) if typed else 0.0

        print(f"\n  {gene}: {resolved}/{typed} resolved ({rate:.1f}%)", file=sys.stderr)

        for note, count in counts.most_common():
            if note in ("resolved", "not_typed"):
                continue
            print(f"      {note:24} {count}", file=sys.stderr)

        if unmatched[gene]:
            print("      unmatched values:", file=sys.stderr)
            for raw, count in unmatched[gene].most_common(arguments.show_unmatched):
                print(f"        {raw!r} x{count}", file=sys.stderr)

    print(f"\nWrote {arguments.normalized} and {arguments.report}", file=sys.stderr)
