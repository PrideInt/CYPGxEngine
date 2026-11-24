"""What a limited panel costs you, measured against GeT-RM.

A real assay interrogates a panel of alleles and reports *1 for everything it
did not look at. So two identical-looking *1/*1 results can mean "sequenced,
reference at every position" or "we checked two positions". This measures the
difference, using diplotypes whose truth is already known.

    python tools/panel_experiment.py --db data/cypgx.db

Writes results/panel_accuracy.csv and results/panel_confusion.csv.
"""

import argparse
import csv
import io
import os
import sys
from collections import Counter

sys.path.insert(0, ".")

from pgx.core.engine import lookup_gene, open_database
from pgx.panel import (
    frequency_ranked_panel,
    is_structural_diplotype,
    mask,
    named_panel,
    panel_names,
    populations
)

# SEVERITY

# Ordered by metabolic capacity, least to most. Used only to say which
# DIRECTION an error went. Calling a poor metabolizer normal gives a patient a
# standard dose they cannot clear; the reverse is usually a lost opportunity
# rather than a hazard. An aggregate accuracy hides the difference.
RANK = {
    "Poor Metabolizer": 0,
    "Likely Poor Metabolizer": 1,
    "Intermediate Metabolizer": 2,
    "Likely Intermediate Metabolizer": 3,
    "Normal Metabolizer": 4,
    "Rapid Metabolizer": 5,
    "Ultrarapid Metabolizer": 6
}

SWEEP = [0, 1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64]

# Cohort

def read_cohort(path):
    with io.open(path, encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))

def samples_for(rows, gene):
    found = []

    for row in rows:
        value = (row.get(gene) or "").strip()

        if value and "/" in value:
            found.append(value)

    return found

# Evaluation

def phenotype_of(connection, gene, diplotype):
    result = lookup_gene(connection, gene, diplotype)
    return result.get("phenotype"), result.get("error")

def evaluate(connection, gene, samples, panel):
    counts = {
        "n": 0, "correct": 0, "overcalled": 0, "undercalled": 0,
        "sideways": 0, "unresolved": 0, "structural_skipped": 0
    }
    confusion = Counter()

    for diplotype in samples:
        # Copy number and hybrids are a different assay problem. Counted, not
        # approximated.
        if is_structural_diplotype(diplotype):
            counts["structural_skipped"] += 1
            continue

        truth, error = phenotype_of(connection, gene, diplotype)

        if error or truth is None:
            counts["unresolved"] += 1
            continue

        called, call_error = phenotype_of(connection, gene, mask(diplotype, panel))

        if call_error or called is None:
            counts["unresolved"] += 1
            continue

        counts["n"] += 1
        confusion[(truth, called)] += 1

        if truth == called:
            counts["correct"] += 1
            continue

        true_rank = RANK.get(truth)
        called_rank = RANK.get(called)

        if true_rank is None or called_rank is None:
            # Indeterminate on one side or the other: wrong, but not on the
            # scale, so it is not called a direction.
            counts["sideways"] += 1
        elif called_rank > true_rank:
            counts["overcalled"] += 1
        else:
            counts["undercalled"] += 1

    return counts, confusion

def accuracy(counts):
    return counts["correct"] / counts["n"] if counts["n"] else 0.0

# OUTPUT

def ensure_directory(path):
    directory = os.path.dirname(os.path.abspath(path))

    if directory:
        os.makedirs(directory, exist_ok=True)

def write_rows(path, header, rows):
    ensure_directory(path)

    with io.open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", default="data/cypgx.db")
    parser.add_argument("--cohort", default="data/getrm_normalized.csv")
    parser.add_argument("--population", default="European",
                        help="Which population's frequencies rank the sweep")
    parser.add_argument("--accuracy-out", default="results/panel_accuracy.csv")
    parser.add_argument("--confusion-out", default="results/panel_confusion.csv")
    arguments = parser.parse_args()

    connection = open_database(arguments.db)
    rows = read_cohort(arguments.cohort)

    genes = [gene for gene in ("CYP2C19", "CYP2D6", "CYP2B6")
             if gene in (rows[0] if rows else {})]

    accuracy_rows = []
    confusion_rows = []

    for gene in genes:
        samples = samples_for(rows, gene)
        available = populations(connection, gene)

        print(f"\n{gene}  —  {len(samples)} samples", file=sys.stderr)
        print(f"{'panel':<22}{'size':>5}{'n':>6}{'acc':>8}"
              f"{'over':>7}{'under':>7}{'skip':>6}", file=sys.stderr)

        runs = [(name, "", named_panel(gene, name)) for name in panel_names(gene)]

        if arguments.population in available:
            for size in SWEEP:
                panel = frequency_ranked_panel(
                    connection, gene, arguments.population, size)

                # Stop repeating the same panel once frequencies run out.
                if size and len(panel) < size:
                    runs.append((f"top{len(panel)}", arguments.population, panel))
                    break
                runs.append((f"top{size}", arguments.population, panel))

        for name, population, panel in runs:
            counts, confusion = evaluate(connection, gene, samples, panel)

            accuracy_rows.append([
                gene, name, population, len(panel), counts["n"], counts["correct"],
                f"{accuracy(counts):.4f}", counts["overcalled"], counts["undercalled"],
                counts["sideways"], counts["unresolved"], counts["structural_skipped"]
            ])

            for (truth, called), count in sorted(confusion.items()):
                confusion_rows.append([gene, name, population, truth, called, count])

            print(f"{name:<22}{len(panel):>5}{counts['n']:>6}"
                  f"{accuracy(counts):>8.3f}{counts['overcalled']:>7}"
                  f"{counts['undercalled']:>7}{counts['structural_skipped']:>6}",
                  file=sys.stderr)

    write_rows(
        arguments.accuracy_out,
        ["gene", "panel", "population", "panel_size", "n", "correct", "accuracy",
         "overcalled", "undercalled", "sideways", "unresolved", "structural_skipped"],
        accuracy_rows
    )
    write_rows(
        arguments.confusion_out,
        ["gene", "panel", "population", "true_phenotype", "called_phenotype", "count"],
        confusion_rows
    )

    print(f"\nWrote {arguments.accuracy_out} and {arguments.confusion_out}",
          file=sys.stderr)
    print("overcalled = reported as MORE metabolically capable than the truth; "
          "that is the direction that gives a standard dose to someone who "
          "cannot clear it.", file=sys.stderr)

if __name__ == "__main__":
    main()
