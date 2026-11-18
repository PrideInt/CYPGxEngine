"""Check the derived phenotype marginal against CPIC's own published figure.

With no panel the posterior reduces to the Hardy-Weinberg marginal implied by
CPIC's allele frequencies. CPIC separately publishes phenotype prevalence per
biogeographic group on its gene_result endpoint. The two are computed from the
same underlying data by different routes, so agreement is evidence the
derivation is sound and disagreement localises where it is not.

    python tools/validate_marginals.py --db data/cypgx.db

Writes results/marginal_validation.csv.

This is the only quantity in the project that can be checked against an
outside number. Everything else is either the guideline text itself or a
measurement against GeT-RM.
"""

import argparse
import csv
import io
import os
import sys

sys.path.insert(0, ".")

from pgx.core.engine import open_database
from pgx.frequency import allele_prior, gene_model, marginal
from pgx.panel import populations

GENES = ("CYP2C19", "CYP2D6", "CYP2B6")

# Agreement to within a percentage point is as much as can be asked: CPIC's
# figure and ours are rounded differently and may come from different cuts.
CLOSE = 0.01


def published(connection, gene, population):
    return {
        phenotype: frequency
        for phenotype, frequency in connection.execute(
            "SELECT phenotype, frequency FROM phenotype_frequency "
            "WHERE gene = ? AND population = ?",
            (gene, population)
        )
    }


def ensure_directory(path):
    directory = os.path.dirname(os.path.abspath(path))

    if directory:
        os.makedirs(directory, exist_ok=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", default="data/cypgx.db")
    parser.add_argument("--out", default="results/marginal_validation.csv")
    arguments = parser.parse_args()

    connection = open_database(arguments.db)

    try:
        connection.execute("SELECT 1 FROM phenotype_frequency LIMIT 1")
    except Exception:
        print("No phenotype_frequency table. Rebuild the knowledge base:",
              file=sys.stderr)
        print("  python -m pgx.query.cpic --genes CYP2C19 CYP2D6 CYP2B6 "
              "--out data/cypgx.db", file=sys.stderr)
        return 2

    rows = []

    for gene in GENES:
        model = gene_model(connection, gene)
        available = populations(connection, gene)

        print(f"\n{gene}", file=sys.stderr)
        print(f"  {'population':<34}{'n':>4}{'max |diff|':>12}"
              f"{'mean |diff|':>13}  worst", file=sys.stderr)

        for population in available:
            theirs = published(connection, gene, population)

            if not theirs:
                print(f"  {population:<34}   —  no published figure", file=sys.stderr)
                continue

            prior = allele_prior(connection, gene, population, model["alleles"])
            ours = marginal(model, prior)

            # Compare on the union: a phenotype one side omits is a real
            # disagreement, not something to drop.
            labels = sorted(set(ours) | set(theirs))
            differences = []

            for phenotype in labels:
                mine = ours.get(phenotype, 0.0)
                cpic = theirs.get(phenotype, 0.0)
                difference = mine - cpic
                differences.append((abs(difference), phenotype, mine, cpic))

                rows.append([
                    gene, population, phenotype,
                    f"{mine:.6f}", f"{cpic:.6f}", f"{difference:+.6f}",
                    "yes" if abs(difference) <= CLOSE else "no",
                    f"{prior['coverage']:.4f}", len(prior["frequencies"])
                ])

            differences.sort(reverse=True)
            worst, worst_phenotype, mine, cpic = differences[0]
            mean = sum(d for d, _, _, _ in differences) / len(differences)

            print(f"  {population:<34}{len(labels):>4}{worst:>12.4f}"
                  f"{mean:>13.4f}  {worst_phenotype} "
                  f"({mine:.3f} vs {cpic:.3f})", file=sys.stderr)

    ensure_directory(arguments.out)

    with io.open(arguments.out, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow([
            "gene", "population", "phenotype", "derived", "cpic_published",
            "difference", "within_1pct", "prior_coverage",
            "alleles_with_frequency"
        ])
        writer.writerows(rows)

    agreeing = sum(1 for row in rows if row[6] == "yes")
    print(f"\nWrote {arguments.out}: {len(rows)} comparisons, "
          f"{agreeing} within {CLOSE:.0%}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
