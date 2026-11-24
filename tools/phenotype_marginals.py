"""Predicted phenotype prevalence, by gene and population.

With no panel at all the posterior reduces to the prior marginal: the
Hardy-Weinberg phenotype distribution implied by CPIC's published allele
frequencies. That is a prediction about a population rather than about a
patient, and it is the one quantity in this project that can be checked
against an outside number — CPIC sets includephenotypefrequencies on the gene
record, so the same figures may be published directly.

    python tools/phenotype_marginals.py --db data/cypgx.db

Writes results/phenotype_marginals.csv.

Until that comparison is made these are a derivation, not a validation. They
inherit every assumption in the prior: Hardy-Weinberg equilibrium, independence
between the two alleles, non-structural diplotypes only, and alleles CPIC has
characterised in that population. `coverage` records the share of the
population's published frequency mass that survived those restrictions.
"""

import argparse
import csv
import io
import os
import sys

sys.path.insert(0, ".")

from pgx.core.engine import open_database
from pgx.frequency import allele_prior, gene_model, marginal, pooled_prior
from pgx.panel import populations

GENES = ("CYP2C19", "CYP2D6", "CYP2B6")

def ensure_directory(path):
    directory = os.path.dirname(os.path.abspath(path))

    if directory:
        os.makedirs(directory, exist_ok=True)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", default="data/cypgx.db")
    parser.add_argument("--out", default="results/phenotype_marginals.csv")
    arguments = parser.parse_args()

    connection = open_database(arguments.db)
    rows = []

    for gene in GENES:
        model = gene_model(connection, gene)
        available = populations(connection, gene)

        if not available:
            print(f"{gene}: no frequencies stored", file=sys.stderr)
            continue

        priors = [("pooled", pooled_prior(connection, gene, model["alleles"], available))]
        priors += [
            (name, allele_prior(connection, gene, name, model["alleles"]))
            for name in available
        ]

        print(f"\n{gene}", file=sys.stderr)
        print(f"{'population':<34}{'phenotype':<34}{'p':>7}{'cover':>8}",
              file=sys.stderr)

        for name, prior in priors:
            distribution = marginal(model, prior)

            for phenotype, probability in sorted(
                    distribution.items(), key=lambda item: -item[1]):
                rows.append([
                    gene, name, phenotype,
                    f"{probability:.6f}", f"{prior['coverage']:.4f}",
                    len(prior["frequencies"])
                ])

            # Only the leading phenotypes on screen; the file carries them all.
            for phenotype, probability in sorted(
                    distribution.items(), key=lambda item: -item[1])[:3]:
                print(f"{name:<34}{phenotype:<34}{probability:>7.3f}"
                      f"{prior['coverage']:>8.3f}", file=sys.stderr)

    ensure_directory(arguments.out)

    with io.open(arguments.out, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow([
            "gene", "population", "phenotype", "probability",
            "prior_coverage", "alleles_with_frequency"
        ])
        writer.writerows(rows)

    print(f"\nWrote {arguments.out}: {len(rows)} rows", file=sys.stderr)
    print("A derivation from CPIC's allele frequencies, not yet checked "
          "against CPIC's own phenotype frequencies.", file=sys.stderr)

if __name__ == "__main__":
    main()
