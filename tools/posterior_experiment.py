"""Does modelling the uncertainty beat pretending there is none?

Part 1 showed that a narrow panel misreports phenotypes. This asks whether a
posterior over phenotypes scores better than the single confident label a lab
would print, on the same samples, against the same known truth.

Two predictions per sample:

    naive      point mass on whatever the masked diplotype resolves to.
               This is what a report says today.
    posterior  P(phenotype | observed, panel, population frequencies).

Scored with Brier, which is bounded — a confidently wrong point call costs
2.0 rather than the infinity log loss would charge it.

    python tools/posterior_experiment.py --db data/cypgx.db

Writes results/posterior_scores.csv.

GeT-RM is enriched for variant alleles by design, so absolute scores here are
not population estimates. The comparison is between two predictions on the
same samples, which the enrichment does not bias.
"""

import argparse
import csv
import io
import os
import sys
from collections import Counter

sys.path.insert(0, ".")

from pgx.core.engine import lookup_gene, open_database
from pgx.frequency import (
    allele_prior,
    brier,
    gene_model,
    log_loss,
    pooled_prior,
    posterior
)
from pgx.panel import (
    frequency_ranked_panel,
    is_structural_diplotype,
    mask,
    named_panel,
    panel_names,
    populations
)

SWEEP = [0, 1, 2, 4, 8, 16, 32]

def read_cohort(path):
    with io.open(path, encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))

def samples_for(rows, gene):
    return [
        (row.get(gene) or "").strip()
        for row in rows
        if (row.get(gene) or "").strip() and "/" in (row.get(gene) or "")
    ]

def phenotype_of(connection, gene, diplotype):
    result = lookup_gene(connection, gene, diplotype)
    return result.get("phenotype") if not result.get("error") else None

def score_panel(connection, gene, model, prior, samples, panel):
    """Both predictions over one cohort, one panel, one prior."""
    cache = {}
    totals = {
        "n": 0,
        "brier_naive": 0.0, "brier_posterior": 0.0,
        "log_naive": 0.0, "log_posterior": 0.0,
        "no_prior": 0, "structural": 0, "unresolved": 0,
        "posterior_better": 0
    }

    for diplotype in samples:
        if is_structural_diplotype(diplotype):
            totals["structural"] += 1
            continue

        truth = phenotype_of(connection, gene, diplotype)

        if truth is None:
            totals["unresolved"] += 1
            continue

        observed = mask(diplotype, panel)
        called = phenotype_of(connection, gene, observed)

        if called is None:
            totals["unresolved"] += 1
            continue

        if observed not in cache:
            cache[observed] = posterior(model, prior, observed, panel)
        distribution = cache[observed]

        # No candidate carried a prior. Not a confident answer and not a bad
        # one — simply outside what this model can speak to.
        if not distribution:
            totals["no_prior"] += 1
            continue

        naive = {called: 1.0}

        naive_brier = brier(naive, truth)
        posterior_brier = brier(distribution, truth)

        totals["n"] += 1
        totals["brier_naive"] += naive_brier
        totals["brier_posterior"] += posterior_brier
        totals["log_naive"] += log_loss(naive, truth)
        totals["log_posterior"] += log_loss(distribution, truth)

        if posterior_brier < naive_brier:
            totals["posterior_better"] += 1

    return totals

def mean(total, n):
    return total / n if n else 0.0

def ensure_directory(path):
    directory = os.path.dirname(os.path.abspath(path))

    if directory:
        os.makedirs(directory, exist_ok=True)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", default="data/cypgx.db")
    parser.add_argument("--cohort", default="data/getrm_normalized.csv")
    parser.add_argument("--out", default="results/posterior_scores.csv")
    parser.add_argument("--rank-population", default="European",
                        help="Whose frequencies order the sweep panels")
    parser.add_argument("--all-populations", action="store_true",
                        help="Score under every population, not just pooled")
    arguments = parser.parse_args()

    connection = open_database(arguments.db)
    rows = read_cohort(arguments.cohort)
    genes = [g for g in ("CYP2C19", "CYP2D6", "CYP2B6") if g in (rows[0] if rows else {})]

    out_rows = []

    for gene in genes:
        model = gene_model(connection, gene)
        samples = samples_for(rows, gene)
        available = populations(connection, gene)

        priors = [("pooled", pooled_prior(connection, gene, model["alleles"], available))]

        if arguments.all_populations:
            priors += [
                (name, allele_prior(connection, gene, name, model["alleles"]))
                for name in available
            ]

        panels = [(name, named_panel(gene, name)) for name in panel_names(gene)]

        # Once the ranking runs out of alleles every larger size returns the
        # same panel, so only the first of each is kept.
        seen = set()

        for size in SWEEP:
            panel = frequency_ranked_panel(
                connection, gene, arguments.rank_population, size)
            key = tuple(sorted(panel))

            if key in seen:
                continue
            seen.add(key)
            panels.append((f"top{len(panel)}", panel))

        print(f"\n{gene}  —  {len(samples)} samples", file=sys.stderr)
        print(f"{'prior':<22}{'panel':<12}{'size':>5}{'n':>6}"
              f"{'naive':>9}{'post':>9}{'gain':>8}{'better':>8}{'nopri':>7}",
              file=sys.stderr)

        for prior_name, prior in priors:
            for panel_name, panel in panels:
                totals = score_panel(
                    connection, gene, model, prior, samples, panel)

                naive = mean(totals["brier_naive"], totals["n"])
                post = mean(totals["brier_posterior"], totals["n"])

                out_rows.append([
                    gene, panel_name, len(panel), prior_name, totals["n"],
                    f"{naive:.4f}", f"{post:.4f}", f"{naive - post:.4f}",
                    f"{mean(totals['log_naive'], totals['n']):.4f}",
                    f"{mean(totals['log_posterior'], totals['n']):.4f}",
                    totals["posterior_better"], totals["no_prior"],
                    totals["structural"], totals["unresolved"],
                    f"{prior['coverage']:.4f}"
                ])

                print(f"{prior_name:<22}{panel_name:<12}{len(panel):>5}"
                      f"{totals['n']:>6}{naive:>9.3f}{post:>9.3f}"
                      f"{naive - post:>8.3f}{totals['posterior_better']:>8}"
                      f"{totals['no_prior']:>7}", file=sys.stderr)

    ensure_directory(arguments.out)

    with io.open(arguments.out, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow([
            "gene", "panel", "panel_size", "prior", "n",
            "brier_naive", "brier_posterior", "brier_gain",
            "logloss_naive", "logloss_posterior",
            "posterior_better", "no_prior", "structural_skipped",
            "unresolved", "prior_coverage"
        ])
        writer.writerows(out_rows)

    print(f"\nWrote {arguments.out}", file=sys.stderr)
    print("gain = naive Brier minus posterior Brier. Positive means modelling "
          "the uncertainty beat asserting none. It should fall to zero as the "
          "panel approaches full, where both become the same point mass.",
          file=sys.stderr)

if __name__ == "__main__":
    main()
