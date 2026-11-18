import argparse
import csv
import sys
from collections import defaultdict

import openpyxl

SHEET = "PGx HLA Genotypes"
HEADER_ROW = 2
SAMPLE_COLUMN = "Coriell ID"

GENES = ["CYP2C19", "CYP2D6", "CYP2B6"]

# Locate the sample and gene columns by header text rather than position
def find_columns(header, genes):
    columns = {}

    for index, label in enumerate(header):
        if label is None:
            continue

        text = str(label).strip()

        if text.startswith("Coriell ID"):
            columns["sample"] = index
            continue

        if text in genes:
            columns[text] = index

    missing = [name for name in ["sample"] + genes if name not in columns]

    if missing:
        raise ValueError(f"Columns not found in the sheet: {', '.join(missing)}")

    return columns

# One sample may appear on several study rows, each filling different genes
def merge_rows(rows, columns, genes):
    merged = defaultdict(dict)
    conflicts = []

    for row in rows:
        sample = row[columns["sample"]]

        if sample in (None, ""):
            continue

        sample = str(sample).strip()

        for gene in genes:
            value = row[columns[gene]]

            if value in (None, ""):
                continue

            value = str(value).strip()
            existing = merged[sample].get(gene)

            if existing is not None and existing != value:
                conflicts.append((sample, gene, existing, value))

            merged[sample][gene] = value

    return merged, conflicts

def extract(workbook_path, genes):
    workbook = openpyxl.load_workbook(workbook_path, read_only=True, data_only=True)

    if SHEET not in workbook.sheetnames:
        raise ValueError(
            f"Sheet {SHEET!r} not found. Available: {', '.join(workbook.sheetnames)}"
        )

    sheet = workbook[SHEET]
    rows = list(sheet.iter_rows(min_row=HEADER_ROW, values_only=True))

    columns = find_columns(rows[0], genes)
    merged, conflicts = merge_rows(rows[1:], columns, genes)

    return merged, conflicts

def write_csv(merged, genes, output_path):
    samples = sorted(merged)

    with open(output_path, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["sample"] + genes)

        for sample in samples:
            writer.writerow([sample] + [merged[sample].get(gene, "") for gene in genes])

    return len(samples)

# The CLI stuff
parser = argparse.ArgumentParser(
    description="Extract GeT-RM consensus diplotypes into a tidy CSV"
)
parser.add_argument("workbook")
parser.add_argument("--out", default="data/getrm_diplotypes.csv")
parser.add_argument("--genes", nargs="+", default=GENES)
parser.add_argument(
    "--complete-only",
    action="store_true",
    help="Keep only samples typed for every requested gene"
)

if __name__ == "__main__":
    arguments = parser.parse_args()
    genes = [gene.upper() for gene in arguments.genes]

    merged, conflicts = extract(arguments.workbook, genes)

    if arguments.complete_only:
        merged = {
            sample: values
            for sample, values in merged.items()
            if len(values) == len(genes)
        }

    written = write_csv(merged, genes, arguments.out)

    print(f"Wrote {arguments.out}: {written} samples", file=sys.stderr)

    for gene in genes:
        typed = sum(1 for values in merged.values() if gene in values)
        print(f"  {gene}: {typed} typed", file=sys.stderr)

    complete = sum(1 for values in merged.values() if len(values) == len(genes))
    print(f"  all {len(genes)} genes: {complete} samples", file=sys.stderr)

    if conflicts:
        print(
            f"  {len(conflicts)} conflicting values across duplicate sample rows",
            file=sys.stderr
        )
        for sample, gene, first, second in conflicts[:5]:
            print(f"    {sample} {gene}: {first!r} vs {second!r}", file=sys.stderr)
