import json
import sqlite3
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from pgx.query.cpic import SCHEMA
from pgx.models import KB_SCHEMA

# CPIC publishes a priority band, an allele function pair and consultation
# prose for every diplotype. The fixture carries a sample for the entries the
# tests assert on; the banding strings are CPIC's own vocabulary.
EXTRAS = {
    ("CYP2C19", "*1/*1"): ("Normal/Routine/Low Risk", ["Normal function", "Normal function"]),
    ("CYP2C19", "*2/*2"): ("Abnormal/Priority/High Risk", ["No function", "No function"]),
    ("CYP2D6", "*4/*4"): ("Abnormal/Priority/High Risk", ["No function", "No function"]),
    ("CYP2D6", "*1/*10"): ("Normal/Routine/Low Risk", ["Normal function", "Decreased function"]),
    ("CYP2B6", "*6/*6"): ("Abnormal/Priority/High Risk", ["No function", "No function"])
}

# Phenotype labels
ULTRARAPID = "Ultrarapid Metabolizer"
RAPID = "Rapid Metabolizer"
NORMAL = "Normal Metabolizer"
INTERMEDIATE = "Intermediate Metabolizer"
POOR = "Poor Metabolizer"
LIKELY_INTERMEDIATE = "Likely Intermediate Metabolizer"
LIKELY_POOR = "Likely Poor Metabolizer"
INDETERMINATE = "Indeterminate"
NO_RESULT = "No Result"

# Diplotype to phenotype mappings
DIPLOTYPES = [
    ("CYP2C19", "*17/*17", ULTRARAPID, None),
    ("CYP2C19", "*1/*17", RAPID, None),
    ("CYP2C19", "*1/*1", NORMAL, None),
    ("CYP2C19", "*1/*9", LIKELY_INTERMEDIATE, None),
    ("CYP2C19", "*9/*9", LIKELY_INTERMEDIATE, None),
    ("CYP2C19", "*1/*2", INTERMEDIATE, None),
    ("CYP2C19", "*1/*3", INTERMEDIATE, None),
    ("CYP2C19", "*2/*17", INTERMEDIATE, None),
    ("CYP2C19", "*2/*9", LIKELY_POOR, None),
    ("CYP2C19", "*2/*2", POOR, None),
    ("CYP2C19", "*2/*3", POOR, None),
    ("CYP2C19", "*1/*12", INDETERMINATE, None),

    ("CYP2D6", "*1/*1xN", ULTRARAPID, "\u22653.0"),
    ("CYP2D6", "*1/*1x2", ULTRARAPID, "3.0"),
    ("CYP2D6", "*2x2/*10", NORMAL, "2.25"),
    ("CYP2D6", "*1/*1", NORMAL, "2.0"),
    ("CYP2D6", "*1/*2", NORMAL, "2.0"),
    ("CYP2D6", "*1/*17", NORMAL, "1.5"),
    ("CYP2D6", "*1/*10", NORMAL, "1.25"),
    ("CYP2D6", "*1/*4", INTERMEDIATE, "1.0"),
    ("CYP2D6", "*1/*5", INTERMEDIATE, "1.0"),
    ("CYP2D6", "*10/*41", INTERMEDIATE, "0.75"),
    ("CYP2D6", "*10/*10", INTERMEDIATE, "0.5"),
    ("CYP2D6", "*4/*10", INTERMEDIATE, "0.25"),
    ("CYP2D6", "*3/*4", POOR, "0.0"),
    ("CYP2D6", "*4/*4", POOR, "0.0"),
    ("CYP2D6", "*5/*5", POOR, "0.0"),
    ("CYP2D6", "*1/*22", INDETERMINATE, "n/a"),

    ("CYP2B6", "*4/*4", ULTRARAPID, None),
    ("CYP2B6", "*1/*4", RAPID, None),
    ("CYP2B6", "*1/*1", NORMAL, None),
    ("CYP2B6", "*1/*6", INTERMEDIATE, None),
    ("CYP2B6", "*4/*18", INTERMEDIATE, None),
    ("CYP2B6", "*6/*6", POOR, None),
    ("CYP2B6", "*6/*18", POOR, None),
    ("CYP2B6", "*3/*6", INDETERMINATE, None)
]

# CYP2D6 recommendations join on activity score, not phenotype
ACTIVITY_SCORE_PHENOTYPES = {
    "\u22653.0": ULTRARAPID,
    "3.0": ULTRARAPID,
    "2.25": NORMAL,
    "2.0": NORMAL,
    "1.5": NORMAL,
    "1.25": NORMAL,
    "1.0": INTERMEDIATE,
    "0.75": INTERMEDIATE,
    "0.5": INTERMEDIATE,
    "0.25": INTERMEDIATE,
    "0.0": POOR,
    "n/a": INDETERMINATE
}

STANDARD_DOSE = "Initiate therapy with recommended starting dose."
NO_RECOMMENDATION = "No recommendation."

# Recommendation tables
RECOMMENDATIONS = [
    ("paroxetine", {"CYP2D6": ULTRARAPID}, "Moderate",
     "Select an alternative not predominantly metabolized by CYP2D6."),
    ("paroxetine", {"CYP2D6": NORMAL}, "Strong", STANDARD_DOSE),
    ("paroxetine", {"CYP2D6": INTERMEDIATE}, "Optional",
     "Consider a lower starting dose and slower titration."),
    ("paroxetine", {"CYP2D6": POOR}, "Moderate",
     "Consider a 50 percent reduction in starting and maintenance dose."),
    ("paroxetine", {"CYP2D6": INDETERMINATE}, "No Recommendation", NO_RECOMMENDATION),

    ("fluvoxamine", {"CYP2D6": ULTRARAPID}, "No Recommendation",
     "No recommendation due to lack of evidence."),
    ("fluvoxamine", {"CYP2D6": NORMAL}, "Strong", STANDARD_DOSE),
    ("fluvoxamine", {"CYP2D6": INTERMEDIATE}, "Moderate", STANDARD_DOSE),
    ("fluvoxamine", {"CYP2D6": POOR}, "Optional",
     "Consider a 25 to 50 percent lower starting dose, or an alternative."),
    ("fluvoxamine", {"CYP2D6": INDETERMINATE}, "No Recommendation", NO_RECOMMENDATION),

    ("venlafaxine", {"CYP2D6": ULTRARAPID}, "No Recommendation",
     "No action recommended. Minimal evidence of clinical impact."),
    ("venlafaxine", {"CYP2D6": NORMAL}, "Strong", STANDARD_DOSE),
    ("venlafaxine", {"CYP2D6": INTERMEDIATE}, "No Recommendation",
     "No action recommended. Minimal evidence of clinical impact."),
    ("venlafaxine", {"CYP2D6": POOR}, "Optional",
     "Consider an alternative not predominantly metabolized by CYP2D6."),
    ("venlafaxine", {"CYP2D6": INDETERMINATE}, "No Recommendation", NO_RECOMMENDATION),

    ("vortioxetine", {"CYP2D6": ULTRARAPID}, "Optional",
     "Select an alternative, or titrate to a higher maintenance dose."),
    ("vortioxetine", {"CYP2D6": NORMAL}, "Strong", STANDARD_DOSE),
    ("vortioxetine", {"CYP2D6": INTERMEDIATE}, "Moderate", STANDARD_DOSE),
    ("vortioxetine", {"CYP2D6": POOR}, "Moderate",
     "Initiate at 50 percent of starting dose, or consider an alternative."),
    ("vortioxetine", {"CYP2D6": INDETERMINATE}, "No Recommendation", NO_RECOMMENDATION),

    ("escitalopram", {"CYP2C19": ULTRARAPID}, "Strong",
     "Consider an alternative not predominantly metabolized by CYP2C19."),
    ("escitalopram", {"CYP2C19": RAPID}, "Optional",
     "Standard starting dose. Consider titration or an alternative if no response."),
    ("escitalopram", {"CYP2C19": NORMAL}, "Strong", STANDARD_DOSE),
    ("escitalopram", {"CYP2C19": INTERMEDIATE}, "Moderate",
     "Standard starting dose. Consider slower titration and lower maintenance."),
    ("escitalopram", {"CYP2C19": LIKELY_INTERMEDIATE}, "Moderate",
     "Standard starting dose. Consider slower titration and lower maintenance."),
    ("escitalopram", {"CYP2C19": LIKELY_POOR}, "Strong",
     "Consider an alternative, or a 50 percent maintenance dose reduction."),
    ("escitalopram", {"CYP2C19": POOR}, "Strong",
     "Consider an alternative, or a 50 percent maintenance dose reduction."),
    ("escitalopram", {"CYP2C19": INDETERMINATE}, "No Recommendation", NO_RECOMMENDATION),

    ("sertraline", {"CYP2C19": NORMAL, "CYP2B6": NO_RESULT}, "Strong",
     STANDARD_DOSE),
    ("sertraline", {"CYP2C19": POOR, "CYP2B6": NO_RESULT}, "Moderate",
     "Consider a 50 percent maintenance dose reduction, or an alternative."),

    ("sertraline", {"CYP2C19": NO_RESULT, "CYP2B6": NORMAL}, "Strong",
     STANDARD_DOSE),
    ("sertraline", {"CYP2C19": NO_RESULT, "CYP2B6": POOR}, "Optional",
     "Consider a 25 percent maintenance dose reduction, or an alternative."),

    ("sertraline", {"CYP2C19": NORMAL, "CYP2B6": NORMAL}, "Strong", STANDARD_DOSE),
    ("sertraline", {"CYP2C19": NORMAL, "CYP2B6": INTERMEDIATE}, "Moderate",
     "Standard starting dose. Consider slower titration and lower maintenance."),
    ("sertraline", {"CYP2C19": NORMAL, "CYP2B6": POOR}, "Optional",
     "Consider a 25 percent maintenance dose reduction, or an alternative."),
    ("sertraline", {"CYP2C19": INTERMEDIATE, "CYP2B6": POOR}, "Optional",
     "Consider a 50 percent maintenance dose reduction."),
    # keyed on a gene this build does not load, so it is not answerable here
    ("warfarin", {"CYP2C9": NORMAL}, "Moderate",
     "Use a published pharmacogenetic dosing algorithm."),

    ("sertraline", {"CYP2C19": POOR, "CYP2B6": POOR}, "Optional",
     "Select an alternative not primarily metabolized by CYP2C19 or CYP2B6.")
]

PAIRS = [
    ("CYP2D6", "paroxetine", "A"),
    ("CYP2B6", "paroxetine", "D"),
    ("CYP2D6", "fluvoxamine", "A"),
    ("CYP2D6", "venlafaxine", "A"),
    ("CYP2D6", "vortioxetine", "A"),
    ("CYP2C19", "escitalopram", "A"),
    ("CYP2C19", "sertraline", "A"),
    ("CYP2B6", "sertraline", "A"),
    ("CYP2C9", "warfarin", "A")
]

# Expand phenotype-keyed CYP2D6 recommendations into score-keyed rows
def expand_lookup_keys(lookup_key):
    if "CYP2D6" not in lookup_key:
        return [lookup_key]

    phenotype = lookup_key["CYP2D6"]
    scores = [
        {"CYP2D6": score}
        for score, mapped in ACTIVITY_SCORE_PHENOTYPES.items()
        if mapped == phenotype
    ]
    return scores or [lookup_key]

# Build the fixture database
def build(database_path="tests/fixture.db"):
    Path(database_path).unlink(missing_ok=True)
    connection = sqlite3.connect(database_path)
    connection.executescript(SCHEMA)

    for gene, diplotype, phenotype, activity_score in DIPLOTYPES:
        if gene == "CYP2D6" and activity_score:
            lookup_key = {gene: activity_score}
        else:
            lookup_key = {gene: phenotype}
        priority, functions = EXTRAS.get((gene, diplotype), (None, []))

        connection.execute(
            "INSERT INTO diplotype ("
            "  gene, diplotype, phenotype, activity_score, total_activity_score,"
            "  function_1, function_2, ehr_priority, description,"
            "  consultation_text, lookupkey"
            ") VALUES (?,?,?,?,?,?,?,?,?,?,?)",
            (
                gene,
                diplotype,
                phenotype,
                activity_score,
                activity_score,
                functions[0] if functions else None,
                functions[1] if len(functions) > 1 else None,
                priority,
                "Fixture description.",
                "Fixture consultation text." if priority else None,
                json.dumps(lookup_key, sort_keys=True)
            )
        )

    # Store recommendations
    recommendation_id = 0
    for drug, lookup_key, classification, text in RECOMMENDATIONS:
        for expanded in expand_lookup_keys(lookup_key):
            recommendation_id += 1
            connection.execute(
                "INSERT INTO recommendation VALUES (?,?,?,?,?,?,?,?,?,?)",
                (
                    str(recommendation_id),
                    drug,
                    drug,
                    json.dumps(expanded, sort_keys=True),
                    classification,
                    text,
                    "Fixture implication text.",
                    None,
                    "adults",
                    "https://example.invalid/guideline"
                )
            )

    # Store pairs and provenance
    for gene, drug, cpic_level in PAIRS:
        connection.execute(
            "INSERT INTO pair VALUES (?,?,?,?)",
            (gene, drug, drug, cpic_level)
        )

    connection.executemany(
        "INSERT INTO meta VALUES (?,?)",
        [
            ("cpic_db", "FIXTURE"),
            ("pharmvar", "FIXTURE"),
            ("retrieved_at", "1970-01-01T00:00:00+00:00"),
            ("kb_schema", KB_SCHEMA)
        ]
    )
    connection.commit()
    connection.close()
    return database_path

if __name__ == "__main__":
    print(build(sys.argv[1] if len(sys.argv) > 1 else "tests/fixture.db"))
