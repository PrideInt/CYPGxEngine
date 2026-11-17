"""Hardy-Weinberg priors over diplotypes, and posteriors given a panel.

An assay reports *1 for every position it did not interrogate, so a narrow
panel turns many different true diplotypes into the same observation. This
inverts that: given what was reported and which alleles the panel carries,
which true diplotypes could have produced it, and with what probability.

Two caveats that belong in any output built on this, not just in a docstring.

The candidate set is restricted to non-structural alleles that appear in the
gene's diplotype table AND have a published frequency in the population being
used. Everything else gets zero prior and can never be proposed. The retained
frequencies are renormalised, so every number here is conditional on the true
diplotype being inside that set.

*1 carries two different meanings that must not be multiplied together. In an
observation it is the panel's catch-all for anything not interrogated. In the
prior it is CPIC's remainder after every characterised allele, which CPIC
derives rather than observes and flags as inferred.
"""

from .core.engine import normalize_diplotype
from .panel import REFERENCE, is_structural, mask

# GENE MODEL

def gene_model(connection, gene):
    """Everything the posterior needs for one gene, read once.

    Structural alleles are dropped: copy number and hybrids are detected by
    different assay chemistry and do not behave like members of an allele pair.
    """
    phenotypes = {}
    alleles = set()

    rows = connection.execute(
        "SELECT diplotype, phenotype FROM diplotype WHERE gene = ?",
        (gene,)
    )

    for diplotype, phenotype in rows:
        phenotypes[diplotype] = phenotype

        parts = diplotype.split("/")
        if not any(is_structural(part) for part in parts):
            alleles.update(parts)

    return {
        "gene": gene,
        "phenotypes": phenotypes,
        "alleles": sorted(alleles)
    }

# PRIORS

def allele_prior(connection, gene, population, alleles):
    """Renormalised allele frequencies over `alleles`, for one population.

    `coverage` is the share of CPIC's published frequency mass that survived
    the restriction. A coverage of 0.94 means 6% of the population's mass sits
    on alleles this model cannot propose — report it rather than letting
    renormalisation quietly absorb it.
    """
    stored = {
        allele: frequency
        for allele, frequency in connection.execute(
            "SELECT allele, frequency FROM allele_frequency "
            "WHERE gene = ? AND population = ?",
            (gene, population)
        )
    }

    kept = {
        allele: stored[allele]
        for allele in alleles
        if allele in stored and stored[allele] > 0
    }

    published = sum(stored.values())
    retained = sum(kept.values())

    return {
        "population": population,
        "frequencies": normalised(kept),
        "coverage": retained / published if published else 0.0,
        "missing": sorted(allele for allele in alleles if allele not in kept)
    }

def pooled_prior(connection, gene, alleles, populations):
    """One prior for the case where ancestry is unknown, which is the usual one.

    An unweighted mean across populations: every group counts the same,
    regardless of how many people it represents or how many subjects CPIC
    characterised it from. That is a choice, and a defensible one only because
    the alternative — weighting by world population — would make the result an
    estimate about humanity rather than about an unknown patient.
    """
    totals = {}
    coverages = []

    for population in populations:
        prior = allele_prior(connection, gene, population, alleles)
        coverages.append(prior["coverage"])

        for allele, frequency in prior["frequencies"].items():
            totals[allele] = totals.get(allele, 0.0) + frequency

    count = len(populations) or 1

    return {
        "population": "pooled",
        "frequencies": normalised({allele: total / count for allele, total in totals.items()}),
        "coverage": sum(coverages) / count if coverages else 0.0,
        "missing": sorted(allele for allele in alleles if allele not in totals)
    }

def normalised(weights):
    total = sum(weights.values())

    if not total:
        return {}
    
    return {key: value / total for key, value in weights.items()}

# POSTERIOR

def candidates_for(observed_allele, alleles, panel):
    """True alleles that an assay would have reported as `observed_allele`.

    `alleles` is already restricted to those carrying a prior, and every
    candidate is checked against it. An allele a panel genuinely interrogates
    can still have no published frequency — CYP2C19 *7 is on most panels and
    has none for Europeans — and such an observation has no posterior at all
    rather than a wrong one.
    """
    carried = set(panel)
    known = set(alleles)

    if observed_allele != REFERENCE and observed_allele in carried:
        return [observed_allele] if observed_allele in known else []

    if observed_allele == REFERENCE:
        return [allele for allele in alleles if allele not in carried]

    # Reported as something the panel does not carry. Not reachable by masking,
    # so nothing could have produced it.
    return []

def posterior(model, prior, observed, panel):
    """P(phenotype | observed, panel, prior).

    Enumerates the true diplotypes consistent with the observation, weights
    them by Hardy-Weinberg, and collapses by phenotype. Pure, and cheap to
    cache on (observed, panel, population): a narrow panel collapses a whole
    cohort onto a handful of distinct observations.

    Returns an EMPTY dict when no candidate carries a prior — an observation
    the model cannot account for. Callers must count those rather than scoring
    them, because an empty distribution is not a confident one.
    """
    frequencies = prior["frequencies"]
    alleles = [allele for allele in model["alleles"] if allele in frequencies]

    first, second = normalize_diplotype(observed).split("/")
    left = candidates_for(first, alleles, panel)
    right = candidates_for(second, alleles, panel)

    weights = {}
    seen = set()

    for a in left:
        for b in right:
            pair = tuple(sorted((a, b)))

            if pair in seen:
                continue
            
            seen.add(pair)

            diplotype = normalize_diplotype(f"{a}/{b}")
            phenotype = model["phenotypes"].get(diplotype)

            if phenotype is None:
                continue

            # The observation has to be reachable from this pair. Masking is
            # not symmetric enough to assume it.
            if mask(diplotype, panel) != normalize_diplotype(observed):
                continue

            weight = frequencies[a] * frequencies[b] * (1.0 if a == b else 2.0)
            weights[phenotype] = weights.get(phenotype, 0.0) + weight

    return normalised(weights)

def marginal(model, prior):
    """The phenotype distribution with no panel at all.

    Written out directly rather than as posterior(..., panel=[]) so the two can
    be checked against each other. It is also the quantity CPIC publishes as
    phenotype frequency, which makes it the one part of this that can be
    validated against an outside number.
    """
    frequencies = prior["frequencies"]
    alleles = [allele for allele in model["alleles"] if allele in frequencies]

    weights = {}

    for index, a in enumerate(alleles):
        for b in alleles[index:]:
            diplotype = normalize_diplotype(f"{a}/{b}")
            phenotype = model["phenotypes"].get(diplotype)

            if phenotype is None:
                continue

            weight = frequencies[a] * frequencies[b] * (1.0 if a == b else 2.0)
            weights[phenotype] = weights.get(phenotype, 0.0) + weight

    return normalised(weights)

# SCORING

def brier(distribution, truth):
    """Multiclass Brier score against a one-hot truth. Lower is better.

    Bounded, unlike log loss, which is what makes it the headline here: a
    point call that is confidently wrong scores 2.0 rather than infinity.
    """
    labels = set(distribution) | {truth}
    return sum(
        (distribution.get(label, 0.0) - (1.0 if label == truth else 0.0)) ** 2
        for label in labels
    )

def log_loss(distribution, truth, floor=1e-6):
    """Negative log probability of the truth, floored.

    The floor is doing real work: a point call assigns zero to everything it
    did not pick, so every miss would otherwise be infinite. Any comparison
    using this is a comparison of how often that floor is hit.
    """
    import math

    return -math.log(max(distribution.get(truth, 0.0), floor))
