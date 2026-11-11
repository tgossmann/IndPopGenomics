#!/usr/bin/env python3
"""Coalescent simulations of a diploid reference donor plus one diploid focal
individual (4 haplotypes), comparing the donor resolved (all 4 haplotypes
observed) with a pseudohaploid (mosaic) reference (donor not resolved).

At every segregating site one of the two donor copies is chosen at random to
form the reference sequence; this is equivalent in expectation to a mosaic
path that switches between the two donor copies.

Outputs: sim_results.csv (one row per replicate and scenario).
Requires: msprime (tested with 1.4.5), numpy.
"""
import msprime, numpy as np, csv, sys

L      = 2_000_000      # sequence length per replicate (bp)
MU     = 1e-8           # mutation rate per bp per generation (infinite sites)
RHO    = 1e-8           # recombination rate per bp per generation
NE     = 10_000
REPS   = 200
SEED   = 20261110

def demography(name):
    d = msprime.Demography()
    if name == "constant":
        d.add_population(name="pop", initial_size=NE)
    elif name == "expansion":          # 20-fold growth 4,000 generations ago
        d.add_population(name="pop", initial_size=20 * NE)
        d.add_population_parameters_change(time=4_000, initial_size=NE)
    elif name == "contraction":        # 5-fold decline 4,000 generations ago
        d.add_population(name="pop", initial_size=NE / 5)
        d.add_population_parameters_change(time=4_000, initial_size=NE)
    return d

def tajimas_d(counts, n):
    """counts: array of derived-allele counts (1..n-1) per segregating site."""
    S = len(counts)
    if S < 2: return np.nan
    pi = np.sum(counts * (n - counts)) / (n * (n - 1) / 2)
    a1 = np.sum(1 / np.arange(1, n)); a2 = np.sum(1 / np.arange(1, n) ** 2)
    b1 = (n + 1) / (3 * (n - 1)); b2 = 2 * (n * n + n + 3) / (9 * n * (n - 1))
    c1 = b1 - 1 / a1; c2 = b2 - (n + 2) / (a1 * n) + a2 / a1 ** 2
    e1 = c1 / a1; e2 = c2 / (a1 ** 2 + a2)
    return (pi - S / a1) / np.sqrt(e1 * S + e2 * S * (S - 1))

def run(scenario, rng, writer):
    ts_iter = msprime.sim_ancestry(
        samples=[msprime.SampleSet(2, ploidy=2)],   # individual 0 = donor, 1 = focal
        demography=demography(scenario), sequence_length=L,
        recombination_rate=RHO, num_replicates=REPS,
        random_seed=rng.integers(1, 2**31))
    for rep, ts in enumerate(ts_iter):
        mts = msprime.sim_mutations(ts, rate=MU, discrete_genome=False,
                                    model=msprime.BinaryMutationModel(),
                                    random_seed=rng.integers(1, 2**31))
        G = mts.genotype_matrix()                   # sites x 4: d1, d2, f1, f2
        donor, focal = G[:, :2], G[:, 2:]
        # donor resolved: all 4 haplotypes
        k4 = G.sum(1)
        m4 = np.minimum(k4, 4 - k4)
        # donor not resolved: reference = one donor copy chosen at random
        pick = rng.integers(0, 2, size=G.shape[0])
        ref = donor[np.arange(G.shape[0]), pick]
        k3 = (focal != ref[:, None]).sum(1)          # focal alleles differing from reference
        seg3 = k3 > 0
        # Tajima's D for the resolved four-haplotype sample (unfolded derived counts).
        # (With three sequences pi equals theta_W at every site, so D is undefined.)
        d4 = tajimas_d(k4, 4)
        writer.writerow(dict(
            scenario=scenario, rep=rep,
            S4=len(k4), M1=int((m4 == 1).sum()), M2=int((m4 == 2).sum()),
            S3=int(seg3.sum()), K1=int((k3 == 1).sum()), K2=int((k3 == 2).sum()),
            lost=int((~seg3).sum()),
            M2_to_K1=int(((m4 == 2) & (k3 == 1)).sum()),
            M1_to_K2=int(((m4 == 1) & (k3 == 2)).sum()),
            M1_lost=int(((m4 == 1) & ~seg3).sum()),
            D4=d4))

if __name__ == "__main__":
    rng = np.random.default_rng(SEED)
    fields = ["scenario","rep","S4","M1","M2","S3","K1","K2","lost",
              "M2_to_K1","M1_to_K2","M1_lost","D4"]
    with open("sim_results.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields); w.writeheader()
        for sc in ["constant", "expansion", "contraction"]:
            run(sc, rng, w); print("done", sc, file=sys.stderr)
