# IndPopGenomics

Code for: Kylmänen A, Javaheri Tehrani S, Vellnow N, Gossmann TI.
*How to Do Population Genetics with One Sample and a Reference Genome.*

The repository contains the coalescent simulations behind Figure 2 and the empirical
analysis of great tit (*Parus major*) genomes behind Figures 3–4, Table 3 and
Supplementary Table S1.

## Contents

| Folder / script | What it does | Manuscript |
|---|---|---|
| `simulations/simulate_mosaic.py` | msprime simulations of one diploid reference donor plus one diploid focal individual; constant size, expansion, contraction; 200 × 2 Mb replicates; pseudohaploid reference made by picking one donor copy at random per site. Writes `sim_results.csv`. | Figure 2, Methods |
| `simulations/summarise_sim.py` | Mean and SD per scenario → `sim_summary.csv` | Figure 2 |
| `simulations/plot_mosaic.py` | Plots Figure 2 from `sim_summary.csv` | Figure 2 |
| `empirical/config.sh` | Paths to the reference assembly and alignments (edit first) | |
| `empirical/01_call_variants_freebayes.sh` | FreeBayes calling of the 11 configurations, BCFtools filtering, strict triallelic records | Methods, Figure 3, Table 3 |
| `empirical/02_popgen_stats.sh` | Biallelic subset; VCFtools windowed π and Tajima's D (50 kb) | Figure 4 |
| `empirical/03_plot_figure3_triallelic.py` | Bar plot of triallelic records | Figure 3 |
| `empirical/04_plot_figure4_configurations.py` | π, θ_W and Tajima's D for five configurations (2, 3, 4, 8, 10 haplotypes) | Figure 4 |
| `empirical/05_pi_with_without_multiallelic.py` | π with and without multiallelic records, chromosome-wide and per window | Supplementary Table S1 |
| `empirical/classify_triallelic_snp_contributors.py` | Genotype classes of triallelic records | Table 3 |
| `empirical/results/` | Small output tables used in the manuscript | |

## Data

- Reference assembly: Parus_major1.1, GCF_001522545.3; analyses restricted to chromosome 24 (NC_031792.1).
- Reads (NCBI SRA): SRR3085443 and SRR3085445 (reference donor, merged), SRR5423293–SRR5423296 (four focal individuals).
- Reads were mapped with BWA-MEM2 v2.0pre2 as part of the pipeline described in Chen et al. (2025), sorted with SAMtools 1.13 and marked for duplicates; CRAM conversion with SAMtools 1.20.
- Mean depth on chromosome 24 after duplicate removal: 37–61× (focal individuals), 21× (reference donor).

## Reproducing the results

```bash
# Figure 2
cd simulations
python3 simulate_mosaic.py      # ~10 s; seeded
python3 summarise_sim.py
python3 plot_mosaic.py

# Figures 3-4, Table 3, Supplementary Table S1
cd ../empirical
# edit config.sh first
bash 01_call_variants_freebayes.sh
bash 02_popgen_stats.sh
python3 03_plot_figure3_triallelic.py realdata_freebayes_highcov_NC_031792.1
python3 04_plot_figure4_configurations.py realdata_freebayes_highcov_NC_031792.1 H1 SRR5423293
python3 05_pi_with_without_multiallelic.py realdata_freebayes_highcov_NC_031792.1
```

## Software

See `environment.yml`. Tested with msprime 1.4.5 and Python 3.
