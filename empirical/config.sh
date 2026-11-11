# Edit these paths before running 01_call_variants_freebayes.sh
CHR="NC_031792.1"     # great tit chromosome 24

# Reference assembly Parus_major1.1 (GCF_001522545.3)
REF="/path/to/GCF_001522545.3_Parus_major1.1_genomic.fna"

# Sorted, duplicate-marked alignments (CRAM or BAM)
REFIND="/path/to/SRR3085443_sorted_dedup.cram"   # reference donor
H1="/path/to/SRR5423293_sorted_dedup.cram"       # focal 1
H2="/path/to/SRR5423294_sorted_dedup.cram"       # focal 2
H3="/path/to/SRR5423295_sorted_dedup.cram"       # focal 3
H4="/path/to/SRR5423296_sorted_dedup.cram"       # focal 4
