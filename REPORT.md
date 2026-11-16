# Genome read mapping: E. coli K-12 MG1655

Reference: `GCF_000005845.2_ASM584v2`, 4,641,652 bp.
Reads: `ERR022075_1.fastq`, 22,720,100 single-end reads of 100 bp. I processed the full dataset.


## My algorithms

- Suffix array of the reference (prefix doubling)
- Binary search over the suffix array
- Pigeonhole seeding
- Banded DP for edit distance
- Reverse complement search
- Difference array for coverage

Each read is cut into k+1 pieces. Every piece is searched exactly in the suffix array, on both strands, and each hit is checked with a banded DP. This is done for k = 3, 6, 9, and a read goes to the next pass only if it is still unmapped.


## What is mapped?

A read is mapped if somewhere in the genome, on either strand, the edit distance between the whole read and a genome substring is at most 9. Nothing is clipped, and `N` counts as an error.

No such read is missed: with at most k edits and k+1 pieces, one piece has no edits and is found exactly, and the DP band of ±k around it covers the whole alignment.

A read is unique if one place has the best distance, and multi-mapped if several do.


## Results

```
genome length: 4641652
threshold: edit distance <= 9 over the whole read
total reads: 22720100
mapped reads: 22458693 (98.85%)
unmapped reads: 261407 (1.15%)
unique: 22041477 (97.01% of all, 98.14% of mapped)
multi-mapped: 417216 (1.84% of all, 1.86% of mapped)
forward / reverse: 11237053 / 11221640

edit distance distribution (mapped reads):
ed = 0: 19387428 (85.33% of all, cumulative 85.33%)
ed = 1: 1927089 (8.48% of all, cumulative 93.81%)
ed = 2: 489377 (2.15% of all, cumulative 95.97%)
ed = 3: 237886 (1.05% of all, cumulative 97.01%)
ed = 4: 142497 (0.63% of all, cumulative 97.64%)
ed = 5: 94300 (0.42% of all, cumulative 98.06%)
ed = 6: 66156 (0.29% of all, cumulative 98.35%)
ed = 7: 49058 (0.22% of all, cumulative 98.56%)
ed = 8: 36533 (0.16% of all, cumulative 98.72%)
ed = 9: 28369 (0.12% of all, cumulative 98.85%)
mean edit distance: 0.26
mean identity: 99.74%
mean base quality: mapped 36.37, unmapped 18.13

coverage (all mapped reads):
covered bases (>= 1x): 4640887 / 4641652 (99.98%)
covered bases (>= 10x): 4638160 (99.92%)
mean depth: 483.77x, median 494x, max 1032x

coverage (unique reads only):
covered bases (>= 1x): 4567691 / 4641652 (98.41%)
covered bases (>= 10x): 4564057 (98.33%)
mean depth: 474.78x, median 491x, max 1032x

index time: 0.82 s
total time: 174.51 s
```


Mapping rate: 98.85% of reads are mapped, 97.01% uniquely and 1.84% to several places.

Alignment quality: 85.33% of reads match exactly, and the mean identity is 99.74%. The last threshold step, from 8 to 9, adds only 0.12% of reads. Unmapped reads have much lower base quality (18.13 against 36.37).

Coverage: 99.98% of the genome is covered, with a mean depth of 483.77x. Unique reads alone cover 98.41%; the rest are repeats.


Speed is 174.51 s for all reads on one thread on MacBook Air M4.