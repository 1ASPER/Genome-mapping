# Genome read mapping: E. coli K-12 MG1655

Reference: `GCF_000005845.2_ASM584v2` (4,641,652 bp).
Reads: `ERR022075_1.fastq`, 22,720,100 single-end reads of 100 bp. The full dataset was processed.


## Algorithms used
- Suffix array [prefix doubling]
- Binary search
- Banded DP [edit distance]
- Reverse complement search
- Difference array
- Pigeonhole seeding

## Definition of "mapped"

A read is mapped if there is a position on the forward or reverse strand where the edit distance (substitutions, insertions, deletions) between the **whole read** and a genome substring is **at most 9**. There is no clipping, and `N` counts as an error.

Confidence in optimality: by the pigeonhole argument the search is exhaustive for this threshold, and the reported distance is the minimum one. On small synthetic genomes the output was compared with a full DP over the whole genome on both strands for 6,000 reads, with no differences.

## Results

| Metric | Value |
|---|---|
| Total reads processed | 22,720,100 |
| Mapped | 22,458,693 (98.85%) |
| Unmapped | 261,407 (1.15%) |
| Uniquely mapped | 22,041,477 (97.01% of all, 98.14% of mapped) |
| Multi-mapped | 417,216 (1.84% of all, 1.86% of mapped) |
| Forward / reverse strand | 11,237,053 / 11,221,640 |

A read is unique if exactly one locus has the best edit distance, and multi-mapped if several non-overlapping loci share it. A multi-mapped read is assigned to one of its best loci (read index modulo number of loci).

## Alignment quality

| Edit distance | Reads | % of all | Cumulative % |
|---|---|---|---|
| 0 | 19,387,428 | 85.33 | 85.33 |
| 1 | 1,927,089 | 8.48 | 93.81 |
| 2 | 489,377 | 2.15 | 95.97 |
| 3 | 237,886 | 1.05 | 97.01 |
| 4 | 142,497 | 0.63 | 97.64 |
| 5 | 94,300 | 0.42 | 98.06 |
| 6 | 66,156 | 0.29 | 98.35 |
| 7 | 49,058 | 0.22 | 98.56 |
| 8 | 36,533 | 0.16 | 98.72 |
| 9 | 28,369 | 0.12 | 98.85 |

- Mean edit distance of mapped reads: 0.26; mean identity: 99.74%.
- The table is also the mapping rate for every threshold from 0 to 9. The last step adds 0.12% of reads, so the threshold sits on the plateau.
- Mean base quality: 36.37 for mapped reads, 18.13 for unmapped reads. Unmapped reads are mostly low-quality reads.

## Genome read coverage

| | All mapped reads | Unique reads only |
|---|---|---|
| Bases covered >= 1x | 4,640,887 (99.98%) | 4,567,691 (98.41%) |
| Bases covered >= 10x | 4,638,160 (99.92%) | 4,564,057 (98.33%) |
| Mean depth | 483.77x | 474.78x |
| Median depth | 494x | 491x |
| Max depth | 1032x | 1032x |

The 1.6% of the genome not covered by unique reads are repeated regions, where a 100 bp read cannot be placed unambiguously.

## Performance

MacBook Air M4, single thread: suffix array built in 0.82 s, total run 174.51 s (CPU time), about 131,000 reads/s.

## Limitations

- No soft clipping: reads with adapter sequence or a badly damaged tail fail the whole-read threshold.
- The genome is circular, but reads spanning the junction between its end and start are not handled.
- Reads whose best alignment needs more than 9 edits are reported as unmapped.

## How to run

Put `A.cpp`, `GCF_000005845.2_ASM584v2_genomic.fna` and `ERR022075_1.fastq` in one folder, then:

```bash
g++ -O2 -std=c++17 A.cpp -o A
./A > report.txt
```

Progress is printed to stderr, the statistics to stdout.
