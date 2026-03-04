# 3' UTR Length Calculation: Genomic Distance vs. Spliced Length


### 1. Genomic Distance (Approximate / Positional)

Formula: `coord (last_exon.end) - cds_end` 
Logic: Measures the total span on the DNA template, including non-coding intervening sequences.
Logic Diagram:
```

Genomic Axis: 100 ------------------------------------------------------> 500

               |  Exon A  |        Intron (Non-coding)    |  Exon B  |

Structure:     [██████████]-------------------------------[██████████]

Content:       |-- CDS --|                                |-- UTR --|

Coordinates:   ^         ^                                          ^

             start    cds_end                                     coord

                      (150)                                       (450)

```

• Result: 300 bp (includes 200 bp Intron).

• Significance: Highlights PolyA Site Choice on the chromosome. In VAE, this amplifies the separation between distal and proximal sites.

---

### 2. Spliced mRNA Length (Precise / Molecular)

Formula: `sum(u.end - u.start + 1 for u in utr_features)` 
Logic: Sums only the segments that remain after splicing (the mature transcript).
Logic Diagram:
```

Mature mRNA:   [██████████][██████████]

Content:       |-- CDS --||-- UTR --|

```

• Result: 100 bp (Exon A tail + Exon B).

• Significance: Reflects the actual Physical Molecule. Critical for calculating miRNA binding density and mRNA decay rates.
### Convergence Case: CDS & PolyA in the Same Exon
In most "housekeeping" scenarios, Genomic Distance equals Spliced Length.
Genomic Axis: 100 --------------------------------------> 300
               | Exon B (last exon) |
Structure:          [████████████████████████]
Content:            |-- CDS --|-- 3' UTR ----|
Coordinates:        ^         ^              ^
                  start    cds_end         coord (250)
                           (150)

calulate: coord (250) - cds_end (150) = 100 bp 
---
Verification Checklist (Action Plan)
IGHG1 (Membrane-bound): Verify if the TM1/TM2 splicing creates a multi-exon 3' UTR.

ACTB/GAPDH: Confirm if these standard controls align with the "Convergence Case."

Data Integrity: Ensure gffutils database includes three_prime_UTR features; otherwise, fall back to Genomic Distance with a warning flag.

Note: IsoDecipher records both metrics to maintain "positional awareness" while respecting "molecular reality."