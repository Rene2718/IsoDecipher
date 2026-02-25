# Algorithmic Optimization: Resolving Group-Level Ambiguity in Isoform Assignment
## Technical Deep Dive: Refining the Isoform Assignment Logic

### 1. The Discovery: Identifying "False Ambiguity"
While reviewing the initial algorithmic output for genes with complex transcriptomic architectures (e.g., BCL2), I identified a logic flaw. The system was flagging high-confidence reads as ambiguous, leading to a loss of usable data. As the developer, I realized that while the code was syntactically "correct," it lacked the biological nuance required for functional grouping.

### 2. Root Cause Analysis (Code Audit)
The initial implementation used a `list` structure for matching features.

**The Conflict:** Since the 3' panel contains multiple transcript entries per gene, a single read often maps to identical coordinates shared by multiple transcripts within the same group. A standard list-based match triggered len(matches) > 1, wrongly categorizing these as ambiguous.

### 3. The Solution: Set-Based Deduplication
I directed the transition from a transcript-level list to a set-based grouping strategy. By leveraging the unique-element property of sets, I ensured that redundancy is resolved at the group level during the scan.

**Optimized Logic:**

```Python
# Transitioned to set-based matching to enforce group-level uniqueness
tier1_matches = set()

for feature in panel[gene]["tier1"]:
    if feature["chrom"] == chrom and feature["strand"] == strand:
        if feature["start"] <= pos <= feature["end"]:
            tier1_matches.add(feature["group"])

# Result: len(tier1_matches) now correctly reflects distinct functional groups
if len(tier1_matches) == 1:
    group = next(iter(tier1_matches))
```
### 4. Impact & Engineering Rigor
Expert Oversight: This fix demonstrates the necessity of domain expertise in AI-assisted development. I identified a logic gap that a general-purpose AI would overlook.

Accuracy: Significantly increased the "Unique Assignment Rate" for 3’ tag data.

Next Phase: I am currently developing unit tests to simulate these edge cases (e.g., overlapping vs. distinct groups) to ensure long-term stability and catch potential regressions.