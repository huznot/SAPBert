# Open questions for Dr. Lix

Two items need PI input before they can be resolved, plus one
scope-changing decision. A previously-flagged discrepancy over the baseline
numbers has since been **resolved** and is kept at the bottom for the
record.

## RESOLVED -- baseline ICD-9-CM -> ICD-10-CA numbers reproduce exactly

**No PI input needed; nothing to action. Recorded because this file
previously asked you to adjudicate it.**

The reference numbers were correct all along. The problem was the parameter
grid: `scripts/02_run_comparison.R` swept top-N {3, 5, 10, 15}, but on the
ICD-9-CM -> ICD-10-CA track the true optimum for both models sits at
**top_n = 30**, outside that range. The search could not reach the optimum,
so it kept reporting a slightly worse combination (F1 0.423/0.523) and the
reference numbers looked unreproducible.

Both tracks now sweep top-N {3, 5, 10, 15, 20, 25, 30} and all four
reference numbers reproduce exactly:

| Track | Model | Reference | Reproduced | At |
|---|---|---|---|---|
| ICD-9 -> ICD-10-CA | ClinicalBERT | 0.427 / 0.271 | 0.427 / 0.271 | thr 0.995, top_n 30, flag 4 |
| ICD-9 -> ICD-10-CA | SapBERT | 0.530 / 0.360 | 0.530 / 0.360 | thr 0.95, top_n 30, flag 4 |
| ICD-9 -> ICDA-8 | ClinicalBERT | 0.716 / 0.557 | 0.716 / 0.557 | thr 0.999, top_n 5, flag 3 |
| ICD-9 -> ICDA-8 | SapBERT | 0.821 / 0.697 | 0.821 / 0.697 | thr 0.99, top_n 3, flag 2 |

Confirmed independently by `scripts/07_full_grid_comparison.R`, which
reaches the same values by a different code path.

**Correcting the record:** an earlier version of this file stated that the
widened range "does not reproduce the claimed numbers." That was wrong. It
was based on a re-run that had been interrupted partway -- it completed only
1 of 4 thresholds for the first of four model/track combinations -- whose
partial output was mistaken for a finished result. The completed sweep says
the opposite. Apologies if that reached you before this correction did.

## 0.5. Should round-trip consistency filter mappings, or just be reported?

`scripts/05_bidirectional_and_roundtrip.R` / `report.Rmd` currently use
round-trip agreement (does a target code that an ICD-9-CM code maps to, map
back to that same ICD-9-CM code?) as a **reported metric only** -- it does
not change which mappings the pipeline outputs. An alternative is to use it
as an active filter (only accept a forward mapping when the reverse
direction agrees), which would change the actual mapping output, not just
its score. Not guessed at since the two options have materially different
consequences. Please confirm which is wanted with Dr. Lix.

## 1. "Removing certain codes because cosine similarity drops"

Not implemented -- skipped per instructions ("do not implement a guess at
this"). Unknowns that need clarifying before this can be built:

- **Which codes?** A specific list, a code range/chapter, or a rule (e.g.
  "codes with fewer than N characters of label text", "codes whose best
  cosine similarity falls below some floor")?
- **"Cosine similarity drops"** -- drops relative to what? Compared across
  the two embedding models (ClinicalBERT vs SapBERT) for the same code?
  Compared to that code's own historical/expected similarity? A drop after
  the filler-word stripping in task 1 (i.e. codes whose top match changes
  once filler words are removed)?
- **Remove from what stage?** From the label text going into the embedding
  model, from the candidate similarity matrix before thresholding, or from
  the manual validation set entirely?
- **Is this related to task 1?** It's plausible this is about filler-word
  stripping making some short labels nearly empty (e.g. a label that's
  almost entirely filler words) and their embeddings becoming degenerate --
  but that's a guess, not something stated.

Please clarify with an example code or two and this can be scoped properly.

## 2. Expanded validation dataset (task 4)

No additional manually-verified ICD-9/ICD-10-CA or ICD-9/ICDA-8 crosswalk
pairs were found anywhere in this repository -- `data/original` contains
exactly the 937-pair and 331-pair sheets already used for the existing
baseline (`Validation_Data .xlsx`, sheets `Validation_ICD9_ICD10` and
`Validaion_ICD9_ICD8`), confirmed by row count. Task 4 is not implemented
because there's no larger validation set to run it against yet.

If Dr. Lix has additional verified pairs (a file, a shared drive location,
a database export), point this pipeline at them and the validation-set
expansion described in the task can be run: same pipeline, same metrics,
comparing old (937/331-pair) vs new results, explicitly reporting any drop
as expected/informative rather than a bug, per the task instructions.
