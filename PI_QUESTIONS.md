# Open questions for Dr. Lix

Three items need PI input before they can be resolved: two task-list items
that would require guessing at specifics the PI didn't give, and one data
discrepancy discovered while verifying the stated baseline. Flagging all
three rather than picking an interpretation.

## 0. Baseline ICD-9-CM -> ICD-10-CA numbers are not reproducible by any code in this repo

The task's stated baseline was F1 0.427/Acc 0.271 (ClinicalBERT) and
F1 0.510-0.530/Acc 0.342-0.360 (SapBERT). The repo's most recent commit
before this work (143b1ff) hand-edited `results/best_by_model.csv` to show
F1 0.427/0.530 and updated `report.Rmd`'s narrative to describe a widened
top-N co-occurrence range ({5,10,15,20,25,30}) as the source -- but
`scripts/02_run_comparison.R` itself was never changed to actually sweep
that wider range; it still used {3,5,10,15} for both tracks.

I fixed the script to actually sweep top-N up to 30 for that track and
re-ran the full grid to check. **The wider range does not reproduce the
claimed numbers.** The best combination is unchanged (F1 0.423/0.268
ClinicalBERT at top_n=15, F1 0.523/0.355 SapBERT at top_n=10) -- none of
the added top-N values (20, 25, 30) beat what {3,5,10,15} already found.
So 0.427/0.271 and 0.510-0.530/0.342-0.360 are not reproducible by any
code currently in this repository, under either grid.

`results/best_by_model.csv` and `report.Rmd` now report the honest,
reproducible numbers (0.423/0.268 and 0.523/0.355) instead of forcing a
match to the unreproducible hand-edited ones, per instruction. This needs
Dr. Lix's input: either the original 0.427/0.530 numbers came from a run
that used different parameters/data not present in this repo (a different
flag rule, a different embedding source, a different validation subset),
or they were a mistake in the earlier hand-edit. Worth asking directly
rather than guessing further.

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
