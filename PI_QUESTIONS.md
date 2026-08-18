# Open questions for Dr. Lix

Two items from the task list are **not implemented** because implementing
them would require guessing at specifics the PI didn't give. Flagging both
here rather than picking an interpretation.

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
