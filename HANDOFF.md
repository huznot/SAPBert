# handoff, icd crosswalk project

paste this into a fresh claude code session on the other machine. it has no
memory of the previous work.

---

## 0. read this part first, it is about how to work

these are rules the user has stated repeatedly. follow them without being asked.

**commit messages**: one short line, all lowercase, like a person in a hurry.
no capitalised subject, no multi paragraph body, no tables of numbers, no
"this is X not Y" phrasing, no em dashes. examples of what is wanted:
```
fix filler word list, was merging two codes
add learning curve, more data doesnt help
cleaning up scripts
```
never add a `Co-Authored-By` trailer. the whole history was rewritten twice to
remove ai sounding messages and trailers, do not reintroduce them.

**code comments**: all lowercase, short, written like a normal programmer.
no essay blocks, no "IMPORTANT:", no "NOTE THAT", no em dashes. the scripts
were all cleaned once already, match the existing style.

**no personal names anywhere** in code, comments, or docs. refer to the gpt-4
paper by what it did, not by author.

**writing style in chat**: the user gets annoyed by long answers. be short, get
to the point, do not over explain. when they ask a direct question give a
direct answer first.

**the user is a high school student** working in a university lab. the repo is
public. everything needs to read as their own work.

---

## 1. what the project is

automating icd diagnosis code crosswalks. when a health system changes coding
standards, every old record becomes hard to use, and someone has to hand build
a mapping from each old code to the new one. this automates most of it.

two tracks, both at the three character code level:
- icd-9-cm to icd-10-ca, 937 manually verified pairs over 345 codes
- icd-9-cm to icda-8, 331 pairs over 302 codes

repo: https://github.com/huznot/SAPBert, branch main

`CHANGELOG.md` is the main writeup now. it has an abstract plus a numbered
list of every change from the original 2024 pipeline to now, with the tables.
the old `report.Rmd` was deleted, the pi asked for a bullet list instead of a
formal report. keep the changelog current, it is what gets shown.

`icd_crosswalk.Rproj` exists so the repo opens as an rstudio project, which
sets the working directory to the repo root.

---

## 2. where it stands, current numbers

all held out, five fold cross validation grouped by icd-9 code.

| | rule based (held out) | two stage system |
|---|---|---|
| icd-9 to icd-10-ca | 0.546 f1 | **0.668** |
| icd-9 to icda-8 | 0.824 f1 | **0.840** |

top-1 accuracy (is the highest scoring target correct): 93.3% and 88.4%.

at 95% precision: 37.5% and 76.1% of the manual crosswalk pairs are built
automatically, covering 86.1% and 84.1% of source codes. review burden is
about 14% of codes on both tracks.

### the main finding

the embedding model was never the bottleneck. every correct answer was already
present in the similarity matrix. the old filtering step kept only candidates
within 5% of each code's best similarity, which let through 37% of correct
mappings, and a filtered candidate can never be recovered. that capped the
whole pipeline at f1 0.770 no matter which model was used. replacing it with
top-k retrieval raised the ceiling to 0.96.

### other findings, all measured

- **filler word stripping** helps clinicalbert, hurts sapbert. not adopted.
  the word list was also merging distinct codes ("other abdominal hernia" and
  "unspecified abdominal hernia" both became "abdominal hernia"). fixed, and
  15_ now fails if any word list merges codes.
- **standard stopword libraries are unsafe here**. the big ones treat single
  letters as stopwords so "hepatitis a" and "hepatitis b" collapse. snowball
  and nltk are clean.
- **mpnet (no medical training) matches sapbert** on icd-10-ca (0.533 vs
  0.534), loses on icda-8 (0.769 vs 0.821). both beat clinicalbert everywhere.
- **more training data would not help**. learning curve plateaus at roughly
  100 to 180 training codes on both tracks, slope slightly negative after.
- **generalises to unseen clinical areas**. holding out whole ccs categories
  instead of random codes costs -0.004 and -0.000, inside run to run noise.
- **ablation**, f1 change when each feature group is dropped:

  | dropped | icd-10-ca | icda-8 |
  |---|---|---|
  | mutual forward+reverse | -0.063 | -0.011 |
  | chapter | -0.051 | -0.003 |
  | co-occurrence | -0.051 | -0.036 |
  | lexical | -0.015 | -0.001 |
  | clinicalbert | -0.012 | +0.003 |
  | rank:ndcg instead of binary | -0.023 | +0.008 |

- **the weak spot is multi target codes**. icd-10 split many icd-9 categories
  into several codes (icd-9 250 diabetes becomes E10, E11, E13, E14). for
  codes with one correct target the system solves 86.7% (icd-10-ca) and 88.2%
  (icda-8). for codes with several it solves 6.0% and 0.0%. the difference
  between the two tracks is entirely composition: 63.2% of icd-10-ca codes are
  multi target vs 7.9% for icda-8.
- **27.7% of icda-8 codes have a perfect similarity match** (identical label
  text), against 0% on icd-10-ca. measured in 25_. this is the cleanest single
  number for why the icda-8 track scores higher.
- **84 icd-9 codes (24%) have no icda-8 co-occurrence data at all**, vs 7 on
  icd-10-ca. the 1972-79 source data is much thinner. also from 25_.
- **69.8% of icda-8 pairs have a target identical to the source code**, so much
  of that track is identity mapping and 0.840 is easier than it looks. this
  belongs in any writeup.
- **code numbers in the embedding text hurt icd-10-ca**. we were embedding
  "250 diabetes mellitus". removing the number took top-1 from 79.7% to 90.1%
  on sapbert and 80.9% to 89.9% on mpnet. it slightly hurt icda-8 (-0.023,
  -0.033) because those code numbers often match exactly.

---

## 3. environment, this will waste an hour if missed

- **rscript is not on PATH**. use the full path. on the old machine it was
  `"/c/Program Files/R/R-4.6.1/bin/Rscript.exe"`. find the equivalent.
- **python for embeddings** was a portable install at `~/Downloads/pyembed/python.exe`
  with torch and transformers. only needed to regenerate embeddings.
- **r packages**: dplyr, tidyr, readxl, writexl, stringr, purrr, ggplot2,
  jsonlite, xgboost, stopwords. xgboost and stopwords
  were installed from cran during the work.
- **never edit an .R file while an Rscript process is running it**. rscript
  parses incrementally and a mid run edit kills it with a fake syntax error.
  this cost a full run once.
- **windows filesystem is case insensitive**. two output filenames differing
  only by case are the same file, and parallel jobs silently clobber each
  other. this corrupted results once. 07_ asserts tags are unique
  case insensitively.
- **bash heredocs mangle backslashes and apostrophes** in this setup. use the
  Write tool for scripts rather than `cat <<EOF`.
- long jobs: run as background processes writing to distinct output files.
  8 cores available.
- **force push is blocked by the sandbox**. the user has to run it themselves.

---

## 4. scripts, what each does

most run from inside `scripts/` and use `../data` paths. the newer ones (21,
24, 25) run from the repo root instead. check the header comment before running
one. worth normalising at some point.

| file | what |
|---|---|
| 01 | generates the original sapbert matrices, rarely needed |
| 02 | original rule based pipeline, full parameter grid |
| 03 | charts from 02 |
| 05 | reverse direction and round trip consistency |
| 06 | superseded by 07, kept for reference |
| 07 | full grid for all 8 embedding conditions, can run conditions in parallel |
| 08 | combines 07 outputs into summary tables |
| 09 | error decomposition, where correct answers are lost, oracle ceiling |
| 10 | candidate generation study, relative threshold vs top-k |
| 11 | builds candidate pool and feature table, writes rerank_features_<track>.rds |
| 12 | trains reranker, grouped cv, evaluates baseline under same protocol |
| 12b | merges per track outputs when 12 is run one track at a time |
| 13 | precision at coverage, triage split |
| 14 | produces the actual crosswalk for codes with no known answer |
| 15 | stopword collision check, fails if a word list merges codes |
| 16 | feature ablation and ranking objective test |
| 17 | retrieval config sensitivity |
| 18 | learning curve |
| 19 | category holdout |
| 20 | top-1 accuracy |
| 21 | performance by single vs multi target code |
| 22 | do multi target sets cluster in blocks |
| 23 | code number in embedding text, help or hurt |
| 24 | worked example of the similarity matrix, for explaining the method |
| 25 | frequency distributions, max similarity and top co-occurrence |
| pipeline_lib.R | shared functions |
| verify_vectorized_equivalence.R | gate on the chapter lookup speed refactor |

typical full rerun:
```
Rscript 11_rerank_features.R
Rscript 12_cv_rerank.R 10_9 &
Rscript 12_cv_rerank.R 8_9 &
Rscript 12b_merge_cv_results.R
Rscript 13_precision_coverage.R
```
about 10 minutes with both tracks in parallel.

---

## 5. next steps, in order

the pi sent a marked up copy of the original methods document on 24 aug 2026
with eight comments. those comments are now the work queue and they override
the older priorities below them. the file is
`Methods_Results_Automatic_Mapping_With Comments_24August2026.docx` in the repo
root, gitignored via `*.docx`, do not commit it.

the document body is byte identical to the 2024 version, verified by diff. all
new content is in the margin comments.

### what the pi settled

- **no bigger code set, ever.** "we cannot do this, because the co-occurrence
  (i.e., empirical) data is limited to co-occurrences of three-digit codes."
  this closes the cihi crosswalk request. do not raise it again. the practical
  consequence: any expansion beyond 354 codes arrives with co-occurrence
  missing, and ablation says co-occurrence is worth -0.051 and -0.036 f1.
- **chronic disease ccs categories only** is deliberate, it builds on a prior
  study. not a limitation to apologise for in a writeup.
- **the odd per track top-n ranges are unexplained on their side too.** she
  asked the previous analyst and got no answer. so our finding that the
  icd-10-ca optimum sits at top_n=30, past where the original grid stopped at
  10, stands unchallenged.

### 1. blocked, needs sign off before running

comments 3 and 4 ask for clinicalbert re-runs in four combinations:

| | code in text | code removed |
|---|---|---|
| stopwords kept | already have this, f1 0.427 / 0.716 | needs running |
| stopwords removed | needs running | needs running |

she explicitly said: "please identify a suitable standardized dictionary and
confirm with me before proceeding to apply it." **do not run these until she
picks a dictionary.**

the evidence for that choice already exists in
`results/stopword_collision_check.csv`:

| dictionary | words | codes merged |
|---|---|---|
| nltk english | 120 | 0 |
| snowball | 175 | 0 |
| smart | 571 | 17 |
| stopwords-iso | 1298 | 20 |

recommend **nltk english**: zero collisions, most widely cited, easy to defend
in a methods section. snowball is equally safe. smart and stopwords-iso are
disqualified, they strip single letters so "hepatitis a" and "hepatitis b"
become the same string.

note this is a different question from the old filler word list. the filler
list was hand written and aimed at redundant clinical phrasing. she wants a
standardized off the shelf dictionary applied to clinicalbert specifically, to
replicate the original result under a documented cleaning rule.

when she approves, the mechanics are already built:
- `generate_embeddings.py` has `--no-code` for dropping the code number
- stopword removal needs a new `--clean stopwords` mode wired to the chosen
  dictionary, alongside the existing `base` and `stripped`
- then run `07_full_grid_comparison.R` on the new conditions and `08_` to
  assemble

also ask her to confirm **354 vs 365**. she wrote "N = 365 ICD-9-CM codes"
twice in comments 5 and 6, but the document body says 354 throughout and 354 is
what is in the data. probably a typo, but if she means a different code set
every denominator changes.

### 2. done, comments 5 and 6

`25_frequency_distributions.R`, run from the repo root. outputs four files to
`results/`. numbers as of 24 aug 2026:

max cosine similarity per icd-9 code, clinicalbert:

| | median | q1 to q3 | range | identical label matches |
|---|---|---|---|---|
| icd-10-ca | 0.930 | 0.914-0.948 | 0.853-0.981 | 0 |
| icda-8 | 0.968 | 0.945-1.000 | 0.865-1.000 | 98 (27.7%) |

most frequent co-occurring target per icd-9 code:

| | median | max | codes with no co-occurrence data |
|---|---|---|---|
| icd-10-ca | 1,296 | 238,500 | 7 |
| icda-8 | 326 | 19,100 | 84 (24%) |

two things in here are worth raising with her unprompted. 98 icda-8 codes have
a perfect similarity match, meaning an identical label, against zero on
icd-10-ca. that is the identity mapping story quantified and it explains most
of why the icda-8 track scores higher. and 84 icd-9 codes have no icda-8
co-occurrence data at all, since the 1972-79 source data is much thinner.

gotcha if this gets rerun or extended: 33 icda-8 similarities come back at
1.0000007 from floating point, so `cut()` with a break at exactly 1 silently
drops them as NA and the percentages stop summing to 100. the script clamps
with `pmin(x, 1)` before binning. keep that.

### 3. label only embeddings through the full pipeline

still the highest certainty win, and it now overlaps with what the pi asked
for in comment 4, so do it in the same pass. the embeddings already exist:
```
data/generated/cosine_similarity_matrices_{10_9,8_9}_{sapbert,mpnet}_base_nocode.xlsx
```
they were only tested at the retrieval stage (23_), where they gained about 10
points of top-1 on icd-10-ca. the reranker has never seen them.

what to do: point the `sim` paths in `11_rerank_features.R` TRACKS at the
`_nocode` files for sapbert and mpnet, rebuild features, rerun 12 and 13,
compare against 0.668 / 0.840.

careful: icda-8 got slightly worse with label only, so either keep the code
number for that track, or better, add "does the source code string equal the
target code string" as a feature so icd-10-ca gets clean labels and icda-8
keeps the identity signal.

### 4. block features for multi target codes

biggest remaining headroom. multi target completion is 6% and everything else
is near its ceiling. 22_ showed 86.2% of multi target sets share a letter
prefix on icd-10-ca, 100% on icda-8.

a crude post hoc expansion rule was tested and it raises f1 (0.594 to 0.687 at
a fixed operating point) but drops precision from 0.951 to 0.847, which breaks
the auto accept guarantee. also the span and floor were chosen by looking at
the results, so that number is selection on test.

do it properly: add features to `11_rerank_features.R` so the model learns when
siblings are valid rather than always assuming it.
- is this candidate in the same letter block as this code's top scoring candidate
- numeric distance to this code's top scoring candidate
- how many pool candidates share that block
- does the label share the block's head term

then rerun 12 and 13 and let the thresholds tune inside folds.

### 5. small cleanups

- drop clinicalbert from the ensemble, the ablation says it earns nothing and
  it slightly helps icda-8 to remove it. note this is about the *ensemble*, the
  pi still wants clinicalbert re-run standalone for comments 3 and 4.
- use `rank:ndcg` on icda-8 only, worth +0.008 there, but it costs -0.023 on
  icd-10-ca so it must be per track.

## 6. the pi requests, resolved

an earlier email asked for two things. both are now answered, see section 5.

1. **the cihi crosswalk / bigger code set: refused, permanently.** the reason is
   that co-occurrence data only exists at three digits. this is not a matter of
   asking again in a different way.
2. **wider co-occurrence coverage: same answer, same reason.** the label files
   have 354 icd-9 codes and the co-occurrence tables have 347, fully
   consistent, so there is no hidden data in the repo to unlock either.

so 354 codes is the permanent ceiling for this project. write the limitations
section accordingly rather than treating it as pending.

### context on the gpt-4 paper

same lab, phd student is first author, pi is senior author. it evaluated gpt-4
on icd-10-ca to icd-9-cm translation across nine prompting strategies and got
48.3% full code accuracy, 84.7% at three characters, and concluded performance
was insufficient for deployment.

**do not present our numbers as beating it.** the direction is reversed, the
reference standard differs, the code subsets differ, and their crosswalk is one
to one where this task averages more than one target per code. a purpose built
system beating zero shot gpt-4 is expected, not a result. the honest framing is
complementary: they showed general llms reach about 85% on broad categories but
cannot say which answers to trust, this adds calibrated confidence and runs
without sending health data to an external api.

they also flagged that public crosswalks may be in gpt-4's training data, so
some of its accuracy may be memorisation. our co-occurrence counts come from
local records and cannot leak that way. that is a real methodological point in
our favour.

the pi was also told not to pursue an llm angle directly, a phd student in the
lab has that project.

## 7. things that were wrong and got corrected

worth knowing so they are not reintroduced.

- an earlier session claimed the baseline numbers were not reproducible. that
  was from a run interrupted after 1 of 4 thresholds. they reproduce exactly
  once the top-n grid goes to 30.
- a "+0.034 free accuracy from fixing retrieval" claim was wrong. the ablation
  script computed the truth set after retrieval, which dropped missed positives
  out of the recall denominator and inflated f1. fixed, and 17_ confirmed there
  is no gain.
- coverage was reported as 43% / 79%, measured against pairs reachable in the
  pool rather than the full crosswalk. correct figures are 37.5% / 76.1%.
- "57% sent for human review" was wrong, it conflated pair incompleteness with
  review load. actual review burden is about 14% of codes.

---

## 8. honest assessment

the work is rigorous for its size: held out grouped cv, ablation, learning
curve, category holdout, oracle ceiling diagnosis, precision at coverage. more
methodological care than a lot of published work in this area.

what limits it: 354 codes from one institution, no external validation, and
retrieve-then-rerank is a standard information retrieval architecture. the
contribution is the diagnosis, the evaluation protocol, and the application,
not the architecture. say that plainly rather than letting a reviewer say it.

f1 0.668 on icd-10-ca is moderate. the top-1 number of 93.3% is the more
flattering and arguably more relevant framing for a crosswalk tool, but the two
measure different things and should not be mixed in one sentence.
