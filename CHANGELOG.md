# Change log

Everything I changed, from the original project to now.

## Abstract

The goal is to automatically match old ICD-9 diagnosis codes to newer ICD-10-CA
codes and older ICDA-8 codes, so hospitals don't have to do it by hand.

The original version worked like this. Take each code's text label, turn it into
a list of numbers using ClinicalBERT (a model that reads medical text), and
compare those numbers to find labels that mean similar things. Keep the closest
matches. Also keep codes that show up together a lot in real patient records.
Throw out any pair from unrelated chapters of the ICD book. Then output whatever
passed one of three rules. Scored 0.427 F1 on ICD-9 to ICD-10-CA and 0.716 on
ICD-9 to ICDA-8. (F1 is a score out of 1 that balances "did we find the right
answers" against "did we also output wrong ones".)

I rebuilt that in R first and checked it gave the same numbers. Then I tried a
better model, SapBERT, which pushed ICD-10-CA from 0.427 to 0.530. But when I
looked at where the mistakes were actually happening, the model wasn't the real
problem. Every correct answer was already sitting in the similarity table. The
filtering step was throwing 37% of them out before anything got ranked, and once
a code is thrown out it can never come back. So no matter how good the model
got, the pipeline was capped at 0.770.

So I widened the filtering step to keep more options, and added a second model
that scores the survivors and picks the good ones. That got it to 0.668 on
ICD-10-CA and 0.840 on ICDA-8, tested on codes the model had never seen. Each
mapping also gets a confidence score now, so 86% can be auto-accepted at 95%
precision and a human only has to check the rest.

The errors that are left are almost all one thing. Some ICD-9 codes have one
correct answer, some have several, because ICD-10 split a lot of the old
categories into finer ones. Codes with one answer get solved 87-88% of the time
on both tracks. Codes with several answers get fully solved 6% of the time on
ICD-10-CA and 0% on ICDA-8. Since 63% of ICD-9 codes have multiple ICD-10-CA
answers but only 8% have multiple ICDA-8 answers, that one problem is basically
the entire difference between the two tracks.

---

## 1. Rebuilding it and checking it still works

- Rewrote the two Python notebooks in R. Same steps in the same order.
- The original scored against the CIHI crosswalk table, which isn't in the
  shared folder. I pointed it at the validation sheets that are shared instead.
- The original tested different numbers of co-occurring codes on each track:
  1 to 10 for ICD-10-CA, 5 to 30 for ICDA-8. I now test 3 to 30 on both. Turns
  out the best setting for ICD-10-CA is 30, which is past where the original
  stopped looking. That's why I couldn't reproduce 0.427 at first.
- Also widened the similarity cutoffs. The original tested 0.99 and up, plus
  1.0, which only keeps exact ties with the best score. I added 0.95.
- After those two fixes the R version reproduces the original numbers exactly:
  0.427 and 0.716.
- Rewrote the chapter lookup to work on the whole list at once instead of one
  row at a time. About 130x faster, same output, checked with a script that
  compares old against new.

## 2. Trying different models

- Wrote a Python script that generates the embeddings and similarity tables for
  any model, so every model gets built the same way and the comparison is fair.
- Tested every combination of settings (112 per track) on eight setups instead
  of one setting per model.

| model | ICD-10-CA F1 | ICDA-8 F1 |
|---|---|---|
| ClinicalBERT (original) | 0.427 | 0.716 |
| ClinicalBERT (rebuilt) | 0.430 | 0.716 |
| ClinicalBERT + filler words removed | 0.448 | 0.722 |
| SapBERT | 0.530 | 0.821 |
| SapBERT + filler words removed | 0.534 | 0.806 |
| all-mpnet-base-v2 | 0.533 | 0.769 |
| all-mpnet-base-v2 + filler words removed | 0.523 | 0.766 |

- SapBERT is the best single model, and it's a big jump over ClinicalBERT
  (+0.10 and +0.11).
- all-mpnet-base-v2 has no medical training at all and it ties SapBERT on
  ICD-10-CA. So being trained on medical text isn't what matters most here.
- Removing filler words helps ClinicalBERT but hurts SapBERT, so it's off by
  default.
- Found a bug where two setups whose names only differed by capital letters were
  writing to the same file on Windows, so jobs running at the same time
  overwrote each other. Fixed and reran all eight.

## 3. Testing filler word lists

- Compared my hand-written filler list against four standard ones (NLTK,
  Snowball, SMART, stopwords-iso).
- The big lists are dangerous here. SMART merges 17 codes and stopwords-iso
  merges 20, because they delete single letters and words like "with" and
  "without". So "hepatitis a" and "hepatitis b" become the same text and no
  model can tell them apart after that.
- Took `other`, `others`, `unspecified`, `with`, `without` out of my own list.
  Those five were merging 10 codes by themselves.
- Added a check that fails if the list ever merges two codes again.

## 4. Finding where the accuracy was actually going

Tracked every correct answer through the pipeline to see where it got lost.

| stage | ICD-10-CA | ICDA-8 |
|---|---|---|
| correct answer exists in the similarity table | 100% | 100% |
| makes it past the similarity cutoff | 37.1% | 79.5% |
| shows up in the co-occurrence list | 55.4% | 71.6% |
| kept as a candidate (either of the above) | 62.6% | 85.2% |
| actually output by the three rules | 42.9% | 78.5% |
| best F1 possible from what's left | 0.770 | 0.920 |

- The model isn't the bottleneck. Every answer is in the table. The cutoff
  throws out 37% of correct ICD-10-CA answers before anything gets ranked, and
  there's no way to get them back later.
- Ran a separate test on different ways of picking candidates. Keeping a fixed
  number of top matches per code instead of a similarity cutoff, dropping the
  chapter filter, and widening the co-occurrence list raises the ceiling from
  0.770 to 0.927 on ICD-10-CA and 0.920 to 0.976 on ICDA-8.
- The chapter filter throws away more right answers than wrong ones. Chapter
  match is now something the model weighs, not a hard rule.

## 5. The new setup

Two stages instead of cutoffs and rules: cast a wide net, then sort.

- **Stage 1, cast a wide net.** Keep the top 10 matches per code from each
  model, plus the top 50 codes it co-occurs with. No chapter filter.
- **Stage 2, sort them.** An xgboost model scores every candidate pair using 52
  pieces of information:
  - how similar the labels are, and what rank that puts the candidate at, for
    each model
  - the same thing in reverse, checking the match from the ICD-10 side too
  - whether the two codes are each other's best match, and how far apart the two
    directions disagree
  - how many other ICD-9 codes are competing for the same target
  - a combined ranking across all three models
  - how often the two codes appear together in patient records
  - plain text overlap: shared words, how many letters you'd have to change to
    turn one label into the other, shared starting words
  - whether the chapters line up, plus how often that particular pair of
    chapters actually matches in the training data
- **Output rule.** Keep a pair if its score is high enough on its own, or close
  enough to the best score for that code, or if it is the best score for that
  code. Both cutoffs are tuned on training data only, never on test data.

Results, tested on codes the model never saw during training:

| | ICD-10-CA | ICDA-8 |
|---|---|---|
| original ClinicalBERT version | 0.427 | 0.716 |
| SapBERT, same old rules | 0.530 | 0.821 |
| SapBERT, same old rules, tested properly | 0.546 | 0.824 |
| wide net + sorting | **0.668** | **0.840** |

## 6. Changes to how it's scored

- The original numbers were measured on the same data used to pick the settings,
  which makes them look better than they are. Everything now splits the codes
  into 5 groups, trains on 4 and tests on the 5th, and rotates. No code is ever
  in both training and testing.
- Every mapping now gets a confidence score. If you set the bar at 95%
  precision:

| | ICD-10-CA | ICDA-8 |
|---|---|---|
| auto-accept | 86.1% | 84.1% |
| human review, has candidates but unsure | 13.6% | 11.9% |
| human review, nothing plausible found | 0.3% | 4.0% |

- Added a check for whether the right answer is even reachable, separate from
  whether the output rule keeps it:

| | right answer ranked 1st | in top 3 | in top 5 |
|---|---|---|---|
| ICD-10-CA | 93.3% | 97.7% | 98.6% |
| ICDA-8 | 88.4% | 93.7% | 94.7% |

## 7. Which parts actually matter

Removed one group of features at a time and re-tested:

| removed | ICD-10-CA | change | ICDA-8 | change |
|---|---|---|---|---|
| nothing | 0.669 | | 0.854 | |
| checking the match in reverse | 0.606 | -0.063 | 0.843 | -0.011 |
| co-occurrence | 0.618 | -0.051 | 0.819 | -0.035 |
| chapter info | 0.618 | -0.051 | 0.851 | -0.003 |
| text overlap | 0.655 | -0.015 | 0.854 | -0.000 |
| ClinicalBERT | 0.657 | -0.012 | 0.857 | +0.003 |
| combined ranking | 0.663 | -0.006 | 0.838 | -0.016 |
| mpnet | 0.664 | -0.005 | 0.835 | -0.019 |

- Checking each match from both directions is the biggest single win on
  ICD-10-CA.
- ClinicalBERT adds nothing once SapBERT and mpnet are in. Could be dropped.
- Tried two other scoring methods for stage 2. Both are worse on ICD-10-CA. One
  is slightly better on ICDA-8 (0.862 vs 0.854).
- **More training data wouldn't help.** Tested this by training on smaller and
  smaller amounts. Performance flattens out at around 100-180 codes on both
  tracks, and we already have more than that.
- **It's not memorizing.** Retested with whole medical categories held out, so
  it has never seen anything similar at test time. Costs 0.004 on ICD-10-CA and
  0.000 on ICDA-8, which is noise.
- **It's not fragile.** Tried 24 different settings for how wide to cast the
  net. Results barely move, so nothing is delicately tuned.

## 8. Taking the code numbers out of the input

The model was being fed the code number and the label together, like
"250 diabetes mellitus". I tested feeding it just the label.

| | right answer ranked 1st | in top 10 | in top 25 | in top 50 |
|---|---|---|---|---|
| SapBERT, code + label, ICD-10-CA | 79.7% | 58.5% | 68.3% | 73.5% |
| SapBERT, label only, ICD-10-CA | **90.1%** | 64.8% | 73.8% | 81.5% |
| mpnet, code + label, ICD-10-CA | 80.9% | 65.5% | 77.0% | 85.1% |
| mpnet, label only, ICD-10-CA | **89.9%** | 67.7% | 77.8% | 85.8% |
| SapBERT, code + label, ICDA-8 | **84.1%** | 93.1% | 95.5% | 97.6% |
| SapBERT, label only, ICDA-8 | 81.8% | 90.6% | 92.5% | 94.0% |

- On ICD-10-CA, dropping the number puts the right answer first 10 points more
  often. The number was just noise, because ICD-9 and ICD-10 numbering have
  nothing to do with each other.
- On ICDA-8 it's slightly worse, because 69.8% of ICDA-8 answers have the same
  number as the ICD-9 code they come from. There the number is a real clue.
- Note these are retrieval numbers (SapBERT alone, before the sorting stage).
  They aren't the same as the 93.3% in section 6, which is the full pipeline.
- The label-only tables are generated but haven't been run through the full
  pipeline yet.

## 9. What's still going wrong

How complete the answers are, per code:

| | ICD-10-CA | ICDA-8 |
|---|---|---|
| found every correct answer | 39.5% | 81.6% |
| found some but not all | 46.2% | 4.4% |
| found none | 14.2% | 14.0% |
| codes with 1 answer, fully solved | 86.7% (n=143) | 88.2% (n=271) |
| codes with 2+ answers, fully solved | 6.0% (n=201) | 0.0% (n=22) |

- ICD-9 codes have 2.72 ICD-10-CA answers on average (up to 13, and 63.2% have
  more than one), but only 1.10 ICDA-8 answers (up to 6, only 7.9% have more
  than one). The original methods document counted the same thing: 127 codes
  with one ICD-10-CA answer vs 218 with several, against 278 vs 24 on ICDA-8.
- Codes with a single answer do equally well on both tracks. So the whole gap
  between tracks is just how many multi-answer codes each one has, not that one
  code system is harder than the other.
- Checked whether the extra answers for a multi-answer code sit next to each
  other in the code list, which would let me grab them as a block. They partly
  do. Haven't built it yet.

## 10. Housekeeping

- Repo made public, MIT license.
- Added a script that shows a real example of the similarity table with the
  correct answers marked, for explaining how this works
  (`24_show_similarity_matrix.R`).
- Deleted the R Markdown report, replaced by this file.

## Things I got wrong and fixed

Listing these because they changed numbers I'd already reported.

- The feature-removal script was calculating the correct answers *after* the
  filtering step instead of before, which quietly dropped the answers filtering
  had already missed and made the scores look better than they were. I'd claimed
  "+0.034 free accuracy" off this. Withdrawn. Fixed, and it now agrees with the
  main results to within 0.001.
- Reported coverage as 43% / 79% at first. That was against reachable pairs.
  Against the full crosswalk it's 37.5% / 76.1%.
- Reported human review as 57% at first. That was the share of pairs, not codes.
  It's 14% of codes.
- Said nine codes had no correct answer and were being silently skipped. Wrong,
  they're on the crosswalk's own exclusion list.

## What's left to do

- Run the label-only version through the whole pipeline.
- Build something specifically for multi-answer codes.
- Bigger validation set. Right now it's 937 ICD-10-CA pairs and 331 ICDA-8 pairs
  across 354 three-digit ICD-9 codes, not the full code list. The CIHI crosswalk
  table the original used would be closer to what the published LLM paper tested
  on.
- More co-occurrence data. It's the second most important feature and it's
  missing for a lot of pairs.
