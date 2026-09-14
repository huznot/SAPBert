# Asked for on 9 September

ICD-9-CM to ICD-10-CA unless a sheet says otherwise. 354 ICD-9-CM codes against
2038 ICD-10-CA codes is 721,452 possible pairs, 937 of them marked correct in
the validation data.

| | |
|---|---|
| [1_threshold_frequencies.xlsx](1_threshold_frequencies.xlsx) | pairs above each cutoff, and TP/FP/FN as counts for all 1008 rules |
| [2_icda8_code_in_label.xlsx](2_icda8_code_in_label.xlsx) | the code number in the label, ICDA-8 against ICD-10-CA |
| [3_similarity_gap_by_model.xlsx](3_similarity_gap_by_model.xlsx) | correct pairs against everything else, per model |
| [4_counts_by_ccs_category.xlsx](4_counts_by_ccs_category.xlsx) | all 130 CCS categories, counts |

Earlier: [the 52 codes with no ICDA-8 match](../review/52_codes_no_icda8.xlsx),
[the nine codes](../review/9_codes_review_UPDATED.xlsx),
[249 and 239](../review/e13_d48_review.xlsx).

Rebuild with `scripts/53_`..`56_`, in that order.
