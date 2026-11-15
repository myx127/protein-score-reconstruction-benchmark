# Canonical input contracts

Configuration paths resolve relative to the JSON file containing them. Defaults for weights, masks, direction units and seeds resolve from the checkout. The examples in `config/` point outside the checkout; adjust them after copying to an external local configuration. Use `data_kind: "governed"` for any real data. The `synthetic` label is reserved for independently generated examples.

## Recovery matrices

One UTF-8 CSV or `.csv.gz` per batch. The first two columns must be `participant_id,record_id`, followed by protein columns. IDs are strings, including leading zeros. The pair of IDs must be unique and complete. Protein column names must match the frozen identifiers **exactly**, including `ACC_TOKEN::` prefixes when present. Column order may differ; matching is by name. All required proteins and finite, already preprocessed values are required; the runner neither imputes recovery matrices nor repeats their normalization.

```text
participant_id,record_id,ACC_TOKEN::EXAMPLE1,ACC_TOKEN::EXAMPLE2,...
SYN_EXAMPLE_001,SYN_RECORD_001,0.12,-0.34,...
```

The example row is schematic, not a real record. The paper matrices have 722/833 rows and 1,222 proteins. Source and target participant sets must not overlap. A generic run needs at least 63 rows per batch so that each outer-fold training set has at least 50 rows; at least 125 rows are needed for five inner folds throughout the five-fold within-target diagnostic. Frozen paper inputs already satisfy these constraints. No target outcome or clinical endpoint is needed.

`metadata/expected_governed_inputs.json` records exact paper matrix hashes, without data. `enforce_paper_fingerprints` defaults to true in governed recovery/timing runs. For an explicitly identified new-data application it may be false, with compatible weights/masks/units/seeds supplied; that run is not an exact paper rerun.

## Shared recovery metadata

| File/config key | Required columns and role |
| --- | --- |
| `weights` | `module,protein,weight`; unique module/protein rows, finite weights and positive absolute weight norm. Extra `source_system` is descriptive. |
| `masks` | `mask_id,protein`; observed proteins for each structural mask, no repeated protein within a mask. |
| `units` | `direction_unit_id,mask_unit_id,module,platform_label,direction,train_batch,test_batch`; one unique recovery setting per row. |
| `seeds` | `stream,replicate,seed,role`; `fold::<direction_unit_id>` at replicate 0, and `<batch>_target` at each bootstrap replicate with role `primary_target_bootstrap`. |

Every mask must include at least two proteins and fewer than the full module. The frozen benchmark has 36 module/mask combinations and both directions, for 72 settings. Public table counts and hashes are checked by tests. Default timing uses the separate 24-unit registry and original randomized order in `metadata/paper_timing_*.csv`.

## Emory source-definition matrix

One CSV, first column a unique source sample identifier (e.g. `sample_id`), then finite protein columns matching **raw `ProteinID`** values from `metadata/emory_module_weights.csv`. This matrix is the approved source representation for CRNP fitting. It is distinct from the batch-standardized recovery matrices; do not infer its scale or create it by concatenation without the original approved preprocessing.

Emory weights use `ModuleColor,ProteinID,kME_ownModule`. Mapping uses `input_id,feature_id,protein_id,conservative_eligible,mapping_route`. Protein/feature IDs are exact strings. Only conservative eligible mappings are used. The mapping is tied to the documented assay annotation version; a new assay version requires a separately audited mapping.

## Six paired assay-fluid matrices

Each branch is a CSV with first column `sample_id`, followed by quantitative assay feature columns. The six `input_id` values are `soma_csf,soma_plasma,olink_csf,olink_plasma,tmt_csf,tmt_plasma`. Platform strings must be `SomaScan`, `Olink`, `TMT-MS`; biofluid strings must be `CSF`, `plasma`. Samples must use the same string identifiers across branches. The paper analysis uses the lexically sorted intersection of **35 paired participants**.

Specify `input_scale` separately for each branch:

| Scale | Required input and transformation performed by the runner |
| --- | --- |
| `raw`, SomaScan | Positive RFU; nonpositive/nonfinite becomes missing; log2 is applied once. |
| `raw`, TMT-MS | **Ratios already formed relative to the 127N reference channel within plex**; nonpositive/nonfinite becomes missing; log2 is applied once. Raw reporter intensities are not accepted as equivalent ratios. |
| `raw`, Olink | NPX as supplied; no additional log transform. |
| `transformed_features` | Feature-level values already on the above log2/NPX scales; no additional log transform. Not a protein-aggregated or previously restandardized matrix. |

The runner retains features with at least 80% finite values (28/35), imputes remaining values with feature medians, standardizes by sample SD, aggregates multiple eligible features per protein by participant median, restandardizes proteins, and excludes zero-variance features/proteins. It then fits source CRNP and the source-calibrated direct comparator and standardizes each branch's scores, following the frozen Emory evaluation. This processing intentionally uses the paired evaluation sample distribution and is not a prospective one-patient deployment procedure.

`strict_paper_design: true` requires 35 paired participants and 52 evaluable concordance cells from 54 planned cells; 27 CSF and 25 plasma cells were evaluable in the paper. Set false only for an explicitly identified different/synthetic design. Bootstrap requires a positive draw count, 10,000 by default. Incomplete/no-valid-cell results must not be silently interpreted as negative evidence.

## Outputs

`RUN_STATUS.json` records data kind, environment and counts. Governed work folders contain participant predictions and R scores and must stay controlled. Recovery's `aggregate/` contains cell metrics, equal-direction point estimates, bootstrap aggregates and intervals. Emory produces concordance and bootstrap aggregates; timing produces interleaved order, warm-ups, all repeat records, workload medians and summaries. Only reviewed aggregate outputs may be selected for a later public release; the runner does not publish them automatically.
