# Repository structure

This document summarizes the organization of the `KNNTrace-Listeria-dataset` repository.

## Top-level directories

```text
KNNTrace-Listeria-dataset/
├── metadata/
├── method_tables/
├── allele_profiles/
├── feature_ranking/
├── scripts/
└── docs/
```

## Directory descriptions

### `metadata/`

Contains isolate-level metadata used in the study.

Expected files:

```text
all_metadata.jsonl.bz2
train_metadata.tsv.bz2
human_metadata.tsv.bz2
seed_metadata.tsv.bz2
```

These files describe the initial NCBI metadata pool, the nine-source food/environmental training set, the external human isolate set, and the seed genomes used for study-specific wgMLST scheme construction.

### `method_tables/`

Contains method-related lookup tables.

Expected file:

```text
source_category_mapping.tsv
```

This table documents the source-category mapping rules used to harmonize NCBI isolation-source text into the predefined categories used in the study.

### `allele_profiles/`

Contains compressed cgMLST, full wgMLST, and Top800 wgMLST allele-profile matrices.

Expected files:

```text
cgMLST_1748_train_matrix.tsv.bz2
cgMLST_1748_human_matrix.tsv.bz2
wgMLST_3586_train_matrix.tsv.bz2
wgMLST_3586_human_matrix.tsv.bz2
Top800_wgMLST_train_matrix.tsv.bz2
Top800_wgMLST_human_matrix.tsv.bz2
```

The training matrices correspond to the nine food/environmental source categories. The human matrices correspond to the external human clinical isolate set.

### `feature_ranking/`

Contains wgMLST locus-ranking and annotation files.

Expected files:

```text
wgMLST_locus_ranking_scores.tsv
Top800_loci.tsv
wgMLST_loci_annotation.tsv
```

These files document the ranked wgMLST loci, the final Top800 feature set, and functional annotation information for the retained loci.

### `scripts/`

Contains the core analysis scripts.

Expected files:

```text
wg_locus_filter.R
knn_locus_rank.R
cg_thresh.R
wg_thresh.R
post_filter_cg.R
post_filter_wg800.R
glmnet_cg.R
glmnet_wg.R
knn_cg.R
knn_wg800.R
rf_cg.R
rf_wg.R
xgb_cg.R
xgb_wg.R
```

These scripts cover wgMLST locus filtering, locus ranking, model benchmarking, source-specific threshold calculation, and post-filtering of external human predictions.

### `docs/`

Contains documentation files for the repository.

Expected files:

```text
file_manifest.tsv
data_dictionary.tsv
repository_structure.md
```

## Files not included

Raw NCBI genome FASTA files, the complete chewBBACA scheme directory, intermediate cache files, model objects, and full model-output directories are not included. The released metadata, allele-profile matrices, feature-ranking files, and analysis scripts provide the core data and code needed to reproduce the main analyses reported in the manuscript.

## Suggested use

1. Inspect `metadata/` and `method_tables/` to understand sample selection and source-category harmonization.
2. Use `allele_profiles/` as model input matrices.
3. Use `feature_ranking/` to identify the retained wgMLST loci and the final Top800 feature set.
4. Use `scripts/` to reproduce locus filtering, feature ranking, model benchmarking, threshold calculation, and external human prediction filtering.
