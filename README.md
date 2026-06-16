# KNNTrace-Listeria-dataset

This repository contains the core dataset files and analysis scripts supporting the wgMLST-based KNNTrace framework for *Listeria monocytogenes* source attribution.

## Repository contents

```text
KNNTrace-Listeria-dataset/
├── metadata/
├── method_tables/
├── allele_profiles/
├── feature_ranking/
├── scripts/
└── docs/
```

## Directory description

### metadata/

This directory contains isolate-level metadata used in the study, including the full NCBI metadata pool, the food/environmental training set, the external human isolate set, and the seed-genome metadata used for wgMLST scheme construction.

### method_tables/

This directory contains method-related lookup tables, including the source-category mapping table used to harmonize NCBI isolation-source metadata into the predefined source categories.

### allele_profiles/

This directory contains the cgMLST, full wgMLST, and Top800 wgMLST allele-profile matrices for the training isolates and external human isolates.

### feature_ranking/

This directory contains wgMLST locus-ranking results, the final Top800 locus list, and functional annotation information for the retained wgMLST loci.

### scripts/

This directory contains the core analysis scripts used for wgMLST locus filtering, locus ranking, model benchmarking, source-specific threshold calculation, and post-filtering of external human predictions.

### docs/

This directory contains documentation files such as the file manifest, data dictionary, and repository-structure description.

## Notes

Raw NCBI genome FASTA files, the complete chewBBACA scheme directory, intermediate cache files, model objects, and full model-output directories are not included in this repository. The released metadata, allele-profile matrices, feature-ranking files, and analysis scripts provide the core data and code required to reproduce the main analyses reported in the manuscript.

## License

This dataset is released under the CC0-1.0 license.
