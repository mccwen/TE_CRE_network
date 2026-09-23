# Example data and required inputs

## Biological question

Transposable elements (TEs) can influence gene regulation by overlapping or
contributing sequence to cis-regulatory elements. This workflow asks whether
co-accessible candidate regulatory elements near a target gene overlap
annotated TEs and, if so, which TE families or genomic TE instances are
represented.

## Example multiome dataset

The example is based on the 10x Genomics
[Fresh Embryonic E18 Mouse Brain 5k single-cell multiome dataset](https://www.10xgenomics.com/datasets/fresh-embryonic-e-18-mouse-brain-5-k-1-standard-1-0-0).
The paired gene-expression and chromatin-accessibility measurements make the
dataset suitable for demonstrating gene-peak correlation and peak
co-accessibility analysis.

The repository does not redistribute the complete 10x dataset. Download the
processed data from 10x Genomics and prepare a matrix directory containing:

| File | Purpose |
| --- | --- |
| `matrix.mtx.gz` | Sparse gene-expression and peak-accessibility matrix |
| `barcodes.tsv.gz` | Cell barcodes |
| `features.tsv.gz` | Feature identifiers, names, types, and chromosome annotation |
| `peaks.bed` | Genomic coordinates of ATAC-seq peaks |

The Cicero script requires a fourth column in `features.tsv.gz` containing the
chromosome for each feature. Standard 10x feature files generally contain only
three columns, so this chromosome annotation must be added during input
preparation.

## Stage 1: gene-peak correlation and Cicero co-accessibility

`Rscript/cicero_monocle3.R` uses the matrix directory to:

1. Calculate correlations between the target gene and peaks on the requested
   chromosome.
2. Construct a Monocle 3 representation of the peak-accessibility matrix.
3. Calculate Cicero co-accessibility scores.
4. Save all connections and a threshold-filtered connection table.

The example filtered table stored in this repository is:

```text
example_data/chr2_coaccess_score_gt0.2.csv
```

## Stage 2: TE intersection and quantification

The TE-quantification tutorial uses:

1. A Cicero CSV containing `Peak1`, `Peak2`, and `coaccess`.
2. A target chromosome and genomic interval.
3. A Dfam mouse TE annotation in `nrph.hits` or `nrph.hits.gz` format.
4. BEDTools for genomic interval intersection.

The mm10 Dfam annotation is available from:

<https://www.dfam.org/releases/Dfam_3.8/annotations/mm10/>

The large Dfam annotation is intentionally excluded from Git. The tutorial
accepts a local annotation in the repository root or a path supplied through
the `TE_ANNOTATION` environment variable.

The tutorial writes TE-family and unique-TE matrices under:

```text
results/TE_quant/
```

## Genome-build consistency

The 10x peak coordinates, feature annotations, chromosome-sizes file, target
interval, and Dfam TE annotations must all use the same genome build. The
documented TE-quantification example uses mm10. Use Cicero's built-in
`mouse.mm9` genome only when every other input also uses mm9 coordinates.
