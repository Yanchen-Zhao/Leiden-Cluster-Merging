# Leiden Cluster Merging

Prototype Python code for testing whether selected Leiden clusters are reasonable candidates to merge in a Scanpy single-cell workflow.

This project currently provides one function:

```python
test_merge()
```

The function is designed to run **after** you have already processed your `AnnData` object, calculated neighbors, generated UMAP, and created Leiden clusters.

It does not create clusters for you. Instead, it helps answer:

> Given two or more existing Leiden clusters, how reasonable is it to merge them?

## Example

For a UMAP like this:

<img height="360" alt="UMAP example" src="https://github.com/user-attachments/assets/663f8c6e-2fbf-455e-af96-281bac94566c" />

The UMAP above comes from the Scanpy clustering tutorial:

https://scanpy.scverse.org/en/stable/tutorials/basics/clustering.html

To test whether clusters 0 and 7 are reasonable to merge:

```python
from leiden_merge_test import test_merge

score = test_merge(
    adata,
    cluster_key="leiden_res_0.50",
    clusters=["0", "7"],
)

score
```

By default, the function returns one number between 0 and 1.

- `0` means the selected clusters look less reasonable to merge
- `1` means the selected clusters look more reasonable to merge

For example, a score of `0.54659` means the evidence is mixed. It is not a final biological decision.

## Basic Use

The default setting is:

```python
simple=True
```

This means the function returns only the final merge-support score:

```python
score = test_merge(
    adata,
    cluster_key="leiden",
    clusters=["0", "2"],
)
```

This mode is faster because it skips p-value-only permutation calculations and detailed output tables.

## Detailed Use

If you want p-values and detailed tables, use:

```python
result = test_merge(
    adata,
    cluster_key="leiden",
    clusters=["0", "2"],
    simple=False,
)
```

Then inspect:

```python
result["p_values"]
```

This table gives p-values for the main statistical tests, such as:

- Pearson correlation test
- Spearman correlation test
- expression similarity permutation test
- marker agreement test
- differential expression combined test
- neighborhood mixing permutation test
- sample composition chi-square test, if `sample_key` is provided

## Inputs

Important arguments:

```python
test_merge(
    adata,
    cluster_key,
    clusters,
    sample_key=None,
    use_rep="X_pca",
    layer=None,
    genes=None,
    max_genes=2000,
    n_permutations=100,
    max_permutation_cells=5000,
    max_de_cells_per_cluster=5000,
    alpha=0.05,
    logfc_threshold=0.25,
    random_state=0,
    weights=None,
    simple=True,
)
```

### Required Inputs

`adata`

Your processed Scanpy `AnnData` object.

`cluster_key`

The column in `adata.obs` containing Leiden cluster labels.

Example:

```python
cluster_key="leiden"
```

or:

```python
cluster_key="leiden_res_0.50"
```

`clusters`

The cluster labels you want to test as one possible merge group.

Example:

```python
clusters=["0", "7"]
```

or:

```python
clusters=["0", "7", "12"]
```

The code supports two clusters and more than two clusters.

### Optional Inputs

`sample_key`

Use this if `adata.obs` has sample or batch labels.

Example:

```python
sample_key="sample"
```

This lets the function test whether the selected clusters have different sample or batch composition.

`use_rep`

The embedding used to compare cluster positions. The default is:

```python
use_rep="X_pca"
```

This uses `adata.obsm["X_pca"]` if it exists.

`layer`

Expression layer to use. If `None`, the function uses `adata.X`.

`genes`

Optional list of genes to test. If not provided, the function uses highly variable genes when available. If highly variable genes are not available, it uses the most variable genes.

`max_genes`

Maximum number of genes to test when `genes=None`.

`n_permutations`

Number of permutations used in detailed mode for permutation p-values.

This is ignored when `simple=True`.

`max_permutation_cells`

Maximum number of selected cells used for permutation tests.

This keeps detailed mode from becoming too slow on large datasets.

`max_de_cells_per_cluster`

Maximum number of cells per cluster used for gene-level differential expression tests.

`simple`

If `True`, return only the final merge-support score.

If `False`, return detailed result tables.

## How It Works

The function compares the selected clusters from several angles.

### 1. Select Cells And Genes

First, the function keeps only cells from the clusters you selected.

For example:

```python
clusters=["3", "7"]
```

means only cells from clusters 3 and 7 are used.

Then the function selects genes.

It uses highly variable genes if they are available in:

```python
adata.var["highly_variable"]
```

If not, it chooses the most variable genes.

### 2. Compare Average Expression Profiles

For each selected cluster, the function calculates the average expression of each selected gene.

This creates one expression profile per cluster.

Then the function compares those profiles using:

- Pearson correlation
- Spearman correlation

High correlation supports merging.

Low correlation argues against merging.

For more than two clusters, every pair is compared. For example:

```python
clusters=["3", "7", "9"]
```

compares:

- 3 vs 7
- 3 vs 9
- 7 vs 9

The function uses the **lowest pairwise Pearson correlation** as one part of the final score. This is important because one dissimilar cluster can make a proposed merge questionable.

### 3. Test Differential Expression

The function tests whether genes are different between the selected clusters.

For exactly two clusters, it uses:

```text
Mann-Whitney U test
```

For three or more clusters, it uses:

```text
Kruskal-Wallis test
```

Because many genes are tested, p-values are adjusted using Benjamini-Hochberg correction.

The function counts:

- number of significant genes
- number of significant genes with large enough log2 fold change

Many strongly different genes argue against merging.

Few strongly different genes support merging.

### 4. Check Marker Agreement

The function checks whether selected clusters rank genes similarly.

For two clusters, it uses Spearman correlation.

For three or more clusters, it uses a Kendall's W-like rank agreement statistic.

High marker agreement supports merging.

Low marker agreement argues against merging.

### 5. Compare PCA Positions

If `adata.obsm["X_pca"]` exists, the function calculates a center point for each selected cluster in PCA space.

Then it measures the distance between cluster centers.

Small distance supports merging.

Large distance argues against merging.

### 6. Check Neighborhood Mixing

If `adata.obsp["connectivities"]` exists, the function checks whether cells from the selected clusters are mixed in the Scanpy neighbor graph.

This matters because Leiden clustering is based on the neighbor graph.

High neighborhood mixing supports merging.

Low neighborhood mixing argues against merging.

In detailed mode, the function also runs a permutation test for neighborhood mixing.

### 7. Check Sample Or Batch Composition

If `sample_key` is provided, the function checks whether the selected clusters have different sample or batch composition.

It uses:

- chi-square test
- Cramer's V

Strong sample or batch association does not automatically mean clusters should or should not merge. It is a warning to inspect the biology carefully.

### 8. Calculate Component Scores

The function converts the evidence into component scores from 0 to 1:

- expression similarity
- marker agreement
- inverse differential-expression evidence
- PCA or embedding similarity
- neighborhood mixing
- sample consistency

A component score closer to 1 supports merging.

A component score closer to 0 argues against merging.

### 9. Calculate Overall Score

The final merge-support score is a weighted average of the component scores.

Default weights:

- expression similarity: 25%
- marker agreement: 20%
- differential-expression evidence: 25%
- PCA/embedding similarity: 10%
- neighborhood mixing: 15%
- sample consistency: 5%

The score is meant to guide inspection. It should not replace biological judgment.

## Detailed Outputs

When `simple=False`, the function returns a dictionary.

Useful outputs:

```python
result["summary"]
```

One-row overview of the proposed merge.

```python
result["p_values"]
```

Easy-to-read table of the main statistical test p-values.

In this table:

- `value` is the statistic or measured quantity
- `p_value` is the p-value for that statistic, when available

For example:

- Pearson correlation test: `value` is Pearson `r`
- Spearman correlation test: `value` is Spearman `r`
- neighborhood mixing permutation test: `value` is the fraction of graph connections crossing between selected clusters
- differential expression combined test: `value` is currently the combined p-value

```python
result["tests"]
```

Detailed test table with short interpretations.

```python
result["pairwise"]
```

Pairwise comparisons between selected clusters.

```python
result["differential_expression"]
```

Gene-level differential-expression results.

```python
result["component_scores"]
```

The 0-to-1 score for each evidence type.

```python
result["overall_score"]
```

The final combined merge-support score.

```python
result["cluster_sizes"]
```

Number of selected cells in each cluster.

```python
result["cluster_means"]
```

Average expression profile for each selected cluster.

```python
result["sample_table"]
```

Cluster-by-sample table, only available if `sample_key` is provided.

```python
result["selected_genes"]
```

Genes used in the test.

## Notes

Very small p-values may appear as `0.0`. This usually means the p-value is smaller than Python can represent, not that the true probability is literally zero.

Permutation p-values are limited by `n_permutations`. For example, with `n_permutations=100`, the smallest possible permutation p-value is about `1 / 101`.

Single-cell datasets often have many cells, so tiny differences can become statistically significant. For this reason, the overall score uses both statistical evidence and effect-size-like summaries.
