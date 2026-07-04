## Leiden Cluster Merging

## Introduction

Prototype tools for quantitatively testing whether selected Leiden clusters are reasonable merge candidates in Scanpy workflows.
After creating leiden clusters and generating UMAP with appropriate parameters such as resolution and min dust size, run this tool to tell you how reasonable it would be to merge any leiden clusters you think may be the same cells.


## How it works

This tool is meant to be used after you already have a normal Scanpy analysis.

For example, you may already have:

- a processed `AnnData` object called `adata`
- Leiden cluster labels in `adata.obs["leiden"]`
- PCA results in `adata.obsm["X_pca"]`
- a neighbor graph from `sc.pp.neighbors`
- a UMAP plot that helps you inspect the clusters

The tool does not create clusters for you. Instead, it helps you test a question like:

> Are clusters 3 and 7 similar enough that merging them is reasonable?

You choose the clusters you want to test, and the function calculates several statistical measurements. Each measurement looks at the proposed merge from a different angle.

### 1. It selects the cells and genes to test

First, the function keeps only the cells from the clusters you selected.

For example, if you test:

```python
clusters=["3", "7"]
```

the function only uses cells from clusters 3 and 7.

Then it chooses genes for testing. If your `adata` has highly variable genes, it uses those. If not, it uses the most variable genes. This keeps the test focused on genes that carry useful biological information.

### 2. It compares average expression profiles

For each selected cluster, the function calculates the average expression of each selected gene.

You can imagine this as making one summary expression profile for each cluster.

Then it compares those profiles using correlation:

- Pearson correlation
- Spearman correlation

High correlation means the clusters have similar gene expression patterns. This supports merging.

Low correlation means the clusters look different. This argues against merging.

For more than two clusters, the function compares every pair. For example, if you test clusters 3, 7, and 9, it compares:

- 3 vs 7
- 3 vs 9
- 7 vs 9

The function reports both the average pairwise similarity and the lowest pairwise similarity. The lowest similarity is important because one very different cluster can make a proposed merge questionable.

### 3. It runs a permutation test for expression similarity

The function also asks:

> Are these clusters more similar than expected by random chance?

To test this, it randomly shuffles the cluster labels many times. Each shuffle creates a random comparison. The function then checks whether the real cluster similarity is stronger than the shuffled results.

This gives a permutation p-value.

A small p-value means the selected clusters are unusually similar compared with random label shuffling.

### 4. It tests for differentially expressed genes

The function checks whether genes are significantly different between the clusters.

If you test exactly two clusters, it uses the Mann-Whitney U test for each gene.

If you test three or more clusters, it uses the Kruskal-Wallis test for each gene.

These tests ask:

> Does this gene show different expression across the selected clusters?

Because many genes are tested, the function adjusts the p-values using Benjamini-Hochberg correction. This helps control false positives.

The function reports:

- gene-level p-values
- adjusted p-values
- effect sizes
- log2 fold changes
- how many genes are significantly different
- how many genes are significantly different with a large enough effect size

This matters because a tiny difference can become statistically significant when there are many cells. The tool therefore looks at both p-values and effect sizes.

If many genes are strongly different, merging is less supported.

If only a few genes are different, merging may be more reasonable.

### 5. It checks marker agreement

The function checks whether the selected clusters rank genes in a similar way.

For two clusters, it uses Spearman correlation.

For three or more clusters, it uses a rank agreement measurement similar to Kendall's W.

This helps answer:

> Do these clusters have similar marker-gene patterns?

High marker agreement supports merging.

Low marker agreement suggests the clusters may represent different cell states or cell types.

### 6. It compares cluster positions in PCA space

If `adata.obsm["X_pca"]` is available, the function calculates the center point of each selected cluster in PCA space.

Then it measures the distance between those cluster centers.

Small distance means the clusters are close in PCA space. This supports merging.

Large distance means the clusters are farther apart. This argues against merging.

### 7. It checks neighborhood mixing

If your `adata` has a neighbor graph in `adata.obsp["connectivities"]`, the function checks whether cells from the selected clusters are mixed together in the graph.

This is useful because Leiden clustering is based on the neighbor graph.

The function measures how often cells connect to cells from another selected cluster.

High neighborhood mixing means the clusters are strongly connected to each other. This supports merging.

Low neighborhood mixing means the clusters are separated in the graph. This argues against merging.

The function also runs a permutation test for this measurement by shuffling cluster labels.

### 8. It checks sample or batch composition

If you provide `sample_key`, the function checks whether the selected clusters have different sample or batch composition.

For example:

```python
sample_key="sample"
```

The function creates a table comparing clusters and samples. Then it runs a chi-square test and calculates Cramer's V.

This helps answer:

> Are these clusters different because of biology, or because of sample/batch effects?

A strong sample or batch association does not automatically mean clusters should merge. It is a warning that the user should inspect the data carefully.

### 9. It creates component scores

The function converts the main results into simple scores from 0 to 1.

The component scores are:

- expression similarity
- marker agreement
- inverse differential-expression evidence
- PCA/embedding similarity
- neighborhood mixing
- sample consistency

A score closer to 1 means stronger support for merging.

A score closer to 0 means weaker support for merging.

### 10. It calculates an overall merge score

Finally, the function combines the component scores into one overall score.

By default, expression similarity and differential-expression evidence are weighted most strongly.

The default weights are:

- expression similarity: 25%
- marker agreement: 20%
- differential-expression evidence: 25%
- PCA/embedding similarity: 10%
- neighborhood mixing: 15%
- sample consistency: 5%

The overall score is not meant to replace biological judgment. It is a quantitative guide.

In simple terms:

- high score: the selected clusters look more reasonable to merge
- low score: the selected clusters look less reasonable to merge

You should still inspect marker genes, UMAP, known cell-type markers, and sample information before making the final decision.

### Main outputs

The function returns a dictionary.

The most useful outputs are:

```python
result["summary"]
```

A one-row overview of the proposed merge.

```python
result["tests"]
```

The main statistical tests and p-values.

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
=======
