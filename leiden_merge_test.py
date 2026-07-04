def test_merge(
    adata,
    cluster_key,
    clusters,
    sample_key=None,
    use_rep="X_pca",
    layer=None,
    genes=None,
    max_genes=2000,
    n_permutations=1000,
    alpha=0.05,
    logfc_threshold=0.25,
    random_state=0,
    weights=None,
):
    """
    Quantitatively test whether selected Leiden clusters are reasonable merge candidates.

    Parameters
    ----------
    adata
        AnnData object after a standard Scanpy workflow.
    cluster_key
        Column in adata.obs containing Leiden cluster labels, for example "leiden".
    clusters
        Cluster labels to test as one proposed merge group. Must contain at least 2 labels.
    sample_key
        Optional column in adata.obs containing sample/batch labels.
    use_rep
        Optional embedding in adata.obsm, usually "X_pca". Set to None to skip.
    layer
        Optional adata.layers key to use for expression. If None, uses adata.X.
    genes
        Optional list of genes to test. If None, uses highly variable genes when available,
        otherwise the most variable genes up to max_genes.
    max_genes
        Maximum number of genes to use when genes is None.
    n_permutations
        Number of permutations for expression-profile and neighborhood-mixing tests.
    alpha
        Adjusted p-value threshold used to count significant genes.
    logfc_threshold
        Absolute log fold-change threshold used to count large-effect genes.
    random_state
        Random seed.
    weights
        Optional dict overriding component weights for the overall score.

    Returns
    -------
    dict
        Dictionary containing summary, tests, pairwise results, component scores,
        overall_score, and selected genes.
    """
    import itertools

    import numpy as np
    import pandas as pd
    from scipy import sparse, stats

    rng = np.random.default_rng(random_state)
    
#Benjamini-Hochberg correction, creates adjusted p-values as each statistical test would generate a p-value.
    def _bh_adjust(p_values):
        p_values = np.asarray(p_values, dtype=float)
        adjusted = np.full(p_values.shape, np.nan, dtype=float)
        valid = np.isfinite(p_values)
        if not valid.any():
            return adjusted

        p = p_values[valid]
        order = np.argsort(p)
        ranked = p[order]
        n = len(ranked)
        q = ranked * n / np.arange(1, n + 1)
        q = np.minimum.accumulate(q[::-1])[::-1]
        q = np.clip(q, 0, 1)

        valid_indices = np.flatnonzero(valid)
        adjusted[valid_indices[order]] = q
        return adjusted

    def _to_dense(x):
        if sparse.issparse(x):
            return x.toarray()
        return np.asarray(x)

    def _safe_mean(x, axis=0):
        if sparse.issparse(x):
            return np.asarray(x.mean(axis=axis)).ravel()
        return np.asarray(x).mean(axis=axis)

    def _safe_var(x, axis=0):
        dense = _to_dense(x)
        return np.var(dense, axis=axis)

    def _safe_pearson(x, y):
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        if x.size < 3 or np.std(x) == 0 or np.std(y) == 0:
            return np.nan, np.nan
        return stats.pearsonr(x, y)

    def _safe_spearman(x, y):
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        if x.size < 3 or np.std(x) == 0 or np.std(y) == 0:
            return np.nan, np.nan
        result = stats.spearmanr(x, y)
        return result.statistic, result.pvalue

    def _combine_pvalues(p_values):
        p_values = np.asarray(p_values, dtype=float)
        p_values = p_values[np.isfinite(p_values)]
        if len(p_values) == 0:
            return np.nan
        p_values = np.clip(p_values, np.nextafter(0, 1), 1.0)
        return stats.combine_pvalues(p_values, method="fisher").pvalue

    def _similarity_score(value, low=0.70, high=0.97):
        if not np.isfinite(value):
            return np.nan
        return float(np.clip((value - low) / (high - low), 0, 1))

    def _inverse_count_score(count, low=0, high=50):
        if not np.isfinite(count):
            return np.nan
        return float(1 - np.clip((count - low) / (high - low), 0, 1))

    def _distance_score(value, low=0.0, high=5.0):
        if not np.isfinite(value):
            return np.nan
        return float(1 - np.clip((value - low) / (high - low), 0, 1))

    def _cramers_v(table):
        chi2, _, _, _ = stats.chi2_contingency(table)
        n = table.to_numpy().sum()
        if n == 0:
            return np.nan
        r, k = table.shape
        denom = n * min(k - 1, r - 1)
        if denom == 0:
            return np.nan
        return float(np.sqrt(chi2 / denom))

    if cluster_key not in adata.obs:
        raise KeyError(f"{cluster_key!r} is not in adata.obs.")

    clusters = [str(cluster) for cluster in clusters]
    if len(clusters) < 2:
        raise ValueError("clusters must contain at least two cluster labels.")

    labels = adata.obs[cluster_key].astype(str)
    available = set(labels.unique())
    missing = [cluster for cluster in clusters if cluster not in available]
    if missing:
        raise ValueError(f"These clusters are not present in {cluster_key!r}: {missing}")

    selected_cells = labels.isin(clusters).to_numpy()
    selected_labels = labels.loc[selected_cells].to_numpy()
    cluster_sizes = pd.Series(selected_labels).value_counts().reindex(clusters).fillna(0).astype(int)

    if layer is None:
        X_all = adata.X
    else:
        if layer not in adata.layers:
            raise KeyError(f"{layer!r} is not in adata.layers.")
        X_all = adata.layers[layer]

    if genes is None:
        if "highly_variable" in adata.var:
            gene_mask = adata.var["highly_variable"].to_numpy(dtype=bool)
            gene_indices = np.flatnonzero(gene_mask)
            if len(gene_indices) > max_genes:
                gene_indices = gene_indices[:max_genes]
        else:
            variances = _safe_var(X_all, axis=0)
            gene_indices = np.argsort(variances)[::-1][:max_genes]
    else:
        gene_to_index = {gene: i for i, gene in enumerate(adata.var_names)}
        missing_genes = [gene for gene in genes if gene not in gene_to_index]
        if missing_genes:
            raise ValueError(f"These genes are not present in adata.var_names: {missing_genes[:10]}")
        gene_indices = np.array([gene_to_index[gene] for gene in genes], dtype=int)

    if len(gene_indices) == 0:
        raise ValueError("No genes were selected for testing.")

    selected_genes = adata.var_names[gene_indices].tolist()
    X = X_all[selected_cells][:, gene_indices]
    X_dense = _to_dense(X)

    cluster_means = {}
    for cluster in clusters:
        cluster_means[cluster] = _safe_mean(X_dense[selected_labels == cluster], axis=0)
    cluster_means_df = pd.DataFrame(cluster_means, index=selected_genes).T

    pairwise_rows = []
    for cluster_a, cluster_b in itertools.combinations(clusters, 2):
        mean_a = cluster_means[cluster_a]
        mean_b = cluster_means[cluster_b]
        pearson_r, pearson_p = _safe_pearson(mean_a, mean_b)
        spearman_r, spearman_p = _safe_spearman(mean_a, mean_b)

        if use_rep is not None and use_rep in adata.obsm:
            rep = np.asarray(adata.obsm[use_rep])[selected_cells]
            centroid_a = rep[selected_labels == cluster_a].mean(axis=0)
            centroid_b = rep[selected_labels == cluster_b].mean(axis=0)
            embedding_distance = float(np.linalg.norm(centroid_a - centroid_b))
        else:
            embedding_distance = np.nan

        pairwise_rows.append(
            {
                "cluster_a": cluster_a,
                "cluster_b": cluster_b,
                "n_a": int(cluster_sizes[cluster_a]),
                "n_b": int(cluster_sizes[cluster_b]),
                "pearson_r": pearson_r,
                "pearson_p_value": pearson_p,
                "spearman_r": spearman_r,
                "spearman_p_value": spearman_p,
                "embedding_centroid_distance": embedding_distance,
            }
        )

    pairwise = pd.DataFrame(pairwise_rows)
    mean_pairwise_pearson = float(pairwise["pearson_r"].mean(skipna=True))
    min_pairwise_pearson = float(pairwise["pearson_r"].min(skipna=True))
    mean_pairwise_spearman = float(pairwise["spearman_r"].mean(skipna=True))
    min_pairwise_spearman = float(pairwise["spearman_r"].min(skipna=True))
    mean_embedding_distance = float(pairwise["embedding_centroid_distance"].mean(skipna=True))

    expression_perm_p = np.nan
    if n_permutations > 0 and len(clusters) > 1:
        observed = mean_pairwise_pearson
        if np.isfinite(observed):
            permuted_stats = []
            for _ in range(n_permutations):
                permuted_labels = rng.permutation(selected_labels)
                perm_means = []
                valid_perm = True
                for cluster in clusters:
                    mask = permuted_labels == cluster
                    if mask.sum() == 0:
                        valid_perm = False
                        break
                    perm_means.append(X_dense[mask].mean(axis=0))
                if not valid_perm:
                    continue
                corrs = []
                for i, j in itertools.combinations(range(len(perm_means)), 2):
                    r, _ = _safe_pearson(perm_means[i], perm_means[j])
                    if np.isfinite(r):
                        corrs.append(r)
                if corrs:
                    permuted_stats.append(np.mean(corrs))
            if permuted_stats:
                permuted_stats = np.asarray(permuted_stats)
                expression_perm_p = float((np.sum(permuted_stats >= observed) + 1) / (len(permuted_stats) + 1))

    de_rows = []
    if len(clusters) == 2:
        cluster_a, cluster_b = clusters
        Xa = X_dense[selected_labels == cluster_a]
        Xb = X_dense[selected_labels == cluster_b]
        for gene_idx, gene in enumerate(selected_genes):
            values_a = Xa[:, gene_idx]
            values_b = Xb[:, gene_idx]
            try:
                u_result = stats.mannwhitneyu(values_a, values_b, alternative="two-sided")
                p_value = float(u_result.pvalue)
                rank_biserial = float((2 * u_result.statistic / (len(values_a) * len(values_b))) - 1)
            except ValueError:
                p_value = np.nan
                rank_biserial = np.nan

            mean_a = float(np.mean(values_a))
            mean_b = float(np.mean(values_b))
            logfc = float(np.log2((mean_a + 1e-9) / (mean_b + 1e-9)))
            de_rows.append(
                {
                    "gene": gene,
                    "test": "mannwhitneyu",
                    "p_value": p_value,
                    "effect_size": rank_biserial,
                    "mean_cluster_a": mean_a,
                    "mean_cluster_b": mean_b,
                    "log2_fold_change": logfc,
                }
            )
    else:
        grouped = [X_dense[selected_labels == cluster] for cluster in clusters]
        for gene_idx, gene in enumerate(selected_genes):
            values = [group[:, gene_idx] for group in grouped]
            try:
                kw_result = stats.kruskal(*values)
                p_value = float(kw_result.pvalue)
                n = sum(len(v) for v in values)
                k = len(values)
                epsilon_squared = float((kw_result.statistic - k + 1) / (n - k)) if n > k else np.nan
            except ValueError:
                p_value = np.nan
                epsilon_squared = np.nan

            means = np.array([np.mean(v) for v in values], dtype=float)
            max_logfc = float(np.log2((means.max() + 1e-9) / (means.min() + 1e-9)))
            de_rows.append(
                {
                    "gene": gene,
                    "test": "kruskal",
                    "p_value": p_value,
                    "effect_size": epsilon_squared,
                    "min_cluster_mean": float(means.min()),
                    "max_cluster_mean": float(means.max()),
                    "max_log2_fold_change": max_logfc,
                }
            )

    de = pd.DataFrame(de_rows)
    de["adjusted_p_value"] = _bh_adjust(de["p_value"].to_numpy())
    if len(clusters) == 2:
        large_effect_mask = de["log2_fold_change"].abs() >= logfc_threshold
    else:
        large_effect_mask = de["max_log2_fold_change"].abs() >= logfc_threshold
    significant_mask = de["adjusted_p_value"] <= alpha
    n_significant_genes = int(significant_mask.sum())
    n_large_effect_genes = int((significant_mask & large_effect_mask).sum())
    de_combined_p = _combine_pvalues(de["p_value"].to_numpy())

    marker_agreement_stat = np.nan
    marker_agreement_p = np.nan
    if len(clusters) == 2:
        marker_agreement_stat = mean_pairwise_spearman
        marker_agreement_p = _combine_pvalues(pairwise["spearman_p_value"].to_numpy())
        marker_test_name = "pairwise_spearman"
    else:
        ranks = cluster_means_df.rank(axis=1, ascending=False).to_numpy()
        n_items = ranks.shape[1]
        n_raters = ranks.shape[0]
        rank_sums = ranks.sum(axis=0)
        s = np.sum((rank_sums - rank_sums.mean()) ** 2)
        if n_items > 1 and n_raters > 1:
            marker_agreement_stat = float(12 * s / (n_raters**2 * (n_items**3 - n_items)))
            chi2_stat = n_raters * (n_items - 1) * marker_agreement_stat
            marker_agreement_p = float(stats.chi2.sf(chi2_stat, df=n_items - 1))
        marker_test_name = "kendall_w"

    neighborhood_mixing = np.nan
    neighborhood_perm_p = np.nan
    if "connectivities" in adata.obsp:
        graph = adata.obsp["connectivities"][selected_cells][:, selected_cells]
        graph = graph.tocsr() if sparse.issparse(graph) else sparse.csr_matrix(graph)
        same_cluster = selected_labels[:, None] == selected_labels[None, :]
        total_weight = graph.sum()
        if total_weight > 0:
            same_weight = graph.multiply(same_cluster).sum()
            neighborhood_mixing = float(1 - same_weight / total_weight)

            if n_permutations > 0:
                permuted_stats = []
                for _ in range(n_permutations):
                    permuted_labels = rng.permutation(selected_labels)
                    perm_same = permuted_labels[:, None] == permuted_labels[None, :]
                    perm_same_weight = graph.multiply(perm_same).sum()
                    permuted_stats.append(1 - perm_same_weight / total_weight)
                permuted_stats = np.asarray(permuted_stats, dtype=float)
                neighborhood_perm_p = float((np.sum(permuted_stats >= neighborhood_mixing) + 1) / (len(permuted_stats) + 1))

    sample_chi2_p = np.nan
    sample_cramers_v = np.nan
    sample_table = None
    if sample_key is not None:
        if sample_key not in adata.obs:
            raise KeyError(f"{sample_key!r} is not in adata.obs.")
        sample_table = pd.crosstab(
            adata.obs.loc[selected_cells, cluster_key].astype(str),
            adata.obs.loc[selected_cells, sample_key].astype(str),
        ).reindex(clusters).fillna(0)
        if sample_table.shape[0] > 1 and sample_table.shape[1] > 1:
            chi2_result = stats.chi2_contingency(sample_table)
            sample_chi2_p = float(chi2_result.pvalue)
            sample_cramers_v = _cramers_v(sample_table)

    component_scores = pd.Series(
        {
            "expression_similarity": _similarity_score(min_pairwise_pearson),
            "marker_agreement": _similarity_score(marker_agreement_stat),
            "de_inverse_evidence": _inverse_count_score(n_large_effect_genes),
            "embedding_similarity": _distance_score(mean_embedding_distance),
            "neighborhood_mixing": neighborhood_mixing,
            "sample_consistency": 1 - sample_cramers_v if np.isfinite(sample_cramers_v) else np.nan,
        },
        dtype=float,
    )

    default_weights = {
        "expression_similarity": 0.25,
        "marker_agreement": 0.20,
        "de_inverse_evidence": 0.25,
        "embedding_similarity": 0.10,
        "neighborhood_mixing": 0.15,
        "sample_consistency": 0.05,
    }
    if weights is not None:
        default_weights.update(weights)

    weight_series = pd.Series(default_weights, dtype=float).reindex(component_scores.index)
    valid_scores = component_scores.notna() & weight_series.notna() & (weight_series > 0)
    if valid_scores.any():
        overall_score = float(
            np.average(component_scores.loc[valid_scores], weights=weight_series.loc[valid_scores])
        )
    else:
        overall_score = np.nan

    tests = pd.DataFrame(
        [
            {
                "test": "expression_profile_similarity",
                "statistic": "mean_pairwise_pearson",
                "value": mean_pairwise_pearson,
                "p_value": expression_perm_p,
                "interpretation": "Higher statistic supports merging; p-value tests unusually high similarity by permutation.",
            },
            {
                "test": "expression_profile_similarity",
                "statistic": "min_pairwise_pearson",
                "value": min_pairwise_pearson,
                "p_value": np.nan,
                "interpretation": "Lowest pairwise expression similarity inside the proposed merge group.",
            },
            {
                "test": "marker_score_agreement",
                "statistic": marker_test_name,
                "value": marker_agreement_stat,
                "p_value": marker_agreement_p,
                "interpretation": "Higher statistic supports merging by marker/rank agreement.",
            },
            {
                "test": "differential_expression",
                "statistic": "combined_gene_p_value",
                "value": de_combined_p,
                "p_value": de_combined_p,
                "interpretation": "Lower p-value indicates stronger evidence that at least one tested gene differs.",
            },
            {
                "test": "differential_expression",
                "statistic": "n_significant_genes",
                "value": n_significant_genes,
                "p_value": np.nan,
                "interpretation": f"Genes with BH-adjusted p-value <= {alpha}.",
            },
            {
                "test": "differential_expression",
                "statistic": "n_large_effect_genes",
                "value": n_large_effect_genes,
                "p_value": np.nan,
                "interpretation": f"Significant genes with absolute log2 fold-change >= {logfc_threshold}.",
            },
            {
                "test": "embedding_similarity",
                "statistic": "mean_centroid_distance",
                "value": mean_embedding_distance,
                "p_value": np.nan,
                "interpretation": f"Lower distance in adata.obsm[{use_rep!r}] supports merging.",
            },
            {
                "test": "neighborhood_mixing",
                "statistic": "cross_cluster_neighbor_fraction",
                "value": neighborhood_mixing,
                "p_value": neighborhood_perm_p,
                "interpretation": "Higher value means selected clusters are more mixed in the neighbor graph.",
            },
            {
                "test": "sample_composition",
                "statistic": "chi_square_p_value",
                "value": sample_chi2_p,
                "p_value": sample_chi2_p,
                "interpretation": "Lower p-value indicates sample/batch composition differs across selected clusters.",
            },
            {
                "test": "sample_composition",
                "statistic": "cramers_v",
                "value": sample_cramers_v,
                "p_value": np.nan,
                "interpretation": "Higher value indicates stronger sample/batch association with cluster labels.",
            },
        ]
    )

    summary = pd.DataFrame(
        [
            {
                "cluster_key": cluster_key,
                "clusters": ",".join(clusters),
                "n_clusters": len(clusters),
                "n_cells": int(selected_cells.sum()),
                "n_genes_tested": len(selected_genes),
                "mean_pairwise_pearson": mean_pairwise_pearson,
                "min_pairwise_pearson": min_pairwise_pearson,
                "mean_pairwise_spearman": mean_pairwise_spearman,
                "min_pairwise_spearman": min_pairwise_spearman,
                "n_significant_genes": n_significant_genes,
                "n_large_effect_genes": n_large_effect_genes,
                "neighborhood_mixing": neighborhood_mixing,
                "sample_cramers_v": sample_cramers_v,
                "overall_score": overall_score,
            }
        ]
    )

    return {
        "summary": summary,
        "tests": tests,
        "pairwise": pairwise,
        "differential_expression": de.sort_values("adjusted_p_value", na_position="last"),
        "component_scores": component_scores,
        "overall_score": overall_score,
        "cluster_sizes": cluster_sizes,
        "cluster_means": cluster_means_df,
        "sample_table": sample_table,
        "selected_genes": selected_genes,
    }
