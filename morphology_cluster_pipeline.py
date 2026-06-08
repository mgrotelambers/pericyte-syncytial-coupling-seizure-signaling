"""
morphology_pipeline.py — Pericyte morphological clustering pipeline.

Classes and functions for morphological feature analysis, K-means clustering
with NaN support, silhouette evaluation, feature selection, and visualisation
(UMAP, PCA, boxplots).

Usage:
    from morphology_pipeline import *

Dependencies:
    numpy, pandas, matplotlib, seaborn, sklearn, umap, KmeansWithNulls
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import seaborn as sns
from itertools import combinations
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import silhouette_score, silhouette_samples
from sklearn.decomposition import PCA
from KmeansWithNulls import KmeansWithNulls
import umap

from fig_utils import save_if_enabled


# =============================================================================
# 1 — EXPLORATORY DATA ANALYSIS
# =============================================================================

def plot_feature_distributions(df, features, kind='box', figsize=None):
    """
    Plot boxplots or histograms for a list of morphological features.

    Parameters
    ----------
    df : DataFrame
        Source data.
    features : list of str
        Column names to plot.
    kind : 'box' or 'hist'
        Plot type.
    figsize : tuple, optional
        Figure size. Defaults to (3*n_features, 5).
    """
    n = len(features)
    figsize = figsize or (3 * n, 5)
    colors = ['lightgreen', 'skyblue', 'lightyellow', 'salmon',
              'gray', 'purple', 'violet', 'orange']

    fig, axes = plt.subplots(1, n, figsize=figsize)
    if n == 1:
        axes = [axes]

    for ax, feat, color in zip(axes, features, colors[:n]):
        data = df[feat].dropna() if feat in df.columns else pd.Series(dtype=float)
        if kind == 'box':
            sns.boxplot(data=data, color=color, ax=ax)
            ax.set_xlabel('')        
            ax.set_ylabel('')
        else:
            sns.histplot(data=data, color=color, ax=ax)
            ax.set_xlabel('')   
            ax.set_ylabel('Count')
        ax.set_title(feat)           

    plt.tight_layout()
    plt.show()


def prepare_derived_features(df):
    """
    Add SurfaceArea/Volume_Ratio and Raw_voxel_count columns in-place.

    Returns the modified DataFrame.
    """
    df['SurfaceArea/Volume_Ratio'] = df['SurfaceArea'] / df['Volume']
    df['Raw_voxel_count'] = (df['# Slab voxels'] +
                             df['# Junction voxels'] +
                             df['# End-point voxels'])
    return df

# =============================================================================
# 2 — OPTIMAL k ANALYSIS
# =============================================================================

def silhouette_analysis(X_scaled, k_list=(2, 3, 4, 5, 6), max_k=12):
    """
    For each k in k_list, plot per-cluster silhouette profiles + PCA scatter.
    Then plot summary elbow + silhouette curves up to max_k.

    Parameters
    ----------
    X_scaled : array
        Standardised feature matrix.
    k_list : tuple/list of int
        Cluster counts for detailed silhouette plots.
    max_k : int
        Upper bound for summary plots.
    """
    pca = PCA(n_components=2)
    X_pca = pca.fit_transform(np.nan_to_num(X_scaled))
    mask = ~np.isnan(X_scaled).any(axis=1)

    # Detailed per-k plots
    for n_clusters in k_list:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 7))
        ax1.set_xlim([-0.1, 1])
        ax1.set_ylim([0, len(X_scaled) + (n_clusters + 1) * 10])

        km = KmeansWithNulls(n_clusters=n_clusters, max_iter=300, random_state=42)
        km.fit(X_scaled)
        labels = km.predict(X_scaled)

        if np.sum(mask) > 1:
            sil_avg = silhouette_score(X_scaled[mask], labels[mask])
            sample_sil = silhouette_samples(X_scaled[mask], labels[mask])
            print(f"k={n_clusters}: avg silhouette = {sil_avg:.3f}")

            y_lower = 10
            for i in range(n_clusters):
                vals = np.sort(sample_sil[labels[mask] == i])
                y_upper = y_lower + len(vals)
                color = cm.nipy_spectral(float(i) / n_clusters)
                ax1.fill_betweenx(np.arange(y_lower, y_upper), 0, vals,
                                  facecolor=color, edgecolor=color, alpha=0.7)
                ax1.text(-0.05, y_lower + 0.5 * len(vals), str(i))
                y_lower = y_upper + 10
            ax1.axvline(x=sil_avg, color='red', ls='--')

        ax1.set_title('Silhouette plot'); ax1.set_xlabel('Silhouette coefficient')
        ax1.set_ylabel('Cluster'); ax1.set_yticks([])

        colors = cm.nipy_spectral(labels.astype(float) / n_clusters)
        ax2.scatter(X_pca[mask, 0], X_pca[mask, 1], marker='.', s=30,
                    lw=0, alpha=0.7, c=colors[mask], edgecolor='k')
        centers_pca = pca.transform(np.nan_to_num(km.centroids))
        ax2.scatter(centers_pca[:, 0], centers_pca[:, 1], marker='o',
                    c='white', s=200, edgecolor='k')
        for i, c in enumerate(centers_pca):
            ax2.scatter(c[0], c[1], marker=f'${i}$', s=50, edgecolor='k')
        ax2.set_title('PCA projection')
        ax2.set_xlabel('PC1'); ax2.set_ylabel('PC2')

        plt.suptitle(f'Silhouette analysis — k = {n_clusters}',
                     fontsize=14, fontweight='bold')
        plt.tight_layout(); plt.show()

    # Summary plots
    def _calc_inertia(data, labels, centroids):
        total = 0.0
        for i in range(len(centroids)):
            pts = data[labels == i]
            valid = pts[~np.isnan(pts).any(axis=1)]
            if len(valid) > 0:
                total += np.sum((valid - centroids[i]) ** 2)
        return total

    k_vals = range(2, max_k + 1)
    inertias, sil_scores = [], []
    for k in k_vals:
        km = KmeansWithNulls(n_clusters=k, random_state=42, max_iter=300)
        km.fit(X_scaled)
        labels = km.predict(X_scaled)
        if np.sum(mask) > 1:
            sil_scores.append(silhouette_score(X_scaled[mask], labels[mask]))
        else:
            sil_scores.append(np.nan)
        inertias.append(_calc_inertia(X_scaled, labels, km.centroids))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    ax1.plot(k_vals, inertias, 'bo-', lw=2, ms=8)
    ax1.set_title('Elbow Method'); ax1.set_xlabel('k'); ax1.set_ylabel('Inertia')
    ax1.grid(True, ls='--', alpha=0.7); ax1.set_xticks(list(k_vals))
    ax2.plot(k_vals, sil_scores, 'ro-', lw=2, ms=8)
    ax2.set_title('Silhouette Score'); ax2.set_xlabel('k'); ax2.set_ylabel('Score')
    ax2.grid(True, ls='--', alpha=0.7); ax2.set_xticks(list(k_vals))
    plt.tight_layout(); plt.show()


# =============================================================================
# 3 — MORPHOLOGY ANALYSER
# =============================================================================

class MorphologyAnalyzer:
    """Cluster pericytes by morphological features, visualise with UMAP/PCA.

    Parameters
    ----------
    features : list of str, optional
        Morphological feature columns. Defaults to the 4-feature base set.
    colors : list of str, optional
        Cluster colour palette.
    umap_kwargs : dict, optional
        Override UMAP hyperparameters.
    """

    DEFAULT_FEATURES = [
        'Longest Shortest Path',
        'Elli.R1/R3',
        'Sphericity',
        '# Junctions',
    ]

    DEFAULT_COLORS = ['#FF7D33', '#d62728', '#00c5cd', '#8C1E95']

    IV_COLORS = {0: '#00008B', 1: '#ADD8E6'}

    def __init__(self, features=None, colors=None, umap_kwargs=None):
        self.FEATURES = features or self.DEFAULT_FEATURES
        self.COLORS = colors or self.DEFAULT_COLORS
        self.scaler = StandardScaler()
        _umap_kw = dict(n_neighbors=25, min_dist=0.0005, spread=2, random_state=42, n_jobs=1)  # n_jobs=1 for reproducibility
        if umap_kwargs:
            _umap_kw.update(umap_kwargs)
        self.umap_reducer = umap.UMAP(**_umap_kw)
        self.pca = PCA(n_components=2)

    def preprocess_data(self, df):
        """Add derived columns if source columns are present, then return standardised feature matrix."""
        required = {'SurfaceArea', 'Volume', '# Slab voxels', '# Junction voxels', '# End-point voxels'}
        if required.issubset(df.columns):
            prepare_derived_features(df)
        return self.scaler.fit_transform(df[self.FEATURES].values)

    def analyze(self, df, n_clusters):
        """Run full morphological cluster analysis.

        Returns dict with features, labels, UMAP embedding, PCA, etc.
        """
        features = self.preprocess_data(df)

        km = KmeansWithNulls(n_clusters=n_clusters, random_state=42)
        km.fit(features)
        labels = km.predict(features)

        df['Morphological_Cluster'] = labels + 1

        clean_mask = ~np.isnan(features).any(axis=1)
        features_clean = features[clean_mask]
        labels_clean = labels[clean_mask]

        np.random.seed(42)  # force full UMAP reproducibility
        embedding = self.umap_reducer.fit_transform(features_clean)
        pca_result = self.pca.fit_transform(features_clean)

        return {
            'features': features,
            'features_clean': features_clean,
            'labels': labels,
            'labels_clean': labels_clean,
            'embedding': embedding,
            'pca_result': pca_result,
            'n_clusters': n_clusters,
            'explained_variance_ratio': self.pca.explained_variance_ratio_,
            'clean_indices': clean_mask,
        }

    # ── Plotting ─────────────────────────────────────────────────────────

    def plot_results(self, results, df, save_path_umap=None, save_path_boxplots=None,
                    cluster_order=None, cluster_labels=None):
        n_clusters = results['n_clusters']
        cluster_names = [f'Cluster {i+1}' for i in range(n_clusters)]
        palette = {name: c for name, c in zip(cluster_names, self.COLORS)}

        # Fall back to 'Cluster 1' etc. if no labels provided
        if cluster_labels is None:
            cluster_labels = cluster_names

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
        for i, color in enumerate(self.COLORS[:n_clusters]):
            m = results['labels_clean'] == i
            ax1.scatter(results['embedding'][m, 0], results['embedding'][m, 1],
                        c=color, label=cluster_labels[i])
            ax2.scatter(results['pca_result'][m, 0], results['pca_result'][m, 1],
                        c=color, label=cluster_labels[i])
        ax1.set_title('UMAP projection')
        ax1.legend(fontsize=8, markerscale=1, loc='lower left')
        ax1.set_xlabel('UMAP1'); ax1.set_ylabel('UMAP2')
        evr = results['explained_variance_ratio'].sum()
        ax2.set_title(f'PCA projection (explained var: {evr:.2f})')
        ax2.legend(fontsize=8, markerscale=1, loc='upper right')
        ax2.set_xlabel('PC1'); ax2.set_ylabel('PC2')
        plt.tight_layout()
        save_if_enabled(save_path_umap)
        plt.show()

        self._plot_feature_boxplots(
            results['features'], results['labels'], n_clusters, palette,
            cluster_order=cluster_order,
            save_path=save_path_boxplots)


    def _plot_feature_boxplots(self, features, labels, n_clusters, palette,
                            cluster_order=None, save_path=None):
        """Boxplots + histograms of standardised features per cluster."""
        if cluster_order is None:
            cluster_order = [f'Cluster {i+1}' for i in range(n_clusters)]

        plot_data = []
        for i, feat in enumerate(self.FEATURES):
            for val, lbl in zip(features[:, i], labels):
                if not np.isnan(val):
                    plot_data.append({'Feature': feat, 'Value': val,
                                    'Cluster': f'Cluster {int(lbl)+1}'})
        plot_df = pd.DataFrame(plot_data)

        n_feat = len(self.FEATURES)

        fig = plt.figure(figsize=(18, 12), constrained_layout=True)
        gs = fig.add_gridspec(2, n_feat, height_ratios=[15, 3])

        for idx, feat in enumerate(self.FEATURES):
            ax = fig.add_subplot(gs[0, idx])
            sns.boxplot(data=plot_df[plot_df['Feature'] == feat],
                        x='Cluster', y='Value', hue='Cluster',
                        order=cluster_order, hue_order=cluster_order,
                        palette=palette, legend=False, ax=ax)
            ax.set_title(feat)
            ax.set_xlabel('Cluster')
            ax.set_xticks(range(len(cluster_order)))
            ax.set_xticklabels([str(i+1) for i in range(len(cluster_order))])
            ax.set_ylabel('Standardized Values (z-scores)')
            ax.grid(True, alpha=0.2)

            ax_h = fig.add_subplot(gs[1, idx])
            sns.histplot(data=plot_df[plot_df['Feature'] == feat]['Value'],
                        bins=30, color='lightgray', edgecolor='black', ax=ax_h)    
            ax_h.set_xlabel('Standardized Values (z-scores)')
            ax_h.set_ylabel('Count')
            ax_h.grid(True, alpha=0.2)

        save_if_enabled(save_path)
        plt.show()

    def plot_results_with_iv(self, results, morpho_df, iv_df, save_path=None):
        """UMAP coloured by IV-clustering class (overlay)."""
        clean_mask = results['clean_indices']
        morpho_clean = morpho_df[clean_mask].copy()

        merged = pd.merge(morpho_clean, iv_df, on='VC_IV',
                          how='left', suffixes=('_morpho', '_iv'))
        unmatched = merged['cluster'].isna()
        if unmatched.any():
            print(f"{unmatched.sum()} cells excluded (no IV match)")

        matched = ~unmatched
        merged_m = merged[matched]
        emb_m = results['embedding'][matched]

        fig, ax = plt.subplots(figsize=(8, 6))
        for iv_class in sorted(self.IV_COLORS):
            m = merged_m['cluster'] == iv_class
            ax.scatter(emb_m[m, 0], emb_m[m, 1],
                    c=self.IV_COLORS[iv_class], label=f'Cluster {iv_class + 1}')
        ax.set_title('UMAP — coloured by IV cluster')
        ax.set_xlabel('UMAP1'); ax.set_ylabel('UMAP2')
        ax.legend()
        plt.tight_layout()
        save_if_enabled(save_path)
        plt.show()

# =============================================================================
# 4 — FEATURE SELECTION
# =============================================================================

def test_feature_additions(df, base_features, all_features, optimal_k):
    """Test every combination of remaining features added to the base set.

    Returns sorted list of result dicts (best first).
    """
    remaining = [f for f in all_features if f not in base_features]

    def _score(feats):
        try:
            X = df[feats].values
            if X.shape[0] < optimal_k or X.shape[1] < 2:
                return -1
            X_sc = StandardScaler().fit_transform(X)
            km = KmeansWithNulls(n_clusters=optimal_k, random_state=42)
            km.fit(X_sc)
            labels = km.predict(X_sc)
            mask = ~np.isnan(X_sc).any(axis=1)
            return silhouette_score(X_sc[mask], labels[mask])
        except Exception as e:
            print(f"  failed for {feats}: {e}")
            return -1

    results = [{'features': base_features,
                'silhouette_score': _score(base_features)}]

    for r in range(1, len(remaining) + 1):
        for combo in combinations(remaining, r):
            feats = base_features + list(combo)
            s = _score(feats)
            if s > -1:
                results.append({'features': feats, 'silhouette_score': s})

    results.sort(key=lambda x: x['silhouette_score'], reverse=True)
    return results


def report_feature_additions(results, base_features):
    """Print silhouette scores for all feature combinations, sorted best first."""
    base_score = next(r['silhouette_score'] for r in results
                      if sorted(r['features']) == sorted(base_features))

    print(f"Base score ({', '.join(base_features)}): {base_score:.4f}\n")
    for r in results:
        if sorted(r['features']) == sorted(base_features):
            print(f"  {base_score:.4f}  —  base features only")
        else:
            added = [f for f in r['features'] if f not in base_features]
            delta = r['silhouette_score'] - base_score
            print(f"  {r['silhouette_score']:.4f} ({delta:+.4f})  —  + {added}")