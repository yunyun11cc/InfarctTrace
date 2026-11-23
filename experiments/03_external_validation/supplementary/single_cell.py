"""在 GSE163465 的三份未给药 CD45+ 文库中描述候选表达。

输入为原始归档、冻结九基因及聚类标签；输出描述性统计、UMAP 数据和运行记录。
"""

import gzip
import tarfile
from hashlib import file_digest
from importlib.metadata import version
import json

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.io import mmread
from scipy import sparse
from sklearn.neighbors import KNeighborsTransformer

from ..config import (
    FROZEN,
    GENES,
    SC_AUDIT_MARKERS,
    SC_HIGHLY_VARIABLE_GENES,
    SC_LABEL_THRESHOLD,
    SC_LEIDEN_ITERATIONS,
    SC_LEIDEN_RESOLUTION,
    SC_LIBRARY_SCALE,
    SC_MARKER_SCORE_FLOOR,
    SC_MARKERS,
    SC_MAX_GENES,
    SC_MAX_MT_PCT,
    SC_MIN_GENES,
    SC_N_NEIGHBORS,
    SC_N_PCS,
    SC_NEIGHBOR_JOBS,
    SEED,
    SINGLE_CELL_ARCHIVE,
    SINGLE_CELL_ARCHIVE_SHA256,
    SINGLE_CELL_LABELS,
    SINGLE_CELL_LIBRARIES as LIBRARIES,
)
from ..tables import read, save


class MatrixKNeighbors(KNeighborsTransformer):
    """为 Scanpy 1.11 的邻居读取器返回 CSR matrix。"""

    def fit_transform(self, X, y=None):
        """返回 Scanpy 读取器所需的 CSR matrix 邻接结果。"""
        return sparse.csr_matrix(super().fit_transform(X, y))


def load_libraries(output):
    """读取三份归档矩阵并通过逐细胞质控，返回合并数据和归档哈希。"""
    path = SINGLE_CELL_ARCHIVE
    with path.open('rb') as file:
        sha = file_digest(file, 'sha256').hexdigest()
    assert sha == SINGLE_CELL_ARCHIVE_SHA256
    datasets, qc = [], []
    with tarfile.open(path) as tar:
        for gsm, tag, time in LIBRARIES:
            prefix = f'{gsm}_{tag}'
            with gzip.open(tar.extractfile(f'{prefix}_features.tsv.gz'), 'rt') as f:
                features = pd.read_csv(
                    f, sep='\t', header=None, names=['id', 'symbol', 'feature_type']
                ).set_index('id')
            with gzip.open(tar.extractfile(f'{prefix}_barcodes.tsv.gz'), 'rt') as f:
                barcodes = [f'{gsm}:{line.strip()}' for line in f]
            with gzip.open(tar.extractfile(f'{prefix}_matrix.mtx.gz'), 'rb') as f:
                counts = sparse.csr_matrix(mmread(f).T, dtype=np.float32)
            cells = ad.AnnData(
                counts, obs=pd.DataFrame({'GSM': gsm, 'Time': time}, index=barcodes), var=features
            )
            cells.var['mt'] = cells.var.symbol.str.lower().str.startswith('mt-')
            sc.pp.calculate_qc_metrics(cells, qc_vars=['mt'], percent_top=None, inplace=True)
            keep = cells.obs.n_genes_by_counts.between(
                SC_MIN_GENES, SC_MAX_GENES
            ) & cells.obs.pct_counts_mt.lt(SC_MAX_MT_PCT)
            qc.append(dict(GSM=gsm, Time=time, Cells_input=cells.n_obs, Cells_retained=int(keep.sum())))
            datasets.append(cells[keep].copy())
    save(pd.DataFrame(qc), 'sc_library_qc', output)
    merged = ad.concat(datasets, join='inner', merge='same')
    merged.var['symbol'] = datasets[0].var.loc[merged.var_names, 'symbol']
    return merged, sha


def plot_location(embedding, clusters, candidate, dot, output):
    """将 UMAP 坐标、标签和表达统计生成定位图与检出率图。"""
    categories = dot.Cell_type.unique().tolist()
    fig, axes = plt.subplots(2, 5, figsize=(14, 6.4), layout='constrained')
    ax = axes.flat[0]
    for i, label in enumerate(categories):
        keep = clusters.eq(label)
        ax.scatter(*embedding[keep].T, s=3, color=plt.colormaps['tab20'](i), label=label, linewidths=0)
    ax.set_title('Marker-based cell labels')
    ax.legend(fontsize=5.8, markerscale=2, loc='best')
    for ax, gene in zip(axes.flat[1:], GENES):
        values = candidate[gene].to_numpy()
        order = np.argsort(values)
        data = ax.scatter(*embedding[order].T, c=values[order], cmap='magma_r', s=3, linewidths=0)
        ax.set_title(gene, fontstyle='italic')
        fig.colorbar(data, ax=ax, shrink=0.68, pad=0.015)
    for ax in axes.flat:
        ax.set_aspect('equal')
        ax.set_axis_off()
    fig.suptitle('GSE163465 | Frozen candidate expression on CD45+ cell UMAP')
    fig.savefig(output / 'figures/single_cell_umap_features.png', dpi=220)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 9.5), layout='constrained')
    row_names = list(dot[['Cell_type', 'Time']].drop_duplicates().itertuples(index=False, name=None))
    row_index = {key: i for i, key in enumerate(row_names)}
    y = [row_index[(r.Cell_type, r.Time)] for r in dot.itertuples()]
    x = dot.Gene.map({gene: i for i, gene in enumerate(GENES)})
    data = ax.scatter(
        x,
        y,
        s=dot.Pct_expressing * 2.1,
        c=dot.Mean_log1p_CP10K,
        cmap='viridis',
        edgecolors='#20242a',
        linewidths=0.25,
    )
    ax.set(
        xticks=range(len(GENES)),
        xticklabels=GENES,
        yticks=range(len(row_names)),
        yticklabels=[f'{cell} | {time}' for cell, time in row_names],
        title='GSE163465 | Candidate expression by cell label and library',
    )
    ax.set_ylim(len(row_names) - 0.5, -0.5)
    ax.set_xlim(-0.5, len(GENES) - 0.5)
    ax.tick_params(axis='x', rotation=35)
    ax.set_axisbelow(True)
    ax.grid(color='#edf0f3')
    legend_dots = [ax.scatter([], [], s=n * 2.1, color='#80959e', label=f'{n}%') for n in (10, 25, 50)]
    ax.legend(
        handles=legend_dots,
        title='Cells detected',
        loc='upper left',
        bbox_to_anchor=(1.02, 0.45),
        fontsize=8,
    )
    fig.colorbar(data, ax=ax, label='Mean log1p(CP10K)', shrink=0.45, anchor=(0, 1))
    fig.savefig(output / 'figures/single_cell_dotplot.png', dpi=220)
    plt.close(fig)


def run(output):
    """执行原 Scanpy 分析，写出每个文库和细胞群的描述性结果。"""
    for directory in ('tables', 'figures', 'reports'):
        (output / directory).mkdir(parents=True, exist_ok=True)
    cells, sha = load_libraries(output)
    frozen = read(FROZEN)
    ids = frozen['基因ID'].tolist()
    assert frozen['基因符号'].tolist() == GENES and set(ids) <= set(cells.var_names)
    sc.pp.normalize_total(cells, target_sum=SC_LIBRARY_SCALE)
    sc.pp.log1p(cells)
    marker_ids = {
        symbol: cells.var_names[cells.var.symbol.eq(symbol)].tolist()
        for symbol in {g for group in SC_MARKERS.values() for g in group}
    }
    assert all(len(found) == 1 for found in marker_ids.values())
    # 九个候选不参与选取聚类特征，也不用于身份注释。
    sc.pp.highly_variable_genes(
        cells, n_top_genes=SC_HIGHLY_VARIABLE_GENES, batch_key='GSM', flavor='seurat'
    )
    hvg = (
        cells.var.highly_variable
        & ~cells.var_names.isin(ids)
        & ~cells.var.symbol.str.lower().str.startswith('mt-')
    )
    sc.pp.pca(cells, n_comps=SC_N_PCS, mask_var=hvg.to_numpy(), random_state=SEED)
    print(f'{cells.n_obs:,} cells after QC. Building the 15-neighbor graph.', flush=True)
    sc.pp.neighbors(
        cells,
        n_neighbors=SC_N_NEIGHBORS,
        n_pcs=SC_N_PCS,
        random_state=SEED,
        transformer=MatrixKNeighbors(
            n_neighbors=SC_N_NEIGHBORS, algorithm='brute', n_jobs=SC_NEIGHBOR_JOBS
        ),
    )
    sc.tl.leiden(
        cells,
        resolution=SC_LEIDEN_RESOLUTION,
        flavor='igraph',
        n_iterations=SC_LEIDEN_ITERATIONS,
        directed=False,
        random_state=SEED,
    )
    print(f'{cells.obs.leiden.nunique()} Leiden clusters. Computing UMAP.', flush=True)
    sc.tl.umap(cells, random_state=SEED)
    cluster = cells.obs.leiden
    ms = pd.DataFrame(
        cells[:, [ids[0] for ids in marker_ids.values()]].X.toarray(),
        columns=list(marker_ids),
        index=cells.obs_names,
    )
    means = ms.groupby(cluster, observed=True).mean()
    extra = {}
    for symbol in SC_AUDIT_MARKERS:
        selection = cells.var.symbol.eq(symbol)
        if selection.sum() == 1:
            extra[symbol] = np.asarray(cells[:, selection.to_numpy()].X.toarray()).ravel()
    extra = pd.DataFrame(extra, index=cells.obs_names).groupby(cluster, observed=True).mean()
    means.join(extra).rename_axis('Cluster').reset_index().pipe(
        save, 'sc_cluster_marker_expression', output
    )
    labels = {}
    audit = []
    # 用相对全体细胞表达的平均标志物偏离量给出宽类群候选，随后按固定映射审核。
    profiles = means.sub(ms.mean()).div(ms.std().clip(lower=SC_MARKER_SCORE_FLOOR))
    for name, values in SC_MARKERS.items():
        profiles[name] = profiles[values].mean(axis=1)
    scoring = profiles[list(SC_MARKERS)]
    annotations = read(SINGLE_CELL_LABELS, dtype={'Cluster': str}).set_index('Cluster')
    assert set(annotations.index) == set(means.index)
    for cid, row in scoring.iterrows():
        ranked = row.sort_values(ascending=False)
        suggested = ranked.index[0] if ranked.iloc[0] >= SC_LABEL_THRESHOLD else 'Unresolved'
        labels[cid] = annotations.loc[cid, 'Cell_type']
        audit.append(
            dict(
                Cluster=cid,
                N_cells=int(cluster.eq(cid).sum()),
                Suggested=suggested,
                Cell_type=labels[cid],
                Best_score=ranked.iloc[0],
                Second_label=ranked.index[1],
                Second_score=ranked.iloc[1],
                Annotation_basis=annotations.loc[cid, 'Annotation_basis'],
            )
        )
    save(pd.DataFrame(audit), 'sc_cluster_labels', output)
    cell_types = cluster.astype(str).map(labels)
    expression = pd.DataFrame(cells[:, ids].X.toarray(), index=cells.obs_names, columns=GENES)
    records = []
    for (gsm, time, label), values in expression.groupby(
        [cells.obs.GSM, cells.obs.Time, cell_types], sort=False, observed=True
    ):
        for gene in GENES:
            records.append(
                dict(
                    GSM=gsm,
                    Time=time,
                    Cell_type=label,
                    Gene=gene,
                    Cells=len(values),
                    Pct_expressing=100 * values[gene].gt(0).mean(),
                    Mean_log1p_CP10K=values[gene].mean(),
                )
            )
    dot = pd.DataFrame(records)
    dot['Time'] = pd.Categorical(dot.Time, ['Sham', 'D3', 'D7'], ordered=True)
    dot = dot.sort_values(['Cell_type', 'Time'])
    save(dot, 'sc_expression_by_group', output)
    # 每时点仅一份公开文库，保存组内描述性变化。
    ref = dot[dot.Time.eq('Sham')][['Cell_type', 'Gene', 'Pct_expressing', 'Mean_log1p_CP10K']]
    changes = dot[dot.Time.ne('Sham')].merge(
        ref, on=['Cell_type', 'Gene'], how='left', suffixes=('', '_Sham')
    )
    changes['Delta_detected_pp'] = changes.Pct_expressing - changes.Pct_expressing_Sham
    changes['Delta_mean_log1p_CP10K'] = changes.Mean_log1p_CP10K - changes.Mean_log1p_CP10K_Sham
    save(changes, 'sc_time_deltas', output)
    coords = pd.DataFrame(cells.obsm['X_umap'], index=cells.obs_names, columns=['UMAP1', 'UMAP2'])
    records = (
        coords.join(cells.obs[['GSM', 'Time', 'leiden']]).assign(Cell_type=cell_types).join(expression)
    )
    records.rename_axis('Cell_ID').reset_index().pipe(save, 'sc_umap_candidate_expression', output)
    plot_location(coords.to_numpy(), cell_types, expression, dot, output)
    (output / 'reports/sc_provenance.json').write_text(
        json.dumps(
            dict(
                archive_sha256=sha,
                scanpy=version('scanpy'),
                seed=SEED,
                retained_cells=cells.n_obs,
                highly_variable_features=int(hvg.sum()),
                Leiden_resolution=SC_LEIDEN_RESOLUTION,
                libraries=LIBRARIES,
                inferential_unit='GSM',
            ),
            indent=2,
        )
        + '\n',
        encoding='utf-8',
    )
    print(f'{cells.n_obs:,} cells, {len(labels)} clusters', flush=True)
