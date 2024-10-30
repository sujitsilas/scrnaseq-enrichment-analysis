# scRNA-seq Enrichment Analysis for Olfr2+ and Olfr2− Macrophages
This repository contains the script used for cross-validation and validation of Olfr2+ expressing macrophages and Olfr2− macrophages in a single-cell RNA-sequencing (scRNA-seq) dataset, as described in [Armstrong Suthahar et al., 2024](https://doi.org/10.1093/cvr/cvae153). The goal of the analysis is to identify cells expressing these specific signature genes within an integrated mouse and human dataset.

We employed the R package escape (v1.8.0) to perform Gene Set Enrichment Analysis (GSEA) using the enrichIt function, which calculates enrichment scores (ES) for each cell. These scores are used to group and plot cells based on their cell-type annotations, enabling the identification of enriched cell types and visualizing the results via feature plots.

# Key Features:
Enrichment analysis using the escape package's enrichIt function.
Identification of enriched Olfr2+ and Olfr2− macrophages.
Visualization of enriched cells through feature plots.
Cross-validation of results within the integrated mouse and human dataset.
For more details on how this method was implemented, refer to the original paper: [Armstrong Suthahar et al., 2024](https://doi.org/10.1093/cvr/cvae153).
