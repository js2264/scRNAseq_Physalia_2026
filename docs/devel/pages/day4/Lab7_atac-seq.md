# 13  Lab 7: Single-cell ATAC-seq analysis workflow

> **NOTE:**
>
> The estimated time for this lab is around 1h10.
>
> This lab relies on the `fragments` file of a 10x Genomics scATAC-seq experiment (~2 GB): its code is not run when the book is built, and no output is shown on this page. Run it in your workshop session!

> **TIP:**
>
> - Import cells from human PBMC scATACseq dataset
> - Perform scATACseq quality controls and checks
> - Filter, normalize and plot PBMC scATACseq dataset
> - Compute gene activity scores, check known markers and annotate scATAC clusters
> - Perform differential peak accessibility analysis

ATAC-seq data may be obtained in isolation using a single-cell ATAC-seq protocol (e.g. [10X scATACseq](https://www.10xgenomics.com/products/single-cell-atac)) or in combination with gene expression data using a single-cell multiome protocole (e.g. [10X multiome](https://www.10xgenomics.com/products/single-cell-multiome-atac-plus-gene-expression) and [SHARE-seq](https://pubmed.ncbi.nlm.nih.gov/33098772/)).

Several packages are currently avaialble to process scATAC-seq data in R. These include [Signac](https://stuartlab.org/signac/) and [ArchR](https://www.archrproject.com/bookdown/index.html). This lab will closely follow the [processing steps](https://stuartlab.org/signac/articles/pbmc_vignette) outlined in `Signac`, which interfaces well with `Seurat` for single-cell analysis.

In this lab, we will process a PBMC single-cell ATAC-seq (scATAC-seq) dataset and perform preliminary analysis to assess quality of these data. The data for this lab comes from [10X Genomics](https://www.10xgenomics.com/datasets/10-k-peripheral-blood-mononuclear-cells-pbm-cs-from-a-healthy-donor-1-standard-1-0-1). The dataset contains roughly ~ 10,000 cells.

## 13.1 Process human PBMC dataset

### 13.1.1 Download data

The files we need are directly available from the 10x Genomics website. We download them in a `~/data/scATAC` folder (this may take a few minutes, the `fragments` file weighs ~2 GB). Note that the `fragments` file and its index (`.tbi`) must sit next to each other, with matching names.

``` r
scatac_dir <- file.path(path.expand("~"), "data", "scATAC")
dir.create(scatac_dir, recursive = TRUE, showWarnings = FALSE)
options(timeout = 3600)
url <- "https://cf.10xgenomics.com/samples/cell-atac/1.0.1/atac_v1_pbmc_10k/"
files <- c(
    "atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5",
    "atac_v1_pbmc_10k_peaks.bed",
    "atac_v1_pbmc_10k_singlecell.csv",
    "atac_v1_pbmc_10k_fragments.tsv.gz",
    "atac_v1_pbmc_10k_fragments.tsv.gz.tbi"
)
for (f in files) {
    if (!file.exists(file.path(scatac_dir, f)))
        download.file(paste0(url, f), destfile = file.path(scatac_dir, f))
}
list.files(scatac_dir)
h5_file <- file.path(scatac_dir, "atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5")
peaks_file <- file.path(scatac_dir, "atac_v1_pbmc_10k_peaks.bed")
metadata_file <- file.path(scatac_dir, "atac_v1_pbmc_10k_singlecell.csv")
fragments_file <- file.path(scatac_dir, "atac_v1_pbmc_10k_fragments.tsv.gz")
```

### 13.1.2 Import data

Notice how the count matrix is in a `.h5` format. We have already encountered this kind of 10x Genomics outputs in [Lab 3](../../pages/day2/Lab3_Rbioc.md). Back then, we imported them with [`DropletUtils::read10xCounts`](https://rdrr.io/pkg/DropletUtils/man/read10xCounts.html).

> **NOTE:**
>
> Does this function work here?
>
> > **NOTE:**
> >
> > ``` r
> > scatac <- DropletUtils::read10xCounts(h5_file)
> > scatac
> > ```

This works because 10X Genomics make sure to distribute files in `.h5` format that are consistent accross single-cell sequencing methodologies.

However, the `SingleCellExperiment` obtained with this approach is not the most convenient, as it cannot natively leverage `fragments` file (see below).

Instead, we can create a `Signac` object, a flavour of `Seurat` objects.

> **NOTE:**
>
> Import counts matrix and feature annotations using an import function provided by `Seurat`.
>
> > **NOTE:**
> >
> > ``` r
> > library(Seurat)
> > library(Signac)
> > library(rtracklayer)
> > library(stringr)
> > cnts <- Read10X_h5(h5_file)
> > features <- import(peaks_file)
> > features$peak <- as.character(features) |> str_replace(':', '-')
> > metadata <- read.csv(
> >   file = metadata_file,
> >   header = TRUE,
> >   row.names = 1
> > )
> > ```
>
> How many accessible genomic segments were found in this dataset?
>
> > **NOTE:**
> >
> > ``` r
> > features
> > length(features)
> > ```

### 13.1.3 Create a `Seurat` object

The next step is to aggregate `counts` and `features` into a `ChromatinAssay`, a scATAC-seq flavour of `Seurat` standard `Assays`. The documentation for [`?CreateChromatinAssay`](https://stuartlab.org/signac/reference/CreateChromatinAssay.html) indicates that the user can provide:

1.  A `fragments` file, corresponding to the full list of all unique fragments mapped across all single cells.
2.  Genomic annotations to the `ChromatinAssay`, corresponding to gene annotations, promoter positions, etc. Such annotations can be generated from Ensembl.

> **NOTE:**
>
> Generate human annotations from Ensembl using a parsing function from `Signac`.
>
> > **NOTE:**
> >
> > ``` r
> > ## - Get human gene annotations (hg19/GRCh37) to feed it into the future `ChromatinAssay`
> > annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75)
> > GenomeInfoDb::seqlevelsStyle(annotations) <- 'UCSC'
> > GenomeInfoDb::genome(annotations) <- "hg19"
> > ```
>
> Create a `ChromatinAssay` using `counts`, `features`, `fragments` and `annotations`.
>
> > **NOTE:**
> >
> > ``` r
> > ## - Create Chromatin Assay
> > assay <- Signac::CreateChromatinAssay(
> >     counts = cnts,
> >     ranges = features,
> >     fragments = fragments_file,
> >     annotation = annotations,
> >     genome = "hg19",
> >     min.cells = 10,
> >     min.features = 10
> > )
> > assay
> > ```
>
> What are the dimensions of this object? Are they comparable to the count matrix? Comment.
>
> > **NOTE:**
> >
> > ``` r
> > dim(cnts)
> > dim(assay)
> > ```

It’s finally time to wrap the `ChromatinAssay` into a `Seurat` standard object! This is done using the `CreateSeuratObject`.

> **NOTE:**
>
> Create a `Seurat` object from the `ChromatinAssay` object.
>
> > **NOTE:**
> >
> > ``` r
> > ## - Create Seurat object
> > PBMC <- Seurat::CreateSeuratObject(
> >     counts = assay,
> >     assay = 'ATAC',
> >     meta.data = metadata[metadata$is__cell_barcode == 1, ]
> > )
> > PBMC
> >
> > PBMC[['ATAC']]
> > granges(PBMC)
> > Annotation(PBMC)
> >
> > ## - Free some memory: the intermediate objects are no longer needed
> > rm(scatac, cnts, assay); invisible(gc())
> > ```

### 13.1.4 Check QCs

#### 13.1.4.1 Cell-based QCs

The fraction of reads in peaks (FRiP) is a good indicator of how well each cell was handled during scATACseq processing.

``` r
PBMC$FRiP <- PBMC$peak_region_fragments / PBMC$passed_filters
PBMC$nCount_ATAC <- colSums(GetAssayData(PBMC, layer = "counts"))
PBMC$nFeature_ATAC <- colSums(GetAssayData(PBMC, layer = "counts") > 0)

quantile(PBMC$FRiP, seq(0, 1, 0.1))
quantile(PBMC$nCount_ATAC, seq(0, 1, 0.1))
quantile(PBMC$nFeature_ATAC, seq(0, 1, 0.1))
```

#### 13.1.4.2 Peaks-based QCs

> **NOTE:**
>
> Which analysis are the `fragments` required for, exactly?
>
> Could we still perform normalization/clustering/annotation without them? And motif enrichment analysis?
>
> Since we do have the `fragments` file at hand, most of the QC steps are available (e.g. `TSSEnrichment`, `NucleosomeSignal` or fragment size distribution). Let’s go through them one by one.
>
> > **NOTE:**
> >
> > ``` r
> > # compute nucleosome signal score per cell
> > PBMC <- NucleosomeSignal(object = PBMC)
> > PBMC
> >
> > # compute TSS enrichment score per cell
> > PBMC <- Signac::TSSEnrichment(object = PBMC)
> > PBMC
> >
> > # plot TSS enrichment vs. number of fragments in peaks
> > DensityScatter(PBMC, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
> > ```

Let’s compare the distributions of the different QC metrics, for cells with a high or a low TSS enrichment score.

``` r
PBMC$high.tss <- ifelse(PBMC$TSS.enrichment > 2.5, 'High', 'Low')
table(PBMC$high.tss)
VlnPlot(
    object = PBMC,
    features = c('nCount_ATAC', 'FRiP', 'TSS.enrichment', 'nucleosome_signal'),
    group.by = 'high.tss',
    pt.size = 0,
    ncol = 4
)
```

> **TIP:**
>
> The [`TSSPlot()`](https://stuartlab.org/signac/reference/TSSPlot.html) function from `Signac` plots the average accessibility profile around TSSs, for groups of cells. It requires running [`TSSEnrichment()`](https://stuartlab.org/signac/reference/TSSEnrichment.html) with `fast = FALSE`, which stores the accessibility of each position around each TSS: this takes a lot of memory (\> 16 GB for this dataset), so we do not run it in the workshop sessions.

The `FragmentHistogram` function from `Signac` can be used to plot the fragment size distribution in peaks with different nucleosome signals.

> **NOTE:**
>
> Group genomic loci by their nucleosome signal (\> 4 or \< 4) and plot the fragment size distribution.
>
> > **NOTE:**
> >
> > ``` r
> > PBMC$nucleosome_group <- ifelse(PBMC$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
> > FragmentHistogram(object = PBMC, group.by = 'nucleosome_group')
> > ```

### 13.1.5 Filter cells and features

Just like standard scRNAseq data, scATACseq data needs to be filtered to remove low-quality cells and peaks.

> **NOTE:**
>
> Filter the `Seurat` object (cells and features) as deemed appropriate.
>
> > **NOTE:**
> >
> > ``` r
> > ## - Filter data
> > PBMC <- subset(PBMC, subset = nCount_ATAC > 3000 & nCount_ATAC < 100000)
> > PBMC <- subset(PBMC, subset = nFeature_ATAC > 1000 & nFeature_ATAC < 20000)
> > PBMC <- subset(PBMC, subset = FRiP > 0.30)
> > PBMC <- subset(PBMC, subset = TSS.enrichment > 2.5 & nucleosome_signal < 4)
> > PBMC
> >
> > ## - Remove peaks with low coverage
> > PBMC <- PBMC[rowSums(GetAssayData(PBMC, layer = "counts")) > 10, ]
> > PBMC <- PBMC[rowSums(GetAssayData(PBMC, layer = "counts") > 0) > 10, ]
> > PBMC
> > ```

### 13.1.6 Dimensionality reduction and clustering

Now that the dataset is filtered, we can proceed to data normalization, dimensionality reduction and clustering. Compared to scRNAseq, scATACseq data is much sparser and has a higher dimensionality. This has important implications for the choice of normalization and dimensionality reduction methods. scATACseq data are generally normalized using a `TF-IDF` normalization, which is implemented in the `RunTFIDF` function in `Signac`, and dimensionality reduction is performed using `LSI` (Latent Semantic Indexing) implemented in the `RunSVD` function.

> **NOTE:**
>
> Normalize then further reduce the dimensionality for visualization purposes.
>
> > **NOTE:**
> >
> > ``` r
> > ## - Normalize data
> > PBMC <- RunTFIDF(PBMC)
> >
> > ## - Reduce dimensionality
> > PBMC <- FindTopFeatures(PBMC, min.cutoff = 'q50')
> > PBMC <- RunSVD(PBMC)
> >
> > ## - The first LSI component is often correlated with sequencing depth
> > DepthCor(PBMC)
> >
> > ## - Embed in 2D UMAP
> > PBMC <- RunUMAP(object = PBMC, reduction = 'lsi', dims = 2:30)
> >
> > ## - Cluster cells
> > PBMC <- FindNeighbors(object = PBMC, reduction = 'lsi', dims = 2:30)
> > PBMC <- FindClusters(object = PBMC, verbose = FALSE, algorithm = 3)
> >
> > ## - Visualize data
> > DimPlot(PBMC, label = TRUE) + NoLegend()
> > ```
>
> What can you observe in the `UMAP` projection of the dataset? Comment on the separation of some cell types into different spatially-resolved clusters.

## 13.2 Compute gene activity scores

Signac’s [`GeneActivity()`](https://stuartlab.org/signac/reference/GeneActivity.html) function require scATACseq fragment information. Since we have them, we can estimate a gene activity score for each gene in the `annotations`.

``` r
gene.activities <- GeneActivity(PBMC)
```

We can now save this new object as an `Assay` in the `PBMC` object and normalize it.

``` r
PBMC[['RNA']] <- CreateAssayObject(counts = gene.activities)

# - Normalize the new RNA assay.
PBMC <- NormalizeData(
    object = PBMC,
    assay = 'RNA',
    normalization.method = 'LogNormalize',
    scale.factor = median(PBMC$nCount_RNA)
)
PBMC
```

One can now visualize expression of individual markers across clusters to infer cluster identity.

> **NOTE:**
>
> Assess the expression of known markers across clusters, projecting the cells in a 2D UMAP.
>
> > **NOTE:**
> >
> > ``` r
> > DefaultAssay(PBMC) <- 'RNA'
> > VlnPlot(
> >   object = PBMC,
> >   features = c(
> >     'MS4A1', # B-cell
> >     'CD3D', # T-cell
> >     'NKG7', # NK cell
> >     'TREM1' # Monocytes
> >   )
> > )
> > FeaturePlot(
> >   object = PBMC,
> >   features = c('MS4A1', 'CD3D', 'NKG7', 'TREM1'),
> >   max.cutoff = 'q95'
> > )
> > ```
> >
> > Based on these markers, we can annotate clusters. **Note that the cluster numbers obtained in your session may differ from the ones below**: adapt the mapping to your own results!
> >
> > ``` r
> > PBMC$annotation <- dplyr::case_when(
> >   Idents(PBMC) %in% c(0, 5, 9, 12) ~ "Monocytes",
> >   Idents(PBMC) %in% c(1, 2, 3, 4, 7, 10) ~ "T-cells",
> >   Idents(PBMC) %in% c(6) ~ "NK-cells",
> >   Idents(PBMC) %in% c(8, 11) ~ "B-cells",
> >   Idents(PBMC) %in% c(13, 14) ~ "DC"
> > )
> > DimPlot(PBMC, group.by = 'annotation', label = TRUE) + NoLegend()
> > ```

## 13.3 Find peaks differentially accessible across clusters

Finally, we can perform differential peak accessibility analysis to identify peaks that are differentially accessible between clusters. This is done using the `FindMarkers` function from `Seurat`, run on the `ATAC` assay. By default, `FindMarkers` uses a (fast) Wilcoxon rank sum test.

> **TIP:**
>
> The `Signac` documentation recommends a logistic regression test, using the total number of fragments per cell as a latent variable to account for differences in sequencing depth (`test.use = 'LR', latent.vars = 'nCount_ATAC'`). This is more rigorous, but much slower (\> 30 minutes for this comparison): try it after the workshop!

``` r
DefaultAssay(PBMC) <- 'ATAC'
Idents(PBMC) <- PBMC$annotation
da_peaks <- FindMarkers(
  object = PBMC,
  ident.1 = "Monocytes",
  ident.2 = "B-cells",
  min.pct = 0.05
)
head(da_peaks)
VlnPlot(
  object = PBMC,
  features = rownames(da_peaks)[1],
  pt.size = 0.1,
  idents = c("Monocytes","B-cells")
)
CoveragePlot(
  object = PBMC,
  region = rownames(da_peaks)[1],
  extend.upstream = 10000,
  extend.downstream = 10000
)
```

## 13.4 Session info

> **NOTE:**
>
>     ##  ─ Session info ────────────────────────────────────────────────────────────
>     ##   setting  value
>     ##   version  R version 4.6.1 (2026-06-24)
>     ##   os       Ubuntu 24.04.5 LTS
>     ##   system   x86_64, linux-gnu
>     ##   ui       X11
>     ##   language (EN)
>     ##   collate  C
>     ##   ctype    en_US.UTF-8
>     ##   tz       Etc/UTC
>     ##   date     2026-11-06
>     ##   pandoc   3.13 @ /usr/bin/ (via rmarkdown)
>     ##   quarto   1.11.5 @ /usr/local/bin/quarto
>     ##  
>     ##  ─ Packages ────────────────────────────────────────────────────────────────
>     ##   package     * version date (UTC) lib source
>     ##   cli           3.6.6   2026-04-09 [2] RSPM (R 4.6.0)
>     ##   digest        0.6.39  2025-11-19 [2] RSPM (R 4.6.0)
>     ##   evaluate      1.0.5   2025-08-27 [2] RSPM (R 4.6.0)
>     ##   fastmap       1.2.0   2024-05-15 [2] RSPM (R 4.6.0)
>     ##   htmltools     0.5.9   2025-12-04 [2] RSPM (R 4.6.0)
>     ##   htmlwidgets   1.6.4   2023-12-06 [2] RSPM (R 4.6.0)
>     ##   jsonlite      2.0.0   2025-03-27 [2] RSPM (R 4.6.0)
>     ##   knitr         1.52    2026-09-06 [2] RSPM (R 4.6.0)
>     ##   otel          0.2.0   2025-08-29 [2] RSPM (R 4.6.0)
>     ##   rlang         1.3.0   2026-07-05 [2] RSPM (R 4.6.0)
>     ##   rmarkdown     2.32    2026-09-01 [2] RSPM (R 4.6.0)
>     ##   sessioninfo   1.2.4   2026-06-04 [2] RSPM (R 4.6.0)
>     ##   xfun          0.62    2026-10-16 [2] RSPM (R 4.6.0)
>     ##  
>     ##   [1] /tmp/Rtmpz8OMVQ/Rinstb1b9be642
>     ##   [2] /usr/local/lib/R/site-library
>     ##   [3] /usr/local/lib/R/library
>     ##  
>     ##  ───────────────────────────────────────────────────────────────────────────
