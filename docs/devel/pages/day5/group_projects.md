# 16  Group projects

> **NOTE:**
>
> - 14:00 - 18:00: project work, in groups
> - 18:00 - 20:00: group presentations (10’ each group, max 5 slides: what / why / where / when / how, conclusions)

> **TIP:**
>
> - Analyse a public single-cell RNA-seq dataset by yourself, from A to Z
> - Make (and justify!) your own analysis choices
> - Present and discuss your results with the other participants

## 16.1 Choose a dataset

Each group picks one dataset. It can be a dataset from your own research, a dataset from a paper you like, or one of the many datasets distributed by the `scRNAseq` package:

``` r
library(scRNAseq)
datasets <- surveyDatasets()
datasets
searchDatasets("pancreas")
sce <- fetchDataset("baron-pancreas-2016", "2023-12-14", path = "human")
```

The [CZ CELLxGENE Discover](https://cellxgene.cziscience.com/datasets) and [Human Cell Atlas](https://data.humancellatlas.org/explore/projects) portals, as well as [GEO](https://www.ncbi.nlm.nih.gov/geo/), are other good places to find processed count matrices.

## 16.2 Analyse it

Here is a checklist of the analysis steps we covered during the week. Not all of them will be relevant to your dataset: decide which ones are, and why.

1.  Import the data into a `SingleCellExperiment` object ([Lab 3](../../pages/day2/Lab3_Rbioc.md), [Lab 4](../../pages/day2/Lab4_data_wrangling_scRNAseq.md))
2.  Quality control and filtering of cells and genes ([Lab 4](../../pages/day2/Lab4_data_wrangling_scRNAseq.md))
3.  Normalization, feature selection and dimensionality reduction ([Lab 5](../../pages/day3/Lab5_clustering.md))
4.  Batch correction, if required ([Lab 6](../../pages/day3/Lab6_batch_correction.md))
5.  Clustering, marker genes and cell annotation ([Lab 5](../../pages/day3/Lab5_clustering.md))
6.  Differential expression analysis between conditions, with pseudobulk ([Lab 6](../../pages/day3/Lab6_batch_correction.md))
7.  Trajectory inference and pseudotime, if relevant ([Lab 8](../../pages/day4/Lab8_pseudotime.md))

## 16.3 Present it

Prepare a short presentation (10 minutes, max 5 slides), covering:

- **What**: what is the biological question, and what is the dataset?
- **Why**: why is single-cell RNA-seq relevant to answer this question?
- **Where / When / How**: which analysis steps did you perform, and which choices did you make?
- **Conclusions**: what did you learn, and what would you do next?
