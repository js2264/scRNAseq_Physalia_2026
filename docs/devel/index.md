# Single-cell RNA-seq analysis with R/Bioconductor

Physalia workshop, 2026 edition

Authors

Jacques Serizay

Orr Ashenberg

Fabrício Almeida-Silva

![](assets/cover.png "Single-cell RNA-seq analysis with R/Bioconductor")

**Package:** scRNAseqPhysalia2026\
**Authors:** Jacques Serizay \[aut, cre\], Orr Ashenberg \[aut\], Fabricio Almeida-Silva \[aut\]\
**Compiled:** 2026-11-06\
**Package version:** 1.0.0\
**R version:** **R version 4.6.1 (2026-06-24)**\
**BioC version:** **3.24**\
**License:** MIT + file LICENSE\

# Welcome

This is the landing page for the **“Single-cell RNA-seq analysis with R/Bioconductor”** workshop organised by [Physalia](https://www.physalia-courses.org/), **2026 edition**.

It gathers all the lectures and hands-on labs of the course. Every lab can be run in the cloud, in your web browser, with nothing to install: see [Getting started](pages/getting-started.md).

## Program

| **Day** | **Time** | **Session** | **Instructor** |
|----|----|----|----|
| **Day 1** | 14:00 - 15:00 | [Lecture 1: Introduction to scRNAseq analysis](pages/day1/Lecture1_introduction.md) | Jacques |
|  | 15:10 - 16:10 | [Lecture 2: From sequencing reads to expression matrices](pages/day1/Lecture2_Bcl2matrix.md) | Jacques |
|  | 16:10 - 17:10 | **Break** |  |
|  | 17:10 - 18:30 | [Lab 1: Familiarizing yourself with the workshop environment](pages/day1/Lab1_Intro.md) | Jacques |
|  | 18:40 - 20:00 | [Lab 2: From .bcl to count matrix](pages/day1/Lab2_processingreads.md) | Fabricio |
| **Day 2** | 14:00 - 15:00 | [Lecture 3: Quality control for scRNAseq data](pages/day2/Lecture3_qualitycontrol.md) | Jacques |
|  | 15:10 - 16:30 | [Lab 3: Introduction to R/Bioconductor](pages/day2/Lab3_Rbioc.md) | Fabricio |
|  | 16:30 - 17:30 | **Break** |  |
|  | 17:30 - 19:00 | [Lab 4: scRNAseq data wrangling](pages/day2/Lab4_data_wrangling_scRNAseq.md) | Fabricio |
|  | 19:10 - 20:00 | Flash talks | All |
| **Day 3** | 14:00 - 15:30 | [Lecture 4: Identifying cell populations](pages/day3/Lecture4_clustering.md) | Jacques |
|  | 15:40 - 17:00 | [Lab 5: Dimension reduction, clustering, and annotations](pages/day3/Lab5_clustering.md) | Fabricio |
|  | 17:00 - 18:00 | **Break** |  |
|  | 18:00 - 19:00 | [Lecture 5: Data integration and batch effect correction](pages/day3/Lecture5_batchcorrection.md) | Jacques |
|  | 19:10 - 20:00 | [Lab 6: Batch correction](pages/day3/Lab6_batch_correction.md) | Jacques |
| **Day 4** | 14:00 - 15:00 | [Lecture 6: Advances in single-cell genomics - The epigenome](pages/day4/Lecture6_ATAC.md) | Jacques |
|  | 15:10 - 16:20 | [Lab 7: scATACseq analysis workflow](pages/day4/Lab7_atac-seq.md) | Jacques |
|  | 16:20 - 17:20 | **Break** |  |
|  | 17:20 - 18:20 | [Lecture 7: Trajectories and pseudotimes](pages/day4/Lecture7_pseudotime.md) | Jacques |
|  | 18:30 - 20:00 | [Lab 8: Pseudotime analyses](pages/day4/Lab8_pseudotime.md) | Fabricio |
| **Day 5** | 14:00 - 18:00 | [Group projects: project work](pages/day5/group_projects.md) | Jacques |
|  | 18:00 - 20:00 | Wrap-Up: Finalize Projects and Discussion |  |

More details about the program can be found [here](pages/program.md). All the times are in Paris time (CET).

## What

This course will introduce biologists and bioinformaticians to the field of single-cell RNA sequencing. We will cover a range of software and analysis workflows that extend over the spectrum from the best practices in the filtering scRNA-seq data to the downstream analysis of cell clusters and temporal ordering. This course will help the attendees gain accurate insights in pre-processing, analysis and interpretation of scRNA-seq data.

We will start by introducing general concepts about single-cell RNA-sequencing. From there, we will then continue to describe the main analysis steps to go from raw sequencing data to processed and usable data. We will present classical analysis workflows, their output and the possible paths to investigate downstream of this.

Throughout the workshop, `bash` tools and `R/Bioconductor` packages will be used to analyse datasets and learn new approaches.

## When

**Dates of the 2026 edition: to be announced.** Classes run from 14:00 to 20:00 (Paris time), over five days.

## Where

This course will be held online. All the hands-on labs run in the cloud, on the [Bioconductor Orchestra platform](https://orchestraplatform.org/): each participant gets their own `RStudio` session, in their web browser, started from the Docker image built from this book.

## How

The course is structured in modules over five days. Each day will include formal lectures covering the key concepts required to understand scRNA-seq analysis. The remainder of each day will consist in practical hands-on sessions focusing on analysis of scRNA-seq data. These sessions will involve a combination of both mirroring exercises with the instructor to demonstrate a skill, as well as applying these skills on your own to complete individual exercises. During and after each exercise, interpretation of results will be discussed as a group.

## Who

The course will be mostly beneficial to those who have, or will shortly have, scRNA-seq data ready to analyse.

The material is suitable both for experimentalists who want to learn more about data-analysis as well as computational biologists who want to learn about scRNA-seq methods.

Examples demonstrated in this course can be applied to any experimental protocol or biological system.

The requirements for this course are:

1.  Working knowledge of Unix / command line interface (managing files, running programs, reading manuals!). Basic bash commands (`cd`, `ls`, `...`) and CLI usage will *not* be covered in this course. We advice attendees to not register if they lack fundamental experience in CLI.
2.  Programming experience in `R` (writing a function, basic I/O operations, variable types, using packages). `Bioconductor` experience is a plus.
3.  Familiarity with next-generation sequencing data and its analyses (using alignment and quantification tools for bulk sequencing data)

See the [Prerequisites](pages/prerequisites.md) page for refresher material.

## Why

At the end of this course, you should be able to:

- Understand the pros/cons of different single-cell RNA-seq methods
- Process and QC of scRNA-seq data
- Normalize scRNA-seq data
- Correct for batch effects
- Visualise the data and applying dimensionality reduction
- Perform cell clustering and annotation
- Perform differential gene expression analysis
- Infer cell trajectory and pseudotime, and perform temporal differential expression

Throughout the course, we will also have a focus on reproducible research, documented content and interactive reports.

## Instructors

[Jacques Serizay](mailto:jacques.serizay@pasteur.fr)

[Orr Ashenberg](mailto:orr@broadinstitute.org)

[Fabrício Almeida-Silva](mailto:fabricio_almeidasilva@hotmail.com)

# Docker image

This book is a [`BiocBook`](https://github.com/js2264/BiocBook): every version of the book comes with a `Docker` image providing `R`, `Bioconductor`, every package and command-line tool used in the labs, and the datasets they use. This is the image the workshop sessions run on Orchestra. It is available here:

👉 [ghcr.io/js2264/scrnaseqphysalia2026](https://github.com/js2264/scRNAseq_Physalia_2026/pkgs/container/scrnaseqphysalia2026) 🐳

> **TIP:**
>
> You can get access to all the packages used in this book in \< 1 minute, using this command in a terminal:
>
> ``` sh
> docker run -it ghcr.io/js2264/scrnaseqphysalia2026:devel R
> ```

# RStudio Server

An RStudio Server instance can be initiated from the `Docker` image as follows:

``` sh
docker run \
    --volume <local_folder>:<destination_folder> \
    -e PASSWORD=physalia \
    -p 8787:8787 \
    ghcr.io/js2264/scrnaseqphysalia2026:devel
```

The initiated RStudio Server instance will be available at <http://localhost:8787> (user: `rstudio`, password: `physalia`). See [Getting started](pages/getting-started.md) for more details.

# Session info

> **NOTE:**
>
> ``` r
> sessioninfo::session_info(
>     installed.packages()[,"Package"],
>     include_base = TRUE
> )
> ##  ─ Session info ────────────────────────────────────────────────────────────
> ##   setting  value
> ##   version  R version 4.6.1 (2026-06-24)
> ##   os       Ubuntu 24.04.5 LTS
> ##   system   x86_64, linux-gnu
> ##   ui       X11
> ##   language (EN)
> ##   collate  C
> ##   ctype    en_US.UTF-8
> ##   tz       Etc/UTC
> ##   date     2026-11-06
> ##   pandoc   3.13 @ /usr/bin/ (via rmarkdown)
> ##   quarto   1.11.5 @ /usr/local/bin/quarto
> ##  
> ##  ─ Packages ────────────────────────────────────────────────────────────────
> ##   package              * version   date (UTC) lib source
> ##   abind                  1.4-8     2024-09-12 [2] RSPM (R 4.6.0)
> ##   airway                 1.33.2    2026-07-21 [2] Bioconductor 3.24 (R 4.6.1)
> ##   alabaster.base         1.13.5    2026-11-04 [2] Bioconductor 3.24 (R 4.6.1)
> ##   alabaster.matrix       1.13.4    2026-09-23 [2] Bioconductor 3.24 (R 4.6.1)
> ##   alabaster.ranges       1.13.1    2026-09-06 [2] Bioconductor 3.24 (R 4.6.1)
> ##   alabaster.sce          1.13.1    2026-09-06 [2] Bioconductor 3.24 (R 4.6.1)
> ##   alabaster.schemas      1.13.0    2026-04-28 [2] Bioconductor 3.24 (R 4.6.1)
> ##   alabaster.se           1.13.1    2026-09-06 [2] Bioconductor 3.24 (R 4.6.1)
> ##   AnnotationDbi          1.75.2    2026-07-21 [2] Bioconductor 3.24 (R 4.6.1)
> ##   AnnotationFilter       1.37.0    2026-04-28 [2] Bioconductor 3.24 (R 4.6.1)
> ##   AnnotationHub          4.3.3     2026-10-23 [2] Bioconductor 3.24 (R 4.6.1)
> ##   askpass                1.2.1     2024-10-04 [2] RSPM (R 4.6.0)
> ##   assorthead             1.7.11    2026-10-08 [2] Bioconductor 3.24 (R 4.6.1)
> ##   backports              1.5.1     2026-04-03 [2] RSPM (R 4.6.0)
> ##   base                 * 4.6.1     2026-10-29 [3] local
> ##   base64enc              0.1-6     2026-02-02 [2] RSPM (R 4.6.0)
> ##   basilisk               1.25.0    2026-04-29 [2] Bioconductor 3.24 (R 4.6.1)
> ##   batchelor              1.29.0    2026-04-28 [2] Bioconductor 3.24 (R 4.6.1)
> ##   beachmat               2.29.3    2026-09-23 [2] Bioconductor 3.24 (R 4.6.1)
> ##   beeswarm               0.4.0     2021-06-01 [2] RSPM (R 4.6.0)
> ##   BH                     1.90.0-1  2025-12-14 [2] RSPM (R 4.6.0)
> ##   Biobase                2.73.2    2026-07-29 [2] Bioconductor 3.24 (R 4.6.1)
> ##   BiocBaseUtils          1.15.1    2026-05-10 [2] Bioconductor 3.24 (R 4.6.1)
> ##   BiocBook               1.11.7    2026-11-04 [2] Bioconductor 3.24 (R 4.6.1)
> ##   BiocFileCache          3.3.2     2026-11-06 [2] xgit (https://git.bioconductor.org/packages/BiocFileCache@709e8bb)
> ##   BiocGenerics           0.59.12   2026-08-11 [2] Bioconductor 3.24 (R 4.6.1)
> ##   BiocIO                 1.23.3    2026-04-29 [2] Bioconductor 3.24 (R 4.6.1)
> ##   biocmake               1.5.1     2026-08-21 [2] Bioconductor 3.24 (R 4.6.1)
> ##   BiocManager            1.30.27   2025-11-14 [2] CRAN (R 4.6.1)
> ##   BiocNeighbors          2.7.3     2026-08-21 [2] Bioconductor 3.24 (R 4.6.1)
> ##   BiocParallel           1.47.0    2026-04-28 [2] Bioconductor 3.24 (R 4.6.1)
> ##   BiocSingular           1.29.2    2026-10-08 [2] Bioconductor 3.24 (R 4.6.1)
> ##   BiocStyle              2.41.0    2026-04-28 [2] Bioconductor 3.24 (R 4.6.1)
> ##   BiocVersion            3.24.0    2026-04-28 [2] Bioconductor 3.24 (R 4.6.1)
> ##   Biostrings             2.81.9    2026-09-06 [2] Bioconductor 3.24 (R 4.6.1)
> ##   biovizBase             1.61.0    2026-04-28 [2] Bioconductor 3.24 (R 4.6.1)
> ##   bit                    4.6.0     2025-03-06 [2] RSPM (R 4.6.0)
> ##   bit64                  4.8.6     2026-09-01 [2] RSPM (R 4.6.0)
> ##   bitops                 1.1-0     2026-07-30 [2] RSPM (R 4.6.0)
> ##   blob                   1.3.0     2026-01-14 [2] RSPM (R 4.6.0)
> ##   bluster                1.23.1    2026-08-25 [2] Bioconductor 3.24 (R 4.6.1)
> ##   bookdown               0.48      2026-08-28 [2] RSPM (R 4.6.0)
> ##   boot                   1.3-32    2025-08-29 [3] CRAN (R 4.6.1)
> ##   brew                   1.0-10    2023-12-16 [2] RSPM (R 4.6.0)
> ##   brio                   1.1.5     2024-04-24 [2] RSPM (R 4.6.0)
> ##   broom                  1.0.13    2026-05-14 [2] RSPM (R 4.6.0)
> ##   BSgenome               1.81.1    2026-07-29 [2] Bioconductor 3.24 (R 4.6.1)
> ##   bslib                  0.12.0    2026-08-04 [2] RSPM (R 4.6.0)
> ##   cachem                 1.1.0     2024-05-16 [2] RSPM (R 4.6.0)
> ##   callr                  3.8.0     2026-06-05 [2] RSPM (R 4.6.0)
> ##   caTools                1.18.4    2026-07-20 [2] RSPM (R 4.6.0)
> ##   celldex                1.23.0    2026-04-30 [2] Bioconductor 3.24 (R 4.6.1)
> ##   cellranger             1.1.0     2016-07-27 [2] RSPM (R 4.6.0)
> ##   checkmate              2.3.4     2026-02-03 [2] RSPM (R 4.6.0)
> ##   cigarillo              1.3.1     2026-07-13 [2] Bioconductor 3.24 (R 4.6.1)
> ##   class                  7.3-24    2026-08-03 [2] RSPM (R 4.6.1)
> ##   cli                    3.6.6     2026-04-09 [2] RSPM (R 4.6.0)
> ##   clipr                  0.8.1     2026-05-25 [2] RSPM (R 4.6.0)
> ##   cluster                2.1.8.3   2026-07-30 [2] RSPM (R 4.6.1)
> ##   clustree               0.5.1     2023-11-05 [2] RSPM (R 4.6.0)
> ##   codetools              0.2-20    2024-03-31 [3] CRAN (R 4.6.1)
> ##   colorspace             2.1-4     2026-10-08 [2] RSPM (R 4.6.0)
> ##   commonmark             2.0.0     2025-07-07 [2] RSPM (R 4.6.0)
> ##   compiler               4.6.1     2026-10-29 [3] local
> ##   conflicted             1.2.0     2023-02-01 [2] RSPM (R 4.6.0)
> ##   cowplot                1.2.0     2025-07-07 [2] RSPM (R 4.6.0)
> ##   cpp11                  0.5.5     2026-05-06 [2] RSPM (R 4.6.0)
> ##   crayon                 1.5.3     2024-06-20 [2] RSPM (R 4.6.0)
> ##   credentials            2.0.3     2025-09-12 [2] RSPM (R 4.6.0)
> ##   crosstalk              1.2.2     2025-08-26 [2] RSPM (R 4.6.0)
> ##   curl                   8.1.0     2026-10-14 [2] RSPM (R 4.6.0)
> ##   data.table             1.18.6.1  2026-08-24 [2] RSPM (R 4.6.0)
> ##   datasets             * 4.6.1     2026-10-29 [3] local
> ##   DBI                    1.3.0     2026-02-25 [2] RSPM (R 4.6.0)
> ##   dbplyr                 2.6.0     2026-06-17 [2] RSPM (R 4.6.0)
> ##   DelayedArray           0.39.8    2026-09-30 [2] Bioconductor 3.24 (R 4.6.1)
> ##   DelayedMatrixStats     1.35.0    2026-04-28 [2] Bioconductor 3.24 (R 4.6.1)
> ##   deldir                 2.0-4     2024-02-28 [2] RSPM (R 4.6.0)
> ##   desc                   1.4.3     2023-12-10 [2] RSPM (R 4.6.0)
> ##   DESeq2                 1.53.7    2026-11-04 [2] Bioconductor 3.24 (R 4.6.1)
> ##   devtools               2.5.2     2026-04-30 [2] RSPM (R 4.6.0)
> ##   dichromat              2.0-1     2026-07-22 [2] RSPM (R 4.6.0)
> ##   diffobj                0.3.9     2026-09-11 [2] RSPM (R 4.6.0)
> ##   digest                 0.6.39    2025-11-19 [2] RSPM (R 4.6.0)
> ##   dir.expiry             1.21.0    2026-04-28 [2] Bioconductor 3.24 (R 4.6.1)
> ##   docopt                 0.7.2     2025-03-25 [2] RSPM (R 4.6.1)
> ##   dotCall64              1.2       2024-10-04 [2] RSPM (R 4.6.0)
> ##   downlit                0.4.5     2025-11-14 [2] RSPM (R 4.6.0)
> ##   dplyr                  1.2.1     2026-04-03 [2] RSPM (R 4.6.0)
> ##   dqrng                  0.4.1     2024-05-28 [2] RSPM (R 4.6.0)
> ##   DropletUtils           1.33.2    2026-09-24 [2] Bioconductor 3.24 (R 4.6.1)
> ##   dtplyr                 1.3.3     2026-02-11 [2] RSPM (R 4.6.0)
> ##   easy.utils             0.1.0     2025-02-19 [2] RSPM (R 4.6.0)
> ##   edgeR                  4.99.12   2026-11-04 [2] Bioconductor 3.24 (R 4.6.1)
> ##   ellipsis               0.3.3     2026-04-04 [2] RSPM (R 4.6.0)
> ##   EnsDb.Hsapiens.v75     2.99.0    2026-07-20 [2] Bioconductor
> ##   ensembldb              2.37.3    2026-06-07 [2] Bioconductor 3.24 (R 4.6.1)
> ##   evaluate               1.0.5     2025-08-27 [2] RSPM (R 4.6.0)
> ##   ExperimentHub          3.3.2     2026-08-12 [2] Bioconductor 3.24 (R 4.6.1)
> ##   fansi                  1.0.7     2025-11-19 [2] RSPM (R 4.6.0)
> ##   farver                 2.1.2     2024-05-13 [2] RSPM (R 4.6.0)
> ##   fastDummies            1.7.6     2026-04-22 [2] RSPM (R 4.6.0)
> ##   fastmap                1.2.0     2024-05-15 [2] RSPM (R 4.6.0)
> ##   fastmatch              1.1-8     2026-01-17 [2] RSPM (R 4.6.0)
> ##   filelock               1.0.3     2023-12-11 [2] RSPM (R 4.6.0)
> ##   fitdistrplus           1.2-6     2026-01-24 [2] RSPM (R 4.6.0)
> ##   FNN                    1.1.4.1   2024-09-22 [2] RSPM (R 4.6.0)
> ##   fontawesome            0.5.3     2024-11-16 [2] RSPM (R 4.6.0)
> ##   forcats                1.0.1     2025-09-25 [2] RSPM (R 4.6.0)
> ##   foreign                0.8-92    2026-10-17 [2] RSPM (R 4.6.1)
> ##   formatR                1.14      2023-01-17 [2] RSPM (R 4.6.0)
> ##   Formula                1.2-6     2026-08-03 [2] RSPM (R 4.6.0)
> ##   fs                     2.1.0     2026-04-18 [2] RSPM (R 4.6.0)
> ##   futile.logger          1.4.9     2025-12-29 [2] RSPM (R 4.6.0)
> ##   futile.options         1.0.1     2018-04-20 [2] RSPM (R 4.6.0)
> ##   future                 1.80.0    2026-10-16 [2] RSPM (R 4.6.0)
> ##   future.apply           1.21.0    2026-10-16 [2] RSPM (R 4.6.0)
> ##   gargle                 1.6.1     2026-01-29 [2] RSPM (R 4.6.0)
> ##   generics               0.1.4     2025-05-09 [2] RSPM (R 4.6.0)
> ##   GenomeInfoDb           1.49.1    2026-05-24 [2] Bioconductor 3.24 (R 4.6.1)
> ##   GenomicAlignments      1.49.2    2026-09-02 [2] Bioconductor 3.24 (R 4.6.1)
> ##   GenomicFeatures        1.65.0    2026-04-28 [2] Bioconductor 3.24 (R 4.6.1)
> ##   GenomicRanges          1.65.4    2026-09-02 [2] Bioconductor 3.24 (R 4.6.1)
> ##   gert                   2.5.0     2026-10-14 [2] RSPM (R 4.6.0)
> ##   ggbeeswarm             0.7.3     2025-11-29 [2] RSPM (R 4.6.0)
> ##   ggforce                0.5.0     2025-06-18 [2] RSPM (R 4.6.0)
> ##   ggplot2                4.0.3     2026-04-22 [2] RSPM (R 4.6.0)
> ##   ggraph                 2.2.2     2025-08-24 [2] RSPM (R 4.6.0)
> ##   ggrepel                0.9.8     2026-03-17 [2] RSPM (R 4.6.0)
> ##   ggridges               0.5.7     2025-08-27 [2] RSPM (R 4.6.0)
> ##   gh                     1.6.1     2026-07-20 [2] RSPM (R 4.6.0)
> ##   gitcreds               0.1.2     2022-09-08 [2] RSPM (R 4.6.0)
> ##   globals                0.20.0    2026-10-11 [2] RSPM (R 4.6.0)
> ##   glue                   1.8.1     2026-04-17 [2] RSPM (R 4.6.0)
> ##   goftest                1.2-3     2021-10-07 [2] RSPM (R 4.6.0)
> ##   googledrive            2.1.2     2025-09-10 [2] RSPM (R 4.6.0)
> ##   googlesheets4          1.1.2     2025-09-03 [2] RSPM (R 4.6.0)
> ##   gplots                 3.3.0     2025-11-30 [2] RSPM (R 4.6.0)
> ##   graphics             * 4.6.1     2026-10-29 [3] local
> ##   graphlayouts           1.2.5     2026-07-22 [2] RSPM (R 4.6.0)
> ##   grDevices            * 4.6.1     2026-10-29 [3] local
> ##   grid                   4.6.1     2026-10-29 [3] local
> ##   gridExtra              2.3.1     2026-06-25 [2] RSPM (R 4.6.0)
> ##   gtable                 0.3.6     2024-10-25 [2] RSPM (R 4.6.0)
> ##   gtools                 3.9.5     2023-11-20 [2] RSPM (R 4.6.0)
> ##   gypsum                 1.9.0     2026-04-28 [2] Bioconductor 3.24 (R 4.6.1)
> ##   h5mread                1.5.3     2026-09-02 [2] Bioconductor 3.24 (R 4.6.1)
> ##   harmony                2.0.5     2026-06-12 [2] RSPM (R 4.6.0)
> ##   haven                  2.5.5     2025-05-30 [2] RSPM (R 4.6.0)
> ##   HDF5Array              1.41.3    2026-08-21 [2] Bioconductor 3.24 (R 4.6.1)
> ##   hdf5r                  1.3.16    2026-09-11 [2] RSPM (R 4.6.0)
> ##   hdf5r.Extra            0.1.0     2024-10-18 [2] RSPM (R 4.6.0)
> ##   here                   1.0.2     2025-09-15 [2] RSPM (R 4.6.0)
> ##   HighFive               3.3.0     2026-01-08 [2] RSPM (R 4.6.0)
> ##   highr                  0.12      2026-03-06 [2] RSPM (R 4.6.0)
> ##   Hmisc                  5.3-0     2026-09-06 [2] RSPM (R 4.6.0)
> ##   hms                    1.1.4     2025-10-17 [2] RSPM (R 4.6.0)
> ##   htmlTable              2.5.0     2026-04-22 [2] RSPM (R 4.6.0)
> ##   htmltools              0.5.9     2025-12-04 [2] RSPM (R 4.6.0)
> ##   htmlwidgets            1.6.4     2023-12-06 [2] RSPM (R 4.6.0)
> ##   httpuv                 1.6.17    2026-03-18 [2] RSPM (R 4.6.0)
> ##   httr                   1.4.9     2026-09-01 [2] RSPM (R 4.6.0)
> ##   httr2                  1.3.0     2026-07-13 [2] RSPM (R 4.6.0)
> ##   ica                    1.0-3     2022-07-08 [2] RSPM (R 4.6.0)
> ##   ids                    1.0.1     2017-05-31 [2] RSPM (R 4.6.0)
> ##   igraph                 2.3.4     2026-09-30 [2] RSPM (R 4.6.0)
> ##   ini                    0.3.1     2018-05-20 [2] RSPM (R 4.6.0)
> ##   IRanges                2.47.5    2026-08-27 [2] Bioconductor 3.24 (R 4.6.1)
> ##   irlba                  2.4.1     2026-10-05 [2] RSPM (R 4.6.0)
> ##   isoband                0.3.0     2025-12-07 [2] RSPM (R 4.6.0)
> ##   jquerylib              0.1.4     2021-04-26 [2] RSPM (R 4.6.0)
> ##   jsonlite               2.0.0     2025-03-27 [2] RSPM (R 4.6.0)
> ##   jsonvalidate           1.5.0     2025-02-07 [2] RSPM (R 4.6.0)
> ##   KEGGREST               1.53.6    2026-07-23 [2] Bioconductor 3.24 (R 4.6.1)
> ##   KernSmooth             2.23-27   2026-08-12 [2] RSPM (R 4.6.1)
> ##   knitr                  1.52      2026-09-06 [2] RSPM (R 4.6.0)
> ##   labeling               0.4.3     2023-08-29 [2] RSPM (R 4.6.0)
> ##   lambda.r               1.2.6     2026-10-13 [2] RSPM (R 4.6.0)
> ##   later                  1.4.8     2026-03-05 [2] RSPM (R 4.6.0)
> ##   lattice                0.23-1    2026-08-12 [2] RSPM (R 4.6.1)
> ##   lazyeval               0.2.3     2026-04-04 [2] RSPM (R 4.6.0)
> ##   leidenAlg              1.1.8     2026-05-31 [2] RSPM (R 4.6.0)
> ##   lifecycle              1.0.5     2026-01-08 [2] RSPM (R 4.6.0)
> ##   limma                  3.99.4    2026-11-04 [2] Bioconductor 3.24 (R 4.6.1)
> ##   listenv                1.1.0     2026-10-04 [2] RSPM (R 4.6.0)
> ##   littler                0.3.23    2026-04-12 [2] RSPM (R 4.6.1)
> ##   lmtest                 0.9-40    2022-03-21 [2] RSPM (R 4.6.0)
> ##   locfit                 1.5-9.12  2025-03-05 [2] RSPM (R 4.6.0)
> ##   LoomExperiment         1.31.0    2026-04-29 [2] Bioconductor 3.24 (R 4.6.1)
> ##   lubridate              1.9.5     2026-02-04 [2] RSPM (R 4.6.0)
> ##   magrittr               2.0.5     2026-04-04 [2] RSPM (R 4.6.0)
> ##   MASS                   7.3-66    2026-07-15 [2] RSPM (R 4.6.1)
> ##   Matrix                 1.7-6     2026-07-25 [2] RSPM (R 4.6.1)
> ##   MatrixGenerics         1.25.0    2026-04-28 [2] Bioconductor 3.24 (R 4.6.1)
> ##   matrixStats            1.5.0     2025-01-07 [2] RSPM (R 4.6.0)
> ##   memoise                2.0.1     2021-11-26 [2] RSPM (R 4.6.0)
> ##   metapod                1.21.0    2026-04-28 [2] Bioconductor 3.24 (R 4.6.1)
> ##   methods              * 4.6.1     2026-10-29 [3] local
> ##   mgcv                   1.9-4     2025-11-07 [3] CRAN (R 4.6.1)
> ##   mime                   0.13      2025-03-17 [2] RSPM (R 4.6.0)
> ##   miniUI                 0.1.2     2025-04-17 [2] RSPM (R 4.6.0)
> ##   modelr                 0.1.11    2023-03-22 [2] RSPM (R 4.6.0)
> ##   nlme                   3.1-171   2026-09-01 [2] RSPM (R 4.6.1)
> ##   nnet                   7.3-21    2026-08-03 [2] RSPM (R 4.6.1)
> ##   openssl                2.4.2     2026-06-09 [2] RSPM (R 4.6.0)
> ##   org.Hs.eg.db           3.24.4    2026-11-06 [2] Bioconductor
> ##   otel                   0.2.0     2025-08-29 [2] RSPM (R 4.6.0)
> ##   pak                    0.11.1    2026-07-22 [2] RSPM (R 4.6.0)
> ##   parallel               4.6.1     2026-10-29 [3] local
> ##   parallelly             1.48.0    2026-06-29 [2] RSPM (R 4.6.0)
> ##   patchwork              1.3.2     2025-08-25 [2] RSPM (R 4.6.0)
> ##   pbapply                1.7-5     2026-09-01 [2] RSPM (R 4.6.0)
> ##   pbmcapply              1.5.1     2022-04-28 [2] RSPM (R 4.6.0)
> ##   pheatmap               1.0.13    2025-06-05 [2] RSPM (R 4.6.0)
> ##   pillar                 1.11.1    2025-09-17 [2] RSPM (R 4.6.0)
> ##   pkgbuild               1.4.8     2025-05-26 [2] RSPM (R 4.6.0)
> ##   pkgconfig              2.0.3     2019-09-22 [2] RSPM (R 4.6.0)
> ##   pkgdown                2.2.1     2026-07-07 [2] RSPM (R 4.6.0)
> ##   pkgload                1.5.3     2026-06-15 [2] RSPM (R 4.6.0)
> ##   plotly                 4.12.1    2026-07-22 [2] RSPM (R 4.6.0)
> ##   plyr                   1.8.9     2023-10-02 [2] RSPM (R 4.6.0)
> ##   plyranges              1.33.2    2026-07-30 [2] Bioconductor 3.24 (R 4.6.1)
> ##   png                    0.1-9     2026-03-15 [2] RSPM (R 4.6.0)
> ##   polyclip               1.10-7    2024-07-23 [2] RSPM (R 4.6.0)
> ##   praise                 1.0.0     2015-08-11 [2] RSPM (R 4.6.0)
> ##   preprocessCore         1.75.1    2026-08-31 [2] Bioconductor 3.24 (R 4.6.1)
> ##   prettyunits            1.2.0     2023-09-24 [2] RSPM (R 4.6.0)
> ##   princurve              2.1.6     2021-01-18 [2] RSPM (R 4.6.0)
> ##   pROC                   1.19.1    2026-08-27 [2] RSPM (R 4.6.0)
> ##   processx               3.9.0     2026-04-22 [2] RSPM (R 4.6.0)
> ##   profvis                0.4.0     2024-09-20 [2] RSPM (R 4.6.0)
> ##   progress               1.2.3     2023-12-06 [2] RSPM (R 4.6.0)
> ##   progressr              1.1.0     2026-10-18 [2] RSPM (R 4.6.0)
> ##   promises               1.5.0     2025-11-01 [2] RSPM (R 4.6.0)
> ##   ProtGenerics           1.45.0    2026-04-28 [2] Bioconductor 3.24 (R 4.6.1)
> ##   ps                     1.9.3     2026-04-20 [2] RSPM (R 4.6.0)
> ##   purrr                  1.2.2     2026-04-10 [2] RSPM (R 4.6.0)
> ##   quarto                 1.5.1     2025-09-04 [2] RSPM (R 4.6.0)
> ##   R.methodsS3            1.8.2     2022-06-13 [2] RSPM (R 4.6.0)
> ##   R.oo                   1.27.1    2025-05-02 [2] RSPM (R 4.6.0)
> ##   R.utils                2.13.0    2025-02-24 [2] RSPM (R 4.6.0)
> ##   R6                     2.6.1     2025-02-15 [2] RSPM (R 4.6.0)
> ##   ragg                   1.5.2     2026-03-23 [2] RSPM (R 4.6.0)
> ##   RANN                   2.6.3     2026-08-24 [2] RSPM (R 4.6.0)
> ##   rappdirs               0.3.4     2026-01-17 [2] RSPM (R 4.6.0)
> ##   rcmdcheck              1.4.0     2021-09-27 [2] RSPM (R 4.6.0)
> ##   RColorBrewer           1.1-3     2022-04-03 [2] RSPM (R 4.6.0)
> ##   Rcpp                   1.1.2     2026-07-05 [2] RSPM (R 4.6.0)
> ##   RcppAnnoy              0.0.23    2026-01-12 [2] RSPM (R 4.6.0)
> ##   RcppArmadillo          15.6.2-1  2026-10-17 [2] RSPM (R 4.6.0)
> ##   RcppEigen              0.3.4.0.2 2024-08-24 [2] RSPM (R 4.6.0)
> ##   RcppHNSW               0.7.0     2026-05-26 [2] RSPM (R 4.6.0)
> ##   RcppML                 0.3.7.1   2026-03-09 [2] RSPM (R 4.6.0)
> ##   RcppPlanc              2.0.16    2026-09-09 [2] RSPM (R 4.6.0)
> ##   RcppProgress           0.4.2     2020-02-06 [2] RSPM (R 4.6.0)
> ##   RcppRoll               0.4.0     2026-09-06 [2] RSPM (R 4.6.0)
> ##   RcppThread             2.4.0     2026-06-11 [2] RSPM (R 4.6.0)
> ##   RcppTOML               0.2.3     2025-03-08 [2] RSPM (R 4.6.0)
> ##   RCurl                  1.98-1.20 2026-08-21 [2] RSPM (R 4.6.0)
> ##   rdtools                0.1.0     2026-07-16 [2] RSPM (R 4.6.0)
> ##   readr                  2.2.0     2026-02-19 [2] RSPM (R 4.6.0)
> ##   readxl                 1.6.0     2026-10-17 [2] RSPM (R 4.6.0)
> ##   rematch                2.0.0     2023-08-30 [2] RSPM (R 4.6.0)
> ##   rematch2               2.1.2     2020-05-01 [2] RSPM (R 4.6.0)
> ##   remotes                2.5.0     2024-03-17 [2] RSPM (R 4.6.0)
> ##   renv                   1.3.1     2026-10-06 [2] RSPM (R 4.6.0)
> ##   reprex                 2.1.1     2024-07-06 [2] RSPM (R 4.6.0)
> ##   reshape2               1.4.5     2025-11-12 [2] RSPM (R 4.6.0)
> ##   ResidualMatrix         1.23.0    2026-04-28 [2] Bioconductor 3.24 (R 4.6.1)
> ##   restfulr               0.0.17    2026-06-11 [2] RSPM (R 4.6.0)
> ##   reticulate             1.47.0    2026-09-03 [2] RSPM (R 4.6.0)
> ##   rhdf5                  2.57.20   2026-10-06 [2] Bioconductor 3.24 (R 4.6.1)
> ##   rhdf5filters           1.25.5    2026-10-28 [2] Bioconductor 3.24 (R 4.6.1)
> ##   Rhdf5lib               2.1.0     2026-04-28 [2] Bioconductor 3.24 (R 4.6.1)
> ##   RhpcBLASctl            0.23-42   2023-02-11 [2] RSPM (R 4.6.0)
> ##   Rhtslib                3.9.0     2026-04-28 [2] Bioconductor 3.24 (R 4.6.1)
> ##   Rigraphlib             1.5.1     2026-05-08 [2] Bioconductor 3.24 (R 4.6.1)
> ##   rjson                  0.2.23    2024-09-16 [2] RSPM (R 4.6.0)
> ##   rlang                  1.3.0     2026-07-05 [2] RSPM (R 4.6.0)
> ##   rliger                 2.2.1     2025-08-26 [2] RSPM (R 4.6.0)
> ##   rmarkdown              2.32      2026-09-01 [2] RSPM (R 4.6.0)
> ##   ROCR                   1.0-12    2026-01-23 [2] RSPM (R 4.6.0)
> ##   roxygen2               8.1.1     2026-10-08 [2] RSPM (R 4.6.0)
> ##   rpart                  4.1.27    2026-03-27 [3] CRAN (R 4.6.1)
> ##   rprojroot              2.1.1     2025-08-26 [2] RSPM (R 4.6.0)
> ##   Rsamtools              2.29.0    2026-04-28 [2] Bioconductor 3.24 (R 4.6.1)
> ##   RSpectra               0.16-2    2024-07-18 [2] RSPM (R 4.6.0)
> ##   RSQLite                3.53.3    2026-06-30 [2] RSPM (R 4.6.0)
> ##   rstudioapi             0.19.0    2026-06-11 [2] RSPM (R 4.6.0)
> ##   rsvd                   1.0.6     2026-10-17 [2] RSPM (R 4.6.0)
> ##   rtracklayer            1.73.0    2026-04-28 [2] Bioconductor 3.24 (R 4.6.1)
> ##   Rtsne                  0.17      2023-12-07 [2] RSPM (R 4.6.0)
> ##   rversions              3.0.0     2025-10-09 [2] RSPM (R 4.6.0)
> ##   rvest                  1.0.5     2025-08-29 [2] RSPM (R 4.6.0)
> ##   S4Arrays               1.13.2    2026-09-30 [2] Bioconductor 3.24 (R 4.6.1)
> ##   S4Vectors              0.51.10   2026-09-16 [2] Bioconductor 3.24 (R 4.6.1)
> ##   S7                     0.2.2     2026-04-22 [2] RSPM (R 4.6.0)
> ##   sass                   0.4.10    2025-04-11 [2] RSPM (R 4.6.0)
> ##   ScaledMatrix           1.21.0    2026-04-28 [2] Bioconductor 3.24 (R 4.6.1)
> ##   scales                 1.4.0     2025-04-24 [2] RSPM (R 4.6.0)
> ##   scater                 1.41.3    2026-10-25 [2] Bioconductor 3.24 (R 4.6.1)
> ##   scattermore            1.2       2023-06-12 [2] RSPM (R 4.6.0)
> ##   sccore                 1.0.7     2026-04-06 [2] RSPM (R 4.6.0)
> ##   scDblFinder            1.27.12   2026-10-14 [2] Bioconductor 3.24 (R 4.6.1)
> ##   scran                  1.41.1    2026-05-20 [2] Bioconductor 3.24 (R 4.6.1)
> ##   scrapper               1.7.11    2026-11-04 [2] Bioconductor 3.24 (R 4.6.1)
> ##   scRNAseq               2.27.0    2026-04-30 [2] Bioconductor 3.24 (R 4.6.1)
> ##   scRNAseqPhysalia2026   1.0.0     2026-11-06 [1] local
> ##   sctransform            0.4.3     2026-01-10 [2] RSPM (R 4.6.0)
> ##   scuttle                1.23.2    2026-08-27 [2] Bioconductor 3.24 (R 4.6.1)
> ##   selectr                0.8-0     2026-09-27 [2] RSPM (R 4.6.0)
> ##   Seqinfo                1.3.2     2026-08-27 [2] Bioconductor 3.24 (R 4.6.1)
> ##   sessioninfo            1.2.4     2026-06-04 [2] RSPM (R 4.6.0)
> ##   Seurat                 5.6.0     2026-10-03 [2] RSPM (R 4.6.0)
> ##   SeuratObject           5.4.0     2026-04-11 [2] RSPM (R 4.6.0)
> ##   shiny                  1.14.0    2026-06-21 [2] RSPM (R 4.6.0)
> ##   Signac                 1.17.1    2026-04-06 [2] RSPM (R 4.6.0)
> ##   SingleCellExperiment   1.35.2    2026-07-16 [2] Bioconductor 3.24 (R 4.6.1)
> ##   SingleR                2.15.4    2026-09-10 [2] Bioconductor 3.24 (R 4.6.1)
> ##   sitmo                  2.0.2     2021-10-13 [2] RSPM (R 4.6.0)
> ##   slingshot              2.21.0    2026-04-28 [2] Bioconductor 3.24 (R 4.6.1)
> ##   snow                   0.4-4     2021-10-27 [2] RSPM (R 4.6.0)
> ##   sourcetools            0.1.7-2   2026-03-28 [2] RSPM (R 4.6.0)
> ##   sp                     2.2-4     2026-10-11 [2] RSPM (R 4.6.0)
> ##   spam                   2.11-4    2026-05-29 [2] RSPM (R 4.6.0)
> ##   SparseArray            1.13.4    2026-09-30 [2] Bioconductor 3.24 (R 4.6.1)
> ##   sparseMatrixStats      1.25.0    2026-04-28 [2] Bioconductor 3.24 (R 4.6.1)
> ##   spatial                7.3-19    2026-08-03 [2] RSPM (R 4.6.1)
> ##   spatstat.data          3.1-9     2025-10-18 [2] RSPM (R 4.6.0)
> ##   spatstat.explore       3.8-3     2026-09-23 [2] RSPM (R 4.6.0)
> ##   spatstat.geom          3.8-3     2026-09-18 [2] RSPM (R 4.6.0)
> ##   spatstat.random        3.5-2     2026-09-22 [2] RSPM (R 4.6.0)
> ##   spatstat.sparse        3.2-0     2026-05-21 [2] RSPM (R 4.6.0)
> ##   spatstat.univar        3.2-0     2026-05-18 [2] RSPM (R 4.6.0)
> ##   spatstat.utils         3.2-5     2026-09-11 [2] RSPM (R 4.6.0)
> ##   splines                4.6.1     2026-10-29 [3] local
> ##   statmod                1.5.2     2026-05-17 [2] RSPM (R 4.6.0)
> ##   stats                * 4.6.1     2026-10-29 [3] local
> ##   stats4                 4.6.1     2026-10-29 [3] local
> ##   stringi                1.8.9     2026-08-04 [2] RSPM (R 4.6.0)
> ##   stringr                1.6.0     2025-11-04 [2] RSPM (R 4.6.0)
> ##   SummarizedExperiment   1.43.0    2026-04-28 [2] Bioconductor 3.24 (R 4.6.1)
> ##   survival               3.8-12    2026-09-09 [2] RSPM (R 4.6.1)
> ##   sys                    3.4.3     2024-10-04 [2] RSPM (R 4.6.0)
> ##   systemfonts            1.3.2     2026-03-05 [2] RSPM (R 4.6.0)
> ##   tcltk                  4.6.1     2026-10-29 [3] local
> ##   tensor                 1.5.1     2025-06-17 [2] RSPM (R 4.6.0)
> ##   TENxPBMCData           1.31.0    2026-04-30 [2] Bioconductor 3.24 (R 4.6.1)
> ##   testthat               3.3.2     2026-01-11 [2] RSPM (R 4.6.0)
> ##   textshaping            1.0.5     2026-03-06 [2] RSPM (R 4.6.0)
> ##   tibble                 3.3.1     2026-01-11 [2] RSPM (R 4.6.0)
> ##   tidygraph              1.3.1     2024-01-30 [2] RSPM (R 4.6.0)
> ##   tidyr                  1.3.2     2025-12-19 [2] RSPM (R 4.6.0)
> ##   tidyselect             1.2.1     2024-03-11 [2] RSPM (R 4.6.0)
> ##   tidyverse              2.0.0     2023-02-22 [2] RSPM (R 4.6.0)
> ##   timechange             0.4.0     2026-01-29 [2] RSPM (R 4.6.0)
> ##   tinytex                0.61      2026-09-17 [2] RSPM (R 4.6.0)
> ##   tools                  4.6.1     2026-10-29 [3] local
> ##   TrajectoryUtils        1.21.0    2026-04-28 [2] Bioconductor 3.24 (R 4.6.1)
> ##   tweenr                 2.0.3     2024-02-26 [2] RSPM (R 4.6.0)
> ##   tzdb                   0.5.0     2025-03-15 [2] RSPM (R 4.6.0)
> ##   UCSC.utils             1.9.0     2026-04-28 [2] Bioconductor 3.24 (R 4.6.1)
> ##   urlchecker             2.0.0     2026-07-08 [2] RSPM (R 4.6.0)
> ##   usethis                3.2.2     2026-09-10 [2] RSPM (R 4.6.0)
> ##   utf8                   1.2.6     2025-06-08 [2] RSPM (R 4.6.0)
> ##   utils                * 4.6.1     2026-10-29 [3] local
> ##   uuid                   1.2-2     2026-01-23 [2] RSPM (R 4.6.0)
> ##   uwot                   0.2.5     2026-08-29 [2] RSPM (R 4.6.0)
> ##   V8                     8.2.1     2026-10-26 [2] RSPM (R 4.6.0)
> ##   VariantAnnotation      1.59.4    2026-07-24 [2] Bioconductor 3.24 (R 4.6.1)
> ##   vctrs                  0.7.3     2026-04-11 [2] RSPM (R 4.6.0)
> ##   velociraptor           1.23.0    2026-04-29 [2] Bioconductor 3.24 (R 4.6.1)
> ##   vipor                  0.4.7     2023-12-18 [2] RSPM (R 4.6.0)
> ##   viridis                0.6.5     2024-01-29 [2] RSPM (R 4.6.0)
> ##   viridisLite            0.4.3     2026-02-04 [2] RSPM (R 4.6.0)
> ##   vroom                  1.7.1     2026-03-31 [2] RSPM (R 4.6.0)
> ##   waldo                  0.6.2     2025-07-11 [2] RSPM (R 4.6.0)
> ##   whisker                0.4.1     2022-12-05 [2] RSPM (R 4.6.0)
> ##   withr                  3.0.3     2026-06-19 [2] RSPM (R 4.6.0)
> ##   xfun                   0.62      2026-10-16 [2] RSPM (R 4.6.0)
> ##   xgboost                3.2.1.1   2026-03-18 [2] RSPM (R 4.6.0)
> ##   XML                    3.99-0.25 2026-09-27 [2] RSPM (R 4.6.0)
> ##   xml2                   1.6.0     2026-06-22 [2] RSPM (R 4.6.0)
> ##   xopen                  1.0.1     2024-04-25 [2] RSPM (R 4.6.0)
> ##   xtable                 1.8-8     2026-02-22 [2] RSPM (R 4.6.0)
> ##   XVector                0.53.0    2026-04-28 [2] Bioconductor 3.24 (R 4.6.1)
> ##   yaml                   2.3.12    2025-12-10 [2] RSPM (R 4.6.0)
> ##   zellkonverter          1.23.0    2026-04-29 [2] Bioconductor 3.24 (R 4.6.1)
> ##   zip                    3.0.2     2026-08-04 [2] RSPM (R 4.6.0)
> ##   zoo                    1.9-1     2026-09-25 [2] RSPM (R 4.6.0)
> ##  
> ##   [1] /tmp/RtmpYDUaaM/Rinstca9466ba
> ##   [2] /usr/local/lib/R/site-library
> ##   [3] /usr/local/lib/R/library
> ##   * ── Packages attached to the search path.
> ##  
> ##  ─ Python configuration ────────────────────────────────────────────────────
> ##   Python is not available
> ##  
> ##  ───────────────────────────────────────────────────────────────────────────
> ```
