<!-- badges: start -->
📦 [Repo](https://github.com/js2264/scRNAseq_Physalia_2026) [![rworkflows](https://img.shields.io/github/actions/workflow/status/js2264/scRNAseq_Physalia_2026/rworkflows.yml?label=Package%20check)](https://github.com/js2264/scRNAseq_Physalia_2026/actions/workflows/rworkflows.yml)   
📖 [Book](https://js2264.github.io/scRNAseq_Physalia_2026/devel) [![deployment](https://img.shields.io/github/actions/workflow/status/js2264/scRNAseq_Physalia_2026/pages/pages-build-deployment?label=Book%20deployment)](https://github.com/js2264/scRNAseq_Physalia_2026/actions/workflows/pages/pages-build-deployment)  
🐳 [Docker](https://github.com/js2264/scRNAseq_Physalia_2026/pkgs/container/scrnaseqphysalia2026) [![biocbook](https://img.shields.io/github/actions/workflow/status/js2264/scRNAseq_Physalia_2026/biocbook.yml?label=Docker%20image)](https://github.com/js2264/scRNAseq_Physalia_2026/actions/workflows/biocbook.yml)  
<!-- badges: end -->

# Single-cell RNA-seq analysis with R/Bioconductor — Physalia 2026

Course material for the 2026 edition of the Physalia workshop
*"Single-cell RNA-seq analysis with R/Bioconductor"*, written as a
[`BiocBook`](https://github.com/js2264/BiocBook): a package-based, versioned
online book with a supporting Docker image.

📖 **Book: <https://js2264.github.io/scRNAseq_Physalia_2026/devel/>**  
🐳 **Docker image: `ghcr.io/js2264/scrnaseqphysalia2026:devel`**

The workshop is run on the [Bioconductor Orchestra platform](https://orchestraplatform.org/),
from the Docker image built from this repository: every participant gets their
own `RStudio` session, in their web browser, with all the packages, tools and
datasets of the labs.

## Content

| Day | Lectures | Labs |
|---|---|---|
| 1 | Introduction to scRNA-seq; from reads to count matrices | Workshop environment, public data (bash); `cellranger` (bash) |
| 2 | Quality control | Introduction to R/Bioconductor; scRNA-seq data wrangling |
| 3 | Clustering; batch correction | Dimensionality reduction, clustering, annotation; batch correction |
| 4 | scATAC-seq; trajectories | scATAC-seq with `Signac`; pseudotime with `slingshot` |
| 5 | — | Group projects |

The lecture slides (PDF) are in `inst/pages/day*/`. The `R` labs are run
(and their outputs shown) when the book is built, with public datasets fetched
from GEO, 10x Genomics, the Human Protein Atlas, `ExperimentHub` and the
`scRNAseq`/`celldex` packages. The labs relying on command-line tools (Labs 1-2,
`cellranger`) and on large scATAC-seq `fragments` files (Lab 7) are displayed
but not run.

The labs follow the current Bioconductor devel (3.24) recommendations: most of
`scuttle` and `scran` is deprecated in favour of `scrapper` (QC, normalization,
HVGs, PCA, clustering, markers, pseudobulk aggregation) and `bluster` (graphs).

## Running the material

**On Orchestra** (during the workshop): open the dashboard link given by the
instructors, sign in with a Google account, launch the workshop and open it.

**With Docker**:

```sh
docker run -e PASSWORD=physalia -p 8787:8787 ghcr.io/js2264/scrnaseqphysalia2026:devel
```

then open <http://localhost:8787> and log in as `rstudio` / `physalia`. The
datasets of the `R` labs are already cached in the image. Give Docker at least
12 GB of RAM (16 GB recommended) to run every lab.

**Locally**: `BiocManager::install("js2264/scRNAseq_Physalia_2026", dependencies = TRUE)`
(R >= 4.6, Bioconductor devel).

## Editing the book

Book pages live in `inst/` (`index.qmd` and `pages/`), and the book settings
(chapters, title, sidebar) in `inst/assets/_book.yml`. Any package used in a
page must be listed in `DESCRIPTION` (`Suggests:`), so that it is installed in
the Docker image. Command-line tools go in `inst/requirements.yml` (conda env).

```r
library(BiocBook)
bb <- BiocBook(".")
check_deps(bb)   # are all packages used in the pages listed in DESCRIPTION?
preview(bb)      # render and preview the book locally
```

Every push to `devel` builds the Docker image, renders the book inside it, then
pushes the image to `ghcr.io/js2264/scrnaseqphysalia2026` and deploys the book to
`gh-pages` (`docs/devel/`), via `.github/workflows/biocbook.yml`.

## Deployment checklist

`BiocBook::init()` normally configures these; this repository was set up by
hand, so they must be done once in the GitHub settings:

1. **Actions permissions**: Settings → Actions → General → Workflow permissions
   → *Read and write permissions*.
2. **GitHub Pages**: after the first successful `biocbook` workflow run (which
   creates the `gh-pages` branch), Settings → Pages → Source: *Deploy from a
   branch*, branch `gh-pages`, folder `/docs`.
3. **⚠️ Make the Docker image public**: after the first push of the image, go to
   <https://github.com/users/js2264/packages/container/scrnaseqphysalia2026/settings>
   → *Danger Zone* → *Change visibility* → **Public**. Orchestra (and anyone
   without a GitHub account) cannot pull a private image.
4. **Orchestra template**: submit the workshop through the
   [Orchestra submission form](https://orchestraplatform.org/docs/contribute/submit/),
   with image `ghcr.io/js2264/scrnaseqphysalia2026:devel`, app port `8787`
   (no env/args needed: Bioconductor images start without login on Orchestra),
   a session duration covering a day of class (e.g. `8h`), and the
   **X-Large** size: Lab 7 (scATAC-seq) peaks at ~12 GB of RAM, Lab 8 at ~8 GB.

## TODO before the workshop

- [ ] Dates of the 2026 edition and list of instructors (`inst/index.qmd`, `DESCRIPTION`, `inst/assets/_book.yml`).
- [ ] Update the lecture slides (`inst/pages/day*/*.pdf`, currently the 2025 ones).
- [ ] Labs 1-2: decide how `cellranger`, `bcl2fastq` and the `cellranger` transcriptome
      reference (~10 GB) are provided on Orchestra. `cellranger` is distributed under a
      10x Genomics license and is not included in the (public) Docker image.
- [ ] Lab 8: provide the `Guo_testis_full-counts.loom` file (spliced/unspliced
      counts) for the RNA velocity section.
- [ ] Remove the temporary `BiocFileCache` fix from the `Dockerfile` once
      BiocFileCache >= 3.3.2 is available from the Bioconductor devel repositories.

## Licence

MIT (see `LICENSE`). Lecture slides © the instructors.
