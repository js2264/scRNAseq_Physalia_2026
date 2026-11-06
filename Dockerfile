# syntax=docker/dockerfile:1
ARG BIOC_VERSION
FROM bioconductor/bioconductor_docker:${BIOC_VERSION}
COPY . /opt/pkg

## Keep R's user cache (the book's `micromamba`, its conda env, and every
## dataset downloaded when rendering the book with `BiocFileCache`,
## `ExperimentHub` or `gypsum`) inside the image, and make RStudio sessions
## use it too (RStudio does not inherit the container's environment
## variables): the datasets of the labs are then readily available to the
## workshop participants
ENV R_USER_CACHE_DIR=/opt/R-cache
RUN echo "R_USER_CACHE_DIR=/opt/R-cache" >> "$(R RHOME)/etc/Renviron.site"

## Optionally pin `quarto`, e.g. to a release that renders `llms.txt` for books
## (quarto >= 1.11). Empty: keep the quarto shipped with the Bioconductor image.
## The download is retried on transient errors (e.g. HTTP 503 from GitHub)
ARG QUARTO_VERSION=""
ARG TARGETARCH
RUN if [ -n "${QUARTO_VERSION}" ]; then \
      curl -fsSL --retry 5 --retry-delay 10 -o /tmp/quarto.tar.gz "https://github.com/quarto-dev/quarto-cli/releases/download/v${QUARTO_VERSION}/quarto-${QUARTO_VERSION}-linux-${TARGETARCH:-amd64}.tar.gz" && \
      mkdir -p "/opt/quarto-${QUARTO_VERSION}" && \
      tar -xzf /tmp/quarto.tar.gz -C "/opt/quarto-${QUARTO_VERSION}" --strip-components=1 && \
      ln -sf "/opt/quarto-${QUARTO_VERSION}/bin/quarto" /usr/local/bin/quarto && \
      rm /tmp/quarto.tar.gz ; \
    fi && \
    quarto --version

## Install book package. GITHUB_PAT arrives as a BuildKit secret, so that
## GitHub-only dependencies do not hit the anonymous API rate limit, and is
## never stored in the image
RUN --mount=type=secret,id=GITHUB_PAT \
    GITHUB_PAT="$(cat /run/secrets/GITHUB_PAT 2>/dev/null || true)" \
    Rscript -e 'install.packages("remotes") ; repos <- BiocManager::repositories() ; remotes::install_local(path = "/opt/pkg/", repos=repos, dependencies=TRUE, build_vignettes=FALSE, upgrade=TRUE) ; sessioninfo::session_info(installed.packages()[,"Package"], include_base = TRUE)'

## Temporary fix (Bioconductor devel, Nov. 2026): BiocFileCache 3.3.1 cannot
## download any ExperimentHub resource (e.g. `TENxPBMCData`). This is fixed in
## BiocFileCache 3.3.2: install it from the Bioconductor git server until it
## reaches the devel repositories. This step does nothing once BiocFileCache
## >= 3.3.2 is installed, and can then be removed
RUN Rscript -e 'if (packageVersion("BiocFileCache") < "3.3.2") remotes::install_git("https://git.bioconductor.org/packages/BiocFileCache", upgrade = "never")' && \
    Rscript -e 'packageVersion("BiocFileCache")'

## Install the book's conda environment (`inst/requirements.yml`), which
## provides the command-line tools used in the labs (`samtools`, `seqkit`,
## `fastqc`, `ffq`), point every R session in the image at it, and append its
## `bin/` folder to the `PATH` of R sessions and of terminals (including the
## `RStudio` terminal), after the system folders
RUN Rscript -e 'BiocBook::setup_python(requirements = "/opt/pkg/inst/requirements.yml"); cat(sprintf("RETICULATE_PYTHON=%s\n", reticulate::py_config()$python), file = file.path(R.home("etc"), "Renviron.site"), append = TRUE)' && \
    ENV_BIN="$(Rscript -e 'cat(dirname(reticulate::py_config()$python))')" && \
    echo "PATH=\${PATH}:${ENV_BIN}" >> "$(R RHOME)/etc/Renviron.site" && \
    echo "export PATH=\"\${PATH}:${ENV_BIN}\"" > /etc/profile.d/scrnaseqphysalia2026.sh && \
    echo "export PATH=\"\${PATH}:${ENV_BIN}\"" >> /etc/bash.bashrc && \
    chmod -R a+rX "${R_USER_CACHE_DIR}" && \
    Rscript -e 'stopifnot(nzchar(Sys.getenv("RETICULATE_PYTHON"))); cat("book python:", reticulate::py_config()$python, "\n")' && \
    bash -lc 'samtools --version | head -n 1 && seqkit version && fastqc --version && ffq --version'

## Build/install using same approach than BBS
RUN R CMD INSTALL /opt/pkg
RUN R CMD build --keep-empty-dirs --no-resave-data /opt/pkg

## Let every user (e.g. `rstudio`) read and update the datasets cached while
## rendering the book
RUN chmod -R a+rwX "${R_USER_CACHE_DIR}"
