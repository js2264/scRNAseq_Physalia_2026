# 3  Lab 1: Familiarizing yourself with the workshop environment

> **NOTE:**
>
> The estimated time for this lab is around 1h15.
>
> This lab only uses the `terminal` and the `R` console interactively: none of the code on this page is run when the book is built.

> **TIP:**
>
> - Launch your own workshop session on the Bioconductor Orchestra platform.
> - Familiarize yourself with the `RStudio` interface.
> - Learn how to use the terminal within `RStudio`.
> - Learn how to find and download public single-cell data from the web.

## 3.1 Connect to RStudio Server

Most of single-cell RNA-seq analysis takes place either in `python` or in `R`. Here, we focus on how to leverage `R` to investigate scRNAseq data. `RStudio` is an `IDE` (Integrated Development Environment, in other words: a nice graphical interface to run `R`-related commands).

For this workshop, every participant gets their own `RStudio Server` session, running in the cloud from the **Docker image built from this book**. Follow the steps described in [Getting started](../../pages/getting-started.md) to launch it from the Orchestra platform: a new tab opens in your browser with `RStudio` already logged in.

Notice how when you log in `RStudio`, there are multiple panels. Familiarize yourself with the different panels.

The interactive `R` console is generally found in the bottom left corner of `RStudio` (though it may be in another corner sometimes). All the rest (history panel, environment panel, directory explorer panel, editor panel) are extra features provided by `RStudio`.

> **TIP:**
>
> Within the `R` console, you can safely use `R`-dedicated commands. Do you know the most common ones? The semantics are a different from the `terminal` commands you may be used to…
>
> ``` r
> getwd() # equivalent of `pwd` in terminal
> dir.create('~/data/') # equivalent of `mkdir ~/data/` in terminal
> setwd("~/data/") # equivalent of `cd ~/data/` in terminal
> list.files("~/data/") # equivalent of `ls` in terminal
> download.file("...") # equivalent of `wget ...` in terminal
> ```

## 3.2 Use a terminal within RStudio

A general issue with bioinformatic analyses stems from the fact that nobody works in the same environment:

- Are you working on Mac? Linux? Windows?
- Do you have a lot of computational power? Perhaps a GPU card?
- Are you connected to the Internet? With a fast connection? Are you working behind a proxy?
- Which versions of `R`, `Bioconductor` and of each package are installed?

To ensure that we are all working in the same environment, we all run the **same Docker image**. A Docker image is a frozen, self-contained Linux system: it ships a given version of `R`, `Bioconductor`, every `R` package and every command-line tool used in this workshop, and even the datasets used in the labs. Each participant runs their own copy of this image (a *container*) in the cloud, so:

- **Everybody gets exactly the same software**, down to package versions;
- **Your session is yours**: you can install, modify or delete files without affecting anyone else;
- **Your session is temporary**: download the files you want to keep before the end of the workshop!

The easiest way to launch **`bash` commands from a terminal** is to do it through `RStudio`: You can open up a `terminal` directly from within `RStudio` as follows: go to `Tools` \> `Terminal` \> `New terminal`. This should open up a new tab in the bottom left corner (next to the `R` console).

> **WARNING:**
>
> From here onwards, be sure you completely understand the difference between **“`R console`”** and **“`terminal` (or `shell`)”**. They are entirely different things, and can be both accessed within `RStudio`. It is crucial you understand the difference between the two to not get confused for the rest of the course.

## 3.3 Basic terminal commands

One can list files, download files, check help pages, …, just like in `R`.

- Check the your present directory

``` sh
pwd
```

- Check history

``` sh
history
```

- put history into a `history.txt` file, in your home folder (`~`)

``` sh
history > ~/history.txt
```

- make a new folder called data

``` sh
mkdir ~/data
```

- Go to the new `data` directory

``` sh
cd ~/data
```

- move `history.txt` file into `data` directory

``` sh
mv ../history.txt ./
```

- check the help of the `wget` command

``` sh
wget --help
```

> **NOTE:**
>
> On most Linux and macOS computers, `man <command>` shows the full manual page of a command (e.g. `man wget`, hit `q` to exit). Like in most Docker images, manual pages are not installed in the workshop environment: use `<command> --help` instead.

- check specific help for `samtools` command and subcommands

``` sh
samtools --help
samtools view --help
```

- redirect the help output into a file called `samtools-view-help.txt`

``` sh
samtools view --help > samtools-view-help.txt
```

- Download a file from Internet with `wget`

``` sh
wget https://cf.10xgenomics.com/supp/cell-exp/cellranger-tiny-bcl-1.2.0.tar.gz
```

- List all files in a folder

``` sh
ls -lh ~/data/
```

> **TIP:**
>
> Download the `git` repository for this course from `GitHub`:
>
> ``` sh
> git clone https://github.com/js2264/scRNAseq_Physalia_2026.git
> ```
>
> This downloads the repository for this course to your home folder in your workshop session. To get it on your local computer (to save the lectures and exercises), go to [the GitHub repo page](https://github.com/js2264/scRNAseq_Physalia_2026), click on the green `Code` button, then `Download ZIP`. Beware, the download may take a significant time based on your internet connection (~100 MB).

## 3.4 Single-cell RNA-seq datasets

> “This is a course about single-cell RNA-seq analysis, after all, so where is my data?”

Ok, **“your”** data is (most likely) yet to be sequenced! Or maybe you’re interested in digging already existing databases!

The [Human Cell Atlas data portal](https://data.humancellatlas.org/explore/projects) and the [CZ CELLxGENE Discover](https://cellxgene.cziscience.com/datasets) portal are probably good places to start digging, if you are interested in mammal-related studies. For instance, let’s say I am interested in epididymis differentiation. Boom: here is an entry from the HCA focusing on epididymis: [link to HCA data portal](https://data.humancellatlas.org/explore/projects/842605c7-375a-47c5-9e2c-a71c2c00fcad).

### 3.4.1 Raw fastq reads from GEO

Here is the link to the actual paper studying epididymis: [An atlas of human proximal epididymis reveals cell-specific functions and distinct roles for CFTR](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7471510/).

> **NOTE:**
>
> Find and check out the corresponding `GEO` entries for this study. What type of sequencing data is available?
>
> > **TIP:**
> >
> > Here is the link to the `GEO` page: [link](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE148963).

> **NOTE:**
>
> Can you find links to download the raw data from this paper?

There are several ways to find this information, e.g. `ffq` command line tool, the web-based `sra-explorer` page ([here](https://sra-explorer.info/)), or the [ENA browser](https://www.ebi.ac.uk/ena/browser/home). You generally will need the `GEO` corresponding ID or SRA project ID (e.g. `SRPxxxxxx`…).

> **NOTE:**
>
> The [`ffq` tool](https://github.com/pachterlab/ffq) from the Pachter lab is already installed in your workshop session. Check `ffq` help and try fetching the metadata of the sequencing run of sample `AXH009` (`GSM4486714`): `SRR11575369`.
>
> > **TIP:**
> >
> > `ffq` is a `python` package, available from `PyPI`:
> >
> > ``` sh
> > pip install ffq
> > ```
>
> > **NOTE:**
> >
> > ``` sh
> > ffq --help
> > ffq SRR11575369 > SRR11575369_ffq.json
> > head -n 30 SRR11575369_ffq.json
> > ```
> >
> > `ffq` returns the metadata of the run in `JSON` format: the experiment, study and sample it belongs to, and the download links of its files.
> >
> > `ffq` can, in principle, also be queried with GEO accessions directly (e.g. `ffq GSE148963`). With the current release of `ffq` (0.3.1), these queries (as well as queries of SRA studies, `SRP...`) fail, as NCBI changed the GEO web pages `ffq` reads. Queries of SRA runs (`SRR...`) work fine.

> **NOTE:**
>
> Can you find the links to raw data associated with the `GSE148963` GEO ID?
>
> > **NOTE:**
> >
> > `ffq` lists the download links of the files of the run, from several mirrors (`ftp`: the ENA, `aws`, `gcp` and `ncbi`). You can use a `grep` command: `grep` returns the lines which match a given pattern (e.g. a link…)!
> >
> > ``` sh
> > grep '"url": "ftp://' SRR11575369_ffq.json
> > ```
> >
> > And with a bit of `sed` magick…
> >
> > ``` sh
> > grep '"url": "ftp://' SRR11575369_ffq.json | sed 's,.*"ftp://,ftp://,' | sed 's,".*,,' > SRR11575369_fastqlist.txt
> > cat SRR11575369_fastqlist.txt
> > ## `ffq` looks through the SRA/ENA databases to find metadata associated with the `SRR11575369` run (including download links)
> > ## grep '"url": "ftp://' recovers the text lines that contain FTP download links
> > ## the `sed` commands clean up the text lines
> > ## the `wget` command downloads locally the links listed in the generated file (SRR11575369_fastqlist.txt)
> > # wget -i SRR11575369_fastqlist.txt ## Do not run, it would take too long (2 x 35 GB)...
> > ```
> >
> > To get the runs of every sample of the study at once, you can directly query the [ENA portal API](https://www.ebi.ac.uk/ena/portal/api/) with the SRA study ID of `GSE148963`: `SRP257542` (listed in the `Relations` section of the GEO page, and in the `study` field returned by `ffq`).
> >
> > ``` sh
> > curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRP257542&result=read_run&fields=run_accession,sample_alias,sample_title,fastq_bytes,fastq_ftp" > SRP257542_runs.tsv
> > cut -f 1-4 SRP257542_runs.tsv
> > cut -f 5 SRP257542_runs.tsv | tail -n +2 | tr ';' '\n' > GSE148963_fastqlist.txt
> > cat GSE148963_fastqlist.txt
> > ```

> **NOTE:**
>
> The full `fastq` files for a single run of this study weigh ~35 GB each… Let’s only download the first 25,000 reads of the `SRR11575369` run (sample `AXH009`), by streaming the file from the ENA and stopping after 100,000 lines (each read spans 4 lines). Have a look at it!
>
> > **NOTE:**
> >
> > ``` sh
> > mkdir -p ~/data && cd ~/data
> > URL=https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR115/069/SRR11575369
> > for READ in 1 2; do
> >     echo "Streaming the first 25,000 reads of SRR11575369_${READ}.fastq.gz..."
> >     curl -s ${URL}/SRR11575369_${READ}.fastq.gz \
> >         | zcat 2>/dev/null \
> >         | head -n 100000 \
> >         | gzip > SRR11575369_${READ}.fq.gz
> > done
> > ls -lh SRR11575369_*.fq.gz
> > zcat SRR11575369_1.fq.gz | head -n 12
> > ```
> >
> > `curl` streams the compressed file, `zcat` decompresses it on the fly, `head` stops the stream after 100,000 lines and `gzip` compresses these lines into a new file. When `head` stops, `curl` and `zcat` complain that they cannot write anymore: `-s` and `2>/dev/null` silence these (expected) errors.
> >
> > Each file takes ~30 seconds to download and weighs ~2 MB. If a file only weighs a few bytes, the download failed (e.g. the ENA server is busy): a copy of both files is available in `~/Share/`.
> >
> > ``` sh
> > cp ~/Share/SRR11575369_*.fq.gz ~/data/
> > ```
>
> > **TIP:**
> >
> > Running these commands on your own Mac? There, `zcat file.gz` looks for a `file.gz.Z` file and fails: use `gzcat file.gz` (or `gunzip -c file.gz`) instead. Reading a stream (`... | zcat`) works as on Linux.

> **NOTE:**
>
> Try and understand the structure of the `fq.gz` file. What is the meaning of each line? What is the difference between the `_1` and `_2` files?
>
> How many reads are there in the `fq.gz` file? And how long are they? Can you get a summary of what is in this file? All of these questions can be quickly answered using the [`seqkit`](https://bioinf.shenwei.me/seqkit/usage/) tool:
>
> > **NOTE:**
> >
> > ``` sh
> > seqkit --help
> > seqkit stats --help
> > seqkit stats ~/data/SRR11575369_*.fq.gz
> > ```

### 3.4.2 Processed count matrices

Many times, researchers will provide a filtered count matrix when they publish scRNAseq experiments (along with mandatory raw `fastq` data, of course). It’s way lighter than `fastq` reads, and you can go ahead with downstream analyses a lot quicker. So how do you get these matrices? Human Cell Atlas Consortium provides many processed datasets. For instance, in our case, the `A single cell atlas of the human proximal epididymis reveals cell-type specific functions and distinct roles for CFTR (Leir et al, 2020)` study is available at the [following link](https://data.humancellatlas.org/explore/projects/842605c7-375a-47c5-9e2c-a71c2c00fcad).

GEO also provides access to both raw and processed files.

> **NOTE:**
>
> Find GEO-hosted processed files for the `Leir et al` study.
>
> > **NOTE:**
> >
> > Check the [Geo Dataset search page](https://www.ncbi.nlm.nih.gov/gds/).
> >
> > [Search for `Leir 2020 CFTF`](https://www.ncbi.nlm.nih.gov/gds/?term=leir+2020+CFTR). You should find the [following entry](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE148963).
> >
> > Scrolling down to the bottom of the page, there is a box labelled “Supplementary data”. Each sample (`GSM...`) also has its own supplementary files: for instance, the processed files for sample `AXH009` are listed on the [GSM4486714 page](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4486714).
>
> Download and check the content of the count matrix, the genes and the barcodes files.
>
> What type of information does each file contain? How is it formatted? is it easily imported in R?
>
> > **NOTE:**
> >
> > ``` sh
> > cd ~/data
> > URL=https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4486nnn/GSM4486714/suppl
> > for FILE in matrix.mtx genes.tsv barcodes.tsv; do
> >     wget ${URL}/GSM4486714_AXH009_${FILE}.gz
> > done
> > ls -lh GSM4486714_AXH009_*
> > zcat GSM4486714_AXH009_genes.tsv.gz | head
> > zcat GSM4486714_AXH009_barcodes.tsv.gz | head
> > zcat GSM4486714_AXH009_matrix.mtx.gz | head
> > ```
> >
> > These three files are the standard output of `cellranger` (a “Matrix Market” sparse matrix, see [Lab 2](../../pages/day1/Lab2_processingreads.md)): the `.mtx` file only lists the non-zero entries of the count matrix, one per line, as `<gene index> <cell index> <count>`, after a 3-line header. The genes and barcodes files give the names of the rows and of the columns of the matrix. We will see how to import them in `R` in [Lab 3](../../pages/day2/Lab3_Rbioc.md).
>
> - How many cells were sequenced? How many genes were counted?
> - Is it easy to interpret the count matrix? Why is it in such format?
> - Comment on the file sizes between processed count matrix files and raw reads.

## 3.5 Bonus

For those of you who are already familiar with the basics, you can fast-forward through this lab and [Lab 2](../../pages/day1/Lab2_processingreads.md), and start working on scRNAseq data in `R` directly with [Lab 3](../../pages/day2/Lab3_Rbioc.md).
