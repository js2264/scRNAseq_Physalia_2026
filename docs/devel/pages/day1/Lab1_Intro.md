# 3  Lab 1: Familiarizing yourself with the workshop environment

> **NOTE:**
>
> The estimated time for this lab is around 1h15.
>
> The `terminal` commands of this lab are run when the book is built, in the same environment as your workshop session: the outputs shown on this page are what you should get.

> **TIP:**
>
> - Launch your own workshop session on the Bioconductor Orchestra platform.
> - Familiarize yourself with the `RStudio` interface.
> - Learn how to use the terminal within `RStudio`.
> - Learn how to find and download public single-cell data from the web.

## 3.1 Connect to RStudio Server on Orchestra

Most of single-cell RNA-seq analysis takes place either in `python` or in `R`. Here, we focus on how to leverage `R` to investigate scRNAseq data. `RStudio` is an `IDE` (Integrated Development Environment, in other words: a nice graphical interface to run `R`-related commands).

For this workshop, every participant gets their own `RStudio Server` session, running in the cloud from the **Docker image built from this book**, on the [Orchestra platform](https://orchestraplatform.org/). Orchestra is an open-source platform, supported by the Bioconductor community, which runs workshop environments in the cloud: each participant gets their own, isolated `RStudio` session, in their web browser. There is nothing to install on your computer, and no password to type.

1.  Open the Orchestra dashboard: **the link will be given by your instructors at the beginning of the workshop**.
2.  Sign in with a Google account (any Google account works, either a personal Gmail address or an institutional account backed by Google).
3.  Pick the **“Single-cell RNA-seq analysis with R/Bioconductor (Physalia 2026)”** workshop from the catalog and click on **Launch**. Starting the session can take a few minutes the first time (the status moves from `Pending` to `Starting` to `Ready`).
4.  When it’s `Ready`, click on **Open**: `RStudio` opens in a new tab of your browser, already logged in.

> **WARNING:**
>
> Every session has a time limit, shown as a countdown on its card in the dashboard (you can extend it by clicking on `+1h`). **When a session expires (or when you terminate it), everything in it is removed**, including the files you created.
>
> Download anything you want to keep (scripts, notebooks, results) before the end of each day: in `RStudio`, tick the files in the `Files` pane, then click on `More` \> `Export...`.
>
> If Orchestra is unavailable, you can run the same environment on your own computer with Docker: see [Getting started](../../pages/getting-started.md).

Notice how, once `RStudio` is open, there are multiple panels. Familiarize yourself with the different panels.

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

> **NOTE:**
>
> On this page, each block of `bash` code is run in a new terminal, which starts in your home folder (`~`): this is why some blocks start with `cd ~/data`. In your own terminal, moving to `~/data` once is enough.

- Check the your present directory

``` bash
pwd
##  /home/rstudio
```

- Check history (in your terminal, this lists the commands you typed; on this page, each block runs in a new terminal, so the history is empty)

``` bash
history
```

- put history into a `history.txt` file, in your home folder (`~`)

``` bash
history > ~/history.txt
```

- make a new folder called data (`-p`: do not complain if it already exists)

``` bash
mkdir -p ~/data
```

- Go to the new `data` directory

``` bash
cd ~/data
pwd
##  /home/rstudio/data
```

- move `history.txt` file into `data` directory

``` bash
cd ~/data
mv ../history.txt ./
ls -l
##  total 0
##  -rw-r--r-- 1 root root 0 Nov 13 15:35 history.txt
```

- check the help of the `wget` command (here, only its first 20 lines, with `head`)

``` bash
wget --help | head -n 20
##  GNU Wget 1.21.4, a non-interactive network retriever.
##  Usage: wget [OPTION]... [URL]...
##  
##  Mandatory arguments to long options are mandatory for short options too.
##  
##  Startup:
##    -V,  --version                   display the version of Wget and exit
##    -h,  --help                      print this help
##    -b,  --background                go to background after startup
##    -e,  --execute=COMMAND           execute a `.wgetrc'-style command
##  
##  Logging and input file:
##    -o,  --output-file=FILE          log messages to FILE
##    -a,  --append-output=FILE        append messages to FILE
##    -d,  --debug                     print lots of debugging information
##    -q,  --quiet                     quiet (no output)
##    -v,  --verbose                   be verbose (this is the default)
##    -nv, --no-verbose                turn off verboseness, without being quiet
##         --report-speed=TYPE         output bandwidth as TYPE.  TYPE can be bits
##    -i,  --input-file=FILE           download URLs found in local or external FILE
```

> **NOTE:**
>
> On most Linux and macOS computers, `man <command>` shows the full manual page of a command (e.g. `man wget`, hit `q` to exit). Like in most Docker images, manual pages are not installed in the workshop environment: use `<command> --help` instead.

- check specific help for `samtools` command and subcommands

``` bash
samtools --help
samtools view --help | head -n 30
##  
##  Program: samtools (Tools for alignments in the SAM format)
##  Version: 1.21 (using htslib 1.24)
##  
##  Usage:   samtools <command> [options]
##  
##  Commands:
##    -- Indexing
##       dict           create a sequence dictionary file
##       faidx          index/extract FASTA
##       fqidx          index/extract FASTQ
##       index          index alignment
##  
##    -- Editing
##       calmd          recalculate MD/NM tags and '=' bases
##       fixmate        fix mate information
##       reheader       replace BAM header
##       targetcut      cut fosmid regions (for fosmid pool only)
##       addreplacerg   adds or replaces RG tags
##       markdup        mark duplicates
##       ampliconclip   clip oligos from the end of reads
##  
##    -- File operations
##       collate        shuffle and group alignments by name
##       cat            concatenate BAMs
##       consensus      produce a consensus Pileup/FASTA/FASTQ
##       merge          merge sorted alignments
##       mpileup        multi-way pileup
##       sort           sort alignment file
##       split          splits a file by read group
##       quickcheck     quickly check if SAM/BAM/CRAM file appears intact
##       fastq          converts a BAM to a FASTQ
##       fasta          converts a BAM to a FASTA
##       import         Converts FASTA or FASTQ files to SAM/BAM/CRAM
##       reference      Generates a reference from aligned data
##       reset          Reverts aligner changes in reads
##  
##    -- Statistics
##       bedcov         read depth per BED region
##       coverage       alignment depth and percent coverage
##       depth          compute the depth
##       flagstat       simple stats
##       idxstats       BAM index stats
##       cram-size      list CRAM Content-ID and Data-Series sizes
##       phase          phase heterozygotes
##       stats          generate stats (former bamcheck)
##       ampliconstats  generate amplicon specific stats
##  
##    -- Viewing
##       flags          explain BAM flags
##       head           header viewer
##       tview          text alignment viewer
##       view           SAM<->BAM<->CRAM conversion
##       depad          convert padded BAM to unpadded BAM
##       samples        list the samples in a set of SAM/BAM/CRAM files
##  
##    -- Misc
##       help [cmd]     display this help message or help for [cmd]
##       version        detailed version information
##  
##  
##  Usage: samtools view [options] <in.bam>|<in.sam>|<in.cram> [region ...]
##  
##  Output options:
##    -b, --bam                  Output BAM
##    -C, --cram                 Output CRAM (requires -T)
##    -1, --fast                 Use fast BAM compression (and default to --bam)
##    -u, --uncompressed         Uncompressed BAM output (and default to --bam)
##    -h, --with-header          Include header in SAM output
##    -H, --header-only          Print SAM header only (no alignments)
##        --no-header            Print SAM alignment records only [default]
##    -c, --count                Print only the count of matching records
##    -o, --output FILE          Write output to FILE [standard output]
##    -U, --unoutput FILE, --output-unselected FILE
##                               Output reads not selected by filters to FILE
##    -p, --unmap                Set flag to UNMAP on reads not selected
##                               then write to output file.
##    -P, --fetch-pairs          Retrieve complete pairs even when outside of region
##  Input options:
##    -t, --fai-reference FILE   FILE listing reference names and lengths
##    -M, --use-index            Use index and multi-region iterator for regions
##        --region[s]-file FILE  Use index to include only reads overlapping FILE
##    -X, --customized-index     Expect extra index file argument after <in.bam>
##  
##  Filtering options (Only include in output reads that...):
##    -L, --target[s]-file FILE  ...overlap (BED) regions in FILE
##    -N, --qname-file [^]FILE   ...whose read name is listed in FILE ("^" negates)
##    -r, --read-group STR       ...are in read group STR
##    -R, --read-group-file [^]FILE
##                               ...are in a read group listed in FILE
```

- redirect the help output into a file called `samtools-view-help.txt`

``` bash
cd ~/data
samtools view --help > samtools-view-help.txt
head -n 5 samtools-view-help.txt
##  
##  Usage: samtools view [options] <in.bam>|<in.sam>|<in.cram> [region ...]
##  
##  Output options:
##    -b, --bam                  Output BAM
```

- Download a file from Internet with `wget` (`-nv`: only print a summary line, rather than a progress bar)

``` bash
cd ~/data
wget -nv https://cf.10xgenomics.com/supp/cell-exp/cellranger-tiny-bcl-1.2.0.tar.gz
##  2026-11-13 15:35:50 URL:https://cf.10xgenomics.com/supp/cell-exp/cellranger-tiny-bcl-1.2.0.tar.gz [136259701/136259701] -> "cellranger-tiny-bcl-1.2.0.tar.gz" [1]
```

- List all files in a folder

``` bash
ls -lh ~/data/
##  total 130M
##  -rw-r--r-- 1 root root 130M Mar  5  2017 cellranger-tiny-bcl-1.2.0.tar.gz
##  -rw-r--r-- 1 root root    0 Nov 13 15:35 history.txt
##  -rw-r--r-- 1 root root 5.8K Nov 13 15:35 samtools-view-help.txt
```

> **TIP:**
>
> Download the `git` repository for this course from `GitHub` (`--depth 1`: without the history of the repository; `--quiet`: without progress messages):
>
> ``` bash
> git clone --quiet --depth 1 https://github.com/js2264/scRNAseq_Physalia_2026.git
> ls scRNAseq_Physalia_2026/inst/pages
> ##  day1
> ##  day2
> ##  day3
> ##  day4
> ##  day5
> ##  extra.qmd
> ##  getting-started.qmd
> ##  images
> ##  prerequisites.qmd
> ##  program.qmd
> ```
>
> This downloads the repository for this course to your home folder in your workshop session. To get it on your local computer (to save the lectures and exercises), go to [the GitHub repo page](https://github.com/js2264/scRNAseq_Physalia_2026), click on the green `Code` button, then `Download ZIP`. Beware, the download may take a significant time based on your internet connection (~270 MB).

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
> > **NOTE:**
> >
> > ``` bash
> > cd ~/data
> > ffq --help
> > ffq SRR11575369 > SRR11575369_ffq.json
> > head -n 30 SRR11575369_ffq.json
> > ##  usage: ffq [-h] [-o OUT] [-l LEVEL] [--ftp] [--aws] [--gcp] [--ncbi] [--split]
> > ##             [--verbose] [--version]
> > ##             IDs [IDs ...]
> > ##  
> > ##  ffq 0.3.1: A command line tool to find sequencing data from SRA / GEO / ENCODE
> > ##  / ENA / EBI-EMBL / DDBJ / Biosample.
> > ##  
> > ##  positional arguments:
> > ##    IDs         One or multiple SRA / GEO / ENCODE / ENA / EBI-EMBL / DDBJ /
> > ##                Biosample accessions, DOIs, or paper titles
> > ##  
> > ##  options:
> > ##    -h, --help  Show this help message and exit
> > ##    -o OUT      Path to write metadata (default: standard out)
> > ##    -l LEVEL    Max depth to fetch data within accession tree
> > ##    --ftp       Return FTP links
> > ##    --aws       Return AWS links
> > ##    --gcp       Return GCP links
> > ##    --ncbi      Return NCBI links
> > ##    --split     Split output into separate files by accession (`-o` is a
> > ##                directory)
> > ##    --verbose   Print debugging information
> > ##    --version   show program's version number and exit
> > ##  [2026-11-13 15:35:59,340]    INFO Parsing run SRR11575369
> > ##  {
> > ##      "SRR11575369": {
> > ##          "accession": "SRR11575369",
> > ##          "experiment": "SRX8143297",
> > ##          "study": "SRP257542",
> > ##          "sample": "SRS6506205",
> > ##          "title": "Illumina NovaSeq 6000 paired end sequencing; GSM4486714: single cell RNA seq of human epididymis epithelial tissue [AXH009]; Homo sapiens; RNA-Seq",
> > ##          "attributes": {
> > ##              "loader": "latf-load",
> > ##              "ENA-STATUS": "PUBLIC",
> > ##              "ENA-FIRST-PUBLIC": "2020-08-27",
> > ##              "ENA-LAST-UPDATE": "2020-08-27",
> > ##              "ENA-BASE-COUNT": 136010047290,
> > ##              "ENA-SPOT-COUNT": 450364395
> > ##          },
> > ##          "files": {
> > ##              "ftp": [
> > ##                  {
> > ##                      "accession": "SRR11575369",
> > ##                      "filename": "SRR11575369_1.fastq.gz",
> > ##                      "filetype": "fastq",
> > ##                      "filesize": 34995132199,
> > ##                      "filenumber": 1,
> > ##                      "md5": "f1c825a9a7ca9ce11a872fcd8777bd5e",
> > ##                      "urltype": "ftp",
> > ##                      "url": "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR115/069/SRR11575369/SRR11575369_1.fastq.gz"
> > ##                  },
> > ##                  {
> > ##                      "accession": "SRR11575369",
> > ##                      "filename": "SRR11575369_2.fastq.gz",
> > ```
> >
> > `ffq` returns the metadata of the run in `JSON` format: the experiment, study and sample it belongs to, and the download links of its files.
> >
> > > **WARNING:**
> > >
> > > `ffq` queries the NCBI and ENA servers, which limit the number of requests they accept from the same place: when all the participants of the workshop run `ffq` at the same time, NCBI may reject some of the queries (`ERROR 429 Client Error: Too Many Requests`), and `ffq` then writes an empty file. Wait a little and try again or, if it keeps failing, use the copy of the result provided in `~/Share/`:
> > >
> > > ``` bash
> > > cd ~/data
> > > ## Use the copy in ~/Share/ if `ffq` failed (i.e. if SRR11575369_ffq.json is empty)
> > > [ -s SRR11575369_ffq.json ] || cp ~/Share/SRR11575369_ffq.json ./
> > > ls -l SRR11575369_ffq.json
> > > ##  -rw-r--r-- 1 root root 3420 Nov 13 15:36 SRR11575369_ffq.json
> > > ```
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
> > ``` bash
> > cd ~/data
> > grep '"url": "ftp://' SRR11575369_ffq.json
> > ##                      "url": "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR115/069/SRR11575369/SRR11575369_1.fastq.gz"
> > ##                      "url": "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR115/069/SRR11575369/SRR11575369_2.fastq.gz"
> > ```
> >
> > And with a bit of `sed` magick…
> >
> > ``` bash
> > cd ~/data
> > grep '"url": "ftp://' SRR11575369_ffq.json | sed 's,.*"ftp://,ftp://,' | sed 's,".*,,' > SRR11575369_fastqlist.txt
> > cat SRR11575369_fastqlist.txt
> > ## `ffq` looks through the SRA/ENA databases to find metadata associated with the `SRR11575369` run (including download links)
> > ## grep '"url": "ftp://' recovers the text lines that contain FTP download links
> > ## the `sed` commands clean up the text lines
> > ## the `wget` command downloads locally the links listed in the generated file (SRR11575369_fastqlist.txt)
> > # wget -i SRR11575369_fastqlist.txt ## Do not run, it would take too long (2 x 35 GB)...
> > ##  ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR115/069/SRR11575369/SRR11575369_1.fastq.gz
> > ##  ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR115/069/SRR11575369/SRR11575369_2.fastq.gz
> > ```

Each of these two `fastq` files weighs ~35 GB: downloading them would take far too long. Let’s pretend we downloaded them with `wget -i SRR11575369_fastqlist.txt`: their first 25,000 reads are available in `~/Share/`.

``` bash
ls -lh ~/Share/SRR11575369_*.fq.gz
##  -rw-r--r-- 1 root root 1.9M Nov  6 18:09 /home/rstudio/Share/SRR11575369_1.fq.gz
##  -rw-r--r-- 1 root root 1.9M Nov  6 18:09 /home/rstudio/Share/SRR11575369_2.fq.gz
```

> **NOTE:**
>
> Try and understand the structure of the `fq.gz` files. What is the meaning of each line? What is the difference between the `_1` and `_2` files?
>
> How many reads are there in each `fq.gz` file? And how long are they? Can you get a summary of what is in these files? All of these questions can be quickly answered using the [`seqkit`](https://bioinf.shenwei.me/seqkit/usage/) tool:
>
> > **NOTE:**
> >
> > ``` bash
> > zcat ~/Share/SRR11575369_1.fq.gz | head -n 8
> > zcat ~/Share/SRR11575369_2.fq.gz | head -n 8
> > ##  @SRR11575369.1 1/1
> > ##  ANCAACAGTGGAATTCGTTTTGATGAAAAAAATAAATTTGTTCTTCAAAGCAGAGTGAATGATGCAGTACGAGCTCTTGCTCTTGAAAACCCATCACAACTTTATAATTTAATAATTAGTGAAAAATTAAAAAAATAATATTTCTTATATT
> > ##  +
> > ##  F#F:FFFFF:FFFFFFFFFF:F,FF,FF::F,::FFFFFFFFFFF:FF:FF:F,FFF,F:F,FFF,FF:FF:F:FFF:FF:FFFF,:F::FFFF:FF,FF:F:FF,,:F:F,::,F:F:FF,,::,FF:,,F,,,,,:FFFF::F:,,,FF
> > ##  @SRR11575369.2 2/1
> > ##  TNGCCAGTCAATACCGCCCTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTGTGTTTTTTTTTTTTTATATAATAAAATTTTTTTTTTTTTTAATAAAATTTTTTTTTTATTTT
> > ##  +
> > ##  F#FFFFFFF::FF:FF::FFFFFFFFFFFFFFFF:FFFFFFFFFFFFF,,FFFFFFFFFF:FFFFFF:FFFFFFFFFFF:F,FFFFF,,,:,,::FFFFFFFF,,:,,,F,,,,,,FFFFFF,F,,FF:,F,:,F,FFFFF,,,:,,F::F
> > ##  @SRR11575369.1 1/2
> > ##  TCCCGCAGACCAGATGATAATCGTATTTTTTTTTTCTATTTTTTTAATATTAATTAAATTCGGAATTACAAATTTTTAAAAAGTTTGAAAATGTAGAGAATTTATAAAAATTATTGTTAAATTTGTTGTGTCGAAATCTAAGAAATATTGT
> > ##  +
> > ##  F,:FFF,,,,F,,:FF,F:F:F,:F:FFF,FFF:F,::F,::FF:,,F,:F:,:,,,,F:,:F,F:F,:,F,FFFFF,F,::FFF:F,::,FFFF:,F,F:FF,F,,,:,,F,FF,FF:,,FFFFFF,FF,,,,F,F,F,,,:,FF:FF,:
> > ##  @SRR11575369.2 2/2
> > ##  GCCCAGAGCCCTGCTATCATCGATTCGGACCCAGTGGATGAGGAAGTGCTGATGTCGCTGGTGGTGGAACTGGGGTTGGACCGAGCCAATGAGCTTCCGGAGCTGTGGCTGGGGCAGAATGAGTTTGACTTCACTGCGGACTTTCCATCTA
> > ##  +
> > ##  FFF:,F,FFFFFFFFFF:FFFF,FFF,FFFFFFFFFFFFFF:FFFFFFFFF,FF,FF:F:FFFFFFFFFFFFFFF,FF:FFFFFFFFFF:FFF:F,FF:FFFF:,FFF:FFFFFFFFFF:FF:FFF:,FFFFFFF,,,F,FFF::FFFFFF
> > ```
> >
> > ``` bash
> > seqkit --help | head -n 30
> > seqkit stats --help | head -n 20
> > seqkit stats --quiet ~/Share/SRR11575369_*.fq.gz
> > ##  SeqKit -- a cross-platform and ultrafast toolkit for FASTA/Q file manipulation
> > ##  
> > ##  Version: 2.9.0
> > ##  
> > ##  Author: Wei Shen <shenwei356@gmail.com>
> > ##  
> > ##  Documents  : http://bioinf.shenwei.me/seqkit
> > ##  Source code: https://github.com/shenwei356/seqkit
> > ##  Please cite:
> > ##    1. https://doi.org/10.1002/imt2.191
> > ##    2. https://doi.org/10.1371/journal.pone.0163962
> > ##  
> > ##  
> > ##  Seqkit utlizies the pgzip (https://github.com/klauspost/pgzip) package to
> > ##  read and write gzip file, and the outputted gzip file would be slighty
> > ##  larger than files generated by GNU gzip.
> > ##  
> > ##  Seqkit writes gzip files very fast, much faster than the multi-threaded pigz,
> > ##  therefore there's no need to pipe the result to gzip/pigz.
> > ##  
> > ##  Seqkit also supports reading and writing xz (.xz) and zstd (.zst) formats since v2.2.0.
> > ##  Bzip2 format is supported since v2.4.0.
> > ##  
> > ##  Compression level:
> > ##    format   range   default  comment
> > ##    gzip     1-9     5        https://github.com/klauspost/pgzip sets 5 as the default value.
> > ##    xz       NA      NA       https://github.com/ulikunitz/xz does not support.
> > ##    zstd     1-4     2        roughly equals to zstd 1, 3, 7, 11, respectively.
> > ##    bzip     1-9     6        https://github.com/dsnet/compress
> > ##  
> > ##  simple statistics of FASTA/Q files
> > ##  
> > ##  Columns:
> > ##  
> > ##    1.  file      input file, "-" for STDIN
> > ##    2.  format    FASTA or FASTQ
> > ##    3.  type      DNA, RNA, Protein or Unlimit
> > ##    4.  num_seqs  number of sequences
> > ##    5.  sum_len   number of bases or residues       , with gaps or spaces counted
> > ##    6.  min_len   minimal sequence length           , with gaps or spaces counted
> > ##    7.  avg_len   average sequence length           , with gaps or spaces counted
> > ##    8.  max_len   miximal sequence length           , with gaps or spaces counted
> > ##    9.  Q1        first quartile of sequence length , with gaps or spaces counted
> > ##    10. Q2        median of sequence length         , with gaps or spaces counted
> > ##    11. Q3        third quartile of sequence length , with gaps or spaces counted
> > ##    12. sum_gap   number of gaps
> > ##    13. N50       N50. https://en.wikipedia.org/wiki/N50,_L50,_and_related_statistics#N50
> > ##    14. N50_num   N50_num or L50. https://en.wikipedia.org/wiki/N50,_L50,_and_related_statistics#L50
> > ##    15. Q20(%)    percentage of bases with the quality score greater than 20
> > ##    16. Q30(%)    percentage of bases with the quality score greater than 30
> > ##  file                                     format  type  num_seqs    sum_len  min_len  avg_len  max_len
> > ##  /home/rstudio/Share/SRR11575369_1.fq.gz  FASTQ   DNA     25,000  3,775,000      151      151      151
> > ##  /home/rstudio/Share/SRR11575369_2.fq.gz  FASTQ   DNA     25,000  3,775,000      151      151      151
> > ```
> >
> > > **TIP:**
> > >
> > > Running these commands on your own Mac? There, `zcat file.gz` looks for a `file.gz.Z` file and fails: use `gzcat file.gz` (or `gunzip -c file.gz`) instead. Reading a stream (`... | zcat`) works as on Linux.

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
> > ``` bash
> > cd ~/data
> > URL=https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4486nnn/GSM4486714/suppl
> > for FILE in matrix.mtx genes.tsv barcodes.tsv; do
> >     wget -nv ${URL}/GSM4486714_AXH009_${FILE}.gz
> > done
> > ls -lh GSM4486714_AXH009_*
> > zcat GSM4486714_AXH009_genes.tsv.gz | head
> > zcat GSM4486714_AXH009_barcodes.tsv.gz | head
> > zcat GSM4486714_AXH009_matrix.mtx.gz | head
> > ##  2026-11-13 15:36:05 URL:https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4486nnn/GSM4486714/suppl/GSM4486714_AXH009_matrix.mtx.gz [30406214/30406214] -> "GSM4486714_AXH009_matrix.mtx.gz" [1]
> > ##  2026-11-13 15:36:06 URL:https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4486nnn/GSM4486714/suppl/GSM4486714_AXH009_genes.tsv.gz [257279/257279] -> "GSM4486714_AXH009_genes.tsv.gz" [1]
> > ##  2026-11-13 15:36:06 URL:https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4486nnn/GSM4486714/suppl/GSM4486714_AXH009_barcodes.tsv.gz [8649/8649] -> "GSM4486714_AXH009_barcodes.tsv.gz" [1]
> > ##  -rw-r--r-- 1 root root 8.5K Apr 20  2020 GSM4486714_AXH009_barcodes.tsv.gz
> > ##  -rw-r--r-- 1 root root 252K Apr 20  2020 GSM4486714_AXH009_genes.tsv.gz
> > ##  -rw-r--r-- 1 root root  29M Apr 20  2020 GSM4486714_AXH009_matrix.mtx.gz
> > ##  ENSG00000243485 MIR1302-10
> > ##  ENSG00000237613 FAM138A
> > ##  ENSG00000186092 OR4F5
> > ##  ENSG00000238009 RP11-34P13.7
> > ##  ENSG00000239945 RP11-34P13.8
> > ##  ENSG00000237683 AL627309.1
> > ##  ENSG00000239906 RP11-34P13.14
> > ##  ENSG00000241599 RP11-34P13.9
> > ##  ENSG00000228463 AP006222.2
> > ##  ENSG00000237094 RP4-669L17.10
> > ##  AAACGGGAGATGGGTC-1
> > ##  AAACGGGAGCTAGTTC-1
> > ##  AAACGGGGTTACCAGT-1
> > ##  AAACGGGGTTGGGACA-1
> > ##  AAACGGGTCCTAGTGA-1
> > ##  AAACGGGTCGGAGGTA-1
> > ##  AAAGATGGTCGAACAG-1
> > ##  AAAGATGGTTAAGACA-1
> > ##  AAAGATGTCACCTTAT-1
> > ##  AAAGCAAAGAAGGTGA-1
> > ##  %%MatrixMarket matrix coordinate integer general
> > ##  %
> > ##  32738 1876 9847254
> > ##  32722 1 1
> > ##  32709 1 42
> > ##  32708 1 1
> > ##  32707 1 2
> > ##  32706 1 26
> > ##  32705 1 2
> > ##  32704 1 28
> > ```
> >
> > These three files are the standard output of `cellranger` (a “Matrix Market” sparse matrix, see [Lab 2](../../pages/day1/Lab2_processingreads.md)): the `.mtx` file only lists the non-zero entries of the count matrix, one per line, as `<gene index> <cell index> <count>`, after a 3-line header. The genes and barcodes files give the names of the rows and of the columns of the matrix. We will see how to import them in `R` in [Lab 3](../../pages/day2/Lab3_Rbioc.md).
>
> - How many cells were sequenced? How many genes were counted?
> - Is it easy to interpret the count matrix? Why is it in such format?
> - Comment on the file sizes between processed count matrix files and raw reads.

## 3.5 Bonus

For those of you who are already familiar with the basics, you can fast-forward through this lab and [Lab 2](../../pages/day1/Lab2_processingreads.md), and start working on scRNAseq data in `R` directly with [Lab 3](../../pages/day2/Lab3_Rbioc.md).
