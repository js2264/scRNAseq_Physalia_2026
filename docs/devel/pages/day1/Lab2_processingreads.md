# 4  Lab 2: From .bcl to count matrix

> **NOTE:**
>
> The estimated time for this lab is around 1h15.
>
> **We will not run `cellranger` in this lab**: it is not available in the workshop environment (see below). The two `cellranger` commands of this lab (`cellranger mkfastq` and `cellranger count`) were run beforehand by the instructors, with `cellranger` 7.1.0, on the same data. Their output is shown below, and **all the files they produced are available in your session**, in the (read-only) `~/Share/data_wrangling/` folder. All the other commands of this lab read these files: run them in your terminal!

> **TIP:**
>
> - Learn how to demultiplex sequencing data using `cellranger mkfastq`.
> - Learn how to generate gene count matrices using `cellranger count`.
> - Explore the outputs of both commands.

> **WARNING:**
>
> [`cellranger`](https://www.10xgenomics.com/support/software/cell-ranger/latest) is distributed by 10x Genomics under its own end-user license agreement, and requires a pre-built transcriptome reference (~10 GB for the mouse reference used here): neither is included in the workshop environment. To run it on your own machine, download `cellranger` from the [10x Genomics website](https://www.10xgenomics.com/support/software/cell-ranger/downloads) and a pre-built reference from [here](https://www.10xgenomics.com/support/software/cell-ranger/downloads#reference-downloads).
>
> `cellranger mkfastq` additionally relies on Illumina’s `bcl2fastq`, and recent `cellranger` releases deprecate it in favour of Illumina’s `BCL Convert`. Check the documentation of the version you use: e.g., since `cellranger` 8.0, `cellranger count` also requires a `--create-bam true` (or `false`) argument.

## 4.1 Demultiplexing sequencing data with `cellranger mkfastq`

Navigate to your terminal in `RStudio`.

Go to the [cellranger mkfastq page](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/inputs/cr-mkfastq) and read the **Overview**.

For this workshop, we are going to look at a toy dataset provided by 10x Genomics: `cellranger-tiny-bcl-1.2.0`, a sequencing run folder containing `.bcl` files. It comes with a samplesheet file, `cellranger-tiny-bcl-simple-1.2.0.csv`, which contains the sample information. Both are available in the `~/Share/data_wrangling/` folder, along with the outputs of `cellranger mkfastq` (`tiny-bcl/`) and of `cellranger count` (`counts/`):

``` sh
ls -l ~/Share/data_wrangling/
##  total 16
##  drwxr-xr-x 5 root root 4096 Nov  6 18:55 cellranger-tiny-bcl-1.2.0
##  -rw-r--r-- 1 root root   41 Nov  6 18:09 cellranger-tiny-bcl-simple-1.2.0.csv
##  drwxr-xr-x 3 root root 4096 Nov  6 18:09 counts
##  drwxr-xr-x 3 root root 4096 Nov  6 18:55 tiny-bcl
```

> **NOTE:**
>
> The run folder was downloaded from 10x Genomics (no need to run this):
>
> ``` sh
> wget https://cf.10xgenomics.com/supp/cell-exp/cellranger-tiny-bcl-1.2.0.tar.gz
> tar -xzf cellranger-tiny-bcl-1.2.0.tar.gz
> ```
>
> Note that the samplesheet currently distributed by 10x Genomics (<https://cf.10xgenomics.com/supp/cell-exp/cellranger-tiny-bcl-simple-1.2.0.csv>) lists a different sample index set (`SI-P03-C9`) from the one used to produce the outputs of this lab (`SI-GA-E3`, in `~/Share/data_wrangling/`).

> **NOTE:**
>
> Go to the Terminal tab in your `RStudio` and take a look at the 10x samplesheet `.csv` file
>
> > **NOTE:**
> >
> > ``` sh
> > cat ~/Share/data_wrangling/cellranger-tiny-bcl-simple-1.2.0.csv
> > ##  Lane,Sample,Index
> > ##  1,test_sample,SI-GA-E3
> > ```
> >
> > The reads of lane `1` tagged with one of the 4 index sequences of the 10x sample index set `SI-GA-E3` (`AGGTATTG`, `CTCCTAGT`, `TCAAGGCC`, `GATGCCAA`) belong to the sample `test_sample`.

> **NOTE:**
>
> Next, explore the contents of the sequencing directory:
>
> > **NOTE:**
> >
> > ``` sh
> > ls -l ~/Share/data_wrangling/cellranger-tiny-bcl-1.2.0
> > cat ~/Share/data_wrangling/cellranger-tiny-bcl-1.2.0/RunInfo.xml
> > ##  total 28
> > ##  drwxr-xr-x 2 root root 4096 Nov  6 18:55 Config
> > ##  drwxr-xr-x 3 root root 4096 Nov  6 18:55 Data
> > ##  drwxr-xr-x 2 root root 4096 Nov  6 18:55 InterOp
> > ##  -rw-r--r-- 1 root root   48 Nov  3  2025 RTAComplete.txt
> > ##  -rw-r--r-- 1 root root  673 Nov  3  2025 RunInfo.xml
> > ##  -rw-r--r-- 1 root root 5366 Nov  3  2025 runParameters.xml
> > ##  <?xml version="1.0"?>
> > ##  <RunInfo xmlns:xsd="http://www.w3.org/2001/XMLSchema" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" Version="2">
> > ##    <Run Id="161020_D00547_0905_BH35KCBCXY" Number="905">
> > ##      <Flowcell>H35KCBCXY</Flowcell>
> > ##      <Instrument>D00547</Instrument>
> > ##      <Date>161020</Date>
> > ##      <Reads>
> > ##        <Read Number="1" NumCycles="26" IsIndexedRead="N" />
> > ##        <Read Number="2" NumCycles="8" IsIndexedRead="Y" />
> > ##        <Read Number="3" NumCycles="98" IsIndexedRead="N" />
> > ##      </Reads>
> > ##      <FlowcellLayout LaneCount="1" SurfaceCount="1" SwathCount="1" TileCount="1" />
> > ##      <AlignToPhiX>
> > ##        <Lane>1</Lane>
> > ##      </AlignToPhiX>
> > ##    </Run>
> > ##  </RunInfo>
> > ```

The `.bcl` files are demultiplexed into `fastq` files with the `cellranger mkfastq` command. Its arguments are listed in the [cellranger mkfastq documentation](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/inputs/cr-mkfastq) (or with `cellranger mkfastq --help`, if `cellranger` is installed).

> **NOTE:**
>
> What is the purpose of `--id`, `--run`, and `--csv` arguments ?
>
> Choose these arguments wisely to demultiplex the sequencing data.
>
> > **NOTE:**
> >
> > - `--id`: the name of the pipeline run, and of the folder created to store its outputs;
> > - `--run`: the path to the Illumina run folder, containing the `.bcl` files;
> > - `--csv`: the path to the samplesheet, assigning sample index(es) to each sample, in each lane.
> >
> > This is the command that was run by the instructors, from the folder containing the run folder and the samplesheet:
> >
> > ``` sh
> > ## Not run: `cellranger` is not available in your session, and its outputs
> > ## are already in ~/Share/data_wrangling/tiny-bcl/
> > cellranger mkfastq \
> >     --id tiny-bcl \
> >     --run cellranger-tiny-bcl-1.2.0 \
> >     --csv cellranger-tiny-bcl-simple-1.2.0.csv
> > ```
> >
> > And here is its output (abridged, and with the paths edited to the location of the files in your session):
> >
> > ``` text
> > cellranger mkfastq (cellranger-7.1.0)
> > Copyright (c) 2021 10x Genomics, Inc.  All rights reserved.
> > -------------------------------------------------------------------------------
> >
> > Martian Runtime - v4.0.10
> > Checking run folder...
> > Checking RunInfo.xml...
> > Checking system environment...
> > Emitting run information...
> > Checking read specification...
> > Checking samplesheet specs...
> > 2025-11-13 18:54:42 [runtime] (ready)           ID.tiny-bcl.MAKE_FASTQS_CS.MAKE_FASTQS.PREPARE_SAMPLESHEET
> > 2025-11-13 18:54:42 [runtime] (run:local)       ID.tiny-bcl.MAKE_FASTQS_CS.MAKE_FASTQS.PREPARE_SAMPLESHEET.fork0.chnk0.main
> > 2025-11-13 18:54:44 [runtime] (chunks_complete) ID.tiny-bcl.MAKE_FASTQS_CS.MAKE_FASTQS.PREPARE_SAMPLESHEET
> > 2025-11-13 18:54:44 [runtime] (ready)           ID.tiny-bcl.MAKE_FASTQS_CS.MAKE_FASTQS.BCL2FASTQ_WITH_SAMPLESHEET
> > ...
> > 2025-11-13 18:54:56 [runtime] (join_complete)   ID.tiny-bcl.MAKE_FASTQS_CS.MAKE_FASTQS.MERGE_FASTQS_BY_LANE_SAMPLE
> >
> > Outputs:
> > - FASTQ output folder:   ~/Share/data_wrangling/tiny-bcl/outs/fastq_path
> > - Interop output folder: ~/Share/data_wrangling/tiny-bcl/outs/interop_path
> > - Input samplesheet:     ~/Share/data_wrangling/tiny-bcl/outs/input_samplesheet.csv
> >
> > Waiting 6 seconds for UI to do final refresh.
> > Pipestance completed successfully!
> >
> > Saving pipestance info to "tiny-bcl/tiny-bcl.mri.tgz"
> > ```

The output folders can be viewed by running the `ls` command:

``` sh
ls -l ~/Share/data_wrangling/tiny-bcl/outs/fastq_path
ls -l ~/Share/data_wrangling/tiny-bcl/outs/fastq_path/H35KCBCXY/test_sample
##  total 12
##  drwxr-xr-x 3 root root 4096 Nov  6 18:55 H35KCBCXY
##  drwxr-xr-x 3 root root 4096 Nov  3  2025 Reports
##  drwxr-xr-x 2 root root 4096 Nov  6 18:55 Stats
##  total 48148
##  -rw-r--r-- 1 root root  3930346 Nov  3  2025 test_sample_S1_L001_I1_001.fastq.gz
##  -rw-r--r-- 1 root root 11093755 Nov  3  2025 test_sample_S1_L001_R1_001.fastq.gz
##  -rw-r--r-- 1 root root 34273874 Nov  3  2025 test_sample_S1_L001_R2_001.fastq.gz
```

`cellranger mkfastq` also writes the reads which could not be assigned to any sample in `Undetermined_*.fastq.gz` files: these were removed from `~/Share/data_wrangling/tiny-bcl/` to save space.

> **NOTE:**
>
> Look at the first reads listed in index read (`I1`), read 1 (`R1`), and read 2 (`R2`) files.
>
> What is the purpose of the index read? What does it contain? And read 1?
>
> > **NOTE:**
> >
> > ``` sh
> > cd ~/Share/data_wrangling/tiny-bcl/outs/fastq_path/H35KCBCXY/test_sample
> > zcat test_sample_S1_L001_I1_001.fastq.gz | head -n 12
> > zcat test_sample_S1_L001_R1_001.fastq.gz | head -n 4
> > zcat test_sample_S1_L001_R2_001.fastq.gz | head -n 4
> > ##  @D00547:905:H35KCBCXY:1:1101:1284:1973 1:N:0:AGGTATTG
> > ##  AGGTATTG
> > ##  +
> > ##  GGGGGIII
> > ##  @D00547:905:H35KCBCXY:1:1101:1729:1981 1:N:0:AGGTATTG
> > ##  AGGTATTG
> > ##  +
> > ##  GAGGGIIG
> > ##  @D00547:905:H35KCBCXY:1:1101:2005:1944 1:N:0:AGGTATTG
> > ##  AGGTATTG
> > ##  +
> > ##  AA<G.A.G
> > ##  @D00547:905:H35KCBCXY:1:1101:1284:1973 1:N:0:AGGTATTG
> > ##  CGGGTCACAAGCGAGTGTTTTTAGTN
> > ##  +
> > ##  GGGGGIIIIIIIIGIGGGIIIIIII#
> > ##  @D00547:905:H35KCBCXY:1:1101:1284:1973 2:N:0:AGGTATTG
> > ##  TCAGCCCTAACCGCTACTGGCTGTGTGCTGCCACAGGCCCCAGCATCAAGATCTGGGATTTAGAGGGAAAGATCATTGTAGATGAACTGAAGCAAGAN
> > ##  +
> > ##  AGGGGIIIGIIIIIGGIIIIIGIGGIIIIIGGIIIIIIIIIGGIIIIGIGIIIGIGGGGGIIIGGGIIIIIIIIIGGGIGIGIGIIIIIIIIIIIIG#
> > ```
> >
> > The index read (`I1`) contains the sample index (here, `AGGTATTG`, one of the 4 sequences of the `SI-GA-E3` index set): it identifies the library, i.e. the sample, each read comes from. Read 1 (`R1`, 26 bp here) contains the cell barcode (first 16 bp) and the UMI (next 10 bp). Read 2 (`R2`) contains the sequence of the transcript (cDNA).

Open the html file `~/Share/data_wrangling/tiny-bcl/outs/fastq_path/Reports/html/index.html` by navigating to the file in RStudio, using the Files Tab. When you click on the file, select the option to View in Web Browser. Take some time to explore the demultiplexed outputs.

## 4.2 Generating gene count matrices with `cellranger count`

Go to the [cellranger count algorithm overview](https://www.10xgenomics.com/support/software/cell-ranger/latest/algorithms-overview/cr-gex-algorithm) and read the section on **Alignment** (Read Trimming, Genome Alignment, MAPQ adjustment, Transcriptome Alignment, UMI Counting).

The `cellranger count` command generates gene count matrices from the demultiplexed reads. Its arguments are described in the [cellranger documentation](https://www.10xgenomics.com/support/software/cell-ranger/latest) (or with `cellranger count --help`, if `cellranger` is installed).

The gene count matrix of the `test_sample` reads was generated with `cellranger count` and the `refdata-gex-mm10-2020-A` mouse transcriptome reference provided by 10x Genomics.

> **NOTE:**
>
> What is the purpose of the `--id`, `--transcriptome`, `--fastqs` and `--sample` arguments?
>
> > **NOTE:**
> >
> > - `--id`: the name of the pipeline run, and of the folder created to store its outputs;
> > - `--transcriptome`: the path to the `cellranger`-compatible transcriptome reference (genome index and gene annotations);
> > - `--fastqs`: the path to the folder containing the `fastq` files;
> > - `--sample`: the prefix of the names of the `fastq` files to use (`test_sample` for `test_sample_S1_L001_R1_001.fastq.gz`, …).
> >
> > This is the command that was run by the instructors (from the same folder as `cellranger mkfastq`):
> >
> > ``` sh
> > ## Not run: `cellranger` is not available in your session, and its outputs
> > ## are already in ~/Share/data_wrangling/counts/
> > cellranger count \
> >     --id counts \
> >     --transcriptome refdata-gex-mm10-2020-A \
> >     --fastqs tiny-bcl/outs/fastq_path/H35KCBCXY/test_sample \
> >     --sample test_sample
> > ```
> >
> > And here is its output (abridged, and with the paths edited to the location of the files in your session):
> >
> > ``` text
> > Martian Runtime - v4.0.10
> > Checking sample info...
> > Checking FASTQ folder...
> > Checking reference...
> > Checking reference_path (refdata-gex-mm10-2020-A)...
> > Checking optional arguments...
> > mro: v4.0.10
> > mrp: v4.0.10
> > Anaconda: Python 3.9.7
> > numpy: 1.22.3
> > scipy: 1.9.1
> > pysam: 0.17.0
> > h5py: 3.6.0
> > pandas: 1.3.5
> > STAR: 2.7.2a
> > samtools: samtools 1.12
> > 2025-11-13 19:17:53 [runtime] (ready)           ID.counts.SC_RNA_COUNTER_CS.WRITE_GENE_INDEX
> > 2025-11-13 19:17:53 [runtime] (run:local)       ID.counts.SC_RNA_COUNTER_CS.WRITE_GENE_INDEX.fork0.chnk0.main
> > ...
> >
> > Outputs:
> > - Run summary HTML:                         ~/Share/data_wrangling/counts/outs/web_summary.html
> > - Run summary CSV:                          ~/Share/data_wrangling/counts/outs/metrics_summary.csv
> > - BAM:                                      ~/Share/data_wrangling/counts/outs/possorted_genome_bam.bam
> > - BAM BAI index:                            ~/Share/data_wrangling/counts/outs/possorted_genome_bam.bam.bai
> > - BAM CSI index:                            null
> > - Filtered feature-barcode matrices MEX:    ~/Share/data_wrangling/counts/outs/filtered_feature_bc_matrix
> > - Filtered feature-barcode matrices HDF5:   ~/Share/data_wrangling/counts/outs/filtered_feature_bc_matrix.h5
> > - Unfiltered feature-barcode matrices MEX:  ~/Share/data_wrangling/counts/outs/raw_feature_bc_matrix
> > - Unfiltered feature-barcode matrices HDF5: ~/Share/data_wrangling/counts/outs/raw_feature_bc_matrix.h5
> > - Secondary analysis output CSV:            ~/Share/data_wrangling/counts/outs/analysis
> > - Per-molecule read information:            ~/Share/data_wrangling/counts/outs/molecule_info.h5
> > - CRISPR-specific analysis:                 null
> > - Antibody aggregate barcodes:              null
> > - Loupe Browser file:                       ~/Share/data_wrangling/counts/outs/cloupe.cloupe
> > - Feature Reference:                        null
> > - Target Panel File:                        null
> > - Probe Set File:                           null
> >
> > Waiting 6 seconds for UI to do final refresh.
> > Pipestance completed successfully!
> > ```
> >
> > With this version of `cellranger`, reads mapped to introns are also counted by default (`--include-introns true`), and a `BAM` file is created by default.

Read about the [format of the feature-barcode matrices](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/outputs/cr-outputs-mex-matrices).

The pipeline outputs can be viewed with `ls`:

``` sh
ls -l ~/Share/data_wrangling/counts/outs
##  total 56900
##  drwxr-xr-x 7 root root     4096 Nov  6 18:55 analysis
##  -rw-r--r-- 1 root root  2782024 Nov  3  2025 cloupe.cloupe
##  drwxr-xr-x 2 root root     4096 Nov  6 18:55 filtered_feature_bc_matrix
##  -rw-r--r-- 1 root root   506767 Nov  3  2025 filtered_feature_bc_matrix.h5
##  -rw-r--r-- 1 root root      627 Nov  3  2025 metrics_summary.csv
##  -rw-r--r-- 1 root root  1806311 Nov  3  2025 molecule_info.h5
##  -rw-r--r-- 1 root root 47960294 Nov  3  2025 possorted_genome_bam.bam
##  -rw-r--r-- 1 root root  1744976 Nov  3  2025 possorted_genome_bam.bam.bai
##  drwxr-xr-x 2 root root     4096 Nov  6 18:55 raw_feature_bc_matrix
##  -rw-r--r-- 1 root root   921409 Nov  3  2025 raw_feature_bc_matrix.h5
##  -rw-r--r-- 1 root root  2513806 Nov  3  2025 web_summary.html
```

> **NOTE:**
>
> Can you locate the feature-barcode matrices? What is the difference between the `raw_feature_bc_matrix` and `filtered_feature_bc_matrix` data types?
>
> > **TIP:**
> >
> > ``` sh
> > cd ~/Share/data_wrangling/counts/outs
> > ls -l raw_feature_bc_matrix.h5 filtered_feature_bc_matrix.h5
> > ls -l raw_feature_bc_matrix/ filtered_feature_bc_matrix/
> > ##  -rw-r--r-- 1 root root 506767 Nov  3  2025 filtered_feature_bc_matrix.h5
> > ##  -rw-r--r-- 1 root root 921409 Nov  3  2025 raw_feature_bc_matrix.h5
> > ##  filtered_feature_bc_matrix/:
> > ##  total 748
> > ##  -rw-r--r-- 1 root root   3981 Nov  3  2025 barcodes.tsv.gz
> > ##  -rw-r--r-- 1 root root 290896 Nov  3  2025 features.tsv.gz
> > ##  -rw-r--r-- 1 root root 465525 Nov  3  2025 matrix.mtx.gz
> > ##  
> > ##  raw_feature_bc_matrix/:
> > ##  total 1100
> > ##  -rw-r--r-- 1 root root 186199 Nov  3  2025 barcodes.tsv.gz
> > ##  -rw-r--r-- 1 root root 260194 Nov  3  2025 features.tsv.gz
> > ##  -rw-r--r-- 1 root root 674205 Nov  3  2025 matrix.mtx.gz
> > ```
>
> How many cells and genes are stored in each matrix? Comment.
>
> > **NOTE:**
> >
> > ``` sh
> > # Number of genes
> > zcat ~/Share/data_wrangling/counts/outs/filtered_feature_bc_matrix/features.tsv.gz | wc -l
> > zcat ~/Share/data_wrangling/counts/outs/raw_feature_bc_matrix/features.tsv.gz | wc -l
> >
> > # Number of cells
> > zcat ~/Share/data_wrangling/counts/outs/filtered_feature_bc_matrix/barcodes.tsv.gz | wc -l
> > zcat ~/Share/data_wrangling/counts/outs/raw_feature_bc_matrix/barcodes.tsv.gz | wc -l
> > ##  32285
> > ##  32285
> > ##  686
> > ##  42543
> > ```
> >
> > Both matrices list the same genes. The raw matrix contains every barcode with at least one UMI (i.e. mostly empty droplets), while the filtered matrix only contains the barcodes `cellranger` called as cells.
>
> Now open the html file `~/Share/data_wrangling/counts/outs/web_summary.html` by navigating to the file in RStudio, using the Files Tab. When you click on the file, select the option to View in Web Browser. Take some time to explore the gene expression matrix outputs.

> **NOTE:**
>
> `cellranger count` also produced a `BAM` file, containing the alignment of every read. Use `samtools` (see [Lab 1](../../pages/day1/Lab1_Intro.md)) to look at one of the reads assigned to a gene. Where are the cell barcode, the UMI and the gene of the read stored?
>
> > **NOTE:**
> >
> > ``` sh
> > samtools view -d GN ~/Share/data_wrangling/counts/outs/possorted_genome_bam.bam | head -n 1 | tr '\t' '\n'
> > ##  D00547:905:H35KCBCXY:1:1101:2865:43071
> > ##  16
> > ##  chr1
> > ##  4774072
> > ##  255
> > ##  98M
> > ##  *
> > ##  0
> > ##  0
> > ##  CAGCTCCAAGCATGTAACGAGCAAGCTAGAGTGTAGTGGTCTGCCTACAGCAAGTTATCACGAAAGACAACCGTCATGGGACTCACACTAACATCAGA
> > ##  GGGGAGIGGGIIIIIIIIIIIIIIIGGIIIIIGIIIIIIGIIIIGIIIIIIIIGGAGGGIIIIIIIIGGGGIGGGG.IIGIIGIGGGIIIGGGAGGGG
> > ##  NH:i:1
> > ##  HI:i:1
> > ##  AS:i:96
> > ##  nM:i:0
> > ##  RG:Z:counts:0:1:H35KCBCXY:1
> > ##  TX:Z:ENSMUST00000130201,+930,98M;ENSMUST00000156816,+3244,98M
> > ##  GX:Z:ENSMUSG00000033845
> > ##  GN:Z:Mrpl15
> > ##  fx:Z:ENSMUSG00000033845
> > ##  RE:A:E
> > ##  xf:i:25
> > ##  CR:Z:GTTAAGCCAGCCTGTG
> > ##  CY:Z:GGGGGIIIIIIIIIII
> > ##  CB:Z:GTTAAGCCAGCCTGTG-1
> > ##  UR:Z:ATCGGTCAGC
> > ##  UY:Z:IIIIIIIIII
> > ##  UB:Z:ATCGGTCAGC
> > ```
> >
> > After the 11 mandatory fields of each alignment, `cellranger` adds tags to each read: `CR`/`CB` are the raw/corrected cell barcode, `UR`/`UB` the raw/corrected UMI, and `GX`/`GN` the ID/name of the gene the read is assigned to. Their description is available in the [cellranger documentation](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/outputs/cr-outputs-bam).
