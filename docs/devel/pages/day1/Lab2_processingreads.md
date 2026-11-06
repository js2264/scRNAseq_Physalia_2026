# 4  Lab 2: From .bcl to count matrix

> **NOTE:**
>
> The estimated time for this lab is around 1h15.
>
> This lab only uses command-line tools (`cellranger`) from the `terminal`: none of the code on this page is run when the book is built. The outputs shown below were obtained by the instructors when preparing the lab.

> **TIP:**
>
> - Learn how to demultiplex sequencing data using `cellranger mkfastq`.
> - Learn how to generate gene count matrices using `cellranger count`.

> **WARNING:**
>
> [`cellranger`](https://www.10xgenomics.com/support/software/cell-ranger/latest) is distributed by 10x Genomics under its own end-user license agreement, and `cellranger mkfastq` additionally relies on Illumina’s `bcl2fastq`. Your instructors will tell you where to find both tools, as well as the pre-built `cellranger` transcriptome reference used in this lab (~10 GB), in your workshop session. On your own machine, `cellranger` can be downloaded from the [10x Genomics website](https://www.10xgenomics.com/support/software/cell-ranger/downloads) and pre-built references are listed [here](https://www.10xgenomics.com/support/software/cell-ranger/downloads#reference-downloads).
>
> Recent `cellranger` releases deprecate `cellranger mkfastq` in favour of Illumina’s `BCL Convert`: check the documentation of the version you use.

## 4.1 Demultiplexing sequencing data with `cellranger mkfastq`

Navigate to your terminal in `RStudio`.

Go to the [cellranger mkfastq page](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/inputs/cr-mkfastq) and read the **Overview**.

For this workshop, we are going to look at a toy dataset provided by 10x Genomics: `cellranger-tiny-bcl-1.2.0`. It comes with a samplesheet file, `cellranger-tiny-bcl-simple-1.2.0.csv`, which contains the sample information.

``` sh
mkdir -p ~/data/data_wrangling && cd ~/data/data_wrangling
wget https://cf.10xgenomics.com/supp/cell-exp/cellranger-tiny-bcl-1.2.0.tar.gz
wget https://cf.10xgenomics.com/supp/cell-exp/cellranger-tiny-bcl-simple-1.2.0.csv
tar -xzf cellranger-tiny-bcl-1.2.0.tar.gz
```

> **NOTE:**
>
> Go to the Terminal tab in your `RStudio` and take a look at the 10x samplesheet `.csv` file
>
> > **NOTE:**
> >
> > ``` sh
> > cat ~/data/data_wrangling/cellranger-tiny-bcl-simple-1.2.0.csv
> > ```
> >
> > ``` text
> > Lane,Sample,Index
> > 1,test_sample,SI-P03-C9
> > ```

> **NOTE:**
>
> Next, explore the contents of the sequencing directory:
>
> > **NOTE:**
> >
> > ``` sh
> > ls -l ~/data/data_wrangling/cellranger-tiny-bcl-1.2.0
> > ```
> >
> > ``` text
> > total 28
> > drwxr-xr-x 2 rstudio rstudio 4096 Nov 12  2016 Config
> > drwxr-xr-x 3 rstudio rstudio 4096 Nov 12  2016 Data
> > drwxr-xr-x 2 rstudio rstudio 4096 Nov 12  2016 InterOp
> > -rw-r--r-- 1 rstudio rstudio   48 Nov 12  2016 RTAComplete.txt
> > -rw-r--r-- 1 rstudio rstudio  673 Nov 12  2016 RunInfo.xml
> > -rw-r--r-- 1 rstudio rstudio 5366 Nov 12  2016 runParameters.xml
> > ```

Now we can demultiplex our bcl files by running the following command in the terminal. To do that, we will use the `cellranger mkfastq` command.

``` sh
cellranger mkfastq --help
```

> **NOTE:**
>
> What is the purpose of `--id`, `--run`, and `--csv` arguments ?
>
> Choose these arguments wisely and execute the command to demultiplex the sequencing data.
>
> > **NOTE:**
> >
> > ``` sh
> > cd ~/data/data_wrangling
> > cellranger mkfastq \
> >     --id tiny-bcl \
> >     --run ~/data/data_wrangling/cellranger-tiny-bcl-1.2.0 \
> >     --csv ~/data/data_wrangling/cellranger-tiny-bcl-simple-1.2.0.csv
> > ```
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
> >
> > ...
> >
> > Outputs:
> > - FASTQ output folder:   ~/data/data_wrangling/tiny-bcl/outs/fastq_path
> > - Interop output folder: ~/data/data_wrangling/tiny-bcl/outs/interop_path
> > - Input samplesheet:     ~/data/data_wrangling/tiny-bcl/outs/input_samplesheet.csv
> >
> > Waiting 6 seconds for UI to do final refresh.
> > Pipestance completed successfully!
> > ```

The output folders can be viewed by running the `ls` command:

``` sh
ls -l ~/data/data_wrangling/tiny-bcl/outs/fastq_path/H35KCBCXY/test_sample
```

``` text
-rw-r--r-- 1 rstudio rstudio test_sample_S1_L001_I1_001.fastq.gz
-rw-r--r-- 1 rstudio rstudio test_sample_S1_L001_R1_001.fastq.gz
-rw-r--r-- 1 rstudio rstudio test_sample_S1_L001_R2_001.fastq.gz
```

> **NOTE:**
>
> Look at the first reads listed in index read (`I1`), read 1 (`R1`), and read (`R2`) files.
>
> What is the purpose of the index read? What does it contain?
>
> > **NOTE:**
> >
> > ``` sh
> > zcat ~/data/data_wrangling/tiny-bcl/outs/fastq_path/H35KCBCXY/test_sample/test_sample_S1_L001_I1_001.fastq.gz | head -n 12
> > ```
> >
> > ``` text
> > @D00547:905:H35KCBCXY:1:1101:1284:1973 1:N:0:AGGTATTG
> > AGGTATTG
> > +
> > GGGGGIII
> > @D00547:905:H35KCBCXY:1:1101:1729:1981 1:N:0:AGGTATTG
> > AGGTATTG
> > +
> > GAGGGIIG
> > @D00547:905:H35KCBCXY:1:1101:2005:1944 1:N:0:AGGTATTG
> > AGGTATTG
> > +
> > AA<G.A.G
> > ```

Open the html file `tiny-bcl/outs/fastq_path/Reports/html/index.html` by navigating to the file in RStudio, using the Files Tab. When you click on the file, select the option to View in Web Browser. Take some time to explore the demultiplexed outputs.

## 4.2 Generating gene count matrices with `cellranger count`

Go to the [cellranger count algorithm overview](https://www.10xgenomics.com/support/software/cell-ranger/latest/algorithms-overview/cr-gex-algorithm) and read the section on **Alignment** (Read Trimming, Genome Alignment, MAPQ adjustment, Transcriptome Alignment, UMI Counting).

In bash, the `cellranger count` command is used to generate gene count matrices from the demultiplexed reads.

``` sh
cellranger count --help
```

Now, we will generate a gene count matrix for the `test_sample` reads using `cellranger count` and the `refdata-gex-mm10-2020-A` transcriptome reference provided by 10x Genomics.

> **NOTE:**
>
> What is the purpose of the `--id`, `--transcriptome`, `--fastqs`, `--sample`, and `--create-bam` arguments?
>
> > **NOTE:**
> >
> > ``` sh
> > ## Change this path to the location of the pre-built reference in your session
> > REF=~/data/refdata-gex-mm10-2020-A
> > cd ~/data/data_wrangling
> > cellranger count \
> >     --id counts \
> >     --transcriptome ${REF} \
> >     --fastqs tiny-bcl/outs/fastq_path/H35KCBCXY/test_sample \
> >     --sample test_sample \
> >     --create-bam true
> > ```
> >
> > ``` text
> > Martian Runtime - v4.0.10
> > Checking sample info...
> > Checking FASTQ folder...
> > Checking reference...
> > Checking optional arguments...
> >
> > ...
> >
> > Outputs:
> > - Run summary HTML:                         ~/data/data_wrangling/counts/outs/web_summary.html
> > - Run summary CSV:                          ~/data/data_wrangling/counts/outs/metrics_summary.csv
> > - BAM:                                      ~/data/data_wrangling/counts/outs/possorted_genome_bam.bam
> > - BAM BAI index:                            ~/data/data_wrangling/counts/outs/possorted_genome_bam.bam.bai
> > - Filtered feature-barcode matrices MEX:    ~/data/data_wrangling/counts/outs/filtered_feature_bc_matrix
> > - Filtered feature-barcode matrices HDF5:   ~/data/data_wrangling/counts/outs/filtered_feature_bc_matrix.h5
> > - Unfiltered feature-barcode matrices MEX:  ~/data/data_wrangling/counts/outs/raw_feature_bc_matrix
> > - Unfiltered feature-barcode matrices HDF5: ~/data/data_wrangling/counts/outs/raw_feature_bc_matrix.h5
> > - Secondary analysis output CSV:            ~/data/data_wrangling/counts/outs/analysis
> > - Per-molecule read information:            ~/data/data_wrangling/counts/outs/molecule_info.h5
> > - Loupe Browser file:                       ~/data/data_wrangling/counts/outs/cloupe.cloupe
> >
> > Waiting 6 seconds for UI to do final refresh.
> > Pipestance completed successfully!
> > ```

While the count command is running, read about the [format of the feature-barcode matrices](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/outputs/cr-outputs-mex-matrices).

Once the count command is finished running, the pipeline outputs can be viewed with `ls`:

``` sh
ls -l ~/data/data_wrangling/counts/outs
```

> **NOTE:**
>
> Can you locate the feature-barcode matrices? What is the difference between the `raw_feature_bc_matrix` and `filtered_feature_bc_matrix` data types?
>
> > **TIP:**
> >
> > ``` sh
> > ls -l ~/data/data_wrangling/counts/outs/raw_feature_bc_matrix.h5
> > ls -l ~/data/data_wrangling/counts/outs/filtered_feature_bc_matrix.h5
> > ```
>
> How many cells and genes are stored in each matrix? Comment.
>
> > **NOTE:**
> >
> > ``` sh
> > # Number of genes
> > zcat ~/data/data_wrangling/counts/outs/filtered_feature_bc_matrix/features.tsv.gz | wc -l
> > zcat ~/data/data_wrangling/counts/outs/raw_feature_bc_matrix/features.tsv.gz | wc -l
> >
> > # Number of cells
> > zcat ~/data/data_wrangling/counts/outs/filtered_feature_bc_matrix/barcodes.tsv.gz | wc -l
> > zcat ~/data/data_wrangling/counts/outs/raw_feature_bc_matrix/barcodes.tsv.gz | wc -l
> > ```
>
> Now open the html file `counts/outs/web_summary.html` by navigating to the file in RStudio, using the Files Tab. When you click on the file, select the option to View in Web Browser. Take some time to explore the gene expression matrix outputs.
