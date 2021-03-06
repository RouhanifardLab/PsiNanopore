# PsiNanopore
You can use this package to calculate the p-value of positions on the genome by comparing direct versus IVT files. The lowest the pvalue, the highest the chance that the position is psuedouridylated. The main input files to our tool are the aligned reads (bam files, please read the step by step guide on how to generate the bam files from "Computational pipeline of DNA/RNA sequencing analysis" section).

## System Requirement

### Hardware Requirements
package requires only a standard computer with enough RAM to support the in-memory operations.

### Software Requirements
#### OS Requirements 
This package is supported for macOS and been tested on the following system:

mcOS:Catalina(10.15.6)

#### R Dependecies
##### R 
  + R (required | Version>/3.6.1 recommended)
  + Rstudio (recommended)
##### R Packages
  + ggplot2(3.3.5), optparse(1.6.6), BiocManager(1.30.16), Biobase(2.52.0), Rsamtools(2.8.0), BSgenome(1.60.0) (**No action is required on your end. These Packages will be installed automatically on your system the first time you run the code (if not already installed).**)
## Installation Guide

First, download the package (click on the green button on top left of this page that says 'code', then click on 'Download ZIP'). Then unzip the compressed file. Next, open terminal and navigate to the directory where you've downloaded the package.

```
cd /path/to/PsiNanopore-main
```

Then use the following command to get the complete instruction on the required inputs:
```
Rscript PsiDetect.R -h
```
*Please note that the first time you run this it might take some time (couple of minutes for our computer), it's because some required packages are being installed*

*Reference genome fasta file can be downloaded from [here](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.36/)*

Here is a sample command you can use to get a list of pseudouridine candidate sites from the sample data (.bam) in "demo" folder. It takes couple of minutes (can be variable depending on the device). 
```
Rscript PsiDetect.R -f ~/Downloads/PsiNanopore-main/data/DirSeq.bam -g ~/Downloads/PsiNanopore-main/data/IVT.bam -k ~/Downloads/PsiNanopore-main/data/kmer_summary.csv -r /PATH/TO/REFERENCE/GENOME/FILE.FA -s 35599541 -e 35641526 -c chr1 -m 0.05 -o ~/Desktop/psi_candidates.csv
```
*The output is a .csv file (psi_candidates.csv in demo folder) that include the targets with the pvalue of lower that specified in -m option.*
# SignalView
*SignalView* is a visualization tool that allows you to visualize the raw ionic signals we obtain from DNA/RNA sequencing.
You can also plot signal intensity values using this tool. Let's run the code on a sample data. First navigate to the path where you've downloaded the package into:
```
cd /path/to/PsiNanopore
```
Then, run the following:
```
Rscript SignalView.R -f data/sample_001.txt  -g data/sample_002.txt -s "+" -p 75814112 -c "chr4" -o "out.pdf" 
```

## Options
  + Run ```Rscript SignalView.R -h```  in your command line to see the description of all the available options.
 
# Computational pipeline of DNA/RNA sequencing analysis
If you have not analyzed your sequencing data yet, here is a review of few tools you can use to generate the desired **txt** files. Here we break down the computational work flow required to work with SignalView to smaller steps and provide further details on each step.

## Step 0: Sequencing and basecalling

### Step 0.1: Sequencing

Once the sequencing is done, MinION outputs the raw signals in **FAST5** format. 

### Step 0.2: Basecalling

A basecaller must then be used to perform the signal segmentation, and to output the bases. Guppy is one of the best tools and can be found [here](https://github.com/nanoporetech/pyguppyclient).

## Step 1: Merging the **FASTQ** files

Once the basecalling is over, we need to merge the **FASTQ** files into a single *fq* file. We can easily do this in terminal by the following command.

```
cat path/to/fastqs/ merged_fastq.fq
```
## Step 2: Indxing the **FASTQ** files
Next step is to index the *fastq* files with their corresponding **FAST5** files. The indexing allows the tools we use in the next step to efficiently connect the basecalled sequence of each read (**FASTQ**) to their raw signal (**FAST5**) information. *Nanopolish* is a great tool for this. A complete documentation is publicly available
[here](https://github.com/jts/nanopolish), but here is the command to use:

```
nanopolish  index -d path/to/fast5s/  merged_fastq.fq
```
<a name="step3"></a>
## Step 3: Aligning the basecalled reads against a reference database
Once the indexing is done, we can proceed to the alignment step. At this step, we must align the DNA or mRNA sequences against a reference data base in **FASTA** (*fa*) format. *minimap2* is an alignment program is a great tool to align the Oxford Nanopore genomic reads to the human genome. A complete documentation for *minimap2* is publicly available [here](https://github.com/lh3/minimap2). An example command for the alignment of ONT mRNA reads is shown below. This command uses the reference (*fa*) file and the merged *fq* file (generated at Step 1), and outputs the alignment results in *sam* format.
```
minimap2 -ax splice -uf -k14 /path/to/ref_name.fa merged_fastq.fq  > alignment.sam
```
An example of the reference file can be found on the [NCBI](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.36/) website.

## Step 4: Converting **SAM** to **BAM** format, sorting and indexing the **BAM** file

Once the alignment is finished, we need to convert the Sequence Alignment Map (**SAM**) file to the Binary Alignment Map (**BAM**) file for more convenient data handling at the next steps. We use *Samtools* for this step; *Samtools* is a hybrid of programs for working with high-throughput sequencing data. This tool is also publicly available [here](https://github.com/samtools/samtools). This format conversion is done in three steps as follows.
<a name="step4.1"></a>
### Step 4.1: Converting **SAM** to **BAM** format

First, we need to convert **SAM** to **BAM** format, the following sample command uses the **SAM** file we generated in [step 3](#step3), and generates the corresponding **BAM** file.

```
samtools view -h -Sb path/to/alignment.sam  > alignment.bam
```
<a name="step4.2"></a>
### Step 4.2: Sorting  the **BAM** file

Once the **BAM** file is ready, we proceed to sorting the **BAM** file we generated at [step 4.1](#step4.1), here is a sample command:
```
samtools sort path/to/alignment.bam  -o sorted_alignment.bam
```
### Step 4.3: Indexing the sorted **BAM** file

After sorting, we then index the sorted **BAM** file we generated at [step 4.2](#step4.2), here is a sample command:
```
samtools index path/to/sorted_alignment.bam/
```

