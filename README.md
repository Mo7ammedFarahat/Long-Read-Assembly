## BYOD workshop 21th-25th Oct 2024, Cape Town - South Africa

Note: Material used to prepare for the workshop was extracted from https://github.com/human-pangenomics/hprc-tutorials/tree/GA-workshop and https://genomicsaotearoa.github.io/long-read-assembly/

Customized by: **Mohammed Farahat**, to be up and running on University of Cape Town HPC (ilifu)

## Table of Contents

* [Ilifu HPC Configuration](#ilifu-hpc-configuration)
     * [Ilifu modules](#ilifu-modules)
     * [Conda environment](#conda-environment)
* [Dataset](#dataset)
* [Part 1: Familiarize Ourselves With The Data](#part-1-familiarize-ourselves-with-the-data)
    * [Subset Our Input Data](#subset-our-input-data)
    * [Data Comparison](#Now-lets-compare-the-data)
    * [Cleaning Data For Assembly](#cleaning-data-for-assembly)
        * [PacBio Adapter Trimming](#pacBio-adapter-trimming)
        * [ONT Read Length Filtering](#ont-read-length-filtering)
    * [Phasing Data](#phasing-data)
        * [Trio data: Meryl](#trio-data-meryl)
        * [Trio data: Yak](#trio-data-yak)

## Ilifu HPC Configuration

Ilifu is a regional node, known as a Tier II node, in a national infrastructure, and partly funded by the Department of Science and Technology (DST) through their Data-Intensive Research Initiative of South Africa (DIRISA).

https://www.ilifu.ac.za/about/


Ilifu is using Slurm as a job scheduling system.


### Ilifu modules
Some of the tools required for this workshop are already installed as modules on ilifu, to display the available modules use `module avail` and to load a module use `module load tool/version`.

### Conda environment
For the tools that are not preinstalled on ilifu, I have created a conda env that contains all the required tools and packages you will need. 

Download the yml file from here [refgraph.yml](https://github.com/Mo7ammedFarahat/Long-Read-Assembly/blob/main/refgraph.yml)  
Or use it directly from the workshop project dir `/cbio/projects/037/mohammed/condaEnv/refgraph.yml`

To create the environment:
```bash
username@slurm-login:~$ conda env create -f refgraph.yml
username@slurm-login:~$ conda activate refgraph
```
Now, you can check one of the tools that we will use later, like `mashmap`:
```bash
username@slurm-login:~$ mashmap --version
3.1.3
```
---
## Dataset

  In this workshop we will be using data from HG002, which is a reference sample from the [Genome In A Bottle (GIAB)](https://www.nist.gov/programs-projects/genome-bottle) consortium. The GIAB project releases benchmark data for genomic characterization, and you may have seen their benchmark variant calls and regions out in the wild. As part of their benchmarking material generation, they release datasets for their reference samples. We will be using those in this workshop.
    
    
    HG002 is actually part of a trio of reference samples. Below is the family listing, also known as the Ashkenazim trio:
    
    * HG002: Son
    * HG003: Father
    * HG004: Mother


You will find this dataset here:
`/cbio/projects/037/mohammed/T2T/HG002`

## Part 1: Familiarize Ourselves With The Data
**Create a directoty**
```bash
username@slurm-login:~$ mkdir ~/lra
username@slurm-login:~$ cd ~/lra
username@slurm-login:~$ mkdir day1_data
username@slurm-login:~$ cd day1_data
```
### Subset Our Input Data
In order to get a feel for the data, we only need a small portion of it. Pull the first few thousand reads of the **HiFi reads** and write them to new files.
```bash
zcat /cbio/projects/037/mohammed/T2T/HG002/LRA/resources/deepconsensus/m64011_190830_220126.Q20.fastq.gz \
    | head -n 200000 \
    | pigz > hifi_50k_reads.fq.gz &
```
Next, downsample the **ONT UL reads**, too.

```bash
samtools fastq -@4 \
    /cbio/projects/037/mohammed/T2T/HG002/LRA/resources/ont_ul/03_08_22_R941_HG002_1_Guppy_6.1.2_5mc_cg_prom_sup.bam \
    | head -n 20000 \
    | pigz > ont_ul_5k_reads.fq.gz &

```
### Now let's compare the data  
We are going to use a tool called NanoComp. This tool can take in multiple FASTQs (or BAMs) and will create summary statistics and nice plots that show things like read length and quality scores. NanoComp has nano in the name, and has some ONT-specific functionality, but it can be used with PacBio data just fine.

```bash
NanoComp --fastq \
    hifi_50k_reads.fq.gz \
    ont_ul_5k_reads.fq.gz \
    --names PacBio_HiFi ONT_UL \
    --outdir nanocomp_hifi_vs_ul
```
Once the run is complete (~2 minutes), navigate in your file browser to the folder that NanoComp just created and then click on the `NanoComp-report.html`

### Cleaning Data For Assembly
#### PacBio Adapter Trimming
PacBio's CCS software attempts to identify adapters and remove them. This process is getting better all the time, but some datasets (especially older ones) can have adapters remaining. If this is the case, adapters can find their way into the assemblies. 

Run CutAdapt to check for adapter sequences in the downsampled data that we are currently using. (The results will print to stdout on your terminal screen.)


```bash
module load cutadapt

cutadapt \
    -b "AAAAAAAAAAAAAAAAAATTAACGGAGGAGGAGGA;min_overlap=35" \
    -b "ATCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGAT;min_overlap=45" \
    --discard-trimmed \
    -o /dev/null \
    hifi_50k_reads.fq.gz \
    -j 0 \
    --revcomp \
    -e 0.05
```
Notice that we are writing output to `/dev/null`. We are working on a subset of these reads so the runtime is reasonable. There is no need to hold onto the reads that we are filtering on, it is just a subset of the data.

#### ONT Read Length Filtering
Hifiasm is often run with ONT data filtered to be over 50kb in length, so let's filter that data now to see how much of the data remains. 
```bash
module load seqkit/2.6.0
seqkit seq \
    -m 50000 \
    ont_ul_5k_reads.fq.gz \
    | pigz > ont_ul_5k_reads.50kb.fq.gz &
```
Now we can quickly check how many reads are retained.
```bash
NanoComp --fastq \
    ont_ul_5k_reads.fq.gz \
    ont_ul_5k_reads.50kb.fq.gz \
    --names ONT_UL ONT_UL_50kb \
    --outdir nanocomp_ul_vs_ul_50kb
```
### Phasing Data
Now that we've introduced the data that creates the graphs, it's time to talk about data types that can phase them in order to produce fully phased diploid assemblies.

At the moment the easiest and most effective way to phase human assemblies is with trio information. Meaning you sequence a sample, and then you also sequence its parents. You then look at which parts of the genome the sample inherited from one parent and not the other. This is done with kmer databases (DBs). In our case, we will use both Meryl (for Verkko) and yak (for hifiasm) so let's take a moment to learn about kmer DBs.

#### Trio data: Meryl
[Meryl](https://github.com/marbl/meryl) is a kmer counter that dates back to Celera. It creates kmer DBs, but it is also a toolset that you can use for finding kmers and manipulating kmer count sets. Meryl is to kmers what BedTools is to genomic regions.

Today we want to use Meryl in the context of creating databases from PCR-free Illumina readsets. These can be used both during the assembly process and during the post-assembly QC. 

**Make sure you are in the right directory**

```bash
cd day1_data
```

**Now create a small file to work with**

```bash
zcat /cbio/projects/037/mohammed/T2T/HG002/LRA/resources/ilmn/pat/HG003_HiSeq30x_subsampled_R1.fastq.gz \
    | head -n 20000000 \
    | pigz > HG003_HiSeq30x_5M_reads_R1.fastq.gz &
```    

**Create a *k*-mer DB from an Illumina read set**


```bash
module load merqury/1.3
module load meryl/1.4.1

meryl count \
    compress \
    k=30 \
    threads=4 \
    memory=8 \
    HG003_HiSeq30x_5M_reads_R1.fastq.gz \
    output paternal_5M_compress.k30.meryl
```
This should be pretty fast because we are just using a small amount of data to get a feel for the program. The output of Meryl is a folder that contains 64 index files and 64 data files. If you try and look at the data files you'll see that they aren't human readable. In order to look at the actual k-mers, you have to use meryl to print them.

**Look at the *k*-mers**


```bash
module load meryl/1.4.1

meryl print \
    greater-than 1 \
    paternal_5M_compress.k30.meryl \
    | head

```
The first column is the *k*-mer and the second column is the count of that *k*-mer in the dataset.

**Take a look at some statistics for the DB**


```bash
module load meryl/1.4.1

meryl statistics \
    paternal_5M_compress.k30.meryl \
    | head -n 20

```

We see a lot of *k*-mers missing and the histogram (frequency column) has a ton of counts at 1. This makes sense for a heavily downsampled dataset. Great. We just got a feel for how to use Meryl in general on subset data. Now let's actually take a look at how to create Meryl DBs for Verkko assemblies.
<details>
<summary>clipboard-question "How would we run Meryl for Verkko?"</summary>

    **Here is what the Slurm script would look like:**
    
    (Don't run this, it is slow! We have made these for you already.)

    
```bash
        #!/bin/bash 

        #SBATCH --job-name='meryl_run'
        #SBATCH --cpus-per-task=32
        #SBATCH --time=12:00:00
        #SBATCH --mem=96G
        #SBATCH --partition=Main
        #SBATCH --output=testjob-%j-stdout.log
        #SBATCH --error=testjob-%j-stderr.log


        module load merqury/1.3
        module load meryl/1.4.1
        export MERQURY=$(dirname $(which merqury.sh))

        ## Create mat/pat/child DBs
        meryl count compress k=30 \
            threads=32 memory=96 \
            maternal.*fastq.gz \
            output maternal_compress.k30.meryl

        meryl count compress k=30 \
            threads=32 memory=96 \
            paternal.*fastq.gz \
            output paternal_compress.k30.meryl

        meryl count compress k=30 \
            threads=32 memory=96    \
            child.*fastq.gz output    \
            child_compress.k30.meryl

        ## Create the hapmer DBs
        $MERQURY/trio/hapmers.sh \
        maternal_compress.k30.meryl \
        paternal_compress.k30.meryl \
            child_compress.k30.meryl
```
</details>

#### Trio data: Yak

Yak (Yet-Another Kmer Analyzer) is the kmer counter that we need for Hifiasm assemblies and to QC assemblies made with either assembler so let's learn about how to make yak dbs. 

**In the Meryl section we subset R1, now subset R2 as well**


```bash
zcat /nesi/nobackup/nesi02659/LRA/resources/ilmn/pat/HG003_HiSeq30x_subsampled_R2.fastq.gz \
    | head -n 20000000 \
    | pigz > HG003_HiSeq30x_5M_reads_R2.fastq.gz &
```  

**Look up yak's github and figure out how to make a count/kmer db for this data**

Yak won't work on our Jupyter instances, so create a slurm script that has 32 cores and 96GB of memory. That way it will work on our subset data and it will also work on full size data -- you'd just have to extend the time variable in slurm.
