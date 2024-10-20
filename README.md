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
