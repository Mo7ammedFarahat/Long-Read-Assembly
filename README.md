## BYOD workshop 21th-25th Oct 2024, Cape Town - South Africa

Note: Material used to prepare for the workshop was extracted from https://github.com/human-pangenomics/hprc-tutorials/tree/GA-workshop and https://genomicsaotearoa.github.io/long-read-assembly/

Customized by: **Mohammed Farahat**, to be up and running on University of Cape Town HPC (ilifu)

**Ilifu HPC Configuration**

Ilifu is a regional node, known as a Tier II node, in a national infrastructure, and partly funded by the Department of Science and Technology (DST) through their Data-Intensive Research Initiative of South Africa (DIRISA).

https://www.ilifu.ac.za/about/


Ilifu is using Slurm as a job scheduling system.

---

Some of the tools required for this workshop are already installed as modules on ilifu, to display the available modules use `module avail` and to load a module use `module load tool/version`.

For the tools that are not preinstalled on ilifu, I have created a conda env that contains all the required tools and packages you will need. 

Download the yml file from here [refgraph.yml](https://github.com/Mo7ammedFarahat/Long-Read-Assembly/blob/main/refgraph.yml)  
Or use it directly from the workshop project dir `/cbio/projects/037/mohammed/condaEnv/refgraph.yml`

To create the environment:
```bash
username@slurm-login:~$conda env create -f refgraph.yml
username@slurm-login:~$conda activate refgraph
```
Now, you can check one of the tools that we will use later, like `mashmap`:
```bash
username@slurm-login:~$mashmap --version
3.1.3
```
---
**Dataset**

  In this workshop we will be using data from HG002, which is a reference sample from the [Genome In A Bottle (GIAB)](https://www.nist.gov/programs-projects/genome-bottle) consortium. The GIAB project releases benchmark data for genomic characterization, and you may have seen their benchmark variant calls and regions out in the wild. As part of their benchmarking material generation, they release datasets for their reference samples. We will be using those in this workshop.
    
    
    HG002 is actually part of a trio of reference samples. Below is the family listing, also known as the Ashkenazim trio:
    
    * HG002: Son
    * HG003: Father
    * HG004: Mother


You will find this dataset here:
`/cbio/projects/037/mohammed/T2T/HG002`




