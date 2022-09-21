# T2T-Assembly

## Telomere-to-telomere assmebly workflow

![Assembly workflow](data/workflow/assembly_workflow.png)


## pipeline


Quick run:

```
snakemake -p  --snakefile  pipeline/T2T_Assembly.snakemake.py --configfile pipeline/T2T_Assembly.snakemake.yaml    --jobs 24 
```


## Tools

Most of the tools listed can be installed via [Anaconda](https://www.anaconda.com/) through searching them in [anaconda](https://anaconda.org/)

The versions listed for thes tools were tested and used in this project.

* snakemake v5.9.1
* fastp v0.20.1
* jellyfish v2.2.10
* genomescope v2.0
* hifiasm v0.15.5-r350
* gfatools v0.4-r214-dirty
* nanoQC v0.8.1
* nextDenovo v2.4.0
* minimap2 v2.17-r941
* racon v1.4.13
* mashmap v2.0
* samtools v1.12
* bwa v0.7.17-r1188
* pilon v1.23
* picard v2.18.29



