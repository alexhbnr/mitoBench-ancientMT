# mitoBench - ancient mtDNA pipeline

This is the pipeline for processing (non-)published ancient mtDNA genomes for
subsequent import into mitoBench.

## Installation

In order to guarantee that all samples are processed reproducibly independent
from which versions of underlying programs are installed on the system, we build
this pipeline upon [Conda](https://anaconda.org/anaconda/conda) as an
environment management system and
[Snakemake](https://snakemake.readthedocs.io/en/stable/) as a workflow
management system. 

For this to work, we assume that **Conda** is installed on the system. If this
is not the case, **Conda** is available for download from
[here](https://docs.conda.io/en/latest/miniconda.html). As **Snakemake** is
relying on Python3, please choose the *Python3 64-bit installer* for your system.

Next to Conda, **Snakemake** with version >=5.0 needs to installed on the
system, too. In order to guarantee this, you can type

```bash
pip install 'snakemake>=5.0'
```

While most programs will be automatically downloaded and installed from the
[Conda YAML](env/mitoBench_bioconda.yaml) using the
[BioConda](https://bioconda.github.io/) channel, there are few programs and
other resources that need to be separately installed. To do so, run

```bash
snakemake -s setup_env.Snakefile --use-conda
```

This will

  - download the commandline version of [HaploGrep](http://haplogrep.uibk.ac.at/)
    for determining the mtDNA haplogroups
  - download the Python package [mixemt](https://github.com/svohr/mixemt) for
    estimating the mtDNA contamination based on read mixtures
  - prepare the revised Cambridge reference sequence (rCRS) as a reference
    genome for the subsequent analyses of the pipeline

## Run pipeline

The pipeline itself is stored in the file `mitoBench_pipeline.Snakefile`. In
order to execute, two config files are required: one that provides the
information about input and output files, and one for providing the information
on the required resources to the scheduler system for cluster usage.

The config file contains the major information on where the input files are
stored, the naming pattern of the samples and input files, and where the output
files are to be written. This file should be adapted for each project and a
template can be found: [mitoBench_pipeline-config.json](mitoBench_pipeline-config.json).
Six different parameters have to specified:

  1. `samplelist`: path to a text file with the sample IDs, one sample ID per
     line. `snakemake` will determine the wildcards for expected output files
     from this filelist.
  2. `bamdir`: the folder in which the BAM files of the sample specified in the
     samplelist are stored. It is expected that all BAM files are stored in the
     same directory.
  3. `projdir`: the final output directory for the results of the pipeline
  4. `tmpdir`: the directory for storing the temporary output during processing
  5. `bamsuffix`: file suffix after the sample ID following the scheme
     `<sampleID>.<file suffix>` (default: bam)
  6. `sampleIDconstraint`: in order to avoid unexpected behaviour, only
      sample IDs with uppercase and lowercase letters, numbers, and periods are
      allowed. Underscores in the sample IDs should be avoided as the pipeline
      uses underscores to separate the file extensions from the sampleIDs. In
      case, sample IDs that do not follow this scheme are to be used, the regex
      constraining the sample ID pattern can be altered (default: [A-Za-z0-9\.]+).

The cluster config file specifies all the computational resource requirements
that each step needs to be run. There are two cluster config files provided in
JSON format for the cluster scheduling systems
[SGE](mitoBench_pipeline-SGE.json) and [SLURM](mitoBench_pipeline-SLURM.json).
These files should set the correct requirements for the expected sequencing
depths for ancient DNA capture data. If you sequenced your samples to a much
deeper depth, higher requirements regarding memory usage might be necessary.

Next to submitting the individual pipeline steps to the cluster, a few steps
will be directly run on the head node, when the run time and memory requirement
of a step are so low that the overhead of sending the step to the cluster is
larger than just quickly running them on the head node itself. In order to
avoid an excessive usage of resources on the head node, we will restrict the
maximum number of head node cores that can be used using the option
`local-cores`.

In order to run the pipeline using **SGE** as the scheduling system run:

```bash
snakemake -s mitoBench_pipeline.Snakefile \
          --configfile mitoBench_pipeline-config.json \
          --cluster-config mitoBench_pipeline-SGE.json \
          --cluster 'qsub -pe smp {threads} -l virtual_free={cluster.vfree},h_vmem={cluster.hvmem},class={cluster.class}' \
          --use-conda \
          --local-cores 8 \
          --cores 100
```

If some of these parameters are not available or not feasible for your local
scheduling system, please alter them accordingly.

Using **SLURM** as the scheduling system:

```bash
snakemake -s mitoBench_pipeline.Snakefile \
          --configfile mitoBench_pipeline-config.json \
          --cluster-config mitoBench_pipeline-SLURM.json \
          --cluster 'sbatch --mem {cluster.mem} -p {cluster.partition} -t {cluster.time}' -o {cluster.out} -e {cluster.err} -n {threads}' \
          --use-conda \
          --local-cores 8 \
          --cores 20
```

In order to run the pipeline on a local computer without using a scheduling
system:

```bash
snakemake -s mitoBench_pipeline.Snakefile \
          --configfile mitoBench_pipeline-config.json \
          --use-conda \
          --cores 20
```

## Overview of the pipeline
