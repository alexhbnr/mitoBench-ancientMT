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
  - unpack the tarball of the R package contamMix
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
          --cluster 'qsub -pe smp {threads} -l virtual_free={cluster.vfree},h_vmem={cluster.hvmem},class={cluster.class} -o {cluster.out} -e {cluster.err}' \
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
          --cluster 'sbatch --mem {cluster.mem} -p {cluster.partition} -t {cluster.time} -o {cluster.out} -e {cluster.err} -n {threads}' \
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

In all three examples, replace the config file `mitoBench_pipeline-config.json`
with the version of the config file in which you specified the paths to the
sample list, the temporary directory, the project directory etc.

## Overview of the pipeline

This pipeline is built to process per-sample BAM files and return a mtDNA
consensus sequence and some additional quality statistics, e.g. coverage across
the mitochondrial genome, DNA damage profile, mtDNA haplogroup etc. The pipeline
assumes that raw sequencing data has been de-multiplexed, adapters have been
removed and overlapping reads have been merged. Additionally, it expects an
alignment against the human reference genome (hg19) and that the output is
sorted BAM file.

The steps of the pipeline are:

  1. Create a softlink of each BAM file specified in the `samplelist` into the
     temporary directory and create an index using `samtools index`.
  2. Extract all reads aligned to the mitochondrial genome (assuming that the
     chromosome name is "MT") regardless of the mapping quality.
  3. Align the extracted reads only against the rCRS reference genome, which is
     extended by pasting the first 1,000 base pairs to the end of the sequence
     to account for the circular nature of the DNA molecule, using `bwa aln`.
     Generate alignment files using `bwa sampe` & `bwa samse`.
  4. Join the separate output files for merged and non-merged paired-end reads
     using `samtools merge`, adjust the alignment coordinates to the original
     length of rCRS (16,569 bp) using circularmapper's `realign` and sort the
     output.
  5. Remove duplicate sequences using `DeDup` and sort DeDup's output file
     again.
  6. Analyse the read length distribution and the DNA damage profile using
     `DamageProfiler`.
  7. Count the number of reads available per sample using `samtools flagstat`
     and calculate the fraction of reads that correspond to 30,000 reads for
     later down-sampling using `samtools view` using seed 0.
  8. Infer the per-base sequencing depth of all samples using `samtools depth`.
  9. Calculate the number of bases to trim from both the 5' and 3' end of a read
     by calculating the mean frequency of C to T substitutions at the 5' end and
     of G to A substitutions on the 3' end, respectively, as well as the
     standard deviation. All bases are trimmed that are above the threshold of
     the mean frequency plus one standard deviation regarding their
     substitutions, in case the threshold is >= 2%, otherwise no bases are
     trimmed. Trim the bases using `bamUtil trimBam`.
  10. Call the consensus sequences using `ANGSD` requiring a minimal base
      quality of 30 and minimal mapping quality of 20 based on both the clipped
      and non-clipped sequencing data.
  11. Infer the haplogroup using `HaploGrep2` using the consensus sequence
      called from the clipped sequencing data.
  12. Infer the haplogroups and their proportion of contribution to the pool of
      sequences using `mixEMT`.
  13. Calculate the proportions of authentic, non-contaminant DNA using
      `contamMix`.

When all these steps are run, the results are summarised in a table called
`summary_table.csv`. This table contains the following information:

  - the number of unique reads that aligned to the rCRS
  - the mean and the standard deviation of the coverage across the genome and
    the number of sites with at least 5-fold coverage
  - the mode of the read length distribution
  - the number of trimmed bases separately calculated for the 5' and the 3' end
  - the percentage of missing information (Ns) in the consensus sequences called
    from the clipped sequencing data
  - the number of differences between the consensus sequences called from the
    clipped and non-clipped sequencing data
  - the inferred haplogroup, the quality of the haplogroup assignment and the
    number of non-found and remaining polymorphisms, respectively
  - the proportion of authentic, non-contaminant DNA
  - the haplogroups and their proportions of contribution to the sequencing pool

Further, the consensus sequences based on the clipped sequencing data will be
copied to the directory `fasta` in the specified project folder and the raw
information that were summarised in the summary table are copied into a
per-sample folder in the directory `sample_stats`.

The temporary folder can later be removed by running the Snakemake pipeline with
the specified rule `clean_tmp`.
