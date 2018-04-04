# hiplexpipe
#Generic_pipe

## A very generic bioinformatics pipeline for variant calling.

Original author of wrapper scripts Khalid Mahmood (kmahmood@unimelb.edu.au).
This version modified heavily from original intended purpose and for use Monash M3.massive cluster by Jason Steen (jason.steen@monash.edu)


genericpipe is based on the [Ruffus](http://www.ruffus.org.uk/) library for writing bioinformatics pipelines. Its features include:

 * Job submission on a cluster using DRMAA (currently only tested with SLURM).
 * Job dependency calculation and checkpointing.
 * Pipeline can be displayed as a flowchart.
 * Re-running a pipeline will start from the most up-to-date stage. It will not redo previously completed tasks.

## License

See LICENSE.txt in source repository.

## Installation

#### External tools dependencies

`genericpipe` depends on the following programs and libraries:

 * [python](https://www.python.org/download/releases/2.7.5/) (version 2.7.5)
 * [java](https://java.com/en/download/) (version 1.8)
 * [DRMAA](http://www.drmaa.org/) for submitting jobs to the cluster (it uses the Python wrapper to do this).
   You need to install your own `libdrama.so` for your local job submission system. There are versions
   available for common schedulers such as Torque/PBS, [SLURM](http://apps.man.poznan.pl/trac/slurm-drmaa) and so on.
 * [SAMtools](http://www.htslib.org/doc/samtools-1.1.html) (version 1.3.1)
 * [bwa](http://bio-bwa.sourceforge.net/) for aligning reads to the reference genome (version 0.7.15)  
 * [GATK](https://software.broadinstitute.org/gatk/) for calling variants and genotyping (version 3.6)

`hiplexpipe` assumes the tools above are installed by the users themselves.

#### Python dependencies

`genericpipe` depends on the following python libraries, tools and wrappers.

* Python 2.7.5+
* [PyVCF](https://pypi.python.org/pypi/PyVCF)  
* [Biopython](https://pypi.python.org/pypi/biopython)
* [pybedtools](https://daler.github.io/pybedtools/)
* [cyvcf2](http://brentp.github.io/cyvcf2/)

We recommend using a python virtual environment. Following is an examples of how to setup a `hiplexpipe` virtual environment ready for analysis:

#### Installation example on Monash clusters

```
/usr/local/python/2.7.12_static/bin/virtualenv --system-site-packages pipeline_env
module load drmaa
export DRMAA_LIBRARY_PATH=/usr/local/drmaa/1.0.7/lib/libdrmaa.so
source pipeline_env/bin/activate
pip install numpy scipy
pip install git+https://github.com/SoutheyLab/genericpipe
genericpipe --config pipeline.config --use_threads --log_file pipeline.log --jobs 10 --verbose 3 --just_print
```


The following lines need to be added to your .bash_rofile so that VEP will function correctly.

```
PERL5LIB=$PERL5LIB:/projects/vh83/reference/annotation_databases/VEP_Modules/.vep/Plugins
export PERL5LIB
```


[to be continued]

