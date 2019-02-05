# AD-LIBS
AD-LIBS is a local ancestry inference tool for low-coverage resequencing data. It is described in [this paper](https://doi.org/10.1186/s12859-017-1613-0). The goal of AD-LIBS is to take a low-coverage genome from one or more hybrid/admixed individuals, along with one or more low-coverage genomes from representatives of two ancestral populations, and produce a map of where in the genome the hybrid individual derives ancestry from each ancestral population. Maps consist of fixed-width windows across the genome, and each window is labeled as homozygous for ancestry from one or the other ancestral population, or heterozygous for ancestry. Because of this, AD-LIBS is only suited for analysis of diploid organisms.

To circumvent problems with identifying heterozygous sites in low-coverage samples, AD-LIBS reads in one FASTA file per individual. FASTA files must be aligned to the same reference genome and should be "pseudo-haploid," where every base is randomly sampled from one or the other parental haplotype, as has been done in studies like [this one](https://doi.org/10.1371/journal.pgen.1003345). To create pseudo-haploid FASTA files, we recommend using the [pu2fa program](https://github.com/Paleogenomics/Chrom-Compare) after mapping reads. 

Results are in [BED format](https://genome.ucsc.edu/FAQ/FAQformat.html#format1), with a fourth column denoting ancestry in a given window (AA for homozygous population A, AB for heterozygous, and BB for homozygous population B). [BEDTools](https://github.com/arq5x/bedtools2) contains many programs helpful for downstream analysis of these files.

## Installation

### Python packages
AD-LIBS requires the following Python packages:

- [SciPy](https://www.scipy.org/)
- [NumPy](http://www.numpy.org/)
- [Matplotlib](http://matplotlib.org/)
- [pomegranate](https://pomegranate.readthedocs.io/en/latest/)

#### Anaconda
An easy way to obtain the first three is to use [Anaconda](https://www.continuum.io/downloads), 
which contains these packages.

If using anaconda and they are not installed, try using `conda install`:

`conda install scipy`

`conda install numpy`

`conda install matplotlib`

Note that pomegranate is not available through Anaconda.

#### Pip
These packages are also installable using [Pip](https://pypi.python.org/pypi/pip):

`pip install scipy`

`pip install numpy`

`pip install matplotlib`

`pip install pomegranate`

If you do not have root access, you install these packages for your user using the --user flag.

### C programs
AD-LIBS equires gcc for compiling C components, as well as [zlib](http://www.zlib.net/) 
for reading compressed [FASTA](https://en.wikipedia.org/wiki/FASTA_format) files.

To compile C programs, simply type `make` in the main project directory.
Binary files will be stored in `/bin` and run by the main Python program.

Note that `src/kseq.h` is an excellent FASTA parser in C, taken from 
[HTSLib](https://github.com/samtools/htslib), part of the 
[SAMTools](https://github.com/samtools) codebase.

## Running AD-LIBS

The Python program bin/adlibs.py is the main program with which users should interact.
AD-LIBS works in pieces; this program interprets input, `bin/adlibs_score` 
reads in FASTA files and prints out scores computed in windows across sequences, 
and `bin/adlibs_hmm.py` uses the Pomegranate package to build and run the 
Viterbi algorithm on the HMM. Additionally, `bin/allele_freq_drift.py` uses the 
diffusion approximation to the genetic drift problem to compute an important 
probability used by the HMM, and `adlibs_math.py` and `adlibs_dists.py` both 
contain code used by `adlibs_hmm.py`.

`bin/adlibs.py` can use command-line input (run with -h for help) or, since 
the many parameters can get cumbersome, a configuration file containing all 
the information needed for one run. Whether or not you provide a configuration 
file, AD-LIBS will create a new one with all given parameters along with your 
output files, in case you want to run AD-LIBS again with all or most of the 
same parameters.

Note that in order for `bin/adlibs.py` to work, you will need to have 
`python2.7` in your `$PATH` environment variable, or you can run it by 
explicitly including the path to `python2.7` as follows:

`/path/to/python2.7 bin/adlibs.py`

## Example input

The example directory contains an example configuration file, along with 
(short) pseudo-haploid FASTA sequences of genomic scaffolds for two brown bears 
(`OFS01_s317.fa` and `Swe_s317.fa`), two polar bears (`PB12_s317.fa` and 
`PB105_s317.fa`), and two admixed ABC Islands brown bears (`ABC01_s317.fa` and 
`ABC02_s317.fa`). `adlibs.config` contains all the parameters needed to run 
AD-LIBS on this data set, with polar bears as ancestral population "A" and 
brown bears as ancestral population "B". All output will be prefixed with 
"example." and the `adlibs.config` file contains comments explaining the 
parameters.

To run AD-LIBS on this example data set, simply type 
`bin/adlibs.py example/adlibs.config` while in the main project directory. 
When complete, you will see 
[BED files](https://genome.ucsc.edu/FAQ/FAQformat#format1) 
containing results in the example directory. Each line is a genomic window, 
with start and end coordinates, along with a label (AA for homozygous 
population A ancestry, AB for heterozygous ancestry, and BB for homozygous 
population B ancestry). You can try altering parameters and re-running as you 
wish to see how it affects results.

## Troubleshooting

If AD-LIBS creates one or more empty BED files, the most likely problem is that one or more of your genomes has too many ambiguous bases (Ns). In this case, you can try to increase the skip threshold parameter -- this is the percent of all bases in each window allowed to be N before the window is skipped -- so setting it to 0.75 means AD-LIBS will still process a window that is up to 75% N bases. If this doesn't work, you can increase the window size to try to increase the number of valid bases in each window. Otherwise, you can re-generate your pseudohaploid FASTA files using lower quality filters to allow more bases to come through.

## Citing AD-LIBS

AD-LIBS was described in the article [AD-LIBS: inferring ancestry across hybrid genomes using low-coverage sequence data](https://doi.org/10.1186/s12859-017-1613-0), NK Schaefer, B Shapiro, and RE Green, BMC Bioinformatics, 2017. We'd appreciate a citation if you use it in your research.

