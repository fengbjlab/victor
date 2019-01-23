# VICTOR

VICTOR = variant interpretation for clinical testing or research.

### How to download

```
git clone https://github.com/fengbjlab/victor.git
```

### How to compile

The following commands will compile the programs and put the binary executable files in victor/.

```
cd victor/src
OS=`uname | tr '[:upper:]' '[:lower:]'`
cp makefile.$OS makefile
make
make clean
```

### How to get data

```
cd victor/data
./get_GRCh37_inst.sh
./get_GRCh37_update.sh
```

### Third-party programs

You need to install the following third-party programs:

|Programs     |Comment    |It will be used for        | License   | URL                                                                  |
|-------------|-----------|---------------------------|-----------|----------------------------------------------------------------------|
|logistf(1.23)|required   |logistic regression        | GPL       | https://cran.r-project.org/web/packages/logistf/logistf.pdf          |
|tabix        |required   |retrieval of file contents | MIT       | https://sourceforge.net/projects/samtools/files/tabix/               |
|GNU parallel |required   |parallel computing         | GPLv3     | https://www.gnu.org/software/parallel/                               |
|PROVEAN      |recommended|PROVEAN score calculation  | GPLv3     | http://provean.jcvi.org/                                             |
|blast        |recommended|PROVEAN                    | Public    | https://blast.ncbi.nlm.nih.gov/                                      |
|CD-HIT 4.5.8 |recommended|PROVEAN                    | GPLv2     | http://weizhongli-lab.org/cd-hit/                                    |
|PLINK 1.9    |recommended|quality control            | GPLv3     | https://www.cog-genomics.org/plink2/                                 |
|KING         |recommended|quality control            | Unknown   | http://people.virginia.edu/~wc9c/KING/                               |
|ShapeIt2     |optional   |phasing                    | Academic  | https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html|
|BEAGLE       |optional   |phasing                    | GPLv3     | https://faculty.washington.edu/browning/beagle/beagle.html           |
|GATK         |optional   |VCF file combining         | BSD       | https://software.broadinstitute.org/gatk/download/                   |
|gnuplot      |optional   |drawing Manhattan plots    | gnuplot   | http://www.gnuplot.info/                                             |

### How to use

Put the victor/ folder in your PATH or module. Please see the doc/ folder for manuals and tutorials. These PDF files were taken from VICTOR's website http://bjfenglab.org/victor/. 
