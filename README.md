# VICTOR

VICTOR (Variant Interpretation for Clinical Testing Or Research) is a pipeline for the analysis of next-generation sequencing data starting from a multi-sample raw VCF file. It can be used for disease gene discovery research or clinical genetic testing. It is designed to be scalable to whole genome sequencing (WGS) of a large sample of individuals that is typical of a research on a complex disease. It is mostly a self-contained package, whereby the requirement for third-party applications is minimal (see below), and all necessary databases are already included. It provides data updates on a monthly basis.

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
# to get GRCh37
cat get_GRCh37_inst.sh | bash
cat get_GRCh37_update.sh | bash
# to get GRCh38
cat get_GRCh38_inst.sh | bash
cat get_GRCh38_update.sh | bash
# to get hg19
cat get_hg19_inst.sh | bash
cat get_hg19_update.sh | bash
```

### Third-party programs

You need to install the following third-party programs:

|Programs     |Comment    |Usage          | License   | URL                                                                  |
|-------------|-----------|---------------|-----------|----------------------------------------------------------------------|
|logistf(1.23)|required   |analysis       | GPL       | https://cran.r-project.org/web/packages/logistf/logistf.pdf          |
|tabix        |required   |read file      | MIT       | https://sourceforge.net/projects/samtools/files/tabix/               |
|GNU parallel |required   |parallelism    | GPLv3     | https://www.gnu.org/software/parallel/                               |
|PROVEAN      |recommended|PROVEAN        | GPLv3     | http://provean.jcvi.org/                                             |
|blast        |recommended|PROVEAN        | Public    | https://blast.ncbi.nlm.nih.gov/                                      |
|CD-HIT 4.5.8 |recommended|PROVEAN        | GPLv2     | http://weizhongli-lab.org/cd-hit/                                    |
|PLINK 1.9    |recommended|QC             | GPLv3     | https://www.cog-genomics.org/plink2/                                 |
|KING         |recommended|QC             | Unknown   | http://people.virginia.edu/~wc9c/KING/                               |
|ShapeIt2     |optional   |phasing        | Academic  | https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html|
|BEAGLE       |optional   |phasing        | GPLv3     | https://faculty.washington.edu/browning/beagle/beagle.html           |
|GATK         |optional   |combine VCF    | BSD       | https://software.broadinstitute.org/gatk/download/                   |
|gnuplot      |optional   |Manhattan plots| gnuplot   | http://www.gnuplot.info/                                             |

### How to use

Put the victor/ folder in your PATH or module. Please see the doc/ folder for manuals and tutorials. These PDF files were taken from VICTOR's website http://bjfenglab.org/victor/. 
