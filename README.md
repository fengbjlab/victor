# victor
VICTOR = variant interpretation for clinical testing or research

### What are included

##### documents in the doc/ folder

##### c++ source codes in the src/ folder

##### precompiled libraries in mac/ and linux/. Libraries include the following:

* bzip2(1.0.6)		http://bzip.org
* zlib(1.2.8)		http://zlib.net
* boost(1.58.0)		http://www.boost.org
* TFT(2019-Jan-03)	https://github.com/fengbjlab/libTFT
* NLopt(2.4.1)		https://nlopt.readthedocs.io/
* Eigen(3.3.3)		http://eigen.tuxfamily.org/
* ptoc

### How to download

```
git clone https://github.com/fengbjlab/victor.git
```

### How to compile

The following commands will compile the programs, create a bin/ folder, and put the binary executable files in the bin/.

```
cd victor/src
OS=`uname | tr '[:upper:]' '[:lower:]'`
cp makefile.$OS makefile
make
make clean
```
