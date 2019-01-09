# victor
VICTOR = variant interpretation for clinical testing or research

### How to compile

1. Install the required libraries listed below. In Mac, install to /opt/local/; in Linux, install to ~/local/. If you like to install other locations, modify the makefile accordingly.

* bzip2(1.0.6)		http://bzip.org
* zlib(1.2.8)		http://zlib.net
* boost(1.58.0)		http://www.boost.org
* TFT(2019-Jan-03)	https://github.com/fengbjlab/libTFT
* NLopt(2.4.1)		https://nlopt.readthedocs.io/

2. Compile

For linux
```
cd src
mv makefile.linux
make
make clean
```

