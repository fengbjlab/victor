# victor
VICTOR = variant interpretation for clinical testing or research

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
