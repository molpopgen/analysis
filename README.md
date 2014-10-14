##The "analysis" package

This is the homepage for __analysis__, which is a set of C++ tools for analyzing population-genetic data

###Warning

I originally wrote these programs with high-quality data in mind (_e.g._ double-pass Sanger sequencing of PCR amplicons).  As a result, they may give incorrect and/or biased results when applied to data from NGS studies.  I still keep the software online for archival purposes, and the tools may yet evolve to the NGS era.

###Revision history

[Here](REVISION_HISTORY.md)

###Source code

The library code is [here](https://github.com/molpopgen/analysis)

###Dependencies 

The package dependes on the following libraries for compilation:

* [libsequence](http://molpopgen.github.io/libsequence)
* [GSL](http://gnu.org/software/gsl)

In addition, you'll need a compiler supporting the C++11 language standard.

###Installation

Simplest case:

```
./configure
make
sudo make install
```

Installing into your home directory:

```
./configure --prefix=$HOME
make
make install
```

Installing when dependencies are in non-standard locations:

(for this example, I'm assuming dependencies are installed into /opt)

```
./configure CXXFLAGS=-I/opt/include LDFLAGS=-L/opt/lib
make 
sudo make install
```

If the depdendencies are in your home folder:

```
./configure CXXFLAGS=-I/$HOME/include LDFLAGS=-L/$HOME/lib
make 
make install
```

The options to ./configure shown above may be mixed and matched.

###Documentation

The following programs have "man" pages that will be available after installation:

* compute
* gestimator
* kimura80
* polydNdS
* sharedPoly
* descPoly
* HBKpermute
* MKtest
* rsq


