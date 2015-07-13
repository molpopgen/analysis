##This is a retired package

As of July 13, 2015, this package is "retired" and will no longer be supported or maintained.

The programs in this package have their origin in the days of using Sanger sequencing to study genetic variation.  The movement towards short-read-based resequencing means that software like this may not be applicable to modern data sets based on modest coverage.

This software did its job for a decade or more, and still has some use in certain applications, but it is time to stop supporting it.

"Retired" means that this code:

* Will no longer be updated to track the development of [libsequence](http://github.com/molpopgen/libsequence).
* Bugs will not be fixed, unless they trace back to libsequence itself.
* Installation issues, etc., will not be dealt with

##The "analysis" package

This is the homepage for __analysis__, which is a set of C++ tools for analyzing population-genetic data

###Warning

I originally wrote these programs with high-quality data in mind (_e.g._ double-pass Sanger sequencing of PCR amplicons).  As a result, they may give incorrect and/or biased results when applied to data from NGS studies.  I still keep the software online for archival purposes, and the tools may yet evolve to the NGS era.

###Revision history

[Here](REVISION_HISTORY.md)

###Source code

The source code is available [here](https://github.com/molpopgen/analysis)

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


