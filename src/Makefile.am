bin_PROGRAMS=rsq compute kimura80  gestimator  polydNdS  HBKpermute MKtest sharedPoly descPoly translateCoding extractCoding Fexact snntest
rsq_SOURCES=rsq.cc \
	Correlations.hpp \
	RandomNumbers.hpp \
	rsq.1
compute_SOURCES=compute.cc \
	compute_functions.hpp \
	compute_functions.cc \
	compute_classes.hpp\
	compute_classes.cc\
	int_handler.hpp \
	int_handler.cc \
	compute.1
kimura80_SOURCES=distance.cc kimura80.1 
gestimator_SOURCES=gestimator.cc \
	int_handler.cc \
	int_handler.hpp \
	gestimator.1
polydNdS_SOURCES = polydNdS.cc \
	polydNdSbase.cc polydNdSbase.hpp \
	tstvbias.cc tstvbias.hpp \
	codingRegionProcessor.cc codingRegionProcessor.hpp \
	cRPimpl.cc cRPimpl.hpp \
	CodonMaker.hpp CodonMaker.cc \
	Warnings.cc Warnings.hpp \
	polydNdS.1
translateCoding_SOURCES = translateCoding.cc \
	polydNdSbase.cc  polydNdSbase.hpp 
MKtest_SOURCES=MKtest.cc \
	MKtestOutput.hpp MKtestOutput.cc\
	polydNdSbase.cc polydNdSbase.hpp \
	codingRegionProcessor.hpp codingRegionProcessor.cc \
	cRPimpl.hpp cRPimpl.cc \
	CodonMaker.hpp CodonMaker.cc \
	Warnings.hpp Warnings.cc \
	chisq.cc chisq.hpp \
	ctest.h fexact.c g2.h g2.c \
	MKtest.1

HBKpermute_SOURCES = HBKpermute.cc HBKpermute.1
sharedPoly_SOURCES = sharedPoly.cc sharedPoly.1
descPoly_SOURCES = descPoly.cc \
	polyCommon.cc polyCommon.hpp \
	descPoly.1
extractCoding_SOURCES =  extractCoding.cc \
	polydNdSbase.cc polydNdSbase.hpp 
Fexact_SOURCES=Fexact_main.cc fexact.c ctest.h

snntest_SOURCES=snntest.cc #getopt.c getopt1.c getopt.h getoptFix.h


# if DEBUG
# OTHERFLAGS = -g
# else
# OTHERFLAGS = -DNDEBUG -O2
# endif

# if HAVE_LIMITS
# LIMITS= -DHAVE_LIMITS
# else
# LIMITS =
# endif

# if HAVE_CFLOAT
# CFLOAT= -DHAVE_CFLOAT
# else
# CFLOAT =
# endif

if HAVE_GSL_HEADERS
MKtest_LDADD=-lm -lgsl -lgslcblas
compute_LDADD=-lgsl -lgslcblas
snntest_LDADD=-lgsl -lgslcblas
rsq_LDADD=-lgsl -lgslcblas
GSL= -DHAVE_GSL
else
GSL=
endif

#polydNdS_LDADD=-lboost_regex

# if HAVE_STRSTREAM
# STRSTREAM = -DHAVE_STRSTREAM
# else
# STRSTREAM =
# endif

# if HAVE_SSTREAM
# SSTREAM = -DHAVE_SSTREAM
# else
# SSTREAM =
# endif

#if HAVEGCC
#STD=-std=c99
#else
#STD=
#endif 

AM_CXXFLAGS= $(OTHERFLAGS)  $(STRSTREAM) $(SSTREAM) $(LIMITS) $(CFLOAT) $(GSL) -Wall -I$(includedir) 
AM_CFLAGS= $(STD) $(OTHERFLAGS) -Wall -I$(includedir) $(GSL) 
CLINK= -shared-libgcc


man_MANS= polydNdS.1 rsq.1 kimura80.1 gestimator.1 compute.1 MKtest.1 sharedPoly.1 descPoly.1 HBKpermute.1


