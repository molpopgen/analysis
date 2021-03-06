AC_PREREQ(2.59)

AC_INIT(analysis, [0.8.8], [https://github.com/molpopgen/analysis])
AC_CONFIG_SRCDIR([src/compute.cc])
AM_INIT_AUTOMAKE
AC_CONFIG_HEADERS([config.h])

AC_CONFIG_MACRO_DIR([m4])

AC_PROG_CC
AC_C_CONST
AC_PROG_CXX
AX_CXX_COMPILE_STDCXX_11([noext],[mandatory])

AM_MAINTAINER_MODE([disable])

AC_LANG(C++)
AC_CONFIG_FILES([Makefile src/Makefile])

AC_ARG_ENABLE(debug,
[  --enable-debug    Turn on debugging],
[case "${enableval}" in
  yes) debug=true ;;
  no)  debug=false ;;
  *) AC_MSG_ERROR(bad value ${enableval} for --enable-debug) ;;
esac],[debug=false])
AM_CONDITIONAL(DEBUG, test x$debug = xtrue)

AC_ARG_ENABLE(profiling,
[  --enable-profiling    Turn on profiling],
[case "${enableval}" in
  yes) profiling=true ;;
  no)  profiling=false ;;
  *) AC_MSG_ERROR(bad value ${enableval} for --enable-profiling) ;;
esac],[profiling=false])
AM_CONDITIONAL(PROFILING, test x$profiling = xtrue)

AC_MSG_CHECKING(for namespace Sequence::coalsim)
AC_TRY_COMPILE(
[
#include <Sequence/Coalescent/Coalescent.hpp>                                                            

using namespace Sequence::coalsim;
],
[
Sequence::coalsim::arg x;
],
[AC_MSG_RESULT(yes)
AC_DEFINE([HAVE_SEQUENCE_COALSIM],[],"Have namespace Sequence::coalsim")],
[AC_MSG_RESULT(no)])

dnl check for things that we need

dnl boost (only needed prior to 1.8.3)
dnl AC_CHECK_HEADER(boost/type_traits.hpp, LIBSEQUENCE_USE_BOOST=1 , [AC_MSG_ERROR([boost::type_traits required - either install boost (see http://www.boost.org) or write a replacement])])
dnl AC_CHECK_HEADER(boost/utility.hpp, LIBSEQUENCE_USE_BOOST=1 , [AC_MSG_ERROR([boost::noncopyable required (boost/utility.hpp) - either install boost (see http://www.boost.org) or write a replacement])])
dnl AC_CHECK_HEADER(boost/static_assert.hpp, LIBSEQUENCE_USE_BOOST=1 , [AC_MSG_ERROR([BOOST_STATIC_ASSERT required - either install boost (see http://www.boost.org) or write a replacement])])
dnl AC_CHECK_HEADER(boost/tuple/tuple.hpp,,[AC_MSG_ERROR([boot::tuple required - either install boost (see http://www.boost.org) or write a replacement])])

dnl gsl.  dependency removed in 1.8.3
dnl AC_CHECK_HEADER(gsl/gsl_rng.h, LIBSEQUENCE_USE_GSL=1 , [AC_MSG_ERROR([gsl/gsl_rng.h not found.  Please install GSL from http://www.gnu.org/software/gsl])])

GSL_HEADERS_FOUND="false"
AC_CHECK_HEADER(gsl/gsl_sf_gamma.h,GSL_HEADERS_FOUND="true",)
AC_CHECK_HEADER(gsl/gsl_cdf.h,GSL_HEADERS_FOUND="true",)
AC_CHECK_HEADER(gsl/gsl_rng.h,GSL_HEADERS_FOUND="true",[AC_MSG_ERROR([gsl/gsl_rng.h required for compilation - please install GSL (www.gnu.org/software/gsl)])])
AC_CHECK_HEADER(gsl/gsl_randist.h,GSL_HEADERS_FOUND="true",[AC_MSG_ERROR([gsl/gsl_rng.h required for compilation - please install GSL (www.gnu.org/software/gsl)])])
AM_CONDITIONAL(HAVE_GSL_HEADERS,test x$GSL_HEADERS_FOUND = xtrue)
AM_CONDITIONAL(HAVEGCC,test x$GCC = xyes)

dnl check for zlib
ZLIB_HEADERS_FOUND="false"
AC_CHECK_HEADER(zlib.h,ZLIB_HEADERS_FOUND="true")
AM_CONDITIONAL(HAVE_ZLIB_HEADERS,test x$ZLIB_HEADERS_FOUND = xtrue)

AC_LANG_SAVE
          AC_LANG_CPLUSPLUS
          AC_REQUIRE_CPP
	  AC_CHECK_LIB(sequence, main,,[echo "libsequence runtime library does not seem to be here";exit 1])	
AC_LANG_RESTORE
dnl AM_CONFIG_HEADER(config.h)

dnl check for libsequence < 1.8.3 or >= 1.8.4. differences
AC_MSG_CHECKING(for Sequence::make_polySiteVector)
AC_TRY_LINK(
[
#include <Sequence/polySiteVector.hpp>                                                            
#include <Sequence/SimData.hpp>
using namespace Sequence;
],
[
SimData d;
auto rd = make_polySiteVector(d);
],
[AC_MSG_RESULT(yes)
AC_DEFINE([HAVE_SEQUENCE_MAKE_POLYSITEVECTOR],[],"Have Sequence::make_polySiteVector")],
[AC_MSG_RESULT(no)])

AC_MSG_CHECKING(for Sequence::rotatePolyTable)
AC_TRY_LINK(
[
#include <Sequence/PolyTableManip.hpp>                                                            
#include <Sequence/SimData.hpp>
using namespace Sequence;
],
[
SimData d;
auto rd = rotatePolyTable(&d);
],
[AC_MSG_RESULT(yes)
AC_DEFINE([HAVE_SEQUENCE_ROTATEPOLYTABLE],[],"Have Sequence::rotatePolyTable")],
[AC_MSG_RESULT(no)])

AC_OUTPUT
