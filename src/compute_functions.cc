#include <iostream>
#include <ctime>
#if defined (__SVR4) && defined (__sun)
#include <ieeefp.h>
#endif

#include <compute_functions.hpp>
#include <compute_classes.hpp>
#include <gsl/gsl_rng.h>
using namespace std;

void process(glob_t *files,  compute_params *args, ostream & ofstr,
	      int  argc,  char** argv)
{
  gsl_rng * rng = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(rng,time(0));
  for(int i=0;i<files->gl_pathc;++i)
    {
      results r(files->gl_pathv[i],args);
      if (i==0)
	{
	  if (! args->suppress_headers)
	    {
	      if(args->outfile != NULL)
		makeheader(argc,argv,args,ofstr);
	      else
		{
		  if (args->pretty == false)
		    makeheader(argc,argv,args,cout);
		}
	    }
	}
      if (r.wasSkipped())
	{
	  cerr << "file "
	       << files->gl_pathv[i]
	       << " was skipped...\n";
	}
      else
	{
	  if (args->pretty)
	    ofstr << r << '\n';
	  else
	    ofstr << r << '\t';
	  if (args->probs)
	    {
	      //get pvals from 10000 coal sims w/o recombination
	      pvals p(rng,r,args);
	      ofstr << p;
	    }
	}
     if (args->pretty)
	ofstr << "\n//\n";
      else
	ofstr << '\n';
    }
}

void makeheader(const int argc, char *argv[],const compute_params * args,ostream &o)
{
  time_t timeval;
  (void)time(&timeval);
  o << '#';
  for(int i = 0 ;i < argc;++i)
    o << argv[i] << ' ';
  o << endl
    << '#' << ctime(&timeval)
    << endl
    << "locus\tnsam\tnsites\tnsites_ug\tS\tSingletons\t";
  if(args->haveOutgroup)
    o << "DerSingletons\t";
  o << "NMut\tnhap\thapdiv\tWallsB\tWallsQ\tThetaW\tThetaPi\t";
  if(args->haveOutgroup)
    o << "ThetaH\t";
  o << "TajD\t";
  if (args->haveOutgroup)
    o << "FuLiD\tFuLiF\t";
  else
    o << "FuLiDStar\tFuLiFStar\t";
  
  o << "Rmin\t"
    << "rho87";

  if (args->probs)
    {
      o<< '\t'
       << "pTajD\t"
       << "pnhap\t"
       << "phapdiv\t"
       << "pWallsB\t"
       << "pWallsQ\t";
      if(args->haveOutgroup)
        {
	  o << "pFuLiD\t"
	    << "pFuLiF\t"
	    << "pFayWuH\t";
        }
      else
        {
	  o << "pFuLiDStar\t"
	    << "pFuLiFStar\t";
        }
    }
  o << endl;
}
