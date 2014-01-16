#include <MKtestOutput.hpp>
#include <iostream>
#include <iterator>
#include <algorithm>
#include <numeric>

//following headers are for calculating test statistics
//ctest.h is from the R project
#include <ctest.h>
#include <chisq.hpp>
#include <g2.h>

using namespace std;
enum {FIXEDA=0,FIXEDS,POLYA,POLYS};

MKtestOutput::MKtestOutput( const vector<unsigned> & contingency_table,
			    const params & args ) :
  _contingency_table(&contingency_table),
  _args(&args)
{
}
 
std::ostream & MKtestOutput::print(std::ostream & s)
{
  //If we have GSL, we do the following tests: Fisher's Exact, G-test, and Chi-squared,
  //else, only do the Fisher's test.

  bool untestable_table =  ( ((*_contingency_table)[0]==0&&(*_contingency_table)[1]==0) ||
			     ((*_contingency_table)[2]==0&&(*_contingency_table)[3]==0) );
  vector<double> contingency_table_dbl(_contingency_table->begin(),
				       _contingency_table->end());
  //Fisher's Exact test
  //BEG: for fexact routine from R
  int a=2,y=2;
  int work=100000;
  int leading = 2;
  double expect,percnt,emin,prt,pre;
  expect=percnt=emin=prt=pre=-1.;
  //END: for fexact routine from R
  fexact(&a,&y,&contingency_table_dbl[0],&leading,&expect,&percnt,&emin,&prt,&pre,&work);

  //check if any of the cell counts are zero
  unsigned zeros = std::count(_contingency_table->begin(),
			      _contingency_table->end(),
			      unsigned(0));

  //get the total number of counts in the table
  unsigned tcounts = accumulate(_contingency_table->begin(),
				_contingency_table->end(),unsigned(0));


  //Chi squared test
  double Csq,CsqYates,pCsq,pCsqYates;
  Csq = CsqYates = pCsq = pCsqYates = 0.;
  //G test
  double Gstat,williams,Gprob,GprobCorr;
  Gstat=williams=Gprob=GprobCorr=0.;
#if defined(HAVE_GSL)
  bool all_tests = true;
  if (untestable_table == false)
    {
      chisq2x2(&contingency_table_dbl[0],4,&Csq,&CsqYates,&pCsq,&pCsqYates);
      if(zeros==0)
	{
	  g2(&contingency_table_dbl[0],4,&Gstat,&williams,&Gprob,&GprobCorr);
	}
    }
#else
  bool all_tests = false;
#endif

  if (_args->parseableOutput == false)
    {
      s << "#Data from file "<<_args->infile<<":\n";
      s << "#\t" << 'A' << '\t' << 'S' << '\n'
	<< "#Fixed\t" << (*_contingency_table)[FIXEDA] << '\t' 
	<<  (*_contingency_table)[FIXEDS] << '\n'
	<< "#Poly\t" << (*_contingency_table)[POLYA] << '\t' 
	<<  (*_contingency_table)[POLYS] << "\n\n";      
      s << "#Statistics:\n";
      if ( untestable_table == false)
	{
	  s << "#Fisher's Exact test (p-value for table, p-value for fixed marginals):\n"
	    << pre << '\t'<< prt <<'\n';
	  if (all_tests)
	    {
	      if(tcounts>20)
		{
		  s << "#Chi-squared test (chi-squared, p-value):\n" 
		    << Csq <<'\t'<<pCsq<<'\n';
		}
	      else
		{
		  s << "#Chi-squared test, Yate's correction (chi-squared, p-value):\n" 
		    << CsqYates <<'\t'<<pCsqYates<<'\n';
		}
	      if(zeros==0)
		{
		  s << "#G test (G,p-value):\n"<< Gstat << '\t' << Gprob << '\n'; 
	      s << "#G test with William's correction (G,p-value):\n"<< Gstat/williams 
		<< '\t'<< GprobCorr <<'\n';
		}
	      else
		{
		  s << "#G-test not calculable for these data\n";
		}
	    }
	}
      else
	{
	  s << "Table untestable: there is at least 1 row filled with zero counts\n";
	}
    }
  else //parseable_args == true
    {
      s << _args->infile << '\t';
      copy(_contingency_table->begin(),_contingency_table->end(),
	   ostream_iterator<unsigned>(s,"\t"));
      if (untestable_table == false)
	{
	  s << pre
	    << '\t'
	    << prt;
	  if (all_tests)
	    {
	      s << '\t';
	      if (tcounts>20)
		{
		  s << Csq <<'\t'<<pCsq<<'\t';
		}
	      else
		{
		  s << CsqYates <<'\t'<<pCsqYates<<'\t';
		}
	      if(zeros==0)
		{
		  s << Gstat 
		    << '\t'
		    << Gprob 
		    << '\t'
		    << Gstat/williams 
		    << '\t' 
		    << GprobCorr;
		}
	      else
		{
		  s << "NA\tNA\tNA\tNA";
		}
	    }
	}
      else
	{
	  s << "untestable";
	}
    }
  return s;
}

