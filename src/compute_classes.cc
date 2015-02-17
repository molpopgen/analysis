#include <compute_classes.hpp>

//libsequence headers
#if __GNUG__ && __GNUC__ >= 3
#include <Sequence/FastaExplicit.hpp>
#else
#include <Sequence/Fasta.hpp>
#include <Sequence/Alignment.hpp>
#endif
#include <Sequence/PolySites.hpp>
#include <Sequence/SimpleSNP.hpp>
#include <Sequence/PolySNP.hpp>
#include <Sequence/PolySIM.hpp>
#include <Sequence/SeqExceptions.hpp>
#include <Sequence/Coalescent/Coalescent.hpp>
#include <Sequence/SeqConstants.hpp>

//standard/unix headers
#include <iostream>
#include <sstream>
#include <fstream>
#include <exception>
#include <limits>
#include <random>
#include <ctime> 

#include <config.h>
using namespace std;
using namespace Sequence;
using namespace Alignment;
#ifdef HAVE_SEQUENCE_COALSIM
using namespace Sequence::coalsim;
#endif
namespace Sequence
{
  MAX_SEG_T MAX_SEGSITES = 200;
  MAX_SEG_T MAX_SEGS_INC = 100;
}

compute_params::compute_params() : infileglob(NULL),
				   outfile(NULL),
				   haveOutgroup(false),
				   outgroup(0),
				   useTotMuts(true),
				   is_table(false),
				   suppress_headers(false),
				   bi_allelic_only(false),
				   no_missing(false),
				   print_and_die(false),
				   probs(false),
				   verbose(false),
				   use_theta(false),
				   pretty(false)
{
}

struct resultsImpl
{
  bool pretty,wasSkipped,haveOutgroup;
  double thetaW,thetaPi,thetaH,tajd,dvh,fld,flf,rho87,wallB,wallQ;
  unsigned nsam,S,n1,nm,nsites,nsites_ug,dvk,rm,outgroup,nexternal;
  void calculate(const PolyTable * snptable,
		 const bool & tmuts);
  void removeFixedInsertions(vector<Fasta> * alignment,
			     unsigned site,
			     const unsigned & outgroup);
  unsigned columnWithChar(const vector<Fasta> & alignment,
			  const char & c);
  explicit resultsImpl (const bool & _p) : pretty(_p),
					   wasSkipped(false)
  {
  }
};

void resultsImpl::calculate(const PolyTable * snpTable,
			    const bool & tmuts)
{
  if (!wasSkipped)
    {
      PolySNP analyze(snpTable,haveOutgroup,outgroup,tmuts);
      
      thetaW = analyze.ThetaW();
      thetaPi = analyze.ThetaPi();
      tajd = analyze.TajimasD();
      dvk = analyze.DandVK();
      dvh = analyze.DandVH();
      rho87 = analyze.HudsonsC();
      wallB = analyze.WallsB();
      wallQ = analyze.WallsQ();
      S = analyze.NumPoly();
      n1 = analyze.NumSingletons();
      nm = analyze.NumMutations();
      rm = analyze.Minrec();
      
      if (haveOutgroup)
	{
	  thetaH = analyze.ThetaH();
	  nexternal = analyze.NumExternalMutations();
	}
      
      fld = (haveOutgroup==true) ? analyze.FuLiD() : analyze.FuLiDStar();
      flf = (haveOutgroup==true) ? analyze.FuLiF() : analyze.FuLiFStar();
    }
}

unsigned resultsImpl::columnWithChar(const vector<Fasta> & alignment,
				     const char & c)
{
  const unsigned ns = alignment.size();
  unsigned count=0;
  for(unsigned site=0;site<alignment[0].length();++site)
    {
      unsigned sc = 0;
      for(unsigned seq=0;seq<alignment.size();++seq)
	{
	  if (::toupper(alignment[seq][site]) == ::toupper(c))
	    ++sc;
	  else
	    seq = alignment.size(); //skip rest of this site
	}
      if(sc==ns) ++count;
    }
  return count;
}

void resultsImpl::removeFixedInsertions(vector<Fasta> * alignment,
					unsigned site,
					const unsigned & outgroup)
{
  const unsigned nsam = alignment->size()-1;
  int s = (*alignment)[0].length()-1;
  for(  ; s >= 0 ; --s )
    {
      unsigned ngap=0;
      for(unsigned ind=0;ind<alignment->size();++ind)
	{
 	  if (ind != outgroup)
	    {
	      ngap += ((*alignment)[ind][s] == '-') ? 1 : 0;
	    }
 	}
      if(ngap==nsam)
 	{
 	  for(unsigned ind=0;ind<alignment->size();++ind)
 	    {
 	      assert (s != (*alignment)[ind].length()-1);
 	      (*alignment)[ind].second.erase(s,1);
 	    }
	}
    }
}

results::results(const char * _infile,
		 compute_params * args)
  : impl(new resultsImpl(args->pretty)),
    infile(_infile)
{
  PolyTable * snptable = NULL;

  impl->haveOutgroup = args->haveOutgroup;
  impl->outgroup = args->outgroup;

  bool fileOK = true;
  if(args->verbose)
    cerr << "processing file " << infile << endl;
  if (args->is_table == false)
    {
      vector<Sequence::Fasta> data;//store the sequence objects

      ifstream in(infile);
      if (!in)
	{
	  fileOK = false;
	  cerr << infile
	       << " could not be opened\n";
	}
      else
	{
	  try 
	    {
	      GetData(data,infile);
	      //remove fixed insertions from the alignment
	      //so that the # sites used is appropriate for the ingroup
	      if (args->haveOutgroup == true)
		{
		  impl->removeFixedInsertions(&data,0,args->outgroup);
		}
	    }
	  catch (SeqException &b)
	    {
	      //badFormat is thrown if
	      //the input is not in the expected format,
	      //such as a FASTA file that does not begin
	      //with '>'
	      cerr <<"error: processing of file " << infile << " threw the following\n";
	      cerr <<"exception: ";
	      cerr << b;
	      cerr << endl;
	      fileOK = 0;
	    }
	  catch (exception & e)
	    {
	      cerr <<"error: processing of file " << infile 
		   << " threw the following exception:\n"
		   << e.what();
	      cerr << "This is an exception from the C++ library.\n"
		   << "This may occur if input is out of range,\n"
		   << "such as the number of the outgroup sequence in the data.\n";
	      fileOK=0;
	    }
	}
      if (! IsAlignment(data) )
	{
	  cerr << "data in file "
	       << infile
	       << " do not appear to be aligned, skipping...\n";
	  fileOK = false;
	}

      if(fileOK)
	{
	  impl->nsites = data[0].second.length();
	  const unsigned nmissingcols = impl->columnWithChar(data,'N');
	  impl->nsites -= nmissingcols;
	  impl->nsites_ug = UnGappedLength(data);
	  impl->nsites_ug -= nmissingcols;
	  snptable = new PolySites(data);
	}
    }
  else
    {
      snptable = new SimpleSNP;
      //set # sites to 1 to simplify output
      impl->nsites=impl->nsites_ug=1;
      ifstream in(infile);
      if (! in)
	{
	  fileOK = false;
	}
      else
	{
	  try {
	    in >> *snptable;
	  }
	  catch (SeqException &e)
	    {
	      cerr << "While processing "
		   << infile
		   << ", the following exception was thrown:\n"
		   << e << '\n';
	      fileOK = false;
	    }
	  catch (exception & e)
	    {
	      cerr << "While processing "
		   << infile
		   << ", the following exception was thrown:\n"
		   << e.what() << '\n';
	      fileOK = false;
	    }
	}

      //need downcast in order to use member fxn of derived type
      //yep, that's ugly
      if (static_cast<SimpleSNP*>(snptable)->outgroup() == true)
	{
	  impl->haveOutgroup=true;
	  args->haveOutgroup=true;
	  impl->outgroup=0;
	}
    }
  if (! validForPolyAnalysis(snptable->begin(),snptable->end()) )
    {
      cerr << infile
	   << " contained characters other than {A,G,C,T,N,-},\n"
	   << "which are not handled by this program.\n"
	   << "Sites containing such characters will be removed from\n"
	   << "the analysis.\n";
    }

  if (impl->haveOutgroup && 
      ( impl->outgroup>=snptable->size()))
    {
      cerr << "outgroup out of range for file "
	   << infile
	   << '\n';
      fileOK=false;
    }
  if (fileOK)
    {
      if (args->bi_allelic_only == true)
	snptable->RemoveMultiHits();
      if (args->no_missing == true)
	snptable->RemoveMissing();
      snptable->RemoveAmbiguous();
      impl->nsam = (impl->haveOutgroup==true) ? snptable->size()-1 : snptable->size();
      if (args->print_and_die)
	{
	  cout << *snptable << endl;
	  exit(0);
	}
      else
	impl->calculate(snptable,args->useTotMuts);
    }
  else
    {
      //file had to be skipped
      if (args->verbose)
	{
	  cerr << "skipping infile, please check data file\n";
	}
      impl->wasSkipped = true;
    }

  delete snptable;
}

std::ostream & results::print(std::ostream & o) const
{
  if (impl->pretty)
    {
      o << "file:\t" << infile << '\n'
	<< "n:\t" << impl->nsam << '\n'
	<< "number of sites:\t" << impl->nsites << '\n'
	<< "number of sites (w/o gaps):\t" << impl->nsites_ug << '\n'
	<< "number of segregating sites:\t" << impl->S << '\n'
	<< "number of singletons:\t" << impl->n1 << '\n';
      if(impl->haveOutgroup)
	o << "number of derived singletons:\t" << impl->nexternal << '\n';
      o << "number of mutations:\t" << impl->nm << '\n'
	<< "number of haplotypes:\t" << impl->dvk << '\n'
	<< "haplotype diversity:\t" << impl->dvh << '\n';
      if(impl->S >= 2)
	{
	  o << "Wall's B:\t"<< impl->wallB << '\n'
	    << "Wall's Q:\t" << impl->wallQ << '\n';
	}
      else
	{
	  o << "Wall's B:\tNA (< 2 segregating sites)\n"
	    << "Wall's Q:\tNA (< 2 segregating sites)\n";
	}
      o << "Watterson's Theta per site:\t" << impl->thetaW/double(impl->nsites_ug) << '\n'
	<< "Pi per site:\t" << impl->thetaPi/double(impl->nsites_ug) << '\n';
      if (impl->haveOutgroup)
	{
	  o << "Theta from site homozygosity (i.e. thetaH)\t"
	    << impl->thetaH/double(impl->nsites_ug) << '\n';
	}
      if (impl->S > 0)
	{
	  o << "Tajima's D:\t"<< impl->tajd << '\n';
	  if (impl->haveOutgroup)
	    {
	      o << "Fu & Li D:\t" << impl->fld << '\n'
		<< "Fu & Li F:\t" << impl->flf << '\n';
	    }
	  else
	    {
	      o << "Fu & Li D*:\t" << impl->fld << '\n'
		<< "Fu & Li F*:\t" << impl->flf << '\n';
	    }
	  o << "min. number of recombination events:\t";
	  if ( impl->rm != SEQMAXUNSIGNED )
	    o << impl->rm << '\n';
	  else
	    o << strtod("NAN",NULL) << '\n';
	  o << "Hudson's (1987) C:\t" << impl->rho87;
	}
      else
	{
	  o << "Tajima's D:\tNA\n";
	  if (impl->haveOutgroup)
	    {
	      o << "Fu & Li D:\tNA\n"
		<< "Fu & Li F:\tNA\n";
	    }
	  else
	    {
	      o << "Fu & Li D*:\tNA\n"
		<< "Fu & Li F*:\tNA\n";
	    }
	  o << "min. number of recombination events:\tNA\n"
	    << "Hudson's (1987) C:\tNA";
	}
    }
  else
    {
      o << infile << '\t'
	<< impl->nsam << '\t'
	<< impl->nsites << '\t'
	<< impl->nsites_ug << '\t';
      if (impl->wasSkipped == false)
	{
	  o << impl->S << '\t'
	    << impl->n1 << '\t';
	  if(impl->haveOutgroup)
	    o << impl->nexternal << '\t';
	  o << impl->nm << '\t'
	    << impl->dvk << '\t'
	    << impl->dvh << '\t';
	  if (impl->S >= 2)
	    {
	      o << impl->wallB << '\t'
		<< impl->wallQ << '\t';
	    }
	  else
	    {
	      o << "NA\tNA\t";
	    }
	  o << impl->thetaW/double(impl->nsites_ug) << '\t'
	    << impl->thetaPi/double(impl->nsites_ug) << '\t';
	  if (impl->haveOutgroup)
	    o << impl->thetaH/double(impl->nsites_ug) << '\t';
	  if(impl->S > 0)
	    {
	      o << impl->tajd << '\t'
		<< impl->fld << '\t'
		<< impl->flf << '\t'
		<< impl->rm << '\t'
		<< impl->rho87;
	    }
	  else
	    {
	      o << "NA\tNA\tNA\tNA\tNA";
	    }
	}
      else //file was skipped, so print a bunch of NA
	{
	  for(unsigned i = 0 ; i < 13 ; ++i)
	    {
	      o << "NA\t";
	    }
	  if (impl->haveOutgroup)
	    {
	      o << "NA\tNA";
	    }
	  else
	    o << "NA";
	}
    }
  return o;
}

results::~results()
{
  delete impl;
}

struct pvalsImpl
{
  double ptajd, pfld,pflf,pdvk,pdvh,
    pfwh,pb,pq;
  bool pretty,haveOutgroup;
  mutable double NA,DEPS;
  std::string stat_as_string(const double & d);
  string as_string(const double & d);

  explicit pvalsImpl(const bool & _p,
		     const bool &_ho) 
    : ptajd(NA), pfld(NA),pflf(NA),pdvk(NA),pdvh(NA),
      pfwh(NA),pb(NA),pq(NA), pretty(_p), haveOutgroup(_ho)
  {
    DEPS = std::numeric_limits<double>::epsilon();
    NA = std::numeric_limits<double>::min();
  }
};

string pvalsImpl::as_string(const double & d)
{
  typedef ostringstream _ostr;
  _ostr o;
  o<<d;
  return o.str();
}

std::string  pvalsImpl::stat_as_string(const double & d)
{
  return ( (d != NA) ? as_string(d) : "NA" );
}

pvals::pvals(const results & r,
	     const compute_params * args)
  : impl(new pvalsImpl(args->pretty,r.impl->haveOutgroup))
{
  if (!r.impl->wasSkipped && r.impl->S > 0)
    {
      //initialize pvals from NA to 0
      impl->ptajd =  impl->pfld = impl->pflf = impl->pdvk = impl->pdvh = 
	impl->pfwh = 0.;
      if (r.impl->S>1)
	{
	  //statistic is only defined for data sets with S>1,
	  //otherwise pvalue is not applicable
	  impl->pb=impl->pq=0.;
	}
      double theta = (args->use_theta) ? r.impl->thetaW : 0.;
      unsigned S = 0;
      if (theta == 0.)
	{
	  S = (args->useTotMuts) ? r.impl->nm : r.impl->S;
	}

      std::mt19937 generator(std::time(0));
      std::uniform_real_distribution<double> uni01(0.,1.);
      std::function<double(const double&,const double&)> uni = [&generator](const double & a, const double & b){ return std::uniform_real_distribution<double>(a,b)(generator); };
      std::function<double(const double&)> expo = [&generator](const double & mean){ return std::exponential_distribution<double>(1./mean)(generator); };
      std::function<double(const double&)> poiss = [&generator](const double & mean){ return std::poisson_distribution<unsigned>(mean)(generator); };

      vector<chromosome> initialized_sample = init_sample( std::vector<int>(1,r.impl->nsam),1 );
      marginal initialized_marginal = init_marginal(r.impl->nsam);
      const unsigned NRUNS = 1000;
      int nlinks = 0;
      unsigned nruns = NRUNS;
      while(nruns--)
	{
	  vector<chromosome> sample(initialized_sample);
#ifdef HAVE_SEQUENCE_COALSIM
	  Sequence::coalsim::arg history(1,initialized_marginal);
#else
	  Sequence::arg history(1,initialized_marginal);
#endif
	  int NSAM = r.impl->nsam;
	  double t = 0.;
	  while(NSAM>1)
	    {
	      double rcoal = double(NSAM*(NSAM-1));
	      t += expo( 1./rcoal );
	      pair<int,int> two = pick2(uni,NSAM);
	      NSAM -= coalesce(t,r.impl->nsam,NSAM,two.first,two.second,1,
			       &nlinks,&sample,&history); 
	    }
	  SimData simsample;
	  if(args->use_theta)
	    {
	      simsample = infinite_sites_sim_data(poiss,uni,1,history,theta);
	    }
	  else
	    {
	      double ttime = total_time(history.begin()->begin(),history.begin()->nsam);
	      simsample = infinite_sites_sim_data(uni,1,history,&ttime,&S);
	    }
	  PolySIM a(&simsample);
	  impl->ptajd += (a.TajimasD() <= r.impl->tajd) ? 1. : 0.;
	  impl->pdvk += (a.DandVK() <= r.impl->dvk) ? 1. : 0.;
	  impl->pdvh += (a.DandVH() <= r.impl->dvh) ? 1. : 0.;
	  if(r.impl->S > 1)
	    {
	      impl->pb += (a.WallsB() >= r.impl->wallB) ? 1. : 0.;
	      impl->pq += (a.WallsQ() >= r.impl->wallQ) ? 1. : 0.;
	    }
	  if (args->haveOutgroup)
	    {
	      impl->pflf += (a.FuLiF() <= r.impl->flf) ? 1. : 0.;
	      impl->pfld += (a.FuLiD() <= r.impl->fld) ? 1. : 0.;
	      impl->pfwh += (a.ThetaPi()-a.ThetaH() <= r.impl->thetaPi-r.impl->thetaH) ? 1. : 0.;
	    }
	  else
	    {
	      impl->pflf += (a.FuLiFStar() <= r.impl->flf) ? 1. : 0.;
	      impl->pfld += (a.FuLiDStar() <= r.impl->fld) ? 1. : 0.;
	    }
	}
      impl->ptajd /= double(NRUNS);
      impl->pdvk /= double(NRUNS);
      impl->pdvh /= double(NRUNS);
      if (r.impl->S > 1)
	{
	  impl->pb /= double(NRUNS);
	  impl->pq /= double(NRUNS);
	}
      impl->pflf /= double(NRUNS);
      impl->pfld /= double(NRUNS);
      impl->pfwh /= double(NRUNS);
    }
}


bool results::wasSkipped() const
{
  return impl->wasSkipped;
}

std::ostream & pvals::print(std::ostream & o) const
{
  if ( impl->pretty )
    {
      o << "P(Taj D <= Taj D obs.):\t" 
	<< impl->stat_as_string(impl->ptajd)
	<< '\n'
	<< "P(Num Haplotypes <= Num Haplotypes Obs.):\t" 
	<< impl->stat_as_string(impl->pdvk)
	<< '\n'
	<< "P(Haplotype Diversity <= Haplotype Diversity Obs.):\t" 
	<< impl->stat_as_string(impl->pdvh)
	<< '\n'
	<< "P(Wall's B >= Wall's B Obs.):\t" 
	<< impl->stat_as_string(impl->pb)
	<< '\n'
	<< "P(Wall's Q >= Wall's Q Obs.):\t" 
	<< impl->stat_as_string(impl->pq) << '\n';
      if (impl->haveOutgroup)
	{
	  o << "P(Fu & Li D <= Fu & Li D Obs.):\t" 
	    << impl->stat_as_string(impl->pfld)
	    << '\n'
	    << "P(Fu & Li F <= Fu & Li F Obs.):\t" 
	    << impl->stat_as_string(impl->pflf)
	    << '\n'
	    << "P(Fay-Wu H <= Fay-Wu H Obs.):\t" 
	    << impl->stat_as_string(impl->pfwh);
	}
      else
	{
	  o << "P(Fu & Li D* <= Fu & Li D* Obs.):\t" 
	    << impl->stat_as_string(impl->pfld)
	    << '\n'
	    << "P(Fu & Li F* <= Fu & Li F* Obs.):\t" 
	    << impl->stat_as_string(impl->pflf)
	    << '\n';
	}
    }
  else
    {
      o << impl->stat_as_string(impl->ptajd) << '\t'
	<< impl->stat_as_string(impl->pdvk) << '\t'
	<< impl->stat_as_string(impl->pdvh) << '\t'
	<< impl->stat_as_string(impl->pb) << '\t'
	<< impl->stat_as_string(impl->pq) << '\t'
	<< impl->stat_as_string(impl->pfld) << '\t'
	<< impl->stat_as_string(impl->pflf) << '\t';
      if (impl->haveOutgroup)
	{
	  o << impl->stat_as_string(impl->pfwh);
	}
    }
  return o;
}

pvals::~pvals() { delete impl; }
