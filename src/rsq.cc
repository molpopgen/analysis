#include <iostream>
#include <vector>
#include <memory>
#include <cmath>
#include <ctime>
#include <Sequence/SimpleSNP.hpp>
#include <Sequence/PolySites.hpp>
#include <Sequence/PolySNP.hpp>
#include <Sequence/Correlations.hpp>
#if defined( __GNUG__ ) && __GNUC__ >= 3
#include <Sequence/FastaExplicit.hpp>
#else
#include <Sequence/Fasta.hpp>
#include <Sequence/Alignment.hpp>
#endif
#include <Sequence/Recombination.hpp>
#include <getopt.h>
#include <limits>
#include <algorithm>
#include <functional>

using namespace std;
using namespace Sequence;
using namespace Sequence::Alignment;

//simple wrapper for us to use
void shuffle_it( std::vector<double>::iterator beg,
		 std::vector<double>::iterator end,
		 std::mt19937 & generator )
{
  std::shuffle(beg,end,generator);
}

int main (int argc, char *argv[])
{
  if (argc == 1)
    {
      cerr << "Usage:\n";
      cerr << "rsq -i <fasta file>\n";
      cerr << "rsq -h <spreadsheet file>\n";
      cerr << " -c <mincount> : apply a frequency cutoff\n";
      cerr << " -t <int> : test significance of correlation of\n";
      cerr << " -m <double> : maximum distance b/w markers\n";
      cerr<< "\tLD with distance with <int> permutations\n";
      exit(1);
    }
  char *infile = NULL;
  int c;
  ifstream in;
  int mincount=1;
  bool isFasta = false;
  bool ish2k1 = false;
  bool test = false;
  double max_marker_distance=numeric_limits<double>::max();
  unsigned int nperms=10000;
  while ((c = getopt (argc, argv, "i:h:c:t:m:")) != -1)
    {
      switch (c)
        {
        case 'i':
          infile = optarg;
          isFasta = true;
          break;
        case 'h':
          infile = optarg;
          ish2k1 = true;
          break;
        case 'c':
          mincount = atoi(optarg);
          break;
        case 't':
          test = true;
          nperms=unsigned(atoi(optarg));
          break;
	case 'm':
	  max_marker_distance = atof(optarg);
	  break;
        }
    }
  vector<Sequence::Fasta > data;
  PolySites *p = NULL;
  auto_ptr<Hudson2001> h(new Hudson2001);
  bool h2k1 = false;
  if(isFasta)
    {
      try
        {
          GetData(data,infile);
          if(!IsAlignment(data))
            {
              cerr << infile <<" does not contain aligned sequences\n";
              exit(1);
            }
          if(Gapped(data))
            RemoveTerminalGaps(data);
          p = new PolySites(data);
        }
      catch (SeqException &e)
        {
          e.print(cerr);
          cerr << endl;
          exit(1);
        }
    }
  else if(ish2k1)
    {
      try
        {
          if (infile != NULL)
            {
              in.open (infile);
	      if (!in)
		{
		  std::cerr << "error:"
			    << infile
			    << " cannot be found\n";
		  exit(10);
		}
              in >> *h;
            }
          else
            cin >> *h;
        }
      catch (SeqException &e)
        {
          e.print(cerr);
          cerr << endl;
          exit(1);
        }
      p = new PolySites (h->GetPositions (), h->GetData ());
      h2k1 = true;
      
    }

  PolySNP *P;
  if( ! h2k1 )
    {
      P = new PolySNP (p,false,false);
    }
  else
    {
      if( h->outgroup() )
	{
	  P = new PolySNP (p,true,0);
	}
      else
	{
	  P = new PolySNP (p,false,false);
	}
    }
  int npoly = P->NumPoly();
  int numsing = P->NumSingletons();
  unsigned i=0,j=1;
  vector<double> LDSTATS(6);
  vector < vector<double> > diseq_vals;// =P->Disequilibrium(mincount,max_marker_distance);
  vector<double> distance,rsq,Dprime;

  cout << "sitei\tsitej\trsq\tD\tDprime\n";
  bool maintest=true;
  //while( Recombination::Disequilibrium(p,LDSTATS,&i,&j,false,0,mincount,max_marker_distance) )
  do
    {
      maintest = Recombination::Disequilibrium(p,LDSTATS,&i,&j,false,0,mincount,max_marker_distance);
      if(! LDSTATS[LDSTATS.size()-1] )
	{
	  if(fabs(LDSTATS[0]-LDSTATS[1]) <= max_marker_distance)
	    {
	      //cout << i << '\t' << j << '\t' 
	      cout << LDSTATS[0] << '\t' 
		   << LDSTATS[1] << '\t' 
		   << LDSTATS[2] << '\t' 
		   << LDSTATS[3] << '\t' 
		   << LDSTATS[4] << endl;
	      if(test)
		{
		  distance.push_back( abs(LDSTATS[0]-LDSTATS[1]) );
		  rsq.push_back(LDSTATS[2]);
		  Dprime.push_back(LDSTATS[4]);
		}
	    }
	}
    } while (maintest);
  delete P;
  if( test )
    {
      std::mt19937 generator(std::time(0));
      std::function<void(std::vector<double>::iterator,
			 std::vector<double>::iterator)> __x = [&generator](std::vector<double>::iterator  a,
									    std::vector<double>::iterator  b) { return shuffle_it(a,b,generator); }; 
      double obs = ProductMoment()(distance.begin(),distance.end(),rsq.begin());
      double p = (obs > 0) ? PermuteCorrelation(distance.begin(),distance.end(),rsq.begin(),
						ProductMoment(),std::greater_equal<double>(),
						__x,
						nperms) :
	PermuteCorrelation(distance.begin(),distance.end(),rsq.begin(),
			   ProductMoment(),std::less_equal<double>(),
			   __x,
			   nperms);

      cout <<"#Product moment correlation between distance and r^2 is " << obs;
      cout << ", p = " << p << '\n';
      obs = ProductMoment()(distance.begin(),distance.end(),Dprime.begin());
      p = (obs > 0) ? PermuteCorrelation(distance.begin(),distance.end(),Dprime.begin(),
					 ProductMoment(),std::greater_equal<double>(),
					 __x,
					 nperms) :
	PermuteCorrelation(distance.begin(),distance.end(),Dprime.begin(),
			   ProductMoment(),std::less_equal<double>(),
			   __x,
			   nperms);

      cout <<"#Product moment correlation between distance and D' is " << obs
	   <<", p = " << p << '\n';
    }
}
