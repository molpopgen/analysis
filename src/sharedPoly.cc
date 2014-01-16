#include <iostream>
#include <fstream>
#include <vector>
#include <utility>
#include <stdexcept>
#include <algorithm>
#include <functional>
#include <iterator>
#include <getopt.h>

#if defined ( __GNUG__ ) && __GNUC__ >= 3
#include <Sequence/FastaExplicit.hpp>
#else
#include <Sequence/Fasta.hpp>
#include <Sequence/Alignment.hpp>
#endif

#include <Sequence/SeqExceptions.hpp>
#include <Sequence/PolySites.hpp>
#include <Sequence/SimpleSNP.hpp>
#include <Sequence/stateCounter.hpp>
#include <Sequence/FST.hpp>

enum dataType {SEQUENCE,TABLE};
enum siteType {FIXED,PRIVA,PRIVB,SHARED,UNKNOWN};

struct params
{
  char *infile;
  dataType type;
  unsigned n1;
  bool Verbose,biAllelic;
};

using namespace std;

void process(const Sequence::PolyTable *data, const params *args,const unsigned & sites_compared);
void streamVector( set<double>::const_iterator beg, 
		   set<double>::const_iterator end);
void parseargs(int argc, char **argv, params *args);
void usage(void);

int main(int argc,char **argv)
{
  params args;
  parseargs(argc,argv,&args);
  if(args.infile == NULL ||
     args.n1 == 0)
    {
      usage();
      exit(1);
    }

  Sequence::PolyTable *data = NULL;
  unsigned sites_compared = 0;
  try 
    {
      if (args.type == SEQUENCE)
	{
	  vector< Sequence::Fasta > seqs;
	  Sequence::Alignment::GetData(seqs,args.infile);
	  data = new Sequence::PolySites(seqs);
	  sites_compared = Sequence::Alignment::UnGappedLength(seqs);
	}
      else if (args.type == TABLE)
	{
	  data = new Sequence::Hudson2001;
	  ifstream in(args.infile);
	  in >> *data;
	  in.close();
	}
    }
  catch (Sequence::SeqException &e)
    {
      cerr << e << endl;
      exit(1);
    }
  catch (exception &std_exception)
    {
      cerr << std_exception.what() << endl;
      exit(1);
    }

  if (args.n1 >= data->size())
    {
      cerr << "error: n must be > 1 and < sample size\n";
      exit(1);
    }
  if (args.biAllelic == true)
    {
      data->RemoveMultiHits();
    }
  process(data,&args,sites_compared);

  delete data;
  exit(0);
}

void process(const Sequence::PolyTable *data, const params *args, const unsigned & sites_compared)
{
  unsigned config[2];
  config[0]=args->n1;
  config[1]=data->size()-args->n1;
  Sequence::FST analysis(data,2,config);
  std::set<double> shared = analysis.shared(0,1);
  std::set<double> fixed = analysis.fixed(0,1);
  std::pair<std::set<double>,std::set<double> > priv = analysis.Private(0,1);

  cout << "Statistics for " << args->infile << endl
       << "Sites compared " << sites_compared << endl
       << "# shared polymorphisms:\t"<< shared.size() << endl
       << "# fixed differences:\t"<< fixed.size() << endl//fixed.nMuts() << endl
       << "# privates in partition 1:\t"<< priv.first.size()<<endl//nMuts() <<endl
       << "# privates in partition 2:\t"<< priv.second.size()<< endl//nMuts() << endl
       << endl;

  if (args->Verbose == true)
    {
      cout << "Sites with fixed differences:\n";
      streamVector(fixed.begin(),fixed.end());
      cout << "Sites with shared polymorphisms:\n";
      streamVector(shared.begin(),shared.end());
      cout << "Sites with mutations private to partition 1:\n";
      streamVector(priv.first.begin(),priv.first.end());
      cout << "Sites with mutations private to partition 2:\n";
      streamVector(priv.second.begin(),priv.second.end());
    }
}

void streamVector( set<double>::const_iterator beg, 
		   set<double>::const_iterator end)
  /*!
    prints out data in space-delimited rows of 10
  */
{
  unsigned i=1;
  while (beg!=end)
    {
      cout << *beg;
      if (i%10 == 0.)
	cout << '\n';
      else
	cout << ' ';
      ++i;
      ++beg;
    }
  cout << endl;
}

void parseargs(int argc, char **argv, params *args)
{
  if(argc == 1)
    {
      usage();
      exit(1);
    }

  args->infile = NULL;
  args->type = SEQUENCE;
  args->n1 = 0;
  args->Verbose = false;
  args->biAllelic=false;
  int c;
  while ((c = getopt (argc, argv, "i:h:n:vb")) != -1)
    {
      switch (c)
	{
	case 'i':
	  args->infile = optarg;
	  args->type = SEQUENCE;
	  break;
	case 'h':
	  args->infile = optarg;
	  args->type = TABLE;
	  break;
	case 'n':
	  args->n1 = atoi(optarg);
	  break;
	case 'v':
	  args->Verbose = true;
	  break;
	case 'b':
	  args->biAllelic = true;
	  break;
	default:
	  usage();
	  exit(1);
	}
    }
}

void usage(void)
{
  cerr << "Usage:\n" 
       << "sharedPoly [options]\n"
       << "options include\n"
       << "-i"
       << '\t'
       << "infile (for sequence data in FASTA format)\n"
       << "-h"
       << '\t'
       << "infile (for SNP data in table format)\n"
       << "-n"
       << '\t'
       << "[integer] (sample size of the first partition in the data)\n"
       << "-v\t(verbose output)\n";
}
