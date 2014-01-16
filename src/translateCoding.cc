#include <vector>
#include <string>
#include <algorithm>
#include <iostream>
#if defined(__GNUG__) && __GNUC__ < 3
#include <strstream>
#else
#include <sstream>
#endif

#if defined( __GNUG__ ) && __GNUC__ >= 3
#include <Sequence/FastaExplicit.hpp>
#else
#include <Sequence/Fasta.hpp>
#include <Sequence/Alignment.hpp>
#endif

#include <Sequence/Translate.hpp>
#include <polydNdSbase.hpp>

using namespace std;
using namespace Sequence;
using namespace Sequence::Alignment;

void usage(void);

void (*usg)(void) = &usage;

int main (int argc, char *argv[])
{
  params args;
  try
    {
      parseargs(argc,argv,&args,usg);
    }
  catch (std::exception &e)
    {
      cerr << e.what() << endl;
      exit(1);
    }
  vector<Fasta > data;
  try
    {
      GetData(data,args.infile);
      if(!IsAlignment(data))
        {
          cerr << "error: data don't seem to be aligned...\n";
          exit(1);
        }
    }
  catch(SeqException &e)
    {
      cerr << "Processing file "<<args.infile<<" resulted in an exception being thrown:"<<endl;
      e.print(cerr);
      cerr <<endl;
      exit(1);
    }

  //reduce intervals vector to just the
  //translatable portion
  size_t radjust=0;
  vector<Fasta> coding_region;
  if (args.intervals.size() > 0)
    {
      if(args.codon_start > 1)
	{
	  switch (args.codon_start)
	    {
	    case 2:
	      //          ladjust = 2;
	      args.intervals[0]+=2;
	      break;
	    case 3:
	      args.intervals[0]+=1;
	      //          ladjust = 1;
	      break;
	    }
	}
      if (args.codon_end < 3)
	{
	  switch (args.codon_end)
	    {
	    case 1:
	      args.intervals[args.intervals.size()-1]-=1;
	      radjust = 1;
	      break;
	    case 2:
	      args.intervals[args.intervals.size()-1]-=2;
	      radjust=2;
	      break;
	    }
	}
      coding_region = Trim(data,args.intervals);
    }
  else 
    {
      coding_region=data;
    }
  vector< Fasta > translations(coding_region.size());
  string::size_type pos;
  for(unsigned i = 0 ; i < coding_region.size() ; ++i)
    {
      //remove the gaps from each sequence individually
      while((pos=coding_region[i].second.find('-'))!=string::npos)
	{
	  coding_region[i].second.erase(pos,1);
	}
      
      std::string translation = Sequence::Translate(coding_region[i].begin(),
						    coding_region[i].end());
#if defined(__GNUG__) && __GNUC__ >= 3
      std::ostringstream name;
#else
      std::ostrstream name;
#endif
      name << coding_region[i].GetName() << " (translated CDS";
      if (std::find(translation.begin(),translation.end(),'*') != translation.end())
	{
	  name << ", stop codons found";
	}
      name << ')';
      //use cast to std::string for compatability with both ostringstream
      //and the (deprecated) ostrstream
      translations[i] = Sequence::Fasta(std::string(name.str()),translation);
    }

  for(unsigned i = 0 ; i < translations.size() ; ++i)
    {
      cout << translations[i] << endl;
    }
}

void usage(void)
{
}
