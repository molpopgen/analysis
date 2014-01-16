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

#include <polydNdSbase.hpp>

using namespace std;
using namespace Sequence;
using namespace Sequence::Alignment;

#if defined(__GNUG__) && __GNUC__ >= 3
typedef std::ostringstream _ostr;
#else
typedef std::ostrstream _ostr;
#endif

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
          cerr << "error: file " <<args.infile<<" doesn't seem to be aligned...\n";
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

  if (args.XtractFullCodons == true)
    {
      switch (args.codon_start)
	{
	case 1:
	  break;
	case 2:
	  args.intervals[0]+=2;
	  break;
	case 3:
	  args.intervals[0]+=1;
	  break;
	}
      switch (args.codon_end)
	{
	case 1:
	  args.intervals[args.intervals.size()-1] -= 1;
	  break;
	case 2:
	  args.intervals[args.intervals.size()-1] -= 2;
	  break;
	case 3:
	  break;
	}
    }

  vector< Fasta > coding_region = Trim(data,args.intervals);

  for (unsigned i = 0 ; i < coding_region.size() ; ++ i)
    {
      _ostr new_name;
      new_name << coding_region[i].GetName()
	       << " (extracted CDS from file "
	       << args.infile
	       << ')';
      cout << Sequence::Fasta(new_name.str(),coding_region[i].GetSeq()) << endl;
    }

  cerr << "1 1" << ' ' << coding_region[0].length() << endl;
  if (args.XtractFullCodons == true)
    {
      cerr << "C 1\nE 3\n";
    }
  else
    {
      cerr << "C " << args.codon_start << endl;
      cerr << "E " << args.codon_end << endl;
    }
}

void usage(void)
{
}
