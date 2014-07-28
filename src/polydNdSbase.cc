#include <polydNdSbase.hpp>
#include <cassert>
#include <iostream>
#include <fstream>
#include <cctype>
#include <stdexcept>
#include <sstream>
#include <getopt.h>

//using declarations are given explicitly here because the code 
//is complicated, so this should help document where certain
//functions come from
using std::vector;
using std::string;
using std::cerr;
using std::endl;
using std::ifstream;


typedef std::ostringstream _ostr;

bool InCoding(int site_index, const vector<int> &intervals)
{
  for(unsigned int p = 0 ; p < intervals.size() ; p+=2)
    {
      if ( site_index >= intervals[p] && site_index <= intervals[p+1] )
        {
          return true;
        }
    }
  return false;
}

int GetCodonPos(int site_index, const vector<int> &intervals)
{
  int len = 0;
  for(unsigned int p = 0 ; p < intervals.size() ; p += 2)
    {
      if ( site_index >= intervals[p] && site_index <= intervals[p+1] )
	{
	  len += site_index - intervals[p] + 1;
	  switch (int(len % 3))
	    {
	    case 1: //first position
	      return 0;
	      break;
	    case 2: //second position
	      return 1;
	      break;
	    case 0: //thirs position
	      return 2;
	      break;
	    }
	}
      else
        len += intervals[p+1]-intervals[p]+1;
    }
  switch (int(len % 3))
    {
    case 1:
      return 0;
      break;
    case 2:
      return 1;
      break;
    case 0:
      return 2;
      break;
    }
  return -1;//error case
}

void MakeCodon(int site_index,int seq, int codon_pos,
               const vector<Fasta >&data,
               const vector<int> &intervals,string &codon)
{
  codon.resize(3);
  unsigned int d_left = codon_pos;
  unsigned int d_right = 2-codon_pos;

  assert (d_left<3 && d_right<3);
  unsigned int j,k;
  for(j=site_index,k=0 ; k<=d_left ; --j)
    {
      if(InCoding(j,intervals))
        {
          codon[codon_pos-k] = data[seq][j];
          ++k;
        }
      else
	{}
    }

  for(j=site_index,k=0 ; j<data[0].length()&&k<=d_right ; ++j)
    {
      if(InCoding(j,intervals))
        {
          codon[codon_pos+k] = data[seq][j];
	  ++k;
        }
    }
}

void parseargs(int argc, char *argv[],params *args,USG usage)
  throw (std::exception)
{
  if(argc==1)
    {
      usage();
      exit(1);
    }

  //assign some defaults
  args->infile = NULL;
  args->haveOutgroup = 0;
  args->useTsTv=0;
  args->outgroup = 0;
  args->print_tables = false;
  args->skipMissing = false;
  args->SuperStrict = false;
  args->parseableOutput = false;//make output less verbose
  args->XtractFullCodons=false;
  args->ApproximateTreatmentOfDivergentCodons=false;
  args->codon_start = 1;
  args->codon_end = 3;
  args->intervals.resize(0);
  args->n1 = 0;
  extern int optind;
  int c;
  int num_intervals = 0,i;
  char *posfile=NULL;

  //LEOPARD PROBLEM?
  while ((c = getopt (argc, argv, "i:O:I:C:E:F:n:kPNSTXA")) != -1)
    {
      switch (c)
        {
        case 'i':
          args->infile = optarg;
          break;
        case 'I':
          num_intervals = atoi(optarg);
          args->intervals.resize(2*num_intervals);
          for( i = 0 ; i < 2*num_intervals ; i+=2)
            {
              //note that these are real positions, not array indices,
              //i.e. they count from one
              args->intervals[i] = atoi(argv[optind+i])-1;//subtract 1 to turn it into an index
              args->intervals[i+1] = atoi(argv[optind+i+1])-1;
            }
	  optind += i;
          break;
        case 'F':
          posfile = optarg;
          break;
        case 'O':
          args->haveOutgroup=1;
          args->outgroup=atoi(optarg)-1;//turn it into an array index
          break;
	case 'A':
	  args->ApproximateTreatmentOfDivergentCodons=true;
	  break;
        case 'C':
          args->codon_start = atoi(optarg);
          break;
        case 'E':
          args->codon_end = atoi(optarg);
          break;
        case 'k':
          args->useTsTv=1;
          break;
        case 'P':
          args->print_tables=true;
          break;
	case 'n':
	  args->n1 = atoi(optarg);
	  break;
	case 'N':
	  args->skipMissing = true;
	  break;
	case 'S':
	  args->SuperStrict = true;
	  break;
	case 'T':
	  args->parseableOutput = true;
	  break;
	case 'X':
	  args->XtractFullCodons = true;
	  break;
        case '?':
          cerr<<"error: unknown argument "<< argv[optind-1]<<endl;
          usage();
          exit(0);
          break;
        }
    }
  if(posfile != NULL)
    {
      ifstream in(posfile);
      if (! in)
	{
	  _ostr o;
	  o << "error : positions file "
	    << posfile
	    << " cannot be opened";
	  throw (std::runtime_error(o.str()));
	}
      in >> num_intervals;

      args->intervals.resize(2*num_intervals);
      for(int i = 0 ; i < 2*num_intervals ; i+=2)
        {
          //note that these are real positions, not array indices,
          //i.e. they count from one
          in >>args->intervals[i];
          in >>args->intervals[i+1];
        }
      for(unsigned i = 0 ; i < args->intervals.size() ; ++i)
	{
	  args->intervals[i] -= 1;//turn into an index
	}
      while (!(in.eof()))
        {
          char ch;
          in >> ch;
          ch = toupper(ch);
          switch (ch)
            {
            case 'E':
              in >> args->codon_end;
              break;
            case 'C':
              in >> args->codon_start;
              break;
            }
        }
      in.close();
    }
  int len=0;
	if(args->intervals.empty())
	{
		ifstream in(args->infile);
		Fasta f;
		in >> f;
		args->intervals.push_back(0);
		args->intervals.push_back(f.length()-1);
		in.close();
	}
  vector<int> interval_copy = args->intervals;
  switch(args->codon_start)
    {
    case 1:
      break;
    case 2:
      interval_copy[0] += 2;
      break;
    case 3:
      interval_copy[0] += 1;
      break;
    default:
      _ostr o;
      o << "error: invalid codon start position entered"<<endl;
      o << "valid values are 1,2, and 3"<<endl;
      throw (std::runtime_error(o.str()));
      break;
    }
  switch(args->codon_end)
    {
    case 1:
      interval_copy[interval_copy.size()-1] -= 1;
      break;
    case 2:
      interval_copy[interval_copy.size()-1] -= 2;
      break;
    case 3:
      break;
    default:
      _ostr o;
      o << "error: invalid codon end position entered"<<endl;
      o << "valid values are 1,2, and 3"<<endl;
      throw (std::runtime_error(o.str()));
      break;
    }

  for(unsigned int i = 0 ; i < interval_copy.size() ; i+=2)
    {
      len += (interval_copy[i+1]-interval_copy[i]+1);
    }

  //assert(len % 3 == 0.0);
  if (args->SuperStrict && len % 3 != 0.)
    {
      _ostr o;
      o << "Error: the coding positions for file "
	<< args->infile
	<< " sum to "<<len<<", which is not a multiple of 3.\n"
	<< "There are "<<interval_copy.size()<<" intervals specified.\n";
      for (unsigned i = 0 ; i < interval_copy.size() ; ++i)
        {
          o << interval_copy[i]<<'\n';
        }
      throw (std::runtime_error(o.str()));
    }
}
