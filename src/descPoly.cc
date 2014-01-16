#include <iostream>
#include <fstream>
#include <vector>
#include <utility>
#include <algorithm>
#include <functional>
#include <stdexcept>
#include <getopt.h>

#include <Sequence/SeqExceptions.hpp>
#include <polyCommon.hpp>

struct params
{
  char *infile;
  bool haveOutgroup;
  unsigned outgroup;
  bool Verbose;
};

using namespace std;
using namespace Sequence;
using namespace Sequence::Alignment;

//prototypes
void parseargs(int argc, char **argv, params *args);
void outputMulti(const vector<stateCounter> &counts);
void outputMissing(const vector<Fasta> &data, const vector< vpuu > & gapInfo);
void outputFreqSpec(const vpuu &freqSpec, const params *args);
void outputAdjacentPolyInfo( const vpuu &adjacents );
void outputGapDist(const vector<Fasta> &data, const vector< vpuu > & gapInfo);
void outputGapSNP(const vector<stateCounter> &counts);
void usage(void);


int main(int argc, char **argv)
{
  params args;
  parseargs(argc,argv,&args);

  vector<Fasta> data;
  try
    {
      GetData(data,args.infile);
      if(! IsAlignment(data))
	{
	  cerr << "Error: data file "<<args.infile<< " does not appear to be aligned.\n";
	  exit(1);
	}
    }
  catch (SeqException &SE)
    {
      cerr << SE << endl;
      exit (1);
    }
  catch (exception &e)
    {
      cerr << e.what() << endl;
      exit(1);
    }

  unsigned datalen = data[0].length();
  unsigned nsam = data.size();

  cout << "#Summary of sequence file: "<< args.infile << endl;
  cout << "#n = "<<nsam<<endl;
  cout << "#Alignment length = "<<datalen<<endl;
  cout << "#Alignment length (excluding gaps) = "<<UnGappedLength(data)<<endl;
  cout << endl;

  //get site by site state info
  vector<stateCounter> stateCounts  = processStateInfo(data,datalen,args.haveOutgroup,args.outgroup);

  //get sequence by sequence gap info
  vector < vpuu > gapInfo(nsam);
  unsigned gPos;
  bool gapsFound = false;
  for (unsigned i=0;i<nsam;++i)
    {
      std::string::iterator beg = data[i].begin(),pos;
      while( (pos = std::find(beg,data[i].end(),'-' )) != data[i].end() )
	{
	  gapsFound = true;
	  gPos = pos-data[i].begin();//array index of the gap character
	  unsigned len = 1;
	  if (gPos < datalen - 1)
	    {
	      for(unsigned j = gPos + 1 ; j < datalen ; ++j)//find out how long the gap is
		{
		  if( data[i][j] == '-' )
		    {
		      ++len;
		    }
		  else
		    {
		      j = datalen; //exit the loop
		    }
		}
	    }
	  gapInfo[i].push_back( puu(pos-data[i].begin()+1,len) );
	  beg=pos + len;//begin the search again at index pos+len
	}
    }
  if (gapsFound == false)
    gapInfo.clear();

  //get sequence by sequence info on missing data
  vector < vpuu > missingInfo(nsam);
  bool missingFound = false;
  for (unsigned i=0;i<nsam;++i)
    {
      std::string::iterator beg = data[i].begin(),pos;
      while( (pos = std::find(beg,data[i].end(),'N' )) != data[i].end() )
	{
	  missingFound = true;
	  gPos = pos-data[i].begin();//array index of the gap character
	  unsigned len = 1;
	  if (gPos < datalen - 1)
	    {
	      for(unsigned j = gPos + 1 ; j < datalen ; ++j)//find out how long the gap is
		{
		  if( data[i][j] == 'N' )
		    {
		      ++len;
		    }
		  else
		    {
		      j = datalen; //exit the loop
		    }
		}
	    }
	  missingInfo[i].push_back( puu(pos-data[i].begin()+1,len) );
	  beg=pos + len;//begin the search again at index pos+len
	}
    }
  if(missingFound == false)
    missingInfo.clear();

  //frequency spectrum info
  outputFreqSpec( getFreqSpectrum(stateCounts,nsam-args.haveOutgroup,args.haveOutgroup),&args );

  //output info on multiple hits
  outputMulti(stateCounts);

  //info on adjacent polymorphisms
  outputAdjacentPolyInfo( checkAdjacentPoly(data,args.haveOutgroup,args.outgroup) );

  //output position and lengths of gaps:
  outputGapDist(data,gapInfo);

  //info on SNPs that overlap gaps
  outputGapSNP(stateCounts);

  //info on missing data
  outputMissing(data,missingInfo);

  exit(0);
}

void outputMulti(const vector<stateCounter> &counts)
{
  vector<stateCounter>::const_iterator vsCi,vsCiBeg;

  //output info on sites with missing data and multiple hits
  bool haveMulti = ( (vsCi = find_if(counts.begin(),counts.end(),multiHit())) 
		       != counts.end() );
 
  if (haveMulti == true)
    {
      cout << "#There are more than 2 alleles present at the following sites\n"
	   <<"#(excluding gaps):\n";
      cout <<"#site\t(states segregating)\n";
      vsCiBeg = counts.begin();
      while ( (vsCi = find_if(vsCiBeg,counts.end(),multiHit()))
	      != counts.end() )
	{
	  //need to add 1 to these values because
	  //they are differences between addresses, so
	  //if a multiple hit (in this case) occurs
	  //at the very first position in the data,
	  //the value of vsCi-counts.begin() is zero,
	  //which we convert from an index to a position
	  //by adding one
	  cout << vsCi - counts.begin() + 1 << '\t';
	  cout << '(';
	  if ( (*vsCi).a > 0 )
	    {
	      cout << " A ";
	    }
	  if ( (*vsCi).g > 0 )
	    {
	      cout << " G ";
	    }
	  if ( (*vsCi).c > 0 )
	    {
	      cout << " C ";
	    }
	  if ( (*vsCi).t > 0 )
	    {
	      cout << " T ";
	    }
	  cout << ')' << endl;
	  vsCiBeg = vsCi+1;
	}
    }
  else
    {
      cout << "#All polymorphic sites (excluding those with gaps) are biallelic\n";
    }

  cout << endl;
}

void outputMissing(const vector<Fasta> &data, const vector< vpuu > &missingInfo)
{
  if (missingInfo.empty() == true)
    {
      cout << "#There are no sites with missing data in the alignment.\n";
    }
  else
    {
      cout << "#Distribution of missing data:\n"
	   << "#Sequence\t (start, length) ... \n";
      for (unsigned i = 0 ; i<missingInfo.size(); ++i)
	{
	  cout << data[i].GetName() << '\t';
	  for (unsigned j = 0 ; j<missingInfo[i].size(); ++j)
	    {
	      cout << '(' << missingInfo[i][j].first << ',' << missingInfo[i][j].second << ") ";
	    }
	  cout << endl;
	}
    }
  cout << endl;
}

void outputFreqSpec(const vpuu &freqSpec, const params *args)
{
  cout << "#Frequency Spectrum From Biallelc Sites:\n#";
  if (args->haveOutgroup == true)
    {
      cout << "(based on derived state)";
    }
  else
    {
      cout << "(based on folded site-frequency spectrum)";
    }
  cout <<endl << "#(excluding sites with gaps)" << endl;
  cout << "#Freq\tNum\n";
  for(unsigned i = 0 ; i < freqSpec.size() ; ++i)
    {
      cout << freqSpec[i].first << '\t' << freqSpec[i].second << '\n';
    }

  cout << endl;
}

void outputAdjacentPolyInfo( const vpuu &adjacents )
{
  if (adjacents.empty())
    {
      cout << "#No adjacent polymorphisms were found on the same background "
	   << "at biallelic sites.\n";
    }
  else
    {
      cout <<"#The following is a list of adjacent polymorphic sites where\n"
	   <<"#the minor allele occurs on the same haplotype at least once at each site.\n";
      cout << "#These sites may deserve manual inspection:\n";
      cout << "#site_i\tsite_j\n";
      for(unsigned i = 0 ; i < adjacents.size() ; ++i)
	{
	  //no need to add 1 since these are assigned from a Sequence::PolySites
	  //object, which already adds the 1 when converting from sequence data
	  //to the polymorphism table
	  cout << adjacents[i].first << '\t' << adjacents[i].second << endl;
	}
    }

  cout << endl;
}

void outputGapDist(const vector<Fasta> &data, const vector< vpuu > & gapInfo)
{
  if (gapInfo.empty() == true)
    {
      cout << "#There are no gaps in the alignment.\n";
    }
  else
    {
      cout << "#Distribution of alignment gaps:\n"
	   << "#Sequence\t (gap start, length) ... \n";
      for (unsigned i = 0 ; i<gapInfo.size(); ++i)
	{
	  cout << data[i].GetName() << '\t';
	  for (unsigned j = 0 ; j<gapInfo[i].size(); ++j)
	    {
	      cout << '(' << gapInfo[i][j].first << ',' << gapInfo[i][j].second << ") ";
	    }
	  cout << endl;
	}
    }
  cout << endl;
}

void outputGapSNP(const vector<stateCounter> &counts)
{
  vector<stateCounter>::const_iterator vsCi,vsCiBeg;
  bool haveSNPsOverGaps = ( (vsCi = find_if(counts.begin(),counts.end(),gapOverSNP()))
			    != counts.end() );
  
  if (haveSNPsOverGaps == false)
    {
      cout << "#No SNPs were found at positions that contain alignment gaps.\n";
    }
  else
    {
      cout << "#The following is a list of sites at which a SNP is present\n"
	   <<"#at a site containing an alignment gap:\n";
      cout << "#site\t(character states)\n";
      vsCiBeg = counts.begin();
      while ( (vsCi = find_if(vsCiBeg,counts.end(),gapOverSNP()))
			    != counts.end() )
	{
	  cout << vsCi - counts.begin() + 1 << '\t';
	  cout << '(';
	  if ( (*vsCi).a > 0 )
	    {
	      cout << " A ";
	    }
	  if ( (*vsCi).g > 0 )
	    {
	      cout << " G ";
	    }
	  if ( (*vsCi).c > 0 )
	    {
	      cout << " C ";
	    }
	  if ( (*vsCi).t > 0 )
	    {
	      cout << " T ";
	    }
	  cout << ')' << endl;
	  vsCiBeg = vsCi+1;
	}
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
  args->haveOutgroup = false;
  args->outgroup = 0;
  args->Verbose = false;
  int c;
  while ((c = getopt (argc, argv, ":i:O:")) != -1)
    {
      switch (c)
        {
        case 'i':
          args->infile = optarg;
          break;
	case 'O':
	  args->haveOutgroup = true;
	  if (atoi(optarg) != ':')
	    {
	      args->outgroup = atoi(optarg)-1;
	    }
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
       << "descPoly [options]\n"
       << "options include\n"
       << "-i"
       << '\t'
       << "infile (for sequence data in FASTA format)\n"
       << "-o [integer]\t"
       << "specify the sequence (i.e. count from one)"
       << " of the outgroup, if present in data\n";
}
