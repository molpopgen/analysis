#include <vector>
#include <iostream>
#include <fstream>
#include <utility>
#include <numeric>
#include <limits>
#include <algorithm>
#include <iterator>
#include <functional>
#include <set>
#if defined (__SVR4) && defined (__sun)
#include <ieeefp.h>
#endif

#include <tuple>

//libsequence stuff
#if defined(__GNUG__) && __GNUC__ >= 3
#include <Sequence/FastaExplicit.hpp>
#else
#include <Sequence/Fasta.hpp>
#include <Sequence/Alignment.hpp>
#endif
#include <Sequence/PolySites.hpp>
#include <Sequence/stateCounter.hpp>
#include <Sequence/shortestPath.hpp>

//base routines for dealing with coding sequence data
#include <polydNdSbase.hpp>
#include <codingRegionProcessor.hpp>
#include <MKtestOutput.hpp>

using namespace std;
using namespace Sequence;
using namespace Alignment;

enum {FIXEDA=0,FIXEDS,POLYA,POLYS};

const double MAXDBL = numeric_limits<double>::max();


struct isPresent : public std::binary_function < 
  iterator_traits<Sequence::PolySites::const_site_iterator>::value_type,
  double,
  bool >
{
  typedef iterator_traits<Sequence::PolyTable::const_site_iterator>::value_type vt;
  inline bool operator()(const vt & v, const double &d) const
  {
    return (v.first == d);
  }
};

vector<string> makePopCodons(const int codon_pos,
			     const int site_index,
			     const vector<Fasta> & alignment,
			     const vector<int> & intervals);

set<string> makeCodSet(vector<string>::const_iterator beg,
		       vector<string>::const_iterator end);

set<shortestPath::pathType> catalogueDiffs(const set<string> &cod1,
					   const set<string> &cod2,
					   const int & cPos);

set<shortestPath::pathType> catalogueDiffsComplex(const set<string> &cod1,
						  const set<string> &cod2,
						  const int & cPos);
unsigned countFixations(const PolySites & table,
			const size_t & n1);

size_t IncrementSNPitr(PolySites::const_site_iterator beg,
		       PolySites::const_site_iterator end,
		       const unsigned &n,
		       const vector<int> &intervals);

void updateFixCells( const shortestPath::pathType & type,
		     unsigned * fixedA,
		     unsigned * fixedS);

tuple<double,double,double> getPositions(PolySites::const_site_iterator beg,
						PolySites::const_site_iterator end,
						const unsigned &n,
						const vector<int> &intervals);

void usage(void);

void (*usg)(void) = &usage;

int main (int argc, char *argv[])
{
  params args;
  try
    {
      parseargs(argc,argv,&args,usg);

      bool Approx = args.ApproximateTreatmentOfDivergentCodons;
      //check that args are valid;
      if(args.infile==NULL
	 || args.n1 == 0)
	{
	  cerr << "usage error.\n";
	  usage();
	  exit(1);
	}


      vector<Fasta> data;
      ifstream in(args.infile);
      if (! in)
	{
	  cerr << "Could not open data file " 
	       << args.infile
	       << "\n";
	  exit(10);
	}
      GetData(data,args.infile);
      if(!IsAlignment(data))
        {
          cerr << "error: data don't seem to be aligned...\n";
          exit(1);
        }

      if (args.intervals.size() == 0)
	{
	  args.intervals.resize(2);
	  args.intervals[0]=0;
	  args.intervals[1]=data[0].length()-1;
	}


      //get Poly data for the 2 genes separately
      vector<Fasta> gene1(data.begin(),data.begin()+args.n1);
      vector<Fasta> gene2(data.begin()+args.n1,data.end());
      vector<stateCounter> states1(data[0].length());
      vector<stateCounter> states2(data[0].length());
      vector<double> pos1,pos2;

      for(unsigned site = 0 ; site < data[0].length() ; ++site)
	{
	  for(unsigned seq=0;seq<data.size();++seq)
	    {
	      if (seq < args.n1)
		{
		  states1[site](data[seq][site]);
		}
	      else
		{
		  states2[site](data[seq][site]);
		}
	    }
	}

      for(unsigned i = 0 ; i < states1.size() ; ++i)
	{
	  if (states1[i].gap == 0 && states2[i].gap ==0)
	    {
	      if (states1[i].nStates() == 2)
		{
		  pos1.push_back(i+1);
		}
	      if (states2[i].nStates() == 2)
		{
		  pos2.push_back(i+1);
		}
	    }
	}

      Warnings w;

      codingRegionProcessor c1(gene1,pos1,args.intervals,
			       args.haveOutgroup,args.outgroup,
			       Approx);
      w = c1.warnings();

      //print the warnings to stderr
      cerr << "Warnings from file "<<args.infile<<":\n";
      std::copy(w.begin(),w.end(),
		std::ostream_iterator<const string>(std::cerr,"\n"));

      codingRegionProcessor c2(gene2,pos2,args.intervals,
			       args.haveOutgroup,args.outgroup,
			       Approx);
      w = c2.warnings();
      std::copy(w.begin(),w.end(),
		std::ostream_iterator<const string>(std::cerr,"\n"));

      PolySites *A1 = new PolySites(c1.replacementTable());
      PolySites *S1 = new PolySites(c1.synonymousTable());
      PolySites *A2 = new PolySites(c2.replacementTable());
      PolySites *S2 = new PolySites(c2.synonymousTable());
      //the following is necessary to handle sites polymorphic in 
      //both genes/species
      set<double> uniqueA,uniqueS;
      uniqueA.insert(A1->pbegin(),A1->pend());
      uniqueA.insert(A2->pbegin(),A2->pend());//now know # AA poly in alignment
      uniqueS.insert(S1->pbegin(),S1->pend());
      uniqueS.insert(S2->pbegin(),S2->pend());

      set<double> complex;
      //Process divergence
      PolySites *allSNP = new PolySites(data);
      PolySites::const_site_iterator sbeg = allSNP->sbegin(),
	send = allSNP->send();

      vector<unsigned> contingency_table(4,0);   //store cell counts
      contingency_table[POLYA] = uniqueA.size(); //repl Poly
      contingency_table[POLYS] = uniqueS.size(); //syn poly
      
      unsigned tot_fixations = 0;
      while (sbeg < send)
	{
	  int indexed_pos = int(sbeg->first)-1;
	  unsigned nPolyPerCodon = 0;
	  if (InCoding(indexed_pos,args.intervals)==true)
	    {
	      // need to assign a temp 
	      // b/c erase-remove bit below
	      // doesn't work w/std::set
	      string temp(sbeg->second.begin(),
			  sbeg->second.begin()+args.n1);
	      //remove missing data
	      temp.erase( remove(temp.begin(),temp.end(),'N'),temp.end() );
	      set<char> one(temp.begin(),temp.end());

	      //same for 2nd aligment partition
	      temp = string(sbeg->second.begin()+args.n1,
			    sbeg->second.end());
	      temp.erase( remove(temp.begin(),temp.end(),'N'),temp.end() );
	      set<char> two(temp.begin(),temp.end());

	      //calculate overlap b/w partitions
	      vector<char> overlap(one.size()+two.size());
	      
	      vector<char>::iterator itr = set_intersection(one.begin(),one.end(),
							    two.begin(),two.end(),
							    overlap.begin());


	      if ( itr-overlap.begin() == 0 && one.size()==1 && two.size()==1) //fixed diff
		{
	      int cPos = GetCodonPos(indexed_pos,args.intervals);//codon pos of mut
		  //The set overlap is of size 0, therefore
		  //no states are shared b/w the two partitions.
		  //This is a fixed difference.
		  vector<string> codons = makePopCodons(cPos,indexed_pos,
							data,args.intervals);
		  PolySites codPol(codons);
		  nPolyPerCodon = codPol.numsites();
		  //number of fixations at this codon
		  unsigned nfixations = countFixations(codPol,args.n1); 
		  tot_fixations += nfixations;
		  set<string> cod1 = makeCodSet(codons.begin(),
						codons.begin()+args.n1);
		  set<string> cod2 = makeCodSet(codons.begin()+args.n1,
						codons.end());

		  if (nfixations == 1)
		    {
		      //this is the easiest case
		      set<shortestPath::pathType> types = catalogueDiffs(cod1,cod2,cPos);
		      if(types.size()==1)
			{
			  switch ( *(types.begin()) )
			    {
			    case shortestPath::pathType::S :
			      contingency_table[FIXEDS]++;
			      break;
			    case shortestPath::pathType::N :
			      contingency_table[FIXEDA]++;
			      break;
			    default:
			      complex.insert(sbeg->first);
			      break;
			    }
			}
		      else
			{
			  complex.insert(sbeg->first);
			}
		    }
		  else if (nfixations > 1)
		    {
		      //harder case
		      set<shortestPath::pathType> types = catalogueDiffsComplex(cod1,cod2,cPos);
		      if (types.size()==1)
			{
			  updateFixCells( *(types.begin()),
					  &contingency_table[FIXEDA],
					  &contingency_table[FIXEDS]);
			}
		      else
			{
			  //fixation occurs at a codon w/complex history,
			  //so we save the positions in the set "complex"
			  tuple<double,double,double> compPos = 
			    getPositions(sbeg,send,nPolyPerCodon,args.intervals);
			  double x = get<0>(compPos);
			  if ( x != MAXDBL)
			    {
			      complex.insert(x);
			    }
			  x = get<1>(compPos);
			  if ( x != MAXDBL)
			    {
			      complex.insert(x);
			    }
			  x = get<2>(compPos);
			  if ( x != MAXDBL)
			    {
			      complex.insert(x);
			    }
			}
		    }
		}
	      else
		{
		  //states are shared, therefore no fixation
		}
	    }
 	  sbeg += IncrementSNPitr(sbeg,send,nPolyPerCodon,
				  args.intervals);
	}
      MKtestOutput output(contingency_table,args);
      cout << output<<endl;
      if (args.parseableOutput==false)
	{
	  cout << "Total number of fixations:\t"<<tot_fixations<<endl;
	  if (!complex.empty()) 
	    {
	      cout << "The following positions contain substitutions with\n"
		   << "complex histories:\n";
	      copy(complex.begin(),complex.end(),ostream_iterator<double>(cout," "));
	      cout << endl;
	    }
	}
      cout << endl;
      exit(0);
    }
  catch (Sequence::SeqException &e)
    {
      cerr << "Processing file "
	   << args.infile
	   << "resulting in the following exception begin thrown:\n";
      cerr << e << endl;
      exit(10);
    }
  catch (std::exception &e)
    {
      cerr << "Processing file "
	   << args.infile
	   << "resulting in the following exception begin thrown:\n";
      cerr << e.what() << endl;
      exit(10);
    }
}

vector<string> makePopCodons(const int codon_pos,
			     const int site_index,
			     const vector<Fasta> & alignment,
			     const vector<int> & intervals)
/*!
  Reduces a data set to just the codon of interest.  Returns
  the codon in a vector of strings
*/
{
  vector<string> rv;

  for(unsigned i=0;i<alignment.size();++i)
    {
      string codon;
      MakeCodon(site_index,i,codon_pos,alignment,intervals,codon);
      rv.push_back(codon);
    }
  return rv;
}

set<string> makeCodSet(vector<string>::const_iterator beg,
		       vector<string>::const_iterator end)
/*!
  Makes a std::set of the strings in the range beg,end.
  Useful if you just want the unique codons in a range
*/
{
  set<string>rv;
  while(beg<end)
    {
      //skip codons w/missing data
      if (beg->find('N')==string::npos)
	rv.insert(*beg);
      ++beg;
    }
  return rv;
}

set<shortestPath::pathType> catalogueDiffs(const set<string> &cod1,
					   const set<string> &cod2,
					   const int & cPos)
/*!
  Calculates the set of shortestPath::pathType differences at 
  codon position cPos between all codons in cod1 and all in cod2.
*/
{
  typedef tuple<shortestPath::pathType,
    shortestPath::pathType,
    shortestPath::pathType> codonTuple;
  set<shortestPath::pathType> types;
  set<string>::const_iterator b1,b2;
  for(b1=cod1.begin() ; b1 != cod1.end() ; ++b1)
    {
      for(b2=cod2.begin() ; b2 != cod2.end() ; ++b2)
	{
	  codonTuple t = diffTypeMulti(*b1,*b2);
	  shortestPath::pathType type;
	  
	  //retrieve type at the site of substitution
	  if (cPos == 0)
	    {
	      type = get<0>(t);
	    }
	  else if (cPos == 1)
	    {
	      type = get<1>(t);
	    }
	  else if (cPos == 2)
	    {
	      type = get<2>(t);
	    }
	  types.insert(type);
	}
    }
  return types;
}

set<shortestPath::pathType> catalogueDiffsComplex(const set<string> &cod1,
						  const set<string> &cod2,
						  const int & cPos)
{
  set<shortestPath::pathType> types;
  set<string>::const_iterator b1,b2;
  for(b1=cod1.begin() ; b1 != cod1.end() ; ++b1)
    {
      for(b2=cod2.begin() ; b2 != cod2.end() ; ++b2)
	{
	  shortestPath sp(*b1,*b2);
	  shortestPath::pathType type = sp.type();
 	  types.insert(type);
	}
    }
  return types;
}

unsigned countFixations(const PolySites & table,const size_t &n1)
/*!
  Counts all the fixations in a snp table
*/
{		  
  PolySites::const_site_iterator beg = table.sbegin(),
    end=table.send();
  unsigned nfix=0;
  while(beg<end)
    {
      set<char> _one(beg->second.begin(),
		     beg->second.begin()+n1);
      set<char> _two(beg->second.begin()+n1,
		     beg->second.end());
      vector<char> overlap(_one.size()+_two.size());
      vector<char>::iterator itr = set_intersection(_one.begin(),_one.end(),
						    _two.begin(),_two.end(),
						    overlap.begin());
      if (itr - overlap.begin() == 0)
	++nfix;
      ++beg;
    }
  return nfix;
}

void updateFixCells( const shortestPath::pathType & type,
		     unsigned * fixedA,
		     unsigned * fixedS )
/*!
  Updates cell entries in MK tables.
*/
{
  switch (type)
    {
    case shortestPath::pathType::S :
      ++(*fixedS);
      break;
    case shortestPath::pathType::N :
      ++(*fixedA);
      break;
    case shortestPath::pathType::SN :
      (*fixedS)++;
      (*fixedA)++;
      break;
    case shortestPath::pathType::SS :
      (*fixedS)+=2;
      break;
    case shortestPath::pathType::NN :
      (*fixedA)+=2;
      
      break;
    case shortestPath::pathType::SSS :
            (*fixedS)+=3;
      
      break;
    case shortestPath::pathType::SSN :
            (*fixedS)+=2;
      
            (*fixedA)++;
      break;
    case shortestPath::pathType::SNN :
      (*fixedS)++;
            (*fixedA)+=2;
      break;
    case shortestPath::pathType::NNN :
            (*fixedA)+=3;
	    //(*fixedA)+=1;
      break;
    default:
      break;
    }
}

size_t IncrementSNPitr(PolySites::const_site_iterator beg,
		       PolySites::const_site_iterator end,
		       const unsigned &n,
		       const vector<int> &intervals)
/*!
  Figures out how much to increment the pointer to site columns.  Used in
  processing fixed differences.
*/
{
  size_t k=1;
  while(beg<end && k<n)
    {
      int index = int(beg->first)-1;
      if (InCoding(index,intervals))
	{
	  ++k;
	}
      ++beg;
    }
  return k;
}

tuple<double,double,double> getPositions(PolySites::const_site_iterator beg,
						PolySites::const_site_iterator end,
						const unsigned &n,
						const vector<int> &intervals)
{
  double rv[3];
  rv[0]=rv[1]=rv[2]=MAXDBL;
  unsigned k = 0;
  while(beg<end && k<n)
    {
      int index = int(beg->first)-1;
      if (InCoding(index,intervals))
	{
	  rv[k] = beg->first;
	  ++k;
	}
      ++beg;
    }
  return make_tuple(rv[0],rv[1],rv[2]);
}

void usage(void)
{
  cerr << "MKtest -i <infile> -I nint i j k l ..."<<endl;
  cerr << "see man MKtest for details"<<endl;
}
