#include <CodonMaker.hpp>
#include <polydNdSbase.hpp>
#include <algorithm>
#include <Sequence/Comparisons.hpp>
#include <Sequence/SeqProperties.hpp>

using Sequence::NumDiffs;

CodonMaker::CodonMaker(const vector<Fasta> &data,const vector<int> &intervals,
		       const int &indexed_pos,const int &codonPos):
  maxdiffs(0),success(false)
{
  success = true;
  std::string tc1,tc2;
  for (unsigned i = 0 ; i < data.size()-1 ; ++i)
    {
      MakeCodon(indexed_pos,i,codonPos,data,intervals,tc1);
      for (unsigned j = i+1 ; j < data.size() ; ++j)
	{
	  MakeCodon(indexed_pos,j,codonPos,data,intervals,tc2);
	  unsigned ndiffs = NumDiffs(tc1,tc2);
	  maxdiffs = (ndiffs > maxdiffs) ? ndiffs : maxdiffs;
	}
    }
  //make sure the codons we return differ by maxdiffs!
  codon1 = tc1;
  codon2 = tc2;
  if (unsigned(NumDiffs(codon1,codon2)) != maxdiffs || 
      (containsAmbiguousBases(codon1)==true && containsAmbiguousBases(codon1)==true) )
    {
      success = false;
      for (unsigned i = 0 ; i < data.size()-1 ; ++i)
	{
	  MakeCodon(indexed_pos,i,codonPos,data,intervals,tc1);
	  for (unsigned j = i+1 ; j < data.size() ; ++j)
	    {
	      MakeCodon(indexed_pos,j,codonPos,data,intervals,tc2);
	      if( unsigned(NumDiffs(tc1,tc2)) == maxdiffs)
		{
		  if (containsAmbiguousBases(tc1) == false &&
		      containsAmbiguousBases(tc2) == false && tc1 != tc2)
		    {
		      codon1 = tc1;
		      codon2 = tc2;
		      success = true;
		      i = j = data.size();
		    }
		}
	    }
	}
    }
}

bool CodonMaker::containsAmbiguousBases(const std::string & codon)
/*! 
  returns true if there is a non-standared DNA character,
  false otherwise
*/
{
  return (std::find_if(codon.begin(),codon.end(),Sequence::ambiguousNucleotide()) != codon.end());
}

 bool CodonMaker::Success(void)
{
  return success;
}

 unsigned CodonMaker::ndiffs(void)
{
  return maxdiffs;
} 

 std::string CodonMaker::c2(void)
{
  return codon2;
}

 std::string CodonMaker::c1(void)
{
  return codon1;
}
    
