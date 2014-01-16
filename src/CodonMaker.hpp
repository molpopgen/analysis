#ifndef __CODON_MAKER_HPP__
#define __CODON_MAKER_HPP__

#include <string>
#include <vector>
#include <Sequence/Fasta.hpp>

using std::vector;
using Sequence::Fasta;

class CodonMaker
{
 private:
  std::string codon1,codon2;
  unsigned maxdiffs;
  bool success;
 public:
  explicit CodonMaker(const vector<Fasta> &data,const vector<int> &intervals,
		      const int &indexed_pos,const int &codonPos);
  std::string c1(void);
  std::string c2(void);
  unsigned ndiffs(void); 
  bool containsAmbiguousBases(const std::string & codon);
  bool Success(void);
  };

#endif
