#ifndef __CRPIMPL_HPP__
#define __CRPIMPL_HPP__

#include <vector>
#include <string>
#include <stdexcept>
#include <Sequence/Fasta.hpp>
//#include <Sequence/PolyTable.hpp>
#include <Warnings.hpp>

struct cRPimpl
{
  enum mutantType {UNDETERMINED,SILENT,REPLACEMENT};
  Warnings w;
  std::vector<double> silent_pos,repl_pos;
  std::vector<std::string> silent_char,repl_char;
  const std::vector< int > intervals;
  const bool haveOutgroup;
  const unsigned outgroup;
  const bool Approximate;
  void process_coding_region(const std::vector< Sequence::Fasta > &data,
			     const std::vector<double> & snp_positions) throw (std::exception);
  bool isStop(const std::string &codon);
  void addToTable(const std::vector<Sequence::Fasta> &data, const std::vector<double> & snp_positions,
		  const int &indexed_pos, std::vector<double> &pos, 
		  std::vector<std::string> &matrix);
  //  unsigned incrementSiteIndex(const unsigned &nDiffsCodon,const unsigned &cPos);
  explicit cRPimpl(const std::vector< Sequence::Fasta > &data,
		   const std::vector< double > & snp_positions,
		   const std::vector< int > & _intervals,
		   const bool & _haveOutgroup,
		   const unsigned &_outgroup,
		   const bool & _Approximate) 
    throw (std::exception);
};
#endif
