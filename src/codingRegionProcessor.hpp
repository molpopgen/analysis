#ifndef __CODING_REGION_PROCESSOR_HPP__
#define __CODING_REGION_PROCESSOR_HPP__

#include <vector>
#include <stdexcept>
#include <Sequence/PolySites.hpp>
#include <Sequence/Fasta.hpp>
#include <Warnings.hpp>

class cRPimpl;
class codingRegionProcessor
{
 private:
  cRPimpl *impl;
 public:
  explicit codingRegionProcessor(const std::vector< Sequence::Fasta > &data,
				 const std::vector<double> & snp_positions,
				 const std::vector< int > & intervals,
				 const bool & haveOutgroup = false,
				 const unsigned &outgroup = 0,
				 const bool & Approximate = false)
    throw (std::exception);
  ~codingRegionProcessor(void);
  Sequence::PolySites synonymousTable(void) const;
  Sequence::PolySites replacementTable(void) const;
  Warnings & warnings(void) const;
};

#endif
