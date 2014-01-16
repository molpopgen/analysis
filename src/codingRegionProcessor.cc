#include <codingRegionProcessor.hpp>
#include <cRPimpl.hpp>

codingRegionProcessor::codingRegionProcessor(const std::vector< Sequence::Fasta > &data,
					     const std::vector<double> & snp_positions,
					     const std::vector< int > & intervals,
					     const bool & haveOutgroup,
					     const unsigned &outgroup,
					     const bool & Approximate)
  throw (std::exception)
/*!
  \param data the raw sequence data
  \param snp_positions positions of SNPs in data (start counting from 1!!)
  \param intervals a vector of integers demarcating the exon boundaries in \a data.  Please be aware
  that this class assumes that \a intervals contains the indexes of the sites of exon boundaries, 
  rather than the actual positions (i.e. count starting from 0, not 1).
  \param haveOutgroup true if there is outgroup sequence present in \a data, false otherwise
  \param outgroup the index of the outgroup in \a data if \a haveOutgroup is true
  \param Approximate if false (the default), then codons that differ at 2 positions
  will be skipped if the 2 alternate pathways between the 2 codons are not complementary
  (by complementary, I mean both pathways must have equal numbers of the same types of changes).
  If true, the most conservative pathway between such ambiguous codons is assumed to be
  the true pathway.
  \exception std::exception
*/
{
  try
    {
      impl = new cRPimpl(data,snp_positions,
			 intervals,
			 haveOutgroup,
			 outgroup,
			 Approximate);
    }
  catch (std::exception &e)
    {
      delete impl;
      throw( e );
    }
}

codingRegionProcessor::~codingRegionProcessor(void)
{
  delete impl;
}

Sequence::PolySites codingRegionProcessor::synonymousTable(void) const
/*!
  \return an object of type Sequence::PolySites representing the synonymous
  changes in the exons
*/
{
  return Sequence::PolySites(impl->silent_pos,impl->silent_char);
}

Sequence::PolySites codingRegionProcessor::replacementTable(void) const
/*!
  \return an object of type Sequence::PolySites representing the nonsynonymous
  changes in the exons
*/
{
  return Sequence::PolySites(impl->repl_pos,impl->repl_char);
}

Warnings & codingRegionProcessor::warnings(void) const
/*!
  \return an object of type warnings containing std::string's describing
  any potential warnings that arose during the calculation
*/
{
  return impl->w;
}
