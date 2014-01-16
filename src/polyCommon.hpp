#ifndef __POLY_COMMON_HPP__
#define __POLY_COMMON_HPP__
#include <vector>
#include <utility>

//libsequence
#include <Sequence/stateCounter.hpp>
#if defined(__GNUG__) && __GNUC__ >= 3
#include <Sequence/FastaExplicit.hpp>
#else
#include <Sequence/Fasta.hpp>
#include <Sequence/Alignment.hpp>
#endif

using std::vector;
using std::pair;
using Sequence::stateCounter;
using Sequence::Fasta;

//typedefs
typedef std::pair<unsigned,unsigned> puu;
typedef std::vector< puu > vpuu;

//prototypes
vector<stateCounter> processStateInfo(const vector<Fasta> &data, const unsigned &datalen,
				      const bool &haveOutgroup, const unsigned &outgroup);
vpuu getFreqSpectrum(const vector<stateCounter> &counts, unsigned sample_size,bool unfolded);
vpuu checkAdjacentPoly(const vector<Fasta> &data, const bool &haveOutgroup, 
		       const unsigned &outgroup);

//predictates
struct findFreq : public std::binary_function< puu,unsigned,bool >
{
  inline bool operator()( const puu &x, const unsigned &y ) const
  {
    return x.first == y;
  }
};

struct sortFreq : public std::binary_function< puu,puu,bool >
{
  inline bool operator()( const puu &x, const puu &y ) const
  {
    return x.first < y.first;
  }
};

struct hasMissing : public std::unary_function< Sequence::stateCounter, bool >
{
  inline bool operator()( const Sequence::stateCounter &counts ) const
  {
    return (counts.n > 0);
  }
};

struct multiHit : public std::unary_function< Sequence::stateCounter, bool >
{
  inline bool operator()( const Sequence::stateCounter &counts ) const
  {
    return (counts.nStates() > 2);
  }
};

struct gapOverSNP : public std::unary_function< Sequence::stateCounter, bool >
{
  inline bool operator()( const Sequence::stateCounter &counts ) const
  {
    return (counts.gap > 0 && counts.nStates() > 1);
  }
};

#endif
