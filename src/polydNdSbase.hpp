#include <vector>
#include <string>
#include <stdexcept>
#include <Sequence/PolySites.hpp>
#include <Sequence/Fasta.hpp>
#include <Warnings.hpp>
// #if defined(__GNUG__) && __GNUC__ < 3
// #include <climits>
// #else
#ifndef __POLYDNDSBASE_HPP__
#define __POLYDNDSBASE_HPP__

using std::string;
using std::vector;
using std::exception;
using Sequence::Fasta;
using Sequence::PolyTable;
using Sequence::PolySites;

// #if defined(__GNUG__) && __GNUC__ < 3
//   const unsigned UINTLOCALMAX = UINT_MAX;
// #else
//   const unsigned UINTLOCALMAX = std::numeric_limits<unsigned>::max();
// #endif

struct params
{
    char *infile;//file containing data
    vector<int> intervals;//intervals of coding sequence.
    bool haveOutgroup;//whether a list of ancestral states is available
    bool useTsTv;//estimate Ts/Tv ratio
    bool print_tables;//print out poly tables
    bool skipMissing;//skip sites with missing data
  bool SuperStrict;//require that annotation length be a multiple of 3bp
    bool parseableOutput;//make output less verbose
    bool XtractFullCodons;//for extractCoding, extract only complete codons
    bool ApproximateTreatmentOfDivergentCodons;//if true, treat codons diverged at 2 sites (with "ambiguous" pathways) and codons divergent at all 3 positions as silent or replacement based simply on translating the codon.
    int codon_start;//position in which first codon begins
    int codon_end;//position in which last codon ends
    unsigned int outgroup;//index of outgroup in data vector
    unsigned n1;//sample size of first population in data vector
  };

typedef void (*USG)(void);
void parseargs(int argc, char *argv[],params *args,USG usage)
  throw (std::exception);
bool InCoding(int site_index, const vector<int> &intervals);
int GetCodonPos(int site_index, const vector<int> &intervals);
void MakeCodon(int site_index,int seq, int codon_pos,
               const vector<Fasta >&data,
               const vector<int> &intervals,string &codon);
#endif
