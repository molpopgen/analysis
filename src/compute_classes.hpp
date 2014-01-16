#ifndef __COMPUTE_CLASSES_HPP__
#define __COMPUTE_CLASSES_HPP__

#include <iosfwd>
#include <gsl/gsl_rng.h>

struct compute_params
{
  char *infileglob;
  char *outfile;
  bool haveOutgroup;
  unsigned int outgroup;
  bool useTotMuts;
  bool is_table;
  bool suppress_headers;
  bool bi_allelic_only;
  bool no_missing;
  bool print_and_die;
  bool probs;
  bool verbose;
  bool use_theta;
  bool pretty;
  explicit compute_params();
};

class resultsImpl;
class results
/*!
  store summary statistics for a data set
*/
{
  friend class pvals;  //yes, I'm cheating...
private:
  resultsImpl * impl;
  const char * infile;
public:
  explicit results(const char * _infile,
		    compute_params * args);
  bool wasSkipped() const;
  std::ostream & print(std::ostream & o) const;
  ~results();
};

class pvalsImpl;
class pvals
/*!
  store pvalues for a dataset, calculated by coalescent
  simulations with no recombination
*/
{
private:
  pvalsImpl * impl;
public:
  explicit pvals(gsl_rng * r,const results & res,
		 const compute_params * args);
  std::ostream & print(std::ostream & o) const;
  ~pvals();
};

inline
std::ostream & operator<<(std::ostream &o,
			  const results & r)
{
  return r.print(o);
}

inline
std::ostream & operator<<(std::ostream &o,
			  const pvals & p)
{
  return p.print(o);
}

#endif
