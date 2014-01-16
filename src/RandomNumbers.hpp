#ifndef __RANDOM_NUMBERS_H__
#define __RANDOM_NUMBERS_H__
/*! \defgroup RandomNum
  These are mostly for convenience.  The templates defined in group
  RandomNumBoost should be preferred.
  @quick and dirty random number routines
 */

/*! \file RandomNumbers.hpp
  \ingroup RandomNum
  Functors related to random numbers.
*/

#include <cmath>
#include <ctime>
#include <functional>  

namespace Numerology
{
  class UniformDeviate:public std::unary_function<int,int>
  /*!
    \ingroup RandomNum
    A functor to return a random number uniformly
    distributed between 0 and N, where N = \c RAND_MAX
    by default.  The random number generation is done
    by a call to the C function \c rand().  Can be
    used as the functor for \c std::random_shuffle()
    @short base of the Random Numbers classes
  */
  {
  private:
    static bool seeded;

  public:
    inline UniformDeviate(bool seedWithTime = true);
    
      int operator()(int N=RAND_MAX) const
      /*!
	This can be used by \c std::random_shuffle
	as a generator.  
	\return \c ((double(rand())/RAND_MAX))*N , 
 	following Stroustrup, 3rd ed., p. 685,
 	which guarantees a return value
 	>=0 and < N-1.
      */
    {
      return int(((double(rand())/RAND_MAX)) * N);
    }
  };

  bool UniformDeviate::seeded = false;

   UniformDeviate::UniformDeviate(bool seedWithTime)
    /*!
      The constructor checks if the RNG has been
      seeded.  If not, it seeds it with the system 
      time
      \param seedWithTime if true, RNG seeds with system time,
      otherwise default seed is used
    */
  {
    if (seeded == false && seedWithTime == true) 
      {
	//seed with system time
	time_t timeval;  
	static_cast<void>(std::time(&timeval));
	srand(long(timeval));
	seeded = true ;
      }
  }

  struct Probability : private UniformDeviate, 
		       public std::unary_function<void,double>
  /*!
    \ingroup RandomNum
    @short Functor to generate a double on the interval [0,1) that is uniformly distributed
  */
  {
  public:
    explicit Probability (bool seedWithTime = true):UniformDeviate(seedWithTime){}
      double operator()(void) const
    {
      return UniformDeviate::operator()(RAND_MAX)/(double(RAND_MAX)+1.);
    }
  };

  struct Exponential : private Probability,
		       public std::unary_function<double,double>
  /*!
    \ingroup RandomNum
    @short Functor to return exponentially-distributed deviates
  */
  {
  public:
    explicit Exponential(bool seedWithTime=true):Probability(seedWithTime)
      /*!
	calls constructor of Numerology::UniformDeviate
	to ensure that RNG is seeded if seedWithTime==true
      */
    {
    }
      double operator()(double mean = 1.0) const
      /*!
	\param mean the mean of an exponential distribution
	\return a value drawn from an exponential distribution with
	mean \a mean
      */
    {
      return -1.0*mean*std::log(1.0 - Probability::operator()());
    }
  };

  struct Gaussian: private Probability, 
		   public std::binary_function<double,double,double>
  /*!
    \ingroup RandomNum
    @short function object to return Gaussian Deviates
  */
  {
  public:
    explicit Gaussian(bool seedWithTime = true):Probability(seedWithTime){}
      double operator()(double mean = 0.0,double variance = 1.0)
      //stolen from Dick Hudson's C code for his "ms" simulation
    {
      static int iset=0;
      static float gset;
      float fac,r,v1,v2;
      
      if  (iset == 0) 
	{
	  do 
	    {
	      v1=2.0*Probability::operator()() - 1.0;
	      v2=2.0*Probability::operator()() - 1.0;
	      r=v1*v1+v2*v2;
	    } 
	  while (r >= 1.0);
	  fac=sqrt(-2.0*log(r)/r);
	  gset= v1*fac;
	  iset=1;
	  return( mean + sqrt(variance)*v2*fac);
	} 
      else 
	{
	  iset=0;
	  return( mean + sqrt(variance)*gset ) ;
	}
    }
  };

  struct Poisson : private Probability,
		   public std::unary_function<double,unsigned>
  /*!
    \ingroup RandomNum
    Functor to return deviates drawn from a Poisson distributions
  */
  {
  public:
    explicit Poisson(bool seedWithTime=true):Probability(seedWithTime){}
    unsigned operator()(double mean = 1., double maxMean = 30.0)
      /*!
	\param mean the mean of the Poisson distribution
	\param maxMean if \a mean > \a maxMean, a gaussian deviate is returned instead
      */
    {
      //poisson's w/large means are well-approximated by Gaussians
      if(mean >= maxMean)
	return ( unsigned(0.5 + Gaussian()(mean,mean) ) );
      double uniform = Probability::operator()();
      double p = exp(-mean);
      if (uniform < p ) return(0);
      double cumProb = p;
      unsigned i = 1;
      
      while( uniform > ( cumProb += (p *= mean/i ) ) )
	{
	  ++i;
	}
      return i;
    }
  };
}

#endif
