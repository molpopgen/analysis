#ifndef __CORRELATIONS_H__
#define __CORRELATIONS_H__

#include <vector>
#include <cassert>
#include <cfloat>
#include <cmath>
#include <algorithm>
#include <iterator>
#include "RandomNumbers.hpp"

/*! \file Correlations.hpp
  \ingroup Numerology
  Generic algorithms for calculating and testing the significance
  of correlation coefficients
*/
using std::pow;
using std::fabs;

namespace Numerology
{
  struct ProductMoment
  /*!
    \ingroup Numerology
    Calculates Pearson's product-moment correlation between two containers, \c x and \c y.
    \f[r_{1,2}=\frac{y_1\times y_2}{\sqrt{\sum_{i=1}^{i=n} y_1^2 \times \sum \sum_{i=1}^{i=n} y_2^2}}\f]
    From Sokal and Rohlf. Biometry, 3rd ed. 1995, Freeman Press, page  571.
    \param beg_x pointer to beginning of \c x
    \param end_x pointer to one past the end of \c x
    \param beg_y pointer to beginning of \c y
    \param end_y pointer to one past the end of \c y
    \note the square of the return value is the well-know coefficient of determination, r^2
    \warning if the ::value_type of the pointers is not a floating-point or rational type,
    this function returns rubbish
    \n \n Example: \n
    \code 
    #include "Numerology.hpp"
    #include <string>
    #include <iostream>
    using namespace std;
    using namespace Numerology;
    int main()
    {
    vector<double> x(12);
    vector<double> y(12);
      
    x[0]=159;
    x[1]=179;
    x[2]=100;
    x[3]=45;
    x[4]=384;
    x[5]=230;
    x[6]=100;
    x[7]=320;
    x[8]=80;
    x[9]=220;
    x[10]=320;
    x[11]=210;
    y[0]=14.4;
    y[1]=15.2;
    y[2]=11.3;
    y[3]=2.5;
    y[4]=22.7;
    y[5]=14.9;
    y[6]=1.41;
    y[7]=15.81;
    y[8]=4.19;
    y[9]=15.39;
    y[10]=17.25;
    y[11]=9.52;
    cout << ProductMoment()( x.begin(),x.end(),y.begin(),y.end() )<<endl;
    }
    \endcode
  */
  {
    template <class ForwardIterator > 
    inline typename std::iterator_traits<ForwardIterator>::value_type
    operator()( ForwardIterator beg_x, 
		ForwardIterator end_x,
		ForwardIterator beg_y,
		ForwardIterator end_y)
    {
      assert ( std::distance(beg_x,end_x)==std::distance(beg_y,end_y) );
      typedef typename std::iterator_traits<ForwardIterator>::value_type rtype;
      unsigned nsam = 0;
      rtype x_bar=rtype();
      rtype y_bar=rtype();
      rtype x_squared=rtype();
      rtype y_squared=rtype();
      rtype xy=rtype();
      
      for ( ; beg_x != end_x && beg_y != end_y  ; ++beg_x,++beg_y)
	{
	  x_bar += *beg_x;
	  x_squared += (*beg_x)*(*beg_x);
	  y_bar += *beg_y;
	  y_squared += (*beg_y)*(*beg_y);	
	  xy += (*beg_x)*(*beg_y);
	  ++nsam;
	}
      x_bar /= rtype(nsam);
      y_bar /= rtype(nsam);
      x_squared /= rtype(nsam);
      y_squared /= rtype(nsam);
      xy -= x_bar*rtype(nsam)*y_bar;
      xy -= y_bar*rtype(nsam)*x_bar;
      xy += rtype(nsam)*x_bar*y_bar;
      rtype sum_x_sq = rtype(nsam)*(x_squared-pow(x_bar,2.0));
      rtype sum_y_sq = rtype(nsam)*(y_squared-pow(y_bar,2.0));
      return(xy/pow((sum_x_sq*sum_y_sq),0.5));
    }
  };

  template<class T> struct PermuteProductMomentNoScramble
  /*! 
    \ingroup Numerology
    This function allows you to test the significance of a Pearson's Product-Moment correlation
    by a permutation test.
    \param beg_x pointer to beginning of \c x
    \param end_x pointer to one past the end of \c x
    \param beg_y pointer to beginning of \c y
    \param end_y pointer to one past the end of \c y
    \param NPERM the number of permutations to do
    \return the two-tailed probability of observing a p-m correlation
    \note This is the logically \c const version of 
    Numerology::PermuteProductMoment<class BidirectionalIterator>.  It does
    not scramble any of the containers pointed to
    \note Implementation involves calls to Numerology::ProductMoment.
    \note This function suffices for a permutation test of the coefficient of
    determination, since its just a transformation of the product-moment correlation
    \warning if the ::value_type of the pointers is not a floating-point or rational type,
    this function returns rubbish
    \n \n Example: \n
    \code 
    #include "Numerology.hpp"
    #include <string>
    #include <iostream>
    using namespace std;
    using namespace Numerology;
    int main()
    {
    vector<double> x(12);
    vector<double> y(12);
      
    x[0]=159;
    x[1]=179;
    x[2]=100;
    x[3]=45;
    x[4]=384;
    x[5]=230;
    x[6]=100;
    x[7]=320;
    x[8]=80;
    x[9]=220;
    x[10]=320;
    x[11]=210;
    y[0]=14.4;
    y[1]=15.2;
    y[2]=11.3;
    y[3]=2.5;
    y[4]=22.7;
    y[5]=14.9;
    y[6]=1.41;
    y[7]=15.81;
    y[8]=4.19;
    y[9]=15.39;
    y[10]=17.25;
    y[11]=9.52;
    cout << ProductMoment()( x.begin(),x.end(),y.begin(),y.end() )<<endl;
    cout << PermuteProductMomentNoScramble()( &x[0],&x[x.size()]
    ,&y[0],&y[y.size()])<<endl;
    }
    \endcode
  */
  {
  public:
    template<class BidirectionalIterator>
    inline typename std::iterator_traits<BidirectionalIterator>::value_type
    operator()(BidirectionalIterator beg_x, 
	       BidirectionalIterator end_x,
	       BidirectionalIterator beg_y,
	       BidirectionalIterator end_y,
	       unsigned NPERM=10000)
    {
      assert ( std::distance(beg_x,end_x)==std::distance(beg_y,end_y) );
      typedef typename std::iterator_traits<BidirectionalIterator>::value_type rtype;
      double _pm = ProductMoment()(beg_x,end_x,beg_y,end_y);
      unsigned _prob=0;
      std::vector<T> copy_x;
      std::vector<T> copy_y;
	
      copy_x.assign(beg_x,end_x);
      copy_y.assign(beg_y,end_y);
      Numerology::UniformDeviate rand;
      for(unsigned i = 0 ; i < NPERM ; ++i)
	{
	  std::random_shuffle(copy_x.begin(),copy_x.end(),rand);
	  if (fabs(ProductMoment()(copy_x.begin(),copy_x.end()
				   ,copy_y.begin(),copy_y.end())) 
	      - fabs(_pm) >= DBL_EPSILON)
	    ++_prob;
	}
      return rtype(_prob)/rtype(NPERM);
    }
  };

  struct PermuteProductMoment
  /*!
    \ingroup Numerology
    This function allows you to test the significance of a Pearson's Product-Moment correlation
    by a permutation test.
    \param beg_x pointer to beginning of \c x
    \param end_x pointer to one past the end of \c x
    \param beg_y pointer to beginning of \c y
    \param end_y pointer to one past the end of \c y
    \param NPERM the number of permutations to do
    \return the two-tailed probability of observing a p-m correlation
    \note implementation involves calls to Numerology::ProductMoment.
    \note This function suffices for a permutation test of the coefficient of
    determination, since its just a transformation of the product-moment correlation
    \warning this function scrambles the container pointed to by \a beg_x and \a end_x. Also, 
    if the ::value_type of the pointers is not a floating-point or rational type,
    this function returns rubbish
    \n \n Example: \n
    \code 
    #include "Numerology.hpp"
    #include <string>
    #include <iostream>
    using namespace std;
    using namespace Numerology;
    int main()
    {
    vector<double> x(12);
    vector<double> y(12);
      
    x[0]=159;
    x[1]=179;
    x[2]=100;
    x[3]=45;
    x[4]=384;
    x[5]=230;
    x[6]=100;
    x[7]=320;
    x[8]=80;
    x[9]=220;
    x[10]=320;
    x[11]=210;
    y[0]=14.4;
    y[1]=15.2;
    y[2]=11.3;
    y[3]=2.5;
    y[4]=22.7;
    y[5]=14.9;
    y[6]=1.41;
    y[7]=15.81;
    y[8]=4.19;
    y[9]=15.39;
    y[10]=17.25;
    y[11]=9.52;
    cout << ProductMoment()( x.begin(),x.end(),y.begin(),y.end() )<<endl;
    cout << PermuteProductMoment()(x.begin(),x.end(),y.begin(),y.end())<<endl;
    }
    \endcode
  */
  {
  public:
    template<class BidirectionalIterator>
    inline typename std::iterator_traits<BidirectionalIterator>::value_type
    operator()(BidirectionalIterator beg_x, 
	       BidirectionalIterator end_x,
	       BidirectionalIterator beg_y,
	       BidirectionalIterator end_y,
	       unsigned NPERM=10000)
    {
      assert ( std::distance(beg_x,end_x)==std::distance(beg_y,end_y) );
      typedef typename std::iterator_traits<BidirectionalIterator>::value_type rtype;
      double _pm = ProductMoment()(beg_x,end_x,beg_y,end_y);
      unsigned _prob=0;
      Numerology::UniformDeviate rand;
      for(unsigned i = 0 ; i < NPERM ; ++i)
	{
	  std::random_shuffle(beg_x,end_x,rand);
	  if (fabs(ProductMoment()(beg_x,end_x,beg_y,end_y)) 
	      - fabs(_pm) >= DBL_EPSILON)
	    ++_prob;
	}
      return rtype(_prob)/rtype(NPERM);
    }
  };

 
}
#endif
