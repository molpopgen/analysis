#include <polyCommon.hpp>

#include<algorithm>
#include <config.h>
#include <Sequence/SeqConstants.hpp>
#include <Sequence/PolySites.hpp>
#if defined(HAVE_SEQUENCE_ROTATEPOLYTABLE) && !defined(HAVE_SEQUENCE_MAKE_POLYSITEVECTOR)
#include <Sequence/PolyTableManip.hpp>
#else
#include <Sequence/polySiteVector.hpp>
#endif 
#include <Sequence/PolySites.hpp>



using namespace std;
using namespace Sequence;

//global constants
const char _alphabet[6] = {'A','G','C','T','0','1'};

const unsigned maxUINT = Sequence::SEQMAXUNSIGNED;

vector<stateCounter> processStateInfo(const vector<Fasta> &data, const unsigned &datalen,
				      const bool &haveOutgroup, const unsigned &outgroup)
{
  vector<stateCounter> stateCounts(datalen,stateCounter('-'));
  unsigned nsam = data.size();
  for(unsigned i = 0 ; i < datalen ; ++i)
    {
      for(unsigned j = 0 ; j < nsam ; ++j)
	{
	  if ( (! haveOutgroup)
	       || (haveOutgroup && j != outgroup 
		   && data[j][i] != data[outgroup][i] && data[outgroup][i]!='-') )
	    {
	      stateCounts[i](data[j][i]);
	    }
	}
    }
  return stateCounts;
}

vpuu getFreqSpectrum(const vector<stateCounter> &counts,unsigned sample_size,bool unfolded)
{

  vpuu freqSpec;
  vpuu::iterator freqSpeqItr;

  unsigned mincount = maxUINT;
  if(!unfolded)
    {
      for(unsigned i = 0 ; i < counts.size() ; ++i)
	{
	  mincount = maxUINT;
	  if(counts[i].nStates() == 2 && counts[i].gap == 0)
	    {
	      mincount = (counts[i].a > 0 && counts[i].a < mincount) ? counts[i].a : mincount;
	      mincount = (counts[i].g > 0 && counts[i].g < mincount) ? counts[i].g : mincount;
	      mincount = (counts[i].c > 0 && counts[i].c < mincount) ? counts[i].c : mincount;
	      mincount = (counts[i].t > 0 && counts[i].t < mincount) ? counts[i].t : mincount;
	      unsigned tofind = min(mincount,sample_size-mincount);
	      if ( (freqSpeqItr = find_if(freqSpec.begin(),freqSpec.end(),
					  bind2nd(findFreq(),tofind))) != freqSpec.end())
		{
		  (*freqSpeqItr).second++;
		}
	      else
		{
		  freqSpec.push_back( puu(tofind,1) );
		}
	    }
	}
    }
  else
    {
      for(unsigned i=0;i<counts.size();++i)
	{
	  mincount = maxUINT;
	  if(counts[i].nStates()==1 && counts[i].gap==0)
	    {
	      mincount = (counts[i].a > 0) ? counts[i].a : maxUINT;
	      mincount = (counts[i].g > 0 && mincount==maxUINT) ? counts[i].g : mincount;
	      mincount = (counts[i].c > 0 && mincount==maxUINT) ? counts[i].c : mincount;
	      mincount = (counts[i].t > 0 && mincount==maxUINT) ? counts[i].t : mincount;
	      unsigned tofind = mincount;
	      if(tofind < sample_size)
	      {
		if ( (freqSpeqItr = find_if(freqSpec.begin(),freqSpec.end(),
					    bind2nd(findFreq(),tofind))) != freqSpec.end())
		  {
		    (*freqSpeqItr).second++;
		  }
		else
		  {
		    freqSpec.push_back( puu(tofind,1) );
		  }
		}
	    }
	}
    }
  sort(freqSpec.begin(),freqSpec.end(),sortFreq());
  return freqSpec;
}


vpuu checkAdjacentPoly(const vector<Fasta> &data, const bool &haveOutgroup, 
		       const unsigned &outgroup)
{
  PolySites poly(data);
  vpuu adjacents;
  if(poly.numsites() > 1)
    {
      //make a 90-degree counterclockwise rotation of the poly table
#if defined(HAVE_SEQUENCE_ROTATEPOLYTABLE) && !defined(HAVE_SEQUENCE_MAKE_POLYSITEVECTOR)
      polySiteVector _Data = rotatePolyTable(&poly);
#else
      polySiteVector _Data = make_polySiteVector(poly);
#endif

      //if the outgroup is present, we need
      //to remove it from the rotated data
      if (haveOutgroup == true)
	{
	  for(unsigned site = 0 ; site < _Data.size() ; ++site)
	    {
	      _Data[site].second.erase(outgroup,1);
	    }
	}
      std::pair<char,char> chars1,chars2;
      for(unsigned site = 1 ; site < poly.numsites() ; ++site)
	{
	  if( poly.position(site)-1. == poly.position(site-1) )
	    {
	      stateCounter countsI,countsJ;
	      for (unsigned seq = 0 ; seq < poly.size() ; ++seq)
		{
		  if( haveOutgroup == false
		      || (haveOutgroup == true && seq !=  outgroup) )
		    {
		      countsI(poly[seq][site]);
		      countsJ(poly[seq][site-1]);
		    }
		}
	      if(countsI.nStates()==2 && countsJ.nStates()==2)
		{
		  chars1.first = chars2.first = 'Z';//Z is a dummy value
		  chars1.second = chars2.second = 'Z';
		  
		  //find out which 2 states segregate in sites site and site-1
		  for(unsigned i = 0 ; i < 6 ; ++i)
		    {
		      if (std::find(_Data[site].second.begin(),
				    _Data[site].second.end(),
				    _alphabet[i]) != _Data[site].second.end())
			{
			  if (chars1.first == 'Z')
			    chars1.first = _alphabet[i];
			  else
			    chars1.second = _alphabet[i];
			}
		      if (std::find(_Data[site-1].second.begin(),
				    _Data[site-1].second.end(),
				    _alphabet[i]) != _Data[site-1].second.end())
			{
			  if (chars2.first == 'Z')
			    chars2.first = _alphabet[i];
			  else
			    chars2.second = _alphabet[i];
			}
		    }
		  //now, which is the minor state and
		  //ensure that .first is the minor character
		  if(std::count(_Data[site].second.begin(),
				_Data[site].second.end(),chars1.first)
		     > std::count(_Data[site].second.begin(),
				  _Data[site].second.end(),chars1.second))
		    {
		      swap(chars1.first,chars1.second);
		    }
		  if(std::count(_Data[site-1].second.begin(),
				_Data[site-1].second.end(),chars2.first)
		     > std::count(_Data[site-1].second.begin(),
				  _Data[site-1].second.end(),chars2.second))
		    {
		      swap(chars2.first,chars2.second);
		    }
		  for (unsigned seq = 0 ; seq < poly.size() ; ++seq)
		    {
		      if( haveOutgroup == false
			  || (haveOutgroup == true && seq !=  outgroup) )
			{
			  if (poly[seq][site]==chars1.first && poly[seq][site-1]==chars2.first)
			    {
			      adjacents.push_back(puu(unsigned(poly.position(site-1)),
								   unsigned(poly.position(site))));
			      break;
			    }
			}
		    }
		}
	    }
	}
    }
  return adjacents;
}
