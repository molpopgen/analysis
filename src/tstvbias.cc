#include <tstvbias.hpp>
#include <Sequence/stateCounter.hpp>

TsTvBias::TsTvBias(const vector<string> &data,unsigned int length)
  :nTs(0),nTv(0),datalen((!data.empty())?data[0].length():0)
{
  TsTvBias::Compute(data);
}

void TsTvBias::Compute(const vector<string> &data)
{
  for(unsigned int i= 0 ; i < datalen ; ++i)
    {
      Sequence::stateCounter Counts;
      for(unsigned int j = 0 ; j <data.size() ; ++j)
        {
          Counts(data[j][i]);
        }
      unsigned nstates = Counts.nStates();
      if (nstates == 2 && Counts.gap == 0)//sum over bi-allelic, ungapped sites
        {
          nTs += (Counts.a > 0 && Counts.g > 0) ? 1 : 0;
          nTs += (Counts.c > 0 && Counts.t > 0) ? 1 : 0;

          nTv += (Counts.a > 0 && (Counts.c > 0 || Counts.t > 0)) ? 1 : 0;
          nTv += (Counts.g > 0 && (Counts.c > 0 || Counts.t > 0)) ? 1 : 0;
          nTv += (Counts.c > 0 && (Counts.a > 0 || Counts.g > 0)) ? 1 : 0;
          nTv += (Counts.t > 0 && (Counts.a > 0 || Counts.g > 0)) ? 1 : 0;
        }
    }
}

double TsTvBias::pTs(void)
{
  return double(nTs)/double(nTs+nTv);
}

double TsTvBias::pTv(void)
{
  return double(nTv)/double(nTs+nTv);
}
