#include<cassert>
#include<vector>
#include<string>

using namespace std;
//template class TsTvBias {
class TsTvBias
  {
  private:
    void Compute(const vector<string> &data);
    unsigned nTs,nTv;
    unsigned int datalen;
  public:
    TsTvBias(const vector<string> &data,unsigned int length);
    double pTs(void);
    double pTv(void);
  };
