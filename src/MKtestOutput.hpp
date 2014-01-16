#include <vector>
#include <polydNdSbase.hpp>
#include <iosfwd>

class MKtestOutput
{
 private:
  const vector<unsigned> * _contingency_table;
  const params * _args;
 public:
  MKtestOutput( const vector<unsigned> & contingency_table,
		const params & args );
  std::ostream & print(std::ostream & s);
};

inline std::ostream &
operator<< (std::ostream & s, MKtestOutput & c)
{
  return c.print(s);
}

