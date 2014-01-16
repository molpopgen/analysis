#include <Warnings.hpp>

void Warnings::add(const std::string &s)
{
  _warnings.push_back(s);
}

Warnings::const_iterator Warnings::begin()
{
  return _warnings.begin();
}

Warnings::const_iterator Warnings::end()
{
  return _warnings.end();
}
