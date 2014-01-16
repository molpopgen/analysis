#ifndef __WARNINGS_HPP__
#define __WARNINGS_HPP__

#include <vector>
#include <string>

class Warnings {
 private:
  std::vector< std::string > _warnings;
 public:
  void add(const std::string &s);
  typedef std::vector< std::string >::const_iterator const_iterator;
  const_iterator begin();
  const_iterator end();
};

#endif
