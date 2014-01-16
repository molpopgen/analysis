#ifndef __COMPUTE__FUNC_H__
#define __COMPUTE__FUNC_H__
#include <iosfwd>
#include <glob.h>

class compute_params;
void process(glob_t *files,  compute_params *args, std::ostream &ofstr,
	      int argc, char ** argv);
void makeheader(const int argc,  char *argv[], const compute_params *args, std::ostream &o);

#endif
