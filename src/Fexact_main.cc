#include <cstdio>
#include <cstdlib>
#include <ctest.h>

using namespace std;

int main( int argc, char **argv)
{
  int a=2,y=2;
  int work=100000;
  int leading = 2;
  double expect,percnt,emin,prt,pre;
  double table[4];
  expect=percnt=emin=prt=pre=-1.;
  if (argc != 5)
    {
      fprintf(stderr,"usage: Fexact i j k l\n");
      exit(1);
    }
  table[0]=atof(argv[1]);
  table[1]=atof(argv[2]);
  table[2]=atof(argv[3]);
  table[3]=atof(argv[4]);
  fexact(&a,&y,table,&leading,&expect,&percnt,&emin,&prt,&pre,&work);
  fprintf(stdout, "%e\t%e\n", pre, prt);
  exit(0);
}
