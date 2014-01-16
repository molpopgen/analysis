#include <iostream>
#include <memory>

#include <getopt.h>

#if defined(__GNUG__) && __GNUC__ >= 3
#include <Sequence/FastaExplicit.hpp>
#else
#include <Sequence/Alignment.hpp>
#include <Sequence/Fasta.hpp>
#endif

#include <Sequence/Kimura80.hpp>
#include <Sequence/SeqExceptions.hpp>

using namespace std;
using namespace Sequence;
using namespace Alignment;

int main(int argc,char *argv[])
{
  char *infile = NULL;
  char *outfile = NULL;
  int c;


  while ((c = getopt (argc, argv, "i:")) != -1)
    {
      switch (c)
        {
        case 'i':
          infile = optarg;
          break;
        case 'o':
          outfile = optarg;
          break;
        }
    }

  vector<Fasta > data;
  try
    {
      if (infile != NULL)
        GetData(data,infile);
      else
        {
          cerr << "error: no infile specified"<<endl;
          exit(10);
        }
    }
  catch (SeqException &e)
    {
      e.print(cerr);
      exit(10);
    }
  //  vector<double> distances;
  for(unsigned int i = 0 ; i < data.size()-1 ; ++i)
    for(unsigned int j = i+1 ; j < data.size() ; ++j)
      {
        try
          {
            auto_ptr<Kimura80> K(new Kimura80(&data[i],&data[j]));
            cout << data[i].GetName() <<'\t'<<data[j].GetName()<<'\t'<<K->K()<<endl;
          }
        catch (SeqException &e)
          {
            cerr << "error: sequences "<<data[i].GetName()<<" and " <<data[j].GetName()<<" are not the same length, skipping.\n";
          }
      }
  return(0);
}
