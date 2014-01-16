#include <iostream>
#include <fstream>
#include <glob.h>
#include <getopt.h>

//for interupt handling
#include <int_handler.hpp>
#include <compute_functions.hpp>
#include <compute_classes.hpp>

/*!
  \example compute.cc
  This program analyzes polymorphism data stored in Fasta files.
  Output is in spreadsheet format.
*/

using namespace std;

void parseargs(int argc, char *argv[],compute_params *args);
void usage(void);

int main(int argc, char *argv[])
{
  glob_t files;
  ofstream ofstr;
  compute_params args;

  parseargs(argc,argv,&args);

  if(args.infileglob == NULL)
    {
      cout << "fatal error: no infile(s) specified\n";
      usage();
      exit(10);
    }

  //filenames are obtained via a pattern passed from
  //the shell, i.e. compute -i '*.fasta', where the
  //pattern in single quotes is what gets processed
  //see man glob(3)
#if defined (__SVR4) && defined (__sun)
  glob(args.infileglob,GLOB_ERR,NULL,&files);
#else
  glob(args.infileglob,GLOB_TILDE|GLOB_ERR,NULL,&files);
#endif

  if (args.outfile != NULL)
    ofstr.open(args.outfile);

//   if (! args.suppress_headers)
//     {
//       if(args.outfile != NULL)
//         makeheader(argc,argv,&args,ofstr);
//       else
// 	{
// 	  if (args.pretty == false)
// 	    makeheader(argc,argv,&args,cout);
// 	}
//     }

  signal(SIGINT,cntrl_c_handler);
  try
    {
      if (args.outfile != NULL)
	process(&files,&args,ofstr,argc,argv);
      else
	process(&files,&args,cout,argc,argv);
    }
  catch (std::exception &e)
  {
    std::cerr << e.what() << std::endl;
    exit(10);
  }

  globfree(&files);
  ofstr.flush();
  ofstr.close();
  exit(0);
}

void parseargs(int argc, char *argv[],compute_params *args)
{
  if(argc==1)
    {
      usage();
      exit(1);
    }
  extern int optind;
  int c;

  while ((c = getopt (argc, argv, "i:h:o:O:nNsbPpvtV")) != -1)
    {
      switch (c)
        {
        case 'i':
          args->infileglob = optarg;
          break;
        case 'h':
          args->infileglob = optarg;
          args->is_table = true;
          break;
        case 'o':
          args->outfile = optarg;
          break;
        case 'n':
          args->useTotMuts=0;
          break;
	case 'N':
	  args->no_missing = true;
	  break;
        case 'O':
          args->haveOutgroup=1;
          args->outgroup=atoi(optarg)-1;
          break;
        case 's':
          args->suppress_headers=1;
          break;
        case 'b':
          args->bi_allelic_only = 1;
          break;
        case 'P':
          args->print_and_die = 1;
          args->suppress_headers=1;
          break;
        case 'p':
          args->probs = 1;
          break;
        case 'v':
          args->verbose=1;
          break;
	case 'V':
	  args->pretty=true;
	  break;
        case 't':
          args->use_theta=1;
          break;
        case '?':
          cerr<<"error: unknown argument "<< argv[optind-1]<<endl;
          usage();
          exit(0);
          break;
        }
    }
}

void usage(void)
{
  cerr<<"usage: ";
  cerr<<"compute -i infile"<<endl;
  cerr<<"\tother options\n";
  cerr<<"\t\t-o outfile : write results to outfile\n";
  cerr<<"\t\t-O n : use the nth sequence in the data as an outgroup\n";
  cerr<<"\t\t(start counting from one!!!)\n";
  cerr<<"\t\t-n : use the total # of segregating sites, rather than # of mutations\n";
  cerr<<"\t\t-v : verbose progress reporting to standard error\n";
  cerr<<"\t\t-V : pretty output\n";
  cerr<<"see manpage for more details"<<endl;
  cerr << endl;
}
