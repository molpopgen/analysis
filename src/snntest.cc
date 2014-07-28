/*! 
  snn.cc - perform the Snn test of Hudson (2000) Genetics 155: 2011-2014
*/
#include <Sequence/SimpleSNP.hpp>
#include <Sequence/SimData.hpp>
#include <Sequence/PolySites.hpp>
#include <Sequence/Snn.hpp>
#include <Sequence/PolyTableFunctions.hpp>
#if __GNUG__ && __GNUC__ >= 3
#include <Sequence/FastaExplicit.hpp>
#else
#include <Sequence/Fasta.hpp>
#include <Sequence/Alignment.hpp>
#endif
#include <getopt.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <vector>
#include <numeric>
#include <iostream>
#include <fstream>
#include <ctime>
#include <cstdlib>

enum FILETYPE { UNKNOWN,FASTA,SIMPLE,MS };
struct params
{
  FILETYPE ft;
  bool biallelic;
  std::string infile,outfile;
  std::vector<unsigned> config;
  unsigned nperms;
  int outgroup;
  params() : ft(UNKNOWN),biallelic(false),
	     infile(std::string()),outfile(std::string()),
	     config(std::vector<unsigned>()),nperms(10000),
	     outgroup(-1)
  {
  }
};

void usage();

params parseargs(int argc, char **argv);

int main(int argc, char **argv)
{
  params p = parseargs(argc,argv);
  const unsigned totsam = std::accumulate(p.config.begin(),p.config.end(),0u);

  unsigned npop = p.config.size();
  //offsets are the size_t for pointer arithmetic
  unsigned * offsets = new unsigned[npop+1];
  offsets[0] = 0;
  for(unsigned i=1;i<=npop;++i)
    {
      offsets[i] = offsets[i-1]+p.config[i-1];
    }

  Sequence::PolySites snpTable;
  bool isms = false;
  try
    {
      if( p.ft == MS )
	{
	  isms = true;
	}
      else
	{
	  std::ifstream in(p.infile.c_str());
	  
	  if(in)
	    {
	      if(p.ft==FASTA)
		{
		  std::vector<Sequence::Fasta> data;
		  if(p.outgroup >= 0)
		    {
		      data.erase(data.begin()+p.outgroup);
		    }
		  Sequence::Alignment::GetData(data,in);
		  snpTable = Sequence::PolySites(data);//,(p.outgroup >= 0),p.outgroup);
		}
	      else if (p.ft == SIMPLE)
		{
		  Sequence::SimpleSNP temp;
		  in >> temp;
		  if(temp.outgroup() == true)
		    {
		      //the analysis requires that
		      //no outgroup sequence be in the data,
		      //so get rid of it.
		      //For class Sequence::SimpleSNP, 
		      //the outgroup is always the 1st
		      //sequence in the data
		      snpTable.assign( &*(temp.pbegin()),temp.numsites(),
				       &*(temp.begin()+1),temp.size()-1);
		      Sequence::RemoveInvariantColumns(&snpTable);
		    }
		  else
		    {
		      snpTable.assign( &*(temp.pbegin()),temp.numsites(),
				       &*(temp.begin()),temp.size());
		    }
		}
	    }
	  else
	    {
	      std::cerr << "fatal error: "
			<< p.infile
			<< " could neither be found nor opened\n";
	      std::exit(1);
	    }
	}
    }
  catch (Sequence::SeqException &e)
    {
      std::cerr << e << '\n';
      std::exit(1);
    }
  catch (std::exception &e)
    {
      std::cerr << e.what() << '\n';
      std::exit(1);
    }
  catch (...)
    {
      std::cerr << "unexpected exception caught!\n";
      std::exit(1);
    }
  if(p.ft != MS && snpTable.size() != totsam)
    {
      std::cerr << "error: sample size in infile does not equal sum of sizes from locales\n";
      std::exit(1);
    }

  if(! Sequence::Alignment::validForPolyAnalysis(snpTable.begin(),snpTable.end()) )
    {
      std::cerr << "error: snp table resulting from data contains characters which are not valid\n";
    }

  if(p.biallelic)
    {
      snpTable.RemoveMultiHits();
    }

  gsl_rng * r = gsl_rng_alloc(gsl_rng_default);
  gsl_rng_set(r,0);

  std::ofstream outfile;
  if( ! p.outfile.empty() )
    {
      outfile.open(p.outfile.c_str());
      std::cout.rdbuf(outfile.rdbuf());
    }
  std::pair<double,double> result;
  std::vector<std::vector<double> > pairwise_result;
  if(isms == false)
    {
      std::cout << '#'<<p.infile<<'\n'
		<<"#n = "<< snpTable.size()
		<< '\t'
		<<"S = "<< snpTable.numsites()
		<< '\n';

      std::cout << "pop1\tpop2\tsnn\tp\n";
      result = Sequence::Snn_test(snpTable,&(p.config[0]),npop,
				  std::bind(gsl_ran_flat,r,0.,std::placeholders::_1),
				  p.nperms);
      std::cout << "nan\tnan\t" << result.first << ' ' << result.second << '\n';
      if(npop>2)
	{
	  pairwise_result = Sequence::Snn_test_pairwise(snpTable,
							&(p.config[0]),npop,
							std::bind(gsl_ran_flat,r,0.,std::placeholders::_1),
							p.nperms);
	  for( unsigned i = 0 ; i < pairwise_result.size() ; ++i )
	    {
	      std::cout << pairwise_result[i][0] << ' '
			<< pairwise_result[i][1] << ' '
			<< pairwise_result[i][2] << ' '
			<< pairwise_result[i][3] << '\n';
	    }
	}
    }
  else
    {
      int rv,rep=0;
      Sequence::SimData d;
      std::cout << "rep\tpop1\tpop2\tsnn\tp\n";
      while( (rv=d.fromfile(stdin)) != EOF )
	{
	  result = Sequence::Snn_test(d,&(p.config[0]),npop,
	    std::bind(gsl_ran_flat,r,0.,std::placeholders::_1),
				      p.nperms);
	  std::cout << rep << "\tnan\tnan\t" << result.first << '\t' << result.second << '\n';
	  if(npop>2)
	    {
	      pairwise_result = Sequence::Snn_test_pairwise(d,&(p.config[0]),npop,
							    std::bind(gsl_ran_flat,r,0.,std::placeholders::_1),
							    p.nperms);
	      
	      for( unsigned i = 0 ; i < pairwise_result.size() ; ++i )
		{
		  std::cout << rep << '\t'
			    << pairwise_result[i][0] << '\t'
			    << pairwise_result[i][1] << '\t'
			    << pairwise_result[i][2] << '\t'
			    << pairwise_result[i][3] << '\n';
		}
	    }
	  rep++;
	}
    }
  if(! p.outfile.empty() )
    {
      outfile.close();
    }
}


params parseargs(int argc, char **argv)
{
  params p;
  int c,npop,dummy;
  if(argc == 1)
    {
      usage();
      exit(1);
    }
  //LEOPARD PROBLEM?
  while ((c = getopt (argc, argv, "f:s:n:c:o:O:bm")) != -1)
    {
      switch (c)
        {
	case 'f':
	  p.ft = FASTA;
	  p.infile = std::string(optarg);
	  break;
	case 's':
	  p.ft = SIMPLE;
	  p.infile = std::string(optarg);
	  break;
	case 'n':
	  p.nperms = atoi(optarg);
	  break;
	case 'c':
	  npop = atoi(optarg);
	  dummy=0;
	  for (int k = optind ; k < optind+npop ; ++k,++dummy)
            {
              p.config.push_back(atoi(argv[k]));
            }
	  optind+=dummy;
          break;
	case 'o':
	  p.outfile = std::string(optarg);
	  break;
	case 'O':
	  p.outgroup = atoi(optarg-1);
	  break;
	case 'b':
	  p.biallelic = true;
	  break;
	case 'm':
	  p.ft = MS;
	  break;
	default:
	  usage();
	  exit(1);
	  break;
	}
    }
  return p;
}

void usage()
{
  std::cerr << "snn [options], where options can be:\n"
	    << "\t-f filename : data are in fasta format in a file called filename\n"
	    << "\t\tfurther options for fasta data:\n"
	    << "\t\t-O outgroup : if there is an outgroup in the file, pass it the number (1 <= outgroup <= n)\n"
	    << "\t-s filename : data are in the format used for Hudson's (2001) programs\n"
	    << "\t-m : read ms-like data from stdin\n"
	    << "(note, only one of -f, -s, or -m are valid!)\n"
	    << "\t-c npop n_1 n_2 ... n_npop : the number of pops followed by list of sample sizes per pop\n"
	    << "\t-o outfilename : write data to outfilename (defaults to writing to the screen)\n"
	    << "\t-n nperms : get p-values from nperms permutations of the data (default=10000)\n"
	    << "\t-b : only analyze bi-allelic sites (only considers the ingroup!)\n";
}
