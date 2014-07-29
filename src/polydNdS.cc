#include <vector>
#include <iostream>
#include <memory>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <iterator>
#if defined (__SVR4) && defined (__sun)
#include <ieeefp.h>
#endif

#if defined( __GNUG__ ) && __GNUC__ >= 3
#include <Sequence/FastaExplicit.hpp>
#else
#include <Sequence/Fasta.hpp>
#include <Sequence/Alignment.hpp>
#endif

#include <Sequence/PolySites.hpp>
#include <Sequence/PolyTableFunctions.hpp>
#include <Sequence/PolySNP.hpp>
#include <Sequence/Comeron95.hpp>
#include <Sequence/RedundancyCom95.hpp>
//#include <Sequence/SeqRegexes.hpp>
#include <codingRegionProcessor.hpp>
#include <tstvbias.hpp>
#include <polydNdSbase.hpp>

using namespace std;
using namespace Sequence;
using namespace Alignment;

void output_statistics(const PolyTable *poly_table,params *args,double n_sites);
void simple_output(const PolyTable *poly_table,params *args,double n_sites);
void GetSiteCounts(const vector<Fasta > &data, const PolyTable * all_data,
		   const params *args,
		   double *n_repl_sites,double *n_silent_sites,
		   double * n_third_sites, double * n_4fold_sites,
		   bool useTsTv=false,double pTs=0.0,double pTv=0.0);
void PrintTable(const PolyTable *table, const string &filename,
		const std::vector<string> &seqnames, const string & outgroup_states,
		bool haveOutgroup=false,
		int outgroup = 0);
void MergeTables(std::vector<double> &newpos, std::vector<string> &newchar,
		 const std::vector<double> &pos1, const std::vector<string> &char1,
		 const std::vector<double> &pos2, const std::vector<string> &char2);
unsigned getNonCodingLen(const std::vector<Sequence::Fasta> &data, 
			 const std::vector<int> &intervals,const params *args);
void usage(void);

void (*usg)(void) = &usage;

int main (int argc, char *argv[])
{
  params args;
  try
    {
      parseargs(argc,argv,&args,usg);
    }
  catch (std::exception &e)
    {
      std::cerr << e.what() << std::endl;
      exit(10);
    }
  vector<int> interval_copy =args.intervals;
  if(args.codon_start > 1)
    {
      switch (args.codon_start)
	{
	case 2:
	  //if it starts at a 2nd position,
	  //move 2 positions to the right
	  //to a 1st position
	  args.intervals[0] += 2;
	  break;
	case 3:
	  //if it starts at a 3rd position,
	  //move 1 position to the right
	  //to a 1st position
	  args.intervals[0] += 1;
	  break;
	}
    }
  if (args.codon_end < 3)
    {
      switch (args.codon_end)
	{
	case 1:
	  //if it ends at a first position, move 1
	  //position left to a 3rd position
	  args.intervals[args.intervals.size()-1] -= 1;
	  break;
	case 2:
	  //if it ends at a second position, move 2
	  //positions left to a 3rd position
	  args.intervals[args.intervals.size()-1] -= 2;
	  break;
	}
    }
  vector<Fasta > data;
  try
    {
      ifstream input(args.infile);
      if (!input)
	{
	  cerr << "error : "
	       << args.infile
	       << " cannot be opened.\n";
	  exit(1);
	}
      GetData(data,input);
      if (! validForPolyAnalysis(data.begin(),data.end()) )
	{
	  throw badFormat("characters other than A,G,C,T,N were encountered in the alignment.\nThis program is only intended to deal with haploid (or phased diploid) data, and no ambiguity symbols.");
	}
      if(!IsAlignment(data))
	{
	  cerr << "error: data don't seem to be aligned...\n";
	  exit(1);
	}
    }
  catch(SeqException &e)
    {
      cerr << "Processing file "<<args.infile<<" resulted in an exception being thrown:"<<endl;
      e.print(cerr);
      cerr <<endl;
      exit(1);
    }
  //now, analyze the polymorphism at silent and replacement sites.
  //save outgroup states if we output poly tables
  //This is a total HACK!
  string outgroup_sequence = (args.haveOutgroup) ? data[args.outgroup] : string();
  bool hadOutgroup = args.haveOutgroup;
  if(args.haveOutgroup)
    {

      vector<Fasta> temp_data(data.begin(),data.begin()+args.outgroup);
      copy(data.begin()+args.outgroup+1,data.end(),back_inserter(temp_data));
      data=temp_data;
      args.haveOutgroup=false;
      args.outgroup=0;
    }

  PolySites * all_data = new PolySites(data,true,true,args.skipMissing);
  unsigned NSITES = all_data->numsites();
  all_data->RemoveAmbiguous();
  cerr << "Warnings from file "<<args.infile<<":\n";
  if(all_data->numsites() != NSITES)
    {
      cerr << (NSITES-all_data->numsites()) << " positions removed due to presence of ambiguous bases\n";
    }
  //skip all codons that have 3 alleles at 1 site
  vector<unsigned> to_remove;
  for(unsigned site=0;site<all_data->numsites();++site)
    {
      int index = int(all_data->position(site))-1;
      if (InCoding(index,args.intervals))
	{
	  int cPos = GetCodonPos(index,args.intervals);
	  vector<int> to_check;
	  switch(cPos)
	    {
	    case 0:
	      to_check.push_back(index+1);
	      to_check.push_back(index+2);
	      break;
	    case 1:
	      to_check.push_back(index-1);
	      to_check.push_back(index+1);
	      break;
	    case 2:
	      to_check.push_back(index-1);
	      to_check.push_back(index-2);
	      break;
	    }
	  for (unsigned ch=0;ch<to_check.size();++ch)
	    {
	      if(InCoding(to_check[ch],args.intervals))
		{
		  stateCounter c;
		  for(unsigned ind=0;ind<data.size();++ind)
		    {
		      c(data[ind][to_check[ch]]);
		    }
		  if(c.nStates()>2)
		    {
		      to_remove.push_back(site);
		      ch = to_check.size();
		    }
		}
	    }
	}
    }
  if (!to_remove.empty())
    {
      vector<double> new_pos;
      vector<string> new_data(all_data->size());
      for(unsigned site = 0;site<all_data->numsites();++site)
	{
	  if (find(to_remove.begin(),to_remove.end(),site) == to_remove.end())
	    {
	      new_pos.push_back(all_data->position(site));
	      for(unsigned ind=0;ind<all_data->size();++ind)
		{
		  new_data[ind]+=(*all_data)[ind][site];
		}
	    }
	}
      all_data->assign(&new_pos[0],new_pos.size(),
		       &new_data[0],new_data.size());
    }
  TsTvBias bias(all_data->GetData(),all_data->numsites());
  vector<double> silent_pos,repl_pos;
  vector<string> silent_char,repl_char;
  Warnings w;
  try
    {
      codingRegionProcessor c(data,all_data->GetPositions(),args.intervals,
			      args.haveOutgroup,args.outgroup,
			      args.ApproximateTreatmentOfDivergentCodons);
      w = c.warnings();
      silent_pos = c.synonymousTable().GetPositions();
      silent_char = c.synonymousTable().GetData();
      repl_pos = c.replacementTable().GetPositions();
      repl_char = c.replacementTable().GetData();
    }
  catch (std::logic_error &e)
    {
      cerr << e.what() << endl;
    }
  catch (std::runtime_error &e)
    {
      cerr << e.what() << endl;
      exit(1);
    }
  catch(...)
    {
      cerr << "whoah!  unexpected exception caught!"<<endl;
      exit(1);
    }

  //print the warnings to stderr
  std::copy(w.begin(),w.end(),
	    std::ostream_iterator<const string>(std::cerr,"\n"));

  //output results
  //get the number of silent and replacement sites in the data
  double n_repl_sites=0.,n_silent_sites=0.,n_third_sites=0.,n_4fold_sites=0.;
  GetSiteCounts(data,all_data,&args,
		&n_repl_sites,&n_silent_sites,&n_third_sites,&n_4fold_sites,
		args.useTsTv,bias.pTs(),bias.pTv());
  double non_coding_len = double(getNonCodingLen(data,interval_copy,&args));
  //1.) do the whole locus
  if (!args.parseableOutput)
    {
      cout << args.infile << '\n'
	   << "Mean # of replacement sites = "<<n_repl_sites<<'\n'
	   << "Mean # of synonymous sites = "<<n_silent_sites<<'\n'
	   << "The # of third positions = " << n_third_sites<<'\n'
	   << "Mean fourfold degenerate sites " << n_4fold_sites<<'\n'
	   << "The # of noncoding sites = " << non_coding_len << '\n'
	   << "Polymorphism in entire rgion:"<<'\n';
    }
  else
    {
      cout << args.infile << '\t';
    }
  double ug_len = double(UnGappedLength(data));
  output_statistics(all_data,&args,ug_len);
  if (!args.parseableOutput)
    cout << endl;

  //  2.) just the exons
  vector<double> exon_pos;
  vector<string> exon_char;
  MergeTables(exon_pos,exon_char,repl_pos,repl_char,silent_pos,silent_char);
  PolySites *exon_table = new PolySites(exon_pos,exon_char);
  if ( args.skipMissing )
    {
      exon_table->RemoveMissing();
    }
  if (!args.parseableOutput)
    cout << "Polymorphism in coding (i.e. translated) region:"<<endl;
  output_statistics(exon_table,&args,n_repl_sites+n_silent_sites);
  if (!args.parseableOutput)
    cout << endl;
  //3.)just the introns and flanking
  vector<double> intron_pos;
  vector<string> intron_char(data.size());
  for(unsigned i = 0 ; i < all_data->numsites() ; ++i)
    {
      //subtract 1 from position to make it the index of the site along the alignment
      if (! InCoding(unsigned(all_data->position(i))-1,interval_copy) )
	{
	  intron_pos.push_back(all_data->position(i));
	  for (unsigned j = 0 ; j < intron_char.size() ; ++j)
	    {
	      intron_char[j] += (*all_data)[j][i];
	    }
	}
    }
  PolySites *intron_table = new PolySites(intron_pos,intron_char);
  if ( args.skipMissing )
    {
      intron_table->RemoveMissing();
    }

  if (!args.parseableOutput)
    cout <<"Noncoding/untranslated polymorphism:"<<endl;
  output_statistics(intron_table,&args,non_coding_len);
  if (!args.parseableOutput)
    cout << endl;
  //4.) just the replacement polymorphisms
  PolySites *repl_table = new PolySites(repl_pos,repl_char);
  if ( args.skipMissing )
    {
      repl_table->RemoveMissing();
    }
  if (!args.parseableOutput)
    cout << "Replacement polymorphisms in coding region:"<<endl;
  output_statistics(repl_table,&args,n_repl_sites);
  if (!args.parseableOutput)
    cout << endl;
  //5.) just the silent polymorphisms
  PolySites *silent_table = new PolySites(silent_pos,silent_char);
  if ( args.skipMissing )
    {
      silent_table->RemoveMissing();
    }
  if (!args.parseableOutput)
    cout << "Silent polymorphisms in coding region:"<<endl;
  output_statistics(silent_table,&args,n_silent_sites);
  if (!args.parseableOutput)
    cout << endl;
  //5a.) just silent polymorphisms at 3rd positions
  vector<double> third_positions;
  vector<string> third_position_characters(silent_char.size());
  for(unsigned p = 0 ; p < silent_pos.size() ; ++p)
    {
      int index = int(silent_pos[p])-1;
      int codon_position = GetCodonPos(index,args.intervals);
      //if it's a 3rd position, add it...
      if(codon_position == 2)
	{
	  third_positions.push_back(silent_pos[p]);
	  for(unsigned seq = 0 ; seq < silent_char.size() ; ++seq)
	    {
	      third_position_characters[seq] += silent_char[seq][p];
	    }
	}
    }
  PolySites * third_positions_table = new PolySites(third_positions,third_position_characters);
  if ( args.skipMissing )
    {
      third_positions_table->RemoveMissing();
    }
  if(!args.parseableOutput)
    {
      cout << "Polymorphisms at third positions\n";
    }
  output_statistics(third_positions_table,&args,n_third_sites);
  if(!args.parseableOutput)
    cout << endl;
  //5b. just silent polymorphisms in 4fold degenerate codons

  vector<double> fourfold_positions;
  vector<string> fourfold_characters(third_position_characters.size());
  RedundancyCom95 redundancy_checker;
  string codon;
  //note--4fold degerate families can only vary at the third position
  for( unsigned p=0 ; p<third_positions.size() ; ++p )
    {
      int index = int(third_positions[p])-1;
      //this loop will return true iff all codons in alignment
      //are 4fold-degenerate
      //bool is4fold=false;
      bool is4fold=true;
      for(unsigned seq=0;is4fold && seq<data.size();++seq)
	{
	  MakeCodon(index,seq,2,data,args.intervals,codon);
	  if( count(codon.begin(),codon.end(),'N') == 0 &&
	      count(codon.begin(),codon.end(),'n') == 0 &&
	      count(codon.begin(),codon.end(),'-') == 0 )
	    {
	      try
		{
		  // if (redundancy_checker.ThirdFour(codon) == 1.)
		  //   {
		  //     is4fold=true;
		  //   }
		  // else
		  //   {
		  //     is4fold=false;
		  //    break;
		  //  }
		  if (redundancy_checker.ThirdFour(codon) != 1.)
		    {
		      is4fold=false;
		    }
		}
	      catch(SeqException & e)
		{
		  cerr << "This codon: " << codon << " raised this exception:\n"
		       << e << '\n';
		  exit(10);
		}
	    }
	}
      if (is4fold)
	{
	  fourfold_positions.push_back(third_positions[p]);
	  for(unsigned seq=0 ; seq < third_position_characters.size() ; ++seq)
	    {
	      //fourfold_characters[seq] += silent_char[seq][p];
	      fourfold_characters[seq] += third_position_characters[seq][p];
	    }
	}
    }
  PolySites * four_fold_sites_table = new PolySites(fourfold_positions,
					fourfold_characters);
  if ( args.skipMissing )
    {
      four_fold_sites_table->RemoveMissing();
    }

  if(!args.parseableOutput)
    cout << "All fourfold-degenerate positions:\n";
  output_statistics(four_fold_sites_table,&args,n_4fold_sites);
  if(!args.parseableOutput)
    cout << endl;

  //6.) silent + intron/flanking polymorphism
  vector<double> silent_intron_flanking;
  vector<string> silent_intron_flanking_chars;
  //merge the positions of the silent and intron/flaking SNPs
  MergeTables(silent_intron_flanking,silent_intron_flanking_chars,
	      silent_pos,silent_char,intron_pos,intron_char);
  //make a new polymorphism table
  PolySites *all_silent_table = new PolySites(silent_intron_flanking,
					      silent_intron_flanking_chars);
  if ( args.skipMissing )
    {
      all_silent_table->RemoveMissing();
    }
  if (!args.parseableOutput)
    cout << "All silent polymorphisms (coding + noncoding):"<<endl;
  output_statistics(all_silent_table,&args,n_silent_sites+non_coding_len);
  if (!args.parseableOutput)
    cout << endl;
  if (args.print_tables == true)
    {
      //get the names of the sequences
      vector<string> seqnames;
      for(unsigned i = 0 ; i < data.size() ; ++i)
	{
	  seqnames.push_back(data[i].GetName());
	}
      string outfile = args.infile;
      outfile += ".replacement";
      PrintTable(repl_table,outfile,seqnames,outgroup_sequence,hadOutgroup);

      outfile = args.infile;
      outfile += ".synonymous";
      PrintTable(silent_table,outfile,seqnames,outgroup_sequence,hadOutgroup);
		
      outfile = args.infile;
      outfile += ".3rdpositions";
      PrintTable(third_positions_table,outfile,seqnames,outgroup_sequence,hadOutgroup);

      outfile = args.infile;
      outfile += ".4fold";
      PrintTable(four_fold_sites_table,outfile,seqnames,outgroup_sequence,hadOutgroup);

      outfile = args.infile;
      outfile += ".exons";
      PrintTable(exon_table,outfile,seqnames,outgroup_sequence,hadOutgroup);

      outfile = args.infile;
      outfile += ".introns_flanking";
      PrintTable(intron_table,outfile,seqnames,outgroup_sequence,hadOutgroup);

      outfile = args.infile;
      outfile += ".all_silent";
      PrintTable(all_silent_table,outfile,seqnames,outgroup_sequence,hadOutgroup);
    }
  if (args.parseableOutput == true)
    cout << endl;
  delete all_data;
  delete exon_table;
  delete intron_table;
  delete repl_table;
  delete silent_table;
  delete four_fold_sites_table;
  delete all_silent_table;
  exit(0);
}

void output_statistics(const PolyTable *poly_table,params *args,double n_sites)
{
  PolySNP *analysis = new PolySNP(poly_table,args->haveOutgroup,args->outgroup);
  double w = analysis->ThetaW()/n_sites;
  if (!isfinite(w))
    w = 0.0;
  double pi =  analysis->ThetaPi()/n_sites;
  if (!isfinite(pi))
    pi=0.0;
  if (args->parseableOutput == false)
    {
      cout << "Segsites\t"
	   << analysis->NumPoly() << endl
	   << "Mutations\t"
	   << analysis->NumMutations() << endl
	   << "Singletons\t"
	   << analysis->NumSingletons() << endl
	   << "Num_Sites\t" 
	   << n_sites<<endl
	   << "ThetaW/site\t";
      cout << w << endl;
      cout << "ThetaPi/site"<<'\t';
      cout << pi <<endl;
    }
  else
    {
      cout << analysis->NumPoly() << '\t'
	   << analysis->NumMutations() << '\t'
	   << analysis->NumSingletons() << '\t'
	   << n_sites << '\t'
	   << w << '\t'
	   << pi << '\t';
    }
  delete analysis;
}

void GetSiteCounts(const vector<Fasta > &data, 
		   const PolyTable * all_data,
		   const params *args,
		   double *n_repl_sites,double *n_silent_sites,
		   double * n_third_sites, double * n_4fold_sites,
		   bool useTsTv,double pTs,double pTv)
{
  vector<Fasta> cds(data.size());
  for(unsigned interval = 0 ; interval < args->intervals.size() ; interval += 2)
    {
      unsigned current = args->intervals[interval];
      unsigned next = args->intervals[interval+1]+1;
      for(unsigned seq=0;seq<data.size();++seq)
	{
	  cds[seq].second += string( data[seq].begin()+current,data[seq].begin()+next);
	}
    }

  double m_L0=0.,m_L2S=0.,m_L2V=0.,m_L4=0.;
  int num_comparisons=0;

  for(unsigned i=0;i<cds.size()-1;++i)
    {
      for(unsigned j=i+1;j<cds.size();++j)
	{
	  try
	    {
	      auto_ptr<Comeron95> distance (new Comeron95(&cds[i],&cds[j]));
	      m_L0  += distance->L0();
	      m_L2S += distance->L2S();
	      m_L2V += distance->L2V();
	      m_L4  += distance->L4();
	    }
	  catch(SeqException &e)
	    {
	      cerr << "Processing file "<<args->infile
		   <<" resulted in an exception being thrown.\n"
		   <<"This is likely due to an out-of-frame indel\n"
		   <<endl;
	      e.print(cerr);
	      cerr << endl;
	      exit(1);
	    }
	  ++num_comparisons;
	}
    }
  m_L0 /= double(num_comparisons);
  m_L2S /= double(num_comparisons);
  m_L2V /= double(num_comparisons);
  m_L4 /= double(num_comparisons);
  *n_repl_sites = m_L0 + m_L2S*(2./3.) + m_L2V*(2./3.);
  *n_silent_sites = m_L4 + m_L2S/3. + m_L2V/3.;
  *n_4fold_sites = m_L4;
  unsigned ngapped = cds[0].length() - Sequence::Alignment::UnGappedLength(cds);
  *n_third_sites = round( (*n_repl_sites+*n_silent_sites)/3. )-ngapped;
}

/*
  Old version (incl. in package 0.8.0).  Replaced with above on 1/6/2012.
*/
/*
void GetSiteCounts(vector<Fasta > &data, 
		   const PolyTable * all_data,
		   const params *args,
		   double *n_repl_sites,double *n_silent_sites,
		   double * n_third_sites, double * n_4fold_sites,
		   bool useTsTv,double pTs,double pTv)
{
  vector<Fasta> cds(data.size());

  //the PolyTable all_data is pre-filtered to remove SNPs from
  //problematic codons (codons at which >= 1 position segregates
  //more than 2 states, for example).
  unsigned ngapped=0;
  for(unsigned site=0;site<data[0].length();)
    {
      if(InCoding(site,args->intervals))
	{
	  int cpos = GetCodonPos(site,args->intervals);
	  vector<string>codons(cds.size());
	  string codon;
	  for(unsigned i=0;i<data.size();++i)
	    {
	      MakeCodon(site,i,cpos,data,args->intervals,codon);
	      codons[i]=codon;
	    }

	  bool skip=false;
	  bool site_in_table = false;
	  unsigned nmuts=0;
	  for(unsigned i=0;i<3;++i)
	    {
	      stateCounter c;
	      for(unsigned j=0;j<codons.size();++j)
		c(codons[j][i]);
	      if(c.nStates()>1)
		++nmuts;
	      if(c.nStates()>2)
		skip=true;
	      if (site_in_table == false)
		{
		  site_in_table = (std::find(all_data->pbegin(),
					     all_data->pend(),
					     double(site+1+i)) == all_data->pend());
		}
	    }
	  if(site_in_table == false && nmuts>0)
	    skip = true;
	  if(!skip)
	    {
	      bool gapped=false;
	      for(unsigned i=0;i<codons.size();++i)
		{
		  cds[i].second+=codons[i];
		  if(cpos == 2 && codons[i][2]=='-')
		    {
		      gapped = true;
		    }
		}
	      if(gapped)++ngapped;
	    }
	  switch (cpos)
	    {
	    case 0:
	      site +=3;
	      break;
	    case 1:
	      site+=2;
	      break;
	    case 2:
	      ++site;
	      break;
	    }
	}
      else
	{
	  ++site;
	}
    }
  double m_L0=0.,m_L2S=0.,m_L2V=0.,m_L4=0.;
  int num_comparisons=0;

  for(unsigned i=0;i<cds.size()-1;++i)
    {
      for(unsigned j=i+1;j<cds.size();++j)
	{
	  try
	    {
	      //	      Fasta x("",cds[i]),y("",cds[j]);
	      auto_ptr<Comeron95> distance (new Comeron95(&cds[i],&cds[j]));
	      m_L0  += distance->L0();
	      m_L2S += distance->L2S();
	      m_L2V += distance->L2V();
	      m_L4  += distance->L4();
	    }
	  catch(SeqException &e)
	    {
	      cerr << "Processing file "<<args->infile
		   <<" resulted in an exception being thrown.\n"
		   <<"This is likely due to an out-of-frame indel\n"
		   <<endl;
	      e.print(cerr);
	      cerr << endl;
	      exit(1);
	    }
	  ++num_comparisons;
	}
    }
  m_L0 /= double(num_comparisons);
  m_L2S /= double(num_comparisons);
  m_L2V /= double(num_comparisons);
  m_L4 /= double(num_comparisons);
  if (useTsTv)
    {
      *n_repl_sites = m_L0 + m_L2S*pTv+m_L2V*pTs;
      *n_silent_sites = m_L4 + m_L2S*(1.-pTv)+m_L2V*(1.-pTs);
      *n_4fold_sites = m_L4;
    }
  else
    {
      *n_repl_sites = m_L0 + m_L2S*(2./3.) + m_L2V*(2./3.);
      *n_silent_sites = m_L4 + m_L2S/3. + m_L2V/3.;
      *n_4fold_sites = m_L4;
    }

  *n_third_sites = round( (*n_repl_sites+*n_silent_sites)/3. )-ngapped;
}
*/

void PrintTable(const PolyTable *table, const string &filename,
		const std::vector<string> &seqnames, const string & outgroup_sequence,bool haveOutgroup,
		int outgroup)
  //print out the table in "Hudson2001" spreadsheet format
  //this allow portability to other programs, etc.
{
  ofstream out;
  out.open(filename.c_str());

  out << table->size() << '\t' << table->numsites() << '\n';

  //output SNP positions
  for(unsigned i = 0 ; i < table->numsites() ; ++i)
    out << table->position(i) << '\t';
  out << '\n';
  //output the outgroup characters, if present
  for(unsigned i = 0 ; i < table->numsites() ; ++i)
    {

      if (haveOutgroup == true)
	{
	  if (i == 0)
	    out << outgroup_sequence[string::size_type(table->position(i)-1)];
	  else
	    out << '\t' << outgroup_sequence[string::size_type(table->position(i)-1)];
	}
      else
	{
	  out << '\t'  << '?';
	}
    }
  out << '\n';

  //output the SNPs and sequence names
  for(unsigned i = 0 ; i < table->size() ; ++i)
    {
      out << seqnames[i];
      for (unsigned j = 0 ; j < table->numsites() ; ++j)
	{
	  out << '\t' << (*table)[i][j];
	}
      out << std::endl;
    }
}

void MergeTables(std::vector<double> &newpos, std::vector<string> &newchar,
		 const std::vector<double> &pos1, const std::vector<string> &char1,
		 const std::vector<double> &pos2, const std::vector<string> &char2)
{
  unsigned int i;

  for(i=0;i<pos1.size();++i)
    newpos.push_back(pos1[i]);

  for(i=0;i<pos2.size();++i)
    newpos.push_back(pos2[i]);

  sort(newpos.begin(),newpos.end());

  newchar.resize(char1.size());

  for(i=0;i<newpos.size();++i)
    {
      for(unsigned j = 0 ; j < pos1.size() ; ++j)
	{
	  if(newpos[i]==pos1[j])
	    {
	      for(unsigned k = 0 ; k < newchar.size() ; ++k)
		{
		  newchar[k] += char1[k][j];
		}
	    }
	}
      for(unsigned j = 0 ; j < pos2.size() ; ++j)
	{
	  if(newpos[i]==pos2[j])
	    {
	      for(unsigned k = 0 ; k < newchar.size() ; ++k)
		{
		  newchar[k] += char2[k][j];
		}
	    }
	}
    }
}

unsigned getNonCodingLen(const std::vector<Sequence::Fasta> &data,
			 const std::vector<int> & intervals,
			 const params *args)
//BUGFIX in pacakge 0.7.4 regarding counting gapped sites
{
  unsigned ugl=0;
  for (unsigned site = 0 ; site < data[0].length() ; ++site)
    {
      if (! InCoding(site,intervals) )
	{
	  unsigned ngapped=0;
	  bool gapped = false;
	  for (unsigned seq = 0 ; seq < data.size() ; ++seq)
	    {
	      if ( args->haveOutgroup == false ||
		   (args->haveOutgroup == true && seq != args->outgroup) )
		{
		  if (data[seq][site] == '-')
		    {
		      ++ngapped;
		      //gapped = true;
		      //seq = data.size();
		    }
		}
	    }
	  gapped = (ngapped == data.size()-args->haveOutgroup);
	  ugl += (gapped == false) ? 1 : 0;
	}
    }
  return ugl;
}

void usage(void)
{
  cerr << "polydNdS -i <infile> -I nint i j k l ..."<<endl;
  cerr << "see man polydNdS for details"<<endl;
}


