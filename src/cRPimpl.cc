#include <Sequence/Translate.hpp>
#include <Sequence/shortestPath.hpp>

#include <cRPimpl.hpp>
#include <polydNdSbase.hpp>
#include <CodonMaker.hpp>

#include <limits>
#if defined(HAVE_SSTREAM)
#include <sstream>
typedef std::ostringstream _ostr;
#elif defined(HAVE_STRSTREAM)
#include <strstream>
typedef std::ostrstream _ostr;
#else
#error
#endif
#include <iostream>
using namespace Sequence;

std::vector<int> nextSNPs(const int & num_next_snps,
			  const unsigned & index,
			  const vector<double> &snp_positions,
			  const vector<int> &intervals)
{
  vector<int> rv;
  for(unsigned i = index+1, count = 0 ; 
      i < snp_positions.size() && count < unsigned(num_next_snps); ++i)
    {
      if (InCoding( int(snp_positions[i]-1), intervals ) == true)
	{
	  rv.push_back(int(snp_positions[i]-1));
	  ++count;
	}
    }
  return rv;
}

    
bool cRPimpl::isStop(const std::string &codon)
/*!
  \param args the arguments passed to the program
  \param codon a string of length 3 representing a codon
  \param indexed_pos the index of the SNP position along the sequence
*/
{
  return  (Translate(codon.begin(),codon.end()) == "*"); 
}
    
void cRPimpl::addToTable(const vector<Fasta> &data, 
			 const std::vector<double> & snp_positions,
			 const int &indexed_pos, vector<double> &pos, 
			 vector<std::string> &matrix)
{
  bool present = false;
  for (unsigned l = 0; l < pos.size(); l++)
    {
      if (pos[l] == indexed_pos+1)
	present = true;
    }
		
  if (!present)
    {
		
			
      pos.push_back(indexed_pos+1);
      for(unsigned k = 0 ;k<data.size();++k)
	{
	  if( !haveOutgroup
	      || (haveOutgroup && k != outgroup) )
	    {
	      if (!haveOutgroup)
		matrix[k]+=data[k][indexed_pos];
	      else
		{
		  if(k>outgroup)
		    matrix[k-1]+=data[k][indexed_pos];
		}
	    
	    }
	}
    }
}
    
/*unsigned cRPimpl::incrementSiteIndex(const unsigned &nDiffsCodon,const unsigned &cPos)
  {
  switch (nDiffsCodon)
  {
  case 1:
  return 0;
  break;
  case 2:
  switch (cPos)
  {
  case 0:
  return 1;
  break;
  case 1:
  return 1;
  break;
  case 2:
  return 0;
  break;
  }
  break;
  case 3:
  switch (cPos)
  {
  case 0:
  //should be the only one that ever occurs
  return 2;
  break;
  case 1:
  return 1;
  break;
  case 2:
  return 0;
  break;
  }
  break;
  }
  return 0;
  }
*/
cRPimpl::cRPimpl(const std::vector< Sequence::Fasta > &data,
		 const std::vector<double> & snp_positions,
		 const std::vector< int > & _intervals,
		 const bool & _haveOutgroup,
		 const unsigned &_outgroup,
		 const bool & _Approximate) 
  throw (std::exception) :
  intervals(_intervals),
  haveOutgroup(_haveOutgroup),
  outgroup(_outgroup),
  Approximate(_Approximate)
{
  process_coding_region(data,snp_positions);
}

void cRPimpl::process_coding_region(const std::vector< Sequence::Fasta > &data,
				    const std::vector<double> & snp_positions) throw (std::exception)
{
  if(haveOutgroup)
    {
      silent_char.resize(data.size()-1);
      repl_char.resize(data.size()-1);
    }
  else
    {
      silent_char.resize(data.size());
      repl_char.resize(data.size());
    }
  for(unsigned int i = 0 ; i < snp_positions.size();++i)
    //iterate over polymorphic sites
    {
      int indexed_pos = int(snp_positions[i]-1);
      if(InCoding(indexed_pos,intervals))
	//then it is a coding polymorphism
	{
	  int cPos = GetCodonPos(indexed_pos,intervals);
	  if (cPos == -1)
	    {
	      _ostr o;
	      o << "The codon position of the SNP could\n"
		<< "not be determined.  This is bad and should\n"
		<< "never happen.  This may be a bug\n";
	      throw (std::runtime_error(o.str()));
	    }

	  CodonMaker process(data,intervals,indexed_pos,cPos);
	  string t1 = process.c1(),t2=process.c2();
	  
	  //unsigned nDiffsCodon = process.ndiffs();
	  if (process.Success()==true)
	    {
	      if( isStop(process.c1()) || isStop(process.c2()) )
		{
		  _ostr o;
		  o << "Found a polymorphic stop codon ("
		    << process.c1()
		    <<", "
		    << process.c2()
		    <<") at or near position " 
		    << indexed_pos+1 << '.';
		  w.add(o.str());
		}
	      // Step 1: figure out positions in alignment of all SNPs in the current codon
	      int indexed_pos_1,indexed_pos_2,indexed_pos_3;
	      //indexed_pos_1 = indexed_pos_2 = indexed_pos_3 = indexed_pos; //was set to min int, but that may have caused problems
	      //does this work?
	      indexed_pos_1 = indexed_pos - cPos;
	      indexed_pos_2 = indexed_pos_1 + 1;
	      indexed_pos_3 = indexed_pos_2 + 1;
	      /* 
		 vector<int> next_snp_indexes;
		 int cPos2;//,cPos3;
		 if (cPos == 0)//current site is a 1st position
		 {
		 switch(nDiffsCodon)
		 {
		 case 1: 
		 //only 1 SNP in codon, the current one
		 indexed_pos_1 = indexed_pos;
		 break;
		 case 2: 
		 //there is an additional SNP in the codon,
		 //so we go find the next SNP that is InCoding()
		 indexed_pos_1 = indexed_pos;
		 next_snp_indexes = nextSNPs(1,i,snp_positions,intervals);
		 cPos2 = GetCodonPos(next_snp_indexes[0],intervals);
		 if (cPos2==1)//a second pos
		 indexed_pos_2 = next_snp_indexes[0];
		 else if (cPos2==2)//a 3rd pos
		 indexed_pos_3 = next_snp_indexes[0];
		 break;
		 case 3:
		 indexed_pos_1 = indexed_pos;
		 next_snp_indexes = nextSNPs(2,i,snp_positions,intervals);
		 indexed_pos_2 = next_snp_indexes[0];
		 indexed_pos_3 = next_snp_indexes[1];
		 break;
		 };
		 }
		 else if (cPos == 1)//current site is a 2nd position
		 {
		 switch(nDiffsCodon)
		 {
		 case 1:
		 indexed_pos_2 = indexed_pos;
		 break;
		 case 2:
		 indexed_pos_2 = indexed_pos;
		 next_snp_indexes = nextSNPs(1,i,snp_positions,intervals);
		 indexed_pos_3 = next_snp_indexes[0];
		 break;
		 case 3:
		 //not possible
				
		 break;
		 };
		 }
		 else if (cPos == 2)//current site is a 3rd position
		 {
		 switch(nDiffsCodon)
		 {
		 case 1:
		 indexed_pos_3 = indexed_pos;
		 break;
		 case 2:
		 //not possible
		 break;
		 case 3:
		 //not possible
		 break;
		 };
		 }
	      */ 
	      shortestPath sp(process.c1(),process.c2());
	      vector<string> pathway(sp.begin(),sp.end());
	      for (unsigned cod_in_path = 1 ; cod_in_path < pathway.size() ; ++cod_in_path)
		{
		  std::pair<unsigned,shortestPath::pathType> pdata = 
		    diffType(pathway[cod_in_path-1],pathway[cod_in_path]);
		  switch (pdata.first) //check codon position of change
		    {
		    case 0: //1st
		      if (pdata.second == shortestPath::N)
			{
			  addToTable(data,snp_positions,indexed_pos_1,
				     repl_pos,repl_char);
			}
		      else if (pdata.second == shortestPath::S)
			{
			  addToTable(data,snp_positions,indexed_pos_1,
				     silent_pos,silent_char);
			}
		      break;
		    case 1: //2nd
		      if (pdata.second == shortestPath::N)
			{
			  addToTable(data,snp_positions,indexed_pos_2,
				     repl_pos,repl_char);
			}
		      else if (pdata.second == shortestPath::S)
			{
			  addToTable(data,snp_positions,indexed_pos_2,
				     silent_pos,silent_char);
			}
		      break;
		    case 2: //3rd
		      if (pdata.second == shortestPath::N)
			{
			  addToTable(data,snp_positions,indexed_pos_3,
				     repl_pos,repl_char);
			}
		      else if (pdata.second == shortestPath::S)
			{
			  addToTable(data,snp_positions,indexed_pos_3,
				     silent_pos,silent_char);
			}
		      break;
		    }
		}
	      // BELOW IS THE OLD IMPLEMENTATION.  EMBRACE IT'S UGLINESS!
	      // 	      if (nDiffsCodon == 1)
	      // 		{
	      // 		  if (Translate(codon1.begin(),codon1.end()) !=
	      // 		      Translate(codon2.begin(),codon2.end()))
	      // 		    {
	      // 		      addToTable(data,snp_positions,indexed_pos,repl_pos,repl_char);
	      // 		    }
	      // 		  else
	      // 		    {
	      // 		      addToTable(data,snp_positions,indexed_pos,silent_pos,silent_char);
	      // 		    }
	      // 		}
	      // 	      else if (nDiffsCodon == 2)
	      // 		{
	      // 		  string intermediates[2];
	      // 		  Intermediates2(intermediates,codon1,codon2);
	      
	      // 		  string refTrans   = Translate(codon1.begin(),codon1.end());
	      // 		  string int1Trans  = Translate(intermediates[0].begin(),intermediates[0].end());
	      // 		  string int2Trans  = Translate(intermediates[1].begin(),intermediates[1].end());
	      // 		  string codonTrans = Translate(codon2.begin(),codon2.end());
	      
	      // 		  //Note: Sequence::Comparisons::NumDiffs(), which was
	      // 		  //used to identify codon2Diffs does not
	      // 		  //allow differences to be counted if the difference
	      // 		  //is the missing data character 'N', so the simplification
	      // 		  //below (i.e. not checking for missing data) is correct
	      // 		  mutantType branch1=UNDETERMINED,branch2=UNDETERMINED,
	      // 		    branch3=UNDETERMINED,branch4=UNDETERMINED;
	      // 		  if ((refTrans != int1Trans) && (refTrans != "X" && int1Trans != "X") )
	      // 		    {
	      // 		      branch1 = REPLACEMENT;
	      // 		    }
	      // 		  else if ( (refTrans != "X" && int1Trans != "X")
	      // 			    && (refTrans == int1Trans) )
	      // 		    {
	      // 		      branch1 = SILENT;
	      // 		    }
	      
	      // 		  if ( (int1Trans != codonTrans) && (int1Trans != "X" && codonTrans != "X") )
	      // 		    {
	      // 		      branch2 = REPLACEMENT;
	      // 		    }
	      // 		  else if ( (int1Trans == codonTrans) && (int1Trans != "X" && codonTrans != "X") )
	      // 		    {
	      // 		      branch2 = SILENT;
	      // 		    }
	      
	      // 		  if ( (refTrans!=int2Trans) && (refTrans!="X" && int2Trans!="X") )
	      // 		    {
	      // 		      branch3 = REPLACEMENT;
	      // 		    }
	      // 		  else if  ( (refTrans==int2Trans) && (refTrans!="X" && int2Trans!="X") )
	      // 		    {
	      // 		      branch3 = SILENT;
	      // 		    }
	      
	      // 		  if ( (int2Trans != codonTrans) && (int2Trans!="X" && codonTrans!="X") )
	      // 		    {
	      // 		      branch4 = REPLACEMENT;
	      // 		    }
	      // 		  else if ( (int2Trans == codonTrans) && (int2Trans!="X" && codonTrans!="X") )
	      // 		    {
	      // 		      branch4 = SILENT;
	      // 		    }
	      
	      // 		  //The definition of Sequence::PathwayHelper::Intermediates2
	      // 		  //guarantees that intermediates[0] differs from codon1
	      // 		  //at the leftmost SNP in the codon, and that intermediates[1]
	      // 		  //differs from codon1 at the rightmost SNP.  Therefore,
	      // 		  //the assignment of SILENT/REPLACEMENT to these SNPS is only
	      // 		  //independent of order (i.e. which SNP occurred 1st, evolutionarily 
	      // 		  //speaking) if branch1==branch4 && branch2==branch3.  Draw a picture,
	      // 		  //it helps...
	      // 		  if(allDetermined(branch1,branch2,branch3,branch4))
	      // 		    {
	      // 		      if ( branch1==branch4 && branch2==branch3 )
	      // 			{
	      // 			  if (branch1 == REPLACEMENT)
	      // 			    {
	      // 			      //the i-th SNP is a replacement poly
	      // 			      addToTable(data,snp_positions,indexed_pos,
	      // 					 repl_pos,repl_char);
	      // 			    } 
	      // 			  else if (branch1 == SILENT)
	      // 			    {
	      // 			      //the i-th SNP is a silent poly
	      // 			      addToTable(data,snp_positions,indexed_pos,
	      // 					 silent_pos,silent_char);
	      // 			    }
		      
	      // 			  int next_indexed_pos = int(snp_positions[i+1])-1;
	      // 			  if (branch2 == REPLACEMENT)
	      // 			    {
	      // 			      //the i+1-th SNP is a replacement poly
	      // 			      addToTable(data,snp_positions,next_indexed_pos,
	      // 					 repl_pos,repl_char);
	      // 			    } 
	      // 			  else if (branch2 == SILENT)
	      // 			    {
	      // 			      //the i+1-th SNP is a silent poly
	      // 			      addToTable(data,snp_positions,next_indexed_pos,
	      // 					 silent_pos,silent_char);
	      // 			    }
	      // 			}
	      // 		      else
	      // 			{
	      // 			  if (Approximate == false)
	      // 			    {
	      // 			      _ostr o;
	      // 			      o << "/////////////////\n"
	      // 				<< " Ambiguous codons: " << codon1 
	      // 				<< " and " << codon2
	      // 				<< " were encountered at or near position "<<indexed_pos+1<<".\n"
	      // 				<< "This means that whether or not the changes\n"
	      // 				<< "are silent or replacement depends on the order\n"
	      // 				<< "in which they occur, and so they are not considered\n"
	      // 				<< "in this analysis.\n"
	      // 				<< "/////////////////";
	      // 			      w.add(o.str());
	      // 			    }
	      // 			  else if (Approximate == true)
	      // 			    {
	      // 			      if ( (branch1 == SILENT || branch2 == SILENT)
	      // 				   && (branch3==REPLACEMENT && branch4==REPLACEMENT) )
	      // 				{
	      // 				  //the treat the first position as a silent change, 
	      // 				  //the second as a replacement
	      // 				  addToTable(data,snp_positions,indexed_pos,
	      // 					     silent_pos,silent_char);
	      // 				  addToTable(data,snp_positions,indexed_pos+1,
	      // 					     repl_pos,repl_char);
	      // 				}
	      // 			      else if  ( (branch3 == SILENT || branch4 == SILENT)
	      // 					 && (branch1==REPLACEMENT && branch2==REPLACEMENT) )
	      // 				{
	      // 				  //the treat the first position as a replacement change, 
	      // 				  //the first as a replacement
	      // 				  addToTable(data,snp_positions,indexed_pos,
	      // 					     repl_pos,repl_char);
	      // 				  addToTable(data,snp_positions,indexed_pos+1,
	      // 					     silent_pos,silent_char);
	      // 				}
	      // 			      else
	      // 				{
	      // 				  //resort to translating the codons
	      // 				  if(Translate(codon1.begin(),codon1.end()) 
	      // 				     != Translate(codon2.begin(),codon2.end()))
	      // 				    {
	      // 				      addToTable(data,snp_positions,indexed_pos,
	      // 						 repl_pos,repl_char);
	      // 				    }
	      // 				  else
	      // 				    {
	      // 				      addToTable(data,snp_positions,indexed_pos,
	      // 						 silent_pos,silent_char);
	      // 				    }
	      // 				}
	      // 			    }
	      // 			}
	      // 		    }//if(allDetermined())
	      // 		}//if (nDiffsCodon == 2)
	      // 	      else if (nDiffsCodon == 3)
	      // 		{
	      // 		  if (Approximate == false)
	      // 		    {
	      // 		      _ostr o;
	      // 		      o << "A codon at or near position " << indexed_pos+1
	      // 			<< " appears to have mutations at all 3 positions\n"
	      // 			<< "and is being skipped.";
	      // 		      w.add(o.str());
	      // 		    }
	      // 		  else if (Approximate == true)
	      // 		    {
	      // 		      //not dealt with yet
	      // 		      _ostr o;
	      // 		      o << "A codon at or near position " << indexed_pos+1
	      // 			<< " appears to have mutations at all 3 positions\n"
	      // 			<< "and is being skipped";
	      // 		      w.add(o.str());
	      // 		    }
	      // 		}

	      //i may need to be incremented by more than
	      //1 if there is > 1 mutation in the codon
	      //this switch statement handles such cases
	      //      i += incrementSiteIndex(nDiffsCodon,cPos);
	    }//if either codon is "NULL"
	  else
	    {
	      _ostr o;
	      o << "The codons "<<process.c1()<<" and "<<process.c2()
		<< " at or near position "<< indexed_pos+1 << ' '
		<< "were entirely ambiguous\n"
		<< "and cannot be processed.";
	      w.add(o.str());
	    }
	}//if (InCoding(...))
    }//iterate over sites
}
