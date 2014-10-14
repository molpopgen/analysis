***analysis README****
Copyright 2002 Kevin Thornton, University of Chicago

NOTE: it appears that there was a bug in polydNdS (likely affecting at least versions 0.7.7-0.8.0) with respect to counting the number of sites in coding regions.  The bug was likely introduced in a relatively recent version, but I cannot track down when.  As a result, the program should not be used until further notics.

* current version 0.8.4 (Nov 9, 2013)  Updated for libsequence 1.7.8, which is required for this version

* version 0.8.3 (Dec. 18, 2012)  Bug fixed:  rsq was outputting too many columns, such that the output didn't match the column headers.

* version 0.8.2 (August 5, 2012)  There was a bug in counting four-fold degenerate sites in polydNdS.  This has now been fixed.  Thanks to Brendan Epstein for discovering the issue.

* version 0.8.1 (Jan 12, 2012).  polydNdS has changed in three ways:  first the site-counting bug alluded to in the above note is now fixed.  Second, the program has zero tolerance for characters other than A,G,C,T, and N (case-insensitive).  This is a reversal of the change made in 0.8.0 (see below), but the program just wasn't written to deal with messy data (it was written to deal with data from Drosophila, where lines are highly inbred).  If you have heterozygous/ambiguous bases in your data, this program is not for you, and this program will not be further developed to support such data.  polydNdS now also requires installation of the boost_regex run-time library.

* version 0.8.0 (Nov 29, 2010).  polydNdS now checks for ambiguous bases (e.g. IUPAC codes), removes them if they are present, and tells you how many SNPs were removed.  

* version 0.7.9 (Sept 2, 2010).  Turns out that version 0.7.8 did not fix the issue with ignoring the -N option.  JJ Emerson figured this out, and I tracked down a second bug that was missed before.  Should work now.

* version 0.7.8 (July 27, 2010).  The option to set max marker distance in rsq should now work.  Some bug fixes involving missing data were taken care of in polydNdS (specifically, a crash when processing 4-fold degenerate sites with missing data, and the -N option to skip alignment columns with missing data should now work).

* version 0.7.7 (June 7, 2010)  Bug fix in "rsq".  For some data sets, the last pair of sites was being skipped.  This should only have affected version 0.7.5 and 0.7.6.

* version 0.7.6 (March 30, 2010)  Removed uncessary printing to STDERR in rsq.

* version 0.7.5 (March 29, 2010).  Several changes in "rsq":  1. is is now much more efficient.  2.  The sign of D is now output  3.  A serious bug was fixed in the permutation test.  Any results based on that should be re-checked.  Requires libsequence 1.6.8 or greater.

* version 0.7.4 (March 24, 2010).  Added an option -m to rsq, which implements a max marker distance. So, -m 1000 will only include markers within 1 kb of each other.  If anyone runs data with positions on a scale from 0 to 1, be sure to adjust the value passed to -m accordingly.

* version 0.7.3 (August 18, 2009).  Fexact now outputs p-values in log scale.  rsq now works correctly when an ancestral state is specified in "snp table" files (i.e., the -h option).  Fasta alignments, however, still cannot have an outgroup.

* version 0.7.2 (Feb 4, 2009)  Updated to compile cleanly on gcc 4.3.x, and on OS X systems where the type "Sequence::arg" conflicts with something in the standard library unless explicit namespaces are used.

* version 0.7.1 (Aug 21, 2008).  Should compile fine on OS X now. There was a problem on some systems where compute_classes.cc failed to compile.

* version 0.7.0 (Jan. 22, 2008).  Changed use of "getopt" so that programs behave properly on Leopard, and are compatible with Tiger and Linux, too.

* version 0.6.9 (Nov. 29, 2007).  Updated code to compile on Apple's "Leopard" operating system.

* version 0.6.8 (Jun 22, 2006).  Snn test now added (snntest).
descPoly now calculates the site-frequency spectrum as either folded or
unfolded, instead of based on the minor allele (which is kind of useless for
the most part...). compute now outputs the number of derived singletons for
data sets with an outgroup.

* version 0.6.7 (Mar 05, 2006) compute now outputs "nan" for undefined
values of Hudson &amp; Kaplan's Rm statistic

* version 0.6.6 (Sept. 21, 2005).  Fixed a bug in writing the output
from compute to a file.  output of compute now more R-friendly.  Requires
libsequence &gt;= 1.5.9

* version 0.6.5 (Jun 6, 2005).  Minor bugfix release.  The program
"compute" in this package would crash on some systems when processing very
large (&gt; 100 kb) alignments with an outgroup.  The problem was a recursive
function.

* 0.6.4 (May 27, 2004) -- updated to libsequence 1.5.4
(required for this release).  This fixes the bugs mentioned above in
calculation of p-value in compute.

* version 0.6.3 (April 29, 2005) -- no changes except that code now
compiles without error under gcc 3.4

* version 0.6.2 (Nov 03, 2004) -- compute now handles outgroups with
value > 1 and doesn't exit with an error.  Also, rsq no longer outputs a
column label "dij"

* version 0.6.1 (Oct 06, 2004) -- compute now reports "NA" for stats with 0
seg sites when appropriate.  Programs dealing with silent and replacement
changes now handle codons with missing data much better.

* version 0.6.0 (Apr 12, 2004) -- fixed bug in compute.  when calculating
pvals, the pval for the Fay-Wu H statistic was incorrect.

* 0.5.9 (Apr 08, 2004) -- updated to latest libsequence (1.3.8).  Made
an improvement to the counting of fixed differences in MKtest. (The change
is just some extra checking, which made no differences in output for the
data sets I looked at.)

* version 0.5.8 (Feb 02, 2004) -- all programs that calculate
polymorphims at silent and replacement sites (polydNdS, MKtest, etc.) are
now based on Sequence::shortestPath, which improved the implementation of
these programs by quite  a bit.  HBKpermute is also updated to a newer,
cleaner version of Sequence::FST.

* version 0.5.7 (Jan 20, 2004) -- The program "MKtest" in versions of
the analysis package up to 0.5.6  generated SNP tables incorrectly.  
Rather amazingly, this was done in a way that had no effect on the outcome for close to 40 data sets (!!).  
This is due to an inherent symmetry of the MKtest, but the bug was still fixed
for version 0.5.7. UPDATE--the bug was a non-issue.  Yes, SNP tables were
made incorrectly, but those tables are just used as guides for processing
the alignment (i.e. actually filling cell counts), which was done correctly.
It was an example of garbage-in, correct results coming out.  Ugh.
	

* version 0.5.6 (Jan 19, 2004) -- configure script now properly detects
all needed headers, including those from GSL.  If GSL headers are not
present, MKtest will not be able to calculate chi-squared or G-test
statistics.

* version 0.5.5 (Jan 12, 2004) -- configure script now checks for
needed headers.  updated to libsequence-1.3.4

* version 0.5.4 (Nov 11, 2003) -- Backend changes to MKtest and polydNdS.
Also, there was a bug in the calculation of the Chi-squared statistic in
MKtest that's now fixed.  The bug affected all previous versions, but that's
ok since everyone uses the Fisher's p-value, right?

* version 0.5.3 (Nov 01, 2003) -- requires libsequence &gt;= 1.3.3.
polydNdS and compute now do more graceful checking to make sure data files
are there.

* 0.5.2 (Aug 22, 2003).  polydNdS now takes a -A option to
provide approximate treatments of codons for the case where there are 2
mutations in the codon and the pathways are ambiguous.  see manpage for
details

* version 0.5.1 (Aug 6, 2003).  sharedPoly now uses Sequence::FST to do
calculations.  manpages updated.  requiress libsequence 1.3.1 or greater to
compile and run.

* vearsion 0.5.0 (Jul 22, 2003) -- sharedPoly now handles sites with &gt;
2 alleles when assigning sites to various classes (fixed, shared, private).

* version 0.4.9 (Jul 15, 2003) -- rewrite of code underlying polydNdS
and MKtest.  The new code base will be much more maintainable in the future,
and will also be easily extensible to handle any problems that should arise.
(Yes, the old code had gotten to be a bit of a mess by 0.4.8) 
MKtest.cc also rewritten and performs much better now.

* version 0.4.8 (Jun 18 2003) -- polydNdSbase.cc was terribly
bug-ridden, in particular for cases where the CDS of the data do not start
and end on a 1st and 3rd position, respectively.  Massive bug-fixing has
taken place, and results have been verified for nearly 40 datasets.  Also,
diagnostic info is output to stderr, so you can see if there are any
potential problems.

* version 0.4.7 (Jun 5 2003) -- fixed bugs in polydNdSbase.cc that
caused SNPs within 1bp of intron/exon boudaries to be mislabelled as coding,
non-coding, or outside the alignment. This would have affected the programs
polydNdS and MKtest in my "analysis" package.

* version 0.4.6 (Jun 3 2003) -- 2 bugs fixed in polydNdS.  The first
bug caused coding sites with &gt; 2 alleles and sites with &gt; 2 SNPs in a codon
to mistakenly be labelled as non-coding SNPs.  The second case affected
samples with only 1 coding SNP--it was omitted from the analysis.  Two new
programs are added: descPoly and sharedPoly

* version 0.4.5 -- compute is now more strict about outputting "NA" in
cases where the data have no segregating sites.  This includes all summary
statistics of the frequency spectrum, and statistics such as Rm.  The reason
for this is to force the user to acknowledge that something special may need
to be done with such numbers.

* version 0.4.4 -- bugfix for "compute".  when there were no seg
sites, and the -p option was used, a segfault occured (due to running sims
with S = 0).  This is now fixed and "NA" is output for p-values for
statistics for such data.  Also, first release of "MKtest," a program to do
the McDonald-Kreitman test

* version 0.4.3 -- bugfix in polydNds.  The bug resulted in first and
second positions being misassigned.  In general, output was correct, but
rare cases (such as a SNP at a redundant first position) resulted in a
silent change being labelled as a replacement change.

* version 0.4.2 -- now compiles on systems without BOOST

* version 0.4.1 -- updated to libsequence >= 1.2.1

* version 0.3.9 -- compatible with libsequence &gt;= 1.1.8

* version 0.3.8 -- now works with g++ 3.1 on OS X

* version 0.3.7 -- added #include directives to compute that fix
compile problems on OS X

* version 0.3.6--updates to new cpptools package, resulting in more
STL-like fun and faster code

* version 0.3.4 -- added D' back to output of rsq, and -c option now
works for rsq.  Requires libsequence &gt;= 1.0.8

* version: 0.3.3 -- removed D' statistic from output of rsq to reflect
changes to libsequence 1.0.7
