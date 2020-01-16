#######################################################
#
# Script name: fastqc_v1
# Language: Perl
# Program: FastQC
# Purpose: Quality control of fastq files
# Author: Arindam Ghosh
# Date: 4 September 2019
#
#######################################################

#!/usr/bin/perl

use Cwd;

$cwd = getcwd;

open (RS, '>', "$cwd/FastqcOut.txt");
open (PL, "$cwd/ProjectList.txt");
@projectlist = <PL>;
close PL;


for($i=0; $i<@projectlist; $i++)
 {
  chomp($projectlist[$i]);
  open (FH, "$cwd/$projectlist[$i]/$projectlist[$i].txt");
  @file = <FH>;
  close FH;
  print RS "\n$projectlist[$i]";

  for($j=1; $j<@file; $j++)
   {
    @col = split(/\t/, "$file[$j]");
    chdir "$cwd/$projectlist[$i]/RawRead/$col[1]";
    $fastqc = 'fastqc';


    $path1 = "$cwd/$projectlist[$i]/RawRead/$col[1]/$col[1]_1.fastq.gz";

    if (-e $path1){
        print "\nPaired End\n";
        $path2 = "$cwd/$projectlist[$i]/RawRead/$col[1]/$col[1]_1_fastqc.html";
        if (-e $path2){
	 print "exists";
         }else{
	system("$fastqc $col[1]_1.fastq.gz");
        print RS "\n\t\t$col[1]_1.fastq.gz";
        }
    
        $path3 = "$cwd/$projectlist[$i]/RawRead/$col[1]/$col[1]_2_fastqc.html";
        if (-e $path3){
	print "exists";
        }else{
	system("$fastqc $col[1]_2.fastq.gz");
        print RS "\n\t\t$$col[1]_2.fastq.gz";
        }

    }else{
        print "\nSingle End\n";
       
        $path4 = "$cwd/$projectlist[$i]/RawRead/$col[1]/$col[1]_fastqc.html";
        if (-e $path4){
	print "exists";
        }else{
	system("$fastqc $col[1].fastq.gz");
        print RS "\n\t\t$col[1].fastq.gz";
        }

    };
chdir $cwd; 
   };
 };

close RS;
