#######################################################
#
# Script name: fastq_screen_v1
# Language: Perl
# Program: FastQ Screen
# Purpose: Quality control of fastq files using FastQ Screen
# Author: Arindam Ghosh
# Date: 6 September 2019
#
#######################################################

#!/usr/bin/perl

use Cwd;

$cwd = getcwd;

open (RS, '>', "$cwd/FastqScreenOut.txt");
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
    $fastqscreen = 'fastq_screen';



    $path = "$cwd/$projectlist[$i]/RawRead/$col[1]/$col[1]_1.fastq.gz";

    if (-e $path){
        print "\nPaired End\n";
	$file1 = "$fastqscreen $col[1]_1.fastq.gz --subset 0 --threads 10 --outdir $cwd/$projectlist[$i]/RawRead/$col[1]/$col[1]_fastqscreen";
        system($file1);
        print RS "\n\t\t$col[1]_1.fastq.gz";
	$file2 = "$fastqscreen $col[1]_2.fastq.gz --subset 0 --threads 10 --outdir $cwd/$projectlist[$i]/RawRead/$col[1]/$col[1]_fastqscreen";
       	system($file2);
        print RS "\n\t\t$col[1]_2.fastq.gz";
    }else{
        print "\nSingle End\n";
	$file = "$fastqscreen $col[1].fastq.gz --subset 0 --threads 10 --outdir $cwd/$projectlist[$i]/RawRead/$col[1]/$col[1]_fastqscreen";       	
	system($file);
        print RS "\n\t\t$col[1].fastq.gz";
    };
chdir $cwd; 
   };
 };

close RS;
