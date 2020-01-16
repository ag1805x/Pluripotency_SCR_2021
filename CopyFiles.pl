#!/usr/bin/perl

use Cwd;

$cwd = getcwd;

open (RS, '>', "$cwd/Copy.txt");
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
    chdir "$cwd/$projectlist[$i]/Trim/$col[1]";
    $html_copy = "cp $cwd/$projectlist[$i]/Trim/$col[1]/*_fastqc.html $cwd/AllSampleReports/Trim_FastQC";
    system($html_copy);
    print "\n$html_copy";
    $zip_copy = "cp $cwd/$projectlist[$i]/Trim/$col[1]/*_fastqc.zip $cwd/AllSampleReports/Trim_FastQC";
    system($zip_copy);
    print "\n$zip_copy\n\n";

   
chdir $cwd; 
   };
 };

close RS;
