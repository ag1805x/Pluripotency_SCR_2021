#######################################################
#
# Script name: data_download_PE_v1
# Language: Perl
# Program: Aspera
# Purpose: Download paired end RNA-seq fastq files from ENA
# Author: Arindam Ghosh
# Date: 4 September 2019
#
#######################################################



#!/usr/bin/perl
use Cwd;

$cwd = getcwd;

open (FH, "/home/dell/Documents/Arindam/WorkB/Samples.txt");

@file = <FH>;

close FH;

for($i=1; $i<@file; $i++)
{
@col = split(/\t/, "$file[$i]");
mkdir $col[1];
chdir "$cwd/$col[1]";

@link = split(/;/, "$col[2]");

$file1 = "ascp -QT -P33001 -k 1 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp\@$link[0] .";
$file2 = "ascp -QT -P33001 -k 1 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp\@$link[1] .";

system($file1);
system($file2);

open (FH2, '>', 'cmd.txt');
print FH2 $file1;
print FH2 $file2;
close FH2;
chdir $cwd;
};


