#######################################################
#
# Script name: data_download_SE_v1
# Language: Perl
# Program: Aspera
# Purpose: Download single end RNA-seq fastq files from ENA
# Author: Arindam Ghosh
# Date: 4 September 2019
#
#######################################################

#!/usr/bin/perl
use Cwd;

$cwd = getcwd;

open (FH, "/home/dell/Documents/Work/PRJNA420980/RawRead/PRJNA420980_cont.txt");

@file = <FH>;

close FH;

for($i=1; $i<@file; $i++)
{
@col = split(/\t/, "$file[$i]");
mkdir $col[1];
chdir "$cwd/$col[1]";


$file = "ascp -QT -P33001 -k 1 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp\@$col[2] .";
system($file);

open (FH2, '>', 'cmd.txt');
print FH2 $file;
close FH2;
chdir $cwd;
};
