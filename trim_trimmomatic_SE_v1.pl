#######################################################
#
# Script name: trim_trimmomatic_SE_v1
# Language: Perl
# Program: Trimmomatic
# Purpose: Trim single end RNA-seq data
# Author: Arindam Ghosh
# Date: 4 September 2019
#
#######################################################

#!/usr/bin/perl
use Cwd;

$cwd = getcwd;



open (FH, "/home/dell/Documents/Arindam/WorkB/PRJEB4879/PRJEB4879.txt");

@file = <FH>;

close FH;

for($i=1; $i<@file; $i++)
{
@col = split(/\t/, "$file[$i]");


mkdir $col[1];
chdir "$cwd/$col[1]";


#TRIMMING USING TRIMMOMATIC
print "\nTRIMMING $col[1]\n";

$trim_cmd = "java -jar ~/Documents/Softwares/Trimmomatic/Trimmomatic-0.36/trimmomatic-0.36.jar SE -trimlog $col[1]_trim.log /home/dell/Documents/Arindam/WorkB/PRJEB4879/RawRead/$col[1]/$col[1].fastq.gz $col[1]_clean.fq.gz ILLUMINACLIP:/home/dell/Documents/Softwares/Trimmomatic/Trimmomatic-0.36/adapters/TruSeq3-SE.fa:2:30:10 TRAILING:20 MINLEN:60";
system ($trim_cmd);
system ("fastqc *.gz");

open (FP, '>', 'trim_cmd.txt');
print FP $trim_cmd;
close FP;


chdir $cwd;
};
