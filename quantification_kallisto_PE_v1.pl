#######################################################
#
# Script name: quantification_kallisto_PE_v1
# Language: Perl
# Program: Kallisto/quant
# Purpose: Quantification using Kallisto
# Author: Arindam Ghosh
# Date: 16 September 2019
#
#######################################################

#!/usr/bin/perl
use Cwd;

$cwd = getcwd;



open (FH, "/home/dell/Documents/Arindam/WorkB/PRJNA268222/PRJNA268222.txt");

@file = <FH>;

close FH;

for($i=1; $i<@file; $i++)
{
@col = split(/\t/, "$file[$i]");


mkdir $col[1];
chdir "$cwd/$col[1]";



print "\nQuantification for $col[1]\n";

$quant_cmd = "kallisto quant -o . -i /home/dell/Documents/Arindam/MSc2019/Kallisto/GRCh38.97/Homo_sapiens.GRCh38.cdna.all.release-97_k31.idx /home/dell/Documents/Arindam/WorkB/PRJNA268222/Trim/$col[1]/$col[1]_clean_1P.fq.gz /home/dell/Documents/Arindam/WorkB/PRJNA268222/Trim/$col[1]/$col[1]_clean_2P.fq.gz";
system ($quant_cmd);

open (FP, '>', 'quant_cmd.txt');
print FP $quant_cmd;
close FP;


chdir $cwd;
};
