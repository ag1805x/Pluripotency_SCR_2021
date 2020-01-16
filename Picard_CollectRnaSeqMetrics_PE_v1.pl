#######################################################
#
# Script name: Picard_CollectRnaSeqMetrics_PE_v1
# Language: Perl
# Program: Picard/CollectRnaSeqMetrics
# Purpose: Quality control of bam files from paired end sequenceing data
# Author: Arindam Ghosh
# Date: 4 September 2019
#
#######################################################

#!/usr/bin/perl
use Cwd;

$cwd = getcwd;



open (FH, "/home/dell/Documents/Arindam/WorkB/PRJEB7132/PRJEB7132.txt");

@file = <FH>;

close FH;

for($i=1; $i<@file; $i++)
{
@col = split(/\t/, "$file[$i]");


mkdir $col[1];
chdir "$cwd/$col[1]";


print "\nCollecting Alignment Metrics for  $col[1]\n";

$picard_cmd = "java -jar /home/dell/Documents/Softwares/Picard/picard.jar CollectRnaSeqMetrics I=/home/dell/Documents/Arindam/WorkB/PRJEB7132/Align/$col[1]/$col[1].bam O=$cwd/$col[1]/$col[1]_Picard.RNAmetrics REF_FLAT=/home/dell/Documents/Arindam/Work/ReferenceGenome/Human_84/refFlat/refFlat.txt.gz STRAND_SPECIFICITY=NONE RIBOSOMAL_INTERVALS=/home/dell/Documents/Arindam/Work/ReferenceGenome/Human_84/Ribosomal_Intervals/GRCh38.rRNA.interval_list";


system($picard_cmd);

open (FP, '>', 'align_cmd.txt');
print FP $picard_cmd;
close FP;

chdir $cwd;
};


