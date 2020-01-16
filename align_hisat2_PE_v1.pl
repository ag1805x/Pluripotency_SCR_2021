#######################################################
#
# Script name: align_hisat2_PE_v1
# Language: Perl
# Program: Hisat2
# Purpose: Align paired end RNA-seq data to reference genome
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



print "\nALIGNING  $col[1]\n";

$hisat2_cmd = "hisat2 -p 12 --dta --rna-strandness RF -x /home/dell/Documents/Arindam/Work/ReferenceGenome/Human_84/hisat2_index/grch38_tran/genome_tran -1 /home/dell/Documents/Arindam/WorkB/PRJEB7132/Trim/$col[1]/$col[1]_clean_1P.fq.gz  -2 /home/dell/Documents/Arindam/WorkB/PRJEB7132/Trim/$col[1]/$col[1]_clean_2P.fq.gz -S $col[1].sam --un-gz $col[1]_unaln.fq.gz --summary-file $col[1]_sum.txt --met-file $col[1]_met.tsv";


system ($hisat2_cmd);

open (FP, '>', 'align_cmd.txt');
print FP $hisat2_cmd;
close FP;


print "\nSAM to BAM conversion\n";
$samtools = 'samtools';
system ("$samtools sort -o $col[1].bam $col[1].sam ");

print "\nQC for BAM\n";
system ("fastqc *.bam");




chdir $cwd;
};
