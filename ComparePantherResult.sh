#######################################################
#
# Script name: ComparePantherResult.sh
# Language: Bash
# Program: Bash
# Purpose: Extract unique terms from two PANTHER output files and prepare for REVIGO
# Author: Arindam Ghosh
# Date: 13 February 2020
#
#######################################################


file1="GO_BPslim_TurquoiseMod.txt"
file2="GO_BPslim_BlueMod.txt"
tail -n +12  $file1 | cut -f1 | while read pat; do grep "^$pat" $file2; done | cut -f1 > Common.txt
tail -n +12 $file1 | grep -vf Common.txt > "Unique_${file1}"
tail -n +12 $file2 | grep -vf Common.txt > "Unique_${file2}"
rm Common.txt
paste <(cut -f1 "Unique_${file1}" | cut -d "(" -f2 | cut -d ")" -f1) <(cut -f7 "Unique_${file1}")  >> "RevigoIn_Unique_${file1}"
paste <(cut -f1 "Unique_${file2}" | cut -d "(" -f2 | cut -d ")" -f1) <(cut -f7 "Unique_${file2}")  >> "RevigoIn_Unique_${file2}"
