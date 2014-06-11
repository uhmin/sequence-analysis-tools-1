Sequence analysis tools 1

Installation
type
make.sh

ORF-finder
==========
Just a small C program that can find ORFs from nucleotide sequences in (multi-)fasta format.

revcmp
==========
Just a small C program that can reverse/complement nucleotide sequences in (multi-)fasta format.

translate
==========
Just a small C program that can translate nucleotide sequence to amino acid sequence.

docono
==========
Just a small C program that can find short sequences from (multi-)fasta nucleotide sequence.

filterN
==========

convert BLAST text formatted result into table format

function : Convert blast textbase output into tabular format.
swhitches
       -i : infile name. stdin if default.
       -o : outfile name. stdout if default.
       -e : errmessege output. stderr if default.
       -v : message level (0-2). 1 if default.
       -Q : number of query  sequence  to show. If 0 show all. 0 if default
       -D : number of matched database to show. If 0 show all. 0 if default
       -A : number of aligments to show.        If 0 show all. 0 if default
       -S : allow Same score [T/F]. F if default. If select T allow more than -D specified number if the result show the Same score.


***** NOTICE ******
This program may find BLAST format err.
If format error were found, this program will stop processing.
If you DO WANT TO FORCE processing, try '-f T' option.
