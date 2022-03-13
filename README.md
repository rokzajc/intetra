# intetra
Command-line program for intragenomic oligonucleotide frequency analysis

INTRODUCTION:
The program splits the nucleotide sequence stored in fasta file to windows of specified lenghts(-f argument) and count the specified oligonucleotide (-n argument: dinucleotide, trinucleotide, tetranucleotide...) occurances of each window. From the windows' counts, specified statistical scores are calculated(-m argument: z-score, zero'th order Markov model, relative oligonucleotide frequncies), which are used to calculate the matrix of correlations between all windows. Windows can also be generated as sliding windows(-s argument), meaning two adjacent windows will have some overlapping sequence. Program can create windows of varius lenghts in one execution using arguments --maxlen and --minlen.
Using --autocorr argument will calculate the correlations between statistical scores of windows and the whole genome.

INSTALLATION:
Linux:
Download ZIP file and extract it anywhere. Open terminal in the directory which was created and run these commands:

chmod +x intetra.py
cp intetra.py ~/.local/bin/intetra
cp programi_args ~/.local/bin/programi_args -r
chmod +x coligo.py
cp coligo.py ~/.local/bin/coligo


REQIREMENTS:
  python3.6
  biopython
  numpy
  pandas
  matplotlib
