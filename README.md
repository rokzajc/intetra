# intetra
Command-line program for intragenomic oligonucleotide frequency analysis
## INTRODUCTION
### intetra
The program splits the nucleotide sequence stored in fasta file to windows of specified lenghts(-f argument) and count the specified oligonucleotide (-n argument: dinucleotide, trinucleotide, tetranucleotide...) occurances of each window. From the windows' counts, specified statistical scores are calculated(-m argument: z-score, zero'th order Markov model, relative oligonucleotide frequncies), which are used to calculate the matrix of correlations between all windows. Windows can also be generated as sliding windows(-s argument), meaning two adjacent windows will have some overlapping sequence. Program can create windows of varius lenghts in one execution using arguments --maxlen and --minlen.
Using --autocorr argument will calculate the correlations between statistical scores of windows and the whole genome.

### coligo
Program compares the oligonucleotide composition of sequences in fasta files(located in current working directory or another directory). The oligonucleotide composition of different sequences is converted into chosen statistical score(-m argument: z-score, zero'th order Markov model, relative oligonucleotide frequncies) which are used for calculation of correlations between them. Using -n argument, the user can choose the lenght of oligonucleotide words that are counted.

## INSTALLATION
### Linux:
Download ZIP file and extract it anywhere. Open terminal in the directory which was created and run these commands:
```
chmod +x intetra.py
cp intetra.py ~/.local/bin/intetra
cp programi_args ~/.local/bin/programi_args -r
chmod +x coligo.py
cp coligo.py ~/.local/bin/coligo
```

After the installation the coligo.py and intetra.py scripts should be executable from any directory using commands "intetra" and "coligo".

## EXAMPLES
`intetra -i <inputfile.fna> -f 5000 -n 2 -m zom`

`intetra -i <inputfile.fna> -o <outputdirectory> -f 5000 -s 0.5 -n 2 4 -m zom --autocoor`

`intetra -i <inputfile.fna> -o <outputdirectory> -f 3000 -s 2000 -n 6 -m zscr zom --maxlen 300000 --minlen 30000 --autocoor --blockfasta`

`coligo -i <inputdirectory> -o <outputfile> -n 4 5 -m zom zscr -t upgma`


## REQIREMENTS
  python3.6
  biopython
  numpy
  pandas
  matplotlib
