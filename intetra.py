#!/usr/bin/env python3

#Author: Rok Zajc

print('INTETRA - Intragenomic oligonucleotide frequency analysis tool\n----------------------------------------------------------\nAuthor: Rok Zajc\n\n')
from math import gcd
import programi_args.frekvence3 as frekvence3, programi_args.new_frames as new_frames, programi_args.z_scores as z_scores, programi_args.coorelations2 as coorelations2, programi_args.autocoorelation as autocoorelation
import sys
import argparse


parser = argparse.ArgumentParser(description='Slice fasta file into smaller windows and calculate correlations between them based on their oligonucleotide frequences. The tool can also create sliding windows by using -s argument. Using --maxlen will create windows of multiple lenghts.')
parser.add_argument('-i', help='input file (must be fasta format)', type=str, dest='input', required=True )
parser.add_argument('-o', help='output folder (folder will be created)', type=str, dest='output', default=f'Output')
parser.add_argument('-f', help='lenght of windows by which genome is sliced', type=int, dest='frame_len', required=True )
parser.add_argument('-s', help='lenght by which a window is slid(Default=frame_len), can be expressed as ratio of -f argument(example: 0.5)', type=str, dest='slide_len', default='')
parser.add_argument('--maxlen', help='maximum lenght of sliding windows(Default=frame_len)', type=int, dest='maxlen', default=-1)
parser.add_argument('--minlen', help='minimum lenght of sliding windows(Default=frame_len)', dest='minlen', type=int, default=0)
parser.add_argument('-m', help='method of calculating frequncies, multiple methods can be chosen (Default=zscr)', dest='method', type=str, nargs='*', default=['zscr'],choices=['zscr', 'zom', 'rof'])
parser.add_argument('-n', help='lenght of nucleotide words whose occurrences in the sequence are counted(default=4)', type=int, dest='nuc_number', nargs='*', default=[4])
parser.add_argument('--increase_slide', help='lenght by which the slide is increased, multiple lenghts can be chosen(default=slidelen)', type=int, dest='inc_sld', default=-1)
parser.add_argument('--blockcorr', help='turn off calculation of correlations between windows', dest='blockcorr', action='store_true')
parser.add_argument('--autocorr', help='turn on autocorrelation',dest='autocorr', action='store_true')
parser.add_argument('--blockfasta', help='new fasta files are not created', dest='blockfasta', action='store_true')
args = parser.parse_args()


if args.slide_len !='' :
    if '.' in args.slide_len:
        slide_len=round(float(args.slide_len)*args.frame_len)
    else:
        try:
            slide_len=int(args.slide_len)
        except:
            print('Error: -s must be an intrger or float')
            sys.exit()
    if slide_len>args.frame_len:
        print('Error: Slide lenght(-s) must be smaller than frame lenght(-f)')
        sys.exit()
    read_by=gcd(args.frame_len,slide_len)
    if read_by<100:
        print('Error: Greatest common divisor between slide lenght(-s) and frame lenght(-f) must be at least 100.')
        sys.exit()
else:
    read_by=args.frame_len
    slide_len=args.frame_len
if (2 in args.nuc_number or 1 in args.nuc_number) and 'zscr'in args.method:
    print('Error: Cannot calculate z-score for -n < 3')
    sys.exit()
if args.maxlen==-1:
    args.maxlen=args.frame_len
if args.inc_sld!=-1:
    inc_sld=int(args.inc_sld/slide_len)*slide_len
else:
    inc_sld=slide_len
if args.minlen==-1:
    args.minlen=args.frame_len
#for i in range(10):
#    frekvence3.main(args.input, args.output, read_by, max(args.nuc_number), args.blockfasta)
#new_frames.main(args.input, args.output, args.frame_len, read_by, slide_len, inc_sld, args.maxlen, args.minlen, max(args.nuc_number), args.blockfasta)#z_scores.main(args.output, args.method)
#z_scores.main(args.output, args.method, args.nuc_number)
#coorelations2.main(args.input, args.output, args.method,args.nuc_number, args.frame_len, slide_len, inc_sld, args.minlen, args.maxlen)
#autocoorelation.main(args.input, args.output, args.frame_len, slide_len, inc_sld, args.maxlen, args.minlen, args.method, args.nuc_number)

print('Please wait...')
frekvence3.main(args.input, args.output, read_by, max(args.nuc_number), args.blockfasta)
new_frames.main(args.input, args.output, args.frame_len, read_by, slide_len, inc_sld, args.maxlen, args.minlen, max(args.nuc_number), args.blockfasta)
z_scores.main(args.output, args.method, args.nuc_number)
if args.blockcorr==False:
    coorelations2.main(args.input, args.output, args.method,args.nuc_number, args.frame_len, slide_len, inc_sld, args.minlen, args.maxlen)
if args.autocorr==True:
    autocoorelation.main(args.input, args.output, args.frame_len, slide_len, inc_sld, args.maxlen, args.minlen, args.method, args.nuc_number)
