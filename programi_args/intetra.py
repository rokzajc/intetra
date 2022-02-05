from math import gcd
from pickle import FALSE
try:
    import frekvence3,new_frames,z_scores,coorelations2,autocoorelation
except ModuleNotFoundError:
    import programi_args.frekvence3 as frekvence3,programi_args.new_frames as new_frames, programi_args.z_scores as z_scores, programi_args.coorelations2 as coorelations2, programi_args.autocoorelation as autocoorelation
import sys,os
import getopt
navodila='python3 intetra.py -i <inputfile> -o <outputfile>  '
try:
    opts, args = getopt.getopt(sys.argv[1:],"hi:o:f:s:m:",['ifile=','ofile=','framelen=',"slidelen=","increase_slide=",'maxlen=','method=','autocoor','minlen=','coorelations=','blockfasta'])
except:
    pass
slide_len=False
frame_len=False
inc_sld=False
autocoor=False
minlen=False
blockfasta=False
coor='fast'
method='zscr'
output=''
maxlen=0
for opt, arg in opts:
    if opt == '-h':
        print (navodila)
        sys.exit(2)
    elif opt in ("-i", "--ifile"):
        input = arg
    elif opt in ("-o", "--ofile"):
        output= os.path.join(os.getcwd(),arg)
    elif opt in ("-s", "--slidelen"):
        slide_len = int(arg)
    elif opt in ("-f", "--framelen"):
        frame_len= int(arg)
    elif opt =="--maxlen":
        maxlen = int(arg)
    elif opt in ("-m", "--method"):
        method= arg
    elif opt == ("--autocoor"):
        autocoor= True    
    elif opt == ("--increase_slide"):
        inc_sld=int(arg)
    elif opt=='--minlen':
        minlen=int(arg)
    elif opt=='--coorelations':
        coor=arg
    elif opt=='--blockfasta':
        blockfasta=True



if frame_len:
    if slide_len:
        read_by=gcd(frame_len,slide_len)
    else:
        read_by=frame_len
        slide_len=frame_len
if output=='':
    output='Output_'+input
if maxlen==False:
    maxlen=frame_len
if inc_sld!=False:
    inc_sld=int(inc_sld/slide_len)*slide_len
else:
    inc_sld=slide_len
if minlen==False:
    minlen=frame_len
#coorelations2.main(input,output,method,coor,frame_len,slide_len,inc_sld)
#new_frames.main(input,output,frame_len,read_by,slide_len,inc_sld,maxlen,minlen,blockfasta)
#autocoorelation.main(input,output,frame_len,slide_len,inc_sld,maxlen,minlen,method,)
try: 
    frekvence3.main(input,output,read_by,blockfasta)
except:
    pass
try:
    new_frames.main(input,output,frame_len,read_by,slide_len,inc_sld,maxlen,minlen,blockfasta)
except:
    pass
try:
    z_scores.main(output,method)
except:
    pass
try:
    coorelations2.main(input,output,method,coor,frame_len,slide_len,inc_sld)
except:
    pass
if autocoor==True:
    autocoorelation.main(input,output,frame_len,slide_len,inc_sld,maxlen,minlen,method,)