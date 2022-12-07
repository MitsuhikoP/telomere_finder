#!/usr/bin/env python3
# copyright (c) 2022 Mitsuhiko P. Sato. All Rights Reserved.
# Mitsuhiko Sato ( E-mail: mitsuhikoevolution@gmail.com )
    
def rev_com(seq):
    return(seq.translate(str.maketrans("ATGCatgc", "TACGtagc"))[::-1])

def find_telomere(seq, telomere, INTERVAL_LIMIT):
    counts=0
    pos=seq.find(telomere)
    first_pos=pos
    spos=0
    while spos < INTERVAL_LIMIT and spos >= 0:
        pos+=spos
        count=0
        while seq[pos:pos+len(telomere)] == telomere:
            count+=1
            pos+=len(telomere)
        spos=seq[pos:].find(telomere)
        counts+=count
    #print(ID,telomere, pos, counts)
    return(first_pos, pos, counts)


def main():
    from argparse import ArgumentParser
    parser=ArgumentParser(description="",usage="python3 telomere_finder.py -i input.fasta", epilog="OUTPUT\nfasta_ID\ttelomere_sequence\tF/R\tDist_from_edge\tnum_telomere\n")
    parser.add_argument("-i",required=True, type=str,metavar="str",help="input file")
    parser.add_argument("-e", type=int,metavar="int",default=100, help="Length of the end that is skipped when non-telomere sequences (default=100)")
    parser.add_argument("-a", type=int,metavar="int",default=28, help="Distance to allow not telomere sequences (default=28)")
    parser.add_argument("-l", type=int,metavar="int",default=100000, help="Minimum contig length (default=100000)")
    parser.add_argument("-t", type=str,metavar="comma separated str",default="CCCTAAA,TTTAGGG", help="telomere sequence. comma separated. (default=CCCTAAA,TTTAGGG for plant)")
    parser.add_argument("-o",type=str,metavar="str",help="output file name of bed format. ")
    args = parser.parse_args()

    EDGE_LIMIT=args.e
    INTERVAL_LIMIT=args.a
    CONTIG_LENGTH_LIMIT=args.l
    TELOMERE=args.t.split(",")

    fhr=open(args.i,"r")
    seqs={}
    ID=""
    for line in fhr:
        line=line.rstrip()
        if line[0]== ">":
            ID = line
            seqs[ID]=[]
        else:
            seqs[ID].append(line)
    fhr.close()
    for id in seqs:
        seqs[id]="".join(seqs[id])
        
    outs=""
    for ID in seqs:
        if len(seqs[ID]) < CONTIG_LENGTH_LIMIT:
            #print(ID[1:]+"\tcontig length\t"+str(len(seqs[ID]))+"\tless than\t"+str(CONTIG_LENGTH_LIMIT))
            continue
        
        for telomere in TELOMERE:
            spos=seqs[ID].find(telomere)
            pos=spos
            if spos > EDGE_LIMIT or spos < 0:
                continue
            tmp=find_telomere(seqs[ID], telomere, INTERVAL_LIMIT)
            print(ID[1:], telomere, "F", tmp[0], tmp[1], tmp[2])
            outs+=ID[1:] +"\t"+ str(tmp[0]) +"\t"+ str(tmp[1]) +"\t"+ telomere +"\t"+ str(tmp[2]) +"\t+\n" 
            
            revseq=rev_com(seqs[ID])
            spos=revseq.find(telomere)
            pos=spos
            if spos > EDGE_LIMIT or spos < 0:
                continue

            tmp=find_telomere(revseq, telomere, INTERVAL_LIMIT)
            print(ID[1:], rev_com(telomere), "R", tmp[0], tmp[1], tmp[2])
            outs+=ID[1:] +"\t"+ str(len(seqs[ID])-tmp[1]) +"\t"+ str(len(seqs[ID]) - tmp[0]) +"\t"+ telomere +"\t"+ str(tmp[2]) +"\t-\n"

    if args.o:
        fhw=open(args.o,"w")
        fhw.write(outs)
        fhw.close()
if __name__ == '__main__': main()

