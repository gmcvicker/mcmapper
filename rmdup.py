from random import choice
import os
import sys
import optparse


parser = optparse.OptionParser()
args= parser.parse_args()

print args[1][0]
print args[1][1]
infile = open(args[1][0],'r')
outfile = open(args[1][1],'w')

line=infile.readline()
chr=""
pos=0
linelistplus=[]
linelistminus=[]
while line:
    read = line.strip().split()
    if read[1]!=chr or read[2]!=pos:
        if(len(linelistplus)>0):
            outfile.write(choice(linelistplus))
        if(len(linelistminus)>0):
            outfile.write(choice(linelistminus))
        chr=read[1]
        pos=read[2]
        linelistplus=[]
        linelistminus=[]
    if(read[3]=="+"):
        linelistplus.append(line)
    if(read[3]=="-"):
        linelistminus.append(line)
    line=infile.readline()
infile.close()
outfile.close()
