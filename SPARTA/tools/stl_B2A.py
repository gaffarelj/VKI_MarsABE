# This script has been slightly modified from the one here:
# https://github.com/IsseiMori/binary-stl-toASCII

def error(str=None):
  if str: print "ERROR:",str
  else: print "Syntax: stl_B2A.py stlfile"
  sys.exit()

# ----------------------------------------------------------------------
# main program

import sys

if len(sys.argv) not in [2, 3] : error()

bin_name = sys.argv[1]

f = 1
if len(sys.argv) == 3 and sys.argv[2] == "-rs": # rescale from mm to m
  f = 1000

# -*- encoding: utf-8 -*-
import struct

infile = open(bin_name) #import file
out_name = bin_name.split(".")
out_name = out_name[0] + "_ASCII" + ".stl"
out = open(out_name, 'w') #export file

data = infile.read()

out.write("solid Body 1\n")

number = data[80] + data[81] + data[82] + data[83] 
faces = struct.unpack('I',number)[0]

for x in range(0,faces):
    out.write("facet normal ")

    xc = data[84+x*50] + data[85+x*50] + data[86+x*50] + data[87+x*50]
    yc = data[88+x*50] + data[89+x*50] + data[90+x*50] + data[91+x*50]
    zc = data[92+x*50] + data[93+x*50] + data[94+x*50] + data[95+x*50]

    out.write(str(struct.unpack('f',xc)[0]) + " ")
    out.write(str(struct.unpack('f',yc)[0]) + " ")
    out.write(str(struct.unpack('f',zc)[0]) + "\n")

    out.write("outer loop\n")

    for y in range(1,4):
        out.write("vertex ")

        xc = data[84+y*12+x*50] + data[85+y*12+x*50] + data[86+y*12+x*50] + data[87+y*12+x*50]
        yc = data[88+y*12+x*50] + data[89+y*12+x*50] + data[90+y*12+x*50] + data[91+y*12+x*50]
        zc = data[92+y*12+x*50] + data[93+y*12+x*50] + data[94+y*12+x*50] + data[95+y*12+x*50]

        out.write(str(struct.unpack('f',xc)[0]/f) + " ")
        out.write(str(struct.unpack('f',yc)[0]/f) + " ")
        out.write(str(struct.unpack('f',zc)[0]/f) + "\n")

    out.write("endloop\n")
    out.write("endfacet\n")

out.write("endsolid Body 1\n")

out.close()
