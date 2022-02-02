#!/usr/bin/env python3

# Convert UBNT antenna pattern files into consolidated .ant format
# UBNT patterns are in 2 files, have 2 different planes of resolution...

import sys 
import os
import csv

def rotate(l, n):
    return l[-n:] + l[:-n]

if len(sys.argv) < 3:
	print("Usage: ubnt2ant.py azimuth.csv elevation.csv")
	quit()
# Build output filename from input filename(s)
out = sys.argv[1].split(".")[0]+".ant"

fout = open(out,"w")
az = []
el = []

# Read rows, take every other
csvfile = csv.reader(open(sys.argv[1]),delimiter=',')
n=0
for row in csvfile:
	if n % 2:
		print(row[1])
		az.append(row[1])
	n=n+1

az = rotate(az, 180)

for v in az:
	fout.write(v+"\n")

# Read rows, take every row
csvfile = csv.reader(open(sys.argv[2]),delimiter=',')
n=0
for row in csvfile:
	if n < 360:
		print(row[1])
		el.append(row[1])
	n=n+1

el = rotate(el, -90)

for v in el:
	fout.write(v+"\n")

fout.close()
print ("Pattern written to %s" % out)
