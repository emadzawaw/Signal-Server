#!/bin/bash
set -e

if [[ $# -eq 0 ]];
then
	echo "Usage:";
	echo "$0 -lat 46.4879 -lon -123.2144 -txh 60 -f 145 -erp 30 -R 80 -o outfile | ./genkmz.sh";
	echo "The last argument must be the output file name.";
else
    # detect file name
	for name; do true; done;

	# This script should contain options common to every run on the target system.

	#time signalserverHD -sdf /mnt/data -pm 1 -rxh 6 -te 3 -cl 5 -pe 2 -dbg $@ 2>&1;
	time signalserver -sdf /mnt/data -pm 1 -rxh 6 -te 3 -cl 5 -pe 2 -dbg $@ 2>&1;
  
	# to resize, add: -resize 7000x7000\>
	convert $name.ppm -transparent white $name.png;
fi
