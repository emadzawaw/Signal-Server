#!/bin/sh
################################################################################
# utils/sdf/convert_sdf.sh
#
# Converts from *.hgt files to *.sdf files
#
# SD quality
#
################################################################################

for file in *.hgt
do
	srtm2sdf -d /dev/null $file
done
