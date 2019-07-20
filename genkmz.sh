#!/bin/bash
while read line
do
	echo $line
	if [[ $line == Writing* ]]
	then
		while IFS='"' read -ra writingline
		do
#		        echo "writingline: ${writingline}"
			filename=${writingline[1]%.*}
#			echo "filename: ${filename}"
		done <<< $line
	fi
	if [[ $line == \|* ]]
	then
		while IFS='|' read -ra coords
		do
			cat << EOF > doc.kml
<?xml version="1.0" encoding="UTF-8"?>
<kml xmlns="http://www.opengis.net/kml/2.2" xmlns:gx="http://www.google.com/kml/ext/2.2" xmlns:kml="http://www.opengis.net/kml/2.2" xmlns:atom="http://www.w3.org/2005/Atom">
<GroundOverlay>
    <name>${filename}</name>
    <Icon>
        <href>${filename}.png</href>
    </Icon>
    <LatLonBox>
        <north>${coords[1]}</north>
        <east>${coords[2]}</east>
        <south>${coords[3]}</south>
        <west>${coords[4]}</west>
    </LatLonBox>
</GroundOverlay>
</kml>
EOF
		done <<< $line
	fi
done
#echo "filename: ${filename}"
zip "${filename}" "${filename}.png" "doc.kml"
mv "${filename}".zip "${filename}".kmz
rm "doc.kml"
echo Generated "${filename}".kmz
