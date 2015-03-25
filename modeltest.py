from math import sin, cos, asin, sqrt, degrees, radians
import os, sys

Earth_radius_km = 6371.0
RADIUS = Earth_radius_km

def dist2degs(lat, lon, km):
    dlat = km / RADIUS
    dlon = asin(sin(dlat) / cos(radians(lat)))
    return degrees(dlon)


lat=50
lon=10

pm=sys.argv[1]
pe=sys.argv[2]
erp=sys.argv[3]
f=int(sys.argv[4])
rad=int(sys.argv[5])

for km in range(1,31):
	rlo = lon+dist2degs(lat,lon,rad);
	cmd = "./signalserver -m -lat "+str(lat)+" -lon "+str(lon)+" -rla "+str(lat)+" -rlo "+str(rlo)+" -txh 30 -rxh 2 -f "+str(f)+" -pm "+pm+" -pe "+pe+" -res 1200"
	out = os.popen(cmd).read().split("\n")
	#print out
	try:
	   	db = out[1]
		dbm = out[2]
		dbuv = out[3]
	except:
		db = 0
		dbm = 0
		dbuv = 0
	#print str(rad)+"km, "+str(f)+"MHz = "+str(db)
	print str(db)+",",
	f=f+200
