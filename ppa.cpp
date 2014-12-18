/****************************************************************************\
*	   PPA: Path Profile Analysis Tool v1.3.1 derived from SPLAT! 1.3           *
******************************************************************************
*	     Project started in 1997 by John A. Magliacane, KD2BD 	     *
*			  Last update: 10-Apr-2009			     *
******************************************************************************
*         Please consult the documentation for a complete list of	     *
*	     individuals who have contributed to this project. 		     *
******************************************************************************
*                                                                            *
*  This program is free software; you can redistribute it and/or modify it   *
*  under the terms of the GNU General Public License as published by the     *
*  Free Software Foundation; either version 2 of the License or any later    *
*  version.								     *
* 									     *
*  This program is distributed in the hope that it will useful, but WITHOUT  *
*  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or     *
*  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License     *
*  for more details.							     *
*									     *
******************************************************************************
* g++ -Wall -O3 -s -lm -fomit-frame-pointer itm.cpp ppa.cpp -o ppa     * 
\****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>

#define GAMMA 2.5
#define BZBUFFER 65536
#define HD_MODE 0
#define MAXPAGES 64
#define ARRAYSIZE 76810
#define IPPD 1200


#ifndef PI
#define PI 3.141592653589793
#endif

#ifndef TWOPI
#define TWOPI 6.283185307179586
#endif

#ifndef HALFPI
#define HALFPI 1.570796326794896
#endif

#define DEG2RAD 1.74532925199e-02
#define EARTHRADIUS 20902230.97
#define	METERS_PER_MILE 1609.344
#define METERS_PER_FOOT 0.3048
#define	KM_PER_MILE 1.609344
#define FOUR_THIRDS 1.3333333333333

char 	string[255], sdf_path[255], opened=0, gpsav=0, ppa_name[10],
	ppa_version[6], dashes[80];

double	earthradius, max_range=0.0, forced_erp=-1.0, dpp, ppd,
	fzone_clearance=0.6, clutter, loss, dBm, field_strength;

int	min_north=90, max_north=-90, min_west=360, max_west=-1, ippd, mpi,
	max_elevation=-32768, min_elevation=32768, bzerror, contour_threshold, output=2;

unsigned char got_elevation_pattern, got_azimuth_pattern, metric=0, dbm=0;



struct site {	double lat;
		double lon;
		float alt;
		char name[50];
		char filename[255];
	    } 	site;

struct path {	double lat[ARRAYSIZE];
		double lon[ARRAYSIZE];
		double elevation[ARRAYSIZE];
		double distance[ARRAYSIZE];
		int length;
	    }	path;

struct dem {	int min_north;
		int max_north;
		int min_west;
		int max_west;
		int max_el;
		int min_el;
		short data[IPPD][IPPD];
		unsigned char mask[IPPD][IPPD];
		unsigned char signal[IPPD][IPPD];
           }	dem[MAXPAGES];

struct LR {	double eps_dielect; 
		double sgm_conductivity; 
		double eno_ns_surfref;
		double frq_mhz; 
		double conf; 
		double rel;
		double erp;
		int radio_climate;  
		int pol;
		float antenna_pattern[361][1001];
          }	LR;

struct region { unsigned char color[32][3];
		int level[32];
		int levels;
	      }	region;

double elev[ARRAYSIZE+10];

void point_to_point(double elev[], double tht_m, double rht_m,
	  double eps_dielect, double sgm_conductivity, double eno_ns_surfref,
	  double frq_mhz, int radio_climate, int pol, double conf,
	  double rel, double &dbloss, char *strmode, int &errnum);

double arccos(double x, double y)
{
	/* This function implements the arc cosine function,
	   returning a value between 0 and TWOPI. */

	double result=0.0;

	if (y>0.0)
		result=acos(x/y);

	if (y<0.0)
		result=PI+acos(x/y);

	return result;
}

int ReduceAngle(double angle)
{
	/* This function normalizes the argument to
	   an integer angle between 0 and 180 degrees */

	double temp;

	temp=acos(cos(angle*DEG2RAD));

	return (int)rint(temp/DEG2RAD);
}

double LonDiff(double lon1, double lon2)
{
	/* This function returns the short path longitudinal
	   difference between longitude1 and longitude2 
	   as an angle between -180.0 and +180.0 degrees.
	   If lon1 is west of lon2, the result is positive.
	   If lon1 is east of lon2, the result is negative. */

	double diff;

	diff=lon1-lon2;

	if (diff<=-180.0)
		diff+=360.0;

	if (diff>=180.0)
		diff-=360.0;

	return diff;
}

char *dec2dms(double decimal)
{
	/* Converts decimal degrees to degrees, minutes, seconds,
	   (DMS) and returns the result as a character string. */

	char	sign;
	int	degrees, minutes, seconds;
	double	a, b, c, d;

	if (decimal<0.0)
	{
		decimal=-decimal;
		sign=-1;
	}

	else
		sign=1;

	a=floor(decimal);
	b=60.0*(decimal-a);
	c=floor(b);
	d=60.0*(b-c);

	degrees=(int)a;
	minutes=(int)c;
	seconds=(int)d;

	if (seconds<0)
		seconds=0;

	if (seconds>59)
		seconds=59;

	string[0]=0;
	snprintf(string,250,"%d%c %d\' %d\"", degrees*sign, 176, minutes, seconds);
	return (string);
}



int OrMask(double lat, double lon, int value)
{
	
	int	x, y, indx;
	char	found;

	for (indx=0, found=0; indx<MAXPAGES && found==0;)
	{
		x=(int)rint(ppd*(lat-dem[indx].min_north));
		y=mpi-(int)rint(ppd*(LonDiff(dem[indx].max_west,lon)));

		if (x>=0 && x<=mpi && y>=0 && y<=mpi)
			found=1;
		else
			indx++;
	}

	if (found)
	{
		dem[indx].mask[x][y]|=value;
		return ((int)dem[indx].mask[x][y]);
	}

	else
		return -1;
}


double GetElevation(struct site location)
{
	/* This function returns the elevation (in feet) of any location
	   represented by the digital elevation model data in memory.
	   Function returns -5000.0 for locations not found in memory. */

	char	found;
	int	x, y, indx;
	double	elevation;

	for (indx=0, found=0; indx<MAXPAGES && found==0;)
	{
		x=(int)rint(ppd*(location.lat-dem[indx].min_north));
		y=mpi-(int)rint(ppd*(LonDiff(dem[indx].max_west,location.lon)));

		if (x>=0 && x<=mpi && y>=0 && y<=mpi)
			found=1;
		else
			indx++;
	}

	if (found)
		elevation=3.28084*dem[indx].data[x][y];
	else
		elevation=-5000.0;
	
			
	return elevation;
}

int AddElevation(double lat, double lon, double height)
{
	/* This function adds a user-defined terrain feature
	   (in meters AGL) to the digital elevation model data
	   in memory.  Does nothing and returns 0 for locations
	   not found in memory. */

	char	found;
	int	x, y, indx;

	for (indx=0, found=0; indx<MAXPAGES && found==0;)
	{
		x=(int)rint(ppd*(lat-dem[indx].min_north));
		y=mpi-(int)rint(ppd*(LonDiff(dem[indx].max_west,lon)));

		if (x>=0 && x<=mpi && y>=0 && y<=mpi)
			found=1;
		else
			indx++;
	}

	if (found)
		dem[indx].data[x][y]+=(short)rint(height);

	return found;
}

double Distance(struct site site1, struct site site2)
{
	/* This function returns the great circle distance
	   in miles between any two site locations. */

	double	lat1, lon1, lat2, lon2, distance;

	lat1=site1.lat*DEG2RAD;
	lon1=site1.lon*DEG2RAD;
	lat2=site2.lat*DEG2RAD;
	lon2=site2.lon*DEG2RAD;

	distance=3959.0*acos(sin(lat1)*sin(lat2)+cos(lat1)*cos(lat2)*cos((lon1)-(lon2)));

	return distance;
}

double Azimuth(struct site source, struct site destination)
{
	/* This function returns the azimuth (in degrees) to the
	   destination as seen from the location of the source. */

	double	dest_lat, dest_lon, src_lat, src_lon,
		beta, azimuth, diff, num, den, fraction;

	dest_lat=destination.lat*DEG2RAD;
	dest_lon=destination.lon*DEG2RAD;

	src_lat=source.lat*DEG2RAD;
	src_lon=source.lon*DEG2RAD;
		
	/* Calculate Surface Distance */

	beta=acos(sin(src_lat)*sin(dest_lat)+cos(src_lat)*cos(dest_lat)*cos(src_lon-dest_lon));

	/* Calculate Azimuth */

	num=sin(dest_lat)-(sin(src_lat)*cos(beta));
	den=cos(src_lat)*sin(beta);
	fraction=num/den;

	/* Trap potential problems in acos() due to rounding */

	if (fraction>=1.0)
		fraction=1.0;

	if (fraction<=-1.0)
		fraction=-1.0;

	/* Calculate azimuth */

	azimuth=acos(fraction);

	/* Reference it to True North */

	diff=dest_lon-src_lon;

	if (diff<=-PI)
		diff+=TWOPI;

	if (diff>=PI)
		diff-=TWOPI;

	if (diff>0.0)
		azimuth=TWOPI-azimuth;

	return (azimuth/DEG2RAD);		
}

double ElevationAngle(struct site source, struct site destination)
{
	/* This function returns the angle of elevation (in degrees)
	   of the destination as seen from the source location.
	   A positive result represents an angle of elevation (uptilt),
	   while a negative result represents an angle of depression
	   (downtilt), as referenced to a normal to the center of
	   the earth. */
	   
	register double a, b, dx;

	a=GetElevation(destination)+destination.alt+earthradius;
	b=GetElevation(source)+source.alt+earthradius;

 	dx=5280.0*Distance(source,destination);

	/* Apply the Law of Cosines */

	return ((180.0*(acos(((b*b)+(dx*dx)-(a*a))/(2.0*b*dx)))/PI)-90.0);
}

void ReadPath(struct site source, struct site destination)
{
	/* This function generates a sequence of latitude and
	   longitude positions between source and destination
	   locations along a great circle path, and stores
	   elevation and distance information for points
	   along that path in the "path" structure. */

	int	c;
	double	azimuth, distance, lat1, lon1, beta, den, num,
		lat2, lon2, total_distance, dx, dy, path_length,
		miles_per_sample, samples_per_radian=68755.0;
	struct	site tempsite;

	lat1=source.lat*DEG2RAD;
	lon1=source.lon*DEG2RAD;

	lat2=destination.lat*DEG2RAD;
	lon2=destination.lon*DEG2RAD;

	if (ppd==1200.0)
		samples_per_radian=68755.0;

	if (ppd==3600.0)
		samples_per_radian=206265.0;

	azimuth=Azimuth(source,destination)*DEG2RAD;

	total_distance=Distance(source,destination);

	if (total_distance>(30.0/ppd))		/* > 0.5 pixel distance */
	{
		dx=samples_per_radian*acos(cos(lon1-lon2));
		dy=samples_per_radian*acos(cos(lat1-lat2));

		path_length=sqrt((dx*dx)+(dy*dy));		/* Total number of samples */

		miles_per_sample=total_distance/path_length;	/* Miles per sample */
	}

	else
	{
		c=0;
		dx=0.0;
		dy=0.0;
		path_length=0.0;
		miles_per_sample=0.0;
		total_distance=0.0;

		lat1=lat1/DEG2RAD;
		lon1=lon1/DEG2RAD;

		path.lat[c]=lat1;
		path.lon[c]=lon1;
		path.elevation[c]=GetElevation(source);
		path.distance[c]=0.0;
	}

	for (distance=0.0, c=0; (total_distance!=0.0 && distance<=total_distance && c<ARRAYSIZE); c++, distance=miles_per_sample*(double)c)
	{
		beta=distance/3959.0;
		lat2=asin(sin(lat1)*cos(beta)+cos(azimuth)*sin(beta)*cos(lat1));
		num=cos(beta)-(sin(lat1)*sin(lat2));
		den=cos(lat1)*cos(lat2);

		if (azimuth==0.0 && (beta>HALFPI-lat1))
			lon2=lon1+PI;

		else if (azimuth==HALFPI && (beta>HALFPI+lat1))
				lon2=lon1+PI;

		else if (fabs(num/den)>1.0)
				lon2=lon1;

		else
		{
			if ((PI-azimuth)>=0.0)
				lon2=lon1-arccos(num,den);
			else
				lon2=lon1+arccos(num,den);
		}
	
		while (lon2<0.0)
			lon2+=TWOPI;

		while (lon2>TWOPI)
			lon2-=TWOPI;
 
		lat2=lat2/DEG2RAD;
		lon2=lon2/DEG2RAD;

		path.lat[c]=lat2;
		path.lon[c]=lon2;
		tempsite.lat=lat2;
		tempsite.lon=lon2;
		path.elevation[c]=GetElevation(tempsite);
		path.distance[c]=distance;
	}

	/* Make sure exact destination point is recorded at path.length-1 */

	if (c<ARRAYSIZE)
	{
		path.lat[c]=destination.lat;
		path.lon[c]=destination.lon;
		path.elevation[c]=GetElevation(destination);
		path.distance[c]=total_distance;
		c++;
	}

	if (c<ARRAYSIZE)
		path.length=c;
	else
		path.length=ARRAYSIZE-1;
}

double ElevationAngle2(struct site source, struct site destination, double er)
{
	/* This function returns the angle of elevation (in degrees)
	   of the destination as seen from the source location, UNLESS
	   the path between the sites is obstructed, in which case, the
	   elevation angle to the first obstruction is returned instead.
	   "er" represents the earth radius. */

	int	x;
	char	block=0;
	double	source_alt, destination_alt, cos_xmtr_angle,
		cos_test_angle, test_alt, elevation, distance,
		source_alt2, first_obstruction_angle=0.0;
	struct	path temp;

	temp=path;

	ReadPath(source,destination);

	distance=5280.0*Distance(source,destination);
	source_alt=er+source.alt+GetElevation(source);
	destination_alt=er+destination.alt+GetElevation(destination);
	source_alt2=source_alt*source_alt;

	/* Calculate the cosine of the elevation angle of the
	   destination (receiver) as seen by the source (transmitter). */

	cos_xmtr_angle=((source_alt2)+(distance*distance)-(destination_alt*destination_alt))/(2.0*source_alt*distance);

	/* Test all points in between source and destination locations to
	   see if the angle to a topographic feature generates a higher
	   elevation angle than that produced by the destination.  Begin
	   at the source since we're interested in identifying the FIRST
	   obstruction along the path between source and destination. */
 
	for (x=2, block=0; x<path.length && block==0; x++)
	{
		distance=5280.0*path.distance[x];

		test_alt=earthradius+(path.elevation[x]==0.0?path.elevation[x]:path.elevation[x]+clutter);

		cos_test_angle=((source_alt2)+(distance*distance)-(test_alt*test_alt))/(2.0*source_alt*distance);

		/* Compare these two angles to determine if
		   an obstruction exists.  Since we're comparing
		   the cosines of these angles rather than
		   the angles themselves, the sense of the
		   following "if" statement is reversed from
		   what it would be if the angles themselves
		   were compared. */

		if (cos_xmtr_angle>=cos_test_angle)
		{
			block=1;
			first_obstruction_angle=((acos(cos_test_angle))/DEG2RAD)-90.0;
		}
	}

	if (block)
		elevation=first_obstruction_angle;

	else
		elevation=((acos(cos_xmtr_angle))/DEG2RAD)-90.0;

	path=temp;

	return elevation;
}

double AverageTerrain(struct site source, double azimuthx, double start_distance, double end_distance)
{
	/* This function returns the average terrain calculated in
	   the direction of "azimuth" (degrees) between "start_distance"
	   and "end_distance" (miles) from the source location.  If
	   the terrain is all water (non-critical error), -5000.0 is
	   returned.  If not enough SDF data has been loaded into
	   memory to complete the survey (critical error), then
	   -9999.0 is returned. */
 
	int	c, samples, endpoint;
	double	beta, lat1, lon1, lat2, lon2, num, den, azimuth, terrain=0.0;
	struct	site destination;

	lat1=source.lat*DEG2RAD;
	lon1=source.lon*DEG2RAD;

	/* Generate a path of elevations between the source
	   location and the remote location provided. */

	beta=end_distance/3959.0;

	azimuth=DEG2RAD*azimuthx;

	lat2=asin(sin(lat1)*cos(beta)+cos(azimuth)*sin(beta)*cos(lat1));
	num=cos(beta)-(sin(lat1)*sin(lat2));
	den=cos(lat1)*cos(lat2);

	if (azimuth==0.0 && (beta>HALFPI-lat1))
		lon2=lon1+PI;

	else if (azimuth==HALFPI && (beta>HALFPI+lat1))
			lon2=lon1+PI;

	else if (fabs(num/den)>1.0)
			lon2=lon1;

	else
	{
		if ((PI-azimuth)>=0.0)
			lon2=lon1-arccos(num,den);
		else
			lon2=lon1+arccos(num,den);
	}
	
	while (lon2<0.0)
		lon2+=TWOPI;

	while (lon2>TWOPI)
		lon2-=TWOPI;
 
	lat2=lat2/DEG2RAD;
	lon2=lon2/DEG2RAD;

	destination.lat=lat2;
	destination.lon=lon2;

	/* If SDF data is missing for the endpoint of
	   the radial, then the average terrain cannot
	   be accurately calculated.  Return -9999.0 */

	if (GetElevation(destination)<-4999.0)
		return (-9999.0);
	else
	{
		ReadPath(source,destination);

		endpoint=path.length;

		/* Shrink the length of the radial if the
		   outermost portion is not over U.S. land. */

		for (c=endpoint-1; c>=0 && path.elevation[c]==0.0; c--);

		endpoint=c+1;

		for (c=0, samples=0; c<endpoint; c++)
		{
			if (path.distance[c]>=start_distance)
			{
				terrain+=(path.elevation[c]==0.0?path.elevation[c]:path.elevation[c]+clutter);
				samples++;
			}
		}

		if (samples==0)
			terrain=-5000.0;  /* No land */
		else
			terrain=(terrain/(double)samples);

		return terrain;
	}
}

double haat(struct site antenna)
{
	/* This function returns the antenna's Height Above Average
	   Terrain (HAAT) based on FCC Part 73.313(d).  If a critical
	   error occurs, such as a lack of SDF data to complete the
	   survey, -5000.0 is returned. */

	int	azi, c;
	char	error=0;
	double	terrain, avg_terrain, haat, sum=0.0;

	/* Calculate the average terrain between 2 and 10 miles
	   from the antenna site at azimuths of 0, 45, 90, 135,
	   180, 225, 270, and 315 degrees. */

	for (c=0, azi=0; azi<=315 && error==0; azi+=45)
	{
		terrain=AverageTerrain(antenna, (double)azi, 2.0, 10.0);

		if (terrain<-9998.0)  /* SDF data is missing */
			error=1;

		if (terrain>-4999.0)  /* It's land, not water */
		{
			sum+=terrain;  /* Sum of averages */
			c++;
		}
	}

	if (error)
		return -5000.0;
	else
	{
		avg_terrain=(sum/(double)c);
		haat=(antenna.alt+GetElevation(antenna))-avg_terrain;
		return haat;
	}
}



double ReadBearing(char *input)
{
	/* This function takes numeric input in the form of a character
	   string, and returns an equivalent bearing in degrees as a
	   decimal number (double).  The input may either be expressed
	   in decimal format (40.139722) or degree, minute, second
	   format (40 08 23).  This function also safely handles
	   extra spaces found either leading, trailing, or
	   embedded within the numbers expressed in the
	   input string.  Decimal seconds are permitted. */
 
	double	seconds, bearing=0.0;
	char	string[20];
	int	a, b, length, degrees, minutes;

	/* Copy "input" to "string", and ignore any extra
	   spaces that might be present in the process. */

	string[0]=0;
	length=strlen(input);

	for (a=0, b=0; a<length && a<18; a++)
	{
		if ((input[a]!=32 && input[a]!='\n') || (input[a]==32 && input[a+1]!=32 && input[a+1]!='\n' && b!=0))
		{
			string[b]=input[a];
			b++;
		}	 
	}

	string[b]=0;

	/* Count number of spaces in the clean string. */

	length=strlen(string);

	for (a=0, b=0; a<length; a++)
		if (string[a]==32)
			b++;

	if (b==0)  /* Decimal Format (40.139722) */
		sscanf(string,"%lf",&bearing);

	if (b==2)  /* Degree, Minute, Second Format (40 08 23.xx) */
	{
		sscanf(string,"%d %d %lf",&degrees, &minutes, &seconds);

		bearing=fabs((double)degrees);
		bearing+=fabs(((double)minutes)/60.0);
		bearing+=fabs(seconds/3600.0);

		if ((degrees<0) || (minutes<0) || (seconds<0.0))
			bearing=-bearing;
	}

	/* Anything else returns a 0.0 */

	if (bearing>360.0 || bearing<-360.0)
		bearing=0.0;

	return bearing;
}

int LoadSDF_SDF(char *name,int winfiles)
{
	
	int	x, y, data, indx, minlat, minlon, maxlat, maxlon;
	char	found, free_page=0, line[20], sdf_file[255],
		path_plus_name[255], *s=NULL;
	FILE	*fd;

	for (x=0; name[x]!='.' && name[x]!=0 && x<250; x++)
		sdf_file[x]=name[x];

	sdf_file[x]=0;

	/* Parse filename for minimum latitude and longitude values */
	if(winfiles==1){
	sscanf(sdf_file,"%d=%d=%d=%d",&minlat,&maxlat,&minlon,&maxlon);
	}else{
	sscanf(sdf_file,"%d:%d:%d:%d",&minlat,&maxlat,&minlon,&maxlon);
	}
	sdf_file[x]='.';
	sdf_file[x+1]='s';
	sdf_file[x+2]='d';
	sdf_file[x+3]='f';
	sdf_file[x+4]=0;

	/* Is it already in memory? */

	for (indx=0, found=0; indx<MAXPAGES && found==0; indx++)
	{
		if (minlat==dem[indx].min_north && minlon==dem[indx].min_west && maxlat==dem[indx].max_north && maxlon==dem[indx].max_west)
			found=1;
	}

	/* Is room available to load it? */

	if (found==0)
	{	
		for (indx=0, free_page=0; indx<MAXPAGES && free_page==0; indx++)
			if (dem[indx].max_north==-90)
				free_page=1;
	}

	indx--;

	if (free_page && found==0 && indx>=0 && indx<MAXPAGES)
	{
		/* Search for SDF file in current working directory first */

		strncpy(path_plus_name,sdf_file,255);

		fd=fopen(path_plus_name,"rb");

		if (fd==NULL)
		{
			
			strncpy(path_plus_name,sdf_path,255);
			strncat(path_plus_name,sdf_file,255);

			fd=fopen(path_plus_name,"rb");
		}

		if (fd!=NULL)
		{
		

			s=fgets(line,19,fd);
			sscanf(line,"%d",&dem[indx].max_west);

			s=fgets(line,19,fd);
			sscanf(line,"%d",&dem[indx].min_north);

			s=fgets(line,19,fd);
			sscanf(line,"%d",&dem[indx].min_west);

			s=fgets(line,19,fd);
			sscanf(line,"%d",&dem[indx].max_north);

			for (x=0; x<ippd; x++)
				for (y=0; y<ippd; y++)
				{
					s=fgets(line,19,fd);
					data=atoi(line);

					dem[indx].data[x][y]=data;
					dem[indx].signal[x][y]=0;
					dem[indx].mask[x][y]=0;

					if (data>dem[indx].max_el)
						dem[indx].max_el=data;

					if (data<dem[indx].min_el)
						dem[indx].min_el=data;
				}

			fclose(fd);

			if (dem[indx].min_el<min_elevation)
				min_elevation=dem[indx].min_el;

			if (dem[indx].max_el>max_elevation)
				max_elevation=dem[indx].max_el;

			if (max_north==-90)
				max_north=dem[indx].max_north;

			else if (dem[indx].max_north>max_north)
					max_north=dem[indx].max_north;

			if (min_north==90)
				min_north=dem[indx].min_north;

			else if (dem[indx].min_north<min_north)
					min_north=dem[indx].min_north;

			if (max_west==-1)
				max_west=dem[indx].max_west;

			else
			{
				if (abs(dem[indx].max_west-max_west)<180)
				{
 					if (dem[indx].max_west>max_west)
						max_west=dem[indx].max_west;
				}

				else
				{
 					if (dem[indx].max_west<max_west)
						max_west=dem[indx].max_west;
				}
			}

			if (min_west==360)
				min_west=dem[indx].min_west;

			else
			{
				if (fabs(dem[indx].min_west-min_west)<180.0)
				{
 					if (dem[indx].min_west<min_west)
						min_west=dem[indx].min_west;
				}

				else
				{
 					if (dem[indx].min_west>min_west)
						min_west=dem[indx].min_west;
				}
			}

			

			return 1;
		}

		else
			return -1;
	}

	else
		return 0;
}

char LoadSDF(char *name, int winfiles)
{
	/* This function loads the requested SDF file from the filesystem.
	   It first tries to invoke the LoadSDF_SDF() function to load an
	   uncompressed SDF file (since uncompressed files load slightly
	   faster).  If that attempt fails, then it tries to load a
	   compressed SDF file by invoking the LoadSDF_BZ() function.
	   If that fails, then we can assume that no elevation data
	   exists for the region requested, and that the region
	   requested must be entirely over water. */

	int	x, y, indx, minlat, minlon, maxlat, maxlon;
	char	found, free_page=0;
	int	return_value=-1;

	/* Try to load an uncompressed SDF first. */

	return_value=LoadSDF_SDF(name,winfiles);


	/* If neither format can be found, then assume the area is water. */

	if (return_value==0 || return_value==-1)
	{
	if(winfiles==1){
		sscanf(name,"%d=%d=%d=%d",&minlat,&maxlat,&minlon,&maxlon);
	}else{
		sscanf(name,"%d:%d:%d:%d",&minlat,&maxlat,&minlon,&maxlon);
	}
		/* Is it already in memory? */

		for (indx=0, found=0; indx<MAXPAGES && found==0; indx++)
		{
			if (minlat==dem[indx].min_north && minlon==dem[indx].min_west && maxlat==dem[indx].max_north && maxlon==dem[indx].max_west)
				found=1;
		}

		/* Is room available to load it? */

		if (found==0)
		{	
			for (indx=0, free_page=0; indx<MAXPAGES && free_page==0; indx++)
				if (dem[indx].max_north==-90)
					free_page=1;
		}

		indx--;

		if (free_page && found==0 && indx>=0 && indx<MAXPAGES)
		{
			fprintf(stdout,"Region  \"%s\" assumed as sea-level into page %d...",name,indx+1);
			fflush(stdout);

			dem[indx].max_west=maxlon;
			dem[indx].min_north=minlat;
			dem[indx].min_west=minlon;
			dem[indx].max_north=maxlat;

			/* Fill DEM with sea-level topography */

			for (x=0; x<ippd; x++)
				for (y=0; y<ippd; y++)
				{
		    			dem[indx].data[x][y]=0;
					dem[indx].signal[x][y]=0;
					dem[indx].mask[x][y]=0;

					if (dem[indx].min_el>0)
						dem[indx].min_el=0;
				}

			if (dem[indx].min_el<min_elevation)
				min_elevation=dem[indx].min_el;

			if (dem[indx].max_el>max_elevation)
				max_elevation=dem[indx].max_el;

			if (max_north==-90)
				max_north=dem[indx].max_north;

			else if (dem[indx].max_north>max_north)
					max_north=dem[indx].max_north;

			if (min_north==90)
				min_north=dem[indx].min_north;

			else if (dem[indx].min_north<min_north)
					min_north=dem[indx].min_north;

			if (max_west==-1)
				max_west=dem[indx].max_west;

			else
			{
				if (abs(dem[indx].max_west-max_west)<180)
				{
 					if (dem[indx].max_west>max_west)
						max_west=dem[indx].max_west;
				}

				else
				{
 					if (dem[indx].max_west<max_west)
						max_west=dem[indx].max_west;
				}
			}

			if (min_west==360)
				min_west=dem[indx].min_west;

			else
			{
				if (abs(dem[indx].min_west-min_west)<180)
				{
 					if (dem[indx].min_west<min_west)
						min_west=dem[indx].min_west;
				}

				else
				{
 					if (dem[indx].min_west>min_west)
						min_west=dem[indx].min_west;
				}
			}

			fprintf(stdout," Done!\n");
			fflush(stdout);

			return_value=1;
		}
	}

	return return_value;
}

int GetMask(double lat, double lon)
{
	/* This function returns the mask bits based on the latitude
	   and longitude given. */

	return (OrMask(lat,lon,0));
}


void PlotPath(struct site source, struct site destination, char mask_value)
{
	/* This function analyzes the path between the source and
	   destination locations.  It determines which points along
	   the path have line-of-sight visibility to the source.
	   Points along with path having line-of-sight visibility
	   to the source at an AGL altitude equal to that of the
	   destination location are stored by setting bit 1 in the
	   mask[][] array, which are displayed in green when PPM
	   maps are later generated by SPLAT!. */

	char block;
	int x, y;
	register double cos_xmtr_angle, cos_test_angle, test_alt;
	double distance, rx_alt, tx_alt;

	ReadPath(source,destination);

	for (y=0; y<path.length; y++)
	{
		/* Test this point only if it hasn't been already
		   tested and found to be free of obstructions. */

		if ((GetMask(path.lat[y],path.lon[y])&mask_value)==0)
		{
			distance=5280.0*path.distance[y];
			tx_alt=earthradius+source.alt+path.elevation[0];
			rx_alt=earthradius+destination.alt+path.elevation[y];

			/* Calculate the cosine of the elevation of the
			   transmitter as seen at the temp rx point. */

			cos_xmtr_angle=((rx_alt*rx_alt)+(distance*distance)-(tx_alt*tx_alt))/(2.0*rx_alt*distance);

			for (x=y, block=0; x>=0 && block==0; x--)
			{
				distance=5280.0*(path.distance[y]-path.distance[x]);
				test_alt=earthradius+(path.elevation[x]==0.0?path.elevation[x]:path.elevation[x]+clutter);

				cos_test_angle=((rx_alt*rx_alt)+(distance*distance)-(test_alt*test_alt))/(2.0*rx_alt*distance);

				/* Compare these two angles to determine if
				   an obstruction exists.  Since we're comparing
				   the cosines of these angles rather than
				   the angles themselves, the following "if"
				   statement is reversed from what it would
				   be if the actual angles were compared. */

				if (cos_xmtr_angle>=cos_test_angle)
					block=1;
			}

			if (block==0)
				OrMask(path.lat[y],path.lon[y],mask_value);
		}
	}
}




void GraphHeight(struct site source, struct site destination, char *name, unsigned char fresnel_plot, unsigned char normalized, int pngwidth, int pngheight)
{
	/* This function invokes gnuplot to generate an appropriate
	   output file indicating the terrain height profile between
	   the source and destination locations referenced to the
	   line-of-sight path between the receive and transmit sites
	   when the -h or -H command line option is used.  "basename"
	   is the name assigned to the output file generated by gnuplot.
	   The filename extension is used to set gnuplot's terminal
	   setting and output file type.  If no extension is found,
	   .png is assumed.  */

	int	x, y, z;
	char	basename[255], term[30], ext[15];
	double	a, b, c, height=0.0, refangle, cangle, maxheight=-100000.0,
		minheight=100000.0, lambda=0.0, f_zone=0.0, fpt6_zone=0.0,
		nm=0.0, nb=0.0, ed=0.0, es=0.0, r=0.0, d=0.0, d1=0.0,
		terrain, azimuth, distance, dheight=0.0, minterrain=100000.0,
		minearth=100000.0, miny, maxy, min2y, max2y;
	struct	site remote;
	FILE	*fd=NULL, *fd1=NULL, *fd2=NULL, *fd3=NULL, *fd4=NULL, *fd5=NULL;

	ReadPath(destination,source);  /* destination=RX, source=TX */
	azimuth=Azimuth(destination,source);
	distance=Distance(destination,source);
	refangle=ElevationAngle(destination,source);
	b=GetElevation(destination)+destination.alt+earthradius;

	/* Wavelength and path distance (great circle) in feet. */

	//LR.frq_mhz = freq;

	if (fresnel_plot)
	{
		lambda=9.8425e8/(LR.frq_mhz*1e6);
		d=5280.0*path.distance[path.length-1];
	}

	
	
	if (normalized)
	{
		ed=GetElevation(destination);
		es=GetElevation(source);
		nb=-destination.alt-ed;
		nm=(-source.alt-es-nb)/(path.distance[path.length-1]);
	}

	fd=fopen("profile.gp","wb");

	if (clutter>0.0)
		fd1=fopen("clutter.gp","wb");

	fd2=fopen("reference.gp","wb");
	fd5=fopen("curvature.gp", "wb");

	if ((LR.frq_mhz>=20.0) && (LR.frq_mhz<=100000.0) && fresnel_plot)
	{
		fd3=fopen("fresnel.gp", "wb");
		fd4=fopen("fresnel_pt_6.gp", "wb");
	}

			
	for (x=0; x<path.length-1; x++)
	{
		remote.lat=path.lat[x];
		remote.lon=path.lon[x];
		remote.alt=0.0;

		terrain=GetElevation(remote);

		if (x==0)
			terrain+=destination.alt;  /* RX antenna spike */

		a=terrain+earthradius;
 		cangle=5280.0*Distance(destination,remote)/earthradius;
		c=b*sin(refangle*DEG2RAD+HALFPI)/sin(HALFPI-refangle*DEG2RAD-cangle);

		height=a-c;

		/* Per Fink and Christiansen, Electronics
		 * Engineers' Handbook, 1989:
		 *
		 *   H = sqrt(lamba * d1 * (d - d1)/d)
		 *
		 * where H is the distance from the LOS
		 * path to the first Fresnel zone boundary.
		 */

		
		if ((LR.frq_mhz>=20.0) && (LR.frq_mhz<=100000.0) && fresnel_plot)
		{
			d1=5280.0*path.distance[x];
			f_zone=-1.0*sqrt(lambda*d1*(d-d1)/d);
			fpt6_zone=f_zone*fzone_clearance;
		}

		if (normalized)
		{
			r=-(nm*path.distance[x])-nb;
			height+=r;

			if ((LR.frq_mhz>=20.0) && (LR.frq_mhz<=100000.0) && fresnel_plot)
			{
				f_zone+=r;
				fpt6_zone+=r;
			}
		}

		else
			r=0.0;

		
	
		if (metric)
		{
		
			//segfault here
			fprintf(fd,"%f\t%f\n",KM_PER_MILE*path.distance[x],METERS_PER_FOOT*height);

			 
	
			if (fd1!=NULL && x>0 && x<path.length-2)
				fprintf(fd1,"%f\t%f\n",KM_PER_MILE*path.distance[x],METERS_PER_FOOT*(terrain==0.0?height:(height+clutter)));

			fprintf(fd2,"%f\t%f\n",KM_PER_MILE*path.distance[x],METERS_PER_FOOT*r);
			fprintf(fd5,"%f\t%f\n",KM_PER_MILE*path.distance[x],METERS_PER_FOOT*(height-terrain));
		
			
		}

		else
		{
			fprintf(fd,"%f\t%f\n",path.distance[x],height);

			if (fd1!=NULL && x>0 && x<path.length-2)
				fprintf(fd1,"%f\t%f\n",path.distance[x],(terrain==0.0?height:(height+clutter)));

			fprintf(fd2,"%f\t%f\n",path.distance[x],r);
			fprintf(fd5,"%f\t%f\n",path.distance[x],height-terrain);
		}

		
	
	
		if ((LR.frq_mhz>=20.0) && (LR.frq_mhz<=100000.0) && fresnel_plot)
		{
			if (metric)
			{
				fprintf(fd3,"%f\t%f\n",KM_PER_MILE*path.distance[x],METERS_PER_FOOT*f_zone);
				fprintf(fd4,"%f\t%f\n",KM_PER_MILE*path.distance[x],METERS_PER_FOOT*fpt6_zone);
			}

			else
			{
				fprintf(fd3,"%f\t%f\n",path.distance[x],f_zone);
				fprintf(fd4,"%f\t%f\n",path.distance[x],fpt6_zone);
			}

			if (f_zone<minheight)
				minheight=f_zone;
		}

		if ((height+clutter)>maxheight)
			maxheight=height+clutter;

		if (height<minheight)
			minheight=height;

		if (r>maxheight)
			maxheight=r;

		if (terrain<minterrain)
			minterrain=terrain;

		if ((height-terrain)<minearth)
			minearth=height-terrain;
	}

	if (normalized)
		r=-(nm*path.distance[path.length-1])-nb;
	else
		r=0.0;

	if (metric)
	{
		fprintf(fd,"%f\t%f\n",KM_PER_MILE*path.distance[path.length-1],METERS_PER_FOOT*r);
		fprintf(fd2,"%f\t%f\n",KM_PER_MILE*path.distance[path.length-1],METERS_PER_FOOT*r);
	}

	else
	{
		fprintf(fd,"%f\t%f\n",path.distance[path.length-1],r);
		fprintf(fd2,"%f\t%f\n",path.distance[path.length-1],r);
	}

	if ((LR.frq_mhz>=20.0) && (LR.frq_mhz<=100000.0) && fresnel_plot)
	{
		if (metric)
		{
			fprintf(fd3,"%f\t%f\n",KM_PER_MILE*path.distance[path.length-1],METERS_PER_FOOT*r);
			fprintf(fd4,"%f\t%f\n",KM_PER_MILE*path.distance[path.length-1],METERS_PER_FOOT*r);
		}

		else
		{
			fprintf(fd3,"%f\t%f\n",path.distance[path.length-1],r);
			fprintf(fd4,"%f\t%f\n",path.distance[path.length-1],r);
		}
	}
	
	if (r>maxheight)
		maxheight=r;

	if (r<minheight)
		minheight=r;

	fclose(fd);

	if (fd1!=NULL)
		fclose(fd1);

	fclose(fd2);
	fclose(fd5);

	if ((LR.frq_mhz>=20.0) && (LR.frq_mhz<=100000.0) && fresnel_plot)
	{
		fclose(fd3);
		fclose(fd4);
	}

	
	
	if (name[0]=='.')
	{
		/* Default filename and output file type */

		strncpy(basename,"profile\0",8);
		strncpy(term,"png\0",4);
		strncpy(ext,"png\0",4);
	}

	else
	{
		/* Extract extension and terminal type from "name" */

		ext[0]=0;
		y=strlen(name);
		strncpy(basename,name,254);

		for (x=y-1; x>0 && name[x]!='.'; x--);

		if (x>0)  /* Extension found */
		{
			for (z=x+1; z<=y && (z-(x+1))<10; z++)
			{
				ext[z-(x+1)]=tolower(name[z]);
				term[z-(x+1)]=name[z];
			}

			ext[z-(x+1)]=0;  /* Ensure an ending 0 */
			term[z-(x+1)]=0;
			basename[x]=0;
		}

		if (ext[0]==0)	/* No extension -- Default is png */
		{
			strncpy(term,"png\0",4);
			strncpy(ext,"png\0",4);
		}
	}

	/* Either .ps or .postscript may be used
	   as an extension for postscript output. */
	
	
	if (strncmp(term,"postscript",10)==0)
		strncpy(ext,"ps\0",3);

	else if (strncmp(ext,"ps",2)==0)
			strncpy(term,"postscript enhanced color\0",26);

	fd=fopen("ppa.gp","w");

	dheight=maxheight-minheight;
	miny=minheight-0.15*dheight;
	maxy=maxheight+0.05*dheight;

	if (maxy<20.0)
		maxy=20.0;

	dheight=maxheight-minheight;
	min2y=miny-minterrain+0.05*dheight;

	if (minearth<min2y)
	{
		miny-=min2y-minearth+0.05*dheight;
		min2y=minearth-0.05*dheight;
	}

	max2y=min2y+maxy-miny;
 
	fprintf(fd,"set grid\n");
	fprintf(fd,"set yrange [%2.3f to %2.3f]\n", metric?miny*METERS_PER_FOOT:miny, metric?maxy*METERS_PER_FOOT:maxy);
	fprintf(fd,"set y2range [%2.3f to %2.3f]\n", metric?min2y*METERS_PER_FOOT:min2y, metric?max2y*METERS_PER_FOOT:max2y);
	fprintf(fd,"set xrange [-0.5 to %2.3f]\n",metric?KM_PER_MILE*rint(distance+0.5):rint(distance+0.5));
	fprintf(fd,"set encoding iso_8859_1\n");
	fprintf(fd,"set term %s size %d, %d\n",term,pngwidth,pngheight);

	if ((LR.frq_mhz>=20.0) && (LR.frq_mhz<=100000.0) && fresnel_plot)
		fprintf(fd,"set title \"Link from %s to %s (%.2f%c)\"\n",source.name, destination.name, azimuth,176);

	else
		fprintf(fd,"set title Height Profile Between %s and %s (%.2f%c)\"\n", source.name, destination.name, azimuth,176);

	if (metric)
		fprintf(fd,"set xlabel \"Distance: %.2f Km Path Loss: %.2f dB Received Power: %.2f dBm\"\n",KM_PER_MILE*Distance(source,destination),loss,dBm);
	else
		fprintf(fd,"set xlabel \"Distance: %.2f Mi Path Loss: %.2f dB Received Power: %.2f dBm\"\n",Distance(source,destination),loss,dBm);

	
		if (metric)
			fprintf(fd,"set ylabel \"Height (meters)\"\n");

		else
			fprintf(fd,"set ylabel \"Height (feet)\"\n");
	

	fprintf(fd,"set output \"%s.%s\"\n",basename,ext);

	if ((LR.frq_mhz>=20.0) && (LR.frq_mhz<=100000.0) && fresnel_plot)
	{
		if (clutter>0.0)
		{
			if (metric)
				fprintf(fd,"plot \"profile.gp\" title \"Terrain Profile\" with filledcurve x1 lt 1 lc rgb \"#cbbc8a\", \"clutter.gp\" title \"Ground Clutter (%.2f meters)\" with lines, \"reference.gp\" title \"Line of Sight\" with lines, \"curvature.gp\" axes x1y2 title \"Earth's Curvature\" with lines, \"fresnel.gp\" axes x1y1 title \"First Fresnel Zone (%.3f MHz)\" with lines, \"fresnel_pt_6.gp\" title \"%.0f%% of Fresnel Zone\" with lines\n",clutter*METERS_PER_FOOT,LR.frq_mhz,fzone_clearance*100.0);
			else
				fprintf(fd,"plot \"profile.gp\" title \"Terrain Profile\" with filledcurve x1 lt 1 lc rgb \"#cbbc8a\", \"clutter.gp\" title \"Ground Clutter (%.2f feet)\" with lines, \"reference.gp\" title \"Line of Sight\" with lines, \"curvature.gp\" axes x1y2 title \"Earth's Curvature\" with lines, \"fresnel.gp\" axes x1y1 title \"First Fresnel Zone (%.3f MHz)\" with lines, \"fresnel_pt_6.gp\" title \"%.0f%% of Fresnel Zone\" with lines\n",clutter,LR.frq_mhz,fzone_clearance*100.0);
		}

		else
			fprintf(fd,"plot \"profile.gp\" title \"Terrain\" with filledcurve x1 lt 1 lc rgb \"#cbbc8a\", \"reference.gp\" title \"Line of sight\" with lines, \"curvature.gp\" axes x1y2 title \"Earth's Curvature\" with lines, \"fresnel.gp\" axes x1y1 title \"Fresnel Zone for %.1f MHz\" with lines, \"fresnel_pt_6.gp\" title \"%.0f%% of Fresnel Zone\" with lines\n",LR.frq_mhz,fzone_clearance*100.0);
	}

	else
	{
		
		

	}

	fclose(fd);

	x=system("gnuplot ppa.gp");

	
	
	if (x!=-1)
	{
		if (gpsav==0)
		{
			//unlink("ppa.gp");
			//unlink("profile.gp");
			//unlink("reference.gp");
			//unlink("curvature.gp");

			if (fd1!=NULL)
				unlink("clutter.gp");

			if ((LR.frq_mhz>=20.0) && (LR.frq_mhz<=100000.0) && fresnel_plot)
			{
				//unlink("fresnel.gp");
				//unlink("fresnel_pt_6.gp");
			}
		}

		
	}

	
	else
		fprintf(stderr,"\n*** ERROR: Error occurred invoking gnuplot!\n");
}

void ObstructionAnalysis(struct site xmtr, struct site rcvr, double f, FILE *outfile)
{
	/* Perform an obstruction analysis along the
	   path between receiver and transmitter. */

	int	x;
	struct	site site_x;
	double	h_r, h_t, h_x, h_r_orig, cos_tx_angle, cos_test_angle,
		cos_tx_angle_f1, cos_tx_angle_fpt6, d_tx, d_x,
		h_r_f1, h_r_fpt6, h_f, h_los, lambda=0.0;
	char	string[255], string_fpt6[255], string_f1[255];

	ReadPath(xmtr,rcvr);
	h_r=GetElevation(rcvr)+rcvr.alt+earthradius;
	h_r_f1=h_r;
	h_r_fpt6=h_r;
	h_r_orig=h_r;
	h_t=GetElevation(xmtr)+xmtr.alt+earthradius;
	d_tx=5280.0*Distance(rcvr,xmtr);
	cos_tx_angle=((h_r*h_r)+(d_tx*d_tx)-(h_t*h_t))/(2.0*h_r*d_tx);
	cos_tx_angle_f1=cos_tx_angle;
	cos_tx_angle_fpt6=cos_tx_angle;

	if (f)
		lambda=9.8425e8/(f*1e6);

	if (clutter>0.0)
	{
		fprintf(outfile,"Terrain has been raised by");

		if (metric)
			fprintf(outfile," %.2f meters",METERS_PER_FOOT*clutter);
		else
			fprintf(outfile," %.2f feet",clutter);

		fprintf(outfile," to account for ground clutter.\n\n");
	}

	/* At each point along the path calculate the cosine
	   of a sort of "inverse elevation angle" at the receiver.
	   From the antenna, 0 deg. looks at the ground, and 90 deg.
	   is parallel to the ground.

	   Start at the receiver.  If this is the lowest antenna,
	   then terrain obstructions will be nearest to it.  (Plus,
	   that's the way ppa!'s original los() did it.)

	   Calculate cosines only.  That's sufficient to compare
	   angles and it saves the extra computational burden of
	   acos().  However, note the inverted comparison: if
	   acos(A) > acos(B), then B > A. */

	for (x=path.length-1; x>0; x--)
	{
		site_x.lat=path.lat[x];
		site_x.lon=path.lon[x];
		site_x.alt=0.0;

		h_x=GetElevation(site_x)+earthradius+clutter;
		d_x=5280.0*Distance(rcvr,site_x);

		/* Deal with the LOS path first. */

		cos_test_angle=((h_r*h_r)+(d_x*d_x)-(h_x*h_x))/(2.0*h_r*d_x);

		if (cos_tx_angle>cos_test_angle)
		{
			if (h_r==h_r_orig)
				fprintf(outfile,"Between %s and %s, %s detected obstructions at:\n\n",rcvr.name,xmtr.name,ppa_name);

			if (site_x.lat>=0.0)
			{
				if (metric)
					fprintf(outfile,"   %8.4f N,%9.4f W, %5.2f kilometers, %6.2f meters AMSL\n",site_x.lat, site_x.lon, KM_PER_MILE*(d_x/5280.0), METERS_PER_FOOT*(h_x-earthradius));
				else
					fprintf(outfile,"   %8.4f N,%9.4f W, %5.2f miles, %6.2f feet AMSL\n",site_x.lat, site_x.lon, d_x/5280.0, h_x-earthradius);
			}

			else
			{
				if (metric)
					fprintf(outfile,"   %8.4f S,%9.4f W, %5.2f kilometers, %6.2f meters AMSL\n",-site_x.lat, site_x.lon, KM_PER_MILE*(d_x/5280.0), METERS_PER_FOOT*(h_x-earthradius));
				else

					fprintf(outfile,"   %8.4f S,%9.4f W, %5.2f miles, %6.2f feet AMSL\n",-site_x.lat, site_x.lon, d_x/5280.0, h_x-earthradius);
			}
		}

		while (cos_tx_angle>cos_test_angle)
		{
			h_r+=1;
			cos_test_angle=((h_r*h_r)+(d_x*d_x)-(h_x*h_x))/(2.0*h_r*d_x);
			cos_tx_angle=((h_r*h_r)+(d_tx*d_tx)-(h_t*h_t))/(2.0*h_r*d_tx);
		}

		if (f)
		{
			/* Now clear the first Fresnel zone... */

			cos_tx_angle_f1=((h_r_f1*h_r_f1)+(d_tx*d_tx)-(h_t*h_t))/(2.0*h_r_f1*d_tx);
			h_los=sqrt(h_r_f1*h_r_f1+d_x*d_x-2*h_r_f1*d_x*cos_tx_angle_f1);
			h_f=h_los-sqrt(lambda*d_x*(d_tx-d_x)/d_tx);

			while (h_f<h_x)
			{
				h_r_f1+=1;
				cos_tx_angle_f1=((h_r_f1*h_r_f1)+(d_tx*d_tx)-(h_t*h_t))/(2.0*h_r_f1*d_tx);
				h_los=sqrt(h_r_f1*h_r_f1+d_x*d_x-2*h_r_f1*d_x*cos_tx_angle_f1);
				h_f=h_los-sqrt(lambda*d_x*(d_tx-d_x)/d_tx);
			}

			/* and clear the 60% F1 zone. */

			cos_tx_angle_fpt6=((h_r_fpt6*h_r_fpt6)+(d_tx*d_tx)-(h_t*h_t))/(2.0*h_r_fpt6*d_tx);
			h_los=sqrt(h_r_fpt6*h_r_fpt6+d_x*d_x-2*h_r_fpt6*d_x*cos_tx_angle_fpt6);
			h_f=h_los-fzone_clearance*sqrt(lambda*d_x*(d_tx-d_x)/d_tx);

			while (h_f<h_x)
			{
				h_r_fpt6+=1;
				cos_tx_angle_fpt6=((h_r_fpt6*h_r_fpt6)+(d_tx*d_tx)-(h_t*h_t))/(2.0*h_r_fpt6*d_tx);
				h_los=sqrt(h_r_fpt6*h_r_fpt6+d_x*d_x-2*h_r_fpt6*d_x*cos_tx_angle_fpt6);
				h_f=h_los-fzone_clearance*sqrt(lambda*d_x*(d_tx-d_x)/d_tx);
			}
		}
	}
		
	if (h_r>h_r_orig)
	{
		if (metric)
			snprintf(string,150,"\nAntenna at %s must be raised to at least %.2f meters AGL\nto clear all obstructions detected by %s.\n",rcvr.name, METERS_PER_FOOT*(h_r-GetElevation(rcvr)-earthradius),ppa_name);
		else
			snprintf(string,150,"\nAntenna at %s must be raised to at least %.2f feet AGL\nto clear all obstructions detected by %s.\n",rcvr.name, h_r-GetElevation(rcvr)-earthradius,ppa_name);
	}

	else
		snprintf(string,150,"\nNo obstructions to LOS path due to terrain were detected by %s\n",ppa_name);

	if (f)
	{
		if (h_r_fpt6>h_r_orig)
		{
			if (metric)
				snprintf(string_fpt6,150,"\nAntenna at %s must be raised to at least %.2f meters AGL\nto clear %.0f%c of the first Fresnel zone.\n",rcvr.name, METERS_PER_FOOT*(h_r_fpt6-GetElevation(rcvr)-earthradius),fzone_clearance*100.0,37);

			else
				snprintf(string_fpt6,150,"\nAntenna at %s must be raised to at least %.2f feet AGL\nto clear %.0f%c of the first Fresnel zone.\n",rcvr.name, h_r_fpt6-GetElevation(rcvr)-earthradius,fzone_clearance*100.0,37);
		}

		else
			snprintf(string_fpt6,150,"\n%.0f%c of the first Fresnel zone is clear.\n",fzone_clearance*100.0,37);
	
		if (h_r_f1>h_r_orig)
		{
			if (metric)
				snprintf(string_f1,150,"\nAntenna at %s must be raised to at least %.2f meters AGL\nto clear the first Fresnel zone.\n",rcvr.name, METERS_PER_FOOT*(h_r_f1-GetElevation(rcvr)-earthradius));

			else			
				snprintf(string_f1,150,"\nAntenna at %s must be raised to at least %.2f feet AGL\nto clear the first Fresnel zone.\n",rcvr.name, h_r_f1-GetElevation(rcvr)-earthradius);

		}

		else
    		    snprintf(string_f1,150,"\nThe first Fresnel zone is clear.\n");
	}

	fprintf(outfile,"%s",string);

	if (f)
	{
		fprintf(outfile,"%s",string_f1);
		fprintf(outfile,"%s",string_fpt6);
	}
	
	
}

void PathReport(struct site source, struct site destination, char *name, char graph_it)
{
	/* This function writes a PPA Path Report (name.txt) to
	   the filesystem.  If (graph_it == 1), then gnuplot is invoked
	   to generate an appropriate output file indicating the Longley-Rice
	   model loss between the source and destination locations.
    	   "filename" is the name assigned to the output file generated
	   by gnuplot.  The filename extension is used to set gnuplot's
	   terminal setting and output file type.  If no extension is
	   found, .png is assumed. */

	int	x, y, z, errnum;
	char	basename[255], term[30], ext[15], strmode[100],
		report_name[80], block=0;
	double	maxloss=-100000.0, minloss=100000.0, haavt,
		angle1, angle2, azimuth, pattern=1.0, patterndB=0.0,
		total_loss=0.0, cos_xmtr_angle, cos_test_angle=0.0,
		source_alt, test_alt, dest_alt, source_alt2, dest_alt2,
		distance, elevation, four_thirds_earth,
		free_space_loss=0.0, eirp=0.0, voltage, rxp, power_density;
	FILE	*fd=NULL, *fd2=NULL;

	//sprintf(report_name,"%s.txt",*name);
	snprintf(report_name,80,"%s.txt%c",name,0);



	four_thirds_earth=FOUR_THIRDS*EARTHRADIUS;

	/*for (x=0; report_name[x]!=0; x++)
		if (report_name[x]==32 || report_name[x]==17 || report_name[x]==92 || report_name[x]==42 || report_name[x]==47)
			report_name[x]='_';	*/

	fd2=fopen(report_name,"w");

	fprintf(fd2,"\n\t\t--==[ %s v%s Path Analysis ]==--\n\n",ppa_name,ppa_version);
	//fprintf(fd2,"%s\n\n",dashes);
	fprintf(fd2,"Transmitter site: %s\n",source.name);

	if (source.lat>=0.0)
	{
		fprintf(fd2,"Site location: %.4f North / %.4f West",source.lat, source.lon);
		//fprintf(fd2, " (%s N / ", source.lat);
	}

	else
	{

		fprintf(fd2,"Site location: %.4f South / %.4f West",-source.lat, source.lon);
		//fprintf(fd2, " (%s S / ", source.lat);
	}
	
	//fprintf(fd2, "%s W)\n", source.lon);

	if (metric)
	{
		fprintf(fd2,"Ground elevation: %.2f meters AMSL\n",METERS_PER_FOOT*GetElevation(source));
		fprintf(fd2,"Antenna height: %.2f meters AGL / %.2f meters AMSL\n",METERS_PER_FOOT*source.alt,METERS_PER_FOOT*(source.alt+GetElevation(source)));
	}

	else
	{
		fprintf(fd2,"Ground elevation: %.2f feet AMSL\n",GetElevation(source));
		fprintf(fd2,"Antenna height: %.2f feet AGL / %.2f feet AMSL\n",source.alt, source.alt+GetElevation(source));
	}

/*
	haavt=haat(source);

	if (haavt>-4999.0)
	{
		if (metric)
			fprintf(fd2,"Antenna height above average terrain: %.2f meters\n",METERS_PER_FOOT*haavt);
		else
			fprintf(fd2,"Antenna height above average terrain: %.2f feet\n",haavt);
	}
*/
	azimuth=Azimuth(source,destination);
	angle1=ElevationAngle(source,destination);
	angle2=ElevationAngle2(source,destination,earthradius);

	if (got_azimuth_pattern || got_elevation_pattern)
	{
		x=(int)rint(10.0*(10.0-angle2));

		if (x>=0 && x<=1000)
			pattern=(double)LR.antenna_pattern[(int)rint(azimuth)][x];

		patterndB=20.0*log10(pattern);
	}

	if (metric)
		fprintf(fd2,"Distance to %s: %.2f kilometers\n",destination.name,KM_PER_MILE*Distance(source,destination));

	else
		fprintf(fd2,"Distance to %s: %.2f miles\n",destination.name,Distance(source,destination));

	fprintf(fd2,"Azimuth to %s: %.2f degrees\n",destination.name,azimuth);

	if (angle1>=0.0)
		fprintf(fd2,"Elevation angle to %s: %+.4f degrees\n",destination.name,angle1);

	else
		fprintf(fd2,"Depression angle to %s: %+.4f degrees\n",destination.name,angle1);

	if ((angle2-angle1)>0.0001)
	{
		if (angle2<0.0)
			fprintf(fd2,"Depression");
		else
			fprintf(fd2,"Elevation");

		//fprintf(fd2," angle to the first obstruction: %+.4f degrees\n",angle2);
	}

	//fprintf(fd2,"\n%s\n\n",dashes);

	/* Receiver */

	fprintf(fd2,"Receiver site: %s\n",destination.name);

	if (destination.lat>=0.0)
	{
		fprintf(fd2,"Site location: %.4f North / %.4f West",destination.lat, destination.lon);
		//fprintf(fd2, " (%s N / ", destination.lat);
	}

	else
	{
		fprintf(fd2,"Site location: %.4f South / %.4f West",-destination.lat, destination.lon);
		//fprintf(fd2, " (%s S / ", destination.lat);
	}

	//fprintf(fd2, "%s W)\n", destination.lon);

	if (metric)
	{
		fprintf(fd2,"Ground elevation: %.2f meters AMSL\n",METERS_PER_FOOT*GetElevation(destination));
		fprintf(fd2,"Antenna height: %.2f meters AGL / %.2f meters AMSL\n",METERS_PER_FOOT*destination.alt, METERS_PER_FOOT*(destination.alt+GetElevation(destination)));
	}

	else
	{
		fprintf(fd2,"Ground elevation: %.2f feet AMSL\n",GetElevation(destination));
		fprintf(fd2,"Antenna height: %.2f feet AGL / %.2f feet AMSL\n",destination.alt, destination.alt+GetElevation(destination));
	}

	/*haavt=haat(destination);

	if (haavt>-4999.0)
	{
		if (metric)
			fprintf(fd2,"Antenna height above average terrain: %.2f meters\n",METERS_PER_FOOT*haavt);
		else
			fprintf(fd2,"Antenna height above average terrain: %.2f feet\n",haavt);
	}*/

	if (metric)
		fprintf(fd2,"Distance to %s: %.2f kilometers\n",source.name,KM_PER_MILE*Distance(source,destination));

	else
		fprintf(fd2,"Distance to %s: %.2f miles\n",source.name,Distance(source,destination));

	azimuth=Azimuth(destination,source);

	angle1=ElevationAngle(destination,source);
	angle2=ElevationAngle2(destination,source,earthradius);

	fprintf(fd2,"Azimuth to %s: %.2f degrees\n",source.name,azimuth);

	if (angle1>=0.0)
		fprintf(fd2,"Elevation angle to %s: %+.4f degrees\n",source.name,angle1);

	else
		fprintf(fd2,"Depression angle to %s: %+.4f degrees\n",source.name,angle1);

	if ((angle2-angle1)>0.0001)
	{
		if (angle2<0.0)
			fprintf(fd2,"Depression");
		else
			fprintf(fd2,"Elevation");

		//fprintf(fd2," angle to the first obstruction: %+.4f degrees\n",angle2);
	}

	//fprintf(fd2,"\n%s\n\n",dashes);

	if (LR.frq_mhz>0.0)
	{
		fprintf(fd2,"Longley-Rice path calculation parameters used in this analysis:\n\n");
		fprintf(fd2,"Earth's Dielectric Constant: %.3lf\n",LR.eps_dielect);
		fprintf(fd2,"Earth's Conductivity: %.3lf Siemens/meter\n",LR.sgm_conductivity);
		fprintf(fd2,"Atmospheric Bending Constant (N-units): %.3lf ppm\n",LR.eno_ns_surfref);
		fprintf(fd2,"Frequency: %.3lf MHz\n",LR.frq_mhz);
		fprintf(fd2,"Radio Climate: %d (",LR.radio_climate);

		switch (LR.radio_climate)
		{
			case 1:
			fprintf(fd2,"Equatorial");
			break;

			case 2:
			fprintf(fd2,"Continental Subtropical");
			break;

			case 3:
			fprintf(fd2,"Maritime Subtropical");
			break;

			case 4:
			fprintf(fd2,"Desert");
			break;

			case 5:
			fprintf(fd2,"Continental Temperate");
			break;

			case 6:
			fprintf(fd2,"Martitime Temperate, Over Land");
			break;

			case 7:
			fprintf(fd2,"Maritime Temperate, Over Sea");
			break;

			default:
			fprintf(fd2,"Unknown");
		}

		fprintf(fd2,")\nPolarisation: %d (",LR.pol);

		if (LR.pol==0)
			fprintf(fd2,"Horizontal");

		if (LR.pol==1)
			fprintf(fd2,"Vertical");

		fprintf(fd2,")\nFraction of Situations: %.1lf%c\n",LR.conf*100.0,37);
		fprintf(fd2,"Fraction of Time: %.1lf%c\n",LR.rel*100.0,37);
	
		if (LR.erp!=0.0)
		{
			fprintf(fd2,"Transmitter ERP: ");

			if (LR.erp<1.0)
				fprintf(fd2,"%.1lf milliwatts",1000.0*LR.erp);

			if (LR.erp>=1.0 && LR.erp<10.0)
				fprintf(fd2,"%.1lf Watts",LR.erp);

			if (LR.erp>=10.0 && LR.erp<10.0e3)
				fprintf(fd2,"%.0lf Watts",LR.erp);

			if (LR.erp>=10.0e3)
				fprintf(fd2,"%.3lf kilowatts",LR.erp/1.0e3);

			dBm=10.0*(log10(LR.erp*1000.0));
			fprintf(fd2," (%+.2f dBm)\n",dBm);

			/* EIRP = ERP + 2.14 dB */

			fprintf(fd2,"Transmitter EIRP: ");

			eirp=LR.erp*1.636816521;

			if (eirp<1.0)
				fprintf(fd2,"%.1lf milliwatts",1000.0*eirp);

			if (eirp>=1.0 && eirp<10.0)
				fprintf(fd2,"%.1lf Watts",eirp);

			if (eirp>=10.0 && eirp<10.0e3)
				fprintf(fd2,"%.0lf Watts",eirp);

			if (eirp>=10.0e3)
				fprintf(fd2,"%.3lf kilowatts",eirp/1.0e3);

			dBm=10.0*(log10(eirp*1000.0));
			fprintf(fd2," (%+.2f dBm)\n",dBm);
		}

		fprintf(fd2,"\n%s\n\n",dashes);

		fprintf(fd2,"Summary for the link between %s and %s:\n\n",source.name, destination.name);

		if (patterndB!=0.0)
			fprintf(fd2,"%s antenna pattern towards %s: %.3f (%.2f dB)\n", source.name, destination.name, pattern, patterndB);

		ReadPath(source, destination);  /* source=TX, destination=RX */

		/* Copy elevations plus clutter along
		   path into the elev[] array. */

		for (x=1; x<path.length-1; x++)
			elev[x+2]=METERS_PER_FOOT*(path.elevation[x]==0.0?path.elevation[x]:(clutter+path.elevation[x]));

		/* Copy ending points without clutter */

		elev[2]=path.elevation[0]*METERS_PER_FOOT;
		elev[path.length+1]=path.elevation[path.length-1]*METERS_PER_FOOT;

		fd=fopen("profile.gp","w");

		azimuth=rint(Azimuth(source,destination));

		for (y=2; y<(path.length-1); y++)  /* path.length-1 avoids LR error */
		{
			distance=5280.0*path.distance[y];
			source_alt=four_thirds_earth+source.alt+path.elevation[0];
			dest_alt=four_thirds_earth+destination.alt+path.elevation[y];
			dest_alt2=dest_alt*dest_alt;
			source_alt2=source_alt*source_alt;

			/* Calculate the cosine of the elevation of
			   the receiver as seen by the transmitter. */

			cos_xmtr_angle=((source_alt2)+(distance*distance)-(dest_alt2))/(2.0*source_alt*distance);

			if (got_elevation_pattern)
			{
				/* If an antenna elevation pattern is available, the
				   following code determines the elevation angle to
			   	   the first obstruction along the path. */

				for (x=2, block=0; x<y && block==0; x++)
				{
					distance=5280.0*(path.distance[y]-path.distance[x]);
					test_alt=four_thirds_earth+path.elevation[x];

					/* Calculate the cosine of the elevation
					   angle of the terrain (test point)
					   as seen by the transmitter. */

					cos_test_angle=((source_alt2)+(distance*distance)-(test_alt*test_alt))/(2.0*source_alt*distance);

					/* Compare these two angles to determine if
					   an obstruction exists.  Since we're comparing
					   the cosines of these angles rather than
					   the angles themselves, the sense of the
					   following "if" statement is reversed from
				   	   what it would be if the angles themselves
				   	   were compared. */

					if (cos_xmtr_angle>=cos_test_angle)
						block=1;
				}

				/* At this point, we have the elevation angle
				   to the first obstruction (if it exists). */
			}

			/* Determine path loss for each point along the
			   path using Longley-Rice's point_to_point mode
		  	   starting at x=2 (number_of_points = 1), the
		  	   shortest distance terrain can play a role in
		  	   path loss. */

			elev[0]=y-1;	/* (number of points - 1) */

			/* Distance between elevation samples */

			elev[1]=METERS_PER_MILE*(path.distance[y]-path.distance[y-1]);

			point_to_point(elev, source.alt*METERS_PER_FOOT, 
			destination.alt*METERS_PER_FOOT, LR.eps_dielect,
			LR.sgm_conductivity, LR.eno_ns_surfref, LR.frq_mhz,
			LR.radio_climate, LR.pol, LR.conf, LR.rel, loss,
			strmode, errnum);

			if (block)
				elevation=((acos(cos_test_angle))/DEG2RAD)-90.0;
			else
				elevation=((acos(cos_xmtr_angle))/DEG2RAD)-90.0;

			/* Integrate the antenna's radiation
			   pattern into the overall path loss. */

			x=(int)rint(10.0*(10.0-elevation));

			if (x>=0 && x<=1000)
			{
				pattern=(double)LR.antenna_pattern[(int)azimuth][x];

				if (pattern!=0.0)
					patterndB=20.0*log10(pattern);
			}

			else
				patterndB=0.0;

			total_loss=loss-patterndB;

			if (metric)
				fprintf(fd,"%f\t%f\n",KM_PER_MILE*(path.distance[path.length-1]-path.distance[y]),total_loss);

			else
				fprintf(fd,"%f\t%f\n",path.distance[path.length-1]-path.distance[y],total_loss);

			if (total_loss>maxloss)
				maxloss=total_loss;

			if (total_loss<minloss)
				minloss=total_loss;
		}

		fclose(fd);

		distance=Distance(source,destination);


		if (distance!=0.0)
		{
			free_space_loss=36.6+(20.0*log10(LR.frq_mhz))+(20.0*log10(distance));

			fprintf(fd2,"Free space path loss: %.2f dB\n",free_space_loss);
		}

		fprintf(fd2,"Longley-Rice path loss: %.2f dB\n",loss);

		if (free_space_loss!=0.0)
			fprintf(fd2,"Attenuation due to terrain shielding: %.2f dB\n",loss-free_space_loss);

		if (patterndB!=0.0)
			fprintf(fd2,"Total path loss including %s antenna pattern: %.2f dB\n",source.name,total_loss);

		if (LR.erp!=0.0)
		{
			field_strength=(139.4+(20.0*log10(LR.frq_mhz))-total_loss)+(10.0*log10(LR.erp/1000.0));

			/* dBm is referenced to EIRP */

			rxp=eirp/(pow(10.0,(total_loss/10.0)));
			dBm=10.0*(log10(rxp*1000.0));
			power_density=(eirp/(pow(10.0,(total_loss-free_space_loss)/10.0)));
			/* divide by 4*PI*distance_in_meters squared */
			power_density/=(4.0*PI*distance*distance*2589988.11);

			fprintf(fd2,"Field strength at %s: %.2f dBuV/meter\n", destination.name,field_strength);
			fprintf(fd2,"Signal power level at %s: %+.2f dBm\n",destination.name,dBm);
			fprintf(fd2,"Signal power density at %s: %+.2f dBW per square meter\n",destination.name,10.0*log10(power_density));
			voltage=1.0e6*sqrt(50.0*(eirp/(pow(10.0,(total_loss-2.14)/10.0))));
			fprintf(fd2,"Voltage across 50 ohm dipole at %s: %.2f uV (%.2f dBuV)\n",destination.name,voltage,20.0*log10(voltage));

			voltage=1.0e6*sqrt(75.0*(eirp/(pow(10.0,(total_loss-2.14)/10.0))));
			fprintf(fd2,"Voltage across 75 ohm dipole at %s: %.2f uV (%.2f dBuV)\n",destination.name,voltage,20.0*log10(voltage));
		}

		fprintf(fd2,"Mode of propagation: %s\n",strmode);
		fprintf(fd2,"Longley-Rice model error number: %d",errnum);

		switch (errnum)
		{
			case 0:
				fprintf(fd2," (No error)\n");
				break;

			case 1:
				fprintf(fd2,"\n  Warning: Some parameters are nearly out of range.\n");
				fprintf(fd2,"  Results should be used with caution.\n");
				break;

			case 2:
				fprintf(fd2,"\n  Note: Default parameters have been substituted for impossible ones.\n");
				break;

			case 3:
				fprintf(fd2,"\n  Warning: A combination of parameters is out of range.\n");
				fprintf(fd2,"  Results are probably invalid.\n");
				break;

			default:
				fprintf(fd2,"\n  Warning: Some parameters are out of range.\n");
				fprintf(fd2,"  Results are probably invalid.\n");
		}

		fprintf(fd2,"\n%s\n\n",dashes);
	}



	ObstructionAnalysis(source, destination, LR.frq_mhz, fd2);

	fclose(fd2);

	
	
	/* Skip plotting the graph if ONLY a path-loss report is needed. */

	if (graph_it)
	{
		if (name[0]=='.')
		{
			/* Default filename and output file type */

			strncpy(basename,"profile\0",8);
			strncpy(term,"png\0",4);
			strncpy(ext,"png\0",4);
		}

		else
		{
			/* Extract extension and terminal type from "name" */

			ext[0]=0;
			y=strlen(name);
			strncpy(basename,name,254);

			for (x=y-1; x>0 && name[x]!='.'; x--);

			if (x>0)  /* Extension found */
			{
				for (z=x+1; z<=y && (z-(x+1))<10; z++)
				{
					ext[z-(x+1)]=tolower(name[z]);
					term[z-(x+1)]=name[z];
				}

				ext[z-(x+1)]=0;  /* Ensure an ending 0 */
				term[z-(x+1)]=0;
				basename[x]=0;
			}
		}

		if (ext[0]==0)	/* No extension -- Default is png */
		{
			strncpy(term,"png\0",4);
			strncpy(ext,"png\0",4);
		}

		/* Either .ps or .postscript may be used
		   as an extension for postscript output. */

		if (strncmp(term,"postscript",10)==0)
			strncpy(ext,"ps\0",3);

		else if (strncmp(ext,"ps",2)==0)
				strncpy(term,"postscript enhanced color\0",26);

		fd=fopen("ppa.gp","w");

		fprintf(fd,"set grid\n");
		fprintf(fd,"set yrange [%2.3f to %2.3f]\n", minloss, maxloss);
		fprintf(fd,"set encoding iso_8859_1\n");
		fprintf(fd,"set term %s\n",term);
		fprintf(fd,"set title \"%s Loss Profile Along Path Between %s and %s (%.2f%c azimuth)\"\n",ppa_name, destination.name, source.name, Azimuth(destination,source),176);

		if (metric)
			fprintf(fd,"set xlabel \"Distance Between %s and %s (%.2f kilometers)\"\n",destination.name,source.name,KM_PER_MILE*Distance(destination,source));
		else
			fprintf(fd,"set xlabel \"Distance Between %s and %s (%.2f miles)\"\n",destination.name,source.name,Distance(destination,source));

		if (got_azimuth_pattern || got_elevation_pattern)
			fprintf(fd,"set ylabel \"Total Path Loss (including TX antenna pattern) (dB)");
		else
			fprintf(fd,"set ylabel \"Longley-Rice Path Loss (dB)");

		fprintf(fd,"\"\nset output \"%s.%s\"\n",basename,ext);
		fprintf(fd,"plot \"profile.gp\" title \"Path Loss\" with lines\n");

		fclose(fd);
			
		x=system("gnuplot ppa.gp");

		if (x!=-1)
		{
			if (gpsav==0)
			{
				//unlink("ppa.gp");
				//unlink("profile.gp");
				//unlink("reference.gp");
			}	

			
		}

		else
			fprintf(stderr,"\n*** ERROR: Error occurred invoking gnuplot!\n");
	}


	
	
}

void SiteReport(struct site xmtr)
{
	char	report_name[80];
	double	terrain;
	int	x, azi;
	FILE	*fd;

	sprintf(report_name,"%s-site_report.txt",xmtr.name);

	for (x=0; report_name[x]!=0; x++)
		if (report_name[x]==32 || report_name[x]==17 || report_name[x]==92 || report_name[x]==42 || report_name[x]==47)
			report_name[x]='_';	

	fd=fopen(report_name,"w");

	fprintf(fd,"\n\t--==[ %s v%s Site Analysis Report For: %s ]==--\n\n",ppa_name, ppa_version, xmtr.name);

	fprintf(fd,"%s\n\n",dashes);

	if (xmtr.lat>=0.0)
	{
		fprintf(fd,"Site location: %.4f North / %.4f West",xmtr.lat, xmtr.lon);
		fprintf(fd, " (%s N / ",xmtr.lat);
	}

	else
	{
		fprintf(fd,"Site location: %.4f South / %.4f West",-xmtr.lat, xmtr.lon);
		fprintf(fd, " (%s S / ",xmtr.lat);
	}

	fprintf(fd, "%s W)\n",xmtr.lon);

	if (metric)
	{
		fprintf(fd,"Ground elevation: %.2f meters AMSL\n",METERS_PER_FOOT*GetElevation(xmtr));
		fprintf(fd,"Antenna height: %.2f meters AGL / %.2f meters AMSL\n",METERS_PER_FOOT*xmtr.alt, METERS_PER_FOOT*(xmtr.alt+GetElevation(xmtr)));
	}

	else
	{
		fprintf(fd,"Ground elevation: %.2f feet AMSL\n",GetElevation(xmtr));
		fprintf(fd,"Antenna height: %.2f feet AGL / %.2f feet AMSL\n",xmtr.alt, xmtr.alt+GetElevation(xmtr));
	}

	terrain=haat(xmtr);

	if (terrain>-4999.0)
	{
		if (metric)
			fprintf(fd,"Antenna height above average terrain: %.2f meters\n\n",METERS_PER_FOOT*terrain);
		else
			fprintf(fd,"Antenna height above average terrain: %.2f feet\n\n",terrain);

		/* Display the average terrain between 2 and 10 miles
		   from the transmitter site at azimuths of 0, 45, 90,
		   135, 180, 225, 270, and 315 degrees. */

		for (azi=0; azi<=315; azi+=45)
		{
			fprintf(fd,"Average terrain at %3d degrees azimuth: ",azi);
			terrain=AverageTerrain(xmtr,(double)azi,2.0,10.0);

			if (terrain>-4999.0)
			{
				if (metric)
					fprintf(fd,"%.2f meters AMSL\n",METERS_PER_FOOT*terrain);
				else
					fprintf(fd,"%.2f feet AMSL\n",terrain);
			}

			else
				fprintf(fd,"No terrain\n");
		}
	}

	fprintf(fd,"\n%s\n\n",dashes);
	fclose(fd);
	
}

void LoadTopoData(int max_lon, int min_lon, int max_lat, int min_lat, int winfiles)
{
	/* This function loads the SDF files required
	   to cover the limits of the region specified. */ 

	int x, y, width, ymin, ymax;

	width=ReduceAngle(max_lon-min_lon);

	if ((max_lon-min_lon)<=180.0)
	{
		for (y=0; y<=width; y++)
			for (x=min_lat; x<=max_lat; x++)
			{
				ymin=(int)(min_lon+(double)y);

				while (ymin<0)
					ymin+=360;

				while (ymin>=360)
					ymin-=360;

				ymax=ymin+1;

				while (ymax<0)
					ymax+=360;

				while (ymax>=360)
					ymax-=360;

				if (winfiles==1){
					if (ippd==3600)
						snprintf(string,19,"%d=%d=%d=%d=hd",x, x+1, ymin, ymax);
					else
						snprintf(string,16,"%d=%d=%d=%d",x, x+1, ymin, ymax);

				}else{
					if (ippd==3600)
						snprintf(string,19,"%d:%d:%d:%d=hd",x, x+1, ymin, ymax);
					else
						snprintf(string,16,"%d:%d:%d:%d",x, x+1, ymin, ymax);
				}  
				LoadSDF(string,winfiles);
			}
	}

	else
	{
		for (y=0; y<=width; y++)
			for (x=min_lat; x<=max_lat; x++)
			{
				ymin=max_lon+y;

				while (ymin<0)
					ymin+=360;

				while (ymin>=360)
					ymin-=360;
					
				ymax=ymin+1;

				while (ymax<0)
					ymax+=360;

				while (ymax>=360)
					ymax-=360;

				if (winfiles==1){
					if (ippd==3600)
						snprintf(string,19,"%d=%d=%d=%d=hd",x, x+1, ymin, ymax);
					else
						snprintf(string,16,"%d=%d=%d=%d",x, x+1, ymin, ymax);

				}else{
					if (ippd==3600)
						snprintf(string,19,"%d:%d:%d:%d=hd",x, x+1, ymin, ymax);
					else
						snprintf(string,16,"%d:%d:%d:%d",x, x+1, ymin, ymax);
				}  
				LoadSDF(string,winfiles);
			}
	}
}


int main(int argc, char *argv[])
{
	int		x, y, z=0, min_lat, min_lon, max_lat, max_lon,
			rxlat, rxlon, txlat, txlon, width=1000, height=350, winfiles=0;

	unsigned char	height_plot=0, norm=0, pt2pt_mode=1,
			max_txsites, nolospath=0, nositereports=0, fresnel_plot=1;
 
	char		header[80], height_file[255], return_height_file[255], 
			longley_file[255], string[255], ext[20];

	double		er_mult,txla=999, txlo=999, rxla=999, rxlo=999, txh=0, rxh=0;

	struct		site tx_site[32], rx_site;

	strncpy(ppa_version,"1.3.1\0",5);

	strncpy(ppa_name,"PPA\0",4);

	strncpy(dashes,"---------------------------------------------------------------------------\0",76);

	if (argc==1)
	{
		fprintf(stdout,"\n\t\t --==[ %s v%s Available Options... ]==--\n\n",ppa_name, ppa_version);

		fprintf(stdout,"       -m Metric units (Default = Imperial)\n");			
		fprintf(stdout,"       -tla Tx latitude)\n");	
		fprintf(stdout,"       -tlo Tx longitude (W) Positive value 0 to 360)\n");		
		fprintf(stdout,"       -th Tx Height)\n");	
		fprintf(stdout,"       -rla Rx latitude)\n");	
		fprintf(stdout,"       -rlo Rx longitude (W) Positive value 0 to 360)\n");		
		fprintf(stdout,"       -rh Rx Height)\n");
		fprintf(stdout,"       -fz Fresnel zone clearance percentage (default = 60)\n");	
		fprintf(stdout,"       -f frequency for Fresnel zone calculation (MHz)\n");
		fprintf(stdout,"       -w Effective radiated power in Watts. Default=0\n");
		fprintf(stdout,"       -p Polarisation. Default=1 (Vertical). 0=Horizontal\n");
		fprintf(stdout,"       -d sdf file directory path (overrides path in ~/.ppa_path file)\n");
		fprintf(stdout,"       -x PNG graph width (default = 800)\n");
		fprintf(stdout,"       -y PNG graph height (default = 200)\n");
		fprintf(stdout,"       -n Normalise graph\n");
		fprintf(stdout,"       -v Output value. 1=Path loss dB, 2=Rxd power dBm (default),3=Field strength dBuV/m\n");
		fprintf(stdout,"       -o PNG filename. Return PNG will be called $filename_R.png\n");
		fprintf(stdout,"       -wf Windows SDF tiles with equals not colons\n");
	
		
	}

	y=argc-1;

	metric=0;
	sdf_path[0]=0;
	path.length=0;
	max_txsites=2;
	fzone_clearance=0.6;
	contour_threshold=0;
	rx_site.lat=91.0;
	rx_site.lon=361.0;
	longley_file[0]=0;
	earthradius=EARTHRADIUS;
	
	LR.eps_dielect=15.0;
	LR.sgm_conductivity=0.005;
	LR.eno_ns_surfref=301.0;
	LR.frq_mhz=300.0; // 
	LR.radio_climate=5;
	LR.pol=1; //
	LR.conf=0.50;
	LR.rel=0.50;
	LR.erp=0.0; //

	ippd=IPPD;		/* pixels per degree (integer) */
	ppd=(double)ippd;	/* pixels per degree (double)  */
	dpp=1.0/ppd;		/* degrees per pixel */
	mpi=ippd-1;		/* maximum pixel index per degree */


	for (x=0; x<4; x++)
	{
		tx_site[x].lat=91.0;
		tx_site[x].lon=361.0;
	}

	for (x=0; x<MAXPAGES; x++)
	{
		dem[x].min_el=32768;
		dem[x].max_el=-32768;
		dem[x].min_north=90;
		dem[x].max_north=-90;
		dem[x].min_west=360;
		dem[x].max_west=-1;
	}

	/* Scan for command line arguments */

	for (x=1; x<=y; x++)
	{
		
		if (strcmp(argv[x],"-m")==0)
		{
			z=x+1;

			if (z<=y && argv[z][0] && argv[z][0]!='-')
			{
				sscanf(argv[z],"%lf",&er_mult);

				if (er_mult<0.1)
					er_mult=1.0;

				if (er_mult>1.0e6)
					er_mult=1.0e6;

				earthradius*=er_mult;
			}			 
		}

		
		if (strcmp(argv[x],"-fz")==0)
		{
			z=x+1;

			if (z<=y && argv[z][0] && argv[z][0]!='-')
			{
				sscanf(argv[z],"%lf",&fzone_clearance);

				if (fzone_clearance<0.0 || fzone_clearance>100.0)
					fzone_clearance=60.0;

				fzone_clearance/=100.0;
			}
		}

		if (strcmp(argv[x],"-x")==0)
		{
			z=x+1;

			if (z<=y && argv[z][0] && argv[z][0]!='-')
			{
				sscanf(argv[z],"%d",&width);

			}
		}
		if (strcmp(argv[x],"-y")==0)
		{
			z=x+1;

			if (z<=y && argv[z][0] && argv[z][0]!='-')
			{
				sscanf(argv[z],"%d",&height);

			}
		}

		


		
			if (strcmp(argv[x],"-o")==0)
			{
				z=x+1;
				if (z<=y && argv[z][0] && argv[z][0]!='-')
				{
				strncpy(height_file,argv[z],253);
				height_plot=1;
				pt2pt_mode=1;
				}
			}
			
			
		

		if (strcmp(argv[x],"-m")==0)
			metric=1;
		if (strcmp(argv[x],"-n")==0)
			norm=1;

		if (strcmp(argv[x],"-d")==0)
		{
			z=x+1;

			if (z<=y && argv[z][0] && argv[z][0]!='-')
				strncpy(sdf_path,argv[z],253);
		}

		if (strcmp(argv[x],"-tla")==0)
		{
			/* Read Transmitter Lat */

			z=x+1;

			if (z<=y && argv[z][0]) 
			{
				//strncpy(txla,argv[z],253);
				sscanf(argv[z],"%lf",&txla);
				
			}

		}

		if (strcmp(argv[x],"-tlo")==0)
                {
                        /* Read Transmitter Lon */

                        z=x+1;

                        if (z<=y && argv[z][0])
                        {
                                //strncpy(txla,argv[z],253);
                                sscanf(argv[z],"%lf",&txlo);

                        }

                }

		if (strcmp(argv[x],"-rla")==0)
                {
                        /* Read Rx Lat */

                        z=x+1;

                        if (z<=y && argv[z][0])
                        {
                                //strncpy(txla,argv[z],253);
                                sscanf(argv[z],"%lf",&rxla);

                        }

                }
		if (strcmp(argv[x],"-rlo")==0)
                {
                        /* Read Rx lon */

                        z=x+1;

                        if (z<=y && argv[z][0])
                        {
                                //strncpy(txla,argv[z],253);
                                sscanf(argv[z],"%lf",&rxlo);

                        }

                }



		if (strcmp(argv[x],"-th")==0)
                {
                        /* Read Transmitter height */

                        z=x+1;

                        if (z<=y && argv[z][0] && argv[z][0]!='-')
                        {
                                //strncpy(txla,argv[z],253);
                                sscanf(argv[z],"%lf",&txh);

                        }


                }

		if (strcmp(argv[x],"-rh")==0)
                {
                        /* Read Rx height */

                        z=x+1;

                        if (z<=y && argv[z][0] && argv[z][0]!='-')
                        {
                               
                                sscanf(argv[z],"%lf",&rxh);

                        }

                }


		
		if (strcmp(argv[x],"-f")==0)
		{
			z=x+1;

			if (z<=y && argv[z][0] && argv[z][0]!='-')
			{
				sscanf(argv[z],"%lf",&LR.frq_mhz);

				if (LR.frq_mhz<20.0)
					LR.frq_mhz=0.0;

				if (LR.frq_mhz>100.0e3)
					LR.frq_mhz=100.0e3;
			}
		}

		if (strcmp(argv[x],"-p")==0)
		{
			z=x+1;

			if (z<=y && argv[z][0] && argv[z][0]!='-')
			{
				sscanf(argv[z],"%d",&LR.pol);

			}
		}
		
		if (strcmp(argv[x],"-w")==0)
		{
			z=x+1;

			if (z<=y && argv[z][0] && argv[z][0]>0)
			{
				sscanf(argv[z],"%lf",&LR.erp);

			}
		}
		
			if (strcmp(argv[x],"-v")==0)
		{
			z=x+1;

			if (z<=y && argv[z][0] && argv[z][0]>0)
			{
				sscanf(argv[z],"%d",&output);

			}
		}
		
		//Windows friendly SDF filenames
			if (strcmp(argv[x],"-wf")==0)
		{
			z=x+1;

			winfiles=1;
		}
 	
}	



if(txla==999 || txlo==999 || rxla==999 || rxlo==999){
fprintf(stdout,"\nERROR: BAD LOCATION GIVEN\n");
exit(1);
}

if (metric)
	{
		txh/=METERS_PER_FOOT;	        /* meters --> feet */
		rxh/=METERS_PER_FOOT;		/* kilometers --> miles */
	}

// populate sites
tx_site[0].lat = txla;
tx_site[0].lon = txlo;
tx_site[0].alt = txh;


                        
						
strncpy(tx_site[0].name,"A",2);

rx_site.lat = rxla;
rx_site.lon = rxlo;
rx_site.alt = rxh;
						
strncpy(rx_site.name,"B",2);


// 100 mile limit
if(Distance(tx_site[0],rx_site) > 100){
fprintf(stdout,"\nToo far!\n");
fflush(stdout);
exit(1);
}
	

	if (sdf_path[0])
	{
		x=strlen(sdf_path);

		if (sdf_path[x-1]!='/' && x!=0)
		{
			sdf_path[x]='/';
			sdf_path[x+1]=0;
		}
	}

	fprintf(stdout,"%s",header);
	fflush(stdout);

	x=0;
	y=0;

	min_lat=90;
	max_lat=-90;

	min_lon=(int)floor(tx_site[0].lon);
	max_lon=(int)floor(tx_site[0].lon);

	
		txlat=(int)floor(tx_site[0].lat);
		txlon=(int)floor(tx_site[0].lon);

		if (txlat<min_lat)
			min_lat=txlat;

		if (txlat>max_lat)
			max_lat=txlat;

		if (LonDiff(txlon,min_lon)<0.0)
			min_lon=txlon;

		if (LonDiff(txlon,max_lon)>=0.0)
			max_lon=txlon;
	
	
		rxlat=(int)floor(rx_site.lat);
		rxlon=(int)floor(rx_site.lon);

		if (rxlat<min_lat)
			min_lat=rxlat;

		if (rxlat>max_lat)
			max_lat=rxlat;

		if (LonDiff(rxlon,min_lon)<0.0)
			min_lon=rxlon;

		if (LonDiff(rxlon,max_lon)>=0.0)
			max_lon=rxlon;
	
	/* Load the required SDF files */ 
		
	LoadTopoData(max_lon, min_lon, max_lat, min_lat, winfiles);

	
		//PlaceMarker(rx_site);

		if (height_plot)
		{
			/* Extract extension (if present)
			   from "height_file" */

			y=strlen(height_file);

			for (x=y-1; x>0 && height_file[x]!='.'; x--);

			if (x>0)  /* Extension found */
			{
				for (z=x+1; z<=y && (z-(x+1))<10; z++)
					ext[z-(x+1)]=tolower(height_file[z]);

				ext[z-(x+1)]=0;    /* Ensure an ending 0 */
				height_file[x]=0;  /* Chop off extension */
			}

			else
				strncpy(ext,"png\0",4);
		}


			//Plot it
			PlotPath(tx_site[0],rx_site,1);
			PathReport(tx_site[0],rx_site,height_file,0);
			//add .png extension to filename
			snprintf(string,250,"%s.%s%c",height_file,ext,0);
			GraphHeight(tx_site[0],rx_site,string,fresnel_plot,norm,width,height);
			
			double outvalue=dBm;
			char* outunit="dBm";
			
			if(output==1){
			outvalue=loss;
			outunit="dB";
			}
			if(output==3){
			outvalue=field_strength;
			outunit="dBuV/m";
			}
			fprintf(stdout,"%.2f %s\n",outvalue,outunit);
			fflush(stdout);
			
			//RETURN PATH
			//Plot it
			PlotPath(rx_site,tx_site[0],1);
			snprintf(return_height_file,250,"%s_R%c",height_file,0);
			PathReport(rx_site,tx_site[0],return_height_file,0);
			//add .png extension to filename_R
			snprintf(string,250,"%s_R.%s%c",height_file,ext,0);
			GraphHeight(rx_site,tx_site[0],string,fresnel_plot,norm,width,height);
			
			outvalue=dBm;
			outunit="dBm";
			if(output==1){
			outvalue=loss;
			outunit="dB";
			}
			if(output==3){
			outvalue=field_strength;
			outunit="dBuV/m";
			}
			fprintf(stdout,"%.2f %s\n",outvalue,outunit);
			fflush(stdout);
	
	return 0;
}

