double version=2.43;
/****************************************************************************\
*       Signal Server: Server optimised SPLAT! by Alex Farrant               *
******************************************************************************
*    SPLAT! Project started in 1997 by John A. Magliacane, KD2BD             *
*                                                                            *
******************************************************************************
*         Please consult the SPLAT! documentation for a complete list of     *
*         individuals who have contributed to this project.                  *
******************************************************************************
*                                                                            *
*  This program is free software; you can redistribute it and/or modify it   *
*  under the terms of the GNU General Public License as published by the     *
*  Free Software Foundation; either version 2 of the License or any later    *
*  version.                                                                  *
*                                                                            *
*  This program is distributed in the hope that it will useful, but WITHOUT  *
*  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or     *
*  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License     *
*  for more details.                                                         *
*                                                                            *
\****************************************************************************/

/*
REQUIRES GCC >= 4.7
90m mode
g++ -Wall -Ofast -s -lm itwom3.0.cpp models.cpp main.cpp -o signalserver
30m HD mode
g++ -Wall -Ofast -s -lm itwom3.0.cpp models.cpp main.cpp -DHD -o signalserverHD
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>


#ifdef HD
// HD mode, 30m res
#define MAXPAGES 9
#define ARRAYSIZE 14844
#define IPPD 3600
#else
// 90m mode (default)
#define MAXPAGES 64
#define ARRAYSIZE 76810
#define IPPD 1200
#endif

#define GAMMA 2.5

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
#define    METERS_PER_MILE 1609.344
#define METERS_PER_FOOT 0.3048
#define    KM_PER_MILE 1.609344
#define FOUR_THIRDS 1.3333333333333

char     string[255], sdf_path[255], udt_file[255], opened=0, gpsav=0, ss_name[16], dashes[80];

double    earthradius, max_range=0.0, forced_erp, dpp, ppd,
    fzone_clearance=0.6, forced_freq, clutter, lat,lon,txh,tercon,terdic,
        north,east,south,west,dBm,loss,field_strength;

int    min_north=90, max_north=-90, min_west=360, max_west=-1, ippd, mpi,
    max_elevation=-32768, min_elevation=32768, bzerror, contour_threshold,
        pred,pblue,pgreen,ter,multiplier=256,debug=0,loops=64,jgets=0, MAXRAD, hottest=10;

unsigned char got_elevation_pattern, got_azimuth_pattern, metric=0, dbm=0;



struct site {    double lat;
        double lon;
        float alt;
        char name[50];
        char filename[255];
        }     site;

struct path {    double lat[ARRAYSIZE];
        double lon[ARRAYSIZE];
        double elevation[ARRAYSIZE];
        double distance[ARRAYSIZE];
        int length;
        }    path;

struct dem {    int min_north;
        int max_north;
        int min_west;
        int max_west;
        int max_el;
        int min_el;
        short data[IPPD][IPPD];
        unsigned char mask[IPPD][IPPD];
        unsigned char signal[IPPD][IPPD];
           }    dem[MAXPAGES];

struct LR {    double eps_dielect;
        double sgm_conductivity;
        double eno_ns_surfref;
        double frq_mhz;
        double conf;
        double rel;
        double erp;
        int radio_climate;
        int pol;
        float antenna_pattern[361][1001];
          }    LR;

struct region { unsigned char color[128][3];
        int level[128];
        int levels;
          }    region;

double elev[ARRAYSIZE+10];

struct        site tx_site[2];

//ITWOM
void point_to_point(double elev[], double tht_m, double rht_m,
      double eps_dielect, double sgm_conductivity, double eno_ns_surfref,
      double frq_mhz, int radio_climate, int pol, double conf,
      double rel, double &dbloss, char *strmode, int &errnum);
//ITM      
void point_to_point_ITM(double elev[], double tht_m, double rht_m,
      double eps_dielect, double sgm_conductivity, double eno_ns_surfref,
      double frq_mhz, int radio_climate, int pol, double conf,
      double rel, double &dbloss, char *strmode, int &errnum);
      
double HATApathLoss(float f,float TxH, float RxH, float d, int mode);

double COST231pathLoss(float f,float TxH, float RxH, float d, int mode);

double FSPLpathLoss(float f, float d);

double SUIpathLoss(float f,float txH, float rxH, float d, int terrain);

double ECC33pathLoss(float f,float TxH, float RxH, float d, int mode);

double EricssonpathLoss(float f,float TxH, float RxH, float d, int mode);

double ked(double freq, double elev[], double rxh, double dkm);

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

    char    sign;
    int    degrees, minutes, seconds;
    double    a, b, c, d;

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

int PutMask(double lat, double lon, int value)
{
    /* Lines, text, markings, and coverage areas are stored in a
       mask that is combined with topology data when topographic
       maps are generated by ss.  This function sets and resets
       bits in the mask based on the latitude and longitude of the
       area pointed to. */

    int    x, y, indx;
    char    found;

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
        dem[indx].mask[x][y]=value;
        return ((int)dem[indx].mask[x][y]);
    }

    else
        return -1;
}

int OrMask(double lat, double lon, int value)
{
    /* Lines, text, markings, and coverage areas are stored in a
       mask that is combined with topology data when topographic
       maps are generated by ss.  This function sets bits in
       the mask based on the latitude and longitude of the area
       pointed to. */

    int    x, y, indx;
    char    found;

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

int GetMask(double lat, double lon)
{
    /* This function returns the mask bits based on the latitude
       and longitude given. */

    return (OrMask(lat,lon,0));
}

int PutSignal(double lat, double lon, unsigned char signal)
{
    /* This function writes a signal level (0-255)
       at the specified location for later recall. */
    char dotfile[255];
    snprintf(dotfile,80,"%s.dot%c",tx_site[0].filename,0);
    snprintf(dotfile,80,"%s.dot%c",tx_site[0].filename,0);
            
    int    x, y, indx;
    char    found;
    
    if(signal > hottest) // dBm, dBuV
        hottest=signal;
        
    //lookup x/y for this co-ord
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
        // Write values to file
        dem[indx].signal[x][y]=signal;
        
        return (dem[indx].signal[x][y]);
    }

    else
        return 0;
}

unsigned char GetSignal(double lat, double lon)
{
    /* This function reads the signal level (0-255) at the
       specified location that was previously written by the
       complimentary PutSignal() function. */

    int    x, y, indx;
    char    found;

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
        return (dem[indx].signal[x][y]);
    else
        return 0;
}

double GetElevation(struct site location)
{
    /* This function returns the elevation (in feet) of any location
       represented by the digital elevation model data in memory.
       Function returns -5000.0 for locations not found in memory. */

    char    found;
    int    x, y, indx;
    double    elevation;

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

    char    found;
    int    x, y, indx;
    

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

    double    lat1, lon1, lat2, lon2, distance;

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

    double    dest_lat, dest_lon, src_lat, src_lon,
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

    int    c;
    double    azimuth, distance, lat1, lon1, beta, den, num,
        lat2, lon2, total_distance, dx, dy, path_length,
        miles_per_sample, samples_per_radian=68755.0;
    struct    site tempsite;
    
    lat1=source.lat*DEG2RAD;
    lon1=source.lon*DEG2RAD;
    lat2=destination.lat*DEG2RAD;
    lon2=destination.lon*DEG2RAD;
    samples_per_radian=ppd*57.295833;
    azimuth=Azimuth(source,destination)*DEG2RAD;

    total_distance=Distance(source,destination);

    if (total_distance>(30.0/ppd))        
    {
        dx=samples_per_radian*acos(cos(lon1-lon2));
        dy=samples_per_radian*acos(cos(lat1-lat2));

        path_length=sqrt((dx*dx)+(dy*dy));        

        miles_per_sample=total_distance/path_length;    
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

    int    x;
    char    block=0;
    double    source_alt, destination_alt, cos_xmtr_angle,
        cos_test_angle, test_alt, elevation, distance,
        source_alt2, first_obstruction_angle=0.0;
    struct    path temp;

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

    double    seconds, bearing=0.0;
    char    string[20];
    int    a, b, length, degrees, minutes;

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

void LoadPAT(char *filename)
{
    /* This function reads and processes antenna pattern (.az
       and .el) files that correspond in name to previously
       loaded ss .lrp files.  */

    int    a, b, w, x, y, z, last_index, next_index, span;
    char    string[255], azfile[255], elfile[255], *pointer=NULL;
    float    az, xx, elevation, amplitude, rotation, valid1, valid2,
        delta, azimuth[361], azimuth_pattern[361], el_pattern[10001],
        elevation_pattern[361][1001], slant_angle[361], tilt,
        mechanical_tilt=0.0, tilt_azimuth, tilt_increment, sum;
    FILE    *fd=NULL;
    unsigned char read_count[10001];

    for (x=0; filename[x]!='.' && filename[x]!=0 && x<250; x++)
    {
        azfile[x]=filename[x];
        elfile[x]=filename[x];
    }

    azfile[x]='.';
    azfile[x+1]='a';
    azfile[x+2]='z';
    azfile[x+3]=0;

    elfile[x]='.';
    elfile[x+1]='e';
    elfile[x+2]='l';
    elfile[x+3]=0;

    rotation=0.0;

    got_azimuth_pattern=0;
    got_elevation_pattern=0;

    /* Load .az antenna pattern file */

    fd=fopen(azfile,"r");

    if (fd!=NULL)
    {
        /* Clear azimuth pattern array */

        for (x=0; x<=360; x++)
        {
            azimuth[x]=0.0;
            read_count[x]=0;
        }


        /* Read azimuth pattern rotation
           in degrees measured clockwise
           from true North. */

        if(fgets(string,254,fd) == NULL){
           //fprintf(stdout,"Azimuth read error\n");
           //exit(0);
        }
        pointer=strchr(string,';');

        if (pointer!=NULL)
            *pointer=0;

        sscanf(string,"%f",&rotation);


        /* Read azimuth (degrees) and corresponding
           normalized field radiation pattern amplitude
           (0.0 to 1.0) until EOF is reached. */

        if(fgets(string,254,fd) == NULL){
           //fprintf(stdout,"Azimuth read error\n");
           //exit(0);
        }
        pointer=strchr(string,';');

        if (pointer!=NULL)
            *pointer=0;

        sscanf(string,"%f %f",&az, &amplitude);

        do
        {
            x=(int)rintf(az);

            if (x>=0 && x<=360 && fd!=NULL)
            {
                azimuth[x]+=amplitude;
                read_count[x]++;
            }

            if(fgets(string,254,fd)== NULL){
                //fprintf(stdout,"Azimuth read error\n");
               // exit(0);
            }
            pointer=strchr(string,';');

            if (pointer!=NULL)
                *pointer=0;

            sscanf(string,"%f %f",&az, &amplitude);

        } while (feof(fd)==0);

        fclose(fd);


        /* Handle 0=360 degree ambiguity */

        if ((read_count[0]==0) && (read_count[360]!=0))
        {
            read_count[0]=read_count[360];
            azimuth[0]=azimuth[360];
        }

        if ((read_count[0]!=0) && (read_count[360]==0))
        {
            read_count[360]=read_count[0];
            azimuth[360]=azimuth[0];
        }

        /* Average pattern values in case more than
            one was read for each degree of azimuth. */

        for (x=0; x<=360; x++)
        {
            if (read_count[x]>1)
                azimuth[x]/=(float)read_count[x];
        }

        /* Interpolate missing azimuths
           to completely fill the array */

        last_index=-1;
        next_index=-1;

        for (x=0; x<=360; x++)
        {
            if (read_count[x]!=0)
            {
                if (last_index==-1)
                    last_index=x;
                else
                    next_index=x;
            }

            if (last_index!=-1 && next_index!=-1)
            {
                valid1=azimuth[last_index];
                valid2=azimuth[next_index];

                span=next_index-last_index;
                delta=(valid2-valid1)/(float)span;

                for (y=last_index+1; y<next_index; y++)
                    azimuth[y]=azimuth[y-1]+delta;

                last_index=y;
                next_index=-1;
            }
        }

        /* Perform azimuth pattern rotation
           and load azimuth_pattern[361] with
           azimuth pattern data in its final form. */

        for (x=0; x<360; x++)
        {
            y=x+(int)rintf(rotation);

            if (y>=360)
                y-=360;

            azimuth_pattern[y]=azimuth[x];
        }

        azimuth_pattern[360]=azimuth_pattern[0];

        got_azimuth_pattern=255;
    }

    /* Read and process .el file */

    fd=fopen(elfile,"r");

    if (fd!=NULL)
    {
        for (x=0; x<=10000; x++)
        {
            el_pattern[x]=0.0;
            read_count[x]=0;
        }

        /* Read mechanical tilt (degrees) and
           tilt azimuth in degrees measured
           clockwise from true North. */

        if(fgets(string,254,fd)==NULL){
            //fprintf(stdout,"Tilt read error\n");
            //exit(0);
        }
        pointer=strchr(string,';');

        if (pointer!=NULL)
            *pointer=0;

        sscanf(string,"%f %f",&mechanical_tilt, &tilt_azimuth);

        /* Read elevation (degrees) and corresponding
           normalized field radiation pattern amplitude
           (0.0 to 1.0) until EOF is reached. */

        if(fgets(string,254,fd)==NULL){
           //fprintf(stdout,"Ant elevation read error\n");
           //exit(0);
        }
        pointer=strchr(string,';');

        if (pointer!=NULL)
            *pointer=0;

        sscanf(string,"%f %f", &elevation, &amplitude);

        while (feof(fd)==0)
        {
            /* Read in normalized radiated field values
               for every 0.01 degrees of elevation between
               -10.0 and +90.0 degrees */

            x=(int)rintf(100.0*(elevation+10.0));

            if (x>=0 && x<=10000)
            {
                el_pattern[x]+=amplitude;
                read_count[x]++;
            }

            if(fgets(string,254,fd)!=NULL){
                pointer=strchr(string,';');
            }
            if (pointer!=NULL)
                *pointer=0;

            sscanf(string,"%f %f", &elevation, &amplitude);
        }

        fclose(fd);

        /* Average the field values in case more than
           one was read for each 0.01 degrees of elevation. */

        for (x=0; x<=10000; x++)
        {
            if (read_count[x]>1)
                el_pattern[x]/=(float)read_count[x];
        }

        /* Interpolate between missing elevations (if
           any) to completely fill the array and provide
           radiated field values for every 0.01 degrees of
           elevation. */

        last_index=-1;
        next_index=-1;

        for (x=0; x<=10000; x++)
        {
            if (read_count[x]!=0)
            {
                if (last_index==-1)
                    last_index=x;
                else
                    next_index=x;
            }

            if (last_index!=-1 && next_index!=-1)
            {
                valid1=el_pattern[last_index];
                valid2=el_pattern[next_index];

                span=next_index-last_index;
                delta=(valid2-valid1)/(float)span;

                for (y=last_index+1; y<next_index; y++)
                    el_pattern[y]=el_pattern[y-1]+delta;

                last_index=y;
                next_index=-1;
            }
        }

        /* Fill slant_angle[] array with offset angles based
           on the antenna's mechanical beam tilt (if any)
           and tilt direction (azimuth). */

        if (mechanical_tilt==0.0)
        {
            for (x=0; x<=360; x++)
                slant_angle[x]=0.0;
        }

        else
        {
            tilt_increment=mechanical_tilt/90.0;

            for (x=0; x<=360; x++)
            {
                xx=(float)x;
                y=(int)rintf(tilt_azimuth+xx);

                while (y>=360)
                    y-=360;

                while (y<0)
                    y+=360;

                if (x<=180)
                    slant_angle[y]=-(tilt_increment*(90.0-xx));

                if (x>180)
                    slant_angle[y]=-(tilt_increment*(xx-270.0));
            }
        }

        slant_angle[360]=slant_angle[0];   /* 360 degree wrap-around */

        for (w=0; w<=360; w++)
        {
            tilt=slant_angle[w];

            /** Convert tilt angle to
                an array index offset **/

            y=(int)rintf(100.0*tilt);

            /* Copy shifted el_pattern[10001] field
               values into elevation_pattern[361][1001]
               at the corresponding azimuth, downsampling
               (averaging) along the way in chunks of 10. */

            for (x=y, z=0; z<=1000; x+=10, z++)
            {
                for (sum=0.0, a=0; a<10; a++)
                {
                    b=a+x;

                    if (b>=0 && b<=10000)
                        sum+=el_pattern[b];
                    if (b<0)
                        sum+=el_pattern[0];
                    if (b>10000)
                        sum+=el_pattern[10000];
                }

                elevation_pattern[w][z]=sum/10.0;
            }
        }

        got_elevation_pattern=255;
    }

    for (x=0; x<=360; x++)
    {
        for (y=0; y<=1000; y++)
        {
            if (got_elevation_pattern)
                elevation=elevation_pattern[x][y];
            else
                elevation=1.0;

            if (got_azimuth_pattern)
                az=azimuth_pattern[x];
            else
                az=1.0;

            LR.antenna_pattern[x][y]=az*elevation;
        }
    }
}

int LoadSDF_SDF(char *name, int winfiles)
{
    /* This function reads uncompressed ss Data Files (.sdf)
       containing digital elevation model data into memory.
       Elevation data, maximum and minimum elevations, and
       quadrangle limits are stored in the first available
       dem[] structure. */

    int    x, y, data, indx, minlat, minlon, maxlat, maxlon,j;
    char    found, free_page=0, line[20], jline[20], sdf_file[255],
        path_plus_name[255], *junk=NULL;
        
        
    FILE    *fd;

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
            /* Next, try loading SDF file from path specified
               in $HOME/.ss_path file or by -d argument */

            strncpy(path_plus_name,sdf_path,255);
            strncat(path_plus_name,sdf_file,255);
            fd=fopen(path_plus_name,"rb");
        }

        if (fd!=NULL)
        {
                    if (debug==1){
            fprintf(stdout,"Loading \"%s\" into page %d...",path_plus_name,indx+1);
            fflush(stdout);
                    }

            if(fgets(line,19,fd)!=NULL){
                sscanf(line,"%d",&dem[indx].max_west);
            }
            
            if(fgets(line,19,fd)!=NULL){
                sscanf(line,"%d",&dem[indx].min_north);
            }
            
            if(fgets(line,19,fd)!=NULL){
                sscanf(line,"%d",&dem[indx].min_west);
            }
            
            if(fgets(line,19,fd)!=NULL){
                sscanf(line,"%d",&dem[indx].max_north);
            }
            /*
            Here X lines of DEM will be read until IPPD is reached.
            Each .sdf tile contains 1200x1200 = 1.44M 'points'
            Each point is sampled for 1200 resolution!
            */
            for (x=0; x<ippd; x++)
                        {
                for (y=0; y<ippd; y++)
                {
                            
                    for (j=0; j<jgets; j++)
                    {
                        junk=fgets(jline,19,fd);   
                    }
                                    
                    if(fgets(line,19,fd)!=NULL){                                                                        
                    data=atoi(line); 
                    }
                                        
                                        dem[indx].data[x][y]=data;
                    dem[indx].signal[x][y]=0;
                    dem[indx].mask[x][y]=0;

                    if (data>dem[indx].max_el)
                        dem[indx].max_el=data;

                    if (data<dem[indx].min_el)
                        dem[indx].min_el=data;
                                       
                } 
                                
                                if (ippd==600){ 
                                for (j=0; j<IPPD; j++)
                                      {
                                        junk=fgets(jline,19,fd);   
                                        }
                                }
                                if (ippd==300){ 
                                for (j=0; j<IPPD; j++)
                                      {
                                        junk=fgets(jline,19,fd); 
                                        junk=fgets(jline,19,fd);
                                        junk=fgets(jline,19,fd);
                                                                           
                                        }
                                }
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

    int    x, y, indx, minlat, minlon, maxlat, maxlon;
    char    found, free_page=0;
    int    return_value=-1;

    return_value=LoadSDF_SDF(name, winfiles);


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
                    if(debug==1){
            fprintf(stdout,"Region  \"%s\" assumed as sea-level into page %d...",name,indx+1);
            fflush(stdout);
                    }

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
                       
            return_value=1;
        }
    }

    return return_value;
}

void PlotLOSPath(struct site source, struct site destination, char mask_value, FILE *fd)
{
    /* This function analyzes the path between the source and
       destination locations.  It determines which points along
       the path have line-of-sight visibility to the source.
       Points along with path having line-of-sight visibility
       to the source at an AGL altitude equal to that of the
       destination location are stored by setting bit 1 in the
       mask[][] array, which are displayed in green when PPM
       maps are later generated by ss. */

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

void PlotPropPath(struct site source, struct site destination, unsigned char mask_value, FILE *fd, int propmodel, int knifeedge, int pmenv)
{

    int    x, y, ifs, ofs, errnum;
    char    block=0, strmode[100];
    double    loss, azimuth, pattern=0.0,
        xmtr_alt, dest_alt, xmtr_alt2, dest_alt2,
        cos_rcvr_angle, cos_test_angle=0.0, test_alt,
        elevation=0.0, distance=0.0, radius=0.0, four_thirds_earth,
        field_strength=0.0, rxp, dBm, txelev, dkm, diffloss;
    struct    site temp;

    radius = Distance(source,destination);

    ReadPath(source,destination);

    four_thirds_earth=FOUR_THIRDS*EARTHRADIUS;


    for (x=1; x<path.length-1; x++)
        elev[x+2]=(path.elevation[x]==0.0?path.elevation[x]*METERS_PER_FOOT:(clutter+path.elevation[x])*METERS_PER_FOOT);

    /* Copy ending points without clutter */

    elev[2]=path.elevation[0]*METERS_PER_FOOT;
    txelev=elev[2]+(source.alt*METERS_PER_FOOT);
    
    elev[path.length+1]=path.elevation[path.length-1]*METERS_PER_FOOT;
    
    /* Since the only energy the Longley-Rice model considers
       reaching the destination is based on what is scattered
       or deflected from the first obstruction along the path,
       we first need to find the location and elevation angle
       of that first obstruction (if it exists).  This is done
       using a 4/3rds Earth radius to match the model used by
       Longley-Rice.  This information is required for properly
       integrating the antenna's elevation pattern into the
       calculation for overall path loss. */

    for (y=2; (y<(path.length-1) && path.distance[y]<=max_range); y++)
    {
        /* Process this point only if it
           has not already been processed. */
        
    
        if ((GetMask(path.lat[y],path.lon[y])&248)!=(mask_value<<3))
        {
            distance=5280.0*path.distance[y];
            xmtr_alt=four_thirds_earth+source.alt+path.elevation[0];
            dest_alt=four_thirds_earth+destination.alt+path.elevation[y];
            dest_alt2=dest_alt*dest_alt;
            xmtr_alt2=xmtr_alt*xmtr_alt;

            /* Calculate the cosine of the elevation of
               the receiver as seen by the transmitter. */

            cos_rcvr_angle=((xmtr_alt2)+(distance*distance)-(dest_alt2))/(2.0*xmtr_alt*distance);

            if (cos_rcvr_angle>1.0)
                cos_rcvr_angle=1.0;

            if (cos_rcvr_angle<-1.0)
                cos_rcvr_angle=-1.0;

            if (got_elevation_pattern || fd!=NULL)
            {
                /* Determine the elevation angle to the first obstruction
                   along the path IF elevation pattern data is available
                   or an output (.ano) file has been designated. */

                for (x=2, block=0; (x<y && block==0); x++)
                {
                    distance=5280.0*path.distance[x];

                    test_alt=four_thirds_earth+(path.elevation[x]==0.0?path.elevation[x]:path.elevation[x]+clutter);

                    /* Calculate the cosine of the elevation
                       angle of the terrain (test point)
                       as seen by the transmitter. */

                    cos_test_angle=((xmtr_alt2)+(distance*distance)-(test_alt*test_alt))/(2.0*xmtr_alt*distance);

                    if (cos_test_angle>1.0)
                        cos_test_angle=1.0;

                    if (cos_test_angle<-1.0)
                        cos_test_angle=-1.0;

                    /* Compare these two angles to determine if
                       an obstruction exists.  Since we're comparing
                       the cosines of these angles rather than
                       the angles themselves, the sense of the
                       following "if" statement is reversed from
                         what it would be if the angles themselves
                       were compared. */

                    if (cos_rcvr_angle>=cos_test_angle)
                        block=1;
                }

                if (block)
                    elevation=((acos(cos_test_angle))/DEG2RAD)-90.0;
                else
                    elevation=((acos(cos_rcvr_angle))/DEG2RAD)-90.0;
            }

            /* Determine attenuation for each point along the
               path using a prop model starting at y=2 (number_of_points = 1), the
               shortest distance terrain can play a role in
               path loss. */

            elev[0]=y-1;  /* (number of points - 1) */

            /* Distance between elevation samples */

            elev[1]=METERS_PER_MILE*(path.distance[y]-path.distance[y-1]);

            if(path.elevation[y] < 1){
            path.elevation[y]=1;
            }
            
            dkm=(elev[1]*elev[0])/1000; // km
    
            switch (propmodel)
            {
              case 1:
                // Longley Rice ITM
                 point_to_point_ITM(elev,source.alt*METERS_PER_FOOT,
                    destination.alt*METERS_PER_FOOT, LR.eps_dielect,
                    LR.sgm_conductivity, LR.eno_ns_surfref, LR.frq_mhz,
                    LR.radio_climate, LR.pol, LR.conf, LR.rel, loss,
                    strmode, errnum);
                 break;
              case 3:
                //HATA 1, 2 & 3
                loss=HATApathLoss(LR.frq_mhz,txelev,path.elevation[y]+(destination.alt*METERS_PER_FOOT),dkm, pmenv);
                break;
              case 4:
                // COST231-HATA
				loss=ECC33pathLoss(LR.frq_mhz,txelev,path.elevation[y]+(destination.alt*METERS_PER_FOOT),dkm, pmenv);
                break;
              case 5:
                // SUI
                loss=SUIpathLoss(LR.frq_mhz,txelev,path.elevation[y]+(destination.alt*METERS_PER_FOOT),dkm, pmenv);
                break;
              case 6:
                loss=COST231pathLoss(LR.frq_mhz,txelev,path.elevation[y]+(destination.alt*METERS_PER_FOOT),dkm, pmenv);
                break;
              case 7:
			  // ITU-R P.525 Free space path loss
                loss=FSPLpathLoss(LR.frq_mhz,dkm);
                break;				
              case 8:
                 // ITWOM 3.0 
                 point_to_point(elev,source.alt*METERS_PER_FOOT,
                    destination.alt*METERS_PER_FOOT, LR.eps_dielect,
                    LR.sgm_conductivity, LR.eno_ns_surfref, LR.frq_mhz,
                    LR.radio_climate, LR.pol, LR.conf, LR.rel, loss,
                    strmode, errnum);
                 break;
              case 9:
			   // Ericsson
                loss=EricssonpathLoss(LR.frq_mhz,txelev,path.elevation[y]+(destination.alt*METERS_PER_FOOT),dkm, pmenv);
                 break;

                
              default:
                 point_to_point_ITM(elev,source.alt*METERS_PER_FOOT,
                    destination.alt*METERS_PER_FOOT, LR.eps_dielect,
                    LR.sgm_conductivity, LR.eno_ns_surfref, LR.frq_mhz,
                    LR.radio_climate, LR.pol, LR.conf, LR.rel, loss,
                    strmode, errnum);
  
            }
            
    
            if(knifeedge==1){
            diffloss = ked(LR.frq_mhz,elev,destination.alt*METERS_PER_FOOT,dkm); 
            loss+=(diffloss); // ;)
            }
            
            //Key stage. Link dB for p2p is returned as 'loss'.
        
            temp.lat=path.lat[y];
            temp.lon=path.lon[y];

            azimuth=(Azimuth(source,temp));

            if (fd!=NULL)
                fprintf(fd,"%.7f, %.7f, %.3f, %.3f, ",path.lat[y], path.lon[y], azimuth, elevation);

            /* If ERP==0, write path loss to alphanumeric
               output file.  Otherwise, write field strength
               or received power level (below), as appropriate. */

            if (fd!=NULL && LR.erp==0.0)
                fprintf(fd,"%.2f",loss);

            /* Integrate the antenna's radiation
               pattern into the overall path loss. */

            x=(int)rint(10.0*(10.0-elevation));

            if (x>=0 && x<=1000)
            {
                azimuth=rint(azimuth);

                pattern=(double)LR.antenna_pattern[(int)azimuth][x];

                if (pattern!=0.0)
                {
                    pattern=20.0*log10(pattern);
                    loss-=pattern;
                }
            }

            if (LR.erp!=0.0)
            {
                if (dbm)
                {
                    /* dBm is based on EIRP (ERP + 2.14) */

                    rxp=LR.erp/(pow(10.0,(loss-2.14)/10.0));

                    dBm=10.0*(log10(rxp*1000.0));

                    if (fd!=NULL)
                        fprintf(fd,"%.3f",dBm);

                    /* Scale roughly between 0 and 255 */

                    ifs=200+(int)rint(dBm);

                    if (ifs<0)
                        ifs=0;

                    if (ifs>255)
                        ifs=255;

                    ofs=GetSignal(path.lat[y],path.lon[y]);

                    if (ofs>ifs)
                        ifs=ofs;

                    
                    PutSignal(path.lat[y],path.lon[y],(unsigned char)ifs);
                    
                }

                else
                {
                    field_strength=(139.4+(20.0*log10(LR.frq_mhz))-loss)+(10.0*log10(LR.erp/1000.0));

                    ifs=100+(int)rint(field_strength);

                    if (ifs<0)
                        ifs=0;

                    if (ifs>255)
                        ifs=255;

                    ofs=GetSignal(path.lat[y],path.lon[y]);

                    if (ofs>ifs)
                        ifs=ofs;

                    PutSignal(path.lat[y],path.lon[y],(unsigned char)ifs);

                    if (fd!=NULL)
                        fprintf(fd,"%.3f",field_strength);
                }
            }

            else
            {
                if (loss>255)
                    ifs=255;
                else
                    ifs=(int)rint(loss);

                ofs=GetSignal(path.lat[y],path.lon[y]);

                if (ofs<ifs && ofs!=0)
                    ifs=ofs;
            
                 PutSignal(path.lat[y],path.lon[y],(unsigned char)ifs);
            }

            if (fd!=NULL)
            {
                if (block)
                    fprintf(fd," *");

                fprintf(fd,"\n");
            }

            /* Mark this point as having been analyzed */

            
            PutMask(path.lat[y],path.lon[y],(GetMask(path.lat[y],path.lon[y])&7)+(mask_value<<3));
        }
    }
    
}

void PlotLOSMap(struct site source, double altitude, char *plo_filename)
{
    /* This function performs a 360 degree sweep around the
       transmitter site (source location), and plots the
       line-of-sight coverage of the transmitter on the ss
       generated topographic map based on a receiver located
       at the specified altitude (in feet AGL).  Results
       are stored in memory, and written out in the form
       of a topographic map when the WritePPM() function
       is later invoked. */

    int y, z;
    struct site edge;
    unsigned char x;
    double lat, lon, minwest, maxnorth, th;
    static unsigned char mask_value=1;
    FILE *fd=NULL;    

    if (plo_filename[0]!=0)
        fd=fopen(plo_filename,"wb");

    if (fd!=NULL)
    {
        fprintf(fd,"%d, %d\t; max_west, min_west\n%d, %d\t; max_north, min_north\n",max_west, min_west, max_north, min_north);
    }
    
    th=ppd/loops;

    z=(int)(th*ReduceAngle(max_west-min_west));

    minwest=dpp+(double)min_west;
    maxnorth=(double)max_north-dpp;


    for (lon=minwest, x=0, y=0; (LonDiff(lon,(double)max_west)<=0.0); y++, lon=minwest+(dpp*(double)y))
    {
        if (lon>=360.0)
            lon-=360.0;

        edge.lat=max_north;
        edge.lon=lon;
        edge.alt=altitude;

        PlotLOSPath(source,edge,mask_value,fd);
    }

    

    z=(int)(th*(double)(max_north-min_north));

    for (lat=maxnorth, x=0, y=0; lat>=(double)min_north; y++, lat=maxnorth-(dpp*(double)y))
    {
        edge.lat=lat;
        edge.lon=min_west;
        edge.alt=altitude;

        PlotLOSPath(source,edge,mask_value,fd);
    
    }



    z=(int)(th*ReduceAngle(max_west-min_west));

    for (lon=minwest, x=0, y=0; (LonDiff(lon,(double)max_west)<=0.0); y++, lon=minwest+(dpp*(double)y))
    {
        if (lon>=360.0)
            lon-=360.0;

        edge.lat=min_north;
        edge.lon=lon;
        edge.alt=altitude;

        PlotLOSPath(source,edge,mask_value,fd);

    }

    
    z=(int)(th*(double)(max_north-min_north));

    for (lat=(double)min_north, x=0, y=0; lat<(double)max_north; y++, lat=(double)min_north+(dpp*(double)y))
    {
        edge.lat=lat;
        edge.lon=max_west;
        edge.alt=altitude;

        PlotLOSPath(source,edge,mask_value,fd);
        

    }

    
    switch (mask_value)
    {
        case 1:
            mask_value=8;
            break;

        case 8:
            mask_value=16;
            break;

        case 16:
            mask_value=32;
    }
}

void PlotPropagation(struct site source, double altitude, char *plo_filename, int propmodel, int knifeedge, int haf, int pmenv)
{
    int y, z, count;
    struct site edge;
    double lat, lon, minwest, maxnorth, th;
    unsigned char x;
    static unsigned char mask_value=1;
    FILE *fd=NULL;

    minwest=dpp+(double)min_west;
    maxnorth=(double)max_north-dpp;

    count=0;
   

    if (LR.erp==0.0 && debug)
        fprintf(stdout,"path loss");
    else
    {
            if(debug){
        if (dbm)
            fprintf(stdout,"signal power level");
        else
            fprintf(stdout,"field strength");
            }
    }
    if (debug){
		fprintf(stdout," contours of \"%s\"\nout to a radius of %.2f %s with Rx antenna(s) at %.2f %s AGL\n",source.name,metric?max_range*KM_PER_MILE:max_range,metric?"kilometers":"miles",metric?altitude*METERS_PER_FOOT:altitude,metric?"meters":"feet");
    }

    if (clutter>0.0 && debug)
        fprintf(stdout,"\nand %.2f %s of ground clutter",metric?clutter*METERS_PER_FOOT:clutter,metric?"meters":"feet");
  
    if(debug){
        fprintf(stdout,"...\n\n 0%c to  25%c ",37,37);
        fflush(stdout);
    }
    
    if (plo_filename[0]!=0)
        fd=fopen(plo_filename,"wb");

    if (fd!=NULL)
    {
        fprintf(fd,"%d, %d\t; max_west, min_west\n%d, %d\t; max_north, min_north\n",max_west, min_west, max_north, min_north);
    }

    th=ppd/loops;

    // Four sections start here
    
    //S1
    if(haf==0 || haf==1){
    z=(int)(th*ReduceAngle(max_west-min_west));

    for (lon=minwest, x=0, y=0; (LonDiff(lon,(double)max_west)<=0.0); y++, lon=minwest+(dpp*(double)y))
    {
        if (lon>=360.0)
            lon-=360.0;

        edge.lat=max_north;
        edge.lon=lon;
        edge.alt=altitude;

        PlotPropPath(source,edge,mask_value,fd,propmodel,knifeedge,pmenv);
        count++;

        if (count==z) 
        {
            count=0;

            if (x==3)
                x=0;
            else
                x++;
        }
    }

    }
    
    //S2
    if(haf==0 || haf==1){
    count=0;
    if(debug){
        fprintf(stdout,"\n25%c to  50%c ",37,37);
        fflush(stdout);
    }
    
    z=(int)(th*(double)(max_north-min_north));

    for (lat=maxnorth, x=0, y=0; lat>=(double)min_north; y++, lat=maxnorth-(dpp*(double)y))
    {
        edge.lat=lat;
        edge.lon=min_west;
        edge.alt=altitude;

        PlotPropPath(source,edge,mask_value,fd,propmodel,knifeedge,pmenv);
        count++;

        if (count==z) 
        {
            count=0;

            if (x==3)
                x=0;
            else
                x++;
        }
    }

    }
    //S3
    if(haf==0 || haf==2){
    count=0;
    if(debug){
        fprintf(stdout,"\n50%c to  75%c ",37,37);
        fflush(stdout);
    }
    
    z=(int)(th*ReduceAngle(max_west-min_west));

    for (lon=minwest, x=0, y=0; (LonDiff(lon,(double)max_west)<=0.0); y++, lon=minwest+(dpp*(double)y))
    {
        if (lon>=360.0)
            lon-=360.0;

        edge.lat=min_north;
        edge.lon=lon;
        edge.alt=altitude;

        PlotPropPath(source,edge,mask_value,fd,propmodel,knifeedge,pmenv);
        count++;
                if (count==z) 
        {
            count=0;

            if (x==3)
                x=0;
            else
                x++;
        }
        
    }

    }
    //S4
    if(haf==0 || haf==2){
    count=0;
    if(debug){
        fprintf(stdout,"\n75%c to 100%c ",37,37);
        fflush(stdout);
    }
    z=(int)(th*(double)(max_north-min_north));

    for (lat=(double)min_north, x=0, y=0; lat<(double)max_north; y++, lat=(double)min_north+(dpp*(double)y))
    {
        edge.lat=lat;
        edge.lon=max_west;
        edge.alt=altitude;

        PlotPropPath(source,edge,mask_value,fd,propmodel,knifeedge,pmenv);
        count++;

        if (count==z) 
        {
        
            count=0;

            if (x==3)
                x=0;
            else
                x++;
        }
    }

    } //S4 
    
    if (fd!=NULL)
        fclose(fd);

    if (mask_value<30)
        mask_value++;
}

void LoadSignalColors(struct site xmtr)
{
    int x, y, ok, val[4];
    char filename[255], string[80], *pointer=NULL, *s=NULL;
    FILE *fd=NULL;

    for (x=0; xmtr.filename[x]!='.' && xmtr.filename[x]!=0 && x<250; x++)
        filename[x]=xmtr.filename[x];

    filename[x]='.';
    filename[x+1]='s';
    filename[x+2]='c';
    filename[x+3]='f';
    filename[x+4]=0;

    /* Default values */

    region.level[0]=128;
    region.color[0][0]=255;
    region.color[0][1]=0;
    region.color[0][2]=0;

    region.level[1]=118;
    region.color[1][0]=255;
    region.color[1][1]=165;
    region.color[1][2]=0;

    region.level[2]=108;
    region.color[2][0]=255;
    region.color[2][1]=206;
    region.color[2][2]=0;

    region.level[3]=98;
    region.color[3][0]=255;
    region.color[3][1]=255;
    region.color[3][2]=0;

    region.level[4]=88;
    region.color[4][0]=184;
    region.color[4][1]=255;
    region.color[4][2]=0;

    region.level[5]=78;
    region.color[5][0]=0;
    region.color[5][1]=255;
    region.color[5][2]=0;

    region.level[6]=68;
    region.color[6][0]=0;
    region.color[6][1]=208;
    region.color[6][2]=0;

    region.level[7]=58;
    region.color[7][0]=0;
    region.color[7][1]=196;
    region.color[7][2]=196;

    region.level[8]=48;
    region.color[8][0]=0;
    region.color[8][1]=148;
    region.color[8][2]=255;

    region.level[9]=38;
    region.color[9][0]=80;
    region.color[9][1]=80;
    region.color[9][2]=255;

    region.level[10]=28;
    region.color[10][0]=0;
    region.color[10][1]=38;
    region.color[10][2]=255;

    region.level[11]=18;
    region.color[11][0]=142;
    region.color[11][1]=63;
    region.color[11][2]=255;

    region.level[12]=8;
    region.color[12][0]=140;
    region.color[12][1]=0;
    region.color[12][2]=128;

    region.levels=13;

    fd=fopen(filename,"r");

    if (fd==NULL)
        fd=fopen(filename,"r");

    if (fd==NULL)
    {
        fd=fopen(filename,"w");

        
        for (x=0; x<region.levels; x++)
            fprintf(fd,"%3d: %3d, %3d, %3d\n",region.level[x], region.color[x][0], region.color[x][1], region.color[x][2]);

        fclose(fd);
    }

    else
    {
        x=0;
        s=fgets(string,80,fd);

        while (x<128 && feof(fd)==0)
        {
            pointer=strchr(string,';');

            if (pointer!=NULL)
                *pointer=0;

            ok=sscanf(string,"%d: %d, %d, %d", &val[0], &val[1], &val[2], &val[3]);

            if (ok==4)
            {
                for (y=0; y<4; y++)
                {
                    if (val[y]>255)
                        val[y]=255;

                    if (val[y]<0)
                        val[y]=0;
                }

                region.level[x]=val[0];
                region.color[x][0]=val[1];
                region.color[x][1]=val[2];
                region.color[x][2]=val[3];
                x++;
            }

            s=fgets(string,80,fd);
        }

        fclose(fd);
        region.levels=x;
    }
}

void LoadLossColors(struct site xmtr)
{
    int x, y, ok, val[4];
    char filename[255], string[80], *pointer=NULL, *s=NULL;
    FILE *fd=NULL;

    for (x=0; xmtr.filename[x]!='.' && xmtr.filename[x]!=0 && x<250; x++)
        filename[x]=xmtr.filename[x];

    filename[x]='.';
    filename[x+1]='l';
    filename[x+2]='c';
    filename[x+3]='f';
    filename[x+4]=0;

    /* Default values */

    region.level[0]=80;
    region.color[0][0]=255;
    region.color[0][1]=0;
    region.color[0][2]=0;

    region.level[1]=90;
    region.color[1][0]=255;
    region.color[1][1]=128;
    region.color[1][2]=0;

    region.level[2]=100;
    region.color[2][0]=255;
    region.color[2][1]=165;
    region.color[2][2]=0;

    region.level[3]=110;
    region.color[3][0]=255;
    region.color[3][1]=206;
    region.color[3][2]=0;

    region.level[4]=120;
    region.color[4][0]=255;
    region.color[4][1]=255;
    region.color[4][2]=0;

    region.level[5]=130;
    region.color[5][0]=184;
    region.color[5][1]=255;
    region.color[5][2]=0;

    region.level[6]=140;
    region.color[6][0]=0;
    region.color[6][1]=255;
    region.color[6][2]=0;

    region.level[7]=150;
    region.color[7][0]=0;
    region.color[7][1]=208;
    region.color[7][2]=0;

    region.level[8]=160;
    region.color[8][0]=0;
    region.color[8][1]=196;
    region.color[8][2]=196;

    region.level[9]=170;
    region.color[9][0]=0;
    region.color[9][1]=148;
    region.color[9][2]=255;

    region.level[10]=180;
    region.color[10][0]=80;
    region.color[10][1]=80;
    region.color[10][2]=255;

    region.level[11]=190;
    region.color[11][0]=0;
    region.color[11][1]=38;
    region.color[11][2]=255;

    region.level[12]=200;
    region.color[12][0]=142;
    region.color[12][1]=63;
    region.color[12][2]=255;

    region.level[13]=210;
    region.color[13][0]=196;
    region.color[13][1]=54;
    region.color[13][2]=255;

    region.level[14]=220;
    region.color[14][0]=255;
    region.color[14][1]=0;
    region.color[14][2]=255;

    region.level[15]=230;
    region.color[15][0]=255;
    region.color[15][1]=194;
    region.color[15][2]=204;

    region.levels=16;

    fd=fopen(filename,"r");

    if (fd==NULL)
        fd=fopen(filename,"r");

    if (fd==NULL)
    {
        fd=fopen(filename,"w");

        

        for (x=0; x<region.levels; x++)
            fprintf(fd,"%3d: %3d, %3d, %3d\n",region.level[x], region.color[x][0], region.color[x][1], region.color[x][2]);

        fclose(fd);
    }

    else
    {
        x=0;
        s=fgets(string,80,fd);

        while (x<128 && feof(fd)==0)
        {
            pointer=strchr(string,';');

            if (pointer!=NULL)
                *pointer=0;

            ok=sscanf(string,"%d: %d, %d, %d", &val[0], &val[1], &val[2], &val[3]);

            if (ok==4)
            {
                for (y=0; y<4; y++)
                {
                    if (val[y]>255)
                        val[y]=255;

                    if (val[y]<0)
                        val[y]=0;
                }

                region.level[x]=val[0];
                region.color[x][0]=val[1];
                region.color[x][1]=val[2];
                region.color[x][2]=val[3];
                x++;
            }

            s=fgets(string,80,fd);
        }

        fclose(fd);
        region.levels=x;
    }
}

void LoadDBMColors(struct site xmtr)
{
    int x, y, ok, val[4];
    char filename[255], string[80], *pointer=NULL, *s=NULL;
    FILE *fd=NULL;

    for (x=0; xmtr.filename[x]!='.' && xmtr.filename[x]!=0 && x<250; x++)
        filename[x]=xmtr.filename[x];

    filename[x]='.';
    filename[x+1]='d';
    filename[x+2]='c';
    filename[x+3]='f';
    filename[x+4]=0;

    /* Default values */

    region.level[0]=0;
    region.color[0][0]=255;
    region.color[0][1]=0;
    region.color[0][2]=0;

    region.level[1]=-10;
    region.color[1][0]=255;
    region.color[1][1]=128;
    region.color[1][2]=0;

    region.level[2]=-20;
    region.color[2][0]=255;
    region.color[2][1]=165;
    region.color[2][2]=0;

    region.level[3]=-30;
    region.color[3][0]=255;
    region.color[3][1]=206;
    region.color[3][2]=0;

    region.level[4]=-40;
    region.color[4][0]=255;
    region.color[4][1]=255;
    region.color[4][2]=0;

    region.level[5]=-50;
    region.color[5][0]=184;
    region.color[5][1]=255;
    region.color[5][2]=0;

    region.level[6]=-60;
    region.color[6][0]=0;
    region.color[6][1]=255;
    region.color[6][2]=0;

    region.level[7]=-70;
    region.color[7][0]=0;
    region.color[7][1]=208;
    region.color[7][2]=0;

    region.level[8]=-80;
    region.color[8][0]=0;
    region.color[8][1]=196;
    region.color[8][2]=196;

    region.level[9]=-90;
    region.color[9][0]=0;
    region.color[9][1]=148;
    region.color[9][2]=255;

    region.level[10]=-100;
    region.color[10][0]=80;
    region.color[10][1]=80;
    region.color[10][2]=255;

    region.level[11]=-110;
    region.color[11][0]=0;
    region.color[11][1]=38;
    region.color[11][2]=255;

    region.level[12]=-120;
    region.color[12][0]=142;
    region.color[12][1]=63;
    region.color[12][2]=255;

    region.level[13]=-130;
    region.color[13][0]=196;
    region.color[13][1]=54;
    region.color[13][2]=255;

    region.level[14]=-140;
    region.color[14][0]=255;
    region.color[14][1]=0;
    region.color[14][2]=255;

    region.level[15]=-150;
    region.color[15][0]=255;
    region.color[15][1]=194;
    region.color[15][2]=204;

    region.levels=16;

    fd=fopen(filename,"r");

    if (fd==NULL)
        fd=fopen(filename,"r");

    if (fd==NULL)
    {
        fd=fopen(filename,"w");

        
        for (x=0; x<region.levels; x++)
            fprintf(fd,"%+4d: %3d, %3d, %3d\n",region.level[x], region.color[x][0], region.color[x][1], region.color[x][2]);

        fclose(fd);
    }

    else
    {
        x=0;
        s=fgets(string,80,fd);

        while (x<128 && feof(fd)==0)
        {
            pointer=strchr(string,';');

            if (pointer!=NULL)
                *pointer=0;

            ok=sscanf(string,"%d: %d, %d, %d", &val[0], &val[1], &val[2], &val[3]);

            if (ok==4)
            {
                if (val[0]<-200)
                    val[0]=-200;

                if (val[0]>+40)
                    val[0]=+40;

                region.level[x]=val[0];

                for (y=1; y<4; y++)
                {
                    if (val[y]>255)
                        val[y]=255;

                    if (val[y]<0)
                        val[y]=0;
                }

                region.color[x][0]=val[1];
                region.color[x][1]=val[2];
                region.color[x][2]=val[3];
                x++;
            }

            s=fgets(string,80,fd);
        }

        fclose(fd);
        region.levels=x;
    }
}


void DoPathLoss(char *filename, unsigned char geo, unsigned char kml, unsigned char ngs, struct site *xmtr, unsigned char txsites)
{
    /* This function generates a topographic map in Portable Pix Map
       (PPM) format based on the content of flags held in the mask[][]
       array (only).  The image created is rotated counter-clockwise
       90 degrees from its representation in dem[][] so that north
       points up and east points right in the image generated. */

    char mapfile[255];
    unsigned width, height, red, green, blue, terrain=0;
    unsigned char found, mask, cityorcounty;
    int indx, x, y, z, x0, y0, loss,match;
    double lat, lon, conversion, one_over_gamma,minwest;
    FILE *fd;

    one_over_gamma=1.0/GAMMA;
    conversion=255.0/pow((double)(max_elevation-min_elevation),one_over_gamma);

    width=(unsigned)(ippd*ReduceAngle(max_west-min_west));
    height=(unsigned)(ippd*ReduceAngle(max_north-min_north));

    LoadLossColors(xmtr[0]);

    if (filename[0]==0)
    {
        strncpy(filename, xmtr[0].filename,254);
        filename[strlen(filename)-4]=0;  /* Remove .qth */
    }

    y=strlen(filename);

    if (y>4)
    {
        if (filename[y-1]=='m' && filename[y-2]=='p' && filename[y-3]=='p' && filename[y-4]=='.')
            y-=4;
    }

    for (x=0; x<y; x++)
    {
        mapfile[x]=filename[x];
    }

    mapfile[x]='.';
    mapfile[x+1]='p';
    mapfile[x+2]='p';
    mapfile[x+3]='m';
    mapfile[x+4]=0;

    minwest=((double)min_west)+dpp;

    if (minwest>360.0)
        minwest-=360.0;

    north=(double)max_north-dpp;

    if (kml || geo)
        south=(double)min_north;    /* No bottom legend */
    else
        south=(double)min_north-(30.0/ppd); /* 30 pixels for bottom legend */

    east=(minwest<180.0?-minwest:360.0-min_west);
    west=(double)(max_west<180?-max_west:360-max_west);


    fd=fopen(mapfile,"wb");

    fprintf(fd,"P6\n%u %u\n255\n",width,(kml?height:height+30));
        if(debug){
    fprintf(stdout,"\nWriting \"%s\" (%ux%u pixmap image)... ",mapfile,width,(kml?height:height+30));
    fflush(stdout);
        }
    for (y=0, lat=north; y<(int)height; y++, lat=north-(dpp*(double)y))
    {
        for (x=0, lon=max_west; x<(int)width; x++, lon=max_west-(dpp*(double)x))
        {
            if (lon<0.0)
                lon+=360.0;

            for (indx=0, found=0; indx<MAXPAGES && found==0;)
            {
                x0=(int)rint(ppd*(lat-(double)dem[indx].min_north));
                y0=mpi-(int)rint(ppd*(LonDiff((double)dem[indx].max_west,lon)));

                if (x0>=0 && x0<=mpi && y0>=0 && y0<=mpi)
                    found=1;
                else
                    indx++;
            }
                        
                        
                        
            if (found)
            {
                mask=dem[indx].mask[x0][y0];
                loss=(dem[indx].signal[x0][y0]);
                cityorcounty=0;

                if(debug){
                    fprintf(stdout,"\n%d\t%d\t%d\t%d",loss,indx,x0,y0);
                    fflush(stdout);
                }
                match=255;

                red=0;
                green=0;
                blue=0;

                if (loss<=region.level[0])
                    match=0;
                else
                {
                    for (z=1; (z<region.levels && match==255); z++)
                    {
                        if (loss>=region.level[z-1] && loss<region.level[z])
                            match=z;
                    }
                }

                if (match<region.levels)
                {
                    red=region.color[match][0];
                    green=region.color[match][1];
                    blue=region.color[match][2];
                }

                 if (mask&2)
                {
                    /* Text Labels: Red or otherwise */

                    if (red>=180 && green<=75 && blue<=75 && loss==0)
                                                fprintf(fd,"%c%c%c",255^red,255^green,255^blue);
                                        else
                                                fprintf(fd,"%c%c%c",255,0,0);

                                        cityorcounty=1;
                }

                else if (mask&4)
                {
                    /* County Boundaries: Black */

                    fprintf(fd,"%c%c%c",0,0,0);

                    cityorcounty=1;
                }

                if (cityorcounty==0)
                {
                    if (loss==0 || (contour_threshold!=0 && loss>abs(contour_threshold)))
                    {
                        if (ngs)  /* No terrain */
                            fprintf(fd,"%c%c%c",255,255,255);
                        else
                        {
                            /* Display land or sea elevation */

                            if (dem[indx].data[x0][y0]==0)
                                fprintf(fd,"%c%c%c",0,0,170);
                            else
                            {
                                terrain=(unsigned)(0.5+pow((double)(dem[indx].data[x0][y0]-min_elevation),one_over_gamma)*conversion);
                                fprintf(fd,"%c%c%c",terrain,terrain,terrain);
                            }
                        }
                    }

                    else
                    {
                        /* Plot path loss in color */

                        if (red!=0 || green!=0 || blue!=0)
                            fprintf(fd,"%c%c%c",red,green,blue);

                        else  /* terrain / sea-level */
                        {
                            if (dem[indx].data[x0][y0]==0)
                                fprintf(fd,"%c%c%c",0,0,170);
                            else
                            {
                                /* Elevation: Greyscale */
                                terrain=(unsigned)(0.5+pow((double)(dem[indx].data[x0][y0]-min_elevation),one_over_gamma)*conversion);
                                fprintf(fd,"%c%c%c",terrain,terrain,terrain);
                            }
                        }
                    }
                }
            }

            else
            {
                /* We should never get here, but if */
                /* we do, display the region as black */

                fprintf(fd,"%c%c%c",0,0,0);
            }
        }
    }



    fclose(fd);

}

void DoSigStr(char *filename, unsigned char geo, unsigned char kml, unsigned char ngs, struct site *xmtr, unsigned char txsites)
{
    /* This function generates a topographic map in Portable Pix Map
       (PPM) format based on the signal strength values held in the
       signal[][] array.  The image created is rotated counter-clockwise
       90 degrees from its representation in dem[][] so that north
       points up and east points right in the image generated. */

    char mapfile[255];
    unsigned width, height, terrain, red, green, blue;
    unsigned char found, mask, cityorcounty;
    int indx, x, y, z=1, x0, y0, signal,match;
    double conversion, one_over_gamma, lat, lon, minwest;
    FILE *fd;

    one_over_gamma=1.0/GAMMA;
    conversion=255.0/pow((double)(max_elevation-min_elevation),one_over_gamma);

    width=(unsigned)(ippd*ReduceAngle(max_west-min_west));
    height=(unsigned)(ippd*ReduceAngle(max_north-min_north));

    LoadSignalColors(xmtr[0]);

    if (filename[0]==0)
    {
        strncpy(filename, xmtr[0].filename,254);
        filename[strlen(filename)-4]=0;  /* Remove .qth */
    }

    y=strlen(filename);

    if (y>4)
    {
        if (filename[y-1]=='m' && filename[y-2]=='p' && filename[y-3]=='p' && filename[y-4]=='.')
            y-=4;
    }

    for (x=0; x<y; x++)
    {
        mapfile[x]=filename[x];
    }

    mapfile[x]='.';
    mapfile[x+1]='p';
    mapfile[x+2]='p';
    mapfile[x+3]='m';
    mapfile[x+4]=0;

    minwest=((double)min_west)+dpp;

    if (minwest>360.0)
        minwest-=360.0;

    north=(double)max_north-dpp;

    
    south=(double)min_north;    /* No bottom legend */
    
    east=(minwest<180.0?-minwest:360.0-min_west);
    west=(double)(max_west<180?-max_west:360-max_west);
    fd=fopen(mapfile,"wb");

    fprintf(fd,"P6\n%u %u\n255\n",width,(kml?height:height+30));
        if(debug){
    fprintf(stdout,"\nWriting \"%s\" (%ux%u pixmap image)... ",mapfile,width,(kml?height:height+30));
    fflush(stdout);
        }
    for (y=0, lat=north; y<(int)height; y++, lat=north-(dpp*(double)y))
    {
        for (x=0, lon=max_west; x<(int)width; x++, lon=max_west-(dpp*(double)x))
        {
            if (lon<0.0)
                lon+=360.0;

            for (indx=0, found=0; indx<MAXPAGES && found==0;)
            {
                x0=(int)rint(ppd*(lat-(double)dem[indx].min_north));
                y0=mpi-(int)rint(ppd*(LonDiff((double)dem[indx].max_west,lon)));

                if (x0>=0 && x0<=mpi && y0>=0 && y0<=mpi)
                    found=1;
                else
                    indx++;
            }

            if (found)
            {
                mask=dem[indx].mask[x0][y0];
                signal=(dem[indx].signal[x0][y0])-100;
                cityorcounty=0;

                if(debug){
                    fprintf(stdout,"\n%d\t%d\t%d\t%d",signal,indx,x0,y0);
                    fflush(stdout);
                }
                    
                match=255;

                red=0;
                green=0;
                blue=0;

                if (signal>=region.level[0])
                    match=0;
                else
                {
                    for (z=1; (z<region.levels && match==255); z++)
                    {
                        if (signal<region.level[z-1] && signal>=region.level[z])
                            match=z;
                    }
                }

                if (match<region.levels)
                {
                    red=region.color[match][0];
                    green=region.color[match][1];
                    blue=region.color[match][2];
                }

                 if (mask&2)
                {
                    /* Text Labels: Red or otherwise */

                    if (red>=180 && green<=75 && blue<=75)
                        fprintf(fd,"%c%c%c",255^red,255^green,255^blue);
                    else
                        fprintf(fd,"%c%c%c",255,0,0);

                    cityorcounty=1;
                }

                else if (mask&4)
                {
                    /* County Boundaries: Black */

                    fprintf(fd,"%c%c%c",0,0,0);

                    cityorcounty=1;
                }

                if (cityorcounty==0)
                {
                    if (contour_threshold!=0 && signal<contour_threshold)
                    {
                        if (ngs)
                            fprintf(fd,"%c%c%c",255,255,255);
                        else
                        {
                            /* Display land or sea elevation */

                            if (dem[indx].data[x0][y0]==0)
                                fprintf(fd,"%c%c%c",0,0,170);
                            else
                            {
                                terrain=(unsigned)(0.5+pow((double)(dem[indx].data[x0][y0]-min_elevation),one_over_gamma)*conversion);
                                fprintf(fd,"%c%c%c",terrain,terrain,terrain);
                            }
                        }
                    }

                    else
                    {
                        /* Plot field strength regions in color */

                        if (red!=0 || green!=0 || blue!=0)
                            fprintf(fd,"%c%c%c",red,green,blue);

                        else  /* terrain / sea-level */
                        {
                            if (ngs)
                                fprintf(fd,"%c%c%c",255,255,255);
                            else
                            {
                                if (dem[indx].data[x0][y0]==0)
                                    fprintf(fd,"%c%c%c",0,0,170);
                                else
                                {
                                    /* Elevation: Greyscale */
                                    terrain=(unsigned)(0.5+pow((double)(dem[indx].data[x0][y0]-min_elevation),one_over_gamma)*conversion);
                                    fprintf(fd,"%c%c%c",terrain,terrain,terrain);
                                }
                            }
                        }
                    }
                }
            }

            else
            {
                /* We should never get here, but if */
                /* we do, display the region as black */

                fprintf(fd,"%c%c%c",0,0,0);
            }
        }
    }


    fclose(fd);
   
}

void DoRxdPwr(char *filename, unsigned char geo, unsigned char kml, unsigned char ngs, struct site *xmtr, unsigned char txsites)
{
    /* This function generates a topographic map in Portable Pix Map
       (PPM) format based on the signal power level values held in the
       signal[][] array.  The image created is rotated counter-clockwise
       90 degrees from its representation in dem[][] so that north
       points up and east points right in the image generated. */

    char mapfile[255];
    unsigned width, height, terrain, red, green, blue;
    unsigned char found, mask, cityorcounty;
    int indx, x, y, z=1, x0, y0, dBm, match;
    double conversion, one_over_gamma, lat, lon, minwest;
    FILE *fd;

    one_over_gamma=1.0/GAMMA;
    conversion=255.0/pow((double)(max_elevation-min_elevation),one_over_gamma);

    width=(unsigned)(ippd*ReduceAngle(max_west-min_west));
    height=(unsigned)(ippd*ReduceAngle(max_north-min_north));

    LoadDBMColors(xmtr[0]);

    if (filename[0]==0)
    {
        strncpy(filename, xmtr[0].filename,254);
        filename[strlen(filename)-4]=0;  /* Remove .qth */
    }

    y=strlen(filename);

    if (y>4)
    {
        if (filename[y-1]=='m' && filename[y-2]=='p' && filename[y-3]=='p' && filename[y-4]=='.')
            y-=4;
    }

    for (x=0; x<y; x++)
    {
        mapfile[x]=filename[x];
    }

    mapfile[x]='.';
    mapfile[x+1]='p';
    mapfile[x+2]='p';
    mapfile[x+3]='m';
    mapfile[x+4]=0;

    minwest=((double)min_west)+dpp;

    if (minwest>360.0)
        minwest-=360.0;

    north=(double)max_north-dpp;

    
    south=(double)min_north;    /* No bottom legend */
    

    east=(minwest<180.0?-minwest:360.0-min_west);
    west=(double)(max_west<180?-max_west:360-max_west);

    fd=fopen(mapfile,"wb");

    fprintf(fd,"P6\n%u %u\n255\n",width,(kml?height:height+30));
        if(debug){
    fprintf(stdout,"\nWriting \"%s\" (%ux%u pixmap image)... ",mapfile,width,(kml?height:height+30));
    fflush(stdout);
        }
        // WriteKML()
        //writeKML(xmtr,filename);
        
    for (y=0, lat=north; y<(int)height; y++, lat=north-(dpp*(double)y))
    {
        for (x=0, lon=max_west; x<(int)width; x++, lon=max_west-(dpp*(double)x))
        {
            if (lon<0.0)
                lon+=360.0;

            for (indx=0, found=0; indx<MAXPAGES && found==0;)
            {
                x0=(int)rint(ppd*(lat-(double)dem[indx].min_north));
                y0=mpi-(int)rint(ppd*(LonDiff((double)dem[indx].max_west,lon)));

                if (x0>=0 && x0<=mpi && y0>=0 && y0<=mpi)
                    found=1;
                else
                    indx++;
            }

            if (found)
            {
                mask=dem[indx].mask[x0][y0];
                dBm=(dem[indx].signal[x0][y0])-200;
                cityorcounty=0;

 
                 if(debug){
                fprintf(stdout,"\n%d\t%d\t%d\t%d",dBm,indx,x0,y0);
                fflush(stdout);
                    }
                
                match=255;

                red=0;
                green=0;
                blue=0;

                if (dBm>=region.level[0])
                    match=0;
                else
                {
                    for (z=1; (z<region.levels && match==255); z++)
                    {
                        if (dBm<region.level[z-1] && dBm>=region.level[z])
                            match=z;
                    }
                }

                if (match<region.levels)
                {
                    red=region.color[match][0];
                    green=region.color[match][1];
                    blue=region.color[match][2];
                }

    
                 if (mask&2)
                {
                    /* Text Labels: Red or otherwise */

                    if (red>=180 && green<=75 && blue<=75 && dBm!=0)
                        fprintf(fd,"%c%c%c",255^red,255^green,255^blue);
                    else
                        fprintf(fd,"%c%c%c",255,0,0);

                    cityorcounty=1;
                }

                else if (mask&4)
                {
                    /* County Boundaries: Black */

                    fprintf(fd,"%c%c%c",0,0,0);

                    cityorcounty=1;
                }

                if (cityorcounty==0)
                {
                    if (contour_threshold!=0 && dBm<contour_threshold)
                    {
                        if (ngs) /* No terrain */
                            fprintf(fd,"%c%c%c",255,255,255);
                        else
                        {
                            /* Display land or sea elevation */

                            if (dem[indx].data[x0][y0]==0)
                                fprintf(fd,"%c%c%c",0,0,170);
                            else
                            {
                                terrain=(unsigned)(0.5+pow((double)(dem[indx].data[x0][y0]-min_elevation),one_over_gamma)*conversion);
                                fprintf(fd,"%c%c%c",terrain,terrain,terrain);
                            }
                        }
                    }

                    else
                    {
                        /* Plot signal power level regions in color */

                        if (red!=0 || green!=0 || blue!=0)
                            fprintf(fd,"%c%c%c",red,green,blue);

                        else  /* terrain / sea-level */
                        {
                            if (ngs)
                                fprintf(fd,"%c%c%c",255,255,255);
                            else
                            {
                                if (dem[indx].data[x0][y0]==0)
                                    fprintf(fd,"%c%c%c",0,0,170);
                                else
                                {
                                    /* Elevation: Greyscale */
                                    terrain=(unsigned)(0.5+pow((double)(dem[indx].data[x0][y0]-min_elevation),one_over_gamma)*conversion);
                                    fprintf(fd,"%c%c%c",terrain,terrain,terrain);
                                }
                            }
                        }
                    }
                }
            }

            else
            {
                /* We should never get here, but if */
                /* we do, display the region as black */

                fprintf(fd,"%c%c%c",0,0,0);
            }
        }
    }

    fclose(fd);
   
}

void DoLOS(char *filename, unsigned char geo, unsigned char kml, unsigned char ngs, struct site *xmtr, unsigned char txsites)
{
    /* This function generates a topographic map in Portable Pix Map
       (PPM) format based on the signal power level values held in the
       signal[][] array.  The image created is rotated counter-clockwise
       90 degrees from its representation in dem[][] so that north
       points up and east points right in the image generated. */

    char mapfile[255];
    unsigned width, height, terrain;
    unsigned char found, mask;
    int indx, x, y, x0, y0;
    double conversion, one_over_gamma, lat, lon, minwest;
    FILE *fd;

    one_over_gamma=1.0/GAMMA;
    conversion=255.0/pow((double)(max_elevation-min_elevation),one_over_gamma);

    width=(unsigned)(ippd*ReduceAngle(max_west-min_west));
    height=(unsigned)(ippd*ReduceAngle(max_north-min_north));

    if (filename[0]==0)
    {
        strncpy(filename, xmtr[0].filename,254);
        filename[strlen(filename)-4]=0;  /* Remove .qth */
    }

    y=strlen(filename);

    if (y>4)
    {
        if (filename[y-1]=='m' && filename[y-2]=='p' && filename[y-3]=='p' && filename[y-4]=='.')
            y-=4;
    }

    for (x=0; x<y; x++)
    {
        mapfile[x]=filename[x];
    }

    mapfile[x]='.';
    mapfile[x+1]='p';
    mapfile[x+2]='p';
    mapfile[x+3]='m';
    mapfile[x+4]=0;

    minwest=((double)min_west)+dpp;

    if (minwest>360.0)
        minwest-=360.0;

    north=(double)max_north-dpp;

    
    south=(double)min_north;    /* No bottom legend */
    

    east=(minwest<180.0?-minwest:360.0-min_west);
    west=(double)(max_west<180?-max_west:360-max_west);

    fd=fopen(mapfile,"wb");

    fprintf(fd,"P6\n%u %u\n255\n",width,(kml?height:height+30));
        if(debug){
    fprintf(stdout,"\nWriting \"%s\" (%ux%u pixmap image)... ",mapfile,width,(kml?height:height+30));
    fflush(stdout);
        }
        // WriteKML()
        //writeKML(xmtr,filename);
        
    for (y=0, lat=north; y<(int)height; y++, lat=north-(dpp*(double)y))
    {
        for (x=0, lon=max_west; x<(int)width; x++, lon=max_west-(dpp*(double)x))
        {
            if (lon<0.0)
                lon+=360.0;

            for (indx=0, found=0; indx<MAXPAGES && found==0;)
            {
                x0=(int)rint(ppd*(lat-(double)dem[indx].min_north));
                y0=mpi-(int)rint(ppd*(LonDiff((double)dem[indx].max_west,lon)));

                if (x0>=0 && x0<=mpi && y0>=0 && y0<=mpi)
                    found=1;
                else
                    indx++;
            }

             if (found)
                     {
                            mask=dem[indx].mask[x0][y0];

                            if (mask&2)
                                   /* Text Labels: Red */
                                   fprintf(fd,"%c%c%c",255,0,0);

                            else if (mask&4)
                                   /* County Boundaries: Light Cyan */
                                   fprintf(fd,"%c%c%c",128,128,255);

                            else switch (mask&57)
                            {
                                   case 1:
                                   /* TX1: Green */
                                   fprintf(fd,"%c%c%c",0,255,0);
                                   break;

                                   case 8:
                                   /* TX2: Cyan */
                                   fprintf(fd,"%c%c%c",0,255,255);
                                   break;

                                   case 9:
                                   /* TX1 + TX2: Yellow */
                                   fprintf(fd,"%c%c%c",255,255,0);
                                   break;

                                   case 16:
                                   /* TX3: Medium Violet */
                                   fprintf(fd,"%c%c%c",147,112,219);
                                   break;

                                   case 17:
                                   /* TX1 + TX3: Pink */
                                   fprintf(fd,"%c%c%c",255,192,203);
                                   break;

                                   case 24:
                                   /* TX2 + TX3: Orange */
                                   fprintf(fd,"%c%c%c",255,165,0);
                                   break;

                                   case 25:
                                   /* TX1 + TX2 + TX3: Dark Green */
                                   fprintf(fd,"%c%c%c",0,100,0);
                                   break;

                                   case 32:
                                   /* TX4: Sienna 1 */
                                   fprintf(fd,"%c%c%c",255,130,71);
                                   break;

                                   case 33:
                                   /* TX1 + TX4: Green Yellow */
                                   fprintf(fd,"%c%c%c",173,255,47);
                                   break;

                                   case 40:
                                   /* TX2 + TX4: Dark Sea Green 1 */
                                   fprintf(fd,"%c%c%c",193,255,193);
                                   break;

                                   case 41:
                                   /* TX1 + TX2 + TX4: Blanched Almond */
                                   fprintf(fd,"%c%c%c",255,235,205);
                                   break;

                                   case 48:
                                   /* TX3 + TX4: Dark Turquoise */
                                   fprintf(fd,"%c%c%c",0,206,209);
                                   break;

                                   case 49:
                                   /* TX1 + TX3 + TX4: Medium Spring Green */
                                   fprintf(fd,"%c%c%c",0,250,154);
                                   break;

                                   case 56:
                                   /* TX2 + TX3 + TX4: Tan */
                                   fprintf(fd,"%c%c%c",210,180,140);
                                   break;

                                   case 57:
                                   /* TX1 + TX2 + TX3 + TX4: Gold2 */
                                   fprintf(fd,"%c%c%c",238,201,0);
                                   break;

                                   default:
                                   if (ngs)  /* No terrain */
                                          fprintf(fd,"%c%c%c",255,255,255);
                                   else
                                   {
                                          /* Sea-level: Medium Blue */
                                          if (dem[indx].data[x0][y0]==0)
                                                 fprintf(fd,"%c%c%c",0,0,170);
                                          else
                                          {
                                                 /* Elevation: Greyscale */
                                                 terrain=(unsigned)(0.5+pow((double)(dem[indx].data[x0][y0]-min_elevation),one_over_gamma)*conversion);
                                                 fprintf(fd,"%c%c%c",terrain,terrain,terrain);
                                          }
                                   }
                            }
                     }

                     else
                     {
                            /* We should never get here, but if */
                            /* we do, display the region as black */

                            fprintf(fd,"%c%c%c",0,0,0);
                     }
              }
       }

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
                        snprintf(string,19,"%d:%d:%d:%d-hd",x, x+1, ymin, ymax);
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
                        snprintf(string,19,"%d:%d:%d:%d-hd",x, x+1, ymin, ymax);
                    else
                        snprintf(string,16,"%d:%d:%d:%d",x, x+1, ymin, ymax);
                }  
            
                LoadSDF(string,winfiles);
            }
    }
}



void LoadUDT(char *filename)
{
    /* This function reads a file containing User-Defined Terrain
       features for their addition to the digital elevation model
       data used by SPLAT!.  Elevations in the UDT file are evaluated
       and then copied into a temporary file under /tmp.  Then the
       contents of the temp file are scanned, and if found to be unique,
       are added to the ground elevations described by the digital
       elevation data already loaded into memory. */

    int    i, x, y, z, ypix, xpix, tempxpix, tempypix, fd=0, n=0, pixelfound=0;
    char    input[80], str[3][80], tempname[15], *pointer=NULL, *s=NULL;
    double    latitude, longitude, height, tempheight;
    FILE    *fd1=NULL, *fd2=NULL;

    strcpy(tempname,"/tmp/XXXXXX\0");

    fd1=fopen(filename,"r");

    if (fd1!=NULL)
    {
        fd=mkstemp(tempname);
        fd2=fopen(tempname,"w");

        s=fgets(input,78,fd1);

        pointer=strchr(input,';');

        if (pointer!=NULL)
            *pointer=0;

        
        while (feof(fd1)==0)
        {
            /* Parse line for latitude, longitude, height */

            for (x=0, y=0, z=0; x<78 && input[x]!=0 && z<3; x++)
            {
                if (input[x]!=',' && y<78)
                {
                    str[z][y]=input[x];
                    y++;
                }

                else
                {
                    str[z][y]=0;
                    z++;
                    y=0;
                }
            }

            latitude=ReadBearing(str[0]);
            longitude=ReadBearing(str[1]);

            if (longitude<0.0)
                longitude+=360; 
            
            /* Remove <CR> and/or <LF> from antenna height string */

            for (i=0; str[2][i]!=13 && str[2][i]!=10 && str[2][i]!=0; i++);

            str[2][i]=0;

            /* The terrain feature may be expressed in either
               feet or meters.  If the letter 'M' or 'm' is
               discovered in the string, then this is an
               indication that the value given is expressed
               in meters.  Otherwise the height is interpreted
               as being expressed in feet.  */

            for (i=0; str[2][i]!='M' && str[2][i]!='m' && str[2][i]!=0 && i<48; i++);

            if (str[2][i]=='M' || str[2][i]=='m')
            {
                str[2][i]=0;
                height=rint(atof(str[2]));
            }

            else
            {
                str[2][i]=0;
                height=rint(METERS_PER_FOOT*atof(str[2]));
            }

            if (height>0.0)
                fprintf(fd2,"%d, %d, %f\n",(int)rint(latitude/dpp), (int)rint(longitude/dpp), height);
        

            s=fgets(input,78,fd1);

            pointer=strchr(input,';');

            if (pointer!=NULL)
                *pointer=0;
        }

        fclose(fd1);
        fclose(fd2);
        close(fd);

        
        fd1=fopen(tempname,"r");
        fd2=fopen(tempname,"r");

        y=0;

        n=fscanf(fd1,"%d, %d, %lf", &xpix, &ypix, &height);

        do
        {
            x=0;
            z=0;

            n=fscanf(fd2,"%d, %d, %lf", &tempxpix, &tempypix, &tempheight);

            do
            {
                if (x>y && xpix==tempxpix && ypix==tempypix)
                {
                        z=1;  /* Dupe! */

                        if (tempheight>height)
                            height=tempheight;
                }

                else
                {
                    n=fscanf(fd2,"%d, %d, %lf", &tempxpix, &tempypix, &tempheight);
                    x++;
                }

            } while (feof(fd2)==0 && z==0);

            if (z==0)  /* No duplicate found */
        
                //fprintf(stdout,"%lf, %lf \n",xpix*dpp, ypix*dpp);
                fflush(stdout);
                pixelfound = AddElevation(xpix*dpp, ypix*dpp, height);
                //fprintf(stdout,"%d \n",pixelfound);
                fflush(stdout);
                
            n=fscanf(fd1,"%d, %d, %lf", &xpix, &ypix, &height);
            y++;

            rewind(fd2);

        } while (feof(fd1)==0);

        fclose(fd1);
        fclose(fd2);
        unlink(tempname);
    }

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

void ObstructionAnalysis(struct site xmtr, struct site rcvr, double f, FILE *outfile)
{
    /* Perform an obstruction analysis along the
       path between receiver and transmitter. */

    int    x;
    struct    site site_x;
    double    h_r, h_t, h_x, h_r_orig, cos_tx_angle, cos_test_angle,
        cos_tx_angle_f1, cos_tx_angle_fpt6, d_tx, d_x,
        h_r_f1, h_r_fpt6, h_f, h_los, lambda=0.0;
    char    string[255], string_fpt6[255], string_f1[255];

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
                fprintf(outfile,"Between %s and %s, obstructions were detected at:\n\n",rcvr.name,xmtr.name);

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
            snprintf(string,150,"\nAntenna at %s must be raised to at least %.2f meters AGL\nto clear all obstructions detected.\n",rcvr.name, METERS_PER_FOOT*(h_r-GetElevation(rcvr)-earthradius));
        else
            snprintf(string,150,"\nAntenna at %s must be raised to at least %.2f feet AGL\nto clear all obstructions detected.\n",rcvr.name, h_r-GetElevation(rcvr)-earthradius);
    }

    else
        snprintf(string,150,"\nNo obstructions to LOS path due to terrain were detected\n");

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


void PathReport(struct site source, struct site destination, char *name, char graph_it, int propmodel, int pmenv)
{
    /* This function writes a PPA Path Report (name.txt) to
       the filesystem.  If (graph_it == 1), then gnuplot is invoked
       to generate an appropriate output file indicating the Longley-Rice
       model loss between the source and destination locations.
           "filename" is the name assigned to the output file generated
       by gnuplot.  The filename extension is used to set gnuplot's
       terminal setting and output file type.  If no extension is
       found, .png is assumed. */

    int    x, y, z, errnum;
    char    basename[255], term[30], ext[15], strmode[100],
        report_name[80], block=0;
    double    maxloss=-100000.0, minloss=100000.0, angle1, angle2, 
        azimuth, pattern=1.0, patterndB=0.0,
        total_loss=0.0, cos_xmtr_angle, cos_test_angle=0.0,
        source_alt, test_alt, dest_alt, source_alt2, dest_alt2,
        distance, elevation, four_thirds_earth,
        free_space_loss=0.0, eirp=0.0, voltage, rxp, power_density, dkm, txelev;
    FILE    *fd=NULL, *fd2=NULL;

    snprintf(report_name,80,"%s.txt%c",name,0);
    four_thirds_earth=FOUR_THIRDS*EARTHRADIUS;

    fd2=fopen(report_name,"w");

    fprintf(fd2,"\n\t\t--==[ Path Profile Analysis ]==--\n\n");
    fprintf(fd2,"Transmitter site: %s\n",source.name);

    if (source.lat>=0.0)
    {
        fprintf(fd2,"Site location: %.4f North / %.4f West\n",source.lat, source.lon);
    }

    else
    {

        fprintf(fd2,"Site location: %.4f South / %.4f West\n",-source.lat, source.lon);
    }
    

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
            fprintf(fd2,"Depression\n");
        else
            fprintf(fd2,"Elevation\n");
    }

    /* Receiver */

    fprintf(fd2,"Receiver site: %s\n",destination.name);

    if (destination.lat>=0.0)
    {
        fprintf(fd2,"Site location: %.4f North / %.4f West\n",destination.lat, destination.lon);
    }

    else
    {
        fprintf(fd2,"Site location: %.4f South / %.4f West\n",-destination.lat, destination.lon);
    }

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
    }

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
            fprintf(fd2,"Maritime Temperate, Over Land");
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

            elev[0]=y-1;    /* (number of points - 1) */

            /* Distance between elevation samples */

            elev[1]=METERS_PER_MILE*(path.distance[y]-path.distance[y-1]);

			/*
            point_to_point(elev, source.alt*METERS_PER_FOOT, 
            destination.alt*METERS_PER_FOOT, LR.eps_dielect,
            LR.sgm_conductivity, LR.eno_ns_surfref, LR.frq_mhz,
            LR.radio_climate, LR.pol, LR.conf, LR.rel, loss,
            strmode, errnum);
			*/
			dkm=(elev[1]*elev[0])/1000; // km
			txelev=elev[2]+(source.alt*METERS_PER_FOOT);
			
			 switch (propmodel)
            {
              case 1:
                // Longley Rice ITM
                 point_to_point_ITM(elev,source.alt*METERS_PER_FOOT,
                    destination.alt*METERS_PER_FOOT, LR.eps_dielect,
                    LR.sgm_conductivity, LR.eno_ns_surfref, LR.frq_mhz,
                    LR.radio_climate, LR.pol, LR.conf, LR.rel, loss,
                    strmode, errnum);
                 break;
              case 3:
                //HATA 1, 2 & 3
                loss=HATApathLoss(LR.frq_mhz,txelev,path.elevation[y]+(destination.alt*METERS_PER_FOOT),dkm, pmenv);
                break;
              case 4:
                // COST231-HATA
				loss=ECC33pathLoss(LR.frq_mhz,txelev,path.elevation[y]+(destination.alt*METERS_PER_FOOT),dkm, pmenv);
                break;
              case 5:
                // SUI
                loss=SUIpathLoss(LR.frq_mhz,txelev,path.elevation[y]+(destination.alt*METERS_PER_FOOT),dkm, pmenv);
                break;
              case 6:
                loss=COST231pathLoss(LR.frq_mhz,txelev,path.elevation[y]+(destination.alt*METERS_PER_FOOT),dkm, pmenv);
                break;
              case 7:
			  // ITU-R P.525 Free space path loss
                loss=FSPLpathLoss(LR.frq_mhz,dkm);
                break;				
              case 8:
                 // ITWOM 3.0 
                 point_to_point(elev,source.alt*METERS_PER_FOOT,
                    destination.alt*METERS_PER_FOOT, LR.eps_dielect,
                    LR.sgm_conductivity, LR.eno_ns_surfref, LR.frq_mhz,
                    LR.radio_climate, LR.pol, LR.conf, LR.rel, loss,
                    strmode, errnum);
                 break;
              case 9:
			   // Ericsson
                loss=EricssonpathLoss(LR.frq_mhz,txelev,path.elevation[y]+(destination.alt*METERS_PER_FOOT),dkm, pmenv);
                 break;

                
              default:
                 point_to_point_ITM(elev,source.alt*METERS_PER_FOOT,
                    destination.alt*METERS_PER_FOOT, LR.eps_dielect,
                    LR.sgm_conductivity, LR.eno_ns_surfref, LR.frq_mhz,
                    LR.radio_climate, LR.pol, LR.conf, LR.rel, loss,
                    strmode, errnum);
  
            }
			
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

   
            if (total_loss>maxloss)
                maxloss=total_loss;

            if (total_loss<minloss)
                minloss=total_loss;
        }

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


    fprintf(stdout,"Path loss (dB), Received Power (dBm), Field strength (dBuV):\n%.1f\n%.1f\n%.1f",loss,dBm,field_strength);

    
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

        if (ext[0]==0)    /* No extension -- Default is png */
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
        fprintf(fd,"set title \"Path Loss Profile Along Path Between %s and %s (%.2f%c azimuth)\"\n",destination.name, source.name, Azimuth(destination,source),176);

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

void SeriesData(struct site source, struct site destination, char *name, unsigned char fresnel_plot, unsigned char normalised)
{

    int    x, y, z;
    char    basename[255], term[30], ext[15], profilename[255], referencename[255],cluttername[255],curvaturename[255],fresnelname[255],fresnel60name[255];
    double    a, b, c, height=0.0, refangle, cangle, maxheight=-100000.0,
        minheight=100000.0, lambda=0.0, f_zone=0.0, fpt6_zone=0.0,
        nm=0.0, nb=0.0, ed=0.0, es=0.0, r=0.0, d=0.0, d1=0.0,
        terrain, azimuth, distance, minterrain=100000.0,
        minearth=100000.0;
    struct    site remote;
    FILE    *fd=NULL, *fd1=NULL, *fd2=NULL, *fd3=NULL, *fd4=NULL, *fd5=NULL;

    ReadPath(destination,source);  
    azimuth=Azimuth(destination,source);
    distance=Distance(destination,source);
    refangle=ElevationAngle(destination,source);
    b=GetElevation(destination)+destination.alt+earthradius;

    if (fresnel_plot)
    {
        lambda=9.8425e8/(LR.frq_mhz*1e6);
        d=5280.0*path.distance[path.length-1];
    }

    if (normalised)
    {
        ed=GetElevation(destination);
        es=GetElevation(source);
        nb=-destination.alt-ed;
        nm=(-source.alt-es-nb)/(path.distance[path.length-1]);
    }

    strcpy(profilename,name);
    strcat(profilename,"_profile\0");
    strcpy(referencename,name);
    strcat(referencename,"_reference\0");
    strcpy(cluttername,name);
    strcat(cluttername,"_clutter\0");
    strcpy(curvaturename,name);
    strcat(curvaturename,"_curvature\0");
    strcpy(fresnelname,name);
    strcat(fresnelname,"_fresnel\0");
    strcpy(fresnel60name,name);
    strcat(fresnel60name,"_fresnel60\0");    
    
    fd=fopen(profilename,"wb");
    if (clutter>0.0)
        fd1=fopen(cluttername,"wb");
    fd2=fopen(referencename,"wb");
    fd5=fopen(curvaturename, "wb");

    if ((LR.frq_mhz>=20.0) && (LR.frq_mhz<=100000.0) && fresnel_plot)
    {
        fd3=fopen(fresnelname, "wb");
        fd4=fopen(fresnel60name, "wb");
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

        if (normalised)
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
            if (METERS_PER_FOOT*height > 0){
                fprintf(fd,"%.3f %.3f\n",KM_PER_MILE*path.distance[x],METERS_PER_FOOT*height);
            }
            
            if (fd1!=NULL && x>0 && x<path.length-2)
                fprintf(fd1,"%.3f %.3f\n",KM_PER_MILE*path.distance[x],METERS_PER_FOOT*(terrain==0.0?height:(height+clutter)));

            fprintf(fd2,"%.3f %.3f\n",KM_PER_MILE*path.distance[x],METERS_PER_FOOT*r);
            fprintf(fd5,"%.3f %.3f\n",KM_PER_MILE*path.distance[x],METERS_PER_FOOT*(height-terrain));
        
        }

        else
        {
            fprintf(fd,"%.3f %.3f\n",path.distance[x],height);

            if (fd1!=NULL && x>0 && x<path.length-2)
                fprintf(fd1,"%.3f %.3f\n",path.distance[x],(terrain==0.0?height:(height+clutter)));

            fprintf(fd2,"%.3f %.3f\n",path.distance[x],r);
            fprintf(fd5,"%.3f %.3f\n",path.distance[x],height-terrain);
        }

    
        if ((LR.frq_mhz>=20.0) && (LR.frq_mhz<=100000.0) && fresnel_plot)
        {
            if (metric)
            {
                fprintf(fd3,"%.3f %.3f\n",KM_PER_MILE*path.distance[x],METERS_PER_FOOT*f_zone);
                fprintf(fd4,"%.3f %.3f\n",KM_PER_MILE*path.distance[x],METERS_PER_FOOT*fpt6_zone);
            }

            else
            {
                fprintf(fd3,"%.3f %.3f\n",path.distance[x],f_zone);
                fprintf(fd4,"%.3f %.3f\n",path.distance[x],fpt6_zone);
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
    } // End of loop

    
    if (normalised)
        r=-(nm*path.distance[path.length-1])-nb;
    else
        r=0.0;

    if (metric)
    {
        fprintf(fd,"%.3f %.3f",KM_PER_MILE*path.distance[path.length-1],METERS_PER_FOOT*r);
        fprintf(fd2,"%.3f %.3f",KM_PER_MILE*path.distance[path.length-1],METERS_PER_FOOT*r);
    }

    else
    {
        fprintf(fd,"%.3f %.3f",path.distance[path.length-1],r);
        fprintf(fd2,"%.3f %.3f",path.distance[path.length-1],r);
    }

    if ((LR.frq_mhz>=20.0) && (LR.frq_mhz<=100000.0) && fresnel_plot)
    {
        if (metric)
        {
            fprintf(fd3,"%.3f %.3f",KM_PER_MILE*path.distance[path.length-1],METERS_PER_FOOT*r);
            fprintf(fd4,"%.3f %.3f",KM_PER_MILE*path.distance[path.length-1],METERS_PER_FOOT*r);
        }

        else
        {
            fprintf(fd3,"%.3f %.3f",path.distance[path.length-1],r);
            fprintf(fd4,"%.3f %.3f",path.distance[path.length-1],r);
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
        strncpy(basename,"profile\0",8);
        strncpy(term,"png\0",4);
        strncpy(ext,"png\0",4);
    }

    else
    {
        
        ext[0]=0;
        y=strlen(name);
        strncpy(basename,name,254);

        for (x=y-1; x>0 && name[x]!='.'; x--);

        if (x>0)  
        {
            for (z=x+1; z<=y && (z-(x+1))<10; z++)
            {
                ext[z-(x+1)]=tolower(name[z]);
                term[z-(x+1)]=name[z];
            }

            ext[z-(x+1)]=0;  
            term[z-(x+1)]=0;
            basename[x]=0;
        }

        if (ext[0]==0)    
        {
            strncpy(term,"png\0",4);
            strncpy(ext,"png\0",4);
        }
    }

    fprintf(stdout,"\n");
    fflush(stdout);
    
}

int main(int argc, char *argv[])
{
    int     x, y, z=0, min_lat, min_lon, max_lat, max_lon,
            rxlat, rxlon, txlat, txlon, west_min, west_max,
            nortRxHin, nortRxHax, propmodel, winfiles,knifeedge=0,ppa=0,normalise=0, haf=0, pmenv=1;

    unsigned char    LRmap=0, txsites=0,topomap=0, geo=0, kml=0, area_mode=0, max_txsites,ngs=0;

    char    mapfile[255], longley_file[255], udt_file[255],ano_filename[255];

    double  altitude=0.0, altitudeLR=0.0, tx_range=0.0,
            rx_range=0.0, deg_range=0.0, deg_limit=0.0,
            deg_range_lon;

    strncpy(ss_name,"Signal Server\0",14);
    
    if (argc==1)
    {
		
                fprintf(stdout,"\n\t\t -- %s %.2f --\n",ss_name, version);
		fprintf(stdout,"\tCompiled for %d tiles at %d pixels/degree\n\n",MAXPAGES,IPPD);
                fprintf(stdout,"     -d Directory containing .sdf tiles\n");
                fprintf(stdout,"     -lat Tx Latitude (decimal degrees) -70/+70\n");
                fprintf(stdout,"     -lon Tx Longitude (decimal degrees) -180/+180\n");
                fprintf(stdout,"     -txh Tx Height (above ground)\n");
                fprintf(stdout,"     -rla (Optional) Rx Latitude for PPA (decimal degrees) -70/+70\n");
                fprintf(stdout,"     -rlo (Optional) Rx Longitude for PPA (decimal degrees) -180/+180\n");
                fprintf(stdout,"     -f Tx Frequency (MHz) 20MHz to 100GHz (LOS after 20GHz)\n");
                fprintf(stdout,"     -erp Tx Effective Radiated Power (Watts)\n");
                fprintf(stdout,"     -rxh Rx Height(s) (optional. Default=0.1)\n");
                fprintf(stdout,"     -rt Rx Threshold (dB / dBm / dBuV/m)\n");
                fprintf(stdout,"     -hp Horizontal Polarisation (default=vertical)\n");
                fprintf(stdout,"     -gc Ground clutter (feet/meters)\n");
                fprintf(stdout,"     -udt User defined terrain filename\n");
                fprintf(stdout,"     -dbm Plot Rxd signal power instead of field strength\n");
                fprintf(stdout,"     -m Metric units of measurement\n");
                fprintf(stdout,"     -te Terrain code 1-6 (optional)\n");
                fprintf(stdout,"     -terdic Terrain dielectric value 2-80 (optional)\n");
                fprintf(stdout,"     -tercon Terrain conductivity 0.01-0.0001 (optional)\n");
                fprintf(stdout,"     -cl Climate code 1-6 (optional)\n");
                fprintf(stdout,"     -o Filename. Required. \n");
                fprintf(stdout,"     -R Radius (miles/kilometers)\n");
                fprintf(stdout,"     -res Pixels per degree. 300/600/1200(default)/3600 (optional)\n");
                fprintf(stdout,"     -t Terrain background\n");
                fprintf(stdout,"     -pm Prop model. 1: ITM, 2: LOS, 3: Hata, 4: ECC33,\n");
                fprintf(stdout,"     		5: SUI, 6: COST-Hata, 7: FSPL, 8: ITWOM, 9: Ericsson\n");
		fprintf(stdout,"     -pe Prop model mode: 1=Urban,2=Suburban,3=Rural\n");
                fprintf(stdout,"     -ked Knife edge diffraction (Default for ITM)\n");
                fprintf(stdout,"     -ng Normalise Path Profile graph\n");
                fprintf(stdout,"     -haf Halve 1 or 2 (optional)\n");
 

        fflush(stdout);

        return 1;
    }

    y=argc-1;

    kml=0;
    geo=0;
    dbm=0;
    gpsav=0;
    metric=0;
    string[0]=0;
    mapfile[0]=0;
    clutter=0.0;
    forced_erp=-1.0;
    forced_freq=0.0;
    sdf_path[0]=0;
    udt_file[0]=0;
    path.length=0;
    max_txsites=30;
    fzone_clearance=0.6;
    contour_threshold=0;
    
    longley_file[0]=0;
    ano_filename[0]=0;
    earthradius=EARTHRADIUS;
    max_range=1.0;
    propmodel=1; //ITM
    winfiles=0;           

    lat=0;
    lon=0;
    txh=0;
    ngs=1; // no terrain background
    kml=1;
    LRmap=1;
    area_mode=1;
    ippd=IPPD; // default resolution

        sscanf("0.1","%lf",&altitudeLR);

        // Defaults
        LR.eps_dielect=15.0; // Farmland
        LR.sgm_conductivity=0.005; // Farmland
        LR.eno_ns_surfref=301.0; 
        LR.frq_mhz=19.0; // Deliberately too low
        LR.radio_climate=5; // continental
        LR.pol=1; // vert
        LR.conf=0.50; 
        LR.rel=0.50; 
        LR.erp=0.0; // will default to Path Loss

        tx_site[0].lat=91.0;
        tx_site[0].lon=361.0;
        tx_site[1].lat=91.0;
        tx_site[1].lon=361.0;

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

            if (strcmp(argv[x],"-res")==0)
        {
                z=x+1;

            if (z<=y && argv[z][0] && argv[z][0]!='-')
            {
                            sscanf(argv[z],"%d",&ippd);
                             
                            switch (ippd)
                                {
                                case 300: 
                                    MAXRAD=300;
                                    jgets=3;
                                break;

                                case 600: 
                                    MAXRAD=150;
                                    jgets=1;
                                break;
                
                                case 3600: 
                                    MAXRAD=100;
                                    ippd=3600;
                                    jgets=0;
                                break;

                                default:
                                    MAXRAD=100;
                                    ippd=1200;
                                    jgets=0;
                                    break;
                        }
                        }
            }
        if (strcmp(argv[x],"-R")==0)
        {
            z=x+1;

            if (z<=y && argv[z][0] && argv[z][0]!='-')
            {
                sscanf(argv[z],"%lf",&max_range);
                
            }
        }


        if (strcmp(argv[x],"-gc")==0)
        {
            z=x+1;

            if (z<=y && argv[z][0] && argv[z][0]!='-')
            {
                sscanf(argv[z],"%lf",&clutter);

                if (clutter<0.0)
                    clutter=0.0;
            }
        }


        if (strcmp(argv[x],"-o")==0)
        {
            z=x+1;

            if (z<=y && argv[z][0] && argv[z][0]!='-')
                        {
                 strncpy(mapfile,argv[z],253);
                                 strncpy (tx_site[0].name,"Tx",2);
                                 strncpy (tx_site[0].filename,argv[z],253);
                                 LoadPAT(argv[z]);
                                                              
                        }
        }

            
        if (strcmp(argv[x],"-rt")==0)
        {
            z=x+1;

            if (z<=y && argv[z][0]) /* A minus argument is legal here */
                sscanf(argv[z],"%d",&contour_threshold);
        }

        if (strcmp(argv[x],"-m")==0)
                {
            metric=1;
                        
                }

        if (strcmp(argv[x],"-t")==0)
        {
        ngs=0; // greyscale background
        }

        if (strcmp(argv[x],"-dbm")==0)
            dbm=1;

        if (strcmp(argv[x],"-d")==0)
        {
            z=x+1;

            if (z<=y && argv[z][0] && argv[z][0]!='-')
                strncpy(sdf_path,argv[z],253);
        }

        if (strcmp(argv[x],"-lat")==0)
        {
            z=x+1;

            if (z<=y && argv[z][0]) 
            {
                    tx_site[0].lat = ReadBearing(argv[z]);
            }
        }
        if (strcmp(argv[x],"-lon")==0)
        {
            z=x+1;
            if (z<=y && argv[z][0]) 
            {
                tx_site[0].lon = ReadBearing(argv[z]);
                tx_site[0].lon*=-1;
                if (tx_site[0].lon<0.0)
                    tx_site[0].lon+=360.0;
                        }
        }
        //Switch to Path Profile Mode if Rx co-ords specified
        if (strcmp(argv[x],"-rla")==0)
        {
            z=x+1;

            if (z<=y && argv[z][0]) 
            {
                    ppa=1;
                    tx_site[1].lat = ReadBearing(argv[z]);
                    
            }
        }
        if (strcmp(argv[x],"-rlo")==0)
        {
            z=x+1;
            if (z<=y && argv[z][0]) 
            {
                tx_site[1].lon = ReadBearing(argv[z]);
                tx_site[1].lon*=-1;
                if (tx_site[1].lon<0.0)
                    tx_site[1].lon+=360.0;
                        }
        }
        
        if (strcmp(argv[x],"-txh")==0)
        {
            z=x+1;

            if (z<=y && argv[z][0] && argv[z][0]!='-')
            {
              sscanf(argv[z],"%f",&tx_site[0].alt);
             
            }
             txsites=1;
        }


        if (strcmp(argv[x],"-rxh")==0)
        {
            z=x+1;

            if (z<=y && argv[z][0] && argv[z][0]!='-')
            {
                sscanf(argv[z],"%lf",&altitudeLR);
                sscanf(argv[z],"%f",&tx_site[1].alt);
            }
        }


        if (strcmp(argv[x],"-f")==0)
        {
            z=x+1;

            if (z<=y && argv[z][0] && argv[z][0]!='-')
            {
                sscanf(argv[z],"%lf",&LR.frq_mhz);
            }
        }

        
        
        if (strcmp(argv[x],"-erp")==0)
        {
            z=x+1;

            if (z<=y && argv[z][0] && argv[z][0]!='-')
            {
                sscanf(argv[z],"%lf",&LR.erp);
            }
        }
               

                if (strcmp(argv[x],"-cl")==0)
        {
            z=x+1;

            if (z<=y && argv[z][0] && argv[z][0]!='-')
            {

                            sscanf(argv[z],"%d",&LR.radio_climate);

            }
        }
                 if (strcmp(argv[x],"-te")==0)
        {
            z=x+1;

            if (z<=y && argv[z][0] && argv[z][0]!='-')
            {

                            sscanf(argv[z],"%d",&ter);

                          switch (ter)
                                {
                                    case 1: // Water
                                    terdic = 80;
                                    tercon = 0.010;
                                    break;

                                    case 2: // Marsh
                                    terdic = 12;
                                    tercon = 0.007;
                                    break;

                                       case 3: // Farmland
                                    terdic = 15;
                                    tercon = 0.005;
                                    break;

                                    case 4: //Mountain
                                    terdic = 13;
                                    tercon = 0.002;
                                    break;
                                     case 5: //Desert
                                    terdic=13;
                                    tercon=0.002;
                                    break;
                                    case 6: //Urban
                                    terdic=5;
                                    tercon=0.001;
                                    break;
                        }
                          LR.eps_dielect = terdic;
                          LR.sgm_conductivity = tercon;

            }
    }

        if (strcmp(argv[x],"-terdic")==0)
        {
            z=x+1;

            if (z<=y && argv[z][0] && argv[z][0]!='-')
            {

                           sscanf(argv[z],"%lf",&terdic);

                          LR.eps_dielect = terdic;
                         
            }
        }

        if (strcmp(argv[x],"-tercon")==0)
        {
            z=x+1;

            if (z<=y && argv[z][0] && argv[z][0]!='-')
            {

                           sscanf(argv[z],"%lf",&tercon);

                          LR.sgm_conductivity = tercon;
                         
            }
        }

        if (strcmp(argv[x],"-hp")==0)
        {
            // Horizontal polarisation (0)
                    LR.pol = 0;
        }

        if (strcmp(argv[x],"-dbg")==0)
        {
                   debug=1;
        }

    
                ppd=(double)ippd;    /* pixels per degree (double)  */
                dpp=1.0/ppd;        /* degrees per pixel */
                mpi=ippd-1;        /* maximum pixel index per degree */


        /*UDT*/
   
        if (strcmp(argv[x],"-udt")==0)
        {
            z=x+1;

            if (z<=y && argv[z][0])
            {
           strncpy(udt_file,argv[z],253);
            }
        }
    
        /*Prop model*/
   
        if (strcmp(argv[x],"-pm")==0)
        {
            z=x+1;

            if (z<=y && argv[z][0])
            {
                sscanf(argv[z],"%d",&propmodel);
            }
        }
 
		// Prop model variant eg. urban/suburban
        if (strcmp(argv[x],"-pe")==0)
        {
            z=x+1;

            if (z<=y && argv[z][0])
            {
                sscanf(argv[z],"%d",&pmenv);
            }
        }
		
        //Knife edge diffraction
        if (strcmp(argv[x],"-ked")==0)
        {
            z=x+1;
            knifeedge=1;
        }
        
        //Windows friendly SDF filenames
        if (strcmp(argv[x],"-wf")==0)
        {
            z=x+1;
            winfiles=1;
        }
        
        //Normalise Path Profile chart
        if (strcmp(argv[x],"-ng")==0)
        {
            z=x+1;
            normalise=1;
        }
        
        //Halve the problem
        if (strcmp(argv[x],"-haf")==0)
        {
            z=x+1;
            if (z<=y && argv[z][0])
            {
                sscanf(argv[z],"%d",&haf);
            }
        }
        
    
    }
    
        /* ERROR DETECTION */
		if (tx_site[0].lat > 90 || tx_site[0].lat < -90)
				{
                      fprintf(stdout,"ERROR: Either the lat was missing or out of range!");
                      exit(0);

                }
                if (tx_site[0].lon > 360 || tx_site[0].lon < 0)
				{
                      fprintf(stdout,"ERROR: Either the lon was missing or out of range!");
                      exit(0);

                }
                if (LR.frq_mhz < 20 || LR.frq_mhz > 100000)
                {
                    fprintf(stdout,"ERROR: Either the Frequency was missing or out of range!");
                      exit(0);
                }
                if (LR.erp>5000000)
                {
                     fprintf(stdout,"ERROR: Power was out of range!");
                      exit(0);

                }
		if (LR.eps_dielect > 80 || LR.eps_dielect < 0.1)
                {
                     fprintf(stdout,"ERROR: Ground Dielectric value out of range!");
                      exit(0);

                }
		if (LR.sgm_conductivity > 0.01 || LR.sgm_conductivity < 0.000001)
                {
                     fprintf(stdout,"ERROR: Ground conductivity out of range!");
                      exit(0);

                }

                if (tx_site[0].alt < 0 || tx_site[0].alt > 60000)
                {
                     fprintf(stdout,"ERROR: Tx altitude above ground was too high: %f", tx_site[0].alt);
                      exit(0);
                }
                if (altitudeLR < 0 || altitudeLR > 60000)
                {
                     fprintf(stdout,"ERROR: Rx altitude above ground was too high!");
                      exit(0);
                }
                if (ippd < 300 || ippd > 3600){
                    fprintf(stdout,"ERROR: resolution out of range!");
                    exit(0);
                }
                if(contour_threshold < -200 || contour_threshold > 200)
                {
                    fprintf(stdout,"ERROR: Receiver threshold out of range (-200 / +200)");
                    exit(0);
                }
                if(propmodel>2 && propmodel<8 && LR.frq_mhz < 150){
                    fprintf(stdout,"ERROR: Frequency too low for Propagation model");
                    exit(0);
                }

         
                
    if (metric)
    {
        altitudeLR/=METERS_PER_FOOT;    /* 10ft * 0.3 = 3.3m */
        max_range/=KM_PER_MILE;        /* 10 / 1.6 = 7.5 */
        altitude/=METERS_PER_FOOT;    
        tx_site[0].alt/=METERS_PER_FOOT;    /* Feet to metres */
        tx_site[1].alt/=METERS_PER_FOOT;    /* Feet to metres */
        clutter/=METERS_PER_FOOT;        /* Feet to metres */
    }

    /* Ensure a trailing '/' is present in sdf_path */

    if (sdf_path[0])
    {
        x=strlen(sdf_path);

        if (sdf_path[x-1]!='/' && x!=0)
        {
            sdf_path[x]='/';
            sdf_path[x+1]=0;
        }
    }

    x=0;
    y=0;

    min_lat=70;
    max_lat=-70;

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
    
    if (ppa==1)
    {
        rxlat=(int)floor(tx_site[1].lat);
        rxlon=(int)floor(tx_site[1].lon);

        if (rxlat<min_lat)
            min_lat=rxlat;

        if (rxlat>max_lat)
            max_lat=rxlat;

        if (LonDiff(rxlon,min_lon)<0.0)
            min_lon=rxlon;

        if (LonDiff(rxlon,max_lon)>=0.0)
            max_lon=rxlon;
    }

    /* Load the required SDF files */

    LoadTopoData(max_lon, min_lon, max_lat, min_lat, winfiles);

    if (area_mode || topomap)
    {
        for (z=0; z<txsites && z<max_txsites; z++)
        {
            /* "Ball park" estimates used to load any additional
               SDF files required to conduct this analysis. */

            tx_range=sqrt(1.5*(tx_site[z].alt+GetElevation(tx_site[z])));

            if (LRmap)
                rx_range=sqrt(1.5*altitudeLR);
            else
                rx_range=sqrt(1.5*altitude);

            /* deg_range determines the maximum
               amount of topo data we read */

            deg_range=(tx_range+rx_range)/57.0;

            /* max_range regulates the size of the
               analysis.  A small, non-zero amount can
               be used to shrink the size of the analysis
               and limit the amount of topo data read by
               ss  A large number will increase the
               width of the analysis and the size of
               the map. */

            if (max_range==0.0)
                max_range=tx_range+rx_range;

            deg_range=max_range/57.0;

            // No more than 8 degs
            deg_limit=3.5;


            if (fabs(tx_site[z].lat)<70.0)
                deg_range_lon=deg_range/cos(DEG2RAD*tx_site[z].lat);
            else
                deg_range_lon=deg_range/cos(DEG2RAD*70.0);

            /* Correct for squares in degrees not being square in miles */

            if (deg_range>deg_limit)
                deg_range=deg_limit;

            if (deg_range_lon>deg_limit)
                deg_range_lon=deg_limit;

            nortRxHin=(int)floor(tx_site[z].lat-deg_range);
            nortRxHax=(int)floor(tx_site[z].lat+deg_range);

            west_min=(int)floor(tx_site[z].lon-deg_range_lon);

            while (west_min<0)
                west_min+=360;

            while (west_min>=360)
                west_min-=360;

            west_max=(int)floor(tx_site[z].lon+deg_range_lon);

            while (west_max<0)
                west_max+=360;

            while (west_max>=360)
                west_max-=360;

            if (nortRxHin<min_lat)
                min_lat=nortRxHin;

            if (nortRxHax>max_lat)
                max_lat=nortRxHax;

            if (LonDiff(west_min,min_lon)<0.0)
                min_lon=west_min;

            if (LonDiff(west_max,max_lon)>=0.0)
                max_lon=west_max;
        }

        /* Load any additional SDF files, if required */

        LoadTopoData(max_lon, min_lon, max_lat, min_lat, winfiles);
    }
 
    // UDT clutter
    LoadUDT   (udt_file);
 
    if(ppa==0){
        if (propmodel==2){
        PlotLOSMap(tx_site[0],altitudeLR,ano_filename);
        DoLOS(mapfile,geo,kml,ngs,tx_site,txsites);
        }
        else{
        // 90% of effort here
        PlotPropagation(tx_site[0],altitudeLR,ano_filename,propmodel,knifeedge,haf,pmenv);
        
        // Near field bugfix
         PutSignal(tx_site[0].lat,tx_site[0].lon,hottest);
        for(lat=tx_site[0].lat-0.002;lat<=tx_site[0].lat+0.002;lat=lat+0.0005){
            for(lon=tx_site[0].lon-0.002;lon<=tx_site[0].lon+0.002;lon=lon+0.0005){
                PutSignal(lat,lon,hottest);
            }
        }
                if (LR.erp==0.0)
                    DoPathLoss(mapfile,geo,kml,ngs,tx_site,txsites);
                else
                    if (dbm)
                        DoRxdPwr(mapfile,geo,kml,ngs,tx_site,txsites);
                    else
                        DoSigStr(mapfile,geo,kml,ngs,tx_site,txsites);
        }
        fprintf(stdout,"|%.5f",north);
        fprintf(stdout,"|%.5f",east);
        fprintf(stdout,"|%.5f",south);
        fprintf(stdout,"|%.5f|",west);
        
    }else{  
        strncpy(tx_site[0].name,"Tx",3);
        strncpy(tx_site[1].name,"Rx",3);
        PlotPath(tx_site[0],tx_site[1],1);
        PathReport(tx_site[0],tx_site[1],tx_site[0].filename,0,propmodel,pmenv);
        SeriesData(tx_site[0],tx_site[1],tx_site[0].filename,1,normalise);
    }    
    fflush(stdout);
    return 0;
}
