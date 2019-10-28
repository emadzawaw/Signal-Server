#include "config.h"

#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>

void print_usage() {
	printf("Version: %s %.2f (Built for %d DEM tiles at %d pixels)\n", ss_name, version,MAXPAGES, IPPD);
	printf("License: GNU General Public License (GPL) version 2\n\n");
	printf("Radio propagation simulator by Alex Farrant QCVS, 2E0TDW\n");
	printf("Contributions by Darcy Buskermolen VA7DBI\n");
	printf("Based upon SPLAT! by John Magliacane, KD2BD\n\n");
	printf("Usage: %s [data options] [input options] [output options] -o outputfile\n\n",argv[0]);
	printf("Data Files:\n");
	printf("     --sdf           Directory containing SRTM derived .sdf DEM tiles\n");
	printf("     --lidar         ASCII grid tile (LIDAR) with dimensions and resolution defined in header\n");
	printf("     --pointclutter  User defined point clutter as decimal co-ordinates: 'latitude,longitude,height'\n");
	printf("     --modisclutter  MODIS 17-class wide area clutter in ASCII grid format\n");
	printf("Input:\n");
	printf("     --latitude      Transmitter Latitude  (decimal degrees -70/+70)\n");
	printf("     --longitude     Transmitter Longitude (decimal degrees -180/+180)\n");
	printf("     --txheight      Transmitter Height above ground  (feet/meters)\n");
	printf("     --frequency     Transmiter Frequency (MHz) 20MHz to 100GHz (LOS after 20GHz)\n");
	printf("     --erp           Tranmsitter Effective Radiated Power including Tx+Rx gain (Watts)\n");
	printf("Optional Inputs:\n");
	printf("     --rxlatitude    Receiver Latitude for Point-To-Point  (decimal degrees -70/+70)\n");
	printf("     --rxlogitude    Receiver Longitude for Point-ToPoint (decimal degrees -180/+180)\n");
	printf("     --rxheight      Receiver Height (feet/meters Default=0.1)\n");
	printf("     --rxgain        Receiver gain in dBi (text report)\n");
	printf("     --randomclutter Random ground clutter (feet/meters)\n");
	printf("     --dieletric     Terrain dielectric value 2-80 (optional)\n");
	printf("     --conductivity  Terrain conductivity 0.01-0.0001 (optional)\n");
	printf("Mondifiers:\n");
	printf("     --metric        Metric units of measurement\n");
	printf("     --horizontal    Horizontal Polarisation (default=vertical)\n");
	printf("     --terrain       Terrain code 1-6 (optional)\n");
	printf("     --resample      Resample Lidar input to specified resolution in meters (optional)\n");
	printf("     --climate       Climate code 1-6 (optional)\n");
	printf("     --reliability   Reliability for ITM model 50 to 99 (optional)\n");
	printf("Output:\n");
	printf("     --output        Output filename. Required. \n");
	printf("     --dbm           Plot Received signal power instead of field strength\n");
	printf("     --rxthreshold   Reciever Threshold (dB / dBm / dBuV/m)\n");
	printf("     --radius        Radius (miles/kilometers)\n");
	printf("     --resolution    Pixels per tile. 300/600/1200/3600 (Optional. LIDAR res is within the tile)\n");
	printf("     --model         Propagation model. 1: ITM, 2: LOS, 3: Hata, 4: ECC33,5: SUI, 6: COST-Hata, \n");
        printf("                               7: FSPL, 8: ITWOM, 9: Ericsson, 10: Plane earth, 11: Egli VHF/UHF\n");
	printf("     --mode          Propagation model mode: 1=Urban,2=Suburban,3=Rural\n");
	printf("     --kinfeedge     Knife edge diffraction (Already on for ITM)\n");
	printf("Debugging:\n");
	printf("     --grayscale     Terrain greyscale background\n");
	printf("     --verbose       Verbose debug messages\n");
	printf("     --normalise     Normalise Path Profile graph\n");
	printf("     --half          Halve 1 or 2 (optional)\n");
	printf("     --nothreads     Turn off threaded processing\n");
	return 1;
}

int main(int argc, char *argv[]) {
    int opt= 0;

    //Specifying the expected options
    //The two options l and b expect numbers as argument
    static struct option long_options[] = {
        {"sdf",		required_argument,	0, 'a'},
        {"lidar",	required_argument,	0, 'p'},
        {"pointclutter",required_argument,      0, 'p'},
	{"modisclutter",required_argument,      0, 'p'},
	{"latitude",	required_argument,      0, 'p'},
	{"longitude",	required_argument,      0, 'p'},
	{"txheight",	required_argument,      0, 'p'},
	{"frequency",	required_argument,      0, 'p'},
	{"erp",		required_argument,      0, 'p'},
	{"rxlatitude",	required_argument,      0, 'p'},
	{"rxlogitude",  required_argument,      0, 'p'},
	{"rxheight",	required_argument,      0, 'p'},
	{"rxgain",	required_argument,      0, 'p'},
	{"randomclutter", required_argument,	0, 'p'},
	{"dieletric",	required_argument,      0, 'p'},
	{"conductivity", required_argument,	0, 'p'},
	{"metric",	no_argument,		0, 'p'},
	{"horizontal",	no_argument,            0, 'p'},
	{"terrain",	required_argument,      0, 'p'},
	{"resample",	required_argument,      0, 'p'},
	{"climate",	required_argument,      0, 'p'},
	{"reliability", required_argument,      0, 'p'},
	{"output",	required_argument,      0, 'p'},
	{"dbm",		no_argument,            0, 'p'},
	{"rxthreshold", required_argument,      0, 'p'},
	{"radius",	required_argument,      0, 'p'},
	{"resolution",	required_argument,      0, 'p'},
	{"model",	required_argument,      0, 'p'},
	{"mode",	required_argument,      0, 'p'},
	{"kinfeedge",	no_argument,            0, 'p'},
	{"grayscale",	no_argument,            0, 'p'},
	{"verbose",	no_argument,            0, 'p'},
	("normalise",	no_argument,            0, 'p'},
	{"half",	required_argument,      0, 'p'},
	{"nothreads",	no_argument,            0, 'p'},
        {0,		0,                 	0,  0 }
    };

    int long_index =0;
    while ((opt = getopt_long(argc, argv,"apl:b:", 
                   long_options, &long_index )) != -1) {
        switch (opt) {
             case 'a':
		area = 0;
                break;
             case 'p':
		perimeter = 0;
                break;
             case 'l':
		length = atoi(optarg); 
                break;
             case 'b':
		breadth = atoi(optarg);
                break;
             default: 
		print_usage(); 
                exit(EXIT_FAILURE);
        }
    }
    if (length == -1 || breadth ==-1) {
        print_usage();
        exit(EXIT_FAILURE);
    }
    return 0;
}
