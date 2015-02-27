/*****************************************************************************
*  RF propagation models for Signal Server by Alex Farrant, CloudRF.com      *  
*                   									                     *
*  This program is free software; you can redistribute it and/or modify it   *
*  under the terms of the GNU General Public License as published by the     *
*  Free Software Foundation; either version 2 of the License or any later    *
*  version.								   								     *
* 									    									 *
*  This program is distributed in the hope that it will useful, but WITHOUT  *
*  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or     *
*  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License     *
*  for more details.													     *
*****************************************************************************/

#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <stdio.h>

using namespace std;

/*
Whilst every effort has been made to ensure the accuracy of the models, their accuracy is not guaranteed.
Finding a reputable paper to source these models from took a while. There was lots of bad copy-pasta out there.
A good paper: http://www.cl.cam.ac.uk/research/dtg/lce-pub/public/vsa23/VTC05_Empirical.pdf
*/

#define PI 3.14159265

/* Acute Angle from Rx point to an obstacle of height (opp) and distance (adj) */ 
double incidenceAngle(double opp, double adj){
	return atan2(opp,adj) * 180 / PI;
}

/* 
Knife edge diffraction:
This is based upon a recognised formula like Huygens, but trades thoroughness for increased speed 
which adds a proportional diffraction effect to obstacles.
*/
double ked(double freq, double elev[], double rxh, double dkm){
double obh,obd,rxobaoi=0,d,dipheight=25;

obh=0; // Obstacle height
obd=0; // Obstacle distance

dkm=dkm*1000; // KM to metres

	// walk along path 
	for(int n=2;n<(dkm/elev[1]);n++){
	
	d = (n-2)*elev[1]; // no of points * delta = km
	
	    //Find dip(s)
		if(elev[n]<(obh+dipheight)){
		
		// Angle from Rx point to obstacle
		rxobaoi = incidenceAngle((obh-(elev[n]+rxh)),d-obd);
		} else{
		// Line of sight or higher
		rxobaoi=0;
		}
	
		//note the highest point
		if(elev[n]>obh){
		obh=elev[n];
		obd=d; 
		}
		
	}
	
if(rxobaoi >= 0){
	return rxobaoi / (300/freq); // Diffraction angle divided by wavelength (m)
}else{
	return 0;
}

}

double HATApathLoss(float f,float TxH, float RxH, float d, int mode){
/*
HATA model for cellular planning
Frequency (MHz) 150 to 1500MHz
Base station height 30-200m
Mobile station height 1-10m
Distance 1-20km
modes 1 = URBAN, 2 = SUBURBAN, 3 = OPEN
*/

	if(f<150 || f>1500){
	    printf("Error: Hata model frequency range 150-1500MHz\n");
		return 0;
	}
	float lRxH = log10(11.75*RxH);
	float C_H = 3.2*lRxH*lRxH-4.97;
	float logf = log10(f);
	float L_u = 69.55 + 26.16*logf - 13.82*log10(TxH) - C_H + (44.9 - 6.55*log10(TxH))*log10(d);

	if(!mode || mode==1){
	return L_u; //URBAN
	}
	
	if(mode==2){ //SUBURBAN
	float logf_28 = log10(f/28);
	return L_u - 2*logf_28*logf_28 - 5.4;
	}

	if(mode==3){ //OPEN
	return L_u - 4.78*logf*logf + 18.33*logf - 40.94;
	}

	return 0;
}

double COST231pathLoss(float f,float TxH, float RxH, float d, int mode){
/*
COST231 extension to HATA model
Frequency 1500 to 2000MHz
TxH = Base station height 30 to 200m
RxH = Mobile station height 1 to 10m
Distance 1-20km
modes 1 = URBAN, 2 = SUBURBAN, 3 = OPEN
*/
	if(f<1500 || f>2000){
	    printf("Error: COST231 Hata model frequency range 1500-2000MHz\n");
		return 0;
	}

	int C = 3; // 3dB for Urban
	if(mode==2){
		C = 0; // Suburban, rural
	}
	if(mode==3){
		C = -3; // Suburban, rural
	}
	float lRxH = log10(11.75*RxH);
	float C_H = 3.2*lRxH*lRxH-4.97;
	float logf = log10(f);
    double dbloss = 46.3 + (33.9 * logf) - (13.82 * log10(TxH)) - C_H + (44.9 - 6.55 * log10(TxH)) * log10(d) + C;  
	return dbloss;
}

double SUIpathLoss(float f,float TxH, float RxH, float d, int mode){
	/*
	f = Frequency (MHz)
	TxH =  Transmitter height (m) 
	RxH = Receiver height (m)
	d = distance (km)
	mode 1 = Hilly + trees
	mode 2 = Flat + trees OR hilly + light foliage
	mode 3 = Flat + light foliage
	*/
	d=d*1000; // km to m
	if(f<1900 || f>11000){
	    printf("Error: SUI model frequency range 1.9-11GHz\n");
		return 0;
	}
	// Terrain mode A is default
	double a = 4.6;
	double b = 0.0075;
	double c = 12.6; 
	double s = 10.6; 
	int XhCF = -10.8;
	if(mode==2){
		a=4.0;
		b=0.0065;
		c=17.1;
		s=9.6;
	}
	if(mode==3){
		a=3.6;
		b=0.005;
		c=20;
		s=8.2;
		XhCF = -20;
	}
	double d0 = 100; 
	double A = 20 * log10((4*M_PI*d0)/(300/f));
	double y = (a - b * TxH) + (c/TxH);
	double Xf = 6 * log10(f/2000);
	double Xh = XhCF * log10(RxH/2);
	return A + (10*y*log10(d/d0)) + Xf + Xh + s;
}

double ECC33pathLoss(float f,float TxH, float RxH, float d, int mode){
	// MHz to GHz
	f=f/1000;

	double Gr = 0.759 * RxH - 1.862; // Big city (1)
	// PL = Afs + Abm - Gb - Gr
	double Afs = 92.4 + 20 * log10(d) + 20 * log10(f);
	double Abm = 20.41 + 9.83 * log10(d) + 7.894 * log10(f) + 9.56 * (log10(f) * log10(f));
	double Gb = log10(TxH/200) * (13.958 + 5.8 * (log10(d) * log10(d)));
	if(mode>1){ // Medium city
		Gr = (42.57 + 13.7 * log10(f)) * (log10(RxH) - 0.585);
	}
	return Afs+Abm-Gb-Gr;
}

double EricssonpathLoss(float f,float TxH, float RxH, float d, int mode){
	double a0=36.2, a1=30.2, a2=12, a3=0.1;
	if(mode==2){ // Med loss
		a0=43.2;
		a1=68.93;
	}
	if(mode==1){ // High loss
		a0=45.95;
		a1=100.6;
	}
	double g1 = (11.75 * RxH) * (11.75 * RxH); 
	double g2 = (44.49 * log10(f)) - 4.78 * ((log10(f) * log10(f))); 
	double PL = a0+ a1 * log10(d) + a2 * log10(TxH) + a3 * log10(TxH) * log10(d) - (3.2 * log10(g1)) + g2;
	return PL;
}

double FSPLpathLoss(float f, float d){
/*
Free Space Path Loss (ITU-R P.525) model
Frequency: Any
Distance: Any 
*/
	//MHz to GHz
	f = f / 1000;
	double dbloss = (20 * log10(d)) + (20 * log10(f)) + 92.45;
	return dbloss;
}
/*
int main(int argc, char* argv[]){
	if(argc<5){
		printf("Need freq,TxH,RxH,dist,terr\n");
		return 0;
	}
	int dis, ter;
	double frq, TxH, RxH;
	
	sscanf(argv[1],"%lf",&frq);
	sscanf(argv[2],"%lf",&TxH);
	sscanf(argv[3],"%lf",&RxH);
	sscanf(argv[4],"%d",&dis);
	sscanf(argv[5],"%d",&ter);
	// ALL are freq in MHz and distances in metres
	printf("FSPL: %.2f dB\n",FSPLpathLoss(frq,dis));
	printf("HATA (%d): %.2f dB\n",ter,HATApathLoss(frq,TxH,RxH,dis,ter));
	printf("COST-HATA (%d): %.2f dB\n",ter,COST231pathLoss(frq,TxH,RxH,dis,ter));
	printf("SUI (%d): %.2f dB\n",ter,SUIpathLoss(frq,TxH,RxH,dis,ter));
	printf("ECC33 (%d): %.2f dB\n",ter,ECC33pathLoss(frq,TxH,RxH,dis,ter));
	printf("Ericsson (%d): %.2f dB\n",ter,EricssonpathLoss(frq,TxH,RxH,dis,ter));
}
*/
