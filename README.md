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

		 -- Signal Server 2.41 --
	Compiled for 64 tiles at 1200 pixels/degree

     -d Directory containing .sdf tiles
     -lat Tx Latitude (decimal degrees) -70/+70
     -lon Tx Longitude (decimal degrees) -180/+180
     -txh Tx Height (above ground)
     -rla (Optional) Rx Latitude for PPA (decimal degrees) -70/+70
     -rlo (Optional) Rx Longitude for PPA (decimal degrees) -180/+180
     -f Tx Frequency (MHz) 20MHz to 100GHz (LOS after 20GHz)
     -erp Tx Effective Radiated Power (Watts)
     -rxh Rx Height(s) (optional. Default=0.1)
     -rt Rx Threshold (dB / dBm / dBuV/m)
     -hp Horizontal Polarisation (default=vertical)
     -gc Ground clutter (feet/meters)
     -udt User defined terrain filename
     -dbm Plot Rxd signal power instead of field strength
     -m Metric units of measurement
     -te Terrain code 1-6 (optional)
     -terdic Terrain dielectric value 2-80 (optional)
     -tercon Terrain conductivity 0.01-0.0001 (optional)
     -cl Climate code 1-6 (optional)
     -o Filename. Required. 
     -R Radius (miles/kilometers)
     -res Pixels per degree. 300/600/1200(default)/3600 (optional)
     -t Terrain background
     -pm Prop model. 1: ITM, 2: LOS, 3: Hata, 4: ECC33,
     		5: SUI, 6: COST-Hata, 7: FSPL, 8: ITWOM, 9: Ericsson
     -pe Prop model mode: 1=Urban,2=Suburban,3=Rural
     -ked Knife edge diffraction (Default for ITM)
     -ng Normalise Path Profile graph
     -haf Halve 1 or 2 (optional)

