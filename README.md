# Signal Server
Server optimised SPLAT! by Alex Farrant, M6ZUJ. 

SPLAT! Project started in 1997 by John A. Magliacane, KD2BD

This server application will generate RF coverage predictions, producing either 2D profile plots (Point-to-Point) or 360 degree polar plots in WGS-84 projection.
For detailed information and historical reference data related to this project see the SPLAT! documentation. Propagation models added to this project have been sourced from reputable academic sources and all efforts have been taken to ensure their accurate implementation. Not all models are ITU ratified and you use them entirely at your own risk.

The accuracy of the output is directly proportional to the accuracy of the inputs and the time taken defining and validating them. 


## Requirements
* Linux
* GCC,G++
* Multicore CPU (optional)
* ~2GB Memory
* SRTM terrain tiles

Signal Server is a very resource intensive multicore application. Only publish it for common use if you know what you are doing and you are advised to wrap it with another script to perform input validation.

Additional programs/scripts will be required to prepare inputs such as .sdf tiles (srtm2sdf.c), 3D antenna patterns (.ant) and user defined clutter (.udt) or manipulate the bitmap output (.ppm). More information can be found in the SPLAT! project.

## Installation
```sh
$ make
```
## Parameters
		 -- Signal Server 2.72 --
	Set for 64 tiles at 1200 pixels/degree

     -sdf Directory containing .sdf tiles
     -lid LIDAR ASCII tile with WGS84 bounds (Dimensions defined in file metadata)
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
     -res Pixels per tile. 300/600/1200/3600 (Optional. LIDAR res is defined within the tile)
     -t Terrain background
     -pm Prop model. 1: ITM, 2: LOS, 3: Hata, 4: ECC33,
     		5: SUI, 6: COST-Hata, 7: FSPL, 8: ITWOM, 9: Ericsson
     -pe Prop model mode: 1=Urban,2=Suburban,3=Rural
     -ked Knife edge diffraction (Default for ITM)
     -ng Normalise Path Profile graph
     -haf Halve 1 or 2 (optional)
     -nothreads Turn off threaded processing (optional)


## Examples
	
### 90m resolution	
- INPUTS: 900MHz tower at 25m AGL with 5W ERP, 30km radius
- OUTPUTS: 1200 resolution, 30km radius, -90dBm receiver threshold, Longley Rice model
- ./signalserver -sdf /data/SRTM3 -lat 51.849 -lon -2.2299 -txh 25 -f 900 -erp 5 -rxh 2 -rt -90 -dbm -m -o test1 -R 30 -res 1200 -pm 1

### 30m resolution
- INPUTS: 450MHz tower at 25f AGL with 20W ERP, 10km radius
- OUTPUTS: 3600 resolution, 30km radius, 10dBuV receiver threshold, Hata model
- ./signalserverHD -sdf /data/SRTM1 -lat 51.849 -lon -2.2299 -txh 25 -f 450 -erp 20 -rxh 2 -rt 10 -o test2 -R 10 -res 3600 -pm 3

### 2m resolution (LIDAR)
- INPUTS: 1800MHz tower at 15m AGL with 1W ERP, 1 km radius
- OUTPUTS: 2m LIDAR resolution, 5km radius, -90dBm receiver threshold, Longley Rice model
- ./signalserverLIDAR -lid /data/LIDAR/Gloucester_2m.asc -lat 51.849 -lon -2.2299 -txh 15 -f 1800 -erp 1 -rxh 2 -rt -90 -dbm -m -o test3 -R 1 -pm 1
