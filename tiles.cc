#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <math.h>
#include "tiles.hh"
#include "common.h"

#define MAX_LINE 25000

/* Computes the distance between two long/lat points */
double haversine_formulaz(double th1, double ph1, double th2, double ph2)
{
	#define TO_RAD (3.1415926536 / 180)
	int R = 6371;
	double dx, dy, dz;
	ph1 -= ph2;
	ph1 *= TO_RAD, th1 *= TO_RAD, th2 *= TO_RAD;
 
	dz = sin(th1) - sin(th2);
	dx = cos(ph1) * cos(th1) - cos(th2);
	dy = sin(ph1) * cos(th1);
	return asin(sqrt(dx * dx + dy * dy + dz * dz) / 2) * 2 * R;
}

int tile_load_lidar(tile_t *tile, char *filename){
	FILE *fd;
	char line[MAX_LINE];
	int nextval;
	char *pch;

	/* Clear the tile data */
	memset(tile,0x00,sizeof(tile_t));

	/* Open the file handle and return on error */
	if ( (fd = fopen(filename,"r")) == NULL )
		return errno;

	/* This is where we read the header data */
	/* The string is split for readability but is parsed as a block */
	if( fscanf(fd,"%*s %d\n" "%*s %d\n" "%*s %lf\n" "%*s %lf\n" "%*s %lf\n" "%*s %d\n",&tile->width,&tile->height,&tile->xll,&tile->yll,&tile->cellsize,&tile->nodata) != 6 ){
		fclose(fd);
		return -1;
	}

	tile->datastart = ftell(fd);

	if(debug){
		fprintf(stderr,"w:%d h:%d s:%lf\n", tile->width, tile->height, tile->cellsize);
		fflush(stderr);
	}

	/* Set the filename */
	tile->filename = strdup(filename);

	/* Perform xur calcs */
	// Degrees with GDAL option: -co "FORCE_CELLSIZE=YES"
	tile->xur = tile->xll+(tile->cellsize*tile->width);
	tile->yur = tile->yll+(tile->cellsize*tile->height);

	if (tile->xur > eastoffset)
		eastoffset = tile->xur;
	if (tile->xll < westoffset)
		westoffset = tile->xll;

	// if (debug)
	// 	fprintf(stderr,"%d, %d, %.7f, %.7f, %.7f, %.7f, %.7f\n",width,height,xll,yll,cellsize,yur,xur);

	// Greenwich straddling hack
	if (tile->xll <= 0 && tile->xur > 0) {
		tile->xll = (tile->xur - tile->xll); // full width
		tile->xur = 0.0; // budge it along so it's west of greenwich
		delta = eastoffset; // add to Tx longitude later
	} else {
		// Transform WGS84 longitudes into 'west' values as society finishes east of Greenwich ;)
		if (tile->xll >= 0)
			tile->xll = 360-tile->xll;
		if(tile->xur >= 0)
			tile->xur = 360-tile->xur;
		if(tile->xll < 0)
			tile->xll = tile->xll * -1;
		if(tile->xur < 0)
			tile->xur = tile->xur * -1;
	}

	if (debug)
		fprintf(stderr, "POST yll %.7f yur %.7f xur %.7f xll %.7f delta %.6f\n", tile->yll, tile->yur, tile->xur, tile->xll, delta);

	/* Read the actual tile data */
	/* Allocate the array for the lidar data */
	if ( (tile->data = (int*) calloc(tile->width * tile->height, sizeof(int))) == NULL ) {
		fclose(fd);
		free(tile->filename);
		return ENOMEM;
	}

	size_t loaded = 0;
	for (size_t h = 0; h < tile->height; h++) {
		if (fgets(line, MAX_LINE, fd) != NULL) {
			pch = strtok(line, " "); // split line into values
			for (size_t w = 0; w < tile->width && pch != NULL; w++) {
				/* If the data is less than a *magic* minimum, normalize it to zero */
				nextval = atoi(pch);
				if (nextval <= -9999)
					nextval = 0;
				tile->data[h*tile->width + w] = nextval;
				loaded++;
				if ( nextval > tile->max_el )
					tile->max_el = nextval;
				if ( nextval < tile->min_el )
					tile->min_el = nextval;
				pch = strtok(NULL, " ");
			}//while
		} else {
			fprintf(stderr, "LIDAR error @ h %zu file %s\n", h, filename);
		}//if
	}

	double current_res_km = haversine_formulaz(tile->max_north, tile->max_west, tile->max_north, tile->min_west);
	tile->resolution = (int) ceil((current_res_km/MAX(tile->width,tile->height))*1000);

	if (debug)
		fprintf(stderr,"Pixels loaded: %zu/%d\n",loaded,tile->width*tile->height);

	/* All done, close the LIDAR file */
	fclose(fd);

	return 0;
}

/*
 *	A positive scale will _increase_ the size of the data
 */
int tile_rescale(tile_t *tile, float scale){
	int *new_data;
	size_t skip_count = 1;
	size_t copy_count = 1;

	size_t new_height = tile->height * scale;
	size_t new_width = tile->width * scale;

	/* Allocate the array for the lidar data */
	if ( (new_data = (int*) calloc(new_height * new_width, sizeof(int))) == NULL ) {
		return ENOMEM;
	}

	tile->max_el = -32768;
	tile->min_el = 32768;

	/* Making the tile data smaller */
	if (scale < 0) {
		skip_count = 1 / scale;
	} else {
		copy_count = (size_t) scale;
	}

	fprintf(stderr,"Skip: %zu Copy: %zu\n", skip_count, copy_count);

	

	/* Nearest neighbour normalization. For each subsample of the original, simply
	 * assign the value in the top left to the new pixel */
	if (scale < 0){
		for (size_t x = 0, i = 0; i < new_width; x += skip_count, i++) {
			for (size_t y = 0, j = 0; j < new_height; y += skip_count, j++){
				new_data[i*new_width+j] = tile->data[x*tile->width+y];
				/* Update local min / max values */
				if (tile->data[x * tile->width + y] > tile->max_el)
					tile->max_el = tile->data[x * tile->width + y];
				if (tile->data[x * tile->width + y] < tile->min_el)
					tile->min_el = tile->data[x * tile->width + y];
			}
		}
	}else{
		for (size_t x = 0; x < tile->width; x++) {
			for (size_t y = 0; y < tile->height; y++){
				/* These are for scaling up the data */
				for (size_t copy_x = 0; copy_x < copy_count; copy_x++) {
					for (size_t copy_y = 0; copy_y < copy_count; copy_y++) {
						size_t new_x = (x * skip_count) + copy_x;
						size_t new_y = (y * skip_count) + copy_y;
						new_data[ new_x * new_width + new_y ] = tile->data[x * tile->width + y];
						// new_data[(x + copy_x) * new_width + (y + copy_y)] = tile->data[x * tile->width + y];
					}
				}
				/* Update local min / max values */
				if (tile->data[x * tile->width + y] > tile->max_el)
					tile->max_el = tile->data[x * tile->width + y];
				if (tile->data[x * tile->width + y] < tile->min_el)
					tile->min_el = tile->data[x * tile->width + y];
			}
		}
	}

	/* Update the date in the tile */
	free(tile->data);
	tile->data = new_data;

	/* Update the height and width values */
	tile->height = new_height;
	tile->width = new_width;

	return 0;
}
