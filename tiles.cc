#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include "tiles.hh"
#include "common.h"

#define MAX_LINE 25000

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

	if (debug)
		fprintf(stderr,"Pixels loaded: %zu/%d\n",loaded,tile->width*tile->height);

	/* All done, close the LIDAR file */
	fclose(fd);

	return 0;
}