#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>

#include "common.h"
#include "main.hh"

int LoadSDF_SDF(char *name, int winfiles)
{
	/* This function reads uncompressed ss Data Files (.sdf)
	   containing digital elevation model data into memory.
	   Elevation data, maximum and minimum elevations, and
	   quadrangle limits are stored in the first available
	   dem[] structure. */

	int x, y, data = 0, indx, minlat, minlon, maxlat, maxlon, j;
	char found, free_page = 0, line[20], jline[20], sdf_file[255],
	    path_plus_name[255], *junk = NULL;

	FILE *fd;

	for (x = 0; name[x] != '.' && name[x] != 0 && x < 250; x++)
		sdf_file[x] = name[x];

	sdf_file[x] = 0;

	/* Parse filename for minimum latitude and longitude values */
	if (winfiles == 1) {
		sscanf(sdf_file, "%d=%d=%d=%d", &minlat, &maxlat, &minlon,
		       &maxlon);
	} else {
		sscanf(sdf_file, "%d:%d:%d:%d", &minlat, &maxlat, &minlon,
		       &maxlon);
	}

	sdf_file[x] = '.';
	sdf_file[x + 1] = 's';
	sdf_file[x + 2] = 'd';
	sdf_file[x + 3] = 'f';
	sdf_file[x + 4] = 0;

	/* Is it already in memory? */

	for (indx = 0, found = 0; indx < MAXPAGES && found == 0; indx++) {
		if (minlat == dem[indx].min_north
		    && minlon == dem[indx].min_west
		    && maxlat == dem[indx].max_north
		    && maxlon == dem[indx].max_west)
			found = 1;
	}

	/* Is room available to load it? */

	if (found == 0) {
		for (indx = 0, free_page = 0; indx < MAXPAGES && free_page == 0;
		     indx++)
			if (dem[indx].max_north == -90)
				free_page = 1;
	}

	indx--;

	if (free_page && found == 0 && indx >= 0 && indx < MAXPAGES) {
		/* Search for SDF file in current working directory first */

		strncpy(path_plus_name, sdf_file, 255);

		fd = fopen(path_plus_name, "rb");

		if (fd == NULL) {
			/* Next, try loading SDF file from path specified
			   in $HOME/.ss_path file or by -d argument */

			strncpy(path_plus_name, sdf_path, 255);
			strncat(path_plus_name, sdf_file, 255);
			fd = fopen(path_plus_name, "rb");
		}

		if (fd != NULL) {
			if (debug == 1) {
				fprintf(stdout,
					"Loading \"%s\" into page %d...",
					path_plus_name, indx + 1);
				fflush(stdout);
			}

			if (fgets(line, 19, fd) != NULL) {
				sscanf(line, "%d", &dem[indx].max_west);
			}

			if (fgets(line, 19, fd) != NULL) {
				sscanf(line, "%d", &dem[indx].min_north);
			}

			if (fgets(line, 19, fd) != NULL) {
				sscanf(line, "%d", &dem[indx].min_west);
			}

			if (fgets(line, 19, fd) != NULL) {
				sscanf(line, "%d", &dem[indx].max_north);
			}
			/*
			   Here X lines of DEM will be read until IPPD is reached.
			   Each .sdf tile contains 1200x1200 = 1.44M 'points'
			   Each point is sampled for 1200 resolution!
			 */
			for (x = 0; x < ippd; x++) {
				for (y = 0; y < ippd; y++) {

					for (j = 0; j < jgets; j++) {
						junk = fgets(jline, 19, fd);
					}

					if (fgets(line, 19, fd) != NULL) {
						data = atoi(line);
					}

					dem[indx].data[x][y] = data;
					dem[indx].signal[x][y] = 0;
					dem[indx].mask[x][y] = 0;

					if (data > dem[indx].max_el)
						dem[indx].max_el = data;

					if (data < dem[indx].min_el)
						dem[indx].min_el = data;

				}

				if (ippd == 600) {
					for (j = 0; j < IPPD; j++) {
						junk = fgets(jline, 19, fd);
					}
				}
				if (ippd == 300) {
					for (j = 0; j < IPPD; j++) {
						junk = fgets(jline, 19, fd);
						junk = fgets(jline, 19, fd);
						junk = fgets(jline, 19, fd);

					}
				}
			}

			fclose(fd);

			if (dem[indx].min_el < min_elevation)
				min_elevation = dem[indx].min_el;

			if (dem[indx].max_el > max_elevation)
				max_elevation = dem[indx].max_el;

			if (max_north == -90)
				max_north = dem[indx].max_north;

			else if (dem[indx].max_north > max_north)
				max_north = dem[indx].max_north;

			if (min_north == 90)
				min_north = dem[indx].min_north;

			else if (dem[indx].min_north < min_north)
				min_north = dem[indx].min_north;

			if (max_west == -1)
				max_west = dem[indx].max_west;

			else {
				if (abs(dem[indx].max_west - max_west) < 180) {
					if (dem[indx].max_west > max_west)
						max_west = dem[indx].max_west;
				}

				else {
					if (dem[indx].max_west < max_west)
						max_west = dem[indx].max_west;
				}
			}

			if (min_west == 360)
				min_west = dem[indx].min_west;

			else {
				if (fabs(dem[indx].min_west - min_west) < 180.0) {
					if (dem[indx].min_west < min_west)
						min_west = dem[indx].min_west;
				}

				else {
					if (dem[indx].min_west > min_west)
						min_west = dem[indx].min_west;
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

	int x, y, indx, minlat, minlon, maxlat, maxlon;
	char found, free_page = 0;
	int return_value = -1;

	return_value = LoadSDF_SDF(name, winfiles);

	/* If neither format can be found, then assume the area is water. */

	if (return_value == 0 || return_value == -1) {

		if (winfiles == 1) {
			sscanf(name, "%d=%d=%d=%d", &minlat, &maxlat, &minlon,
			       &maxlon);
		} else {
			sscanf(name, "%d:%d:%d:%d", &minlat, &maxlat, &minlon,
			       &maxlon);
		}
		/* Is it already in memory? */

		for (indx = 0, found = 0; indx < MAXPAGES && found == 0; indx++) {
			if (minlat == dem[indx].min_north
			    && minlon == dem[indx].min_west
			    && maxlat == dem[indx].max_north
			    && maxlon == dem[indx].max_west)
				found = 1;
		}

		/* Is room available to load it? */

		if (found == 0) {
			for (indx = 0, free_page = 0;
			     indx < MAXPAGES && free_page == 0; indx++)
				if (dem[indx].max_north == -90)
					free_page = 1;
		}

		indx--;

		if (free_page && found == 0 && indx >= 0 && indx < MAXPAGES) {
			if (debug == 1) {
				fprintf(stdout,
					"Region  \"%s\" assumed as sea-level into page %d...",
					name, indx + 1);
				fflush(stdout);
			}

			dem[indx].max_west = maxlon;
			dem[indx].min_north = minlat;
			dem[indx].min_west = minlon;
			dem[indx].max_north = maxlat;

			/* Fill DEM with sea-level topography */

			for (x = 0; x < ippd; x++)
				for (y = 0; y < ippd; y++) {
					dem[indx].data[x][y] = 0;
					dem[indx].signal[x][y] = 0;
					dem[indx].mask[x][y] = 0;

					if (dem[indx].min_el > 0)
						dem[indx].min_el = 0;
				}

			if (dem[indx].min_el < min_elevation)
				min_elevation = dem[indx].min_el;

			if (dem[indx].max_el > max_elevation)
				max_elevation = dem[indx].max_el;

			if (max_north == -90)
				max_north = dem[indx].max_north;

			else if (dem[indx].max_north > max_north)
				max_north = dem[indx].max_north;

			if (min_north == 90)
				min_north = dem[indx].min_north;

			else if (dem[indx].min_north < min_north)
				min_north = dem[indx].min_north;

			if (max_west == -1)
				max_west = dem[indx].max_west;

			else {
				if (abs(dem[indx].max_west - max_west) < 180) {
					if (dem[indx].max_west > max_west)
						max_west = dem[indx].max_west;
				}

				else {
					if (dem[indx].max_west < max_west)
						max_west = dem[indx].max_west;
				}
			}

			if (min_west == 360)
				min_west = dem[indx].min_west;

			else {
				if (abs(dem[indx].min_west - min_west) < 180) {
					if (dem[indx].min_west < min_west)
						min_west = dem[indx].min_west;
				}

				else {
					if (dem[indx].min_west > min_west)
						min_west = dem[indx].min_west;
				}
			}

			return_value = 1;
		}
	}

	return return_value;
}

void LoadPAT(char *filename)
{
	/* This function reads and processes antenna pattern (.az
	   and .el) files that correspond in name to previously
	   loaded ss .lrp files.  */

	int a, b, w, x, y, z, last_index, next_index, span;
	char string[255], azfile[255], elfile[255], *pointer = NULL;
	float az, xx, elevation, amplitude, rotation, valid1, valid2,
	    delta, azimuth[361], azimuth_pattern[361], el_pattern[10001],
	    elevation_pattern[361][1001], slant_angle[361], tilt,
	    mechanical_tilt = 0.0, tilt_azimuth, tilt_increment, sum;
	FILE *fd = NULL;
	unsigned char read_count[10001];

	for (x = 0; filename[x] != '.' && filename[x] != 0 && x < 250; x++) {
		azfile[x] = filename[x];
		elfile[x] = filename[x];
	}

	azfile[x] = '.';
	azfile[x + 1] = 'a';
	azfile[x + 2] = 'z';
	azfile[x + 3] = 0;

	elfile[x] = '.';
	elfile[x + 1] = 'e';
	elfile[x + 2] = 'l';
	elfile[x + 3] = 0;

	rotation = 0.0;

	got_azimuth_pattern = 0;
	got_elevation_pattern = 0;

	/* Load .az antenna pattern file */

	fd = fopen(azfile, "r");

	if (fd != NULL) {
		/* Clear azimuth pattern array */

		for (x = 0; x <= 360; x++) {
			azimuth[x] = 0.0;
			read_count[x] = 0;
		}

		/* Read azimuth pattern rotation
		   in degrees measured clockwise
		   from true North. */

		if (fgets(string, 254, fd) == NULL) {
			//fprintf(stdout,"Azimuth read error\n");
			//exit(0);
		}
		pointer = strchr(string, ';');

		if (pointer != NULL)
			*pointer = 0;

		sscanf(string, "%f", &rotation);

		/* Read azimuth (degrees) and corresponding
		   normalized field radiation pattern amplitude
		   (0.0 to 1.0) until EOF is reached. */

		if (fgets(string, 254, fd) == NULL) {
			//fprintf(stdout,"Azimuth read error\n");
			//exit(0);
		}
		pointer = strchr(string, ';');

		if (pointer != NULL)
			*pointer = 0;

		sscanf(string, "%f %f", &az, &amplitude);

		do {
			x = (int)rintf(az);

			if (x >= 0 && x <= 360 && fd != NULL) {
				azimuth[x] += amplitude;
				read_count[x]++;
			}

			if (fgets(string, 254, fd) == NULL) {
				//fprintf(stdout,"Azimuth read error\n");
				// exit(0);
			}
			pointer = strchr(string, ';');

			if (pointer != NULL)
				*pointer = 0;

			sscanf(string, "%f %f", &az, &amplitude);

		} while (feof(fd) == 0);

		fclose(fd);

		/* Handle 0=360 degree ambiguity */

		if ((read_count[0] == 0) && (read_count[360] != 0)) {
			read_count[0] = read_count[360];
			azimuth[0] = azimuth[360];
		}

		if ((read_count[0] != 0) && (read_count[360] == 0)) {
			read_count[360] = read_count[0];
			azimuth[360] = azimuth[0];
		}

		/* Average pattern values in case more than
		   one was read for each degree of azimuth. */

		for (x = 0; x <= 360; x++) {
			if (read_count[x] > 1)
				azimuth[x] /= (float)read_count[x];
		}

		/* Interpolate missing azimuths
		   to completely fill the array */

		last_index = -1;
		next_index = -1;

		for (x = 0; x <= 360; x++) {
			if (read_count[x] != 0) {
				if (last_index == -1)
					last_index = x;
				else
					next_index = x;
			}

			if (last_index != -1 && next_index != -1) {
				valid1 = azimuth[last_index];
				valid2 = azimuth[next_index];

				span = next_index - last_index;
				delta = (valid2 - valid1) / (float)span;

				for (y = last_index + 1; y < next_index; y++)
					azimuth[y] = azimuth[y - 1] + delta;

				last_index = y;
				next_index = -1;
			}
		}

		/* Perform azimuth pattern rotation
		   and load azimuth_pattern[361] with
		   azimuth pattern data in its final form. */

		for (x = 0; x < 360; x++) {
			y = x + (int)rintf(rotation);

			if (y >= 360)
				y -= 360;

			azimuth_pattern[y] = azimuth[x];
		}

		azimuth_pattern[360] = azimuth_pattern[0];

		got_azimuth_pattern = 255;
	}

	/* Read and process .el file */

	fd = fopen(elfile, "r");

	if (fd != NULL) {
		for (x = 0; x <= 10000; x++) {
			el_pattern[x] = 0.0;
			read_count[x] = 0;
		}

		/* Read mechanical tilt (degrees) and
		   tilt azimuth in degrees measured
		   clockwise from true North. */

		if (fgets(string, 254, fd) == NULL) {
			//fprintf(stdout,"Tilt read error\n");
			//exit(0);
		}
		pointer = strchr(string, ';');

		if (pointer != NULL)
			*pointer = 0;

		sscanf(string, "%f %f", &mechanical_tilt, &tilt_azimuth);

		/* Read elevation (degrees) and corresponding
		   normalized field radiation pattern amplitude
		   (0.0 to 1.0) until EOF is reached. */

		if (fgets(string, 254, fd) == NULL) {
			//fprintf(stdout,"Ant elevation read error\n");
			//exit(0);
		}
		pointer = strchr(string, ';');

		if (pointer != NULL)
			*pointer = 0;

		sscanf(string, "%f %f", &elevation, &amplitude);

		while (feof(fd) == 0) {
			/* Read in normalized radiated field values
			   for every 0.01 degrees of elevation between
			   -10.0 and +90.0 degrees */

			x = (int)rintf(100.0 * (elevation + 10.0));

			if (x >= 0 && x <= 10000) {
				el_pattern[x] += amplitude;
				read_count[x]++;
			}

			if (fgets(string, 254, fd) != NULL) {
				pointer = strchr(string, ';');
			}
			if (pointer != NULL)
				*pointer = 0;

			sscanf(string, "%f %f", &elevation, &amplitude);
		}

		fclose(fd);

		/* Average the field values in case more than
		   one was read for each 0.01 degrees of elevation. */

		for (x = 0; x <= 10000; x++) {
			if (read_count[x] > 1)
				el_pattern[x] /= (float)read_count[x];
		}

		/* Interpolate between missing elevations (if
		   any) to completely fill the array and provide
		   radiated field values for every 0.01 degrees of
		   elevation. */

		last_index = -1;
		next_index = -1;

		for (x = 0; x <= 10000; x++) {
			if (read_count[x] != 0) {
				if (last_index == -1)
					last_index = x;
				else
					next_index = x;
			}

			if (last_index != -1 && next_index != -1) {
				valid1 = el_pattern[last_index];
				valid2 = el_pattern[next_index];

				span = next_index - last_index;
				delta = (valid2 - valid1) / (float)span;

				for (y = last_index + 1; y < next_index; y++)
					el_pattern[y] =
					    el_pattern[y - 1] + delta;

				last_index = y;
				next_index = -1;
			}
		}

		/* Fill slant_angle[] array with offset angles based
		   on the antenna's mechanical beam tilt (if any)
		   and tilt direction (azimuth). */

		if (mechanical_tilt == 0.0) {
			for (x = 0; x <= 360; x++)
				slant_angle[x] = 0.0;
		}

		else {
			tilt_increment = mechanical_tilt / 90.0;

			for (x = 0; x <= 360; x++) {
				xx = (float)x;
				y = (int)rintf(tilt_azimuth + xx);

				while (y >= 360)
					y -= 360;

				while (y < 0)
					y += 360;

				if (x <= 180)
					slant_angle[y] =
					    -(tilt_increment * (90.0 - xx));

				if (x > 180)
					slant_angle[y] =
					    -(tilt_increment * (xx - 270.0));
			}
		}

		slant_angle[360] = slant_angle[0];	/* 360 degree wrap-around */

		for (w = 0; w <= 360; w++) {
			tilt = slant_angle[w];

	    /** Convert tilt angle to
                an array index offset **/

			y = (int)rintf(100.0 * tilt);

			/* Copy shifted el_pattern[10001] field
			   values into elevation_pattern[361][1001]
			   at the corresponding azimuth, downsampling
			   (averaging) along the way in chunks of 10. */

			for (x = y, z = 0; z <= 1000; x += 10, z++) {
				for (sum = 0.0, a = 0; a < 10; a++) {
					b = a + x;

					if (b >= 0 && b <= 10000)
						sum += el_pattern[b];
					if (b < 0)
						sum += el_pattern[0];
					if (b > 10000)
						sum += el_pattern[10000];
				}

				elevation_pattern[w][z] = sum / 10.0;
			}
		}

		got_elevation_pattern = 255;
	}

	for (x = 0; x <= 360; x++) {
		for (y = 0; y <= 1000; y++) {
			if (got_elevation_pattern)
				elevation = elevation_pattern[x][y];
			else
				elevation = 1.0;

			if (got_azimuth_pattern)
				az = azimuth_pattern[x];
			else
				az = 1.0;

			LR.antenna_pattern[x][y] = az * elevation;
		}
	}
}

void LoadSignalColors(struct site xmtr)
{
	int x, y, ok, val[4];
	char filename[255], string[80], *pointer = NULL, *s = NULL;
	FILE *fd = NULL;

	for (x = 0; xmtr.filename[x] != '.' && xmtr.filename[x] != 0 && x < 250;
	     x++)
		filename[x] = xmtr.filename[x];

	filename[x] = '.';
	filename[x + 1] = 's';
	filename[x + 2] = 'c';
	filename[x + 3] = 'f';
	filename[x + 4] = 0;

	/* Default values */

	region.level[0] = 128;
	region.color[0][0] = 255;
	region.color[0][1] = 0;
	region.color[0][2] = 0;

	region.level[1] = 118;
	region.color[1][0] = 255;
	region.color[1][1] = 165;
	region.color[1][2] = 0;

	region.level[2] = 108;
	region.color[2][0] = 255;
	region.color[2][1] = 206;
	region.color[2][2] = 0;

	region.level[3] = 98;
	region.color[3][0] = 255;
	region.color[3][1] = 255;
	region.color[3][2] = 0;

	region.level[4] = 88;
	region.color[4][0] = 184;
	region.color[4][1] = 255;
	region.color[4][2] = 0;

	region.level[5] = 78;
	region.color[5][0] = 0;
	region.color[5][1] = 255;
	region.color[5][2] = 0;

	region.level[6] = 68;
	region.color[6][0] = 0;
	region.color[6][1] = 208;
	region.color[6][2] = 0;

	region.level[7] = 58;
	region.color[7][0] = 0;
	region.color[7][1] = 196;
	region.color[7][2] = 196;

	region.level[8] = 48;
	region.color[8][0] = 0;
	region.color[8][1] = 148;
	region.color[8][2] = 255;

	region.level[9] = 38;
	region.color[9][0] = 80;
	region.color[9][1] = 80;
	region.color[9][2] = 255;

	region.level[10] = 28;
	region.color[10][0] = 0;
	region.color[10][1] = 38;
	region.color[10][2] = 255;

	region.level[11] = 18;
	region.color[11][0] = 142;
	region.color[11][1] = 63;
	region.color[11][2] = 255;

	region.level[12] = 8;
	region.color[12][0] = 140;
	region.color[12][1] = 0;
	region.color[12][2] = 128;

	region.levels = 13;

	fd = fopen(filename, "r");

	if (fd == NULL)
		fd = fopen(filename, "r");

	if (fd == NULL) {
		fd = fopen(filename, "w");

		for (x = 0; x < region.levels; x++)
			fprintf(fd, "%3d: %3d, %3d, %3d\n", region.level[x],
				region.color[x][0], region.color[x][1],
				region.color[x][2]);

		fclose(fd);
	}

	else {
		x = 0;
		s = fgets(string, 80, fd);

		while (x < 128 && feof(fd) == 0) {
			pointer = strchr(string, ';');

			if (pointer != NULL)
				*pointer = 0;

			ok = sscanf(string, "%d: %d, %d, %d", &val[0], &val[1],
				    &val[2], &val[3]);

			if (ok == 4) {
				for (y = 0; y < 4; y++) {
					if (val[y] > 255)
						val[y] = 255;

					if (val[y] < 0)
						val[y] = 0;
				}

				region.level[x] = val[0];
				region.color[x][0] = val[1];
				region.color[x][1] = val[2];
				region.color[x][2] = val[3];
				x++;
			}

			s = fgets(string, 80, fd);
		}

		fclose(fd);
		region.levels = x;
	}
}

void LoadLossColors(struct site xmtr)
{
	int x, y, ok, val[4];
	char filename[255], string[80], *pointer = NULL, *s = NULL;
	FILE *fd = NULL;

	for (x = 0; xmtr.filename[x] != '.' && xmtr.filename[x] != 0 && x < 250;
	     x++)
		filename[x] = xmtr.filename[x];

	filename[x] = '.';
	filename[x + 1] = 'l';
	filename[x + 2] = 'c';
	filename[x + 3] = 'f';
	filename[x + 4] = 0;

	/* Default values */

	region.level[0] = 80;
	region.color[0][0] = 255;
	region.color[0][1] = 0;
	region.color[0][2] = 0;

	region.level[1] = 90;
	region.color[1][0] = 255;
	region.color[1][1] = 128;
	region.color[1][2] = 0;

	region.level[2] = 100;
	region.color[2][0] = 255;
	region.color[2][1] = 165;
	region.color[2][2] = 0;

	region.level[3] = 110;
	region.color[3][0] = 255;
	region.color[3][1] = 206;
	region.color[3][2] = 0;

	region.level[4] = 120;
	region.color[4][0] = 255;
	region.color[4][1] = 255;
	region.color[4][2] = 0;

	region.level[5] = 130;
	region.color[5][0] = 184;
	region.color[5][1] = 255;
	region.color[5][2] = 0;

	region.level[6] = 140;
	region.color[6][0] = 0;
	region.color[6][1] = 255;
	region.color[6][2] = 0;

	region.level[7] = 150;
	region.color[7][0] = 0;
	region.color[7][1] = 208;
	region.color[7][2] = 0;

	region.level[8] = 160;
	region.color[8][0] = 0;
	region.color[8][1] = 196;
	region.color[8][2] = 196;

	region.level[9] = 170;
	region.color[9][0] = 0;
	region.color[9][1] = 148;
	region.color[9][2] = 255;

	region.level[10] = 180;
	region.color[10][0] = 80;
	region.color[10][1] = 80;
	region.color[10][2] = 255;

	region.level[11] = 190;
	region.color[11][0] = 0;
	region.color[11][1] = 38;
	region.color[11][2] = 255;

	region.level[12] = 200;
	region.color[12][0] = 142;
	region.color[12][1] = 63;
	region.color[12][2] = 255;

	region.level[13] = 210;
	region.color[13][0] = 196;
	region.color[13][1] = 54;
	region.color[13][2] = 255;

	region.level[14] = 220;
	region.color[14][0] = 255;
	region.color[14][1] = 0;
	region.color[14][2] = 255;

	region.level[15] = 230;
	region.color[15][0] = 255;
	region.color[15][1] = 194;
	region.color[15][2] = 204;

	region.levels = 16;

	fd = fopen(filename, "r");

	if (fd == NULL)
		fd = fopen(filename, "r");

	if (fd == NULL) {
		fd = fopen(filename, "w");

		for (x = 0; x < region.levels; x++)
			fprintf(fd, "%3d: %3d, %3d, %3d\n", region.level[x],
				region.color[x][0], region.color[x][1],
				region.color[x][2]);

		fclose(fd);
	}

	else {
		x = 0;
		s = fgets(string, 80, fd);

		while (x < 128 && feof(fd) == 0) {
			pointer = strchr(string, ';');

			if (pointer != NULL)
				*pointer = 0;

			ok = sscanf(string, "%d: %d, %d, %d", &val[0], &val[1],
				    &val[2], &val[3]);

			if (ok == 4) {
				for (y = 0; y < 4; y++) {
					if (val[y] > 255)
						val[y] = 255;

					if (val[y] < 0)
						val[y] = 0;
				}

				region.level[x] = val[0];
				region.color[x][0] = val[1];
				region.color[x][1] = val[2];
				region.color[x][2] = val[3];
				x++;
			}

			s = fgets(string, 80, fd);
		}

		fclose(fd);
		region.levels = x;
	}
}

void LoadDBMColors(struct site xmtr)
{
	int x, y, ok, val[4];
	char filename[255], string[80], *pointer = NULL, *s = NULL;
	FILE *fd = NULL;

	for (x = 0; xmtr.filename[x] != '.' && xmtr.filename[x] != 0 && x < 250;
	     x++)
		filename[x] = xmtr.filename[x];

	filename[x] = '.';
	filename[x + 1] = 'd';
	filename[x + 2] = 'c';
	filename[x + 3] = 'f';
	filename[x + 4] = 0;

	/* Default values */

	region.level[0] = 0;
	region.color[0][0] = 255;
	region.color[0][1] = 0;
	region.color[0][2] = 0;

	region.level[1] = -10;
	region.color[1][0] = 255;
	region.color[1][1] = 128;
	region.color[1][2] = 0;

	region.level[2] = -20;
	region.color[2][0] = 255;
	region.color[2][1] = 165;
	region.color[2][2] = 0;

	region.level[3] = -30;
	region.color[3][0] = 255;
	region.color[3][1] = 206;
	region.color[3][2] = 0;

	region.level[4] = -40;
	region.color[4][0] = 255;
	region.color[4][1] = 255;
	region.color[4][2] = 0;

	region.level[5] = -50;
	region.color[5][0] = 184;
	region.color[5][1] = 255;
	region.color[5][2] = 0;

	region.level[6] = -60;
	region.color[6][0] = 0;
	region.color[6][1] = 255;
	region.color[6][2] = 0;

	region.level[7] = -70;
	region.color[7][0] = 0;
	region.color[7][1] = 208;
	region.color[7][2] = 0;

	region.level[8] = -80;
	region.color[8][0] = 0;
	region.color[8][1] = 196;
	region.color[8][2] = 196;

	region.level[9] = -90;
	region.color[9][0] = 0;
	region.color[9][1] = 148;
	region.color[9][2] = 255;

	region.level[10] = -100;
	region.color[10][0] = 80;
	region.color[10][1] = 80;
	region.color[10][2] = 255;

	region.level[11] = -110;
	region.color[11][0] = 0;
	region.color[11][1] = 38;
	region.color[11][2] = 255;

	region.level[12] = -120;
	region.color[12][0] = 142;
	region.color[12][1] = 63;
	region.color[12][2] = 255;

	region.level[13] = -130;
	region.color[13][0] = 196;
	region.color[13][1] = 54;
	region.color[13][2] = 255;

	region.level[14] = -140;
	region.color[14][0] = 255;
	region.color[14][1] = 0;
	region.color[14][2] = 255;

	region.level[15] = -150;
	region.color[15][0] = 255;
	region.color[15][1] = 194;
	region.color[15][2] = 204;

	region.levels = 16;

	fd = fopen(filename, "r");

	if (fd == NULL)
		fd = fopen(filename, "r");

	if (fd == NULL) {
		fd = fopen(filename, "w");

		for (x = 0; x < region.levels; x++)
			fprintf(fd, "%+4d: %3d, %3d, %3d\n", region.level[x],
				region.color[x][0], region.color[x][1],
				region.color[x][2]);

		fclose(fd);
	}

	else {
		x = 0;
		s = fgets(string, 80, fd);

		while (x < 128 && feof(fd) == 0) {
			pointer = strchr(string, ';');

			if (pointer != NULL)
				*pointer = 0;

			ok = sscanf(string, "%d: %d, %d, %d", &val[0], &val[1],
				    &val[2], &val[3]);

			if (ok == 4) {
				if (val[0] < -200)
					val[0] = -200;

				if (val[0] > +40)
					val[0] = +40;

				region.level[x] = val[0];

				for (y = 1; y < 4; y++) {
					if (val[y] > 255)
						val[y] = 255;

					if (val[y] < 0)
						val[y] = 0;
				}

				region.color[x][0] = val[1];
				region.color[x][1] = val[2];
				region.color[x][2] = val[3];
				x++;
			}

			s = fgets(string, 80, fd);
		}

		fclose(fd);
		region.levels = x;
	}
}

void LoadTopoData(int max_lon, int min_lon, int max_lat, int min_lat,
		  int winfiles)
{
	/* This function loads the SDF files required
	   to cover the limits of the region specified. */

	int x, y, width, ymin, ymax;

	width = ReduceAngle(max_lon - min_lon);

	if ((max_lon - min_lon) <= 180.0) {
		for (y = 0; y <= width; y++)
			for (x = min_lat; x <= max_lat; x++) {
				ymin = (int)(min_lon + (double)y);

				while (ymin < 0)
					ymin += 360;

				while (ymin >= 360)
					ymin -= 360;

				ymax = ymin + 1;

				while (ymax < 0)
					ymax += 360;

				while (ymax >= 360)
					ymax -= 360;

				if (winfiles == 1) {
					if (ippd == 3600)
						snprintf(string, 19,
							 "%d=%d=%d=%d=hd", x,
							 x + 1, ymin, ymax);
					else
						snprintf(string, 16,
							 "%d=%d=%d=%d", x,
							 x + 1, ymin, ymax);

				} else {
					if (ippd == 3600)
						snprintf(string, 19,
							 "%d:%d:%d:%d-hd", x,
							 x + 1, ymin, ymax);
					else
						snprintf(string, 16,
							 "%d:%d:%d:%d", x,
							 x + 1, ymin, ymax);
				}

				LoadSDF(string, winfiles);
			}
	}

	else {
		for (y = 0; y <= width; y++)
			for (x = min_lat; x <= max_lat; x++) {
				ymin = max_lon + y;

				while (ymin < 0)
					ymin += 360;

				while (ymin >= 360)
					ymin -= 360;

				ymax = ymin + 1;

				while (ymax < 0)
					ymax += 360;

				while (ymax >= 360)
					ymax -= 360;

				if (winfiles == 1) {
					if (ippd == 3600)
						snprintf(string, 19,
							 "%d=%d=%d=%d=hd", x,
							 x + 1, ymin, ymax);
					else
						snprintf(string, 16,
							 "%d=%d=%d=%d", x,
							 x + 1, ymin, ymax);

				} else {
					if (ippd == 3600)
						snprintf(string, 19,
							 "%d:%d:%d:%d-hd", x,
							 x + 1, ymin, ymax);
					else
						snprintf(string, 16,
							 "%d:%d:%d:%d", x,
							 x + 1, ymin, ymax);
				}

				LoadSDF(string, winfiles);
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

	int i, x, y, z, ypix, xpix, tempxpix, tempypix, fd = 0, n = 0;
	char input[80], str[3][80], tempname[15], *pointer = NULL, *s = NULL;
	double latitude, longitude, height, tempheight;
	FILE *fd1 = NULL, *fd2 = NULL;

	strcpy(tempname, "/tmp/XXXXXX\0");

	fd1 = fopen(filename, "r");

	if (fd1 != NULL) {
		fd = mkstemp(tempname);
		fd2 = fopen(tempname, "w");

		s = fgets(input, 78, fd1);

		pointer = strchr(input, ';');

		if (pointer != NULL)
			*pointer = 0;

		while (feof(fd1) == 0) {
			/* Parse line for latitude, longitude, height */

			for (x = 0, y = 0, z = 0;
			     x < 78 && input[x] != 0 && z < 3; x++) {
				if (input[x] != ',' && y < 78) {
					str[z][y] = input[x];
					y++;
				}

				else {
					str[z][y] = 0;
					z++;
					y = 0;
				}
			}

			latitude = ReadBearing(str[0]);
			longitude = ReadBearing(str[1]);

			if (longitude < 0.0)
				longitude += 360;

			/* Remove <CR> and/or <LF> from antenna height string */

			for (i = 0;
			     str[2][i] != 13 && str[2][i] != 10
			     && str[2][i] != 0; i++) ;

			str[2][i] = 0;

			/* The terrain feature may be expressed in either
			   feet or meters.  If the letter 'M' or 'm' is
			   discovered in the string, then this is an
			   indication that the value given is expressed
			   in meters.  Otherwise the height is interpreted
			   as being expressed in feet.  */

			for (i = 0;
			     str[2][i] != 'M' && str[2][i] != 'm'
			     && str[2][i] != 0 && i < 48; i++) ;

			if (str[2][i] == 'M' || str[2][i] == 'm') {
				str[2][i] = 0;
				height = rint(atof(str[2]));
			}

			else {
				str[2][i] = 0;
				height = rint(METERS_PER_FOOT * atof(str[2]));
			}

			if (height > 0.0)
				fprintf(fd2, "%d, %d, %f\n",
					(int)rint(latitude / dpp),
					(int)rint(longitude / dpp), height);

			s = fgets(input, 78, fd1);

			pointer = strchr(input, ';');

			if (pointer != NULL)
				*pointer = 0;
		}

		fclose(fd1);
		fclose(fd2);
		close(fd);

		fd1 = fopen(tempname, "r");
		fd2 = fopen(tempname, "r");

		y = 0;

		n = fscanf(fd1, "%d, %d, %lf", &xpix, &ypix, &height);

		do {
			x = 0;
			z = 0;

			n = fscanf(fd2, "%d, %d, %lf", &tempxpix, &tempypix,
				   &tempheight);

			do {
				if (x > y && xpix == tempxpix
				    && ypix == tempypix) {
					z = 1;	/* Dupe! */

					if (tempheight > height)
						height = tempheight;
				}

				else {
					n = fscanf(fd2, "%d, %d, %lf",
						   &tempxpix, &tempypix,
						   &tempheight);
					x++;
				}

			} while (feof(fd2) == 0 && z == 0);

			if (z == 0)
				/* No duplicate found */
				//fprintf(stdout,"%lf, %lf \n",xpix*dpp, ypix*dpp);
				fflush(stdout);
			    AddElevation(xpix * dpp, ypix * dpp, height);
			fflush(stdout);

			n = fscanf(fd1, "%d, %d, %lf", &xpix, &ypix, &height);
			y++;

			rewind(fd2);

		} while (feof(fd1) == 0);

		fclose(fd1);
		fclose(fd2);
		unlink(tempname);
	}

}
