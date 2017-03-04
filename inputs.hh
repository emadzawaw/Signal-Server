#ifndef _INPUTS_HH_
#define _INPUTS_HH_

#include "common.h"

int LoadSDF_SDF(char *name, int winfiles);
int LoadSDF(char *name, int winfiles);
int LoadPAT(char *az_filename, char *el_filename);
int LoadSignalColors(struct site xmtr);
int LoadLossColors(struct site xmtr);
int LoadDBMColors(struct site xmtr);
int LoadTopoData(int max_lon, int min_lon, int max_lat, int min_lat);
void LoadUDT(char *filename);
int loadLIDAR(char *filename);
int loadClutter(char *filename, double radius, struct site tx);

static const char AZ_FILE_SUFFIX[] = ".az";
static const char EL_FILE_SUFFIX[] = ".el"; 

#endif /* _INPUTS_HH_ */
