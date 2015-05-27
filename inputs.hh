#ifndef _INPUTS_HH_
#define _INPUTS_HH_

#include "common.h"

int LoadSDF_SDF(char *name, int winfiles);
char LoadSDF(char *name, int winfiles);
void LoadPAT(char *filename);
void LoadSignalColors(struct site xmtr);
void LoadLossColors(struct site xmtr);
void LoadDBMColors(struct site xmtr);
void LoadTopoData(int max_lon, int min_lon, int max_lat, int min_lat,
		  int winfiles);
void LoadUDT(char *filename);

#endif /* _INPUTS_HH_ */
