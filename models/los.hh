#ifndef _LOS_HH_
#define _LOS_HH_

#include <stdio.h>

#include "../common.h"

void PlotLOSPath(struct site source, struct site destination, char mask_value,
		 FILE *fd);
void PlotPropPath(struct site source, struct site destination,
		  unsigned char mask_value, FILE * fd, int propmodel,
		  int knifeedge, int pmenv);
void PlotLOSMap(struct site source, double altitude, char *plo_filename);
void PlotPropagation(struct site source, double altitude, char *plo_filename,
		     int propmodel, int knifeedge, int haf, int pmenv);
void PlotPath(struct site source, struct site destination, char mask_value);

#endif /* _LOS_HH_ */
