#ifndef _TILES_HH_
#define _TILES_HH_

typedef struct _tile_t{
	char	*filename;
	union{
		int	cols;
		int	width;
	};
	union{
		int	rows;
		int	height;
	};
	double	xll;
	double	yll;
	double	xur;
	double	yur;
	double	cellsize;
	long long datastart;
	int		nodata;
	int 	max_el;
	int		min_el;
	int		*data;
} tile_t, *ptile_t;

int tile_load_lidar(tile_t*, char *, int);

#endif