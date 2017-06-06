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
	union{
		double	xll;
		double	max_west;
	};
	union{
		double	yll;
		double	min_north;
	};
	union{
		double	xur;
		double	min_west;
	};
	union{
		double	yur;
		double	max_north;
	};
	double	cellsize;
	long long datastart;
	int		nodata;
	int 	max_el;
	int		min_el;
	int		*data;
} tile_t, *ptile_t;

int tile_load_lidar(tile_t*, char *);

#endif