#ifndef _IMAGE_HH_
#define _IMAGE_HH_

#include <stdint.h>

#define RGB_SIZE  3
#define RGBA_SIZE 4

typedef enum _IMAGE_FORMAT{	IMAGE_DEFAULT = 0, \
							IMAGE_PPM, \
							IMAGE_LIBRARY, \
							IMAGE_FORMAT_MAX \
						} IMAGE_FORMAT;

typedef enum _IMAGE_MODEL{	IMAGE_RGB, \
							IMAGE_RGBA, \
							IMAGE_MODEL_MAX
						} IMAGE_MODEL;

typedef struct _IMAGE_CTX{
	size_t width;
	size_t height;
	IMAGE_MODEL model;
	IMAGE_FORMAT format;
	uint8_t *canvas;
	uint8_t *next_pixel;
	uint32_t initialized;
	char *extension;
	void *_dt;
} IMAGE_CTX, *PIMAGE_CTX;

int image_set_format(IMAGE_FORMAT);
int image_init(PIMAGE_CTX, const size_t, const size_t, const IMAGE_MODEL, const IMAGE_FORMAT);
int image_add_pixel(PIMAGE_CTX ctx, const uint8_t, const uint8_t, const uint8_t, const uint8_t);
int image_set_pixel(PIMAGE_CTX ctx, const size_t, const size_t, const uint8_t, const uint8_t, const uint8_t, const uint8_t);
int image_get_pixel(PIMAGE_CTX ctx,const size_t,const size_t, uint8_t const*, uint8_t const*, uint8_t const*, uint8_t const*);
int image_get_filename(PIMAGE_CTX, char*, size_t, char*);
int image_write(PIMAGE_CTX, FILE*);
void image_free(PIMAGE_CTX);
int image_set_library(char*);

#define ADD_PIXEL(ctx,r,g,b) image_add_pixel((ctx),(r),(g),(b),0xff)
#define ADD_PIXELA(ctx,r,g,b,a) image_add_pixel((ctx),(r),(g),(b),(a))

#endif