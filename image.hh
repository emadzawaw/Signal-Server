#ifndef _IMAGE_HH_
#define _IMAGE_HH_

#include <stdint.h>

#define RGB_SIZE  3
#define RGBA_SIZE 4

typedef enum _IMAGE_FORMAT{	IMAGE_PPM, \
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

int image_init(PIMAGE_CTX ctx, const size_t width, const size_t height, const IMAGE_MODEL model);
int image_add_pixel(PIMAGE_CTX ctx, const uint8_t r, const uint8_t g, const uint8_t b, const uint8_t a);
int image_set_pixel(PIMAGE_CTX ctx, const size_t x, const size_t y, const uint8_t r, const uint8_t g, const uint8_t b, const uint8_t a);
int image_get_pixel(PIMAGE_CTX ctx,const size_t x,const size_t y, uint8_t const *r, uint8_t const *g, uint8_t const *b, uint8_t const *a);
int image_get_filename(PIMAGE_CTX ctx, char* out, size_t len, char* in);
int image_write(PIMAGE_CTX ctx, FILE *fd);
int image_set_library(char *library);

#define ADD_PIXEL(ctx,r,g,b) image_add_pixel((ctx),(r),(g),(b),0)
#define ADD_PIXELA(ctx,r,g,b,a) image_add_pixel((ctx),(r),(g),(b),(a))

#endif