#ifndef _IMAGE_PPM_HH
#define _IMAGE_PPM_HH

#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>

int ppm_init(PIMAGE_CTX ctx);
int ppm_add_pixel(PIMAGE_CTX ctx,const uint8_t r,const uint8_t g,const uint8_t b,const uint8_t a);
int ppm_get_pixel(PIMAGE_CTX ctx,const size_t x,const size_t y,const uint8_t *r,const uint8_t *g,const uint8_t *b,const uint8_t *a);
int ppm_write(PIMAGE_CTX ctx, FILE* fd);

#endif