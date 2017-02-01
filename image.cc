#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdio.h>
#include <errno.h>
#include <dlfcn.h>
#include "image.hh"
#include "image-ppm.hh"

typedef int _init(PIMAGE_CTX);
typedef int _add_pixel(PIMAGE_CTX,const uint8_t,const uint8_t,const uint8_t,const uint8_t);
typedef int _get_pixel(PIMAGE_CTX,const size_t,const size_t,const uint8_t*,const uint8_t*,const uint8_t*,const uint8_t*);
typedef int _write(PIMAGE_CTX,FILE*);

struct image_dispatch_table{
	_init *init;
	_add_pixel *add_pixel;
	_get_pixel *get_pixel;
	_write *write;
};

static IMAGE_FORMAT format = IMAGE_PPM;
char *dynamic_backend = NULL;

int load_library(struct image_dispatch_table *dt){
	void *hndl;
	int success = 0;

	if(dynamic_backend == NULL){
		fprintf(stderr,"Custom image processor requested without specification\n");
		return EINVAL;
	}
	
	hndl = dlopen(dynamic_backend,RTLD_LAZY);
	if(hndl == NULL){
		fprintf(stderr,"Error loading shared object\n");
		return EINVAL;
	}
	
	dt->init = (_init*)dlsym(hndl,"lib_init");
	dt->add_pixel = (_add_pixel*)dlsym(hndl,"lib_add_pixel");
	dt->get_pixel = (_get_pixel*)dlsym(hndl,"lib_get_pixel");
	dt->write = (_write*)dlsym(hndl,"lib_write");

	if(dt->init == NULL ||
		dt->add_pixel == NULL ||
		dt->get_pixel == NULL ||
		dt->write == NULL){
		fprintf(stderr,"Invalid image processing module specified\n");
		success = dlclose(hndl);
	}

	return success;
}

int get_dt(struct image_dispatch_table *dt, IMAGE_FORMAT format){
	int success = 0;

	memset((void*)dt,0x00,sizeof(struct image_dispatch_table));
	switch(format){
		case IMAGE_PPM:
			dt->init = ppm_init;
			dt->add_pixel = ppm_add_pixel;
			dt->get_pixel = ppm_get_pixel;
			dt->write = ppm_write;
			break;
		case IMAGE_LIBRARY:
			success = load_library(dt);
			break;
		default:
			success = EINVAL;
	}
	return success;
}

int image_init(PIMAGE_CTX ctx, \
				const size_t width, \
				const size_t height, \
				const IMAGE_MODEL model) {

	int success = 0;
	struct image_dispatch_table *dt;

	/* Perform some sanity checking on provided arguments */
	if(ctx == NULL)
		return EINVAL;
	if(width == 0 || height == 0)
		return EINVAL;
	if(model < 0 || model > IMAGE_MODEL_MAX)
		return EINVAL;
	if(format < 0 || format > IMAGE_FORMAT_MAX)
		return EINVAL;

	memset(ctx,0x00,sizeof(IMAGE_CTX));

	/* Assign the initialize values to the processing context */
	ctx->width = width;
	ctx->height = height;
	ctx->model = model;
	ctx->format = format;

	/* Get the dispatch table for this image format */
	dt = (struct image_dispatch_table*) calloc(1,sizeof(struct image_dispatch_table));
	if(dt == NULL){
		fprintf(stderr,"Error allocating dispatch table\n");
		return ENOMEM;
	}
	success = get_dt(dt,format);
	if(success != 0){
		fprintf(stderr,"Error locating dispatch table\n");
		free(dt);
		return success;
	}
	ctx->_dt = (void*)dt;

	/* Call the format-specific initialization function */
	success = dt->init(ctx);
	if(success != 0){
		fprintf(stderr,"Error initializing image context\n");
		free(dt);
		return success;
	}

	ctx->initialized = 1;

	return success;
}
int image_add_pixel(PIMAGE_CTX ctx, const uint8_t r, const uint8_t g, const uint8_t b, const uint8_t a){
	struct image_dispatch_table *dt = (struct image_dispatch_table*)ctx->_dt;
	return dt->add_pixel(ctx,r,g,b,a);
}
int image_set_pixel(PIMAGE_CTX ctx, const size_t x, const size_t y, const uint8_t r, const uint8_t g, const uint8_t b, const uint8_t a){
	size_t block_size;
	block_size = ctx->model == IMAGE_RGB ? RGB_SIZE : RGBA_SIZE;
	ctx->next_pixel = ctx->canvas + (x * block_size) + (ctx->width * block_size * y);
	struct image_dispatch_table *dt = (struct image_dispatch_table*)ctx->_dt;
	return dt->add_pixel(ctx,r,g,b,a);
}
int image_get_pixel(PIMAGE_CTX ctx, const size_t x, const size_t y, const uint8_t *r, const uint8_t *g, const uint8_t *b, const uint8_t *a){
	struct image_dispatch_table *dt = (struct image_dispatch_table*)ctx->_dt;
	return dt->get_pixel(ctx,x,y,r,g,b,a);
}
int image_get_filename(PIMAGE_CTX ctx, char *out, size_t len_out, char *in){
	size_t len_src;
	size_t len_ext;
	int success = 0;

	if(ctx->initialized != 1){
		fprintf(stderr,"Called get filename before initialization\n");
		return EINVAL;
	}

	len_src = strlen(in);
	len_ext = strlen(ctx->extension);

	if(len_src == 0){
		in = "output";
		len_src = 6;
	}

	if(len_src > len_ext && strcmp(in+len_src-len_ext,ctx->extension) == 0){
		/* Already has correct extension and fits in buffer */
		if(len_src < len_out)
			strncpy(in,out,len_out);
		else
			success = ENOMEM;
	}else if(len_src > len_ext){
		/* Doesn't have correct extension and fits */
		if(len_src + len_ext < len_out){
			strncpy(out,in,len_out);
			strncat(out,ctx->extension,len_out);
		}else
			success = ENOMEM;
	}else{
		/* The input buffer plus an extension cannot fit in the output buffer */
		fprintf(stderr,"Error building image output filename\n");
		success = ENOMEM;
	}
	return success;
}
int image_write(PIMAGE_CTX ctx, FILE *fd){
	struct image_dispatch_table *dt = (struct image_dispatch_table*)ctx->_dt;
	return dt->write(ctx,fd);
}
int image_set_library(char *library){
	char *libname;
	size_t length;

	length = strlen(library) + 1;
	libname = (char*)calloc(length,sizeof(char));
	if(libname == NULL)
		return ENOMEM;
	strncpy(libname,library,length);

	dynamic_backend = libname;
	format = IMAGE_LIBRARY;
	return 0;
}

