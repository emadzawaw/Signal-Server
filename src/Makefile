<<<<<<< HEAD
# The PostgreSQL make files exploit features of GNU make that other
# makes do not have. Because it is a common mistake for users to try
# to build Postgres with a different make, we have this make file
# that, as a service, will look for a GNU make and invoke it, or show
# an error message if none could be found.

# If the user were using GNU make now, this file would not get used
# because GNU make uses a make file named "GNUmakefile" in preference
# to "Makefile" if it exists. PostgreSQL is shipped with a
# "GNUmakefile". If the user hasn't run the configure script yet, the
# GNUmakefile won't exist yet, so we catch that case as well.


all check install clean signalserver utils:
	@if [ ! -f GNUmakefile ] ; then \
	   echo "ACK!!!!! Your repo is not complete, you need to fix that"; \
	   false ; \
	 fi
	@IFS=':' ; \
	 for dir in $$PATH; do \
	   for prog in gmake gnumake make; do \
	     if [ -f $$dir/$$prog ] && ( $$dir/$$prog -f /dev/null --version 2>/dev/null | grep GNU >/dev/null 2>&1 ) ; then \
	       GMAKE=$$dir/$$prog; \
	       break 2; \
	     fi; \
	   done; \
	 done; \
	\
	 if [ x"$${GMAKE+set}" = xset ]; then \
	   echo "Using GNU make found at $${GMAKE}"; \
	   unset MAKEFLAGS; unset MAKELEVEL; \
	   $${GMAKE} $@ ; \
	 else \
	   echo "You must use GNU make to build PostgreSQL." ; \
	   false; \
	 fi
=======
SHELL		= /bin/sh

<<<<<<< HEAD:src/GNUmakefile
CC		= cc
CXX		= c++
CFLAGS		= -Wall -O3 -ffast-math -DCROPPING
CXXFLAGS	= -Wall -O3 -ffast-math
LIBS		= -lm -lpthread -ldl -lz -lbz2
=======
CC		= gcc
CXX		= g++
CFLAGS		= -Wall -O3 -s -ffast-math -DCROPPING
CXXFLAGS	= -Wall -O3 -s -ffast-math
LIBS		= -lm -lpthread -ldl -lbz2 -lz
>>>>>>> 77bb415059863568992ffa8fc82143df8c7b3ebf:src/Makefile

INSTALL_PATH	= ..
VPATH		= models
objects 	= main.o cost.o ecc33.o ericsson.o fspl.o hata.o itwom3.0.o \
		  los.o sui.o pel.o egli.o soil.o inputs.o outputs.o image.o \
		  image-ppm.o tiles.o

GCC_MAJOR	:= $(shell $(CXX) -dumpversion 2>&1 | cut -d . -f 1)
GCC_MINOR	:= $(shell $(CXX) -dumpversion 2>&1 | cut -d . -f 2)
GCC_VER_OK	:= $(shell test $(GCC_MAJOR) -ge 4 && \
			   test $(GCC_MINOR) -ge 7 && \
			   echo 1)

#ifneq "$(GCC_VER_OK)" "1"
#error:
#	@echo "Requires GCC version >= 4.7"
#	@exit
#endif

%.o : %.cc
	@echo -e "    CXX\t$@"
	@$ $(CXX) $(CXXFLAGS) -c $<

%.o : %.c
	@echo -e "    CC\t$@"
	@$ $(CC) $(CFLAGS) -c $<

signalserver: $(objects)
	@echo -e "    LNK\t$@"
	@$(CXX) $(objects) -o $@ ${LIBS}
	@echo -e " SYMLNK\tsignalserverHD -> $@"
	@ln -sf $@ signalserverHD
	@echo -e " SYMLNK\tsignalserverLIDAR -> $@"
	@ln -sf $@ signalserverLIDAR

main.o: main.cc common.h inputs.hh outputs.hh itwom3.0.hh los.hh pel.hh \
	image.hh

inputs.o: inputs.cc common.h main.hh tiles.hh

outputs.o: outputs.cc common.h main.hh inputs.hh cost.hh ecc33.hh ericsson.hh \
	   fspl.hh hata.hh itwom3.0.hh sui.hh image.hh

image.o: image.cc image-ppm.o image.hh image-ppm.hh

image-ppm.o: image-ppm.cc image.hh

los.o: los.cc los.hh main.hh cost.hh ecc33.hh ericsson.hh fspl.hh hata.hh \
       itwom3.0.hh sui.hh pel.hh egli.hh soil.hh

tiles.o: tiles.cc tiles.hh common.h


testmodels.o: testmodels.cc cost.hh ecc33.hh ericsson.hh fspl.hh hata.hh \
	egli.hh sui.hh itwom3.0.hh

testmodels: testmodels.o cost.o ecc33.o ericsson.o fspl.o hata.o egli.o sui.o itwom3.0.o
	@echo -e "    LNK\t$@"
	@$(CXX) testmodels.o cost.o ecc33.o ericsson.o fspl.o hata.o egli.o \
	sui.o itwom3.0.o -o $@ ${LIBS}

.PHONY: clean
clean:
<<<<<<< HEAD:src/GNUmakefile
	rm -f $(objects) signalserver signalserverHD signalserverLIDAR

all:	signalserver
=======
	rm -f $(objects) signalserver signalserverHD signalserverLIDAR testmodels

all: signalserver

install: all
	cp -a signalserver signalserverHD signalserverLIDAR $(INSTALL_PATH)

>>>>>>> 77bb415059863568992ffa8fc82143df8c7b3ebf:src/Makefile
>>>>>>> 77bb415059863568992ffa8fc82143df8c7b3ebf
