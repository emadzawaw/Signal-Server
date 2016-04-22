SHELL		= /bin/sh

CC		= gcc
CXX		= g++
CFLAGS		= -Wall -O3 -s -ffast-math
CXXFLAGS	= -Wall -O3 -s -ffast-math
LIBS		= -lm -lpthread

VPATH		= models
objects 	= main.o cost.o ecc33.o ericsson.o fspl.o hata.o itwom3.0.o \
		  los.o sui.o inputs.o outputs.o

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
	
main.o: main.cc common.h inputs.hh outputs.hh itwom3.0.hh los.hh

inputs.o: inputs.cc common.h main.hh

outputs.o: outputs.cc common.h inputs.hh main.hh cost.hh ecc33.hh ericsson.hh \
	   fspl.hh hata.hh itwom3.0.hh sui.hh

los.o: los.cc common.h main.hh cost.hh ecc33.hh ericsson.hh fspl.hh hata.hh \
       itwom3.0.hh sui.hh

.PHONY: clean
clean:
	rm -f $(objects) signalserver signalserverHD signalserverLIDAR
