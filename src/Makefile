INSTALL_PREFIX = /usr/local
DCAMPROF_VERSION = $(shell cat VERSION)
CC = gcc
CFLAGS = -Wall -std=c99 -g -O2 -D_XOPEN_SOURCE=600 -D__EXTENSIONS__ -fopenmp
LDFLAGS = -fopenmp
INCLUDE = -I.
DEFINE =

PROGS = dcamprof

DCAMPROF_OBJS = dcamprof.o \
nmsimplex.o \
profio.o \
argyllio.o \
tps.o \
colmath.o \
matopt.o \
lut.o \
lut3d.o \
look.o \
dngref.o \
gamut.o \
spectrum.o \
jsonio.o \
cJSON.o \
interp.o \
bisection.o \
icclut.o \
dnglut.o \
glare.o \
target.o \
xyz2spec.o \
observers.o \
spectraldb.o \
tifio.o \
elog.o \
wcompat.o \
strbuf.o \
spectraldb_munsell.o \
spectraldb_cc24.o \
spectraldb_kuopio_natural.o

all: dcamprof

%.o: %.c
	$(CC) -o $@ -c $(INCLUDE) $(CFLAGS) $(DEFINE) -DDCAMPROF_VERSION=\"$(DCAMPROF_VERSION)\" $<

dcamprof: $(DCAMPROF_OBJS)
	$(CC) -o $@ $(LDFLAGS) $(DCAMPROF_OBJS) -llcms2 -ltiff -lm

install: $(PROGS)
	install -d $(INSTALL_PREFIX)/bin
	install $(PROGS) $(INSTALL_PREFIX)/bin

clean:
	rm -f $(DCAMPROF_OBJS) ${PROGS}
