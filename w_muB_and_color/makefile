SHELL=/bin/sh

SRCS= \
HydrodynamicalNoiseChecks.cpp \
gauss_quadrature.cpp \
lib.cpp

HDRS= \
gauss_quadrature.h \
defs.h \
lib.h

MAKEFILE=makefile

COMMAND=HNC

OBJS= $(addsuffix .o, $(basename $(SRCS)))

CC= g++
CFLAGS=  -pg
WARNFLAGS= -ansi -pedantic -Wall -W \
   -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align \
   -Wwrite-strings -fshort-enums -fno-common 
LDFLAGS= -lgsl -lgslcblas 
LIBS= -L/sw/lib -I/sw/include

$(COMMAND): $(OBJS) $(HDRS) $(MAKEFILE) 
	$(CC) $(CFLAGS) -o $(COMMAND) $(OBJS) $(LDFLAGS) $(LIBS)

clean:
	rm -f $(OBJS)
 
tarz:
	tar zcf - $(MAKEFILE) $(SRCS) $(HDRS) > $(COMMAND).tar.gz
 
