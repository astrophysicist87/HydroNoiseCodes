SHELL=/bin/sh

SRCS= \
ElectricChargeFluctuations.cpp \
gauss_quadrature.cpp \
lib.cpp

HDRS= \
gauss_quadrature.h \
defs.h \
lib.h

MAKEFILE=makefile

COMMAND=ECF

OBJS= $(addsuffix .o, $(basename $(SRCS)))

CC= g++
CFLAGS=  -g
WARNFLAGS= -ansi -pedantic -Wall -W \
   -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align \
   -Wwrite-strings -fshort-enums -fno-common 
LDFLAGS= $(CFLAGS) -lgsl -lgslcblas 
LIBS= -L/sw/lib -I/sw/include

 
$(COMMAND): $(OBJS) $(HDRS) $(MAKEFILE) 
	$(CC) -o $(COMMAND) $(OBJS) $(LDFLAGS) $(LIBS)

ElectricChargeFluctuations.o : ElectricChargeFluctuations.cpp defs.h
	$(CC) $(CFLAGS) $(WARNFLAGS) $(LIBS) -c ElectricChargeFluctuations.cpp -o ElectricChargeFluctuations.o

clean:
	rm -f $(OBJS)
 
tarz:
	tar zcf - $(MAKEFILE) $(SRCS) $(HDRS) > $(COMMAND).tar.gz
 
