
# Author: Nicholas DeCicco <nsd.cicco@gmail.com>
#                          <deciccon0@students.rowan.edu>

CC = gcc
CFLAGS = -Wall -g
LIBS = -lnetcdf
AR = ar

STATIC_LIB = libgeomag.a

SHARED_LIB = libgeomag.so.1.0.0

OBJS = geomag.o

all: $(STATIC_LIB) #$(SHARED_LIB)

$(OBJS): %.o : %.c
	$(CC) -fPIC $(CFLAGS) -c $< -o $@

$(STATIC_LIB): $(OBJS)
	$(AR) rcs $@ $<

$(SHARED_LIB): $(OBJS)
	$(CC) -shared -Wl,-soname,libgeomag.so.1 -o $@ $<

clean:
	-rm $(OBJS) $(SHARED_LIB) $(STATIC_LIB)
