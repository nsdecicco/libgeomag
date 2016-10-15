
# Author: Nicholas DeCicco <nsd.cicco@gmail.com>
#                          <deciccon0@students.rowan.edu>

CC = gcc
CFLAGS = -Wall -g
LIBS = -lnetcdf
AR = ar

STATIC_LIB = libgeomag.a

SHARED_LIB = libgeomag.so.1.0.0

OBJS = geomag.o

UNIT_TESTS = test_geomag

all: $(UNIT_TESTS) $(STATIC_LIB) #$(SHARED_LIB)

$(OBJS): %.o : %.c
	$(CC) -fPIC $(CFLAGS) -c $< -o $@

$(STATIC_LIB): $(OBJS)
	$(AR) rcs $@ $<

$(SHARED_LIB): $(OBJS)
	$(CC) -shared -Wl,-soname,libgeomag.so.1 -o $@ $<

$(UNIT_TESTS): % : %.c $(STATIC_LIB)
	$(CC) -L. $(CFLAGS) $< -Bstatic -lgeomag -Bdynamic $(LIBS) -o $@

clean:
	-rm $(UNIT_TESTS) $(OBJS) $(SHARED_LIB) $(STATIC_LIB)
