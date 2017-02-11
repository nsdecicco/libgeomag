
# Author: Nicholas DeCicco <nsd.cicco@gmail.com>
#                          <deciccon0@students.rowan.edu>

CC = gcc
CFLAGS = -Wall -g
LIBS = -lnetcdf
AR = ar

STATIC_LIB = libgeomag.a

SHARED_LIB = libgeomag.so.1.0.0

EMBEDDED_LIB = libgeomag-embedded.a

OBJS = geomag.o

UNIT_TESTS = test_geomag print_model

all: $(UNIT_TESTS) $(STATIC_LIB) #$(SHARED_LIB)

$(OBJS): %.o : %.c
	$(CC) -fPIC $(CFLAGS) -c $< -o $@

$(STATIC_LIB): $(OBJS)
	$(AR) rcs $@ $<

geomag-embedded.o: geomag.c
	$(CC) -fPIC $(CFLAGS) -DTARGET_EMBEDDED -c $< -o $@

$(EMBEDDED_LIB): geomag-embedded.o igrf.o
	$(AR) rcs $@ $^

$(SHARED_LIB): $(OBJS)
	$(CC) -shared -Wl,-soname,libgeomag.so.1 -o $@ $<

$(UNIT_TESTS): % : %.c $(STATIC_LIB)
	$(CC) -L. $(CFLAGS) $< -Bstatic -lgeomag -Bdynamic $(LIBS) -o $@

igrf.c: print_model
	./print_model

igrf.o: igrf.c
	$(CC) $(CFLAGS) -c $< -o $@

test_geomag_embedded: test_geomag.c libgeomag-embedded.a
	$(CC) -L. $(CFLAGS) -DTEST_EMBEDDED $< -Bstatic -lgeomag-embedded -Bdynamic $(LIBS) -o $@

clean:
	-rm $(UNIT_TESTS) $(OBJS) $(SHARED_LIB) $(STATIC_LIB)
	-rm libgeomag-embedded.a igrf.o geomag-embedded.o test_geomag_embedded
