IDIR=
CC=g++

LIB1=.
CFLAGS=-c -g -Wall `root-config --cflags` -I${LIB1}
LDFLAGS=`root-config --glibs`
LDFLAGS2=${ROOTSYS}/lib
SOURCES=dnde_loglog.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=$(SOURCES:.cpp=.exe)
HEADERS=$(SOURCES:.cpp=.h)

#.cpp.o:
#$(CC) $(CFLAGS) $< -o $@

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) ${LDFLAGS2}/libMathMore.so $(OBJECTS) -o $@

%.o:%.cpp %.h
	 $(CC)  $(CFLAGS) -c $<


.PHONY: clean

clean:
	rm -f *.exe *.o



