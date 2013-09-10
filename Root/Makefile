IDIR=
CC=g++

LIB1=.
CFLAGS=-c -g -Wall `root-config --cflags` -I${LIB1}
LDFLAGS=`root-config --glibs`
LDFLAGS2=${ROOTSYS}/lib
SOURCES=DATA.cpp DataExp.cpp read.cpp main.cpp
OBJECTS=$(SOURCES:.cpp=.o)
#EXECUTABLE=$(SOURCES:.cpp=.exe)
HEADERS=$(SOURCES:.cpp=.h)

exes = main.exe

#.cpp.o:
#$(CC) $(CFLAGS) $< -o $@

#all: $(SOURCES) $(EXECUTABLE)

$(exes): $(OBJECTS)
	$(CC) -o $@ $^ $(LDFLAGS) ${LDFLAGS2}/libMathMore.so 

%.o:%.cpp %.h
	 $(CC)  $(CFLAGS) -c $<


.PHONY: clean

clean:
	rm -f *.exe *.o



