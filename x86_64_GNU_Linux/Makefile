SHELL = /bin/sh

CC = gcc
CXX = g++

CFLAGS = -static --std=c99 -O3 -g -march=native -W -Wall -pedantic
CXXFLAGS = -static-libstdc++ -O3 -g


LDFLAGS = 
LIBS = -lpthread -lm

SRC = ibs.c ibs_readdata.c ibs_trees.c ibs_markers.c \
ibs_obs.c ibs_kin.c ibs_founders.c

OBJ = ../main.o ../perm.o
PROGS = diffExprPermutation

all: $(PROGS)

diffExprPermutation: $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJ) $(LDFLAGS) $(LIBS)

clean: 
	rm -f *~ *.o *.a *.bak a.out core seedfile depend

distclean: clean
	rm -f diffExprPermutation
